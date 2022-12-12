#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

#include "FlowField.hpp"
#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/MaxNuStencil.hpp"
#include "Stencils/TurbulentVTKStencil.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, TurbulentFlowField& flowField):
  Simulation(parameters, flowField),
  turbulentFlowField_(flowField),
  maxNuStencil_(parameters),
  maxNuFieldIterator_(flowField, parameters, maxNuStencil_),
  turbulentViscosityStencil_(parameters),
  turbulentViscosityIterator_(flowField, parameters, turbulentViscosityStencil_),
  turbulentfghStencil_(parameters),
  turbulentfghIterator_(turbulentFlowField_, parameters, turbulentfghStencil_),
  parallel_manager_(parameters, flowField) {}

void TurbulentSimulation::initializeFlowField() {

  if (parameters_.simulation.scenario == "taylor-green") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitTaylorGreenFlowFieldStencil stencil(parameters_);
    FieldIterator<FlowField>                  iterator(turbulentFlowField_, parameters_, stencil);
    iterator.iterate();
  } else if (parameters_.simulation.scenario == "channel") {
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(turbulentFlowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
    wallVelocityIterator_.iterate();
  } else if (parameters_.simulation.scenario == "pressure-channel") {
    // Set pressure boundaries here for left wall
    const RealType value = parameters_.walls.scalarLeft;
    ScalarField&   rhs   = flowField_.getRHS();

    if (parameters_.geometry.dim == 2) {
      const int sizey = flowField_.getNy();
      for (int i = 0; i < sizey + 3; i++) {
        rhs.getScalar(0, i) = value;
      }
    } else {
      const int sizey = flowField_.getNy();
      const int sizez = flowField_.getNz();
      for (int i = 0; i < sizey + 3; i++) {
        for (int j = 0; j < sizez + 3; j++) {
          rhs.getScalar(0, i, j) = value;
        }
      }
    }

    // Do same procedure for domain flagging as for regular channel
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(turbulentFlowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
  }
  // Adding Nearest Wall distance for turbulence
  // This initialisation is only required to calculate Nearest Wall distance for turbulence.
  Stencils::InitWallDistanceStencil wallDistancestencil(parameters_);
  FieldIterator<TurbulentFlowField> wallDistanceiterator(turbulentFlowField_, parameters_, wallDistancestencil);
  wallDistanceiterator.iterate();

  Stencils::InitBoundaryLayerThickness boundaryLayerThicknessstencil(parameters_);
  FieldIterator<TurbulentFlowField>    boundaryLayerThicknesssiterator(
    turbulentFlowField_, parameters_, boundaryLayerThicknessstencil
  );
  boundaryLayerThicknesssiterator.iterate();

  solver_->reInitMatrix();
}

void TurbulentSimulation::setTimeStep() {

  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());
  // Determine maximum velocity
  maxUStencil_.reset();
  maxNuStencil_.reset_Nu();
  maxNuFieldIterator_.iterate();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + EPSILON);
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + EPSILON);
  }
  // std::cout<<"MAX Nu value is: "<<maxNuStencil_.getMaxNuValues()<<std::endl;
  //  localMin = std::min(parameters_.timestep.dt, std::min(std::min(parameters_.flow.Re/(2 * factor), 1.0 /
  //  maxUStencil_.getMaxValues()[0]), 1.0 / maxUStencil_.getMaxValues()[1]));
  localMin = std::min(
    1 / ((1 / parameters_.flow.Re) + maxNuStencil_.getMaxNuValues()) / (2 * factor),
    std::min(
      parameters_.timestep.dt,
      std::min(1 / (maxUStencil_.getMaxValues()[0] + EPSILON), 1 / (maxUStencil_.getMaxValues()[1] + EPSILON))
    )
  );

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}

void TurbulentSimulation::solveTimestep() {
  turbulentViscosityIterator_.iterate();
#ifndef NDEBUG

  feclearexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  if (fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT))
    raise(SIGFPE);
#endif

  // Determine and set max. timestep which is allowed in this simulation
  parallel_manager_.communicateViscosity();
  setTimeStep();
  // Compute FGH
  turbulentfghIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // TODO WS1: compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  parallel_manager_.communicatePressure();
  // TODO WS2: communicate pressure values
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  parallel_manager_.communicateVelocity();
  // TODO WS2: communicate velocity values
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::TurbulentVTKStencil     vtkStencil(parameters_);
  FieldIterator<TurbulentFlowField> vtkIterator(turbulentFlowField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(turbulentFlowField_, timeStep, simulationTime);
}
