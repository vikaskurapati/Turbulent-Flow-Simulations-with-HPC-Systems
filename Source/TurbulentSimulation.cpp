#include "TurbulentSimulation.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, FlowField& flowField):
    Simulation(Parameters& parameters, FlowField& flowField) {}

void TurbulentSimulation::initializeFlowField() {

  if (parameters_.simulation.scenario == "taylor-green") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitTaylorGreenFlowFieldStencil stencil(parameters_);
    FieldIterator<FlowField>                  iterator(flowField_, parameters_, stencil);
    iterator.iterate();
  } else if (parameters_.simulation.scenario == "channel") {
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
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
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
  }

  // Adding Nearest Wall distance for turbulence
    // This initialisation is only required to calculate Nearest Wall distance for turbulence.
    Stencils::InitWallDistanceStencil wallDistancestencil(parameters_);
    FieldIterator<FlowField>          wallDistanceiterator(flowField_, parameters_, wallDistancestencil);
    wallDistanceiterator.iterate();

    Stencils::InitBoundaryLayerThickness boundaryLayerThicknessstencil(parameters_);
    FieldIterator<FlowField>          boundaryLayerThicknesssiterator(flowField_, parameters_, boundaryLayerThicknessstencil);
    boundaryLayerThicknesssiterator.iterate();

  solver_->reInitMatrix();
}

void Simulation::setTimeStep() {

  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());
  // Determine maximum velocity
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2] + EPSILON);
  } else {
    parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0] + EPSILON);
  }

  // localMin = std::min(parameters_.timestep.dt, std::min(std::min(parameters_.flow.Re/(2 * factor), 1.0 /
  // maxUStencil_.getMaxValues()[0]), 1.0 / maxUStencil_.getMaxValues()[1]));
  localMin = std::min(
    parameters_.flow.Re / (2 * factor),
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