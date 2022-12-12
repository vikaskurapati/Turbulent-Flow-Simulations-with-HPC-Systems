#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"

#include "Solvers/LinearSolver.hpp"
#include "Stencils/BFInputStencils.hpp"
#include "Stencils/BFStepInitStencil.hpp"
#include "Stencils/FGHStencil.hpp"
#include "Stencils/RHSStencil.hpp"
#include "Stencils/InitTaylorGreenFlowFieldStencil.hpp"
#include "Stencils/InitWallDistanceStencil.hpp"
#include "Stencils/InitBoundaryLayerThickness.hpp"
#include "Stencils/MaxUStencil.hpp"
#include "Stencils/MovingWallStencils.hpp"
#include "Stencils/NeumannBoundaryStencils.hpp"
#include "Stencils/ObstacleStencil.hpp"
#include "Stencils/PeriodicBoundaryStencils.hpp"
#include "Stencils/VelocityStencil.hpp"
#include "Stencils/VTKStencil.hpp"
#include "./ParallelManagers/PetscParallelManager.hpp"
#include "Stencils/TurbulentViscosityStencil.hpp"

class Simulation {
protected:
  Parameters& parameters_;
  FlowField& flowField_;

  Stencils::MaxUStencil             maxUStencil_;
  FieldIterator<FlowField>          maxUFieldIterator_;
  GlobalBoundaryIterator<FlowField> maxUBoundaryIterator_;

  // Set up the boundary conditions
  GlobalBoundaryFactory             globalBoundaryFactory_;
  GlobalBoundaryIterator<FlowField> wallVelocityIterator_;
  GlobalBoundaryIterator<FlowField> wallFGHIterator_;
  
  Stencils::FGHStencil     fghStencil_;
  FieldIterator<FlowField> fghIterator_;
  FieldIterator<FlowField> rhsIterator_;
  Stencils::RHSStencil     rhsStencil_;

  Stencils::VelocityStencil velocityStencil_;
  Stencils::ObstacleStencil obstacleStencil_;
  FieldIterator<FlowField>  velocityIterator_;
  FieldIterator<FlowField>  obstacleIterator_;

  std::unique_ptr<Solvers::LinearSolver> solver_;

  ParallelManagers::PetscParallelManager parallel_manager_;

  virtual void setTimeStep();

  //ParallelManagers::PetscParallelManager comm_;
  //ParallelBoundaryIterator<FlowField> velocityParallelBoundaryIterator_;
  // ParallelBoundaryIterator<FlowField> fghParallelBoundaryIterator_;

  int rank_;

public:
  Simulation(Parameters& parameters, FlowField& flowField);
  virtual ~Simulation() = default;

  /** Initialises the flow field according to the scenario */
  virtual void initializeFlowField();

  virtual void solveTimestep();

  /** Plots the flow field */
  virtual void plotVTK(int timeStep, RealType simulationTime);
};
