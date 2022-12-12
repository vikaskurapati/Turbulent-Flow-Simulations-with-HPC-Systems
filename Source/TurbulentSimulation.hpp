#include "StdAfx.hpp"

#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/MaxNuStencil.hpp"

class TurbulentSimulation: public Simulation {

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& turbulentFlowField);
  virtual ~TurbulentSimulation() = default;

private:
  TurbulentFlowField                  turbulentFlowField_;
  Stencils::MaxNuStencil              maxNuStencil_;
  FieldIterator<TurbulentFlowField>   maxNuFieldIterator_;
  Stencils::TurbulentViscosityStencil turbulentViscosityStencil_;
  FieldIterator<TurbulentFlowField>   turbulentViscosityIterator_;

  void initializeFlowField() override;

  Stencils::TurbulentFGHStencil     turbulentfghStencil_;
  FieldIterator<TurbulentFlowField> turbulentfghIterator_;

  void solveTimestep() override;

  void plotVTK(int timeStep, RealType simulationTime) override;
  void setTimeStep() override;
  ParallelManagers::PetscTurbulentParallelManager parallel_manager_;
};