#include "StdAfx.hpp"

#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/MaxNuStencil.hpp"

class TurbulentSimulation: public Simulation {

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& turbulentFlowField);
  ~TurbulentSimulation() = default;

private:
  TurbulentFlowField turbulentFlowField_;
  Stencils::MaxNuStencil            maxNuStencil_;
  FieldIterator<TurbulentFlowField> maxNuFieldIterator_;
  void                              setTimeStep() override;
  Stencils::TurbulentViscosityStencil turbulentViscosityStencil_;
  FieldIterator<TurbulentFlowField>  turbulentViscosityIterator_;

  void initializeFlowField() override;

  void solveTimestep() override;
};