#include "StdAfx.hpp"

#include "Simulation.hpp"

#include "Stencils/MaxNuStencil.hpp"

class TurbulentSimulation: public Simulation {

public:
  TurbulentSimulation(Parameters& parameters, FlowField& flowField);
  ~TurbulentSimulation() = default;

private:
  Stencils::MaxNuStencil maxNuStencil_;
  FieldIterator<FlowField>          maxNuFieldIterator_;
  void                   setTimeStep() override;

  void initializeFlowField() override;
};