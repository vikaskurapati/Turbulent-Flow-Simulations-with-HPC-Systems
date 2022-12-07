#include "StdAfx.hpp"
#include "Stencils/MaxNuStencil.hpp"
#include "Simulation.hpp"

class TurbulentSimulation: public Simulation{

    Stencils::MaxNuStencil maxNuStencil_;
    
    TurbulentSimulation(Parameters& parameters, FlowField& flowField);
    ~TurbulentSimulation();

    void setTimeStep() override;

    void initializeFlowField() override;
};