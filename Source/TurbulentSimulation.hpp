#include "StdAfx.hpp"

#include "Simulation.hpp"

class TurbulentSimulation: public Simulation{
    TurbulentSimulation(Parameters& parameters, FlowField& flowField);
    ~TurbulentSimulation();

    void setTimeStep() override;

    void initializeFlowField() override;
};