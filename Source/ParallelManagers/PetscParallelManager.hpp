#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"

#include "../FlowField.hpp"
#include "../Iterators.hpp"
#include "../Stencils/PressureBufferFillStencil.hpp"
#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferFillStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

namespace ParallelManagers {
  class PetscParallelManager {

  private:
    Parameters& parameters_;
    FlowField&  flowfield_;

    Stencils::VelocityBufferFillStencil fillVelocityStencil;
    Stencils::VelocityBufferReadStencil readVelocityStencil;

    Stencils::PressureBufferFillStencil fillPressureStencil;
    Stencils::PressureBufferReadStencil readPressureStencil;

    ParallelBoundaryIterator<FlowField> velocityfillIterator;
    ParallelBoundaryIterator<FlowField> velocityreadIterator;
    ParallelBoundaryIterator<FlowField> pressurereadIterator;
    ParallelBoundaryIterator<FlowField> pressurefillIterator;

  public:

    PetscParallelManager(Parameters& parameters, FlowField& flowfield);
    void communicateVelocity();
    void communicatePressure();
    
    ~PetscParallelManager() {}
  };
} // namespace ParallelManagers