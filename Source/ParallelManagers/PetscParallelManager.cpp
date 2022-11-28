#include "PetscParallelManager.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters &parameters, FlowField &flowfield):
  parameters_(parameters),
  flowfield_(flowfield),
  fillVelocityStencil(parameters),
  readVelocityStencil(parameters),
  velocityfillIterator(flowfield, parameters, fillVelocityStencil, 1, -1),
  velocityreadIterator(flowfield, parameters, readVelocityStencil, 1, -1),
  fillPressureStencil(parameters),
  readPressureStencil(parameters),
  pressurefillIterator(flowfield, parameters, fillPressureStencil, 1, -1),
  pressurereadIterator(flowfield, parameters, readPressureStencil, 1, -1) {}

