#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftPressureFillBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightPressureFillBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomPressureFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topPressureFillBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontPressureFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backPressureFillBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftPressureFillBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightPressureFillBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomPressureFillBuffer = std::make_unique<RealType[]>(localSize[0]);
    topPressureFillBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
}
// For 2D Cases
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    *(leftPressureFillBuffer.get() + (j - 2)) = (flowField.getPressure().getScalar(i + 2, j));
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    // Need to verify indices
    *(rightPressureFillBuffer.get() + (j - 2)) = (flowField.getPressure().getScalar(i - 1, j));
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(bottomPressureFillBuffer.get() + (i - 2)) = (flowField.getPressure().getScalar(i, j + 2));
  }
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(topPressureFillBuffer.get() + (i - 2)) = (flowField.getPressure().getScalar(i, j - 1));
  }
}

// For 3D Cases
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(leftPressureFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getPressure().getScalar(i + 2, j, k)
    );
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(rightPressureFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]
    ) = (flowField.getPressure().getScalar(i - 1, j, k));
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(bottomPressureFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)
    ) = (flowField.getPressure().getScalar(i, j + 2, k));
  }
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(topPressureFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getPressure().getScalar(i, j - 1, k)
    );
  }
}

void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(frontPressureFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)
    ) = (flowField.getPressure().getScalar(i, j, k + 2));
  }
}

void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(backPressureFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getPressure().getScalar(i, j, k - 1)
    );
  }
}