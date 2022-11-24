#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {
    
    leftPressureBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightPressureBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomPressureBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topPressureBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontPressureBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backPressureBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  } 
  
  else {
    leftPressureBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightPressureBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomPressureBuffer = std::make_unique<RealType[]>(localSize[0]);
    topPressureBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
    if (j >= 2) {
    *(leftPressureBuffer.get() + (j - 2)) = (flowField.getPressure().getScalar(i + 1, j));
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
    if (j >= 2) {
    *(rightPressureBuffer.get() + (j - 2)) = (flowField.getPressure().getScalar(i, j));
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
      if ((i >= 2)) {
    *(bottomPressureBuffer.get() + (i - 2) * localSize[2]) = (flowField.getPressure().getScalar(i, j + 1));
  }
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
      if ((i >= 2)) {
    *(topPressureBuffer.get() + (i - 2) * localSize[2]) = (flowField.getPressure().getScalar(i, j));
  }
}

void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2)) {
    *(leftPressureBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getPressure().getScalar(i + 1, j, k));
  }
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2)) {
    *(rightPressureBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getPressure().getScalar(i, j, k));
  }
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2)) {
    *(bottomPressureBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getPressure().getScalar(i, j + 1, k));
  }
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2)) {
    *(topPressureBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getPressure().getScalar(i, j, k));
  }
}

void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2)) {
    *(frontPressureBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getPressure().getScalar(i, j, k + 1));
  }
}

void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2)) {
    *(backPressureBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getPressure().getScalar(i, j, k));
  }
}