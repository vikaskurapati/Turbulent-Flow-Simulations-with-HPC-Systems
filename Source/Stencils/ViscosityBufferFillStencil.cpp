#include "StdAfx.hpp"

#include "ViscosityBufferFillStencil.hpp"

#include "Definitions.hpp"

Stencils::ViscosityBufferFillStencil::ViscosityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0]);
    topViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
}
// For 2D Cases
void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    *(leftViscosityFillBuffer.get() + (j - 2)) = (flowField.getTurbulentViscosity().getScalar(i + 2, j));
  }
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    // Need to verify indices
    *(rightViscosityFillBuffer.get() + (j - 2)) = (flowField.getTurbulentViscosity().getScalar(i - 1, j));
  }
}

void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(bottomViscosityFillBuffer.get() + (i - 2)) = (flowField.getTurbulentViscosity().getScalar(i, j + 2));
  }
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(topViscosityFillBuffer.get() + (i - 2)) = (flowField.getTurbulentViscosity().getScalar(i, j - 1));
  }
}

// For 3D Cases
void Stencils::ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(leftViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getTurbulentViscosity().getScalar(i + 2, j, k)
    );
  }
}

void Stencils::ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(rightViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]
    ) = (flowField.getTurbulentViscosity().getScalar(i - 1, j, k));
  }
}

void Stencils::ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(bottomViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)
    ) = (flowField.getTurbulentViscosity().getScalar(i, j + 2, k));
  }
}

void Stencils::ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(topViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getTurbulentViscosity().getScalar(i, j - 1, k)
    );
  }
}

void Stencils::ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(frontViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)
    ) = (flowField.getTurbulentViscosity().getScalar(i, j, k + 2));
  }
}

void Stencils::ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(backViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getTurbulentViscosity().getScalar(i, j, k - 1)
    );
  }
}