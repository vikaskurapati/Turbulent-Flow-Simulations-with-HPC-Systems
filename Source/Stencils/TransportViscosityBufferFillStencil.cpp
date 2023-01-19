#include "StdAfx.hpp"

#include "TransportViscosityBufferFillStencil.hpp"

#include "Definitions.hpp"

Stencils::TransportViscosityBufferFillStencil::TransportViscosityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftTransportViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightTransportViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomTransportViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topTransportViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontTransportViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backTransportViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftTransportViscosityFillBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightTransportViscosityFillBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomTransportViscosityFillBuffer = std::make_unique<RealType[]>(localSize[0]);
    topTransportViscosityFillBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
}
// For 2D Cases
void Stencils::TransportViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    *(leftTransportViscosityFillBuffer.get() + (j - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 2, j));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (j >= 2 && j <= (localSize[1] + 1)) {
    // Need to verify indices
    *(rightTransportViscosityFillBuffer.get() + (j - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(bottomTransportViscosityFillBuffer.get() + (i - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j + 2));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if ((i >= 2 && i <= localSize[0] + 1)) {
    *(topTransportViscosityFillBuffer.get() + (i - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j - 1));
  }
}

// For 3D Cases
void Stencils::TransportViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(leftTransportViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 2, j, k)
    );
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
    *(rightTransportViscosityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]
    ) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j, k));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(bottomTransportViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)
    ) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j + 2, k));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
    *(topTransportViscosityFillBuffer.get() + (i - 2) * localSize[2] + (k - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j - 1, k)
    );
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(frontTransportViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)
    ) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k + 2));
  }
}

void Stencils::TransportViscosityBufferFillStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) {
  if ((i >= 2) && (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
    *(backTransportViscosityFillBuffer.get() + (i - 2) * localSize[1] + (j - 2)) = (flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k - 1)
    );
  }
}