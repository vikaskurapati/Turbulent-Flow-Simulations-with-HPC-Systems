#include "StdAfx.hpp"

#include "TransportViscosityBufferReadStencil.hpp"

#include "Definitions.hpp"

Stencils::TransportViscosityBufferReadStencil::TransportViscosityBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftTransportViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightTransportViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomTransportViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topTransportViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontTransportViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backTransportViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftTransportViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightTransportViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomTransportViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0]);
    topTransportViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
} // End of constructor

// For 2D Cases
void Stencils::TransportViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.leftNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 1, j) = *(leftTransportViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.rightNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) = *(rightTransportViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j + 1) = *(bottomTransportViscosityReadBuffer.get() + (i - 2));
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) = *(topTransportViscosityReadBuffer.get() + (i - 2));
    }
  }
}

// For 3D Cases
void Stencils::TransportViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.leftNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(
        i + 1, j, k
      ) = *(leftTransportViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.rightNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(i, j, k) = *(rightTransportViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(
        i, j + 1, k
      ) = *(bottomTransportViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(i, j, k) = *(topTransportViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.frontNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(
        i, j, k + 1
      ) = *(frontTransportViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}

void Stencils::TransportViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.backNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getCurrentTurbulentViscosityTransport().getScalar(i, j, k) = *(backTransportViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}
