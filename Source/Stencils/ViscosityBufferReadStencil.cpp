#include "StdAfx.hpp"

#include "ViscosityBufferReadStencil.hpp"

#include "Definitions.hpp"

Stencils::ViscosityBufferReadStencil::ViscosityBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<TurbulentFlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftViscosityReadBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightViscosityReadBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomViscosityReadBuffer = std::make_unique<RealType[]>(localSize[0]);
    topViscosityReadBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
} // End of constructor

// For 2D Cases
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.leftNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getTurbulentViscosity().getScalar(i + 1, j) = *(leftViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.rightNb >= 0) {
    if (j >= 2 && j <= (localSize[1] + 1)) {
      flowField.getTurbulentViscosity().getScalar(i, j) = *(rightViscosityReadBuffer.get() + (j - 2));
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getTurbulentViscosity().getScalar(i, j + 1) = *(bottomViscosityReadBuffer.get() + (i - 2));
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2 && i <= localSize[0] + 1)) {
      flowField.getTurbulentViscosity().getScalar(i, j) = *(topViscosityReadBuffer.get() + (i - 2));
    }
  }
}

// For 3D Cases
void Stencils::ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.leftNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(
        i + 1, j, k
      ) = *(leftViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.rightNb >= 0) {
    if ((j >= 2) && (k >= 2) && j <= (localSize[1] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(i, j, k) = *(rightViscosityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]);
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.bottomNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(
        i, j + 1, k
      ) = *(bottomViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.topNb >= 0) {
    if ((i >= 2) && (k >= 2) && i <= (localSize[0] + 1) && k <= (localSize[2] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(i, j, k) = *(topViscosityReadBuffer.get() + (k - 2) + (i - 2) * localSize[2]);
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.frontNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(
        i, j, k + 1
      ) = *(frontViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}

void Stencils::ViscosityBufferReadStencil::applyBackWall(TurbulentFlowField& flowfield, int i, int j, int k) {
  if (parameters_.parallel.backNb >= 0) {
    if ((i >= 2) & (j >= 2) && i <= (localSize[0] + 1) && j <= (localSize[1] + 1)) {
      flowfield.getTurbulentViscosity().getScalar(i, j, k) = *(backViscosityReadBuffer.get() + (j - 2) + (i - 2) * localSize[1]);
    }
  }
}
