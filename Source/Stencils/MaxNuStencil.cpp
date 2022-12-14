#include "StdAfx.hpp"

#include "MaxNuStencil.hpp"

#include "TurbulentFlowField.hpp"

Stencils::MaxNuStencil::MaxNuStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {
  reset_Nu();
}

void Stencils::MaxNuStencil::apply(TurbulentFlowField& flowField, int i, int j) { cellMaxNuValue(flowField, i, j); }

void Stencils::MaxNuStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  cellMaxNuValue(flowField, i, j, k);
}

// Finding max nu in 2D
void Stencils::MaxNuStencil::cellMaxNuValue(TurbulentFlowField& flowField, int i, int j) {
  RealType localViscosity = flowField.getTurbulentViscosity().getScalar(i, j);
  if (localViscosity > maxNuValue_) {
    maxNuValue_ = localViscosity;
  }
}
// Finding max nu in 3D
void Stencils::MaxNuStencil::cellMaxNuValue(TurbulentFlowField& flowField, int i, int j, int k) {
  RealType localViscosity = flowField.getTurbulentViscosity().getScalar(i, j, k);

  if (localViscosity > maxNuValue_) {
    maxNuValue_ = localViscosity;
  }
}

void Stencils::MaxNuStencil::reset_Nu() { maxNuValue_ = 0; }

RealType Stencils::MaxNuStencil::getMaxNuValues() const { return maxNuValue_; }
