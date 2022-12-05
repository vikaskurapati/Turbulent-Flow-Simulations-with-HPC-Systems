#include "StdAfx.hpp"

#include "InitBoundaryLayerThickness.hpp"

Stencils::InitBoundaryLayerThickness::InitBoundaryLayerThickness(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::InitBoundaryLayerThickness::apply(FlowField& flowField, int i, int j) {
  if (parameters_.turbulence.boundaryLayerType == "inviscid") {
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.0;
  } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
    // Need to ADD Reynolds number at different x positions
    // =========================================
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 4.91 * (parameters_.meshsize->getDx(i, j) * (parameters_.parallel.firstCorner[0] + i))
                                              / pow((parameters_.flow.Re), 0.5);
  } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.382 * (parameters_.meshsize->getDx(i, j) * (parameters_.parallel.firstCorner[0] + i))
                                              / pow((parameters_.flow.Re), 0.2);
    // =========================================
  }
}

void Stencils::InitBoundaryLayerThickness::apply(FlowField& flowField, int i, int j, int k) {
  if (parameters_.turbulence.boundaryLayerType == "inviscid") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.0;
  } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
    // Need to ADD Reynolds number at different x positions
    // =========================================
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 4.91 * (parameters_.meshsize->getDx(i, j, k) * (parameters_.parallel.firstCorner[0] + i))
                                              / pow((parameters_.flow.Re), 0.5);
  } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.382 * (parameters_.meshsize->getDx(i, j, k) * (parameters_.parallel.firstCorner[0] + i))
                                              / pow((parameters_.flow.Re), 0.2);
    // =========================================
  }
}
