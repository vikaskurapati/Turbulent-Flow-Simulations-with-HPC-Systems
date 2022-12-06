#include "StdAfx.hpp"

#include "InitBoundaryLayerThickness.hpp"

Stencils::InitBoundaryLayerThickness::InitBoundaryLayerThickness(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::InitBoundaryLayerThickness::apply(FlowField& flowField, int i, int j) {
  RealType x    = parameters_.meshsize->getDx(i, j) * (parameters_.parallel.firstCorner[0] + i);
  RealType Re_x = parameters_.walls.vectorLeft[0] * x / parameters_.geometry.lengthX;
  if (parameters_.turbulence.boundaryLayerType == "inviscid") {
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.0;
  } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 4.91 * (x) / pow(Re_x, 0.5);
    //std::cout<<flowField.getBoundaryLayerThickness().getScalar(i, j)<<std::endl;
    //exit(0);
  } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
    flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.382 * (x) / pow(Re_x, 0.2);
  }
}

void Stencils::InitBoundaryLayerThickness::apply(FlowField& flowField, int i, int j, int k) {
  RealType x    = parameters_.meshsize->getDx(i, j,k) * (parameters_.parallel.firstCorner[0] + i);
  RealType Re_x = parameters_.walls.vectorLeft[0] * x / parameters_.geometry.lengthX;
  if (parameters_.turbulence.boundaryLayerType == "inviscid") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.0;
  } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 4.91 * (x) / pow(Re_x, 0.5);
  } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.382 * (x) / pow(Re_x, 0.2);
  }
}
