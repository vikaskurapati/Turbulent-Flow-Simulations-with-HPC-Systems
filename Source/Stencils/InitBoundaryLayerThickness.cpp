#include "StdAfx.hpp"

#include "InitBoundaryLayerThickness.hpp"

Stencils::InitBoundaryLayerThickness::InitBoundaryLayerThickness(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::InitBoundaryLayerThickness::apply(TurbulentFlowField& flowField, int i, int j) {

  if (i > 1 && j > 1 && i <= (parameters_.geometry.sizeX + 1) && j <= (parameters_.geometry.sizeY + 1)) {
    RealType x    = parameters_.meshsize->getPosX(i, j) + (0.5 * parameters_.meshsize->getDx(i, j));
    RealType Re_x = parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re; 

    if (parameters_.turbulence.boundaryLayerType == "inviscid") {
      flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.0;
    } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
      flowField.getBoundaryLayerThickness().getScalar(i, j) = 4.91 * (x) / pow(Re_x, 0.5);
    } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
      flowField.getBoundaryLayerThickness().getScalar(i, j) = 0.382 * (x) / pow(Re_x, 0.2);
    }
  }
}

void Stencils::InitBoundaryLayerThickness::apply(TurbulentFlowField& flowField, int i, int j, int k) {

  if (i > 1 && j>1 && k>1 && i<=(parameters_.geometry.sizeX+1) && j<=(parameters_.geometry.sizeY+1) && k<=(parameters_.geometry.sizeZ+1) ){
  RealType x    = parameters_.meshsize->getPosX(i, j, k) + (0.5 * parameters_.meshsize->getDx(i, j, k));
  RealType Re_x = parameters_.walls.vectorLeft[0] * x * parameters_.flow.Re;
  if (parameters_.turbulence.boundaryLayerType == "inviscid") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.0;
  } else if (parameters_.turbulence.boundaryLayerType == "laminar") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 4.91 * (x) / pow(Re_x, 0.5);
  } else if (parameters_.turbulence.boundaryLayerType == "turbulence") {
    flowField.getBoundaryLayerThickness().getScalar(i, j, k) = 0.382 * (x) / pow(Re_x, 0.2);
  }
  }
}
