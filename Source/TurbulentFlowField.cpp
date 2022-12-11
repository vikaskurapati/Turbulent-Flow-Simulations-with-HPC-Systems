#include "StdAfx.hpp"

#include "FlowField.hpp"
#include "TurbulentFlowField.hpp"

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny):
   // Pressure field doesn't need to have an extra layer, but this allows to address the same
  // positions with the same iterator for both pressures and velocities.
  FlowField(Nx, Ny), boundaryLayerThickness_(ScalarField(Nx + 3, Ny + 3)),
  wallDistance_(ScalarField(Nx + 3, Ny + 3)),
  turbulentViscosity_(ScalarField(Nx + 3, Ny + 3))
{
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny, int Nz):
  FlowField(Nx, Ny, Nz), boundaryLayerThickness_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  wallDistance_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  turbulentViscosity_(ScalarField(Nx + 3, Ny + 3, Nz + 3))
{
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
  ASSERTION(Nz > 0);
}

TurbulentFlowField::TurbulentFlowField(const Parameters& parameters):
FlowField(parameters),
  boundaryLayerThickness_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  wallDistance_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  turbulentViscosity_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ) {}

ScalarField& TurbulentFlowField::getBoundaryLayerThickness() { return boundaryLayerThickness_; }

ScalarField& TurbulentFlowField::getWallDistance() { return wallDistance_; }

ScalarField& TurbulentFlowField::getTurbulentViscosity() { return turbulentViscosity_; }

void TurbulentFlowField::getViscosity(RealType& viscosity, int i, int j) {
  viscosity = getTurbulentViscosity().getScalar(i, j);
}

void TurbulentFlowField::getViscosity(RealType& viscosity, int i, int j, int k) {
  viscosity = getTurbulentViscosity().getScalar(i, j, k);
}

void TurbulentFlowField::getH(RealType& h, int i, int j) {
  h = getWallDistance().getScalar(i, j);
}

void TurbulentFlowField::getH(RealType& h, int i, int j, int k) {
  h = getWallDistance().getScalar(i, j, k);
}

void TurbulentFlowField::getDelta(RealType& delta, int i, int j) {
  delta = getBoundaryLayerThickness().getScalar(i, j);
}

void TurbulentFlowField::getDelta(RealType& delta, int i, int j, int k) {
  delta  = getBoundaryLayerThickness().getScalar(i, j, k);
}