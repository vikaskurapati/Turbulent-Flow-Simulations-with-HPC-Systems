#include "StdAfx.hpp"
#include <algorithm>

#include "TurbulentFlowField.hpp"

#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny):
  // Pressure field doesn't need to have an extra layer, but this allows to address the same
  // positions with the same iterator for both pressures and velocities.
  FlowField(Nx, Ny),
  boundaryLayerThickness_(ScalarField(Nx + 3, Ny + 3)),
  wallDistance_(ScalarField(Nx + 3, Ny + 3)),
  turbulentViscosity_(ScalarField(Nx + 3, Ny + 3)),
  currentTurbulentViscosityTransport_(ScalarField(Nx + 3, Ny + 3)),
  previousTurbulentViscosityTransport_(ScalarField(Nx + 3, Ny + 3)) {
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny, int Nz):
  FlowField(Nx, Ny, Nz),
  boundaryLayerThickness_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  wallDistance_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  turbulentViscosity_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  currentTurbulentViscosityTransport_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  previousTurbulentViscosityTransport_(ScalarField(Nx + 3, Ny + 3, Nz + 3)) {
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
  ASSERTION(Nz > 0);
}

TurbulentFlowField::TurbulentFlowField(const Parameters& parameters):
  FlowField(parameters),
  boundaryLayerThickness_(
    parameters.geometry.dim == 2
      ? ScalarField(getNx() + 3, getNy() + 3)
      : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  wallDistance_(
    parameters.geometry.dim == 2
      ? ScalarField(getNx() + 3, getNy() + 3)
      : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  turbulentViscosity_(
    parameters.geometry.dim == 2
      ? ScalarField(getNx() + 3, getNy() + 3)
      : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  currentTurbulentViscosityTransport_(
    parameters.geometry.dim == 2
      ? ScalarField(getNx() + 3, getNy() + 3, 0.0)
      : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3, 0.0)
  ),
  previousTurbulentViscosityTransport_(
    parameters.geometry.dim == 2
      ? ScalarField(getNx() + 3, getNy() + 3, (1/parameters.flow.Re)*5)
      : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3, (1/parameters.flow.Re)*5)
  ),
  method_(parameters.simulation.type) {

    //BOUNDARY CONDITIONS FOR nu_tilda

    if(parameters.geometry.dim == 2){

      for (int i =0; i<getNx()+3; i++)
      {
        getCurrentTurbulentViscosityTransport().getScalar(i,0)=0.0;
        getCurrentTurbulentViscosityTransport().getScalar(i,getNy()+2)=0.0;
      }

      for (int j =0; j<getNy()+3; j++)
      {
        getCurrentTurbulentViscosityTransport().getScalar(0,j)=0.0;
        getCurrentTurbulentViscosityTransport().getScalar(getNx()+2,j)=0.0;
      }
    }
  }

ScalarField& TurbulentFlowField::getBoundaryLayerThickness() { return boundaryLayerThickness_; }

ScalarField& TurbulentFlowField::getWallDistance() { return wallDistance_; }

ScalarField& TurbulentFlowField::getTurbulentViscosity() { return turbulentViscosity_; }

ScalarField& TurbulentFlowField::getCurrentTurbulentViscosityTransport() { return currentTurbulentViscosityTransport_; }

ScalarField& TurbulentFlowField::getPreviousTurbulentViscosityTransport() { return previousTurbulentViscosityTransport_; }


void TurbulentFlowField::getViscosity(RealType& viscosity, int i, int j) {
  viscosity = getTurbulentViscosity().getScalar(i, j);
}

void TurbulentFlowField::getViscosity(RealType& viscosity, int i, int j, int k) {
  viscosity = getTurbulentViscosity().getScalar(i, j, k);
}

void TurbulentFlowField::getCurrentViscosityTransport(RealType& viscosity, int i, int j) {
  viscosity = getCurrentTurbulentViscosityTransport().getScalar(i, j);
}

void TurbulentFlowField::getCurrentViscosityTransport(RealType& viscosity, int i, int j, int k) {
  viscosity = getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
}

void TurbulentFlowField::getPreviousViscosityTransport(RealType& viscosity, int i, int j) {
  viscosity = getPreviousTurbulentViscosityTransport().getScalar(i, j);
}

void TurbulentFlowField::getPreviousViscosityTransport(RealType& viscosity, int i, int j, int k) {
  viscosity = getPreviousTurbulentViscosityTransport().getScalar(i, j, k);
}

void TurbulentFlowField::setPreviousViscosityTransport(ScalarField& current_visc_field ){
  //make the previous field equal to the current field
  previousTurbulentViscosityTransport_ = current_visc_field;
}

void TurbulentFlowField::getH(RealType& h, int i, int j) { h = getWallDistance().getScalar(i, j); }

void TurbulentFlowField::getH(RealType& h, int i, int j, int k) { h = getWallDistance().getScalar(i, j, k); }

void TurbulentFlowField::getDelta(RealType& delta, int i, int j) {
  delta = getBoundaryLayerThickness().getScalar(i, j);
}

void TurbulentFlowField::getDelta(RealType& delta, int i, int j, int k) {
  delta = getBoundaryLayerThickness().getScalar(i, j, k);
}