#include "StdAfx.hpp"

#include "InitWallDistanceStencil.hpp"

Stencils::InitWallDistanceStencil::InitWallDistanceStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xObsCells(parameters.bfStep.xRatio * parameters.geometry.sizeX),
  yObsCells(parameters.bfStep.yRatio * parameters.geometry.sizeY) {}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j) {

  const int obstacle = flowField.getFlags().getValue(i, j);
  
  RealType distance = 0;

  RealType temp = 0;

  if ((obstacle & OBSTACLE_SELF) == 0) {
    distance = (parameters_.geometry.sizeY - (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)))
                        * parameters_.meshsize->getDy(i, j);
    temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) * parameters_.meshsize->getDy(i, j);
    if (temp < distance) {
      distance = temp;
    }

    //For BFS case:
    if (xObsCells > 0 && yObsCells > 0) { 
      // For cells above the step
      if ((parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) > yObsCells && (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) < xObsCells){              
        temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5) - yObsCells) * parameters_.meshsize->getDy(i, j);
        if (temp < distance) {
          distance = temp;
        }
      }

      // For cells on the rightside of the step
      if ((parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) < yObsCells && (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) > xObsCells){
        temp = (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5) - xObsCells) * parameters_.meshsize->getDx(i, j);
        if (temp < distance) {
          distance = temp;
        }
      }
    }

    flowField.getWallDistance().getScalar(i, j) = distance;
  }
}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {

  const int obstacle = flowField.getFlags().getValue(i, j, k);

  RealType distance = 0;

  RealType temp = 0;

  if ((obstacle & OBSTACLE_SELF) == 0) {
    distance = (parameters_.geometry.sizeY - (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)))
                        * parameters_.meshsize->getDy(i, j, k);
    temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) * parameters_.meshsize->getDy(i, j, k);
    if (temp < distance) {
      distance = temp;
    }
    temp = (parameters_.parallel.firstCorner[2] + (k - 2 + 0.5)) * parameters_.meshsize->getDz(i, j, k);
    if (temp < distance) {
      distance = temp;
    }
    temp = (parameters_.geometry.sizeZ - (parameters_.parallel.firstCorner[2] + (k - 2 + 0.5)))
           * parameters_.meshsize->getDz(i, j, k);
    if (temp < distance) {
      distance = temp;
    }

    //For BFS case:
    if (xObsCells > 0 && yObsCells > 0) {
      // For cells above the step
      if ((parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) > yObsCells && (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) < xObsCells){              
        temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5) - yObsCells) * parameters_.meshsize->getDy(i, j, k);
        if (temp < distance) {
          distance = temp;
        }
      }

      // For cells on the rightside of the step
      if ((parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) < yObsCells && (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) > xObsCells){
        temp = (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5) - xObsCells) * parameters_.meshsize->getDx(i, j, k);
        if (temp < distance) {
          distance = temp;
        }
      }
    }

    flowField.getWallDistance().getScalar(i, j, k) = distance;
  }
}
