#include "StdAfx.hpp"

#include "InitWallDistanceStencil.hpp"

Stencils::InitWallDistanceStencil::InitWallDistanceStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xObsCells(parameters.bfStep.xRatio * parameters.geometry.sizeX),
  yObsCells(parameters.bfStep.yRatio * parameters.geometry.sizeY) {}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j) {}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {

  const int obstacle = flowField.getFlags().getValue(i, j, k);

  if ((obstacle & OBSTACLE_SELF) == 0) {
    RealType distance = (parameters_.geometry.sizeY - (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)))
                        * parameters_.meshsize->getDy(i, j, k);
    RealType temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) * parameters_.meshsize->getDy(i, j, k);
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

    if (parameters_.simulation.scenario != "channel") {
      temp = (parameters_.geometry.sizeX - (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)))
             * parameters_.meshsize->getDx(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
      temp = (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) * parameters_.meshsize->getDx(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
    }

    flowField.getWallDistance().getScalar(i, j, k) = distance;
  }

  const int* firstcorner = parameters_.parallel.firstCorner;
  // For BFS case, considering two cases.
  if (flowField.getFlags().getValue(firstcorner[0] + 2, firstcorner[1] + 2, firstcorner[2] + 2) == OBSTACLE_SELF) {
    int wallOffset[3];
    // Case 1: when fluid is above the step
    if (i + parameters_.parallel.firstCorner[0] - 2 < xObsCells) {
      wallOffset[1]  = parameters_.parallel.firstCorner[1] + yObsCells;
      RealType distance = (parameters_.geometry.sizeY - (wallOffset[1] + (j - 2 + 0.5))) * parameters_.meshsize->getDy(i, j, k);
      RealType temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) * wallOffset[1];
      if (temp < distance) {
        distance = temp;
      }
      temp = (parameters_.parallel.firstCorner[2] + (k - 2 + 0.5)) * parameters_.meshsize->getDz(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
      temp = ( parameters_.geometry.sizeZ
                        - ( parameters_.parallel.firstCorner[2] + ( k - 2 + 0.5 ) ) )
                        * parameters_.meshsize->getDz(i, j, k);
      if (temp < distance ) {
        distance = temp;
      }

      // Need to verify if this if condition is required
      if (parameters_.simulation.scenario != "channel") {
        temp = ( parameters_.geometry.sizeX
                            - ( parameters_.parallel.firstCorner[0] + ( i - 2 + 0.5 ) ) )
                            * parameters_.meshsize->getDx(i, j, k);
        if ( temp < distance ) {
          distance = temp;
        }
        temp = (parameters_.parallel.firstCorner[0] + (i - 2 + 0.5)) * parameters_.meshsize->getDx(i, j, k);
        if (temp < distance) {
          distance = temp;
        }
      }

      flowField.getWallDistance().getScalar(i, j, k) = distance;
    }
    // Case 2: When fluid is below the step and infront of the step
    if (i + parameters_.parallel.firstCorner[0] - 2 > xObsCells && j + parameters_.parallel.firstCorner[1] - 2 < yObsCells) {
      wallOffset[0]  = parameters_.parallel.firstCorner[0] + xObsCells;
      RealType distance = (parameters_.geometry.sizeY - (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)))
                       * parameters_.meshsize->getDy(i, j, k);
      RealType temp = (parameters_.parallel.firstCorner[1] + (j - 2 + 0.5)) * parameters_.meshsize->getDy(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
      temp = (parameters_.parallel.firstCorner[2] + (k - 2 + 0.5)) * parameters_.meshsize->getDz(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
      temp = ( parameters_.geometry.sizeZ
                        - ( parameters_.parallel.firstCorner[2] + ( k - 2 + 0.5 ) ) )
                        * parameters_.meshsize->getDz(i, j, k);
      if ( temp < distance ) {
        distance = temp;
      }
      temp = (parameters_.geometry.sizeX - (wallOffset[0] + (i - 2 + 0.5))) * parameters_.meshsize->getDx(i, j, k);
      if (temp < distance) {
        distance = temp;
      }
      temp = (wallOffset[0] + (i - 2 + 0.5)) * parameters_.meshsize->getDx(i, j, k);
      if (temp < distance) {
        distance = temp;
      }

      flowField.getWallDistance().getScalar(i, j, k) = distance;
    }
  }
}
