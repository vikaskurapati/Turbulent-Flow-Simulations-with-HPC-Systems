#include "StdAfx.hpp"

#include "InitWallDistanceStencil.hpp"

Stencils::InitWallDistanceStencil::InitWallDistanceStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  xObsCells(parameters.bfStep.xRatio * parameters.geometry.sizeX),
  yObsCells(parameters.bfStep.yRatio * parameters.geometry.sizeY) {}

void Stencils::InitWallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j) {

  // Location of center
  RealType xPos = parameters_.meshsize->getPosX(i, j) + (0.5 * parameters_.meshsize->getDx(i, j));
  RealType yPos = parameters_.meshsize->getPosY(i, j) + (0.5 * parameters_.meshsize->getDy(i, j));

  // For Channel Flow
  if (xObsCells == 0 && yObsCells == 0) {
    flowField.getWallDistance().getScalar(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos));
  }
  // For BFS
  else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.getWallDistance().getScalar(i, j) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY), (parameters_.geometry.lengthY - yPos)
      );
    } else if (yPos < parameters_.bfStep.yRatio * parameters_.geometry.lengthY) {
      flowField.getWallDistance().getScalar(i, j) = std::min(
        xPos - (parameters_.bfStep.xRatio * parameters_.geometry.lengthX),
        (std::min(yPos, (parameters_.geometry.lengthY - yPos)))
      );
    } else {
      flowField.getWallDistance().getScalar(i, j) = std::min(yPos, (parameters_.geometry.lengthY - yPos));
    }
  }
}

void Stencils::InitWallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {

  RealType xPos = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
  RealType yPos = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i, j, k);
  RealType zPos = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i, j, k);

  if (parameters_.bfStep.xRatio < 0 || parameters_.bfStep.yRatio < 0) {
    flowField.getWallDistance().getScalar(i, j, k) = std::min(
      yPos, std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, parameters_.geometry.lengthZ - zPos))
    );
    // *Remove x pos for cells near to inlet*
  } else {
    if (xPos <= parameters_.bfStep.xRatio * parameters_.geometry.lengthX) {
      flowField.getWallDistance().getScalar(i, j, k) = std::min(
        yPos - (parameters_.bfStep.yRatio * parameters_.geometry.lengthY),
        std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, (parameters_.geometry.lengthZ - zPos)))
      );
    } else if (yPos < parameters_.bfStep.yRatio * parameters_.geometry.lengthY) {
      flowField.getWallDistance().getScalar(i, j, k) = std::min(
        xPos - (parameters_.bfStep.xRatio * parameters_.geometry.lengthX),
        std::min(
          yPos, std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, (parameters_.geometry.lengthZ - zPos)))
        )
      );
    } else {
      flowField.getWallDistance().getScalar(i, j, k) = std::min(
        yPos, std::min((parameters_.geometry.lengthY - yPos), std::min(zPos, parameters_.geometry.lengthZ - zPos))
      );
    }
  }
}
