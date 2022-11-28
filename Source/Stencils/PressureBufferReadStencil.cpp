#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Definitions.hpp"

Stencils::PressureBufferReadStencil::PressureBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {

    leftPressureReadBuffer  = std::make_unique<RealType[]>(localSize[1] * localSize[2]);
    rightPressureReadBuffer = std::make_unique<RealType[]>(localSize[1] * localSize[2]);

    bottomPressureReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[2]);
    topPressureReadBuffer    = std::make_unique<RealType[]>(localSize[0] * localSize[2]);

    frontPressureReadBuffer = std::make_unique<RealType[]>(localSize[0] * localSize[1]);
    backPressureReadBuffer  = std::make_unique<RealType[]>(localSize[0] * localSize[1]);

  }

  else {
    leftPressureReadBuffer  = std::make_unique<RealType[]>(localSize[1]);
    rightPressureReadBuffer = std::make_unique<RealType[]>(localSize[1]);

    bottomPressureReadBuffer = std::make_unique<RealType[]>(localSize[0]);
    topPressureReadBuffer    = std::make_unique<RealType[]>(localSize[0]);
  }
} // End of constructor

// For 2D Cases
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
    if(parameters_.parallel.leftNb >=0 ){
    if (j >= 2) {
    flowField.getPressure().getScalar(i + 1, j) = *(leftPressureReadBuffer.get() + (j - 2));
  }
    }
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
    if(parameters_.parallel.rightNb >=0 ){
    if (j >= 2) {
    //Need to verify indices
     flowField.getPressure().getScalar(i, j) = *(rightPressureReadBuffer.get() + (j - 2));
  }
    }
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
      if(parameters_.parallel.bottomNb >=0 ){
      if ((i >= 2)) {
     flowField.getPressure().getScalar(i, j + 1) = *(bottomPressureReadBuffer.get() + (i - 2) * localSize[2]);
  }
      }
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
      if(parameters_.parallel.topNb >=0 ){
      if ((i >= 2)) {
        //Need to verify indices
     flowField.getPressure().getScalar(i, j)= *(topPressureReadBuffer.get() + (i - 2) * localSize[2]);
  }
      }
}

//For 3D Cases
void Stencils::PressureBufferReadStencil::applyLeftWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.leftNb >=0 ){
        if((j>=2) && (k>=2)){
            //Check indices
            flowfield.getPressure().getScalar(i,j,k) = *(leftPressureReadBuffer.get() + (j-2) + (k-2)*localSize[1]);
        }
    }
}

void Stencils::PressureBufferReadStencil::applyRightWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.rightNb >=0 ){
        if((j>=2) && (k>=2)){
            flowfield.getPressure().getScalar(i + 1,j,k) = *(rightPressureReadBuffer.get() + (j-2) + (k-2)*localSize[1]);
        }
    }
}

void Stencils::PressureBufferReadStencil::applyBottomWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.bottomNb >=0 ){
        if((i>=2) && (k>=2)){
            //Check indices
            flowfield.getPressure().getScalar(i,j,k) = *(bottomPressureReadBuffer.get() + (k-2) + (i-2)*localSize[2]);
        }
    }
}

void Stencils::PressureBufferReadStencil::applyTopWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.topNb >=0 ){
        if((i>=2) && (k>=2)){
            flowfield.getPressure().getScalar(i,j + 1,k) = *(topPressureReadBuffer.get() + (k-2) + (i-2)*localSize[2]);
        }
    }
}

void Stencils::PressureBufferReadStencil::applyFrontWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.frontNb >=0 ){
        if((i>=2) & (j>=2)){
            //Check Indices
            flowfield.getPressure().getScalar(i,j,k) = *(frontPressureReadBuffer.get() + (j-2) + (i-2)*localSize[1]);
        }
    }
}

void Stencils::PressureBufferReadStencil::applyBackWall(FlowField & flowfield, int i, int j , int k){
    if(parameters_.parallel.backNb >=0 ){
        if((i>=2) & (j>=2)){
            flowfield.getPressure().getScalar(i,j,k +1) = *(backPressureReadBuffer.get() + (j-2) + (i-2)*localSize[1]);
        }
    }
}

