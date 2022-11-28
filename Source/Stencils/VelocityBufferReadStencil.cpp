#include "StdAfx.hpp" 
#include "VelocityBufferReadStencil.hpp"
#include "Definitions.hpp"

Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {
    leftVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1)
    );
    rightVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1)
    );
    bottomVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1)
    );
    topVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1)
    );
    frontVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1)
    );
    backVelocityReadBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1)
    );

  }

  else {
    leftVelocityReadBuffer   = std::make_unique<RealType[]>(localSize[1] + (localSize[1] + 1));
    rightVelocityReadBuffer  = std::make_unique<RealType[]>(localSize[1] + (localSize[1] + 1));
    bottomVelocityReadBuffer = std::make_unique<RealType[]>(localSize[0] + (localSize[0] + 1));
    topVelocityReadBuffer    = std::make_unique<RealType[]>(localSize[0] + (localSize[0] + 1));
  }
} // End of constructor

//For 2D
void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
if( parameters_.parallel.leftNb >= 0){
  if(j>=2){
   (flowField.getVelocity().getVector(i-1,j))[0] =  *(leftVelocityReadBuffer.get() + (j - 2));
   (flowField.getVelocity().getVector(i,j))[1] = *(leftVelocityReadBuffer.get() + localSize[1] + (j-1));
  }
  else if(j==1){
   (flowField.getVelocity().getVector(i,j))[1] = *(leftVelocityReadBuffer.get() + localSize[1] + (j-1));
  }
}
}

void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
if( parameters_.parallel.rightNb >= 0){
  if(j>=2){
    (flowField.getVelocity().getVector(i+1,j))[0] = *(rightVelocityReadBuffer.get() + (j - 2));
    (flowField.getVelocity().getVector(i+1,j))[1] = *(rightVelocityReadBuffer.get() + localSize[1] + (j-1));
  }
  else if(j==1){
   (flowField.getVelocity().getVector(i+1,j))[1] = *(rightVelocityReadBuffer.get() + localSize[1] + (j-1));
  }
}
}

void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
if( parameters_.parallel.topNb >= 0){
if ((i >= 2) ) {
	(flowField.getVelocity().getVector(i, j+1))[0] = *(topVelocityReadBuffer.get() + (i - 1) );
    (flowField.getVelocity().getVector(i, j+1))[1] = *(topVelocityReadBuffer.get() + localSize[0] + (i - 1));
}
else if  ((i == 1))  {
    (flowField.getVelocity().getVector(i, j+1))[0] = *(topVelocityReadBuffer.get() + (i - 1));
}   

}
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
if( parameters_.parallel.bottomNb >= 0){
    if((i >= 2) ) {
	(flowField.getVelocity().getVector(i, j))[0] = *(bottomVelocityReadBuffer.get() + (i - 1) );
    (flowField.getVelocity().getVector(i, j-1))[1] = *(bottomVelocityReadBuffer.get() + localSize[0] + (i - 1));
}
    else if  ((i == 1))  {
    (flowField.getVelocity().getVector(i, j))[0] = *(bottomVelocityReadBuffer.get() + (i - 1));
}

}
}

//For 3D
void Stencils::VelocityBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) { 
if( parameters_.parallel.leftNb >= 0){
  if(j>=2 && k>=2){
    (flowField.getVelocity().getVector(i-1,j,k))[0] = *(leftVelocityReadBuffer.get() + (j - 2) + (k - 2) * localSize[1]) ;
    (flowField.getVelocity().getVector(i,j,k))[1] = *(leftVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2]) + (k-2));
    (flowField.getVelocity().getVector(i,j,k))[2] = *(leftVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((localSize[1]+1)*localSize[2]) + (j-2) + ((k-1)*localSize[1]));
  }
  else if(j==1 && k>=2){
    (flowField.getVelocity().getVector(i,j,k))[1] = *(leftVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2])+ (k-2));
  }
  else if(k==1 && j>=2){
    (flowField.getVelocity().getVector(i,j,k))[2] = *(leftVelocityReadBuffer.get()+ (localSize[1] * localSize[2])+ ((localSize[1] + 1) * (localSize[2])) + (j - 2)+ ((k - 1) * localSize[1]));
  }

}
}

void Stencils::VelocityBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) { 
if( parameters_.parallel.rightNb >= 0){
if ((j >= 2) && (k >= 2)) {
	(flowField.getVelocity().getVector(i+1,j,k))[0] = *(rightVelocityReadBuffer.get() +(j - 2) + (k - 2) * localSize[1]);
    (flowField.getVelocity().getVector(i+1,j,k))[1] = *(rightVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2]) + (k-2));
    (flowField.getVelocity().getVector(i+1,j,k))[2] = *(rightVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((localSize[1]+1)*localSize[2]) + (j-2) + ((k-1)*localSize[1]));
}
else if ((j == 1) && (k >= 2)) { 
  (flowField.getVelocity().getVector(i+1,j,k))[1] = *(rightVelocityReadBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2])+ (k-2));
 
 }
else if ((k == 1) && (j >= 2)) {
    (flowField.getVelocity().getVector(i+1,j,k))[2] = *(rightVelocityReadBuffer.get()+ (localSize[1] * localSize[2])+ ((localSize[1] + 1) * (localSize[2])) + (j - 2)+ ((k - 1) * localSize[1]));
  }
}
}

void Stencils::VelocityBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) { 
if( parameters_.parallel.topNb >= 0){

if ((i >= 2) && (k >= 2)) {
	(flowField.getVelocity().getVector(i, j+1, k))[0] = *(topVelocityReadBuffer.get() + ((i - 1) * localSize[2]) + (k - 2));
    (flowField.getVelocity().getVector(i, j+1, k))[1] = *(topVelocityReadBuffer.get() + ((localSize[0] + 1) * localSize[2]) + ((i - 2) * localSize[2]) + (k - 2));
    (flowField.getVelocity().getVector(i, j+1, k))[2] = *(topVelocityReadBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) );
    }

else if  ((i == 1) && (k >= 2))  {
    (flowField.getVelocity().getVector(i, j+1, k))[0] = *(topVelocityReadBuffer.get() + ((i - 1) * localSize[2]) + (k - 2));
    }

else if ((k == 1) && (i >= 2)){
  (flowField.getVelocity().getVector(i, j+1, k))[2] = *(topVelocityReadBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) );
    }
    }
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) { 
if( parameters_.parallel.bottomNb >= 0){
    if ((i >= 2) && (k >= 2)) {
	(flowField.getVelocity().getVector(i, j  , k))[0] = *(bottomVelocityReadBuffer.get() + ((i - 1) * localSize[2]) + (k - 2));
    (flowField.getVelocity().getVector(i, j-1, k))[1] = *(bottomVelocityReadBuffer.get() + ((localSize[0] + 1) * localSize[2]) + ((i - 2) * localSize[2]) + (k - 2));
    (flowField.getVelocity().getVector(i, j, k))[2] = *(bottomVelocityReadBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) );
    }

    else if  ((i == 1) && (k >= 2))  {
    (flowField.getVelocity().getVector(i, j, k))[0] = *(bottomVelocityReadBuffer.get() + ((i - 1) * localSize[2]) + (k - 2));
    }

    else if ((k == 1) && (i >= 2)){
  (flowField.getVelocity().getVector(i, j, k))[2] = *(bottomVelocityReadBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) );
    }

    }
}

void Stencils::VelocityBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
if( parameters_.parallel.frontNb >= 0){
   if ((i >= 2) & (j >= 2)) {
    (flowField.getVelocity().getVector(i, j, k  ))[0] = *(frontVelocityReadBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  );
    (flowField.getVelocity().getVector(i, j, k  ))[1] = *(frontVelocityReadBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2));
    (flowField.getVelocity().getVector(i, j, k-1))[2] = *(frontVelocityReadBuffer.get() +  (localSize[0]*(localSize[1]+1)) + ((localSize[0]+1)*localSize[1]) + ((i-2)*localSize[1]) + (j-2));
  }
  else if ((i == 1) & (j >= 2)) {
    (flowField.getVelocity().getVector(i, j, k  ))[0] = *(frontVelocityReadBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  );
  }
  else if ((j == 1) & (i >= 2)) {
    (flowField.getVelocity().getVector(i, j, k  ))[1] = *(frontVelocityReadBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2));
  }

}
}

void Stencils::VelocityBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
if( parameters_.parallel.backNb >= 0){
  if ((i >= 2) & (j >= 2)) {
    (flowField.getVelocity().getVector(i, j, k+1))[0] = *(backVelocityReadBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  );
    (flowField.getVelocity().getVector(i, j, k+1))[1] = *(backVelocityReadBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2));
    (flowField.getVelocity().getVector(i, j, k+1))[2] = *(backVelocityReadBuffer.get() +  (localSize[0]*(localSize[1]+1)) + ((localSize[0]+1)*localSize[1]) + ((i-2)*localSize[1]) + (j-2));
  }
  else if ((i == 1) & (j >= 2)) {
    (flowField.getVelocity().getVector(i, j, k+1))[0] = *(backVelocityReadBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  );
  }
  else if ((j == 1) & (i >= 2)) {
    (flowField.getVelocity().getVector(i, j, k+1))[1] = *(backVelocityReadBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2));
  }

}
}



