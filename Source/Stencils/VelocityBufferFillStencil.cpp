#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  localSize(parameters.parallel.localSize) {

  if (parameters.geometry.dim == 3) {
    leftVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1)
    );
    rightVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1)
    );
    bottomVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1)
    );
    topVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[2] + (localSize[0] + 1) * localSize[2] + localSize[0] * (localSize[2] + 1)
    );
    frontVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1)
    );

    backVelocityFillBuffer = std::make_unique<RealType[]>(
      localSize[0] * localSize[1] + (localSize[0] + 1) * localSize[1] + localSize[0] * (localSize[1] + 1)
    );

  }

  else {
    leftVelocityFillBuffer   = std::make_unique<RealType[]>(localSize[1] + (localSize[1] + 1));
    rightVelocityFillBuffer  = std::make_unique<RealType[]>(localSize[1] + (localSize[1] + 1));
    bottomVelocityFillBuffer = std::make_unique<RealType[]>(localSize[0] + (localSize[0] + 1));
    topVelocityFillBuffer    = std::make_unique<RealType[]>(localSize[0] + (localSize[0] + 1));
  }

} // End of constructor

// For 2D   (WHAT ARE i and j, start from there)
void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {  
  //In the domain where we want to exchange the velocities, say we have 4 cells, we need to exchange 4 u , 5 v values
  //Where are i and j starting from? if it is starting from the SW corner(including the ghost cells) Shouldn't it be then i+2 ??
  if(j>=2){
    *(leftVelocityFillBuffer.get() + (j - 2)) = (flowField.getVelocity().getVector(i+1,j))[0];
    *(leftVelocityFillBuffer.get() + localSize[1] + (j-1)) = (flowField.getVelocity().getVector(i+1,j))[1];
  }
  else if(j==1){
    *(leftVelocityFillBuffer.get() + localSize[1] + (j-1)) =(flowField.getVelocity().getVector(i+1,j))[1];
  }
  }

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  if(j>=2){
    *(rightVelocityFillBuffer.get() + (j - 2)) = (flowField.getVelocity().getVector(i-1,j))[0];
    *(rightVelocityFillBuffer.get() + localSize[1] + (j-1)) = (flowField.getVelocity().getVector(i,j))[1];
  }
  else if(j==1){
    *(rightVelocityFillBuffer.get() + localSize[1] + (j-1)) =(flowField.getVelocity().getVector(i,j))[1];
  }
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  if ((i >= 2) ) {
		*(topVelocityFillBuffer.get() + (i - 1) ) = (flowField.getVelocity().getVector(i, j))[0];
    *(topVelocityFillBuffer.get() + localSize[0] + (i - 1)) = (flowField.getVelocity().getVector(i, j-1))[1];
}
else if  ((i == 1))  {
    *(topVelocityFillBuffer.get() + (i - 1)) = (flowField.getVelocity().getVector(i, j))[0];
}
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  if ((i >= 2) ) {
		*(bottomVelocityFillBuffer.get() + (i - 1) ) = (flowField.getVelocity().getVector(i, j+1))[0];
    *(bottomVelocityFillBuffer.get() + localSize[0] + (i - 1)) = (flowField.getVelocity().getVector(i, j+1))[1];
}
else if  ((i == 1))  {
    *(bottomVelocityFillBuffer.get() + (i - 1)) = (flowField.getVelocity().getVector(i, j+1))[0];
}
}
// End of 2D apply functions

// For 3D
void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  
  if(j>=2 && k>=2){
    *(leftVelocityFillBuffer.get() + (j - 2) + (k - 2) * localSize[1]) = (flowField.getVelocity().getVector(i+1,j,k))[0];
    *(leftVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2]) + (k-2)) = (flowField.getVelocity().getVector(i+1,j,k))[1];
    *(leftVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + ((localSize[1]+1)*localSize[2]) + (j-2) + ((k-1)*localSize[1])) = (flowField.getVelocity().getVector(i+1,j,k))[2];
  }
  else if(j==1 && k>=2){
    *(leftVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2])+ (k-2)) =(flowField.getVelocity().getVector(i+1,j,k))[1] ;
  }
  else if(k==1 && j>=2){
    *(leftVelocityFillBuffer.get()+ (localSize[1] * localSize[2])+ ((localSize[1] + 1) * (localSize[2])) + (j - 2)+ ((k - 1) * localSize[1])) = (flowField.getVelocity().getVector(i + 1, j, k))[2];
  }
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  if ((j >= 2) && (k >= 2)) {
		*(rightVelocityFillBuffer.get() +(j - 2) + ((k - 2) * localSize[1])) =(flowField.getVelocity().getVector(i - 1, j, k))[0];
    //check if i-1 in RHS is correct or not for the below two (modified back to --)
    *(rightVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2]) + (k-2)) = (flowField.getVelocity().getVector(i,j,k))[1];
    *(rightVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + (localSize[1]+1)*localSize[2] + (j-2) + (k-1)*localSize[1]) = (flowField.getVelocity().getVector(i,j,k))[2];
}
 else if ((j == 1) && (k >= 2)) {
  //check if i-1 in RHS is correct or not for the below (modified back to --)
  *(rightVelocityFillBuffer.get() + (localSize[1]*localSize[2]) + ((j-1)*localSize[2])+ (k-2)) =(flowField.getVelocity().getVector(i,j,k))[1];
 
 }
  else if ((k == 1) && (j >= 2)) {
    //check if i-1 is correct or not for the below (modified back to --)
    *(rightVelocityFillBuffer.get()+ (localSize[1] * localSize[2])+ ((localSize[1] + 1) * (localSize[2])) + (j - 2)+ ((k - 1) * localSize[1])) = (flowField.getVelocity().getVector(i , j, k))[2];
  }
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
if ((i >= 2) && (k >= 2)) {
		*(topVelocityFillBuffer.get() + ((i - 1) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j, k))[0];
    *(topVelocityFillBuffer.get() + ((localSize[0] + 1) * localSize[2]) + ((i - 2) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j-1, k))[1];
    *(topVelocityFillBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) ) =(flowField.getVelocity().getVector(i, j, k))[2];
}

else if  ((i == 1) && (k >= 2))  {
    *(topVelocityFillBuffer.get() + ((i - 1) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j, k))[0];
}

else if ((k == 1) && (i >= 2)){
  *(topVelocityFillBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) ) =(flowField.getVelocity().getVector(i, j, k))[2];
}
    
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {

if ((i >= 2) && (k >= 2)) {
		*(bottomVelocityFillBuffer.get() + ((i - 1) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j+1, k))[0];
    *(bottomVelocityFillBuffer.get() + ((localSize[0] + 1) * localSize[2]) + ((i - 2) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j+1, k))[1];
    *(bottomVelocityFillBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) ) =(flowField.getVelocity().getVector(i, j+1, k))[2];
}

else if  ((i == 1) && (k >= 2))  {
    *(bottomVelocityFillBuffer.get() + ((i - 1) * localSize[2]) + (k - 2)) = (flowField.getVelocity().getVector(i, j+1, k))[0];
}

else if ((k == 1) && (i >= 2)){
  *(bottomVelocityFillBuffer.get() + ((localSize[0] + 1) *localSize[2]) + (localSize[0] * localSize[2]) + ((k - 1) * localSize[0])+ (i - 2) ) =(flowField.getVelocity().getVector(i, j+1, k))[2];
}
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {

  if ((i >= 2) & (j >= 2)) {
    *(frontVelocityFillBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  )= (flowField.getVelocity().getVector(i, j, k + 1))[0];
    *(frontVelocityFillBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2)) = (flowField.getVelocity().getVector(i, j, k + 1))[1];
    *(frontVelocityFillBuffer.get() +  (localSize[0]*(localSize[1]+1)) + ((localSize[0]+1)*localSize[1]) + ((i-2)*localSize[1]) + (j-2)) = (flowField.getVelocity().getVector(i, j, k + 1))[2];
  }
  else if ((i == 1) & (j >= 2)) {
    *(frontVelocityFillBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  )= (flowField.getVelocity().getVector(i, j, k + 1))[0];
  }
  else if ((j == 1) & (i >= 2)) {
    *(frontVelocityFillBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2)) = (flowField.getVelocity().getVector(i, j, k + 1))[1];
  }

}

void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  if ((i >= 2) & (j >= 2)) {
    *(backVelocityFillBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  )= (flowField.getVelocity().getVector(i, j, k))[0];
    *(backVelocityFillBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2)) = (flowField.getVelocity().getVector(i, j, k))[1];
    *(backVelocityFillBuffer.get() +  (localSize[0]*(localSize[1]+1)) + ((localSize[0]+1)*localSize[1]) + ((i-2)*localSize[1]) + (j-2)) = (flowField.getVelocity().getVector(i, j, k - 1))[2];
  }
  else if ((i == 1) & (j >= 2)) {
    *(backVelocityFillBuffer.get() + ((i - 1) * localSize[1]) + (j - 2)  )= (flowField.getVelocity().getVector(i, j, k))[0];
  }
  else if ((j == 1) & (i >= 2)) {
    *(backVelocityFillBuffer.get() +  ((localSize[0] + 1) *localSize[1])+ ((j - 1) * localSize[0]) + (i - 2)) = (flowField.getVelocity().getVector(i, j, k))[1];
  }

}

// End of 3D apply functions