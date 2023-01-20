#pragma once
#include "PetscParallelManager.hpp"

#include <mpi.h>
#include <petsclog.h>

#include "FlowField.hpp"
#include "Parameters.hpp"

#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/TransportViscosityBufferFillStencil.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowfield):
  parameters_(parameters),
  flowfield_(flowfield),
  fillVelocityStencil(parameters),
  readVelocityStencil(parameters),
  velocityfillIterator(flowfield, parameters, fillVelocityStencil, 0, 0),
  velocityreadIterator(flowfield, parameters, readVelocityStencil, 0, 0),
  fillPressureStencil(parameters),
  readPressureStencil(parameters),
  pressurefillIterator(flowfield, parameters, fillPressureStencil, 0, 0),
  pressurereadIterator(flowfield, parameters, readPressureStencil, 0, 0) {}

ParallelManagers::PetscTurbulentParallelManager::PetscTurbulentParallelManager(
  Parameters& parameters, TurbulentFlowField& flowfield
):
  PetscParallelManager(parameters, flowfield),
  flowfield_(flowfield),
  fillViscosityStencil(parameters),
  readViscosityStencil(parameters),
  fillTransportViscosityStencil(parameters),
  readTransportViscosityStencil(parameters),
  viscosityfillIterator(flowfield, parameters, fillViscosityStencil, 0, 0),
  viscosityreadIterator(flowfield, parameters, readViscosityStencil, 0, 0),
  transportViscosityfillIterator(flowfield, parameters, fillTransportViscosityStencil, 0, 0),
  transportViscosityreadIterator(flowfield, parameters, readTransportViscosityStencil, 0, 0)
  {}

void ParallelManagers::PetscParallelManager::communicatePressure() {
  pressurefillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];
    const int*  localSize = fillPressureStencil.localSize;

    MPI_Isend(
      // left to right
      fillPressureStencil.leftPressureFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );
    MPI_Isend(
      // front to back
      fillPressureStencil.frontPressureFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      205,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readPressureStencil.backPressureReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      205,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillPressureStencil.backPressureFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      206,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readPressureStencil.frontPressureReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      206,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  if (parameters_.geometry.dim == 2) {
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillPressureStencil.localSize;

    MPI_Isend(
      // left to right
      fillPressureStencil.leftPressureFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );
    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
  pressurereadIterator.iterate();
}

void ParallelManagers::PetscTurbulentParallelManager::communicateViscosity() {
  viscosityfillIterator.iterate();
//****************************************
// For 3D 
//***************************************
  if (parameters_.geometry.dim == 3) {
    
    // for eddy visc
    MPI_Request request[12];
    MPI_Status  status[12];

    const int*  localSize = fillViscosityStencil.localSize; // localSize same for transport visc

    // left to right eddy visc
    MPI_Isend(
      fillViscosityStencil.leftViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      301,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readViscosityStencil.rightViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      301,
      PETSC_COMM_WORLD,
      &request[1]
    );

    // right to left eddy visc
    MPI_Isend(
      fillViscosityStencil.rightViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      302,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readViscosityStencil.leftViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      302,
      PETSC_COMM_WORLD,
      &request[3]
    );

    // top to bottom eddy visc
    MPI_Isend(
      fillViscosityStencil.topViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      303,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readViscosityStencil.bottomViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      303,
      PETSC_COMM_WORLD,
      &request[5]
    );

    // bottom to top eddy visc
    MPI_Isend(
      fillViscosityStencil.bottomViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      304,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readViscosityStencil.topViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      304,
      PETSC_COMM_WORLD,
      &request[7]
    );

    // front to back eddy visc
    MPI_Isend(
      fillViscosityStencil.frontViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      305,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readViscosityStencil.backViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      305,
      PETSC_COMM_WORLD,
      &request[9]
    );

    // back to front eddy visc
    MPI_Isend(
      fillViscosityStencil.backViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      306,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readViscosityStencil.frontViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      306,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
  //********************************************************
  //for 2D
  //********************************************************
  if (parameters_.geometry.dim == 2) {
    
    // for eddy visc
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillViscosityStencil.localSize;

    // left to right eddy visc
    MPI_Isend(
      fillViscosityStencil.leftViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      301,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readViscosityStencil.rightViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      301,
      PETSC_COMM_WORLD,
      &request[1]
    );

    // right to left eddy visc
    MPI_Isend(
      fillViscosityStencil.rightViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      302,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readViscosityStencil.leftViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      302,
      PETSC_COMM_WORLD,
      &request[3]
    );

    // top to bottom eddy visc
    MPI_Isend(
      fillViscosityStencil.topViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      303,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readViscosityStencil.bottomViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      303,
      PETSC_COMM_WORLD,
      &request[5]
    );

    // bottom to top eddy visc
    MPI_Isend(
      fillViscosityStencil.bottomViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      304,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readViscosityStencil.topViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      304,
      PETSC_COMM_WORLD,
      &request[7]
    );

    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
  
  viscosityreadIterator.iterate();
}

void ParallelManagers::PetscTurbulentParallelManager::communicateTransportViscosity() {
  transportViscosityfillIterator.iterate();
//****************************************
// For 3D 
//***************************************
  if (parameters_.geometry.dim == 3) {
    // for transport visc
    MPI_Request transport_request[12];
    MPI_Status  transport_status[12];

    const int*  localSize = fillTransportViscosityStencil.localSize; // localSize same for transport visc

    //left to right transport visc (IF CONDITION REQ FOR SA)
    MPI_Isend(
      fillTransportViscosityStencil.leftTransportViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      3001,
      PETSC_COMM_WORLD,
      &transport_request[0]
    );
    MPI_Irecv(
      readTransportViscosityStencil.rightTransportViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      3001,
      PETSC_COMM_WORLD,
      &transport_request[1]
    );

    // right to left transport visc
    MPI_Isend(
      fillTransportViscosityStencil.rightTransportViscosityFillBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      3002,
      PETSC_COMM_WORLD,
      &transport_request[2]
    );
    MPI_Irecv(
      readTransportViscosityStencil.leftTransportViscosityReadBuffer.get(),
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      3002,
      PETSC_COMM_WORLD,
      &transport_request[3]
    );

    //top to bottom transport visc
    MPI_Isend(
      fillTransportViscosityStencil.topTransportViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3003,
      PETSC_COMM_WORLD,
      &transport_request[4]
    );
    MPI_Irecv(
      readTransportViscosityStencil.bottomTransportViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3003,
      PETSC_COMM_WORLD,
      &transport_request[5]
    );

    // bottom to top transport visc
    MPI_Isend(
      fillTransportViscosityStencil.bottomTransportViscosityFillBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3004,
      PETSC_COMM_WORLD,
      &transport_request[6]
    );
    MPI_Irecv(
      readTransportViscosityStencil.topTransportViscosityReadBuffer.get(),
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3004,
      PETSC_COMM_WORLD,
      &transport_request[7]
    );

    // front to back transport visc
    MPI_Isend(
      fillTransportViscosityStencil.frontTransportViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      3005,
      PETSC_COMM_WORLD,
      &transport_request[8]
    );
    MPI_Irecv(
      readTransportViscosityStencil.backTransportViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      3005,
      PETSC_COMM_WORLD,
      &transport_request[9]
    );

    // back to front transport visc
    MPI_Isend(
      fillTransportViscosityStencil.backTransportViscosityFillBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      3006,
      PETSC_COMM_WORLD,
      &transport_request[10]
    );
    MPI_Irecv(
      readTransportViscosityStencil.frontTransportViscosityReadBuffer.get(),
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      3006,
      PETSC_COMM_WORLD,
      &transport_request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&transport_request[i], &transport_status[i]);

    }
  }
  //********************************************************
  //for 2D
  //********************************************************
  if (parameters_.geometry.dim == 2) {
    
    // for transport visc
    MPI_Request transport_request[8];
    MPI_Status  transport_status[8];

    const int* localSize = fillTransportViscosityStencil.localSize;

    // left to right transport visc
    MPI_Isend(
      fillTransportViscosityStencil.leftTransportViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      3001,
      PETSC_COMM_WORLD,
      &transport_request[0]
    );
    MPI_Irecv(
      readTransportViscosityStencil.rightTransportViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      3001,
      PETSC_COMM_WORLD,
      &transport_request[1]
    );

    // right to left transport visc
    MPI_Isend(
      fillTransportViscosityStencil.rightTransportViscosityFillBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      3002,
      PETSC_COMM_WORLD,
      &transport_request[2]
    );
    MPI_Irecv(
      readTransportViscosityStencil.leftTransportViscosityReadBuffer.get(),
      localSize[1],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      3002,
      PETSC_COMM_WORLD,
      &transport_request[3]
    );

    //top to bottom transport viscs
    MPI_Isend(
      fillTransportViscosityStencil.topTransportViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3003,
      PETSC_COMM_WORLD,
      &transport_request[4]
    );
    MPI_Irecv(
      readTransportViscosityStencil.bottomTransportViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3003,
      PETSC_COMM_WORLD,
      &transport_request[5]
    );

    // bottom to top transport Visc
        MPI_Isend(
      fillTransportViscosityStencil.bottomTransportViscosityFillBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3004,
      PETSC_COMM_WORLD,
      &transport_request[6]
    );
    MPI_Irecv(
      readTransportViscosityStencil.topTransportViscosityReadBuffer.get(),
      localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3004,
      PETSC_COMM_WORLD,
      &transport_request[7]
    );
    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&transport_request[i], &transport_status[i]);

    }
  }
    transportViscosityreadIterator.iterate();



}

void ParallelManagers::PetscParallelManager::communicateVelocity() {
  velocityfillIterator.iterate();
  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];

    const int* localSize = fillVelocityStencil.localSize;

    MPI_Isend(
      // left to right
      fillVelocityStencil.leftVelocityFillBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      101,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readVelocityStencil.rightVelocityReadBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      101,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillVelocityStencil.rightVelocityFillBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      102,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readVelocityStencil.leftVelocityReadBuffer.get(),
      localSize[1] * localSize[2] + (localSize[1] + 1) * localSize[2] + localSize[1] * (localSize[2] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      102,
      PETSC_COMM_WORLD,
      &request[3]
    );

    MPI_Isend(
      // top to bottom
      fillVelocityStencil.topVelocityFillBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      103,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readVelocityStencil.bottomVelocityReadBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      103,
      PETSC_COMM_WORLD,
      &request[5]
    );

    MPI_Isend(
      // bottom to top
      fillVelocityStencil.bottomVelocityFillBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      104,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readVelocityStencil.topVelocityReadBuffer.get(),
      (localSize[0] * localSize[2]) + localSize[0] * (localSize[2] + 1) + (localSize[0] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      104,
      PETSC_COMM_WORLD,
      &request[7]
    );

    MPI_Isend(
      // front to back
      fillVelocityStencil.frontVelocityFillBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      105,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readVelocityStencil.backVelocityReadBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      105,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillVelocityStencil.backVelocityFillBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      106,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readVelocityStencil.frontVelocityReadBuffer.get(),
      (localSize[1] * localSize[0]) + localSize[1] * (localSize[0] + 1) + (localSize[1] + 1) * localSize[0],
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      106,
      PETSC_COMM_WORLD,
      &request[11]
    );

    for (size_t i = 0; i < 12; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  if (parameters_.geometry.dim == 2) {
    MPI_Request request[8];
    MPI_Status  status[8];

    const int* localSize = fillVelocityStencil.localSize;

    MPI_Isend(
      // left to right
      fillVelocityStencil.leftVelocityFillBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      101,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readVelocityStencil.rightVelocityReadBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      101,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillVelocityStencil.rightVelocityFillBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      102,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readVelocityStencil.leftVelocityReadBuffer.get(),
      (localSize[1]) + (localSize[1] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      102,
      PETSC_COMM_WORLD,
      &request[3]
    );

    MPI_Isend(
      // top to bottom
      fillVelocityStencil.topVelocityFillBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      103,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readVelocityStencil.bottomVelocityReadBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      103,
      PETSC_COMM_WORLD,
      &request[5]
    );

    MPI_Isend(
      // bottom to top
      fillVelocityStencil.bottomVelocityFillBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      104,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readVelocityStencil.topVelocityReadBuffer.get(),
      (localSize[0]) + (localSize[0] + 1),
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      104,
      PETSC_COMM_WORLD,
      &request[7]
    );

    for (size_t i = 0; i < 8; i++) {
      MPI_Wait(&request[i], &status[i]);
    }
  }

  velocityreadIterator.iterate();
}