#include "PetscParallelManager.hpp"

#include <mpi.h>
#include <petsclog.h>

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowfield):
  parameters_(parameters),
  flowfield_(flowfield),
  fillVelocityStencil(parameters),
  readVelocityStencil(parameters),
  velocityfillIterator(flowfield, parameters, fillVelocityStencil, 1, -1),
  velocityreadIterator(flowfield, parameters, readVelocityStencil, 1, -1),
  fillPressureStencil(parameters),
  readPressureStencil(parameters),
  pressurefillIterator(flowfield, parameters, fillPressureStencil, 1, -1),
  pressurereadIterator(flowfield, parameters, readPressureStencil, 1, -1) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {
  pressurefillIterator.iterate();

  if (parameters_.geometry.dim == 3) {
    MPI_Request request[12];
    MPI_Status  status[12];
    const int*  localSize = fillPressureStencil.localSize;

    MPI_Isend(
      // left to right
      fillPressureStencil.leftPressureFillBuffer,
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer,
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer,
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer,
      localSize[1] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer,
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer,
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer,
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer,
      localSize[0] * localSize[2],
      MY_MPI_FLOAT,
      _parameters.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );
    MPI_Isend(
      // front to back
      fillPressureStencil.frontPressureFillBuffer,
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.frontNb,
      205,
      PETSC_COMM_WORLD,
      &request[8]
    );
    MPI_Irecv(
      readPressureStencil.backPressureReadBuffer,
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.backNb,
      205,
      PETSC_COMM_WORLD,
      &request[9]
    );
    MPI_Isend(
      // back to front
      fillPressureStencil.backPressureFillBuffer,
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.backNb,
      206,
      PETSC_COMM_WORLD,
      &request[10]
    );
    MPI_Irecv(
      readPressureStencil.frontPressureReadBuffer,
      localSize[0] * localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.frontNb,
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
      fillPressureStencil.leftPressureFillBuffer,
      localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.leftNb,
      201,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readPressureStencil.rightPressureReadBuffer,
      localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.rightNb,
      201,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillPressureStencil.rightPressureFillBuffer,
      localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.rightNb,
      202,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readPressureStencil.leftPressureReadBuffer,
      localSize[1],
      MY_MPI_FLOAT,
      _parameters.parallel.leftNb,
      202,
      PETSC_COMM_WORLD,
      &request[3]
    );
    MPI_Isend(
      // top to bottom
      fillPressureStencil.topPressureFillBuffer,
      localSize[0],
      MY_MPI_FLOAT,
      _parameters.parallel.topNb,
      203,
      PETSC_COMM_WORLD,
      &request[4]
    );
    MPI_Irecv(
      readPressureStencil.bottomPressureReadBuffer,
      localSize[0],
      MY_MPI_FLOAT,
      _parameters.parallel.bottomNb,
      203,
      PETSC_COMM_WORLD,
      &request[5]
    );
    MPI_Isend(
      // bottom to top
      fillPressureStencil.bottomPressureFillBuffer,
      localSize[0],
      MY_MPI_FLOAT,
      _parameters.parallel.bottomNb,
      204,
      PETSC_COMM_WORLD,
      &request[6]
    );
    MPI_Irecv(
      readPressureStencil.topPressureReadBuffer,
      localSize[0],
      MY_MPI_FLOAT,
      _parameters.parallel.topNb,
      204,
      PETSC_COMM_WORLD,
      &request[7]
    );
  }
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
      (localSize[1] * localSize[2]) + localSize[1] * (localSize[2] + 1) + (localSize[1] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      101,
      PETSC_COMM_WORLD,
      &request[0]
    );
    MPI_Irecv(
      readVelocityStencil.rightVelocityReadBuffer.get(),
      (localSize[1] * localSize[2]) + localSize[1] * (localSize[2] + 1) + (localSize[1] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      101,
      PETSC_COMM_WORLD,
      &request[1]
    );

    MPI_Isend(
      // right to left
      fillVelocityStencil.rightVelocityFillBuffer.get(),
      (localSize[1] * localSize[2]) + localSize[1] * (localSize[2] + 1) + (localSize[1] + 1) * localSize[2],
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      102,
      PETSC_COMM_WORLD,
      &request[2]
    );
    MPI_Irecv(
      readVelocityStencil.leftVelocityReadBuffer.get(),
      (localSize[1] * localSize[2]) + localSize[1] * (localSize[2] + 1) + (localSize[1] + 1) * localSize[2],
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