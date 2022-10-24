#include "StdAfx.hpp"

#include "DataStructures.hpp"

constexpr auto SIZE_X = 10;
constexpr auto SIZE_Y = 10;
constexpr auto SIZE_Z = 10;

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  spdlog::info("Testing scalar fields");

  ScalarField sfield2D(SIZE_X, SIZE_Y);
  ScalarField sfield3D(SIZE_X, SIZE_Y, SIZE_Z);

  RealType d2counter = 1.0;
  RealType d3counter = 1.0;

  // Fill the fields completely with stuff
  for (int i = 0; i < SIZE_X; i++) {
    for (int j = 0; j < SIZE_Y; j++) {
      sfield2D.getScalar(i, j) = d2counter;
      d2counter += 0.5;

      for (int k = 0; k < SIZE_Z; k++) {
        sfield3D.getScalar(i, j, k) = d3counter;
        d3counter += 0.5;
      }
    }
  }

  // Now read and see if the same values come out
  d2counter = 1.0;
  d3counter = 1.0;

  for (int i = 0; i < SIZE_X; i++) {
    for (int j = 0; j < SIZE_Y; j++) {
      if (sfield2D.getScalar(i, j) != d2counter) {
        spdlog::error("Error while reading the scalar values for 2D field");
        return EXIT_FAILURE;
      }
      d2counter += 0.5;

      for (int k = 0; k < SIZE_Z; k++) {
        if (sfield3D.getScalar(i, j, k) != d3counter) {
          spdlog::error("Error while reading the scalar values for 3D field");
          return EXIT_FAILURE;
        }
        d3counter += 0.5;
      }
    }
  }

  spdlog::info("Test for scalar fields completed successfully");

  return EXIT_SUCCESS;
}
