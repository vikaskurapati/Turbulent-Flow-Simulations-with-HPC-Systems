#include "StdAfx.hpp"

#include "DataStructures.hpp"
#include "FlowField.hpp"

constexpr auto SIZE_X = 20;
constexpr auto SIZE_Y = 25;

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
  spdlog::info("Testing flow field");

  FlowField field(SIZE_X, SIZE_Y);

  if (field.getPressure().getScalar(10, 10) != 0) {
    return EXIT_FAILURE;
  }

  field.getPressure().getScalar(10, 10) = 10;
  if (field.getPressure().getScalar(10, 10) / 2 != 5) {
    return EXIT_FAILURE;
  }

  if (field.getFlags().getNx() != 23) {
    return EXIT_FAILURE;
  }

  if (field.getFlags().getNy() != 28) {
    return EXIT_FAILURE;
  }

  field.getFlags().getValue(10, 10) = 7;
  if (field.getFlags().getValue(10, 10) / 3 != 2) {
    return EXIT_FAILURE;
  }

  spdlog::info("Test for flow field completed successfully");

  return EXIT_SUCCESS;
}
