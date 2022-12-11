#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /** Stencil for initialising the wall distance
   */
  class InitWallDistanceStencil: public FieldStencil<TurbulentFlowField> {
  private:
    const int xObsCells;  // Number of obstacle cells in Backwards Facing step in X direction
    const int yObsCells;  // Number of obstacle cells in Backwards Facing step in Y direction
    // The step extends in Z direction completely
  public:
    InitWallDistanceStencil(const Parameters& parameters);

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
