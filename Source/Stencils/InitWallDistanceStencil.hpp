#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil for initialising the wall distance
   */
  class InitWallDistanceStencil: public FieldStencil<FlowField> {
  public:
    InitWallDistanceStencil(const Parameters& parameters);

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
