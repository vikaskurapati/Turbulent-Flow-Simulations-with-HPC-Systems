#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil for initialising the Boundary layer thickness
   */
  class InitBoundaryLayerThickness: public FieldStencil<FlowField> {
  public:
    InitBoundaryLayerThickness(const Parameters& parameters);

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils