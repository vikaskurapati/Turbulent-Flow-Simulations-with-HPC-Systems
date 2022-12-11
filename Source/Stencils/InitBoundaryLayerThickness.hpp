#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /** Stencil for initialising the Boundary layer thickness
   */
  class InitBoundaryLayerThickness: public FieldStencil<TurbulentFlowField> {
  public:
    InitBoundaryLayerThickness(const Parameters& parameters);

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
