#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

/**
 * @brief Stencil to initialise the boundary layer thickness
 * 
 */
  class InitBoundaryLayerThickness: public FieldStencil<TurbulentFlowField> {
  public:
  /**
   * @brief Construct a new Init Boundary Layer Thickness object
   * 
   * @param parameters parameters of the flow simulation
   */
    InitBoundaryLayerThickness(const Parameters& parameters);
    /**
     * @brief function to initialise the boundary layer thickness
     * 
     * @param flowField flowField which holds the flow quantities
     * @param i index in x
     * @param j index in y
     */

    void apply(TurbulentFlowField& flowField, int i, int j) override;
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
