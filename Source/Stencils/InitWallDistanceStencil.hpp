#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /**
   * @brief Stencil for initialising the wall distance
   *
   */
  class InitWallDistanceStencil: public FieldStencil<TurbulentFlowField> {
  private:
    const int xObsCells; // Number of obstacle cells in Backwards Facing step in X direction
    const int yObsCells; // Number of obstacle cells in Backwards Facing step in Y direction
    // The step extends in Z direction completely
  public:
    /**
     * @brief Construct a new Init Wall Distance Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    InitWallDistanceStencil(const Parameters& parameters);
    /**
     * @brief function to initialise the wall distance
     *
     * @param flowField data structure which holds the flow quantities
     * @param i index in x
     * @param j index in y
     */
    void apply(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief function to initialise the wall distance
     *
     * @param flowField data structure which holds the flow quantities
     * @param i index in x

     * @param j index in y
     * @param k index in z
     */
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
