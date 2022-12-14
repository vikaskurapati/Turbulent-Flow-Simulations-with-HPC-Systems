#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  class FGHStencil: public FieldStencil<FlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    /**
     * @brief Construct a new FGHStencil object
     *
     * @param parameters parameters of the flow simulation
     */
    FGHStencil(const Parameters& parameters);
    /**
     * @brief Destroy the FGHStencil object
     *
     */
    ~FGHStencil() override = default;

    /**
     * @brief apply function on the cell of the flow field
     *
     * @param flowField data structure holding all the flow field quantities
     * @param i index in x
     * @param j index in y
     */

    void apply(FlowField& flowField, int i, int j) override;
    /**
     * @brief apply function on the cell of the flow field
     *
     * @param flowField data structure holding all the flow field quantities
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

  class TurbulentFGHStencil: public FieldStencil<TurbulentFlowField> {
    /**
     * @brief FGH Stencil which takes care of turbulence
     *
     */
  private:
    /**
     * @brief local velocity variable to hold values of self and neighbours to calculate the required derivatives
     *
     */
    RealType localVelocity_[27 * 3];
    /**
     * @brief local mesh size variable to hold values of self and neighbours to calculate the required derivatives
     *
     */
    RealType localMeshsize_[27 * 3];

  public:
    /**
     * @brief Construct a new Turbulent F G H Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    TurbulentFGHStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Turbulent F G H Stencil object
     *
     */
    ~TurbulentFGHStencil() override = default;
    /**
     * @brief apply function on the cell of the flow field
     *
     * @param flowField data structure holding all the flow field quantities
     * @param i index in x
     * @param j index in y
     */
    void apply(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief apply function on the cell of the flow field
     *
     * @param flowField data structure holding all the flow field quantities
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
