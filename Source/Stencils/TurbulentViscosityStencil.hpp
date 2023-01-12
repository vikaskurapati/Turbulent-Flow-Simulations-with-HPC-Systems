#pragma once

#include <string>
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  /**
   * @brief Stencil for Turbulent Viscosity
   *
   */

  class TurbulentViscosityStencil: public FieldStencil<TurbulentFlowField> {
    private:
     std::string method_; //! Method for turbulence - turbulence or turbulence-sa

  public:
    /**
     * @brief Construct a new Turbulent Viscosity Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    TurbulentViscosityStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Turbulent Viscosity Stencil object
     *
     */
    ~TurbulentViscosityStencil() override = default;
    /**
     * @brief function to calculate the turbulent viscosity in 2D
     *
     * @param flowfield data structure holding the flowfield quantities
     * @param i index in x
     * @param j index in y
     */
    void apply(TurbulentFlowField& flowfield, int i, int j) override;
    /**
     * @brief function to calculate the turbulent viscosity in 3D
     *
     * @param flowfield data structure holding the flowfield quantities
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */
    void apply(TurbulentFlowField& flowfield, int i, int j, int k) override;
  };

} // namespace Stencils
