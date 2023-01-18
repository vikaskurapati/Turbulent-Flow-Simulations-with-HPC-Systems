#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
/**
 * @brief Stencil to calculate the maximum turbulent viscosity
 * 
 */
  class MaxNuStencil: public FieldStencil<TurbulentFlowField> {
  private:
    RealType maxNuValue_; //! Stores the maximum module of every component

    /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
     *
     * 2D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     */
    void cellMaxNuValue(TurbulentFlowField& flowField, int i, int j);

    /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
     *
     * 3D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     * @param k Position in the Z direction.
     */
    void cellMaxNuValue(TurbulentFlowField& flowField, int i, int j, int k);

    public:
    /**
     * @brief Construct a new Max Nu Stencil object
     * 
     * @param parameters parameters of the flow simulation
     */
    MaxNuStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Max Nu Stencil object
     * 
     */
    ~MaxNuStencil() override = default;
    /**
     * @brief function which calculates the maximum turbulent viscosity in 2D
     * 
     * @param flowField datastructure holding all flow quantities
     * @param i index in x
     * @param j index in y
     */
    void apply(TurbulentFlowField& flowField, int i, int j) override;

    /**
     * @brief function which calculates the maximum turbulent viscosity in 3D
     * 
     * @param flowField datastructure holding all flow quantities
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    /** Resets the maximum values to zero before computing the timestep.
     */
    void reset_Nu();

    /** Returns the array with the maximum modules of the components of the velocity,
     *  divided by the respective local meshsize.
     */
     RealType getMaxNuValues() const;

  }; //End of class

} //End of Namespace


