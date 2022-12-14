#pragma once

#include <algorithm>
#include <memory>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  /**
   * @brief Stencil which fills the Viscosity values from the domain into the buffer
   *
   */

  class ViscosityBufferFillStencil: public BoundaryStencil<TurbulentFlowField> {
  public:
    /**
     * @brief Construct a new Viscosity Buffer Fill Stencil object
     *
     * @param parameters parameters of the flow field
     */
    ViscosityBufferFillStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Viscosity Buffer Fill Stencil object
     *
     */
    ~ViscosityBufferFillStencil() override = default;

    /**
     * @brief Fill values from left wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from right wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from bottom wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from top wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from left wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from right wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from bottom wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from top wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from front wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from back wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief pointer to the fillbuffer to the left wall
     *
     */

    std::unique_ptr<RealType[]> leftViscosityFillBuffer;
        /**
     * @brief pointer to the fillbuffer to the right wall
     *
     */

    std::unique_ptr<RealType[]> rightViscosityFillBuffer;
        /**
     * @brief pointer to the fillbuffer to the top wall
     *
     */

    std::unique_ptr<RealType[]> topViscosityFillBuffer;
        /**
     * @brief pointer to the fillbuffer to the bottom wall
     *
     */

    std::unique_ptr<RealType[]> bottomViscosityFillBuffer;
        /**
     * @brief pointer to the fillbuffer to the front wall
     *
     */

    std::unique_ptr<RealType[]> frontViscosityFillBuffer;
        /**
     * @brief pointer to the fillbuffer to the back wall
     *
     */

    std::unique_ptr<RealType[]> backViscosityFillBuffer;

    const int* localSize;
  };
} // namespace Stencils