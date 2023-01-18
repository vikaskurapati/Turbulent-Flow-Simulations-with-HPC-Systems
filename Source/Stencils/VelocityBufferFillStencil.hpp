#pragma once
#include <algorithm>
#include <memory>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {
  /**
   * @brief Stencil which fills the Velocity values from the domain into the buffer. Holds all us, all vs, all ws in
   * that order
   *
   */
  class VelocityBufferFillStencil: public BoundaryStencil<FlowField> {

  public:
    /**
     * @brief Construct a new Velocity Buffer Fill Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    VelocityBufferFillStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Velocity Buffer Fill Stencil object
     *
     */
    ~VelocityBufferFillStencil() override = default;

    /**
     * @brief Fill values from left wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from right wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyRightWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from bottom wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyBottomWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from top wall into the buffer in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyTopWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from left wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from right wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from bottom wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from top wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from front wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from back wall into the buffer in 3D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

    /**
     * @brief pointer to the fillbuffer to the left wall
     *
     */
    std::unique_ptr<RealType[]> leftVelocityFillBuffer;
    /**
     * @brief pointer to the fillbuffer to the right wall
     *
     */

    std::unique_ptr<RealType[]> rightVelocityFillBuffer;
    /**
     * @brief pointer to the fillbuffer to the top wall
     *
     */

    std::unique_ptr<RealType[]> topVelocityFillBuffer;
    /**
     * @brief pointer to the fillbuffer to the bottom wall
     *
     */

    std::unique_ptr<RealType[]> bottomVelocityFillBuffer;
    /**
     * @brief pointer to the fillbuffer to the front wall
     *
     */

    std::unique_ptr<RealType[]> frontVelocityFillBuffer;
    /**
     * @brief pointer to the fillbuffer to the back wall
     *
     */

    std::unique_ptr<RealType[]> backVelocityFillBuffer;
    /**
     * @brief local size of the domain
     *
     */
    const int* localSize;

  }; // end of class

} // namespace Stencils
