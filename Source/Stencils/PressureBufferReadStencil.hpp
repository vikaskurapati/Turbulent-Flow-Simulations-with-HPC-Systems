#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {
  /**
   * @brief Stencil which fills the Pressure value into the domain from the buffer
   *
   */
  class PressureBufferReadStencil: public BoundaryStencil<FlowField> {
  public:
    /**
     * @brief Construct a new Pressure Buffer Read Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    PressureBufferReadStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Pressure Buffer Read Stencil object
     *
     */
    ~PressureBufferReadStencil() override = default;
    /**
     * @brief Fill values from buffer back into the domains left wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains right wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyRightWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains bottom wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyBottomWall(FlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains top wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyTopWall(FlowField& flowField, int i, int j) override;

    /**
     * @brief Fill values from buffer back into the domains left wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains right wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains bottom wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains top wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains front wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains back wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBackWall(FlowField& flowField, int i, int j, int k) override;
    /**
     * @brief pointer to the readbuffer to the left wall
     *
     */
    std::unique_ptr<RealType[]> leftPressureReadBuffer;
    /**
     * @brief pointer to the readbuffer to the right wall
     *
     */

    std::unique_ptr<RealType[]> rightPressureReadBuffer;
    /**
     * @brief pointer to the readbuffer to the top wall
     *
     */

    std::unique_ptr<RealType[]> topPressureReadBuffer;
    /**
     * @brief pointer to the readbuffer to the bottom wall
     *
     */

    std::unique_ptr<RealType[]> bottomPressureReadBuffer;
    /**
     * @brief pointer to the readbuffer to the front wall
     *
     */

    std::unique_ptr<RealType[]> frontPressureReadBuffer;
    /**
     * @brief pointer to the readbuffer to the back wall
     *
     */

    std::unique_ptr<RealType[]> backPressureReadBuffer;
    /**
     * @brief local size of the domain
     *
     */
    const int* localSize;
  };
} // namespace Stencils
