#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {
  /**
   * @brief Stencil which fills the Viscosity values into the domain from the buffer
   *
   */

  class ViscosityBufferReadStencil: public BoundaryStencil<TurbulentFlowField> {
  public:
    /**
     * @brief Construct a new Viscosity Buffer Read Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    ViscosityBufferReadStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Viscosity Buffer Read Stencil object
     *
     */
    ~ViscosityBufferReadStencil() override = default;
    /**
     * @brief Fill values from buffer back into the domains left wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains right wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains bottom wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains top wall in 2D
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     */

    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;
    /**
     * @brief Fill values from buffer back into the domains left wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains right wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains bottom wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains top wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains front wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief Fill values from buffer back into the domains back wall in 3d
     *
     * @param flowField datastructure holding flow quanitites
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */

    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    /**
     * @brief pointer to the readbuffer to the left wall
     *
     */

    std::unique_ptr<RealType[]> leftViscosityReadBuffer;
    /**
     * @brief pointer to the readbuffer to the right wall
     *
     */

    std::unique_ptr<RealType[]> rightViscosityReadBuffer;
    /**
     * @brief pointer to the readbuffer to the top wall
     *
     */

    std::unique_ptr<RealType[]> topViscosityReadBuffer;
    /**
     * @brief pointer to the readbuffer to the bottom wall
     *
     */

    std::unique_ptr<RealType[]> bottomViscosityReadBuffer;
    /**
     * @brief pointer to the readbuffer to the front wall
     *
     */

    std::unique_ptr<RealType[]> frontViscosityReadBuffer;
    /**
     * @brief pointer to the readbuffer to the back wall
     *
     */

    std::unique_ptr<RealType[]> backViscosityReadBuffer;
    /**
     * @brief local size of the domain
     *
     */

    const int* localSize;
  };
} // namespace Stencils
