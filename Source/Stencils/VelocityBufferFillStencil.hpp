#pragma once
#include <algorithm>
#include <memory>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class VelocityBufferFillStencil: public BoundaryStencil<FlowField> {

  public:
    VelocityBufferFillStencil(const Parameters& parameters);
    ~VelocityBufferFillStencil() override = default;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

    std::unique_ptr<RealType[]> leftVelocityFillBuffer;
    std::unique_ptr<RealType[]> rightVelocityFillBuffer;
    std::unique_ptr<RealType[]> topVelocityFillBuffer;
    std::unique_ptr<RealType[]> bottomVelocityFillBuffer;
    std::unique_ptr<RealType[]> frontVelocityFillBuffer;
    std::unique_ptr<RealType[]> backVelocityFillBuffer;

    const int* localSize;

  }; // end of class

} // namespace Stencils
