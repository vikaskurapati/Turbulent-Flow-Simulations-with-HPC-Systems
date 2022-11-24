#pragma once

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace Stencils {

  class PressureBufferFillStencil: public BoundaryStencil<FlowField> {
  public:
    PressureBufferFillStencil(const Parameters& parameters);
    ~PressureBufferFillStencil() override = default;

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

    std::unique_ptr<RealType[]> leftPressureBuffer;
    std::unique_ptr<RealType[]> rightPressureBuffer;
    std::unique_ptr<RealType[]> topPressureBuffer;
    std::unique_ptr<RealType[]> bottomPressureBuffer;
    std::unique_ptr<RealType[]> frontPressureBuffer;
    std::unique_ptr<RealType[]> backPressureBuffer;
    
    const int* localSize;

  };
} // namespace Stencils
