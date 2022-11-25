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

  class PressureBufferReadStencil: public BoundaryStencil<FlowField> {
  public:

    PressureBufferReadStencil(const Parameters& parameters);
    ~PressureBufferReadStencil() override = default;

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

    std::unique_ptr<RealType[]> leftPressureReadBuffer;
    std::unique_ptr<RealType[]> rightPressureReadBuffer;
    std::unique_ptr<RealType[]> topPressureReadBuffer;
    std::unique_ptr<RealType[]> bottomPressureReadBuffer;
    std::unique_ptr<RealType[]> frontPressureReadBuffer;
    std::unique_ptr<RealType[]> backPressureReadBuffer;

    const int* localSize; 
    };
} // namespace Stencils
