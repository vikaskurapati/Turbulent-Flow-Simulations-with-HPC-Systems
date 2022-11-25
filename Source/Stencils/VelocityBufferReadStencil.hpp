#pragma once

#include <algorithm>
#include <memory>

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils{
      class VelocityBufferReadStencil: public BoundaryStencil<FlowField> {
  public:

    VelocityBufferReadStencil(const Parameters& parameters);
    ~VelocityBufferReadStencil() override = default;

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

    std::unique_ptr<RealType[]> leftVelocityReadBuffer;
    std::unique_ptr<RealType[]> rightVelocityReadBuffer;
    std::unique_ptr<RealType[]> topVelocityReadBuffer;
    std::unique_ptr<RealType[]> bottomVelocityReadBuffer;
    std::unique_ptr<RealType[]> frontVelocityReadBuffer;
    std::unique_ptr<RealType[]> backVelocityReadBuffer;

    const int* localSize; 

      }

}