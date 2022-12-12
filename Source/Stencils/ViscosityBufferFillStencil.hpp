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
  class ViscosityBufferFillStencil: public BoundaryStencil<TurbulentFlowField> {
  public:
    ViscosityBufferFillStencil(const Parameters& parameters);
    ~ViscosityBufferFillStencil() override = default;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j) override;

    void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) override;
    void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) override;

    std::unique_ptr<RealType[]> leftViscosityFillBuffer;
    std::unique_ptr<RealType[]> rightViscosityFillBuffer;
    std::unique_ptr<RealType[]> topViscosityFillBuffer;
    std::unique_ptr<RealType[]> bottomViscosityFillBuffer;
    std::unique_ptr<RealType[]> frontViscosityFillBuffer;
    std::unique_ptr<RealType[]> backViscosityFillBuffer;

    const int* localSize;
  };
} // namespace Stencils