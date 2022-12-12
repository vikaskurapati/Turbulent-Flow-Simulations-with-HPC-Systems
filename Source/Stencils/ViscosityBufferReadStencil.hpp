#pragma once

#include "BoundaryStencil.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace Stencils {

  class ViscosityBufferReadStencil: public BoundaryStencil<TurbulentFlowField> {
  public:

    ViscosityBufferReadStencil(const Parameters& parameters);
    ~ViscosityBufferReadStencil() override = default;

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

    std::unique_ptr<RealType[]> leftViscosityReadBuffer;
    std::unique_ptr<RealType[]> rightViscosityReadBuffer;
    std::unique_ptr<RealType[]> topViscosityReadBuffer;
    std::unique_ptr<RealType[]> bottomViscosityReadBuffer;
    std::unique_ptr<RealType[]> frontViscosityReadBuffer;
    std::unique_ptr<RealType[]> backViscosityReadBuffer;

    const int* localSize; 
    };
} // namespace Stencils
