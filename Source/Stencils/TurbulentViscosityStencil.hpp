#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  class TurbulentViscosityStencil: public FieldStencil<TurbulentFlowField>{
    public:
    TurbulentViscosityStencil(const Parameters& parameters);
    ~TurbulentViscosityStencil() override = default;

    void apply(TurbulentFlowField& flowfield, int i, int j) override;
    void apply(TurbulentFlowField& flowfield, int i, int j, int k) override;
  };

} // namespace Stencils
