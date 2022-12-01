#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class TurbulentViscosityStencil: public FieldStencil<FlowField>{
    public:
    TurbulentViscosityStencil(const Parameters& parameters);
    ~TurbulentViscosityStencil() override = default;

    void apply(FlowField& flowfield, int i, int j) override;
    void apply(FlowField& flowfield, int i, int j, int k) override;
  };

} // namespace Stencils
