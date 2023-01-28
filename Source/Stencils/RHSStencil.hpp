#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {
  class RHSStencil:public FieldStencil<FlowField>{
    public:
    RHSStencil(const Parameters& parameters);
    ~RHSStencil() override = default;

    void apply(FlowField& flowfield, int i, int j) override;
    void apply(FlowField& flowfield, int i, int j, int k) override;
  };

} // namespace Stencils
