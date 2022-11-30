#include "StdAfx.hpp"

#include "InitWallDistanceStencil.hpp"

Stencils::InitWallDistanceStencil::InitWallDistanceStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters){}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j) {
 
}

void Stencils::InitWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {

}

