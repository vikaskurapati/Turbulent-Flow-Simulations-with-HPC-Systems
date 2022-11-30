#include "StdAfx.hpp"

#include "TurbulentViscosityStencil.hpp"

Stencils::TurbulentViscosityStencil::TurbulentViscosityStencil(const Parameters& parameters):
FieldStencil<FlowField>(parameters){}

void Stencils::TurbulentViscosityStencil::apply(FlowField& flowField, int i, int j)
{

}

void Stencils::TurbulentViscosityStencil::apply(FlowField& flowField, int i, int j, int k)
{
   
}