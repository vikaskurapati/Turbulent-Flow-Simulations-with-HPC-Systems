#include "StdAfx.hpp"

#include "RHSStencil.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
FieldStencil<FlowField>(parameters){}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j)
{
    RealType& value = flowField.getRHS().getScalar(i, j);

    // value = //calculation for RHS
    auto dx = 0.5*(parameters_.meshsize->getDx(i,j)+parameters_.meshsize->getDx(i-1,j));
    auto dy = 0.5*(parameters_.meshsize->getDy(i,j)+parameters_.meshsize->getDy(i,j-1));
    auto dt = parameters_.timestep.dt;

    const RealType* fgh = flowField.getFGH().getVector(i, j);
    const RealType* fgh_x = flowField.getFGH().getVector(i-1,j);
    const RealType* fgh_y = flowField.getFGH().getVector(i, j-1);

    value = (((fgh[0]-fgh_x[0])/dx) + ((fgh[1]-fgh_y[1])/dy))/dt;

}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k)
{
    RealType& value = flowField.getRHS().getScalar(i, j, k);

    auto dx = 0.5*(parameters_.meshsize->getDx(i, j, k)+parameters_.meshsize->getDx(i-1, j, k));
    auto dy = 0.5*(parameters_.meshsize->getDy(i, j, k)+parameters_.meshsize->getDy(i, j-1, k));
    auto dz = 0.5*(parameters_.meshsize->getDz(i, j, k)+parameters_.meshsize->getDx(i, j, k-1));
    auto dt = parameters_.timestep.dt;

    const RealType* fgh = flowField.getFGH().getVector(i,j,k);
    const RealType* fgh_x = flowField.getFGH().getVector(i-1, j, k);
    const RealType* fgh_y = flowField.getFGH().getVector(i, j-1, k);
    const RealType* fgh_z = flowField.getFGH().getVector(i, j, k-1);

    value = (((fgh[0]-fgh_x[0])/dx) + ((fgh[1]-fgh_y[1])/dy) + ((fgh[2]-fgh_z[2])/dz))/dt;
}