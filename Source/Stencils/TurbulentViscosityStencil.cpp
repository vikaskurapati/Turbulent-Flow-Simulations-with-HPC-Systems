#include "StdAfx.hpp"

#include "TurbulentViscosityStencil.hpp"

#include <algorithm>
#include <cmath>

#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"

Stencils::TurbulentViscosityStencil::TurbulentViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  method_(parameters.simulation.type) {}

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  //std::cout<<i<<"    "<<j<<"     "<<flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)<<"   "<<flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j)<<std::endl;
  // Do it for fluid cells only

  if ((obstacle && OBSTACLE_SELF) == 0) {

    if (method_ == "turbulence-sa") {

      //std::cout<<flowField.getTurbulentViscosityTransport().getScalar(i, j)<<std::endl;

      RealType W_12 = 0.5*(((flowField.getVelocity().getVector(i, j+1)[0] - flowField.getVelocity().getVector(i, j)[0])/(0.5*(parameters_.meshsize->getDy(i,j+1)+parameters_.meshsize->getDy(i,j)))) - ((flowField.getVelocity().getVector(i+1, j)[1] - flowField.getVelocity().getVector(i, j)[1])/(0.5*(parameters_.meshsize->getDx(i+1,j)+parameters_.meshsize->getDx(i,j)))));

      RealType chi = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) * parameters_.flow.Re;

      RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

      RealType f_v1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(7.1, 3.0));

      RealType f_v2 = 1.0 - (chi / (1 + (chi * f_v1)));

      RealType S_hat = 2.0 * std::sqrt((W_12 * W_12)) + (flowField.getPreviousTurbulentViscosityTransport().getScalar(i,j)*f_v2)/((parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j)*flowField.getWallDistance().getScalar(i,j)));

      RealType r = std::min(10.0,(flowField.getPreviousTurbulentViscosityTransport().getScalar(i,j))/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j)*flowField.getWallDistance().getScalar(i,j))));

      //std::cout<<i<<"   "<<j<<"   "<<S_hat<<"   "<<(2.0 * std::sqrt((W_12 * W_12))*0.3)<<std::endl;

      RealType g = r + 0.3 * (std::pow(r, 6.0) - r);

      RealType f_w = g * std::pow(((1 + std::pow(2.0, 6.0)) / (std::pow(g, 6.0) + std::pow(2, 6.0))), 1.0 / 6.0);
          // std::cout << "Here:   " <<r<<"   "<<(flowField.getPreviousTurbulentViscosityTransport().getScalar(i,j))/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j)*flowField.getWallDistance().getScalar(i,j)))<<std::endl;

      //Term 1
      RealType dx1 = 0.5 * (parameters_.meshsize->getDx(i, j) + parameters_.meshsize->getDx(i + 1, j));
      RealType dx0 = 0.5 * (parameters_.meshsize->getDx(i - 1, j) + parameters_.meshsize->getDx(i, j));

      RealType term1
        = 0.5 * (1 / dx1) * (flowField.getVelocity().getVector(i, j)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j));

      term1 -=  0.5*(1 / dx0) * (flowField.getVelocity().getVector(i-1, j)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i-1, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dx1) * std::fabs(flowField.getVelocity().getVector(i, j)[0])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) - flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dx0) * std::fabs(flowField.getVelocity().getVector(i - 1, j)[0])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
               ));

      RealType dy1 = 0.5 * (parameters_.meshsize->getDy(i, j) + parameters_.meshsize->getDy(i, j + 1));
      RealType dy0 = 0.5 * (parameters_.meshsize->getDy(i, j - 1) + parameters_.meshsize->getDy(i, j));

      term1
        += 0.5 * (1 / dy1) * (flowField.getVelocity().getVector(i, j)[1])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1));

      term1 = term1 - 0.5*(1 / dy0) * (flowField.getVelocity().getVector(i, j-1)[1])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dy1) * std::fabs(flowField.getVelocity().getVector(i, j)[1])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dy0) * std::fabs(flowField.getVelocity().getVector(i, j - 1)[1])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j - 1)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
               ));
      
      //Term 2         
      RealType term2 = 0.1355 * (1 - f_t2) * S_hat * flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j);

      RealType C_w1 = (0.1355) / (parameters_.turbulence.kappa * parameters_.turbulence.kappa)
                      + 
                      ((1 + 0.622) / (2.0 / 3.0));

      //Term 3
      RealType term3
        = (C_w1 * f_w - (0.1355 * f_t2 / (parameters_.turbulence.kappa * parameters_.turbulence.kappa)))
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / flowField.getWallDistance().getScalar(i, j))
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / flowField.getWallDistance().getScalar(i, j));

    // Term 4
      dx1 =  parameters_.meshsize->getDx(i + 1, j);
      dx0 =  parameters_.meshsize->getDx(i, j);

      dy1 =  parameters_.meshsize->getDy(i, j + 1);
      dy0 =  parameters_.meshsize->getDy(i, j);
      
      RealType viscosity_laplacian = flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j) / (dx1 * (dx1 + dx0));

      viscosity_laplacian -= flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / (dx1 * dx0);

      viscosity_laplacian += flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j) / (dx0 * (dx1 + dx0));

      viscosity_laplacian += flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1) / (dy1 * (dy1 + dy0));

      viscosity_laplacian -= flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / (dy1 * dy0);

      viscosity_laplacian += flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1) / (dy0 * (dy1 + dy0));

      viscosity_laplacian = 2.0 * viscosity_laplacian;

      RealType term4 = ((1 / parameters_.flow.Re)+flowField.getPreviousTurbulentViscosityTransport().getScalar(i,j))* viscosity_laplacian;

      //CHECK THE indices i,j in the DIVISION BY dx
      RealType viscgradsquare = ((flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
                                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j))
                                / (parameters_.meshsize->getDx(i, j)))
                                *
                                ((flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
                                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j))
                                / (parameters_.meshsize->getDx(i, j)));

      viscgradsquare += ((flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
                                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1))
                                / (parameters_.meshsize->getDy(i, j)))
                                *
                                ((flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
                                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1))
                                / (parameters_.meshsize->getDy(i, j)));

      // std::cout << "Here:   " <<viscosity_laplacian <<"   "<<viscgradsquare <<std::endl;

      term4 = term4 + (1 + 0.622) * viscgradsquare;

      term4 = term4 / (2.0 / 3.0);

      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) +  parameters_.timestep.dt * (term2 - term3 + term4 - term1);

      //BOundary conditions for additional ghost layer on bottom and left
      // if(i==1 || j==1){
      //   flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) = 0.0;
      // }

      flowField.getTurbulentViscosity().getScalar(i, j) = f_v1 * flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);

      // if(i==2 && j==10){
      // //std::cout << "Here:   " <<term1<<"   "<<term2<<"   "<<term3<<"    "<<term4 <<std::endl;
      // }//BOUNDARY CONDITION 
      
      // std::cout << "Here:   " <<i<<"   "<<j<<"    "<<term1<<"   "<<term2<<"   "<<term3<<"    "<<term4 <<"   "<<flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j)<<std::endl;
      
    }
    // This is for the mixing length model (WS-2)
    else {
      // calculating Mixing length
      RealType mixing_length = 0;
      mixing_length          = std::min(
        parameters_.turbulence.kappa * flowField.getWallDistance().getScalar(i, j),
        0.09 * flowField.getBoundaryLayerThickness().getScalar(i, j)
      );
      // calculating shear stress term Sij*Sij
      RealType Sij_e2, S11, S22, S12;

      S11 = (flowField.getVelocity().getVector(i + 1, j)[0] - flowField.getVelocity().getVector(i, j)[0])
            / (0.5 * (parameters_.meshsize->getDx(i + 1, j) + parameters_.meshsize->getDx(i, j)));

      S22 = (flowField.getVelocity().getVector(i, j + 1)[1] - flowField.getVelocity().getVector(i, j)[1])
            / (0.5 * (parameters_.meshsize->getDy(i, j + 1) + parameters_.meshsize->getDy(i, j)));

      RealType u_avg_jm1, u_avg_j, u_avg_jp1, dudy_avg;
      u_avg_jm1 = (flowField.getVelocity().getVector(i, j - 1)[0] + flowField.getVelocity().getVector(i - 1, j - 1)[0])
                  / 2;
      u_avg_j   = (flowField.getVelocity().getVector(i, j)[0] + flowField.getVelocity().getVector(i - 1, j)[0]) / 2;
      u_avg_jp1 = (flowField.getVelocity().getVector(i, j + 1)[0] + flowField.getVelocity().getVector(i - 1, j + 1)[0])
                  / 2;

      dudy_avg =  0.5 * (((u_avg_jp1-u_avg_j)/((parameters_.meshsize->getDy(i, j+1)/2)+(parameters_.meshsize->getDy(i, j)/2)))+((u_avg_j - u_avg_jm1 )/((parameters_.meshsize->getDy(i, j)/2) + (parameters_.meshsize->getDy(i, j-1)/2))));

      RealType v_avg_im1, v_avg_i, v_avg_ip1, dvdx_avg;
      v_avg_im1 = (flowField.getVelocity().getVector(i - 1, j)[1] + flowField.getVelocity().getVector(i - 1, j - 1)[1])
                  / 2;
      v_avg_i   = (flowField.getVelocity().getVector(i, j)[1] + flowField.getVelocity().getVector(i, j - 1)[1]) / 2;
      v_avg_ip1 = (flowField.getVelocity().getVector(i + 1, j)[1] + flowField.getVelocity().getVector(i + 1, j - 1)[1])
                  / 2;

      dvdx_avg =   0.5 * (((v_avg_ip1-v_avg_i)/((parameters_.meshsize->getDx(i+1, j)/2)+(parameters_.meshsize->getDx(i, j)/2)))+((v_avg_i - v_avg_im1 )/((parameters_.meshsize->getDx(i, j)/2) +(parameters_.meshsize->getDx(i-1, j)/2) )));

      S12 = 0.5 * (dudy_avg + dvdx_avg);

      Sij_e2 = ((S11 * S11) + (S22 * S22)) + (2 * (S12 * S12));

      flowField.getTurbulentViscosity().getScalar(i, j) = mixing_length * mixing_length * sqrt(2 * Sij_e2);
    }
  }
}

//***************************************************************
//  FOR 3D
//***************************************************************
void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  // Do it for fluid cells only
  if ((obstacle && OBSTACLE_SELF) == 0) {
          
      if (method_ == "turbulence-sa") {
       /*RealType omega_12, omega_13, omega_23;

      // CHECK THESE
      omega_12 = 0.5*(((flowField.getVelocity().getVector(i, j+1,k)[0] - flowField.getVelocity().getVector(i, j,k)[0])/(0.5*(parameters_.meshsize->getDy(i,j+1,k)+parameters_.meshsize->getDy(i,j,k)))) - ((flowField.getVelocity().getVector(i+1, j,k)[1] - flowField.getVelocity().getVector(i, j,k)[1])/(0.5*(parameters_.meshsize->getDx(i+1,j,k)+parameters_.meshsize->getDx(i,j,k)))));
      omega_13 = 0.5*(((flowField.getVelocity().getVector(i, j,k+1)[0] - flowField.getVelocity().getVector(i, j,k)[0])/(0.5*(parameters_.meshsize->getDz(i,j, k+1)+parameters_.meshsize->getDz(i,j,k)))) - ((flowField.getVelocity().getVector(i+1, j,k)[2] - flowField.getVelocity().getVector(i, j,k)[2])/(0.5*(parameters_.meshsize->getDx(i+1,j,k)+parameters_.meshsize->getDx(i,j,k)))));
      omega_23 = 0.5*(((flowField.getVelocity().getVector(i, j,k+1)[1] - flowField.getVelocity().getVector(i, j,k)[1])/(0.5*(parameters_.meshsize->getDz(i,j, k+1)+parameters_.meshsize->getDz(i,j,k)))) - ((flowField.getVelocity().getVector(i, j+1,k)[2] - flowField.getVelocity().getVector(i, j,k)[2])/(0.5*(parameters_.meshsize->getDy(i,j+1,k)+parameters_.meshsize->getDy(i,j,k)))));

      // RealType S = 2.0 * std::sqrt((omega_12 * omega_12 + omega_13 * omega_13 + omega_23 * omega_23));

      RealType chi = flowField.getTurbulentViscosityTransport().getScalar(i, j, k) * parameters_.flow.Re;

      RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

      RealType f_v1 = std::pow(chi, 3) / (std::pow(chi, 3) + std::pow(7.1, 3));

      RealType f_v2 = 1 - (chi / (1 + chi * f_v1));

      RealType S_hat = 2.0 * std::sqrt((omega_12 * omega_12 + omega_13 * omega_13 + omega_23 * omega_23)) + (flowField.getTurbulentViscosityTransport().getScalar(i,j,k)*f_v2)/((parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j,k)*flowField.getWallDistance().getScalar(i,j,k)));

      RealType r = std::min(RealType(10.0),(flowField.getTurbulentViscosityTransport().getScalar(i,j,k)*f_v2)/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j,k)*flowField.getWallDistance().getScalar(i,j,k))));

      RealType g = r + 0.3 * (std::pow(r, 6) - r);

      RealType f_w = g * pow((1 + std::pow(2, 6)) / (std::pow(g, 6) + std::pow(2, 6)), 1 / 6); */

    }

    // This is for the mixing length model (WS-2)
    else {
      // calculating Mixing length
      RealType mixing_length = 0;
      mixing_length          = std::min(
        parameters_.turbulence.kappa * flowField.getWallDistance().getScalar(i, j, k),
        0.09 * flowField.getBoundaryLayerThickness().getScalar(i, j, k)
      );

      // calculating shear stress term Sij*Sij
      RealType Sij_e2, S11, S22, S33, S12, S13, S23;

      S11 = (flowField.getVelocity().getVector(i + 1, j, k)[0] - flowField.getVelocity().getVector(i, j, k)[0])
            / (parameters_.meshsize->getDx(i, j, k));
      S22 = (flowField.getVelocity().getVector(i, j + 1, k)[1] - flowField.getVelocity().getVector(i, j, k)[1])
            / (parameters_.meshsize->getDy(i, j, k));
      S33 = (flowField.getVelocity().getVector(i, j, k + 1)[2] - flowField.getVelocity().getVector(i, j, k)[2])
            / (parameters_.meshsize->getDz(i, j, k));

      RealType u_avg_jm1, u_avg_j, u_avg_jp1, dudy_avg;
      u_avg_jm1 = (flowField.getVelocity().getVector(i, j - 1, k)[0]
                   + flowField.getVelocity().getVector(i - 1, j - 1, k)[0])
                  / 2;
      u_avg_j = (flowField.getVelocity().getVector(i, j, k)[0] + flowField.getVelocity().getVector(i - 1, j, k)[0]) / 2;
      u_avg_jp1 = (flowField.getVelocity().getVector(i, j + 1, k)[0]
                   + flowField.getVelocity().getVector(i - 1, j + 1, k)[0])
                  / 2;

      dudy_avg =  0.5 * (((u_avg_jp1-u_avg_j)/parameters_.meshsize->getDy(i, j,k))+((u_avg_j - u_avg_jm1 )/parameters_.meshsize->getDy(i, j,k)));

      RealType u_avg_km1, u_avg_k, u_avg_kp1, dudz_avg;
      u_avg_km1 = (flowField.getVelocity().getVector(i, j, k - 1)[0]
                   + flowField.getVelocity().getVector(i - 1, j, k - 1)[0])
                  / 2;
      u_avg_k = (flowField.getVelocity().getVector(i, j, k)[0] + flowField.getVelocity().getVector(i - 1, j, k)[0]) / 2;
      u_avg_kp1 = (flowField.getVelocity().getVector(i, j, k + 1)[0]
                   + flowField.getVelocity().getVector(i - 1, j, k + 1)[0])
                  / 2;

      dudz_avg =  0.5 * (((u_avg_kp1-u_avg_k)/parameters_.meshsize->getDz(i, j,k))+((u_avg_k - u_avg_km1 )/parameters_.meshsize->getDz(i, j,k)));

      RealType v_avg_im1, v_avg_i, v_avg_ip1, dvdx_avg;
      v_avg_im1 = (flowField.getVelocity().getVector(i - 1, j, k)[1]
                   + flowField.getVelocity().getVector(i - 1, j - 1, k)[1])
                  / 2;
      v_avg_i = (flowField.getVelocity().getVector(i, j, k)[1] + flowField.getVelocity().getVector(i, j - 1, k)[1]) / 2;
      v_avg_ip1 = (flowField.getVelocity().getVector(i + 1, j, k)[1]
                   + flowField.getVelocity().getVector(i + 1, j - 1, k)[1])
                  / 2;

      dvdx_avg =   0.5 * (((v_avg_ip1-v_avg_i)/parameters_.meshsize->getDx(i, j,k))+((v_avg_i - v_avg_im1 )/parameters_.meshsize->getDx(i, j,k)));

      RealType v_avg_km1, v_avg_k, v_avg_kp1, dvdz_avg;
      v_avg_km1 = (flowField.getVelocity().getVector(i, j, k - 1)[1]
                   + flowField.getVelocity().getVector(i, j - 1, k - 1)[1])
                  / 2;
      v_avg_k = (flowField.getVelocity().getVector(i, j, k)[1] + flowField.getVelocity().getVector(i, j - 1, k)[1]) / 2;
      v_avg_kp1 = (flowField.getVelocity().getVector(i, j, k + 1)[1]
                   + flowField.getVelocity().getVector(i, j - 1, k + 1)[1])
                  / 2;

      dvdz_avg =   0.5 * (((v_avg_kp1-v_avg_k)/parameters_.meshsize->getDz(i, j,k))+((v_avg_k - v_avg_km1 )/parameters_.meshsize->getDz(i, j,k)));

      RealType w_avg_im1, w_avg_i, w_avg_ip1, dwdx_avg;
      w_avg_im1 = (flowField.getVelocity().getVector(i - 1, j, k)[2]
                   + flowField.getVelocity().getVector(i - 1, j, k - 1)[2])
                  / 2;
      w_avg_i = (flowField.getVelocity().getVector(i, j, k)[2] + flowField.getVelocity().getVector(i, j, k - 1)[2]) / 2;
      w_avg_ip1 = (flowField.getVelocity().getVector(i + 1, j, k)[2]
                   + flowField.getVelocity().getVector(i + 1, j, k - 1)[2])
                  / 2;

      dwdx_avg =   0.5 * (((w_avg_ip1-w_avg_i)/parameters_.meshsize->getDx(i, j,k))+((w_avg_i - w_avg_im1 )/parameters_.meshsize->getDx(i, j,k)));

      RealType w_avg_jm1, w_avg_j, w_avg_jp1, dwdy_avg;
      w_avg_jm1 = (flowField.getVelocity().getVector(i, j - 1, k)[2]
                   + flowField.getVelocity().getVector(i, j - 1, k - 1)[2])
                  / 2;
      w_avg_j = (flowField.getVelocity().getVector(i, j, k)[2] + flowField.getVelocity().getVector(i, j, k - 1)[2]) / 2;
      w_avg_jp1 = (flowField.getVelocity().getVector(i, j + 1, k)[2]
                   + flowField.getVelocity().getVector(i, j + 1, k - 1)[2])
                  / 2;

      dwdy_avg =  0.5 * (((w_avg_jp1-w_avg_j)/parameters_.meshsize->getDy(i, j,k))+((w_avg_j - w_avg_jm1 )/parameters_.meshsize->getDy(i, j,k)));

      S12 = 0.5 * (dudy_avg + dvdx_avg);
      S13 = 0.5 * (dudz_avg + dwdx_avg);
      S23 = 0.5 * (dvdz_avg + dwdy_avg);

      Sij_e2 = (S11 * S11 + S22 * S22 + S33 * S33) + 2 * (S12 * S12 + S13 * S13 + S23 * S23);

      flowField.getTurbulentViscosity().getScalar(i, j, k) = mixing_length * mixing_length * sqrt(2 * Sij_e2);
    }
  }
}