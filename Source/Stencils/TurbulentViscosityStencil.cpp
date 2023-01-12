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
  // Do it for fluid cells only
  if ((obstacle && OBSTACLE_SELF) == 0) {

    if (method_ == "turbulence-sa") {

      RealType omega_12 = 0.5*(((flowField.getVelocity().getVector(i, j+1)[0] - flowField.getVelocity().getVector(i, j)[0])/(0.5*(parameters_.meshsize->getDy(i,j+1)+parameters_.meshsize->getDy(i,j)))) - ((flowField.getVelocity().getVector(i+1, j)[1] - flowField.getVelocity().getVector(i, j)[1])/(0.5*(parameters_.meshsize->getDx(i+1,j)+parameters_.meshsize->getDx(i,j)))));

      RealType chi = flowField.getTurbulentViscosityTransport().getScalar(i, j) * parameters_.flow.Re;

      RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

      RealType f_v1 = std::pow(chi, 3) / (std::pow(chi, 3) + std::pow(7.1, 3));

      RealType f_v2 = 1 - (chi / (1 + chi * f_v1));

      RealType S_hat = 2.0 * std::sqrt((omega_12 * omega_12)) + (flowField.getTurbulentViscosityTransport().getScalar(i,j)*f_v2)/((parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j)*flowField.getWallDistance().getScalar(i,j)));

      RealType r = std::min(10,(flowField.getTurbulentViscosityTransport().getScalar(i,j)*fv2)/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j)*flowField.getWallDistance().getScalar(i,j))));

      RealType g = r + 0.3 * (std::pow(r, 6) - r);

      RealType f_w = g * pow((1 + std::pow(2, 6)) / (std::pow(g, 6) + std::pow(2, 6)), 1 / 6);

      

      // RealType trip = flowField.getVelocity().getVector(i, j)[0] * flowField.getVelocity().getVector(i, j)[0]
      //                 + flowField.getVelocity().getVector(i, j)[1] * flowField.getVelocity().getVector(i, j)[1];

      // RealType omega_12_wall;
      // RealType omega_t, gt;

      // if (j < parameters_.geometry.sizeY / 2) // need to think of a better one for BFS
      // {
      //   omega_12_wall = (((flowField.getVelocity().getVector(i, 1)[0] - flowField.getVelocity().getVector(i,
      //   0)[0])/(0.5*(parameters_.meshsize->getDy(i,1)+parameters_.meshsize->getDy(i,0)))) -
      //   ((flowField.getVelocity().getVector(i+1, j)[1] - flowField.getVelocity().getVector(i,
      //   j)[1])/(0.5*(parameters_.meshsize->getDy(i+1,j)+parameters_.meshsize->getDy(i,j))))); omega_t =
      //   omega_12_wall; // be careful for 3d, it is summation of three squares and a square root there
      // } else {
      //   omega_12_wall = (((flowField.getVelocity().getVector(i, parameters_.geometry.sizeY-1)[0] -
      //   flowField.getVelocity().getVector(i,
      //   parameters_.geometry.sizeY-2)[0])/(0.5*(parameters_.meshsize->getDy(i,parameters_.geometry.sizeY-1)+parameters_.meshsize->getDy(i,parameters_.geometry.sizeY-2))))
      //   - ((flowField.getVelocity().getVector(i+1, parameters_.geometry.sizeY-1)[1] -
      //   flowField.getVelocity().getVector(i,
      //   parameters_.geometry.sizeY-1)[1])/(0.5*(parameters_.meshsize->getDy(i+1,parameters_.geometry.sizeY-1)+parameters_.meshsize->getDy(i,parameters_.geometry.sizeY-1)))));
      //   omega_t = omega_12_wall; // be careful for 3d, it is summation of three squares and a square root there
      // }
      // gt = std::min(0.1, std::min(flowField.getVelocity().getVector(i, j)[0] / (omega_t *
      // parameters_.meshsize->getDx(i, 0)), flowField.getVelocity().getVector(i,
      // j)[1]/(omega_t*parameters_.meshsize->getDy(i,0))));

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

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  // Do it for fluid cells only
  if ((obstacle && OBSTACLE_SELF) == 0) {
    if (method_ == "turbulence-sa") {
      RealType omega_12, omega_13, omega_23;

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

      RealType r = std::min(10,(flowField.getTurbulentViscosityTransport().getScalar(i,j,k)*f_v2)/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j,k)*flowField.getWallDistance().getScalar(i,j,k))));

      RealType g = r + 0.3 * (std::pow(r, 6) - r);

      RealType f_w = g * pow((1 + std::pow(2, 6)) / (std::pow(g, 6) + std::pow(2, 6)), 1 / 6);

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