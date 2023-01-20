#include "StdAfx.hpp"

#include "TurbulentViscosityStencil.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

#include "BoundaryType.hpp"
#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/StencilFunctions.hpp"

Stencils::TurbulentViscosityStencil::TurbulentViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  method_(parameters.simulation.type) {}

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  //  Do it for fluid cells only

  // if ((obstacle && OBSTACLE_SELF) == 0) {

  // Strictly iterating only in Fluid cells

  if ((i >= 2) && (j >= 2) && i < (parameters_.geometry.sizeX + 2) && j < (parameters_.geometry.sizeY + 2)) {

    if (method_ == "turbulence-sa") {

      RealType dudy = (flowField.getVelocity().getVector(i, j + 1)[0] - flowField.getVelocity().getVector(i, j - 1)[0])
                      / (0.5*parameters_.meshsize->getDy(i, j-1) + parameters_.meshsize->getDy(i,j) + 0.5*parameters_.meshsize->getDy(i,j+1));

      RealType dvdx = (flowField.getVelocity().getVector(i + 1, j)[1] - flowField.getVelocity().getVector(i - 1, j)[1])
                      / (0.5*parameters_.meshsize->getDx(i-1, j) + parameters_.meshsize->getDx(i,j) + 0.5*parameters_.meshsize->getDx(i+1,j));

      RealType W_12 = 0.5 * (dudy - dvdx);

      RealType chi = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / (1 / parameters_.flow.Re);

      RealType c_v2 = 0.7;
      RealType c_v3 = 0.9;

      RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

      RealType temp1 = std::pow(7.1, 3.0); // c_v1 = 7.1
      RealType temp2 = std::pow(chi, 3.0);

      RealType f_v1 = temp2 / (temp2 + temp1);

      RealType f_v2 = 1.0 - (chi / (1 + (chi * f_v1)));

      RealType temp3 = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)/(parameters_.turbulence.kappa*parameters_.turbulence.kappa*(flowField.getWallDistance().getScalar(i, j)+1e-6)*(flowField.getWallDistance().getScalar(i, j)+1e-6));

      RealType S = 2.0 * std::sqrt((W_12 * W_12));

      RealType S_bar = temp3 * f_v2;

      // Three options for S_hat
      // RealType S_hat = std::max(0.0, S + S_bar);

      RealType S_hat = std::max(0.3 * S, S + S_bar);

      // RealType S_hat;

      // if (S_bar >= -c_v2 * S) {
      //   // std::cout << "first term: " << S << " second term: " << S_bar << std::endl;
      //   S_hat = S + S_bar;
      // } else {
      //   // std::cout << "first term: " << S << " second term: " << (S*(c_v2*c_v2*S + c_v3*S_bar))/((c_v3- 2.0*c_v2)*S
      //   -
      //   // S_bar) << std::endl;
      //   S_hat = S + (S * (c_v2 * c_v2 * S + c_v3 * S_bar)) / ((c_v3 - 2.0 * c_v2) * S - S_bar);
      // }

      RealType r = std::min(10.0, temp3 / (S_hat + 1e-6));

      RealType g = r + 0.3 * (std::pow(r, 6.0) - r);

      RealType cw3_pow6 = 64.0;

      RealType f_w = g * std::pow(((1 + cw3_pow6) / (std::pow(g, 6.0) + cw3_pow6)), 1.0 / 6.0);

      // ********************************************************************************
      // Term 1 (ADVECTION)
      // *********************************************************************************

      RealType dx1 = parameters_.meshsize->getDx(i, j);
      RealType dx0 = parameters_.meshsize->getDx(i, j);

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

      RealType dy1 = parameters_.meshsize->getDy(i, j);
      RealType dy0 = parameters_.meshsize->getDy(i, j);

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

      //**********************************************************************
      // Term 2 (PRODUCTION/SOURCE Term)
      //**********************************************************************
      RealType term2 = 0.1355 * (1 - f_t2) * S_hat * flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j);
      //**********************************************************************
      // Term 3 (Wall DESTRUCTION Source Term)
      //*************************************************************************
      RealType C_w1 = 3.239067817;

      RealType temp4
        = (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / (flowField.getWallDistance().getScalar(i, j)+1e-6));

      RealType term3 = (C_w1 * f_w - (0.1355 * f_t2 / (parameters_.turbulence.kappa * parameters_.turbulence.kappa)))
                       * temp4 * temp4;
      //************************************************************
      // Term 4 (DIFFUSION Term)
      //************************************************************

      RealType dx = parameters_.meshsize->getDx(i, j);
      RealType dy = parameters_.meshsize->getDy(i, j);

      RealType nu_ij       = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j);
      RealType nu_iplus1j  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j);
      RealType nu_iminus1j = flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j);
      RealType nu_ijplus1  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1);
      RealType nu_ijminus1 = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j - 1);
      RealType nu          = 1 / parameters_.flow.Re;

      RealType viscosity_laplacian_x = (nu + 0.5 * (nu_ij + nu_iplus1j)) * (nu_iplus1j - nu_ij) / dx;
      viscosity_laplacian_x -= (nu + 0.5 * (nu_iminus1j + nu_ij)) * (nu_ij - nu_iminus1j) / dx;

      viscosity_laplacian_x = viscosity_laplacian_x / dx;

      RealType viscosity_laplacian_y = (nu + 0.5 * (nu_ij + nu_ijplus1)) * (nu_ijplus1 - nu_ij) / dy;
      viscosity_laplacian_y -= (nu + 0.5 * (nu_ijminus1 + nu_ij)) * (nu_ij - nu_ijminus1) / dy;

      viscosity_laplacian_y = viscosity_laplacian_y / dy;

      RealType viscosity_laplacian = viscosity_laplacian_x + viscosity_laplacian_y;

      RealType dnudx = 0.5 * (nu_iplus1j - nu_iminus1j) / dx;
      RealType dnudy = 0.5 * (nu_ijplus1 - nu_ijminus1) / dy;

      RealType advec_term = 0.622 * (dnudx * dnudx + dnudy * dnudy);

      RealType term4 = (1.5) * (viscosity_laplacian + advec_term);

      flowField.getCurrentTurbulentViscosityTransport().getScalar(
        i, j
      ) = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
          + parameters_.timestep.dt * (term2 - term3 + term4 - term1);

      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) = std::max(
        0.0, flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j)
      );

      //***********************************************
      // Boundary Conditions for nu_transport:
      //***********************************************
      // std::cout << parameters_.bfStep.xRatio << ":   " << parameters_.bfStep.yRatio << std::endl;

      // int rrrank;
      // MPI_Comm_rank(PETSC_COMM_WORLD, &rrrank);

      // if (rrrank == 1) {
      //   std::cout << parameters_.geometry.sizeX << "  " << parameters_.geometry.sizeY << std::endl;
      // }

      // Boundary condition for channel flow
      if (parameters_.bfStep.xRatio * parameters_.geometry.sizeX <= 0 && parameters_.bfStep.xRatio * parameters_.geometry.sizeX <= 0) {

        // at the lower wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == 2) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j - 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }

        // at the top wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == parameters_.geometry.sizeY + 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j + 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }

        // at the inlet of the channel
        if (((i + parameters_.parallel.firstCorner[0]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j) = 3.0 / parameters_.flow.Re;
        }

        // //at the outlet of the channel
        if (((i + parameters_.parallel.firstCorner[0]) == parameters_.geometry.sizeX + 1) && ((j + parameters_.parallel.firstCorner[1]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i + 1, j
          ) = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
          // std::cout<<"here"<<std::endl;
        }

      }

      // Boundary condition for BFS
      else {

        if ((obstacle & OBSTACLE_SELF) == 1) {
          // If top cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_TOP) == 0) {
              
            
            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j + 1);
          }
          // Same for bottom
          if ((obstacle & OBSTACLE_BOTTOM) == 0) {
            
            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j - 1);
          }
          // If right cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_RIGHT) == 0) {

            if ((obstacle & OBSTACLE_TOP) == 0) {


              flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) = -0.5*(flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j+1) + flowField.getCurrentTurbulentViscosityTransport().getScalar(i+1, j) );
            }

            else {

              flowField.getCurrentTurbulentViscosityTransport().getScalar(
                i, j
              ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 1, j);
            }
          }
          // Same for left
          if ((obstacle & OBSTACLE_LEFT) == 0) {

            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j);
          }
        }
        // at the lower wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == 2) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j - 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }

        // at the top wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == parameters_.geometry.sizeY + 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j + 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }
        // at the inlet of the BFS channel
        if ((i + parameters_.parallel.firstCorner[0]== 2) && (j+parameters_.parallel.firstCorner[1] > 1) && ((obstacle & OBSTACLE_SELF) == 0)) {                                       
          flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j) = 3.0 / parameters_.flow.Re;
        }

        // at the outlet of the BFS channel
        if ((i+parameters_.parallel.firstCorner[0] == parameters_.geometry.sizeX + 1) && (j+parameters_.parallel.firstCorner[1] > 1)) {
                                      //std::cout<<i<<"  "<<j<<std::endl;

          flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 1, j) = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }

        // std::cout << "here " << i << "  :  " << j << " : " << parameters_.geometry.sizeX << std::endl;

        // at  the outlet of the channel

        // // Set normal velocity to zero if right neighbour is not obstacle
        // if ((obstacle & OBSTACLE_RIGHT) == 0) {
        //   velocity.getVector(i, j)[0] = 0.0;
        // }

        // // Set normal velocity to zero if top neighbour is not obstacle
        // if ((obstacle & OBSTACLE_TOP) == 0) {
        //   velocity.getVector(i, j)[1] = 0.0;
        // }
      }

      chi  = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) * parameters_.flow.Re;
      f_v1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(7.1, 3.0));
      flowField.getTurbulentViscosity().getScalar(
        i, j
      ) = f_v1 * flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
    }

    //****************************************************************
    // This is for the mixing length model (WS-2)
    //****************************************************************
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
     omega_12 = 0.5*(((flowField.getVelocity().getVector(i, j+1,k)[0] - flowField.getVelocity().getVector(i,
     j,k)[0])/(0.5*(parameters_.meshsize->getDy(i,j+1,k)+parameters_.meshsize->getDy(i,j,k)))) -
     ((flowField.getVelocity().getVector(i+1, j,k)[1] - flowField.getVelocity().getVector(i,
     j,k)[1])/(0.5*(parameters_.meshsize->getDx(i+1,j,k)+parameters_.meshsize->getDx(i,j,k))))); omega_13 =
     0.5*(((flowField.getVelocity().getVector(i, j,k+1)[0] - flowField.getVelocity().getVector(i,
     j,k)[0])/(0.5*(parameters_.meshsize->getDz(i,j, k+1)+parameters_.meshsize->getDz(i,j,k)))) -
     ((flowField.getVelocity().getVector(i+1, j,k)[2] - flowField.getVelocity().getVector(i,
     j,k)[2])/(0.5*(parameters_.meshsize->getDx(i+1,j,k)+parameters_.meshsize->getDx(i,j,k))))); omega_23 =
     0.5*(((flowField.getVelocity().getVector(i, j,k+1)[1] - flowField.getVelocity().getVector(i,
     j,k)[1])/(0.5*(parameters_.meshsize->getDz(i,j, k+1)+parameters_.meshsize->getDz(i,j,k)))) -
     ((flowField.getVelocity().getVector(i, j+1,k)[2] - flowField.getVelocity().getVector(i,
     j,k)[2])/(0.5*(parameters_.meshsize->getDy(i,j+1,k)+parameters_.meshsize->getDy(i,j,k)))));

     // RealType S = 2.0 * std::sqrt((omega_12 * omega_12 + omega_13 * omega_13 + omega_23 * omega_23));

     RealType chi = flowField.getTurbulentViscosityTransport().getScalar(i, j, k) * parameters_.flow.Re;

     RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

     RealType f_v1 = std::pow(chi, 3) / (std::pow(chi, 3) + std::pow(7.1, 3));

     RealType f_v2 = 1 - (chi / (1 + chi * f_v1));

     RealType S_hat = 2.0 * std::sqrt((omega_12 * omega_12 + omega_13 * omega_13 + omega_23 * omega_23)) +
     (flowField.getTurbulentViscosityTransport().getScalar(i,j,k)*f_v2)/((parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j,k)*flowField.getWallDistance().getScalar(i,j,k)));

     RealType r =
     std::min(RealType(10.0),(flowField.getTurbulentViscosityTransport().getScalar(i,j,k)*f_v2)/(S_hat*(parameters_.turbulence.kappa*parameters_.turbulence.kappa)*(flowField.getWallDistance().getScalar(i,j,k)*flowField.getWallDistance().getScalar(i,j,k))));

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