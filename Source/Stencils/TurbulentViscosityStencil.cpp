#include "StdAfx.hpp"

#include "TurbulentViscosityStencil.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <petscvec.h>

#include "BoundaryType.hpp"
#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/StencilFunctions.hpp"

Stencils::TurbulentViscosityStencil::TurbulentViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  method_(parameters.simulation.type) {}

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

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

      RealType S_hat = std::max(0.3 * S, S + S_bar);

      RealType r = std::min(10.0, temp3 / (S_hat + 1e-6));

      RealType g = r + 0.3 * (std::pow(r, 6.0) - r);

      RealType cw3_pow6 = 64.0;

      RealType f_w = g * std::pow(((1 + cw3_pow6) / (std::pow(g, 6.0) + cw3_pow6)), 1.0 / 6.0);

      // ********************************************************************************
      // Term 1 (ADVECTION)
      // *********************************************************************************
      // Delete either of the below two variables
      RealType dx0 = parameters_.meshsize->getDx(i, j);

      RealType term1
        = 0.5 * (1 / dx0) * (flowField.getVelocity().getVector(i, j)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j));

      term1 -=  0.5*(1 / dx0) * (flowField.getVelocity().getVector(i-1, j)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i-1, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dx0) * std::fabs(flowField.getVelocity().getVector(i, j)[0])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) - flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dx0) * std::fabs(flowField.getVelocity().getVector(i - 1, j)[0])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j)
               ));

      RealType dy0 = parameters_.meshsize->getDy(i, j);

      term1
        += 0.5 * (1 / dy0) * (flowField.getVelocity().getVector(i, j)[1])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1));

      term1 = term1 - 0.5*(1 / dy0) * (flowField.getVelocity().getVector(i, j-1)[1])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dy0) * std::fabs(flowField.getVelocity().getVector(i, j)[1])
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
        = (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j) / (flowField.getWallDistance().getScalar(i, j) + 1e-6));

      RealType term3 = (C_w1 * f_w - (0.1355 * f_t2 / (parameters_.turbulence.kappa * parameters_.turbulence.kappa)))
                       * temp4 * temp4;
      //************************************************************
      // Term 4 (DIFFUSION Term)
      //************************************************************

      RealType dx          = parameters_.meshsize->getDx(i, j);
      RealType dx_iplus1j  = parameters_.meshsize->getDx(i + 1, j);
      RealType dx_iminus1j = parameters_.meshsize->getDx(i - 1, j);

      RealType dy          = parameters_.meshsize->getDy(i, j);
      RealType dy_ijplus1  = parameters_.meshsize->getDy(i, j + 1);
      RealType dy_ijminus1 = parameters_.meshsize->getDy(i, j - 1);

      RealType nu_ij       = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j);
      RealType nu_iplus1j  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j);
      RealType nu_iminus1j = flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j);
      RealType nu_ijplus1  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1);
      RealType nu_ijminus1 = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j - 1);
      RealType nu          = 1 / parameters_.flow.Re;

      RealType viscosity_laplacian_x = (nu + 0.5 * (nu_ij + nu_iplus1j))
                                       * ((nu_iplus1j - nu_ij) / (0.5 * (dx_iplus1j + dx)));
      viscosity_laplacian_x -= (nu + 0.5 * (nu_iminus1j + nu_ij)) * ((nu_ij - nu_iminus1j) / (0.5 * (dx_iminus1j + dx)));

      viscosity_laplacian_x = viscosity_laplacian_x / dx;

      RealType viscosity_laplacian_y = (nu + 0.5 * (nu_ij + nu_ijplus1))
                                       * ((nu_ijplus1 - nu_ij) / (0.5 * (dy_ijplus1 + dy)));
      viscosity_laplacian_y -= (nu + 0.5 * (nu_ijminus1 + nu_ij)) * ((nu_ij - nu_ijminus1) / (0.5 * (dy_ijminus1 + dy)));

      viscosity_laplacian_y = viscosity_laplacian_y / dy;

      RealType viscosity_laplacian = viscosity_laplacian_x + viscosity_laplacian_y;

      RealType dnudx = (nu_iplus1j - nu_iminus1j) / ((0.5 * dx_iplus1j) + dx + (0.5 * dx_iminus1j));
      RealType dnudy = (nu_ijplus1 - nu_ijminus1) / ((0.5 * dy_ijplus1) + dy + (0.5 * dy_ijminus1));

      RealType advec_term = 0.622 * ((dnudx * dnudx) + (dnudy * dnudy));

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

      // Boundary condition for ONLY channel flow
      if ((parameters_.bfStep.xRatio * parameters_.geometry.sizeX <= 0) && (parameters_.bfStep.yRatio * parameters_.geometry.sizeY <= 0)) {
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
        // if the current cell is an obstacle
        if ((obstacle & OBSTACLE_SELF) == 1) {
          // If top cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_TOP) == 0) {
            const RealType dy_t         = parameters_.meshsize->getDy(i, j + 1);
            const RealType dy_t1        = parameters_.meshsize->getDy(i, j);
            velocity.getVector(i, j)[0] = -dy_t1 / dy_t * velocity.getVector(i, j + 1)[0];
          }
          // Same for bottom
          if ((obstacle & OBSTACLE_BOTTOM) == 0) {
            const RealType dy_b         = parameters_.meshsize->getDy(i, j - 1);
            const RealType dy_b1        = parameters_.meshsize->getDy(i, j);
            velocity.getVector(i, j)[0] = -dy_b1 / dy_b * velocity.getVector(i, j - 1)[0];
          }
          // If right cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_RIGHT) == 0) {
            const RealType dx_r         = parameters_.meshsize->getDx(i + 1, j);
            const RealType dx_r1        = parameters_.meshsize->getDx(i, j);
            velocity.getVector(i, j)[1] = -dx_r1 / dx_r * velocity.getVector(i + 1, j)[1];
          }
          // Same for left
          if ((obstacle & OBSTACLE_LEFT) == 0) {
            const RealType dx_l         = parameters_.meshsize->getDx(i - 1, j);
            const RealType dx_l1        = parameters_.meshsize->getDx(i, j);
            velocity.getVector(i, j)[1] = -dx_l1 / dx_l * velocity.getVector(i - 1, j)[1];
          }

          // Set normal velocity to zero if right neighbour is not obstacle
          if ((obstacle & OBSTACLE_RIGHT) == 0) {
            velocity.getVector(i, j)[0] = 0.0;
          }

          // Set normal velocity to zero if top neighbour is not obstacle
          if ((obstacle & OBSTACLE_TOP) == 0) {
            velocity.getVector(i, j)[1] = 0.0;
          }
        }

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
            } else {
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
        if ((i + parameters_.parallel.firstCorner[0] == 2) && (j + parameters_.parallel.firstCorner[1] > 1) && ((obstacle & OBSTACLE_SELF) == 0)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j) = 3.0 / parameters_.flow.Re;
        }

        // at the outlet of the BFS channel
        if ((i + parameters_.parallel.firstCorner[0] == parameters_.geometry.sizeX + 1) && (j + parameters_.parallel.firstCorner[1] > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i + 1, j
          ) = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j);
        }
      }

      chi  = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j) * parameters_.flow.Re;
      f_v1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(7.1, 3.0));
      // calculating eddy viscosity from nu_tilda
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
  if ((i >= 2) && (j >= 2) && (k>=2) && i < (parameters_.geometry.sizeX + 2) && j < (parameters_.geometry.sizeY + 2) && k< (parameters_.geometry.sizeZ +2)) {

    if (method_ == "turbulence-sa") {

      RealType dudy = (flowField.getVelocity().getVector(i, j + 1,k)[0] - flowField.getVelocity().getVector(i, j - 1,k)[0])
                      / (0.5*parameters_.meshsize->getDy(i, j-1,k) + parameters_.meshsize->getDy(i,j,k) + 0.5*parameters_.meshsize->getDy(i,j+1,k));

      RealType dvdx = (flowField.getVelocity().getVector(i + 1, j,k)[1] - flowField.getVelocity().getVector(i - 1, j,k)[1])
                      / (0.5*parameters_.meshsize->getDx(i-1, j,k) + parameters_.meshsize->getDx(i,j,k) + 0.5*parameters_.meshsize->getDx(i+1,j,k));

      RealType W_12 = 0.5 * (dudy - dvdx);

      RealType dwdx = (flowField.getVelocity().getVector(i+1, j ,k)[2] - flowField.getVelocity().getVector(i-1, j ,k)[2])
                      / (0.5*parameters_.meshsize->getDx(i-1, j,k) + parameters_.meshsize->getDx(i,j,k) + 0.5*parameters_.meshsize->getDx(i+1,j,k));

      RealType dudz =  (flowField.getVelocity().getVector(i, j ,k+1)[0] - flowField.getVelocity().getVector(i, j ,k-1)[0])
                      / (0.5*parameters_.meshsize->getDz(i, j,k-1) + parameters_.meshsize->getDz(i,j,k) + 0.5*parameters_.meshsize->getDz(i,j,k+1));

      RealType W_13 = 0.5 * (dwdx - dudz);

      RealType dwdy = (flowField.getVelocity().getVector(i, j + 1,k)[2] - flowField.getVelocity().getVector(i, j - 1,k)[2])
                      / (0.5*parameters_.meshsize->getDy(i, j-1,k) + parameters_.meshsize->getDy(i,j,k) + 0.5*parameters_.meshsize->getDy(i,j+1,k));

      RealType dvdz = (flowField.getVelocity().getVector(i, j ,k+1)[1] - flowField.getVelocity().getVector(i, j ,k-1)[1])
                      / (0.5*parameters_.meshsize->getDz(i, j,k-1) + parameters_.meshsize->getDz(i,j,k) + 0.5*parameters_.meshsize->getDz(i,j,k+1));

      RealType W_23 = 0.5 * (dwdy - dvdz);

      RealType chi = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k) / (1 / parameters_.flow.Re);

      RealType f_t2 = 1.2 * std::exp(-0.5 * chi * chi);

      RealType temp1 = std::pow(7.1, 3.0); // c_v1 = 7.1
      RealType temp2 = std::pow(chi, 3.0);

      RealType f_v1 = temp2 / (temp2 + temp1);

      RealType f_v2 = 1.0 - (chi / (1 + (chi * f_v1)));

      RealType temp3 = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k)/(parameters_.turbulence.kappa*parameters_.turbulence.kappa*(flowField.getWallDistance().getScalar(i, j,k)+1e-6)*(flowField.getWallDistance().getScalar(i, j,k)+1e-6));

      RealType S = 2.0 * std::sqrt((W_12 * W_12) + (W_13 * W_13) + (W_23 * W_23));

      RealType S_bar = temp3 * f_v2;

      RealType S_hat = std::max(0.3 * S, S + S_bar);

      RealType r = std::min(10.0, temp3 / (S_hat + 1e-6));

      RealType g = r + 0.3 * (std::pow(r, 6.0) - r);

      RealType cw3_pow6 = 64.0;

      RealType f_w = g * std::pow(((1 + cw3_pow6) / (std::pow(g, 6.0) + cw3_pow6)), 1.0 / 6.0);

      // ********************************************************************************
      // Term 1 (ADVECTION)
      // *********************************************************************************

      RealType dx = parameters_.meshsize->getDx(i, j, k);

      RealType term1
        = 0.5 * (1 / dx) * (flowField.getVelocity().getVector(i, j,k)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j,k));

      term1 -=  0.5*(1 / dx) * (flowField.getVelocity().getVector(i-1, j,k)[0])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i-1, j,k) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dx) * std::fabs(flowField.getVelocity().getVector(i - 1, j, k)[0])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j, k)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k)
               ));

      RealType dy = parameters_.meshsize->getDy(i, j);

      term1
        += 0.5 * (1 / dy) * (flowField.getVelocity().getVector(i, j,k)[1])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1,k));

      term1 = term1 - 0.5*(1 / dy) * (flowField.getVelocity().getVector(i, j-1,k)[1])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j-1,k) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dy) * std::fabs(flowField.getVelocity().getVector(i, j,k)[1])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k) - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1,k));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dy) * std::fabs(flowField.getVelocity().getVector(i, j - 1, k)[1])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j - 1, k)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k)
               ));

      RealType dz = parameters_.meshsize->getDz(i, j, k);

      term1
        += 0.5 * (1 / dz) * (flowField.getVelocity().getVector(i, j,k)[2])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j ,k+1));

      term1 = term1 - 0.5*(1 / dz) * (flowField.getVelocity().getVector(i, j,k-1)[2])
          * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k-1) + flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k));

      term1
        += 0.5 * parameters_.solver.gamma * (1 / dz) * std::fabs(flowField.getVelocity().getVector(i, j,k)[2])
           * (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j,k) - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j ,k+1));

      term1 -= 0.5 * parameters_.solver.gamma * (1 / dz) * std::fabs(flowField.getVelocity().getVector(i, j, k - 1)[2])
               * ((
                 flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k - 1)
                 - flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k)
               ));

      //**********************************************************************
      // Term 2 (PRODUCTION/SOURCE Term)
      //**********************************************************************
      RealType term2 = 0.1355 * (1 - f_t2) * S_hat
                       * flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k);
      //**********************************************************************
      // Term 3 (Wall DESTRUCTION Source Term)
      //*************************************************************************
      RealType C_w1 = 3.239067817;

      RealType temp4
        = (flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k) / (flowField.getWallDistance().getScalar(i, j, k)));

      RealType term3 = (C_w1 * f_w - (0.1355 * f_t2 / (parameters_.turbulence.kappa * parameters_.turbulence.kappa)))
                       * temp4 * temp4;
      //************************************************************
      // Term 4 (DIFFUSION Term)
      //************************************************************

      RealType dx_ijk       = parameters_.meshsize->getDx(i, j, k);
      RealType dx_iplus1jk  = parameters_.meshsize->getDx(i + 1, j, k);
      RealType dx_iminus1jk = parameters_.meshsize->getDx(i - 1, j, k);

      RealType dy_ijk       = parameters_.meshsize->getDy(i, j, k);
      RealType dy_ijplus1k  = parameters_.meshsize->getDy(i, j + 1, k);
      RealType dy_ijminus1k = parameters_.meshsize->getDy(i, j - 1, k);

      RealType dz_ijk       = parameters_.meshsize->getDz(i, j, k);
      RealType dz_ijkplus1  = parameters_.meshsize->getDz(i, j, k + 1);
      RealType dz_ijkminus1 = parameters_.meshsize->getDz(i, j, k - 1);

      RealType nu_ijk = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k);
      RealType nu     = (1 / parameters_.flow.Re);

      RealType nu_iplus1jk  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i + 1, j, k);
      RealType nu_iminus1jk = flowField.getPreviousTurbulentViscosityTransport().getScalar(i - 1, j, k);

      RealType nu_ijplus1k  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j + 1, k);
      RealType nu_ijminus1k = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j - 1, k);

      RealType nu_ijkplus1  = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k + 1);
      RealType nu_ijkminus1 = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k - 1);

      RealType viscosity_laplacian_x = (nu + 0.5 * (nu_ijk + nu_iplus1jk))
                                       * ((nu_iplus1jk - nu_ijk) / (0.5 * (dx_iplus1jk + dx_ijk)));
      viscosity_laplacian_x -= (nu + 0.5 * (nu_iminus1jk + nu_ijk))
                               * ((nu_ijk - nu_iminus1jk) / (0.5 * (dx_iminus1jk + dx_ijk)));

      viscosity_laplacian_x = viscosity_laplacian_x / dx_ijk;

      RealType viscosity_laplacian_y = (nu + 0.5 * (nu_ijk + nu_ijplus1k))
                                       * ((nu_ijplus1k - nu_ijk) / (0.5 * (dy_ijplus1k + dy_ijk)));
      viscosity_laplacian_y -= (nu + 0.5 * (nu_ijminus1k + nu_ijk))
                               * ((nu_ijk - nu_ijminus1k) / (0.5 * (dy_ijminus1k + dy_ijk)));

      viscosity_laplacian_y = viscosity_laplacian_y / dy_ijk;

      RealType viscosity_laplacian_z = (nu + 0.5 * (nu_ijk + nu_ijkplus1))
                                       * ((nu_ijkplus1 - nu_ijk) / (0.5 * (dz_ijkplus1 + dz_ijk)));
      viscosity_laplacian_z -= (nu + 0.5 * (nu_ijkminus1 + nu_ijk))
                               * ((nu_ijk - nu_ijkminus1) / (0.5 * (dz_ijkminus1 + dz_ijk)));

      viscosity_laplacian_z = viscosity_laplacian_z / dz_ijk;

      RealType viscosity_laplacian = viscosity_laplacian_x + viscosity_laplacian_y + viscosity_laplacian_z;

      RealType dnudx = (nu_iplus1jk - nu_iminus1jk) / ((0.5 * dx_iplus1jk) + dx_ijk + (0.5 * dx_iminus1jk));
      RealType dnudy = (nu_ijplus1k - nu_ijminus1k) / ((0.5 * dy_ijplus1k) + dy_ijk + (0.5 * dy_ijminus1k));
      RealType dnudz = (nu_ijkplus1 - nu_ijkminus1) / ((0.5 * dz_ijkplus1) + dz_ijk + (0.5 * dz_ijkminus1));

      RealType advec_term = 0.622 * ((dnudx * dnudx) + (dnudy * dnudy) + (dnudz * dnudz));

      RealType term4 = (1.5) * (viscosity_laplacian + advec_term);

      flowField.getCurrentTurbulentViscosityTransport().getScalar(
        i, j, k
      ) = flowField.getPreviousTurbulentViscosityTransport().getScalar(i, j, k)
          + parameters_.timestep.dt * (term2 - term3 + term4 - term1);

      flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k) = std::max(
        0.0, flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k)
      );

      //***********************************************
      // Boundary Conditions for nu_transport:
      //***********************************************
      // Boundary condition for CHANNEL FLOW 3D
      if ((parameters_.bfStep.xRatio * parameters_.geometry.sizeX <= 0) && (parameters_.bfStep.yRatio * parameters_.geometry.sizeY <= 0)  ) {

        // at the lower wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == 2) && ((k + parameters_.parallel.firstCorner[2]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j - 1, k
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

        // at the top wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == parameters_.geometry.sizeY + 1) && ((k + parameters_.parallel.firstCorner[2]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j + 1, k
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

        // at the front wall of the channel (no slip)
        if (((k + parameters_.parallel.firstCorner[2]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j, k - 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

        // at the back wall of the channel (no slip)
        if (((k + parameters_.parallel.firstCorner[2]) == parameters_.geometry.sizeZ + 1) && ((j + parameters_.parallel.firstCorner[0]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j, k + 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

        // at the inlet of the channel generating constant turbulent transport visc throughout the duration of the
        // simulation
        if (((i + parameters_.parallel.firstCorner[0]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((k + parameters_.parallel.firstCorner[2]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j, k) = 3.0 / parameters_.flow.Re;
        }

        // at the outlet of the channel
        if (((i + parameters_.parallel.firstCorner[0]) == parameters_.geometry.sizeX + 1) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((k + parameters_.parallel.firstCorner[2]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i + 1, j, k
          ) = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

      } // end of BC for 3D Channel

      // Boundary conditions for BFS
      else {
        // if the current cell is an obstacle
        if ((obstacle & OBSTACLE_SELF) == 1) {
          // If top cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_TOP) == 0) {
            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j, k
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j + 1, k);
          }
          // Same for bottom
          if ((obstacle & OBSTACLE_BOTTOM) == 0) {

            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j, k
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j - 1, k);
          }
          // If right cell is fluid, then the no-slip boundary has to be enforced
          if ((obstacle & OBSTACLE_RIGHT) == 0) {
            if ((obstacle & OBSTACLE_TOP) == 0) {
              // no need for k+1 becasue the step exists completely in the z direction, so no fluid neighbour for the
              // corner element in the z direction
              flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j,k) = -0.5*(flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j+1,k) + flowField.getCurrentTurbulentViscosityTransport().getScalar(i+1, j,k) );
            } else {
              flowField.getCurrentTurbulentViscosityTransport().getScalar(
                i, j, k
              ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i + 1, j, k);
            }
          }
          // Same for left
          if ((obstacle & OBSTACLE_LEFT) == 0) {
            flowField.getCurrentTurbulentViscosityTransport().getScalar(
              i, j, k
            ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j, k);
          }
        }
        // at the lower wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == 2) && ((k + parameters_.parallel.firstCorner[2]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j - 1, k
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }
        // at the top wall of the channel (no slip)
        if (((j + parameters_.parallel.firstCorner[1]) == parameters_.geometry.sizeY + 1) && ((k + parameters_.parallel.firstCorner[2]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j + 1, k
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }
        // at the front wall of the channel (no slip)
        if (((k + parameters_.parallel.firstCorner[2]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j, k - 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }
        // at the back wall of the channel (no slip)
        if (((k + parameters_.parallel.firstCorner[2]) == parameters_.geometry.sizeZ + 1) && ((j + parameters_.parallel.firstCorner[0]) > 1) && ((i + parameters_.parallel.firstCorner[0]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i, j, k + 1
          ) = -flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }
        // at the inlet of the channel generating constant turbulent transport visc throughout the duration of the
        // simulation
        if (((i + parameters_.parallel.firstCorner[0]) == 2) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((k + parameters_.parallel.firstCorner[2]) > 1) && ((obstacle & OBSTACLE_SELF) == 0)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(i - 1, j, k) = 3.0 / parameters_.flow.Re;
        }
        // at the outlet of the channel
        if (((i + parameters_.parallel.firstCorner[0]) == parameters_.geometry.sizeX + 1) && ((j + parameters_.parallel.firstCorner[1]) > 1) && ((k + parameters_.parallel.firstCorner[2]) > 1)) {
          flowField.getCurrentTurbulentViscosityTransport().getScalar(
            i + 1, j, k
          ) = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);
        }

      } // end of BC for BFS

      chi = flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k) * parameters_.flow.Re;

      f_v1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(7.1, 3.0));

      // calculating eddy viscosity from nu_tilda
      flowField.getTurbulentViscosity().getScalar(
        i, j, k
      ) = f_v1 * flowField.getCurrentTurbulentViscosityTransport().getScalar(i, j, k);

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