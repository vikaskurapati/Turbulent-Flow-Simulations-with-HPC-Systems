#include "StdAfx.hpp"

#include "TurbulentViscosityStencil.hpp"
#include "TurbulentFlowField.hpp"

Stencils::TurbulentViscosityStencil::TurbulentViscosityStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters) {}

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j) {
  const int obstacle = flowField.getFlags().getValue(i, j);
  // Do it for fluid cells only
  if ((obstacle && OBSTACLE_SELF) == 0) {
    // calculating Mixing length
    RealType mixing_length = 0;
    mixing_length  = std::min(parameters_.turbulence.kappa * flowField.getWallDistance().getScalar(i, j), 0.09 * flowField.getBoundaryLayerThickness().getScalar(i, j));

      // calculating shear stress term Sij*Sij
    RealType Sij_e2 ,S11, S22, S12;
    
    // S11 = (flowField.getVelocity().getVector(i + 1, j)[0] - flowField.getVelocity().getVector(i - 1, j)[0])
    //       / (2 * parameters_.meshsize->getDx(i, j));
    // S22 = (flowField.getVelocity().getVector(i, j + 1)[1] - flowField.getVelocity().getVector(i, j - 1)[1])
    //       / (2 * parameters_.meshsize->getDy(i, j));

    // S12 = 0.5*( (( flowField.getVelocity().getVector( i, j + 1)[0]- flowField.getVelocity().getVector( i, j - 1)[0] )
    //                 / ( 2 * parameters_.meshsize->getDy(i, j))) 
    //             + (( flowField.getVelocity().getVector( i + 1, j)[1]- flowField.getVelocity().getVector( i - 1, j)[1] )
    //                 / ( 2 * parameters_.meshsize->getDx(i, j))) ) ;
    

    S11 = (flowField.getVelocity().getVector(i +1, j)[0] - flowField.getVelocity().getVector(i , j)[0])
          / (parameters_.meshsize->getDx(i, j));

    S22 = (flowField.getVelocity().getVector(i , j+1)[1] - flowField.getVelocity().getVector(i , j)[1])
          / (parameters_.meshsize->getDy(i, j));

    // S12= 0.5*(
    //      0.25*((flowField.getVelocity().getVector(i,j+1)[0] - flowField.getVelocity().getVector(i,j)[0] )/(parameters_.meshsize->getDy(i, j))+
    //      (flowField.getVelocity().getVector(i,j)[0] - flowField.getVelocity().getVector(i,j-1)[0] )/(parameters_.meshsize->getDy(i, j))+
    //      (flowField.getVelocity().getVector(i-1,j)[0] - flowField.getVelocity().getVector(i-1,j-1)[0] )/(parameters_.meshsize->getDy(i, j))+
    //      (flowField.getVelocity().getVector(i-1,j+1)[0] - flowField.getVelocity().getVector(i-1,j)[0] )/(parameters_.meshsize->getDy(i, j)) ) 
    //      +
    //      0.25*((flowField.getVelocity().getVector(i+1,j)[1] - flowField.getVelocity().getVector(i,j)[1] )/(parameters_.meshsize->getDx(i, j))+
    //      (flowField.getVelocity().getVector(i,j)[1] - flowField.getVelocity().getVector(i-1,j)[1] )/(parameters_.meshsize->getDx(i, j))+
    //      (flowField.getVelocity().getVector(i+1,j-1)[1] - flowField.getVelocity().getVector(i,j-1)[1] )/(parameters_.meshsize->getDx(i, j))+
    //      (flowField.getVelocity().getVector(i,j-1)[1] - flowField.getVelocity().getVector(i-1,j-1)[1] )/(parameters_.meshsize->getDx(i, j)) ) 
    //      ) ;
        RealType u_avg_jm1, u_avg_j, u_avg_jp1, dudy_avg;
         u_avg_jm1 =  (flowField.getVelocity().getVector(i,j-1)[0] - flowField.getVelocity().getVector(i-1,j-1)[0] )/ 2;
         u_avg_j =(flowField.getVelocity().getVector(i,j)[0] - flowField.getVelocity().getVector(i-1,j)[0] )/ 2;
         u_avg_jp1 = (flowField.getVelocity().getVector(i,j+1)[0] - flowField.getVelocity().getVector(i-1,j+1)[0] )/ 2;

         dudy_avg =  0.5 * (((u_avg_jp1-u_avg_j)/parameters_.meshsize->getDy(i, j))+((u_avg_j - u_avg_jm1 )/parameters_.meshsize->getDy(i, j)));

        RealType v_avg_im1, v_avg_i, v_avg_ip1, dvdx_avg;
         v_avg_im1 = (flowField.getVelocity().getVector(i-1,j)[1] - flowField.getVelocity().getVector(i-1,j-1)[1] )/ 2 ; //At i-1
         v_avg_i = (flowField.getVelocity().getVector(i,j)[1] - flowField.getVelocity().getVector(i,j-1)[1] )/ 2 ; 
         v_avg_ip1 =(flowField.getVelocity().getVector(i+1,j)[1] - flowField.getVelocity().getVector(i+1,j-1)[1] )/ 2 ;

         dvdx_avg =   0.5 * (((v_avg_ip1-v_avg_i)/parameters_.meshsize->getDx(i, j))+((v_avg_i - v_avg_im1 )/parameters_.meshsize->getDx(i, j)));

         S12 = 0.5*(dudy_avg + dvdx_avg);

    Sij_e2 = ((S11 * S11) + (S22 * S22)) + (2 * ( S12 * S12));
  
    flowField.getTurbulentViscosity().getScalar( i, j) = mixing_length * mixing_length * sqrt( 2 * Sij_e2 );
  }
}

void Stencils::TurbulentViscosityStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  const int obstacle = flowField.getFlags().getValue(i, j, k);
  // Do it for fluid cells only
  if ((obstacle && OBSTACLE_SELF) == 0) {
    // calculating Mixing length
    RealType mixing_length = 0;
    mixing_length  = std::min(parameters_.turbulence.kappa * flowField.getWallDistance().getScalar(i, j, k), 0.09 * flowField.getBoundaryLayerThickness().getScalar(i, j, k));

      // calculating shear stress term Sij*Sij
    RealType Sij_e2 ,S11, S22, S33, S12, S13, S23;

    S11 = (flowField.getVelocity().getVector(i + 1, j, k)[0] - flowField.getVelocity().getVector(i - 1, j, k)[0])
          / (2 * parameters_.meshsize->getDx(i, j, k));
    S22 = (flowField.getVelocity().getVector(i, j + 1, k)[1] - flowField.getVelocity().getVector(i, j - 1, k)[1])
          / (2 * parameters_.meshsize->getDy(i, j, k));
    S33 = (flowField.getVelocity().getVector(i, j, k + 1)[2] - flowField.getVelocity().getVector(i, j, k - 1)[2])
          / (2 * parameters_.meshsize->getDz(i, j, k));

    S12 = 0.5*( (( flowField.getVelocity().getVector( i, j + 1, k )[0]- flowField.getVelocity().getVector( i, j - 1, k )[0] )
                    / ( 2 * parameters_.meshsize->getDy(i, j, k))) 
                + (( flowField.getVelocity().getVector( i + 1, j, k )[1]- flowField.getVelocity().getVector( i - 1, j, k )[1] )
                    / ( 2 * parameters_.meshsize->getDx(i, j, k))) ) ;

    S13 = 0.5* ( (( flowField.getVelocity().getVector( i, j, k + 1 )[0]- flowField.getVelocity().getVector( i, j, k - 1 )[0] )
                    / ( 2 *  parameters_.meshsize->getDz(i, j, k) ))
                + (( flowField.getVelocity().getVector( i + 1, j, k )[2]- flowField.getVelocity().getVector( i - 1, j, k )[2] )
                    / ( 2 *  parameters_.meshsize->getDx(i, j, k))) );

    S23 = 0.5
            * ( (( flowField.getVelocity().getVector( i, j, k + 1 )[1]- flowField.getVelocity().getVector( i, j, k - 1 )[1] )
                / ( 2 * parameters_.meshsize->getDz(i, j, k) ))
            +   (( flowField.getVelocity().getVector( i, j + 1, k )[2]- flowField.getVelocity().getVector( i, j - 1, k )[2] )
                / ( 2 * parameters_.meshsize->getDy(i, j, k))) );
    Sij_e2 = (S11 * S11 + S22 * S22 + S33 * S33) + 2 * ( S12 * S12 + S13 * S13 + S23 * S23 );
  
    flowField.getTurbulentViscosity().getScalar( i, j, k ) = mixing_length * mixing_length * sqrt( 2 * Sij_e2 );
  
  }
}