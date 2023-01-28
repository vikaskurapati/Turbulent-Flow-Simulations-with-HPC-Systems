#include "StdAfx.hpp"

#include "TurbulentVTKStencil.hpp"

#include "TurbulentFlowField.hpp"

Stencils::TurbulentVTKStencil::TurbulentVTKStencil(const Parameters& parameters):
  FieldStencil<TurbulentFlowField>(parameters),
  written_(false),
  prefix_(parameters.vtk.prefix) {

  if (parameters_.parallel.rank == 0) {
    const std::string outputFolder = "Output/" + prefix_;

#ifdef _MSC_VER
    const int success = _mkdir(outputFolder.data());
    if (success != 0) {
      if (errno != EEXIST) {
        spdlog::error("Cannot create folder {}", prefix_);
        throw std::runtime_error("Error while creating folder for VTK output");
      }
    }
#else
    char str_mkdir[256] = "mkdir -p ";
    strcat(str_mkdir, outputFolder.data());
    const int success = system(str_mkdir);
    if (success != 0) {
      spdlog::error("Cannot create folder {}", prefix_);
      throw std::runtime_error("Error while creating folder for VTK output");
    }
#endif
  }

  std::ofstream gitignore;
  gitignore.open(prefix_ + "/.gitignore");
  gitignore << "*";
  gitignore.close();
}

void Stencils::TurbulentVTKStencil::writeVTKHeader(std::ostream& file) const {
  file << "# vtk DataFile Version 2.0" << std::endl << "NS-EOF" << std::endl << "ASCII" << std::endl << std::endl;
}

void Stencils::TurbulentVTKStencil::writePoints(std::ostream& file, RealType simulationTime) const {
  // Number of points in every direction
  int px = parameters_.parallel.localSize[0] + 1;
  int py = parameters_.parallel.localSize[1] + 1;
  int pz = parameters_.geometry.dim == 2 ? 1 : parameters_.parallel.localSize[2] + 1;

  std::string grid;
  char        buffer[256];

  grid.reserve((file.precision() + 6) * px * py * pz * 3);

  sprintf(
    buffer,
    "DATASET STRUCTURED_GRID\nFIELD FieldData 1\nTIME 1 1 double\n%f\nDIMENSIONS %d %d %d\nPOINTS %d float\n",
    simulationTime,
    px,
    py,
    pz,
    px * py * pz
  );
  grid.append(buffer);

  if (parameters_.geometry.dim == 3) {
    for (int k = 2; k < 2 + pz; k++) {
      for (int j = 2; j < 2 + py; j++) {
        for (int i = 2; i < 2 + px; i++) {
          // Determine positions of grid points at lower/left/front corner of the respective grid cell (i,j,k) -> use
          // meshsize-ptr
          sprintf(
            buffer,
            "%f %f %f\n",
            parameters_.meshsize->getPosX(i, j, k),
            parameters_.meshsize->getPosY(i, j, k),
            parameters_.meshsize->getPosZ(i, j, k)
          );
          grid.append(buffer);
        }
      }
    }
  } else {
    for (int j = 2; j < 2 + py; j++) {
      for (int i = 2; i < 2 + px; i++) {
        sprintf(buffer, "%f %f 0.0\n", parameters_.meshsize->getPosX(i, j), parameters_.meshsize->getPosY(i, j));
        grid.append(buffer);
      }
    }
  }
  grid.append("\n");
  file << grid;
}

void Stencils::TurbulentVTKStencil::openFile(int timeStep, RealType simulationTime, std::string parameter) {
  written_ = false;
  std::stringstream namestream;
  std::string       name;
  namestream.precision(4);
  namestream
    << "Output"
    << "/" << prefix_ << "/" << prefix_ << "." << parameter << "." << parameters_.parallel.rank << "." << timeStep
    << ".vtk";
  name = namestream.str();
  ofile_.open(name.c_str());
  namestream.str("");

  writeVTKHeader(ofile_);
  writePoints(ofile_, simulationTime);
}

void Stencils::TurbulentVTKStencil::apply(TurbulentFlowField& flowField, int i, int j) {
#ifndef NDEBUG
  RealType h     = 0.0;
  RealType delta = 0.0;
#endif

  ASSERTION(FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 2);
  // std::cout<< i << "  "<<j<<std::endl;
  RealType pressure    = 0.0;
  RealType velocity[2] = {0.0, 0.0};
  RealType viscosity   = 0.0;

  //*********************************************************************
  // Shear Stress (tau_net) Calculation
  // //******************************************************************

  RealType tau_net, dudy_avg;
  RealType u_avg_im1_top, u_avg_im1_bottom, u_avg_i_top, u_avg_i_bottom, u_avg_top, u_avg_bottom;

  u_avg_im1_top = 0.5
                  * (flowField.getVelocity().getVector(i - 1, j + 1)[0] + flowField.getVelocity().getVector(i - 1, j)[0]);
  u_avg_im1_bottom
    = 0.5 * (flowField.getVelocity().getVector(i - 1, j)[0] + flowField.getVelocity().getVector(i - 1, j - 1)[0]);

  u_avg_i_top    = 0.5 * (flowField.getVelocity().getVector(i, j + 1)[0] + flowField.getVelocity().getVector(i, j)[0]);
  u_avg_i_bottom = 0.5 * (flowField.getVelocity().getVector(i, j)[0] + flowField.getVelocity().getVector(i, j - 1)[0]);

  u_avg_top    = 0.5 * (u_avg_i_top + u_avg_im1_top);
  u_avg_bottom = 0.5 * (u_avg_i_bottom + u_avg_im1_bottom);

  dudy_avg = (u_avg_top - u_avg_bottom) / parameters_.meshsize->getDy(i, j);

  // RealType nu_top, nu_bottom, nu_avg;
  // nu_top=  0.5*(flowField.getTurbulentViscosity().getScalar(i,j+1) +
  // flowField.getTurbulentViscosity().getScalar(i,j)); nu_bottom =
  // 0.5*(flowField.getTurbulentViscosity().getScalar(i,j) + flowField.getTurbulentViscosity().getScalar(i,j-1));

  // nu_avg = 0.5*(nu_top + nu_bottom);

  // tau_net = std::fabs(dudy_avg)*((1 / parameters_.flow.Re) + 0.5*(flowField.getTurbulentViscosity().getScalar(i,j+1)
  // + flowField.getTurbulentViscosity().getScalar(i,j-1)));
  tau_net = std::fabs(dudy_avg) * ((1 / parameters_.flow.Re) + flowField.getTurbulentViscosity().getScalar(i, j));

  //*********************************************************************
  // Wall Shear Stress Calculation
  // //******************************************************************
  RealType u_avg_im1_wall_top, u_avg_i_wall_top, u_avg_wall_top, dudy_wall_avg, tau_wall, tau;
  u_avg_im1_wall_top = 0.5
                       * (flowField.getVelocity().getVector(i - 1, 3)[0] + flowField.getVelocity().getVector(i - 1, 2)[0]);
  u_avg_i_wall_top = 0.5 * (flowField.getVelocity().getVector(i, 3)[0] + flowField.getVelocity().getVector(i, 2)[0]);

  u_avg_wall_top = 0.5 * (u_avg_i_wall_top + u_avg_im1_wall_top);

  dudy_wall_avg = 0.5 * (u_avg_wall_top) / parameters_.meshsize->getDy(i, 2);

  // std::cout<<i<<"   "<<j<<"   "<<flowField.getVelocity().getVector(i, 3)[0] <<"
  // "<<flowField.getVelocity().getVector(i, 2)[0] <<"  " << flowField.getVelocity().getVector(i, 1)[0] <<"  " <<
  // flowField.getVelocity().getVector(i, 0)[0]<<std::endl;

  // dudy_wall_avg = 2*(flowField.getVelocity().getVector(i, 2)[0]  ) / 0.5*(parameters_.meshsize->getDy(i,2) +
  // parameters_.meshsize->getDy(i,1) );

  tau_wall = std::fabs(dudy_wall_avg) * ((1 / parameters_.flow.Re));

  RealType u_tau, l_plus, y_plus, u_plus;

  // wall shear velocity
  u_tau = std::sqrt(tau_wall);
  // wall_unit
  l_plus = 1.0 / (parameters_.flow.Re * u_tau);
  // wall shear reynolds number
  y_plus = parameters_.meshsize->getPosY(i, j) / l_plus;

  u_plus = (0.5 * (flowField.getVelocity().getVector(i - 1, j)[0] + flowField.getVelocity().getVector(i, j)[0]))
           / (u_tau);

  tau = tau_net / 1.0;

  // if not an obstacle, write the data
  if ((flowField.getFlags().getValue(i, j) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " 0" << std::endl;
    flowField.getViscosity(viscosity, i, j);
    viscosityStream_ << viscosity << std::endl;
    tauStream_ << tau << std::endl;
    u_plusStream_ << u_plus << std::endl;
    y_plusStream_ << y_plus << std::endl;

#ifndef NDEBUG
    flowField.getH(h, i, j);
    flowField.getDelta(delta, i, j);
    hStream << h << std::endl;
    deltaStream << delta << std::endl;

#endif
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
    viscosityStream_ << "0.0" << std::endl;
    tauStream_ << "0.0" << std::endl;
    u_plusStream_ << "0.0" << std::endl;
    y_plusStream_ << "0.0" << std::endl;

#ifndef NDEBUG
    hStream << "0.0" << std::endl;
    deltaStream << "0.0" << std::endl;

#endif
  }
}

void Stencils::TurbulentVTKStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
  ASSERTION(FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 3);

  RealType pressure    = 0.0;
  RealType velocity[3] = {0.0, 0.0, 0.0};
  RealType viscosity   = 0.0;

  RealType u_avg_jm1, u_avg_j, u_avg_jp1, dudy_avg, tau;
  u_avg_jm1
    = 0.5 * (flowField.getVelocity().getVector(i, j - 1, k)[0] + flowField.getVelocity().getVector(i - 1, j - 1, k)[0]);
  u_avg_j = 0.5 * (flowField.getVelocity().getVector(i, j, k)[0] + flowField.getVelocity().getVector(i - 1, j, k)[0]);
  u_avg_jp1
    = 0.5 * (flowField.getVelocity().getVector(i, j + 1, k)[0] + flowField.getVelocity().getVector(i - 1, j + 1, k)[0]);

  dudy_avg
    = 0.5
      * (((u_avg_jp1 - u_avg_j) / parameters_.meshsize->getDy(i, j, k)) + ((u_avg_j - u_avg_jm1) / parameters_.meshsize->getDy(i, j, k)));
  tau = std::fabs(dudy_avg) * ((1 / parameters_.flow.Re) + flowField.getTurbulentViscosity().getScalar(i, j));

#ifndef NDEBUG
  RealType h     = 0.0;
  RealType delta = 0.0;
#endif
  if ((flowField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j, k);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
    flowField.getViscosity(viscosity, i, j, k);
    viscosityStream_ << viscosity << std::endl;
    tauStream_ << tau << std::endl;

#ifndef NDEBUG
    flowField.getH(h, i, j, k);
    flowField.getDelta(delta, i, j, k);
    hStream << h << std::endl;
    deltaStream << delta << std::endl;

#endif
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
    viscosityStream_ << "0.0" << std::endl;
    tauStream_ << "0.0" << std::endl;
#ifndef NDEBUG
    hStream << "0.0" << std::endl;
    deltaStream << "0.0" << std::endl;
#endif
  }
}

void Stencils::TurbulentVTKStencil::write(
  TurbulentFlowField& flowField, int timeStep, RealType simulationTime, std::string parameter
) {

  openFile(timeStep, simulationTime, parameter);

  if (FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 2) {
    // Write pressure
    ofile_
      << "CELL_DATA " << flowField.getNx() * flowField.getNy() << std::endl
      << "SCALARS pressure float 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");

    // Write viscosity
    ofile_ << "SCALARS viscosity float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << viscosityStream_.str() << std::endl;
    viscosityStream_.str("");

    // Write shear stress tau
    ofile_ << "SCALARS shear_stress float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << tauStream_.str() << std::endl;
    tauStream_.str("");

    // Write u_plus
    ofile_ << "SCALARS u_plus float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << u_plusStream_.str() << std::endl;
    u_plusStream_.str("");

    // Write y_plus
    ofile_ << "SCALARS y_plus float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << y_plusStream_.str() << std::endl;
    y_plusStream_.str("");

// Write nearest wall thickness(h) and boundary layer thickness(delta)
#ifndef NDEBUG
    ofile_ << "SCALARS h float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << hStream.str() << std::endl;
    hStream.str("");
    ofile_ << "SCALARS delta float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << deltaStream.str() << std::endl;
    deltaStream.str("");

#endif
  }

  if (FieldStencil<TurbulentFlowField>::parameters_.geometry.dim == 3) {
    // Write pressure
    ofile_
      << "CELL_DATA " << flowField.getNx() * flowField.getNy() * flowField.getNz() << std::endl
      << "SCALARS pressure float 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");

    // Write viscosity
    ofile_ << "SCALARS viscosity float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << viscosityStream_.str() << std::endl;
    viscosityStream_.str("");

    // Write shear stress
    ofile_ << "SCALARS shear_stress float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << tauStream_.str() << std::endl;
    tauStream_.str("");

#ifndef NDEBUG
    // Write , nearest wall distance(h), boundary layer thickness(delta)
    ofile_ << "SCALARS h float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << hStream.str() << std::endl;
    hStream.str("");
    ofile_ << "SCALARS delta float 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
    ofile_ << deltaStream.str() << std::endl;
    deltaStream.str("");

#endif
  }

  written_ = true;
  closeFile();
}

void Stencils::TurbulentVTKStencil::closeFile() { ofile_.close(); }
