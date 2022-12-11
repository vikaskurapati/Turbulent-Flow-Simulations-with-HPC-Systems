#pragma once

#include "DataStructures.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"

/** Turbulent Flow field
 *
 * Class intended to contain the state of the domain for turbulent flows.
 */
class TurbulentFlowField: public FlowField {

private:
  ScalarField boundaryLayerThickness_; //! Scalar field representing the Local boundary layer thickness
  ScalarField wallDistance_;           //! Scalar field representing the nearest wall distance
  ScalarField turbulentViscosity_;     //! Scalar field representing the turbulent viscosity
  
public:
  /** Constructor for the 2D Turbulent flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   */
  TurbulentFlowField(int Nx, int Ny);

  /** Constructor for the 3D Turbulent flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
   */
  TurbulentFlowField(int Nx, int Ny, int Nz);

  /** Constructs a turbulent field from parameters object
   *
   * Constructs a turbulent field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  TurbulentFlowField(const Parameters& parameters);

  virtual ~TurbulentFlowField() = default;

  ScalarField& getBoundaryLayerThickness();
  ScalarField& getWallDistance();
  ScalarField& getTurbulentViscosity();

  void getViscosity(RealType& viscosity, int i, int j);
  void getViscosity(RealType& viscosity, int i, int j, int k);

  void getH(RealType& h, int i, int j);
  void getH(RealType& h, int i, int j, int k);

  void getDelta(RealType& delta, int i, int j);
  void getDelta(RealType& delta, int i, int j, int k);
};
