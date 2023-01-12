#pragma once

#include <string>

#include "DataStructures.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Parameters.hpp"

/** Turbulent Flow field
 *
 * Class intended to contain the state of the domain for turbulent flows.
 */
class TurbulentFlowField: public FlowField {

private:
  ScalarField boundaryLayerThickness_;      //! Scalar field representing the Local boundary layer thickness
  ScalarField wallDistance_;                //! Scalar field representing the nearest wall distance
  ScalarField turbulentViscosity_;          //! Scalar field representing the turbulent viscosity
  ScalarField turbulentViscosityTransport_; //! Scalar field representing the variable which is being solved with the
                                            //! viscosity transport equation
  std::string method_;                      //! Method for Turbulence

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
  /**
   * @brief Destroy the Turbulent Flow Field object
   *
   */
  virtual ~TurbulentFlowField() = default;

  /**
   * @brief Get the Boundary Layer Thickness object
   *
   * @return ScalarField& reference of boundary layer thickness of the cell
   */
  ScalarField& getBoundaryLayerThickness();
  /**
   * @brief Get the Wall Distance object
   *
   * @return ScalarField& reference of wall distance of the cell
   */
  ScalarField& getWallDistance();
  /**
   * @brief Get the Turbulent Viscosity object
   *
   * @return ScalarField& reference of the turbulent viscosity of the cell
   */
  ScalarField& getTurbulentViscosity();

  /**
   * @brief Get the Turbulent Viscosity Transport object
   *
   * @return ScalarField& reference of the Scalar Field
   */
  ScalarField& getTurbulentViscosityTransport();

  /**
   * @brief Get the Viscosity object in 2D
   *
   * @param viscosity reference of viscosity to store viscosity of the flow field
   * @param i index in x
   * @param j index in y
   */
  void getViscosity(RealType& viscosity, int i, int j);
  /**
   * @brief Get the Viscosity object in 3D
   *
   * @param viscosity reference of viscosity to store viscosity of the flow field
   * @param i index in x
   * @param j index in y
   * @param k index in z
   */
  void getViscosity(RealType& viscosity, int i, int j, int k);
  /**
   * @brief Get the Turbulent Viscosity Transport object in 2D
   * 
   * @param viscosity reference of transport viscosity to store transport viscosity of the flow field
   * @param i index in x
   * @param j index in y
   */
  void getViscosityTransport(RealType& viscosity, int i, int j);
  /**
   * @brief Get the Turbulent Viscosity Transport object in 3D
   * 
   * @param viscosity reference of transport viscosity to store transport viscosity of the flow field
   * @param i index in x
   * @param j index in y
   * @param k index in z
   */
  void getViscosityTransport(RealType& viscosity, int i, int j, int k);
  /**
   * @brief get Wall Distance in 2D
   *
   * @param h reference to store the wall distance
   * @param i index in x
   * @param j index in y
   */
  void getH(RealType& h, int i, int j);
  /**
   * @brief get Wall Distance in 3D
   *
   * @param h reference to store the wall distance
   * @param i index in x
   * @param j index in y
   * @param k index in z
   */
  void getH(RealType& h, int i, int j, int k);
  /**
   * @brief Get the Delta object in 2D
   *
   * @param delta reference to store boundary layer thickness
   * @param i index in x
   * @param j index in y
   */
  void getDelta(RealType& delta, int i, int j);
  /**
   * @brief Get the Delta object
   *
   * @param delta reference to store boundary layer thickness
   * @param i index in x
   * @param j index in y
   * @param k index in z
   */
  void getDelta(RealType& delta, int i, int j, int k);
};
