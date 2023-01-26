#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

namespace Stencils {

  /** Stencil for writting VTK files
   *
   * When iterated with, creates a VTK file.
   */
  class TurbulentVTKStencil: public FieldStencil<TurbulentFlowField> {
  private:
    bool          written_; //! Whether the file has already been written
    std::string   prefix_;  //! Prefix to be attached to the vtk files
    std::ofstream ofile_;   //! Output file stream

    std::stringstream pressureStream_;  //! Stream for the pressure data
    std::stringstream velocityStream_;  //! Stream for the velocity data
    std::stringstream viscosityStream_; //! Stream for the viscosity data
    std::stringstream hStream;          //! Stream for the nearest neighbour data
    std::stringstream deltaStream;      //! Stream for the boundary layer thickness
    std::stringstream tauStream_;      //! Stream for the shear stress
    std::stringstream u_plusStream_;   //! Stream for the u_plus
    std::stringstream y_plusStream_;   //! Stream for y_plus


    void writeVTKHeader(std::ostream& file) const;
    void writePoints(std::ostream& file, RealType simulationTime) const;

    /** Open output file
     * Opens the output file and prepares for writing.
     * It also writes the header of the file and the points of the grid.
     * @param timeStep Current time step of the simulation
     * @param simulationTime Current simulation time in double format
     */
    void openFile(int timeStep, RealType simulationTime);

    /** Finish writing. Must be called once the file has been written.
     *
     * Stores all the streams and closes the file.
     */
    void closeFile();

  public:
    /**
     * @brief Construct a new Turbulent V T K Stencil object
     *
     * @param parameters parameters of the flow simulation
     */
    TurbulentVTKStencil(const Parameters& parameters);
    /**
     * @brief Destroy the Turbulent V T K Stencil object
     *
     */
    ~TurbulentVTKStencil() override = default;
    /**
     * @brief add the respective streams in 2D
     *
     * @param turbulentFlowField data structure holding the flow quantities
     * @param i index in x
     * @param j index in y
     */
    void apply(TurbulentFlowField& turbulentFlowField, int i, int j) override;
    /**
     * @brief add the respective streams in 3D
     *
     * @param turbulentFlowField data structure holding the flow quantities
     * @param i index in x
     * @param j index in y
     * @param k index in z
     */
    void apply(TurbulentFlowField& turbulentFlowField, int i, int j, int k) override;
    /**
     * @brief function to write the vtk
     *
     * @param turbulentFlowField data structure holding the flow quantities
     * @param timeStep timestep of the global domain
     * @param simulationTime simulation time elapsed
     */
    void write(TurbulentFlowField& turbulentFlowField, int timeStep, RealType simulationTime);
  };

} // namespace Stencils
