#include "StdAfx.hpp"
#include <string>

#include "Simulation.hpp"
#include "TurbulentFlowField.hpp"

#include "Stencils/MaxNuStencil.hpp"

class TurbulentSimulation: public Simulation {
  /**
   * @brief Inherited Simulation class for Turbulent cases
   * 
   */

public:
/**
 * @brief Construct a new Turbulent Simulation object
 * 
 * @param parameters parameters of the flow simulation
 * @param turbulentFlowField flowField holding turbulent flow quantities
 */
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& turbulentFlowField);
  /**
   * @brief Destroy the Turbulent Simulation object
   * 
   */
  virtual ~TurbulentSimulation() = default;

private:
/**
 * @brief flowField holding turbulent flow quantities
 * 
 */
  TurbulentFlowField                  turbulentFlowField_;
  /**
   * @brief Stencil to calculate maximum turbulent viscosity
   * 
   */
  Stencils::MaxNuStencil              maxNuStencil_;
  /**
   * @brief iterator to operate on maxNuStencil_
   * 
   */
  FieldIterator<TurbulentFlowField>   maxNuFieldIterator_;
  /**
   * @brief Stencil for turbulent viscosity calculations
   * 
   */
  Stencils::TurbulentViscosityStencil turbulentViscosityStencil_;
  /**
   * @brief Iterator to operate on turbulentViscosityStencil_
   * 
   */
  FieldIterator<TurbulentFlowField>   turbulentViscosityIterator_;
/**
 * @brief Function to initialise the flow field
 * 
 */
  void initializeFlowField() override;
/**
 * @brief FGHStencil to incorportate turbulent functionalities
 * 
 */
  Stencils::TurbulentFGHStencil     turbulentfghStencil_;
  /**
   * @brief iterator to operate on turbulentfghStencil_
   * 
   */
  FieldIterator<TurbulentFlowField> turbulentfghIterator_;
/**
 * @brief method to solve one time step of the simulation
 * 
 */
  void solveTimestep() override;
/**
 * @brief method to write the VTK
 * 
 * @param timeStep timestep of the global domain
 * @param simulationTime time of the simulation elapsed
 */
  void plotVTK(int timeStep, RealType simulationTime) override;
  /**
   * @brief Set the Time Step object
   * 
   */
  void setTimeStep() override;
  /**
   * @brief parallel manager to incorporate communication of turbulent attributes
   * 
   */
  ParallelManagers::PetscTurbulentParallelManager parallel_manager_;

  std::string method_;
};