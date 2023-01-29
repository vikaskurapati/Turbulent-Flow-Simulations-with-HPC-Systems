#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"

#include "FlowField.hpp"
#include "Iterators.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"
#include "Stencils/ViscosityBufferFillStencil.hpp"
#include "Stencils/ViscosityBufferReadStencil.hpp"
#include "Stencils/TransportViscosityBufferReadStencil.hpp"
#include "Stencils/TransportViscosityBufferFillStencil.hpp"

namespace ParallelManagers {
  class PetscParallelManager {
    /**
     * @brief PetscParallelManager class to define all communication related attributes and methods for MPI
     * parallelisation
     *
     */

  protected:
    /**
     * @brief parameters for flow
     *
     */
    Parameters& parameters_;
    /**
     * @brief Flowfield for all the flow quantities
     *
     */
    FlowField& flowfield_;

    /**
     * @brief Stencil which fills the velocities during communication
     *
     */

    Stencils::VelocityBufferFillStencil fillVelocityStencil;
    /**
     * @brief Stencil which reads the velocities during communication
     *
     */
    Stencils::VelocityBufferReadStencil readVelocityStencil;
    /**
     * @brief Stencil which fills the pressures during communication
     *
     */
    Stencils::PressureBufferFillStencil fillPressureStencil;
    /**
     * @brief Stencil which reads the pressures during communication
     *
     */
    Stencils::PressureBufferReadStencil readPressureStencil;
    /**
     * @brief Iterator which operates on the velocityFillStencil
     *
     */
    ParallelBoundaryIterator<FlowField> velocityfillIterator;
    /**
     * @brief Iterator which operates on the velocityReadStencil
     *
     */
    ParallelBoundaryIterator<FlowField> velocityreadIterator;
    /**
     * @brief Iterator which operates on the pressureReadStencil
     *
     */
    ParallelBoundaryIterator<FlowField> pressurereadIterator;
    /**
     * @brief Iteratore which operates on the pressureFillStencil
     *
     */
    ParallelBoundaryIterator<FlowField> pressurefillIterator;

  public:
    /**
     * @brief Construct a new Petsc Parallel Manager object
     *
     * @param parameters parameters of the flow simulation
     * @param flowfield flowfield to store the flow quantities
     */
    PetscParallelManager(Parameters& parameters, FlowField& flowfield);
    /**
     * @brief method which communicates velocities using MPI
     *
     */
    void communicateVelocity();
    /**
     * @brief method which communicates pressures using MPI
     *
     */
    void communicatePressure();
    /**
     * @brief Destroy the Petsc Parallel Manager object
     *
     */
    virtual ~PetscParallelManager() = default;
  };

  class PetscTurbulentParallelManager: public PetscParallelManager {
    /**
     * @brief ParallelManger inherited from the default one to handle turbulent scenarios
     *
     */
  private:
    /**
     * @brief Turbulent flowField which holds flow quantities
     *
     */
    TurbulentFlowField& flowfield_;
    /**
     * @brief Stencil which fills the viscosities during communication
     *
     */
    Stencils::ViscosityBufferFillStencil fillViscosityStencil;
    /**
     * @brief Stencil which reads the viscosities during communication
     *
     */
    Stencils::ViscosityBufferReadStencil readViscosityStencil;
    /**
     * @brief Stencil which fills the transport viscosities during communication
     *
     */
    Stencils::TransportViscosityBufferFillStencil fillTransportViscosityStencil;
    /**
     * @brief Stencil which reads the transport viscosities during communication
     *
     */
    Stencils::TransportViscosityBufferReadStencil readTransportViscosityStencil;
    /**
     * @brief Iterator which operates on viscosityFillStencil
     *
     */
    ParallelBoundaryIterator<TurbulentFlowField> viscosityfillIterator;
    /**
     * @brief Iterator which operates on viscosityReadStencil
     *
     */
    ParallelBoundaryIterator<TurbulentFlowField> viscosityreadIterator;
        /**
     * @brief Iterator which operates on fillTransportViscosityStencil
     *
     */
    ParallelBoundaryIterator<TurbulentFlowField> transportViscosityfillIterator;
    /**
     * @brief Iterator which operates on readTransportViscosityStencil
     *
     */
    ParallelBoundaryIterator<TurbulentFlowField> transportViscosityreadIterator;

  public:
    /**
     * @brief Construct a new Petsc Turbulent Parallel Manager object
     *
     * @param parameters parameters of the flow simulation
     * @param flowfield Turbulent Flow field to hold the flow quantities
     */
    PetscTurbulentParallelManager(Parameters& parameters, TurbulentFlowField& flowfield);
    /**
     * @brief Method to communicate viscosity in a parallel simulation using MPI
     *
     */
    void communicateViscosity();
    /**
     * @brief Method to communicate transport viscosity in a parallel simulation using MPI
     *
     */
    void communicateTransportViscosity();
    /**
     * @brief Destroy the Petsc Turbulent Parallel Manager object
     *
     */
    virtual ~PetscTurbulentParallelManager() = default;
  };

} // namespace ParallelManagers