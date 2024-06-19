
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares parameters needed during simulation time
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_COLVARSIMULATIONSPARAMETERS_H
#define GMX_APPLIED_FORCES_COLVARSIMULATIONSPARAMETERS_H

#include <memory>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/logger.h"

struct gmx_mtop_t;

namespace gmx
{

/*! \internal
 * \brief Collect colvars parameters only available during simulation setup.
 *
 * To build the colvars force provider during simulation setup,
 * one needs access to parameters that become available only during simulation setup.
 *
 * This class collects these parameters via MdModuleNotifications in the
 * simulation setup phase and provides a check if all necessary parameters have
 * been provided.
 */
class ColvarsSimulationsParameters
{
public:
    ColvarsSimulationsParameters() = default;

    //! Set the local atom set Manager for colvars.
    void setLocalAtomSetManager(LocalAtomSetManager* localAtomSetManager);
    //! Get the local atom set Manager for colvars.
    LocalAtomSetManager* localAtomSetManager() const;


    /*! \brief Construct the topology of the system.
     *
     * \param[in] mtop is the pointer to the global topology struct
     */
    void setTopology(const gmx_mtop_t& mtop);

    //! Get the topology
    t_atoms topology() const;

    /*! \brief Set the periodic boundary condition via MdModuleNotifier.
     *
     * The pbc type is wrapped in PeriodicBoundaryConditionType to
     * allow the MdModuleNotifier to statically distinguish the callback
     * function type from other 'int' function callbacks.
     *
     * \param[in] pbcType enumerates the periodic boundary condition.
     */
    void setPeriodicBoundaryConditionType(const PbcType& pbcType);

    //! Get the periodic boundary conditions
    PbcType periodicBoundaryConditionType();

    //! Set the simulation time step
    void setSimulationTimeStep(double timeStep);
    //! Return the simulation time step
    double simulationTimeStep() const;

    //! Set the communicator
    void setComm(const t_commrec& cr);
    //! Return the communicator
    const t_commrec* comm() const;

    /*! \brief Set the logger for QMMM during mdrun
     * \param[in] logger Logger instance to be used for output
     */
    void setLogger(const MDLogger& logger);

    //! Get the logger instance
    const MDLogger* logger() const;

private:
    //! The LocalAtomSetManager
    LocalAtomSetManager* localAtomSetManager_;
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;
    //! The simulation time step
    double simulationTimeStep_ = 1;
    //! The topology
    t_atoms gmxAtoms_;
    //! The communicator
    const t_commrec* cr_;
    //! MDLogger for notifications during mdrun
    const MDLogger* logger_ = nullptr;


    // GMX_DISALLOW_COPY_AND_ASSIGN(ColvarsSimulationsParameters);
};

} // namespace gmx

#endif
