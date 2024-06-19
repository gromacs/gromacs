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
 * Implements the class holding parameters needed during simulation time
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "colvarssimulationsparameters.h"

#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/exceptions.h"

struct gmx_mtop_t;

namespace gmx
{

void ColvarsSimulationsParameters::setLocalAtomSetManager(LocalAtomSetManager* localAtomSetManager)
{
    localAtomSetManager_ = localAtomSetManager;
}

LocalAtomSetManager* ColvarsSimulationsParameters::localAtomSetManager() const
{
    if (localAtomSetManager_ == nullptr)
    {
        GMX_THROW(InternalError("Local atom manager not set for Colvars simulation."));
    }
    return localAtomSetManager_;
}

void ColvarsSimulationsParameters::setTopology(const gmx_mtop_t& mtop)
{
    gmxAtoms_ = gmx_mtop_global_atoms(mtop);
}

t_atoms ColvarsSimulationsParameters::topology() const
{
    return gmxAtoms_;
}


void ColvarsSimulationsParameters::setPeriodicBoundaryConditionType(const PbcType& pbcType)
{
    pbcType_ = std::make_unique<PbcType>(pbcType);
}

PbcType ColvarsSimulationsParameters::periodicBoundaryConditionType()
{
    if (pbcType_ == nullptr)
    {
        GMX_THROW(
                InternalError("Periodic boundary condition enum not set for Colvars simulation."));
    }
    return *pbcType_;
}


void ColvarsSimulationsParameters::setSimulationTimeStep(double timeStep)
{
    simulationTimeStep_ = timeStep;
}

double ColvarsSimulationsParameters::simulationTimeStep() const
{
    return simulationTimeStep_;
}


void ColvarsSimulationsParameters::setComm(const t_commrec& cr)
{
    cr_ = &cr;
}

const t_commrec* ColvarsSimulationsParameters::comm() const
{
    if (cr_ == nullptr)
    {
        GMX_THROW(InternalError("Communication record not set for Colvars simulation."));
    }
    return cr_;
}


void ColvarsSimulationsParameters::setLogger(const MDLogger& logger)
{
    logger_ = &logger;
}


const MDLogger* ColvarsSimulationsParameters::logger() const
{
    if (logger_ == nullptr)
    {
        GMX_THROW(InternalError("Logger not set for Colvars simulation."));
    }
    return logger_;
}

} // namespace gmx
