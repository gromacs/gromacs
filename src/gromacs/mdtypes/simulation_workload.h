/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief Declares step, domain-lifetime, and run workload managers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDTYPES_SIMULATION_WORKLOAD_H
#define GMX_MDTYPES_SIMULATION_WORKLOAD_H

namespace gmx
{

/*! \libinternal
 * \brief Data structure that describes work that can change per-step.
 *
 * Note that the contents of an object of this type has a lifetime
 * of a single step and it is expected to be set at the beginning each step.
 *
 * The initial set of flags map the legacy force flags to boolean flags;
 * these have the role of directing per-step compute tasks undertaken by a PP rank.
 *
 */
class StepWorkload
{
    public:
        //! Whether the state has changed, always set unless TPI is used.
        bool stateChanged = false;
        //! Whether the box might have changed
        bool haveDynamicBox = false;
        //! Whether neighbor searching needs to be done this step
        bool doNeighborSearch = false;
        //! Whether virial needs to be computed this step
        bool computeVirial = false;
        //! Whether energies need to be computed this step this step
        bool computeEnergy = false;
        //! Whether (any) forces need to be computed this step, not only energies
        bool computeForces = false;
        //! Whether nonbonded forces need to be computed this step
        bool computeNonbondedForces = false;
        //! Whether listed forces need to be computed this step
        bool computeListedForces = false;
        //! Whether this step DHDL needs to be computed
        bool computeDhdl = false;
};

/*! \libinternal
 * \brief Manage computational work that has the lifetime of decomposition.
 *
 * An object of this type is updated every decomposition step
 * (i.e. domain decomposition / neighbour search)
 * reflecting what work is required during the lifetime of a domain.
 * e.g. whether there are bonded interactions in this PP task.
 *
 * This will remove the desire for inline getters from modules that
 * describe whether they have work to do, because that can be set up
 * once per simulation or neighborlist lifetime and not changed
 * thereafter.
 */
class DomainLifetimeWorkload
{
    public:
        //! Whether the current nstlist step-range has bonded work to run on a GPU.
        bool haveGpuBondedWork = false;
        //! Whether the current nstlist step-range has bonded work to run on he CPU.
        bool haveCpuBondedWork = false;
        //! Whether the current nstlist step-range has restraints work to run on he CPU.
        bool haveRestraintsWork = false;
        //! Whether the current nstlist step-range has listed forces work to run on he CPU.
        //  Note: currently this is haveCpuBondedWork | haveRestraintsWork
        bool haveCpuListedForceWork = false;
        //! Whether the current nstlist step-range has special forces on the CPU.
        bool haveSpecialForces = false;

        // TODO
        //! Whether the current nstlist step-range Free energy work on the CPU.
        bool haveFreeEnergyWork = false;
};

/*! \libinternal
 * \brief Manage what computation is required during the simulation.
 *
 * Holds information on the type of workload constant for the entire
 * simulation.
 *
 * An object of this type is constructed at the beginning of the
 * simulation and is expected to not change.
 */
class SimulationWorkload
{
};

class MdrunScheduleWorkload
{
    public:
        //! Workload descriptor for information constant for an entire run
        gmx::SimulationWorkload     simulationWork;

        //! Workload descriptor for information constant for an nstlist range of steps
        gmx::DomainLifetimeWorkload domainWork;

        //! Workload descriptor for information that may change per-step
        gmx::StepWorkload           stepWork;
};

}      // namespace gmx

#endif // GMX_MDTYPES_SIMULATION_WORKLOAD_H
