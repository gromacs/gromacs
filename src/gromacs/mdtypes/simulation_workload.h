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
 * \brief Describes work done on this domain that may change per-step.
 *
 * This work description is based on the SimulationWorkload in the context of the
 * current particle interactions assigned to this domain as well as other
 * factors that may change during the lifetime of a domain.
 *
 * Note that the contents of an object of this type is valid for
 * a single step and it is expected to be set at the beginning each step.
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
    /*! \brief Whether coordinate buffer ops are done on the GPU this step
     * \note This technically belongs to DomainLifetimeWorkload but due
     * to needing the flag before DomainLifetimeWorkload is built we keep
     * it here for now.
     */
    bool useGpuXBufferOps = false;
    //! Whether force buffer ops are done on the GPU this step
    bool useGpuFBufferOps = false;
    //! Whether PME forces are reduced with other contributions on the GPU this step
    bool useGpuPmeFReduction = false; // TODO: add this flag to the internal PME GPU data structures too
};

/*! \libinternal
 * \brief Describes work done on this domain on every step of its lifetime,
 * but which might change after the next domain paritioning.
 *
 * This work description is based on the SimulationWorkload in the context of the
 * current particle interactions assigned to this domain. The latter might change
 * after the next domain partitioning.
 *
 * An object of this type is updated every domain decomposition / neighbour search step
 * and reflects what work is required during the lifetime of a domain;
 * e.g. whether there are bonded interactions in this PP task.
 *
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
    //! Whether there are currently any local forces to be computed on the CPU
    bool haveCpuLocalForceWork = false;

    //! Whether the current nstlist step-range Free energy work on the CPU.
    bool haveFreeEnergyWork = false;
};

/*! \libinternal
 * \brief Manage what computation is required during the simulation.
 *
 * Holds information on the type of workload constant for the entire
 * simulation, and independent of the particle interactions handled
 * on any specific domain.
 *
 * An object of this type is constructed at the beginning of the
 * simulation and is expected to not change.
 * Additionally, the initialization is uniform across ranks of a
 * simulation, even with MPMD decomposition and separate PME ranks.
 */
class SimulationWorkload
{
public:
    //! If we have calculation of short range nonbondeds on CPU
    bool useCpuNonbonded = false;
    //! If we have calculation of short range nonbondeds on GPU
    bool useGpuNonbonded = false;
    //! If we have calculation of long range PME in GPU
    bool useCpuPme = false;
    //! If we have calculation of long range PME in GPU
    bool useGpuPme = false;
    //! If PME FFT solving is done on GPU.
    bool useGpuPmeFft = false;
    //! If bonded interactions are calculated on GPU.
    bool useGpuBonded = false;
    //! If update and constraint solving is performed on GPU.
    bool useGpuUpdate = false;
    //! If buffer operations are performed on GPU.
    bool useGpuBufferOps = false;
    //! If domain decomposition halo exchange is performed on GPU.
    bool useGpuHaloExchange = false;
    //! If direct PP-PME communication between GPU is used.
    bool useGpuPmePpCommunication = false;
    //! If direct GPU-GPU communication is enabled.
    bool useGpuDirectCommunication = false;
    //! If there is an Ewald surface (dipole) term to compute
    bool haveEwaldSurfaceContribution = false;
};

class MdrunScheduleWorkload
{
public:
    //! Workload descriptor for information constant for an entire run
    SimulationWorkload simulationWork;

    //! Workload descriptor for information constant for an nstlist range of steps
    DomainLifetimeWorkload domainWork;

    //! Workload descriptor for information that may change per-step
    StepWorkload stepWork;
};

} // namespace gmx

#endif // GMX_MDTYPES_SIMULATION_WORKLOAD_H
