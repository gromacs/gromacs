/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Describes work done on this domain by the current rank that may change per-step.
 *
 * This work description is based on the SimulationWorkload in the context of the
 * current particle interactions assigned to this domain as well as other
 * factors that may change during the lifetime of a domain.
 *
 * Note that unlike the other workload descriptors, these flags are also used on
 * dedicated PME ranks, hence the content is rank-specific (at least when it
 * comes to flags related to PME).
 *
 * Note that the contents of an object of this type is valid for
 * a single step and it is expected to be set at the beginning each step.
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
    //! Whether the slow forces need to be computed this step (in addition to the faster forces)
    bool computeSlowForces = false;
    //! Whether virial needs to be computed this step
    bool computeVirial = false;
    //! Whether energies need to be computed this step this step
    bool computeEnergy = false;
    //! Whether (any) forces need to be computed this step, not only energies
    bool computeForces = false;
    //! Whether only the MTS combined force buffers are needed and not the separate normal force buffer.
    bool useOnlyMtsCombinedForceBuffer = false;
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
    //! Whether GPU coordinates halo exchange is active this step
    bool useGpuXHalo = false;
    //! Whether GPU forces halo exchange is active this step
    bool useGpuFHalo = false;
    //! Whether GPU PME work is computed on the current rank this step (can be false on PP-only ranks or on fast steps with MTS)
    bool haveGpuPmeOnThisRank = false;
    //! Whether a separate PME rank has any work this step
    bool computePmeOnSeparateRank = false;
    //! Whether to combine the forces for multiple time stepping before the halo exchange
    bool combineMtsForcesBeforeHaloExchange = false;
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
    //! Whether the current nstlist step-range has bonded work to run on the CPU.
    bool haveCpuBondedWork = false;
    //! Whether the current nstlist step-range has listed (bonded + restraints) forces work to run on the CPU.
    bool haveCpuListedForceWork = false;
    //! Whether the current nstlist step-range has special forces on the CPU.
    bool haveSpecialForces = false;
    //! Whether there are currently any local forces to be computed on the CPU.
    bool haveCpuLocalForceWork = false;
    //! Whether there are currently any non-local forces to be computed on the CPU and, with GPU update and DD, later reduced on the GPU.
    bool haveCpuNonLocalForceWork = false;
    //! Whether the current nstlist step-range Free energy work on the CPU.
    bool haveFreeEnergyWork = false;
    //! Whether the CPU force buffer has contributions to local atoms that need to be reduced on the GPU (with DD).
    // This depends on whether there are CPU-based force tasks
    // or when DD is active the halo exchange has resulted in contributions
    // from the non-local part.
    bool haveLocalForceContribInCpuBuffer = false;
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
    //! Whether to compute nonbonded pair interactions
    bool computeNonbonded = false;
    //! Whether nonbonded pair forces are to be computed at slow MTS steps only
    bool computeNonbondedAtMtsLevel1 = false;
    //! Whether total dipole needs to be computed
    bool computeMuTot = false;
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
    //! If X buffer operations are performed on GPU.
    bool useGpuXBufferOps = false;
    //! If F buffer operations are performed on GPU.
    bool useGpuFBufferOps = false;
    //! If PP domain decomposition is active.
    bool havePpDomainDecomposition = false;
    //! If domain decomposition halo exchange is performed on CPU (in CPU-only runs or with staged GPU communication).
    bool useCpuHaloExchange = false;
    //! If domain decomposition halo exchange is performed on GPU.
    bool useGpuHaloExchange = false;
    //! If separate PME rank(s) are used.
    bool haveSeparatePmeRank = false;
    //! If PP-PME communication is done purely on CPU (in CPU-only runs or with staged GPU communication).
    bool useCpuPmePpCommunication = false;
    //! If direct PP-PME communication between GPU is used.
    bool useGpuPmePpCommunication = false;
    //! If direct GPU-GPU communication is enabled.
    bool useGpuDirectCommunication = false;
    //! If GPU PME decomposition is enabled.
    bool useGpuPmeDecomposition = false;
    //! If there is an Ewald surface (dipole) term to compute
    bool haveEwaldSurfaceContribution = false;
    //! Whether to use multiple time stepping
    bool useMts = false;
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
