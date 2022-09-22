/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#ifndef GMX_MDLIB_CALC_VERLETBUF_H
#define GMX_MDLIB_CALC_VERLETBUF_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_inputrec;

namespace gmx
{
template<typename>
class ArrayRef;
class RangePartitioning;
} // namespace gmx

namespace Nbnxm
{
enum class KernelType;
} // namespace Nbnxm


struct VerletbufListSetup
{
    int cluster_size_i; /* Cluster pair-list i-cluster size atom count */
    int cluster_size_j; /* Cluster pair-list j-cluster size atom count */
};


/* Add a 5% and 10% rlist buffer for simulations without dynamics (EM, NM, ...)
 * and NVE simulations with zero initial temperature, respectively.
 * 10% should be enough for any NVE simulation with PME and nstlist=10,
 * for other settings it might not be enough, but then it's difficult
 * to come up with any reasonable (not crazily expensive) value
 * and grompp will notify the user when using the 10% buffer.
 */
static const real verlet_buffer_ratio_nodynamics = 0.05;
static const real verlet_buffer_ratio_NVE_T0     = 0.10;


/* Returns the pair-list setup for the given nbnxn kernel type.
 */
VerletbufListSetup verletbufGetListSetup(Nbnxm::KernelType nbnxnKernelType);

//! \brief Chance target to use in minCellSizeForAtomDisplacement()
enum class ChanceTarget
{
    Atom,  //<! Indicates that the chance passed is applied to each atom individually
    System //<! Indicates that the chance passed is for the max displacement over all atoms
};

/* Enum for choosing the list type for verletbufGetSafeListSetup() */
enum class ListSetupType
{
    CpuNoSimd,            /* CPU Plain-C 4x4 list */
    CpuSimdWhenSupported, /* CPU 4xN list, where N=4 when the binary doesn't support SIMD or the smallest N supported by SIMD in this binary */
    Gpu                   /* GPU (8x2x)8x4 list */
};

/* Returns the pair-list setup assumed for the current Gromacs configuration.
 * The setup with smallest cluster sizes is returned, such that the Verlet
 * buffer size estimated with this setup will be conservative.
 */
VerletbufListSetup verletbufGetSafeListSetup(ListSetupType listType);

/* Returns the atom density weighted over cells of size \p cutoff
 *
 * A grid is put over the unit cell with a grid size of, approximately, cutoff.
 * The atom count is determined for each set. The effective density is then
 * calculated as the cell density weighted by atom count.
 * This in intended for passing to \p calcVerletBufferSize().
 *
 * \note This is an expensive function as it needs to loop over all atoms
 *       in the system.
 * \note When communicator!=MPI_COMM_NULL computes on rank 0 and broadcasts to the other ranks.
 *
 * \param[in] coordinates   The coordinates of all atoms, can be empty on non-main ranks
 * \param[in] box           The simulation unit cell
 * \param[in] cutoff        The maximum non-bonded interaction cut-off distance
 * \param[in] communicator  MPI communicator, can be MPI_COMM_NULL when not called in parallel
 */
real computeEffectiveAtomDensity(gmx::ArrayRef<const gmx::RVec> coordinates,
                                 const matrix                   box,
                                 real                           cutoff,
                                 MPI_Comm                       communicator);

/* Returns the non-bonded pair-list radius including computed buffer
 *
 * Calculate the non-bonded pair-list buffer size for the Verlet list
 * based on the particle masses, temperature, LJ types, charges
 * and constraints as well as the non-bonded force behavior at the cut-off.
 * The pair list update frequency and the list lifetime, which is nstlist-1
 * for normal pair-list buffering, are passed separately, as in some cases
 * we want an estimate for different values than the ones set in the inputrec.
 * If referenceTemperature < 0, the maximum coupling temperature will be used.
 * The target is a maximum average energy jump per atom of
 * inputrec.verletbuf_tol*nstlist*inputrec.delta_t over the lifetime of the list.
 *
 * \note For non-linear virtual sites it can be problematic to determine their
 *       contribution to the drift exaclty, so we approximate.
 *
 * \param[in] mtop          The system topology
 * \param[in] effectiveAtomDensity  The effective atom density, use computeEffectiveAtomDensity()
 * \param[in] inputrec      The input record
 * \param[in] nstlist       The pair list update frequency in steps (is not taken from \p inputrec)
 * \param[in] listLifetime  The lifetime of the pair-list, usually nstlist-1, but could be different
 *                          for dynamic pruning
 * \param[in] referenceTemperature  The reference temperature for the ensemble
 * \param[in] listSetup     The pair-list setup
 * \returns The computed pair-list radius including buffer
 */
real calcVerletBufferSize(const gmx_mtop_t&         mtop,
                          real                      effectiveAtomDensity,
                          const t_inputrec&         inputrec,
                          int                       nstlist,
                          int                       listLifetime,
                          real                      referenceTemperature,
                          const VerletbufListSetup& listSetup);

/* Convenience type */
using PartitioningPerMoltype = gmx::ArrayRef<const gmx::RangePartitioning>;

/*! \brief Determines the minimum cell size based on atom displacement
 *
 * The value returned is the minimum size for which the chance that
 * each atom (with \p chanceTarget=ChangeTarget::Atom) or any atom or
 * update group (with \p chanceTarget=ChangeTarget::System) crosses into
 * a non nearest-neighbor cell is <= chanceRequested within ir.nstlist steps.
 * Update groups are used when !updateGrouping.empty(). In that case
 * the displacement of the center of geometry of the update group is considered.
 * Without T-coupling, SD or BD, we can not estimate atom displacements
 * and fall back to the, crude, estimate of using the pair list buffer size.
 *
 * Note: Like the Verlet buffer estimate, this estimate is based on
 *       non-interacting atoms and constrained atom-pairs. Therefore for
 *       any system that is not an ideal gas, this will be an overestimate.
 *
 * Note: With \p chanceTarget=ChangeTarget::System this size increases
 *       (very slowly) with the number of atoms in the system.
 *
 * \param[in] mtop            The system topology
 * \param[in] ir              The input record
 * \param[in] updateGrouping  The update grouping within each molecule type
 * \param[in] chanceRequested The requested chance
 * \param[in] chanceTarget    Whether \p chance refors to a displacement per-atom or maximum over all atoms
 *
 * \returns  The minimum cell size given \p chanceRequested.
 */
real minCellSizeForAtomDisplacement(const gmx_mtop_t&      mtop,
                                    const t_inputrec&      ir,
                                    PartitioningPerMoltype updateGrouping,
                                    real                   chanceRequested,
                                    ChanceTarget           chanceTarget);

/* Struct for unique atom type for calculating the energy drift.
 * The atom displacement depends on mass and constraints.
 * The energy jump for given distance depend on LJ type and q.
 */
struct atom_nonbonded_kinetic_prop_t
{
    real mass     = 0;     /* mass */
    int  type     = 0;     /* type (used for LJ parameters) */
    real q        = 0;     /* charge */
    bool bConstr  = false; /* constrained, if TRUE, use #DOF=2 iso 3 */
    real con_mass = 0;     /* mass of heaviest atom connected by constraints */
    real con_len  = 0;     /* constraint length to the heaviest atom */
};

/* This function computes two components of the estimate of the variance
 * in the displacement of one atom in a system of two constrained atoms.
 * Returns in sigma2_2d the variance due to rotation of the constrained
 * atom around the atom to which it constrained.
 * Returns in sigma2_3d the variance due to displacement of the COM
 * of the whole system of the two constrained atoms.
 *
 * Only exposed here for testing purposes.
 */
void constrained_atom_sigma2(real kT_fac, const atom_nonbonded_kinetic_prop_t* prop, real* sigma2_2d, real* sigma2_3d);

#endif
