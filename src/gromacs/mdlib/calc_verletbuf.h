/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_MDLIB_CALC_VERLETBUF_H
#define GMX_MDLIB_CALC_VERLETBUF_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_inputrec;

namespace gmx
{
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
 * \param[in] boxVolume     The volume of the unit cell
 * \param[in] inputrec      The input record
 * \param[in] nstlist       The pair list update frequency in steps (is not taken from \p inputrec)
 * \param[in] listLifetime  The lifetime of the pair-list, usually nstlist-1, but could be different
 * for dynamic pruning \param[in] referenceTemperature  The reference temperature for the ensemble
 * \param[in] listSetup     The pair-list setup
 * \returns The computed pair-list radius including buffer
 */
real calcVerletBufferSize(const gmx_mtop_t&         mtop,
                          real                      boxVolume,
                          const t_inputrec&         inputrec,
                          int                       nstlist,
                          int                       listLifetime,
                          real                      referenceTemperature,
                          const VerletbufListSetup& listSetup);

/* Convenience type */
using PartitioningPerMoltype = gmx::ArrayRef<const gmx::RangePartitioning>;

/* Determines the mininum cell size based on atom displacement
 *
 * The value returned is the minimum size for which the chance that
 * an atom or update group crosses to non nearest-neighbor cells
 * is <= chanceRequested within ir.nstlist steps.
 * Update groups are used when !updateGrouping.empty().
 * Without T-coupling, SD or BD, we can not estimate atom displacements
 * and fall back to the, crude, estimate of using the pairlist buffer size.
 *
 * Note: Like the Verlet buffer estimate, this estimate is based on
 *       non-interacting atoms and constrained atom-pairs. Therefore for
 *       any system that is not an ideal gas, this will be an overestimate.
 *
 * Note: This size increases (very slowly) with system size.
 */
real minCellSizeForAtomDisplacement(const gmx_mtop_t&      mtop,
                                    const t_inputrec&      ir,
                                    PartitioningPerMoltype updateGrouping,
                                    real                   chanceRequested);

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
