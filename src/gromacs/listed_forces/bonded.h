/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *
 * \brief This file contains declarations necessary for low-level
 * functions for computing energies and forces for bonded interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_listed_forces
 */

#ifndef GMX_LISTED_FORCES_BONDED_H
#define GMX_LISTED_FORCES_BONDED_H

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_cmap_t;
struct t_fcdata;
struct t_mdatom;
struct t_nrnb;
struct t_pbc;

/*! \brief Calculate bond-angle. No PBC is taken into account (use mol-shift) */
real bond_angle(const rvec          xi,
                const rvec          xj,
                const rvec          xk,
                const struct t_pbc* pbc,
                rvec                r_ij,
                rvec                r_kj,
                real*               costh,
                int*                t1,
                int*                t2); /* out */

/*! \brief Calculate dihedral-angle. No PBC is taken into account (use mol-shift) */
real dih_angle(const rvec          xi,
               const rvec          xj,
               const rvec          xk,
               const rvec          xl,
               const struct t_pbc* pbc,
               rvec                r_ij,
               rvec                r_kj,
               rvec                r_kl,
               rvec                m,
               rvec                n, /* out */
               int*                t1,
               int*                t2,
               int*                t3);

/*! \brief Do an update of the forces for dihedral potentials */
void do_dih_fup(int                 i,
                int                 j,
                int                 k,
                int                 l,
                real                ddphi,
                rvec                r_ij,
                rvec                r_kj,
                rvec                r_kl,
                rvec                m,
                rvec                n,
                rvec4               f[],
                rvec                fshift[],
                const struct t_pbc* pbc,
                const rvec*         x,
                int                 t1,
                int                 t2,
                int                 t3);

/*! \brief Make a dihedral fall in the range (-pi,pi) */
void make_dp_periodic(real* dp);

/*! \brief Compute CMAP dihedral energies and forces */
real cmap_dihs(int                 nbonds,
               const t_iatom       forceatoms[],
               const t_iparams     forceparams[],
               const gmx_cmap_t*   cmap_grid,
               const rvec          x[],
               rvec4               f[],
               rvec                fshift[],
               const struct t_pbc* pbc,
               real gmx_unused lambda,
               real gmx_unused* dvdlambda,
               const t_mdatoms gmx_unused* md,
               t_fcdata gmx_unused* fcd,
               int gmx_unused* global_atom_index);

/*! \brief For selecting which flavor of bonded kernel is used for simple bonded types */
enum class BondedKernelFlavor
{
    ForcesSimdWhenAvailable, //!< Compute only forces, use SIMD when available; should not be used with perturbed parameters
    ForcesNoSimd,             //!< Compute only forces, do not use SIMD
    ForcesAndVirialAndEnergy, //!< Compute forces, virial and energy (no SIMD)
    ForcesAndEnergy,          //!< Compute forces and energy (no SIMD)
    Count                     //!< The number of flavors
};

/*! \brief Returns whether the energy should be computed */
static constexpr inline bool computeEnergy(const BondedKernelFlavor flavor)
{
    return (flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy
            || flavor == BondedKernelFlavor::ForcesAndEnergy);
}

/*! \brief Returns whether the virial should be computed */
static constexpr inline bool computeVirial(const BondedKernelFlavor flavor)
{
    return (flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy);
}

/*! \brief Returns whether the energy and/or virial should be computed */
static constexpr inline bool computeEnergyOrVirial(const BondedKernelFlavor flavor)
{
    return (flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy
            || flavor == BondedKernelFlavor::ForcesAndEnergy);
}

/*! \brief Calculates bonded interactions for simple bonded types
 *
 * Exits with an error when the bonded type is not simple
 * All pointers should be non-null, except for pbc and g which can be nullptr.
 * \returns the energy or 0 when \p bondedKernelFlavor did not request the energy.
 */
real calculateSimpleBond(int                 ftype,
                         int                 numForceatoms,
                         const t_iatom       forceatoms[],
                         const t_iparams     forceparams[],
                         const rvec          x[],
                         rvec4               f[],
                         rvec                fshift[],
                         const struct t_pbc* pbc,
                         real                lambda,
                         real*               dvdlambda,
                         const t_mdatoms*    md,
                         t_fcdata*           fcd,
                         int gmx_unused*    global_atom_index,
                         BondedKernelFlavor bondedKernelFlavor);

//! Getter for finding the flop count for an \c ftype interaction.
int nrnbIndex(int ftype);

#endif
