/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * \ingroup module_listed-forces
 */

#ifndef GMX_LISTED_FORCES_BONDED_H
#define GMX_LISTED_FORCES_BONDED_H

#include <stdio.h>

#include "gromacs/legacyheaders/types/fcdata.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/mdatom.h"
#include "gromacs/legacyheaders/types/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc-simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_wallcycle;
struct t_graph;
struct t_pbc;

/*! \brief Calculate bond-angle. No PBC is taken into account (use mol-shift) */
real bond_angle(const rvec xi, const rvec xj, const rvec xk,
                const struct t_pbc *pbc,
                rvec r_ij, rvec r_kj, real *costh,
                int *t1, int *t2);  /* out */

/*! \brief Calculate dihedral-angle. No PBC is taken into account (use mol-shift) */
real dih_angle(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
               const struct t_pbc *pbc,
               rvec r_ij, rvec r_kj, rvec r_kl, rvec m, rvec n, /* out */
               real *sign,
               int *t1, int *t2, int *t3);

/*! \brief Do an update of the forces for dihedral potentials */
void do_dih_fup(int i, int j, int k, int l, real ddphi,
                rvec r_ij, rvec r_kj, rvec r_kl,
                rvec m, rvec n, rvec f[], rvec fshift[],
                const struct t_pbc *pbc, const struct t_graph *g,
                const rvec *x, int t1, int t2, int t3);

/*! \brief Make a dihedral fall in the range (-pi,pi) */
void make_dp_periodic(real *dp);

/*! \brief Compute CMAP dihedral energies and forces */
real
    cmap_dihs(int nbonds,
              const t_iatom forceatoms[], const t_iparams forceparams[],
              const gmx_cmap_t *cmap_grid,
              const rvec x[], rvec f[], rvec fshift[],
              const struct t_pbc *pbc, const struct t_graph *g,
              real gmx_unused lambda, real gmx_unused *dvdlambda,
              const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
              int  gmx_unused *global_atom_index);

//! \cond
/*************************************************************************
 *
 *  Bonded force functions
 *
 *************************************************************************/
t_ifunc bonds, g96bonds, morse_bonds, cubic_bonds, FENE_bonds, restraint_bonds;
t_ifunc angles, g96angles, cross_bond_bond, cross_bond_angle, urey_bradley, quartic_angles, linear_angles;
t_ifunc restrangles;
t_ifunc pdihs, idihs, rbdihs;
t_ifunc restrdihs, cbtdihs;
t_ifunc tab_bonds, tab_angles, tab_dihs;
t_ifunc polarize, anharm_polarize, water_pol, thole_pol, angres, angresz, dihres, unimplemented;

/* As pdihs(), but without calculating energies and shift forces */
void
    pdihs_noener(int nbonds,
                 const t_iatom forceatoms[], const t_iparams forceparams[],
                 const rvec x[], rvec f[],
                 const struct t_pbc gmx_unused *pbc,
                 const struct t_graph gmx_unused *g,
                 real lambda,
                 const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                 int gmx_unused *global_atom_index);

#ifdef GMX_SIMD_HAVE_REAL

/* As angles(), but using SIMD to calculate many angles at once.
 * This routines does not calculate energies and shift forces.
 */
void
    angles_noener_simd(int nbonds,
                       const t_iatom forceatoms[], const t_iparams forceparams[],
                       const rvec x[], rvec f[],
                       const struct t_pbc *pbc,
                       const struct t_graph gmx_unused *g,
                       real gmx_unused lambda,
                       const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                       int gmx_unused *global_atom_index);

/* As pdihs_noener(), but using SIMD to calculate many dihedrals at once. */
void
    pdihs_noener_simd(int nbonds,
                      const t_iatom forceatoms[], const t_iparams forceparams[],
                      const rvec x[], rvec f[],
                      const struct t_pbc *pbc,
                      const struct t_graph gmx_unused *g,
                      real gmx_unused lambda,
                      const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                      int gmx_unused *global_atom_index);

/* As rbdihs(), when not needing energy or shift force, using SIMD to calculate many dihedrals at once. */
void
    rbdihs_noener_simd(int nbonds,
                       const t_iatom forceatoms[], const t_iparams forceparams[],
                       const rvec x[], rvec f[],
                       const struct t_pbc *pbc,
                       const struct t_graph gmx_unused *g,
                       real gmx_unused lambda,
                       const t_mdatoms gmx_unused *md, t_fcdata gmx_unused *fcd,
                       int gmx_unused *global_atom_index);

#endif

//! \endcond

#ifdef __cplusplus
}
#endif

#endif
