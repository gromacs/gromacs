/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_listed-forces Interactions between lists of particles
 * \ingroup group_mdrun
 *
 * \brief Computes energies and forces for interactions between a
 * small number of particles, e.g bonds.
 *
 * More functionality will move into this module shortly.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 */

/*! \file
 *
 * \brief
 * This file contains function declarations necessary for mdrun and tools
 * to compute energies and forces for bonded interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_listed-forces
 */

#ifndef GMX_LISTED_FORCES_BONDED_H
#define GMX_LISTED_FORCES_BONDED_H

#include <stdio.h>

#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_graph;
struct t_pbc;

/*! \brief Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When global_atom_index=NULL returns i+1.
 */
int glatnr(int *global_atom_index, int i);

/*! \brief Return whether this is an interaction that actually
 * calculates a potential and works on multiple atoms (not e.g. a
 * connection or a position restraint).
 *
 * \todo This function could go away when idef is not a big bucket of
 * everything. */
gmx_bool
ftype_is_bonded_potential(int ftype);

/*! \brief Calculates all bonded force interactions. */
void calc_bonds(const gmx_multisim_t *ms,
                const t_idef *idef,
                const rvec x[], history_t *hist,
                rvec f[], t_forcerec *fr,
                const struct t_pbc *pbc, const struct t_graph *g,
                gmx_enerdata_t *enerd, t_nrnb *nrnb, real *lambda,
                const t_mdatoms *md,
                t_fcdata *fcd, int *ddgatindex,
                int force_flags);

/*! \brief As calc_bonds, but only determines the potential energy
 * for the perturbed interactions.
 * The shift forces in fr are not affected. */
void calc_bonds_lambda(const t_idef *idef,
                       const rvec x[],
                       t_forcerec *fr,
                       const struct t_pbc *pbc, const struct t_graph *g,
                       gmx_grppairener_t *grpp, real *epot, t_nrnb *nrnb,
                       real *lambda,
                       const t_mdatoms *md,
                       t_fcdata *fcd, int *global_atom_index);

/*! \brief Position restraints require a different pbc treatment from other bondeds */
real posres(int nbonds,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec f[], rvec vir_diag,
            struct t_pbc *pbc,
            real lambda, real *dvdlambda,
            int refcoord_scaling, int ePBC, rvec comA, rvec comB);

/*! \brief Flat-bottom posres. Same PBC treatment as in normal position restraints */
real fbposres(int nbonds,
              const t_iatom forceatoms[], const t_iparams forceparams[],
              const rvec x[], rvec f[], rvec vir_diag,
              struct t_pbc *pbc, int refcoord_scaling, int ePBC, rvec com);

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
//! \endcond

#ifdef __cplusplus
}
#endif

#endif
