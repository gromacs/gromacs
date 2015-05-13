/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * \brief Handles computing energies and forces for listed
 * interactions.
 *
 * Located here is the the code for
 * - computing energies and forces for interactions between a small
     number of particles, e.g bonds, position restraints and listed
     non-bonded interactions (e.g. 1-4).
 * - high-level functions used by mdrun for computing a set of such
     quantities
 * - managing thread-wise decomposition, thread-local buffer output,
     and reduction of output data across threads.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 */
/*! \file
 *
 * \brief This file contains declarations of high-level functions used
 * by mdrun to compute energies and forces for listed interactions.
 *
 * Clients of libgromacs that want to evaluate listed interactions
 * should call functions declared here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inpublicapi
 * \ingroup module_listed-forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_H
#define GMX_LISTED_FORCES_LISTED_FORCES_H

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Return whether this is an interaction that actually
 * calculates a potential and works on multiple atoms (not e.g. a
 * connection or a position restraint).
 *
 * \todo This function could go away when idef is not a big bucket of
 * everything. */
gmx_bool
ftype_is_bonded_potential(int ftype);

/*! \brief Calculates all listed force interactions.
 *
 * Note that pbc_full is used only for position restraints, and is
 * not initialized if there are none. */
void calc_listed(const gmx_multisim_t *ms,
                 struct gmx_wallcycle *wcycle,
                 const t_idef *idef,
                 const rvec x[], history_t *hist,
                 rvec f[], t_forcerec *fr,
                 const struct t_pbc *pbc, const struct t_pbc *pbc_full,
                 const struct t_graph *g,
                 gmx_enerdata_t *enerd, t_nrnb *nrnb, real *lambda,
                 const t_mdatoms *md,
                 t_fcdata *fcd, int *ddgatindex,
                 int force_flags);

/*! \brief As calc_listed(), but only determines the potential energy
 * for the perturbed interactions.
 *
 * The shift forces in fr are not affected. */
void calc_listed_lambda(const t_idef *idef,
                        const rvec x[],
                        t_forcerec *fr,
                        const struct t_pbc *pbc, const struct t_graph *g,
                        gmx_grppairener_t *grpp, real *epot, t_nrnb *nrnb,
                        real *lambda,
                        const t_mdatoms *md,
                        t_fcdata *fcd, int *global_atom_index);

/*! \brief Do all aspects of energy and force calculations for mdrun
 * on the set of listed interactions */
void
do_force_listed(struct gmx_wallcycle     *wcycle,
                matrix                    box,
                const t_lambda           *fepvals,
                const gmx_multisim_t     *ms,
                const t_idef             *idef,
                const rvec                x[],
                history_t                *hist,
                rvec                      f[],
                t_forcerec               *fr,
                const struct t_pbc       *pbc,
                const struct t_graph     *graph,
                gmx_enerdata_t           *enerd,
                t_nrnb                   *nrnb,
                real                     *lambda,
                const t_mdatoms          *md,
                t_fcdata                 *fcd,
                int                      *global_atom_index,
                int                       flags);

#ifdef __cplusplus
}
#endif

#endif
