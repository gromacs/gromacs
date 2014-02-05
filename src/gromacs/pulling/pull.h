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

/*! \libinternal \file
 *
 *
 * \brief
 * This file contains datatypes and function declarations necessary
   for mdrun to interface with the pull code.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_H
#define GMX_PULLING_PULL_H

#include "typedefs.h"
#include "../fileio/filenm.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Get the distance to the reference and deviation for pull coord coord_ind.
 *
 * \param[in]  pull         The pull group.
 * \param[in]  coord_ind    Number of the pull coordinate.
 * \param[in]  pbc          Information structure about periodicity.
 * \param[in]  t            Time.
 * \param[out] dr           The pull coordinate difference vector.
 * \param[out] dev          The deviation from the reference distance.
 */
void get_pull_coord_distance(const t_pull *pull,
                             int coord_ind,
                             const t_pbc *pbc, double t,
                             dvec dr, double *dev);


/*! \brief Set the all the pull forces to zero.
 *
 * \param pull              The pull group.
 */
void clear_pull_forces(t_pull *pull);


/*! \brief Determine the COM pull forces and add them to f, return the potential
 *
 * \param[in] ePull      Enum defining the type of pulling: umbrella, const force, ...
 * \param[in] pull       The pull group.
 * \param[in] md         All atoms.
 * \param[in] pbc        Information struct about periodicity.
 * \param[in] cr         Struct for communication info.
 * \param[in] t          Time.
 * \param[in] lambda     The value of lambda in FEP calculations.
 * \param[in] x          Positions.
 * \param[in] f          Forces.
 * \param[in,out] vir    The virial, which, if != NULL, gets a pull correction.
 * \param[out] dvdlambda Pull contribution to dV/d(lambda).
 *
 * \returns The pull potential energy.
 */
real pull_potential(int ePull, t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, rvec *f, tensor vir, real *dvdlambda);


/*! \brief Constrain the coordinates xp in the directions in x
 * and also constrain v when v != NULL.
 *
 * \param[in] pull       The pull group.
 * \param[in] md         All atoms.
 * \param[in] pbc        Information struct about periodicity.
 * \param[in] cr         Struct for communication info.
 * \param[in] dt         The time step length.
 * \param[in] t          The time.
 * \param[in] x          Positions.
 * \param[in,out] xp     Updated x, can be NULL.
 * \param[in,out] v      Velocities, which may get a pull correction.
 * \param[in,out] vir    The virial, which, if != NULL, gets a pull correction.
 */
void pull_constraint(t_pull *pull, t_mdatoms *md, t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir);


/*! \brief Make a selection of the home atoms for all pull groups.
 * Should be called at every domain decomposition.
 *
 * \param dd             Structure containing domain decomposition data.
 * \param pull           The pull group.
 * \param md             All atoms.
 */
void dd_make_local_pull_groups(gmx_domdec_t *dd,
                               t_pull *pull, t_mdatoms *md);


/*! \brief Get memory and initialize the fields of pull that still need it, and
 * do runtype specific initialization.
 *
 * \param fplog      General output file, normally md.log.
 * \param ir         The inputrec.
 * \param nfile      Number of files.
 * \param fnm        Standard filename struct.
 * \param mtop       The topology of the whole system.
 * \param cr         Struct for communication info.
 * \param oenv       Output options.
 * \param lambda     FEP lambda.
 * \param bOutFile   Open output files?
 * \param Flags      Flags passed over from main, used to determine
 *                   whether or not we are appending.
 */
void init_pull(FILE              *fplog,
               t_inputrec        *ir,
               int                nfile,
               const t_filenm     fnm[],
               gmx_mtop_t        *mtop,
               t_commrec        * cr,
               const output_env_t oenv,
               real               lambda,
               gmx_bool           bOutFile,
               unsigned long      Flags);


/*! \brief Close the pull output files.
 *
 * \param pull       The pull group.
 */
void finish_pull(t_pull *pull);


/*! \brief Print the pull output (x and/or f)
 *
 * \param pull     The pull group.
 * \param step     Time step number.
 * \param time     Time.
 */
void pull_print_output(t_pull *pull, gmx_int64_t step, double time);


/*! \brief Calculates centers of mass all pull groups.
 *
 * \param[in] cr       Struct for communication info.
 * \param[in] pull     The pull group.
 * \param[in] md       All atoms.
 * \param[in] pbc      Information struct about periodicity.
 * \param[in] t        Time, only used for cylinder ref.
 * \param[in] x        The local positions.
 * \param[in,out] xp   Updated x, can be NULL.
 *
 */
void pull_calc_coms(t_commrec *cr,
                    t_pull    *pull,
                    t_mdatoms *md,
                    t_pbc     *pbc,
                    double     t,
                    rvec       x[],
                    rvec      *xp);

#ifdef __cplusplus
}
#endif

#endif
