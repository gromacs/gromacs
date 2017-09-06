/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \brief
 * Declares functions to enforce rotational motion upon a group of particles.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_ROTATION_H
#define GMX_PULLING_PULL_ROTATION_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct MdrunOptions;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_rot;
class t_state;

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Initializes the enforced rotation groups.
 *
 * This routine does the memory allocation for various helper arrays, opens
 * the output files etc.
 *
 * \param fplog    General output file, normally md.log.
 * \param ir       Struct containing MD input parameters, among those
 *                 also the enforced rotation parameters.
 * \param nfile    Number of entries in the fnm structure.
 * \param fnm      The filenames struct containing also the names
 *                 of the rotation output files.
 * \param cr       Pointer to MPI communication data.
 * \param globalState  The global state, only used on the master rank.
 * \param mtop     Molecular topology.
 * \param oenv     Needed to open the rotation output xvgr file.
 * \param mdrunOptions  Options for mdrun.
 */
extern void init_rot(FILE *fplog, t_inputrec *ir, int nfile, const t_filenm fnm[],
                     struct t_commrec *cr, const t_state *globalState, gmx_mtop_t *mtop, const gmx_output_env_t *oenv,
                     const MdrunOptions &mdrunOptions);


/*! \brief Make a selection of the home atoms for all enforced rotation groups.
 *
 * This routine is similar to \ref dd_make_local_pull_groups, but works only with
 * domain decomposition. It should be called at every domain decomposition.
 *
 * \param dd      Structure containing domain decomposition data.
 * \param rot     Pointer to all the enforced rotation data.
 */
extern void dd_make_local_rotation_groups(struct gmx_domdec_t *dd, t_rot *rot);


/*! \brief Calculates the enforced rotation potential(s).
 *
 * This is the main enforced rotation module which is called during every time
 * step. Here the rotation potential as well as the resulting forces are
 * calculated.
 *
 * \param cr      Pointer to MPI communication data.
 * \param ir      Struct containing MD input parameters, among those
 * \param box     Simulation box, needed to make group whole.
 * \param x       The positions of all the local particles.
 * \param t       Time.
 * \param step    The time step.
 * \param bNS     After domain decomposition / neighbor searching several
 *                local arrays have to be updated (masses, shifts)
 */
extern void do_rotation(struct t_commrec *cr, t_inputrec *ir, matrix box, rvec x[], real t,
                        gmx_int64_t step, gmx_bool bNS);


/*! \brief Add the enforced rotation forces to the official force array.
 *
 * Adds the forces from enforced rotation potential to the local forces and
 * sums up the contributions to the rotation potential from all the nodes. Since
 * this needs communication, this routine should be called after the short range
 * forces have been evaluated (in order not to spoil cycle counts).
 * This routine also outputs data to the rotation output files (e.g.
 * the potential, the angle of the group(s), and torques).
 *
 * \param rot     Pointer to all the enforced rotation data.
 * \param f       The local forces to which the rotational forces have
 *                to be added.
 * \param cr      Pointer to MPI communication data.
 * \param step    The time step, used for output.
 * \param t       Time, used for output.
 * \returns       The potential energy of the rotation potentials.
 */
extern real add_rot_forces(t_rot *rot, rvec f[], struct t_commrec *cr, gmx_int64_t step, real t);


/*! \brief Close the enforced rotation output files.
 *
 * \param rot     Pointer to all the enforced rotation data.
 */
extern void finish_rot(t_rot *rot);


#ifdef __cplusplus
}
#endif


#endif
