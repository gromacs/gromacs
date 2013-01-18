/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _edsam_h
#define _edsam_h
#include "visibility.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void do_edsam(t_inputrec *ir, gmx_large_int_t step,
              t_commrec *cr, rvec xs[], rvec v[], matrix box, gmx_edsam_t ed);
/* Essential dynamics constraints, called from constrain() */

GMX_LIBMD_EXPORT
gmx_edsam_t ed_open(int natoms, edsamstate_t *EDstate, int nfile, const t_filenm fnm[],
                    unsigned long Flags, const output_env_t oenv, t_commrec *cr);
/* Sets the ED input/output filenames, opens output file */

void init_edsam(gmx_mtop_t *mtop, t_inputrec *ir, t_commrec *cr,
                gmx_edsam_t ed, rvec x[], matrix box, edsamstate_t *edsamstate);
/* Init routine for ED and flooding. Calls init_edi in a loop for every .edi-cycle
 * contained in the input file, creates a NULL terminated list of t_edpar structures */

void dd_make_local_ed_indices(gmx_domdec_t *dd, gmx_edsam_t ed);
/* Make a selection of the home atoms for the ED groups.
 * Should be called at every domain decomposition. */

void do_flood(t_commrec *cr, t_inputrec *ir, rvec x[], rvec force[], gmx_edsam_t ed,
              matrix box, gmx_large_int_t step, gmx_bool bNS);
/* Flooding - called from do_force() */

#ifdef __cplusplus
}
#endif

#endif  /* _edsam_h */
