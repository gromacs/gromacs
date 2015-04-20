/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef _orires_h
#define _orires_h

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_pbc;

void init_orires(FILE *fplog, const gmx_mtop_t *mtop,
                 rvec x[],
                 const t_inputrec *ir,
                 const t_commrec *cr, t_oriresdata *od,
                 t_state *state);
/* Decides whether orientation restraints can work, and initializes
   all the orientation restraint stuff in *od (and assumes *od is
   already allocated. */

real calc_orires_dev(const gmx_multisim_t *ms,
                     int nfa, const t_iatom fa[], const t_iparams ip[],
                     const t_mdatoms *md, const rvec x[],
                     const struct t_pbc *pbc, t_fcdata *fcd, history_t *hist);
/*
 * Calculates the time averaged D matrices, the S matrix for each experiment.
 * Returns the weighted RMS deviation of the orientation restraints.
 */

void diagonalize_orires_tensors(t_oriresdata *od);
/*
 * Diagonalizes the order tensor(s) of the orienation restraints.
 * For each experiment eig containts first 3 eigenvalues and then
 * the 3 eigenvectors. The eigenvalues are ordered on magnitude.
 */

void print_orires_log(FILE *log, t_oriresdata *od);
/* Print order parameter, eigenvalues and eigenvectors to the log file */

t_ifunc orires;
/* Does only the orientation restraint force calculation */

void update_orires_history(t_fcdata *fcd, history_t *hist);
/* Copy the new time averages that have been calculated in calc_orires_dev */

#ifdef __cplusplus
}
#endif

#endif  /* _orires_h */
