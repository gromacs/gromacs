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

#ifndef _calc_verletbuf_h
#define _calc_verletbuf_h
#include "visibility.h"
#include "typedefs.h"

typedef struct
{
    int  cluster_size_i;  /* Cluster pair-list i-cluster size atom count */
    int  cluster_size_j;  /* Cluster pair-list j-cluster size atom count */
} verletbuf_list_setup_t;


/* Sets the pair-list setup assumed for the current Gromacs configuration.
 * The setup with smallest cluster sizes is return, such that the Verlet
 * buffer size estimated with this setup will be conservative.
 */
GMX_LIBGMXPREPROCESS_EXPORT
void verletbuf_get_list_setup(gmx_bool                bGPU,
                              verletbuf_list_setup_t *list_setup);


/* Calculate the non-bonded pair-list buffer size for the Verlet list
 * based on the particle masses, temperature, LJ types, charges
 * and constraints as well as the non-bonded force behavior at the cut-off.
 * The target is a maximum energy drift.
 * Returns the number of non-linear virtual sites. For these it's difficult
 * to determine their contribution to the drift exaclty, so we approximate.
 * Returns the pair-list cut-off.
 */
GMX_LIBGMXPREPROCESS_EXPORT
void calc_verlet_buffer_size(const gmx_mtop_t *mtop, real boxvol,
                             const t_inputrec *ir, real drift_target,
                             const verletbuf_list_setup_t *list_setup,
                             int *n_nonlin_vsite,
                             real *rlist);

#endif  /* _calc_verletbuf_h */
