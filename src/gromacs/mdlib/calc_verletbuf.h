/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int  cluster_size_i;  /* Cluster pair-list i-cluster size atom count */
    int  cluster_size_j;  /* Cluster pair-list j-cluster size atom count */
} verletbuf_list_setup_t;


/* Add a 5% and 10% rlist buffer for simulations without dynamics (EM, NM, ...)
 * and NVE simulations with zero initial temperature, respectively.
 * 10% should be enough for any NVE simulation with PME and nstlist=10,
 * for other settings it might not be enough, but then it's difficult
 * to come up with any reasonable (not crazily expensive) value
 * and grompp will notify the user when using the 10% buffer.
 */
static const real verlet_buffer_ratio_nodynamics = 0.05;
static const real verlet_buffer_ratio_NVE_T0     = 0.10;


/* Sets the pair-list setup assumed for the current Gromacs configuration.
 * The setup with smallest cluster sizes is return, such that the Verlet
 * buffer size estimated with this setup will be conservative.
 * bSIMD tells if to take into account SIMD, when supported.
 * bGPU tells to estimate for GPU kernels (bSIMD is ignored with bGPU=TRUE)
 */
void verletbuf_get_list_setup(gmx_bool                bSIMD,
                              gmx_bool                bGPU,
                              verletbuf_list_setup_t *list_setup);


/* Calculate the non-bonded pair-list buffer size for the Verlet list
 * based on the particle masses, temperature, LJ types, charges
 * and constraints as well as the non-bonded force behavior at the cut-off.
 * If reference_temperature < 0, the maximum coupling temperature will be used.
 * The target is a maximum energy drift of ir->verletbuf_tol.
 * Returns the number of non-linear virtual sites. For these it's difficult
 * to determine their contribution to the drift exaclty, so we approximate.
 * Returns the pair-list cut-off.
 */
void calc_verlet_buffer_size(const gmx_mtop_t *mtop, real boxvol,
                             const t_inputrec *ir,
                             real reference_temperature,
                             const verletbuf_list_setup_t *list_setup,
                             int *n_nonlin_vsite,
                             real *rlist);

#ifdef __cplusplus
}
#endif

#endif
