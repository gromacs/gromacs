/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
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

#include "nbnxn_kernel_simd_4xn_ewald_twin_comb_geom_ener.h"

#ifdef GMX_NBNXN_SIMD_4XN

#define CALC_COUL_EWALD
#define VDW_CUTOFF_CHECK /* Use twin-range cut-off */
#define LJ_COMB_GEOM
#define CALC_ENERGIES

#include "nbnxn_kernel_simd_4xn_preparation.h"

void
nbnxn_kernel_simd_4xn_ewald_twin_comb_geom_ener(const nbnxn_pairlist_t     *nbl,
                                                const nbnxn_atomdata_t     *nbat,
                                                const interaction_const_t  *ic,
                                                rvec                       *shift_vec,
                                                real                       *f,
                                                real                       *fshift,
                                                real                       *Vvdw,
                                                real                       *Vc)
{
#include "nbnxn_kernel_simd_4xn_outer.h"
}

#endif
