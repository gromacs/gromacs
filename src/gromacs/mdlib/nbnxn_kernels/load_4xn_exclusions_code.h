/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
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
#ifndef _load_4xn_exclusions_code_h
#define _load_4xn_exclusions_code_h

#include "load_4xn_exclusions.h"
#include "gromacs/simd/types.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

static inline void
gmx_load_4xn_exclusions(unsigned int excl,
                        gmx_exclmask mask_S0,
                        gmx_exclmask mask_S1,
                        gmx_exclmask mask_S2,
                        gmx_exclmask mask_S3,
                        gmx_mm_pb   *interact_S0,
                        gmx_mm_pb   *interact_S1,
                        gmx_mm_pb   *interact_S2,
                        gmx_mm_pb   *interact_S3)
{
    /* Load integer interaction mask */
    gmx_exclmask mask_S = gmx_load1_exclmask(excl);
    *interact_S0  = gmx_checkbitmask_pb(mask_S, mask_S0);
    *interact_S1  = gmx_checkbitmask_pb(mask_S, mask_S1);
    *interact_S2  = gmx_checkbitmask_pb(mask_S, mask_S2);
    *interact_S3  = gmx_checkbitmask_pb(mask_S, mask_S3);
}

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
