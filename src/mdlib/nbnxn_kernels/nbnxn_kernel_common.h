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

#ifndef _nbnxn_kernel_common_h
#define _nbnxn_kernel_common_h

#include "typedefs.h"
#include "types/force_flags.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \brief Macro to make defining kernel functions in headers easy
 *
 * C doesn't let you declare or define functions using a typedef,
 * because that's impossible to parse correctly in all cases. So these
 * definitions support changing the NBNXN function signature at some
 * future time.
 */
#define DECLARE_NBNXN_KERNEL_ENER(kernel_name) \
void \
kernel_name(const nbnxn_pairlist_t     *nbl,\
            const nbnxn_atomdata_t     *nbat,\
            const interaction_const_t  *ic,\
            rvec                       *shift_vec,\
            real                       *f,\
            real                       *fshift,\
            real                       *Vvdw,\
            real                       *Vc)

#define DECLARE_NBNXN_KERNEL_NOENER(kernel_name) \
void \
kernel_name(const nbnxn_pairlist_t     *nbl,\
            const nbnxn_atomdata_t     *nbat,\
            const interaction_const_t  *ic,\
            rvec                       *shift_vec,\
            real                       *f,\
            real                       *fshift)

/*! \brief Typedefs for declaring lookup tables of kernel functions.
 */
typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                                const nbnxn_atomdata_t     *nbat,
                                const interaction_const_t  *ic,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vvdw,
                                real                       *Vc);

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                                  const nbnxn_atomdata_t     *nbat,
                                  const interaction_const_t  *ic,
                                  rvec                       *shift_vec,
                                  real                       *f,
                                  real                       *fshift);

/* Clear the force buffer f. Either the whole buffer or only the parts
 * used by the current thread when nbat->bUseBufferFlags is set.
 * In the latter case output_index is the task/thread list/buffer index.
 */
void
clear_f(const nbnxn_atomdata_t *nbat, int output_index, real *f);

/* Clear the shift forces */
void
clear_fshift(real *fshift);

/* Reduce the collected energy terms over the pair-lists/threads */
void
reduce_energies_over_lists(const nbnxn_atomdata_t     *nbat,
                           int                         nlist,
                           real                       *Vvdw,
                           real                       *Vc);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
