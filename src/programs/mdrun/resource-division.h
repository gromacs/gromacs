/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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

#ifndef GMX_RESOURCE_DIVISION_H
#define GMX_RESOURCE_DIVISION_H

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Return the number of threads to use for thread-MPI based on how many
 * were requested, which algorithms we're using,
 * and how many particles there are.
 * At the point we have already called check_and_update_hw_opt.
 * Thus all options should be internally consistent and consistent
 * with the hardware, except that ntmpi could be larger than #GPU.
 */
int get_nthreads_mpi(const gmx_hw_info_t *hwinfo,
                     const gmx_hw_opt_t  *hw_opt,
                     const t_inputrec    *inputrec,
                     const gmx_mtop_t    *mtop,
                     const t_commrec     *cr,
                     FILE                *fplog);

/* Checks we can do when we don't (yet) know the cut-off scheme */
void check_and_update_hw_opt_1(gmx_hw_opt_t *hw_opt,
                               gmx_bool      bIsSimMaster);

/* Checks we can do when we know the cut-off scheme */
void check_and_update_hw_opt_2(gmx_hw_opt_t *hw_opt,
                               int           cutoff_scheme);

#ifdef __cplusplus
}
#endif

#endif /* GMX_RESOURCE_DIVISION_H */
