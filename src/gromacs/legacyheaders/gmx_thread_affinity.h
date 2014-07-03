/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
#ifndef GMX_THREAD_AFFINITY_H_
#define GMX_THREAD_AFFINITY_H_

#include <stdio.h>

#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

struct t_commrec;

/* Sets the thread affinity using the requested setting stored in hw_opt.
 * The hardware topologu is requested from hwinfo, when present.
 */
void
gmx_set_thread_affinity(FILE                       *fplog,
                        const struct t_commrec     *cr,
                        gmx_hw_opt_t               *hw_opt,
                        const gmx_hw_info_t        *hwinfo);

/* Check the process affinity mask and if it is found to be non-zero,
 * will honor it and disable mdrun internal affinity setting.
 * This function should be called first before the OpenMP library gets
 * initialized with the last argument FALSE (which will detect affinity
 * set by external tools like taskset), and later, after the OpenMP
 * initialization, with the last argument TRUE to detect affinity changes
 * made by the OpenMP library.
 *
 * Note that this will only work on Linux as we use a GNU feature.
 * With bAfterOpenmpInit false, it will also detect whether OpenMP environment
 * variables for setting the affinity are set.
 */
void
gmx_check_thread_affinity_set(FILE *fplog, const struct t_commrec *cr,
                              gmx_hw_opt_t *hw_opt, int ncpus,
                              gmx_bool bAfterOpenmpInit);

#ifdef __cplusplus
}
#endif

#endif /* GMX_THREAD_AFFINITY_H_ */
