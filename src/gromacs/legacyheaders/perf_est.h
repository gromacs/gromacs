/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014,2015, by the GROMACS development team, led by
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
#ifndef _perf_est_h
#define _perf_est_h

#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_mtop_t;

void count_bonded_distances(struct gmx_mtop_t *mtop, const t_inputrec *ir,
                            double *ndistance_c, double *ndistance_simd);
/* Count the number of distance calculations in bonded interactions,
 * separately for plain-C and SIMD bonded functions.
 * The computational cost is nearly proportional to the numbers.
 * It is allowed to pass NULL for the last two arguments.
 */

float pme_load_estimate(struct gmx_mtop_t *mtop, t_inputrec *ir, matrix box);
/* Returns an estimate for the relative load of the PME mesh calculation
 * in the total force calculation.
 * This estimate is reasonable for recent Intel and AMD x86_64 CPUs.
 */

#ifdef __cplusplus
}
#endif

#endif  /* _perf_est_h */
