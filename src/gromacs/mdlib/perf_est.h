/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MDLIB_PERF_EST_H
#define GMX_MDLIB_PERF_EST_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/vectypes.h"

struct gmx_mtop_t;
struct t_inputrec;

void count_bonded_distances(const gmx_mtop_t& mtop, const t_inputrec& ir, double* ndistance_c, double* ndistance_simd);
/* Count the number of distance calculations in bonded interactions,
 * separately for plain-C and SIMD bonded functions.
 * The computational cost is nearly proportional to the numbers.
 * It is allowed to pass NULL for the last two arguments.
 */

float pme_load_estimate(const gmx_mtop_t& mtop, const t_inputrec& ir, const matrix box);
/* Returns an estimate for the relative load of the PME mesh calculation
 * in the total force calculation.
 * This estimate is reasonable for recent Intel and AMD x86_64 CPUs.
 */

#endif
