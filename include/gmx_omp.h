/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * GROup of MAchos and Cynical Suckers
 */

#ifndef GMX_OMP_H
#define GMX_OMP_H

#include "types/commrec.h"
#include "mdrun.h"

/* This module defines wrappers for OpenMP API functions and enables compiling
 * code even when OpenMP is turned off in the build system.
 * Therefore, OpenMP API functions should always be used through these wrappers
 * and omp.h should never be directly included. Instead, this header should be
 * used whenever OpenMP API functions are needed.
 */

/*! Returns an integer equal to or greater than the number of threads
 *  that would be available if a parallel region without num_threads were
 *  defined at that point in the code. Acts as a wrapper for omp_set_num_threads(void). */
int  gmx_omp_get_max_threads(void);

/*! Returns the number of processors available when the function is called.
 *  Acts as a wrapper around omp_get_num_procs() */
int gmx_omp_get_num_procs(void);

/*! Returns the thread number of the thread executing within its thread team.
 *  Acts as a wrapper for omp_get_thread_num(void). */
GMX_LIBGMX_EXPORT int  gmx_omp_get_thread_num(void);

/*! Sets the number of threads in subsequent parallel regions, unless overridden
 *  by a num_threads clause. Acts as a wrapper for omp_get_max_threads(void). */
void gmx_omp_set_num_threads(int num_threads);

/*! Check for externally set thread affinity to avoid conflicts with GROMACS internal setting. */
GMX_LIBGMX_EXPORT void gmx_omp_check_thread_affinity(FILE *fplog, const t_commrec *cr,
                                   gmx_hw_opt_t *hw_opt);

#endif /* GMX_OMP_H */
