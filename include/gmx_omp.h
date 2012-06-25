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

/* This module defines wrappers for OpenMP API functions and enables compiling
 * code even when OpenMP is turned off in the build system.
 * Therfore, OpenMP API functions should always be used through these wrappers
 * and omp.h should never be directly included. Instead, this header should be
 * used whnever OpenMP API functions are needed.
 */

/*! Sets the number of threads in subsequent parallel regions, unless overridden
 *  by a num_threads clause. Acts as a wrapper for omp_get_max_threads(void). */
int  gmx_omp_get_max_threads(void);

/*! Returns the thread number of the thread executing within its thread team.
 *  Acts as a warpper for omp_get_thread_num(void). */
int  gmx_omp_get_thread_num(void);

/*! Returns an integer that is equal to or greater than the number of threads
 * that would be available if a parallel region without num_threads were
 * defined at that point in the code. Acts as a wapepr for omp_set_num_threads(void). */
void gmx_omp_set_num_threads(int num_threads);

#endif /* GMX_OMP_H */
