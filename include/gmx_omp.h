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

/* This header file should always be included whenever OpenMP API functions are
 * needed. It is meant enable compiling code even when OpenMP is turned off in
 * the build system (and without #ifdef-decoration).
 */
#ifdef GMX_OPENMP
#include <omp.h>
#else
static inline int  omp_get_max_threads() { return 1; }
static inline int  omp_get_thread_num()  { return 0; }
static inline void omp_set_num_threads(int dummy) { return; }
#endif /* GMX_OPENMP */

#endif /* GMX_OMP_H */
