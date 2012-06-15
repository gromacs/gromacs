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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifdef GMX_OPENMP
#include <omp.h>
#endif

#include "gmx_omp.h"

int gmx_omp_get_max_threads(void)
{
#ifdef GMX_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}


int gmx_omp_get_thread_num(void)
{
#ifdef GMX_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

void gmx_omp_set_num_threads(int num_threads)
{
#ifdef GMX_OPENMP
    return omp_set_num_threads(num_threads);
#else
    return;
#endif
}
