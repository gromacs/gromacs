/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2010
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _GMX_INPUTREC_H_
#define _GMX_INPUTREC_H_

/** @file gmx_sort.h
 *
 *  @brief Portable implementation of thread-safe sort routines.
 *
 *
 *  This module provides a Gromacs version of the qsort() routine defined.
 *  It is not highly optimized, but it is thread safe, i.e. multiple threads
 *  can simultaneously call gmx_qsort with different data.
 */

#include <stdlib.h>

#include "types/inputrec.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif



int ir_optimal_nstcalcenergy(const t_inputrec *ir);

int tcouple_min_integration_steps(int etc);

int ir_optimal_nsttcouple(const t_inputrec *ir);

int pcouple_min_integration_steps(int epc);

int ir_optimal_nstpcouple(const t_inputrec *ir);

#ifdef __cplusplus
}
#endif


#endif /* _GMX_INPUTREC_H_ */
