/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * This file is part of Gromacs        Copyright (c) 1991-2010
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
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
#include "visibility.h"
#include "types/inputrec.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif



GMX_LIBGMX_EXPORT
int ir_optimal_nstcalcenergy(const t_inputrec *ir);

GMX_LIBGMX_EXPORT
int tcouple_min_integration_steps(int etc);

GMX_LIBGMX_EXPORT
int ir_optimal_nsttcouple(const t_inputrec *ir);

GMX_LIBGMX_EXPORT
int pcouple_min_integration_steps(int epc);

GMX_LIBGMX_EXPORT
int ir_optimal_nstpcouple(const t_inputrec *ir);

#ifdef __cplusplus
}
#endif


#endif /* _GMX_INPUTREC_H_ */
