/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef _sortwater_h
#define _sortwater_h

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/random/random.h"

#ifdef __cplusplus
extern "C" {
#endif

void randwater(int astart, int nwater, int nwatom,
               rvec x[], rvec v[], gmx_rng_t rng);
/* Randomize the order of nwater molecules of length nwatom, the
 * first atom of which is at astart.
 * If v is not NULL it will be shuffled along
 * IS NOT THREAD SAFE
 */


void sortwater(int astart, int nwater, int nwatom, rvec x[], rvec v[]);
/* Sort the order of nwater molecules of length nwatom on X coordinate
 * If v is not NULL it will be shuffled along
 * IS NOT THREAD SAFE
 */

void mkcompact(int astart, int nwater, int nwatom, rvec x[], rvec v[],
               int nnode, matrix box);
/* Make compact subboxes
 * IS NOT THREAD SAFE  */

#ifdef __cplusplus
}
#endif

#endif
