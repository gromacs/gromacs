/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef GMX_LEGACYHEADERS_TYPES_FORCETABLE_H
#define GMX_LEGACYHEADERS_TYPES_FORCETABLE_H

#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Structure describing the data in a single table */
typedef struct
{
    enum gmx_table_interaction  interaction; /* Types of interactions stored in this table */
    enum gmx_table_format       format;      /* Interpolation type and data format */

    real                        r;           /* range of the table */
    int                         n;           /* n+1 is the number of table points */
    real                        scale;       /* distance (nm) between two table points */
    real                        scale_exp;   /* distance for exponential part of VdW table, not always used */
    real *                      data;        /* the actual table data */

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it makes it
     * much easier to access the tables in the nonbonded kernels when we can set the data from variables.
     * It is always true that stride = formatsize*ninteractions
     */
    int                         formatsize;    /* Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.) */
    int                         ninteractions; /* Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp. */
    int                         stride;        /* Distance to next table point (number of fp variables per table point in total) */
} t_forcetable;

#ifdef __cplusplus
}
#endif
#endif
