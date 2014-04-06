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

#ifndef _simple_h
#define _simple_h

#include "../../utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#define XX      0                 /* Defines for indexing in */
#define YY      1                 /* vectors                 */
#define ZZ      2
#define DIM     3                 /* Dimension of vectors    */

typedef int         atom_id;      /* To indicate an atoms id         */
#define NO_ATID     (atom_id)(~0) /* Use this to indicate invalid atid */

/*! \brief Double precision accuracy */
#define GMX_DOUBLE_EPS   1.11022302E-16

/*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define GMX_DOUBLE_MAX   1.79769312E+308

/*! \brief Minimum double precision value */
#define GMX_DOUBLE_MIN   2.22507386E-308

/*! \brief Single precision accuracy */
#define GMX_FLOAT_EPS    5.96046448E-08

/*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define GMX_FLOAT_MAX    3.40282346E+38

/*! \brief Minimum single precision value */
#define GMX_FLOAT_MIN    1.17549435E-38


/* Check whether we already have a real type! */
#ifdef GMX_DOUBLE

#ifndef HAVE_REAL
typedef double      real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL    MPI_DOUBLE
#define GMX_REAL_EPS    GMX_DOUBLE_EPS
#define GMX_REAL_MIN    GMX_DOUBLE_MIN
#define GMX_REAL_MAX    GMX_DOUBLE_MAX
#define gmx_real_fullprecision_pfmt "%21.14e"
#else

#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL    MPI_FLOAT
#define GMX_REAL_EPS    GMX_FLOAT_EPS
#define GMX_REAL_MIN    GMX_FLOAT_MIN
#define GMX_REAL_MAX    GMX_FLOAT_MAX
#define gmx_real_fullprecision_pfmt "%14.7e"
#endif

typedef real            rvec[DIM];

typedef double          dvec[DIM];

typedef real            matrix[DIM][DIM];

typedef real            tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];

/* Standard sizes for char* string buffers */
#define STRLEN 4096

#ifdef __cplusplus
}
#endif

#endif
