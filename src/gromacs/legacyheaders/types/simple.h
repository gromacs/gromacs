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

#include "../../math/vectypes.h"
#include "../../utility/basedefinitions.h"
#include "../../utility/real.h"

typedef int         atom_id;      /* To indicate an atoms id         */
#define NO_ATID     (atom_id)(~0) /* Use this to indicate invalid atid */

#ifdef __PGI
/* The portland group x86 C/C++ compilers do not treat negative zero initializers
 * correctly, but "optimizes" them to positive zero, so we implement it explicitly.
 * These constructs are optimized to simple loads at compile time. If you want to
 * use them on other compilers those have to support gcc preprocessor extensions.
 * Note: These initializers might be sensitive to the endianness (which can
 * be different for byte and word order), so check that it works for your platform
 * and add a separate section if necessary before adding to the ifdef above.
 */
#    define GMX_DOUBLE_NEGZERO  ({ const union { int  di[2]; double d; } _gmx_dzero = {0, -2147483648}; _gmx_dzero.d; })
#    define GMX_FLOAT_NEGZERO   ({ const union { int  fi; float f; } _gmx_fzero = {-2147483648}; _gmx_fzero.f; })
#else
/*! \brief Negative zero in double */
#    define GMX_DOUBLE_NEGZERO  (-0.0)

/*! \brief Negative zero in float */
#    define GMX_FLOAT_NEGZERO   (-0.0f)
#endif

#endif
