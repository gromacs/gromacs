/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#include <assert.h>

#include "gromacs/utility/basedefinitions.h"

/* We only use the C interface of ThreeFry and r123array2x64. This file is a
   replacment for the original from the Random123 distribution. It sets all
   defines (they all start with R123_), which are used by those parts of
   Random123 being used. Instead of determining values based on the Compiler
   (Name, Version, ..) we set values based on our assumptions in Gromacs, and
   the defines used in Gromacs */

/* Random123 isn't used from Cuda thus this can always be empty */
#define R123_CUDA_DEVICE
/* For "inline" use the Gromacs own gmx_inline */
#define R123_STATIC_INLINE static gmx_inline
/* force_inline isn't used in Gromacs - if it matters for a compiler it probably
   not only matters here and should be defined in basedefinitions.h */
#define R123_FORCE_INLINE(decl) decl
/* We assume in Gromacs that assert is available outside of Cuda */
#define R123_ASSERT assert
/* Not used  (only used by C++ interface of ThreeFry) */
#define R123_STATIC_ASSERT(expr, msg)
