/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2017,2018, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Implements partial least squares regression
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 * \inlibraryapi
 */

#ifndef GMX_LINEARALGEBRA_PLS_H
#define GMX_LINEARALGEBRA_PLS_H

#include "gromacs/simd/simd_setup.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/data_structures/irreg_array_2d.h"

/** \brief The matrix datatype to be used in pls_denham */
using Pls2dArray = gmx::IrregArray2D<real, typename gmx::SimdSetup<real>::allocator_type>;
/** \brief The vector datatype to be used in pls_denham */
using Pls1dArray = std::vector<real, typename gmx::SimdSetup<real>::allocator_type>;

/*! \brief Partial least squares (PLS) regression
 *
 * Implemented as described in:
 *
 * Denham, M. C. "Implementing Partial Least Squares."
 * Statistics and Computing 5, no. 3 (September 1995): 191–202. doi:10.1007/BF00142661.
 *
 * The algoritm aims to find a vector t to minimize
 *
 * abs(y-X*t)
 *
 * where the vector t=w*q is a linear combination of vectors.
 *
 * \param [in]  X  contains the centered (n,p)-matrix X.
 * \param [in]  y  contains the centered (n)-vector y.
 * \param [in]  n  the number of rows of the matrix X to use in analysis.
 * \param [in]  p  the number of columns of the matrix X to use in analysis.
 * \param [in]  k  the number of PLS factors to include in the regression of y on the matrix X.
 * \param [out] W  a (k,p)-matrix containing the coefficient vectors stored by column of
 *                 the a PLS regression factors obtained from the matrix X.
 * \param [out] q  a (k)-vector containing the least squares regression coefficients of
 *                 the ordinary least squares regression of y on the PLS factor matrix t.
 *
 */
ptrdiff_t pls_denham(const Pls2dArray &X,
                     const Pls1dArray &y,
                     const ptrdiff_t n, const ptrdiff_t p, const ptrdiff_t k,
                     Pls2dArray &W,
                     Pls1dArray &q);
#endif
