/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 *
 *
 * \brief
 * This file contains math functions used by the AWH code.
 *
 * \author Viveca Lindahl
 * \inlibraryapi
 */

#ifndef GMX_AWH_MATH_H
#define GMX_AWH_MATH_H

#include "gromacs/utility/basedefinitions.h"

/*! \brief Generate a sample from a discrete probability distribution defined on [0, ndistr - 1].
 *
 * \param[in] distr      Normalized probability distribution to generate a sample from.
 * \param[in] ndistr     Length of distr.
 * \param[in] step   Time step needed by the random number generator.
 * \param[in] seed   Random seed needed by the random number generator.
 * \param[in] iseed  Second random seed.
 * \returns a sample in [0, ndistr - 1]
 */
int get_sample_from_distribution(const double *distr, int ndistr, gmx_int64_t step, gmx_int64_t seed, int iseed);

/*! \brief Returns the exponent c where exp(c) = exp(a) + exp(b).
 *
 * \param[in] a     First exponent.
 * \param[in] b     Second exponent.
 * \returns c.
 */
double expsum(double a, double b);


/*! \brief
 * Returns an approximation of the geometry factor used for initializing the AWH update size.
 *
 * The geometry factor is defined as the following sum of Gaussians:
 * sum_{k!=0} exp(-0.5*(k*pi*x)^2)/(pi*k)^2,
 * where k is a dim-dimensional integer vector with k_i in {0,1,..}.
 *
 * \param[in] xarray  dim-dimensional array to evaluate.
 * \param[in] dim   The number of dimension.
 * \returns the geometry factor.
 */
double gaussian_geometry_factor(const double *xarray, int ndim);

/*! \brief Returns the the minimum power of 2 exponent k >= 0 such that x <= 2^k.
 * \param[in] x  real value.
 * \returns ceil(log2(x)) for x >= 1 and 0 for x <= 1, i.e. ceil(log2(max(1, x))).
 */
int ceil_log2(double x);

#endif
