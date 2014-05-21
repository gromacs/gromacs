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
/*! \libinternal
 * \file
 * \brief
 * Declares routines for integrating a data set
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_INTEGRATE_H
#define GMX_INTEGRATE_H

#include <stdio.h>

#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Integrate the equispaced data in c[] from 0 to n using trapezium rule.
 * If fit != NULL the fit is written as well.
 * \param[in] fp File pointer to write to (maybe NULL)
 * \param[in] n Number of data points
 * \param[in] dt The time step between data points
 * \param[in] c The data set
 * \param[in] fit Fit to the function that is printed too if not a NULL pointer is passed.
 * \param[in] nskip Determines whether all elements are written to the output file
 * (written when i % nskip == 0)
 * \return The integral
 */
real print_and_integrate(FILE *fp, int n, real dt, const real c[], const real *fit, int nskip);

/*! \brief
 * Integrate data in y using the trapezium rule, and, if given, use dy as weighting
 *
 * \param[in] n The number of data points
 * \param[in] x The x coordinate
 * \param[in] y The y data (function values)
 * \param[in] dy The uncertainties (can be NULL)
 * \param[in] aver_start should be set to a value where the function has
 * converged to 0.
 * \param[out] stddev The standard deviation in the integral
 * \return the integral
 */
real evaluate_integral(int n, const real x[], const real y[], const real dy[], real aver_start,
                       real *stddev);

#ifdef __cplusplus
}
#endif

#endif
