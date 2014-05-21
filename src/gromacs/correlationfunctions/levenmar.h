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
 * Declares routines for non-linear fitting a data set to a function
 *
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_LEVENMAR_H
#define GMX_LEVENMAR_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Do non-linear curve fitting of user supplied function to a data set
 * using the Levenberg Marquardt algorithm.
 *
 * The Levenberg-Marquardt method attempts to reduce the value Chi^2
 * of a fit between a set of data points x[1..ndata], y[1..ndata]
 * with individual standard deviations sig[1..ndata], and a nonlinear
 * function dependent on ma coefficients a[1..ma]. The input array
 * ia[1..ma] indicates by nonzero entries those components of a that
 * should be fitted for, and by zero entries those components that
 * should be held fixed at their input values. The program returns
 * current best-fit values for the parameters a[1..ma], and
 * Chi^2 = chisq. The arrays covar[1..ma][1..ma], alpha[1..ma][1..ma]
 * are used as working space during most iterations. Supply a routine
 * funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit,
 * and its derivatives dyda[1..ma] with respect to the fitting
 * parameters a at x. On the first call provide an initial guess for
 * the parameters a, and set alamda < 0 for initialization (which then
 * sets alamda=.001). If a step succeeds chisq becomes smaller and
 * alamda de-creases by a factor of 10. If a step fails alamda grows by
 * a factor of 10. You must call this routine repeatedly until
 * convergence is achieved. Then, make one final call with alamda=0,
 * so that covar[1..ma][1..ma] returns the covariance matrix, and alpha
 * the curvature matrix.
 * (Parameters held fixed will return zero covariances.)
 *
 * \param[in] x The x-axis (time series)
 * \param[in] y The y-axis (the data)
 * \param[in] sig The uncertainty in y
 * \param[in] ndata The number of data points
 * \param[in] a
 * \param[in] ia
 * \param[in] ma
 * \param[out] covar The covariance matrix
 * \param[in] alpha
 * \param[out] chisq The Chi^2 value of the final fit
 * \param[in] funcs A C-function to be called to fit the data
 * \param[out] alamda
 */
gmx_bool mrqmin(real x[], real y[], real sig[], int ndata, real a[],
                int ia[], int ma, real **covar, real **alpha, real *chisq,
                void (*funcs)(real, real [], real *, real []),
                real *alamda);

#ifdef __cplusplus
}
#endif

#endif
