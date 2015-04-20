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
 * Declares routine for computing a cross correlation between two data sets
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_CROSSCORR_H
#define GMX_CROSSCORR_H

#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * fft cross correlation algorithm.
 * Computes corr = f (.) g
 *
 * \param[in] n number of data point
 * \param[in] f First function
 * \param[in] g Second function
 * \param[out] corr The cross correlation
 */
void cross_corr(int n, real f[], real g[], real corr[]);

/*! \brief
 * fft cross correlation algorithm.
 *
 * Computes corr[n] = f[n][i] (.) g[n][i], that is for nFunc
 * pairs of arrays n the cross correlation is computed in parallel
 * using OpenMP.
 *
 * \param[in] nFunc nuber of function to crosscorrelate
 * \param[in] nData number of data point in eatch function
 * \param[in] f 2D array of first function to crosscorrelate
 * \param[in] g 2D array of second function to crosscorrelate
 * \param[out] corr 2D array of the cross correlations
 */
void many_cross_corr(int nFunc, int * nData, real ** f, real ** g, real ** corr);

#ifdef __cplusplus
}
#endif

#endif
