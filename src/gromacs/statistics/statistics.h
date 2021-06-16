/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2019,2021, by the GROMACS development team, led by
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
 * Declares simple statistics toolbox
 *
 * \authors David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 */
#ifndef GMX_STATISTICS_H
#define GMX_STATISTICS_H

#include <cstdio>

#include <tuple>

#include "gromacs/utility/real.h"

//! Abstract container type
typedef struct gmx_stats* gmx_stats_t;

//! Enum for statistical weights
enum
{
    elsqWEIGHT_NONE,
    elsqWEIGHT_X,
    elsqWEIGHT_Y,
    elsqWEIGHT_XY,
    elsqWEIGHT_NR
};

/*! \brief
 * Initiate a data structure
 * \return the data structure
 */
gmx_stats_t gmx_stats_init();

/*! \brief
 * Destroy a data structure
 * \param stats The data structure
 */
void gmx_stats_free(gmx_stats_t stats);

/*! \brief
 * Add a point to the data set
 * \param[in] stats The data structure
 * \param[in] x   The x value
 * \param[in] y   The y value
 * \param[in] dx  The error in the x value
 * \param[in] dy  The error in the y value
 */
void gmx_stats_add_point(gmx_stats_t stats, double x, double y, double dx, double dy);

/*! \brief
 * Fit the data to y = ax + b, possibly weighted, if uncertainties
 * have been input. da and db may be NULL.
 * \param[in] stats The data structure
 * \param[in] weight type of weighting
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] da sigma in a
 * \param[out] db sigma in b
 * \param[out] chi2 normalized quality of fit
 * \param[out] Rfit correlation coefficient
 */
void gmx_stats_get_ab(gmx_stats_t stats, int weight, real* a, real* b, real* da, real* db, real* chi2, real* Rfit);

/*! \brief
 * Computes and returns the average value.
 * \param[in]  stats The data structure
 * \return Average value
 * \throws  InconsistentInputError if given no points to average
 */
real gmx_stats_get_average(gmx_stats_t stats);

/*! \brief
 * Pointers may be null, in which case no assignment will be done.
 * \param[in]  stats The data structure
 * \return Tuple of (average value, its standard deviation, its standard error)
 * \throws  InconsistentInputError if given no points to analyze
 */
std::tuple<real, real, real> gmx_stats_get_ase(gmx_stats_t stats);

/****************************************************
 * Some statistics utilities for convenience: useful when a complete data
 * set is available already from another source, e.g. an xvg file.
 ****************************************************/

/*! \brief
 * Fit a straight line y=ax+b thru the n data points x, y.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] r correlation coefficient
 * \param[out] chi2 quality of fit
 *
 * \throws  InconsistentInputError if given no points to fit
 */
void lsq_y_ax_b(int n, real x[], real y[], real* a, real* b, real* r, real* chi2);

/*! \copydoc lsq_y_ax_b
 * Suits cases where x is already always computed in double precision
 * even in a mixed-precision build configuration.
 */
void lsq_y_ax_b_xdouble(int n, double x[], real y[], real* a, real* b, real* r, real* chi2);

/*! \brief
 * Fit a straight line y=ax+b thru the n data points x, y.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[in] dy uncertainty in data point y
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] da error in slope
 * \param[out] db error in intercept
 * \param[out] r correlation coefficient
 * \param[out] chi2 quality of fit
 *
 * \throws  InconsistentInputError if given no points to fit
 */
void lsq_y_ax_b_error(int n, real x[], real y[], real dy[], real* a, real* b, real* da, real* db, real* r, real* chi2);

#endif
