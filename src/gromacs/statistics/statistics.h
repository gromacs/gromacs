/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef GMX_STATISTICS_H
#define GMX_STATISTICS_H

#include <stdio.h>

#include "gromacs/utility/real.h"

/*! \libinternal \file
 *
 * \brief
 * Declares simple statistics toolbox
 *
 * \authors David van der Spoel <david.vanderspoel@icm.uu.se>
 *
 * \inlibraryapi
 */
#ifdef __cplusplus
extern "C" {
#endif

//! Abstract container type
typedef struct gmx_stats *gmx_stats_t;

//! Error codes returned by the routines
enum {
    estatsOK, estatsNO_POINTS, estatsNO_MEMORY, estatsERROR,
    estatsINVALID_INPUT, estatsNOT_IMPLEMENTED, estatsNR
};

//! Enum for statistical weights
enum {
    elsqWEIGHT_NONE, elsqWEIGHT_X, elsqWEIGHT_Y,
    elsqWEIGHT_XY, elsqWEIGHT_NR
};

//! Enum determining which coordinate to histogram
enum {
    ehistoX, ehistoY, ehistoNR
};

/*! \brief
 * Initiate a data structure
 * \return the data structure
 */
gmx_stats_t gmx_stats_init();

/*! \brief
 * Destroy a data structure
 * \param stats The data structure
 * \return error code
 */
int gmx_stats_done(gmx_stats_t stats);

/*! \brief
 * Remove outliers from a straight line, where level in units of
 * sigma. Level needs to be larger than one obviously.
 * \param[in] stats The data structure
 * \param[in] level The sigma level
 * \return error code
 */
int gmx_stats_remove_outliers(gmx_stats_t stats, double level);

/*! \brief
 * Add a point to the data set
 * \param[in] stats The data structure
 * \param[in] x   The x value
 * \param[in] y   The y value
 * \param[in] dx  The error in the x value
 * \param[in] dy  The error in the y value
 * \return error code
 */
int gmx_stats_add_point(gmx_stats_t stats, double x, double y,
                        double dx, double dy);

/*! \brief
 * Add a series of datapoints at once. The arrays dx and dy may
 * be NULL in that case zero uncertainties will be assumed.
 *
 * \param[in] stats The data structure
 * \param[in] n   Number of points
 * \param[in] x   The array of x values
 * \param[in] y   The array of y values
 * \param[in] dx  The error in the x value
 * \param[in] dy  The error in the y value
 * \return error code
 */
int gmx_stats_add_points(gmx_stats_t stats, int n, real *x, real *y,
                         real *dx, real *dy);

/*! \brief
 * Delivers data points from the statistics.
 *
 * Should be used in a while loop. Variables for either
 * pointer may be NULL, in which case the routine can be used as an
 * expensive point counter.
 * Return the data points one by one. Return estatsOK while there are
 *  more points, and returns estatsNOPOINTS when the last point has
 *  been returned.
 *  If level > 0 then the outliers outside level*sigma are reported
 * only.
 * \param[in] stats The data structure
 * \param[out] x   The array of x values
 * \param[out] y   The array of y values
 * \param[out] dx  The error in the x value
 * \param[out] dy  The error in the y value
 * \param[in]  level sigma level (see above)
 * \return error code
 */
int gmx_stats_get_point(gmx_stats_t stats, real *x, real *y,
                        real *dx, real *dy, real level);

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
 * \return error code
 */
int gmx_stats_get_ab(gmx_stats_t stats, int weight,
                     real *a, real *b,
                     real *da, real *db, real *chi2, real *Rfit);

/*! \brief
 * Fit the data to y = ax, possibly weighted, if uncertainties have
 * have been input. da and db may be NULL.
 * \param[in] stats The data structure
 * \param[in] weight type of weighting
 * \param[out] a slope
 * \param[out] da sigma in a
 * \param[out] chi2 normalized quality of fit
 * \param[out] Rfit correlation coefficient
 * \return error code
 */
int gmx_stats_get_a(gmx_stats_t stats, int weight,
                    real *a, real *da, real *chi2, real *Rfit);

/*! \brief
 * Get the correlation coefficient.
 * \param[in]  stats The data structure
 * \param[out] R the correlation coefficient between the data (x and y) as input to the structure.
 * \return error code
 */
int gmx_stats_get_corr_coeff(gmx_stats_t stats, real *R);

/*! \brief
 * Get the root mean square deviation.
 * \param[in]  stats The data structure
 * \param[out] rmsd  the root mean square deviation between x and y values.
 * \return error code
 */
int gmx_stats_get_rmsd(gmx_stats_t stats, real *rmsd);

/*! \brief
 * Get the number of points.
 * \param[in]  stats The data structure
 * \param[out] N     number of data points
 * \return error code
 */
int gmx_stats_get_npoints(gmx_stats_t stats, int *N);

/*! \brief
 * Computes and returns the average value.
 * \param[in]  stats The data structure
 * \param[out] aver  Average value
 * \return error code
 */
int gmx_stats_get_average(gmx_stats_t stats, real *aver);

/*! \brief
 * Computes and returns the standard deviation.
 * \param[in]  stats The data structure
 * \param[out] sigma  Standard deviation
 * \return error code
 */
int gmx_stats_get_sigma(gmx_stats_t stats, real *sigma);

/*! \brief
 * Computes and returns the standard error.
 * \param[in]  stats The data structure
 * \param[out] error Standard error
 * \return error code
 */
int gmx_stats_get_error(gmx_stats_t stats, real *error);

/*! \brief
 * Pointers may be null, in which case no assignment will be done.
 * \param[in]  stats The data structure
 * \param[out] aver  Average value
 * \param[out] sigma  Standard deviation
 * \param[out] error Standard error
 * \return error code
 */
int gmx_stats_get_ase(gmx_stats_t stats, real *aver, real *sigma, real *error);

/*! \brief
 * Dump the x, y, dx, dy data to a text file
 * \param[in]  stats The data structure
 * \param[in] fp  File pointer
 * \return error code
 */
int gmx_stats_dump_xy(gmx_stats_t stats, FILE *fp);

/*! \brief
 * Make a histogram of the data present.
 *
 * Uses either binwidth to
 * determine the number of bins, or nbins to determine the binwidth,
 * therefore one of these should be zero, but not the other. If *nbins = 0
 * the number of bins will be returned in this variable. ehisto should be one of
 * ehistoX or ehistoY. If
 * normalized not equal to zero, the integral of the histogram will be
 * normalized to one. The output is in two arrays, *x and *y, to which
 * you should pass a pointer. Memory for the arrays will be allocated
 * as needed. Function returns one of the estats codes.
 * \param[in]  stats The data structure
 * \param[in] binwidth For the histogram
 * \param[in] nbins    Number of bins
 * \param[in] ehisto   Type (see enum above)
 * \param[in] normalized see above
 * \param[out] x see above
 * \param[out] y see above
 * \return error code
 */
int gmx_stats_make_histogram(gmx_stats_t stats, real binwidth, int *nbins,
                             int ehisto,
                             int normalized, real **x, real **y);

/*! \brief
 * Return message belonging to error code
 * \param[in] estats error code
 */
const char *gmx_stats_message(int estats);

/****************************************************
 * Some statistics utilities for convenience: useful when a complete data
 * set is available already from another source, e.g. an xvg file.
 ****************************************************/
/*! \brief
 * Fit a straight line y=ax thru the n data points x, y, return the
 * slope in *a.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[out] a slope
 * \return error code
 */
int lsq_y_ax(int n, real x[], real y[], real *a);

/*! \brief
 * Fit a straight line y=ax+b thru the n data points x, y.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] r correlation coefficient
 * \param[out] chi2 quality of fit
 * \return error code
 */
int lsq_y_ax_b(int n, real x[], real y[], real *a, real *b, real *r,
               real *chi2);

/*! \copydoc lsq_y_ax_b
 */
int lsq_y_ax_b_xdouble(int n, double x[], real y[],
                       real *a, real *b, real *r, real *chi2);

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
 * \return error code
 */
int lsq_y_ax_b_error(int n, real x[], real y[], real dy[],
                     real *a, real *b, real *da, real *db,
                     real *r, real *chi2);

#ifdef __cplusplus
}
#endif

#endif
