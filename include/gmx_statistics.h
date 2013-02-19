/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _GMX_STATS_H
#define _GMX_STATS_H
#include "visibility.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "typedefs.h"

typedef struct gmx_stats *gmx_stats_t;

/* Error codes returned by the routines */
enum {
    estatsOK, estatsNO_POINTS, estatsNO_MEMORY, estatsERROR,
    estatsINVALID_INPUT, estatsNOT_IMPLEMENTED, estatsNR
};

enum {
    elsqWEIGHT_NONE, elsqWEIGHT_X, elsqWEIGHT_Y,
    elsqWEIGHT_XY, elsqWEIGHT_NR
};

enum {
    ehistoX, ehistoY, ehistoNR
};

GMX_LIBGMX_EXPORT
gmx_stats_t gmx_stats_init();

GMX_LIBGMX_EXPORT
int gmx_stats_done(gmx_stats_t stats);

/* Remove outliers from a straight line, where level in units of
   sigma. Level needs to be larger than one obviously. */
int gmx_stats_remove_outliers(gmx_stats_t stats, double level);

GMX_LIBGMX_EXPORT
int gmx_stats_add_point(gmx_stats_t stats, double x, double y,
                        double dx, double dy);

/* The arrays dx and dy may be NULL if no uncertainties are available,
   in that case zero uncertainties will be assumed. */
int gmx_stats_add_points(gmx_stats_t stats, int n, real *x, real *y,
                         real *dx, real *dy);

/* Return the data points one by one. Return estatsOK while there are
   more points, and returns estatsNOPOINTS when the last point has
   been returned. Should be used in a while loop. Variables for either
   pointer may be NULL, in which case the routine can be used as an
   expensive point counter. */
GMX_LIBGMX_EXPORT
int gmx_stats_get_point(gmx_stats_t stats, real *x, real *y,
                        real *dx, real *dy);

/* Fit the data to y = ax + b, possibly weighted, if uncertainties
   have been input. Returns slope in *a and intercept in b, *return
   sigmas in *da and *db respectively. Returns normalized *quality of
   fit in *chi2 and correlation of fit with data in Rfit. chi2, Rfit,
   da and db may be NULL. */
GMX_LIBGMX_EXPORT
int gmx_stats_get_ab(gmx_stats_t stats, int weight,
                     real *a, real *b,
                     real *da, real *db, real *chi2, real *Rfit);

/* Fit the data to y = ax, possibly weighted, if uncertainties have
   been input. Returns slope in *a, sigma in a in *da, and normalized
   quality of fit in *chi2 and correlation of fit with data in
   Rfit. chi2, Rfit and da may be NULL. */
int gmx_stats_get_a(gmx_stats_t stats, int weight,
                    real *a, real *da, real *chi2, real *Rfit);

/* Return the correlation coefficient between the data (x and y) as
   input to the structure. */
int gmx_stats_get_corr_coeff(gmx_stats_t stats, real *R);

/* Returns the root mean square deviation between x and y values. */
int gmx_stats_get_rmsd(gmx_stats_t gstats, real *rmsd);

GMX_LIBGMX_EXPORT
int gmx_stats_get_npoints(gmx_stats_t stats, int *N);

GMX_LIBGMX_EXPORT
int gmx_stats_get_average(gmx_stats_t stats, real *aver);

int gmx_stats_get_sigma(gmx_stats_t stats, real *sigma);

int gmx_stats_get_error(gmx_stats_t stats, real *error);

/* Get all three of the above. Pointers may be null, in which case no
   assignment will be done. */
GMX_LIBGMX_EXPORT
int gmx_stats_get_ase(gmx_stats_t gstats, real *aver, real *sigma, real *error);

/* Dump the x, y, dx, dy data to a text file */
int gmx_stats_dump_xy(gmx_stats_t gstats, FILE *fp);

/* Make a histogram of the data present. Uses either bindwith to
   determine the number of bins, or nbins to determine the binwidth,
   therefore one of these should be zero, but not the other. If *nbins = 0
   the number of bins will be returned in this variable. ehisto should be one of
   ehistoX or ehistoY. If
   normalized not equal to zero, the integral of the histogram will be
   normalized to one. The output is in two arrays, *x and *y, to which
   you should pass a pointer. Memory for the arrays will be allocated
   as needed. Function returns one of the estats codes. */
int gmx_stats_make_histogram(gmx_stats_t gstats, real binwidth, int *nbins,
                             int ehisto,
                             int normalized, real **x, real **y);

/* Return message belonging to error code */
GMX_LIBGMX_EXPORT
const char *gmx_stats_message(int estats);

/****************************************************
 * Some statistics utilities for convenience: useful when a complete data
 * set is available already from another source, e.g. an xvg file.
 ****************************************************/
int lsq_y_ax(int n, real x[], real y[], real *a);
/* Fit a straight line y=ax thru the n data points x, y, return the
   slope in *a. Return value can be estatsOK, or something else. */

GMX_LIBGMX_EXPORT
int lsq_y_ax_b(int n, real x[], real y[], real *a, real *b, real *r,
               real *chi2);
/* Fit a straight line y=ax+b thru the n data points x,y.
 * Returns the "fit quality" sigma = sqrt(chi^2/(n-2)).
 * The correlation coefficient is returned in r.
 */

int lsq_y_ax_b_xdouble(int n, double x[], real y[],
                       real *a, real *b, real *r, real *chi2);
/* As lsq_y_ax_b, but with x in double precision.
 */

GMX_LIBGMX_EXPORT
int lsq_y_ax_b_error(int n, real x[], real y[], real dy[],
                     real *a, real *b, real *da, real *db,
                     real *r, real *chi2);
/* Fit a straight line y=ax+b thru the n data points x,y, with sigma dy
 * Returns the "fit quality" sigma = sqrt(chi^2/(n-2)).
 * The correlation coefficient is returned in r.
 */

#ifdef __cplusplus
}
#endif

#endif
