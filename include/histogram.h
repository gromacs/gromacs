/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
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
/*! \file
 * \brief API for calculation of histograms with error estimates.
 *
 * The API is documented in more detail on a separate page:
 * \ref histograms
 *
 * The functions within this file can be used and developed independently of
 * the other parts of the library.
 * Other parts of the library do not reference these functions.
 */
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Type of histogram. */
typedef enum
{
    /*! \brief
     * Simple histogram.
     *
     * Use gmx_histogram_increment() or gmx_histogram_increment_bin()
     * to sample.
     */
    HIST_SIMPLE,
    /*! \brief
     * Weighted histogram where different points contribute different amounts.
     *
     * Use gmx_histogram_add() or gmx_histogram_add_to_bin() to sample.
     */
    HIST_WEIGHT,
    /*! \brief
     * Calculate averages within each bin.
     *
     * Use gmx_histogram_add_item() or gmx_histogram_add_item_to_bin() to sample.
     */
    HIST_BINAVER
} e_histogram_t;

/** Whether bins are centered at integer values. */
#define HIST_INTEGERBINS  1
/** Whether the values outside the range should be included in the histogram. */
#define HIST_ALL          2
/** Whether the initialization used binwidths. */
#define HIST_INITBW       128

/** Stores data for a histogram. */
typedef struct gmx_histogram_t gmx_histogram_t;

/*! \name Initialization functions
 */
/*@{*/
/** Initialize calculation of a histogram. */
int
gmx_histogram_create(gmx_histogram_t **h, e_histogram_t type, int nbins);
/** Initialize calculation of a histogram for a range. */
int
gmx_histogram_create_range(gmx_histogram_t **h, e_histogram_t type,
                           real start, real end, real binw, gmx_bool bIntegerBins);
/** Clears the bins in the histogram. */
void
gmx_histogram_clear(gmx_histogram_t *h);
/** Frees the memory allocated for a histogram. */
void
gmx_histogram_free(gmx_histogram_t *h);
/** Sets histogram range using a starting point and a bin width. */
int
gmx_histogram_set_binwidth(gmx_histogram_t *h, real start, real binw);
/** Sets histogram range using endpoint values. */
int
gmx_histogram_set_range(gmx_histogram_t *h, real start, real end);
/** Sets histogram bins to center at integer values. */
void
gmx_histogram_set_integerbins(gmx_histogram_t *h, gmx_bool bIntegerBins);
/** Sets histogram to include outlying values in the bins at the edges. */
void
gmx_histogram_set_all(gmx_histogram_t *h, gmx_bool bAll);
/** Sets block size for histogram averaging. */
int
gmx_histogram_set_blocksize(gmx_histogram_t *h, int bsize);
/** Sets output file for block histograms. */
int
gmx_histogram_set_block_output(gmx_histogram_t *h, FILE *fp);
/*@}*/

/*! \name Access functions
 */
/*@{*/
/** Finds the histogram bin corresponding to a value. */
int
gmx_histogram_find_bin(gmx_histogram_t *h, real pos);
/** Returns the number of bins in a histogram. */
int
gmx_histogram_get_nbins(gmx_histogram_t *h);
/** Returns the bin width of a histogram. */
real
gmx_histogram_get_binwidth(gmx_histogram_t *h);
/** Returns the value of the histogram at a certain position. */
void
gmx_histogram_get_value(gmx_histogram_t *h, real pos, double *val, double *err);
/** Returns the value of the histogram in a certain bin. */
void
gmx_histogram_get_bin_value(gmx_histogram_t *h, int bin, double *val, double *err);
/** Returns an array of values for the histogram. */
double *
gmx_histogram_get_values(gmx_histogram_t *h);
/** Returns an array of error values for the histogram. */
double *
gmx_histogram_get_errors(gmx_histogram_t *h);
/*@}*/

/*! \name Sampling functions
 */
/*@{*/
/** Increments the count in a histogram bin corresponding to \p pos. */
void
gmx_histogram_increment(gmx_histogram_t *h, real pos);
/** Increments the count in a histogram bin. */
void
gmx_histogram_increment_bin(gmx_histogram_t *h, int bin);

/** Adds a value to a histogram bin corresponding to \p pos. */
void
gmx_histogram_add(gmx_histogram_t *h, real pos, double value);
/** Adds a value to a histogram bin. */
void
gmx_histogram_add_to_bin(gmx_histogram_t *h, int bin, double value);

/** Adds a value to a histogram bin corresponding to \p pos. */
void
gmx_histogram_add_item(gmx_histogram_t *h, real pos, double value);
/** Adds a value to a histogram bin. */
void
gmx_histogram_add_item_to_bin(gmx_histogram_t *h, int bin, double value);

/** Finishes histogram sampling for a frame. */
void
gmx_histogram_finish_frame(gmx_histogram_t *h);
/** Normalizes a histogram. */
void
gmx_histogram_finish(gmx_histogram_t *h);
/*@}*/


/*! \name Post-processing functions
 */
/*@{*/
/** Creates a new histogram with double the binwidth. */
void
gmx_histogram_resample_dblbw(gmx_histogram_t **dest, gmx_histogram_t *src,
                             gmx_bool bIntegerBins);
/** Makes a clone of a histogram. */
void
gmx_histogram_clone(gmx_histogram_t **dest, gmx_histogram_t *src);
/** Normalizes a histogram to a probability distribution. */
void
gmx_histogram_normalize_prob(gmx_histogram_t *h);
/** Scales a histogram with a custom normalization factor. */
void
gmx_histogram_scale(gmx_histogram_t *h, real norm);
/** Scales a histogram with a custom non-uniform normalization factor. */
void
gmx_histogram_scale_vec(gmx_histogram_t *h, real norm[]);
/** Writes a single histogram to a file. */
void
gmx_histogram_write(FILE *fp, gmx_histogram_t *h, gmx_bool bErrors);
/** Writes a set of histograms to a file. */
void
gmx_histogram_write_array(FILE *fp, int n, gmx_histogram_t *h[],
                          gmx_bool bValue, gmx_bool bErrors);
/** Writes a set of cumulative histograms to a file. */
void
gmx_histogram_write_cum_array(FILE *fp, int n, gmx_histogram_t *h[],
                              gmx_bool bValue, gmx_bool bErrors);
/*@}*/

#ifdef __cplusplus
}
#endif

#endif
