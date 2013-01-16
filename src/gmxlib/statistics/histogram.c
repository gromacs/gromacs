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
/*! \page histograms Histogram calculation
 *
 * Functions to calculate histograms are defined in histogram.h.
 * The principle is simple: histogram calculation is set up with
 * gmx_histogram_create() and gmx_histogram_set_*(), the histogram is sampled
 * using the provided functions, and gmx_histogram_finish() is called to
 * finish the sampling.
 * Various post-processing functions can then be used to normalize the
 * histogram in various ways and to write it into a file.
 * gmx_histogram_get_*() can be used to access data from the histogram.
 * gmx_histogram_free() can be used to free the memory allocated for a
 * histogram.
 *
 *
 * \section histogram_sampling Initialization and sampling
 *
 * A histogram calculation is initialized by calling gmx_histogram_create()
 * or gmx_histogram_create_range().
 * The first function takes the type of the histogram (see \ref e_histogram_t)
 * and the number of bins, and allocates the necessary memory.
 * The bin locations can be initialized with gmx_histogram_set_integerbins()
 * and one of gmx_histogram_set_binwidth() and gmx_histogram_set_range().
 * Only uniformly spaced bins are currently supported.
 * The second initialization function takes the beginning and end of the range
 * and the bin width.
 * The treatment of values that fall outside the histogram range can be
 * set with gmx_histogram_set_all().
 * gmx_histogram_set_blocksize() can be used to specify a block size for error
 * estimates. If not called, no error estimates are calculated.
 * Otherwise, a histogram is calculated for each block of given number of
 * frames, and the error estimated as the error of the mean of these block
 * histograms.
 * If the number of frames is not divisible by the block size,
 * the last block is discarded.
 * gmx_histogram_set_block_output() can be used to write the individual
 * block histograms to a file (it is currently not very flexible).
 *
 * After initialization, three sets of sampling functions are provided
 * to sample the differen types of histograms:
 *  - Use gmx_histogram_increment() or gmx_histogram_increment_bin() to
 *    calculate simple histograms (\ref HIST_SIMPLE histograms).
 *  - Use gmx_histogram_add() or gmx_histogram_add_to_bin() to calculate
 *    histograms where different points can have different weight
 *    (\ref HIST_WEIGHT histograms).
 *  - Use gmx_histogram_add_item() or gmx_histogram_add_item_to_bin() to
 *    calculate averages within each bin
 *    (\ref HIST_BINAVER histograms).
 *
 * The functions aboge that take bin indices directly do no range checking,
 * i.e., one must ensure that the bin number is between 0 and the number of
 * bins in the histogram. In contrast, the functions that use real values to
 * locate the correct bin use range checking, and silently discard values
 * outside the histogram (unless gmx_histogram_set_all() has been called).
 * gmx_histogram_find_bin() can be used to find the bin number corresponding
 * to a value, adhering to the gmx_histogram_set_all() setting.
 *
 * After a frame is sampled, gmx_histogram_finish_frame() needs to be
 * called.
 * After all data is sampled, gmx_histogram_finish() should be called.
 * After these calls, the histogram will be normalized such that the
 * sum over all bins equals the average number of samples per frame, and
 * the error field contains the error of the mean at each bin.
 *
 *
 * \section histogram_processing Post-processing
 *
 * The following functions can be used to process the histograms:
 *  - gmx_histogram_clone() makes a deep copy.
 *  - gmx_histogram_resample_dblbw() also makes a deep copy,
 *    but uses a binwidth that is double of the input one
 *    (see below for usage).
 *  - gmx_histogram_normalize_prob() normalizes a histogram such that its
 *    integral becomes one.
 *  - gmx_histogram_scale() and gmx_histogram_scale_vec() can be used
 *    to scale a histogram uniformly or non-uniformly, respectively.
 *
 * There are also functions to write histograms to a file:
 *  - gmx_histogram_write() writes a single histogram, optionally with
 *    errors.
 *  - gmx_histogram_write_array() writes multiple histograms with identical
 *    bins. Values and/or errors can be written.
 *  - gmx_histogram_write_cum_array() writes cumulative histograms for
 *    multiple histograms with identical bins.
 *
 * gmx_histogram_resample_dblbw() is useful if a histogram and its
 * cumulative sum are required at same points.
 * One can then sample the histogram with half the desired binwidth, and
 * then use gmx_histogram_resample_dblbw() to construct two different
 * versions of the histogram (by changing \p bIntegerBins), one suitable for
 * writing out the histogram and the other for writing the cumulative
 * distribution.
 *
 * Access to the histogram data is also provided through the following
 * functions:
 *  - gmx_histogram_get_nbins()
 *  - gmx_histogram_get_binwidth()
 *  - gmx_histogram_get_value()
 *  - gmx_histogram_get_bin_value()
 *  - gmx_histogram_get_values()
 *  - gmx_histogram_get_errors()
 *
 * gmx_histogram_free() can be used to free memory associated with a
 * histogram when it is no longer needed.
 */
/*! \internal \file
 * \brief Implementation of functions in histogram.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include <smalloc.h>
#include <vec.h>

#include <histogram.h>

/*! \internal \brief
 * Stores data for a histogram.
 */
struct gmx_histogram_t
{
    /** The left edge of the first bin. */
    real                 start;
    /** The right edge of the last bin. */
    real                 end;
    /** Bin width. */
    real                 binwidth;
    /** Number of bins. */
    int                  nbins;
    /** The final histogram. */
    double              *hist;
    /** Standard deviation of each bin (calculated from block histograms). */
    double              *histerr;

    /** Histogram type. */
    e_histogram_t        type;
    /** Histogram flags; */
    int                  flags;
    /** Block size for averaging. */
    int                  bsize;
    /** Output file for block histograms (can be NULL). */
    FILE                *blockfp;

    /** Inverse bin width. */
    real                 invbw;
    /** Current histogram value. */
    double              *chist;
    /** Current histogram count. */
    int                 *cn;
    /** Number of frames read for the current block so far. */
    int                  nframes;
    /** Number of blocks averaged so far. */
    int                  nblocks;
};

/*!
 * \param[out] hp      Histogram to create.
 * \param[in]  type    Histogram type to create;
 * \param[in]  nbins   Number of bins.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Initialized the histogram structure \p h with the provided values,
 * allocating the necessary memory.
 */
int
gmx_histogram_create(gmx_histogram_t **hp, e_histogram_t type, int nbins)
{
    gmx_histogram_t *h;

    if (nbins <= 0)
    {
        *hp = NULL;
        gmx_incons("number of histogram bins <= 0");
        return EINVAL;
    }

    snew(h, 1);
    h->start    = 0;
    h->end      = 0;
    h->binwidth = 0;
    h->nbins    = nbins;
    h->hist     = NULL;
    h->histerr  = NULL;

    h->type     = type;
    h->flags    = 0;

    h->bsize    = 0;
    h->blockfp  = NULL;
    h->invbw    = 0;
    h->chist    = NULL;
    h->cn       = NULL;
    h->nframes  = 0;
    h->nblocks  = 0;

    snew(h->hist,    nbins + 1);
    snew(h->histerr, nbins + 1);
    snew(h->cn,      nbins + 1);
    if (type != HIST_SIMPLE)
    {
        snew(h->chist, nbins + 1);
    }
    gmx_histogram_clear(h);

    *hp = h;
    return 0;
}

/*!
 * \param[out] hp      Histogram to create.
 * \param[in]  type    Histogram type to create;
 * \param[in]  start   Left edge/center of the first bin
 *   (depending on \p bIntegerBins).
 * \param[in]  end     Right edge/center of the last bin
 *   (depending on \p bIntegerBins).
 * \param[in]  binw    Bin width.
 * \param[in]  bIntegerBins If true, histogram bins are centered at
 *   \c start+n*binwidth, otherwise the centers are at
 *   \c start+(n+0.5)*binwidth.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Initialized the histogram structure \p h with the provided values,
 * allocating the necessary memory.
 */
int
gmx_histogram_create_range(gmx_histogram_t **hp, e_histogram_t type,
                           real start, real end, real binw, gmx_bool bIntegerBins)
{
    gmx_histogram_t *h;
    int              nbins;
    int              rc;

    *hp = NULL;
    if (start >= end || binw <= 0)
    {
        gmx_incons("histogram left edge larger than right edge or bin width <= 0");
        return EINVAL;
    }

    /* Adjust the end points and calculate the number of bins */
    if (bIntegerBins)
    {
        nbins = ceil((end - start) / binw) + 1;
        end   = start + (nbins - 1) * binw;
    }
    else
    {
        start = binw * floor(start / binw);
        end   = binw * ceil(end / binw);
        if (start != 0)
        {
            start -= binw;
        }
        end  += binw;
        nbins = (int)((end - start) / binw + 0.5);
    }
    /* Create the histogram */
    rc = gmx_histogram_create(&h, type, nbins);
    if (rc != 0)
    {
        return rc;
    }
    /* Set it up */
    gmx_histogram_set_integerbins(h, bIntegerBins);
    gmx_histogram_set_binwidth(h, start, binw);

    *hp = h;
    return 0;
}

/*!
 * \param[in,out] h       Histogram to clear.
 *
 * Histograms are automatically cleared when initialized; you only need to
 * call this function if you want to reuse a histogram structure that has
 * already been used for sampling.
 */
void
gmx_histogram_clear(gmx_histogram_t *h)
{
    int              i;

    if (h->nbins <= 0)
    {
        return;
    }
    for (i = 0; i <= h->nbins; ++i)
    {
        h->hist[i]    = 0;
        h->histerr[i] = 0;
        if (h->chist)
        {
            h->chist[i] = 0;
        }
        h->cn[i]      = 0;
    }
    h->nframes      = 0;
    h->nblocks      = 0;
}

/*!
 * \param[in] h  Histogram to free.
 *
 * The pointer \p h is invalid after the call.
 */
void
gmx_histogram_free(gmx_histogram_t *h)
{
    sfree(h->chist);
    sfree(h->cn);
    sfree(h->hist);
    sfree(h->histerr);
    sfree(h);
}

/*!
 * \param[in,out] h      Histogram data structure.
 * \param[in]     start  Left edge/center of the first bin
 *   (depending on gmx_histogram_set_integerbins()).
 * \param[in]     binw   Bin width.
 * \returns       0 on success, a non-zero error code on error.
 */
int
gmx_histogram_set_binwidth(gmx_histogram_t *h, real start, real binw)
{
    if (binw <= 0)
    {
        gmx_incons("histogram binwidth <= 0");
        return EINVAL;
    }
    if (h->flags & HIST_INTEGERBINS)
    {
        start   -= 0.5*binw;
    }
    h->start     = start;
    h->binwidth  = binw;
    h->end       = start + h->nbins * binw;
    h->invbw     = 1.0/binw;
    h->flags    |= HIST_INITBW;
    return 0;
}

/*!
 * \param[in,out] h      Histogram data structure.
 * \param[in]     start  Left edge/center of the first bin
 *   (depending on gmx_histogram_set_integerbins()).
 * \param[in]     end    Right edge/center of the last bin
 *   (depending on gmx_histogram_set_integerbins()).
 * \returns       0 on success, a non-zero error code on error.
 */
int
gmx_histogram_set_range(gmx_histogram_t *h, real start, real end)
{
    if (start >= end)
    {
        gmx_incons("histogram left edge larger than right edge");
        return EINVAL;
    }
    h->start         = start;
    h->end           = end;
    if (h->flags & HIST_INTEGERBINS)
    {
        h->binwidth  = (end - start) / (h->nbins - 1);
        start       -= 0.5*h->binwidth;
        end         += 0.5*h->binwidth;
    }
    else
    {
        h->binwidth  = (end - start) / h->nbins;
    }
    h->invbw         = 1.0/h->binwidth;
    h->flags        &= ~HIST_INITBW;
    return 0;
}

/*!
 * \param[in,out] h      Histogram data structure.
 * \param[in]     bIntegerBins If true, histogram bins are centered at
 *   \c start+n*binwidth, otherwise the centers are at
 *   \c start+(n+0.5)*binwidth.
 */
void
gmx_histogram_set_integerbins(gmx_histogram_t *h, gmx_bool bIntegerBins)
{
    /* Adjust the ranges if they have been initialized */
    if (h->start < h->end)
    {
        if (h->flags & HIST_INTEGERBINS)
        {
            h->start   += 0.5*h->binwidth;
            if (h->flags & HIST_INITBW)
            {
                h->end += 0.5*h->binwidth;
            }
            else
            {
                h->end -= 0.5*h->binwidth;
            }
        }
        if (bIntegerBins)
        {
            h->start   -= 0.5*h->binwidth;
            if (h->flags & HIST_INITBW)
            {
                h->end -= 0.5*h->binwidth;
            }
            else
            {
                h->end += 0.5*h->binwidth;
            }
        }
    }
    if (bIntegerBins)
    {
        h->flags |= HIST_INTEGERBINS;
    }
    else
    {
        h->flags &= ~HIST_INTEGERBINS;
    }
}

/*!
 * \param[in,out] h     Histogram data structure.
 * \param[in]     bAll  If true, values outside the histogram edges are added
 *   to the bins at the edges.
 *
 * \p bAll can be used to avoid rounding errors in cases where the histogram
 * spans the full range of possible values. If set, the values that are at
 * the exact maximum are still correctly included.
 */
void
gmx_histogram_set_all(gmx_histogram_t *h, gmx_bool bAll)
{
    if (bAll)
    {
        h->flags |= HIST_ALL;
    }
    else
    {
        h->flags &= ~HIST_ALL;
    }
}

/*!
 * \param[in,out] h      Histogram data structure.
 * \param[in]     bsize  Block size for error estimates.
 * \returns       0 on success, a non-zero error code on error.
 *
 * If \p bsize is zero, no block averaging or error estimates are done.
 * This is also the case if this function is not called.
 */
int
gmx_histogram_set_blocksize(gmx_histogram_t *h, int bsize)
{
    if (bsize < 0)
    {
        gmx_incons("histogram block size < 0");
        return EINVAL;
    }
    h->bsize = bsize;
    return 0;
}

/*!
 * \param[in,out] h  Histogram data structure.
 * \param[in]     fp File for block output.
 * \returns       0 on success, a non-zero error code on error.
 *
 * Sets a file into which each block histogram (the histogram after each
 * \c gmx_histogram_t::bsize frames) is written.
 * All histograms are written to the same file, separated by empty lines.
 * If this function is not called, the block histograms are only used for
 * error estimation, not written to a file.
 */
int
gmx_histogram_set_block_output(gmx_histogram_t *h, FILE *fp)
{
    if (h->bsize <= 0)
    {
        gmx_incons("histogram block size not set but output initialized");
        return EINVAL;
    }
    h->blockfp = fp;
    return 0;
}

/*!
 * \param[in] h     Histogram data.
 * \param[in] pos   Position.
 * \returns   Bin index in \p h corresponding to \p pos, or -1 if there is no
 *    such bin.
 *
 * If gmx_histogram_set_all() has been called with TRUE, values outside the
 * histogram are mapped to the bins at the edges.
 */
int
gmx_histogram_find_bin(gmx_histogram_t *h, real pos)
{
    if (pos < h->start)
    {
        if (h->flags & HIST_ALL)
        {
            return 0;
        }
        else
        {
            return -1;
        }
    }
    else if (pos >= h->end)
    {
        if (h->flags & HIST_ALL)
        {
            return h->nbins - 1;
        }
        else
        {
            return -1;
        }
    }
    /* Calculate the bin index if we are inside the box */
    return (int)((pos - h->start)*h->invbw);
}

/*!
 * \param[in] h     Histogram data.
 * \returns   Number of bins in \p h.
 */
int
gmx_histogram_get_nbins(gmx_histogram_t *h)
{
    return h->nbins;
}

/*!
 * \param[in] h     Histogram data.
 * \returns   Bin width of \p h, 0 if no binwidth has been set.
 */
real
gmx_histogram_get_binwidth(gmx_histogram_t *h)
{
    return h->binwidth;
}

/*!
 * \param[in]  h     Histogram data.
 * \param[in]  pos   Position.
 * \param[out] value Pointer to receive the value (can be NULL).
 * \param[out] err   Pointer to receive the value (can be NULL).
 *
 * If \p pos is outside the range of the histogram, zeros are returned.
 */
void
gmx_histogram_get_value(gmx_histogram_t *h, real pos, double *value, double *err)
{
    int     bin;
    double  v, e;

    if (pos < h->start || pos > h->end)
    {
        v = e = 0;
    }
    else
    {
        bin = gmx_histogram_find_bin(h, pos);
        if (bin < 0)
        {
            v = e = 0;
        }
        else
        {
            v = h->hist[bin];
            e = h->histerr[bin];
        }
    }
    if (value)
    {
        *value = v;
    }
    if (err)
    {
        *err = e;
    }
}

/*!
 * \param[in]  h     Histogram data.
 * \param[in]  bin   Bin number.
 * \param[out] value Pointer to receive the value (can be NULL).
 * \param[out] err   Pointer to receive the value (can be NULL).
 *
 * If \p bin is outside the valid range, zeros are returned.
 */
void
gmx_histogram_get_bin_value(gmx_histogram_t *h, int bin, double *value, double *err)
{
    double  v, e;

    if (bin < 0 || bin >= h->nbins)
    {
        v = e = 0;
    }
    else
    {
        v = h->hist[bin];
        e = h->histerr[bin];
    }
    if (value)
    {
        *value = v;
    }
    if (err)
    {
        *err = e;
    }
}

/*!
 * \param[in]  h     Histogram data.
 * \returns    Pointer to an array of values for \p h.
 *
 * The returned array has one element for each bin of \p h.
 * The returned pointer should not be freed.
 */
double *
gmx_histogram_get_values(gmx_histogram_t *h)
{
    return h->hist;
}

/*!
 * \param[in]  h     Histogram data.
 * \returns    Pointer to an array of errors for \p h.
 *
 * The returned array has one element for each bin of \p h.
 * The returned pointer should not be freed.
 */
double *
gmx_histogram_get_errors(gmx_histogram_t *h)
{
    return h->histerr;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     pos   Position.
 */
void
gmx_histogram_increment(gmx_histogram_t *h, real pos)
{
    int bin = gmx_histogram_find_bin(h, pos);
    if (bin < 0)
    {
        return;
    }
    h->cn[bin]++;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     pos   Position.
 * \param[in]     value Value to add.
 */
void
gmx_histogram_add(gmx_histogram_t *h, real pos, double value)
{
    int bin = gmx_histogram_find_bin(h, pos);
    if (bin < 0)
    {
        return;
    }
    h->chist[bin] += value;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     pos   Position.
 * \param[in]     value Value to add.
 */
void
gmx_histogram_add_item(gmx_histogram_t *h, real pos, double value)
{
    int bin = gmx_histogram_find_bin(h, pos);
    if (bin < 0)
    {
        return;
    }
    h->chist[bin] += value;
    h->cn[bin]++;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     bin   Bin number.
 *
 * No checks for out-of-bound errors are done: \p bin should be >=0 and
 * < \p h->nbins.
 */
void
gmx_histogram_increment_bin(gmx_histogram_t *h, int bin)
{
    h->cn[bin]++;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     bin   Bin number.
 * \param[in]     value Value to add.
 *
 * No checks for out-of-bound errors are done: \p bin should be >=0 and
 * < \p h->nbins.
 */
void
gmx_histogram_add_to_bin(gmx_histogram_t *h, int bin, double value)
{
    h->chist[bin] += value;
}

/*!
 * \param[in,out] h     Histogram data.
 * \param[in]     bin   Bin number.
 * \param[in]     value Value to add.
 *
 * No checks for out-of-bound errors are done: \p bin should be >=0 and
 * < \p h->nbins.
 */
void
gmx_histogram_add_item_to_bin(gmx_histogram_t *h, int bin, double value)
{
    h->chist[bin] += value;
    h->cn[bin]++;
}

/*! \brief
 * Processes a sampled block histogram.
 *
 * \param[in,out] h    Histogram data structure.
 */
static void
finish_histogram_block(gmx_histogram_t *h)
{
    int   i;
    real  v;

    if (h->nframes == 0)
    {
        return;
    }

    if (h->flags & HIST_ALL)
    {
        if (h->chist)
        {
            h->chist[h->nbins-1] += h->chist[h->nbins];
        }
        h->cn[h->nbins-1] += h->cn[h->nbins];
    }

    for (i = 0; i <= h->nbins; ++i)
    {
        if (h->chist)
        {
            v = h->chist[i] / (h->cn[i] > 0 ? h->cn[i] : h->nframes);
        }
        else
        {
            v = ((real)h->cn[i]) / h->nframes;
        }
        if (h->blockfp)
        {
            fprintf(h->blockfp, "%10g %10g\n", h->start + (i+0.5)*h->binwidth, v);
        }
        h->hist[i]    += v;
        h->histerr[i] += sqr(v);
        if (h->chist)
        {
            h->chist[i] = 0;
        }
        h->cn[i] = 0;
    }
    if (h->blockfp)
    {
        fprintf(h->blockfp, "\n");
    }
    h->nblocks++;
    h->nframes = 0;
}

/*!
 * \param[in,out] h    Histogram data structure.
 *
 * Should be called after each frame of data.
 *
 * \see gmx_histogram_finish()
 */
void
gmx_histogram_finish_frame(gmx_histogram_t *h)
{
    h->nframes++;
    if (h->nframes == h->bsize)
    {
        finish_histogram_block(h);
    }
}

/*!
 * \param[in,out] h    Histogram data structure.
 *
 * Normalizes the histogram.
 * Should be called after all frames have been sampled and before any
 * post-processing of the histogram.
 * If block size has not been set with gmx_histogram_set_blocksize(),
 * this function does no normalization, but it should still be called,
 * otherwise the \c gmx_histogram_t::hist array will contain only zeros.
 *
 * \see gmx_histogram_finish_frame()
 */
void
gmx_histogram_finish(gmx_histogram_t *h)
{
    int  i;

    if (h->nframes > 0 || h->bsize == 0)
    {
        if (h->nframes < h->bsize)
        {
            fprintf(stderr, "Last block smaller (%d frames) than the target size (%d frames) skipped \n",
                    h->nframes, h->bsize);
        }
        else
        {
            finish_histogram_block(h);
        }
    }
    if (h->nblocks == 0)
    {
        return;
    }

    for (i = 0; i <= h->nbins; ++i)
    {
        h->hist[i]    /= h->nblocks;
        h->histerr[i] /= h->nblocks;
        h->histerr[i]  = sqrt((h->histerr[i] - sqr(h->hist[i])) / h->nblocks);
    }
}

/*!
 * \param[out] destp        Destination histogram.
 * \param[in]  src          Source histogram.
 * \param[in]  bIntegerBins Control bin center position.
 *
 * Uses \p src to create a new histogram that has double the binwidth.
 * If \p bIntegerBins is TRUE, the first new bin will be centered at
 * \c src->start, otherwise the left edge will be at \c src->start.
 *
 * This function is mostly useful if one needs to sample both the
 * histogram and the cumulative histogram at same points.
 * To achieve this, first sample a histogram with half the desired
 * binwidth, and then use this function to create two different verions
 * of it with different values of \p bIntegerBins.
 * gmx_histogram_write_array() and gmx_histogram_write_cum_array() can
 * now used to write out the values correctly at identical values.
 *
 * Should be called only after gmx_histogram_finish() has been called.
 * \p src should have been sampled without gmx_histogram_set_integerbins().
 */
void
gmx_histogram_resample_dblbw(gmx_histogram_t **destp, gmx_histogram_t *src,
                             gmx_bool bIntegerBins)
{
    gmx_histogram_t *dest;
    int              i, j;
    real             v, ve;

    gmx_histogram_create(destp, src->type, src->nbins / 2);
    dest = *destp;
    sfree(dest->chist); dest->chist = NULL;
    sfree(dest->cn);    dest->cn    = NULL;
    gmx_histogram_set_binwidth(dest, src->start, src->binwidth * 2);
    gmx_histogram_set_integerbins(dest, bIntegerBins);

    for (i = j = 0; i < dest->nbins; ++i)
    {
        if (bIntegerBins && i == 0)
        {
            v  = src->hist[0];
            ve = sqr(src->histerr[0]);
            ++j;
        }
        else
        {
            v  = src->hist[j]         + src->hist[j+1];
            ve = sqr(src->histerr[j]) + sqr(src->histerr[j+1]);
            j += 2;
        }
        ve               = sqrt(ve);
        dest->hist[i]    = v;
        dest->histerr[i] = ve;
    }
    dest->hist[dest->nbins]    = 0;
    dest->histerr[dest->nbins] = 0;
}

/*!
 * \param[out] destp Destination histogram.
 * \param[in]  src   Source histogram.
 *
 * Makes a clone of \p src for post-processing.
 */
void
gmx_histogram_clone(gmx_histogram_t **destp, gmx_histogram_t *src)
{
    gmx_histogram_t *dest;

    snew(dest, 1);
    memcpy(dest, src, sizeof(*dest));

    /* These are not needed in post-processing */
    dest->blockfp = NULL;
    dest->chist   = NULL;
    dest->cn      = NULL;

    /* Make a deep copy of the actual histograms */
    snew(dest->hist,    src->nbins+1);
    snew(dest->histerr, src->nbins+1);
    memcpy(dest->hist,    src->hist,    (src->nbins+1)*sizeof(real));
    memcpy(dest->histerr, src->histerr, (src->nbins+1)*sizeof(real));

    *destp = dest;
}

/*!
 * \param[in,out] h  Histogram to normalize.
 *
 * Normalizes the histogram such that its integral equals one.
 */
void
gmx_histogram_normalize_prob(gmx_histogram_t *h)
{
    int  i;
    real sum;
    real normfac;

    sum = 0;
    for (i = 0; i <= h->nbins; ++i)
    {
        sum += h->hist[i];
    }

    normfac = h->invbw / sum;
    gmx_histogram_scale(h, normfac);
}

/*!
 * \param[in,out] h    Histogram to normalize.
 * \param[in]     norm Scaling factor.
 *
 * All bin values are multiplied by \p norm.
 */
void
gmx_histogram_scale(gmx_histogram_t *h, real norm)
{
    int  i;

    for (i = 0; i <= h->nbins; ++i)
    {
        h->hist[i]    *= norm;
        h->histerr[i] *= norm;
    }
}

/*!
 * \param[in,out] h    Histogram to normalize.
 * \param[in]     norm Scaling vector.
 *
 * The i'th bin is multiplied by \p norm[i].
 */
void
gmx_histogram_scale_vec(gmx_histogram_t *h, real norm[])
{
    int  i;

    for (i = 0; i < h->nbins; ++i)
    {
        h->hist[i]    *= norm[i];
        h->histerr[i] *= norm[i];
    }
    h->hist[h->nbins]    *= norm[h->nbins-1];
    h->histerr[h->nbins] *= norm[h->nbins-1];
}

/*! \brief
 * Makes some checks on output histograms and finds the maximum number of bins.
 *
 * \param[in]  n     Number of histograms in \p h.
 * \param[in]  h     Array of histograms.
 * \param[out] nbins Pointer to a value that will receive the maximum number
 *   of bins in \p h.
 */
static void
prepare_output(int n, gmx_histogram_t *h[], int *nbins)
{
    int  j;

    *nbins = 0;
    for (j = 0; j < n; ++j)
    {
        if (!gmx_within_tol(h[j]->start, h[0]->start, GMX_REAL_EPS))
        {
            fprintf(stderr, "gmx_ana_histogram_write: histogram start values not identical\n");
        }
        if (!gmx_within_tol(h[j]->binwidth, h[0]->binwidth, GMX_REAL_EPS))
        {
            fprintf(stderr, "gmx_ana_histogram_write: bin widths not identical\n");
        }
        if (*nbins < h[j]->nbins)
        {
            *nbins = h[j]->nbins;
        }
    }
}

/*!
 * \param[in] fp           Output file.
 * \param[in] h            Array of histograms to write.
 * \param[in] bErrors      If TRUE, histogram errors are written.
 *
 * Convenience wrapper for gmx_histogram_write_array() for writing a single
 * histogram.
 *
 * \see gmx_histogram_write_array()
 */
void
gmx_histogram_write(FILE *fp, gmx_histogram_t *h, gmx_bool bErrors)
{
    gmx_histogram_write_array(fp, 1, &h, TRUE, bErrors);
}

/*!
 * \param[in] fp           Output file.
 * \param[in] n            Number of histograms in \p h.
 * \param[in] h            Array of histograms to write.
 * \param[in] bValue       If TRUE, histogram values are written.
 * \param[in] bErrors      If TRUE, histogram errors are written.
 *
 * All the histograms in the array \p h should have the same bin widths and
 * left edges, otherwise the behavior is undefined.
 * The output format is one bin per line. The first number on the line is the
 * center of the bin, followed by one or two numbers for each histogram
 * (depending on \p bValue and \p bErrors).
 * If both \p bValue and \p bErrors are both TRUE, the values are written
 * before the errors.
 */
void
gmx_histogram_write_array(FILE *fp, int n, gmx_histogram_t *h[],
                          gmx_bool bValue, gmx_bool bErrors)
{
    int           i, j, nbins;

    prepare_output(n, h, &nbins);
    for (i = 0; i < nbins; ++i)
    {
        fprintf(fp, "%10g", h[0]->start + (i+0.5)*h[0]->binwidth);
        for (j = 0; j < n; ++j)
        {
            if (bValue)
            {
                fprintf(fp, " %10g", i < h[j]->nbins ? h[j]->hist[i] : 0.0);
            }
            if (bErrors)
            {
                fprintf(fp, " %10g", i < h[j]->nbins ? h[j]->histerr[i] : 0.0);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

/*!
 * \param[in] fp           Output file.
 * \param[in] n            Number of histograms in \p h.
 * \param[in] h            Array of histograms to write.
 * \param[in] bValue       If TRUE, histogram values are written.
 * \param[in] bErrors      If TRUE, histogram errors are written.
 *
 * Works as gmx_histogram_write_array(), but writes the cumulative
 * histograms.
 * The first column in output will be the right edges of the bins.
 *
 * \note
 * Error output is not currently supported (zeros will be written if asked).
 *
 * \see gmx_histogram_write_array()
 */
void
gmx_histogram_write_cum_array(FILE *fp, int n, gmx_histogram_t *h[],
                              gmx_bool bValue, gmx_bool bErrors)
{
    int           i, j, nbins;
    double       *sum;

    prepare_output(n, h, &nbins);
    snew(sum, n);

    fprintf(fp, "%10g", h[0]->start);
    for (j = 0; j < n; ++j)
    {
        if (bValue)
        {
            fprintf(fp, " %10g", 0.0);
        }
        if (bErrors)
        {
            fprintf(fp, " %10g", 0.0);
        }
    }
    fprintf(fp, "\n");
    for (i = 0; i < nbins; ++i)
    {
        fprintf(fp, "%10g", h[0]->start + (i+1)*h[0]->binwidth);
        for (j = 0; j < n; ++j)
        {
            sum[j] += i < h[j]->nbins ? h[j]->hist[i] : 0.0;
            if (bValue)
            {
                fprintf(fp, " %10g", sum[j]);
            }
            /* TODO: Error output not implemented */
            if (bErrors)
            {
                fprintf(fp, " %10g", 0.0);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    sfree(sum);
}
