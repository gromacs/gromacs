/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \file
 * \brief
 * Declares analysis data modules for calculating histograms.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_HISTOGRAM_H
#define GMX_ANALYSISDATA_MODULES_HISTOGRAM_H

#include "../analysisdata.h"
#include "../arraydata.h"
#include "../datamodule.h"

namespace gmx
{

class HistogramAverageModule;

/*! \brief
 * Abstract base class for per-frame histogram modules.
 *
 * \ingroup module_analysisdata
 */
class AbstractHistogramModule : public AbstractAnalysisDataStored,
                                public AnalysisDataModuleInterface
{
    public:
        virtual ~AbstractHistogramModule();

        /*! \brief
         * Initializes the histogram using bin width and the number of bins.
         */
        int initNBins(real miny, real binw, int nbins,
                      bool bIntegerBins = false);
        /*! \brief
         * Initializes the histogram using a range and a bin width.
         */
        int initRange(real miny, real maxy, real binw,
                      bool bIntegerBins = false);
        /*! \brief
         * Sets the histogram to match all values.
         *
         * If \p bAll is true, the histogram behaves as if the bins at the ends
         * extended to +-infinity.
         */
        void setAll(bool bAll);

        /*! \brief
         * Returns the average histogram over all frames.
         *
         * Can be called already before the histogram is calculated to
         * customize the way the average histogram is calculated.
         *
         * \see HistogramAverageModule
         */
        HistogramAverageModule *averager();

        //! Returns the number of bins in the histogram.
        int nbins() const { return _nbins; }
        //! Returns the width of a bin in the histogram.
        real binwidth() const { return _binwidth; }
        //! Returns a zero-based bin index for a value, or -1 if not in range.
        int findBin(real y) const;

        virtual int flags() const;

        virtual int dataStarted(AbstractAnalysisData *data);
        virtual int frameStarted(real x, real dx);
        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present) = 0;
        virtual int frameFinished();
        virtual int dataFinished();

    protected:
        AbstractHistogramModule();

        //! Actual histogram data.
        real                   *_hist;

    private:
        void createAverager();

        HistogramAverageModule *_averager;
        int                     _nbins;
        real                    _miny;
        real                    _maxy;
        real                    _binwidth;
        real                    _invbw;
        bool                    _bAll;

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for averaging histograms over frames.
 *
 * The averaging module for a per-frame histogram is always created by the
 * AbstractHistogramModule class, and can be accessed using
 * AbstractHistogramModule::averager().
 * The user can alter some properties of the average histogram directly, but
 * the main use of the object is to postprocess the histogram once the
 * calculation is finished.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class HistogramAverageModule : public AbstractAnalysisArrayData,
                               public AnalysisDataModuleInterface
{
    public:
        virtual int flags() const;

        virtual int dataStarted(AbstractAnalysisData *data);
        virtual int frameStarted(real x, real dx);
        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present);
        virtual int frameFinished();
        virtual int dataFinished();

        /*! \brief
         * Sets the averager to ignore missing values.
         */
        void setIgnoreMissing(bool bIgnoreMissing);

        /*! \brief
         * Creates a copy of the histogram with double the bin width.
         */
        HistogramAverageModule *resampleDoubleBinWidth(bool bIntegerBins) const;
        //! Creates a deep copy of the histogram.
        HistogramAverageModule *clone() const;
        //! Normalizes the histogram such that the integral over it is one.
        void normalizeProbability();
        //! Scales the value of each bin by an uniform scaling factor.
        void scale(real norm);
        //! Scales the value of each bin by a different scaling factor.
        void scaleVector(real norm[]);
        /*! \brief
         * Notifies attached modules of the histogram data.
         *
         * After this function has been called, it is no longer possible to
         * alter the histogram.
         */
        int done() { return AbstractAnalysisArrayData::valuesReady(); }

    private:
        HistogramAverageModule();

        int                     _nframes;
        bool                    _bIgnoreMissing;

        friend class AbstractHistogramModule;

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for per-frame histograms.
 *
 * Output data contains the same number of frames as the input data.
 * Each frame contains the histogram for the points in that frame.
 * All input columns are averaged into the same histogram.
 * The number of columns equals the number of bins in the histogram.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataSimpleHistogramModule : public AbstractHistogramModule
{
    public:
        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present);

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for per-frame weighted histograms.
 *
 * Output data contains the same number of frames as the input data.
 * Each frame contains the histogram for the points in that frame, interpreted
 * such that the first column passed to pointsAdded() determines the bin and
 * the rest give weights to be added to that bin (input data should have at
 * least two colums, and at least two columns should be added at the same time).
 * All input columns are averaged into the same histogram.
 * The number of columns equals the number of bins in the histogram.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataWeightedHistogramModule : public AbstractHistogramModule
{
    public:
        virtual int flags() const;

        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present);

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for per-frame bin averages.
 *
 * Output data contains the same number of frames as the input data.
 * Each frame contains the average for the points in that frame within each bin.
 * The input data is interpreted such that the first column passed to
 * pointsAdded() determines the bin and the rest give values to be added to
 * that bin (input data should have at least two colums, and at least two
 * columns should be added at the same time).
 * All input columns are averaged into the same histogram.
 * The number of columns equals the number of bins.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataBinAverageModule : public AbstractHistogramModule
{
    public:
        AnalysisDataBinAverageModule();
        virtual ~AnalysisDataBinAverageModule();

        //! Ignore missing bins in the average histogram.
        void setIgnoreMissing(bool bIgnoreMissing);

        virtual int flags() const;

        virtual int dataStarted(AbstractAnalysisData *data);
        virtual int frameStarted(real x, real dx);
        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present);
        virtual int frameFinished();

    private:
        int                    *_n;
        bool                   *_present;
        bool                    _bIgnoreMissing;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
