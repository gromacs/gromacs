/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares analysis data modules for calculating histograms.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_MODULES_HISTOGRAM_H
#define GMX_ANALYSISDATA_MODULES_HISTOGRAM_H

#include <boost/shared_ptr.hpp>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/arraydata.h"
#include "gromacs/analysisdata/datamodule.h"

namespace gmx
{

class AnalysisHistogramSettings;

/*! \brief
 * Provides "named parameter" idiom for constructing histograms.
 *
 * \see histogramFromBins()
 * \see histogramFromRange()
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisHistogramSettingsInitializer
{
    public:
        /*! \brief
         * Creates an empty initializer.
         *
         * Should not be called directly, but histogramFromRange() or
         * histogramFromBins() should be used instead.
         */
        AnalysisHistogramSettingsInitializer();

        /*! \brief
         * Sets the first bin location.
         *
         * Typically should not be called directly, but through
         * histogramFromBins().
         */
        AnalysisHistogramSettingsInitializer &start(real min)
        { min_ = min; return *this; }
        /*! \brief
         * Sets the number of bins in the histogram.
         *
         * If only the first bin location is specified, this value is required
         * (and automatically provided if histogramFromBins() is used).
         * If both the first and last bins are specified, either this value or
         * binWidth() is required.
         */
        AnalysisHistogramSettingsInitializer &binCount(int binCount)
        { binCount_ = binCount; return *this; }
        /*! \brief
         * Sets the first and last bin locations.
         *
         * Typically should not be called directly, but through
         * histogramFromRange().
         */
        AnalysisHistogramSettingsInitializer &range(real min, real max)
        { min_ = min; max_ = max; return *this; }
        /*! \brief
         * Sets the bin width of the histogram.
         *
         * If only the first bin location is specified, this value is required
         * (and automatically provided if histogramFromBins() is used).
         * If both the first and last bins are specified, either this value or
         * binCount() is required.
         * If a bin width is provided with both first and last bin locations,
         * and the given bin width does not divide the range exactly, the last
         * bin location is adjusted to match.
         */
        AnalysisHistogramSettingsInitializer &binWidth(real binWidth)
        { binWidth_ = binWidth; return *this; }
        /*! \brief
         * Indicate that first and last bin locations to specify bin centers.
         *
         * If set, the first and last bin locations are interpreted as bin
         * centers.
         * If not set (the default), the first and last bin locations are
         * interpreted as the edges of the whole histogram.
         *
         * Cannot be specified together with roundRange().
         */
        AnalysisHistogramSettingsInitializer &integerBins(bool enabled = true)
        { bIntegerBins_ = enabled; return *this; }
        /*! \brief
         * Round first and last bin locations.
         *
         * If set, the resulting histogram will cover the range specified, but
         * the actual bin locations will be rounded such that the edges fall
         * on multiples of the bin width.
         * Only implemented when both first and last bin location and bin width
         * are defined.
         * Cannot be specified together with integerBins() or with binCount().
         */
        AnalysisHistogramSettingsInitializer &roundRange(bool enabled = true)
        { bRoundRange_ = enabled; return *this; }
        /*! \brief
         * Sets the histogram to match all values.
         *
         * If set, the histogram behaves as if the bins at the ends extended to
         * +-infinity.
         */
        AnalysisHistogramSettingsInitializer &includeAll(bool enabled = true)
        { bIncludeAll_ = enabled; return *this; }

    private:
        real                    min_;
        real                    max_;
        real                    binWidth_;
        int                     binCount_;
        bool                    bIntegerBins_;
        bool                    bRoundRange_;
        bool                    bIncludeAll_;

        friend class AnalysisHistogramSettings;
};

/*! \brief
 * Initializes a histogram using a range and a bin width.
 *
 * Does not throw.
 *
 * \inpublicapi
 */
inline AnalysisHistogramSettingsInitializer
histogramFromRange(real min, real max)
{
    return AnalysisHistogramSettingsInitializer().range(min, max);
}

/*! \brief
 * Initializes a histogram using bin width and the number of bins.
 *
 * Does not throw.
 *
 * \inpublicapi
 */
inline AnalysisHistogramSettingsInitializer
histogramFromBins(real start, int nbins, real binwidth)
{
    return AnalysisHistogramSettingsInitializer()
               .start(start).binCount(nbins).binWidth(binwidth);
}


/*! \brief
 * Contains parameters that specify histogram bin locations.
 *
 * Methods in this class do not throw.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisHistogramSettings
{
    public:
        //! Initializes undefined parameters.
        AnalysisHistogramSettings();
        /*! \brief
         * Initializes parameters based on a named parameter object.
         *
         * This constructor is not explicit to allow initialization of
         * histograms directly from AnalysisHistogramSettingsInitializer:
         * \code
           gmx::AnalysisDataSimpleHistogramModule *hist =
                   new gmx::AnalysisDataSimpleHistogramModule(
                           histogramFromRange(0.0, 5.0).binWidth(0.5));
         * \endcode
         */
        AnalysisHistogramSettings(const AnalysisHistogramSettingsInitializer &settings);

        //! Returns the left edge of the first bin.
        real firstEdge() const { return firstEdge_; }
        //! Returns the right edge of the first bin.
        real lastEdge() const { return lastEdge_; }
        //! Returns the number of bins in the histogram.
        int binCount() const { return binCount_; }
        //! Returns the width of a bin in the histogram.
        real binWidth() const { return binWidth_; }
        //! Whether values beyond the edges are mapped to the edge bins.
        bool includeAll() const { return bAll_; }
        //! Returns a zero-based bin index for a value, or -1 if not in range.
        int findBin(real y) const;

    private:
        real                    firstEdge_;
        real                    lastEdge_;
        real                    binWidth_;
        real                    inverseBinWidth_;
        int                     binCount_;
        bool                    bAll_;
};


class AbstractAverageHistogram;

//! Smart pointer to manage an AbstractAverageHistogram object.
typedef boost::shared_ptr<AbstractAverageHistogram>
    AverageHistogramPointer;

/*! \brief
 * Base class for representing histograms averaged over frames.
 *
 * The averaging module for a per-frame histogram is always created by the
 * histogram module class (e.g., AnalysisDataSimpleHistogramModule), and can be
 * accessed using, e.g., AnalysisDataSimpleHistogramModule::averager().
 * The user can alter some properties of the average histogram directly, but
 * the main use of the object is to postprocess the histogram once the
 * calculation is finished.
 *
 * This class can represent multiple histograms in one object: each column in
 * the data is an independent histogram.
 * The X values correspond to center of the bins, except for a cumulative
 * histogram made with makeCumulative().
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AbstractAverageHistogram : public AbstractAnalysisArrayData
{
    public:
        virtual ~AbstractAverageHistogram();

        //! Returns bin properties for the histogram.
        const AnalysisHistogramSettings &settings() const { return settings_; }

        /*! \brief
         * Creates a copy of the histogram with double the bin width.
         *
         * \param[in] bIntegerBins If `true`, the first bin in the result will
         *     cover the first bin from the source. Otherwise, the first bin
         *     will cover first two bins from the source.
         * \throws std::bad_alloc if out of memory.
         *
         * The caller is responsible of deleting the returned object.
         */
        AverageHistogramPointer resampleDoubleBinWidth(bool bIntegerBins) const;
        /*! \brief
         * Creates a deep copy of the histogram.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * The returned histogram is not necessarily of the same dynamic type
         * as the original object, but contains the same data from the point of
         * view of the AbstractAverageHistogram interface.
         *
         * The caller is responsible of deleting the returned object.
         */
        AverageHistogramPointer clone() const;
        //! Normalizes the histogram such that the integral over it is one.
        void normalizeProbability();
        /*! \brief
         * Makes the histograms cumulative by summing up each bin to all bins
         * after it.
         *
         * The X values in the data are adjusted such that they match the right
         * edges of bins instead of bin centers.
         */
        void makeCumulative();
        //! Scales a single histogram by a uniform scaling factor.
        void scaleSingle(int index, real factor);
        //! Scales all histograms by a uniform scaling factor.
        void scaleAll(real factor);
        //! Scales the value of each bin by a different scaling factor.
        void scaleAllByVector(real factor[]);
        /*! \brief
         * Notifies attached modules of the histogram data.
         *
         * After this function has been called, it is no longer possible to
         * alter the histogram.
         */
        void done() { AbstractAnalysisArrayData::valuesReady(); }

    protected:
        /*! \brief
         * Creates a histogram module with undefined bins.
         *
         * Bin parameters must be defined with init() before data input is
         * started.
         */
        AbstractAverageHistogram();
        //! Creates a histogram module with defined bin parameters.
        explicit AbstractAverageHistogram(const AnalysisHistogramSettings &settings);

        /*! \brief
         * (Re)initializes the histogram from settings.
         */
        void init(const AnalysisHistogramSettings &settings);

    private:
        AnalysisHistogramSettings  settings_;

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for per-frame histograms.
 *
 * Output data contains the same number of frames and data sets as the input
 * data.  Each frame contains the histogram(s) for the points in that frame.
 * Each input data set is processed independently into the corresponding output
 * data set.  Missing values are ignored.
 * All input columns for a data set are averaged into the same histogram.
 * The number of columns for all data sets equals the number of bins in the
 * histogram.
 *
 * The histograms are accumulated as 64-bit integers within a frame and summed
 * in double precision across frames, even if the output data is in single
 * precision.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataSimpleHistogramModule : public AbstractAnalysisData,
                                          public AnalysisDataModuleParallel
{
    public:
        /*! \brief
         * Creates a histogram module with undefined bins.
         *
         * Bin parameters must be defined with init() before data input is
         * started.
         */
        AnalysisDataSimpleHistogramModule();
        //! Creates a histogram module with defined bin parameters.
        explicit AnalysisDataSimpleHistogramModule(const AnalysisHistogramSettings &settings);
        virtual ~AnalysisDataSimpleHistogramModule();

        /*! \brief
         * (Re)initializes the histogram from settings.
         */
        void init(const AnalysisHistogramSettings &settings);

        /*! \brief
         * Returns the average histogram over all frames.
         *
         * Can be called already before the histogram is calculated to
         * customize the way the average histogram is calculated.
         *
         * \see AbstractAverageHistogram
         */
        AbstractAverageHistogram &averager();

        //! Returns bin properties for the histogram.
        const AnalysisHistogramSettings &settings() const;

        virtual int frameCount() const;

        virtual int flags() const;

        virtual bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void frameFinishedSerial(int frameIndex);
        virtual void dataFinished();

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        class Impl;

        PrivateImplPointer<Impl> impl_;

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for per-frame weighted histograms.
 *
 * Output data contains the same number of frames and data sets as the input
 * data.  Each frame contains the histogram(s) for the points in that frame,
 * interpreted such that the first column passed to pointsAdded() determines
 * the bin and the rest give weights to be added to that bin (input data should
 * have at least two colums, and at least two columns should be added at the
 * same time).
 * Each input data set is processed independently into the corresponding output
 * data set.
 * All input columns for a data set are averaged into the same histogram.
 * The number of columns for all data sets equals the number of bins in the
 * histogram.
 *
 * The histograms are accumulated in double precision, even if the output data
 * is in single precision.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataWeightedHistogramModule : public AbstractAnalysisData,
                                            public AnalysisDataModuleParallel
{
    public:
        //! \copydoc AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule()
        AnalysisDataWeightedHistogramModule();
        //! \copydoc AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule(const AnalysisHistogramSettings &)
        explicit AnalysisDataWeightedHistogramModule(const AnalysisHistogramSettings &settings);
        virtual ~AnalysisDataWeightedHistogramModule();

        //! \copydoc AnalysisDataSimpleHistogramModule::init()
        void init(const AnalysisHistogramSettings &settings);

        //! \copydoc AnalysisDataSimpleHistogramModule::averager()
        AbstractAverageHistogram &averager();

        //! \copydoc AnalysisDataSimpleHistogramModule::settings()
        const AnalysisHistogramSettings &settings() const;

        virtual int frameCount() const;

        virtual int flags() const;

        virtual bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void frameFinishedSerial(int frameIndex);
        virtual void dataFinished();

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        class Impl;

        PrivateImplPointer<Impl> impl_;

        // Copy and assign disallowed by base.
};


/*! \brief
 * Data module for bin averages.
 *
 * Output data contains one row for each bin; see AbstractAverageHistogram.
 * Output data contains one column for each input data set.
 * The value in a column is the average over all frames of that data set for
 * that bin.
 * The input data is interpreted such that the first column passed to
 * pointsAdded() determines the bin and the rest give values to be added to
 * that bin (input data should have at least two colums, and at least two
 * columns should be added at the same time).
 * All input columns for a data set are averaged into the same histogram.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataBinAverageModule : public AbstractAnalysisArrayData,
                                     public AnalysisDataModuleSerial
{
    public:
        //! \copydoc AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule()
        AnalysisDataBinAverageModule();
        //! \copydoc AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule(const AnalysisHistogramSettings &)
        explicit AnalysisDataBinAverageModule(const AnalysisHistogramSettings &settings);
        virtual ~AnalysisDataBinAverageModule();

        //! \copydoc AnalysisDataSimpleHistogramModule::init()
        void init(const AnalysisHistogramSettings &settings);

        //! \copydoc AnalysisDataSimpleHistogramModule::settings()
        const AnalysisHistogramSettings &settings() const;

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        class Impl;

        PrivateImplPointer<Impl>   impl_;

        // Copy and assign disallowed by base.
};

//! Smart pointer to manage an AnalysisDataSimpleHistogramModule object.
typedef boost::shared_ptr<AnalysisDataSimpleHistogramModule>
    AnalysisDataSimpleHistogramModulePointer;
//! Smart pointer to manage an AnalysisDataWeightedHistogramModule object.
typedef boost::shared_ptr<AnalysisDataWeightedHistogramModule>
    AnalysisDataWeightedHistogramModulePointer;
//! Smart pointer to manage an AnalysisDataBinAverageModule object.
typedef boost::shared_ptr<AnalysisDataBinAverageModule>
    AnalysisDataBinAverageModulePointer;

} // namespace gmx

#endif
