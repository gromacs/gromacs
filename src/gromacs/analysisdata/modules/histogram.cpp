/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements classes in histogram.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/histogram.h"

#include <cmath>
#include <cstdint>

#include <limits>
#include <memory>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/analysisdata/framelocaldata.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "frameaverager.h"

namespace gmx
{
class AnalysisDataParallelOptions;
} // namespace gmx

namespace
{

//! Value used to signify that a real-valued histogram setting is not set.
const real UNDEFINED = std::numeric_limits<real>::max();
//! Checks whether \p value is defined.
bool isDefined(real value)
{
    return value != UNDEFINED;
}

} // namespace

namespace gmx
{

/********************************************************************
 * AnalysisHistogramSettingsInitializer
 */

AnalysisHistogramSettingsInitializer::AnalysisHistogramSettingsInitializer() :
    min_(UNDEFINED),
    max_(UNDEFINED),
    binWidth_(UNDEFINED),
    binCount_(0),
    bIntegerBins_(false),
    bRoundRange_(false),
    bIncludeAll_(false)
{
}


/********************************************************************
 * AnalysisHistogramSettings
 */

AnalysisHistogramSettings::AnalysisHistogramSettings() :
    firstEdge_(0.0), lastEdge_(0.0), binWidth_(0.0), inverseBinWidth_(0.0), binCount_(0), bAll_(false)
{
}


AnalysisHistogramSettings::AnalysisHistogramSettings(const AnalysisHistogramSettingsInitializer& settings)
{
    GMX_RELEASE_ASSERT(isDefined(settings.min_), "Histogram start value must be defined");
    GMX_RELEASE_ASSERT(!isDefined(settings.max_) || settings.max_ > settings.min_,
                       "Histogram end value must be larger than start value");
    GMX_RELEASE_ASSERT(!isDefined(settings.binWidth_) || settings.binWidth_ > 0.0,
                       "Histogram bin width must be positive");
    GMX_RELEASE_ASSERT(settings.binCount_ >= 0, "Histogram bin count must be positive");

    if (!isDefined(settings.max_))
    {
        GMX_RELEASE_ASSERT(isDefined(settings.binWidth_) && settings.binCount_ > 0,
                           "Not all required values provided");
        GMX_RELEASE_ASSERT(!settings.bRoundRange_, "Rounding only supported for min/max ranges");

        firstEdge_ = settings.min_;
        binCount_  = settings.binCount_;
        binWidth_  = settings.binWidth_;
        if (settings.bIntegerBins_)
        {
            firstEdge_ -= 0.5 * binWidth_;
        }
        lastEdge_ = firstEdge_ + binCount_ * binWidth_;
    }
    else
    {
        GMX_RELEASE_ASSERT(!(isDefined(settings.binWidth_) && settings.binCount_ > 0),
                           "Conflicting histogram bin specifications");
        GMX_RELEASE_ASSERT(isDefined(settings.binWidth_) || settings.binCount_ > 0,
                           "Not all required values provided");

        if (settings.bRoundRange_)
        {
            GMX_RELEASE_ASSERT(!settings.bIntegerBins_,
                               "Rounding and integer bins cannot be combined");
            GMX_RELEASE_ASSERT(isDefined(settings.binWidth_),
                               "Rounding only makes sense with defined binwidth");
            binWidth_  = settings.binWidth_;
            firstEdge_ = binWidth_ * std::floor(settings.min_ / binWidth_);
            lastEdge_  = binWidth_ * std::ceil(settings.max_ / binWidth_);
            binCount_  = gmx::roundToInt((lastEdge_ - firstEdge_) / binWidth_);
        }
        else
        {
            firstEdge_ = settings.min_;
            lastEdge_  = settings.max_;
            if (settings.binCount_ > 0)
            {
                binCount_ = settings.binCount_;
                if (settings.bIntegerBins_)
                {
                    GMX_RELEASE_ASSERT(settings.binCount_ > 1,
                                       "Bin count must be at least two with integer bins");
                    binWidth_ = (lastEdge_ - firstEdge_) / (binCount_ - 1);
                    firstEdge_ -= 0.5 * binWidth_;
                    lastEdge_ += 0.5 * binWidth_;
                }
                else
                {
                    binWidth_ = (lastEdge_ - firstEdge_) / binCount_;
                }
            }
            else
            {
                binWidth_ = settings.binWidth_;
                binCount_ = gmx::roundToInt((lastEdge_ - firstEdge_) / binWidth_);
                if (settings.bIntegerBins_)
                {
                    firstEdge_ -= 0.5 * binWidth_;
                    ++binCount_;
                }
                lastEdge_ = firstEdge_ + binCount_ * binWidth_;
            }
        }
    }

    inverseBinWidth_ = 1.0 / binWidth_;
    bAll_            = settings.bIncludeAll_;
}


int AnalysisHistogramSettings::findBin(real y) const
{
    if (y < firstEdge_)
    {
        return bAll_ ? 0 : -1;
    }
    int bin = static_cast<int>((y - firstEdge_) * inverseBinWidth_);
    if (bin >= binCount_)
    {
        return bAll_ ? binCount_ - 1 : -1;
    }
    return bin;
}


/********************************************************************
 * StaticAverageHistogram
 */

namespace
{

/*! \brief
 * Represents copies of average histograms.
 *
 * Methods in AbstractAverageHistogram that return new histogram instances
 * return objects of this class.
 * Initialization of values is handled in those methods.
 *
 * \ingroup module_analysisdata
 */
class StaticAverageHistogram : public AbstractAverageHistogram
{
public:
    StaticAverageHistogram();
    //! Creates an average histogram module with defined bin parameters.
    explicit StaticAverageHistogram(const AnalysisHistogramSettings& settings);

    // Copy and assign disallowed by base.
};

StaticAverageHistogram::StaticAverageHistogram() {}


StaticAverageHistogram::StaticAverageHistogram(const AnalysisHistogramSettings& settings) :
    AbstractAverageHistogram(settings)
{
}

} // namespace


/********************************************************************
 * AbstractAverageHistogram
 */

AbstractAverageHistogram::AbstractAverageHistogram() {}


AbstractAverageHistogram::AbstractAverageHistogram(const AnalysisHistogramSettings& settings) :
    settings_(settings)
{
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(), settings.binWidth());
}


AbstractAverageHistogram::~AbstractAverageHistogram() {}


void AbstractAverageHistogram::init(const AnalysisHistogramSettings& settings)
{
    settings_ = settings;
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(), settings.binWidth());
}


AverageHistogramPointer AbstractAverageHistogram::resampleDoubleBinWidth(bool bIntegerBins) const
{
    const int nbins = bIntegerBins ? (rowCount() + 1) / 2 : rowCount() / 2;

    AverageHistogramPointer dest(new StaticAverageHistogram(
            histogramFromBins(settings().firstEdge(), nbins, 2 * xstep()).integerBins(bIntegerBins)));
    dest->setColumnCount(columnCount());
    dest->allocateValues();

    for (int i = 0, j = 0; i < nbins; ++i)
    {
        const bool bFirstHalfBin = (bIntegerBins && i == 0);
        for (int c = 0; c < columnCount(); ++c)
        {
            const real v1 = bFirstHalfBin ? value(0, c).value() : value(j, c).value();
            const real v2 = bFirstHalfBin ? 0 : value(j + 1, c).value();
            const real e1 = bFirstHalfBin ? value(0, c).error() : value(j, c).error();
            const real e2 = bFirstHalfBin ? 0 : value(j + 1, c).error();
            dest->value(i, c).setValue(v1 + v2, std::sqrt(e1 * e1 + e2 * e2));
        }
        if (bFirstHalfBin)
        {
            ++j;
        }
        else
        {
            j += 2;
        }
    }
    return dest;
}


AverageHistogramPointer AbstractAverageHistogram::clone() const
{
    AverageHistogramPointer dest(new StaticAverageHistogram());
    copyContents(this, dest.get());
    dest->settings_ = settings_;
    return dest;
}


void AbstractAverageHistogram::normalizeProbability()
{
    for (int c = 0; c < columnCount(); ++c)
    {
        double sum = 0;
        for (int i = 0; i < rowCount(); ++i)
        {
            sum += value(i, c).value();
        }
        if (sum > 0.0)
        {
            scaleSingle(c, 1.0 / (sum * xstep()));
        }
    }
}

void AbstractAverageHistogram::makeCumulative()
{
    for (int c = 0; c < columnCount(); ++c)
    {
        double sum = 0;
        for (int i = 0; i < rowCount(); ++i)
        {
            sum += value(i, c).value();
            // Clear the error, as we don't cumulate that.
            value(i, c).clear();
            value(i, c).setValue(sum);
        }
    }
    setXAxis(settings().firstEdge() + settings().binWidth(), settings().binWidth());
}


void AbstractAverageHistogram::scaleSingle(int index, real factor)
{
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, index).value() *= factor;
        value(i, index).error() *= factor;
    }
}


void AbstractAverageHistogram::scaleAll(real factor)
{
    for (int i = 0; i < columnCount(); ++i)
    {
        scaleSingle(i, factor);
    }
}


void AbstractAverageHistogram::scaleAllByVector(const real factor[])
{
    for (int c = 0; c < columnCount(); ++c)
    {
        for (int i = 0; i < rowCount(); ++i)
        {
            value(i, c).value() *= factor[i];
            value(i, c).error() *= factor[i];
        }
    }
}


/********************************************************************
 * BasicAverageHistogramModule
 */

namespace internal
{

/*! \internal
 * \brief
 * Implements average histogram module that averages per-frame histograms.
 *
 * This class is used for accumulating average histograms in per-frame
 * histogram modules (those that use BasicHistogramImpl as their implementation
 * class).
 * There are two columns, first for the average and second for standard
 * deviation.
 *
 * \ingroup module_analysisdata
 */
class BasicAverageHistogramModule : public AbstractAverageHistogram, public AnalysisDataModuleSerial
{
public:
    BasicAverageHistogramModule();
    //! Creates an average histogram module with defined bin parameters.
    explicit BasicAverageHistogramModule(const AnalysisHistogramSettings& settings);

    using AbstractAverageHistogram::init;

    int flags() const override;

    void dataStarted(AbstractAnalysisData* data) override;
    void frameStarted(const AnalysisDataFrameHeader& header) override;
    void pointsAdded(const AnalysisDataPointSetRef& points) override;
    void frameFinished(const AnalysisDataFrameHeader& header) override;
    void dataFinished() override;

private:
    //! Averaging helper objects for each input data set.
    std::vector<AnalysisDataFrameAverager> averagers_;

    // Copy and assign disallowed by base.
};

BasicAverageHistogramModule::BasicAverageHistogramModule() {}


BasicAverageHistogramModule::BasicAverageHistogramModule(const AnalysisHistogramSettings& settings) :
    AbstractAverageHistogram(settings)
{
}


int BasicAverageHistogramModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipleDataSets;
}


void BasicAverageHistogramModule::dataStarted(AbstractAnalysisData* data)
{
    setColumnCount(data->dataSetCount());
    averagers_.resize(data->dataSetCount());
    for (int i = 0; i < data->dataSetCount(); ++i)
    {
        GMX_RELEASE_ASSERT(rowCount() == data->columnCount(i),
                           "Inconsistent data sizes, something is wrong in the initialization");
        averagers_[i].setColumnCount(data->columnCount(i));
    }
}


void BasicAverageHistogramModule::frameStarted(const AnalysisDataFrameHeader& /*header*/) {}


void BasicAverageHistogramModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    averagers_[points.dataSetIndex()].addPoints(points);
}


void BasicAverageHistogramModule::frameFinished(const AnalysisDataFrameHeader& /*header*/) {}


void BasicAverageHistogramModule::dataFinished()
{
    allocateValues();
    for (int i = 0; i < columnCount(); ++i)
    {
        averagers_[i].finish();
        for (int j = 0; j < rowCount(); ++j)
        {
            value(j, i).setValue(averagers_[i].average(j), std::sqrt(averagers_[i].variance(j)));
        }
    }
}


/********************************************************************
 * BasicHistogramImpl
 */

/*! \internal
 * \brief
 * Base class for private implementation classes for histogram modules.
 *
 * Actual implementation classes are derived from this and add an accumulation
 * data member that is specific to the histogram type in question.
 * This is done like this to keep implementation details out of the header, and
 * to not unnecessarily duplicate code.
 *
 * \ingroup module_analysisdata
 */
class BasicHistogramImpl
{
public:
    //! Smart pointer to manage an BasicAverageHistogramModule object.
    typedef std::shared_ptr<BasicAverageHistogramModule> BasicAverageHistogramModulePointer;

    BasicHistogramImpl();
    //! Creates an histogram impl with defined bin parameters.
    explicit BasicHistogramImpl(const AnalysisHistogramSettings& settings);
    // Virtual only for simplicity.
    virtual ~BasicHistogramImpl();

    /*! \brief
     * (Re)initializes the histogram from settings.
     */
    void init(const AnalysisHistogramSettings& settings);

    //! Storage implementation object.
    AnalysisDataStorage storage_;
    //! Settings for the histogram object.
    AnalysisHistogramSettings settings_;
    //! Averager module.
    BasicAverageHistogramModulePointer averager_;
};

BasicHistogramImpl::BasicHistogramImpl() : averager_(new BasicAverageHistogramModule()) {}


BasicHistogramImpl::BasicHistogramImpl(const AnalysisHistogramSettings& settings) :
    settings_(settings), averager_(new BasicAverageHistogramModule(settings))
{
}


BasicHistogramImpl::~BasicHistogramImpl() {}


void BasicHistogramImpl::init(const AnalysisHistogramSettings& settings)
{
    settings_ = settings;
    averager_->init(settings);
}

} // namespace internal


/********************************************************************
 * AnalysisDataSimpleHistogramModule
 */

/*! \internal \brief
 * Private implementation class for AnalysisDataSimpleHistogramModule.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataSimpleHistogramModule::Impl : public internal::BasicHistogramImpl
{
public:
    //! Shorthand for the per-frame accumulation data structure type.
    typedef AnalysisDataFrameLocalData<int64_t> FrameLocalData;

    Impl() {}
    //! Creates an histogram impl with defined bin parameters.
    explicit Impl(const AnalysisHistogramSettings& settings) : BasicHistogramImpl(settings) {}

    //! Accumulates the histogram within a frame.
    FrameLocalData accumulator_;
};

AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule() : impl_(new Impl()) {}


AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule(const AnalysisHistogramSettings& settings) :
    impl_(new Impl(settings))
{
}


AnalysisDataSimpleHistogramModule::~AnalysisDataSimpleHistogramModule() {}


void AnalysisDataSimpleHistogramModule::init(const AnalysisHistogramSettings& settings)
{
    impl_->init(settings);
}


AbstractAverageHistogram& AnalysisDataSimpleHistogramModule::averager()
{
    return *impl_->averager_;
}


const AnalysisHistogramSettings& AnalysisDataSimpleHistogramModule::settings() const
{
    return impl_->settings_;
}


int AnalysisDataSimpleHistogramModule::frameCount() const
{
    return impl_->storage_.frameCount();
}


int AnalysisDataSimpleHistogramModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint | efAllowMissing | efAllowMultipleDataSets;
}


bool AnalysisDataSimpleHistogramModule::parallelDataStarted(AbstractAnalysisData*              data,
                                                            const AnalysisDataParallelOptions& options)
{
    addModule(impl_->averager_);
    const int dataSetCount = data->dataSetCount();
    const int columnCount  = settings().binCount();
    setDataSetCount(dataSetCount);
    impl_->accumulator_.setDataSetCount(dataSetCount);
    for (int i = 0; i < dataSetCount; ++i)
    {
        setColumnCount(i, columnCount);
        impl_->accumulator_.setColumnCount(i, columnCount);
    }
    impl_->accumulator_.init(options);
    impl_->storage_.startParallelDataStorage(this, &moduleManager(), options);
    return true;
}


void AnalysisDataSimpleHistogramModule::frameStarted(const AnalysisDataFrameHeader& header)
{
    impl_->accumulator_.frameData(header.index()).clear();
}


void AnalysisDataSimpleHistogramModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    Impl::FrameLocalData::DataSetHandle handle =
            impl_->accumulator_.frameDataSet(points.frameIndex(), points.dataSetIndex());
    for (int i = 0; i < points.columnCount(); ++i)
    {
        if (points.present(i))
        {
            const int bin = settings().findBin(points.y(i));
            if (bin != -1)
            {
                handle.value(bin) += 1;
            }
        }
    }
}


void AnalysisDataSimpleHistogramModule::frameFinished(const AnalysisDataFrameHeader& header)
{
    Impl::FrameLocalData::FrameHandle handle      = impl_->accumulator_.frameData(header.index());
    AnalysisDataStorageFrame&         frame       = impl_->storage_.startFrame(header);
    const int                         columnCount = settings().binCount();
    for (int s = 0; s < dataSetCount(); ++s)
    {
        Impl::FrameLocalData::DataSetHandle dataSet = handle.dataSet(s);
        frame.selectDataSet(s);
        for (int i = 0; i < columnCount; ++i)
        {
            frame.setValue(i, dataSet.value(i));
        }
    }
    frame.finishFrame();
}


void AnalysisDataSimpleHistogramModule::frameFinishedSerial(int frameIndex)
{
    impl_->storage_.finishFrameSerial(frameIndex);
}


void AnalysisDataSimpleHistogramModule::dataFinished()
{
    impl_->storage_.finishDataStorage();
}


AnalysisDataFrameRef AnalysisDataSimpleHistogramModule::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool AnalysisDataSimpleHistogramModule::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataWeightedHistogramModule
 */

/*! \internal \brief
 * Private implementation class for AnalysisDataWeightedHistogramModule.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataWeightedHistogramModule::Impl : public internal::BasicHistogramImpl
{
public:
    //! Shorthand for the per-frame accumulation data structure type.
    typedef AnalysisDataFrameLocalData<double> FrameLocalData;

    Impl() {}
    //! Creates an histogram impl with defined bin parameters.
    explicit Impl(const AnalysisHistogramSettings& settings) : BasicHistogramImpl(settings) {}

    //! Accumulates the histogram within a frame.
    FrameLocalData accumulator_;
};

AnalysisDataWeightedHistogramModule::AnalysisDataWeightedHistogramModule() : impl_(new Impl()) {}


AnalysisDataWeightedHistogramModule::AnalysisDataWeightedHistogramModule(const AnalysisHistogramSettings& settings) :
    impl_(new Impl(settings))
{
}


AnalysisDataWeightedHistogramModule::~AnalysisDataWeightedHistogramModule() {}


void AnalysisDataWeightedHistogramModule::init(const AnalysisHistogramSettings& settings)
{
    impl_->init(settings);
}


AbstractAverageHistogram& AnalysisDataWeightedHistogramModule::averager()
{
    return *impl_->averager_;
}


const AnalysisHistogramSettings& AnalysisDataWeightedHistogramModule::settings() const
{
    return impl_->settings_;
}


int AnalysisDataWeightedHistogramModule::frameCount() const
{
    return impl_->storage_.frameCount();
}


int AnalysisDataWeightedHistogramModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;
}


bool AnalysisDataWeightedHistogramModule::parallelDataStarted(AbstractAnalysisData* data,
                                                              const AnalysisDataParallelOptions& options)
{
    addModule(impl_->averager_);
    const int dataSetCount = data->dataSetCount();
    const int columnCount  = settings().binCount();
    setDataSetCount(dataSetCount);
    impl_->accumulator_.setDataSetCount(dataSetCount);
    for (int i = 0; i < dataSetCount; ++i)
    {
        setColumnCount(i, columnCount);
        impl_->accumulator_.setColumnCount(i, columnCount);
    }
    impl_->accumulator_.init(options);
    impl_->storage_.startParallelDataStorage(this, &moduleManager(), options);
    return true;
}


void AnalysisDataWeightedHistogramModule::frameStarted(const AnalysisDataFrameHeader& header)
{
    impl_->accumulator_.frameData(header.index()).clear();
}


void AnalysisDataWeightedHistogramModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (points.firstColumn() != 0 || points.columnCount() < 2)
    {
        GMX_THROW(APIError("Invalid data layout"));
    }
    int bin = settings().findBin(points.y(0));
    if (bin != -1)
    {
        Impl::FrameLocalData::DataSetHandle handle =
                impl_->accumulator_.frameDataSet(points.frameIndex(), points.dataSetIndex());
        for (int i = 1; i < points.columnCount(); ++i)
        {
            handle.value(bin) += points.y(i);
        }
    }
}


void AnalysisDataWeightedHistogramModule::frameFinished(const AnalysisDataFrameHeader& header)
{
    Impl::FrameLocalData::FrameHandle handle      = impl_->accumulator_.frameData(header.index());
    AnalysisDataStorageFrame&         frame       = impl_->storage_.startFrame(header);
    const int                         columnCount = settings().binCount();
    for (int s = 0; s < dataSetCount(); ++s)
    {
        Impl::FrameLocalData::DataSetHandle dataSet = handle.dataSet(s);
        frame.selectDataSet(s);
        for (int i = 0; i < columnCount; ++i)
        {
            frame.setValue(i, dataSet.value(i));
        }
    }
    frame.finishFrame();
}


void AnalysisDataWeightedHistogramModule::frameFinishedSerial(int frameIndex)
{
    impl_->storage_.finishFrameSerial(frameIndex);
}


void AnalysisDataWeightedHistogramModule::dataFinished()
{
    impl_->storage_.finishDataStorage();
}


AnalysisDataFrameRef AnalysisDataWeightedHistogramModule::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool AnalysisDataWeightedHistogramModule::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataBinAverageModule
 */

class AnalysisDataBinAverageModule::Impl
{
public:
    Impl() {}
    explicit Impl(const AnalysisHistogramSettings& settings) : settings_(settings) {}

    //! Histogram settings.
    AnalysisHistogramSettings settings_;
    //! Averaging helper objects for each input data set.
    std::vector<AnalysisDataFrameAverager> averagers_;
};

AnalysisDataBinAverageModule::AnalysisDataBinAverageModule() : impl_(new Impl())
{
    setColumnCount(3);
}


AnalysisDataBinAverageModule::AnalysisDataBinAverageModule(const AnalysisHistogramSettings& settings) :
    impl_(new Impl(settings))
{
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(), settings.binWidth());
}


AnalysisDataBinAverageModule::~AnalysisDataBinAverageModule() {}


void AnalysisDataBinAverageModule::init(const AnalysisHistogramSettings& settings)
{
    impl_->settings_ = settings;
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(), settings.binWidth());
}


const AnalysisHistogramSettings& AnalysisDataBinAverageModule::settings() const
{
    return impl_->settings_;
}


int AnalysisDataBinAverageModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;
}


void AnalysisDataBinAverageModule::dataStarted(AbstractAnalysisData* data)
{
    setColumnCount(data->dataSetCount());
    impl_->averagers_.resize(data->dataSetCount());
    for (int i = 0; i < data->dataSetCount(); ++i)
    {
        impl_->averagers_[i].setColumnCount(rowCount());
    }
}


void AnalysisDataBinAverageModule::frameStarted(const AnalysisDataFrameHeader& /*header*/) {}


void AnalysisDataBinAverageModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (points.firstColumn() != 0 || points.columnCount() < 2)
    {
        GMX_THROW(APIError("Invalid data layout"));
    }
    int bin = settings().findBin(points.y(0));
    if (bin != -1)
    {
        AnalysisDataFrameAverager& averager = impl_->averagers_[points.dataSetIndex()];
        for (int i = 1; i < points.columnCount(); ++i)
        {
            averager.addValue(bin, points.y(i));
        }
    }
}


void AnalysisDataBinAverageModule::frameFinished(const AnalysisDataFrameHeader& /*header*/) {}


void AnalysisDataBinAverageModule::dataFinished()
{
    allocateValues();
    for (int i = 0; i < columnCount(); ++i)
    {
        AnalysisDataFrameAverager& averager = impl_->averagers_[i];
        averager.finish();
        for (int j = 0; j < rowCount(); ++j)
        {
            value(j, i).setValue(averager.average(j), std::sqrt(averager.variance(j)));
        }
    }
    valuesReady();
}

} // namespace gmx
