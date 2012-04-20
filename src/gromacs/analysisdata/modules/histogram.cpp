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
/*! \internal \file
 * \brief
 * Implements classes in histogram.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/modules/histogram.h"

#include <cmath>

#include <limits>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "histogram-impl.h"

static const real UNDEFINED = std::numeric_limits<real>::max();
static bool isDefined(real value)
{
    return value != UNDEFINED;
}

namespace gmx
{

/********************************************************************
 * AnalysisHistogramSettingsInitializer
 */

AnalysisHistogramSettingsInitializer::AnalysisHistogramSettingsInitializer()
    : min_(UNDEFINED), max_(UNDEFINED), binWidth_(UNDEFINED),
      binCount_(0), bIntegerBins_(false), bRoundRange_(false),
      bIncludeAll_(false)
{
}


/********************************************************************
 * AnalysisHistogramSettings
 */

AnalysisHistogramSettings::AnalysisHistogramSettings()
    : firstEdge_(0.0), lastEdge_(0.0), binWidth_(0.0), inverseBinWidth_(0.0),
      binCount_(0), bAll_(false)
{
}


AnalysisHistogramSettings::AnalysisHistogramSettings(
        const AnalysisHistogramSettingsInitializer &settings)
{
    GMX_RELEASE_ASSERT(isDefined(settings.min_),
                       "Histogram start value must be defined");
    GMX_RELEASE_ASSERT(!isDefined(settings.max_) || settings.max_ > settings.min_,
                       "Histogram end value must be larger than start value");
    GMX_RELEASE_ASSERT(!isDefined(settings.binWidth_) || settings.binWidth_ > 0.0,
                       "Histogram bin width must be positive");
    GMX_RELEASE_ASSERT(settings.binCount_ >= 0,
                       "Histogram bin count must be positive");

    if (!isDefined(settings.max_))
    {
        GMX_RELEASE_ASSERT(isDefined(settings.binWidth_) && settings.binCount_ > 0,
                           "Not all required values provided");
        GMX_RELEASE_ASSERT(!settings.bRoundRange_,
                           "Rounding only supported for min/max ranges");

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
            firstEdge_ = binWidth_ * floor(settings.min_ / binWidth_);
            lastEdge_  = binWidth_ * ceil(settings.max_ / binWidth_);
            binCount_  = static_cast<int>((lastEdge_ - firstEdge_) / binWidth_ + 0.5);
        }
        else
        {
            firstEdge_     = settings.min_;
            lastEdge_     = settings.max_;
            if (settings.binCount_ > 0)
            {
                binCount_ = settings.binCount_;
                if (settings.bIntegerBins_)
                {
                    GMX_RELEASE_ASSERT(settings.binCount_ > 1,
                                       "Bin count must be at least two with integer bins");
                    binWidth_   = (lastEdge_ - firstEdge_) / (binCount_ - 1);
                    firstEdge_ -= 0.5 * binWidth_;
                    lastEdge_  += 0.5 * binWidth_;
                }
                else
                {
                    binWidth_ = (lastEdge_ - firstEdge_) / binCount_;
                }
            }
            else
            {
                binWidth_ = settings.binWidth_;
                binCount_ = static_cast<int>((lastEdge_ - firstEdge_) / binWidth_ + 0.5);
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


int
AnalysisHistogramSettings::findBin(real y) const
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
 * AbstractAverageHistogram
 */

AbstractAverageHistogram::AbstractAverageHistogram()
{
}


AbstractAverageHistogram::AbstractAverageHistogram(
        const AnalysisHistogramSettings &settings)
    : settings_(settings)
{
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(),
             settings.binWidth());
}


AbstractAverageHistogram::~AbstractAverageHistogram()
{
}


void
AbstractAverageHistogram::init(const AnalysisHistogramSettings &settings)
{
    settings_ = settings;
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(),
             settings.binWidth());
}


AverageHistogramPointer
AbstractAverageHistogram::resampleDoubleBinWidth(bool bIntegerBins) const
{
    int nbins;
    if (bIntegerBins)
    {
        nbins = (rowCount() + 1) / 2;
    }
    else
    {
        nbins = rowCount() / 2;
    }

    AverageHistogramPointer dest(
        new internal::StaticAverageHistogram(
            histogramFromBins(xstart(), nbins, 2*xstep())
                .integerBins(bIntegerBins)));
    dest->setColumnCount(columnCount());
    dest->allocateValues();

    int  i, j;
    for (i = j = 0; i < nbins; ++i)
    {
        const bool bFirstHalfBin = (bIntegerBins && i == 0);
        for (int c = 0; c < columnCount(); ++c)
        {
            real  v1, v2;
            if (bFirstHalfBin)
            {
                v1 = value(0, c);
                v2 = 0;
            }
            else
            {
                v1 = value(j, c);
                v2 = value(j + 1, c);
            }
            if (c == 1)
            {
                dest->setValue(i, c, sqrt(v1 * v1 + v2 * v2));
            }
            else
            {
                dest->setValue(i, c, v1 + v2);
            }
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


AverageHistogramPointer
AbstractAverageHistogram::clone() const
{
    AverageHistogramPointer dest(new internal::StaticAverageHistogram());
    copyContents(this, dest.get());
    return dest;
}


void
AbstractAverageHistogram::normalizeProbability()
{
    real sum = 0;
    for (int i = 0; i < rowCount(); ++i)
    {
        sum += value(i, 0);
    }
    scale(1.0 / (sum * xstep()));
}


void
AbstractAverageHistogram::scale(real norm)
{
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, 0) *= norm;
        value(i, 1) *= norm;
    }
}


void
AbstractAverageHistogram::scaleVector(real norm[])
{
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, 0) *= norm[i];
        value(i, 1) *= norm[i];
    }
}


/********************************************************************
 * StaticAverageHistogram
 */

namespace internal
{

StaticAverageHistogram::StaticAverageHistogram()
{
}


StaticAverageHistogram::StaticAverageHistogram(
        const AnalysisHistogramSettings &settings)
    : AbstractAverageHistogram(settings)
{
}


/********************************************************************
 * BasicAverageHistogramModule
 */

BasicAverageHistogramModule::BasicAverageHistogramModule()
    : frameCount_(0)
{
    setColumnCount(2);
}


BasicAverageHistogramModule::BasicAverageHistogramModule(
        const AnalysisHistogramSettings &settings)
    : AbstractAverageHistogram(settings), frameCount_(0)
{
    setColumnCount(2);
}


int
BasicAverageHistogramModule::flags() const
{
    return efAllowMulticolumn;
}


void
BasicAverageHistogramModule::dataStarted(AbstractAnalysisData *data)
{
    GMX_RELEASE_ASSERT(rowCount() == data->columnCount(),
                       "Inconsistent data sizes, something is wrong in the initialization");
    allocateValues();
}


void
BasicAverageHistogramModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
}


void
BasicAverageHistogramModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    int firstcol = points.firstColumn();
    for (int i = 0; i < points.columnCount(); ++i)
    {
        real y = points.y(i);
        value(firstcol + i, 0) += y;
        value(firstcol + i, 1) += y * y;
    }
}


void
BasicAverageHistogramModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
    ++frameCount_;
}


void
BasicAverageHistogramModule::dataFinished()
{
    for (int i = 0; i < rowCount(); ++i)
    {
        real ave = value(i, 0) / frameCount_;
        real std = sqrt(value(i, 1) / frameCount_ - ave * ave);
        setValue(i, 0, ave);
        setValue(i, 1, std);
    }
}


/********************************************************************
 * BasicHistogramImpl
 */

BasicHistogramImpl::BasicHistogramImpl()
    : averager_(new BasicAverageHistogramModule())
{
}


BasicHistogramImpl::BasicHistogramImpl(const AnalysisHistogramSettings &settings)
    : settings_(settings), averager_(new BasicAverageHistogramModule(settings))
{
}


BasicHistogramImpl::~BasicHistogramImpl()
{
}


void BasicHistogramImpl::init(const AnalysisHistogramSettings &settings)
{
    settings_ = settings;
    averager_->init(settings);
}


void
BasicHistogramImpl::initFrame(AnalysisDataStorageFrame *frame)
{
    for (int i = 0; i < frame->columnCount(); ++i)
    {
        frame->setValue(i, 0.0);
    }
}

} // namespace internal


/********************************************************************
 * AnalysisDataSimpleHistogramModule
 */

AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule()
    : impl_(new internal::BasicHistogramImpl())
{
}


AnalysisDataSimpleHistogramModule::AnalysisDataSimpleHistogramModule(
        const AnalysisHistogramSettings &settings)
    : impl_(new internal::BasicHistogramImpl(settings))
{
}


AnalysisDataSimpleHistogramModule::~AnalysisDataSimpleHistogramModule()
{
}


void AnalysisDataSimpleHistogramModule::init(const AnalysisHistogramSettings &settings)
{
    impl_->init(settings);
}


AbstractAverageHistogram &
AnalysisDataSimpleHistogramModule::averager()
{
    return *impl_->averager_;
}


const AnalysisHistogramSettings &
AnalysisDataSimpleHistogramModule::settings() const
{
    return impl_->settings_;
}


int
AnalysisDataSimpleHistogramModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void
AnalysisDataSimpleHistogramModule::dataStarted(AbstractAnalysisData *data)
{
    addModule(impl_->averager_);
    setColumnCount(settings().binCount());
    notifyDataStart();
    impl_->storage_.startDataStorage(this);
}


void
AnalysisDataSimpleHistogramModule::frameStarted(const AnalysisDataFrameHeader &header)
{
    AnalysisDataStorageFrame &frame = impl_->storage_.startFrame(header);
    impl_->initFrame(&frame);
}


void
AnalysisDataSimpleHistogramModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    AnalysisDataStorageFrame &frame =
        impl_->storage_.currentFrame(points.frameIndex());
    for (int i = 0; i < points.columnCount(); ++i)
    {
        int bin = settings().findBin(points.y(i));
        if (bin != -1)
        {
            frame.value(bin) += 1;
        }
    }
}


void
AnalysisDataSimpleHistogramModule::frameFinished(const AnalysisDataFrameHeader &header)
{
    impl_->storage_.finishFrame(header.index());
}


void
AnalysisDataSimpleHistogramModule::dataFinished()
{
    notifyDataFinish();
}


AnalysisDataFrameRef
AnalysisDataSimpleHistogramModule::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool
AnalysisDataSimpleHistogramModule::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataWeightedHistogramModule
 */

AnalysisDataWeightedHistogramModule::AnalysisDataWeightedHistogramModule()
    : impl_(new internal::BasicHistogramImpl())
{
}


AnalysisDataWeightedHistogramModule::AnalysisDataWeightedHistogramModule(
        const AnalysisHistogramSettings &settings)
    : impl_(new internal::BasicHistogramImpl(settings))
{
}


AnalysisDataWeightedHistogramModule::~AnalysisDataWeightedHistogramModule()
{
}


void AnalysisDataWeightedHistogramModule::init(const AnalysisHistogramSettings &settings)
{
    impl_->init(settings);
}


AbstractAverageHistogram &
AnalysisDataWeightedHistogramModule::averager()
{
    return *impl_->averager_;
}


const AnalysisHistogramSettings &
AnalysisDataWeightedHistogramModule::settings() const
{
    return impl_->settings_;
}


int
AnalysisDataWeightedHistogramModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void
AnalysisDataWeightedHistogramModule::dataStarted(AbstractAnalysisData *data)
{
    addModule(impl_->averager_);
    setColumnCount(settings().binCount());
    notifyDataStart();
    impl_->storage_.startDataStorage(this);
}


void
AnalysisDataWeightedHistogramModule::frameStarted(const AnalysisDataFrameHeader &header)
{
    AnalysisDataStorageFrame &frame = impl_->storage_.startFrame(header);
    impl_->initFrame(&frame);
}


void
AnalysisDataWeightedHistogramModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (points.firstColumn() != 0 || points.columnCount() < 2)
    {
        GMX_THROW(APIError("Invalid data layout"));
    }
    int bin = settings().findBin(points.y(0));
    if (bin != -1)
    {
        AnalysisDataStorageFrame &frame =
            impl_->storage_.currentFrame(points.frameIndex());
        for (int i = 1; i < points.columnCount(); ++i)
        {
            frame.value(bin) += points.y(i);
        }
    }
}


void
AnalysisDataWeightedHistogramModule::frameFinished(const AnalysisDataFrameHeader &header)
{
    impl_->storage_.finishFrame(header.index());
}


void
AnalysisDataWeightedHistogramModule::dataFinished()
{
    notifyDataFinish();
}


AnalysisDataFrameRef
AnalysisDataWeightedHistogramModule::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}


bool
AnalysisDataWeightedHistogramModule::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}


/********************************************************************
 * AnalysisDataBinAverageModule
 */

AnalysisDataBinAverageModule::AnalysisDataBinAverageModule()
{
    setColumnCount(3);
}


AnalysisDataBinAverageModule::AnalysisDataBinAverageModule(
        const AnalysisHistogramSettings &settings)
    : settings_(settings)
{
    setColumnCount(3);
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(),
             settings.binWidth());
}


AnalysisDataBinAverageModule::~AnalysisDataBinAverageModule()
{
}


void
AnalysisDataBinAverageModule::init(const AnalysisHistogramSettings &settings)
{
    settings_ = settings;
    setRowCount(settings.binCount());
    setXAxis(settings.firstEdge() + 0.5 * settings.binWidth(),
             settings.binWidth());
}


int
AnalysisDataBinAverageModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void
AnalysisDataBinAverageModule::dataStarted(AbstractAnalysisData * /*data*/)
{
    allocateValues();
}


void
AnalysisDataBinAverageModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
}


void
AnalysisDataBinAverageModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (points.firstColumn() != 0 || points.columnCount() < 2)
    {
        GMX_THROW(APIError("Invalid data layout"));
    }
    int bin = settings().findBin(points.y(0));
    if (bin != -1)
    {
        for (int i = 1; i < points.columnCount(); ++i)
        {
            real y = points.y(i);
            value(bin, 0) += y;
            value(bin, 1) += y * y;
        }
        value(bin, 2) += points.columnCount() - 1;
    }
}


void
AnalysisDataBinAverageModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
}


void
AnalysisDataBinAverageModule::dataFinished()
{
    for (int i = 0; i < settings().binCount(); ++i)
    {
        real n = value(i, 2);
        if (n > 0)
        {
            real ave = value(i, 0) / n;
            real std = sqrt(value(i, 1) / n - ave * ave);
            setValue(i, 0, ave);
            setValue(i, 1, std);
        }
    }
    valuesReady();
}

} // namespace gmx
