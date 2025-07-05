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
 * Implements gmx::AnalysisDataAverageModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/average.h"

#include <cmath>

#include <utility>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "frameaverager.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataAverageModule
 */

class AnalysisDataAverageModule::Impl
{
public:
    Impl() : bDataSets_(false) {}

    //! Averaging helper objects for each input data set.
    std::vector<AnalysisDataFrameAverager> averagers_;
    //! Whether to average all columns in a data set into a single value.
    bool bDataSets_;
};

AnalysisDataAverageModule::AnalysisDataAverageModule() : impl_(new Impl()) {}

AnalysisDataAverageModule::~AnalysisDataAverageModule() {}

void AnalysisDataAverageModule::setAverageDataSets(bool bDataSets)
{
    impl_->bDataSets_ = bDataSets;
}

int AnalysisDataAverageModule::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing | efAllowMultipleDataSets;
}

void AnalysisDataAverageModule::dataStarted(AbstractAnalysisData* data)
{
    if (impl_->bDataSets_)
    {
        setColumnCount(1);
        setRowCount(data->dataSetCount());
        impl_->averagers_.resize(1);
        impl_->averagers_[0].setColumnCount(data->dataSetCount());
    }
    else
    {
        setColumnCount(data->dataSetCount());
        impl_->averagers_.resize(data->dataSetCount());
        size_t rowCount = 0;
        for (size_t i = 0; i < data->dataSetCount(); ++i)
        {
            impl_->averagers_[i].setColumnCount(data->columnCount(i));
            rowCount = std::max(rowCount, data->columnCount(i));
        }
        setRowCount(rowCount);
    }
}

void AnalysisDataAverageModule::frameStarted(const AnalysisDataFrameHeader& /*header*/) {}

void AnalysisDataAverageModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (impl_->bDataSets_)
    {
        const int dataSet = points.dataSetIndex();
        for (int i = 0; i < points.columnCount(); ++i)
        {
            if (points.present(i))
            {
                impl_->averagers_[0].addValue(dataSet, points.y(i));
            }
        }
    }
    else
    {
        impl_->averagers_[points.dataSetIndex()].addPoints(points);
    }
}

void AnalysisDataAverageModule::frameFinished(const AnalysisDataFrameHeader& /*header*/) {}

void AnalysisDataAverageModule::dataFinished()
{
    allocateValues();
    for (size_t i = 0; i < columnCount(); ++i)
    {
        impl_->averagers_[i].finish();
        size_t j = 0;
        for (; j < impl_->averagers_[i].columnCount(); ++j)
        {
            value(j, i).setValue(impl_->averagers_[i].average(j),
                                 std::sqrt(impl_->averagers_[i].variance(j)));
        }
        for (; j < rowCount(); ++j)
        {
            value(j, i).setValue(0.0, 0.0, false);
        }
    }
    valuesReady();
}

real AnalysisDataAverageModule::average(int dataSet, int column) const
{
    if (impl_->bDataSets_)
    {
        GMX_ASSERT(column == 0, "Column should be zero with setAverageDataSets(true)");
        std::swap(dataSet, column);
    }
    return value(column, dataSet).value();
}

real AnalysisDataAverageModule::standardDeviation(int dataSet, int column) const
{
    if (impl_->bDataSets_)
    {
        GMX_ASSERT(column == 0, "Column should be zero with setAverageDataSets(true)");
        std::swap(dataSet, column);
    }
    return value(column, dataSet).error();
}

int AnalysisDataAverageModule::sampleCount(int dataSet, int column) const
{
    if (impl_->bDataSets_)
    {
        GMX_ASSERT(column == 0, "Column should be zero with setAverageDataSets(true)");
        std::swap(dataSet, column);
    }
    return impl_->averagers_[dataSet].sampleCount(column);
}


/********************************************************************
 * AnalysisDataFrameAverageModule
 */

class AnalysisDataFrameAverageModule::Impl
{
public:
    //! Storage implementation object.
    AnalysisDataStorage storage_;
    //! Number of samples in a frame for each data set.
    std::vector<int> sampleCount_;
};

AnalysisDataFrameAverageModule::AnalysisDataFrameAverageModule() : impl_(new Impl()) {}

AnalysisDataFrameAverageModule::~AnalysisDataFrameAverageModule() {}

size_t AnalysisDataFrameAverageModule::frameCount() const
{
    return impl_->storage_.frameCount();
}

int AnalysisDataFrameAverageModule::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing | efAllowMultipleDataSets;
}

void AnalysisDataFrameAverageModule::dataStarted(AbstractAnalysisData* data)
{
    setColumnCount(0, data->dataSetCount());
    impl_->sampleCount_.resize(data->dataSetCount());
    impl_->storage_.startDataStorage(this, &moduleManager());
}

void AnalysisDataFrameAverageModule::frameStarted(const AnalysisDataFrameHeader& header)
{
    AnalysisDataStorageFrame& frame = impl_->storage_.startFrame(header);
    for (size_t i = 0; i < columnCount(); ++i)
    {
        impl_->sampleCount_[i] = 0;
        frame.setValue(i, 0.0);
    }
}

void AnalysisDataFrameAverageModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    const int                 dataSet = points.dataSetIndex();
    AnalysisDataStorageFrame& frame   = impl_->storage_.currentFrame(points.frameIndex());
    for (int i = 0; i < points.columnCount(); ++i)
    {
        if (points.present(i))
        {
            // TODO: Consider using AnalysisDataFrameAverager
            const real y     = points.y(i);
            const real delta = y - frame.value(dataSet);
            impl_->sampleCount_[dataSet] += 1;
            frame.value(dataSet) += delta / impl_->sampleCount_[dataSet];
        }
    }
}

void AnalysisDataFrameAverageModule::frameFinished(const AnalysisDataFrameHeader& header)
{
    impl_->storage_.finishFrame(header.index());
}

void AnalysisDataFrameAverageModule::dataFinished()
{
    impl_->storage_.finishDataStorage();
}

AnalysisDataFrameRef AnalysisDataFrameAverageModule::tryGetDataFrameInternal(size_t index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}

bool AnalysisDataFrameAverageModule::requestStorageInternal(size_t nframes)
{
    return impl_->storage_.requestStorage(nframes);
}

} // namespace gmx
