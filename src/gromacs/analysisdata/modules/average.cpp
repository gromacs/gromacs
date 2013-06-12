/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements gmx::AnalysisDataAverageModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "average.h"

#include <cmath>

#include <vector>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"

#include "frameaverager.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataAverageModule
 */

class AnalysisDataAverageModule::Impl
{
    public:
        //! Averaging helper object.
        AnalysisDataFrameAverager averager_;
};

AnalysisDataAverageModule::AnalysisDataAverageModule()
    : impl_(new Impl())
{
    setColumnCount(2);
}

AnalysisDataAverageModule::~AnalysisDataAverageModule()
{
}

int
AnalysisDataAverageModule::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing;
}

void
AnalysisDataAverageModule::dataStarted(AbstractAnalysisData *data)
{
    const int valueCount = data->columnCount();
    impl_->averager_.setColumnCount(valueCount);
    setRowCount(valueCount);
}

void
AnalysisDataAverageModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
}

void
AnalysisDataAverageModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    impl_->averager_.addPoints(points);
}

void
AnalysisDataAverageModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
}

void
AnalysisDataAverageModule::dataFinished()
{
    impl_->averager_.finish();
    allocateValues();
    for (int i = 0; i < rowCount(); ++i)
    {
        value(i, 0).setValue(impl_->averager_.average(i));
        value(i, 1).setValue(std::sqrt(impl_->averager_.variance(i)));
    }
    valuesReady();
}

real
AnalysisDataAverageModule::average(int index) const
{
    return value(index, 0).value();
}

real
AnalysisDataAverageModule::stddev(int index) const
{
    return value(index, 1).value();
}


/********************************************************************
 * AnalysisDataFrameAverageModule
 */

class AnalysisDataFrameAverageModule::Impl
{
    public:
        //! Storage implementation object.
        AnalysisDataStorage     storage_;
        //! Number of samples in a frame for each data set.
        std::vector<int>        sampleCount_;
};

AnalysisDataFrameAverageModule::AnalysisDataFrameAverageModule()
    : impl_(new Impl())
{
}

AnalysisDataFrameAverageModule::~AnalysisDataFrameAverageModule()
{
}

int
AnalysisDataFrameAverageModule::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing
           | efAllowMultipleDataSets;
}

void
AnalysisDataFrameAverageModule::dataStarted(AbstractAnalysisData *data)
{
    setColumnCount(0, data->dataSetCount());
    impl_->sampleCount_.resize(data->dataSetCount());
    notifyDataStart();
    impl_->storage_.startDataStorage(this);
}

void
AnalysisDataFrameAverageModule::frameStarted(const AnalysisDataFrameHeader &header)
{
    AnalysisDataStorageFrame &frame = impl_->storage_.startFrame(header);
    for (int i = 0; i < columnCount(); ++i)
    {
        impl_->sampleCount_[i] = 0;
        frame.setValue(i, 0.0);
    }
}

void
AnalysisDataFrameAverageModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    const int                 dataSet = points.dataSetIndex();
    AnalysisDataStorageFrame &frame   =
        impl_->storage_.currentFrame(points.frameIndex());
    for (int i = 0; i < points.columnCount(); ++i)
    {
        if (points.present(i))
        {
            // TODO: Consider using AnalysisDataFrameAverager
            const real y     = points.y(i);
            const real delta = y - frame.value(dataSet);
            impl_->sampleCount_[dataSet] += 1;
            frame.value(dataSet)         += delta / impl_->sampleCount_[dataSet];
        }
    }
}

void
AnalysisDataFrameAverageModule::frameFinished(const AnalysisDataFrameHeader &header)
{
    impl_->storage_.finishFrame(header.index());
}

void
AnalysisDataFrameAverageModule::dataFinished()
{
    notifyDataFinish();
}

AnalysisDataFrameRef
AnalysisDataFrameAverageModule::tryGetDataFrameInternal(int index) const
{
    return impl_->storage_.tryGetDataFrame(index);
}

bool
AnalysisDataFrameAverageModule::requestStorageInternal(int nframes)
{
    return impl_->storage_.requestStorage(nframes);
}

} // namespace gmx
