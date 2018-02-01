/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes in write.h.
 *
 * \ingroup module_analysisdata
 * \author
 */
#include "gmxpre.h"

#include "write.h"

#include <cstdio>
#include <cstring>

#include <string>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datastorage.h"
#include "gromacs/analysisdata/framelocaldata.h"
#include "gromacs/analysisdata/modules/filehandler.h"
#include "gromacs/analysisdata/modules/framehandler.h"
#include "gromacs/analysisdata/modules/settings.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace
{

} // namespace

namespace gmx
{

/********************************************************************
 * AbstractWriteModule::Impl
 */

class AbstractWriteModule::Impl
{
    public:
        explicit Impl(const TrajectoryDataWriteSettings &settings);
        ~Impl();

        void init(const TrajectoryDataWriteSettings &settings);
        void closeAllFiles();

        TrajectoryDataWriteSettings        settings_;
        const Selection                   *sel_;

        typedef AnalysisDataFrameLocalData<t_trxframe> FrameLocalData;
        FrameLocalData            coordinate_;

        std::vector<Filehandler>  filehandler_;
        std::vector<Framehandler> framehandler_;


};

AbstractWriteModule::Impl::Impl(const TrajectoryDataWriteSettings &settings)
    : settings_(settings)
{
}

AbstractWriteModule::Impl::~Impl()
{
    // finish the handling of filehandler objects?
    closeAllFiles();
}

void
AbstractWriteModule::Impl::init(const TrajectoryDataWriteSettings &settings)
{
    settings_ = settings;
}

void
AbstractWriteModule::Impl::closeAllFiles()
{
    for (int i = 0; i < static_cast<int>(filehandler_.size()); ++i)
    {
        filehandler_[i].closeFile();
    }
}

/********************************************************************
 * AbstractWriteModule
 */
/*! \cond libapi */
AbstractWriteModule::AbstractWriteModule()
    : impl_(new Impl(TrajectoryDataWriteSettings()))
{
}

AbstractWriteModule::AbstractWriteModule(const TrajectoryDataWriteSettings &settings)
    : impl_(new Impl(settings))
{
}
//! \endcond

AbstractWriteModule::~AbstractWriteModule()
{
}


void
AbstractWriteModule::setSettings(const TrajectoryDataWriteSettings &settings)
{
    impl_->settings_ = settings;
}

void
AbstractWriteModule::setExternal(const Selection *sel, std::string name, const gmx_mtop_t *mtop, const t_topology *top)
{
    impl_->settings_.setName(name);
    impl_->settings_.setInputSel(sel);
    impl_->settings_.setTopology(mtop, top);
}

bool
AbstractWriteModule::parallelDataStarted(AbstractAnalysisData              *data,
                                         const AnalysisDataParallelOptions &options)
{
    setColumnCount(data->dataSetCount());
    // for each data set, prepare an object for the frame data changing
    impl_->filehandler_.resize(data->dataSetCount());
    impl_->framehandler_.resize(data->dataSetCount());
    impl_->coordinate_.setDataSetCount(data->dataSetCount());
    // reset info in those new objects
    for (int i = 0; i < data->dataSetCount(); ++i)
    {
        impl_->filehandler_[i].setSettings(&impl_->settings_);
        impl_->filehandler_[i].openFile();
        // openFile needs to throw so that no checks are needed later
        impl_->framehandler_[i].setSettings(&impl_->settings_);
        impl_->coordinate_.setColumnCount(i, 1);
    }
    impl_->coordinate_.init(options);
    return true;
}

void
AbstractWriteModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
    // initialize empty storage for coordinate frames

}

void
AbstractWriteModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    Impl::FrameLocalData::DataSetHandle handle
        = impl_->coordinate_.frameDataSet(points.frameIndex(), points.dataSetIndex());
    for (int i = 0; i < points.columnCount(); i++)
    {
        if (points.present(i))
        {
            t_trxframe local;
            clear_trxframe(&local, true);
            t_trxframe old = points.values()[i].valueAsVariant().cast<t_trxframe>();
            impl_->framehandler_[points.dataSetIndex()].modifyFrame(&local, &old);
            handle.value(i) = local;
        }
    }
}

void
AbstractWriteModule::frameFinished(const AnalysisDataFrameHeader &header)
{
    Impl::FrameLocalData::FrameHandle handle = impl_->coordinate_.frameData(header.index());
    for (int i = 0; i < dataSetCount(); ++i)
    {
        Impl::FrameLocalData::DataSetHandle dataSet = handle.dataSet(i);
        if (impl_->filehandler_[i].isFileOpen())
        {
            impl_->filehandler_[i].writeValue(dataSet.value(0));
        }
    }
}

void
AbstractWriteModule::frameFinishedSerial(int /*frameIndex*/)
{
}


void
AbstractWriteModule::dataFinished()
{

    impl_->closeAllFiles();
}


int
AbstractWriteModule::flags() const
{
    return efAllowMissing | efAllowMulticolumn | efAllowMultipoint
           | efAllowMultipleDataSets;
}


//! \endcond

/********************************************************************
 * DataTrjWriteModule
 */

TrajectoryDataWriteModule::TrajectoryDataWriteModule()
{
}

TrajectoryDataWriteModule::TrajectoryDataWriteModule(
        const TrajectoryDataWriteSettings &settings)
    : AbstractWriteModule(settings)
{
}

} // namespace gmx
