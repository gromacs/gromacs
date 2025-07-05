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
 * Implements gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/abstractdata.h"

#include <memory>
#include <vector>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/datamodulemanager.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "dataproxy.h"

namespace gmx
{

/********************************************************************
 * AbstractAnalysisData::Impl
 */

/*! \internal \brief
 * Private implementation class for AbstractAnalysisData.
 *
 * \ingroup module_analysisdata
 */
class AbstractAnalysisData::Impl
{
public:
    Impl();

    //! Column counts for each data set in the data.
    std::vector<int> columnCounts_;
    //! Whether the data is multipoint.
    bool bMultipoint_;
    //! Manager for the added modules.
    AnalysisDataModuleManager modules_;
};

AbstractAnalysisData::Impl::Impl() : bMultipoint_(false)
{
    columnCounts_.push_back(0);
}


/********************************************************************
 * AbstractAnalysisData
 */
/*! \cond libapi */
AbstractAnalysisData::AbstractAnalysisData() : impl_(new Impl()) {}
//! \endcond

AbstractAnalysisData::~AbstractAnalysisData() {}

bool AbstractAnalysisData::isMultipoint() const
{
    return impl_->bMultipoint_;
}

size_t AbstractAnalysisData::dataSetCount() const
{
    return impl_->columnCounts_.size();
}

size_t AbstractAnalysisData::columnCount(size_t dataSet) const
{
    GMX_ASSERT(dataSet < dataSetCount(), "Out of range data set index");
    return impl_->columnCounts_[dataSet];
}

size_t AbstractAnalysisData::columnCount() const
{
    GMX_ASSERT(dataSetCount() == 1, "Convenience method not available for multiple data sets");
    return columnCount(0);
}


AnalysisDataFrameRef AbstractAnalysisData::tryGetDataFrame(size_t index) const
{
    if (index >= frameCount())
    {
        return AnalysisDataFrameRef();
    }
    return tryGetDataFrameInternal(index);
}


AnalysisDataFrameRef AbstractAnalysisData::getDataFrame(size_t index) const
{
    AnalysisDataFrameRef frame = tryGetDataFrame(index);
    if (!frame.isValid())
    {
        GMX_THROW(APIError("Invalid frame accessed"));
    }
    return frame;
}


bool AbstractAnalysisData::requestStorage(int nframes)
{
    GMX_RELEASE_ASSERT(nframes >= -1, "Invalid number of frames requested");
    if (nframes == 0)
    {
        return true;
    }
    return requestStorageInternal(nframes);
}


void AbstractAnalysisData::addModule(const AnalysisDataModulePointer& module)
{
    impl_->modules_.addModule(this, module);
}


void AbstractAnalysisData::addColumnModule(size_t col, size_t span, const AnalysisDataModulePointer& module)
{
    GMX_RELEASE_ASSERT(span >= 1, "Invalid columns specified for a column module");
    std::shared_ptr<AnalysisDataProxy> proxy(new AnalysisDataProxy(col, span, this));
    proxy->addModule(module);
    addModule(proxy);
}


void AbstractAnalysisData::applyModule(IAnalysisDataModule* module)
{
    impl_->modules_.applyModule(this, module);
}

/*! \cond libapi */
void AbstractAnalysisData::setDataSetCount(size_t dataSetCount)
{
    GMX_RELEASE_ASSERT(dataSetCount > 0, "Invalid data column count");
    impl_->modules_.dataPropertyAboutToChange(AnalysisDataModuleManager::eMultipleDataSets,
                                              dataSetCount > 1);
    impl_->columnCounts_.resize(dataSetCount);
}

void AbstractAnalysisData::setColumnCount(size_t dataSet, size_t columnCount)
{
    GMX_RELEASE_ASSERT(dataSet < dataSetCount(), "Out of range data set index");
    GMX_RELEASE_ASSERT(columnCount > 0, "Invalid data column count");

    bool bMultipleColumns = columnCount > 1;
    for (size_t i = 0; i < dataSetCount() && !bMultipleColumns; ++i)
    {
        if (i != dataSet && this->columnCount(i) > 1)
        {
            bMultipleColumns = true;
        }
    }
    impl_->modules_.dataPropertyAboutToChange(AnalysisDataModuleManager::eMultipleColumns, bMultipleColumns);
    impl_->columnCounts_[dataSet] = columnCount;
}

void AbstractAnalysisData::setMultipoint(bool bMultipoint)
{
    impl_->modules_.dataPropertyAboutToChange(AnalysisDataModuleManager::eMultipoint, bMultipoint);
    impl_->bMultipoint_ = bMultipoint;
}

AnalysisDataModuleManager& AbstractAnalysisData::moduleManager()
{
    return impl_->modules_;
}

const AnalysisDataModuleManager& AbstractAnalysisData::moduleManager() const
{
    return impl_->modules_;
}
//! \endcond

} // namespace gmx
