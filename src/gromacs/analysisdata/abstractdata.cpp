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
 * Implements gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/abstractdata.h"

#include <vector>

#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/uniqueptr.h"

#include "dataframe.h"
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
        //! Shorthand for list of modules added to the data.
        typedef std::vector<AnalysisDataModulePointer> ModuleList;

        Impl();

        //! Returns whether any data set has more than one column.
        bool isMultiColumn() const;

        /*! \brief
         * Checks whether a module is compatible with the data properties.
         *
         * \param[in] module Module to check.
         * \throws    APIError if \p module is not compatible with the data.
         *
         * Does not check the actual data (e.g., missing values), but only the
         * dimensionality and other preset properties of the data.
         */
        void checkModuleProperties(const AnalysisDataModuleInterface &module) const;

        /*! \brief
         * Present data already added to the data object to a module.
         *
         * \param[in] data   Data object to read data from.
         * \param[in] module Module to present the data to.
         * \throws    APIError if \p module is not compatible with the data.
         * \throws    APIError if all data is not available through
         *      getDataFrame().
         * \throws    unspecified Any exception thrown by \p module in its data
         *      notification methods.
         *
         * Uses getDataFrame() in \p data to access all data in the object, and
         * calls the notification functions in \p module as if the module had
         * been registered to the data object when the data was added.
         */
        void presentData(AbstractAnalysisData        *data,
                         AnalysisDataModuleInterface *module);

        //! Column counts for each data set in the data.
        std::vector<int>        columnCounts_;
        //! Whether the data is multipoint.
        bool                    bMultipoint_;
        //! List of modules added to the data.
        ModuleList              modules_;
        //! true if all modules support missing data.
        bool                    bAllowMissing_;
        //! Whether notifyDataStart() has been called.
        mutable bool            bDataStart_;
        //! Whether new data is being added.
        mutable bool            bInData_;
        //! Whether data for a frame is being added.
        mutable bool            bInFrame_;
        //! Index of the currently active frame.
        mutable int             currIndex_;
        /*! \brief
         * Total number of frames in the data.
         *
         * The counter is incremented in notifyFrameFinish().
         */
        int                     nframes_;
};

AbstractAnalysisData::Impl::Impl()
    : bMultipoint_(false), bAllowMissing_(true),
      bDataStart_(false), bInData_(false), bInFrame_(false),
      currIndex_(-1), nframes_(0)
{
    columnCounts_.push_back(0);
}

bool
AbstractAnalysisData::Impl::isMultiColumn() const
{
    std::vector<int>::const_iterator i;
    for (i = columnCounts_.begin(); i != columnCounts_.end(); ++i)
    {
        if (*i > 1)
        {
            return true;
        }
    }
    return false;
}

//! Helper macro for testing module flags.
#define TEST_MODULE_FLAG(flags, flagname) \
    ((flags) & AnalysisDataModuleInterface::flagname)
void
AbstractAnalysisData::Impl::checkModuleProperties(
        const AnalysisDataModuleInterface &module) const
{
    const int flags = module.flags();
    if ((!TEST_MODULE_FLAG(flags, efAllowMulticolumn) && isMultiColumn()) ||
        (!TEST_MODULE_FLAG(flags, efAllowMultipoint)  && bMultipoint_) ||
        ( TEST_MODULE_FLAG(flags, efOnlyMultipoint)   && !bMultipoint_) ||
        (!TEST_MODULE_FLAG(flags, efAllowMultipleDataSets)
         && columnCounts_.size() > 1U))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }
}
#undef TEST_MODULE_FLAGS

void
AbstractAnalysisData::Impl::presentData(AbstractAnalysisData        *data,
                                        AnalysisDataModuleInterface *module)
{
    module->dataStarted(data);
    bool bCheckMissing = bAllowMissing_
        && !(module->flags() & AnalysisDataModuleInterface::efAllowMissing);
    for (int i = 0; i < data->frameCount(); ++i)
    {
        AnalysisDataFrameRef frame = data->getDataFrame(i);
        GMX_RELEASE_ASSERT(frame.isValid(), "Invalid data frame returned");
        // TODO: Check all frames before doing anything for slightly better
        // exception behavior.
        if (bCheckMissing && !frame.allPresent())
        {
            GMX_THROW(APIError("Missing data not supported by a module"));
        }
        module->frameStarted(frame.header());
        for (int j = 0; j < frame.pointSetCount(); ++j)
        {
            module->pointsAdded(frame.pointSet(j));
        }
        module->frameFinished(frame.header());
    }
    if (!bInData_)
    {
        module->dataFinished();
    }
}


/********************************************************************
 * AbstractAnalysisData
 */
/*! \cond libapi */
AbstractAnalysisData::AbstractAnalysisData()
    : impl_(new Impl())
{
}
//! \endcond

AbstractAnalysisData::~AbstractAnalysisData()
{
}

bool
AbstractAnalysisData::isMultipoint() const
{
    return impl_->bMultipoint_;
}

int
AbstractAnalysisData::dataSetCount() const
{
    return impl_->columnCounts_.size();
}

int
AbstractAnalysisData::columnCount(int dataSet) const
{
    GMX_ASSERT(dataSet >= 0 && dataSet < dataSetCount(),
               "Out of range data set index");
    return impl_->columnCounts_[dataSet];
}

int
AbstractAnalysisData::columnCount() const
{
    GMX_ASSERT(dataSetCount() == 1,
               "Convenience method not available for multiple data sets");
    return columnCount(0);
}

int
AbstractAnalysisData::frameCount() const
{
    return impl_->nframes_;
}


AnalysisDataFrameRef
AbstractAnalysisData::tryGetDataFrame(int index) const
{
    if (index < 0 || index >= frameCount())
    {
        return AnalysisDataFrameRef();
    }
    return tryGetDataFrameInternal(index);
}


AnalysisDataFrameRef
AbstractAnalysisData::getDataFrame(int index) const
{
    AnalysisDataFrameRef frame = tryGetDataFrame(index);
    if (!frame.isValid())
    {
        GMX_THROW(APIError("Invalid frame accessed"));
    }
    return frame;
}


bool
AbstractAnalysisData::requestStorage(int nframes)
{
    GMX_RELEASE_ASSERT(nframes >= -1, "Invalid number of frames requested");
    if (nframes == 0)
    {
        return true;
    }
    return requestStorageInternal(nframes);
}


void
AbstractAnalysisData::addModule(AnalysisDataModulePointer module)
{
    impl_->checkModuleProperties(*module);

    if (impl_->bDataStart_)
    {
        GMX_RELEASE_ASSERT(!impl_->bInFrame_,
                           "Cannot add data modules in mid-frame");
        impl_->presentData(this, module.get());
    }
    if (!(module->flags() & AnalysisDataModuleInterface::efAllowMissing))
    {
        impl_->bAllowMissing_ = false;
    }
    impl_->modules_.push_back(module);
}


void
AbstractAnalysisData::addColumnModule(int col, int span,
                                      AnalysisDataModulePointer module)
{
    GMX_RELEASE_ASSERT(col >= 0 && span >= 1,
                       "Invalid columns specified for a column module");
    if (impl_->bDataStart_)
    {
        GMX_THROW(NotImplementedError("Cannot add column modules after data"));
    }

    boost::shared_ptr<AnalysisDataProxy> proxy(
            new AnalysisDataProxy(col, span, this));
    proxy->addModule(module);
    addModule(proxy);
}


void
AbstractAnalysisData::applyModule(AnalysisDataModuleInterface *module)
{
    impl_->checkModuleProperties(*module);
    GMX_RELEASE_ASSERT(impl_->bDataStart_ && !impl_->bInData_,
                       "Data module can only be applied to ready data");

    impl_->presentData(this, module);
}

/*! \cond libapi */
void
AbstractAnalysisData::setDataSetCount(int dataSetCount)
{
    GMX_RELEASE_ASSERT(dataSetCount > 0, "Invalid data column count");
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "Data set count cannot be changed after data has been added");
    impl_->columnCounts_.resize(dataSetCount);
}

void
AbstractAnalysisData::setColumnCount(int dataSet, int columnCount)
{
    GMX_RELEASE_ASSERT(dataSet >= 0 && dataSet < dataSetCount(),
                       "Out of range data set index");
    GMX_RELEASE_ASSERT(columnCount > 0, "Invalid data column count");
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "Data column count cannot be changed after data has been added");
    impl_->columnCounts_[dataSet] = columnCount;
}

void
AbstractAnalysisData::setMultipoint(bool multipoint)
{
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "Data type cannot be changed after data has been added");
    impl_->bMultipoint_ = multipoint;
}


/*! \internal
 * This method is not const because the dataStarted() methods of the attached
 * modules can request storage of the data.
 */
void
AbstractAnalysisData::notifyDataStart()
{
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "notifyDataStart() called more than once");
    for (int d = 0; d < dataSetCount(); ++d)
    {
        GMX_RELEASE_ASSERT(columnCount(d) > 0,
                           "Data column count is not set");
    }
    impl_->bDataStart_ = impl_->bInData_ = true;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        impl_->checkModuleProperties(**i);
        (*i)->dataStarted(this);
    }
}


void
AbstractAnalysisData::notifyFrameStart(const AnalysisDataFrameHeader &header) const
{
    GMX_ASSERT(impl_->bInData_, "notifyDataStart() not called");
    GMX_ASSERT(!impl_->bInFrame_,
               "notifyFrameStart() called while inside a frame");
    GMX_ASSERT(header.index() == impl_->nframes_,
               "Out of order frames");
    impl_->bInFrame_  = true;
    impl_->currIndex_ = header.index();

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        (*i)->frameStarted(header);
    }
}


void
AbstractAnalysisData::notifyPointsAdd(const AnalysisDataPointSetRef &points) const
{
    GMX_ASSERT(impl_->bInData_, "notifyDataStart() not called");
    GMX_ASSERT(impl_->bInFrame_, "notifyFrameStart() not called");
    GMX_ASSERT(points.lastColumn() < columnCount(points.dataSetIndex()),
               "Invalid columns");
    GMX_ASSERT(points.frameIndex() == impl_->currIndex_,
               "Points do not correspond to current frame");
    if (!impl_->bAllowMissing_ && !points.allPresent())
    {
        GMX_THROW(APIError("Missing data not supported by a module"));
    }

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        (*i)->pointsAdded(points);
    }
}


void
AbstractAnalysisData::notifyFrameFinish(const AnalysisDataFrameHeader &header)
{
    GMX_ASSERT(impl_->bInData_, "notifyDataStart() not called");
    GMX_ASSERT(impl_->bInFrame_, "notifyFrameStart() not called");
    GMX_ASSERT(header.index() == impl_->currIndex_,
               "Header does not correspond to current frame");
    impl_->bInFrame_  = false;
    impl_->currIndex_ = -1;

    // Increment the counter before notifications to allow frame access from
    // modules.
    ++impl_->nframes_;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        (*i)->frameFinished(header);
    }
}


void
AbstractAnalysisData::notifyDataFinish() const
{
    GMX_RELEASE_ASSERT(impl_->bInData_, "notifyDataStart() not called");
    GMX_RELEASE_ASSERT(!impl_->bInFrame_,
                       "notifyDataFinish() called while inside a frame");
    impl_->bInData_ = false;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        (*i)->dataFinished();
    }
}
//! \endcond

} // namespace gmx
