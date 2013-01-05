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
 * Implements gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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

        /*! \brief
         * Present data already added to the data object to a module.
         *
         * \param[in] data   Data object to read data from.
         * \param[in] module Module to present the data to.
         * \throws    APIError if \p module is not compatible with the data
         *      object.
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
    : bAllowMissing_(true), bDataStart_(false), bInData_(false), bInFrame_(false),
      currIndex_(-1), nframes_(0)
{
}

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
        module->pointsAdded(frame.points());
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
    : impl_(new Impl()), columnCount_(0), bMultiPoint_(false)
{
}
//! \endcond

AbstractAnalysisData::~AbstractAnalysisData()
{
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
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }

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
    GMX_RELEASE_ASSERT(col >= 0 && span >= 1 && col + span <= columnCount_,
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
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }
    GMX_RELEASE_ASSERT(impl_->bDataStart_ && !impl_->bInData_,
                       "Data module can only be applied to ready data");

    impl_->presentData(this, module);
}

/*! \cond libapi */
void
AbstractAnalysisData::setColumnCount(int columnCount)
{
    GMX_RELEASE_ASSERT(columnCount > 0, "Invalid data column count");
    GMX_RELEASE_ASSERT(columnCount_ == 0 || impl_->modules_.empty(),
                       "Data column count cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "Data column count cannot be changed after data has been added");
    columnCount_ = columnCount;
}


void
AbstractAnalysisData::setMultipoint(bool multipoint)
{
    GMX_RELEASE_ASSERT(impl_->modules_.empty(),
                       "Data type cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!impl_->bDataStart_,
                       "Data type cannot be changed after data has been added");
    bMultiPoint_ = multipoint;
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
    GMX_RELEASE_ASSERT(columnCount_ > 0, "Data column count is not set");
    impl_->bDataStart_ = impl_->bInData_ = true;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        if (columnCount_ > 1 && !((*i)->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        {
            GMX_THROW(APIError("Data module not compatible with data object properties"));
        }
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
    GMX_ASSERT(points.lastColumn() < columnCount(), "Invalid columns");
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
