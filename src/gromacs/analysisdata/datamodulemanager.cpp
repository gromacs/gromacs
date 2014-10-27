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
/*! \internal \file
 * \brief
 * Implements gmx::AnalysisDataModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#include "gmxpre.h"

#include "datamodulemanager.h"

#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataModuleManager::Impl
 */

/*! \internal \brief
 * Private implementation class for AnalysisDataModuleManager.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataModuleManager::Impl
{
    public:
        //! Stores information about an attached module.
        struct ModuleInfo
        {
            //! Initializes the module information.
            explicit ModuleInfo(AnalysisDataModulePointer module)
                : module(module), bParallel(false)
            {
            }

            //! Pointer to the actual module.
            AnalysisDataModulePointer   module;
            //! Whether the module supports parallel processing.
            bool                        bParallel;
        };

        //! Shorthand for list of modules added to the data.
        typedef std::vector<ModuleInfo> ModuleList;

        //! Describes the current state of the notification methods.
        enum State
        {
            eNotStarted, //!< Initial state (nothing called).
            eInData,     //!< notifyDataStart() called, no frame in progress.
            eInFrame,    //!< notifyFrameStart() called, but notifyFrameFinish() not.
            eFinished    //!< notifyDataFinish() called.
        };

        Impl();

        /*! \brief
         * Checks whether a module is compatible with a given data property.
         *
         * \param[in] module   Module to check.
         * \param[in] property Property to check.
         * \param[in] bSet     Value of the property to check against.
         * \throws    APIError if \p module is not compatible with the data.
         */
        void checkModuleProperty(const AnalysisDataModuleInterface &module,
                                 DataProperty property, bool bSet) const;
        /*! \brief
         * Checks whether a module is compatible with the data properties.
         *
         * \param[in] module Module to check.
         * \throws    APIError if \p module is not compatible with the data.
         *
         * Does not currently check the actual data (e.g., missing values), but
         * only the dimensionality and other preset properties of the data.
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

        //! List of modules added to the data.
        ModuleList              modules_;
        //! Properties of the owning data for module checking.
        bool                    bDataProperty_[eDataPropertyNR];
        //! true if all modules support missing data.
        bool                    bAllowMissing_;
        //! true if there are modules that do not support parallel processing.
        bool                    bSerialModules_;
        //! true if there are modules that support parallel processing.
        bool                    bParallelModules_;

        /*! \brief
         * Current state of the notification methods.
         *
         * This is used together with \a currIndex_ for sanity checks on the
         * input data; invalid call sequences trigger asserts.
         * The state of these variables does not otherwise affect the behavior
         * of this class; this is the reason they can be changed in const
         * methods.
         */
        //! Whether notifyDataStart() has been called.
        mutable State           state_;
        //! Index of currently active frame or the next frame if not in frame.
        mutable int             currIndex_;
};

AnalysisDataModuleManager::Impl::Impl()
    : bAllowMissing_(true), bSerialModules_(false), bParallelModules_(false),
      state_(eNotStarted), currIndex_(0)
{
    // This must be in sync with how AbstractAnalysisData is actually
    // initialized.
    for (int i = 0; i < eDataPropertyNR; ++i)
    {
        bDataProperty_[i] = false;
    }
}

void
AnalysisDataModuleManager::Impl::checkModuleProperty(
        const AnalysisDataModuleInterface &module,
        DataProperty property, bool bSet) const
{
    bool      bOk   = true;
    const int flags = module.flags();
    switch (property)
    {
        case eMultipleDataSets:
            if (bSet && !(flags & AnalysisDataModuleInterface::efAllowMultipleDataSets))
            {
                bOk = false;
            }
            break;
        case eMultipleColumns:
            if (bSet && !(flags & AnalysisDataModuleInterface::efAllowMulticolumn))
            {
                bOk = false;
            }
            break;
        case eMultipoint:
            if ((bSet && !(flags & AnalysisDataModuleInterface::efAllowMultipoint))
                || (!bSet && (flags & AnalysisDataModuleInterface::efOnlyMultipoint)))
            {
                bOk = false;
            }
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Invalid data property enumeration");
    }
    if (!bOk)
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }
}

void
AnalysisDataModuleManager::Impl::checkModuleProperties(
        const AnalysisDataModuleInterface &module) const
{
    for (int i = 0; i < eDataPropertyNR; ++i)
    {
        checkModuleProperty(module, static_cast<DataProperty>(i), bDataProperty_[i]);
    }
}

void
AnalysisDataModuleManager::Impl::presentData(AbstractAnalysisData        *data,
                                             AnalysisDataModuleInterface *module)
{
    if (state_ == eNotStarted)
    {
        return;
    }
    GMX_RELEASE_ASSERT(state_ != eInFrame,
                       "Cannot apply a modules in mid-frame");
    module->dataStarted(data);
    const bool bCheckMissing = bAllowMissing_
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
        module->frameFinishedSerial(frame.header().index());
    }
    if (state_ == eFinished)
    {
        module->dataFinished();
    }
}

/********************************************************************
 * AnalysisDataModuleManager
 */

AnalysisDataModuleManager::AnalysisDataModuleManager()
    : impl_(new Impl())
{
}

AnalysisDataModuleManager::~AnalysisDataModuleManager()
{
}

void
AnalysisDataModuleManager::dataPropertyAboutToChange(DataProperty property, bool bSet)
{
    GMX_RELEASE_ASSERT(impl_->state_ == Impl::eNotStarted,
                       "Cannot change data properties after data has been started");
    if (impl_->bDataProperty_[property] != bSet)
    {
        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            impl_->checkModuleProperty(*i->module, property, bSet);
        }
        impl_->bDataProperty_[property] = bSet;
    }
}

void
AnalysisDataModuleManager::addModule(AbstractAnalysisData      *data,
                                     AnalysisDataModulePointer  module)
{
    impl_->checkModuleProperties(*module);
    // TODO: Ensure that the system does not end up in an inconsistent state by
    // adding a module in mid-data during parallel processing (probably best to
    // prevent alltogether).
    GMX_RELEASE_ASSERT(impl_->state_ != Impl::eInFrame,
                       "Cannot add a data module in mid-frame");
    impl_->presentData(data, module.get());

    if (!(module->flags() & AnalysisDataModuleInterface::efAllowMissing))
    {
        impl_->bAllowMissing_ = false;
    }
    impl_->modules_.push_back(Impl::ModuleInfo(module));
}

void
AnalysisDataModuleManager::applyModule(AbstractAnalysisData        *data,
                                       AnalysisDataModuleInterface *module)
{
    impl_->checkModuleProperties(*module);
    GMX_RELEASE_ASSERT(impl_->state_ == Impl::eFinished,
                       "Data module can only be applied to ready data");
    impl_->presentData(data, module);
}


bool
AnalysisDataModuleManager::hasSerialModules() const
{
    GMX_ASSERT(impl_->state_ != Impl::eNotStarted,
               "Module state not accessible before data is started");
    return impl_->bSerialModules_;
}


void
AnalysisDataModuleManager::notifyDataStart(AbstractAnalysisData *data)
{
    GMX_RELEASE_ASSERT(impl_->state_ == Impl::eNotStarted,
                       "notifyDataStart() called more than once");
    for (int d = 0; d < data->dataSetCount(); ++d)
    {
        GMX_RELEASE_ASSERT(data->columnCount(d) > 0,
                           "Data column count is not set");
    }
    impl_->state_            = Impl::eInData;
    impl_->bSerialModules_   = !impl_->modules_.empty();
    impl_->bParallelModules_ = false;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        // This should not fail, since addModule() and
        // dataPropertyAboutToChange() already do the checks, but kept here to
        // catch potential bugs (perhaps it would be best to assert on failure).
        impl_->checkModuleProperties(*i->module);
        i->module->dataStarted(data);
    }
}


void
AnalysisDataModuleManager::notifyParallelDataStart(
        AbstractAnalysisData              *data,
        const AnalysisDataParallelOptions &options)
{
    GMX_RELEASE_ASSERT(impl_->state_ == Impl::eNotStarted,
                       "notifyDataStart() called more than once");
    for (int d = 0; d < data->dataSetCount(); ++d)
    {
        GMX_RELEASE_ASSERT(data->columnCount(d) > 0,
                           "Data column count is not set");
    }
    impl_->state_            = Impl::eInData;
    impl_->bSerialModules_   = false;
    impl_->bParallelModules_ = false;

    Impl::ModuleList::iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        // This should not fail, since addModule() and
        // dataPropertyAboutToChange() already do the checks, but kept here to
        // catch potential bugs (perhaps it would be best to assert on failure).
        impl_->checkModuleProperties(*i->module);
        i->bParallel = i->module->parallelDataStarted(data, options);
        if (i->bParallel)
        {
            impl_->bParallelModules_ = true;
        }
        else
        {
            impl_->bSerialModules_ = true;
        }
    }
}


void
AnalysisDataModuleManager::notifyFrameStart(const AnalysisDataFrameHeader &header) const
{
    GMX_ASSERT(impl_->state_ == Impl::eInData, "Invalid call sequence");
    GMX_ASSERT(header.index() == impl_->currIndex_, "Out of order frames");
    impl_->state_     = Impl::eInFrame;

    if (impl_->bSerialModules_)
    {
        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (!i->bParallel)
            {
                i->module->frameStarted(header);
            }
        }
    }
}

void
AnalysisDataModuleManager::notifyParallelFrameStart(
        const AnalysisDataFrameHeader &header) const
{
    if (impl_->bParallelModules_)
    {
        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (i->bParallel)
            {
                i->module->frameStarted(header);
            }
        }
    }
}


void
AnalysisDataModuleManager::notifyPointsAdd(const AnalysisDataPointSetRef &points) const
{
    GMX_ASSERT(impl_->state_ == Impl::eInFrame, "notifyFrameStart() not called");
    // TODO: Add checks for column spans (requires passing the information
    // about the column counts from somewhere).
    //GMX_ASSERT(points.lastColumn() < columnCount(points.dataSetIndex()),
    //           "Invalid columns");
    GMX_ASSERT(points.frameIndex() == impl_->currIndex_,
               "Points do not correspond to current frame");
    if (impl_->bSerialModules_)
    {
        if (!impl_->bAllowMissing_ && !points.allPresent())
        {
            GMX_THROW(APIError("Missing data not supported by a module"));
        }

        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (!i->bParallel)
            {
                i->module->pointsAdded(points);
            }
        }
    }
}


void
AnalysisDataModuleManager::notifyParallelPointsAdd(
        const AnalysisDataPointSetRef &points) const
{
    // TODO: Add checks for column spans (requires passing the information
    // about the column counts from somewhere).
    //GMX_ASSERT(points.lastColumn() < columnCount(points.dataSetIndex()),
    //           "Invalid columns");
    if (impl_->bParallelModules_)
    {
        if (!impl_->bAllowMissing_ && !points.allPresent())
        {
            GMX_THROW(APIError("Missing data not supported by a module"));
        }

        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (i->bParallel)
            {
                i->module->pointsAdded(points);
            }
        }
    }
}


void
AnalysisDataModuleManager::notifyFrameFinish(const AnalysisDataFrameHeader &header) const
{
    GMX_ASSERT(impl_->state_ == Impl::eInFrame, "notifyFrameStart() not called");
    GMX_ASSERT(header.index() == impl_->currIndex_,
               "Header does not correspond to current frame");
    // TODO: Add a check for the frame count in the source data including this
    // frame.
    impl_->state_ = Impl::eInData;
    ++impl_->currIndex_;

    if (impl_->bSerialModules_)
    {
        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (!i->bParallel)
            {
                i->module->frameFinished(header);
            }
        }
    }
    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        i->module->frameFinishedSerial(header.index());
    }
}


void
AnalysisDataModuleManager::notifyParallelFrameFinish(
        const AnalysisDataFrameHeader &header) const
{
    if (impl_->bParallelModules_)
    {
        Impl::ModuleList::const_iterator i;
        for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
        {
            if (i->bParallel)
            {
                i->module->frameFinished(header);
            }
        }
    }
}


void
AnalysisDataModuleManager::notifyDataFinish() const
{
    GMX_RELEASE_ASSERT(impl_->state_ == Impl::eInData, "Invalid call sequence");
    impl_->state_ = Impl::eFinished;

    Impl::ModuleList::const_iterator i;
    for (i = impl_->modules_.begin(); i != impl_->modules_.end(); ++i)
    {
        i->module->dataFinished();
    }
}

} // namespace gmx
