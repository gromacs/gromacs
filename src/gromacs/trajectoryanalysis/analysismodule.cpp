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
 * Implements classes in analysismodule.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "analysismodule.h"

#include <utility>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * TrajectoryAnalysisModule::Impl
 */

/*! \internal \brief
 * Private implementation class for TrajectoryAnalysisModule.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisModule::Impl
{
    public:
        //! Container that associates a data set with its name.
        typedef std::map<std::string, AbstractAnalysisData *> DatasetContainer;
        //! Container that associates a AnalysisData object with its name.
        typedef std::map<std::string, AnalysisData *> AnalysisDatasetContainer;

        //! Initializes analysis module data with given name and description.
        Impl(const char *name, const char *description)
            : name_(name), description_(description)
        {
        }

        //! Name of the module.
        std::string                     name_;
        //! Description of the module.
        std::string                     description_;
        //! List of registered data set names.
        std::vector<std::string>        datasetNames_;
        /*! \brief
         * Keeps all registered data sets.
         *
         * This container also includes datasets from \a analysisDatasets_.
         */
        DatasetContainer                datasets_;
        //! Keeps registered AnalysisData objects.
        AnalysisDatasetContainer        analysisDatasets_;
};

/********************************************************************
 * TrajectoryAnalysisModuleData::Impl
 */

/*! \internal \brief
 * Private implementation class for TrajectoryAnalysisModuleData.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisModuleData::Impl
{
    public:
        //! Container that associates a data handle to its AnalysisData object.
        typedef std::map<const AnalysisData *, AnalysisDataHandle>
            HandleContainer;

        //! \copydoc TrajectoryAnalysisModuleData::TrajectoryAnalysisModuleData()
        Impl(TrajectoryAnalysisModule          *module,
             const AnalysisDataParallelOptions &opt,
             const SelectionCollection         &selections);

        //! Checks whether the given AnalysisData has been initialized.
        bool isInitialized(const AnalysisData &data) const;

        //! Keeps a data handle for each AnalysisData object.
        HandleContainer            handles_;
        //! Stores thread-local selections.
        const SelectionCollection &selections_;
};

TrajectoryAnalysisModuleData::Impl::Impl(
        TrajectoryAnalysisModule          *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection         &selections)
    : selections_(selections)
{
    TrajectoryAnalysisModule::Impl::AnalysisDatasetContainer::const_iterator i;
    for (i = module->impl_->analysisDatasets_.begin();
         i != module->impl_->analysisDatasets_.end(); ++i)
    {
        AnalysisDataHandle handle;
        if (isInitialized(*i->second))
        {
            handle = i->second->startData(opt);
        }
        handles_.insert(std::make_pair(i->second, handle));
    }
}

bool TrajectoryAnalysisModuleData::Impl::isInitialized(
        const AnalysisData &data) const
{
    for (int i = 0; i < data.dataSetCount(); ++i)
    {
        if (data.columnCount(i) > 0)
        {
            // If not all of the column counts are set, startData() in the
            // constructor asserts, so that does not need to be checked here.
            return true;
        }
    }
    return false;
}


/********************************************************************
 * TrajectoryAnalysisModuleData
 */

TrajectoryAnalysisModuleData::TrajectoryAnalysisModuleData(
        TrajectoryAnalysisModule          *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection         &selections)
    : impl_(new Impl(module, opt, selections))
{
}


TrajectoryAnalysisModuleData::~TrajectoryAnalysisModuleData()
{
}


void TrajectoryAnalysisModuleData::finishDataHandles()
{
    // FIXME: Call finishData() for all handles even if one throws
    Impl::HandleContainer::iterator i;
    for (i = impl_->handles_.begin(); i != impl_->handles_.end(); ++i)
    {
        if (i->second.isValid())
        {
            i->second.finishData();
        }
    }
    impl_->handles_.clear();
}


AnalysisDataHandle
TrajectoryAnalysisModuleData::dataHandle(const AnalysisData &data)
{
    Impl::HandleContainer::const_iterator i = impl_->handles_.find(&data);
    GMX_RELEASE_ASSERT(i != impl_->handles_.end(),
                       "Data handle requested on unknown dataset");
    return i->second;
}


Selection TrajectoryAnalysisModuleData::parallelSelection(const Selection &selection)
{
    // TODO: Implement properly.
    return selection;
}


SelectionList
TrajectoryAnalysisModuleData::parallelSelections(const SelectionList &selections)
{
    // TODO: Consider an implementation that does not allocate memory every time.
    SelectionList                 newSelections;
    newSelections.reserve(selections.size());
    SelectionList::const_iterator i = selections.begin();
    for (; i != selections.end(); ++i)
    {
        newSelections.push_back(parallelSelection(*i));
    }
    return newSelections;
}


/********************************************************************
 * TrajectoryAnalysisModuleDataBasic
 */

namespace
{

/*! \brief
 * Basic thread-local trajectory analysis data storage class.
 *
 * Most simple tools should only require data handles and selections to be
 * thread-local, so this class implements just that.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisModuleDataBasic : public TrajectoryAnalysisModuleData
{
    public:
        /*! \brief
         * Initializes thread-local storage for data handles and selections.
         *
         * \param[in] module     Analysis module to use for data objects.
         * \param[in] opt        Data parallelization options.
         * \param[in] selections Thread-local selection collection.
         */
        TrajectoryAnalysisModuleDataBasic(TrajectoryAnalysisModule          *module,
                                          const AnalysisDataParallelOptions &opt,
                                          const SelectionCollection         &selections);

        virtual void finish();
};

TrajectoryAnalysisModuleDataBasic::TrajectoryAnalysisModuleDataBasic(
        TrajectoryAnalysisModule          *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection         &selections)
    : TrajectoryAnalysisModuleData(module, opt, selections)
{
}


void
TrajectoryAnalysisModuleDataBasic::finish()
{
    finishDataHandles();
}

}   // namespace


/********************************************************************
 * TrajectoryAnalysisModule
 */

TrajectoryAnalysisModule::TrajectoryAnalysisModule(const char *name,
                                                   const char *description)
    : impl_(new Impl(name, description))
{
}


TrajectoryAnalysisModule::~TrajectoryAnalysisModule()
{
}


void TrajectoryAnalysisModule::optionsFinished(
        Options                    * /*options*/,
        TrajectoryAnalysisSettings * /*settings*/)
{
}


void TrajectoryAnalysisModule::initAfterFirstFrame(
        const TrajectoryAnalysisSettings & /*settings*/,
        const t_trxframe                 & /*fr*/)
{
}


TrajectoryAnalysisModuleDataPointer
TrajectoryAnalysisModule::startFrames(const AnalysisDataParallelOptions &opt,
                                      const SelectionCollection         &selections)
{
    return TrajectoryAnalysisModuleDataPointer(
            new TrajectoryAnalysisModuleDataBasic(this, opt, selections));
}


void TrajectoryAnalysisModule::finishFrames(TrajectoryAnalysisModuleData * /*pdata*/)
{
}


const char *TrajectoryAnalysisModule::name() const
{
    return impl_->name_.c_str();
}


const char *TrajectoryAnalysisModule::description() const
{
    return impl_->description_.c_str();
}


int TrajectoryAnalysisModule::datasetCount() const
{
    return impl_->datasetNames_.size();
}


const std::vector<std::string> &TrajectoryAnalysisModule::datasetNames() const
{
    return impl_->datasetNames_;
}


AbstractAnalysisData &TrajectoryAnalysisModule::datasetFromIndex(int index) const
{
    if (index < 0 || index >= datasetCount())
    {
        GMX_THROW(APIError("Out of range data set index"));
    }
    Impl::DatasetContainer::const_iterator item
        = impl_->datasets_.find(impl_->datasetNames_[index]);
    GMX_RELEASE_ASSERT(item != impl_->datasets_.end(),
                       "Inconsistent data set names");
    return *item->second;
}


AbstractAnalysisData &TrajectoryAnalysisModule::datasetFromName(const char *name) const
{
    Impl::DatasetContainer::const_iterator item = impl_->datasets_.find(name);
    if (item == impl_->datasets_.end())
    {
        GMX_THROW(APIError("Unknown data set name"));
    }
    return *item->second;
}


void TrajectoryAnalysisModule::registerBasicDataset(AbstractAnalysisData *data,
                                                    const char           *name)
{
    GMX_RELEASE_ASSERT(data != NULL, "Attempting to register NULL data");
    // TODO: Strong exception safety should be possible to implement.
    GMX_RELEASE_ASSERT(impl_->datasets_.find(name) == impl_->datasets_.end(),
                       "Duplicate data set name registered");
    impl_->datasets_[name] = data;
    impl_->datasetNames_.push_back(name);
}


void TrajectoryAnalysisModule::registerAnalysisDataset(AnalysisData *data,
                                                       const char   *name)
{
    // TODO: Strong exception safety should be possible to implement.
    registerBasicDataset(data, name);
    impl_->analysisDatasets_[name] = data;
}


void TrajectoryAnalysisModule::finishFrameSerial(int frameIndex)
{
    Impl::AnalysisDatasetContainer::const_iterator data;
    for (data = impl_->analysisDatasets_.begin();
         data != impl_->analysisDatasets_.end();
         ++data)
    {
        data->second->finishFrameSerial(frameIndex);
    }
}

} // namespace gmx
