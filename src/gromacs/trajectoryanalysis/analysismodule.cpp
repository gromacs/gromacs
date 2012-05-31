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
 * Implements classes in analysismodule.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gromacs/trajectoryanalysis/analysismodule.h"

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

        //! List of registered data set names.
        std::vector<std::string>        _datasetNames;
        /*! \brief
         * Keeps all registered data sets.
         *
         * This container also includes datasets from \a _analysisDatasets.
         */
        DatasetContainer                _datasets;
        //! Keeps registered AnalysisData objects.
        AnalysisDatasetContainer        _analysisDatasets;
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
        Impl(TrajectoryAnalysisModule *module,
             const AnalysisDataParallelOptions &opt,
             const SelectionCollection &selections);

        //! Keeps a data handle for each AnalysisData object.
        HandleContainer         _handles;
        //! Stores thread-local selections.
        const SelectionCollection &_selections;
};

TrajectoryAnalysisModuleData::Impl::Impl(
        TrajectoryAnalysisModule *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection &selections)
    : _selections(selections)
{
    TrajectoryAnalysisModule::Impl::AnalysisDatasetContainer::const_iterator i;
    for (i = module->_impl->_analysisDatasets.begin();
         i != module->_impl->_analysisDatasets.end(); ++i)
    {
        _handles.insert(std::make_pair(i->second, i->second->startData(opt)));
    }
}


/********************************************************************
 * TrajectoryAnalysisModuleData
 */

TrajectoryAnalysisModuleData::TrajectoryAnalysisModuleData(
        TrajectoryAnalysisModule *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection &selections)
    : _impl(new Impl(module, opt, selections))
{
}


TrajectoryAnalysisModuleData::~TrajectoryAnalysisModuleData()
{
}


void TrajectoryAnalysisModuleData::finishDataHandles()
{
    // FIXME: Call finishData() for all handles even if one throws
    Impl::HandleContainer::iterator i;
    for (i = _impl->_handles.begin(); i != _impl->_handles.end(); ++i)
    {
        i->second.finishData();
    }
    _impl->_handles.clear();
}


AnalysisDataHandle
TrajectoryAnalysisModuleData::dataHandle(const AnalysisData &data)
{
    Impl::HandleContainer::const_iterator i = _impl->_handles.find(&data);
    GMX_RELEASE_ASSERT(i != _impl->_handles.end(),
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
    SelectionList newSelections;
    newSelections.reserve(selections.size());
    SelectionList::const_iterator i = selections.begin();
    for ( ; i != selections.end(); ++i)
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

/*! \internal \brief
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
        TrajectoryAnalysisModuleDataBasic(TrajectoryAnalysisModule *module,
                                          const AnalysisDataParallelOptions &opt,
                                          const SelectionCollection &selections);

        virtual void finish();
};

TrajectoryAnalysisModuleDataBasic::TrajectoryAnalysisModuleDataBasic(
        TrajectoryAnalysisModule *module,
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection &selections)
    : TrajectoryAnalysisModuleData(module, opt, selections)
{
}


void
TrajectoryAnalysisModuleDataBasic::finish()
{
    finishDataHandles();
}

} // namespace


/********************************************************************
 * TrajectoryAnalysisModule
 */

TrajectoryAnalysisModule::TrajectoryAnalysisModule()
    : _impl(new Impl)
{
}


TrajectoryAnalysisModule::~TrajectoryAnalysisModule()
{
}


void TrajectoryAnalysisModule::initOptionsDone(TrajectoryAnalysisSettings * /*settings*/)
{
}


void TrajectoryAnalysisModule::initAfterFirstFrame(const t_trxframe &/*fr*/)
{
}


TrajectoryAnalysisModuleDataPointer
TrajectoryAnalysisModule::startFrames(const AnalysisDataParallelOptions &opt,
                                      const SelectionCollection &selections)
{
    return TrajectoryAnalysisModuleDataPointer(
            new TrajectoryAnalysisModuleDataBasic(this, opt, selections));
}


void TrajectoryAnalysisModule::finishFrames(TrajectoryAnalysisModuleData * /*pdata*/)
{
}


int TrajectoryAnalysisModule::datasetCount() const
{
    return _impl->_datasetNames.size();
}


const std::vector<std::string> &TrajectoryAnalysisModule::datasetNames() const
{
    return _impl->_datasetNames;
}


AbstractAnalysisData &TrajectoryAnalysisModule::datasetFromIndex(int index) const
{
    if (index < 0 || index >= datasetCount())
    {
        GMX_THROW(APIError("Out of range data set index"));
    }
    Impl::DatasetContainer::const_iterator item
        = _impl->_datasets.find(_impl->_datasetNames[index]);
    GMX_RELEASE_ASSERT(item != _impl->_datasets.end(),
                       "Inconsistent data set names");
    return *item->second;
}


AbstractAnalysisData &TrajectoryAnalysisModule::datasetFromName(const char *name) const
{
    Impl::DatasetContainer::const_iterator item = _impl->_datasets.find(name);
    if (item == _impl->_datasets.end())
    {
        GMX_THROW(APIError("Unknown data set name"));
    }
    return *item->second;
}


void TrajectoryAnalysisModule::registerBasicDataset(AbstractAnalysisData *data,
                                                    const char *name)
{
    // TODO: Strong exception safety should be possible to implement.
    GMX_RELEASE_ASSERT(_impl->_datasets.find(name) == _impl->_datasets.end(),
                       "Duplicate data set name registered");
    _impl->_datasets[name] = data;
    _impl->_datasetNames.push_back(name);
}


void TrajectoryAnalysisModule::registerAnalysisDataset(AnalysisData *data,
                                                       const char *name)
{
    // TODO: Strong exception safety should be possible to implement.
    registerBasicDataset(data, name);
    _impl->_analysisDatasets[name] = data;
}

} // namespace gmx
