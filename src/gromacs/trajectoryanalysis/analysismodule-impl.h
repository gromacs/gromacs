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
 * Declares private implementation classes for gmx::TrajectoryAnalysisModule
 * and gmx::TrajectoryAnalysisModuleData.
 *
 * \ingroup module_trajectoryanalysis
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_IMPL_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_IMPL_H

#include <map>
#include <string>
#include <vector>

#include "analysismodule.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;

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
        ~Impl();

        //! Keeps a data handle for each AnalysisData object.
        HandleContainer         _handles;
        //! Stores thread-local selections.
        const SelectionCollection &_selections;
};

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

} // namespace gmx

#endif
