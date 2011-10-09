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

class TrajectoryAnalysisModuleData::Impl
{
    public:
        typedef std::map<std::string, AnalysisDataHandle *> HandleContainer;

        Impl(TrajectoryAnalysisModule *module,
             /*AnalysisDataParallelOptions*/ void* opt,
             const SelectionCollection &selections);
        ~Impl();

        void finishHandles();

        HandleContainer         _handles;
        const SelectionCollection &_selections;
};

class TrajectoryAnalysisModule::Impl
{
    public:
        typedef std::map<std::string, AbstractAnalysisData *> DatasetContainer;
        typedef std::map<std::string, AnalysisData *> AnalysisDatasetContainer;

        std::vector<std::string>        _datasetNames;
        DatasetContainer                _datasets;
        AnalysisDatasetContainer        _analysisDatasets;
};

/*! \internal \brief
 * Basic thread-local trajectory analysis data storage class.
 *
 * Most simple tools should only require data handles and selections to be
 * thread-local, so this class implements just that.
 */
class TrajectoryAnalysisModuleDataBasic : public TrajectoryAnalysisModuleData
{
    public:
        TrajectoryAnalysisModuleDataBasic(TrajectoryAnalysisModule *module,
                                          /*AnalysisDataParallelOptions*/ void* opt,
                                          const SelectionCollection &selections);

        virtual void finish();
};

} // namespace gmx

#endif
