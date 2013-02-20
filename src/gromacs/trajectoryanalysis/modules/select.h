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
 * Declares trajectory analysis module for basic selection information.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_SELECT_H
#define GMX_TRAJECTORYANALYSIS_MODULES_SELECT_H

#include <string>
#include <vector>

#include "../analysismodule.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/selection/selection.h"

namespace gmx
{

class SelectionOptionInfo;

namespace analysismodules
{

class Select : public TrajectoryAnalysisModule
{
    public:
        static const char name[];
        static const char shortDescription[];

        Select();
        virtual ~Select();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void optionsFinished(Options                    *options,
                                     TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        SelectionList                    sel_;
        SelectionOptionInfo             *selOpt_;

        std::string                      fnSize_;
        std::string                      fnFrac_;
        std::string                      fnIndex_;
        std::string                      fnNdx_;
        std::string                      fnMask_;
        std::string                      fnOccupancy_;
        std::string                      fnPDB_;
        bool                             bDump_;
        bool                             bTotNorm_;
        bool                             bFracNorm_;
        bool                             bResInd_;
        std::string                      resNumberType_;
        std::string                      pdbAtoms_;

        const TopologyInformation       *top_;
        std::vector<int>                 totsize_;
        AnalysisData                     sdata_;
        AnalysisData                     cdata_;
        AnalysisData                     idata_;
        AnalysisData                     mdata_;
        AnalysisDataAverageModulePointer occupancyModule_;
};

} // namespace analysismodules

} // namespace gmx

#endif
