/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012, by the GROMACS development team, led by
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
