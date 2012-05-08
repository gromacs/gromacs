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
 * Declares trajectory analysis module for angle calculations.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_ANGLE_H
#define GMX_TRAJECTORYANALYSIS_MODULES_ANGLE_H

#include <string>

#include "../analysismodule.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"

namespace gmx
{

class SelectionOptionInfo;

namespace analysismodules
{

class Angle : public TrajectoryAnalysisModule
{
    public:
        static const char name[];
        static const char shortDescription[];

        Angle();
        virtual ~Angle();

        virtual Options &initOptions(TrajectoryAnalysisSettings *settings);
        virtual void initOptionsDone(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        void checkSelections(const SelectionList &sel1,
                             const SelectionList &sel2) const;

        Options                 _options;

        SelectionList           _sel1;
        SelectionList           _sel2;
        SelectionOptionInfo    *_sel1info;
        SelectionOptionInfo    *_sel2info;
        std::string             _fnAngle;
        std::string             _fnDump;

        std::string             _g1type;
        std::string             _g2type;
        bool                    _bSplit1;
        bool                    _bSplit2;
        bool                    _bMulti;
        bool                    _bAll;
        bool                    _bDumpDist;

        AnalysisData            _data;
        int                     _natoms1;
        int                     _natoms2;
        rvec                   *_vt0;

        // Copy and assign disallowed by base.
};

} // namespace analysismodules

} // namespace gmxana

#endif
