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
 * Declares private implementation class for gmx::TrajectoryAnalysisSettings.
 *
 * \ingroup module_trajectoryanalysis
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_IMPL_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISSETTINGS_IMPL_H

#include "analysissettings.h"

#include "../analysisdata/modules/plot.h"
#include "../options/timeunitmanager.h"

namespace gmx
{

/*! \internal \brief
 * Private implementation class for TrajectoryAnalysisSettings.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisSettings::Impl
{
    public:
        //! Initializes the default values for the settings object.
        Impl() : flags(0), frflags(0), bRmPBC(true), bPBC(true) {}

        //! Global time unit setting for the analysis module.
        TimeUnitManager          timeUnitManager;
        //! Global plotting settings for the analysis module.
        AnalysisDataPlotSettings plotSettings;
        //! Flags for the analysis module.
        unsigned long            flags;
        //! Frame reading flags for the analysis module.
        int                      frflags;

        //! Whether to make molecules whole for each frame.
        bool                 bRmPBC;
        //! Whether to pass PBC information to the analysis module.
        bool                 bPBC;
};

} // namespace gmx

#endif
