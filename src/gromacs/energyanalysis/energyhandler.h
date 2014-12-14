/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Declares gmx::EnergyHandler
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_HANDLER_H
#define GMX_ENERGYANALYSIS_HANDLER_H

#include <string>
#include <vector>

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/options/timeunitmanager.h"
#include "analysismodule.h"

namespace gmx
{
/*! \brief
 * Class doing the actual reading of an energy file, and passing the
 * information to different energy analysis tools.
 */
class EnergyHandler
{
    private:
        //! The energy files
        std::vector<std::string>                 fnEnergy_;
        //! Start time of the analysis
        double                                   t0_;
        //! End time of the analysis
        double                                   t1_;
        //! Skipping time of the analysis
        double                                   tDelta_;
        //! Do we want to view the output?
        bool                                     bView_;
        //! Do we want verbose output?
        bool                                     bVerbose_;
        //! Time management: even that is solved by GROMACS
        TimeUnitManager                          timeUnitManager_;
        //! Output environment for xvg writing etcetera
        output_env_t                             oenv_;
        //! Module that does all the work
        EnergyAnalysisModulePointer              module_;
        //! Check whether time is within range
        int checkTime(double t);
    public:
        //! Constructor
        EnergyHandler(EnergyAnalysisModulePointer module);

        //! Destructor
        ~EnergyHandler() {}

        /*! \brief
         * Prepare the reading and processing options.
         *
         * Collects additional options from all registered analysis
         * tools, and then processes all of those.
         */
        void initOptions(Options *options);

        void optionsFinished(Options *options);

        /*! \brief
         * Read the files and call the tools to analyze them.
         *
         * The files are read
         * sequentially and tools have to be able to deal with this.
         * \return 0 if everything went smoothly
         */
        int run();
};

}
#endif
