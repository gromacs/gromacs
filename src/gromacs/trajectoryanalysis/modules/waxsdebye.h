/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Declares trajectory analysis module for distance calculations.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_WAXSDEBYE_H
#define GMX_TRAJECTORYANALYSIS_MODULES_WAXSDEBYE_H

#include <string>
#include "gromacs/legacyheaders/gmx_random.h"
#include "../analysismodule.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/selection/selection.h"
#include "gromacs/waxsdebye/waxs_debye_force.h"

namespace gmx
{

namespace analysismodules
{

/*! \brief
 * Class used to compute waxs scatering in a simulations box.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 * Does not implement any new functionality.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class WaxsDebye : public TrajectoryAnalysisModule
{
    public:
        //! Name of the tool
        static const char name[];

        //! One line description
        static const char shortDescription[];

        //! Constructor
        WaxsDebye();

        //! Destructor
        virtual ~WaxsDebye();

        //! Set the options and setting
        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);

        //! First routine called by the analysis frame work
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        //! Call for each frame of the trajectory
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        //! Last routine called by the analysis frame work
        virtual void finishAnalysis(int nframes);

        //! Routine to write output, that is additional over the built-in
        virtual void writeOutput();

    private:
        std::string                       fnSfactor_, fnSqref_, fnSqdiff_, fnSqcalc_, fnEner_, fnAlpha_;
        Selection                         sel_;
        AnalysisData                      data_;
        AnalysisDataAverageModulePointer  adata_;
        WaxsDebyeForce                   *wdf_;

        // Copy and assign disallowed by base.
};

} // namespace analysismodules

} // namespace gmx

#endif
