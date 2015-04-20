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
 * Declares gmx::TrajectoryAnalysisRunnerCommon.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_RUNNERCOMMON_H
#define GMX_TRAJECTORYANALYSIS_RUNNERCOMMON_H

#include "gromacs/utility/classhelpers.h"

struct t_trxframe;

namespace gmx
{

class Options;
class SelectionCollection;
class TopologyInformation;
class TrajectoryAnalysisSettings;

/*! \internal
 * \brief
 * Implements common trajectory analysis runner functionality.
 *
 * As there is currently only one runner (TrajectoryAnalysisCommandLineRunner),
 * the division of responsibilities is not yet very clear.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisRunnerCommon
{
    public:
        /*! \brief
         * Initializes a new runner helper.
         *
         * \param    settings  Settings object to use.
         */
        explicit TrajectoryAnalysisRunnerCommon(TrajectoryAnalysisSettings *settings);
        ~TrajectoryAnalysisRunnerCommon();

        /*! \brief
         * Initializes common options for trajectory analysis.
         *
         * \param[in,out] options  Options object to add the options to.
         */
        void initOptions(Options *options);
        //! Scales time option values according to the time unit set.
        void scaleTimeOptions(Options *options);
        /*! \brief
         * Processes common option values after they have been parsed.
         *
         * \param[in,out] options Options object in which options are stored.
         */
        void optionsFinished(Options *options);
        //! Initialize index groups for selections.
        void initIndexGroups(SelectionCollection *selections,
                             bool                 bUseDefaults);
        //! Free memory allocated for index groups.
        void doneIndexGroups(SelectionCollection *selections);
        //! Load topology information if provided and/or required.
        void initTopology(SelectionCollection *selections);
        /*! \brief
         * Reads the first frame from the trajectory.
         *
         * After this call, frame() returns the first frame.
         */
        void initFirstFrame();
        /*! \brief
         * Reads the next frame from the trajectory.
         *
         * \returns false if there were no more frames.
         *
         * After this call, frame() returns the newly loaded frame.
         */
        bool readNextFrame();
        /*! \brief
         * Performs common initialization for the currently loaded frame.
         *
         * Currently, makes molecules whole if requested.
         */
        void initFrame();

        //! Returns true if input data comes from a trajectory.
        bool hasTrajectory() const;
        //! Returns the topology information object.
        const TopologyInformation &topologyInformation() const;
        //! Returns the currently loaded frame.
        t_trxframe                &frame() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
