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
 * Declares gmx::TrajectoryAnalysisRunnerCommon.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_RUNNERCOMMON_H
#define GMX_TRAJECTORYANALYSIS_RUNNERCOMMON_H

#include "../utility/common.h"

namespace gmx
{

class Options;
class SelectionCollection;
class TopologyInformation;
class TrajectoryAnalysisSettings;

/*! \internal \brief
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
         * Flags that define what kind of help should be printed.
         */
        enum HelpFlag
        {
            efHelpShowOptions           = 1<<0, //!< Show options.
            efHelpShowHidden            = 1<<1, //!< Show hidden options.
            efHelpShowDescriptions      = 1<<2  //!< Show detailed description.
        };
        //! Combination of \ref HelpFlag values.
        typedef unsigned long HelpFlags;

        /*! \brief
         * Initializes a new runner helper.
         *
         * \param    settings  Settings object to use.
         */
        explicit TrajectoryAnalysisRunnerCommon(TrajectoryAnalysisSettings *settings);
        ~TrajectoryAnalysisRunnerCommon();

        //! Initializes common options for trajectory analysis.
        Options &initOptions();
        //! Scales time option values according to the time unit set.
        void scaleTimeOptions(Options *options);
        /*! \brief
         * Processes common option values after they have been parsed.
         *
         * \returns false if the tool should exit after printing help.
         */
        bool initOptionsDone();
        //! Initialize index groups for selections.
        void initIndexGroups(SelectionCollection *selections);
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

        //! Returns flags for help printing.
        HelpFlags helpFlags() const;
        //! Returns true if input data comes from a trajectory.
        bool hasTrajectory() const;
        //! Returns the topology information object.
        const TopologyInformation &topologyInformation() const;
        //! Returns the currently loaded frame.
        t_trxframe &frame() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
