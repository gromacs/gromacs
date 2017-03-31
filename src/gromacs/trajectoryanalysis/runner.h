/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

/*! \libinternal \file
 * \brief
 * Declares gmx::trajectoryanalysis::Runner.
 *
 * \ingroup module_trajectoryanalysis
 */

#ifndef GMX_TRAJECTORYANALYSIS_RUNNER_H
#define GMX_TRAJECTORYANALYSIS_RUNNER_H

#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

#include "runnercommon.h"


namespace gmx
{

class IOptionsContainer;

namespace trajectoryanalysis
{

/*! \brief Generic runner for modules in Trajectory Analysis Framework
 *
 * The protocol described in the Trajectory Analysis Framework for a runner to
 * drive analysis tools is not implemented generically and publicly for runner
 * classes to use.
 *
 * The ultimate goal of this class is to provide a base behavior class for
 * a runner implementing the protocol described in the Trajectory Analysis
 * Framework to drive an analysis tool chain. The command-line runner and
 * command-line analysis tool template can use the common behavior in this
 * class. It is intended to incorporate behavior of the RunnerModule class (defined in an
 * anonymous namespace in gromacs/trajectoryanalysis/cmdlinerunner.cpp)
 * into a possible replacement for TrajectoryAnalysisRunnerCommon, which will
 * initially simply be used.
 *
 * First revision handles only a single module and no parallelism.
 *
 * Implements the Trajectory Analysis Framework protocol to run analysis modules
 * derived from gmx::TrajectoryAnalysisModule, registered with addModule().
 *
 * Usage: The runner is constructed, a module bound, the runner initialized with
 * options, and then the user may either call next() to step one input frame or
 * run() to process all remaining input. The analysis is finalized when the
 * runner is destroyed, releasing the module(s).
 *
 * Runner objects are not assignable or copyable.
 *
 * Example:
 * ```
 *     // Assume Module is a class derived from TrajectoryAnalysisModule
 *     auto mymodule = std::make_shared<Module>(Module());
 *     auto runner = Runner();
 *     runner.addModule(mymodule);
 *     runner.registerOptions(options)
 *     runner.initialize(options);
 *     bool end_of_data = runner.next();
 *     runner.run();
 *     assert(runner.next() == false);
 *     del runner;
 *     do_something(mymodule);
 *     del mymodule;
 * ```
 * \libinternal \ingroup module_trajectoryanalysis
 */
class Runner
{
    public:
        /// Default constructor
        Runner();

        /// Clean up
        virtual ~Runner();

        /// Disallow copy construction
        Runner(const Runner &runner)            = delete;
        /// Disallow copy assignment
        const Runner &operator=(const Runner &runner) = delete;

        /// Populate an options object provided by the caller.
        //! \param options created and owned by the caller
        void registerOptions(gmx::Options &options);

        /// Finalize options and read first frame.
        //! \param options created and owned by the caller
        void initialize(gmx::Options &options);

        /*! \brief Registers a TrajectoryAnalysisModule with the runner.
         *
         * \param module Currently can only handle a single runner.
         * \return another shared pointer to the module added.
         * Could be set to null if unsuccessful. This exists in case we want to
         * allow an overload that takes ownership of a module that is not already
         * managed.
         */
        gmx::TrajectoryAnalysisModuleSharedPointer addModule(gmx::TrajectoryAnalysisModuleSharedPointer module);

        /*! \brief Advance one frame.
         *
         * \return false if there are no more frames.
         */
        bool next();

        /// Process all remaining available frames
        int run();

    private:
        /// Handle to bound module
        gmx::TrajectoryAnalysisModuleSharedPointer module_;

        TrajectoryAnalysisSettings                 settings_;
        TrajectoryAnalysisRunnerCommon             common_;
        SelectionCollection selections_;

        // data object returned by module_->startFrames()
        std::unique_ptr<TrajectoryAnalysisModuleData> pdata_;

        /// Number of frames processed
        unsigned int nframes_;
        /// True if modules have been initialized and first frame read.
        bool         is_initialized_; // TODO: more complete taf_state enum
        // Actual state machine can be hidden in implementation class.

        /// Indicate we have reached the end of input.
        bool                    end_of_frames_;

        SelectionOptionBehavior selectionOptionBehavior_;
};

}      // end namespace trajectoryanalysis
}      // end namespace gmx
#endif // GMX_TRAJECTORYANALYSIS_RUNNER_H
