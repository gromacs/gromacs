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
/*! \internal \file
 * \brief
 * Defines gmx::trajectoryanalysis::Runner.
 *
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "runner.h"

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

#include "runnercommon.h"


//#include "gromacs/commandline/cmdlineoptionsmodule.h"

namespace gmx
{
namespace trajectoryanalysis
{

using gmx::TrajectoryAnalysisModule;
using gmx::TrajectoryAnalysisModuleSharedPointer;
using gmx::IOptionsContainer;
using gmx::SelectionOptionBehavior;

using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

// Initializes an empty TrajectoryAnalysisSettings and uses it to
// initialize the TrajectoryAnalysisRunnerCommon member.
Runner::Runner() :
    module_(nullptr),
    settings_(),
    common_(&settings_),
    selections_(),
    pdata_(nullptr),
    nframes_(0),
    is_initialized_(false),
    end_of_frames_(false),
    selectionOptionBehavior_(&selections_, common_.topologyProvider())
{
    // TODO: Uncrustify wants to make my brace initialization ugly, at least
    // when it hits certain arguments...
}

Runner::~Runner()
{
    if (module_)
    {
        module_->finishFrames(pdata_.get());
    }
    if (pdata_.get() != nullptr)
    {
        pdata_->finish();
    }
    pdata_.reset();
}

// The Runner will share ownership of the module argument.
// TODO: clarify requirements and behavior for order of invocation of addModule()
// and initialize()
TrajectoryAnalysisModuleSharedPointer Runner::addModule(TrajectoryAnalysisModuleSharedPointer module)
{
    if (is_initialized_)
    {
        // Alert that we can't add modules anymore.
        return nullptr;
    }
    module_ = module;

    return module_;
}

/// Negotiate the options once module(s) have been registered.
/*! Receive an IOptionsContainer from the caller. It is assumed to already
 * have a FileNameOptionManager in place. After calling this method, the
 * caller should parse/process its input using the configured OptionsContainer.
 */
void Runner::registerOptions(Options &options)
{

    /* Set up IOptionsContainer and behaviors as in RunnerModule::initOptions(). Note that a caller-provided
       Options object is a IOptionsContainerWithSections is a
       IOptionsContainer.
       TODO: make this relationship clearer in doxygen.
     */

    // A TimeUnitBehavior is an IOptionsBehavior
    const unique_ptr<TimeUnitBehavior> time_unit_behavior(new TimeUnitBehavior());

    // TODO: extract options behaviors from ICommandLineOptionsModuleSettings and ICommandLineOptionsModuleSettings from TrajectoryAnalysisSettings
    /*
       // Retain settings for module(s)
       ICommandLineOptionsModuleSettings cli_settings{};
       cli_settings.addOptionsBehavior(selectionOptionBehavior);
       cli_settings.addOptionsBehavior(time_unit_behavior);
       settings_.setOptionsModuleSettings(cli_settings);
     */

    // This is where RunnerModule would call CommandLineOptionsModuleSettings::addOptionsBehavior(...), which causes OptionsBehaviorCollection::addBehavior(behavior), which causes behavior->initBehavior(options);
    selectionOptionBehavior_.initBehavior(&options);

    // Let the common_ add its options to the IOptionsContainer
    // and create for it a TimeUnitBehavior.
    IOptionsContainer &common_options = options.addGroup();
    common_.initOptions(&common_options, time_unit_behavior.get());
    // TODO: need methods to report this upstream...
    // Note no one owns the TimeUnitBehavior after this returns.
    // If we need it, we'll have to find a place for it to live
    // or it will be destroyed.

    // Let the module(s) register its options and receive settings.
    IOptionsContainer &module_options = options.addGroup();
    module_->initOptions(&module_options, &settings_);

    selectionOptionBehavior_.initOptions(&common_options);
}

// Parse user parameters, parse selections, and initialized selection variables in module
// parse(&options)
//options.finish();

/// Prepare the runner and modules to start iterating over frames.
/*! Part of initialization is to read the first frame with knowledge of what
 * information is needed by the modules. Thus, modules cannot be added without
 * reinitializing afterwards. For the moment, require forward progress...
 */
void Runner::initialize(Options &options)
{
    // Finalize options and check for errors. Module may adjust selections settings
    module_->optionsFinished(&settings_);

    // Selections are checked and module variables updated.
    // Finalize and compile selections.
    selectionOptionBehavior_.optionsFinishing(&options);
    selectionOptionBehavior_.optionsFinished();

    // finalize options for helper class and check for errors.
    common_.optionsFinished();

    //
    // Perform first frame initialization as in RunnerModule::run()
    //
    common_.initTopology(); // Need some error handling
    const TopologyInformation &topology = common_.topologyInformation();
    module_->initAnalysis(settings_, topology);

    // Load first frame.
    common_.initFirstFrame();      // need some error handling
    common_.initFrameIndexGroup(); // what is this?
    module_->initAfterFirstFrame(settings_, common_.frame());

    AnalysisDataParallelOptions dataOptions;
    pdata_ = decltype           (pdata_)(module_->startFrames(std::move(dataOptions), selections_));

    is_initialized_ = true;
}

bool Runner::next()
{
    if (end_of_frames_)
    {
        // There are no available input frames.
        return false;
    }
    if (!is_initialized_)
    {
        // Can't run if not initialized...
        // TODO: raise APIError?
        return false;
    }
    common_.initFrame();

    // Why isn't this const?
    t_trxframe                &frame = common_.frame();

    const TopologyInformation &topology = common_.topologyInformation();

    std::unique_ptr<t_pbc>     ppbc_ {
        nullptr
    };
    if (settings_.hasPBC())
    {
        // Need to preallocate memory for pbc
        t_pbc* pbc = new t_pbc;
        set_pbc(pbc, topology.ePBC(), frame.box);
        // Take ownership of any memory now pointed to by pbc
        ppbc_ = std::unique_ptr<t_pbc>(pbc);
        //TODO: update set_pbc function and/or documentation.
    }

    // TODO: convert the next two functions not to need non-const pointers to t_pbc
    selections_.evaluate(&frame, ppbc_.get());
    module_->analyzeFrame(nframes_, frame, ppbc_.get(), pdata_.get());
    module_->finishFrameSerial(nframes_);

    ++nframes_;

    if (!common_.readNextFrame())
    {
        // There are no more input frames
        end_of_frames_ = true;
    }
    return true;
}

int Runner::run()
{
    while (next())
    {
        // handle error and return 1;
    }
    ;
    return 0;
}

} // end namespace trajectoryanalysis
} // end namespace gmx
