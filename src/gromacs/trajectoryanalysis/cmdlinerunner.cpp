/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/cmdlinerunner.h"

#include <cstdio>

#include <memory>
#include <utility>

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/filestream.h"

#include "runnercommon.h"

namespace gmx
{
class CommandLineModuleSettings;

namespace
{

/********************************************************************
 * RunnerModule
 */

class RunnerModule : public ICommandLineOptionsModule
{
public:
    explicit RunnerModule(TrajectoryAnalysisModulePointer module) :
        module_(std::move(module)), common_(&settings_)
    {
    }

    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

    TrajectoryAnalysisModulePointer module_;
    TrajectoryAnalysisSettings      settings_;
    TrajectoryAnalysisRunnerCommon  common_;
    SelectionCollection             selections_;
};

void RunnerModule::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    std::shared_ptr<TimeUnitBehavior>        timeUnitBehavior(new TimeUnitBehavior());
    std::shared_ptr<SelectionOptionBehavior> selectionOptionBehavior(
            new SelectionOptionBehavior(&selections_, common_.topologyProvider()));
    settings->addOptionsBehavior(timeUnitBehavior);
    settings->addOptionsBehavior(selectionOptionBehavior);
    IOptionsContainer& commonOptions = options->addGroup();
    IOptionsContainer& moduleOptions = options->addGroup();

    settings_.setOptionsModuleSettings(settings);
    module_->initOptions(&moduleOptions, &settings_);
    settings_.setOptionsModuleSettings(nullptr);
    common_.initOptions(&commonOptions, timeUnitBehavior.get());
    selectionOptionBehavior->initOptions(&commonOptions);
}

void RunnerModule::optionsFinished()
{
    common_.optionsFinished();
    module_->optionsFinished(&settings_);
}

int RunnerModule::run()
{
    common_.initTopology();
    const TopologyInformation& topology = common_.topologyInformation();
    module_->initAnalysis(settings_, topology);

    // Load first frame.
    common_.initFirstFrame();
    common_.initFrameIndexGroup();
    module_->initAfterFirstFrame(settings_, common_.frame());

    t_pbc  pbc;
    t_pbc* ppbc = settings_.hasPBC() ? &pbc : nullptr;

    int                                 nframes = 0;
    AnalysisDataParallelOptions         dataOptions;
    TrajectoryAnalysisModuleDataPointer pdata(module_->startFrames(dataOptions, selections_));
    do
    {
        common_.initFrame();
        t_trxframe& frame = common_.frame();
        if (ppbc != nullptr)
        {
            set_pbc(ppbc, topology.pbcType(), frame.box);
        }

        selections_.evaluate(&frame, ppbc);
        module_->analyzeFrame(nframes, frame, ppbc, pdata.get());
        module_->finishFrameSerial(nframes);

        ++nframes;
    } while (common_.readNextFrame());
    module_->finishFrames(pdata.get());
    if (pdata.get() != nullptr)
    {
        pdata->finish();
    }
    pdata.reset();

    if (common_.hasTrajectory())
    {
        fprintf(stderr, "Analyzed %d frames, last time %.3f\n", nframes, common_.frame().time);
    }
    else
    {
        fprintf(stderr, "Analyzed topology coordinates\n");
    }

    // Restore the maximal groups for dynamic selections.
    selections_.evaluateFinal(nframes);

    module_->finishAnalysis(nframes);
    module_->writeOutput();

    return 0;
}

} // namespace

/********************************************************************
 * TrajectoryAnalysisCommandLineRunner
 */

// static
int TrajectoryAnalysisCommandLineRunner::runAsMain(int argc, char* argv[], const ModuleFactoryMethod& factory)
{
    auto runnerFactory = [factory] { return createModule(factory()); };
    return ICommandLineOptionsModule::runAsMain(argc, argv, nullptr, nullptr, runnerFactory);
}

// static
void TrajectoryAnalysisCommandLineRunner::registerModule(CommandLineModuleManager*  manager,
                                                         const char*                name,
                                                         const char*                description,
                                                         const ModuleFactoryMethod& factory)
{
    auto runnerFactory = [factory] { return createModule(factory()); };
    ICommandLineOptionsModule::registerModuleFactory(manager, name, description, runnerFactory);
}

// static
std::unique_ptr<ICommandLineOptionsModule>
TrajectoryAnalysisCommandLineRunner::createModule(TrajectoryAnalysisModulePointer module)
{
    return ICommandLineOptionsModulePointer(new RunnerModule(std::move(module)));
}

} // namespace gmx
