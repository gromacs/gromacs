/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "cmdlinerunner.h"

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"

#include "runnercommon.h"

namespace gmx
{

/********************************************************************
 * TrajectoryAnalysisCommandLineRunner::Impl
 */

class TrajectoryAnalysisCommandLineRunner::Impl
{
    public:
        explicit Impl(TrajectoryAnalysisModulePointer module)
            : module_(module), common_(&settings_),
              bUseDefaultGroups_(true)
        {
        }

        TrajectoryAnalysisModulePointer module_;
        TrajectoryAnalysisSettings      settings_;
        TrajectoryAnalysisRunnerCommon  common_;
        SelectionCollection             selections_;
        bool                            bUseDefaultGroups_;
};

/********************************************************************
 * TrajectoryAnalysisCommandLineRunner
 */

TrajectoryAnalysisCommandLineRunner::TrajectoryAnalysisCommandLineRunner(
        TrajectoryAnalysisModulePointer module)
    : impl_(new Impl(module))
{
}


TrajectoryAnalysisCommandLineRunner::~TrajectoryAnalysisCommandLineRunner()
{
}


void
TrajectoryAnalysisCommandLineRunner::setUseDefaultGroups(bool bUseDefaults)
{
    impl_->bUseDefaultGroups_ = bUseDefaults;
}


void TrajectoryAnalysisCommandLineRunner::init(
        CommandLineModuleSettings * /*settings*/)
{
}

void
TrajectoryAnalysisCommandLineRunner::initOptions(
        IOptionsContainer *options, ICommandLineOptionsModuleSettings *settings)
{
    boost::shared_ptr<TimeUnitBehavior>        timeUnitBehavior(
            new TimeUnitBehavior());
    boost::shared_ptr<SelectionOptionBehavior> selectionOptionBehavior(
            new SelectionOptionBehavior(&impl_->selections_,
                                        impl_->common_.topologyProvider(),
                                        impl_->bUseDefaultGroups_));
    settings->addOptionsBehavior(timeUnitBehavior);
    settings->addOptionsBehavior(selectionOptionBehavior);
    IOptionsContainer &commonOptions = options->addGroup();
    IOptionsContainer &moduleOptions = options->addGroup();

    impl_->module_->initOptions(&moduleOptions, &impl_->settings_);
    impl_->common_.initOptions(&commonOptions, timeUnitBehavior.get());
    selectionOptionBehavior->initOptions(&commonOptions);
}

void TrajectoryAnalysisCommandLineRunner::optionsFinished()
{
    impl_->common_.optionsFinished();
    impl_->module_->optionsFinished(&impl_->settings_);
}

int
TrajectoryAnalysisCommandLineRunner::run()
{
    TrajectoryAnalysisModule       *module     = impl_->module_.get();
    TrajectoryAnalysisSettings     &settings   = impl_->settings_;
    TrajectoryAnalysisRunnerCommon &common     = impl_->common_;
    SelectionCollection            &selections = impl_->selections_;

    common.initTopology();
    const TopologyInformation      &topology = common.topologyInformation();
    module->initAnalysis(settings, topology);

    // Load first frame.
    common.initFirstFrame();
    module->initAfterFirstFrame(settings, common.frame());

    t_pbc  pbc;
    t_pbc *ppbc = settings.hasPBC() ? &pbc : NULL;

    int    nframes = 0;
    AnalysisDataParallelOptions         dataOptions;
    TrajectoryAnalysisModuleDataPointer pdata(
            module->startFrames(dataOptions, selections));
    do
    {
        common.initFrame();
        t_trxframe &frame = common.frame();
        if (ppbc != NULL)
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        selections.evaluate(&frame, ppbc);
        module->analyzeFrame(nframes, frame, ppbc, pdata.get());
        module->finishFrameSerial(nframes);

        ++nframes;
    }
    while (common.readNextFrame());
    module->finishFrames(pdata.get());
    if (pdata.get() != NULL)
    {
        pdata->finish();
    }
    pdata.reset();

    if (common.hasTrajectory())
    {
        fprintf(stderr, "Analyzed %d frames, last time %.3f\n",
                nframes, common.frame().time);
    }
    else
    {
        fprintf(stderr, "Analyzed topology coordinates\n");
    }

    // Restore the maximal groups for dynamic selections.
    selections.evaluateFinal(nframes);

    module->finishAnalysis(nframes);
    module->writeOutput();

    return 0;
}

// static
int
TrajectoryAnalysisCommandLineRunner::runAsMain(
        int argc, char *argv[], ModuleFactoryMethod factory)
{
    auto runnerFactory = [factory]
    {
        return ICommandLineOptionsModulePointer(
                new TrajectoryAnalysisCommandLineRunner(factory()));
    };
    return ICommandLineOptionsModule::runAsMain(argc, argv, NULL, NULL, runnerFactory);
}

// static
void
TrajectoryAnalysisCommandLineRunner::registerModule(
        CommandLineModuleManager *manager, const char *name,
        const char *description, ModuleFactoryMethod factory)
{
    auto runnerFactory = [factory]
    {
        return ICommandLineOptionsModulePointer(
                new TrajectoryAnalysisCommandLineRunner(factory()));
    };
    ICommandLineOptionsModule::registerModule(manager, name, description, runnerFactory);
}

} // namespace gmx
