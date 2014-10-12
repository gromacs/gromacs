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
#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionmanager.h"
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
        Impl(TrajectoryAnalysisModulePointer module);
        ~Impl();

        TrajectoryAnalysisModulePointer module_;
        bool                            bUseDefaultGroups_;
        int                             debugLevel_;
        SelectionCollection             selections_;
        SelectionOptionManager          selOptManager_;
        TrajectoryAnalysisSettings      settings_;
        TrajectoryAnalysisRunnerCommon  common_;
};


TrajectoryAnalysisCommandLineRunner::Impl::Impl(
        TrajectoryAnalysisModulePointer module)
    : module_(module), bUseDefaultGroups_(true), debugLevel_(0),
      selOptManager_(&selections_), common_(&settings_)
{
}


TrajectoryAnalysisCommandLineRunner::Impl::~Impl()
{
    module_.reset();
}


/********************************************************************
 * TrajectoryAnalysisCommandLineRunner
 */

TrajectoryAnalysisCommandLineRunner::TrajectoryAnalysisCommandLineRunner(
        TrajectoryAnalysisCommandLineRunner::ModuleFactoryMethod factory)
    : impl_(new Impl(factory()))
{
}

TrajectoryAnalysisCommandLineRunner::TrajectoryAnalysisCommandLineRunner(
        TrajectoryAnalysisCommandLineRunner::ModuleFactoryFunctor *functor)
    : impl_(new Impl((*functor)()))
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


void
TrajectoryAnalysisCommandLineRunner::setSelectionDebugLevel(int debuglevel)
{
    impl_->debugLevel_ = debuglevel;
}

void TrajectoryAnalysisCommandLineRunner::init(CommandLineModuleSettings * /*settings*/)
{
}

void
TrajectoryAnalysisCommandLineRunner::initOptions(IOptionsContainer                 *options,
                                                 ICommandLineOptionsModuleSettings *settings)
{
    impl_->selections_.setDebugLevel(impl_->debugLevel_);

    options->addManager(&impl_->selOptManager_);

    IOptionsContainer &commonOptions = options->addGroup();
    IOptionsContainer &moduleOptions = options->addGroup();

    impl_->module_->initOptions(&moduleOptions, &impl_->settings_);
    impl_->common_.initOptions(&commonOptions);
    impl_->selections_.initOptions(&commonOptions);

    const char *help[] = { impl_->settings_.helpText().data() };
    settings->setHelpText(help);
}

void
TrajectoryAnalysisCommandLineRunner::optionsFinished()
{
    // TODO: scaleTimeOptions? (maybe pass options here too, with default arg)
    impl_->common_.optionsFinished();
    impl_->module_->optionsFinished(&impl_->settings_);

    impl_->common_.initIndexGroups(&impl_->selections_, impl_->bUseDefaultGroups_);

    const bool bInteractive = StandardInputStream::instance().isInteractive();
    impl_->selOptManager_.parseRequestedFromStdin(bInteractive);
    impl_->common_.doneIndexGroups(&impl_->selections_);

    impl_->common_.initTopology(&impl_->selections_);
    impl_->selections_.compile();
}


int
TrajectoryAnalysisCommandLineRunner::run()
{
    const TopologyInformation &topology = impl_->common_.topologyInformation();
    impl_->module_->initAnalysis(impl_->settings_, topology);

    // Load first frame.
    impl_->common_.initFirstFrame();
    impl_->module_->initAfterFirstFrame(impl_->settings_, impl_->common_.frame());

    t_pbc  pbc;
    t_pbc *ppbc = impl_->settings_.hasPBC() ? &pbc : NULL;

    int    nframes = 0;
    AnalysisDataParallelOptions         dataOptions;
    TrajectoryAnalysisModuleDataPointer pdata(
            impl_->module_->startFrames(dataOptions, impl_->selections_));
    do
    {
        impl_->common_.initFrame();
        t_trxframe &frame = impl_->common_.frame();
        if (ppbc != NULL)
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        impl_->selections_.evaluate(&frame, ppbc);
        impl_->module_->analyzeFrame(nframes, frame, ppbc, pdata.get());
        impl_->module_->finishFrameSerial(nframes);

        ++nframes;
    }
    while (impl_->common_.readNextFrame());
    impl_->module_->finishFrames(pdata.get());
    if (pdata.get() != NULL)
    {
        pdata->finish();
    }
    pdata.reset();

    if (impl_->common_.hasTrajectory())
    {
        fprintf(stderr, "Analyzed %d frames, last time %.3f\n",
                nframes, impl_->common_.frame().time);
    }
    else
    {
        fprintf(stderr, "Analyzed topology coordinates\n");
    }

    // Restore the maximal groups for dynamic selections.
    impl_->selections_.evaluateFinal(nframes);

    impl_->module_->finishAnalysis(nframes);
    impl_->module_->writeOutput();
    return 0;
}

// static
int
TrajectoryAnalysisCommandLineRunner::runAsMain(
        int argc, char *argv[], ModuleFactoryMethod factory)
{
    TrajectoryAnalysisCommandLineRunner *module = new TrajectoryAnalysisCommandLineRunner(factory);
    return ICommandLineOptionsModule::runAsMain(argc, argv, module);
}

// static
int
TrajectoryAnalysisCommandLineRunner::runAsMain(
        int argc, char *argv[], ModuleFactoryFunctor *factory)
{
    TrajectoryAnalysisCommandLineRunner *module = new TrajectoryAnalysisCommandLineRunner(factory);
    return ICommandLineOptionsModule::runAsMain(argc, argv, module);
}

// static
void
TrajectoryAnalysisCommandLineRunner::registerModule(
        CommandLineModuleManager *manager, const char *name,
        const char *description, ModuleFactoryMethod factory)
{
    ICommandLineOptionsModule::registerModule(manager, name, description,
                                              new TrajectoryAnalysisCommandLineRunner(factory));
}

} // namespace gmx
