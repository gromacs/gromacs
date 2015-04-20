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
#include "gromacs/commandline/cmdlinemodule.h"
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
#include "gromacs/utility/file.h"
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
        class RunnerCommandLineModule;

        Impl(TrajectoryAnalysisModule *module);
        ~Impl();

        void parseOptions(TrajectoryAnalysisSettings *settings,
                          TrajectoryAnalysisRunnerCommon *common,
                          SelectionCollection *selections,
                          int *argc, char *argv[]);

        TrajectoryAnalysisModule *module_;
        bool                      bUseDefaultGroups_;
        int                       debugLevel_;
};


TrajectoryAnalysisCommandLineRunner::Impl::Impl(
        TrajectoryAnalysisModule *module)
    : module_(module), bUseDefaultGroups_(true), debugLevel_(0)
{
}


TrajectoryAnalysisCommandLineRunner::Impl::~Impl()
{
}


void
TrajectoryAnalysisCommandLineRunner::Impl::parseOptions(
        TrajectoryAnalysisSettings *settings,
        TrajectoryAnalysisRunnerCommon *common,
        SelectionCollection *selections,
        int *argc, char *argv[])
{
    FileNameOptionManager  fileoptManager;
    SelectionOptionManager seloptManager(selections);
    Options                options(NULL, NULL);
    Options                moduleOptions(module_->name(), module_->description());
    Options                commonOptions("common", "Common analysis control");
    Options                selectionOptions("selection", "Common selection control");

    options.addManager(&fileoptManager);
    options.addManager(&seloptManager);
    options.addSubSection(&commonOptions);
    options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    module_->initOptions(&moduleOptions, settings);
    common->initOptions(&commonOptions);
    selections->initOptions(&selectionOptions);

    {
        CommandLineParser  parser(&options);
        // TODO: Print the help if user provides an invalid option?
        // Or just add a message advising the user to invoke the help?
        parser.parse(argc, argv);
        common->scaleTimeOptions(&options);
        options.finish();
    }

    common->optionsFinished(&commonOptions);
    module_->optionsFinished(&moduleOptions, settings);

    common->initIndexGroups(selections, bUseDefaultGroups_);

    const bool bInteractive = File::standardInput().isInteractive();
    seloptManager.parseRequestedFromStdin(bInteractive);
    common->doneIndexGroups(selections);

    common->initTopology(selections);
    selections->compile();
}


/********************************************************************
 * TrajectoryAnalysisCommandLineRunner
 */

TrajectoryAnalysisCommandLineRunner::TrajectoryAnalysisCommandLineRunner(
        TrajectoryAnalysisModule *module)
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


void
TrajectoryAnalysisCommandLineRunner::setSelectionDebugLevel(int debuglevel)
{
    impl_->debugLevel_ = debuglevel;
}


int
TrajectoryAnalysisCommandLineRunner::run(int argc, char *argv[])
{
    TrajectoryAnalysisModule *module = impl_->module_;

    SelectionCollection       selections;
    selections.setDebugLevel(impl_->debugLevel_);

    TrajectoryAnalysisSettings      settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    impl_->parseOptions(&settings, &common, &selections, &argc, argv);

    const TopologyInformation &topology = common.topologyInformation();
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


void
TrajectoryAnalysisCommandLineRunner::writeHelp(const CommandLineHelpContext &context)
{
    // TODO: This method duplicates some code from run().
    // See how to best refactor it to share the common code.
    SelectionCollection             selections;
    TrajectoryAnalysisSettings      settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    SelectionOptionManager          seloptManager(&selections);
    Options                         options(NULL, NULL);
    Options                         moduleOptions(impl_->module_->name(), impl_->module_->description());
    Options                         commonOptions("common", "Common analysis control");
    Options                         selectionOptions("selection", "Common selection control");

    options.addManager(&seloptManager);
    options.addSubSection(&commonOptions);
    options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    impl_->module_->initOptions(&moduleOptions, &settings);
    common.initOptions(&commonOptions);
    selections.initOptions(&selectionOptions);

    CommandLineHelpWriter(options)
        .setShowDescriptions(true)
        .setTimeUnitString(settings.timeUnitManager().timeUnitAsString())
        .writeHelp(context);
}


/*! \internal \brief
 * Command line module for a trajectory analysis module.
 *
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisCommandLineRunner::Impl::RunnerCommandLineModule
    : public CommandLineModuleInterface
{
    public:
        /*! \brief
         * Constructs a module.
         *
         * \param[in] name         Name for the module.
         * \param[in] description  One-line description for the module.
         * \param[in] factory      Factory method to create the analysis module.
         *
         * Does not throw.  This is important for correct implementation of
         * runAsMain().
         */
        RunnerCommandLineModule(const char *name, const char *description,
                                ModuleFactoryMethod factory)
            : name_(name), description_(description), factory_(factory)
        {
        }

        virtual const char *name() const { return name_; }
        virtual const char *shortDescription() const { return description_; };

        virtual void init(CommandLineModuleSettings *settings);
        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const CommandLineHelpContext &context) const;

    private:
        const char             *name_;
        const char             *description_;
        ModuleFactoryMethod     factory_;

        GMX_DISALLOW_COPY_AND_ASSIGN(RunnerCommandLineModule);
};

void TrajectoryAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::init(
        CommandLineModuleSettings * /*settings*/)
{
}

int TrajectoryAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::run(
        int argc, char *argv[])
{
    TrajectoryAnalysisModulePointer     module(factory_());
    TrajectoryAnalysisCommandLineRunner runner(module.get());
    return runner.run(argc, argv);
}

void TrajectoryAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::writeHelp(
        const CommandLineHelpContext &context) const
{
    TrajectoryAnalysisModulePointer     module(factory_());
    TrajectoryAnalysisCommandLineRunner runner(module.get());
    runner.writeHelp(context);
}

// static
int
TrajectoryAnalysisCommandLineRunner::runAsMain(
        int argc, char *argv[], ModuleFactoryMethod factory)
{
    Impl::RunnerCommandLineModule module(NULL, NULL, factory);
    return CommandLineModuleManager::runAsMainSingleModule(argc, argv, &module);
}

// static
void
TrajectoryAnalysisCommandLineRunner::registerModule(
        CommandLineModuleManager *manager, const char *name,
        const char *description, ModuleFactoryMethod factory)
{
    CommandLineModulePointer module(
            new Impl::RunnerCommandLineModule(name, description, factory));
    manager->addModule(move(module));
}

} // namespace gmx
