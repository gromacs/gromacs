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
 * Implements gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "copyrite.h"
#include "pbc.h"
#include "rmpbc.h"
#include "statutil.h"

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptioninfo.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/runnercommon.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programinfo.h"
#include "external/mpp/gmxmpp.h"

namespace gmx
{

/********************************************************************
 * TrajectoryAnalysisCommandLineRunner::Impl
 */

class TrajectoryAnalysisCommandLineRunner::Impl
{
    public:
        Impl(TrajectoryAnalysisModule *module);
        ~Impl();

        void printHelp(const Options &options,
                       const TrajectoryAnalysisSettings &settings,
                       const TrajectoryAnalysisRunnerCommon &common);
        bool parseOptions(TrajectoryAnalysisSettings *settings,
                          TrajectoryAnalysisRunnerCommon *common,
                          SelectionCollection *selections,
                          int *argc, char *argv[]);

        TrajectoryAnalysisModule *module_;
        int                     debugLevel_;
        bool                    bPrintCopyright_;
};


TrajectoryAnalysisCommandLineRunner::Impl::Impl(
        TrajectoryAnalysisModule *module)
    : module_(module), debugLevel_(0), bPrintCopyright_(true)
{
}


TrajectoryAnalysisCommandLineRunner::Impl::~Impl()
{
}


void
TrajectoryAnalysisCommandLineRunner::Impl::printHelp(
        const Options &options,
        const TrajectoryAnalysisSettings &settings,
        const TrajectoryAnalysisRunnerCommon &common)
{
    if (mpi::isMaster())
    {
        TrajectoryAnalysisRunnerCommon::HelpFlags flags = common.helpFlags();
        if (flags != 0)
        {
            HelpWriterContext context(&File::standardError(),
                    eHelpOutputFormat_Console);
            CommandLineHelpWriter(options)
            .setShowDescriptions(flags & TrajectoryAnalysisRunnerCommon::efHelpShowDescriptions)
            .setShowHidden(flags & TrajectoryAnalysisRunnerCommon::efHelpShowHidden)
            .setTimeUnitString(settings.timeUnitManager().timeUnitAsString())
            .writeHelp(context);
        }
    }
}


bool
TrajectoryAnalysisCommandLineRunner::Impl::parseOptions(
        TrajectoryAnalysisSettings *settings,
        TrajectoryAnalysisRunnerCommon *common,
        SelectionCollection *selections,
        int *argc, char *argv[])
{
    Options options(NULL, NULL);
    Options moduleOptions(module_->name(), module_->description());
    Options commonOptions("common", "Common analysis control");
    Options selectionOptions("selection", "Common selection control");
    module_->initOptions(&moduleOptions, settings);
    common->initOptions(&commonOptions);
    selections->initOptions(&selectionOptions);

    options.addSubSection(&commonOptions);
    options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    SelectionOptionManager seloptManager(selections);
    setManagerForSelectionOptions(&options, &seloptManager);

    {
        CommandLineParser  parser(&options);
        try
        {
            parser.parse(argc, argv);
        }
        catch (const UserInputError &ex)
        {
            printHelp(options, *settings, *common);
            throw;
        }
        printHelp(options, *settings, *common);
        common->scaleTimeOptions(&options);
        options.finish();
    }

    if (!common->optionsFinished(&commonOptions))
    {
        return false;
    }
    module_->optionsFinished(&moduleOptions, settings);

    common->initIndexGroups(selections);

    // TODO: Check whether the input is a pipe.
    bool bInteractive = true;
    seloptManager.parseRequestedFromStdin(bInteractive);
    common->doneIndexGroups(selections);

    return true;
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
TrajectoryAnalysisCommandLineRunner::setPrintCopyright(bool bPrint)
{
    impl_->bPrintCopyright_ = bPrint;
}


void
TrajectoryAnalysisCommandLineRunner::setSelectionDebugLevel(int debuglevel)
{
    impl_->debugLevel_ = 1;
}


int
TrajectoryAnalysisCommandLineRunner::run(int argc, char *argv[])
{
    TrajectoryAnalysisModule *module = impl_->module_;

    if (impl_->bPrintCopyright_)
    {
        CopyRight(stderr, argv[0]);
    }

    SelectionCollection  selections;
    selections.setDebugLevel(impl_->debugLevel_);

    TrajectoryAnalysisSettings  settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    if (!impl_->parseOptions(&settings, &common, &selections, &argc, argv))
    {
        return 0;
    }

    common.initTopology(&selections);
    selections.compile();

    const TopologyInformation &topology = common.topologyInformation();
    module->initAnalysis(settings, topology);

    // Load first frame.
    common.initFirstFrame();
    module->initAfterFirstFrame(common.frame());
    t_pbc  pbc;
    t_pbc *ppbc = settings.hasPBC() ? &pbc : NULL;

    int nframes = 0;
    AnalysisDataParallelOptions dataOptions;
    TrajectoryAnalysisModuleDataPointer pdata(
            module->startFrames(dataOptions, selections));
    do
    {
        size_t globalIndex = common.getCurrentFrameNumber();
        common.initFrame();
        t_trxframe &frame = common.frame();
        if (ppbc != NULL)
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        selections.evaluate(&frame, ppbc);
        module->analyzeFrame(globalIndex, frame, ppbc, pdata.get());

        nframes++;
    }
    while (common.readNextFrame());
    module->finishFrames(pdata.get());
    if (pdata.get() != NULL)
    {
        pdata->finish();
    }
    pdata.reset();
#ifdef GMX_LIB_MPI
    //TODO: Replace the streams with MPP send/receive calls
    //TODO: This statement incorrectly assumes the last core always has the last frame.
    std::vector<int> allframes(mpi::comm::world.size());
    if (mpi::comm::world.rank()==mpi::comm::world.size()-1)
    {
        mpi::comm::world(0) << common.frame().time;
    }
    else if (mpi::comm::world.isMaster())
    {
        mpi::comm::world(mpi::comm::world.size()-1) >> common.frame().time;
    }
    mpi::comm::world.gather(nframes, allframes, 0);
    for (int i=1; i<mpi::comm::world.size(); i++)
    {
        nframes += allframes[i];
    }
#endif
    if (mpi::isMaster())
    {
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
    }

    return 0;
}


void
TrajectoryAnalysisCommandLineRunner::writeHelp(const HelpWriterContext &context)
{
    // TODO: This method duplicates some code from run() and Impl::printHelp().
    // See how to best refactor it to share the common code.
    SelectionCollection             selections;
    TrajectoryAnalysisSettings      settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    Options options(NULL, NULL);
    Options moduleOptions(impl_->module_->name(), impl_->module_->description());
    Options commonOptions("common", "Common analysis control");
    Options selectionOptions("selection", "Common selection control");

    impl_->module_->initOptions(&moduleOptions, &settings);
    common.initOptions(&commonOptions);
    selections.initOptions(&selectionOptions);

    options.addSubSection(&commonOptions);
    options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    SelectionOptionManager seloptManager(&selections);
    setManagerForSelectionOptions(&options, &seloptManager);

    CommandLineHelpWriter(options)
        .setShowDescriptions(true)
        .setTimeUnitString(settings.timeUnitManager().timeUnitAsString())
        .writeHelp(context);
}

} // namespace gmx
