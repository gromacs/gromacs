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

#include <copyrite.h>
#include <pbc.h>
#include <rmpbc.h>
#include <statutil.h>

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/options/asciihelpwriter.h"
#include "gromacs/options/cmdlineparser.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptioninfo.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/runnercommon.h"

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
                       const TrajectoryAnalysisRunnerCommon &common);
        bool parseOptions(TrajectoryAnalysisSettings *settings,
                          TrajectoryAnalysisRunnerCommon *common,
                          SelectionCollection *selections,
                          Options *options,
                          int *argc, char *argv[]);

        TrajectoryAnalysisModule *_module;
        int                     _debugLevel;
};


TrajectoryAnalysisCommandLineRunner::Impl::Impl(
        TrajectoryAnalysisModule *module)
    : _module(module), _debugLevel(0)
{
}


TrajectoryAnalysisCommandLineRunner::Impl::~Impl()
{
}


void
TrajectoryAnalysisCommandLineRunner::Impl::printHelp(
        const Options &options,
        const TrajectoryAnalysisRunnerCommon &common)
{
    TrajectoryAnalysisRunnerCommon::HelpFlags flags = common.helpFlags();
    if (flags != 0)
    {
        AsciiHelpWriter(options)
            .setShowDescriptions(flags & TrajectoryAnalysisRunnerCommon::efHelpShowDescriptions)
            .setShowHidden(flags & TrajectoryAnalysisRunnerCommon::efHelpShowHidden)
            .writeHelp(stderr);
    }
}


bool
TrajectoryAnalysisCommandLineRunner::Impl::parseOptions(
        TrajectoryAnalysisSettings *settings,
        TrajectoryAnalysisRunnerCommon *common,
        SelectionCollection *selections,
        Options *options,
        int *argc, char *argv[])
{
    Options &moduleOptions = _module->initOptions(settings);
    Options &commonOptions = common->initOptions();
    Options &selectionOptions = selections->initOptions();

    options->addSubSection(&commonOptions);
    options->addSubSection(&selectionOptions);
    options->addSubSection(&moduleOptions);

    setSelectionCollectionForOptions(options, selections);

    {
        CommandLineParser  parser(options);
        try
        {
            parser.parse(argc, argv);
        }
        catch (const UserInputError &ex)
        {
            printHelp(*options, *common);
            throw;
        }
        printHelp(*options, *common);
        common->scaleTimeOptions(options);
        options->finish();
    }

    if (!common->initOptionsDone())
    {
        return false;
    }
    _module->initOptionsDone(settings);

    common->initIndexGroups(selections);

    // TODO: Check whether the input is a pipe.
    bool bInteractive = true;
    selections->parseRequestedFromStdin(bInteractive);
    common->doneIndexGroups(selections);

    return true;
}


/********************************************************************
 * TrajectoryAnalysisCommandLineRunner
 */

TrajectoryAnalysisCommandLineRunner::TrajectoryAnalysisCommandLineRunner(
        TrajectoryAnalysisModule *module)
    : _impl(new Impl(module))
{
}


TrajectoryAnalysisCommandLineRunner::~TrajectoryAnalysisCommandLineRunner()
{
}


void
TrajectoryAnalysisCommandLineRunner::setSelectionDebugLevel(int debuglevel)
{
    _impl->_debugLevel = 1;
}


int
TrajectoryAnalysisCommandLineRunner::run(int argc, char *argv[])
{
    TrajectoryAnalysisModule *module = _impl->_module;

    CopyRight(stderr, argv[0]);

    SelectionCollection  selections(NULL);
    selections.setDebugLevel(_impl->_debugLevel);

    TrajectoryAnalysisSettings  settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    Options  options(NULL, NULL);
    if (!_impl->parseOptions(&settings, &common, &selections, &options,
                             &argc, argv))
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
        common.initFrame();
        t_trxframe &frame = common.frame();
        if (ppbc != NULL)
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        selections.evaluate(&frame, ppbc);
        module->analyzeFrame(nframes, frame, ppbc, pdata.get());

        nframes++;
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

} // namespace gmx
