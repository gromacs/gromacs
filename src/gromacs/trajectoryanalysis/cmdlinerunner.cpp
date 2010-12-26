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
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <copyrite.h>
#include <pbc.h>
#include <rmpbc.h>
#include <statutil.h>

#include "gromacs/errorreporting/standarderrorreporter.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/asciihelpwriter.h"
#include "gromacs/options/cmdlineparser.h"
#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selectioncollection.h"
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
        int parseOptions(TrajectoryAnalysisSettings *settings,
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
    delete _module;
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


int
TrajectoryAnalysisCommandLineRunner::Impl::parseOptions(
        TrajectoryAnalysisSettings *settings,
        TrajectoryAnalysisRunnerCommon *common,
        SelectionCollection *selections,
        Options *options,
        int *argc, char *argv[])
{
    StandardErrorReporter  errors;
    int rc;

    Options *moduleOptions = _module->initOptions(settings);
    if (moduleOptions == NULL)
    {
        GMX_ERROR(eeOutOfMemory,
                  "Could not allocate memory for option storage");
    }

    Options *commonOptions = common->initOptions();
    if (moduleOptions == NULL)
    {
        GMX_ERROR(eeOutOfMemory,
                  "Could not allocate memory for option storage");
    }

    Options *selectionOptions = selections->initOptions();
    if (selectionOptions == NULL)
    {
        GMX_ERROR(eeOutOfMemory,
                  "Could not allocate memory for option storage");
    }

    options->addSubSection(commonOptions);
    options->addSubSection(selectionOptions);
    options->addSubSection(moduleOptions);

    options->globalProperties().setSelectionCollection(selections);
    commonOptions->addDefaultOptions();

    {
        CommandLineParser  parser(options, &errors);
        rc = parser.parse(argc, argv);
        printHelp(*options, *common);
        if (rc != 0)
        {
            GMX_ERROR(rc, "Command-line option parsing failed, "
                          "see higher up for detailed error messages");
        }
        rc = options->finish(&errors);
        if (rc != 0)
        {
            GMX_ERROR(rc, "Command-line option parsing failed, "
                          "see higher up for detailed error messages");
        }
    }

    rc = common->initOptionsDone();
    if (rc != 0)
    {
        return rc;
    }
    rc = _module->initOptionsDone(settings);
    if (rc != 0)
    {
        return rc;
    }

    rc = common->initIndexGroups(selections);
    if (rc != 0)
    {
        return rc;
    }

    // TODO: Check whether the input is a pipe.
    bool bInteractive = true;
    rc = selections->parseRequestedFromStdin(bInteractive, &errors);
    common->doneIndexGroups(selections);
    return rc;
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
    delete _impl;
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
    int                       rc;

    CopyRight(stderr, argv[0]);

    SelectionCollection  selections(NULL);
    rc = selections.init();
    if (rc != 0)
    {
        return rc;
    }
    selections.setDebugLevel(_impl->_debugLevel);

    TrajectoryAnalysisSettings  settings;
    TrajectoryAnalysisRunnerCommon  common(&settings);

    Options  options(NULL, NULL);
    rc = _impl->parseOptions(&settings, &common, &selections, &options,
                             &argc, argv);
    if (rc != 0)
    {
        return rc;
    }

    rc = common.initTopology(&selections);
    if (rc != 0)
    {
        return rc;
    }
    rc = selections.compile();
    if (rc != 0)
    {
        return rc;
    }

    const TopologyInformation &topology = common.topologyInformation();
    rc = module->initAnalysis(topology);
    if (rc != 0)
    {
        return rc;
    }

    // Load first frame.
    rc = common.initFirstFrame();
    if (rc != 0)
    {
        return rc;
    }
    rc = module->initAfterFirstFrame(common.frame());
    if (rc != 0)
    {
        return rc;
    }

    t_pbc  pbc;
    t_pbc *ppbc = settings.hasPBC() ? &pbc : 0;

    int nframes = 0;
    TrajectoryAnalysisModuleData *pdata = NULL;
    rc = module->startFrames(NULL, selections, &pdata);
    if (rc != 0)
    {
        return rc;
    }
    do
    {
        rc = common.initFrame();
        if (rc != 0)
        {
            return rc;
        }
        t_trxframe &frame = common.frame();
        if (ppbc)
        {
            set_pbc(ppbc, topology.ePBC(), frame.box);
        }

        rc = selections.evaluate(&frame, ppbc);
        if (rc != 0)
        {
            return rc;
        }
        rc = module->analyzeFrame(nframes, frame, ppbc, pdata);
        if (rc != 0)
        {
            return rc;
        }

        nframes++;
    }
    while (common.readNextFrame());
    rc = module->finishFrames(pdata);
    if (rc != 0)
    {
        return rc;
    }
    if (pdata)
    {
        rc = pdata->finish();
        delete pdata;
        if (rc != 0)
        {
            return rc;
        }
    }

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
    rc = selections.evaluateFinal(nframes);
    if (rc != 0)
    {
        return rc;
    }

    rc = module->finishAnalysis(nframes);
    if (rc == 0)
    {
        rc = module->writeOutput();
    }

    return rc;
}

} // namespace gmx
