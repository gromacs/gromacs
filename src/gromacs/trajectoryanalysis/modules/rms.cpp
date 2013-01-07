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
 * AGWIP
 * Implements gmx::analysismodules::Rdf.
 *
 * \ingroup module_trajectoryanalysis
 */
#include "rms.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pbc.h"
#include "vec.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

const char Rms::name[] = "rms";
const char Rms::shortDescription[] =
    "Calculate (best-fit) root mean square deviations.";

Rms::Rms()
    : TrajectoryAnalysisModule(name, shortDescription),
        avem_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(1);
    registerAnalysisDataset(&data_, "rmsd_data");
}


Rms::~Rms()
{
}


void
Rms::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "rms tool is rms."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnRms_).defaultBasename("rms")
                           .description("RMS over time"));

    options->addOption(SelectionOption("select").required()
                       .valueCount(1)
                       .description("Selection for rms calculation.")
                       .store(sel_));
}

void
Rms::optionsFinished(Options * options, TrajectoryAnalysisSettings * /*settings*/)
{
}

void
Rms::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation &top)
{
    // here need to:
    // * setup reference conf.
    // * save reference to topology for later calls to fitting functions.

#ifdef GMX_LIB_MPI
    if (!fnRms_.empty() && mpi::isMaster())
#else
    if (!fnRms_.empty())
#endif
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnRms_);
        plotm->setTitle("RMS");
        plotm->setXLabel("time");
        avem_.addModule(plotm);
    }

#ifdef GMX_LIB_MPI
    if (mpi::isMaster())
    {
        fprintf(stderr, "master process: %d\n", getpid());
    }
    else
    {
        fprintf(stderr, "other process: %d\n", getpid());
    }
#endif

}

void
Rdf::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_[0]);
    const int           pos = sel.posCount();

    rvec                dx;
    int                 dsamples;

}


void
Rdf::finishAnalysis(int /*nframes*/)
{
}


void
Rdf::writeOutput()
{
    fprintf(stderr, "Writing some output to stderr!\n");
}

} // namespace analysismodules

} // namespace gmx
