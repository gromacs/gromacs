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

#include "zleep.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

const char Zleep::name[] = "zleep";
const char Zleep::shortDescription[] =
    "Short description: Zleep...";

Zleep::Zleep()
    : TrajectoryAnalysisModule(name, shortDescription)
{
    data_.setColumnCount(1);
    registerAnalysisDataset(&data_, "flag");
}

Zleep::~Zleep()
{
}


void
Zleep::initOptions(Options *options,
                   TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Description: Zleep just zleeps a little while! :)"
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(DoubleOption("z").required()
                        .store(&zleep_time_).timeValue()
                        .description("Time spent zleeping. (us)"));
}


void
Zleep::initAnalysis(const TrajectoryAnalysisSettings &settings,
                    const TopologyInformation & /*top*/)
{
}


void
Zleep::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);

    dh.startFrame(frnr, fr.time);
        usleep(zleep_time_);           // dont do any work! HAHA!!!
        dh.setPoint(0, 1.0);    // just a flag.
    dh.finishFrame();
}


void
Zleep::finishAnalysis(int /*nframes*/)
{
}


void
Zleep::writeOutput()
{
    fprintf(stderr, "Output: Done zleeping!\n");
}

} // namespace analysismodules

} // namespace gmx
