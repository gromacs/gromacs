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
#include "distance.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <pbc.h>
#include <vec.h>

// FIXME: This kind of hackery should not be necessary
#undef min
#undef max
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"

namespace gmx
{

namespace analysismodules
{

Distance::Distance()
    : _options("distance", "Distance calculation")
{
    _sel[0] = _sel[1] = NULL;
}


Distance::~Distance()
{
}


Options *
Distance::initOptions(TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "g_dist can calculate the distance between two positions as",
        "a function of time. The total distance and its",
        "x, y and z components are plotted.",
        NULL
    };

    _options.setDescription(desc);

    _options.addOption(FileNameOption("o").filetype(eftPlot).writeOnly()
                           .store(&_fnDist).defaultValue("dist"));
    _options.addOption(SelectionOption("select").required().valueCount(2)
                           .store(_sel));
    return &_options;
}


void
Distance::initAnalysis(const TopologyInformation & /*top*/)
{
    if (_sel[0]->posCount() != 1)
    {
        GMX_THROW(InvalidInputError("The first selection does not define a single position"));
    }
    if (_sel[1]->posCount() != 1)
    {
        GMX_THROW(InvalidInputError("The second selection does not define a single position"));
    }
    _data.setColumns(4);
    registerAnalysisDataset(&_data, "distance");

    _avem = new AnalysisDataAverageModule();
    _data.addModule(_avem);

    _plotm = new AnalysisDataPlotModule(_options);
    _plotm->setFileName(_fnDist);
    _plotm->setTitle("Distance");
    _plotm->setXLabel("Time [ps]");
    _plotm->setYLabel("Distance [nm]");
    _data.addModule(_plotm);
}


void
Distance::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle *dh = pdata->dataHandle("distance");
    Selection          *sel1 = pdata->parallelSelection(_sel[0]);
    Selection          *sel2 = pdata->parallelSelection(_sel[1]);
    rvec                dx;
    real                r;

    if (pbc != NULL)
    {
        pbc_dx(pbc, sel1->x(0), sel2->x(0), dx);
    }
    else
    {
        rvec_sub(sel1->x(0), sel2->x(0), dx);
    }
    r = norm(dx);
    dh->startFrame(frnr, fr.time);
    dh->addPoint(0, r);
    dh->addPoints(1, 3, dx);
    dh->finishFrame();
}


void
Distance::finishAnalysis(int /*nframes*/)
{
}


void
Distance::writeOutput()
{
    const real *ave;

    _avem->getData(0, NULL, &ave, NULL);
    fprintf(stderr, "Average distance: %f\n", ave[0]);
    fprintf(stderr, "Std. deviation:   %f\n", ave[1]);
}


TrajectoryAnalysisModule *
Distance::create()
{
    return new Distance();
}

} // namespace analysismodules

} // namespace gmx
