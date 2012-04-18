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
 * Implements gmx::analysismodules::Distance.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "distance.h"

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

namespace gmx
{

namespace analysismodules
{

Distance::Distance()
    : _options("distance", "Distance calculation"), _avem(new AnalysisDataAverageModule())
{
}


Distance::~Distance()
{
}


Options &
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
    return _options;
}


void
Distance::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation & /*top*/)
{
    if (_sel[0].posCount() != 1)
    {
        GMX_THROW(InvalidInputError("The first selection does not define a single position"));
    }
    if (_sel[1].posCount() != 1)
    {
        GMX_THROW(InvalidInputError("The second selection does not define a single position"));
    }
    _data.setColumnCount(4);
    registerAnalysisDataset(&_data, "distance");

    _data.addModule(_avem);
    AnalysisDataPlotModulePointer _plotm(new AnalysisDataPlotModule());
    _plotm->setSettings(settings.plotSettings());
    _plotm->setFileName(_fnDist);
    _plotm->setTitle("Distance");
    _plotm->setXAxisIsTime();
    _plotm->setYLabel("Distance (nm)");
    _data.addModule(_plotm);
}


void
Distance::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(_data);
    const Selection    &sel1 = pdata->parallelSelection(_sel[0]);
    const Selection    &sel2 = pdata->parallelSelection(_sel[1]);
    rvec                dx;
    real                r;
    const SelectionPosition &p1 = sel1.position(0);
    const SelectionPosition &p2 = sel2.position(0);

    if (pbc != NULL)
    {
        pbc_dx(pbc, p1.x(), p2.x(), dx);
    }
    else
    {
        rvec_sub(p1.x(), p2.x(), dx);
    }
    r = norm(dx);
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, r);
    dh.setPoints(1, 3, dx);
    dh.finishFrame();
}


void
Distance::finishAnalysis(int /*nframes*/)
{
}


void
Distance::writeOutput()
{
    fprintf(stderr, "Average distance: %f\n", _avem->average(0));
    fprintf(stderr, "Std. deviation:   %f\n", _avem->stddev(0));
}


TrajectoryAnalysisModulePointer
Distance::create()
{
    return TrajectoryAnalysisModulePointer(new Distance());
}

} // namespace analysismodules

} // namespace gmx
