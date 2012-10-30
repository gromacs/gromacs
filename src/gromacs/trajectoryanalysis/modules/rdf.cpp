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
#include "rdf.h"

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

const char Rdf::name[] = "rdf";
const char Rdf::shortDescription[] =
    "Calculate radial distribution function";

Rdf::Rdf()
    : TrajectoryAnalysisModule(name, shortDescription),
      avem_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(1);
    registerAnalysisDataset(&data_, "rdf");
}


Rdf::~Rdf()
{
}


void
Rdf::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "This tool calculates the radial distribution function (RDF)",
        "for a set of output positions relative a set of reference",
        "positions.\n"
        "The RDF is defined such that the quantity\n",
        "\t(N/V)*RDF(r)*4pi*r^2*dr\n",
        "is the average number of output positions found in a",
        "spherical shell of width dr at a distance r from a reference",
        "position in the system, where N/V is the average output",
        "position density in the system.\n",
        "When the output and reference position sets are equal, it is",
        "sufficient to specify only the former."
        #ifdef GMX_LIB_MPI
        ,"\n Interactive selection is disabled for MPI builds."
        #endif
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(SelectionOption("select").required().valueCount(1)
                           .store(sel_));
    options->addOption(SelectionOption("reference").valueCount(1)
                       .description("Reference to calculate distances from")
                           .store(refsel_));
}


void
Rdf::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation & /*top*/)
{
    data_.addModule(avem_);
}


void
Rdf::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_[0]);
    const int           pos = sel.posCount();
    rvec                dx;
    real                r;
    std::vector<real>   dij((pos-1)*pos/2);
    double              sum = 0.0;
    int                 index = 0;


    for (int i = 0; i < pos-1; i++)
    {
        const SelectionPosition &pi = sel.position(i);
        for (int j = i+1; j < pos; j++)
        {
                const SelectionPosition &pj = sel.position(j);
                if (pbc != NULL)
                {
                    pbc_dx(pbc, pi.x(), pj.x(), dx);
                }
                else
                {
                    rvec_sub(pi.x(), pj.x(), dx);
                }
                dij[index] = norm(dx);
                sum += (double) dij[index];
                index++;
        }
    }

    r = (real) sum/dij.size();

    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, r);

    dh.finishFrame();
}


void
Rdf::finishAnalysis(int /*nframes*/)
{
}


void
Rdf::writeOutput()
{
    fprintf(stderr, "Average distance: %f\n", avem_->average(0));
    fprintf(stderr, "Std. deviation:   %f\n", avem_->stddev(0));
}

} // namespace analysismodules

} // namespace gmx
