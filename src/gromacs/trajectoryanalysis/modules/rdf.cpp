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
        // avem_(new AnalysisDataAverageModule()),
        histm_(new AnalysisDataSimpleHistogramModule(histogramFromRange(0.0, 5.0).binCount(10)))
{
    data_.setColumnCount(1);
    data_.setMultipoint(true);
    registerAnalysisDataset(&data_, "dx_data");
}


Rdf::~Rdf()
{
}


void
Rdf::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "The rdf tool calculates the radial distribution function (RDF)",
        "for a set of output positions relative a set of reference",
        "positions.\n"
        "The RDF is defined such that the quantity\n",
        "\t(N/V)*RDF(r)*4pi*r^2*dr\n",
        "is the average number of output positions found in a",
        "spherical shell of width dr at a distance r from a reference",
        "position in the system, where N/V is the average output",
        "position density in the system.\n",
        "When the output and reference position sets are equal, it is",
        "sufficient to specify only the former. Self-distances are",
        "removed from histogram data."
        #ifdef GMX_LIB_MPI
        ,"\n Interactive selection is disabled for MPI builds."
        #endif
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnHist_).defaultBasename("hist")
                           .description("Computed histogram"));

    options->addOption(SelectionOption("select").required()
                       .valueCount(1)
                       .description("Selection for output.")
                       .store(sel_));

    options->addOption(SelectionOption("refsel")
                       .valueCount(1)
                       .description("Selection for reference.")
                       .store(refsel_));

    options->addOption(BooleanOption("surf")
                       .description("RDF relative surface of reference position set.")
                       .store(&surfref_));
}

void
Rdf::optionsFinished(Options * options, TrajectoryAnalysisSettings * /*settings*/)
{
    // check if a reference selection has been given.
    bRefSelectionSet_ = options->isSet("refsel");
}

void
Rdf::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation & /*top*/)
{
    // determine rdf mode
    if (surfref_)
    {
        rdfmode_ = SURFREF;
    }
    else if (bRefSelectionSet_)
    {
        rdfmode_ = REFSEL;
    }
    else
    {
        rdfmode_ = INTRA;
    }

    // check for empty selections
    if (sel_[0].posCount() == 0)
    {
        GMX_THROW(InvalidInputError("Selection does not define any positions."));
    }

    if ((rdfmode_==(REFSEL|SURFREF)) && !bRefSelectionSet_)
    {
        GMX_THROW(InvalidInputError("No reference selection given."));
    }

    if (bRefSelectionSet_ && refsel_[0].posCount() == 0)
    {
        GMX_THROW(InvalidInputError("Reference selection does not define any positions."));
    }

    data_.addModule(histm_);

#ifdef GMX_LIB_MPI
    if (!fnHist_.empty() && mpi::isMaster())
#else
    if (!fnHist_.empty())
#endif
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnHist_);
        plotm->setTitle("Histogram");
        plotm->setYLabel("counts");
        plotm->setXLabel("radius");
        histm_->averager().addModule(plotm);
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

    switch (rdfmode_)
    {
        case INTRA:
        {
            // dsamples = (pos-1)*pos/2;
            dh.startFrame(frnr, fr.time);

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
                    dh.setPoint(0, norm(dx));
                    dh.finishPointSet();
                }
            }
            dh.finishFrame();
        } break;

        case REFSEL:
        {
            const Selection    &refsel = pdata->parallelSelection(refsel_[0]);
            const int           rpos = refsel.posCount();
            real                s;

            dsamples = pos*rpos;
            dh.startFrame(frnr, fr.time);

            for (int i = 0; i < rpos; i++)
            {
                const SelectionPosition &pi = refsel.position(i);
                for (int j = 0; j < pos; j++)
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

                        s = norm(dx);
                        if (s > GMX_REAL_EPS)   // remove selection aliasing
                        {
                            dh.setPoint(0, s);
                            dh.finishPointSet();
                        }
                        else
                        {
                            dsamples--;
                        }
                }
            }
            dh.finishFrame();
        } break;

        case SURFREF:
        {
            const Selection    &refsel = pdata->parallelSelection(refsel_[0]);
            const int           rpos = refsel.posCount();
            real                s, min_s;

            dsamples = pos;
            dh.startFrame(frnr, fr.time);

            for (int j = 0; j < pos; j++)
            {
                const SelectionPosition &pj = sel.position(j);
                min_s = GMX_REAL_MAX;
                for (int i = 0; i < rpos; i++)
                {
                        const SelectionPosition &pi = refsel.position(i);
                        if (pbc != NULL)
                        {
                            pbc_dx(pbc, pi.x(), pj.x(), dx);
                        }
                        else
                        {
                            rvec_sub(pi.x(), pj.x(), dx);
                        }

                        s = norm(dx);
                        if (s < min_s && s > GMX_REAL_EPS)  // remove selection aliasing
                        {
                            min_s = s;
                        }
                }
                if (min_s < GMX_REAL_MAX)
                {
                    dh.setPoint(0, min_s);
                    dh.finishPointSet();
                }
                else
                {
                    dsamples--;
                }
            }
            dh.finishFrame();
        } break;
    }
}


void
Rdf::finishAnalysis(int /*nframes*/)
{
    histm_->averager().done();
}


void
Rdf::writeOutput()
{
    fprintf(stderr, "Writing some output!\n");
}

} // namespace analysismodules

} // namespace gmx
