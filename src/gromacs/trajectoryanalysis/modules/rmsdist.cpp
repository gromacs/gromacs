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
 * Implements gmx::analysismodules::RmsDist.
 *
 * \ingroup module_trajectoryanalysis
 */
#include "rmsdist.h"

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

#include "smalloc.h"
#include "rmpbc.h"


namespace gmx
{

namespace analysismodules
{

const char RmsDist::name[] = "rmsdist";
const char RmsDist::shortDescription[] =
    "Calculate RMS distance differences.";

RmsDist::RmsDist()
    : TrajectoryAnalysisModule(name, shortDescription),
        avem_(new AnalysisDataAverageModule())
{
    data_.setColumnCount(1);
    registerAnalysisDataset(&data_, "rmsd_data");
}


RmsDist::~RmsDist()
{
}


void
RmsDist::initOptions(Options *options, TrajectoryAnalysisSettings * settings)
{
    static const char *const desc[] = {
        "rmsdist calculates RMS deviation of (intra-molecular) atom distances relative a reference structure."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o")
                       .description("Plot of RMSD over time.")
                       .filetype(eftPlot).outputFile()
                       .store(&fnRmsDist_).defaultBasename("rms"));

    options->addOption(SelectionOption("select")
                       .description("Selection for rms calculation.")
                       .required()
                       .valueCount(1)
                       // .onlyStatic()
                       .onlyAtoms()
                       .store(&sel_));

    options->addOption(BooleanOption("uw")
                       .description("put weights to 1.")
                       .defaultValue(TRUE)
                       .store(&bUnitWeights_));

    // topology with coords must be provided for use as reference in RMS calc.
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);

    // // turn off PBC distances to calc intra-mol distances when having whole molecules
    // if (settings->hasRmPCB())
    // {
    //     settings->setPCB(false);
    // }
}

void
RmsDist::optionsFinished(Options * options, TrajectoryAnalysisSettings * /*settings*/)
{
}

void
RmsDist::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation &topInfo)
{
    // here need to:
    // * setup reference conformation from top file.
    // * save reference to topology for later calls to fitting functions.
    matrix          box;


    pRefTop_ = topInfo.topology();
    topAtoms_ = pRefTop_->atoms.nr;
    pRefX_ = NULL;
    // pRefM_ = NULL;

    // get coords for reference structure
    (void) topInfo.getTopologyConf(&pRefX_, box);

    // make reference structure whole
    if (settings.hasRmPBC())
    {
        gmx_rmpbc_t     gpbc = NULL;

        gpbc = gmx_rmpbc_init(&pRefTop_->idef, topInfo.ePBC(), pRefTop_->atoms.nr, box);
        (void) gmx_rmpbc(gpbc, pRefTop_->atoms.nr, box, pRefX_);
        (void) gmx_rmpbc_done(gpbc);
    }

    // // setup weights
    // snew(pRefM_, pRefTop_->atoms.nr);
    // bool bMass = FALSE;
    // for(int i=0; i<pRefTop_->atoms.nr; i++)
    // {
    //     if (bUnitWeights_)
    //         pRefM_[i] = 1.0;
    //     else
    //         pRefM_[i] = pRefTop_->atoms.atom[i].m;
    //     bMass = bMass || (pRefTop_->atoms.atom[i].m != 0);
    // }
    // if (!bMass)
    // {
    //     fprintf(stderr,"All masses in the fit group are 0, using masses of 1\n");
    //     for(int i=0; i<pRefTop_->atoms.nr; i++)
    //     {
    //         pRefM_[i] = 1.0;
    //     }
    // }

/* setup plot output file */
#ifdef GMX_LIB_MPI
    if (!fnRmsDist_.empty() && mpi::isMaster())
#else
    if (!fnRmsDist_.empty())
#endif
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnRmsDist_);
        plotm->setTitle("RMS");
        plotm->setXAxisIsTime();
        data_.addModule(plotm);
    }

    // also calc average
    data_.addModule(avem_);

// #ifdef GMX_LIB_MPI
//     if (mpi::isMaster())
//     {
//         fprintf(stderr, "master process: %d\n", getpid());
//     }
//     else
//     {
//         fprintf(stderr, "other process: %d\n", getpid());
//     }
// #endif

}

void
RmsDist::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_);

    const int                   selatoms = sel.atomCount();
    const ConstArrayRef< int >  selAtomIndsArray = sel.atomIndices();

    // atom_id *   selind;     // indicies of selected atoms (ref. to topology indexing)
    rvec *      x;          // coordinates of selected atoms
    rvec *      xp;         // coordinates of reference atoms corresponding to selected atoms
    real *      m;          // masses of selected atoms

    rvec        dx, dxp;
    real        r, rp;
    real        rms_val = 0.0;  // for accumulating RSMD

    const real  avging_prefactor = 2.0 / (selatoms * (selatoms-1));
    real        avging_massprefactor = 0.0;

    /* first put selection stuff in legacy format
        so that we can use the old routines             */

    // // selection indecies
    // (void) snew(selind, selatoms);
    // for(int i=0; i<selatoms; i++)
    // {
    //     selind[i] = (atom_id) selAtomIndsArray.at(i);
    // }

    // coords and masses
    (void) snew(x, selatoms);
    (void) snew(xp, selatoms);
    (void) snew(m, selatoms);
    for (int i=0; i<selatoms; i++)
    {
        const SelectionPosition &spi = sel.position(i);
        const rvec &spix = spi.x();
        const atom_id &selindi = (atom_id) selAtomIndsArray.at(i);
        for (int m=0; m<DIM; m++)
        {
            x[i][m]     = spix[m];
            xp[i][m]    = pRefX_[selindi][m];
        }

        if (bUnitWeights_)
            m[i] = 1.0;
        else
            m[i] = spi.mass();
    }

    for (int i=0; i<selatoms-1; i++)
    {
        const rvec &xi = x[i];
        const rvec &xpi = xp[i];
        const real &mi = m[i];
        for (int j=i+1; j<selatoms; j++)
        {
            const rvec &xj = x[j];
            const rvec &xpj = xp[j];
            const real &mj = m[j];
            if (pbc != NULL)
            {
                pbc_dx(pbc, xi, xj, dx);
                pbc_dx(pbc, xpi, xpj, dxp);
            }
            else
            {
                rvec_sub(xi, xj, dx);
                rvec_sub(xpi, xpj, dxp);
            }
            r = norm(dx);
            rp = norm(dxp);

            rms_val += avging_prefactor * mi*mj * (r-rp)*(r-rp);
            avging_massprefactor += mi*mj;
        }
    }

    //  correct mass weight averaging
    rms_val = rms_val * 1.0/(avging_prefactor * avging_massprefactor);

    /* write the result */
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, rms_val);
    dh.finishFrame();

    // release temp frame stuff
    // (void) sfree(selind);
    (void) sfree(x);
    (void) sfree(xp);
    (void) sfree(m);
}


void
RmsDist::finishAnalysis(int /*nframes*/)
{
}


void
RmsDist::writeOutput()
{
    fprintf(stderr, "Average RMS: %f\n", avem_->average(0));
    fprintf(stderr, "Std. RMS:   %f\n", avem_->stddev(0));
}

} // namespace analysismodules

} // namespace gmx
