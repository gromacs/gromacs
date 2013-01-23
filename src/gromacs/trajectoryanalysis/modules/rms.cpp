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
 * Implements gmx::analysismodules::Rms.
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

#include "smalloc.h"
#include "do_fit.h"

namespace gmx
{

namespace analysismodules
{

const char Rms::name[] = "rms";
const char Rms::shortDescription[] =
    "Calculate (best-fit) RMS deviations.";

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
Rms::initOptions(Options *options, TrajectoryAnalysisSettings * settings)
{
    static const char *const desc[] = {
        "rms calculates RMS deviation of atomic positions relative a reference structure."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o")
                       .description("Plot of RMS over time.")
                       .filetype(eftPlot).outputFile()
                       .store(&fnRms_).defaultBasename("rms"));

    options->addOption(SelectionOption("select")
                       .description("Selection for rms calculation.")
                       .required()
                       .valueCount(1)
                       // .onlyStatic()
                       .onlyAtoms()
                       .store(&sel_));

    options->addOption(BooleanOption("mw")
                       .description("Use mass weighted RMS.")
                       .defaultValue(true)
                       .store(&bRMSUseMassWeights_));

    options->addOption(BooleanOption("fit")
                       .description("Do translational and rotational least-squares fit.")
                       .defaultValue(true)
                       .store(&bDoFit_));

    options->addOption(SelectionOption("fitsel")
                       .description("Selection for least-squares fitting.")
                       .valueCount(1)
                       // .onlyStatic()
                       .onlyAtoms()
                       .store(&fitsel_));

    options->addOption(BooleanOption("mwfit")
                       .description("Use mass weighted fit.")
                       .defaultValue(true)
                       .store(&bFitUseMassWeights_));

    // topology with coords must be provided for use as reference in RMS calc.
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}

void
Rms::optionsFinished(Options * options, TrajectoryAnalysisSettings * /*settings*/)
{
}

void
Rms::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation &topInfo)
{
    // here need to:
    // * setup reference conformation from top file.
    // * save reference to topology for later calls to fitting functions.
    matrix          box;

    pRefTop_ = topInfo.topology();
    topAtoms_ = pRefTop_->atoms.nr;
        fprintf(stderr, "topAtoms_: %d\n", topAtoms_);
    pRefX_ = NULL;
    pRefM_ = NULL;
    pRefMU_ = NULL;

    // get coords for reference structure
    (void) topInfo.getTopologyConf(&pRefX_, box);

    // make reference structure whole. Q: does this make sense?!
    if (settings.hasRmPBC())
    {
        gmx_rmpbc_t     gpbc = NULL;

        gpbc = gmx_rmpbc_init(&pRefTop_->idef, topInfo.ePBC(), pRefTop_->atoms.nr, box);
        (void) gmx_rmpbc(gpbc, pRefTop_->atoms.nr, box, pRefX_);
        (void) gmx_rmpbc_done(gpbc);
    }

    // setup topology weights
    snew(pRefM_, topAtoms_);
    snew(pRefMU_, topAtoms_);

    bool bMass = FALSE;
    for(int i=0; i<topAtoms_; i++)
    {
        pRefM_[i] = pRefTop_->atoms.atom[i].m;
        pRefMU_[i] = 1.0;
        bMass = bMass || (pRefTop_->atoms.atom[i].m != 0);
    }
    if (!bMass)
    {
        fprintf(stderr,"All masses in the fit group are 0, using masses of 1\n");
        for(int i=0; i<topAtoms_; i++)
        {
            pRefM_[i] = 1.0;
        }
    }

    /* setup plot output file */
    if (!fnRms_.empty() && mpi::isMaster())
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnRms_);
        plotm->setTitle("RMS");
        plotm->setXAxisIsTime();
        plotm->setYFormat(8, 6, 'f');   // y output fmt: "%8.6f"
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
Rms::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_);
    const Selection    &fitsel = pdata->parallelSelection(fitsel_);

    const int                   selatoms = sel.atomCount();
    const int                   fitselatoms = fitsel.atomCount();
    const ConstArrayRef< int >  selAtomIndsArray = sel.atomIndices();
    const ConstArrayRef< int >  fitselAtomIndsArray = fitsel.atomIndices();

    // atom_id *   selind;          // indicies of selected atoms (ref. to topology indexing)
    rvec    *x,      *x_fit;           // coordinates of output selection and fit selection
    rvec    *xp,     *xp_fit;         // coordinates of reference corresponding to selection and fit selection
    real    *w_rms,  *w_fit;       // weights (masses) of output and reference selection atoms

    real      rms_val;            // variable for accumulating data point value (RMSD)


    /* put stuff from selections in format
        so that we can use the legacy routines */

    // // selection indecies
    // (void) snew(selind, selatoms);
    // for(int i=0; i<selatoms; i++)
    // {
    //     selind[i] = (atom_id) selAtomIndsArray.at(i);
    // }

    // coords and masses
    (void) snew(x, selatoms);
    (void) snew(xp, selatoms);
    (void) snew(w_rms, selatoms);
    (void) snew(x_fit, fitselatoms);
    (void) snew(xp_fit, fitselatoms);
    (void) snew(w_fit, fitselatoms);

    // setup output selection arrays
    for (int i=0; i<selatoms; i++)
    {
        const SelectionPosition &spi = sel.position(i);
        // const rvec &spix = spi.x();
        const atom_id &selindi = (atom_id) selAtomIndsArray.at(i);
        // for (int m=0; m<DIM; m++)
        // {
        //     x[i][m]     = spix[m];
        //     xp[i][m]    = pRefX_[selindi][m];
        // }
        (void) copy_rvec(spi.x(), x[i]);
        (void) copy_rvec(pRefX_[selindi], xp[i]);

        if (bRMSUseMassWeights_)
            w_rms[i] = spi.mass();
        else
            w_rms[i] = 1.0;
    }

    // perform fitting of structures
    if (bDoFit_)
    {
        // setup fitting selection arrays
        for (int i=0; i<fitselatoms; i++)
        {
            const SelectionPosition &fitspi = fitsel.position(i);
            // const rvec &spix = spi.x();
            const atom_id &fitselindi = (atom_id) fitselAtomIndsArray.at(i);

            (void) copy_rvec(fitspi.x(), x_fit[i]);
            (void) copy_rvec(pRefX_[fitselindi], xp_fit[i]);

            if (bFitUseMassWeights_)
                w_fit[i] = fitspi.mass();
            else
                w_fit[i] = 1.0;
        }

        /* translate (ie "reset") COM of reference and selected atoms to origin */
        (void) reset_x(selatoms, NULL,      // calc COM from same atoms as in selection
                       selatoms, NULL,
                       xp_fit, w_fit);
        (void) reset_x(selatoms, NULL,      // calc COM from all atoms in selection
                       selatoms, NULL,
                       x_fit, w_fit);

        /* rotate selection for least squares fit */
        (void) do_fit(selatoms, w_fit, xp, x);
    }

    rms_val = rmsdev(selatoms, w_rms, x, xp);   // calc_similar_ind(FALSE, selatoms, NULL, m, x, xp);

    /* write the result */
    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, rms_val);
    dh.finishFrame();

    // release temp frame stuff
    // (void) sfree(selind);
    (void) sfree(x);
    (void) sfree(xp);
    (void) sfree(w_rms);
    (void) sfree(w_fit);
}


void
Rms::finishAnalysis(int /*nframes*/)
{
    sfree(pRefM_);
    sfree(pRefMU_);
}


void
Rms::writeOutput()
{
    fprintf(stderr, "Average RMS: %f\n", avem_->average(0));
    fprintf(stderr, "Std. RMS:   %f\n", avem_->stddev(0));
}

} // namespace analysismodules

} // namespace gmx
