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
#include "smalloc.h"
#include "rmpbc.h"

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
        "rmsdist calculates RMS deviation of (intra-molecular) pair-wise atom distances relative a reference structure."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o")
                       .description("Plot of RMSD over time.")
                       .filetype(eftPlot).outputFile()
                       .store(&fnRmsDist_).defaultBasename("rmsdist"));

    options->addOption(SelectionOption("select")
                       .description("Selection for RMS calculation.")
                       .required()
                       .valueCount(1)
                       .onlyAtoms()
                       .store(&sel_));

    options->addOption(BooleanOption("mw")
                       .description("Use mass weighted RMS.")
                       .defaultValue(false)
                       .store(&bUseMassWeights_));

    options->addOption(BooleanOption("cache")
                       .description("Put reference distances in memory.")
                       .defaultValue(true)
                       .store(&bDoCache_));

    // topology with coords MUST be provided for use as reference in RMS calculations!
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}

void
RmsDist::optionsFinished(Options * options, TrajectoryAnalysisSettings * /*settings*/)
{
}

void
RmsDist::initAnalysis(const TrajectoryAnalysisSettings &settings,
                       const TopologyInformation &topInfo)
{
    matrix      box;

    pRefTop_ = topInfo.topology();
    (void) topInfo.getTopologyConf(&pRefX_, box);

    // cacheing is not implemented for dynamic selections!
    if (bDoCache_ && sel_.isDynamic())
    {
        if (mpi::isMaster())
        {
            fprintf(stderr, "Cacheing not supported for dynamic selections - option turned off!");
        }
        bDoCache_ = false;
    }

    // NOTE: num of cache elements = selCount * (selCount-1) / 2 .
    if (bDoCache_)
    {
        const int                   selCount            = sel_.atomCount();
        const int                   cacheRows           = selCount - 1;
        const ConstArrayRef< int >  selAtomIndsArray    = sel_.atomIndices();

        rvec        dx;
        t_pbc       pbc;


        t_pbc  *ppbc = settings.hasPBC() ? &pbc : NULL;
        if (ppbc != NULL)
        {
            set_pbc(ppbc, topInfo.ePBC(), box);
        }

        // setup cache structure
        snew(pRefDCache_, cacheRows);
        for (int i=0; i<cacheRows; i++)
        {
            const int       rowLength = selCount - (1+i);
            const atom_id  &selindi = (atom_id) selAtomIndsArray[i];
            const rvec     &xi = pRefX_[selindi];

            // fill cache
            snew(pRefDCache_[i], rowLength);
            for (int j=0; j<rowLength; j++)
            {
                const int       ndxj = i + j + 1;
                const atom_id  &selindj = (atom_id) selAtomIndsArray[ndxj];
                const rvec     &xj = pRefX_[selindj];

                if (ppbc != NULL)
                    pbc_dx(ppbc, xi, xj, dx);
                else
                    rvec_sub(xi, xj, dx);
                pRefDCache_[i][j] = norm(dx);
            }
        }
    }

    // setup plot output file
    if (!fnRmsDist_.empty() && mpi::isMaster())
    {
        AnalysisDataPlotModulePointer plotm(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnRmsDist_);
        plotm->setTitle("RMS dist");
        plotm->setXAxisIsTime();
        plotm->setYFormat(8, 6, 'f');   // gives output fmt: "%8.6f"
        data_.addModule(plotm);
    }

    // add average module
    data_.addModule(avem_);
}

void
RmsDist::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(data_);
    const Selection    &sel = pdata->parallelSelection(sel_);

    const int                   selatoms = sel.atomCount();
    const ConstArrayRef< int >  selAtomIndsArray = sel.atomIndices();

    const real  avging_prefactor = 2.0 / (selatoms * (selatoms-1));
    real        inv_massprefactor = 0.0;    // for mass weight averaging

    rvec        dx, dxp;
    real        r, rp;
    real        rms_val = 0.0;              // for RSMD value


    for (int i=0; i<selatoms-1; i++)
    {
        const SelectionPosition     &spi = sel.position(i);
        const rvec                  &spix = spi.x();
        const real                  &spim = spi.mass();

        const atom_id               &selindi = (atom_id) selAtomIndsArray[i];
        const rvec                  &xpi = pRefX_[selindi];

        for (int j=i+1; j<selatoms; j++)
        {
            const SelectionPosition     &spj = sel.position(j);
            const rvec                  &spjx = spj.x();
            const real                  &spjm = spj.mass();

            // selection i-to-j distance
            if (pbc != NULL)
                pbc_dx(pbc, spix, spjx, dx);
            else
                rvec_sub(spix, spjx, dx);
            r = norm(dx);

            // reference i-to-j distance
            if (bDoCache_)
            {
                rp = pRefDCache_[i][j-(i+1)];
            }
            else
            {
                const atom_id       &selindj = (atom_id) selAtomIndsArray[j];
                const rvec          &xpj = pRefX_[selindj];

                if (pbc != NULL)
                    pbc_dx(pbc, xpi, xpj, dxp);
                else
                    rvec_sub(xpi, xpj, dxp);
                rp = norm(dxp);
            }

            // accumulate squared RMS value
            if (bUseMassWeights_)
            {
                const real m_pf = sqrt(spim * spjm);

                rms_val += m_pf * (r-rp)*(r-rp);
                inv_massprefactor += m_pf;
            }
            else
            {
                rms_val += (r-rp)*(r-rp);
            }

        }
    }

    //  apply correct averaging prefactor
    if (bUseMassWeights_)
    {
        rms_val = sqrt(rms_val / inv_massprefactor);
    }
    else
    {
        rms_val = sqrt(rms_val * avging_prefactor);
    }

    dh.startFrame(frnr, fr.time);
    dh.setPoint(0, rms_val);
    dh.finishFrame();
}


void
RmsDist::finishAnalysis(int /*nframes*/)
{
    // release cache
    if (bDoCache_)
    {
        const int selCount      = sel_.atomCount();
        const int cacheRows     = selCount - 1;

        for (int i=0; i<cacheRows; i++)
        {
            sfree(pRefDCache_[i]);
        }

        sfree(pRefDCache_);
    }
}


void
RmsDist::writeOutput()
{
    fprintf(stderr, "Average RMS: %f\n", avem_->average(0));
    fprintf(stderr, "Std. RMS:   %f\n", avem_->stddev(0));
}

} // namespace analysismodules

} // namespace gmx
