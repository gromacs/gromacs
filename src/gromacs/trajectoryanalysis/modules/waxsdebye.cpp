/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Freevolume.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "waxsdebye.h"

#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/atomprop.h"
#include "gromacs/legacyheaders/copyrite.h"

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

const char WaxsDebye::name[]             = "waxsdebye";
const char WaxsDebye::shortDescription[] =
    "Calculate scattering curves";

// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.
WaxsDebye::WaxsDebye()
    : TrajectoryAnalysisModule(name, shortDescription),
      adata_(new AnalysisDataAverageModule())
{
    // We only compute two numbers per frame
    data_.setColumnCount(0, 2);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "waxsdebye");
    
    wdf_ = NULL;
}


WaxsDebye::~WaxsDebye()
{
    // Destroy C structures where there is no automatic memory release
    // C++ takes care of memory in classes (hopefully)
    
    if (NULL != wdf_) 
    {
        delete wdf_;
    }
}


void
WaxsDebye::initOptions(Options                    *options,
                       TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "gmx waxsdebye can compute WAXS/SAXS scattering curves",
        "for each frame in a trajectory, and plot those as",
        "line graphs in an xvg file. It can also compute the deviation",
        "from a reference curve if a reference curve is given.",
        "The code is based on the same data structures that are used",
        "for WAXS/SAXS refinement in mdrun and hence the results",
        "are compatible with the refinement. The value of alpha is computed",
        "using a bisection algorithm for alpha that minimizes the energy,",
        "in effect minimizing the energy with respect to alpha, within",
        "the bounds given in the [TT]tpr[tt] file. Finally the tool can",
        "do some data management, like scaling a reference (e.g. experimental)",
        "data file to match the expected S(0) based on the Debye formula." 
    };

    // Add the descriptive text (program help text) to the options
    options->setDescription(concatenateStrings(desc));

    // Add option for optional input file
    options->addOption(FileNameOption("sfac").filetype(eftGenericData).inputFile()
                       .store(&fnSfactor_).defaultBasename("sq")
                       .description("Structure factors"));

    // Add option for optional input file
    options->addOption(FileNameOption("waxs_ref").filetype(eftGenericData).inputFile()
                       .store(&fnSqref_).defaultBasename("sq")
                       .description("Reference S(q)"));

    options->addOption(FileNameOption("waxs_diff").filetype(eftGenericData).inputFile()
                       .store(&fnSqdiff_).defaultBasename("sq")
                       .description("Difference S(q)"));

    // Add option for optional output file
    options->addOption(FileNameOption("waxs_out").filetype(eftPlot).outputFile()
                       .store(&fnSqcalc_).defaultBasename("sqt")
                       .description("Waxs scattering as function of time"));
    
    // Add option for optional output file
    options->addOption(FileNameOption("was_ener").filetype(eftPlot).outputFile()
                       .store(&fnEner_).defaultBasename("ener")
                       .description("Waxs-Debye energy (see manual)"));

    // Add option for optional output file
    options->addOption(FileNameOption("waxs_alpha").filetype(eftPlot).outputFile()
                       .store(&fnAlpha_).defaultBasename("alpha")
                       .description("Excited state population (see manual)"));

    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select").required().valueCount(1)
                       .store(&sel_)
                       .onlyAtoms());

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);
}

void
WaxsDebye::initAnalysis(const TrajectoryAnalysisSettings &settings,
                        const TopologyInformation        &top)
{
    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);

    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data et.
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(fnSqcalc_);
    plotm_->setTitle("S(q,t)");
    plotm_->setXLabel("q (1/nm)");
    plotm_->setYLabel("S(q) (UNIT)");

    wdf_ = new WaxsDebyeForce(stdout,
                              fnSfactor_.c_str(), 
                              fnSqref_.c_str(),
                              fnSqdiff_.c_str(),
                              fnSqcalc_.c_str(),
                              fnAlpha_.c_str(),
                              NULL, 
                              NULL);
    
    data_.addModule(plotm_);
}

void
WaxsDebye::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                        TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle       dh   = pdata->dataHandle(data_);
    const Selection         &sel  = pdata->parallelSelection(sel_);

    GMX_RELEASE_ASSERT(NULL != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Compute volume and number of insertions to perform
    
    // Magic
    dh.finishFrame();
}


void
WaxsDebye::finishAnalysis(int nframes)
{
    please_cite(stdout, "Bjorling2014a");
}

void
WaxsDebye::writeOutput()
{
    // Final results come from statistics module in analysis framework
}

} // namespace analysismodules

} // namespace gmx
