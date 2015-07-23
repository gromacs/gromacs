/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * Implements classes in simple.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "simple.h"

#include <stdio.h>
#include <string.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

#include "energyhandler.h"
#include "energyhelper.h"
#include "energyinfo.h"
#include "select_energy.h"

namespace gmx
{

namespace energyanalysis
{

SimpleEnergyModule::SimpleEnergyModule()
{
    fc_        = NULL;
    fp_        = NULL;
    bSum_      = false;
    ehelper_   = new(EnergyHelper);
}

void SimpleEnergyModule::initOptions(Options                           *options,
                                     ICommandLineOptionsModuleSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] extracts energy components",
        "data from an energy file. The user is prompted to interactively",
        "select the desired energy terms.[PAR]",

        "Average, RMSD, and drift are calculated with full precision from the",
        "simulation (see printed manual). Drift is calculated by performing",
        "a least-squares fit of the data to a straight line. The reported total drift",
        "is the difference of the fit at the first and last point.",
        "An error estimate of the average is given based on a block averages",
        "over 5 blocks using the full-precision averages. The error estimate",
        "can be performed over multiple block lengths with the options",
        "[TT]-nbmin[tt] and [TT]-nbmax[tt].",
        "[BB]Note[bb] that in most cases the energy files contains averages over all",
        "MD steps, or over many more points than the number of frames in",
        "energy file. This makes the [THISMODULE] statistics output more accurate",
        "than the [TT].xvg[tt] output. When exact averages are not present in the energy",
        "file, the statistics mentioned above are simply over the single, per-frame",
        "energy values.[PAR]"
    };
    settings->setHelpText(desc);

    // Add option for output files
    setOutputFile("energy.xvg");
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnEnergy_).defaultBasename("energy")
                           .description("Energy terms as a function of time"));
    options->addOption(FileNameOption("convergence")
                           .filetype(eftPlot)
                           .outputFile()
                           .store(&fnFluctConv_)
                           .defaultBasename("fluct_conv")
                           .description("Convergence of fluctuation properties"));
    options->addOption(BooleanOption("sum")
                           .store(&bSum_)
                           .description("Sum the energy terms selected rather than display them all"));
    options->addOption(StringOption("term")
                           .store(&term_)
                           .description("List the energy terms within quotes on the command line rather than interactively, e.g. '3 4' or 'Coulomb Angle'"));
    ehelper_->initOptions(options);
}

bool SimpleEnergyModule::initAnalysis(std::vector<std::string> eName,
                                      std::vector<std::string> eUnit)
{
    std::vector<int> set;

    if (ehelper_->nDataSet() == 0)
    {
        printf("No data sets in initAnalysis\n");
        return false;
    }
    // Select which energies to use
    bool bVerbose = (output_env_get_verbosity(ehelper_->getOutputEnvironment()) > 0);
    select_by_name(eName, set, bVerbose, term_);

    for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
    {
        eti->reset();
    }
    for (unsigned int i = 0; (i < set.size()); i++)
    {
        EnergyTerm et(set[i], ehelper_->getStoreData(),
                      eName[set[i]], eUnit[set[i]]);
        ehelper_->addEnergyTerm(et);
    }

    if (getFluctConvFile().size() > 0)
    {
        fc_ = xvgropen(getFluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", ehelper_->getOutputEnvironment());
    }
    if (getOutputFile().size() > 0)
    {
        fp_ = xvgropen(getOutputFile().c_str(),
                       "Energy", "Time (ps)",
                       "", ehelper_->getOutputEnvironment());
        ehelper_->printXvgLegend(fp_);
    }

    return true;
}

bool SimpleEnergyModule::addDataSet(std::string name)
{
    ehelper_->addDataSet(name);

    return true;
}

bool SimpleEnergyModule::addAnalysisFrame(t_enxframe *fr)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (NULL != fc_)
        {
            fprintf(fc_, "%10g", fr->t);
        }
        if (NULL != fp_)
        {
            fprintf(fp_, "%10g", fr->t);
        }
        for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
        {
            unsigned int findex = eti->getIndex();
            if (findex < (unsigned int)fr->nre)
            {
                eti->addData(fr->t,
                             fr->step,
                             fr->nsum,
                             fr->ener[findex].esum,
                             fr->ener[findex].eav,
                             fr->ener[findex].e);
                if (NULL != fc_)
                {
                    fprintf(fc_, "  %10g", eti->standardDeviation());
                }
                if (NULL != fp_)
                {
                    double ee = fr->ener[findex].e;
                    if (eti->isEner())
                    {
                        ee /= ehelper_->getNmol();
                    }
                    fprintf(fp_, "  %10g", ee);
                }
            }
        }
        if (NULL != fc_)
        {
            fprintf(fc_, "\n");
        }
        if (NULL != fp_)
        {
            fprintf(fp_, "\n");
        }
    }
    return true;
}

void SimpleEnergyModule::sumEnergies()
{
    if (ehelper_->getStoreData())
    {
        EnergyTerm etsum(0, true, "Total", "kJ/mol");
        etsum.setStoreData(ehelper_->getStoreData());

        gmx_int64_t nEner = ehelper_->etBegin()->nEnergy();

        for (gmx_int64_t i = 0; (i < nEner); i++)
        {
            double sum = 0;
            for (EnergyTermIterator eti = ehelper_->etBegin(); (eti < ehelper_->etEnd()); ++eti)
            {
                sum += eti->searchEF(i)->getE();
            }
            // It may be possible to sum the variance and sum of energy terms
            // Please check! Does not look like it, but leaving the comment in
            // anyway. DvdS 2014-05-08.
            etsum.addData(ehelper_->etBegin()->searchEF(i)->getT(),
                          i, 1, 0.0, 0.0, sum);
        }
        ehelper_->addEnergyTerm(etsum);
    }
    else
    {
        fprintf(stderr, "WARNING: Energies not stored, so I can not sum them\n");
    }
}

bool SimpleEnergyModule::finalizeAnalysis()
{
    if (NULL != fc_)
    {
        xvgrclose(fc_);
    }
    if (NULL != fp_)
    {
        xvgrclose(fp_);
    }
    if (getSumming())
    {
        sumEnergies();
        ehelper_->printEnergies(getOutputFile());
    }
    ehelper_->printStatistics(stdout);
    viewOutput();

    return true;
}

void SimpleEnergyModule::viewOutput()
{
    if (getOutputFile().size() > 0)
    {
        do_view(ehelper_->getOutputEnvironment(), getOutputFile().c_str(), "-nxy");
    }
    if (getFluctConvFile().size() > 0)
    {
        do_view(ehelper_->getOutputEnvironment(), getFluctConvFile().c_str(), "-nxy");
    }
}

}   // namespace energyanalysis

const char SimpleInfo::name[] = "energy";

const char SimpleInfo::shortDescription[] = "Extract energy from data files and print graphs";

ICommandLineOptionsModule *SimpleInfo::create()
{
    energyanalysis::EnergyAnalysisModulePointer eamp(new energyanalysis::SimpleEnergyModule);
    return new EnergyInfo(eamp);
}

} // namespace gmx
