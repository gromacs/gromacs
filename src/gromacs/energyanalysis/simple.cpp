/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <stdio.h>
#include <string.h>

#include "gromacs/energyanalysis/simple.h"
#include "gromacs/math/units.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/energyanalysis/select.h"

namespace gmx
{

SimpleEnergy::SimpleEnergy()
{
    fc_        = NULL;
    fp_        = NULL;
    bSum_      = false;
}

void SimpleEnergy::initOptions(Options *options)
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
    options->setDescription(desc);

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
}

bool SimpleEnergy::initAnalysis(std::vector<std::string> eName,
                                std::vector<std::string> eUnit)
{
    std::vector<int> set;

    if (helper()->nDataSet() == 0)
    {
        return false;
    }
    // Select which energies to use
    select_by_name(eName, set);

    for (EnergyTermIterator eti = helper()->etBegin(); (eti < helper()->etEnd()); ++eti)
    {
        eti->reset();
    }
    for (unsigned int i = 0; (i < set.size()); i++)
    {
        EnergyTerm et(set[i], eName[set[i]], eUnit[set[i]]);
        helper()->addEnergyTerm(et);
    }

    if (getFluctConvFile().size() > 0)
    {
        fc_ = xvgropen(getFluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", helper()->getOutputEnvironment());
    }
    if (!helper()->getStoreData())
    {
        fp_ = xvgropen(getOutputFile().c_str(),
                       "Energy", "Time (ps)",
                       "", helper()->getOutputEnvironment());
        helper()->printXvgLegend(fp_);
    }

    return true;
}

bool SimpleEnergy::addDataSet(std::string name)
{
    helper()->addDataSet(name);

    return true;
}

bool SimpleEnergy::addAnalysisFrame(t_enxframe *fr)
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
        for (EnergyTermIterator eti = helper()->etBegin(); (eti < helper()->etEnd()); ++eti)
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
                        ee /= helper()->getNmol();
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

void SimpleEnergy::sumEnergies()
{
    if (helper()->getStoreData())
    {
        EnergyTerm etsum(0, "Total", "kJ/mol");
        etsum.setStoreData(helper()->getStoreData());

        gmx_int64_t nEner = helper()->etBegin()->nEnergy();

        for (gmx_int64_t i = 0; (i < nEner); i++)
        {
            double sum = 0;
            for (EnergyTermIterator eti = helper()->etBegin(); (eti < helper()->etEnd()); ++eti)
            {
                sum += eti->searchEF(i)->getE();
            }
            // It may be possible to sum the variance and sum of energy terms
            // Please check! Does not look like it, but leaving the comment in
            // anyway. DvdS 2014-05-08.
            etsum.addData(helper()->etBegin()->searchEF(i)->getT(),
                          i, 1, 0.0, 0.0, sum);
        }
        helper()->addEnergyTerm(etsum);
    }
    else
    {
        fprintf(stderr, "WARNING: Energies not stored, so I can not sum them\n");
    }
}

bool SimpleEnergy::finalizeAnalysis()
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
        helper()->printEnergies(getOutputFile());
    }
    helper()->printStatistics(stdout);
    viewOutput();

    return true;
}

void SimpleEnergy::viewOutput()
{
    if (getOutputFile().size() > 0)
    {
        do_view(helper()->getOutputEnvironment(), getOutputFile().c_str(), "-nxy");
    }
    if (getFluctConvFile().size() > 0)
    {
        do_view(helper()->getOutputEnvironment(), getFluctConvFile().c_str(), "-nxy");
    }
}

const char SimpleInfo::name[] = "energy";

const char SimpleInfo::shortDescription[] = "Extract energy from data files and print graphs";

EnergyAnalysisPointer SimpleInfo::create()
{
    return EnergyAnalysisPointer(new SimpleEnergy);
}

}
