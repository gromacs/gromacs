/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * Implements classes in energytermcontainer.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "energytermcontainer.h"

#include <cmath>
#include <cstring>

#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

namespace gmx
{

namespace energyanalysis
{

void EnergyTermContainer::initOptions(IOptionsContainer *options)
{
    options->addOption(IntegerOption("nmol")
                           .store(&nMol_)
                           .description("Number of molecules in the system"));
    options->addOption(IntegerOption("nblocks")
                           .store(&nB_)
                           .description("Number of blocks for error analysis"));
    options->addOption(BooleanOption("dp")
                           .store(&useDouble_)
                           .description("Write output with more digits"));
}

void EnergyTermContainer::setStoreData(bool storeData)
{
    storeData_ = storeData;
    for (auto &eti : et_)
    {
        eti.setStoreData(storeData);
    }
}

void EnergyTermContainer::addFrame(t_enxframe *fr)
{
    for (auto &eti : et_)
    {
        unsigned int findex = eti.fileIndex();
        if (findex < (unsigned int)fr->nre)
        {
            eti.addData(fr->t,
                        fr->step,
                        fr->nsum,
                        fr->ener[findex].esum,
                        fr->ener[findex].eav,
                        fr->ener[findex].e);
        }
    }
}

bool EnergyTermContainer::energyTerm(std::string term, double *e, double *stddev)
{
    for (auto &eti : *this)
    {
        if (eti.energyTerm() == term)
        {
            *e      = eti.average();
            *stddev = eti.standardDeviation();
            return true;
        }
    }
    return false;
}

EnergyTermIterator EnergyTermContainer::etSearch(unsigned int findex)
{
    for (EnergyTermIterator eti = begin(); eti < end(); ++eti)
    {
        if (eti->fileIndex() == findex)
        {
            return eti;
        }
    }
    return end();
}

EnergyTermIterator EnergyTermContainer::etSearch(std::string eTerm)
{
    EnergyTermIterator eti;
    for (eti = begin(); (eti < end()); ++eti)
    {
        if (eti->energyTerm().compare(eTerm) == 0)
        {
            break;
        }
    }
    return eti;
}

bool EnergyTermContainer::energyTerm(unsigned int ftype, double *e, double *stddev)
{
    return energyTerm(interaction_function[ftype].longname, e, stddev);
}

void EnergyTermContainer::printStatistics(FILE *fp)
{
    if (nEnergyTerm() > 0)
    {
        char               buf[256];
        EnergyTermIterator eti = begin();

        fprintf(fp, "\nStatistics over %s steps [ %.4f through %.4f ps ], %u data sets\n",
                gmx_step_str(eti->nSteps(), buf), eti->timeBegin(), eti->timeEnd(), nEnergyTerm());
        unsigned int nb = nBlocks();
        if (nb > 1)
        {
            fprintf(fp, "Error estimate based on averaging over %u blocks of %g ps.\n",
                    nb, eti->timeSpan()/nb);
        }
        else
        {
            fprintf(fp, "Specify number of blocks in order to provide an error estimate.\n");
        }
        fprintf(fp, "%-24s %10s %10s %10s %10s\n",
                "Energy", "Average", "Err.Est.", "RMSD", "Tot-Drift");
        fprintf(fp, "--------------------------------------------------------------------\n");
        for (std::vector<EnergyTerm>::iterator eti = et_.begin(); (eti < et_.end()); ++eti)
        {
            char drift[32];
            char errorEstimate[32];
            if (eti->calculateDrift())
            {
                snprintf(drift, sizeof(drift), "%10g", eti->driftA() * eti->timeSpan());
            }
            else
            {
                snprintf(drift, sizeof(drift), "N/A");
            }
            if (eti->storeData() && (nb > 1))
            {
                snprintf(errorEstimate, sizeof(errorEstimate), "%10g", eti->errorEstimate(nb));
            }
            else
            {
                snprintf(errorEstimate, sizeof(errorEstimate), "N/A");
            }

            fprintf(fp, "%-24s %10g %10s %10g %10s (%s)\n",
                    eti->energyTerm().c_str(),
                    eti->average(),
                    errorEstimate,
                    eti->standardDeviation(),
                    drift,
                    eti->energyUnit().c_str());
        }
    }
    else
    {
        fprintf(fp, "There are no energy terms to be printed.\n");
    }
}

void EnergyTermContainer::printXvgLegend(FILE *fp)
{
    const char **leg;
    char       **kkk;
    // TODO: Cleanup this mess. Needs overhaul of xvgr code.
    snew(leg, nEnergyTerm());
    snew(kkk, nEnergyTerm());
    int          i = 0;
    for (auto &eti : *this)
    {
        kkk[i] = strdup(eti.energyTerm().c_str());
        leg[i] = kkk[i];
        i++;
    }
    xvgr_legend(fp, nEnergyTerm(), leg, outputEnvironment());
    for (int j = 0; (j < i); j++)
    {
        sfree(kkk[j]);
    }
    sfree(kkk);
    sfree(leg);
}

void EnergyTermContainer::printEnergies(std::string outputFile)
{
    gmx_int64_t nEner = begin()->nEnergy();
    if (storeData() && (nEner > 0))
    {
        FILE *fp = xvgropen(outputFile.c_str(), "Energy", "Time (ps)",
                            "Unit", outputEnvironment());
        printXvgLegend(fp);
        for (auto &eti : *this)
        {
            fprintf(fp, "@type xy\n");
            for (gmx_int64_t i = 0; (i < nEner); i++)
            {
                if (doublePrecision())
                {
                    fprintf(fp, "%15e  %15e\n", eti.findFrame(i)->time(), eti.findFrame(i)->energy());
                }
                else
                {
                    fprintf(fp, "%10g  %10g\n", eti.findFrame(i)->time(), eti.findFrame(i)->energy());
                }
            }
            fprintf(fp, "&\n");
        }
        xvgrclose(fp);
    }
    else
    {
        fprintf(stderr, "WARNING: Energies not stored, so I can not print them\n");
    }
}

}

}
