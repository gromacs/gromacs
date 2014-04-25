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
/*! \internal \file
 * \brief
 * Implements classes in energyanalysis.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#include <math.h>
#include <string.h>
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/statistics/statistics.h"
#include "energyanalysis.h"

namespace gmx
{

void EnergyAnalysis::setStoreData(bool bStoreData)
{
    bStoreData_ = bStoreData;
    for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
    {
        eti->setStoreData(bStoreData);
    }
}

bool EnergyAnalysis::getEnergyTerm(const char *term, double *e, double *stddev)
{
    for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
    {
        if (gmx_strcasecmp(eti->getEterm().c_str(), term) == 0)
        {
            *e      = eti->average();
            *stddev = eti->standardDeviation();
            return true;
        }
    }
    return false;
}

EnergyTermIterator EnergyAnalysis::etSearch(unsigned int findex)
{
    EnergyTermIterator eti;
    for (eti = etBegin(); (eti < etEnd()); ++eti)
    {
        if (eti->getIndex() == findex)
        {
            break;
        }
    }
    return eti;
}

EnergyTermIterator EnergyAnalysis::etSearch(std::string eTerm)
{
    EnergyTermIterator eti;
    for (eti = etBegin(); (eti < etEnd()); ++eti)
    {
        if (eti->getEterm() == eTerm)
        {
            break;
        }
    }
    return eti;
}

bool EnergyAnalysis::getEnergyTerm(unsigned int ftype, double *e, double *stddev)
{
    return getEnergyTerm(interaction_function[ftype].longname, e, stddev);
}

void EnergyAnalysis::printXvgLegend(FILE *fp)
{
    const char **leg;
    snew(leg, nEnergyTerm());
    int          i = 0;
    for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
    {
        leg[i++] = eti->getEterm().c_str();
    }
    xvgr_legend(fp, nEnergyTerm(), leg, getOutputEnvironment());
    sfree(leg);
}

void EnergyAnalysis::printEnergies()
{
    gmx_int64_t nEner = etBegin()->nEnergy();
    if (getStoreData() && (nEner > 0))
    {
        FILE *fp = xvgropen(getOutputFile().c_str(), "Energy", "Time (ps)",
                            "Unit", getOutputEnvironment());
        printXvgLegend(fp);
        for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
        {
            fprintf(fp, "@type xy\n");
            for (gmx_int64_t i = 0; (i < nEner); i++)
            {
                fprintf(fp, "%10g  %10g\n",
                        eti->searchEF(i)->getT(),
                        eti->searchEF(i)->getE());
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

void EnergyAnalysis::printStatistics(FILE *fp)
{
    if (nEnergyTerm() > 0)
    {
        char               buf[256];
        EnergyTermIterator eti = etBegin();

        fprintf(fp, "\nStatistics over %s steps [ %.4f through %.4f ps ], %u data sets\n",
                gmx_step_str(eti->nSteps(), buf), eti->timeBegin(), eti->timeEnd(), nEnergyTerm());
        unsigned int nb = getNblocks();
        fprintf(fp, "Error estimate based on averaging over %u block%s of %g ps.\n",
                nb, (nb > 1) ? "s" : "", eti->timeSpan()/nb);
        fprintf(fp, "%-24s %10s %10s %10s %10s\n",
                "Energy", "Average", "Err.Est.", "RMSD", "Tot-Drift");
        fprintf(fp, "--------------------------------------------------------------------\n");
        for (eti = etBegin(); (eti < etEnd()); ++eti)
        {
            double drift = 0;
            if (eti->calculateDrift())
            {
                drift = eti->driftA() * eti->timeSpan();
            }
            fprintf(fp, "%-24s %10g %10g %10g %10g (%s)\n",
                    eti->getEterm().c_str(),
                    eti->average(),
                    eti->errorEstimate(nb),
                    eti->standardDeviation(),
                    drift,
                    eti->getUnit().c_str());
        }
    }
    else
    {
        fprintf(fp, "There are no energy terms to be printed.\n");
    }
}

void EnergyAnalysis::viewOutput()
{
    do_view(getOutputEnvironment(), getOutputFile().c_str(), "-nxy");
}

}
