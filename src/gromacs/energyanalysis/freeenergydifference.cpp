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
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/legacyheaders/macros.h"
#include "freeenergydifference.h"

namespace gmx
{

FreeEnergyDifference::FreeEnergyDifference()
{
    refTemp_   = 0;
    fileIndex_ = -1;
}

bool FreeEnergyDifference::initAnalysis(int nre, gmx_enxnm_t enm[])
{
    if ((fileIndex_ < 2) && SimpleEnergy::initAnalysis(nre, enm))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FreeEnergyDifference::addDataSet(std::string name)
{
    if (fileIndex_ < 1)
    {
        fileIndex_++;
        helper()->addDataSet(name);

        return true;
    }
    return false;
}

bool FreeEnergyDifference::addAnalysisFrame(t_enxframe *fr)
{
    if ((fileIndex_ < 0) || (fileIndex_ > 1))
    {
        return false;
    }
    else
    {
        return SimpleEnergy::addAnalysisFrame(fr);
    }
}

bool FreeEnergyDifference::finalizeAnalysis()
{
    const char * ravgleg[] = {
        "\\8D\\4E = E\\sB\\N-E\\sA\\N",
        "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N"
    };
    FILE        *fp;
    double       beta;
    bool         bOK;

    /* check */
    if (eh_[0].nEnergyTerm() != eh_[1].nEnergyTerm())
    {
        fprintf(stderr, "Number of energy terms extracted from files differs\n");
        return false;
    }

    // 1/T Boltzmann constant
    beta = 1.0/(BOLTZ*refTemp_);

    /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
    fp = NULL;
    if (getRunAverFile().size() > 0)
    {
        fp = xvgropen(getRunAverFile().c_str(),
                      "Running average free energy difference",
                      "Time (" unit_time ")", "\\8D\\4E (" unit_energy ")",
                      helper()->getOutputEnvironment());
        xvgr_legend(fp, asize(ravgleg), ravgleg, helper()->getOutputEnvironment());
    }
    fprintf(stdout, "\n%-24s %10s\n",
            "Energy", "dF = -kT ln < exp(-(EB-EA)/kT) >A");

    bOK = true;
    for (EnergyTermIterator et0 = eh_[0].etBegin();
         bOK && (et0 < eh_[0].etEnd()); ++et0)
    {
        EnergyTermIterator et1 = eh_[1].etSearch(et0->getEterm());
        if (eh_[1].etEnd() == et1)
        {
            fprintf(stderr, "Can not find %s in second file.\n",
                    et0->getEterm().c_str());
            bOK = false;
        }
        if (bOK && (et0->nEnergy() != et1->nEnergy()))
        {
            fprintf(stderr, "Not same length of energy files\n");
            bOK = false;
        }
        double sum = 0;
        for (int i = 0; (i < et0->nEnergy()); i++)
        {
            EnergyFrameIterator ef0  = et0->searchEF(i);
            EnergyFrameIterator ef1  = et1->searchEF(i);
            double              dE   = ef0->getE() - ef1->getE();
            sum        += exp(-dE*beta);
            if (NULL != fp)
            {
                fprintf(fp, "%10g %10g %10g\n",
                        ef0->getT(), dE, -BOLTZ*refTemp_*log(sum/(i+1)) );
            }
        }
        double aver = -BOLTZ*refTemp_*log(sum/et0->nEnergy());
        fprintf(stdout, "%-24s %10g\n", et0->getEterm().c_str(), aver);
    }
    if (fp)
    {
        gmx_ffclose(fp);
    }

    return bOK;
}

}
