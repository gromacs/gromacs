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
#include "gromacs/legacyheaders/physics.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/xvgr.h"
#include "simple.h"
#include "select.h"

namespace gmx
{

bool SimpleEnergy::initAnalysis(int nre, gmx_enxnm_t enm[])
{
    int nset, *set;

    // Select which energies to use
    set = select_by_name(nre, enm, &nset);

    for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
    {
        eti->reset();
    }
    for (int i = 0; (i < nset); i++)
    {
        EnergyTerm et(set[i], enm[set[i]].name, enm[set[i]].unit);
        addEnergyTerm(et);
    }
    sfree(set);
    if (getFluctConvFile().size() > 0)
    {
        fc_ = xvgropen(getFluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", getOutputEnvironment());
    }
    if (!getStoreData())
    {
        fp_ = xvgropen(getOutputFile().c_str(),
                       "Energy", "Time (ps)",
                       "", getOutputEnvironment());
        printXvgLegend(fp_);
    }

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
        for (EnergyTermIterator eti = etBegin(); (eti < etEnd()); ++eti)
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
                        ee /= getNmol();
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
    return true;
}

void SimpleEnergy::viewOutput()
{
    EnergyAnalysis::viewOutput();
    if (getFluctConvFile().size() > 0)
    {
        do_view(getOutputEnvironment(), getFluctConvFile().c_str(), "-nxy");
    }
}

}
