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
#include "gromacs/legacyheaders/physics.h"
#include "gromacs/legacyheaders/macros.h"
#include "freeenergyestimate.h"

namespace gmx
{

bool FreeEnergyEstimate::initAnalysis(int nre, gmx_enxnm_t enm[])
{
    return SimpleEnergy::initAnalysis(nre, enm);
}

bool FreeEnergyEstimate::addAnalysisFrame(t_enxframe *fr)
{
    return SimpleEnergy::addAnalysisFrame(fr);
}

bool FreeEnergyEstimate::finalizeAnalysis()
{
    const char * ravgleg[] = {
        "\\8D\\4E = E\\sB\\N-E\\sA\\N",
        "<e\\S-\\8D\\4E/kT\\N>\\s0..t\\N"
    };
    FILE        *fp;
    ener_file_t  enx;
    int          timecheck, nenergy, nenergy2, maxenergy;
    int          i, j;
    gmx_bool     bCont;
    real         aver, beta;
    real       **eneset2;
    double       dE, sum;
    gmx_enxnm_t *enm = NULL;
    t_enxframe  *fr;
    char         buf[22];

    /* read second energy file */
    snew(fr, 1);
    enm = NULL;
    enx = open_enx(getEneFile2().c_str(), "r");
    do_enxnms(enx, &(fr->nre), &enm);

    snew(eneset2, nset+1);
    nenergy2  = 0;
    maxenergy = 0;
    timecheck = 0;
    do
    {
        /* This loop searches for the first frame (when -b option is given),
         * or when this has been found it reads just one energy frame
         */
        do
        {
            bCont = do_enx(enx, fr);

            if (bCont)
            {
                timecheck = check_times(fr->t);
            }

        }
        while (bCont && (timecheck < 0));

        /* Store energies for analysis afterwards... */
        if ((timecheck == 0) && bCont)
        {
            if (fr->nre > 0)
            {
                if (nenergy2 >= maxenergy)
                {
                    maxenergy += 1000;
                    for (i = 0; i <= nset; i++)
                    {
                        srenew(eneset2[i], maxenergy);
                    }
                }
                if (fr->t != time[nenergy2])
                {
                    fprintf(stderr, "\nWARNING time mismatch %g!=%g at frame %s\n",
                            fr->t, time[nenergy2], gmx_step_str(fr->step, buf));
                }
                for (i = 0; i < nset; i++)
                {
                    eneset2[i][nenergy2] = fr->ener[set[i]].e;
                }
                nenergy2++;
            }
        }
    }
    while (bCont && (timecheck == 0));

    /* check */
    if (edat->nframes != nenergy2)
    {
        fprintf(stderr, "\nWARNING file length mismatch %d!=%d\n",
                edat->nframes, nenergy2);
    }
    nenergy = std::min(edat->nframes, nenergy2);

    /* calculate fe difference dF = -kT ln < exp(-(E_B-E_A)/kT) >_A */
    fp = NULL;
    if (getRunAverFile().size() > 0)
    {
        fp = xvgropen(getRunAverFile().c_str(),
                      "Running average free energy difference",
                      "Time (" unit_time ")", "\\8D\\4E (" unit_energy ")",
                      getOutputEnvironment());
        xvgr_legend(fp, asize(ravgleg), ravgleg, getOutputEnvironment());
    }
    fprintf(stdout, "\n%-24s %10s\n",
            "Energy", "dF = -kT ln < exp(-(EB-EA)/kT) >A");
    sum  = 0;
    beta = 1.0/(BOLTZ*reftemp_);
    for (i = 0; i < nset; i++)
    {
        if (gmx_strcasecmp(leg[i], enm[set[i]].name) != 0)
        {
            fprintf(stderr, "\nWARNING energy set name mismatch %s!=%s\n",
                    leg[i], enm[set[i]].name);
        }
        for (j = 0; j < nenergy; j++)
        {
            dE   = eneset2[i][j] - edat->s[i].ener[j];
            sum += exp(-dE*beta);
            if (fp)
            {
                fprintf(fp, "%10g %10g %10g\n",
                        time[j], dE, -BOLTZ*reftemp_*log(sum/(j+1)) );
            }
        }
        aver = -BOLTZ*reftemp_*log(sum/nenergy);
        fprintf(stdout, "%-24s %10g\n", leg[i], aver);
    }
    if (fp)
    {
        gmx_ffclose(fp);
    }
    sfree(fr);

    return true;
}

}
