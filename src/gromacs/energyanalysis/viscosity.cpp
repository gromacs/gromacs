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
#include <string.h>
#include "gromacs/math/vec.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/gmxana/gstat.h"
#include "viscosity.h"

namespace gmx
{

enum {
    vPresXX      = 0, vPresXY = 1, vPresXZ = 2,
    vPresYX      = 3, vPresYY = 4, vPresYZ = 5,
    vPresZX      = 6, vPresZY = 7, vPresZZ = 8,
    vTemperature = 9, vVolume = 10, vPressure = 11,
    vNR          = 12
};

static const char *setnm[] = {
    "Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
    "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ", "Temperature",
    "Volume",  "Pressure"
};

bool Viscosity::readInit(int nre, gmx_enxnm_t enm[])
{
    bool bFound[vNR];
    int  nfound = 0;
    for (int j = 0; (j < vNR); j++)
    {
        bFound[j] = false;
    }
    for (int i = 0; (i < nre); i++)
    {
        for (int j = 0; (j < vNR); j++)
        {
            if (strcasecmp(enm[i].name, setnm[j]) == 0)
            {
                EnergyTerm et(i, enm[i].name, enm[i].unit);
                addEnergyTerm(et);
                bFound[j] = true;
                nfound++;
            }
        }
    }
    setStoreData(true);
    for (int j = 0; (j < vNR); j++)
    {
        if (!bFound[j])
        {
            fprintf(stderr, "Energy file %s does not contain %s. Can not compute viscosity.\n",
                    getEneFile().c_str(), setnm[j]);
        }
    }
    return (nfound == vNR);
}

bool Viscosity::readFrame(t_enxframe *fr)
{
    return SimpleEnergy::readFrame(fr);
}

bool Viscosity::readFinalize()
{
    int         nout;
    char        buf[256];
    double      eAver[vNR], eStddev[vNR];
    double      Dt;
    gmx_int64_t nsteps, nframes;

    for (int i = 0; (i < vNR); i++)
    {
        if (!getEnergyTerm(setnm[i], &eAver[i], &eStddev[i]))
        {
            fprintf(stderr, "Could not extract %s from %s.\n",
                    setnm[i], getEneFile().c_str());
            return false;
        }
    }
    nsteps  = etBegin()->nSteps();
    if (nsteps < 2)
    {
        char buf[256];
        fprintf(stderr, "Only %s steps, that's not enough.\n", gmx_step_str(nsteps, buf));
    }
    nframes = etBegin()->nEnergy();
    Dt      = etBegin()->timeSpan() / (nframes-1);

    // Symmetrise tensor! (and store in first three elements)
    real      **eneset;
    snew(eneset, 4);
    for (int i = 0; i < 4; i++)
    {
        snew(eneset[i], nframes);
    }
    /*
       real      **enesum;
       snew(enesum, 3);
       for (int i = 0; i < 3; i++)
       {
        snew(enesum[i], nframes);
        }*/
    EnergyTermIterator etxyz[vNR];
    for (int i = 0; (i < vNR); i++)
    {
        etxyz[i] = etSearch(setnm[i]);
    }
    for (gmx_int64_t i = 0; (i < nframes); i++)
    {
        eneset[0][i] = 0.5*(etxyz[vPresXY]->searchEF(i)->getE() - eAver[vPresXY] +
                            etxyz[vPresYX]->searchEF(i)->getE() - eAver[vPresYX]);
        eneset[1][i] = 0.5*(etxyz[vPresXZ]->searchEF(i)->getE() - eAver[vPresXZ] +
                            etxyz[vPresZX]->searchEF(i)->getE() - eAver[vPresZX]);
        eneset[2][i] = 0.5*(etxyz[vPresYZ]->searchEF(i)->getE() - eAver[vPresYZ] +
                            etxyz[vPresZY]->searchEF(i)->getE() - eAver[vPresZY]);
        // Subtract average pressure.
        eneset[3][i] = etxyz[vPressure]->searchEF(i)->getE() - eAver[vPressure];
        /*
           enesum[0][i]   = 0.5*(edat->s[vPresXY].es[i].sum+edat->s[vPresYX].es[i].sum);
           enesum[1][i]   = 0.5*(edat->s[vPresXZ].es[i].sum+edat->s[vPresZX].es[i].sum);
           enesum[2][i]   = 0.5*(edat->s[vPresYZ].es[i].sum+edat->s[vPresZY].es[i].sum);
         */
    }
    /*
       einstein_visco("evisco.xvg", "eviscoi.xvg",
                   3, nframes, enesum, Volume, Temp, nsteps, time, oenv);
     */
    /* Do it for shear viscosity */
    const char* leg[] = { "Shear", "Bulk" };
    low_do_autocorr(NULL, getOutputEnvironment(), leg[0], nframes, 3,
                    (nframes+1)/2, eneset, Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

    /* Now for bulk viscosity */
    strcpy(buf, "Bulk Viscosity");
    low_do_autocorr(NULL, getOutputEnvironment(), leg[1], nframes, 1,
                    (nframes+1)/2, &(eneset[3]), Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

    double factor = (eAver[vVolume]*1e-26/(BOLTZMANN*eAver[vTemperature]))*Dt;

    FILE  *fp      = xvgropen(getOutputFile().c_str(), "Viscosity", "Time (ps)",
                              "\\8h\\4 (cp)", getOutputEnvironment());
    xvgr_legend(fp, asize(leg), leg, getOutputEnvironment());

    /* Use trapezium rule for integration */
    double integral = 0;
    double intBulk  = 0;
    nout     = get_acfnout();
    if ((nout < 2) || (nout >= nframes/2))
    {
        nout = nframes/2;
    }
    for (int i = 1; (i < nout); i++)
    {
        integral += 0.5*(eneset[0][i-1]  + eneset[0][i])*factor;
        intBulk  += 0.5*(eneset[3][i-1] + eneset[3][i])*factor;
        fprintf(fp, "%10g  %10g  %10g\n", (i*Dt), integral, intBulk);
    }
    for (int i = 0; (i < 4); i++)
    {
        sfree(eneset[i]);
    }
    sfree(eneset);
    xvgrclose(fp);

    return true;
}


void einstein_visco(const char *fn, const char *fni, int nsets,
                    int nframes, real **sum,
                    real V, real T, int nsteps, double time[],
                    const output_env_t oenv)
{
    FILE  *fp0, *fp1;
    double av[4], avold[4];
    double fac, dt, di;
    int    i, j, m, nf4;

    if (nframes < 1)
    {
        return;
    }

    dt  = (time[1]-time[0]);
    nf4 = nframes/4+1;

    for (i = 0; i <= nsets; i++)
    {
        avold[i] = 0;
    }
    fp0 = xvgropen(fni, "Shear viscosity integral",
                   "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)", oenv);
    fp1 = xvgropen(fn, "Shear viscosity using Einstein relation",
                   "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)", oenv);
    for (i = 1; i < nf4; i++)
    {
        fac = dt*nframes/nsteps;
        for (m = 0; m <= nsets; m++)
        {
            av[m] = 0;
        }
        for (j = 0; j < nframes-i; j++)
        {
            for (m = 0; m < nsets; m++)
            {
                di   = sqr(fac*(sum[m][j+i]-sum[m][j]));

                av[m]     += di;
                av[nsets] += di/nsets;
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        fprintf(fp0, "%10g", time[i]-time[0]);
        for (m = 0; (m <= nsets); m++)
        {
            av[m] = fac*av[m];
            fprintf(fp0, "  %10g", av[m]);
        }
        fprintf(fp0, "\n");
        fprintf(fp1, "%10g", 0.5*(time[i]+time[i-1])-time[0]);
        for (m = 0; (m <= nsets); m++)
        {
            fprintf(fp1, "  %10g", (av[m]-avold[m])/dt);
            avold[m] = av[m];
        }
        fprintf(fp1, "\n");
    }
    gmx_ffclose(fp0);
    gmx_ffclose(fp1);
}

}
