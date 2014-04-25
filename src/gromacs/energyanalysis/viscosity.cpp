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
#include "gromacs/math/vec.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/cstringutil.h"
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

/*! \brief
 * Strings that should appear in the energy file in order to
 * do viscosity analysis
 */
static const char *setnm[] = {
    "Pres-XX",     "Pres-XY", "Pres-XZ",
    "Pres-YX",     "Pres-YY", "Pres-YZ",
    "Pres-ZX",     "Pres-ZY", "Pres-ZZ",
    "Temperature", "Volume",  "Pressure"
};

bool Viscosity::initAnalysis(int nre, gmx_enxnm_t enm[])
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
            if (gmx_strcasecmp(enm[i].name, setnm[j]) == 0)
            {
                EnergyTerm et(i, enm[i].name, enm[i].unit);
                helper()->addEnergyTerm(et);
                bFound[j] = true;
                nfound++;
            }
        }
    }
    helper()->setStoreData(true);
    for (int j = 0; (j < vNR); j++)
    {
        if (!bFound[j])
        {
            fprintf(stderr, "Energy file does not contain %s. Can not compute viscosity.\n",
                    setnm[j]);
        }
    }
    return (nfound == vNR);
}

bool Viscosity::addAnalysisFrame(t_enxframe *fr)
{
    return SimpleEnergy::addAnalysisFrame(fr);
}

bool Viscosity::finalizeAnalysis()
{
    double      eAver[vNR], eStddev[vNR];
    double      Dt;
    gmx_int64_t nsteps, nframes;
    int         nout;

    for (int i = 0; (i < vNR); i++)
    {
        if (!helper()->getEnergyTerm(setnm[i], &eAver[i], &eStddev[i]))
        {
            fprintf(stderr, "Could not extract %s from energy file.\n",
                    setnm[i]);
            return false;
        }
    }
    nsteps  = helper()->etBegin()->nSteps();
    if (nsteps < 2)
    {
        char buf[256];
        fprintf(stderr, "Only %s steps, that's not enough.\n",
                gmx_step_str(nsteps, buf));
    }
    nframes = helper()->etBegin()->nEnergy();
    Dt      = helper()->etBegin()->timeSpan() / (nframes-1);

    // Symmetrise tensor! (and store in first three elements)
    real      **eneset;
    snew(eneset, 4);
    for (int i = 0; i < 4; i++)
    {
        snew(eneset[i], nframes);
    }
    real      **eneint;
    snew(eneint, 3);
    for (int i = 0; i < 3; i++)
    {
        snew(eneint[i], nframes);
    }
    EnergyTermIterator etxyz[vNR];
    for (int i = 0; (i < vNR); i++)
    {
        etxyz[i] = helper()->etSearch(setnm[i]);
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
    }
    for (gmx_int64_t i = 0; (i < nframes-1); i++)
    {
        // Assume all components are store equally often
        int nSum = etxyz[vPresXY]->searchEF(i)->getNsum();
        if (0 >= nSum)
        {
            gmx_fatal(FARGS, "nSum = %d", nSum);
        }
        double fac = 0.5*Dt/nSum;

        eneint[0][i+1] = eneint[0][i] +
            fac*(etxyz[vPresXY]->searchEF(i)->getEsum()+
                 etxyz[vPresYX]->searchEF(i)->getEsum());
        eneint[1][i+1] = eneint[1][i] +
            fac*(etxyz[vPresXZ]->searchEF(i)->getEsum()+
                 etxyz[vPresZX]->searchEF(i)->getEsum());
        eneint[2][i+1] = eneint[2][i] +
            fac*(etxyz[vPresYZ]->searchEF(i)->getEsum()+
                 etxyz[vPresZY]->searchEF(i)->getEsum());
    }
    doEinstein(3, eneint,
               etxyz[vVolume]->average(),
               etxyz[vTemperature]->average());

    /* Do it for shear viscosity */
    const char* leg[] = { "Shear", "Bulk" };
    nout = nframes/2+1;
    low_do_autocorr(NULL, helper()->getOutputEnvironment(), leg[0], nframes, 3,
                    nout, eneset, Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

    /* Now for bulk viscosity */
    low_do_autocorr(NULL, helper()->getOutputEnvironment(), leg[1], nframes, 1,
                    nout, &(eneset[3]), Dt,
                    eacNormal, 1, TRUE, FALSE, FALSE, 0.0, 0.0, 0);

    double factor = (eAver[vVolume]*1e-26/(BOLTZMANN*eAver[vTemperature]))*Dt;

    FILE * fp = NULL;
    if (helper()->getOutputFile().size() > 0)
    {
        fp = xvgropen(helper()->getOutputFile().c_str(), "Viscosity", "Time (ps)",
                      "\\8h\\4 (cp)", helper()->getOutputEnvironment());
        xvgr_legend(fp, asize(leg), leg, helper()->getOutputEnvironment());
    }
    /* Use trapezium rule for integration */
    double integral = 0;
    double intBulk  = 0;
    for (int i = 1; (i < nout); i++)
    {
        integral += 0.5*(eneset[0][i-1]  + eneset[0][i])*factor;
        intBulk  += 0.5*(eneset[3][i-1] + eneset[3][i])*factor;
        if (NULL != fp)
        {
            fprintf(fp, "%10g  %10g  %10g\n", (i*Dt), integral, intBulk);
        }
    }
    for (int i = 0; (i < 4); i++)
    {
        sfree(eneset[i]);
    }
    sfree(eneset);
    for (int i = 0; (i < 3); i++)
    {
        sfree(eneint[i]);
    }
    sfree(eneint);
    if (NULL != fp)
    {
        xvgrclose(fp);
    }

    return true;
}


void Viscosity::doEinstein(int nsets, real **sum,
                           real V, real T)
{
    double      av[4], avold[4];
    double      fac, di;
    int         i, j, m, nf4;

    gmx_int64_t nframes = helper()->etBegin()->nEnergy();
    double      t0      = helper()->etBegin()->timeBegin();
    double      Dt      = helper()->etBegin()->timeSpan() / (nframes-1);

    if (nframes < 1)
    {
        return;
    }

    nf4 = nframes/4+1;

    for (i = 0; i <= nsets; i++)
    {
        avold[i] = 0;
    }
    std::string fni = getFnEinsteinIntegral();
    FILE       *fpi = NULL;
    if (fni.size() > 0)
    {
        fpi = xvgropen(fni.c_str(), "Shear viscosity integral",
                       "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N ps)",
                       helper()->getOutputEnvironment());
    }
    std::string fn = getFnEinstein();
    FILE       *fp = NULL;
    if (fn.size() > 0)
    {
        fp = xvgropen(fn.c_str(), "Shear viscosity using Einstein relation",
                      "Time (ps)", "(kg m\\S-1\\N s\\S-1\\N)",
                      helper()->getOutputEnvironment());
    }
    for (i = 1; i < nf4; i++)
    {
        for (m = 0; m <= nsets; m++)
        {
            av[m] = 0;
        }
        for (j = 0; j < nframes-i; j++)
        {
            for (m = 0; m < nsets; m++)
            {
                di   = sqr(sum[m][j+i]-sum[m][j]);

                av[m]     += di;
                av[nsets] += di/nsets;
            }
        }
        /* Convert to SI for the viscosity */
        fac = (V*NANO*NANO*NANO*PICO*1e10)/(2*BOLTZMANN*T)/(nframes-i);
        if (NULL != fpi)
        {
            fprintf(fpi, "%10g", i*Dt);
        }
        for (m = 0; (m <= nsets); m++)
        {
            av[m] = fac*av[m];
            if (NULL != fpi)
            {
                fprintf(fpi, "  %10g", av[m]);
            }
        }
        if (NULL != fpi)
        {
            fprintf(fpi, "\n");
        }
        if (NULL != fp)
        {
            fprintf(fp, "%10g", (i+0.5)*Dt-t0);
        }
        for (m = 0; (m <= nsets); m++)
        {
            if (NULL != fp)
            {
                fprintf(fp, "  %10g", (av[m]-avold[m])/Dt);
            }
            avold[m] = av[m];
        }
        if (NULL != fp)
        {
            fprintf(fp, "\n");
        }
    }
    if (NULL != fpi)
    {
        xvgrclose(fpi);
    }
    if (NULL != fp)
    {
        xvgrclose(fp);
    }
}

}
