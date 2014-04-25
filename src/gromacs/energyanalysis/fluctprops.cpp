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
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/xvgr.h"
#include "fluctprops.h"

namespace gmx
{

//! String describing the product of Enthalpy and Volume
static const char *HVstring = "HV";

FluctProps::FluctProps()
{
    fc_        = NULL;
    iEnthalpy_ = INT_MAX;
    iVolume_   = INT_MAX;
    iTemp_     = INT_MAX;
    iTotal_    = INT_MAX;
    iHV_       = INT_MAX;
};

unsigned int FluctProps::searchTerm(int nre, gmx_enxnm_t enm[], const char *term)
{
    for (int i = 0; (i < nre); i++)
    {
        if (strcmp(term, enm[i].name) == 0)
        {
            EnergyTerm et(i, term, enm[i].unit);
            helper()->addEnergyTerm(et);
            return i;
        }
    }
    return INT_MAX;
}

bool FluctProps::addDataSet(std::string name)
{
    helper()->addDataSet(name);
    return true;
}

bool FluctProps::initAnalysis(int nre, gmx_enxnm_t enm[])
{
    // Select which energies to use
    for (EnergyTermIterator eti = helper()->etBegin(); (eti < helper()->etEnd()); ++eti)
    {
        eti->reset();
    }
    iTemp_     = searchTerm(nre, enm, "Temperature");
    iVolume_   = searchTerm(nre, enm, "Volume");
    iTotal_    = searchTerm(nre, enm, "Total Energy");
    iEnthalpy_ = searchTerm(nre, enm, "Enthalpy");
    iHV_       = searchTerm(nre, enm, HVstring);

    if (((iVolume_ < INT_MAX) && (iEnthalpy_ < INT_MAX)) &&
        (iHV_ == INT_MAX))
    {
        iHV_ = nre;
        EnergyTerm et(iHV_, HVstring, "kJ/mol");
        helper()->addEnergyTerm(et);
    }
    if (getFluctConvFile().size() > 0)
    {
        fc_ = xvgropen(getFluctConvFile().c_str(),
                       "Standard Deviation", "Time (ps)",
                       "", helper()->getOutputEnvironment());
        helper()->printXvgLegend(fc_);
    }
    // We need to store all the energy data for later processing.
    helper()->setStoreData(true);

    return true;
}

bool FluctProps::addAnalysisFrame(t_enxframe *fr)
{
    /* We read a valid frame, so we can use it */
    if (fr->nre > 0)
    {
        if (NULL != fc_)
        {
            fprintf(fc_, "%10g", fr->t);
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
            }
            else if ((iHV_ == findex) && (iEnthalpy_ < INT_MAX) && (iVolume_ < INT_MAX))
            {
                /* Special case for old energy files where HV is not stored. */
                double hv = fr->ener[iEnthalpy_].e * fr->ener[iVolume_].e;
                eti->addData(fr->t, fr->step, 1, 0, 0, hv);
                if (NULL != fc_)
                {
                    fprintf(fc_, "  %10g", eti->standardDeviation());
                }
            }
        }
        if (NULL != fc_)
        {
            fprintf(fc_, "\n");
        }
    }
    return true;
}

bool FluctProps::finalizeAnalysis()
{
    bool   bTemp, bTotal, bEnthalpy, bVolume, bHV;
    double Temp, TempStddev;
    double Total, TotalStddev;
    double Enthalpy, EnthalpyStddev;
    double Volume, VolumeStddev;
    double HV, HVStddev;
    double kappa, cp, dcp, cv;

    /* Compute it all! */
    kappa = cp = dcp = cv = NOTSET;

    if (helper()->getNmol() < 2)
    {
        printf("\nWARNING: #molecules = %d, this may not be what you want.\n",
               helper()->getNmol());
    }
    printf("\nTemperature dependent fluctuation properties.\n");
    printf("\nHeat capacities obtained from fluctuations do *not* include\n");
    printf("quantum corrections. If you want to get a more accurate estimate\n");
    printf("please use the gmx dos tool.\n\n");
    printf("WARNING: Please verify that your simulations are converged and perform\n"
           "a block-averaging error analysis (not implemented in g_energy yet)\n\n");

    /* Temperature */
    if (true == (bTemp = helper()->getEnergyTerm(F_TEMP, &Temp, &TempStddev)))
    {
        printf("Temperature                              = %10g +/- %10g K\n",
               Temp, TempStddev);
    }

    if (true == (bEnthalpy = helper()->getEnergyTerm("Enthalpy", &Enthalpy, &EnthalpyStddev)))
    {
        printf("Enthalpy                                 = %10g +/- %10g kJ/mol\n",
               Enthalpy/helper()->getNmol(), EnthalpyStddev/helper()->getNmol());
        bTotal = false;
    }
    else
    {
        if (true == (bTotal = helper()->getEnergyTerm(F_ETOT, &Total, &TotalStddev)))
        {
            printf("Total energy                         = %10g +/- %10g kJ/mol\n",
                   Total/helper()->getNmol(), TotalStddev/helper()->getNmol());
        }
        else
        {
            fprintf(stderr, "WARNING: Neither Enthalpy nor Total-Energy in your file.\n");
        }
    }

    if (true == (bVolume = helper()->getEnergyTerm("Volume", &Volume, &VolumeStddev)))
    {
        printf("Volume                                   = %10g +/- %10g nm^3/molecule\n",
               Volume/helper()->getNmol(), VolumeStddev/helper()->getNmol());
    }

    bHV = helper()->getEnergyTerm(HVstring, &HV, &HVStddev);

    if (bTemp)
    {
        if (bVolume)
        {
            /* Compute kappa */
            double varv    = dsqr(VolumeStddev);
            kappa   = pow(NANO, 3.0)*(varv/Volume)/(BOLTZMANN*Temp);
            printf("Isothermal Compressibility Kappa         = %10g (1/GPa)\n",
                   kappa/1.0e-9);
            printf("Adiabatic bulk modulus                   = %10g (GPa)\n",
                   1.0e-9/kappa);
        }
        if (bEnthalpy)
        {
            /* Compute cP */
            double varh = dsqr(EnthalpyStddev)/helper()->getNmol();
            cp          = KILO*(varh/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant pressure Cp    = %10g J/mol K\n", cp);
        }
        else if (bTotal)
        {
            /* Total energy.
             * Only compute cv in constant volume runs, which we can test
             * by checking whether the enthalpy was computed.
             */
            double varet = dsqr(TotalStddev);
            cv    = KILO*((varet/helper()->getNmol())/(BOLTZ*Temp*Temp));
            printf("Heat capacity at constant volume Cv      = %10g J/mol K\n", cv);
        }
        /* Alpha, dcp */
        if (bHV)
        {
            double alpha   = ((KILO/AVOGADRO)*(HV-Volume*Enthalpy)/
                              (Volume*BOLTZMANN*Temp*Temp));
            printf("Coefficient of Thermal Expansion Alpha_P = %10g (1/K)\n",
                   alpha);
            dcp     = ((pow(NANO, 3.0)*Volume*AVOGADRO/helper()->getNmol())*Temp*
                       dsqr(alpha)/(kappa));
            printf("Cp-Cv                                    =  %10g J/mol K\n",
                   dcp);
        }

        please_cite(stdout, "Allen1987a");
    }
    else
    {
        fprintf(stderr, "WARNING: No temperature found in your energy file.\n");
    }
    if (NULL != fc_)
    {
        xvgrclose(fc_);
    }
    helper()->printStatistics(stdout);
    helper()->viewOutput();
    return true;
}


}
