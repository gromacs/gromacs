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
#include "gromacs/utility/cstringutil.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/statistics/statistics.h"
#include "energyterm.h"

namespace gmx
{

EnergyTerm::EnergyTerm(unsigned int findex,
                       std::string  eTerm,
                       std::string  eUnit)
{
    reset();
    eTerm_      = eTerm;
    eUnit_      = eUnit;
    findex_     = findex;
    bDrift_     = false;
    bIsEner_    = false;
    bStoreData_ = false;
    bExactData_ = false;
    for (int j = 0; (j <= F_ETOT); j++)
    {
        bIsEner_ = bIsEner_ ||
            (gmx_strcasecmp_min(interaction_function[j].longname, eTerm.c_str()) == 0);
    }
}


void EnergyTerm::reset()
{
    nesum_     = 0;
    energy_    = 0;
    stddev_    = 0;
    findex_    = 0;
    esumTot_   = 0;
    evarTot_   = 0;
    t0_        = NOTSET;
    t1_        = NOTSET;
    step0_     = NOTSET;
    a_         = 0;
    b_         = 0;
    chi2_      = 0;
    R_         = 0;
}

double EnergyTerm::standardDeviation()
{
    return stddev_;
}

double EnergyTerm::average()
{
    return energy_;
}

void EnergyTerm::addData(double t, gmx_int64_t step, int nsum,
                         double esum, double evar, double e)
{
    if (NOTSET == t0_)
    {
        t0_ = t;
    }
    t1_      = t;
    if (NOTSET == step0_)
    {
        step0_ = step;
    }
    step1_   = step;
    if (bStoreData_)
    {
        if (0 == nsum)
        {
            nsum = 1;
        }
        if (0 == evar)
        {
            // We have limited data only!
            esum = nsum*e;
        }
        EnergyFrame ef(t, step, e, nsum, esum, evar);
        ef_.push_back(ef);
    }
    // Equations from Appendix 2.1 in manual
    double m = nesum_;
    double k = nsum;
    evarTot_ += evar;
    if (m > 0)
    {
        evarTot_ += dsqr(esumTot_/m - (esumTot_+esum)/(m+k))*(m*(m+k)/k);
    }
    esumTot_ += esum;
    nesum_   += nsum;
    if (0)
    {
        char buf[256];
        printf("nesum = %s, esumTot = %g, evarTot = %g evar = %g\n",
               gmx_step_str(nesum_, buf), esumTot_, evarTot_, evar);
    }
    // Keep "output" variables up to date.
    energy_ = esumTot_/nesum_;
    stddev_ = sqrt(evarTot_/nesum_);
}

EnergyFrameIterator EnergyTerm::searchEF(gmx_int64_t iframe)
{
    if (getStoreData())
    {
        if ((iframe < nEnergy()) && (iframe >= 0))
        {
            return beginEF() + iframe;
        }
        else if (iframe != nEnergy())
        {
            char buf1[256], buf2[256];
            fprintf(stderr, "WARNING: step %s out of range (0 <= step < %s)\n",
                    gmx_step_str(iframe, buf1), gmx_step_str(nEnergy(), buf2));
        }
    }
    else
    {
        fprintf(stderr, "WARNING: energy frames not stored.\n");
    }
    return endEF();
}

bool EnergyTerm::calculateDrift()
{
    bDrift_ = false;
    if (nEnergy() > 2)
    {
        real *x, *y;

        snew(x, nEnergy());
        snew(y, nEnergy());
        int i = 0;
        for (EnergyFrameIterator efi = beginEF(); (efi < endEF()); ++efi)
        {
            x[i] = efi->getT();
            y[i] = efi->getE();
            i++;
        }
        GMX_RELEASE_ASSERT(i == nEnergy(), "Number of steps in calculateDrift");
        lsq_y_ax_b(i, x, y, &a_, &b_, &R_, &chi2_);
        sfree(y);
        sfree(x);
        bDrift_ = true;
    }
    return bDrift_;
}

bool EnergyTerm::removeDrift()
{
    bool bRemove_ = false;
    if (!bDrift_)
    {
        bRemove_ = calculateDrift();
    }
    fprintf(stderr, "WARNING: removeDrift not implemented yet!\n");
    for (gmx_int64_t i = 0; (i < nEnergy()); i++)
    {

    }
    return bRemove_;
}

double EnergyTerm::errorEstimate(unsigned int nb)
{
    double bSum  = 0;
    double bSum2 = 0;
    for (unsigned int b = 0; (b < nb); b++)
    {
        // There are nb blocks in this analysis
        EnergyFrameIterator f0              = searchEF(b*nEnergy()/nb);
        EnergyFrameIterator f1              = searchEF((b+1)*nEnergy()/nb);
        double              sum             = 0;
        gmx_int64_t         np              = 0;
        for (EnergyFrameIterator f = f0; (f < f1); ++f)
        {
            sum  += f->getEsum();
            np   += f->getNsum();
        }
        sum   /= np;
        bSum  += sum;
        bSum2 += sum*sum;
    }
    double est = sqrt(bSum2/nb - dsqr(bSum/nb));

    return est;
}

}
