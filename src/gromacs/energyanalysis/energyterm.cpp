/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements classes in energyanalysis.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "energyterm.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <optional>

#include "gromacs/commandline/viewit.h"
#include "gromacs/energyanalysis/energyanalysisframe.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

EnergyTerm::EnergyTerm(unsigned int       indexWithinEnergyFile,
                       bool               bStoreData,
                       const std::string& energyTerm,
                       const std::string& energyUnit) :
    energyTerm_(energyTerm),
    energyUnit_(energyUnit),
    indexWithinEnergyFile_(indexWithinEnergyFile),
    storeData_(bStoreData)
{

    for (int j = 0; (j <= F_ETOT); j++)
    {
        termIsEnergy_ = termIsEnergy_ || equalIgnoreDash(interaction_function[j].longname, energyTerm);
    }
}

void EnergyTerm::addFrame(double  time,
                          int64_t step,
                          int     numIntermediateStepsSum,
                          double  energySumOverNumSteps,
                          double  energyVarianceOverNumSteps,
                          double  energyAtTime)
{
    if (!firstFrameRead_)
    {
        startTime_      = time;
        firstStep_      = step;
        firstFrameRead_ = true;
    }
    else
    {
        endTime_  = time;
        lastStep_ = step;
    }
    if (storeData_)
    {
        if (0 == numIntermediateStepsSum)
        {
            numIntermediateStepsSum = 1;
        }
        if (0 == energyVarianceOverNumSteps)
        {
            // We have limited data only!
            energySumOverNumSteps = numIntermediateStepsSum * energyAtTime;
        }
        energyAnalysisFrames_.emplace_back(
                time, step, energyAtTime, numIntermediateStepsSum, energySumOverNumSteps, energyVarianceOverNumSteps);
    }
    // Equations from Appendix 2.1 in manual
    double m = numberOfEnergyTerms_;
    double k = numIntermediateStepsSum;
    totalVarianceOfEnergy_ += energyVarianceOverNumSteps;
    if (m > 0)
    {
        totalVarianceOfEnergy_ +=
                square(totalSumOfEnergy_ / m - (totalSumOfEnergy_ + energySumOverNumSteps) / (m + k))
                * (m * (m + k) / k);
    }
    totalSumOfEnergy_ += energySumOverNumSteps;
    numberOfEnergyTerms_ += numIntermediateStepsSum;
    // Keep "output" variables up to date.
    if (numberOfEnergyTerms_ > 0)
    {
        average_           = totalSumOfEnergy_ / numberOfEnergyTerms_;
        standardDeviation_ = std::sqrt(totalVarianceOfEnergy_ / numberOfEnergyTerms_);
    }
}

EnergyAnalysisFrameIterator EnergyTerm::findFrame(int64_t frameIndex) const
{
    if (storeData())
    {
        if ((frameIndex < numFrames()) && (frameIndex >= 0))
        {
            return begin() + frameIndex;
        }
        else if (frameIndex != numFrames())
        {
            char buf1[256], buf2[256];
            fprintf(stderr,
                    "WARNING: frame %s out of range (0 <= frame < %s)\n",
                    gmx_step_str(frameIndex, buf1),
                    gmx_step_str(numFrames(), buf2));
        }
    }
    else
    {
        fprintf(stderr, "WARNING: energy frames not stored.\n");
    }
    return end();
}

std::optional<real> EnergyTerm::slopeOfLinearFit() const
{
    if (numFrames() > 2)
    {
        real              a;
        std::vector<real> x, y;

        x.resize(numFrames(), 0);
        y.resize(numFrames(), 0);
        int i = 0;
        for (const auto& efi : *this)
        {
            x[i] = efi.time();
            y[i] = efi.energyAtTime();
            i++;
        }
        GMX_RELEASE_ASSERT(i == numFrames(), "Number of steps in drift() is too large");
        real b, R, chi2;
        lsq_y_ax_b(i, x.data(), y.data(), &a, &b, &R, &chi2);
        return std::optional(a);
    }
    return std::nullopt;
}

std::optional<real> EnergyTerm::errorEstimate(unsigned int numBlocks) const
{
    if (!storeData())
    {
        return std::nullopt;
    }

    double bSum  = 0;
    double bSum2 = 0;
    for (unsigned int b = 0; (b < numBlocks); b++)
    {
        // There are nb blocks in this analysis
        EnergyAnalysisFrameIterator f0  = findFrame(b * numFrames() / numBlocks);
        EnergyAnalysisFrameIterator f1  = findFrame((b + 1) * numFrames() / numBlocks);
        double                      sum = 0;
        int64_t                     np  = 0;
        for (EnergyAnalysisFrameIterator f = f0; (f < f1); ++f)
        {
            sum += f->energySumOverNumSteps();
            np += f->numSteps();
        }
        if (np > 0)
        {
            sum /= np;
            bSum += sum;
            bSum2 += sum * sum;
        }
    }
    return (numBlocks > 0) ? std::optional(std::sqrt(bSum2 / numBlocks - square(bSum / numBlocks)))
                           : std::nullopt;
}

} // namespace gmx
