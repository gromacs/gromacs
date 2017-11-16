/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Implements the initialization of the BiasParams class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biasparams.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "grid.h"
#include "math.h"

namespace gmx
{

bool BiasParams::isCheckStep(std::size_t numPointsInHistogram,
                             gmx_int64_t step) const
{
    int numStepsUpdateFreeEnergy = numSamplesUpdateFreeEnergy_*numStepsSampleCoord_;
    int numStepsCheck            = (1 + numPointsInHistogram/numSamplesUpdateFreeEnergy_)*numStepsUpdateFreeEnergy;

    if (step > 0 && step % numStepsCheck == 0)
    {
        GMX_ASSERT(isUpdateFreeEnergyStep(step), "We should only check at free-energy update steps");

        return true;
    }
    else
    {
        return false;
    }
}

namespace
{

/*! \brief Determines the interval for updating the target distribution.
 *
 * The interval value is based on the target distrbution type
 * (this could be made a user-option but there is most likely
 * no big need for tweaking this for most users).
 *
 * \param[in] awhParams      AWH parameters.
 * \param[in] awhBiasParams  Bias parameters.
 * \returns the target update interval in steps.
 */
gmx_int64_t calcTargetUpdateInterval(const AwhParams     &awhParams,
                                     const AwhBiasParams &awhBiasParams)
{
    gmx_int64_t numStepsUpdateTarget = 0;
    /* Set the target update frequency based on the target distrbution type
     * (this could be made a user-option but there is most likely no big need
     * for tweaking this for most users).
     */
    switch (awhBiasParams.eTarget)
    {
        case eawhtargetCONSTANT:
            numStepsUpdateTarget = 0;
            break;
        case eawhtargetCUTOFF:
        case eawhtargetBOLTZMANN:
            /* Updating the target generally requires updating the whole grid so to keep the cost down
               we generally update the target less often than the free energy (unless the free energy
               update step is set to > 100 samples). */
            numStepsUpdateTarget = std::max(100 % awhParams.numSamplesUpdateFreeEnergy,
                                            awhParams.numSamplesUpdateFreeEnergy)*awhParams.nstSampleCoord;
            break;
        case eawhtargetLOCALBOLTZMANN:
            /* The target distribution is set equal to the reference histogram which is updated every free energy update.
               So the target has to be updated at the same time. This leads to a global update each time because it is
               assumed that a target distribution update can take any form. This is a bit unfortunate for a "local"
               target distribution. One could avoid the global update by making a local target update function (and
               postponing target updates for non-local points as for the free energy update). We avoid such additions
               for now and accept that this target type always does global updates. */
            numStepsUpdateTarget = awhParams.numSamplesUpdateFreeEnergy*awhParams.nstSampleCoord;
            break;
        default:
            gmx_incons("Unknown AWH target type");
    }

    return numStepsUpdateTarget;
}

/*! \brief
 * Returns an approximation of the geometry factor used for initializing the AWH update size.
 *
 * The geometry factor is defined as the following sum of Gaussians:
 * sum_{k!=0} exp(-0.5*(k*pi*x)^2)/(pi*k)^2,
 * where k is a xArray.size()-dimensional integer vector with k_i in {0,1,..}.
 *
 * \param[in] xArray  Array to evaluate.
 * \returns the geometry factor.
 */
double gaussianGeometryFactor(gmx::ArrayRef<const double> xArray)
{
    /* For convenience we give the geometry factor function a name: zeta(x) */
    constexpr size_t                    tableSize   = 5;
    std::array<const double, tableSize> xTabulated  =
    { {1e-5, 1e-4, 1e-3, 1e-2, 1e-1} };
    std::array<const double, tableSize> zetaTable1d =
    { { 0.166536811948, 0.16653116886, 0.166250075882, 0.162701098306, 0.129272430287 } };
    std::array<const double, tableSize> zetaTable2d =
    { { 2.31985974274, 1.86307292523, 1.38159772648, 0.897554759158, 0.405578211115 } };

    gmx::ArrayRef<const double> zetaTable;

    if (xArray.size() == 1)
    {
        zetaTable = zetaTable1d;
    }
    else if (xArray.size() == 2)
    {
        zetaTable = zetaTable2d;
    }
    else
    {
        /* TODO... but this is anyway a rough estimate and > 2 dimensions is not so popular. */
        zetaTable = zetaTable2d;
    }

    /* TODO. Really zeta is a function of an ndim-dimensional vector x and we shoudl have a ndim-dimensional lookup-table.
       Here we take the geometric average of the components of x which is ok if the x-components are not very different. */
    double xScalar = 1;
    for (const double &x : xArray)
    {
        xScalar *= x;
    }

    GMX_ASSERT(xArray.size() > 0, "We should have a non-empty input array");
    xScalar = std::pow(xScalar, 1.0/xArray.size());

    /* Look up zeta(x) */
    size_t xIndex = 0;
    while ((xIndex < xTabulated.size()) && (xScalar > xTabulated[xIndex]))
    {
        xIndex++;
    }

    double zEstimate;
    if (xIndex == xTabulated.size())
    {
        /* Take last value */
        zEstimate = zetaTable[xTabulated.size() - 1];
    }
    else if (xIndex == 0)
    {
        zEstimate = zetaTable[xIndex];
    }
    else
    {
        /* Interpolate */
        double x0 = xTabulated[xIndex - 1];
        double x1 = xTabulated[xIndex];
        double w  = (xScalar - x0)/(x1 - x0);
        zEstimate = w*zetaTable[xIndex - 1] + (1 - w)*zetaTable[xIndex];
    }

    return zEstimate;
}

/*! \brief
 * Estimate a reasonable initial reference weight histogram size.
 *
 * \param[in] dimParams         Parameters for the dimensions of the coordinate.
 * \param[in] awhBiasParams     Bias parameters.
 * \param[in] gridAxis          The Grid axes.
 * \param[in] beta              1/(k_B T).
 * \param[in] samplingTimestep  Sampling frequency of probability weights.
 * \returns estimate of initial histogram size.
 */
double getInitialHistSizeEstimate(const std::vector<DimParams> &dimParams,
                                  const AwhBiasParams          &awhBiasParams,
                                  const std::vector<GridAxis>  &gridAxis,
                                  double                        beta,
                                  double                        samplingTimestep)
{
    /* Get diffusion factor */
    double              L2invD = 0.;
    std::vector<double> x;
    for (size_t d = 0; d < gridAxis.size(); d++)
    {
        double L = gridAxis[d].length();
        if (L > 0)
        {
            L2invD   += awhBiasParams.dimParams[d].diffusion/(L*L);
            double s  = 1./std::sqrt(dimParams[d].betak);
            x.push_back(s/L);
        }
    }
    GMX_RELEASE_ASSERT(L2invD > 0, "We need at least one dimension with non-zero length");
    L2invD                  = 1./L2invD;
    double errorInitialInKT = beta*awhBiasParams.errorInitial;
    double histSize         = L2invD*gaussianGeometryFactor(x)/(gmx::square(errorInitialInKT)*samplingTimestep);

    return histSize;
}

/*! \brief Returns the number of simulations sharing bias updates.
 *
 * \param[in] awhBiasParams          Bias parameters.
 * \param[in] numSharingSimulations  The number of simulations to share the bias across.
 * \returns the number of shared updates.
 */
int getNumSharedUpdate(const AwhBiasParams &awhBiasParams,
                       int                  numSharingSimulations)
{
    GMX_RELEASE_ASSERT(numSharingSimulations >= 1, "We should ''share'' at least with ourselves");

    int numShared = 1;

    if (awhBiasParams.shareGroup > 0)
    {
        /* We do not yet support sharing within a simulation */
        int numSharedWithinSim = 1;
        numShared              = numSharingSimulations*numSharedWithinSim;
    }

    return numShared;
}

}   // namespace

/* Constructor.
 *
 * The local Boltzmann target distibution is defined by
 * 1) Adding the sampled weights instead of the target weights to the reference weight histogram.
 * 2) Scaling the weights of these samples by the beta scaling factor.
 * 3) Setting the target distribution equal the reference weight histogram.
 * This requires the following special update settings:
 *   localWeightScaling      = targetParam
 *   idealWeighthistUpdate   = false
 * Note: these variables could in principle be set to something else also for other target distribution types.
 * However, localWeightScaling < 1  is in general expected to give lower efficiency and, except for local Boltzmann,
 * idealWeightHistUpdate = false gives (in my experience) unstable, non-converging results.
 */
BiasParams::BiasParams(const AwhParams              &awhParams,
                       const AwhBiasParams          &awhBiasParams,
                       const std::vector<DimParams> &dimParams,
                       double                        beta,
                       double                        mdTimeStep,
                       DisableUpdateSkips            disableUpdateSkips,
                       int                           numSharingSimulations,
                       const std::vector<GridAxis>  &gridAxis,
                       int                           biasIndex) :
    invBeta(beta > 0 ? 1/beta : 0),
    numStepsSampleCoord_(awhParams.nstSampleCoord),
    numSamplesUpdateFreeEnergy_(awhParams.numSamplesUpdateFreeEnergy),
    numStepsUpdateTarget_(calcTargetUpdateInterval(awhParams, awhBiasParams)),
    eTarget(awhBiasParams.eTarget),
    freeEnergyCutoff(beta*awhBiasParams.targetCutoff),
    temperatureScaleFactor(awhBiasParams.targetBetaScaling),
    idealWeighthistUpdate(eTarget != eawhtargetLOCALBOLTZMANN),
    numSharedUpdate(getNumSharedUpdate(awhBiasParams, numSharingSimulations)),
    updateWeight(numSamplesUpdateFreeEnergy_*numSharedUpdate),
    localWeightScaling(eTarget == eawhtargetLOCALBOLTZMANN ? temperatureScaleFactor : 1),
    // Estimate and initialize histSizeInitial, depends on the grid.
    histSizeInitial(getInitialHistSizeEstimate(dimParams, awhBiasParams, gridAxis, beta, numStepsSampleCoord_*mdTimeStep)),
    convolveForce(awhParams.ePotential == eawhpotentialCONVOLVED),
    biasIndex(biasIndex),
    disableUpdateSkips_(disableUpdateSkips == DisableUpdateSkips::yes)
{
    if (beta <= 0)
    {
        gmx_fatal(FARGS, "To use AWH, the beta=1/(k_B T) should be > 0");
    }

    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        double coverRadiusInNm = 0.5*awhBiasParams.dimParams[d].coverDiameter;
        double spacing         = gridAxis[d].spacing();
        coverRadius_[d]        = spacing > 0 ?  static_cast<int>(std::round(coverRadiusInNm/spacing)) : 0;
    }
}

} // namespace gmx
