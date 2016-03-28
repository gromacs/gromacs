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

#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "grid.h"
#include "math.h"

namespace gmx
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
static int calcTargetUpdateInterval(const AwhParams     &awhParams,
                                    const AwhBiasParams &awhBiasParams)
{
    int nstUpdateTarget = 0;
    /* Set the target update frequency based on the target distrbution type
     * (this could be made a user-option but there is most likely no big need
     * for tweaking this for most users).
     */
    switch (awhBiasParams.eTarget)
    {
        case eawhtargetCONSTANT:
            nstUpdateTarget = 0;
            break;
        case eawhtargetCUTOFF:
        case eawhtargetBOLTZMANN:
            /* Updating the target generally requires updating the whole grid so to keep the cost down
               we generally update the target less often than the free energy (unless the free energy
               update step is set to > 100 samples). */
            nstUpdateTarget = std::max(100 % awhParams.numSamplesUpdateFreeEnergy,
                                       awhParams.numSamplesUpdateFreeEnergy)*awhParams.nstSampleCoord;
            break;
        case eawhtargetLOCALBOLTZMANN:
            /* The target distribution is set equal to the reference histogram which is updated every free energy update.
               So the target has to be updated at the same time. This leads to a global update each time because it is
               assumed that a target distribution update can take any form. This is a bit unfortunate for a "local"
               target distribution. One could avoid the global update by making a local target update function (and
               postponing target updates for non-local points as for the free energy update). We avoid such additions
               for now and accept that this target type always does global updates. */
            nstUpdateTarget = awhParams.numSamplesUpdateFreeEnergy*awhParams.nstSampleCoord;
            break;
        default:
            gmx_incons("Unknown AWH target type");
    }

    return nstUpdateTarget;
}

/*! \brief Returns the target parameter, depending on the target type.
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \returns the target update interval in steps.
 */
static double getTargetParameter(const AwhBiasParams &awhBiasParams)
{
    double targetParam = 0;

    switch (awhBiasParams.eTarget)
    {
        case eawhtargetCUTOFF:
            targetParam = awhBiasParams.targetCutoff;
            break;
        case eawhtargetBOLTZMANN:
        case eawhtargetLOCALBOLTZMANN:
            targetParam = awhBiasParams.targetBetaScaling;
            break;
        default:
            targetParam = 0;
            break;
    }

    return targetParam;
}

/*! \brief
 * Estimate a reasonable initial reference weight histogram size.
 *
 * \param[in] dimParams         Parameters for the dimensions of the coordinate.
 * \param[in] awhBiasParams     Bias parameters.
 * \param[in] gridAxis          The Grid axes.
 * \param[in] samplingTimestep  Sampling frequency of probability weights.
 * \returns estimate of initial histogram size.
 */
static double getInitialHistSizeEstimate(const std::vector<DimParams> &dimParams, const AwhBiasParams &awhBiasParams,
                                         const std::vector<GridAxis> &gridAxis, double samplingTimestep)
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
    L2invD               = 1./L2invD;
    double errorInitial2 = awhBiasParams.errorInitial*awhBiasParams.errorInitial;
    double histSize      = L2invD*Math::gaussianGeometryFactor(x)/(errorInitial2*samplingTimestep);

    return histSize;
}

/*! \brief Returns the number of simulations sharing bias updates.
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \param[in] cr             Struct for communication, can be nullptr.
 * \returns the number of shared updates.
 */
static int getNumSharedUpdate(const AwhBiasParams &awhBiasParams,
                              const t_commrec     *cr)
{
    int numMultiSims = ((cr != nullptr) && MULTISIM(cr)) ? cr->ms->nsim : 1;

    return (awhBiasParams.bShare ? numMultiSims : 1);
}

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
                       const t_commrec              *cr,
                       const std::vector<GridAxis>  &gridAxis,
                       int                           biasIndex) :
    invBeta(beta > 0 ? 1/beta : 0),
    numStepsSampleCoord(awhParams.nstSampleCoord),
    numSamplesUpdateFreeEnergy(awhParams.numSamplesUpdateFreeEnergy),
    nstUpdateTarget(calcTargetUpdateInterval(awhParams, awhBiasParams)),
    eTarget(awhBiasParams.eTarget),
    targetParam(getTargetParameter(awhBiasParams)),
    idealWeighthistUpdate(eTarget != eawhtargetLOCALBOLTZMANN),
    numSharedUpdate(getNumSharedUpdate(awhBiasParams, cr)),
    updateWeight(numSamplesUpdateFreeEnergy*numSharedUpdate),
    localWeightScaling(eTarget == eawhtargetLOCALBOLTZMANN ? targetParam : 1),
    // Estimate and initialize histSizeInitial, depends on the grid.
    histSizeInitial(getInitialHistSizeEstimate(dimParams, awhBiasParams, gridAxis, numStepsSampleCoord*mdTimeStep)),
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
