/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "gmxpre.h"

#include "read_params.h"

#include <algorithm>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/seed.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/stringutil.h"

#include "biasparams.h"
#include "biassharing.h"

namespace gmx
{

const char* enumValueToString(AwhTargetType enumValue)
{
    constexpr gmx::EnumerationArray<AwhTargetType, const char*> awhTargetTypeNames = {
        "constant", "cutoff", "boltzmann", "local-boltzmann"
    };
    return awhTargetTypeNames[enumValue];
}

const char* enumValueToString(AwhHistogramGrowthType enumValue)
{
    constexpr gmx::EnumerationArray<AwhHistogramGrowthType, const char*> awhHistogramGrowthTypeNames = {
        "exp-linear", "linear"
    };
    return awhHistogramGrowthTypeNames[enumValue];
}

const char* enumValueToString(AwhPotentialType enumValue)
{
    constexpr gmx::EnumerationArray<AwhPotentialType, const char*> awhPotentialTypeNames = {
        "convolved", "umbrella"
    };
    return awhPotentialTypeNames[enumValue];
}

const char* enumValueToString(AwhCoordinateProviderType enumValue)
{
    constexpr gmx::EnumerationArray<AwhCoordinateProviderType, const char*> awhCoordinateProviderTypeNames = {
        "pull", "fep-lambda"
    };
    return awhCoordinateProviderTypeNames[enumValue];
}

namespace
{

/*! \brief
 * Check multiple time-stepping consistency between AWH and pull and/or FEP
 *
 * \param[in] inputrec  Inputput parameter struct.
 * \param[in,out] wi    Struct for bookeeping warnings.
 */
void checkMtsConsistency(const t_inputrec& inputrec, WarningHandler* wi)
{
    if (!inputrec.useMts)
    {
        return;
    }

    GMX_RELEASE_ASSERT(inputrec.mtsLevels.size() == 2, "Only 2 MTS levels supported here");

    bool usesPull = false;
    bool usesFep  = false;
    for (const auto& awhBiasParam : inputrec.awhParams->awhBiasParams())
    {
        for (const auto& dimParam : awhBiasParam.dimParams())
        {
            switch (dimParam.coordinateProvider())
            {
                case AwhCoordinateProviderType::Pull: usesPull = true; break;
                case AwhCoordinateProviderType::FreeEnergyLambda: usesFep = true; break;
                default: GMX_RELEASE_ASSERT(false, "Unsupported coord provider");
            }
        }
    }
    const int awhMtsLevel = forceGroupMtsLevel(inputrec.mtsLevels, MtsForceGroups::Awh);
    if (usesPull && forceGroupMtsLevel(inputrec.mtsLevels, MtsForceGroups::Pull) != awhMtsLevel)
    {
        wi->addError(
                "When AWH is applied to pull coordinates, pull and AWH should be computed at "
                "the same MTS level");
    }
    if (usesFep && awhMtsLevel != ssize(inputrec.mtsLevels) - 1)
    {
        wi->addError(
                "When AWH is applied to the free-energy lambda with MTS, AWH should be "
                "computed at the slow MTS level");
    }

    if (inputrec.awhParams->nstSampleCoord() % inputrec.mtsLevels[awhMtsLevel].stepFactor != 0)
    {
        wi->addError("With MTS applied to AWH, awh-nstsample should be a multiple of mts-factor");
    }
}

/*! \brief
 * Check the parameters of an AWH bias pull dimension.
 *
 * \param[in] prefix         Prefix for dimension parameters.
 * \param[in,out] dimParams  AWH dimensional parameters.
 * \param[in] pull_params    Pull parameters.
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
void checkPullDimParams(const std::string&   prefix,
                        const AwhDimParams&  dimParams,
                        const pull_params_t& pull_params,
                        WarningHandler*      wi)
{
    if (dimParams.coordinateIndex() < 0)
    {
        wi->addError(
                gmx::formatString("Failed to read a valid coordinate index for %s-coord-index. "
                                  "Note that the pull coordinate indexing starts at 1.",
                                  prefix.c_str()));
    }
    if (dimParams.coordinateIndex() >= pull_params.ncoord)
    {
        wi->addError(gmx::formatString(
                "The given AWH coordinate index (%d) is larger than the number of pull "
                "coordinates (%d)",
                dimParams.coordinateIndex() + 1,
                pull_params.ncoord));
    }
    if (pull_params.coord[dimParams.coordinateIndex()].rate != 0)
    {
        auto message = formatString(
                "Setting pull-coord%d-rate (%g) is incompatible with AWH biasing this coordinate",
                dimParams.coordinateIndex() + 1,
                pull_params.coord[dimParams.coordinateIndex()].rate);
        wi->addError(message);
    }

    if (gmx_within_tol(dimParams.end() - dimParams.origin(), 0, GMX_REAL_EPS))
    {
        auto message = formatString(
                "The given interval length given by %s-start (%g) and %s-end (%g) is zero. "
                "This will result in only one point along this axis in the coordinate value grid.",
                prefix.c_str(),
                dimParams.origin(),
                prefix.c_str(),
                dimParams.end());
        wi->addWarning(message);
    }

    if (dimParams.forceConstant() <= 0)
    {
        wi->addError("The force AWH bias force constant should be > 0");
    }

    /* Grid params for each axis */
    PullGroupGeometry eGeom = pull_params.coord[dimParams.coordinateIndex()].eGeom;

    /* Check that the requested interval is in allowed range */
    if (eGeom == PullGroupGeometry::Distance)
    {
        if (dimParams.origin() < 0 || dimParams.end() < 0)
        {
            wi->addError(gmx::formatString(
                    "%s-start (%g) or %s-end (%g) set to a negative value. With pull "
                    "geometry distance coordinate values are non-negative. "
                    "Perhaps you want to use geometry %s instead?",
                    prefix.c_str(),
                    dimParams.origin(),
                    prefix.c_str(),
                    dimParams.end(),
                    enumValueToString(PullGroupGeometry::Direction)));
        }
    }
    else if (eGeom == PullGroupGeometry::Angle || eGeom == PullGroupGeometry::AngleAxis)
    {
        if (dimParams.origin() < 0 || dimParams.end() > 180)
        {
            wi->addError(gmx::formatString(
                    "%s-start (%g) and %s-end (%g) are outside of the allowed range "
                    "0 to 180 deg for pull geometries %s and %s ",
                    prefix.c_str(),
                    dimParams.origin(),
                    prefix.c_str(),
                    dimParams.end(),
                    enumValueToString(PullGroupGeometry::Angle),
                    enumValueToString(PullGroupGeometry::AngleAxis)));
        }
    }
    else if (eGeom == PullGroupGeometry::Dihedral)
    {
        if (dimParams.origin() < -180 || dimParams.end() > 180)
        {
            wi->addError(gmx::formatString(
                    "%s-start (%g) and %s-end (%g) are outside of the allowed range "
                    "-180 to 180 deg for pull geometry %s. ",
                    prefix.c_str(),
                    dimParams.origin(),
                    prefix.c_str(),
                    dimParams.end(),
                    enumValueToString(PullGroupGeometry::Dihedral)));
        }
    }
}

/*! \brief
 * Check parameters of an AWH bias in a free energy lambda state dimension.
 *
 * \param[in] prefix         Prefix for dimension parameters.
 * \param[in,out] dimParams  AWH dimensional parameters.
 * \param[in] lambdaParams   The free energy lambda related parameters.
 * \param[in] efep           This is the type of FEP calculation (efep enumerator).
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
void checkFepLambdaDimParams(const std::string&               prefix,
                             const AwhDimParams&              dimParams,
                             const t_lambda*                  lambdaParams,
                             const FreeEnergyPerturbationType efep,
                             WarningHandler*                  wi)
{
    std::string opt;

    if (!lambdaParams)
    {
        wi->addError(
                "There must be free energy input if using AWH to steer the free energy lambda "
                "state.");
        return;
    }

    if (lambdaParams->lambda_neighbors != -1)
    {
        wi->addError(gmx::formatString(
                "When running AWH coupled to the free energy lambda state all lambda states "
                "should be used as neighbors in order to get correct probabilities, i.e. "
                "calc-lambda-neighbors (%d) must be %d.",
                lambdaParams->lambda_neighbors,
                -1));
    }

    if (efep == FreeEnergyPerturbationType::SlowGrowth || lambdaParams->delta_lambda != 0)
    {
        wi->addError(
                "AWH coupled to the free energy lambda state is not compatible with slow-growth "
                "and delta-lambda must be 0.");
    }

    if (efep == FreeEnergyPerturbationType::Expanded)
    {
        wi->addError(
                "AWH is not treated like other expanded ensemble methods. Do not use expanded.");
    }

    if (dimParams.origin() < 0)
    {
        opt = prefix + "-start";
        wi->addError(gmx::formatString(
                "When running AWH coupled to the free energy lambda state the lower lambda state "
                "for AWH, %s (%.0f), must be >= 0.",
                opt.c_str(),
                dimParams.origin()));
    }
    if (dimParams.end() >= lambdaParams->n_lambda)
    {
        opt = prefix + "-end";
        wi->addError(gmx::formatString(
                "When running AWH coupled to the free energy lambda state the upper lambda state "
                "for AWH, %s (%.0f), must be < n_lambda (%d).",
                opt.c_str(),
                dimParams.origin(),
                lambdaParams->n_lambda));
    }
    if (gmx_within_tol(dimParams.end() - dimParams.origin(), 0, GMX_REAL_EPS))
    {
        auto message = formatString(
                "The given interval length given by %s-start (%g) and %s-end (%g) is zero. "
                "This will result in only one lambda point along this free energy lambda state "
                "axis in the coordinate value grid.",
                prefix.c_str(),
                dimParams.origin(),
                prefix.c_str(),
                dimParams.end());
        wi->addWarning(message);
    }

    if (dimParams.forceConstant() != 0)
    {
        wi->addError(

                "The force AWH bias force constant is not used with free energy lambda state as "
                "coordinate provider.");
    }

    if (!fepLambdasChangeAtSameRate(lambdaParams->all_lambda))
    {
        wi->addWarning(
                "Some free-energy lambda components change at different rates over some lambda "
                "interval. The friction metric will be incorrect for such intervals.");
    }
}

/*! \brief
 * Check that AWH FEP is not combined with incompatible decoupling.
 *
 * \param[in] mtop      System topology.
 * \param[in,out] wi    Struct for bookeeping warnings.
 */
void checkFepLambdaDimDecouplingConsistency(const gmx_mtop_t& mtop, WarningHandler* wi)
{
    if (haveFepPerturbedMasses(mtop))
    {
        wi->addWarning(
                "Masses may not be perturbed when using the free energy lambda state as AWH "
                "coordinate provider. If you are using fep-lambdas to specify lambda states "
                "make sure that you also specify mass-lambdas without perturbation.");
    }
    if (havePerturbedConstraints(mtop))
    {
        wi->addWarning(
                "Constraints may not be perturbed when using the free energy lambda state as AWH "
                "coordinate provider. If you are using fep-lambdas to specify lambda states "
                "make sure that you also specify mass-lambdas without perturbation.");
    }
}

/*! \brief
 * Check the parameters of an AWH bias dimension.
 *
 * \param[in] prefix         Prefix for dimension parameters.
 * \param[in,out] dimParams  AWH dimensional parameters.
 * \param[in] ir             Input parameter struct.
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
void checkDimParams(const std::string& prefix, const AwhDimParams& dimParams, const t_inputrec& ir, WarningHandler* wi)
{
    if (dimParams.coordinateProvider() == AwhCoordinateProviderType::Pull)
    {
        if (!ir.bPull)
        {
            wi->addError(
                    "AWH biasing along a pull dimension is only compatible with COM pulling "
                    "turned on");
        }
        checkPullDimParams(prefix, dimParams, *ir.pull, wi);
    }
    else if (dimParams.coordinateProvider() == AwhCoordinateProviderType::FreeEnergyLambda)
    {
        if (ir.efep == FreeEnergyPerturbationType::No)
        {
            wi->addError(
                    "AWH biasing along a free energy lambda state dimension is only compatible "
                    "with free energy turned on");
        }
        checkFepLambdaDimParams(prefix, dimParams, ir.fepvals.get(), ir.efep, wi);
    }
    else
    {
        wi->addError(
                "AWH biasing can only be  applied to pull and free energy lambda state "
                "coordinates");
    }
}

/*! \brief
 * Check the parameters of an AWH bias.
 *
 * \param[in]     awhBiasParams  AWH dimensional parameters.
 * \param[in]     prefix         Prefix for bias parameters.
 * \param[in]     ir             Input parameter struct.
 * \param[in,out] wi             Struct for bookeeping warnings.
 */
void checkBiasParams(const AwhBiasParams& awhBiasParams,
                     const std::string&   prefix,
                     const t_inputrec&    ir,
                     WarningHandler*      wi)
{
    std::string opt = prefix + "-error-init";
    if (awhBiasParams.initialErrorEstimate() <= 0)
    {
        wi->addError(gmx::formatString("%s needs to be > 0.", opt.c_str()));
    }

    if (awhBiasParams.growthFactor() <= 1)
    {
        opt = prefix + "-growth-factor";
        wi->addError(gmx::formatString("%s needs to be > 1.", opt.c_str()));
    }

    opt = prefix + "-equilibrate-histogram";
    if (awhBiasParams.equilibrateHistogram()
        && awhBiasParams.growthType() != AwhHistogramGrowthType::ExponentialLinear)
    {
        auto message =
                formatString("Option %s will only have an effect for histogram growth type '%s'.",
                             opt.c_str(),
                             enumValueToString(AwhHistogramGrowthType::ExponentialLinear));
        wi->addWarning(message);
    }

    if ((awhBiasParams.targetDistribution() == AwhTargetType::LocalBoltzmann)
        && (awhBiasParams.growthType() == AwhHistogramGrowthType::ExponentialLinear))
    {
        auto message = formatString(
                "Target type '%s' combined with histogram growth type '%s' is not "
                "expected to give stable bias updates. You probably want to use growth type "
                "'%s' instead.",
                enumValueToString(AwhTargetType::LocalBoltzmann),
                enumValueToString(AwhHistogramGrowthType::ExponentialLinear),
                enumValueToString(AwhHistogramGrowthType::Linear));
        wi->addWarning(message);
    }

    opt = prefix + "-target-beta-scaling";
    switch (awhBiasParams.targetDistribution())
    {
        case AwhTargetType::Boltzmann:
        case AwhTargetType::LocalBoltzmann:
            if (awhBiasParams.targetBetaScaling() < 0 || awhBiasParams.targetBetaScaling() > 1)
            {
                wi->addError(gmx::formatString("%s = %g is not useful for target type %s.",
                                               opt.c_str(),
                                               awhBiasParams.targetBetaScaling(),
                                               enumValueToString(awhBiasParams.targetDistribution())));
            }
            break;
        default:
            if (awhBiasParams.targetBetaScaling() != 0)
            {
                wi->addError(gmx::formatString(
                        "Value for %s (%g) set explicitly but will not be used for target type %s.",
                        opt.c_str(),
                        awhBiasParams.targetBetaScaling(),
                        enumValueToString(awhBiasParams.targetDistribution())));
            }
            break;
    }

    opt = prefix + "-target-cutoff";
    switch (awhBiasParams.targetDistribution())
    {
        case AwhTargetType::Cutoff:
            if (awhBiasParams.targetCutoff() <= 0)
            {
                wi->addError(gmx::formatString("%s = %g is not useful for target type %s.",
                                               opt.c_str(),
                                               awhBiasParams.targetCutoff(),
                                               enumValueToString(awhBiasParams.targetDistribution())));
            }
            break;
        default:
            if (awhBiasParams.targetCutoff() != 0)
            {
                wi->addError(gmx::formatString(

                        "Value for %s (%g) set explicitly but will not be used for target type %s.",
                        opt.c_str(),
                        awhBiasParams.targetCutoff(),
                        enumValueToString(awhBiasParams.targetDistribution())));
            }
            break;
    }

    opt = prefix + "-share-group";
    if (awhBiasParams.shareGroup() < 0)
    {
        wi->addError("AWH bias share-group should be >= 0");
    }

    opt = prefix + "-ndim";
    if (awhBiasParams.ndim() <= 0 || awhBiasParams.ndim() > c_biasMaxNumDim)
    {
        wi->addError(gmx::formatString(
                "%s (%d) needs to be > 0 and at most %d\n", opt.c_str(), awhBiasParams.ndim(), c_biasMaxNumDim));
    }
    if (awhBiasParams.ndim() > 2)
    {
        wi->addNote(
                "For awh-dim > 2 the estimate based on the diffusion and the initial error is "
                "currently only a rough guideline."
                " You should verify its usefulness for your system before production runs!");
    }
    for (int d = 0; d < awhBiasParams.ndim(); d++)
    {
        std::string prefixdim = prefix + formatString("-dim%d", d + 1);
        checkDimParams(prefixdim, awhBiasParams.dimParams(d), ir, wi);
    }
}

/*! \brief
 * Check consistency of input at the AWH bias level.
 *
 * \param[in]     awhBiasParams  AWH bias parameters.
 * \param[in,out] wi             Struct for bookkeeping warnings.
 */
void checkInputConsistencyAwhBias(const AwhBiasParams& awhBiasParams, WarningHandler* wi)
{
    /* Covering diameter and sharing warning. */
    const auto awhBiasDimensionParams = awhBiasParams.dimParams();
    for (const auto& dimensionParam : awhBiasDimensionParams)
    {
        double coverDiameter = dimensionParam.coverDiameter();
        if (awhBiasParams.shareGroup() <= 0 && coverDiameter > 0)
        {
            wi->addWarning(
                    "The AWH covering diameter is only relevant to set for bias sharing "
                    "simulations.");
        }
        if (awhBiasParams.shareGroup() > 0 && coverDiameter == 0)
        {
            wi->addWarning(
                    "When simulations share an AWH bias, it is strongly recommended to use a "
                    "non-zero covering diameter");
        }
    }
}

/*! \brief
 * Check consistency of input at the AWH level.
 *
 * \param[in]     awhParams  AWH parameters.
 * \param[in,out] wi         Struct for bookkeeping warnings.
 */
void checkInputConsistencyAwh(const AwhParams& awhParams, WarningHandler* wi)
{
    /* Each pull coord can map to at most 1 AWH coord.
     * Check that we have a shared bias when requesting multisim sharing.
     */
    bool       haveSharedBias = false;
    const auto awhBiasParams  = awhParams.awhBiasParams();
    for (int k1 = 0; k1 < awhParams.numBias(); k1++)
    {
        const AwhBiasParams& awhBiasParams1 = awhBiasParams[k1];

        if (awhBiasParams1.shareGroup() > 0)
        {
            haveSharedBias = true;
        }

        /* k1 is the reference AWH, k2 is the AWH we compare with (can be equal to k1) */
        for (int k2 = k1; k2 < awhParams.numBias(); k2++)
        {
            const AwhBiasParams& awhBiasParams2 = awhBiasParams[k2];
            const auto&          dimParams1     = awhBiasParams1.dimParams();
            const auto&          dimParams2     = awhBiasParams2.dimParams();
            for (int d1 = 0; d1 < gmx::ssize(dimParams1); d1++)
            {
                if (dimParams1[d1].coordinateProvider() == AwhCoordinateProviderType::FreeEnergyLambda)
                {
                    continue;
                }
                /* d1 is the reference dimension of the reference AWH. d2 is the dim index of the AWH to compare with. */
                for (int d2 = 0; d2 < gmx::ssize(dimParams2); d2++)
                {
                    if (dimParams2[d2].coordinateProvider() == AwhCoordinateProviderType::FreeEnergyLambda)
                    {
                        continue;
                    }
                    /* Give an error if (d1, k1) is different from (d2, k2) but the pull coordinate is the same */
                    if ((d1 != d2 || k1 != k2)
                        && (dimParams1[d1].coordinateIndex() == dimParams2[d2].coordinateIndex()))
                    {
                        wi->addError(gmx::formatString(
                                "One pull coordinate (%d) cannot be mapped to two separate AWH "
                                "dimensions (awh%d-dim%d and awh%d-dim%d). "
                                "If this is really what you want to do you will have to duplicate "
                                "this pull coordinate.",
                                dimParams1[d1].coordinateIndex() + 1,
                                k1 + 1,
                                d1 + 1,
                                k2 + 1,
                                d2 + 1));
                    }
                }
            }
        }
    }

    if (awhParams.shareBiasMultisim() && !haveSharedBias)
    {
        wi->addWarning(
                "Sharing of biases over multiple simulations is requested, but no bias is marked "
                "as shared (share-group > 0)");
    }

    /* mdrun does not support this (yet), but will check again */
    if (haveBiasSharingWithinSimulation(awhParams))
    {
        wi->addWarning(
                "You have shared biases within a single simulation, but mdrun does not support "
                "this (yet)");
    }
}

} // namespace

AwhDimParams::AwhDimParams(std::vector<t_inpfile>* inp, const std::string& prefix, WarningHandler* wi, bool bComment)
{
    std::string opt;
    if (bComment)
    {
        printStringNoNewline(
                inp,
                "The provider of the reaction coordinate, "
                "currently only 'pull' and 'fep-lambda' (free energy lambda state) is supported");
    }

    opt             = prefix + "-coord-provider";
    eCoordProvider_ = getEnum<AwhCoordinateProviderType>(inp, opt.c_str(), wi);

    if (bComment)
    {
        printStringNoNewline(inp, "The coordinate index for this dimension");
    }
    opt = prefix + "-coord-index";
    int coordIndexInput;
    coordIndexInput = get_eint(inp, opt, 1, wi);
    if (coordIndexInput < 1)
    {
        wi->addError(
                gmx::formatString("Failed to read a valid coordinate index for %s. "
                                  "Note that the pull coordinate indexing starts at 1.",
                                  opt.c_str()));
    }

    /* The pull coordinate indices start at 1 in the input file, at 0 internally */
    coordIndex_ = coordIndexInput - 1;

    if (bComment)
    {
        printStringNoNewline(inp, "Start and end values for each coordinate dimension");
    }

    opt     = prefix + "-start";
    origin_ = get_ereal(inp, opt, 0., wi);

    opt  = prefix + "-end";
    end_ = get_ereal(inp, opt, 0., wi);

    if (bComment)
    {
        printStringNoNewline(
                inp, "The force constant for this coordinate (kJ/mol/nm^2 or kJ/mol/rad^2)");
    }
    opt            = prefix + "-force-constant";
    forceConstant_ = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Estimated diffusion constant (nm^2/ps or rad^2/ps or ps^-1)");
    }
    opt                   = prefix + "-diffusion";
    double diffusionValue = get_ereal(inp, opt, 0, wi);
    if (diffusionValue <= 0)
    {
        const double diffusion_default = 1e-5;
        auto         message           = formatString(
                "%s not explicitly set by user. You can choose to use a default "
                "value (%g nm^2/ps or rad^2/ps) but this may very well be "
                "non-optimal for your system!",
                opt.c_str(),
                diffusion_default);
        wi->addNote(message);
        diffusionValue = diffusion_default;
    }
    diffusion_ = diffusionValue;

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Diameter that needs to be sampled around a point before it is "
                             "considered covered. In FEP dimensions the cover diameter is "
                             "specified in lambda states.");
    }
    opt            = prefix + "-cover-diameter";
    coverDiameter_ = get_ereal(inp, opt, 0, wi);
}

AwhDimParams::AwhDimParams(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(serializer->reading(),
                       "Can not use writing serializer for creating datastructure");
    serializer->doEnumAsInt(&eCoordProvider_);
    serializer->doInt(&coordIndex_);
    serializer->doDouble(&origin_);
    serializer->doDouble(&end_);
    serializer->doDouble(&period_);
    serializer->doDouble(&forceConstant_);
    serializer->doDouble(&diffusion_);
    serializer->doDouble(&coordValueInit_);
    serializer->doDouble(&coverDiameter_);
}

void AwhDimParams::serialize(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Can not use reading serializer for writing datastructure");
    serializer->doEnumAsInt(&eCoordProvider_);
    serializer->doInt(&coordIndex_);
    serializer->doDouble(&origin_);
    serializer->doDouble(&end_);
    serializer->doDouble(&period_);
    serializer->doDouble(&forceConstant_);
    serializer->doDouble(&diffusion_);
    serializer->doDouble(&coordValueInit_);
    serializer->doDouble(&coverDiameter_);
}

AwhBiasParams::AwhBiasParams(std::vector<t_inpfile>* inp, const std::string& prefix, WarningHandler* wi, bool bComment)
{
    if (bComment)
    {
        printStringNoNewline(inp, "Estimated initial PMF error (kJ/mol)");
    }

    std::string opt = prefix + "-error-init";
    errorInitial_   = get_ereal(inp, opt, 10, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Growth rate of the reference histogram determining the bias update "
                             "size: exp-linear or linear");
    }
    opt      = prefix + "-growth";
    eGrowth_ = getEnum<AwhHistogramGrowthType>(inp, opt.c_str(), wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Growth factor during the exponential growth phase");
    }
    opt           = prefix + "-growth-factor";
    growthFactor_ = get_ereal(inp, opt, 2, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Start the simulation by equilibrating histogram towards the target "
                             "distribution: no or yes");
    }
    opt                   = prefix + "-equilibrate-histogram";
    equilibrateHistogram_ = (getEnum<Boolean>(inp, opt.c_str(), wi) != Boolean::No);

    if (bComment)
    {
        printStringNoNewline(
                inp, "Target distribution type: constant, cutoff, boltzmann or local-boltzmann");
    }
    opt      = prefix + "-target";
    eTarget_ = getEnum<AwhTargetType>(inp, opt.c_str(), wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Boltzmann beta scaling factor for target distribution types "
                             "'boltzmann' and 'boltzmann-local'");
    }
    opt                = prefix + "-target-beta-scaling";
    targetBetaScaling_ = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Free energy cutoff value for target distribution type 'cutoff'");
    }
    opt           = prefix + "-target-cutoff";
    targetCutoff_ = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Initialize PMF and target with user data: no or yes");
    }
    opt        = prefix + "-user-data";
    bUserData_ = getEnum<Boolean>(inp, opt.c_str(), wi) != Boolean::No;

    if (bComment)
    {
        printStringNoNewline(inp, "Group index to share the bias with, 0 means not shared");
    }
    opt         = prefix + "-share-group";
    shareGroup_ = get_eint(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Scale the target distribution (can be used to modify any target "
                             "distribution type and can be combined with user data) based on "
                             "the AWH friction metric: no or yes");
    }
    opt                  = prefix + "-target-metric-scaling";
    scaleTargetByMetric_ = getEnum<Boolean>(inp, opt.c_str(), wi) != Boolean::No;

    if (scaleTargetByMetric_
        && (eTarget_ == AwhTargetType::Boltzmann || eTarget_ == AwhTargetType::LocalBoltzmann))
    {
        auto message = formatString(
                "Combining a %s target distribution with scaling the target distribution "
                "by the friction metric (%s) might result in a feedback loop between the two "
                "adaptive update mechanisms.",
                enumValueToString(eTarget_),
                opt.c_str());
        wi->addWarning(message);
    }

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Maximum factor when scaling the target distribution based on the "
                             "friction metric. The inverse of the value is used as the lower "
                             "limit for the scaling.");
    }
    opt                             = prefix + "-target-metric-scaling-limit";
    double targetMetricScalingLimit = get_ereal(inp, opt, 10, wi);
    if (scaleTargetByMetric_ && targetMetricScalingLimit <= 1)
    {
        auto message = formatString(
                "%s (%g) must be > 1. Setting it to the default value 10. This may "
                "not be optimal for your system.",
                opt.c_str(),
                targetMetricScalingLimit);
        wi->addWarning(message);
        targetMetricScalingLimit = 10;
    }
    targetMetricScalingLimit_ = targetMetricScalingLimit;

    if (bComment)
    {
        printStringNoNewline(inp, "Dimensionality of the coordinate");
    }
    opt      = prefix + "-ndim";
    int ndim = get_eint(inp, opt, 0, wi);

    /* Check this before starting to read the AWH dimension parameters. */
    if (ndim <= 0 || ndim > c_biasMaxNumDim)
    {
        gmx_fatal(FARGS, "%s (%d) needs to be > 0 and at most %d\n", opt.c_str(), ndim, c_biasMaxNumDim);
    }
    for (int d = 0; d < ndim; d++)
    {
        bComment              = bComment && d == 0;
        std::string prefixdim = prefix + formatString("-dim%d", d + 1);
        dimParams_.emplace_back(inp, prefixdim, wi, bComment);
    }
}

AwhBiasParams::AwhBiasParams(ISerializer* serializer,
                             const bool   tprWithoutGrowthFactor,
                             const bool   tprWithoutTargetMetricScaling)
{
    GMX_RELEASE_ASSERT(serializer->reading(),
                       "Can not use writing serializer to create datastructure");
    serializer->doEnumAsInt(&eTarget_);
    serializer->doDouble(&targetBetaScaling_);
    serializer->doDouble(&targetCutoff_);
    serializer->doEnumAsInt(&eGrowth_);
    if (tprWithoutGrowthFactor)
    {
        // A factor 3 was the old, hard-coded value
        growthFactor_ = 3;
    }
    else
    {
        serializer->doDouble(&growthFactor_);
    }
    int temp = 0;
    serializer->doInt(&temp);
    bUserData_ = static_cast<bool>(temp);
    if (tprWithoutTargetMetricScaling)
    {
        scaleTargetByMetric_      = false;
        targetMetricScalingLimit_ = 10;
    }
    else
    {
        serializer->doBool(&scaleTargetByMetric_);
        serializer->doDouble(&targetMetricScalingLimit_);
    }
    serializer->doDouble(&errorInitial_);
    int numDimensions = dimParams_.size();
    serializer->doInt(&numDimensions);
    serializer->doInt(&shareGroup_);
    serializer->doBool(&equilibrateHistogram_);

    for (int k = 0; k < numDimensions; k++)
    {
        dimParams_.emplace_back(serializer);
    }
}

void AwhBiasParams::serialize(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Can not use reading serializer to write datastructure");
    serializer->doEnumAsInt(&eTarget_);
    serializer->doDouble(&targetBetaScaling_);
    serializer->doDouble(&targetCutoff_);
    serializer->doEnumAsInt(&eGrowth_);
    serializer->doDouble(&growthFactor_);
    int temp = static_cast<int>(bUserData_);
    serializer->doInt(&temp);
    serializer->doBool(&scaleTargetByMetric_);
    serializer->doDouble(&targetMetricScalingLimit_);
    serializer->doDouble(&errorInitial_);
    int numDimensions = ndim();
    serializer->doInt(&numDimensions);
    serializer->doInt(&shareGroup_);
    serializer->doBool(&equilibrateHistogram_);

    for (int k = 0; k < numDimensions; k++)
    {
        dimParams_[k].serialize(serializer);
    }
}

AwhParams::AwhParams(std::vector<t_inpfile>* inp, WarningHandler* wi)
{
    std::string opt;

    /* Parameters common for all biases */

    printStringNoNewline(inp, "The way to apply the biasing potential: convolved or umbrella");
    opt            = "awh-potential";
    potentialEnum_ = getEnum<AwhPotentialType>(inp, opt.c_str(), wi);

    printStringNoNewline(inp,
                         "The random seed used for sampling the umbrella center in the case of "
                         "umbrella type potential");
    opt   = "awh-seed";
    seed_ = get_eint(inp, opt, -1, wi);
    if (seed_ == -1)
    {
        seed_ = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Setting the AWH bias MC random seed to %" PRId64 "\n", seed_);
    }

    printStringNoNewline(inp, "Data output interval in number of steps");
    opt     = "awh-nstout";
    nstOut_ = get_eint(inp, opt, 100000, wi);

    printStringNoNewline(inp, "Coordinate sampling interval in number of steps");
    opt             = "awh-nstsample";
    nstSampleCoord_ = get_eint(inp, opt, 10, wi);

    printStringNoNewline(inp, "Free energy and bias update interval in number of samples");
    opt                         = "awh-nsamples-update";
    numSamplesUpdateFreeEnergy_ = get_eint(inp, opt, 10, wi);

    printStringNoNewline(
            inp, "When true, biases with share-group>0 are shared between multiple simulations");
    opt                = "awh-share-multisim";
    shareBiasMultisim_ = (getEnum<Boolean>(inp, opt.c_str(), wi) != Boolean::No);

    printStringNoNewline(inp, "The number of independent AWH biases");
    opt         = "awh-nbias";
    int numBias = get_eint(inp, opt, 1, wi);
    /* Check this before starting to read the AWH biases. */
    if (numBias <= 0)
    {
        gmx_fatal(FARGS, "%s needs to be an integer > 0", opt.c_str());
    }

    /* Read the parameters specific to each AWH bias */
    for (int k = 0; k < numBias; k++)
    {
        bool        bComment  = (k == 0);
        std::string prefixawh = formatString("awh%d", k + 1);
        awhBiasParams_.emplace_back(inp, prefixawh, wi, bComment);
        checkInputConsistencyAwhBias(awhBiasParams_.back(), wi);
    }
    checkInputConsistencyAwh(*this, wi);
}

AwhParams::AwhParams(ISerializer* serializer, const bool tprWithoutGrowthFactor, const bool tprWithoutTargetMetricScaling)
{
    GMX_RELEASE_ASSERT(serializer->reading(),
                       "Can not use writing serializer to read AWH parameters");
    int numberOfBiases = awhBiasParams_.size();
    serializer->doInt(&numberOfBiases);
    serializer->doInt(&nstOut_);
    serializer->doInt64(&seed_);
    serializer->doInt(&nstSampleCoord_);
    serializer->doInt(&numSamplesUpdateFreeEnergy_);
    serializer->doEnumAsInt(&potentialEnum_);
    serializer->doBool(&shareBiasMultisim_);

    if (numberOfBiases > 0)
    {
        for (int k = 0; k < numberOfBiases; k++)
        {
            awhBiasParams_.emplace_back(serializer, tprWithoutGrowthFactor, tprWithoutTargetMetricScaling);
        }
    }
}

void AwhParams::serialize(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(),
                       "Can not use reading serializer to write AWH parameters");
    int numberOfBiases = numBias();
    serializer->doInt(&numberOfBiases);
    serializer->doInt(&nstOut_);
    serializer->doInt64(&seed_);
    serializer->doInt(&nstSampleCoord_);
    serializer->doInt(&numSamplesUpdateFreeEnergy_);
    serializer->doEnumAsInt(&potentialEnum_);
    serializer->doBool(&shareBiasMultisim_);

    if (numberOfBiases > 0)
    {
        for (int k = 0; k < numberOfBiases; k++)
        {
            awhBiasParams_[k].serialize(serializer);
        }
    }
}

/*! \brief
 * Gets the period of a pull coordinate.
 *
 * \param[in] pullCoordParams   The parameters for the pull coordinate.
 * \param[in] pbc               The PBC setup
 * \param[in] intervalLength    The length of the AWH interval for this pull coordinate
 * \returns the period (or 0 if not periodic).
 */
static double get_pull_coord_period(const t_pull_coord& pullCoordParams, const t_pbc& pbc, const real intervalLength)
{
    double period = 0;

    if (pullCoordParams.eGeom == PullGroupGeometry::Direction)
    {
        const real margin = 0.001;
        // Make dims periodic when the interval covers > 95%
        const real periodicFraction = 0.95;

        // Check if the pull direction is along a box vector
        for (int dim = 0; dim < pbc.ndim_ePBC; dim++)
        {
            const real boxLength    = norm(pbc.box[dim]);
            const real innerProduct = iprod(pullCoordParams.vec, pbc.box[dim]);
            if (innerProduct >= (1 - margin) * boxLength && innerProduct <= (1 + margin) * boxLength)
            {
                if (intervalLength > (1 + margin) * boxLength)
                {
                    gmx_fatal(FARGS,
                              "The AWH interval (%f nm) for a pull coordinate is larger than the "
                              "box size (%f nm)",
                              intervalLength,
                              boxLength);
                }

                if (intervalLength > periodicFraction * boxLength)
                {
                    period = boxLength;
                }
            }
        }
    }
    else if (pullCoordParams.eGeom == PullGroupGeometry::Dihedral)
    {
        /* The dihedral angle is periodic in -180 to 180 deg */
        period = 360;
    }

    return period;
}

/*! \brief
 * Checks if the given interval is defined in the correct periodic interval.
 *
 * \param[in] origin      Start value of interval.
 * \param[in] end         End value of interval.
 * \param[in] period      Period (or 0 if not periodic).
 * \returns true if the end point values are in the correct periodic interval.
 */
static bool intervalIsInPeriodicInterval(double origin, double end, double period)
{
    return (period == 0) || (std::fabs(origin) <= 0.5 * period && std::fabs(end) <= 0.5 * period);
}

/*! \brief
 * Checks if a value is within an interval.
 *
 * \param[in] origin      Start value of interval.
 * \param[in] end         End value of interval.
 * \param[in] period      Period (or 0 if not periodic).
 * \param[in] value       Value to check.
 * \returns true if the value is within the interval.
 */
static bool valueIsInInterval(double origin, double end, double period, double value)
{
    bool bIn_interval;

    if (period > 0)
    {
        if (origin < end)
        {
            /* The interval closes within the periodic interval */
            bIn_interval = (value >= origin) && (value <= end);
        }
        else
        {
            /* The interval wraps around the periodic boundary */
            bIn_interval = ((value >= origin) && (value <= 0.5 * period))
                           || ((value >= -0.5 * period) && (value <= end));
        }
    }
    else
    {
        bIn_interval = (value >= origin) && (value <= end);
    }

    return bIn_interval;
}

/*! \brief
 * Check if the starting configuration is consistent with the given interval.
 *
 * \param[in]     awhParams  AWH parameters.
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
static void checkInputConsistencyInterval(const AwhParams& awhParams, WarningHandler* wi)
{
    const auto& awhBiasParams = awhParams.awhBiasParams();
    for (int k = 0; k < gmx::ssize(awhBiasParams); k++)
    {
        const auto& dimParams = awhBiasParams[k].dimParams();
        for (int d = 0; d < gmx::ssize(dimParams); d++)
        {
            int    coordIndex = dimParams[d].coordinateIndex();
            double origin = dimParams[d].origin(), end = dimParams[d].end(),
                   period         = dimParams[d].period();
            double coordValueInit = dimParams[d].initialCoordinate();

            if ((period == 0) && (origin > end))
            {
                gmx_fatal(FARGS,
                          "For the non-periodic pull coordinates awh%d-dim%d-start (%f) cannot be "
                          "larger than awh%d-dim%d-end (%f)",
                          k + 1,
                          d + 1,
                          origin,
                          k + 1,
                          d + 1,
                          end);
            }

            /* Currently we assume symmetric periodic intervals, meaning we use [-period/2, period/2] as the reference interval.
               Make sure the AWH interval is within this reference interval.

               Note: we could fairly simply allow using a  more general interval (e.g. [x, x + period]) but it complicates
               things slightly and I don't see that there is a great need for it. It would also mean that the interval would
               depend on AWH input. Also, for dihedral angles you would always want the reference interval to be -180, +180,
               independent of AWH parameters.
             */
            if (!intervalIsInPeriodicInterval(origin, end, period))
            {
                gmx_fatal(FARGS,
                          "When using AWH with periodic pull coordinate geometries "
                          "awh%d-dim%d-start (%.8g) and "
                          "awh%d-dim%d-end (%.8g) should cover at most one period (%.8g) and take "
                          "values in between "
                          "minus half a period and plus half a period, i.e. in the interval [%.8g, "
                          "%.8g].",
                          k + 1,
                          d + 1,
                          origin,
                          k + 1,
                          d + 1,
                          end,
                          period,
                          -0.5 * period,
                          0.5 * period);
            }

            /* Warn if the pull initial coordinate value is not in the grid */
            if (!valueIsInInterval(origin, end, period, coordValueInit))
            {
                auto message = formatString(
                        "The initial coordinate value (%.8g) for pull coordinate index %d falls "
                        "outside "
                        "of the sampling nterval awh%d-dim%d-start (%.8g) to awh%d-dim%d-end "
                        "(%.8g). "
                        "This can lead to large initial forces pulling the coordinate towards the "
                        "sampling interval.",
                        coordValueInit,
                        coordIndex + 1,
                        k + 1,
                        d + 1,
                        origin,
                        k + 1,
                        d + 1,
                        end);
                wi->addWarning(message);
            }
        }
    }
}

/*! \brief
 * Sets AWH parameters, for one AWH pull dimension.
 *
 * \param[in,out] dimParams          AWH dimension parameters.
 * \param[in]     biasIndex             The index of the bias containing this AWH pull dimension.
 * \param[in]     dimIndex              The index of this AWH pull dimension.
 * \param[in]     pull_params           Pull parameters.
 * \param[in,out] pull_work             Pull working struct to register AWH bias in.
 * \param[in]     pbc                   A pbc information structure.
 * \param[in]     compressibility       Compressibility matrix for pressure coupling,
 * pass all 0 without pressure coupling.
 * \param[in,out] wi                    Struct for bookeeping warnings.
 *
 */
static void setStateDependentAwhPullDimParams(AwhDimParams*        dimParams,
                                              const int            biasIndex,
                                              const int            dimIndex,
                                              const pull_params_t& pull_params,
                                              pull_t*              pull_work,
                                              const t_pbc&         pbc,
                                              const tensor&        compressibility,
                                              WarningHandler*      wi)
{
    const t_pull_coord& pullCoordParams = pull_params.coord[dimParams->coordinateIndex()];

    if (pullCoordParams.eGeom == PullGroupGeometry::DirectionPBC)
    {
        gmx_fatal(FARGS,
                  "AWH does not support pull geometry '%s'. "
                  "If the maximum distance between the groups is always "
                  "less than half the box size, "
                  "you can use geometry '%s' instead.",
                  enumValueToString(PullGroupGeometry::DirectionPBC),
                  enumValueToString(PullGroupGeometry::Direction));
    }

    dimParams->setPeriod(
            get_pull_coord_period(pullCoordParams, pbc, dimParams->end() - dimParams->origin()));
    // We would like to check for scaling, but we don't have the full inputrec available here
    if (dimParams->period() > 0
        && !(pullCoordParams.eGeom == PullGroupGeometry::Angle
             || pullCoordParams.eGeom == PullGroupGeometry::Dihedral))
    {
        bool coordIsScaled = false;
        for (int d2 = 0; d2 < DIM; d2++)
        {
            if (pullCoordParams.vec[d2] != 0 && norm2(compressibility[d2]) != 0)
            {
                coordIsScaled = true;
            }
        }
        if (coordIsScaled)
        {
            std::string mesg = gmx::formatString(
                    "AWH dimension %d of bias %d is periodic with pull geometry '%s', "
                    "while you should be applying pressure scaling to the "
                    "corresponding box vector, this is not supported.",
                    dimIndex + 1,
                    biasIndex + 1,
                    enumValueToString(pullCoordParams.eGeom));
            wi->addWarning(mesg);
        }
    }

    /* The initial coordinate value, converted to external user units.
     * Note that we do not pass time here, as time is not allowed as
     * a variable for transformation coordinates when AWH is in use.
     */
    double initialCoordinate = get_pull_coord_value(pull_work, dimParams->coordinateIndex(), pbc);
    initialCoordinate *= pull_conversion_factor_internal2userinput(pullCoordParams);
    dimParams->setInitialCoordinate(initialCoordinate);
}

void setStateDependentAwhParams(AwhParams*           awhParams,
                                const pull_params_t& pull_params,
                                pull_t*              pull_work,
                                const matrix         box,
                                PbcType              pbcType,
                                const tensor&        compressibility,
                                const t_inputrec&    inputrec,
                                const real           initLambda,
                                const gmx_mtop_t&    mtop,
                                WarningHandler*      wi)
{
    /* The temperature is not really state depenendent but is not known
     * when read_awhParams is called (in get ir).
     * It is known first after do_index has been called in grompp.cpp.
     */
    GMX_RELEASE_ASSERT(haveConstantEnsembleTemperature(inputrec),
                       "AWH requires a constant ensemble temperaure");
    if (constantEnsembleTemperature(inputrec) <= 0)
    {
        gmx_fatal(FARGS, "AWH biasing is only supported for temperatures > 0");
    }

    t_pbc pbc;
    set_pbc(&pbc, pbcType, box);

    auto awhBiasParams = awhParams->awhBiasParams();
    for (int k = 0; k < awhParams->numBias(); k++)
    {
        auto awhBiasDimensionParams = awhBiasParams[k].dimParams();
        for (int d = 0; d < awhBiasParams[k].ndim(); d++)
        {
            AwhDimParams* dimParams = &awhBiasDimensionParams[d];
            if (dimParams->coordinateProvider() == AwhCoordinateProviderType::Pull)
            {
                const t_pull_coord& pullCoordParams = pull_params.coord[dimParams->coordinateIndex()];
                if (dimParams->coverDiameter() != 0.0
                    && (pullCoordParams.eGeom == PullGroupGeometry::Angle
                        || pullCoordParams.eGeom == PullGroupGeometry::Dihedral
                        || pullCoordParams.eGeom == PullGroupGeometry::AngleAxis))
                {
                    wi->addNote(

                            "Note that the unit of the AWH cover-diameter parameter for angle and "
                            "dihedral pull coordinates has recently changed from radian to "
                            "degrees");
                }

                setStateDependentAwhPullDimParams(
                        dimParams, k, d, pull_params, pull_work, pbc, compressibility, wi);
            }
            else
            {
                dimParams->setInitialCoordinate(initLambda);
                checkFepLambdaDimDecouplingConsistency(mtop, wi);
            }
        }
    }
    checkInputConsistencyInterval(*awhParams, wi);

    /* Register AWH as external potential with pull (for AWH dimensions that use pull as
     * reaction coordinate provider) to check consistency. */
    Awh::registerAwhWithPull(*awhParams, pull_work);
}

void checkAwhParams(const AwhParams& awhParams, const t_inputrec& ir, WarningHandler* wi)
{
    std::string opt;
    checkMtsConsistency(ir, wi);

    opt = "awh-nstout";
    if (awhParams.nstout() <= 0)
    {
        auto message = formatString("Not writing AWH output with AWH (%s = %d) does not make sense",
                                    opt.c_str(),
                                    awhParams.nstout());
        wi->addError(message);
    }
    /* This restriction can be removed by changing a flag of print_ebin() */
    if (ir.nstenergy == 0 || awhParams.nstout() % ir.nstenergy != 0)
    {
        auto message = formatString(
                "%s (%d) should be a multiple of nstenergy (%d)", opt.c_str(), awhParams.nstout(), ir.nstenergy);
        wi->addError(message);
    }

    opt = "awh-nsamples-update";
    if (awhParams.numSamplesUpdateFreeEnergy() <= 0)
    {
        wi->addError(opt + " needs to be an integer > 0");
    }

    bool       haveFepLambdaDim = false;
    const auto awhBiasParams    = awhParams.awhBiasParams();
    for (int k = 0; k < awhParams.numBias() && !haveFepLambdaDim; k++)
    {
        std::string prefixawh = formatString("awh%d", k + 1);
        checkBiasParams(awhBiasParams[k], prefixawh, ir, wi);
        /* Check if there is a FEP lambda dimension. */
        const auto dimParams = awhBiasParams[k].dimParams();
        haveFepLambdaDim = std::any_of(dimParams.begin(), dimParams.end(), [](const auto& dimParam) {
            return dimParam.coordinateProvider() == AwhCoordinateProviderType::FreeEnergyLambda;
        });
    }

    if (haveFepLambdaDim)
    {
        if (awhParams.nstSampleCoord() % ir.nstcalcenergy != 0)
        {
            opt          = "awh-nstsample";
            auto message = formatString(
                    "%s (%d) should be a multiple of nstcalcenergy (%d) when using AWH for "
                    "sampling an FEP lambda dimension",
                    opt.c_str(),
                    awhParams.nstSampleCoord(),
                    ir.nstcalcenergy);
            wi->addError(message);
        }
        if (awhParams.potential() != AwhPotentialType::Umbrella)
        {
            opt          = "awh-potential";
            auto message = formatString(
                    "%s (%s) must be set to %s when using AWH for sampling an FEP lambda dimension",
                    opt.c_str(),
                    enumValueToString(awhParams.potential()),
                    enumValueToString(AwhPotentialType::Umbrella));
            wi->addError(message);
        }
    }

    if (ir.init_step != 0)
    {
        wi->addError("With AWH init-step should be 0");
    }
}

bool awhHasFepLambdaDimension(const AwhParams& awhParams)
{
    for (const auto& biasParams : awhParams.awhBiasParams())
    {
        for (const auto& dimParams : biasParams.dimParams())
        {
            if (dimParams.coordinateProvider() == AwhCoordinateProviderType::FreeEnergyLambda)
            {
                return true;
            }
        }
    }

    return false;
}

} // namespace gmx
