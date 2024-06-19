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
 * Implements force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingforceprovider.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/densityfittingforce.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/strconvert.h"

#include "densityfittingamplitudelookup.h"
#include "densityfittingparameters.h"

namespace gmx
{

namespace
{

/*! \internal \brief Generate the spread kernel from Gaussian parameters.
 *
 * \param[in] sigma the width of the Gaussian to be spread
 * \param[in] nSigma the range of the Gaussian in multiples of sigma
 * \param[in] scaleToLattice the coordinate transformation into the spreading lattice
 * \returns A Gauss-transform kernel shape
 */
GaussianSpreadKernelParameters::Shape makeSpreadKernel(real sigma, real nSigma, const ScaleCoordinates& scaleToLattice)
{
    RVec sigmaInLatticeCoordinates{ sigma, sigma, sigma };
    scaleToLattice(&sigmaInLatticeCoordinates);
    return { DVec{ sigmaInLatticeCoordinates[XX], sigmaInLatticeCoordinates[YY], sigmaInLatticeCoordinates[ZZ] },
             nSigma };
}

} // namespace

/********************************************************************
 * DensityFittingForceProviderState
 */

const std::string DensityFittingForceProviderState::adaptiveForceConstantScaleName_ =
        "adaptiveForceConstantScale";

const std::string DensityFittingForceProviderState::exponentialMovingAverageStateName_ =
        "exponentialMovingAverageState";

const std::string DensityFittingForceProviderState::stepsSinceLastCalculationName_ =
        "stepsSinceLastCalculation";

void DensityFittingForceProviderState::writeState(KeyValueTreeObjectBuilder kvtBuilder,
                                                  const std::string&        identifier) const
{
    writeKvtCheckpointValue(
            stepsSinceLastCalculation_, stepsSinceLastCalculationName_, identifier, kvtBuilder);
    writeKvtCheckpointValue(
            adaptiveForceConstantScale_, adaptiveForceConstantScaleName_, identifier, kvtBuilder);

    KeyValueTreeObjectBuilder exponentialMovingAverageKvtEntry =
            kvtBuilder.addObject(identifier + "-" + exponentialMovingAverageStateName_);
    exponentialMovingAverageStateAsKeyValueTree(exponentialMovingAverageKvtEntry,
                                                exponentialMovingAverageState_);
}

void DensityFittingForceProviderState::readState(const KeyValueTreeObject& kvtData,
                                                 const std::string&        identifier)
{
    readKvtCheckpointValue(compat::make_not_null(&stepsSinceLastCalculation_),
                           stepsSinceLastCalculationName_,
                           identifier,
                           kvtData);
    readKvtCheckpointValue(compat::make_not_null(&adaptiveForceConstantScale_),
                           adaptiveForceConstantScaleName_,
                           identifier,
                           kvtData);

    if (kvtData.keyExists(identifier + "-" + exponentialMovingAverageStateName_))
    {
        exponentialMovingAverageState_ = exponentialMovingAverageStateFromKeyValueTree(
                kvtData[identifier + "-" + exponentialMovingAverageStateName_].asObject());
    }
}

void DensityFittingForceProviderState::broadcastState(MPI_Comm communicator, bool isParallelRun)
{
    if (isParallelRun)
    {
        block_bc(communicator, stepsSinceLastCalculation_);
        block_bc(communicator, adaptiveForceConstantScale_);
        block_bc(communicator, exponentialMovingAverageState_);
    }
}


/********************************************************************
 * DensityFittingForceProvider::Impl
 */

class DensityFittingForceProvider::Impl
{
public:
    //! \copydoc DensityFittingForceProvider::DensityFittingForceProvider
    Impl(const DensityFittingParameters&             parameters,
         basic_mdspan<const float, dynamicExtents3D> referenceDensity,
         const TranslateAndScale&                    transformationToDensityLattice,
         const LocalAtomSet&                         localAtomSet,
         PbcType                                     pbcType,
         double                                      simulationTimeStep,
         const DensityFittingForceProviderState&     state);
    ~Impl();
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput);
    const DensityFittingForceProviderState& stateToCheckpoint();

private:
    DensityFittingForceProviderState state();
    const DensityFittingParameters&  parameters_;
    DensityFittingForceProviderState state_;
    DensityFittingForceProviderState stateToCheckpoint_;
    LocalAtomSet                     localAtomSet_;

    GaussianSpreadKernelParameters::Shape spreadKernel_;
    GaussTransform3D                      gaussTransform_;
    DensitySimilarityMeasure              measure_;
    DensityFittingForce                   densityFittingForce_;
    //! the local atom coordinates transformed into the grid coordinate system
    std::vector<RVec>             transformedCoordinates_;
    std::vector<RVec>             forces_;
    DensityFittingAmplitudeLookup amplitudeLookup_;
    TranslateAndScale             transformationToDensityLattice_;
    RVec                          referenceDensityCenter_;
    PbcType                       pbcType_;

    //! Optionally scale the force according to a moving average of the similarity
    std::optional<ExponentialMovingAverage> expAverageSimilarity_;

    //! Optionally translate the structure
    std::optional<AffineTransformation> affineTransformation_;
};

DensityFittingForceProvider::Impl::~Impl() = default;

DensityFittingForceProvider::Impl::Impl(const DensityFittingParameters&             parameters,
                                        basic_mdspan<const float, dynamicExtents3D> referenceDensity,
                                        const TranslateAndScale& transformationToDensityLattice,
                                        const LocalAtomSet&      localAtomSet,
                                        PbcType                  pbcType,
                                        double                   simulationTimeStep,
                                        const DensityFittingForceProviderState& state) :
    parameters_(parameters),
    state_(state),
    localAtomSet_(localAtomSet),
    spreadKernel_(makeSpreadKernel(parameters_.gaussianTransformSpreadingWidth_,
                                   parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_,
                                   transformationToDensityLattice.scaleOperationOnly())),
    gaussTransform_(referenceDensity.extents(), spreadKernel_),
    measure_(parameters.similarityMeasureMethod_, referenceDensity),
    densityFittingForce_(spreadKernel_),
    transformedCoordinates_(localAtomSet_.numAtomsLocal()),
    amplitudeLookup_(parameters_.amplitudeLookupMethod_),
    transformationToDensityLattice_(transformationToDensityLattice),
    pbcType_(pbcType),
    expAverageSimilarity_(std::nullopt)
{
    if (parameters_.adaptiveForceScaling_)
    {
        GMX_ASSERT(simulationTimeStep > 0,
                   "Simulation time step must be larger than zero for adaptive for scaling.");
        expAverageSimilarity_.emplace(ExponentialMovingAverage(
                parameters_.adaptiveForceScalingTimeConstant_
                        / (simulationTimeStep * parameters_.calculationIntervalInSteps_),
                state.exponentialMovingAverageState_));
    }

    // set up optional coordinate translation if the translation string contains a vector
    const std::optional<std::array<real, 3>> translationParametersAsArray =
            parsedArrayFromInputString<real, 3>(parameters_.translationString_);
    // set up optional coordinate transformation if the transformation string contains data
    const std::optional<std::array<real, 9>> transformationMatrixParametersAsArray =
            parsedArrayFromInputString<real, 9>(parameters_.transformationMatrixString_);
    if (translationParametersAsArray || transformationMatrixParametersAsArray)
    {
        Matrix3x3 translationMatrix = transformationMatrixParametersAsArray.has_value()
                                              ? *transformationMatrixParametersAsArray
                                              : identityMatrix<real, 3>();
        RVec      translationVector = translationParametersAsArray.has_value()
                                              ? RVec((*translationParametersAsArray)[XX],
                                                (*translationParametersAsArray)[YY],
                                                (*translationParametersAsArray)[ZZ])
                                              : RVec(0, 0, 0);
        affineTransformation_.emplace(translationMatrix.asConstView(), translationVector);
    }

    referenceDensityCenter_ = { real(referenceDensity.extent(XX)) / 2,
                                real(referenceDensity.extent(YY)) / 2,
                                real(referenceDensity.extent(ZZ)) / 2 };
    transformationToDensityLattice_.scaleOperationOnly().inverseIgnoringZeroScale(&referenceDensityCenter_);
    // correct the reference density center for a shift
    // if the reference density does not have its origin at (0,0,0)
    RVec referenceDensityOriginShift(0, 0, 0);
    transformationToDensityLattice_(&referenceDensityOriginShift);
    transformationToDensityLattice_.scaleOperationOnly().inverseIgnoringZeroScale(&referenceDensityOriginShift);
    referenceDensityCenter_ -= referenceDensityOriginShift;
}

void DensityFittingForceProvider::Impl::calculateForces(const ForceProviderInput& forceProviderInput,
                                                        ForceProviderOutput* forceProviderOutput)
{
    // TODO change if checkpointing moves to the start of the md loop
    stateToCheckpoint_ = state();
    // do nothing but count number of steps when not in density fitting step
    if (state_.stepsSinceLastCalculation_ % parameters_.calculationIntervalInSteps_ != 0)
    {
        ++(state_.stepsSinceLastCalculation_);
        return;
    }

    state_.stepsSinceLastCalculation_ = 1;

    transformedCoordinates_.resize(localAtomSet_.numAtomsLocal());
    // pick and copy atom coordinates
    std::transform(std::cbegin(localAtomSet_.localIndex()),
                   std::cend(localAtomSet_.localIndex()),
                   std::begin(transformedCoordinates_),
                   [&forceProviderInput](int index) { return forceProviderInput.x_[index]; });

    // apply additional structure transformations
    if (affineTransformation_)
    {
        (*affineTransformation_)(transformedCoordinates_);
    }

    // pick periodic image that is closest to the center of the reference density
    {
        t_pbc pbc;
        set_pbc(&pbc, pbcType_, forceProviderInput.box_);
        for (RVec& x : transformedCoordinates_)
        {
            rvec dx;
            pbc_dx(&pbc, x, referenceDensityCenter_, dx);
            x = referenceDensityCenter_ + dx;
        }
    }

    // transform local atom coordinates to density grid coordinates
    transformationToDensityLattice_(transformedCoordinates_);

    // spread atoms on grid
    gaussTransform_.setZero();

    std::vector<real> amplitudes = amplitudeLookup_(
            forceProviderInput.chargeA_, forceProviderInput.massT_, localAtomSet_.localIndex());

    if (parameters_.normalizeDensities_)
    {
        real sum = std::accumulate(std::begin(amplitudes), std::end(amplitudes), 0.);
        if (havePPDomainDecomposition(&forceProviderInput.cr_))
        {
            gmx_sum(1, &sum, &forceProviderInput.cr_);
        }
        for (real& amplitude : amplitudes)
        {
            amplitude /= sum;
        }
    }

    auto amplitudeIterator = amplitudes.cbegin();

    for (const auto& r : transformedCoordinates_)
    {
        gaussTransform_.add({ r, *amplitudeIterator });
        ++amplitudeIterator;
    }

    // communicate grid
    if (havePPDomainDecomposition(&forceProviderInput.cr_))
    {
        // \todo update to real once GaussTransform class returns real
        gmx_sumf(gaussTransform_.view().mapping().required_span_size(),
                 gaussTransform_.view().data(),
                 &forceProviderInput.cr_);
    }

    // calculate grid derivative
    const DensitySimilarityMeasure::density& densityDerivative =
            measure_.gradient(gaussTransform_.constView());
    // calculate forces
    forces_.resize(localAtomSet_.numAtomsLocal());
    std::transform(
            std::begin(transformedCoordinates_),
            std::end(transformedCoordinates_),
            std::begin(amplitudes),
            std::begin(forces_),
            [&densityDerivative, this](const RVec r, real amplitude) {
                return densityFittingForce_.evaluateForce({ r, amplitude }, densityDerivative);
            });

    // correct forces for coordinate transformations with chain rule
    // F = -k d U(transform(x)) / d x =
    //    k * -d U(transform(x)) / d transform(x) * d(transform(x)) / d x
    //        --------- calculated above --------   ---correction below---

    // correction for coordinate transformation into density lattice
    transformationToDensityLattice_.scaleOperationOnly()(forces_);
    // correction for affine coordinate transformation
    if (affineTransformation_)
    {
        const Matrix3x3 gradient = affineTransformation_->gradient();
        for (RVec currentForce : forces_)
        {
            matrixVectorMultiply(gradient, &currentForce);
        }
    }

    // multiply with the current force constant
    auto       densityForceIterator = forces_.cbegin();
    const real effectiveForceConstant = state_.adaptiveForceConstantScale_ * parameters_.calculationIntervalInSteps_
                                        * parameters_.forceConstant_;
    for (const auto localAtomIndex : localAtomSet_.localIndex())
    {
        forceProviderOutput->forceWithVirial_.force_[localAtomIndex] +=
                effectiveForceConstant * *densityForceIterator;
        ++densityForceIterator;
    }

    const float similarity = measure_.similarity(gaussTransform_.constView());
    if (MAIN(&(forceProviderInput.cr_)))
    {
        // calculate corresponding potential energy
        const real energy = -similarity * parameters_.forceConstant_ * state_.adaptiveForceConstantScale_;
        forceProviderOutput->enerd_.term[F_DENSITYFITTING] += energy;
    }

    if (expAverageSimilarity_.has_value())
    {
        expAverageSimilarity_->updateWithDataPoint(similarity);

        if (expAverageSimilarity_->increasing())
        {
            state_.adaptiveForceConstantScale_ /= 1._real + expAverageSimilarity_->inverseTimeConstant();
        }
        else
        {
            state_.adaptiveForceConstantScale_ *=
                    1._real + 2 * expAverageSimilarity_->inverseTimeConstant();
        }
    }
}

DensityFittingForceProviderState DensityFittingForceProvider::Impl::state()
{
    if (expAverageSimilarity_.has_value())
    {
        state_.exponentialMovingAverageState_ = expAverageSimilarity_->state();
    }
    return state_;
}

const DensityFittingForceProviderState& DensityFittingForceProvider::Impl::stateToCheckpoint()
{
    return stateToCheckpoint_;
}
/********************************************************************
 * DensityFittingForceProvider
 */

DensityFittingForceProvider::~DensityFittingForceProvider() = default;

DensityFittingForceProvider::DensityFittingForceProvider(const DensityFittingParameters& parameters,
                                                         basic_mdspan<const float, dynamicExtents3D> referenceDensity,
                                                         const TranslateAndScale& transformationToDensityLattice,
                                                         const LocalAtomSet& localAtomSet,
                                                         PbcType             pbcType,
                                                         double              simulationTimeStep,
                                                         const DensityFittingForceProviderState& state) :
    impl_(new Impl(parameters, referenceDensity, transformationToDensityLattice, localAtomSet, pbcType, simulationTimeStep, state))
{
}

void DensityFittingForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                                  ForceProviderOutput*      forceProviderOutput)
{
    impl_->calculateForces(forceProviderInput, forceProviderOutput);
}

void DensityFittingForceProvider::writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting,
                                                      const std::string&           moduleName)
{
    impl_->stateToCheckpoint().writeState(checkpointWriting.builder_, moduleName);
}

} // namespace gmx
