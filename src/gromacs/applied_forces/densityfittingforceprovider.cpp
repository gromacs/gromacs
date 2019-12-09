/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingforceprovider.h"

#include <numeric>

#include "gromacs/compat/optional.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/densityfittingforce.h"
#include "gromacs/math/exponentialmovingaverage.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"

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
    scaleToLattice({ &sigmaInLatticeCoordinates, &sigmaInLatticeCoordinates + 1 });
    return { DVec{ sigmaInLatticeCoordinates[XX], sigmaInLatticeCoordinates[YY],
                   sigmaInLatticeCoordinates[ZZ] },
             nSigma };
}

} // namespace

/********************************************************************
 * DensityFittingForceProvider::Impl
 */

class DensityFittingForceProvider::Impl
{
public:
    //! \copydoc DensityFittingForceProvider(const DensityFittingParameters &parameters)
    Impl(const DensityFittingParameters&             parameters,
         basic_mdspan<const float, dynamicExtents3D> referenceDensity,
         const TranslateAndScale&                    transformationToDensityLattice,
         const LocalAtomSet&                         localAtomSet,
         int                                         pbcType,
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
    int                           pbcType_;

    //! Optionally scale the force according to a moving average of the similarity
    compat::optional<ExponentialMovingAverage> expAverageSimilarity_;
};

DensityFittingForceProvider::Impl::~Impl() = default;

DensityFittingForceProvider::Impl::Impl(const DensityFittingParameters&             parameters,
                                        basic_mdspan<const float, dynamicExtents3D> referenceDensity,
                                        const TranslateAndScale& transformationToDensityLattice,
                                        const LocalAtomSet&      localAtomSet,
                                        int                      pbcType,
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
    expAverageSimilarity_(compat::nullopt)
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
    referenceDensityCenter_ = { real(referenceDensity.extent(XX)) / 2,
                                real(referenceDensity.extent(YY)) / 2,
                                real(referenceDensity.extent(ZZ)) / 2 };
    transformationToDensityLattice_.scaleOperationOnly().inverseIgnoringZeroScale(
            { &referenceDensityCenter_, &referenceDensityCenter_ + 1 });
    // correct the reference density center for a shift
    // if the reference density does not have its origin at (0,0,0)
    RVec referenceDensityOriginShift(0, 0, 0);
    transformationToDensityLattice_({ &referenceDensityOriginShift, &referenceDensityOriginShift + 1 });
    transformationToDensityLattice_.scaleOperationOnly().inverseIgnoringZeroScale(
            { &referenceDensityOriginShift, &referenceDensityOriginShift + 1 });
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
    std::transform(std::cbegin(localAtomSet_.localIndex()), std::cend(localAtomSet_.localIndex()),
                   std::begin(transformedCoordinates_),
                   [&forceProviderInput](int index) { return forceProviderInput.x_[index]; });

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

    std::vector<real> amplitudes =
            amplitudeLookup_(forceProviderInput.mdatoms_, localAtomSet_.localIndex());

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
                 gaussTransform_.view().data(), &forceProviderInput.cr_);
    }

    // calculate grid derivative
    const DensitySimilarityMeasure::density& densityDerivative =
            measure_.gradient(gaussTransform_.constView());
    // calculate forces
    forces_.resize(localAtomSet_.numAtomsLocal());
    std::transform(
            std::begin(transformedCoordinates_), std::end(transformedCoordinates_), std::begin(amplitudes),
            std::begin(forces_), [&densityDerivative, this](const RVec r, real amplitude) {
                return densityFittingForce_.evaluateForce({ r, amplitude }, densityDerivative);
            });

    transformationToDensityLattice_.scaleOperationOnly().inverseIgnoringZeroScale(forces_);

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
    if (MASTER(&(forceProviderInput.cr_)))
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
                                                         int                 pbcType,
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

const DensityFittingForceProviderState& DensityFittingForceProvider::stateToCheckpoint()
{
    return impl_->stateToCheckpoint();
}

} // namespace gmx
