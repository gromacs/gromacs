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
 * Declares force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H

#include <memory>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/exponentialmovingaverage.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

struct DensityFittingParameters;

/*! \internal
 * \brief Parameters defining the internal density fitting force provider state.
 */
struct DensityFittingForceProviderState
{
    /*! \brief The steps since the last force calculation.
     *  Used if density fitting is to be calculated every N steps.
     */
    std::int64_t stepsSinceLastCalculation_ = 0;
    //! The state of the exponential moving average of the similarity measure
    ExponentialMovingAverageState exponentialMovingAverageState_ = {};
    //! An additional factor scaling the force for adaptive force scaling
    real adaptiveForceConstantScale_ = 1.0_real;
};

/*! \internal \brief
 * Implements IForceProvider for density-fitting forces.
 */
class DensityFittingForceProvider final : public IForceProvider
{
public:
    //! Construct force provider for density fitting from its parameters
    DensityFittingForceProvider(const DensityFittingParameters&             parameters,
                                basic_mdspan<const float, dynamicExtents3D> referenceDensity,
                                const TranslateAndScale& transformationToDensityLattice,
                                const LocalAtomSet&      localAtomSet,
                                int                      pbcType,
                                double                   simulationTimeStep,
                                const DensityFittingForceProviderState& state);
    ~DensityFittingForceProvider();
    /*!\brief Calculate forces that maximise goodness-of-fit with a reference density map.
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    /*! \brief Return the state of the forceprovider to be checkpointed
     * TODO update this routine if checkpointing is moved to the beginning of
     *      the md loop
     */
    const DensityFittingForceProviderState& stateToCheckpoint();

private:
    class Impl;
    PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H
