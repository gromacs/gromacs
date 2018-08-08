/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Implements the class MetropolisStepMehlig.
 *
 * \author Sebastian Wingbermuehle
 * \ingroup module_hybridMCMD
 */

#include "gmxpre.h"

#include "metropolis.h"

#include <cmath>

#include "gromacs/math/units.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/topology/idef.h"

namespace gmx
{

/* all constructor arguments receive their values from inputrec, which is available on all PP nodes => can be instanciated with same data on all PP nodes */
MetropolisStepMehlig::MetropolisStepMehlig(const real    temperatureEnsemble,
                                           const real    temperatureVelocities,
                                           const int64_t seed)                  : oldPotentialEnergy_(0),
                                                                                  initialKineticEnergy_(-1),
                                                                                  temperatureEnsemble_(temperatureEnsemble),
                                                                                  temperatureVelocities_(temperatureVelocities),
                                                                                  seed_(seed),
                                                                                  isFirstCall_(true)
{
}

double MetropolisStepMehlig::calculateMetropolisCriterion(const double newPotentialEnergy, const double newKineticEnergy) const
{
    double deltaPotentialEnergy = oldPotentialEnergy_ - newPotentialEnergy;
    double deltaKineticEnergy   = initialKineticEnergy_ - newKineticEnergy;
    real   betaEnsemble         = 1.0/(BOLTZ * temperatureEnsemble_);
    real   betaVelocities       = 1.0/(BOLTZ * temperatureVelocities_);
    double exponent             = betaEnsemble * deltaPotentialEnergy + betaVelocities * deltaKineticEnergy;

    /* we don't bother evaluating if events are more rare than exp(-100) = 3.7x10^-44 */
    constexpr double minExponent = -100;

    double           result;

    if (exponent < minExponent)
    {
        result = 0.0;
    }
    else if (exponent > 0.0)
    {
        result = 1.0;
    }
    else
    {
        result = exp(exponent);
    }

    return result;
}

double MetropolisStepMehlig::drawRandomNumber(const int64_t step)
{
    ThreeFry2x64<0>               rng(seed_, RandomDomain::HybridMCMD);
    rng.restart(step, 0);
    UniformRealDistribution<real> uniformRealDist;
    uniformRealDist.reset();
    return uniformRealDist(rng);
}

void MetropolisStepMehlig::updatePotentialEnergy(const double newPotentialEnergy)
{
    oldPotentialEnergy_ = newPotentialEnergy;
}

void MetropolisStepMehlig::setInitialKineticEnergy(const double initialKineticEnergy)
{
    initialKineticEnergy_ = initialKineticEnergy;
}

bool MetropolisStepMehlig::accept(const int64_t step, const gmx_enerdata_t *enerd)
{
    bool accepted;
    /* We must always accept the starting structure of the simulation and every structure
     * that is read in from a checkpoint file.
     * In both cases, MetropolisStepMehlig has just been constructed and accept is called for the first time.
     */
    if (!isFirstCall_)
    {
        double metropolisCriterion = calculateMetropolisCriterion(enerd->term[F_EPOT], enerd->term[F_EKIN]);
        double randomNumber        = drawRandomNumber(step);
        accepted = (randomNumber < metropolisCriterion);
        if (accepted)
        {
            updatePotentialEnergy(enerd->term[F_EPOT]);
        }
    }
    /* first call (only back-up potential energy) */
    else
    {
        /* isFirstCall_ has to be set to false during first call or structure will always be accepted */
        isFirstCall_ = false;
        accepted     = true;
        updatePotentialEnergy(enerd->term[F_EPOT]);
    }
    return accepted;
}

} // namespace
