/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/topology/idef.h"

#define PROBABILITYCUTOFF (-100)
/* we don't bother evaluating if events are more rare than exp(-100) = 3.7x10^-44 (from repl_ex.cpp) */

namespace gmx
{

// all constructor arguments receive their values from inputrec, which is available on all PP nodes => can be instanciated with same data on all PP nodes
MetropolisStepMehlig::MetropolisStepMehlig(const real temperatureEnsemble, const real temperatureVelocities, const int seed, const int initStep, const bool startingFromCheckpoint) : oldPotentialEnergy_(0),
                                                                                                                                                                                      initialKineticEnergy_(-1),
                                                                                                                                                                                      temperatureEnsemble_(temperatureEnsemble),
                                                                                                                                                                                      temperatureVelocities_(temperatureVelocities),
                                                                                                                                                                                      initStep_(initStep),
                                                                                                                                                                                      startingFromCheckpoint_(startingFromCheckpoint),
                                                                                                                                                                                      rng_(seed, RandomDomain::HybridMCMD) // Is this ok?
{
}

// internal function: only used by MetropolisStepMehlig::accept
double MetropolisStepMehlig::calculateMetropolisCriterion(const double newPotentialEnergy, const double newKineticEnergy) const
{
    double deltaPotentialEnergy = oldPotentialEnergy_ - newPotentialEnergy;
    double deltaKineticEnergy   = initialKineticEnergy_ - newKineticEnergy;
    real   betaEnsemble         = BOLTZ * temperatureEnsemble_;
    real   betaVelocities       = BOLTZ * temperatureVelocities_;
    double exponent             = betaEnsemble * deltaPotentialEnergy + betaVelocities * deltaKineticEnergy;
    double result;
    if (exponent < PROBABILITYCUTOFF)
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

// internal function: only used by MetropolisStepMehlig::accept
double MetropolisStepMehlig::drawRandomNumber(const int64_t step)
{
    // compare mdlib/expanded.cpp and mdlib/repl_ex.cpp
    rng_.restart(step, 0);
    uniformRealDist_.reset();
    return uniformRealDist_(rng_);
}

// internal function: only used by MetropolisStepMehlig::accept
void MetropolisStepMehlig::updatePotentialEnergy(const double newPotentialEnergy)
{
    oldPotentialEnergy_ = newPotentialEnergy;
}

//! \brief Setter to be used by HybridMCMDVelocities (The Metropolis criterion needs the kinetic energy immediately after AcceptOrRewind has been called.)
void MetropolisStepMehlig::setInitialKineticEnergy(const double initialKineticEnergy)
{
    initialKineticEnergy_ = initialKineticEnergy;
}

//! \brief Returns a boolean on whether the Metropolis step has accepted or rejected the proposed configuration; is used by AcceptOrRewind.
bool MetropolisStepMehlig::accept(const int64_t step, const gmx_enerdata_t *enerd)
{
    bool accepted;
    if ((step != 0) && (step != initStep_) && (!startingFromCheckpoint_))
    {
        // Let's check whether the first back-up of potential and kinetic energy has worked
        GMX_ASSERT(oldPotentialEnergy_ != 0, "Unless you are simulating a system of non-interacting particles, the first back-up of the potential energy has not worked.");
        GMX_ASSERT(initialKineticEnergy_ >= 0, "The initial kinetic energy has not been set correctly. Please check whether the interaction with HybridMCMDVelocities has been set up correctly.");

        double metropolisCriterion = calculateMetropolisCriterion(enerd->term[F_EPOT], enerd->term[F_EKIN]);
        double randomNumber        = drawRandomNumber(step);
        accepted = (randomNumber < metropolisCriterion);
        if (accepted)
        {
            updatePotentialEnergy(enerd->term[F_EPOT]);
        }
    }
    // first call (only back-up energies at t = 0 ps)
    else
    {
        // startingFromCheckpoint_ has to be set to false during first call or structure will always be accepted
        if (startingFromCheckpoint_)
        {
            startingFromCheckpoint_ = false;
        }
        updatePotentialEnergy(enerd->term[F_EPOT]);
        accepted = true;
    }
    return accepted;
}

} // namespace
