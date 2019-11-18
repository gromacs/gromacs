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
 * Implements routines to calculate an exponential moving average.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "exponentialmovingaverage.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"

namespace gmx
{

//! Convert the exponential moving average state as key-value-tree object
void exponentialMovingAverageStateAsKeyValueTree(KeyValueTreeObjectBuilder            builder,
                                                 const ExponentialMovingAverageState& state)
{
    builder.addValue<real>("weighted-sum", state.weightedSum_);
    builder.addValue<real>("weighted-count", state.weightedCount_);
    builder.addValue<bool>("increasing", state.increasing_);
}

//! Sets the exponential moving average state from a key-value-tree object
ExponentialMovingAverageState exponentialMovingAverageStateFromKeyValueTree(const KeyValueTreeObject& object)
{
    const real weightedSum   = object["weighted-sum"].cast<real>();
    const real weightedCount = object["weighted-count"].cast<real>();
    const bool increasing    = object["increasing"].cast<bool>();
    return { weightedSum, weightedCount, increasing };
}

ExponentialMovingAverage::ExponentialMovingAverage(real timeConstant,
                                                   const ExponentialMovingAverageState& state) :
    state_(state)
{
    if (timeConstant < 1)
    {
        GMX_THROW(InconsistentInputError(
                "Lag time may not be negative or zero for exponential moving averages."));
    }
    inverseTimeConstant_ = 1. / timeConstant;
}

void ExponentialMovingAverage::updateWithDataPoint(real dataPoint)
{
    state_.weightedSum_   = dataPoint + (1 - inverseTimeConstant_) * state_.weightedSum_;
    state_.weightedCount_ = 1 + (1 - inverseTimeConstant_) * state_.weightedCount_;

    state_.increasing_ = dataPoint * state_.weightedCount_ > state_.weightedSum_;
}

const ExponentialMovingAverageState& ExponentialMovingAverage::state() const
{
    return state_;
}

real ExponentialMovingAverage::biasCorrectedAverage() const
{
    return state_.weightedSum_ / state_.weightedCount_;
}

bool ExponentialMovingAverage::increasing() const
{
    return state_.increasing_;
}

real ExponentialMovingAverage::inverseTimeConstant() const
{
    return inverseTimeConstant_;
}

} // namespace gmx
