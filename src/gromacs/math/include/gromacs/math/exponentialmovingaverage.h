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
/*! \libinternal \file
 * \brief
 * Declares an exponential moving average class.
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_EXPONENTIALMOVINGAVERAGE_H
#define GMX_MATH_EXPONENTIALMOVINGAVERAGE_H

#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \libinternal \brief Store the state of exponential moving averages.
 */
struct ExponentialMovingAverageState
{
    //! The weighted sum
    real weightedSum_ = 0;

    //! The weighted count, used for bias correction
    real weightedCount_ = 0;

    //! Remember if adding the latest data point increased the average
    bool increasing_ = false;
};

//! Convert the exponential moving average state as key-value-tree object
void exponentialMovingAverageStateAsKeyValueTree(KeyValueTreeObjectBuilder            builder,
                                                 const ExponentialMovingAverageState& state);

//! Sets the expoential moving average state from a key-value-tree object
ExponentialMovingAverageState exponentialMovingAverageStateFromKeyValueTree(const KeyValueTreeObject& object);

/*! \libinternal
 * \brief Evaluate the exponential moving average with bias correction.
 *
 * The exponential moving average at the 0th data point \f$Y_0\f$ is
 * \f$ S_0 = Y_0 \f$ and at the n-th data point \f$Y_n\f$ with n>0 it is
 * \f$ S_n = \alpha Y_n + (1-\alpha) S_{n-1} \f$, where the smoothing factor
 * \f$\alpha=1/t\f$ is determined via a time constant \f$t\f$.
 *
 * To avoid large impact of the first data point in a "burn-in" phase, the
 * weight of points are unbiased by substituting for \f$S_{n-1}\f$ above,
 * \f$\hat{S}_{n-1} = S_{n-1} / (1-\alpha^{n})\f$.
 */
class ExponentialMovingAverage
{
public:
    /*! \brief Construct by setting the time constant and state.
     * Allows reinitiating with data from memory.
     * \param[in] timeConstant time in number of data points
     * \param[in] state of the exponential moving average
     * \throws InconsistentInputError if timeConstant < 1
     */
    ExponentialMovingAverage(real timeConstant, const ExponentialMovingAverageState& state = {});

    //! Update the moving average with a data point
    void updateWithDataPoint(real dataPoint);

    //! The exponential weighted average with bias correction
    real biasCorrectedAverage() const;

    //! Returns true if last added data point increased the average
    bool increasing() const;

    //! Return the current state of the exponential moving average
    const ExponentialMovingAverageState& state() const;

    //! The inverse time constant for the exponential moving average
    real inverseTimeConstant() const;

private:
    //! The current state of the exponential moving average
    ExponentialMovingAverageState state_;

    //! The inverse time constant for the exponential moving average
    real inverseTimeConstant_;
};

} // namespace gmx

#endif
