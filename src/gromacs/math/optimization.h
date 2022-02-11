/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *
 * \brief Declare function optimization routines.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */

#ifndef GMX_MATH_OPTIMIZATION_H
#define GMX_MATH_OPTIMIZATION_H

#include <functional>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \internal
 *  \brief Compiles results of an a function optimisation.
 */
struct OptimisationResult
{
    //! The coordinates at which the optimal function value has been found
    std::vector<real> coordinates_;
    //! The value of the function at the optimal coordinates
    real functionValue_;
};

/*! \brief Derivative-free downhill simplex optimisation.
 *
 * Find a local minimum of an N-dimensional mapping
 * \f$\mathbb{R}^N\to\mathbb{R}\f$ using the downhill simplex algorithm by
 * Nelder and Mead as described in
 *
 *       Saša Singer and John Nelder (2009), Scholarpedia, 4(7):2928.
 *       doi:10.4249/scholarpedia.2928
 *
 * Stops when the oriented simplex length is less than a constant factor times the
 * initial lengths or when a maximum step size is reached.
 *
 * For best performance pre-condition problem magnitudes to 1.
 *
 * The following algorithm is implemented in this function
 *  1 Define the N+1 vertices of the initial simplex
 *      The inital simplex is constructed from the initial guess and N
 *      additional vertices by adding 0.05 to the initial guess (or 0.0025 if
 *      the initial guess is the null vector) from the initial vertex (in line
 *      with usual implementations).
 *
 *  1a Sort vertices according to function value with the lowest function value
 *    first in order to minimize the function.
 *
 *  2 Calculate the centroid of the simplex as arithmetic mean of all vertices
 *    except the worst, \f$ x_c = \frac1n \sum_{i=1}{N} x_i\f$.
 *
 *  3 Reflect the worst simplex vertex (the one with the highest function value)
 *    at the centroid to obtain a reflection point
 *    \f$ x_r = x_c + \alpha (x_c - x_{N+1}) \f$ which lies outside the vertex.
 *
 *  3a Replace worst point with reflection point if reflection point function
 *    value is better than second worst point, but not better than best and go
 *    to 1a.
 *
 *  4 If the reflection point is better than all other points so far, attempt
 *    an expansion by calculating the expansion point at
 *    \f$ x_e = x_c + \gamma (x_r - x_c) \f$. Swap out the worst point in the
 *    vertex with the expansion point if better than reflection point, otherwise
 *    use the reflection point and go to 1a.
 *
 *  5 Attempt contraction, because reflection was not successful;
 *    \f$ x_t = x_c + \rho (x_{N+1} - x_c) \f$. If the contraction point is
 *    better than the worst point, swap out worst point with contracted point
 *    and go to 1a.
 *
 *  6 Shrink the vertex. Replace all points except the best one with
 *    \f$ x_i = x_1 + \sigma (x_i - x_1) \f$ and go to 1a.
 *
 * \param[in] functionToMinimize function to be minimized
 * \param[in] initialGuess of coordinates
 * \param[in] minimumRelativeSimplexLength minimal oriented simplex length with
 *                                         respect to initial simplex
 * \param[in] maxSteps to run algorithm for
 *
 * \returns the lowest found function value and corresponding coordinates.
 */
OptimisationResult nelderMead(const std::function<real(ArrayRef<const real>)>& functionToMinimize,
                              ArrayRef<const real>                             initialGuess,
                              real minimumRelativeSimplexLength = 1e-8,
                              int  maxSteps                     = 10'000);

} // namespace gmx
#endif // GMX_MATH_OPTIMIZATION_H
