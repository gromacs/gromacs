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
/*! \libinternal \file
 *
 * \brief Declare classes to aid Nelder-Mead downhill simplex optimisation.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */

#ifndef GMX_MATH_NELDERMEAD_H
#define GMX_MATH_NELDERMEAD_H

#include <functional>
#include <list>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \internal
 * \brief Tie together coordinate and function value at this coordinate.
 */
struct RealFunctionvalueAtCoordinate
{
    //! Vertex coordinate
    std::vector<real> coordinate_;
    //! Function value at this coordinate
    real value_;
};

/*! \internal
 * \brief The simplex for the Nelder-Mead algorithm.
 *
 * Contains the N+1 simplex N-dimensional coordinates and its function values.
 * Allows for simplex manipulations as needed for the Nelder-Mead algorithm.
 *
 * \note Keeps the simplex sorted according to function values with the simplex
 *       at the lowest function value first.
 */
class NelderMeadSimplex
{
public:
    /*! \brief Set up Nelder-Mead simplex from an initial guess.
     *
     * \note Triggers N+1 function evaluations at all simplex points.
     *
     * \param[in] f the function to be evaluated
     * \param[in] initalGuess initial guess of the coordinates.
     *
     */
    NelderMeadSimplex(const std::function<real(ArrayRef<const real>)>& f, ArrayRef<const real> initalGuess);

    //! Return the vertex with the lowest function value at any of the simplex vertices.
    const RealFunctionvalueAtCoordinate& bestVertex() const;

    //! Return the vertex of the simplex with the highest (worst) function value.
    const RealFunctionvalueAtCoordinate& worstVertex() const;

    //! Return the second largest function value at any of the simplex vertices.
    real secondWorstValue() const;

    //! Return the reflection point and the evaluated function value at this point.
    RealFunctionvalueAtCoordinate
    evaluateReflectionPoint(const std::function<real(ArrayRef<const real>)>& f) const;

    //! Evaluate and return the expansion point and function value.
    RealFunctionvalueAtCoordinate evaluateExpansionPoint(const std::function<real(ArrayRef<const real>)>& f) const;

    //! Evaluate and return the contraction point and function value.
    RealFunctionvalueAtCoordinate
    evaluateContractionPoint(const std::function<real(ArrayRef<const real>)>& f) const;

    /*! \brief Replace the simplex vertex with the largest function value.
     *
     * \param[in] newVertex to replace the worst vertex with
     * \note keeps the simplex list sorted and reevaluates the reflection point
     */
    void swapOutWorst(const RealFunctionvalueAtCoordinate& newVertex);

    /*! \brief Shrink the simplex.
     *
     * All points move closer to the best point by a factor \f$\sigma\f$.
     *
     * Replace all point coordinates, except the best, with
     * \f$x_i = x_{\mathrm{best}} + \sigma (x_i - x_{\mathrm{best}})\f$
     */
    void shrinkSimplexPointsExceptBest(const std::function<real(ArrayRef<const real>)>& f);

    /*! \brief The oriented length of the vertex.
     *
     * The oriented length of the simplex is defined as the largest distance
     * between the first simplex vertex coordinate (with the lowest, best function
     * value) and any other simplex coordinate.
     *
     * The oriented length is used as a computationally fast and simple
     * convergence criterion because it is proven that
     * orientedLegnth < simplex_diameter < 2 * orientedLength
     *
     */
    real orientedLength() const;

private:
    /*! \brief Update centroid and reflection point.
     *
     * The arithmetic mean of all vertex coordinates expect the one at the
     * highest (worst) function value.
     *
     */
    void updateCentroidAndReflectionPoint();

    /*! \brief The points of the simplex with the function values.
     * \note This list stays sorted according to function value during the
     *       life-time of this object.
     */
    std::list<RealFunctionvalueAtCoordinate> simplex_;

    //! The centroid of the simplex, skipping the worst point is updated once the simplex changes
    std::vector<real> centroidWithoutWorstPoint_;

    //! The reflection point and its function value is updated once the simplex changes
    std::vector<real> reflectionPointCoordinates_;
};


} // namespace gmx

#endif // GMX_MATH_NELDERMEAD_H
