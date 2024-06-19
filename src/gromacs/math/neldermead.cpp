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
 * \brief Implements routines in neldermead.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/math/neldermead.h"

#include <cmath>

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <numeric>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

/*! \brief Evaluate the linear combination of two vectors a and b.
 *
 * \param[in] alpha scaling factor for a
 * \param[in] a vector to be scaled
 * \param[in] beta scaling factor for b
 * \param[in] b vector to be scaled
 *
 * \returns alpha * a + beta * b.
 */
std::vector<real> linearCombination(real alpha, ArrayRef<const real> a, real beta, ArrayRef<const real> b)
{
    GMX_ASSERT(a.size() == b.size(),
               "Input vectors have to have the same size to evaluate their linear combination.");
    std::vector<real> result(a.size());
    std::transform(std::begin(a),
                   std::end(a),
                   std::begin(b),
                   std::begin(result),
                   [alpha, beta](auto elemA, auto elemB) { return alpha * elemA + beta * elemB; });
    return result;
}

/*! \internal
 *  \brief The parameters for a Nelder-Mead optimisation.
 */
struct NelderMeadParameters
{
    //! Factor to evaluate the reflection point
    real alpha_ = 1;
    //! Factor to evaluate the expansion point
    real gamma_ = 2;
    //! Factor to evaluate the contraction point
    real rho_ = 0.5;
    //! Factor to evaluate the simplex shrinkage
    real sigma_ = 0.5;
};

constexpr NelderMeadParameters defaultNelderMeadParameters = { 1, 2, 0.5, 0.5 };

} // namespace


NelderMeadSimplex::NelderMeadSimplex(const std::function<real(ArrayRef<const real>)>& f,
                                     ArrayRef<const real>                             initalGuess)
{
    // initial simplex contains the initally guessed vertex
    std::vector<real> initalVertex = copyOf(initalGuess);
    simplex_.push_back({ initalVertex, f(initalVertex) });
    // create the missing verticies by moving 0.05 or 0.0025 if null
    // from the initial vertex dimension after dimension
    for (auto& v : initalVertex)
    {
        const auto oldValue = v;
        if (v == 0)
        {
            v = 0.0025;
        }
        else
        {
            v += 0.05;
        }
        simplex_.push_back({ initalVertex, f(initalVertex) });
        v = oldValue;
    }
    simplex_.sort([](const RealFunctionvalueAtCoordinate& lhs,
                     const RealFunctionvalueAtCoordinate& rhs) { return lhs.value_ < rhs.value_; });
    updateCentroidAndReflectionPoint();
}

RealFunctionvalueAtCoordinate
NelderMeadSimplex::evaluateReflectionPoint(const std::function<real(ArrayRef<const real>)>& f) const
{
    return { reflectionPointCoordinates_, f(reflectionPointCoordinates_) };
}

const RealFunctionvalueAtCoordinate& NelderMeadSimplex::bestVertex() const
{
    return simplex_.front();
}

const RealFunctionvalueAtCoordinate& NelderMeadSimplex::worstVertex() const
{
    return simplex_.back();
}

real NelderMeadSimplex::secondWorstValue() const
{
    // go backwards one step from the end of the sorted simplex list
    // and look at the vertex value
    return std::next(std::rbegin(simplex_))->value_;
}

RealFunctionvalueAtCoordinate
NelderMeadSimplex::evaluateExpansionPoint(const std::function<real(ArrayRef<const real>)>& f) const
{
    const std::vector<real> expansionPointCoordinate =
            linearCombination(1 - defaultNelderMeadParameters.gamma_,
                              centroidWithoutWorstPoint_,
                              defaultNelderMeadParameters.gamma_,
                              reflectionPointCoordinates_);
    return { expansionPointCoordinate, f(expansionPointCoordinate) };
}

RealFunctionvalueAtCoordinate
NelderMeadSimplex::evaluateContractionPoint(const std::function<real(ArrayRef<const real>)>& f) const
{
    std::vector<real> contractionPoint = linearCombination(1 - defaultNelderMeadParameters.rho_,
                                                           centroidWithoutWorstPoint_,
                                                           defaultNelderMeadParameters.rho_,
                                                           worstVertex().coordinate_);
    return { contractionPoint, f(contractionPoint) };
}

void NelderMeadSimplex::swapOutWorst(const RealFunctionvalueAtCoordinate& newVertex)
{
    // drop the worst point - we know it's at the back of the simplex list, because
    // we kept the list sorted
    simplex_.pop_back();
    // find the point to insert the new vertex, so that the simplex vertices
    // keep being sorted according to function value
    const auto insertionPoint = std::lower_bound(
            std::begin(simplex_),
            std::end(simplex_),
            newVertex.value_,
            [](const RealFunctionvalueAtCoordinate& lhs, real value) { return lhs.value_ < value; });
    simplex_.insert(insertionPoint, newVertex);
    // now that the simplex has changed, it has a new centroid and reflection point
    updateCentroidAndReflectionPoint();
}

void NelderMeadSimplex::shrinkSimplexPointsExceptBest(const std::function<real(ArrayRef<const real>)>& f)
{
    std::vector<real> bestPointCoordinate = simplex_.front().coordinate_;
    // skipping over the first simplex vertex, pull points closer to the best
    // vertex
    std::transform(std::next(std::begin(simplex_)),
                   std::end(simplex_),
                   std::next(std::begin(simplex_)),
                   [bestPointCoordinate, f](const RealFunctionvalueAtCoordinate& d) -> RealFunctionvalueAtCoordinate {
                       const std::vector<real> shrinkPoint =
                               linearCombination(defaultNelderMeadParameters.sigma_,
                                                 d.coordinate_,
                                                 1 - defaultNelderMeadParameters.sigma_,
                                                 bestPointCoordinate);
                       return { shrinkPoint, f(shrinkPoint) };
                   });

    simplex_.sort([](const RealFunctionvalueAtCoordinate& lhs,
                     const RealFunctionvalueAtCoordinate& rhs) { return lhs.value_ < rhs.value_; });

    // now that the simplex has changed, it has a new centroid and reflection point
    updateCentroidAndReflectionPoint();
}

real NelderMeadSimplex::orientedLength() const
{
    real                    result                       = 0;
    const std::vector<real> firstSimplexVertexCoordinate = simplex_.front().coordinate_;
    // find out which vertex coordinate has the largest distance to the first simplex vertex.
    for (const auto& simplexVertex : simplex_)
    {
        const std::vector<real> differenceVector =
                linearCombination(1, firstSimplexVertexCoordinate, -1, simplexVertex.coordinate_);
        const real thisLength = std::accumulate(
                std::begin(differenceVector), std::end(differenceVector), 0., [](real sum, real value) {
                    return sum + value * value;
                });
        result = std::max(result, thisLength);
    }
    return std::sqrt(result);
}

void NelderMeadSimplex::updateCentroidAndReflectionPoint()
{
    // initialize with first vertex, then add up all other vertex coordinates
    // expect last one
    centroidWithoutWorstPoint_ =
            std::accumulate(std::next(std::begin(simplex_)),
                            std::prev(std::end(simplex_)),
                            simplex_.front().coordinate_,
                            [](std::vector<real> sum, const RealFunctionvalueAtCoordinate& x) {
                                std::transform(std::begin(sum),
                                               std::end(sum),
                                               std::begin(x.coordinate_),
                                               std::begin(sum),
                                               std::plus<>());
                                return sum;
                            });

    // divide the summed up coordinates by N (the simplex has N+1 vertices)
    std::transform(std::begin(centroidWithoutWorstPoint_),
                   std::end(centroidWithoutWorstPoint_),
                   std::begin(centroidWithoutWorstPoint_),
                   [n = simplex_.size() - 1](const auto& x) { return x / n; });

    // now, that we have evaluated the centroid, update the reflection points
    reflectionPointCoordinates_ = linearCombination(
            defaultNelderMeadParameters.alpha_ + 1, centroidWithoutWorstPoint_, -1, worstVertex().coordinate_);
}

} // namespace gmx
