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
 * \brief Implements routines in optimization.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/math/optimization.h"

#include <functional>

#include "gromacs/math/neldermead.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{

OptimisationResult nelderMead(const std::function<real(ArrayRef<const real>)>& functionToMinimize,
                              ArrayRef<const real>                             initalGuess,
                              real minimumRelativeSimplexLength,
                              int  maxSteps)
{
    // Set up the initial simplex, sorting vertices according to function value
    NelderMeadSimplex nelderMeadSimplex(functionToMinimize, initalGuess);

    // Run until maximum step size reached or algorithm is converged, e.g.,
    // the oriented simplex length is smaller or equal a given number.
    const real minimumSimplexLength = minimumRelativeSimplexLength * nelderMeadSimplex.orientedLength();
    for (int currentStep = 0;
         nelderMeadSimplex.orientedLength() > minimumSimplexLength && currentStep < maxSteps;
         ++currentStep)
    {

        // see if simplex can by improved by reflecing the worst vertex at the centroid
        const RealFunctionvalueAtCoordinate& reflectionPoint =
                nelderMeadSimplex.evaluateReflectionPoint(functionToMinimize);

        // Reflection point is not better than best simplex vertex so far
        // but better than second worst
        if ((nelderMeadSimplex.bestVertex().value_ <= reflectionPoint.value_)
            && (reflectionPoint.value_ < nelderMeadSimplex.secondWorstValue()))
        {
            nelderMeadSimplex.swapOutWorst(reflectionPoint);
            continue;
        }

        // If the reflection point is better than the best one see if simplex
        // can be further improved by continuing going in that direction
        if (reflectionPoint.value_ < nelderMeadSimplex.bestVertex().value_)
        {
            RealFunctionvalueAtCoordinate expansionPoint =
                    nelderMeadSimplex.evaluateExpansionPoint(functionToMinimize);
            if (expansionPoint.value_ < reflectionPoint.value_)
            {
                nelderMeadSimplex.swapOutWorst(expansionPoint);
            }
            else
            {
                nelderMeadSimplex.swapOutWorst(reflectionPoint);
            }
            continue;
        }

        // The reflection point was a poor choice, try contracting the
        // worst point coordinates using the centroid instead
        RealFunctionvalueAtCoordinate contractionPoint =
                nelderMeadSimplex.evaluateContractionPoint(functionToMinimize);
        if (contractionPoint.value_ < nelderMeadSimplex.worstVertex().value_)
        {
            nelderMeadSimplex.swapOutWorst(contractionPoint);
            continue;
        }

        // If neither expansion nor contraction of the worst point give a
        // good result shrink the whole simplex
        nelderMeadSimplex.shrinkSimplexPointsExceptBest(functionToMinimize);
    }

    return { nelderMeadSimplex.bestVertex().coordinate_, nelderMeadSimplex.bestVertex().value_ };
}

} // namespace gmx
