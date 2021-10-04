/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "transformationcoordinate.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "pull_internal.h"
#include "pullcoordexpressionparser.h"

namespace gmx
{

namespace
{

//! Calculates the value a for transformation pull coordinate
double getTransformationPullCoordinateValue(pull_coord_work_t* coord)
{
    const int transformationPullCoordinateIndex = coord->params.coordIndex;
    GMX_ASSERT(ssize(coord->transformationVariables) == transformationPullCoordinateIndex,
               "We need as many variables as the transformation pull coordinate index");
    double result = 0;
    try
    {
        result = coord->expressionParser.evaluate(coord->transformationVariables);
    }
#if HAVE_MUPARSER
    catch (mu::Parser::exception_type& e)
    {
        GMX_THROW(InconsistentInputError(
                formatString("failed to evaluate expression for transformation pull-coord%d: %s\n",
                             transformationPullCoordinateIndex + 1,
                             e.GetMsg().c_str())));
    }
#endif
    catch (std::exception& e)
    {
        GMX_THROW(InconsistentInputError(
                formatString("failed to evaluate expression for transformation pull-coord%d.\n"
                             "Last variable pull-coord-index: %d.\n"
                             "Message:  %s\n",
                             transformationPullCoordinateIndex + 1,
                             transformationPullCoordinateIndex + 1,
                             e.what())));
    }
    return result;
}

} // namespace

double getTransformationPullCoordinateValue(pull_coord_work_t*                coord,
                                            ArrayRef<const pull_coord_work_t> variableCoords)
{
    GMX_ASSERT(ssize(variableCoords) == coord->params.coordIndex,
               "We need as many variables as the transformation pull coordinate index");
    int coordIndex = 0;
    for (const auto& variableCoord : variableCoords)
    {
        coord->transformationVariables[coordIndex++] = variableCoord.spatialData.value;
    }

    return getTransformationPullCoordinateValue(coord);
}

/*! \brief Calculates and returns the derivative of a transformation pull coordinate from a dependent coordinate
 *
 * Note #1: this requires that getTransformationPullCoordinateValue() has been called
 * before with the current coordinates.
 *
 * Note #2: this method will not compute inner derivates. That is taken care of in the regular pull code
 *
 * \param[in] coord  The (transformation) coordinate to compute the value for
 * \param[in] variablePcrdIndex Pull coordinate index of a variable.
 */
static double computeDerivativeForTransformationPullCoord(pull_coord_work_t* coord, const int variablePcrdIndex)
{
    GMX_ASSERT(variablePcrdIndex >= 0 && variablePcrdIndex < coord->params.coordIndex,
               "The variable index should be in range of the transformation coordinate");

    // epsilon for numerical differentiation.
    const double transformationPcrdValue = coord->spatialData.value;
    // Perform numerical differentiation of 1st order
    const double valueBackup = coord->transformationVariables[variablePcrdIndex];
    double       dx          = coord->params.dx;
    coord->transformationVariables[variablePcrdIndex] += dx;
    double transformationPcrdValueEps = getTransformationPullCoordinateValue(coord);
    double derivative                 = (transformationPcrdValueEps - transformationPcrdValue) / dx;
    // reset pull coordinate value
    coord->transformationVariables[variablePcrdIndex] = valueBackup;
    return derivative;
}


/*!
 * \brief Distributes the force from a transformation pull coordiante to the dependent pull
 * coordinates by computing the inner derivatives
 *
 * \param pcrd The transformation pull coord
 * \param variableCoords The dependent pull coordinates
 * \param transformationCoordForce The force to distribute
 */
static void distributeTransformationPullCoordForce(pull_coord_work_t*               pcrd,
                                                   gmx::ArrayRef<pull_coord_work_t> variableCoords,
                                                   const double transformationCoordForce)
{
    if (std::abs(transformationCoordForce) < 1e-9)
    {
        // the force is effectively 0. Don't proceed and distribute it recursively
        return;
    }
    GMX_ASSERT(pcrd->params.eGeom == PullGroupGeometry::Transformation,
               "We shouldn't end up here when not using a transformation pull coordinate.");
    GMX_ASSERT(ssize(variableCoords) == pcrd->params.coordIndex,
               "We should have as many variable coords as the coord index of the transformation "
               "coordinate");

    // pcrd->scalarForce += transformationCoordForce;
    for (auto& variableCoord : variableCoords)
    {
        const double derivative =
                computeDerivativeForTransformationPullCoord(pcrd, variableCoord.params.coordIndex);
        const double variablePcrdForce = transformationCoordForce * derivative;
        /* Since we loop over all pull coordinates with smaller index, there can be ones
         * that are not referenced by the transformation coordinate. Avoid apply forces
         * on those by skipping application of zero force.
         */
        if (variablePcrdForce != 0)
        {
            if (debug)
            {
                fprintf(debug,
                        "Distributing force %4.4f for transformation coordinate %d to coordinate "
                        "%d with "
                        "force "
                        "%4.4f\n",
                        transformationCoordForce,
                        pcrd->params.coordIndex,
                        variableCoord.params.coordIndex,
                        variablePcrdForce);
            }
            // Note that we add to the force here, in case multiple biases act on the same pull
            // coord (although that is not recommended it should still work)
            variableCoord.scalarForce += variablePcrdForce;
            if (variableCoord.params.eGeom == PullGroupGeometry::Transformation)
            {
                /*
                 * We can have a transformation pull coordinate depend on another transformation pull coordinate
                 * which in turn leads to inner derivatives between pull coordinates.
                 * Here we redistribute the force via the inner product
                 *
                 * Note that this only works properly if the lower ranked transformation pull coordinate has it's scalarForce set to zero
                 */
                distributeTransformationPullCoordForce(
                        &variableCoord,
                        variableCoords.subArray(0, variableCoord.params.coordIndex),
                        variablePcrdForce);
            }
        }
    }
}

void applyTransformationPullCoordForce(pull_coord_work_t*               pcrd,
                                       gmx::ArrayRef<pull_coord_work_t> variableCoords,
                                       const double                     transformationCoordForce)
{
    pcrd->scalarForce = transformationCoordForce;
    // Note on why we need to call another method here:
    // applyTransformationPullCoordForce is the method that should be called by the rest of the pull code.
    // It's non-recursive and called exactly once for every transformation coordinate for every timestep.
    // In it, we set the force on the transformation coordinate,
    // then pass the force on to the other pull coordinates via the method distributeTransformationPullCoordForce.
    // The latter method is recursive to account for inner derivatives.
    // Note that we don't set the force on the top-level transformation coordinate in distributeTransformationPullCoordForce,
    // we only add to the force, which is why it can be recursive.
    distributeTransformationPullCoordForce(pcrd, variableCoords, transformationCoordForce);
}


} // namespace gmx
