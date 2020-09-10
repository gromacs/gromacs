/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 *
 * \brief Implements routines in freeenergyparameters.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "freeenergyparameters.h"

#include "gromacs/mdtypes/inputrec.h"

namespace gmx
{

namespace
{

std::array<real, efptNR> lambdasAtState(const int stateIndex, double** const lambdaArray, const int lambdaArrayExtent)
{
    std::array<real, efptNR> lambda;
    // set lambda from an fep state index from stateIndex, if stateIndex was defined (> -1)
    if (stateIndex >= 0 && stateIndex < lambdaArrayExtent)
    {
        for (int i = 0; i < efptNR; i++)
        {
            lambda[i] = lambdaArray[i][stateIndex];
        }
    }
    return lambda;
}

/*! \brief Evaluates where in the lambda arrays we are at currently.
 *
 * \param[in] step the current simulation step
 * \param[in] deltaLambdaPerStep the change of the overall controlling lambda parameter per step
 * \param[in] initialFEPStateIndex the FEP state at the start of the simulation. -1 if not set.
 * \param[in] initialLambda the lambda at the start of the simulation . -1 if not set.
 * \param[in] lambdaArrayExtent how many lambda values we have
 * \returns a number that reflects where in the lambda arrays we are at the moment
 */
double currentGlobalLambda(const int64_t step,
                           const double  deltaLambdaPerStep,
                           const int     initialFEPStateIndex,
                           const double  initialLambda,
                           const int     lambdaArrayExtent)
{
    const real fracSimulationLambda = step * deltaLambdaPerStep;

    // Set initial lambda value for the simulation either from initialFEPStateIndex or,
    // if not set, from the initial lambda.
    double initialGlobalLambda = 0;
    if (initialFEPStateIndex > -1)
    {
        if (lambdaArrayExtent > 1)
        {
            initialGlobalLambda = static_cast<double>(initialFEPStateIndex) / (lambdaArrayExtent - 1);
        }
    }
    else
    {
        if (initialLambda > -1)
        {
            initialGlobalLambda = initialLambda;
        }
    }

    return initialGlobalLambda + fracSimulationLambda;
}

/*! \brief Returns an array of lambda values from linear interpolation of a lambda value matrix.
 *
 * \note If there is nothing to interpolate, fills the array with the global current lambda.
 * \note Returns array boundary values if currentGlobalLambda <0 or >1 .
 *
 * \param[in] currentGlobalLambda determines at which position in the lambda array to interpolate
 * \param[in] lambdaArray array of lambda values
 * \param[in] lambdaArrayExtent number of lambda values
 */
std::array<real, efptNR> interpolatedLambdas(const double   currentGlobalLambda,
                                             double** const lambdaArray,
                                             const int      lambdaArrayExtent)
{
    std::array<real, efptNR> lambda;
    // when there is no lambda value array, set all lambdas to steps * deltaLambdaPerStep
    if (lambdaArrayExtent <= 0)
    {
        std::fill(std::begin(lambda), std::end(lambda), currentGlobalLambda);
        return lambda;
    }

    // if we run over the boundary of the lambda array, return the boundary array values
    if (currentGlobalLambda <= 0)
    {
        for (int i = 0; i < efptNR; i++)
        {
            lambda[i] = lambdaArray[i][0];
        }
        return lambda;
    }
    if (currentGlobalLambda >= 1)
    {
        for (int i = 0; i < efptNR; i++)
        {
            lambda[i] = lambdaArray[i][lambdaArrayExtent - 1];
        }
        return lambda;
    }

    // find out between which two value lambda array elements to interpolate
    const int fepStateLeft = static_cast<int>(std::floor(currentGlobalLambda * (lambdaArrayExtent - 1)));
    const int fepStateRight = fepStateLeft + 1;
    // interpolate between this state and the next
    const double fracBetween = currentGlobalLambda * (lambdaArrayExtent - 1) - fepStateLeft;
    for (int i = 0; i < efptNR; i++)
    {
        lambda[i] = lambdaArray[i][fepStateLeft]
                    + fracBetween * (lambdaArray[i][fepStateRight] - lambdaArray[i][fepStateLeft]);
    }
    return lambda;
}

} // namespace

std::array<real, efptNR> currentLambdas(const int64_t step, const t_lambda& fepvals, const int currentLambdaState)
{
    if (fepvals.delta_lambda == 0)
    {
        if (currentLambdaState > -1)
        {
            return lambdasAtState(currentLambdaState, fepvals.all_lambda, fepvals.n_lambda);
        }

        if (fepvals.init_fep_state > -1)
        {
            return lambdasAtState(fepvals.init_fep_state, fepvals.all_lambda, fepvals.n_lambda);
        }

        std::array<real, efptNR> lambdas;
        std::fill(std::begin(lambdas), std::end(lambdas), fepvals.init_lambda);
        return lambdas;
    }
    const double globalLambda = currentGlobalLambda(step, fepvals.delta_lambda, fepvals.init_fep_state,
                                                    fepvals.init_lambda, fepvals.n_lambda);
    return interpolatedLambdas(globalLambda, fepvals.all_lambda, fepvals.n_lambda);
}

} // namespace gmx
