/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * \brief
 * Declares the DimParams struct and AWH vector types.
 *
 * This class holds the physical information for a dimension
 * of the bias reaction-coordinate grid.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_DIMPARAMS_H
#define GMX_AWH_DIMPARAMS_H

#include <variant>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

//! The maximum dimensionality of the AWH coordinate.
static const int c_biasMaxNumDim = 4;

//! A real vector in AWH coordinate space.
typedef double awh_dvec[c_biasMaxNumDim];

//! An integer vector in AWH coordinate space.
typedef int awh_ivec[c_biasMaxNumDim];

/*! \internal \brief Constant parameters for each dimension of the coordinate.
 */
struct DimParams
{
    /*! \internal \brief Type for storing dimension parameters for pull type dimensions
     */
    struct PullDimParams
    {
        const double k;     /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
        const double betak; /**< Inverse variance (1/nm^2) for each coordinate dimension. */
        const double userCoordUnitsToInternal; /**< Conversion factor coordinate units. */
    };

    /*! \internal \brief Type for storing dimension parameters for free-energy lambda type dimensions
     */
    struct FepDimParams
    {
        const double beta;               /**< 1/(k_B T). */
        const int    numFepLambdaStates; /**< Number of lambda points in this dimension. */
    };

private:
    /*! \brief
     * Private constructor called by public builder functions for PullDimParams and FepLambdaDimParams.
     */
    DimParams(double conversionFactor, std::variant<PullDimParams, FepDimParams> dimParams) :
        dimParams_(std::move(dimParams)), userCoordUnitsToInternal_(conversionFactor)
    {
    }

public:
    /*! \brief
     * Builder function for pull dimension parameters.
     *
     * \param[in] conversionFactor  Conversion factor from user coordinate units to bias internal
     * units (=c_deg2Rad for angles).
     * \param[in] forceConstant     The harmonic force constant.
     * \param[in] beta              1/(k_B T).
     */
    static DimParams pullDimParams(double conversionFactor, double forceConstant, double beta)
    {
        PullDimParams pullDimParams = { forceConstant, forceConstant * beta };

        return DimParams(conversionFactor, pullDimParams);
    }

    /*! \brief
     * Builder function for FEP lambda dimension parameters.
     *
     * \param[in] numFepLambdaStates  Number of lambda states in the system.
     * \param[in] beta                1/(k_B T).
     */
    static DimParams fepLambdaDimParams(int numFepLambdaStates, double beta)
    {
        FepDimParams fepDimParams = { beta, numFepLambdaStates };

        return DimParams(1.0, fepDimParams);
    }

    //! Returns whether this dimension is coupled to a pull coordinate.
    bool isPullDimension() const { return std::holds_alternative<PullDimParams>(dimParams_); }

    //! Returns whether this dimension has lambda states and thereby is a dimension coupled to lambda.
    bool isFepLambdaDimension() const { return std::holds_alternative<FepDimParams>(dimParams_); }

    //! Returns pull dimension parameters, only call for pull dimensions
    const PullDimParams& pullDimParams() const { return std::get<PullDimParams>(dimParams_); }

    //! Returns FEP dimension parameters, only call for FEP dimensions
    const FepDimParams& fepDimParams() const { return std::get<FepDimParams>(dimParams_); }

    /*! \brief Convert internal coordinate units to external, user coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleInternalToUserInput(double value) const
    {
        return value / userCoordUnitsToInternal_;
    }

    /*! \brief Convert external, user coordinate units to internal coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleUserInputToInternal(double value) const
    {
        return value * userCoordUnitsToInternal_;
    }

    //! Parameters for pull dimensions, either type pull or free-energy lambda
    const std::variant<PullDimParams, FepDimParams> dimParams_;
    //! Conversion factor for ordinate units
    const double userCoordUnitsToInternal_;
};

} // namespace gmx

#endif /* GMX_AWH_DIMPARAMS_H */
