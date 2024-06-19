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
/*! \internal \file
 * \brief
 * Implement mrc/ccp4-file metadata.
 *
 * \author Christian Blau <blau@kth.se>
 *
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "mrcdensitymapheader.h"

#include <cstdint>

#include <algorithm>
#include <iterator>

#include "gromacs/mdspan/extents.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

//! Returns true if any of the argument values is smaller than zero.
template<typename Container>
bool anySmallerZero(const Container& values)
{
    return std::any_of(std::begin(values), std::end(values), [](auto v) { return v < 0; });
}

//! Returns true if any of the argument values is larger than a given boundary value.
template<typename Container>
bool anyLargerThanValue(const Container& values, typename Container::value_type boundaryValue)
{
    return std::any_of(std::begin(values), std::end(values), [boundaryValue](auto v) {
        return v > boundaryValue;
    });
}

} // namespace

size_t numberOfExpectedDataItems(const MrcDensityMapHeader& header)
{
    if (anySmallerZero(header.numColumnRowSection_))
    {
        GMX_THROW(
                InternalError("Cannot determine data size, because the mrc "
                              "density map header is invalid (Negative number "
                              "describing data extent in at least one dimension)."));
    }

    return header.numColumnRowSection_[XX] * header.numColumnRowSection_[YY]
           * header.numColumnRowSection_[ZZ];
}

TranslateAndScale getCoordinateTransformationToLattice(const MrcDensityMapHeader& header)
{
    constexpr real c_AAtoNmConversion = 0.1;

    RVec       scale = { header.extent_[XX] / (header.cellLength_[XX] * c_AAtoNmConversion),
                   header.extent_[YY] / (header.cellLength_[YY] * c_AAtoNmConversion),
                   header.extent_[ZZ] / (header.cellLength_[ZZ] * c_AAtoNmConversion) };
    const RVec emdbOrigin{ header.userDefinedFloat_[12],
                           header.userDefinedFloat_[13],
                           header.userDefinedFloat_[14] };
    RVec       translation;
    if (emdbOrigin[XX] == 0. && emdbOrigin[YY] == 0. && emdbOrigin[ZZ] == 0.)
    {
        translation = RVec{ -header.columnRowSectionStart_[XX] / scale[XX],
                            -header.columnRowSectionStart_[YY] / scale[YY],
                            -header.columnRowSectionStart_[ZZ] / scale[ZZ] };
    }
    else
    {
        translation = { -emdbOrigin[XX] * c_AAtoNmConversion,
                        -emdbOrigin[YY] * c_AAtoNmConversion,
                        -emdbOrigin[ZZ] * c_AAtoNmConversion };
    }
    return { scale, translation };
}

dynamicExtents3D getDynamicExtents3D(const MrcDensityMapHeader& header)
{
    return { header.numColumnRowSection_[ZZ],
             header.numColumnRowSection_[YY],
             header.numColumnRowSection_[XX] };
};

bool mrcHeaderIsSane(const MrcDensityMapHeader& header)
{
    // Make sure all numbers of columns, row sections, extents and cell angles
    // are positive
    if (anySmallerZero(header.numColumnRowSection_) || anySmallerZero(header.cellAngles_)
        || anySmallerZero(header.extent_))
    {
        return false;
    }

    // The maximum integer number in an mrc header to be considered sane
    constexpr std::int32_t c_maxIntegerNumber = 100'000;
    if (anyLargerThanValue(header.numColumnRowSection_, c_maxIntegerNumber)
        || anyLargerThanValue(header.extent_, c_maxIntegerNumber))
    {
        return false;
    }

    constexpr std::int32_t c_maxCellAngle = 360;
    if (anyLargerThanValue(header.cellAngles_, c_maxCellAngle))
    {
        return false; //NOLINT(readability-simplify-boolean-expr)
    }

    return true;
}

} // namespace gmx
