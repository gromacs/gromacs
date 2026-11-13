/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Definitions of utilities for H5MD compound data types.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "h5md_compound_type.h"

#include <hdf5.h>

#include <cstdint>
#include <cstring>

#include <type_traits>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

hid_t createCompoundDataType(const size_t size, const std::initializer_list<CompoundMember> members)
{
    auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Tcreate(H5T_COMPOUND, size));
    GMX_H5MD_THROW_UPON_INVALID_HID(dataType, "Could not create a compound data type");
    for (const CompoundMember& member : members)
    {
        GMX_H5MD_THROW_UPON_ERROR(
                H5Tinsert(dataType, member.name, member.offset, member.dataType) < 0,
                formatString("Could not add member %s to a compound data type", member.name));
    }
    // The caller takes ownership of the data type
    dataTypeGuard.dismiss();
    return dataType;
}

template<typename ValueType>
bool checkCompoundMember(const hid_t compoundDataType, const int memberIndex, const char* expectedName)
{
    const auto [memberDataType, memberDataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_member_type(compoundDataType, memberIndex));
    GMX_H5MD_THROW_UPON_INVALID_HID(memberDataType,
                                    "Could not access member type from compound data type");

    char* memberName = H5Tget_member_name(compoundDataType, memberIndex);
    GMX_H5MD_THROW_UPON_ERROR(
            !memberName,
            formatString("Could not access compound type member variable %s at %d", expectedName, memberIndex));
    // HDF5 wants us to free allocated memory using its own functions,
    // use a scope guard to do so when we return
    const auto memberNameGuard = sg::make_scope_guard([&]() { H5free_memory(memberName); });

    if (strcmp(expectedName, memberName) != 0)
    {
        return false;
    }
    if constexpr (std::is_same_v<ValueType, FixedLengthString>)
    {
        return H5Tget_class(memberDataType) == H5T_STRING && H5Tis_variable_str(memberDataType) == 0;
    }
    else if constexpr (std::is_arithmetic_v<ValueType>)
    {
        return valueTypeIsSameClassAsDataType<ValueType>(memberDataType);
    }
    else
    {
        return valueTypeIsDataType<ValueType>(memberDataType);
    }
}

template bool checkCompoundMember<int32_t>(hid_t, int, const char*);
template bool checkCompoundMember<int64_t>(hid_t, int, const char*);
template bool checkCompoundMember<float>(hid_t, int, const char*);
template bool checkCompoundMember<double>(hid_t, int, const char*);
template bool checkCompoundMember<FixedLengthString>(hid_t, int, const char*);

} // namespace gmx
