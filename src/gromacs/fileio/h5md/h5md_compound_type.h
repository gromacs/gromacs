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
 * \brief Declarations of utilities for H5MD compound data types.
 *
 * A compound type is a struct with
 * - `static hid_t constructHdf5DataType()`, which returns the HDF5 data type of the struct,
 * - `static bool matches(hid_t dataType)`, which returns whether a stored data type can be
 *   read into the struct, and
 * - `static constexpr char sc_typeName[]`, its name for messages.
 *
 * A struct with fixed-length string members may also have
 * `static hid_t constructHdf5DataType(size_t stringSize)`, which returns the data type to store
 * with those members sized to \c stringSize bytes. H5mdDataSetBuilder uses it when given a
 * maximum string length, and HDF5 converts the members between the two layouts.
 *
 * A struct whose data set stores physical quantities also has
 * `static constexpr std::array<const char*, N> sc_units` with one unit per member, in the
 * order of `sc_names`, and an empty string for a dimensionless member. H5mdDataSetBuilder
 * writes it as the `units` attribute of the data set if and only if the struct declares
 * `sc_units`, whatever its values: a struct whose members are all dimensionless still gets
 * `units` (as the specification requires for `coulomb_14`), and a struct without `sc_units`
 * (such as MoleculeBlock) gets none.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_H5MD_COMPOUND_TYPE_H
#define GMX_FILEIO_H5MD_COMPOUND_TYPE_H

#include <hdf5.h>

#include <cstddef>

#include <initializer_list>
#include <type_traits>
#include <utility>

namespace gmx
{

namespace internal
{

//! \brief Whether \p ValueType has a static constructHdf5DataType().
template<typename ValueType, typename = void>
struct HasConstructHdf5DataType : std::false_type
{
};

//! \copydoc HasConstructHdf5DataType
template<typename ValueType>
struct HasConstructHdf5DataType<ValueType, std::void_t<decltype(ValueType::constructHdf5DataType())>> :
    std::true_type
{
};

//! \brief Whether \p ValueType declares the units of its members in sc_units.
template<typename ValueType, typename = void>
struct HasUnits : std::false_type
{
};

//! \copydoc HasUnits
template<typename ValueType>
struct HasUnits<ValueType, std::void_t<decltype(ValueType::sc_units)>> : std::true_type
{
};

//! \brief Whether \p ValueType has a static constructHdf5DataType(size_t) for sized strings.
template<typename ValueType, typename = void>
struct HasSizedConstructHdf5DataType : std::false_type
{
};

//! \copydoc HasSizedConstructHdf5DataType
template<typename ValueType>
struct HasSizedConstructHdf5DataType<ValueType, std::void_t<decltype(ValueType::constructHdf5DataType(std::declval<size_t>()))>> :
    std::true_type
{
};

} // namespace internal

/*! \brief Return whether the templated \p ValueType is a compound type for HDF5.
 *
 * \todo Replace this with a concept once we are in c++20.
 */
template<typename ValueType>
constexpr bool valueTypeIsHdf5CompoundType() noexcept
{
    return internal::HasConstructHdf5DataType<ValueType>::value;
}

/*! \brief Return whether the compound type \p ValueType can be stored with sized strings.
 *
 * \todo Replace this with a concept once we are in c++20.
 */
template<typename ValueType>
constexpr bool valueTypeHasSizedStrings() noexcept
{
    return internal::HasSizedConstructHdf5DataType<ValueType>::value;
}

/*! \brief Return whether the compound type \p ValueType declares the units of its members.
 *
 * If it does, H5mdDataSetBuilder writes them as the `units` attribute of the data set.
 *
 * \todo Replace this with a concept once we are in c++20.
 */
template<typename ValueType>
constexpr bool valueTypeHasUnits() noexcept
{
    return internal::HasUnits<ValueType>::value;
}

//! \brief Member type for checkCompoundMember(): a fixed-length string of any size.
struct FixedLengthString
{
};

//! \brief Name, offset in its struct and data type of a member of a compound type.
struct CompoundMember
{
    //! Name of the member in the HDF5 compound type
    const char* name;
    //! Offset of the member in its struct, from HOFFSET()
    size_t offset;
    //! HDF5 data type of the member, not owned
    hid_t dataType;
};

/*! \brief Create an HDF5 compound data type for a struct of \p size bytes with \p members.
 *
 * \note The returned handle must be closed with a call to H5Tclose to avoid leaking resources.
 *
 * \throws H5mdError if the data type cannot be created.
 */
hid_t createCompoundDataType(size_t size, std::initializer_list<CompoundMember> members);

/*! \brief Return whether a compound type member is correct.
 *
 * \tparam    ValueType        Type that the compound member should match.
 * \param[in] compoundDataType Compound data type.
 * \param[in] memberIndex      Index of member to check.
 * \param[in] expectedName     Expected name of member.
 *
 * An arithmetic member matches a stored member of the same type class (integer or floating
 * point) of any width, which HDF5 converts when reading, since the specification defines only
 * the class. A FixedLengthString member matches a stored fixed-length string of any size.
 *
 * Instantiated for int32_t, int64_t, float, double and FixedLengthString.
 *
 * \returns True if the member type and name is correct, otherwise false.
 *
 * \throws H5mdError if an error occurred during any type checking operations.
 */
template<typename ValueType>
bool checkCompoundMember(hid_t compoundDataType, int memberIndex, const char* expectedName);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_COMPOUND_TYPE_H
