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
 * \brief Compound HDF5 data types for the H5MD topology.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_H5MD_TOPOLOGYCOMPOUND_H
#define GMX_FILEIO_H5MD_TOPOLOGYCOMPOUND_H

#include <hdf5.h>

#include <cstdint>

#include <array>
#include <string_view>

#include "gromacs/utility/cstringutil.h"

namespace gmx
{

/*! \brief Molecule block entry data.
 *
 * In memory the name has a fixed capacity that holds any molecule type name, and is stored
 * as a fixed-length string. Give H5mdDataSetBuilder::withMaxStringLength() the length of the
 * longest name plus one to store the names in no more space than they need.
 *
 * A fixed-length string is used because it is stored, compressed and read as plain bytes of
 * the record, without the per-string heap allocations of a variable-length string that every
 * reader must free and that many tools, such as NumPy structured arrays in h5py, cannot map.
 */
struct MoleculeBlock
{
    //! \brief Name of this struct, for messages
    static constexpr char sc_typeName[] = "MoleculeBlock";
    //! \brief Number of elements in this struct
    static constexpr int sc_nelems = 2;
    //! \brief Names for elements in this struct to store in the HDF5 compound type
    static constexpr std::array<const char*, sc_nelems> sc_names = { "name", "count" };
    //! \brief Capacity of \c name including the terminating null, as grompp reads the names
    static constexpr size_t sc_nameCapacity = STRLEN;

    /*! \brief Construct the HDF5 data type of the struct in memory and return its handle.
     *
     * \note The returned handle must be closed with a call to H5Tclose to avoid leaking resources.
     *
     * \throws H5mdError if the data type cannot be created.
     */
    static hid_t constructHdf5DataType();
    /*! \brief Construct the HDF5 data type that stores names of \p nameSize bytes and return
     * its handle.
     *
     * \p nameSize includes the terminating null. The members are packed without padding.
     *
     * \note The returned handle must be closed with a call to H5Tclose to avoid leaking resources.
     *
     * \throws H5mdError if \p nameSize is zero or larger than sc_nameCapacity, or if the data
     * type cannot be created.
     */
    static hid_t constructHdf5DataType(size_t nameSize);
    /*! \brief Return true if the HDF5 \p dataType matches this struct.
     *
     * The stored name may have any size up to sc_nameCapacity.
     *
     * \throws H5mdError if the members of \p dataType cannot be read.
     */
    static bool matches(hid_t dataType);

    //! \brief Default constructor, with an empty name and a count of zero.
    MoleculeBlock() = default;
    /*! \brief Constructor, copies \p inputName into \c name and sets \c count.
     *
     * \throws InvalidInputError if \p inputName does not fit in sc_nameCapacity.
     */
    MoleculeBlock(std::string_view inputName, int32_t inputCount);

    std::array<char, sc_nameCapacity> name{};    //!< Name of the molecule type, null-terminated
    int32_t                           count = 0; //!< Number of molecules in block
};

} // namespace gmx

#endif // GMX_FILEIO_H5MD_TOPOLOGYCOMPOUND_H
