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
 * \brief Definitions of the compound HDF5 data types for the H5MD topology.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "h5md_topologycompound.h"

#include <hdf5.h>

#include <algorithm>
#include <string>
#include <string_view>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_compound_type.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

MoleculeBlock::MoleculeBlock(std::string_view inputName, const int32_t inputCount) :
    count{ inputCount }
{
    if (inputName.size() >= sc_nameCapacity)
    {
        GMX_THROW(
                InvalidInputError(formatString("The name of molecule type %s is longer than the "
                                               "%zu characters that H5MD can store",
                                               std::string(inputName).c_str(),
                                               sc_nameCapacity - 1)));
    }
    std::copy(inputName.begin(), inputName.end(), name.begin());
}

namespace
{

//! \brief Return a molecule block data type of \p size bytes with its members at the given offsets.
hid_t createMoleculeBlockDataType(const size_t size,
                                  const size_t nameOffset,
                                  const size_t nameSize,
                                  const size_t countOffset)
{
    // The compound type holds a copy of the member types, so this one can be closed after
    const auto [nameDataType, nameDataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedLengthString(nameSize));
    return createCompoundDataType(
            size,
            { { MoleculeBlock::sc_names[0], nameOffset, nameDataType },
              { MoleculeBlock::sc_names[1], countOffset, hdf5DataTypeFor<int32_t>() } });
}

} // namespace

hid_t MoleculeBlock::constructHdf5DataType()
{
    return createMoleculeBlockDataType(sizeof(MoleculeBlock),
                                       HOFFSET(MoleculeBlock, name),
                                       sc_nameCapacity,
                                       HOFFSET(MoleculeBlock, count));
}

hid_t MoleculeBlock::constructHdf5DataType(const size_t nameSize)
{
    GMX_H5MD_THROW_UPON_ERROR(
            nameSize == 0 || nameSize > sc_nameCapacity,
            formatString("The names of molecule blocks must be stored with 1 to %zu bytes", sc_nameCapacity));
    return createMoleculeBlockDataType(nameSize + sizeof(int32_t), 0, nameSize, nameSize);
}

bool MoleculeBlock::matches(const hid_t dataType)
{
    if (H5Tget_class(dataType) != H5T_COMPOUND || H5Tget_nmembers(dataType) != sc_nelems
        || !checkCompoundMember<FixedLengthString>(dataType, 0, sc_names[0])
        || !checkCompoundMember<int32_t>(dataType, 1, sc_names[1]))
    {
        return false;
    }
    // A longer stored name would be cut when reading it into the struct
    const auto [nameDataType, nameDataTypeGuard] = makeH5mdTypeGuard(H5Tget_member_type(dataType, 0));
    return H5Tget_size(nameDataType) <= sc_nameCapacity;
}

} // namespace gmx
