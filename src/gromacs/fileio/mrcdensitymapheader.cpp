/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Implements mrc file.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "mrcdensitymapheader.h"

namespace gmx
{

namespace
{
//! Default values for mrc headers.
struct MrcHeaderDefaults
{
    //! format identifier
    static constexpr std::array<char, MrcDensityMapHeader::c_formatIdentifierLength>
    c_formatIdentifier = {{'M', 'A', 'P', ' '}};

    //! cell length
    static constexpr std::array<float, 3> c_cellLength = {{1, 1, 1}};

    //! space group
    static constexpr SpaceGroup c_spaceGroup = SpaceGroup::P1;

    //! data mode
    static constexpr MrcDataMode c_dataMode = MrcDataMode::float32;

    //! cell angles
    static constexpr std::array<float, 3> c_cellAngles = {{90, 90, 90}};

    //! crs to xyz mapping
    static constexpr std::array<MrcDensityMapHeader::mrc_integer, 3> c_crsToXyz = {{0, 1, 2}};
    //! crs to xyz mapping
    static constexpr std::array<MrcDensityMapHeader::mrc_integer, 3> c_crsStart = {{0, 0, 0}};
};

}   // namespace

MrcDensityMapHeader::MrcDensityMapHeader() : spaceGroup_(MrcHeaderDefaults::c_spaceGroup),
                                             dataMode_(MrcHeaderDefaults::c_dataMode),
#ifdef GMX_INTEGER_BIG_ENDIAN
                                             machineStamp_(MachineStamp::bigEndian),
#else
                                             machineStamp_(MachineStamp::smallEndian),
#endif
                                             formatIdentifier_(MrcHeaderDefaults::c_formatIdentifier),
                                             userDefinedFloat_({}
                                                               ),
                                             numUsedLabels_(1),
                                             cellLength_(MrcHeaderDefaults::c_cellLength),
                                             cellAngles_(MrcHeaderDefaults::c_cellAngles),
                                             crsToXyz_(MrcHeaderDefaults::c_crsToXyz),
                                             numCrs_({}),
                                             crsStart_(MrcHeaderDefaults::c_crsStart),
                                             extent_({}),
                                             dataStats_({}),
                                             extendedHeader_({})
{
    skewData_.valid_       = false;
    skewData_.matrix_      = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
    skewData_.translation_ = {{0, 0, 0}};

    static constexpr char c_emdbCustomLabel[MrcDensityMapHeader::c_labelSize + 1] =
        "::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention        ::::";
    labels_[0] = std::string(c_emdbCustomLabel, MrcDensityMapHeader::c_labelSize);
    for (int i = 1; i < c_numCrystallographicLabels; ++i)
    {
        labels_[i] = std::string(MrcDensityMapHeader::c_labelSize, ' ');
    }
}

} // namespace gmx
