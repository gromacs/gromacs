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

#include <set>

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

//! Check if cell angles (in deg) have valid values for a crystallographic cell.
bool validUnitCellAngles(const std::array<float, 3> &angles)
{
    bool valid = true;
    // avoid extremely acute angled cells
    if (angles[0] * angles[1] * angles[2] < 1.f)
    {
        valid = false;
    }
    const auto angleIsConvex = [](float angle) {
            return angle > 0.f && angle < 180.f;
        };
    if (!std::all_of(std::begin(angles), std::end(angles), angleIsConvex))
    {
        valid = false;
    }
    return valid;
}

//! True if map contains all elements, 0,1,2, which stand for x-, y- and z-dimension.
bool bijectiveMappingToXyz(std::array<int, 3> map)
{
    return std::set<int>{
               std::begin(map), std::end(map)
    } == std::set<int>{
               0, 1, 2
    };
}

}   // namespace

MrcDensityMapHeader::MrcDensityMapHeader() : spaceGroup(SpaceGroup::P1),
                                             dataMode(MrcDataMode::float32),
#ifdef GMX_INTEGER_BIG_ENDIAN
                                             machineStamp(MachineStamp::bigEndian),
#else
                                             machineStamp(MachineStamp::smallEndian),
#endif
                                             formatIdentifier(MrcHeaderDefaults::c_formatIdentifier),
                                             userDefinedFloat({}
                                                              ),
                                             numUsedLabels(1),
                                             cellLength(MrcHeaderDefaults::c_cellLength),
                                             cellAngles(MrcHeaderDefaults::c_cellAngles),
                                             crsToXyz(MrcHeaderDefaults::c_crsToXyz),
                                             numCrs({}),
                                             crsStart(MrcHeaderDefaults::c_crsStart),
                                             extent({}),
                                             dataStats({}),
                                             extendedHeader({})
{
    skewData.valid       = false;
    skewData.matrix      = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
    skewData.translation = {{0, 0, 0}};

    static constexpr char c_emdbCustomLabel[MrcDensityMapHeader::c_labelSize + 1] =
        "::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention        ::::";
    labels[0] = std::string(c_emdbCustomLabel, MrcDensityMapHeader::c_labelSize);
    for (int i = 1; i < c_numCrystallographicLabels; ++i)
    {
        labels[i] = std::string(MrcDensityMapHeader::c_labelSize, ' ');
    }
}

MrcDiagnostic mrcDiagnosticAndFixes(const MrcDensityMapHeader &mrcFile)
{
    MrcDiagnostic diagnostic;

    diagnostic.fixedHeader = mrcFile;

    std::string diagnosticMessage = {};
    if (mrcFile.dataMode != MrcDataMode::float32)
    {
        diagnostic.report += "Expected data mode 2, corresponding to float32, but "
            " have an unsupported one.\n";

        diagnostic.fixedHeader.dataMode = MrcHeaderDefaults::c_dataMode;
    }

    if (mrcFile.extendedHeader.size() % MrcDensityMapHeader::c_labelSize != 0)
    {
        diagnostic.report += "Expected multiple of "
            + std::to_string(MrcDensityMapHeader::c_labelSize)
            + " bytes in symbol table , but have "
            + std::to_string(mrcFile.extendedHeader.size()) + ".\n";

        diagnostic.fixedHeader.extendedHeader.resize(
                (mrcFile.extendedHeader.size() / MrcDensityMapHeader::c_labelSize + 1)
                * MrcDensityMapHeader::c_labelSize);
    }

    for (const auto &label : mrcFile.labels)
    {
        if (label.size() != MrcDensityMapHeader::c_labelSize)
        {
            diagnostic.report += "Expected crystallographic label length "
                + std::to_string(MrcDensityMapHeader::c_labelSize)
                + ", but is " + std::to_string(label.size()) + ".\n";
        }
    }

    if (mrcFile.formatIdentifier != decltype(mrcFile.formatIdentifier){
            {'M', 'A', 'P', ' '}
        })
    {
        diagnostic.report += "Expected \'MAP \' as format identifier.\n";

        diagnostic.fixedHeader.formatIdentifier = MrcHeaderDefaults::c_formatIdentifier;
    }

    // By convention, unset cell angles (all 0) are interpreted as 90 deg.
    if (!validUnitCellAngles(mrcFile.cellAngles))
    {
        diagnostic.report += "Expected all cell angles larger 0 and smaller 180 degrees, but are (" +
            std::to_string(mrcFile.cellAngles[0]) + "," +
            std::to_string(mrcFile.cellAngles[1]) + "," +
            std::to_string(mrcFile.cellAngles[2]) + ").\n";

        diagnostic.fixedHeader.cellAngles = MrcHeaderDefaults::c_cellAngles;
    }

    //! Check that each column, row section has a corresponding x, y, z dimension and vice verca
    if (!bijectiveMappingToXyz(mrcFile.crsToXyz))
    {
        diagnostic.report += "Expected column row section to x, y, z assignment to be a permutation of "
            " (0,1,2), but is (" +
            std::to_string(mrcFile.crsToXyz[0]) + "," +
            std::to_string(mrcFile.crsToXyz[1]) + "," +
            std::to_string(mrcFile.crsToXyz[2]) + ").\n";

        diagnostic.fixedHeader.crsToXyz = MrcHeaderDefaults::c_crsToXyz;
    }

    return diagnostic;
}

} // namespace gmx
