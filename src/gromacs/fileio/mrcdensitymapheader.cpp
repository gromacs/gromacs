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
#include "gmxpre.h"

#include "mrcdensitymapheader.h"

#include "config.h"

namespace gmx
{

MrcDensityMapHeader::MrcDensityMapHeader() : spaceGroup_(SpaceGroup::P1),
                                             dataMode_(MrcDataMode::float32),
#ifdef GMX_INTEGER_BIG_ENDIAN
                                             machineStamp_(MachineStamp::bigEndian),
#else
                                             machineStamp_(MachineStamp::smallEndian),
#endif
                                             formatIdentifier_({ 'M', 'A', 'P', ' ' }
                                                               ),
                                             userDefinedFloat_(),
                                             numUsedLabels_( 1 ),
                                             cellLength_({ 1, 1, 1 }),
                                             cellAngles_({ 90, 90, 90 }),
                                             crsToXyz_({ 0, 1, 2 }),
                                             numCrs_(),
                                             crsStart_({ 0, 0, 0 }),
                                             extent_(),
                                             dataStats_(),
                                             extendedHeader_()
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
