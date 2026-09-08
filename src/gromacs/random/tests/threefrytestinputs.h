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
 * \brief Shared counter/key inputs for ThreeFry known-answer tests.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \ingroup module_random
 */
#ifndef GMX_RANDOM_TESTS_THREEFRYTESTINPUTS_H
#define GMX_RANDOM_TESTS_THREEFRYTESTINPUTS_H

#include <cstdint>

#include <array>
#include <string>
#include <tuple>

namespace gmx
{
namespace test
{

//! Input configuration for ThreeFry known-answer tests (name, ctr0, ctr1, key0, key1).
using ThreeFryKnownAnswerInput = std::tuple<std::string, uint64_t, uint64_t, uint64_t, uint64_t>;

/*! \brief Reference counter and key inputs for known-answer tests.
 *
 * The 2x64 flavors of ThreeFry64 use ctr0, ctr1, key0, key1.
 */
inline const std::array<ThreeFryKnownAnswerInput, 3> sc_threefryKnownAnswerInputs = { {
        { "AllZero", 0, 0, 0, 0 },
        { "AllOne", 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL },
        { "Pi", 0x243f6a8885a308d3ULL, 0x13198a2e03707344ULL, 0xa4093822299f31d0ULL, 0x082efa98ec4e6c89ULL },
} };

} // namespace test
} // namespace gmx

#endif // GMX_RANDOM_TESTS_THREEFRYTESTINPUTS_H
