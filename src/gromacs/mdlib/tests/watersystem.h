/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include "gromacs/math/vectypes.h"

namespace gmx
{

namespace test
{

//! Database of 51 water atom input positions (taken from spc216.gro) for use as test inputs.
static const std::array<gmx::RVec, 51> c_waterPositions(
        { { { .130, -.041, -.291 },  { .120, -.056, -.192 },  { .044, -.005, -.327 },
            { -.854, -.406, .477 },  { -.900, -.334, .425 },  { -.858, -.386, .575 },
            { .351, -.061, .853 },   { .401, -.147, .859 },   { .416, .016, .850 },
            { -.067, -.796, .873 },  { -.129, -.811, .797 },  { -.119, -.785, .958 },
            { -.635, -.312, -.356 }, { -.629, -.389, -.292 }, { -.687, -.338, -.436 },
            { .321, -.919, .242 },   { .403, -.880, .200 },   { .294, -1.001, .193 },
            { -.404, .735, .728 },   { -.409, .670, .803 },   { -.324, .794, .741 },
            { .461, -.596, -.135 },  { .411, -.595, -.221 },  { .398, -.614, -.059 },
            { -.751, -.086, .237 },  { -.811, -.148, .287 },  { -.720, -.130, .152 },
            { .202, .285, -.364 },   { .122, .345, -.377 },   { .192, .236, -.278 },
            { -.230, -.485, .081 },  { -.262, -.391, .071 },  { -.306, -.548, .069 },
            { .464, -.119, .323 },   { .497, -.080, .409 },   { .540, -.126, .258 },
            { -.462, .107, .426 },   { -.486, .070, .336 },   { -.363, .123, .430 },
            { .249, -.077, -.621 },  { .306, -.142, -.571 },  { .233, -.110, -.714 },
            { -.922, -.164, .904 },  { -.842, -.221, .925 },  { -.971, -.204, .827 },
            { .382, .700, .480 },    { .427, .610, .477 },    { .288, .689, .513 },
            { .781, .264, -.113 },   { .848, .203, -.070 },   { .708, .283, -.048 } } });

} // namespace test

} // namespace gmx
