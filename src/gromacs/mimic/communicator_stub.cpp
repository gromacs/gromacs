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
#include "gmxpre.h"

#include <cstdint>

#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

#include "communicator.h"

struct gmx_mtop_t;

namespace gmx
{

CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

void MimicCommunicator::init()
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

void MimicCommunicator::sendInitData(gmx_mtop_t* /*mtop*/, ArrayRef<const RVec> /*coords*/)
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

int64_t MimicCommunicator::getStepNumber()
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

void MimicCommunicator::getCoords(ArrayRef<RVec> /*x*/, const int /*natoms*/)
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

void MimicCommunicator::sendEnergies(real /*energy*/)
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

void MimicCommunicator::sendForces(ArrayRef<RVec> /*forces*/, int /*natoms*/)
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

void MimicCommunicator::finalize()
{
    GMX_THROW(InternalError(
            "GROMACS is compiled without MiMiC support! Please, reconfigure with -DGMX_MIMIC=ON"));
}

CLANG_DIAGNOSTIC_RESET

} // namespace gmx
