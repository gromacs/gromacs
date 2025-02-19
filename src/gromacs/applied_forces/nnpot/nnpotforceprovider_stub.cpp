/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Stub implementation of the NNPot Force Provider class.
 * Compiled in case Libtorch/NN backend is not linked.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

#include "nnpotforceprovider.h"
#include "nnpotoptions.h"

namespace gmx
{

CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

NNPotForceProvider::NNPotForceProvider(const NNPotParameters& nnpotParameters, const MDLogger* logger) :
    params_(nnpotParameters),
    positions_(params_.numAtoms_, RVec({ 0.0, 0.0, 0.0 })),
    atomNumbers_(params_.numAtoms_, -1),
    idxLookup_(params_.numAtoms_, -1),
    box_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } },
    logger_(logger),
    cr_(params_.cr_)
{
    GMX_THROW(InternalError(
            "Libtorch/NN backend is not linked into GROMACS, NNPot simulation is not possible."
            " Please, reconfigure GROMACS with -DGMX_NNPOT=TORCH\n"));
}

NNPotForceProvider::~NNPotForceProvider() {}

void NNPotForceProvider::calculateForces(const ForceProviderInput& /*fInput*/, ForceProviderOutput* /*fOutput*/)
{
    GMX_THROW(InternalError(
            "Libtorch/NN backend is not linked into GROMACS, NNPot simulation is not possible."
            " Please, reconfigure GROMACS with -DGMX_NNPOT=TORCH\n"));
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void NNPotForceProvider::gatherAtomNumbersIndices()
{
    GMX_THROW(InternalError(
            "Libtorch/NN backend is not linked into GROMACS, NNPot simulation is not possible."
            " Please, reconfigure GROMACS with -DGMX_NNPOT=TORCH\n"));
}

CLANG_DIAGNOSTIC_RESET

} // namespace gmx
