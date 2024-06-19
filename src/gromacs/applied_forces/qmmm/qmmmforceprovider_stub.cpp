/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Stub implementation of QMMMForceProvider
 * Compiled in case if CP2K is not linked
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

#include "qmmmforceprovider.h"

enum class PbcType : int;
struct t_commrec;

namespace gmx
{
class ForceProviderInput;
class ForceProviderOutput;
class LocalAtomSet;
class MDLogger;
struct QMMMParameters;

CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

QMMMForceProvider::QMMMForceProvider(const QMMMParameters& parameters,
                                     const LocalAtomSet&   localQMAtomSet,
                                     const LocalAtomSet&   localMMAtomSet,
                                     PbcType               pbcType,
                                     const MDLogger&       logger) :
    parameters_(parameters),
    qmAtoms_(localQMAtomSet),
    mmAtoms_(localMMAtomSet),
    pbcType_(pbcType),
    logger_(logger),
    box_{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }
{
    GMX_THROW(
            InternalError("CP2K has not been linked into GROMACS, QMMM simulation is not "
                          "possible.\nPlease, reconfigure GROMACS with -DGMX_CP2K=ON\n"));
}

QMMMForceProvider::~QMMMForceProvider() {}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool QMMMForceProvider::isQMAtom(Index /*globalAtomIndex*/)
{
    GMX_THROW(
            InternalError("CP2K has not been linked into GROMACS, QMMM simulation is not "
                          "possible.\nPlease, reconfigure GROMACS with -DGMX_CP2K=ON\n"));
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void QMMMForceProvider::appendLog(const std::string& /*msg*/)
{
    GMX_THROW(
            InternalError("CP2K has not been linked into GROMACS, QMMM simulation is not "
                          "possible.\nPlease, reconfigure GROMACS with -DGMX_CP2K=ON\n"));
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void QMMMForceProvider::initCP2KForceEnvironment(const t_commrec& /*cr*/)
{
    GMX_THROW(
            InternalError("CP2K has not been linked into GROMACS, QMMM simulation is not "
                          "possible.\nPlease, reconfigure GROMACS with -DGMX_CP2K=ON\n"));
}

void QMMMForceProvider::calculateForces(const ForceProviderInput& /*fInput*/, ForceProviderOutput* /*fOutput*/)
{
    GMX_THROW(
            InternalError("CP2K has not been linked into GROMACS, QMMM simulation is not "
                          "possible.\nPlease, reconfigure GROMACS with -DGMX_CP2K=ON\n"));
};

CLANG_DIAGNOSTIC_RESET

} // namespace gmx
