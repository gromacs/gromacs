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
 * Implements QMMM force provider state for checkpointing
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "qmmmforceproviderstate.h"

#include <string_view>

#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"


namespace gmx
{

namespace
{

/*! \brief This tags for parameters which will be stored/read to/from checkpoint file
 */
//! \{
const std::string c_qmTransTag_ = "qmtrans";
//! \}

} // namespace

void QMMMForceProviderState::writeState(KeyValueTreeObjectBuilder kvtBuilder, std::string_view identifier) const
{
    // Write QM Translation vector
    auto doubleArrayAdder =
            kvtBuilder.addUniformArray<double>(std::string(identifier) + "-" + c_qmTransTag_);
    for (int i = 0; i < DIM; i++)
    {
        doubleArrayAdder.addValue(static_cast<double>(qmTrans_[i]));
    }
}

void QMMMForceProviderState::readState(const KeyValueTreeObject& kvtData, std::string_view identifier)
{
    stateRead_ = true;

    // Try to read QM translation vector
    if (kvtData.keyExists(std::string(identifier) + "-" + c_qmTransTag_))
    {
        auto kvtDoubleArray = kvtData[std::string(identifier) + "-" + c_qmTransTag_].asArray().values();
        for (int i = 0; i < DIM; i++)
        {
            qmTrans_[i] = static_cast<real>(kvtDoubleArray[i].cast<double>());
        }
    }
    else
    {
        stateRead_ = false;
        return;
    }
}

} // namespace gmx
