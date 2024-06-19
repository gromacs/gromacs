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

/*! \internal \file
 * \brief Implementions of related classes for tests that want to
 * inspect energies produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "energyframe.h"

#include <cinttypes>

#include <map>
#include <string>
#include <utility>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

EnergyFrame::EnergyFrame(const t_enxframe& enxframe, const std::map<std::string, int>& indicesOfEnergyFields) :
    step_(enxframe.step), time_(enxframe.t)
{
    for (const auto& index : indicesOfEnergyFields)
    {
        if (index.second >= enxframe.nre)
        {
            GMX_THROW(InternalError(formatString(
                    "Index %d for energy %s not present in energy frame with %d energies",
                    index.second,
                    index.first.c_str(),
                    enxframe.nre)));
        }
        values_[index.first] = enxframe.ener[index.second].e;
    }
}

std::string EnergyFrame::frameName() const
{
    return formatString("Time %f Step %" PRId64, time_, step_);
}

const real& EnergyFrame::at(const std::string& name) const
{
    auto valueIterator = values_.find(name);
    if (valueIterator == values_.end())
    {
        GMX_THROW(APIError("Cannot get energy value " + name
                           + " unless previously registered when constructing EnergyFrameReader"));
    }
    return valueIterator->second;
}

EnergyFrame::MapConstIterator EnergyFrame::begin() const
{
    return values_.begin();
}

EnergyFrame::MapConstIterator EnergyFrame::end() const
{
    return values_.end();
}

EnergyFrame::MapConstIterator EnergyFrame::find(const std::string& key) const
{
    return values_.find(key);
}

} // namespace gmx
