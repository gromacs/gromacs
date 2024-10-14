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
 * Defines options for Plumed. This class handles parameters set during
 * pre-processing time.
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#include "plumedOptions.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

void PlumedOptionProvider::setTopology(const gmx_mtop_t& mtop)
{
    opts_.natoms_ = mtop.natoms;
}

void PlumedOptionProvider::setEnsembleTemperature(const EnsembleTemperature& temp)
{
    if (temp.constantEnsembleTemperature_)
    {
        opts_.ensembleTemperature_ = temp.constantEnsembleTemperature_;
    }
}

void PlumedOptionProvider::setPlumedFile(const std::optional<std::string>& fname)
{
    if (fname.has_value())
    {
        opts_.active_     = true;
        opts_.plumedFile_ = fname.value();
    }
}
const PlumedOptions& PlumedOptionProvider::options() const
{
    return opts_;
}

void PlumedOptionProvider::setSimulationTimeStep(double timeStep)
{
    opts_.simulationTimeStep_ = timeStep;
}

void PlumedOptionProvider::setStartingBehavior(const StartingBehavior& behavior)
{
    opts_.startingBehavior_ = behavior;
}

void PlumedOptionProvider::setComm(const t_commrec& cr)
{
    opts_.cr_ = &cr;
}

bool PlumedOptionProvider::active() const
{
    return opts_.active_;
}


} // namespace gmx
