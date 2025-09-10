/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Helper functions for topology generation.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

#include "topologyhelpers.h"

#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

namespace gmx
{
namespace test
{

void addNWaterMolecules(gmx_mtop_t* mtop, int numWaters)
{
    // Target distance between oxygen and hydrogens
    constexpr real dOH = 0.09572;
    // Target distance between hydrogens
    constexpr real dHH = 0.15139;

    gmx_moltype_t moltype;
    moltype.atoms.nr             = NRAL(F_SETTLE);
    std::vector<int>& iatoms     = moltype.ilist[F_SETTLE].iatoms;
    const int         settleType = 0;
    iatoms.push_back(settleType);
    iatoms.push_back(0);
    iatoms.push_back(1);
    iatoms.push_back(2);
    int moleculeTypeIndex = mtop->moltype.size();
    mtop->moltype.push_back(moltype);
    init_t_atoms(&mtop->moltype[0].atoms, NRAL(F_SETTLE), false);
    for (int i = 0; i < NRAL(F_SETTLE); ++i)
    {
        mtop->moltype[0].atoms.atom[i].m = (i % 3 == 0) ? 16 : 1;
    }

    mtop->molblock.emplace_back();
    mtop->molblock.back().type = moleculeTypeIndex;
    mtop->molblock.back().nmol = numWaters;
    mtop->natoms               = moltype.atoms.nr * mtop->molblock.back().nmol;

    // Set up the SETTLE parameters.
    t_iparams iparams;
    iparams.settle.doh = dOH;
    iparams.settle.dhh = dHH;
    mtop->ffparams.iparams.push_back(iparams);
}

} // namespace test
} // namespace gmx
