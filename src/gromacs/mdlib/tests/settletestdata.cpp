/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief Defines the class that accumulates SETTLE test data.
 *
 * \author Mark Abraham  <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settletestdata.h"

#include <cstdlib>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/mdlib/tests/watersystem.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

SettleTestData::SettleTestData(int numSettles) :
    x_(c_waterPositions.size()), xPrime_(c_waterPositions.size()), v_(c_waterPositions.size())
{
    // Initialize coordinates and velocities from the constant set of coordinates
    std::copy(c_waterPositions.begin(), c_waterPositions.end(), x_.begin());
    std::copy(c_waterPositions.begin(), c_waterPositions.end(), xPrime_.begin());

    // Perturb the atom positions, to appear like an
    // "update," and where there is definitely constraining
    // work to do.
    const real deltas[] = { 0.01, -0.01, +0.02, -0.02 };
    int        i        = 0;
    for (auto& xPrime : xPrime_)
    {
        xPrime[XX] += deltas[i % 4];
        ++i;
        xPrime[YY] += deltas[i % 4];
        ++i;
        xPrime[ZZ] += deltas[i % 4];
        ++i;
    }
    std::fill(v_.begin(), v_.end(), RVec{ 0.0, 0.0, 0.0 });

    // Set up the topology.
    const int settleType = 0;
    mtop_.moltype.resize(1);
    mtop_.molblock.resize(1);
    mtop_.molblock[0].type   = 0;
    std::vector<int>& iatoms = mtop_.moltype[0].ilist[F_SETTLE].iatoms;
    for (int i = 0; i < numSettles; ++i)
    {
        iatoms.push_back(settleType);
        iatoms.push_back(i * atomsPerSettle_ + 0);
        iatoms.push_back(i * atomsPerSettle_ + 1);
        iatoms.push_back(i * atomsPerSettle_ + 2);
    }

    // Set up the SETTLE parameters.
    t_iparams iparams;
    iparams.settle.doh = dOH_;
    iparams.settle.dhh = dHH_;
    mtop_.ffparams.iparams.push_back(iparams);

    // Set up the masses.
    mtop_.moltype[0].atoms.atom =
            static_cast<t_atom*>(calloc(numSettles * atomsPerSettle_, sizeof(t_atom)));
    numAtoms_ = numSettles * atomsPerSettle_;
    masses_.resize(numAtoms_);
    inverseMasses_.resize(numAtoms_);
    for (int i = 0; i < numSettles; ++i)
    {
        masses_[i * atomsPerSettle_ + 0] = oxygenMass_;
        masses_[i * atomsPerSettle_ + 1] = hydrogenMass_;
        masses_[i * atomsPerSettle_ + 2] = hydrogenMass_;

        inverseMasses_[i * atomsPerSettle_ + 0] = 1.0 / oxygenMass_;
        inverseMasses_[i * atomsPerSettle_ + 1] = 1.0 / hydrogenMass_;
        inverseMasses_[i * atomsPerSettle_ + 2] = 1.0 / hydrogenMass_;

        mtop_.moltype[0].atoms.atom[i * atomsPerSettle_ + 0].m = oxygenMass_;
        mtop_.moltype[0].atoms.atom[i * atomsPerSettle_ + 1].m = hydrogenMass_;
        mtop_.moltype[0].atoms.atom[i * atomsPerSettle_ + 2].m = hydrogenMass_;
    }

    idef_               = std::make_unique<InteractionDefinitions>(mtop_.ffparams);
    idef_->il[F_SETTLE] = mtop_.moltype[0].ilist[F_SETTLE];
}

SettleTestData::~SettleTestData() {}

} // namespace test
} // namespace gmx
