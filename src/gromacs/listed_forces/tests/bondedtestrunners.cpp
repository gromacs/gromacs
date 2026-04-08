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
 * \brief CPU implementation of the bonded force test runner.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "bondedtestrunners.h"

#include "gromacs/listed_forces/bonded.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{
namespace test
{

void BondedHostTestRunner::run(const iListInput&           input,
                               const PaddedVector<RVec>&   x,
                               const t_pbc&                pbc,
                               real                        lambda,
                               const std::vector<t_iatom>& iatoms,
                               BondedKernelFlavor          flavor,
                               OutputQuantities*           output)
{
    std::vector<int> ddgatindex = { 0, 1, 2, 3 };

    output->energy = calculateSimpleBond(input.ftype.value(),
                                         iatoms.size(),
                                         iatoms.data(),
                                         &input.iparams,
                                         as_rvec_array(x.data()),
                                         output->f,
                                         output->fshift,
                                         &pbc,
                                         lambda,
                                         &output->dvdlambda,
                                         c_bondedTestCharges,
                                         nullptr,
                                         nullptr,
                                         nullptr,
                                         ddgatindex.data(),
                                         flavor);
}

bool BondedHostTestRunner::supportsFlavor(BondedKernelFlavor flavor, const iListInput& input, real lambda) const
{
    // CPU supports both flavors when not FEP or when lambda==0 (force output matches non-FEP)
    if (input.fep && lambda != 0.0)
    {
        return flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy;
    }
    return (flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy
            || flavor == BondedKernelFlavor::ForcesSimdWhenAvailable);
}

} // namespace test
} // namespace gmx
