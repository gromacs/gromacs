/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Implements class for SANS Debye Scattering
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "scattering-debye-sans.h"

#include <utility>
#include <vector>

#include "gromacs/trajectoryanalysis/modules/isotope.h"
#include "gromacs/trajectoryanalysis/modules/scatteringfactors.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/strdb.h"

namespace gmx
{

SansDebye::SansDebye(std::vector<Isotope> isotopes) : isotopes_(std::move(isotopes))
{
    std::vector<AtomicStructureFactor> atomicScatterFactors = readAtomicStructureFactors();
    std::vector<double>                scatterFactors(static_cast<int>(Isotope::Count));
    for (const auto& scatter : atomicScatterFactors)
    {
        Isotope isotope                              = getIsotopeFromString(scatter.isotope);
        scatterFactors.at(static_cast<int>(isotope)) = scatter.neutronIsotropicScatteringLength;
    }
    scatterFactors_ = scatterFactors;
}


double SansDebye::getScatteringLength(int i, double gmx_unused q)
{
    return scatterFactors_[static_cast<int>(isotopes_[i])];
}


} // namespace gmx
