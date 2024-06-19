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
 * Implements helper functions for reading structure factors from datafile
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "scatteringfactors.h"

#include <cstdio>

#include <filesystem>

#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strdb.h"

#include "isotope.h"

namespace gmx
{

std::vector<AtomicStructureFactor> readAtomicStructureFactors()
{
    std::vector<AtomicStructureFactor> atomicStructureFactors;
    gmx::FilePtr                       fptr = gmx::openLibraryFile("scatteringfactors.dat");
    char                               line[1000];
    // loop over all non header lines
    while (get_a_line(fptr.get(), line, 1000))
    {
        int    proton;
        char   currentAtomType[8];
        double cohb, cma1, cma2, cma3, cma4, cmb1, cmb2, cmb3, cmb4, cmc;

        if (sscanf(line,
                   "%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   currentAtomType,
                   &proton,
                   &cohb,
                   &cma1,
                   &cma2,
                   &cma3,
                   &cma4,
                   &cmb1,
                   &cmb2,
                   &cmb3,
                   &cmb4,
                   &cmc)
            == 12)
        {
            std::array<double, 4> cma        = { cma1, cma2, cma3, cma4 };
            std::array<double, 4> cmb        = { cmb1, cmb2, cmb3, cmb4 };
            CromerMannParameters  cromerMann = { cma, cmb, cmc };
            AtomicStructureFactor asf        = { currentAtomType, proton, cohb, cromerMann };
            atomicStructureFactors.push_back(asf);
        }
    }
    return atomicStructureFactors;
}


} // namespace gmx
