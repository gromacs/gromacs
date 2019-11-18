/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 * \brief
 * This file defines functions for setting up a benchmark system
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "bench_system.h"

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"

#include "bench_coords.h"

namespace gmx
{

namespace
{

// A 3-site water model
//! The number of atoms in a molecule
constexpr int numAtomsInMolecule = 3;
//! The atom type of the oxygen atom
constexpr int typeOxygen = 0;
//! The atom type of the hydrogen atom
constexpr int typeHydrogen = 1;
//! The charge of the oxygen atom
constexpr real chargeOxygen = -0.8476;
//! The charge of the hydrogen atom
constexpr real chargeHydrogen = 0.4238;
//! The LJ C6 parameter of the Oxygen atom
constexpr real c6Oxygen = 0.0026173456;
//! The LJ C12 parameter of the Oxygen atom
constexpr real c12Oxygen = 2.634129e-06;
// Note that the hydrogen has LJ parameters all zero

} // namespace

//! Generates coordinates and a box for the base system scaled by \p multiplicationFactor
//
// The parameter \p multiplicationFactor should be a power of 2.
// A fatal error is generated when this is not the case.
static void generateCoordinates(int multiplicationFactor, std::vector<gmx::RVec>* coordinates, matrix box)
{
    if (multiplicationFactor < 1 || (multiplicationFactor & (multiplicationFactor - 1)) != 0)
    {
        gmx_fatal(FARGS, "The size factor has to be a power of 2");
    }

    if (multiplicationFactor == 1)
    {
        *coordinates = coordinates1000;
        copy_mat(box1000, box);

        return;
    }

    ivec factors = { 1, 1, 1 };

    int dim = 0;
    while (multiplicationFactor > 1)
    {
        factors[dim] *= 2;
        multiplicationFactor /= 2;
        dim++;
        if (dim == DIM)
        {
            dim = 0;
        }
    }
    printf("Stacking a box of %zu atoms %d x %d x %d times\n", coordinates1000.size(), factors[XX],
           factors[YY], factors[ZZ]);

    coordinates->resize(factors[XX] * factors[YY] * factors[ZZ] * coordinates1000.size());

    int       i = 0;
    gmx::RVec shift;
    for (int x = 0; x < factors[XX]; x++)
    {
        shift[XX] = x * box1000[XX][XX];
        for (int y = 0; y < factors[YY]; y++)
        {
            shift[YY] = y * box1000[YY][YY];
            for (int z = 0; z < factors[ZZ]; z++)
            {
                shift[ZZ] = z * box1000[ZZ][ZZ];

                for (const gmx::RVec& coordOrig : coordinates1000)
                {
                    (*coordinates)[i] = coordOrig + shift;
                    i++;
                }
            }
        }
    }

    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            box[d1][d2] = factors[d1] * box1000[d1][d2];
        }
    }
}

BenchmarkSystem::BenchmarkSystem(const int multiplicationFactor)
{
    numAtomTypes = 2;
    nonbondedParameters.resize(numAtomTypes * numAtomTypes * 2, 0);
    nonbondedParameters[0] = c6Oxygen;
    nonbondedParameters[1] = c12Oxygen;

    generateCoordinates(multiplicationFactor, &coordinates, box);
    put_atoms_in_box(epbcXYZ, box, coordinates);

    int numAtoms = coordinates.size();
    GMX_RELEASE_ASSERT(numAtoms % numAtomsInMolecule == 0,
                       "Coordinates should match whole molecules");

    atomTypes.resize(numAtoms);
    charges.resize(numAtoms);
    atomInfoAllVdw.resize(numAtoms);
    atomInfoOxygenVdw.resize(numAtoms);
    snew(excls.index, numAtoms + 1);
    snew(excls.a, numAtoms * numAtomsInMolecule);
    excls.index[0] = 0;

    for (int a = 0; a < numAtoms; a++)
    {
        if (a % numAtomsInMolecule == 0)
        {
            // Oxgygen
            atomTypes[a] = typeOxygen;
            charges[a]   = chargeOxygen;
            SET_CGINFO_HAS_VDW(atomInfoAllVdw[a]);
            SET_CGINFO_HAS_VDW(atomInfoOxygenVdw[a]);
        }
        else
        {
            // Hydrogen
            atomTypes[a] = typeHydrogen;
            charges[a]   = chargeHydrogen;
            SET_CGINFO_HAS_VDW(atomInfoAllVdw[a]);
        }
        SET_CGINFO_HAS_Q(atomInfoAllVdw[a]);
        SET_CGINFO_HAS_Q(atomInfoOxygenVdw[a]);

        const int firstAtomInMolecule = a - (a % numAtomsInMolecule);
        for (int aj = 0; aj < numAtomsInMolecule; aj++)
        {
            excls.a[a * numAtomsInMolecule + aj] = firstAtomInMolecule + aj;
        }
        excls.index[a + 1] = (a + 1) * numAtomsInMolecule;
    }

    forceRec.ntype = numAtomTypes;
    forceRec.nbfp  = nonbondedParameters;
    snew(forceRec.shift_vec, SHIFTS);
    calc_shifts(box, forceRec.shift_vec);
}

} // namespace gmx
