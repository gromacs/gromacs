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
 * This file defines functions for setting up a nonbonded system
 *
 * \author Berk Hess <hess@kth.se>
 * \author Joe Jordan <ejjordan@kth.se>
 *
 */

#include "gmxpre.h"

#include "system.h"

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"

#include "coords.h"

namespace
{

//! The number of atoms in a molecule
constexpr int  numAtomsInMolecule = 1;
//! The atom type of the atom
constexpr int  atomType         = 0;
//! The charge of the atom
constexpr real atomCharge       = 0.0;

}   // namespace

//! Generates coordinates and a box for the base system scaled by \p multiplicationFactor
//
// The parameter \p multiplicationFactor should be a power of 2.
// A fatal error is generated when this is not the case.
static void
generateCoordinates(int                     multiplicationFactor,
                    std::vector<gmx::RVec> *coordinates,
                    matrix                  box)
{
    if (multiplicationFactor < 1 ||
        (multiplicationFactor & (multiplicationFactor - 1)) != 0)
    {
        gmx_fatal(FARGS, "The size factor has to be a power of 2");
    }

    if (multiplicationFactor == 1)
    {
        *coordinates = coordinates12;
        copy_mat(box12, box);

        return;
    }

    ivec factors = { 1, 1, 1 };

    int  dim = 0;
    while (multiplicationFactor > 1)
    {
        factors[dim]         *= 2;
        multiplicationFactor /= 2;
        dim++;
        if (dim == DIM)
        {
            dim = 0;
        }
    }
    printf("Stacking a box of %zu atoms %d x %d x %d times\n",
           coordinates12.size(), factors[XX], factors[YY], factors[ZZ]);

    coordinates->resize(factors[XX]*factors[YY]*factors[ZZ]*coordinates12.size());

    int       i = 0;
    gmx::RVec shift;
    for (int x = 0; x < factors[XX]; x++)
    {
        shift[XX] = x*box12[XX][XX];
        for (int y = 0; y < factors[YY]; y++)
        {
            shift[YY] = y*box12[YY][YY];
            for (int z = 0; z < factors[ZZ]; z++)
            {
                shift[ZZ] = z*box12[ZZ][ZZ];

                for (const gmx::RVec &coordOrig : coordinates12)
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
            box[d1][d2] = factors[d1]*box12[d1][d2];
        }
    }
}

NBKernelSystem::NBKernelSystem(const int multiplicationFactor)
{
    nonbondedParameters.resize(numAtomTypes*numAtomTypes*2, 0);
    nonbondedParameters[0] = c6Param;
    nonbondedParameters[1] = c12Param;

    generateCoordinates(multiplicationFactor, &coordinates, box);
    put_atoms_in_box(epbcXYZ, box, coordinates);

    int numAtoms = coordinates.size();
    GMX_RELEASE_ASSERT(numAtoms % numAtomsInMolecule == 0, "Coordinates should match whole molecules");

    atomTypes.resize(numAtoms);
    charges.resize(numAtoms);
    atomInfoAllVdw.resize(numAtoms);
    snew(excls.index, numAtoms + 1);
    snew(excls.a, numAtoms*numAtomsInMolecule);
    excls.index[0] = 0;

    for (int atom = 0; atom < numAtoms; atom++)
    {
        atomTypes[atom] = atomType;
        charges[atom]   = atomCharge;
        SET_CGINFO_HAS_VDW(atomInfoAllVdw[atom]);
        SET_CGINFO_HAS_Q(atomInfoAllVdw[atom]);

        const int firstAtomInMolecule = atom - (atom % numAtomsInMolecule);
        for (int atomJ = 0; atomJ < numAtomsInMolecule; atomJ++)
        {
            excls.a[atom*numAtomsInMolecule + atomJ] = firstAtomInMolecule + atomJ;
        }
        excls.index[atom + 1] = (atom + 1)*numAtomsInMolecule;
    }
}

namespace nblib {
SimState::SimState(const std::vector<gmx::RVec> &coord, Box box, Topology &topo,
                   const std::vector<gmx::RVec> &vel) : box_(box), topo_(topo)
{
    coord_ = coord;
    vel_ = vel;
}

SimState::SimState(const SimState &simState)
    : box_(simState.box_), topo_(simState.topo_)
{
    coord_ = simState.coord_;
    vel_   = simState.vel_;
}

SimState &SimState::operator=(const SimState &simState)
{
    coord_ = simState.coord_;
    vel_   = simState.vel_;
    box_   = simState.box_;
    topo_  = simState.topo_;
}

SimState::SimState(SimState &&simState) noexcept
    : box_(simState.box_), topo_(std::move(simState.topo_))
{
    coord_ = std::move(simState.coord_);
    vel_   = std::move(simState.vel_);
}

SimState& SimState::operator=(nblib::SimState &&simState) noexcept
{
    coord_ = std::move(simState.coord_);
    vel_   = std::move(simState.vel_);
    box_   = simState.box_;
    topo_  = std::move(simState.topo_);
}

} // namespace nblib
