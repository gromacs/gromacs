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
 * Implements nblib utility functions
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "util.h"

#include "gromacs/math/units.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/fatalerror.h"

namespace nblib {

static std::vector<gmx::RVec> low_mspeed(real tempi,  std::vector<real> const& masses, gmx::ThreeFry2x64<>* rng)
{
    int                                    nrdf;
    real                                   boltz;
    real                                   ekin, temp;
    gmx::TabulatedNormalDistribution<real> normalDist;

    std::vector<gmx::RVec> velocities(masses.size());

    boltz = BOLTZ * tempi;
    ekin  = 0.0;
    nrdf  = 0;
    for(size_t i = 0; i < masses.size(); i++) {
        real mass  = masses[i];
        if (mass > 0)
        {
            rng->restart(i, 0);
            real sd = std::sqrt(boltz / mass);
            for (int m = 0; (m < DIM); m++)
            {
                velocities[i][m] = sd * normalDist(*rng);
                ekin += 0.5 * mass * velocities[i][m] * velocities[i][m];
            }
            nrdf += DIM;
        }

    }
    temp = (2.0 * ekin) / (nrdf * BOLTZ);
    if (temp > 0)
    {
        real scal = std::sqrt(tempi / temp);
        for(auto &vel : velocities)
        {
            for (int m = 0; (m < DIM); m++)
            {
                vel[m] *= scal;
            }
        }
    }
    fprintf(stderr, "Velocities were taken from a Maxwell distribution at %g K\n", tempi);
    if (debug)
    {
        fprintf(debug,
                "Velocities were taken from a Maxwell distribution\n"
                "Initial generated temperature: %12.5e (scaled to: %12.5e)\n",
                temp, tempi);
    }

    return velocities;
}

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
std::vector<gmx::RVec> generateVelocity(real tempi, unsigned int seed, std::vector<real> const& masses)
{

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Using random seed %u for generating velocities\n", seed);
    }
    gmx::ThreeFry2x64<> rng(seed, gmx::RandomDomain::MaxwellVelocities);

    return low_mspeed(tempi, masses, &rng);
}

bool checkNumericValues(const std::vector<gmx::RVec> &values)
{
    for (auto val : values)
    {
        for (int m = 0; (m < DIM); m++)
        {
            if(std::isnan(val[m]) or std::isinf(val[m]))
            {
                return false;
            }
        }
    }
    return true;
}

} // namespace nblib
