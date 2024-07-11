/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#pragma once

#include <iostream>
#include <random>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vectypes.h"

namespace gmx
{

class RandomSphericalDirectionGenerator
{
public:
    RandomSphericalDirectionGenerator(int64_t seed, bool use_old_angle_dist = false) :
        engine(seed),
        dist(0.0, 1.0),
        use_old_angle_dist(use_old_angle_dist)
    {
        if (use_old_angle_dist)
        {
            std::cout << "==== RAMD ==== Warning: Old angle distribution is used." << std::endl;
        }
    }

    DVec operator()()
    {
        // azimuth angle
        real theta = 2 * M_PI * dist(engine);

        // polar angle
        real psi;
        if (use_old_angle_dist)
        {
            psi = M_PI * dist(engine);
        }
        else
        {
            psi = std::acos(1.0 - 2 * dist(engine));
        }

        DVec direction;
        direction[0] = std::cos(theta) * std::sin(psi);
        direction[1] = std::sin(theta) * std::sin(psi);
        direction[2] = std::cos(psi);

        return direction;
    }

private:
    /// Random number generator
    std::default_random_engine engine;

    /// Random number distribution
    std::uniform_real_distribution<> dist;

    /// For backward compa
    bool use_old_angle_dist;
};

} // namespace gmx
