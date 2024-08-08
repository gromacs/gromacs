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

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/// Parameters for RAMD group
struct RAMDGroup
{
    real force;      ///< Force to be applied in kcal/mol/Angstrom
    real max_dist;   ///< Specifies the distance in Angstrom between the COMs of the ligand
                     ///  and the receptor when the simulation is stopped
    real r_min_dist; ///< Specifies the minimum distance in Angstrom
                     ///  to be traveled by the ligand in one RAMD step
};

/// Parameters for Random Acceleration Molecular Dynamics (RAMD)
struct RAMDParams
{
    int64_t    seed;              ///< Initialization number for pseudo random number generator
    int        ngroup;            ///< Number of RAMD groups
    std::vector<RAMDGroup> group; ///< List of RAMD receptor-ligand pairs
    int        eval_freq;         ///< Number of MD steps in one RAMD step
    int force_out_freq;           ///< Every 'force_out_freq' steps detailed output of forces will be written
    gmx_bool old_angle_dist;      ///< Use old angle distribution
    gmx_bool connected_ligands;   ///< Behavior of re-entering ligands into the dissociation radius
};

} // namespace gmx
