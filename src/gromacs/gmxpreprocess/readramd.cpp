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
#include "gmxpre.h"

#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/ramd_params.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

void read_ramdparams(std::vector<t_inpfile>* inp, gmx::RAMDParams* ramdparams, WarningHandler* wi)
{
    ramdparams->seed   = get_eint(inp, "ramd-seed", 1234, wi);
    ramdparams->ngroup = get_eint(inp, "ramd-ngroups", 1, wi);
    ramdparams->group.resize(ramdparams->ngroup);

    for (int i = 0; i < ramdparams->ngroup; i++)
    {
        auto* ramdgrp      = &ramdparams->group[i];
        auto  ramd_prefix  = std::string("ramd-group") + std::to_string(i + 1);
        auto  pull1_prefix = std::string("pull-group") + std::to_string(i * 2 + 1);
        auto  pull2_prefix = std::string("pull-group") + std::to_string(i * 2 + 2);

        inp->emplace_back(0,
                          1,
                          false,
                          false,
                          false,
                          pull1_prefix + "-name",
                          get_estr(inp, ramd_prefix + "-receptor", "Protein"));
        inp->emplace_back(0,
                          1,
                          false,
                          false,
                          false,
                          pull1_prefix + "-pbcatom",
                          get_estr(inp, ramd_prefix + "-receptor-pbcatom", "0"));

        inp->emplace_back(0,
                          1,
                          false,
                          false,
                          false,
                          pull2_prefix + "-name",
                          get_estr(inp, ramd_prefix + "-ligand", "INH"));
        inp->emplace_back(0,
                          1,
                          false,
                          false,
                          false,
                          pull2_prefix + "-pbcatom",
                          get_estr(inp, ramd_prefix + "-ligand-pbcatom", "0"));

        ramdgrp->force      = get_ereal(inp, ramd_prefix + "-force", 600, wi);
        ramdgrp->r_min_dist = get_ereal(inp, ramd_prefix + "-r-min-dist", 0.0025, wi);
        ramdgrp->max_dist   = get_ereal(inp, ramd_prefix + "-max-dist", 4.0, wi);
    }

    inp->emplace_back(0,
                      1,
                      false,
                      false,
                      false,
                      "pull-pbc-ref-prev-step-com",
                      get_estr(inp, "ramd-pbc-ref-prev-step-com", "yes"));

    ramdparams->eval_freq      = get_eint(inp, "ramd-eval-freq", 50, wi);
    ramdparams->force_out_freq = get_eint(inp, "ramd-force-out-freq", 100, wi);
    ramdparams->old_angle_dist = getEnum<Boolean>(inp, "ramd-old-angle-dist", wi) != Boolean::No;

    inp->emplace_back(0, 1, false, false, false, "pull-ngroups", std::to_string(ramdparams->ngroup * 2));
    inp->emplace_back(
            0, 1, false, false, false, "pull-nstxout", std::to_string(ramdparams->force_out_freq));
    inp->emplace_back(
            0, 1, false, false, false, "pull-nstfout", std::to_string(ramdparams->force_out_freq));
    inp->emplace_back(0, 1, false, false, false, "pull-ncoords", std::to_string(ramdparams->ngroup * 3));

    std::vector<std::string> v{ "1 0 0", "0 1 0", "0 0 1" };
    for (int i = 0; i < ramdparams->ngroup; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            auto prefix = std::string("pull-coord") + std::to_string(i * 3 + d + 1);
            inp->emplace_back(0,
                              1,
                              false,
                              false,
                              false,
                              prefix + "-groups",
                              std::to_string(i * 2 + 1) + " " + std::to_string(i * 2 + 2));
            inp->emplace_back(0, 1, false, false, false, prefix + "-type", "external-potential");
            inp->emplace_back(0, 1, false, false, false, prefix + "-potential-provider", "RAMD");
            inp->emplace_back(0, 1, false, false, false, prefix + "-geometry", "direction");
            inp->emplace_back(0, 1, false, false, false, prefix + "-vec", v[d]);
        }
    }

    ramdparams->connected_ligands = getEnum<Boolean>(inp, "ramd-connected-ligands", wi) != Boolean::No;
}
