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
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

#include "read_params.h"
#include "ramdparameters.h"

namespace gmx
{

void read_params(std::vector<t_inpfile>* inp, gmx::RAMDParameters* ramdparams, WarningHandler* wi)
{
    ramdparams->seed_ = get_eint(inp, "ramd-seed", 1234, wi);
    ramdparams->ngroups_ = get_eint(inp, "ramd-ngroups", 1, wi);
    ramdparams->groups_.resize(ramdparams->ngroups_);

    for (int i = 0; i < ramdparams->ngroups_; i++)
    {
        auto* ramdgrp = &ramdparams->groups_[i];
        auto  ramd_prefix = std::string("ramd-group") + std::to_string(i + 1);
        ramdgrp->force_ = get_ereal(inp, ramd_prefix + "-force", 600, wi);
        ramdgrp->r_min_dist_ = get_ereal(inp, ramd_prefix + "-r-min-dist", 0.0025, wi);
        ramdgrp->max_dist_ = get_ereal(inp, ramd_prefix + "-max-dist", 4.0, wi);
        ramdgrp->receptor_group_ = get_estr(inp, ramd_prefix + "-receptor", "Protein");
        ramdgrp->ligand_group_ = get_estr(inp, ramd_prefix + "-ligand", "Ligand");
        ramdgrp->receptor_group_pbcatom_ = get_eint(inp, ramd_prefix + "-receptor-pbcatom", 0, wi);
        ramdgrp->ligand_group_pbcatom_ = get_eint(inp, ramd_prefix + "-ligand-pbcatom", 0, wi);
    }

    ramdparams->eval_freq_ = get_eint(inp, "ramd-eval-freq", 50, wi);
    ramdparams->force_out_freq_ = get_eint(inp, "ramd-force-out-freq", 100, wi);
    ramdparams->old_angle_dist_ = getEnum<Boolean>(inp, "ramd-old-angle-dist", wi) != Boolean::No;
    ramdparams->connected_ligands_ = getEnum<Boolean>(inp, "ramd-connected-ligands", wi) != Boolean::No;
}

} // namespace gmx
