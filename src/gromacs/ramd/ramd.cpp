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

#include "ramd.h"

#include <cassert>
#include <numeric>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

RAMD::RAMD(const RAMDParams&           params,
           pull_t*                     pull,
           const gmx::StartingBehavior startingBehavior,
           const t_commrec*            cr,
           int                         nfile,
           const t_filenm              fnm[],
           const gmx_output_env_t*     oenv,
           FILE*                       log)
  : params(params),
    pull(pull),
    random_spherical_direction_generator(params.seed, params.old_angle_dist),
    direction(params.ngroup),
    com_rec_prev(params.ngroup),
    com_lig_prev(params.ngroup),
    out(nullptr),
    cr(cr),
    ligand_exited(params.ngroup, 0),
    write_trajectory(false),
    log(log)
{
    for (int g = 0; g < params.ngroup; ++g)
    {
        register_external_pull_potential(pull, g * 3    , "RAMD");
        register_external_pull_potential(pull, g * 3 + 1, "RAMD");
        register_external_pull_potential(pull, g * 3 + 2, "RAMD");
        direction[g] = random_spherical_direction_generator();
    }

    if (MAIN(cr) and opt2bSet("-ramd", nfile, fnm))
    {
        auto filename = std::string(opt2fn("-ramd", nfile, fnm));

        if (startingBehavior == gmx::StartingBehavior::RestartWithAppending)
        {
            out = gmx_fio_fopen(filename.c_str(), "a+");
        }
        else
        {
            out = gmx_fio_fopen(filename.c_str(), "w+");
            xvgr_header(out, "RAMD Ligand-receptor COM distance", "Time (ps)", "Distance (nm)", exvggtXNY, oenv);

            std::vector<std::string> setnames;
            for (int g = 1; g <= params.ngroup; ++g)
            {
                setnames.push_back(std::to_string(g));
            }
            xvgrLegend(out, setnames, oenv);
        }
        fflush(out);
    }
}

void RAMD::calculateForces(const ForceProviderInput& forceProviderInput,
                           [[maybe_unused]] ForceProviderOutput* forceProviderOutput)
{
    t_pbc pbc;
    set_pbc(&pbc, pull->pbcType, forceProviderInput.box_);

    auto step = forceProviderInput.step_;
    if (MAIN(cr) and out and (step % params.eval_freq == 0))
    {
        fprintf(out, "%.4f", forceProviderInput.t_);
    }

    if (step == 0)
    {
        // Store COM positions for first evaluation
        for (int g = 0; g < params.ngroup; ++g)
        {
            com_rec_prev[g] = pull->group[g * 2 + 1].x;
            com_lig_prev[g] = pull->group[g * 2 + 2].x;

            if (MAIN(cr) and out)
            {
                DVec curr_dist_vect;
                pbc_dx_d(&pbc, com_lig_prev[g], com_rec_prev[g], curr_dist_vect);
                auto dist = std::sqrt(curr_dist_vect.norm2());
                fprintf(out, "\t%g", dist);
            }
        }
    }
    else if (step % params.eval_freq == 0)
    {
        if (MAIN(cr) and debug)
        {
            fprintf(debug, "==== RAMD ==== evaluation %ld\n", step);
        }
        for (int g = 0; g < params.ngroup; ++g)
        {
            DVec com_rec_curr = pull->group[g * 2 + 1].x;
            DVec com_lig_curr = pull->group[g * 2 + 2].x;
            DVec curr_dist_vect;
            pbc_dx_d(&pbc, com_lig_curr, com_rec_curr, curr_dist_vect);
            auto curr_dist = std::sqrt(curr_dist_vect.norm2());

            if (MAIN(cr) and debug)
            {
                fprintf(debug, "==== RAMD ==== group %d\n", g);
                fprintf(debug, "==== RAMD ==== COM ligand position at [%g, %g, %g]\n",
                        com_lig_curr[0], com_lig_curr[1], com_lig_curr[2]);
                fprintf(debug, "==== RAMD ==== COM receptor position at [%g, %g, %g]\n",
                        com_rec_curr[0], com_rec_curr[1], com_rec_curr[2]);
                fprintf(debug,
                        "==== RAMD ==== Distance between COM of receptor and COM of ligand is %g\n",
                        curr_dist);
            }

            if (MAIN(cr) and out)
            {
                fprintf(out, "\t%g", curr_dist);
            }

            if (curr_dist >= params.group[g].max_dist)
            {
                ligand_exited[g] = 1;
                if (!params.connected_ligands) {
                    direction[g] = DVec(0.0, 0.0, 0.0);
                }
                if (MAIN(cr))
                {
                    fprintf(this->log, "==== RAMD ==== RAMD group %d has exited the binding site in step %ld\n",
                            g, step);
                }
            } else if (ligand_exited[g] == 1) {
                ligand_exited[g] = 0;
            }

            // difference of the COM ligand-receptor distance between current and the last evaluation step
            DVec walk_dist_vect;
            pbc_dx_d(&pbc, com_lig_curr - com_rec_curr, com_lig_prev[g] - com_rec_prev[g], walk_dist_vect);
            auto walk_dist = std::sqrt(walk_dist_vect.norm2());

            if (MAIN(cr) and debug)
            {
                fprintf(debug, "==== RAMD ==== Previous COM ligand position at [%g, %g, %g]\n",
                        com_lig_prev[g][0], com_lig_prev[g][1], com_lig_prev[g][2]);
                fprintf(debug, "==== RAMD ==== Previous COM receptor position at [%g, %g, %g]\n",
                        com_rec_prev[g][0], com_rec_prev[g][1], com_rec_prev[g][2]);
                fprintf(debug,
                        "==== RAMD ==== Change in receptor-ligand distance since last RAMD evaluation "
                        "is %g\n",
                        walk_dist);
            }

            if (walk_dist < params.group[0].r_min_dist)
            {
                direction[g] = random_spherical_direction_generator();
                if (MAIN(cr) and debug)
                {
                    fprintf(debug, "==== RAMD ==== New random direction is [%g, %g, %g]\n",
                            direction[g][0], direction[g][1], direction[g][2]);
                }
            }

            com_lig_prev[g] = com_lig_curr;
            com_rec_prev[g] = com_rec_curr;
        }
    }

    if (step % params.eval_freq == 0)
    {
        if (MAIN(cr) and out)
        {
            fprintf(out, "\n");
            fflush(out);
        }

        // Exit if all ligand-receptor COM distances are larger than max_dist
        if (std::accumulate(ligand_exited.begin(), ligand_exited.end(), 0) == params.ngroup)
        {
            if (MAIN(cr))
            {
                fprintf(this->log, "==== RAMD ==== GROMACS will be stopped after %ld steps.\n", step);
            }
            this->write_trajectory = true;
            gmx_set_stop_condition(StopCondition::Next);
        }
    }

    for (int g = 0; g < params.ngroup; ++g)
    {
        for (int i = 0; i < 3; ++i)
        {
            get_pull_coord_value(pull, g * 3 + i, pbc);
            apply_external_pull_coord_force(pull, g * 3 + i, direction[g][i] * params.group[g].force);
        }
    }
}

} // namespace gmx
