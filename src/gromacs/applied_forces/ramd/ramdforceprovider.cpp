/*! \internal \file
 * \brief
 * Implements force provider for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "ramdforceprovider.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

RAMDForceProvider::RAMDForceProvider(const RAMDParameters& parameters,
                                     PbcType               pbcType,
                                     const MDLogger&       logger,
                                     RAMDOutputProvider&   ramdOutputProvider) :
    parameters_(parameters),
    pbcType_(pbcType),
    logger_(logger),
    ramdOutputProvider_(ramdOutputProvider),    
    random_spherical_direction_generator(parameters.seed_, parameters.old_angle_dist_),
    direction(parameters.ngroups_),
    com_rec_prev(parameters.ngroups_),
    com_lig_prev(parameters.ngroups_)
{}

RAMDForceProvider::~RAMDForceProvider()
{}

void RAMDForceProvider::calculateForces(const ForceProviderInput& fInput,
                                        [[maybe_unused]] ForceProviderOutput* fOutput)
{
    t_pbc pbc;
    set_pbc(&pbc, this->pbcType_, fInput.box_);

    if (fInput.mpiComm_.isMainRank())
    {
        if (fInput.step_ == 0)
        {
            GMX_LOG(logger_.warning).appendText("==== RAMD ==== Initial COM calculation");
        }
        if (fInput.step_ % parameters_.eval_freq_ == 0)
        {
            ramdOutputProvider_.addTimePoint(fInput.t_);
        }
    }

    // Store COM positions for first evaluation
    if (fInput.step_ == 0)
    {
        for (int g = 0; g < parameters_.ngroups_; ++g)
        {
            // com_rec_prev[g] = pull->group[g * 2 + 1].x;
            // com_lig_prev[g] = pull->group[g * 2 + 2].x;
            com_rec_prev[g] = DVec(0.0, 0.0, 0.0);
            com_lig_prev[g] = DVec(0.0, 0.0, 0.0);

            if (fInput.mpiComm_.isMainRank())
            {
                DVec curr_dist_vect;
                pbc_dx_d(&pbc, com_lig_prev[g], com_rec_prev[g], curr_dist_vect);
                auto dist = std::sqrt(curr_dist_vect.norm2());
                ramdOutputProvider_.addCOMDistance(dist);
            }
        }
    }
    // else if (step % params.eval_freq == 0)
    // {
    //     if (mpiComm.isMainRank() and debug)
    //     {
    //         fprintf(debug, "==== RAMD ==== evaluation %ld\n", step);
    //     }
    //     for (int g = 0; g < params.ngroup; ++g)
    //     {
    //         DVec com_rec_curr = pull->group[g * 2 + 1].x;
    //         DVec com_lig_curr = pull->group[g * 2 + 2].x;
    //         DVec curr_dist_vect;
    //         pbc_dx_d(&pbc, com_lig_curr, com_rec_curr, curr_dist_vect);
    //         auto curr_dist = std::sqrt(curr_dist_vect.norm2());

    //         if (mpiComm.isMainRank() and debug)
    //         {
    //             fprintf(debug, "==== RAMD ==== group %d\n", g);
    //             fprintf(debug,
    //                     "==== RAMD ==== COM ligand position at [%g, %g, %g]\n",
    //                     com_lig_curr[0],
    //                     com_lig_curr[1],
    //                     com_lig_curr[2]);
    //             fprintf(debug,
    //                     "==== RAMD ==== COM receptor position at [%g, %g, %g]\n",
    //                     com_rec_curr[0],
    //                     com_rec_curr[1],
    //                     com_rec_curr[2]);
    //             fprintf(debug,
    //                     "==== RAMD ==== Distance between COM of receptor and COM of ligand is %g\n",
    //                     curr_dist);
    //         }

    //         if (mpiComm.isMainRank() and out)
    //         {
    //             fprintf(out, "\t%g", curr_dist);
    //         }

    //         if (curr_dist >= params.group[g].max_dist)
    //         {
    //             ligand_exited[g] = 1;
    //             if (!params.connected_ligands)
    //             {
    //                 direction[g] = DVec(0.0, 0.0, 0.0);
    //             }
    //             if (mpiComm.isMainRank())
    //             {
    //                 fprintf(this->log,
    //                         "==== RAMD ==== RAMD group %d has exited the binding site in step "
    //                         "%ld\n",
    //                         g,
    //                         step);
    //             }
    //         }
    //         else if (ligand_exited[g] == 1)
    //         {
    //             ligand_exited[g] = 0;
    //         }

    //         // difference of the COM ligand-receptor distance between current and the last evaluation step
    //         DVec walk_dist_vect;
    //         pbc_dx_d(&pbc, com_lig_curr - com_rec_curr, com_lig_prev[g] - com_rec_prev[g], walk_dist_vect);
    //         auto walk_dist = std::sqrt(walk_dist_vect.norm2());

    //         if (mpiComm.isMainRank() and debug)
    //         {
    //             fprintf(debug,
    //                     "==== RAMD ==== Previous COM ligand position at [%g, %g, %g]\n",
    //                     com_lig_prev[g][0],
    //                     com_lig_prev[g][1],
    //                     com_lig_prev[g][2]);
    //             fprintf(debug,
    //                     "==== RAMD ==== Previous COM receptor position at [%g, %g, %g]\n",
    //                     com_rec_prev[g][0],
    //                     com_rec_prev[g][1],
    //                     com_rec_prev[g][2]);
    //             fprintf(debug,
    //                     "==== RAMD ==== Change in receptor-ligand distance since last RAMD "
    //                     "evaluation "
    //                     "is %g\n",
    //                     walk_dist);
    //         }

    //         if (walk_dist < params.group[0].r_min_dist)
    //         {
    //             direction[g] = random_spherical_direction_generator();
    //             if (mpiComm.isMainRank() and debug)
    //             {
    //                 fprintf(debug,
    //                         "==== RAMD ==== New random direction is [%g, %g, %g]\n",
    //                         direction[g][0],
    //                         direction[g][1],
    //                         direction[g][2]);
    //             }
    //         }

    //         com_lig_prev[g] = com_lig_curr;
    //         com_rec_prev[g] = com_rec_curr;
    //     }
    // }

    // if (step % params.eval_freq == 0)
    // {
    //     if (mpiComm.isMainRank() and out)
    //     {
    //         fprintf(out, "\n");
    //         fflush(out);
    //     }

    //     // Exit if all ligand-receptor COM distances are larger than max_dist
    //     if (std::accumulate(ligand_exited.begin(), ligand_exited.end(), 0) == params.ngroup)
    //     {
    //         if (mpiComm.isMainRank())
    //         {
    //             fprintf(this->log, "==== RAMD ==== GROMACS will be stopped after %ld steps.\n", step);
    //         }
    //         this->write_trajectory = true;
    //         gmx_set_stop_condition(StopCondition::Next);
    //     }
    // }

    // for (int g = 0; g < params.ngroup; ++g)
    // {
    //     for (int i = 0; i < 3; ++i)
    //     {
    //         get_pull_coord_value(pull, g * 3 + i, pbc);
    //         apply_external_pull_coord_force(pull, g * 3 + i, direction[g][i] * params.group[g].force);
    //     }
    // }
}

} // namespace gmx
