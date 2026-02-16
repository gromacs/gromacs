/*! \internal \file
 * \brief
 * Implements force provider for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "ramdforceprovider.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

RAMDForceProvider::RAMDForceProvider(const RAMDParameters& parameters,
                                     const std::vector<std::unique_ptr<LocalAtomSet>>& localAtoms,
                                     const gmx_mtop_t& topology,
                                     PbcType pbcType,
                                     const MDLogger& logger,
                                     RAMDOutputProvider& ramdOutputProvider) :
    parameters_(parameters),
    localAtoms_(localAtoms),
    topology_(topology),
    pbcType_(pbcType),
    logger_(logger),
    ramdOutputProvider_(ramdOutputProvider),
    random_spherical_direction_generator(parameters.seed_, parameters.old_angle_dist_),
    direction_(parameters.groups_.size()),
    com_rec_prev_(parameters.groups_.size()),
    com_lig_prev_(parameters.groups_.size()),
    ligand_exited_(parameters.groups_.size(), 0),
    write_trajectory_(false),
    mTopLookUp_(topology)
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
            ramdOutputProvider_.addTime(fInput.t_);
        }
    }
    
    // Evaluate RAMD every eval_freq steps
    if (fInput.step_ % parameters_.eval_freq_ == 0)
    {
        GMX_LOG(logger_.warning).appendText("==== RAMD ==== evaluation ").appendText(std::to_string(fInput.step_));
        for (int g = 0; g < parameters_.ngroups_; ++g)
        {
            std::string logPrefix = "==== RAMD group " + std::to_string(g) + " ====";
            DVec com_rec_curr = calc_com(fInput.x_, parameters_.groups_[g].receptor_indices_);
            DVec com_lig_curr = calc_com(fInput.x_, parameters_.groups_[g].ligand_indices_);
            DVec curr_dist_vect;
            pbc_dx_d(&pbc, com_lig_curr, com_rec_curr, curr_dist_vect);
            real curr_dist = std::sqrt(curr_dist_vect.norm2());
            ramdOutputProvider_.addDistance(curr_dist);

            GMX_LOG(logger_.info).appendText(logPrefix + "Current COM ligand position at [" +
                                            std::to_string(com_lig_curr[0]) + ", " +
                                            std::to_string(com_lig_curr[1]) + ", " +
                                            std::to_string(com_lig_curr[2]) + "]");
            GMX_LOG(logger_.info).appendText(logPrefix + "Current COM receptor position at [" +
                                            std::to_string(com_rec_curr[0]) + ", " +
                                            std::to_string(com_rec_curr[1]) + ", " +
                                            std::to_string(com_rec_curr[2]) + "]");
            GMX_LOG(logger_.info).appendText(logPrefix + "Distance between COM of receptor and COM of ligand is "
                + std::to_string(curr_dist) + "\n");

            if (curr_dist >= parameters_.groups_[g].max_dist_)
            {
                ligand_exited_[g] = 1;
                if (!parameters_.connected_ligands_)
                {
                    direction_[g] = DVec(0.0, 0.0, 0.0);
                }
                GMX_LOG(logger_.info).appendText(logPrefix + "group " + std::to_string(g) + " has exited the binding site in step " + std::to_string(fInput.step_));
            }
            else if (ligand_exited_[g] == 1)
            {
                ligand_exited_[g] = 0;
            }

            real walk_dist = 0.0;
            if (fInput.step_ != 0)
            {
                DVec walk_dist_vect;
                // difference of the COM ligand-receptor distance between current and the last evaluation step
                pbc_dx_d(&pbc, com_lig_curr - com_rec_curr, com_lig_prev_[g] - com_rec_prev_[g], walk_dist_vect);
                walk_dist = std::sqrt(walk_dist_vect.norm2());

                GMX_LOG(logger_.info).appendText(logPrefix + "Previous COM ligand position at [" +
                                                std::to_string(com_lig_prev_[g][0]) + ", " +
                                                std::to_string(com_lig_prev_[g][1]) + ", " +
                                                std::to_string(com_lig_prev_[g][2]) + "]");
                GMX_LOG(logger_.info).appendText(logPrefix + "Previous COM receptor position at [" +
                                                std::to_string(com_rec_prev_[g][0]) + ", " +
                                                std::to_string(com_rec_prev_[g][1]) + ", " +
                                                std::to_string(com_rec_prev_[g][2]) + "]");
                GMX_LOG(logger_.info).appendText(logPrefix + "Change in receptor-ligand"
                    " distance since last RAMD evaluation is " + std::to_string(walk_dist) + "\n");
            }

            if (walk_dist < parameters_.groups_[0].r_min_dist_)
            {
                direction_[g] = random_spherical_direction_generator();
                GMX_LOG(logger_.warning).appendText("==== RAMD ==== New random direction is [" +
                                                std::to_string(direction_[g][0]) + ", " +
                                                std::to_string(direction_[g][1]) + ", " +
                                                std::to_string(direction_[g][2]) + "]");
            }

            com_lig_prev_[g] = com_lig_curr;
            com_rec_prev_[g] = com_rec_curr;
        }

        // Finish the RAMD output line
        ramdOutputProvider_.newLine();
        ramdOutputProvider_.flush();

        // Exit if all ligand-receptor COM distances are larger than max_dist
        if (std::accumulate(ligand_exited_.begin(), ligand_exited_.end(), 0) == parameters_.ngroups_)
        {
            GMX_LOG(logger_.warning).appendTextFormatted("==== RAMD ==== GROMACS will be stopped after %ld steps.\n", fInput.step_);
            write_trajectory_ = true;
            gmx_set_stop_condition(StopCondition::Next);
        }
    }

    // Apply forces to ligand atoms
    for (size_t g = 0; g < parameters_.groups_.size(); ++g)
    {
        for (size_t i = 0; i < localAtoms_[g]->numAtomsLocal(); ++i)
        {
            fOutput->forceWithVirial_.force_[localAtoms_[g]->localIndex()[i]][XX] += direction_[g][XX] * parameters_.groups_[g].force_;
            fOutput->forceWithVirial_.force_[localAtoms_[g]->localIndex()[i]][YY] += direction_[g][YY] * parameters_.groups_[g].force_;
            fOutput->forceWithVirial_.force_[localAtoms_[g]->localIndex()[i]][ZZ] += direction_[g][ZZ] * parameters_.groups_[g].force_;
        }
    }
}

} // namespace gmx
