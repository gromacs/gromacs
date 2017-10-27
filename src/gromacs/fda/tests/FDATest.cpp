/*
 * FDA.cpp
 *
 *  Created on: Sep 15, 2017
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <gtest/gtest.h>
#include "gromacs/fda/FDA.h"

namespace fda {

TEST(FDATest, add_angle)
{
    FDASettings fda_settings;
    fda_settings.atom_based_result_type = ResultType::NO;
    fda_settings.residue_based_result_type = ResultType::PUNCTUAL_STRESS;
    fda_settings.one_pair = OnePair::SUMMED;
    fda_settings.v2s = Vector2Scalar::NORM;
    fda_settings.residues_renumber = ResiduesRenumber::AUTO;
    fda_settings.no_end_zeros = false;
    fda_settings.syslen_atoms = 0;
    fda_settings.syslen_residues = 0;
    fda_settings.time_averaging_period = 1;
    fda_settings.type = InteractionType_NONE;
    fda_settings.nonbonded_exclusion_on = true;
    fda_settings.bonded_exclusion_on = true;
    fda_settings.index_group1 = -1;
    fda_settings.index_group2 = -1;
    fda_settings.groups = nullptr;
    fda_settings.groupnames = nullptr;

    FDA fda(fda_settings);
    rvec f_i{0.0, 0.0, 0.0};
    rvec f_j{0.0, 0.0, 0.0};
    rvec f_k{0.0, 0.0, 0.0};
    fda.add_angle(0, 1, 2, f_i, f_j, f_k);
}

} // namespace fda
