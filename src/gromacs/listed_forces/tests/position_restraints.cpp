/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \internal \file
 * \brief Implements position restraint tests.
 *
 * \author Kevin Boyd <kevin44boyd@gmail.com>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "gromacs/listed_forces/position_restraints.h"

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Tolerance for float evaluation
constexpr float c_precisionTolerance = 1e-6;


class PositionRestraintsTest : public ::testing::TestWithParam<std::tuple<RefCoordScaling, PbcType>>
{
protected:
    std::vector<RVec> x_;
    std::vector<RVec> f_;

    matrix  box_;
    t_pbc   pbc_;
    PbcType pbcType_;

    t_nrnb                           nrnb_;
    t_forcerec                       fr_;
    InteractionDefinitions           idef_;
    gmx_enerdata_t                   enerd_;
    std::unique_ptr<ForceWithVirial> forceWithVirial_;
    RefCoordScaling                  refCoordScaling_;

    TestReferenceData    refData_;
    TestReferenceChecker checker_;

    PositionRestraintsTest() : idef_({}), enerd_(1, nullptr), checker_(refData_.rootChecker())
    {
        refCoordScaling_ = std::get<0>(GetParam());
        pbcType_         = std::get<1>(GetParam());

        clear_mat(box_);
        box_[0][0] = 0.9;
        box_[1][1] = 1.0;
        box_[2][2] = 1.1;
        set_pbc(&pbc_, pbcType_, box_);

        FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(1.0, c_precisionTolerance);
        checker_.setDefaultTolerance(tolerance);

        fr_.rc_scaling    = refCoordScaling_;
        fr_.pbcType       = pbcType_;
        fr_.posres_com[1] = 0.5;
    }

    //! Prepares the test with the coordinate and force constant input.
    void setValues(gmx::ArrayRef<const RVec> positions,
                   gmx::ArrayRef<const RVec> referencePositions,
                   gmx::ArrayRef<const RVec> forceConstants)
    {
        x_.resize(positions.size());
        std::copy(positions.begin(), positions.end(), x_.begin());

        for (size_t i = 0; i < positions.size(); i++)
        {
            // First item is "type" - each atom will have a different forceparam type
            // Second item is index - we'll just go from 0.
            idef_.il[F_POSRES].iatoms.push_back(i);
            idef_.il[F_POSRES].iatoms.push_back(i);

            auto& entry = idef_.iparams_posres.emplace_back();
            copy_rvec(referencePositions[i], entry.posres.pos0A);
            copy_rvec(forceConstants[i], entry.posres.fcA);
            clear_rvec(entry.posres.pos0B);
            clear_rvec(entry.posres.fcB);
        }
        f_.resize(x_.size(), { 0, 0, 0 });
        forceWithVirial_ = std::make_unique<ForceWithVirial>(f_, /*computeVirial=*/true);
    }
};

std::array<real, static_cast<size_t>(FreeEnergyPerturbationCouplingType::Count)> c_emptyLambdas = { { 0 } };

TEST_P(PositionRestraintsTest, BasicPosResNoFreeEnergy)
{
    SCOPED_TRACE(formatString("Testing PBC type: %s, refcoord type: %s",
                              c_pbcTypeNames[pbcType_].c_str(),
                              enumValueToString(refCoordScaling_)));
    const std::vector<RVec> positions          = { { 0.0, 0.0, 0.0 }, { 0.4, 0.5, 0.6 } };
    const std::vector<RVec> referencePositions = { { 0.0, 0.0, 0.0 }, { 0.5, 0.6, 0.0 } };
    const std::vector<RVec> forceConstants     = { { 1000, 500, 250 }, { 0, 200, 400 } };
    setValues(positions, referencePositions, forceConstants);
    posres_wrapper(&nrnb_,
                   idef_,
                   &pbc_,
                   as_rvec_array(x_.data()),
                   &enerd_,
                   c_emptyLambdas,
                   &fr_,
                   forceWithVirial_.get());
    checker_.checkSequence(
            std::begin(forceWithVirial_->force_), std::end(forceWithVirial_->force_), "Forces");
    checker_.checkSequenceArray(3, forceWithVirial_->getVirial(), "Virial contribution");
    checker_.checkReal(enerd_.term[F_POSRES], "Potential energy");
}

//! PBC values for testing
std::vector<PbcType> c_pbcForTests = { PbcType::No, PbcType::XY, PbcType::Xyz };
//! Reference Coordinate Scaling values for testing
std::vector<RefCoordScaling> c_refCoordScalingForTests = { RefCoordScaling::No,
                                                           RefCoordScaling::Com,
                                                           RefCoordScaling::All };

INSTANTIATE_TEST_SUITE_P(PosResBasicTest,
                         PositionRestraintsTest,
                         ::testing::Combine(::testing::ValuesIn(c_refCoordScalingForTests),
                                            ::testing::ValuesIn(c_pbcForTests)));

} // namespace

} // namespace test

} // namespace gmx
