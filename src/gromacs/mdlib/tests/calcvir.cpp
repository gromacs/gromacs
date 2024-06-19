/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Tests for virial calculation.
 *
 * \author Joe Jordan <ejjordan@kth.se>
 */
#include "gmxpre.h"

#include "gromacs/mdlib/calcvir.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{


class CalcvirTest : public ::testing::Test
{
public:
    TestReferenceData      refData_;
    TestReferenceChecker   checker_;
    std::vector<gmx::RVec> coordinates_;
    std::vector<gmx::RVec> forces_;
    int                    numVirialAtoms_;
    tensor                 virial_{ { 0 } };
    FloatingPointTolerance tolerances =
            FloatingPointTolerance(1e-16, 1.0e-16, 1e-16, 1.0e-7, 1000, 100, false);

    CalcvirTest() :
        checker_(refData_.rootChecker()),
        coordinates_{ { 1, 2, 3 },
                      {
                              4,
                              5,
                              6,
                      },
                      { 7, 8, 9 },
                      { 1, 5, 9 } },
        forces_{ { 0.9, 0.1, 0.3 }, { 0.4, 0.7, 0.2 }, { 0.5, 1, 0.6 }, { 0.9, 0.7, 0.6 } },
        numVirialAtoms_(coordinates_.size())
    {
        checker_.setDefaultTolerance(tolerances);
    }

private:
};

TEST_F(CalcvirTest, CanCalculateVirialAllAtomsInBox)
{
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, false, nullptr);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

TEST_F(CalcvirTest, CanCalculateVirialAllAtomsInBoxScrew)
{

    const matrix box = { { 10, 0, 0 }, { 0, 12, 0 }, { 0, 0, 14 } };
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, true, box);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

TEST_F(CalcvirTest, CanCalculateVirialAtomsOutOfBoxScrewX)
{

    const matrix box = { { 5, 0, 0 }, { 0, 10, 0 }, { 0, 0, 12 } };
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, true, box);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

TEST_F(CalcvirTest, CanCalculateVirialAtomsOutOfBoxScrewY)
{

    const matrix box = { { 10, 0, 0 }, { 0, 5, 0 }, { 0, 0, 12 } };
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, true, box);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

TEST_F(CalcvirTest, CanCalculateVirialAtomsOutOfBoxScrewZ)
{

    const matrix box = { { 10, 0, 0 }, { 0, 12, 0 }, { 0, 0, 5 } };
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, true, box);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

TEST_F(CalcvirTest, CanCalculateVirialAtomsOutOfBoxScrewXYZ)
{

    const matrix box = { { 4, 0, 0 }, { 0, 5, 0 }, { 0, 0, 6 } };
    calc_vir(numVirialAtoms_, as_rvec_array(coordinates_.data()), as_rvec_array(forces_.data()), virial_, true, box);

    checker_.checkVector(virial_[0], "Virial x");
    checker_.checkVector(virial_[1], "Virial y");
    checker_.checkVector(virial_[2], "Virial z");
}

} // namespace
} // namespace test
} // namespace gmx
