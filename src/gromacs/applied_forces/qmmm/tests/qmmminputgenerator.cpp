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
 * \brief
 * Tests for QMMMInputGenerator class for QMMM MDModule
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/qmmm/qmmminputgenerator.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/qmmm/qmmmtypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

class QMMMInputGeneratorTest : public ::testing::Test
{
public:
    //! Reference input 2x TIP3P waters (first one QM) and no Link atoms
    void initParameters2WatersNoLink()
    {
        parameters_.qmIndices_   = { 0, 1, 2 };
        parameters_.mmIndices_   = { 3, 4, 5 };
        parameters_.atomNumbers_ = { 8, 1, 1, 8, 1, 1 };
        q_                       = { 0.0, 0.0, 0.0, -0.834, 0.417, 0.417 };
        pbc_                     = PbcType::Xyz;
        matrix box1_             = { { 2.0, 0.0, 0.0 }, { 0.0, 2.0, 0.0 }, { 0.0, 0.0, 2.0 } };
        copy_mat(box1_, box_);
        coords_ = { RVec(0.005, 0.600, 0.244), RVec(-0.017, 0.690, 0.270),
                    RVec(0.051, 0.610, 0.161), RVec(0.975, 0.631, 1.569),
                    RVec(0.951, 0.615, 1.661), RVec(0.916, 0.701, 1.541) };
    }

    //! Reference input 2x TIP3P waters (first two atoms are QM) and one link should be generated
    void initParameters2WatersWithLink()
    {
        parameters_.qmIndices_   = { 0, 1 };
        parameters_.mmIndices_   = { 2, 3, 4, 5 };
        parameters_.atomNumbers_ = { 8, 1, 1, 8, 1, 1 };
        parameters_.link_        = { { 0, 2 } };
        q_                       = { 0.0, 0.0, 0.417, -0.834, 0.417, 0.417 };
        pbc_                     = PbcType::Xyz;
        matrix box1_             = { { 2.0, 0.0, 0.0 }, { 0.0, 2.0, 0.0 }, { 0.0, 0.0, 2.0 } };
        copy_mat(box1_, box_);
        coords_ = { RVec(0.005, 0.600, 0.244), RVec(-0.017, 0.690, 0.270),
                    RVec(0.051, 0.610, 0.161), RVec(0.975, 0.631, 1.569),
                    RVec(0.951, 0.615, 1.661), RVec(0.916, 0.701, 1.541) };
    }

protected:
    QMMMParameters parameters_;
    PbcType        pbc_;
    matrix         box_ = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };

    std::vector<real> q_;
    std::vector<RVec> coords_;
};


TEST_F(QMMMInputGeneratorTest, CanConstruct)
{
    initParameters2WatersNoLink();
    EXPECT_NO_THROW(QMMMInputGenerator inpGen(parameters_, pbc_, box_, q_, coords_));
}

TEST_F(QMMMInputGeneratorTest, TwoWatersPBENoLink)
{
    // Create Input Generator
    initParameters2WatersNoLink();
    QMMMInputGenerator inpGen(parameters_, pbc_, box_, q_, coords_);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // Tolerance of all coordinates and vectors should be 1E-3 (as in gro or pdb files)
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    checker.checkVector(inpGen.qmTrans(), "QM Translation");
    checker.checkVector(inpGen.qmBox()[0], "QM Box Vector 1");
    checker.checkVector(inpGen.qmBox()[1], "QM Box Vector 2");
    checker.checkVector(inpGen.qmBox()[2], "QM Box Vector 3");
    checker.checkString(inpGen.generateCP2KInput(), "Input");
    checker.checkString(inpGen.generateCP2KPdb(), "PDB");
}

TEST_F(QMMMInputGeneratorTest, TwoWatersPBEWithLink)
{
    // Create Input Generator
    initParameters2WatersWithLink();
    QMMMInputGenerator inpGen(parameters_, pbc_, box_, q_, coords_);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // Tolerance of all coordinates and vectors should be 1E-3 (as in gro or pdb files)
    checker.setDefaultTolerance(gmx::test::absoluteTolerance(0.001));
    checker.checkVector(inpGen.qmTrans(), "QM Translation");
    checker.checkVector(inpGen.qmBox()[0], "QM Box Vector 1");
    checker.checkVector(inpGen.qmBox()[1], "QM Box Vector 2");
    checker.checkVector(inpGen.qmBox()[2], "QM Box Vector 3");
    checker.checkString(inpGen.generateCP2KInput(), "Input");
    checker.checkString(inpGen.generateCP2KPdb(), "PDB");
}

} // namespace gmx
