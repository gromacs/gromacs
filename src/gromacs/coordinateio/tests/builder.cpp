/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*!\internal
 * \file
 * \brief
 * Tests for outputmanager
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "config.h"

#include <string>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/coordinatefileenums.h"
#include "gromacs/coordinateio/requirements.h"
#include "gromacs/coordinateio/tests/coordinate_test.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Test fixture to test different file types are supported by the CoordinateFile.
 */
class TrajectoryFrameWriterTest : public ModuleTest
{
public:
    /*! \brief
     * Test basic behaviour without special requirements.
     *
     * \param[in] filename Name for output file.
     */
    void basicTest(const char* filename)
    {
        addTopology();

        OutputRequirements requirements;

        runTest(filename, requirements);
    }
    /*! \brief
     * Test with extra requirements.
     *
     * \param[in] filename Name for output file.
     * \param[in] requirements Specify extra reqs for output.
     */
    void testWithRequirements(const char* filename, const OutputRequirements& requirements)
    {
        addTopology();
        runTest(filename, requirements);
    }
};

TEST_P(TrajectoryFrameWriterTest, WorksWithFormats)
{
    EXPECT_NO_THROW(basicTest(GetParam()));
}

TEST_F(TrajectoryFrameWriterTest, RejectsWrongFiletype)
{
    EXPECT_THROW(basicTest("test.xvg"), InvalidInputError);
}

TEST_F(TrajectoryFrameWriterTest, BuilderFailsWithPdbAndNoAtoms)
{
    OutputRequirements requirements;
    requirements.atoms = ChangeAtomsType::Never;
    EXPECT_THROW(testWithRequirements("test.pdb", requirements), InconsistentInputError);
}

TEST_F(TrajectoryFrameWriterTest, BuilderFailsWithGroAndNoAtoms)
{
    OutputRequirements requirements;
    requirements.atoms = ChangeAtomsType::Never;
    EXPECT_THROW(testWithRequirements("test.gro", requirements), InconsistentInputError);
}

TEST_F(TrajectoryFrameWriterTest, BuilderImplictlyAddsAtoms)
{
    OutputRequirements requirements;
    requirements.atoms = ChangeAtomsType::PreservedIfPresent;
    {
        EXPECT_NO_THROW(testWithRequirements("test.pdb", requirements));
    }
    {
        EXPECT_NO_THROW(testWithRequirements("test.gro", requirements));
    }
}

#if GMX_USE_TNG
TEST_F(TrajectoryFrameWriterTest, TNGOutputWorks)
{
    OutputRequirements requirements;
    runTest("test.tng", requirements);
}
#endif

/*!\brief
 * Character array of different file names to test.
 */
const char* const trajectoryFileNames[] = { "spc2-traj.trr",
#if GMX_USE_TNG
                                            "spc2-traj.tng",
#endif
                                            "spc2-traj.xtc", "spc2-traj.pdb",
                                            "spc2-traj.gro", "spc2-traj.g96" };

INSTANTIATE_TEST_SUITE_P(CoordinateFileFileFormats,
                         TrajectoryFrameWriterTest,
                         ::testing::ValuesIn(trajectoryFileNames));

} // namespace test

} // namespace gmx
