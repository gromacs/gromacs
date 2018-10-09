/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*!\file
 * \internal
 * \brief
 * Tests for outputmanager
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */


#include "gmxpre.h"

#include "gromacs/coordinateio/builder.h"

#include "config.h"

#include <utility>

#include "gromacs/coordinateio/tests/coordinate_test.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Test fixture to test different file types are supported by the OutputManager.
 */
class OutputManagerBuilderTest : public ModuleTest
{
    public:
        /*! \brief
         * Test basic behaviour without special requirements.
         *
         * \param[in] filename Name for output file.
         */
        void basicTest(const char *filename)
        {
            addTopology();

            OutputRequirements       requirements;

            runTest(filename, requirements);
        }
        /*! \brief
         * Test with extra requirements.
         *
         * \param[in] filename Name for output file.
         * \param[in] requirements Specify extra reqs for output.
         */
        void testWithRequirements(const char        *filename,
                                  OutputRequirements requirements)
        {
            addTopology();
            runTest(filename, std::move(requirements));
        }

};

TEST_P(OutputManagerBuilderTest, WorksWithFormats)
{
    EXPECT_NO_THROW(basicTest(GetParam()));
}

TEST_F(OutputManagerBuilderTest, RejectsWrongFiletype)
{
    EXPECT_THROW(basicTest("test.xvg"), InvalidInputError);
}

TEST_F(OutputManagerBuilderTest, BuilderFailsWithPdbAndNoAtoms)
{
    OutputRequirements       requirements;
    requirements.atoms = ChangeAtomsType::efUserNo;
    EXPECT_THROW(testWithRequirements("test.pdb",
                                      requirements),
                 InconsistentInputError);
    EXPECT_EQ(output_.get(), nullptr);
}

TEST_F(OutputManagerBuilderTest, BuilderFailsWithGroAndNoAtoms)
{
    OutputRequirements       requirements;
    requirements.atoms = ChangeAtomsType::efUserNo;
    EXPECT_THROW(testWithRequirements("test.gro",
                                      requirements),
                 InconsistentInputError);
    EXPECT_EQ(output_.get(), nullptr);
}

TEST_F(OutputManagerBuilderTest, BuilderImplictlyAddsAtoms)
{
    OutputRequirements       requirements;
    requirements.atoms = ChangeAtomsType::efUnchanged;
    {
        EXPECT_NO_THROW(testWithRequirements("test.pdb",
                                             requirements));
        EXPECT_NE(output_.get(), nullptr);
    }
    {
        EXPECT_NO_THROW(testWithRequirements("test.gro",
                                             requirements));
        EXPECT_NE(output_.get(), nullptr);
    }
}

TEST_F(OutputManagerBuilderTest, TNGNeedsTopology)
{
    OutputRequirements requirements;
    EXPECT_THROW(runTest("test.tng", requirements),
                 InconsistentInputError);
}

/*!\brief
 * Character array of different file names to test.
 */
const char *const trajectoryFileNames[] = {
    "spc2-traj.trr",
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.xtc",
    "spc2-traj.pdb",
    "spc2-traj.gro",
    "spc2-traj.g96"
};

INSTANTIATE_TEST_CASE_P(OutputManagerFileFormats,
                        OutputManagerBuilderTest, ::testing::ValuesIn(trajectoryFileNames));

} // namespace test

} // namespace gmx
