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

#include "config.h"

#include "gromacs/coordinateio/tests/coordinate_test.h"

namespace gmx
{

namespace test
{

/*! \brief
 * Method to test the addition of a dummy module to the OutputManager.
 *
 * Allocates a new minimal outputadapter module with a specified \p flag
 * and tries to add it to the OutputManager. Used for testing that adding
 * works and that only matching modules are accepted.
 *
 * \param[in] filename                      Name of file used to create dummy output manager.
 * \param[in] dummyModuleRequirementsFlag   Flag setting propagated to the dummy module to check
 *                                          that mismatched flags are rejected and matching
 *                                          flags accepted.
 * \param[in] dummyModuleId                 Flag representing id to test that multipel modules work.
 */
static void addEmptyModuleToOutputManager(std::string   filename,
                                          unsigned long dummyModuleRequirementsFlag,
                                          unsigned long dummyModuleId)
{
    gmx_mtop_t               dummyTopology;

    OutputAdapters           adapters;
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(dummyModuleRequirementsFlag, dummyModuleId));

    OutputManagerPointer output = createMinimalOutputManager(filename,
                                                             &dummyTopology,
                                                             std::move(adapters));
}

/*!\brief
 * Create minimal OutputManager using the provided builder.
 *
 * \param[in] filename Name of file to create OutputManager for.
 * \param[in] dummyTopology Pointer to input top or null.
 * \throws InconsistentInputError When builder can not create the OutputManager.
 * \returns unique_ptr to new OutputManager object.
 */
OutputManagerPointer createMinimalOutputManager(const std::string       &filename,
                                                const gmx_mtop_t        *dummyTopology,
                                                OutputAdapters           adapters)
{
    Selection                dummySelection;

    return createOutputManager(dummyTopology,
                               dummySelection,
                               filename,
                               std::move(adapters));
}

/*!\brief
 * Test fixture to test different file types are supported by the OutputManager.
 */
class OutputManagerBuilderTest : public gmx::test::CommandLineTestBase,
                                 public ::testing::WithParamInterface<const char *>
{
    public:
        void runTest(const char *filename)
        {
            //! Pointer to new OutputManager object.
            OutputManagerPointer     output;
            //! Dummy topology to use to create OutputManager.
            gmx_mtop_t               dummyTopology;
            OutputAdapters           adapters;
            EXPECT_NO_THROW(output = createMinimalOutputManager(filename,
                                                                &dummyTopology,
                                                                std::move(adapters)));
        }
};

TEST(OutputManagerTest, CanAddFittingModule)
{
    unsigned long moduleId               = efDummyModule;
    unsigned long moduleRequirementsFlag = efBaseOutputManager;
    EXPECT_NO_THROW(addEmptyModuleToOutputManager(
                            "test.pdb",
                            moduleRequirementsFlag,
                            moduleId));
}

TEST(OutputManagerTest, CannotAddMismatchedModule)
{
    unsigned long moduleId               = efDummyModule;
    unsigned long moduleRequirementsFlag = efChangeVelocityModule;
    EXPECT_THROW(addEmptyModuleToOutputManager(
                         "test.pdb",
                         moduleRequirementsFlag,
                         moduleId),
                 InconsistentInputError);
}

TEST(OutputManagerTest, CanAddMultipleModules)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    std::string              filename = "test.tng";
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeAtomInformationModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeVelocityModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeForceModule));
    EXPECT_NO_THROW(output = createMinimalOutputManager(filename,
                                                        &dummyTopology,
                                                        std::move(adapters)));
    EXPECT_NE(output.get(), nullptr);
}

TEST(OutputManagerTest, CannotAddSameModuleTwice1)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    std::string              filename = "test.tng";
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efDummyModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efDummyModule));
    EXPECT_THROW(output = createMinimalOutputManager(filename,
                                                     &dummyTopology,
                                                     std::move(adapters)),
                 InternalError);
    EXPECT_EQ(output.get(), nullptr);
}

TEST(OutputManagerTest, CannotAddSameModuleTwice2)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    std::string              filename = "test.tng";
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efDummyModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeVelocityModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efDummyModule));
    EXPECT_THROW(output = createMinimalOutputManager(filename,
                                                     &dummyTopology,
                                                     std::move(adapters)),
                 InternalError);
    EXPECT_EQ(output.get(), nullptr);
}

TEST(OutputManagerTest, CannotViolateModuleOrder)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    std::string              filename = "test.tng";
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeCoordinateSelectionModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeVelocityModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efBaseOutputManager, efChangeAtomInformationModule));
    EXPECT_THROW(output = createMinimalOutputManager(filename,
                                                     &dummyTopology,
                                                     std::move(adapters)),
                 InternalError);
    EXPECT_EQ(output.get(), nullptr);
}

TEST(OutputManagerTest, CannotHideMismatchedModule)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    std::string              filename = "test.gro";
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efChangeAtomInformationModule, efChangeAtomInformationModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efChangeForceModule, efChangeForceModule));
    adapters.emplace_back(compat::make_unique<DummyOutputModule>(
                                  efChangeVelocityModule, efChangeVelocityModule));
    EXPECT_THROW(output = createMinimalOutputManager(filename,
                                                     &dummyTopology,
                                                     std::move(adapters)),
                 InconsistentInputError);
    EXPECT_EQ(output.get(), nullptr);
}


TEST_P(OutputManagerBuilderTest, WorksWithAllFormats)
{
    runTest(GetParam());
}

TEST(OutputManagerTest, BuilderRejectsWrongFiletype)
{
    OutputManagerPointer     output;
    gmx_mtop_t               dummyTopology;
    OutputAdapters           adapters;
    EXPECT_THROW(output = createMinimalOutputManager("test.xvg",
                                                     &dummyTopology,
                                                     std::move(adapters)),
                 InvalidInputError);
    EXPECT_EQ(output.get(), nullptr);
}

TEST(OutputManagerTest, BuilderFailsWhenNeedingTopology)
{
    OutputManagerPointer     output;
    gmx_mtop_t              *dummyTopology = nullptr;
    std::vector<std::string> filenames     = {"test.pdb", "test.gro", "test.tng"};
    for (const auto &name : filenames)
    {
        OutputAdapters adapters;
        EXPECT_THROW(output = createMinimalOutputManager(name,
                                                         dummyTopology,
                                                         std::move(adapters)),
                     InconsistentInputError);
        EXPECT_EQ(output.get(), nullptr);
    }
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
    "spc2-traj.gro",
    "spc2-traj.pdb",
    "spc2-traj.g96"
};

INSTANTIATE_TEST_CASE_P(OutputManagerFileFormats,
                        OutputManagerBuilderTest, ::testing::ValuesIn(trajectoryFileNames));

} // namespace test

} // namespace gmx
