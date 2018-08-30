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
 * \libinternal
 * \brief
 * Helpers and data for outputadapter module tests.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_MODULE_H
#define GMX_COORDINATEIO_TESTS_MODULE_H

#include "gmxpre.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/coordinateio/modules/outputselector.h"
#include "gromacs/coordinateio/modules/setatoms.h"
#include "gromacs/coordinateio/modules/setbox.h"
#include "gromacs/coordinateio/modules/setforces.h"
#include "gromacs/coordinateio/modules/setprecision.h"
#include "gromacs/coordinateio/modules/settime.h"
#include "gromacs/coordinateio/modules/setvelocities.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include "gromacs/coordinateio/tests/coordinate_test.h"
#include "testutils/testfilemanager.h"


namespace gmx
{

namespace test
{


/*! \libinternal \brief
 * Test fixture to test matching file types for modules.
 */
class ModuleTest : public gmx::test::CommandLineTestBase,
                   public ::testing::WithParamInterface<const char *>
{
    public:
        /*! \brief
         * Run the builder to create an OutputManager during tests.
         *
         * \param[in] filename Name for output file, to determine filetype during construction.
         * \param[in] adapters Container of outputadapters to add to the OutputManager.
         */
        void runTest(const char *filename, CoordinateOutputAdapters adapters)
        {
            output_ = createMinimalOutputManager(filename,
                                                 &dummyTopology_,
                                                 std::move(adapters));
        }
        //! Pointer to new OutputManager object.
        OutputManagerPointer output_;
        //! Dummy topology to use to create OutputManager.
        gmx_mtop_t           dummyTopology_;
};

/*! \libinternal \brief
 * Helper class for tests that need an initialized selection.
 */
class ModuleSelection
{
    public:
        ModuleSelection() : manager_(&sc_)
        {
            options_.addManager(&manager_);
            sc_.setReferencePosType("atom");
            sc_.setOutputPosType("atom");
            top_.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.gro"));
        }

        /*! \brief
         * Method to add a valid selection option to the Module, or an invalid
         * one for testing.
         *
         * \param[in] mtop Topology to use.
         * \param[in] sel Selection to add option to.
         * \param[in] useValid Set if the added selection should be valid for the module.
         */
        void addOptionForSelection(gmx_mtop_t *mtop, Selection *sel, bool useValid);

    private:
        //! Local selection collection.
        SelectionCollection    sc_;
        //! Local selection manager.
        SelectionOptionManager manager_;
        //! Local options manager.
        Options                options_;
        //! Local topology information.
        TopologyInformation    top_;

};

inline void
ModuleSelection::addOptionForSelection(gmx_mtop_t *mtop, Selection *sel, bool useValid)
{
    mtop = top_.mtop();

    sc_.setTopology(mtop, 0);

    if (useValid)
    {
        options_.addOption(SelectionOption("sel").store(sel).onlyAtoms());
    }
    else
    {
        options_.addOption(SelectionOption("sel").store(sel).dynamicMask());
    }

    OptionsAssigner assigner(&options_);
    assigner.start();
    assigner.startOption("sel");
    if (useValid)
    {
        assigner.appendValue("all");
    }
    else
    {
        assigner.appendValue("res_cog of all");
    }
    assigner.finishOption();
    assigner.finish();

    ASSERT_TRUE(sel->isValid());
    EXPECT_NO_THROW(sc_.compile());
}

/*!\libinternal \brief  Helper to test supported file names. */
class SetAtomsSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Local atoms
            t_atoms                 *atoms = nullptr;

            adapters.emplace_back(compat::make_unique<SetAtoms>(ChangeSettingType::efUserNo, atoms));

            EXPECT_NO_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetAtomsUnSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Local atoms
            t_atoms                 *atoms = nullptr;

            adapters.emplace_back(compat::make_unique<SetAtoms>(ChangeSettingType::efUserNo, atoms));

            EXPECT_ANY_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class AnyOutputSupportedFiles : public ModuleTest,
                                public ModuleSelection
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Local atoms
            Selection                sel;
            //! Local box
            matrix                   box;
            //! Value for startTime
            real                     startTime = 0;
            //! Value for timeStep
            real                     timeStep = 1;

            clear_mat(box);

            addOptionForSelection(&dummyTopology_, &sel, true);

            adapters.emplace_back(compat::make_unique<OutputSelector>(sel));
            adapters.emplace_back(compat::make_unique<SetBox>(box));
            adapters.emplace_back(compat::make_unique<SetTime>(startTime, timeStep, ChangeFrameTimeType::efUnchanged));

            EXPECT_NO_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test that invalid selection is rejected */
class OutputSelectorDeathTest : public ModuleTest,
                                public ModuleSelection
{
    public:
        void prepareTest()
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Local atoms
            Selection                sel;

            addOptionForSelection(&dummyTopology_, &sel, false);

            ASSERT_DEATH_IF_SUPPORTED(
                    adapters.emplace_back(compat::make_unique<OutputSelector>(sel)),
                    "Need a valid selection out of simple atom indices");
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetVelocitySupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;

            adapters.emplace_back(compat::make_unique<SetVelocities>(ChangeSettingType::efUserNo));

            EXPECT_NO_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetVelocityUnSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;

            adapters.emplace_back(compat::make_unique<SetVelocities>(ChangeSettingType::efUserNo));

            EXPECT_ANY_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetForceSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;

            adapters.emplace_back(compat::make_unique<SetForces>(ChangeSettingType::efUserNo));

            EXPECT_NO_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetForceUnSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;

            adapters.emplace_back(compat::make_unique<SetForces>(ChangeSettingType::efUserNo));

            EXPECT_ANY_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetPrecisionSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Value for new precision.
            int precision = 3;

            adapters.emplace_back(compat::make_unique<SetPrecision>(precision));

            EXPECT_NO_THROW(runTest(filename, std::move(adapters)));
        }
};

/*!\libinternal \brief  Helper to test supported file names. */
class SetPrecisionUnSupportedFiles : public ModuleTest
{
    public:
        void prepareTest(const char *filename)
        {
            //! Storage for frameadapters.
            CoordinateOutputAdapters adapters;
            //! Value for new precision.
            int precision = 3;

            adapters.emplace_back(compat::make_unique<SetPrecision>(precision));

            EXPECT_ANY_THROW(runTest(filename, std::move(adapters)));
        }
};

//! Names here work for setAtoms module
const char *const setAtomsSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.gro",
    "spc2-traj.pdb",
};

//! Names here don't work for setAtoms module
const char *const setAtomsUnSupported[] = {
    "spc2-traj.trr",
    "spc2-traj.xtc",
    "spc2-traj.g96"
};

//! Names here work for stuff that has no specific requirements.
const char *const anySupported[] = {
    "spc2-traj.trr",
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.xtc",
    "spc2-traj.gro",
    "spc2-traj.pdb",
    "spc2-traj.g96"
};

//! Names here work for setVelocity module
const char *const setVelocitySupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.gro",
    "spc2-traj.trr",
};

//! Names here don't work for setVelocity module
const char *const setVelocityUnSupported[] = {
    "spc2-traj.xtc",
    "spc2-traj.pdb",
    "spc2-traj.g96"
};

//! Names here work for setForce module
const char *const setForceSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.trr",
};

//! Names here don't work for setForce module
const char *const setForceUnSupported[] = {
    "spc2-traj.xtc",
    "spc2-traj.pdb",
    "spc2-traj.gro",
    "spc2-traj.g96"
};

//! Names here work for setPrecision module
const char *const setPrecisionSupported[] = {
#if GMX_USE_TNG
    "spc2-traj.tng",
#endif
    "spc2-traj.xtc",
};

//! Names here don't work for setPrecision module
const char *const setPrecisionUnSupported[] = {
    "spc2-traj.trr",
    "spc2-traj.pdb",
    "spc2-traj.gro",
    "spc2-traj.g96"
};


} // namespace test

} // namespace gmx

#endif
