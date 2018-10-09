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
 * Helper classes for outputmanager and coordinateio tests
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_COORDINATEIO_H
#define GMX_COORDINATEIO_TESTS_COORDINATEIO_H

#include <utility>

#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/coordinateio/builder.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/outputmanager.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Create minimal OutputManager using the provided builder.
 *
 * \param[in] filename      Name of file to create OutputManager for.
 * \param[in] topology      Reference to input top.
 * \param[in] selection     Reference to input selection.
 * \param[in] requirements  Requirements for constructing OutputManagar.
 * \throws InconsistentInputError When builder can not create the OutputManager.
 * \returns unique_ptr to new OutputManager object.
 */
inline OutputManagerPointer createMinimalOutputManager(const std::string         &filename,
                                                       const TopologyInformation &topology,
                                                       const Selection           &selection,
                                                       OutputRequirements         requirements)
{
    requirements.isValid = true;
    return createOutputManager(topology.hasTopology() ? topology.mtop() : nullptr,
                               selection,
                               filename,
                               topology.copyAtoms(),
                               std::move(requirements));
}
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
            sc_.setTopology(top_.mtop(), 0);
        }

        /*! \brief
         * Method to add a valid selection option to the Module, or an invalid
         * one for testing.
         *
         * \param[in] sel Selection to add option to.
         * \param[in] useValid Set if the added selection should be valid for the module.
         */
        void addOptionForSelection(Selection *sel, bool useValid);

        /*! \brief
         * Set the actual values for the selection.
         *
         * \param[in] options Option to set values for.
         * \param[in] sel Selection to use.
         * \param[in] useValid Set if the added selection should be valid for the module.
         */
        void setSelectionOptionValues(Options *options, Selection *sel, bool useValid);

        /*! \brief
         * Get pointer to options to set values.
         *
         * \returns Pointer to options.
         */
        Options *getOption() { return &options_; }

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
ModuleSelection::addOptionForSelection(Selection *sel, bool useValid)
{
    if (useValid)
    {
        options_.addOption(SelectionOption("sel").store(sel).onlyAtoms());
    }
    else
    {
        options_.addOption(SelectionOption("sel").store(sel).dynamicMask());
    }
}

inline void
ModuleSelection::setSelectionOptionValues(Options *options, Selection *sel, bool useValid)
{
    OptionsAssigner assigner(options);
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
         * \param[in] requirements Requirements for adding to the OutputManager.
         */
        void runTest(const char *filename, OutputRequirements requirements)
        {
            output_ = createMinimalOutputManager(filename,
                                                 dummyTopology_,
                                                 dummySelection_,
                                                 std::move(requirements));
        }
        //! Add topology information to test if needed.
        void addTopology()
        {
            dummyTopology_.fillFromInputFile(
                    TestFileManager::getInputFilePath("lysozyme.gro"));
        }
        //! Pointer to new OutputManager object.
        OutputManagerPointer output_;
        //! Dummy topology to use to create OutputManager.
        TopologyInformation  dummyTopology_;
        //! Dummy selection.
        Selection            dummySelection_;
};

} // namespace test

} // namespace gmx

#endif
