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
/*!\file
 * \libinternal
 * \brief
 * Helper classes for coordinatefile and coordinateio tests
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_COORDINATEIO_H
#define GMX_COORDINATEIO_TESTS_COORDINATEIO_H

#include <utility>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/ioutputadapter.h"
#include "gromacs/coordinateio/requirements.h"
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
 * Create minimal TrajectoryFrameWriter using the provided builder.
 *
 * \param[in] filename      Name of file to create object for.
 * \param[in] topology      Reference to input top.
 * \param[in] selection     Reference to input selection.
 * \param[in] requirements  Requirements for constructing OutputManagar.
 * \throws InconsistentInputError When builder can not create the CoordinateFile.
 * \returns unique_ptr to new CoordinateFile object.
 */
inline TrajectoryFrameWriterPointer createMinimalTrajectoryFrameWriter(const std::string& filename,
                                                                       const TopologyInformation& topology,
                                                                       const Selection& selection,
                                                                       OutputRequirements requirements)
{
    return createTrajectoryFrameWriter(topology.mtop(),
                                       selection,
                                       filename,
                                       topology.hasTopology() ? topology.copyAtoms() : nullptr,
                                       requirements);
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
        top_.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.pdb").string());
        sc_.setTopology(top_.mtop(), 0);
    }

    /*! \brief
     * Method to add a valid selection option to the Module, or an invalid
     * one for testing.
     *
     * \param[in] sel Selection to add option to.
     * \param[in] useValid Set if the added selection should be valid for the module.
     */
    void addOptionForSelection(Selection* sel, bool useValid);

    /*! \brief
     * Set the actual values for the selection.
     *
     * \param[in] options Option to set values for.
     * \param[in] sel Selection to use.
     * \param[in] useValid Set if the added selection should be valid for the module.
     */
    void setSelectionOptionValues(Options* options, Selection* sel, bool useValid);

    /*! \brief
     * Get pointer to options to set values.
     *
     * \returns Pointer to options.
     */
    Options* getOption() { return &options_; }

private:
    //! Selection collection used for handling test selection.
    SelectionCollection sc_;
    //! Selection manager for test selection.
    SelectionOptionManager manager_;
    //! Options manager for test selection input.
    Options options_;
    //! Topology information needed for test selection atoms.
    TopologyInformation top_;
};

inline void ModuleSelection::addOptionForSelection(Selection* sel, bool useValid)
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

inline void ModuleSelection::setSelectionOptionValues(Options* options, Selection* sel, bool useValid)
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
class ModuleTest : public gmx::test::CommandLineTestBase, public ::testing::WithParamInterface<const char*>
{
public:
    /*! \brief
     * Run the builder to create an TrajectoryFrameWriter during tests.
     *
     * \param[in] filename Name for output file, to determine filetype during construction.
     * \param[in] requirements Requirements for adding to the object.
     * \returns The newly created object.
     */
    TrajectoryFrameWriterPointer runTest(const char* filename, const OutputRequirements& requirements) const
    {
        return createMinimalTrajectoryFrameWriter(filename, dummyTopology_, dummySelection_, requirements);
    }
    //! Add topology information to test if needed.
    void addTopology()
    {
        dummyTopology_.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.pdb").string());
    }
    //! Dummy topology to use to create CoordinateFile.
    TopologyInformation dummyTopology_;
    //! Dummy selection.
    Selection dummySelection_;
};

} // namespace test

} // namespace gmx

#endif
