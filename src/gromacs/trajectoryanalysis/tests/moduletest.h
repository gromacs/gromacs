/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares test fixture for TrajectoryAnalysisModule subclasses.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_TESTS_MODULETEST_H
#define GMX_TRAJECTORYANALYSIS_TESTS_MODULETEST_H

#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/utility/common.h"

namespace gmx
{

class TrajectoryAnalysisModule;

namespace test
{

class CommandLine;
class FloatingPointTolerance;

/*! \internal
 * \brief
 * Test fixture for trajectory analysis modules.
 *
 * This class implements common logic for all tests for trajectory analysis
 * modules.  The tests simply need to specify the type of the module to test,
 * input files and any additional options, names of datasets to test (if not
 * all), and possible output files to explicitly test by calling the different
 * methods in this class.  runTest() then runs the specified module with the
 * given options and performs all the requested tests against reference data.
 *
 * The actual module to be tested is constructed in the pure virtual
 * createModule() method, which should be implemented in a subclass.
 * Typically, the TrajectoryAnalysisModuleTestFixture template can be used.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_trajectoryanalysis
 */
class AbstractTrajectoryAnalysisModuleTestFixture : public ::testing::Test
{
    public:
        AbstractTrajectoryAnalysisModuleTestFixture();
        virtual ~AbstractTrajectoryAnalysisModuleTestFixture();

        /*! \brief
         * Sets the topology file to use for the test.
         *
         * \param[in] filename  Name of input topology file.
         *
         * \p filename is interpreted relative to the test input data
         * directory, see getTestDataPath().
         *
         * Must be called at most once.  Either this method or setTrajectory()
         * must be called before runTest().
         */
        void setTopology(const char *filename);
        /*! \brief
         * Sets the trajectory file to use for the test.
         *
         * \param[in] filename  Name of input trajectory file.
         *
         * \p filename is interpreted relative to the test input data
         * directory, see getTestDataPath().
         *
         * Must be called at most once.  Either this method or setTopology()
         * must be called before runTest().
         */
        void setTrajectory(const char *filename);
        /*! \brief
         * Sets an output file to use for the test.
         *
         * \param[in] option    Option to set.
         * \param[in] filename  Name of the output file.
         *
         * This method:
         *  - Sets \p option in the tested module to a temporary file name
         *    constructed from \p filename.
         *  - Makes runTest() to check the contents of the file against
         *    reference data after running the module.
         *  - Marks the temporary file for removal at test teardown.
         *
         * \p filename is given to TestTemporaryFileManager to make a unique
         * filename for the temporary file, but is not otherwise used.
         *
         * Currently, this method should not be called for an XVG file, because
         * the comments in the beginning of the file contain timestamps and
         * other variable information, causing the test to fail.  Best used
         * only for custom data formats.  For numeric data, testing the
         * underlying dataset is typically sufficient.
         */
        void setOutputFile(const char *option, const char *filename);
        /*! \brief
         * Sets an output file parameter needed for the test.
         *
         * \param[in] option    Option to set.
         * \param[in] extension Extension for the file to create.
         *
         * This method:
         *  - Sets \p option in the tested module to a temporary file name.
         *  - Marks the temporary file for removal at test teardown.
         *
         * This method provides the mechanism to set output files that are
         * required to trigger computation of values that are required for
         * the test.  The contents of the output file are not tested.
         *
         * \todo
         * Adding facilities to AnalysisData to check whether there are any
         * output modules attached to the data object (directly or indirectly),
         * marking the mocks as output modules, and using these checks in the
         * tools instead of or in addition to the output file presence would be
         * a superior.
         * Also, the extension should be deducible from the option.
         */
        void setOutputFileNoTest(const char *option, const char *extension);
        /*! \brief
         * Includes only specified dataset for the test.
         *
         * \param[in] name  Name of dataset to include.
         *
         * If this method is not called, all datasets are tested by default.
         * If called once, only the specified dataset is tested.
         * If called more than once, also the additional datasets are tested.
         *
         * \p name should be one of the names registered for the tested module
         * using TrajectoryAnalysisModule::registerBasicDataset() or
         * TrajectoryAnalysisModule::registerAnalysisDataset().
         */
        void includeDataset(const char *name);
        /*! \brief
         * Excludes specified dataset from the test.
         *
         * \param[in] name  Name of dataset to exclude.
         *
         * \p name should be one of the names registered for the tested module
         * using TrajectoryAnalysisModule::registerBasicDataset() or
         * TrajectoryAnalysisModule::registerAnalysisDataset().
         */
        void excludeDataset(const char *name);
        /*! \brief
         * Sets a custom tolerance for checking a dataset.
         *
         * \param[in] name       Name of dataset to set the tolerance for.
         * \param[in] tolerance  Tolerance used when verifying the data.
         *
         * \p name should be one of the names registered for the tested module
         * using TrajectoryAnalysisModule::registerBasicDataset() or
         * TrajectoryAnalysisModule::registerAnalysisDataset().
         */
        void setDatasetTolerance(const char                   *name,
                                 const FloatingPointTolerance &tolerance);

        /*! \brief
         * Runs the analysis module with the given additional options.
         *
         * \param[in] args  Options to provide to the module.
         *
         * \p args should be formatted as command-line options, and contain the
         * name of the module as the first argument (the latter requirement is
         * for clarity only).  They are passed to the module in addition to
         * those specified using other methods in this class.
         *
         * All other methods should be called before calling this method.
         *
         * Exceptions thrown by the module are caught by this method.
         */
        void runTest(const CommandLine &args);

    protected:
        /*! \brief
         * Constructs the analysis module to be tested.
         */
        virtual TrajectoryAnalysisModulePointer createModule() = 0;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \internal \brief
 * Test fixture for a trajectory analysis module.
 *
 * \tparam ModuleInfo  Info object for the analysis module to test.
 *
 * \p ModuleInfo should provide a static
 * `TrajectoryAnalysisModulePointer create()` method.
 *
 * \ingroup module_trajectoryanalysis
 */
template <class ModuleInfo>
class TrajectoryAnalysisModuleTestFixture
    : public AbstractTrajectoryAnalysisModuleTestFixture
{
    protected:
        virtual TrajectoryAnalysisModulePointer createModule()
        {
            return ModuleInfo::create();
        }
};

} // namespace test
} // namespace gmx

#endif
