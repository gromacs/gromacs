/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \libinternal \file
 * \brief
 * Helper classes for testing classes that derive from AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_TESTS_DATATEST_H
#define GMX_ANALYSISDATA_TESTS_DATATEST_H

#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/types/simple.h"

#include "testutils/refdata.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;

namespace test
{

class MockAnalysisModule;

//! Constant to use to signify end of data for AnalysisDataTestInput.
const real END_OF_DATA = std::numeric_limits<real>::max();
//! Constant to use to signify end of one data frame for AnalysisDataTestInput.
const real END_OF_FRAME = std::numeric_limits<real>::min();

/*! \libinternal \brief
 * Represents a single frame in AnalysisDataTestInput structure.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInputFrame
{
    public:
        AnalysisDataTestInputFrame();

        real x() const { return x_; }
        real dx() const { return 0.0; }
        real y(int i) const { return y_[i]; }
        const std::vector<real> &yvector() const { return y_; }
        const real *yptr() const { return &y_[0]; }
        const real *dyptr() const { return NULL; }
        const bool *presentptr() const { return NULL; }

    private:
        real                    x_;
        std::vector<real>       y_;

        friend class AnalysisDataTestInput;
};

/*! \libinternal \brief
 * Represents static input data for AbstractAnalysisData tests.
 *
 * Used to construct structured test input data from a static array of reals,
 * and then typically used as input to methods in AnalysisDataTestFixture.
 *
 * \see AnalysisDataTestFixture
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInput
{
    public:
        /*! \brief
         * Constructs data representation from a simple array.
         *
         * \param[in] data  Array to construct data from.
         *
         * The input array should consist of a set of frames, separated by a
         * END_OF_FRAME marker.  The first value for a frame is the X value,
         * all following values are Y values.
         * The final element in the array (after the last END_OF_FRAME) should
         * be END_OF_DATA.
         */
        explicit AnalysisDataTestInput(const real *data);
        ~AnalysisDataTestInput();

        int frameCount() const { return frames_.size(); }
        int columnCount() const { return columnCount_; }
        const AnalysisDataTestInputFrame &frame(int index) const;

    private:
        int                     columnCount_;
        std::vector<AnalysisDataTestInputFrame> frames_;
};

/*! \libinternal \brief
 * Test fixture for AbstractAnalysisData testing.
 *
 * This test fixture is designed to help writing tests for objects that
 * derive from the AbstractAnalysisData class.  Typical flow in such tests is
 * that first the test creates the required data objects, then call static
 * methods in this class to add mock modules (using
 * AbstractAnalysisData::addModule()) to the data objects to check that they
 * produce the correct data, and then invokes methods in the data object to
 * produce the data to be checked.  Static methods are also provided for
 * pushing data from an AnalysisDataTestInput object to some generic types
 * derived from AbstractAnalysisData.
 *
 * Methods addStaticCheckerModule(), addStaticColumnCheckerModule() and
 * addStaticStorageCheckerModule() create and add mock modules that check the
 * data against a given AnalysisDataTestInput instance.
 *
 * Method addReferenceCheckerModule() creates and adds a mock module that
 * checks the output against reference data produced by a previous test
 * execution (see TestReferenceData).  Two versions are provided, a static
 * method to be used with any TestReferenceChecker, and a non-static method
 * that uses the protected \p data_ member.
 *
 * presentAllData() and presentDataFrame() are provided to push data from an
 * AnalysisDataTestInput into an AnalysisData object.  In typical tests, most
 * checks are done during the these methods, by the added mock modules.
 *
 * \todo
 * Support for errors and for multipoint data.
 *
 * \see AnalysisDataTestInput
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestFixture : public ::testing::Test
{
    public:
        AnalysisDataTestFixture();

        /*! \brief
         * Adds all data from AnalysisDataTestInput into an AnalysisData.
         */
        static void presentAllData(const AnalysisDataTestInput &input,
                                   AnalysisData *data);
        /*! \brief
         * Adds a single frame from AnalysisDataTestInput into an AnalysisData.
         */
        static void presentDataFrame(const AnalysisDataTestInput &input, int row,
                                     AnalysisDataHandle *handle);

        /*! \brief
         * Adds a mock module that verifies output against
         * AnalysisDataTestInput.
         *
         * \param[in]  data     Data to compare against.
         * \param      source   Data object to verify.
         *
         * Creates a mock module that verifies that the
         * AnalysisDataModuleInterface methods are called correctly by
         * \p source.  Parameters for the calls are verified against \p data.
         * Adds the created module to \p source using \p data->addModule().
         * Any exceptions from the called functions should be caught by the
         * caller.
         *
         * \see AbstractAnalysisData::addModule()
         */
        static void addStaticCheckerModule(const AnalysisDataTestInput &data,
                                           AbstractAnalysisData *source);
        /*! \brief
         * Adds a column mock module that verifies output against
         * AnalysisDataTestInput.
         *
         * \param[in]  data     Data to compare against.
         * \param[in]  firstcol First column to check.
         * \param[in]  n        Number of columns to check.
         * \param      source   Data object to verify.
         *
         * Creates a mock module that verifies that the
         * AnalysisDataModuleInterface methods are called correctly by
         * \p source.  Parameters for the calls are verified against \p data.
         * Adds the created module to \p source using
         * \p data->addColumnModule().
         * Any exceptions from the called functions should be caught by the
         * caller.
         *
         * \see AbstractAnalysisData::addColumnModule()
         */
        static void addStaticColumnCheckerModule(const AnalysisDataTestInput &data,
                                                 int firstcol, int n,
                                                 AbstractAnalysisData *source);
        /*! \brief
         * Adds a mock module that verifies output and storage against
         * AnalysisDataTestInput.
         *
         * \param[in]  data     Data to compare against.
         * \param[in]  storageCount  Number of previous frames to check
         *                      (-1 = all).
         * \param      source   Data object to verify.
         *
         * Works like addStaticCheckerModule(), except that in addition, for
         * each frame, the mock module also checks that previous frames can be
         * accessed using AbstractAnalysisData::getDataWErr().  In the
         * AnalysisDataModuleInterface::dataStarted() callback, the mock module
         * calls AbstractAnalysisData::requestStorage() with \p storageCount as
         * the parameter.
         */
        static void addStaticStorageCheckerModule(const AnalysisDataTestInput &data,
                                                  int storageCount,
                                                  AbstractAnalysisData *source);
        /*! \brief
         * Adds a mock module that verifies output against reference data.
         *
         * \param[in]  checker  Reference data checker to use for comparison.
         * \param      source   Data object to verify.
         *
         * Creates a mock module that verifies that the
         * AnalysisDataModuleInterface methods are called correctly by
         * \p source.  Parameters for the calls are verified against reference
         * data using \p checker.
         * Adds the created module to \p source using \p data->addModule().
         * Any exceptions from the called functions should be caught by the
         * caller.
         *
         * \see TestReferenceData
         */
        static void addReferenceCheckerModule(const TestReferenceChecker &checker,
                                              AbstractAnalysisData *source);

        /*! \brief
         * Adds a mock module that verifies output against reference data.
         *
         * \param[in]  id       Identifier for reference data compound to use.
         * \param      source   Data object to verify.
         *
         * Creates a reference checker module using a compound checker with id
         * \p id at the root level of \p data_.
         *
         * See the static overload for other details.
         */
        void addReferenceCheckerModule(const char *id,
                                       AbstractAnalysisData *source);

    protected:
        gmx::test::TestReferenceData  data_;
};

} // namespace test

} // namespace gmx

#endif
