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
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_DATATEST_H
#define GMX_TESTUTILS_DATATEST_H

#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/types/simple.h"

#include "gromacs/utility/gmxassert.h"

#include "testutils/refdata.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;

namespace test
{

//! Constant to use to signify end of one data frame for AnalysisDataTestInput.
const real END_OF_FRAME = std::numeric_limits<real>::max();
//! Constant to use to signify end of one multipoint set for AnalysisDataTestInput.
const real MPSTOP = -std::numeric_limits<real>::max();

/*! \libinternal \brief
 * Represents a single set of points in AnalysisDataTestInputFrame structure.
 *
 * If the data is not multipoint, each frame contains exactly one set of
 * points.  If there is more than one set of points, each of these sets is
 * passed separately and AnalysisDataHandle::finishPointSet() called in
 * between.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class AnalysisDataTestInputPointSet
{
    public:
        //! Returns the number of columns in the point set.
        int size() const { return y_.size(); }
        //! Returns the value in column \p i.
        real y(int i) const { return y_[i]; }
        //! Returns the error in column \p i.
        real dy(int i) const { return 0.0; }
        //! Returns whether the value in column \p i is present.
        real present(int i) const { return true; }
        //! Returns a vector of values for all columns.
        const std::vector<real> &yvector() const { return y_; }

    private:
        //! Creates an empty point set.
        AnalysisDataTestInputPointSet();

        std::vector<real> y_;

        friend class AnalysisDataTestInput;
};

/*! \libinternal \brief
 * Represents a single frame in AnalysisDataTestInput structure.
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class AnalysisDataTestInputFrame
{
    public:
        //! Returns zero-based index for the frame.
        int index() const { return index_; }
        //! Returns x coordinate for the frame.
        real x() const { return x_; }
        //! Returns error in the x coordinate for the frame.
        real dx() const { return 0.0; }

        //! Number of individual point sets in the frame.
        int pointSetCount() const { return points_.size(); }
        //! Returns a point set object for a given point set.
        const AnalysisDataTestInputPointSet &points(int index = 0) const
        {
            GMX_ASSERT(index >= 0 && static_cast<size_t>(index) < points_.size(),
                       "Point set index out of range");
            return points_[index];
        }

    private:
        //! Constructs a new frame object with the given values.
        AnalysisDataTestInputFrame(int index, real x);

        int  index_;
        real x_;
        std::vector<AnalysisDataTestInputPointSet> points_;

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
 * \inlibraryapi
 * \ingroup module_testutils
 */
class AnalysisDataTestInput
{
    public:
        /*! \brief
         * Constructs data representation from a simple array.
         *
         * \tparam    count Number of elements in the array.
         * \param[in] data  Array to construct data from.
         *
         * The input array should consist of a set of frames, separated by a
         * END_OF_FRAME marker.  The first value for a frame is the X value,
         * all following values are Y values.
         * For multipoint data, one frame can contain several point sets,
         * separated by MPSTOP markers.  There should be no MPSTOP marker after
         * the last point set, only an END_OF_FRAME marker.  All point sets are
         * assumed to start from column zero, but the sets may contain
         * different number of columns.  For non-multipoint data, all frames
         * must containt the same number of columns.
         * The final element in the array should be an END_OF_FRAME.
         */
        template <size_t count>
        explicit AnalysisDataTestInput(const real (&data)[count])
            : columnCount_(0), bMultipoint_(false)
        {
            initFromArray(data, count);
        }
        ~AnalysisDataTestInput();

        //! Returns the number of frames in the input data.
        int frameCount() const { return frames_.size(); }
        //! Returns the number of columns in the input data.
        int columnCount() const { return columnCount_; }
        //! Whether the input data is multipoint.
        bool isMultipoint() const { return bMultipoint_; }
        //! Returns a frame object for the given input frame.
        const AnalysisDataTestInputFrame &frame(int index) const;

    private:
        void initFromArray(const real *data, size_t count);

        int  columnCount_;
        bool bMultipoint_;
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
 * setupArrayData() performs the same function for classes derived from
 * AbstractAnalysisArrayData.  In that case, the test should separately ensure
 * that AbstractAnalysisArrayData::valuesReady() gets called.
 *
 * \todo
 * Support for errors and for arbitrary multipoint data.
 *
 * \see AnalysisDataTestInput
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class AnalysisDataTestFixture : public ::testing::Test
{
    public:
        AnalysisDataTestFixture();

        /*! \brief
         * Adds all data from AnalysisDataTestInput into an AnalysisData.
         */
        static void presentAllData(const AnalysisDataTestInput &input,
                                   AnalysisData                *data);
        /*! \brief
         * Adds a single frame from AnalysisDataTestInput into an AnalysisData.
         */
        static void presentDataFrame(const AnalysisDataTestInput &input, int row,
                                     AnalysisDataHandle handle);
        /*! \brief
         * Initializes an array data object from AnalysisDataTestInput.
         *
         * \tparam ArrayData  Class derived from AbstractAnalysisArrayData.
         *
         * The ArrayData class should expose the setter methods
         * (setColumnCount(), setRowCount(), allocateValues(), setValue())
         * publicly or declare the fixture class as a friend.
         * The X axis in \p data must be configured to match \p input before
         * calling this method.
         *
         * Does not call AbstractAnalysisArrayData::valuesReady().
         * The test must ensure that this method gets called, otherwise the
         * mock modules never get called.
         */
        template <class ArrayData>
        static void setupArrayData(const AnalysisDataTestInput &input,
                                   ArrayData                   *data);

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
                                           AbstractAnalysisData        *source);
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
         * accessed using AbstractAnalysisData::getDataFrame().  In the
         * AnalysisDataModuleInterface::dataStarted() callback, the mock module
         * calls AbstractAnalysisData::requestStorage() with \p storageCount as
         * the parameter.
         */
        static void addStaticStorageCheckerModule(const AnalysisDataTestInput &data,
                                                  int                          storageCount,
                                                  AbstractAnalysisData        *source);
        /*! \brief
         * Adds a mock module that verifies output against reference data.
         *
         * \param[in]  checker  Reference data checker to use for comparison.
         * \param[in]  id       Identifier for reference data compound to use.
         * \param      source   Data object to verify.
         *
         * Creates a mock module that verifies that the
         * AnalysisDataModuleInterface methods are called correctly by
         * \p source.  Parameters for the calls are verified against reference
         * data using a child compound \p id of \p checker.
         * Adds the created module to \p source using \p data->addModule().
         * Any exceptions from the called functions should be caught by the
         * caller.
         *
         * \see TestReferenceData
         */
        static void addReferenceCheckerModule(TestReferenceChecker  checker,
                                              const char           *id,
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
        void addReferenceCheckerModule(const char           *id,
                                       AbstractAnalysisData *source);

    protected:
        /*! \brief
         * Reference data object used for the reference checker modules.
         *
         * Tests can use the data object also for their own purposes if needed.
         */
        gmx::test::TestReferenceData data_;
};


template <class ArrayData>
void AnalysisDataTestFixture::setupArrayData(const AnalysisDataTestInput &input,
                                             ArrayData                   *data)
{
    GMX_RELEASE_ASSERT(!input.isMultipoint(),
                       "Array data cannot be initialized from multipoint data");
    GMX_RELEASE_ASSERT(data->columnCount() == 0 || data->columnCount() == input.columnCount(),
                       "Mismatching input and target data");
    GMX_RELEASE_ASSERT(data->rowCount() == 0 || data->rowCount() == input.frameCount(),
                       "Mismatching input and target data");
    data->setColumnCount(input.columnCount());
    data->setRowCount(input.frameCount());
    data->allocateValues();
    for (int row = 0; row < input.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame    &frame = input.frame(row);
        EXPECT_FLOAT_EQ(frame.x(), data->xvalue(row));
        const AnalysisDataTestInputPointSet &points = frame.points();
        for (int column = 0; column < input.columnCount(); ++column)
        {
            data->setValue(row, column, points.y(column));
        }
    }
}

} // namespace test

} // namespace gmx

#endif
