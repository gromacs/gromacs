/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Helper classes for testing classes that derive from AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_TESTS_DATATEST_H
#define GMX_ANALYSISDATA_TESTS_DATATEST_H

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"

// currently the bug manifests itself only in AbstractAnalysisData testing
#if (defined __ICL && __ICL >= 1400) || (defined __ICC && __ICC >= 1400) || (defined __PATHSCALE__)
#define STATIC_ANON_NAMESPACE_BUG //see #1558 for details
#endif

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;

namespace test
{

class FloatingPointTolerance;

/*! \libinternal \brief
 * Represents a single set of points in AnalysisDataTestInputFrame structure.
 *
 * If the data is not multipoint, each frame contains exactly one set of
 * points.  If there is more than one set of points, each of these sets is
 * passed separately and AnalysisDataHandle::finishPointSet() called in
 * between.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInputPointSet
{
    public:
        //! Returns zero-based index of this point set in its frame.
        int index() const { return index_; }
        //! Returns zero-based index of the data set of this point set.
        int dataSetIndex() const { return dataSetIndex_; }
        //! Returns zero-based index of the first column in this point set.
        int firstColumn() const { return firstColumn_; }
        //! Returns zero-based index of the last column in this point set.
        int lastColumn() const { return firstColumn_ + size() - 1; }
        //! Returns the number of columns in the point set.
        int size() const { return values_.size(); }
        //! Returns the value in column \p i.
        real y(int i) const { return values_[i].y; }
        //! Returns whether the error is present for column \p i.
        bool hasError(int i) const { return values_[i].bError; }
        //! Returns the error in column \p i.
        real error(int i) const { return values_[i].error; }
        //! Returns whether the value in column \p i is present.
        bool present(int /*i*/) const { return true; }
        //! Returns an AnalysisDataValue for column \p i.
        AnalysisDataValue value(int i) const
        {
            AnalysisDataValue result;
            result.setValue(values_[i].y);
            if (values_[i].bError)
            {
                result.setError(values_[i].error);
            }
            return result;
        }

        //! Appends a value to this point set.
        void addValue(real y) { values_.push_back(Value(y)); }
        //! Appends a value with an error estimate to this point set.
        void addValueWithError(real y, real error)
        {
            values_.push_back(Value(y, error));
        }

    private:
        //! Creates an empty point set.
        AnalysisDataTestInputPointSet(int index, int dataSetIndex,
                                      int firstColumn);

        struct Value
        {
            Value() : y(0.0), error(0.0), bError(false) {}
            explicit Value(real y) : y(y), error(0.0), bError(false) {}
            Value(real y, real error) : y(y), error(error), bError(true) {}

            real                y;
            real                error;
            bool                bError;
        };

        int                     index_;
        int                     dataSetIndex_;
        int                     firstColumn_;
        std::vector<Value>      values_;

        //! For constructing new point sets.
        friend class AnalysisDataTestInputFrame;
};

/*! \libinternal \brief
 * Represents a single frame in AnalysisDataTestInput structure.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
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
        int pointSetCount() const { return pointSets_.size(); }
        //! Returns a point set object for a given point set.
        const AnalysisDataTestInputPointSet &pointSet(int index) const
        {
            GMX_ASSERT(index >= 0 && static_cast<size_t>(index) < pointSets_.size(),
                       "Point set index out of range");
            return pointSets_[index];
        }

        //! Appends an empty point set to this frame.
        AnalysisDataTestInputPointSet &addPointSet(int dataSet, int firstColumn);
        //! Adds a point set with given values to this frame.
        void addPointSetWithValues(int dataSet, int firstColumn, real y1);
        //! Adds a point set with given values to this frame.
        void addPointSetWithValues(int dataSet, int firstColumn,
                                   real y1, real y2);
        //! Adds a point set with given values to this frame.
        void addPointSetWithValues(int dataSet, int firstColumn,
                                   real y1, real y2, real y3);
        //! Adds a point set with given values to this frame.
        void addPointSetWithValueAndError(int dataSet, int firstColumn,
                                          real y1, real e1);

    private:
        //! Constructs a new frame object with the given values.
        AnalysisDataTestInputFrame(int index, real x);

        int                                         index_;
        real                                        x_;
        std::vector<AnalysisDataTestInputPointSet>  pointSets_;

        //! For constructing new frames.
        friend class AnalysisDataTestInput;
};

/*! \libinternal \brief
 * Represents static input data for AbstractAnalysisData tests.
 *
 * Used to construct structured test input data for analysis data unit tests.
 * Typically used as input to methods in AnalysisDataTestFixture.
 *
 * \see AnalysisDataTestFixture
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInput
{
    public:
        /*! \brief
         * Constructs empty input data.
         *
         * \param[in] dataSetCount Number of data sets in the data.
         * \param[in] bMultipoint  Whether the data will be multipoint.
         *
         * The column count for each data set must be set with
         * setColumnCount().
         */
        AnalysisDataTestInput(int dataSetCount, bool bMultipoint);
        ~AnalysisDataTestInput();

        //! Whether the input data is multipoint.
        bool isMultipoint() const { return bMultipoint_; }
        //! Returns the number of data sets in the input data.
        int dataSetCount() const { return columnCounts_.size(); }
        //! Returns the number of columns in a given data set.
        int columnCount(int dataSet) const { return columnCounts_[dataSet]; }
        //! Returns the number of frames in the input data.
        int frameCount() const { return frames_.size(); }
        //! Returns a frame object for the given input frame.
        const AnalysisDataTestInputFrame &frame(int index) const;

        //! Sets the number of columns in a data set.
        void setColumnCount(int dataSet, int columnCount);
        //! Appends an empty frame to this data.
        AnalysisDataTestInputFrame &addFrame(real x);
        //! Adds a frame with a single point set and the given values.
        void addFrameWithValues(real x, real y1);
        //! Adds a frame with a single point set and the given values.
        void addFrameWithValues(real x, real y1, real y2);
        //! Adds a frame with a single point set and the given values.
        void addFrameWithValues(real x, real y1, real y2, real y3);
        //! Adds a frame with a single point set and the given values.
        void addFrameWithValueAndError(real x, real y1, real e1);

    private:
        std::vector<int>                        columnCounts_;
        bool                                    bMultipoint_;
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
 * checks are done during these methods, by the added mock modules.
 * setupArrayData() performs the same function for classes derived from
 * AbstractAnalysisArrayData.  In that case, the test should separately ensure
 * that AbstractAnalysisArrayData::valuesReady() gets called.
 *
 * \todo
 * Support for arbitrary AnalysisDataValues (errors and missing values).
 *
 * \see AnalysisDataTestInput
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataTestFixture : public ::testing::Test
{
    public:
        AnalysisDataTestFixture();

        /*! \brief
         * Initializes an AnalysisData object from input data.
         *
         * Sets the column count and other properties based on the input data.
         */
        static void setupDataObject(const AnalysisDataTestInput &input,
                                    AnalysisData                *data);

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
         * Adds a mock module that verifies parallel output against
         * AnalysisDataTestInput.
         *
         * \param[in]  data     Data to compare against.
         * \param      source   Data object to verify.
         *
         * Creates a parallel mock module that verifies that the
         * AnalysisDataModuleInterface methods are called correctly by
         * \p source.  Parameters for the calls are verified against \p data.
         * Adds the created module to \p source using \p data->addModule().
         * Any exceptions from the called functions should be caught by the
         * caller.
         *
         * Differs from addStaticCheckerModule() in that the created mock
         * module reports that it accepts parallel input data, and accepts and
         * verifies notification calls following the parallel pattern.
         *
         * \see AbstractAnalysisData::addModule()
         */
        static void addStaticParallelCheckerModule(
            const AnalysisDataTestInput &data,
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
         * \param[in]  checker   Reference data checker to use for comparison.
         * \param[in]  id        Identifier for reference data compound to use.
         * \param      source    Data object to verify.
         * \param[in]  tolerance Tolerance to use for comparison.
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
        static void addReferenceCheckerModule(TestReferenceChecker          checker,
                                              const char                   *id,
                                              AbstractAnalysisData         *source,
                                              const FloatingPointTolerance &tolerance);

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
        gmx::test::TestReferenceData  data_;
};


template <class ArrayData>
void AnalysisDataTestFixture::setupArrayData(const AnalysisDataTestInput &input,
                                             ArrayData                   *data)
{
    GMX_RELEASE_ASSERT(!input.isMultipoint(),
                       "Array data cannot be initialized from multipoint data");
    GMX_RELEASE_ASSERT(input.dataSetCount() == 1,
                       "Array data cannot be initialized from multiple data sets");
    GMX_RELEASE_ASSERT(data->columnCount() == 0 || data->columnCount() == input.columnCount(0),
                       "Mismatching input and target data");
    GMX_RELEASE_ASSERT(data->rowCount() == 0 || data->rowCount() == input.frameCount(),
                       "Mismatching input and target data");
    data->setColumnCount(input.columnCount(0));
    data->setRowCount(input.frameCount());
    data->allocateValues();
    for (int row = 0; row < input.frameCount(); ++row)
    {
        const AnalysisDataTestInputFrame    &frame = input.frame(row);
        EXPECT_FLOAT_EQ(frame.x(), data->xvalue(row));
        GMX_RELEASE_ASSERT(frame.pointSetCount() == 1,
                           "Multiple point sets not supported by array data");
        const AnalysisDataTestInputPointSet &points = frame.pointSet(0);
        for (int column = 0; column < points.size(); ++column)
        {
            data->value(row, column + points.firstColumn()) = points.value(column);
        }
    }
}

} // namespace test

} // namespace gmx

#endif
