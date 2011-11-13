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
#include "gromacs/fatalerror/gmxassert.h"

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
//! Constant to use to signify end of one multipoint set for AnalysisDataTestInput.
const real MPSTOP = 0.999*std::numeric_limits<real>::max();

/*! \libinternal \brief
 * Represents a single set of points in AnalysisDataTestInputFrame structure.
 *
 * If the data is not multipoint, each frame contains exactly one set of
 * points.  If there is more than one set of points, each of these sets is
 * passed separately to AnalysisDataHandle::addPoints().
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInputPointSet
{
    public:
        AnalysisDataTestInputPointSet();

        int size() const { return y_.size(); }
        real y(int i) const { return y_[i]; }
        const std::vector<real> &yvector() const { return y_; }
        const real *yptr() const { return &y_[0]; }
        const real *dyptr() const { return NULL; }
        const bool *presentptr() const { return NULL; }

    private:
        std::vector<real>       y_;

        friend class AnalysisDataTestInput;
};

/*! \libinternal \brief
 * Represents a single frame in AnalysisDataTestInput structure.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInputFrame
{
    public:
        AnalysisDataTestInputFrame();

        bool isMultipoint() const { return points_.size() > 1; }

        real x() const { return x_; }
        real dx() const { return 0.0; }

        int pointSetCount() const { return points_.size(); }
        const AnalysisDataTestInputPointSet &points(int index = 0) const
        {
            GMX_ASSERT(index >= 0 && static_cast<size_t>(index) < points_.size(),
                       "Point set index out of range");
            return points_[index];
        }

    private:
        real                    x_;
        std::vector<AnalysisDataTestInputPointSet>  points_;

        friend class AnalysisDataTestInput;
};

/*! \libinternal \brief
 * Represents static input data for AbstractAnalysisData tests.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataTestInput
{
    public:
        /*! \brief
         * Constructs data representation from a simple array.
         */
        explicit AnalysisDataTestInput(const real *data);
        ~AnalysisDataTestInput();

        int frameCount() const { return frames_.size(); }
        int columnCount() const { return columnCount_; }
        bool isMultipoint() const { return bMultipoint_; }
        const AnalysisDataTestInputFrame &frame(int index) const;

    private:
        int                     columnCount_;
        bool                    bMultipoint_;
        std::vector<AnalysisDataTestInputFrame> frames_;
};

/*! \libinternal \brief
 * Test fixture for AbstractAnalysisData testing.
 *
 * \todo
 * Support for errors and for arbitrary multipoint data.
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
         * Initializes an array data object from AnalysisDataTestInput.
         *
         * \tparam ArrayData  Class derived from AbstractAnalysisArrayData.
         *
         * The ArrayData class should expose the setter methods
         * (setColumnCount(), setRowCount(), allocateValues(), setValue())
         * publicly or declare the fixture class as a friend.
         * The X axis in \p data must be configured to match \p input before
         * calling this method.
         */
        template <class ArrayData>
        static void setupArrayData(const AnalysisDataTestInput &input,
                                   ArrayData *data);

        /*! \brief
         * Adds a mock module that verifies output against
         * AnalysisDataTestInput.
         *
         * \param[in]  data     Data to compare against.
         * \param      source   Data object to verify.
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
         * For each frame, the mock module also checks that previous frames can
         * be accessed using AbstractAnalysisData::getDataWErr().  In the
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
         */
        static void addReferenceCheckerModule(const TestReferenceChecker &checker,
                                              AbstractAnalysisData *source);

        /*! \brief
         * Adds a mock module that verifies output against reference data.
         *
         * \param[in]  id       Identifier for reference data compound to use.
         * \param      source   Data object to verify.
         */
        void addReferenceCheckerModule(const char *id,
                                       AbstractAnalysisData *source);

    protected:
        gmx::test::TestReferenceData  data_;
};


template <class ArrayData>
void AnalysisDataTestFixture::setupArrayData(const AnalysisDataTestInput &input,
                                             ArrayData *data)
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
        const AnalysisDataTestInputFrame &frame = input.frame(row);
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
