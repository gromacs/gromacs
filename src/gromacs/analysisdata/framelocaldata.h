/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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
 * Defines gmx::AnalysisDataFrameLocalData and supporting types.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_FRAMELOCALDATA_H
#define GMX_ANALYSISDATA_FRAMELOCALDATA_H

#include <algorithm>
#include <numeric>
#include <vector>

#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

//! \addtogroup module_analysisdata
//! \{

/*! \internal
 * \brief
 * Handle to a single data set within frame-local data array.
 *
 * Methods in this class do not throw.
 *
 * \see AnalysisDataFrameLocalData
 */
template<typename ValueType>
class AnalysisDataFrameLocalDataSetHandle
{
    public:
        //! Constructs a handle from an array of values.
        explicit AnalysisDataFrameLocalDataSetHandle(ArrayRef<ValueType> values)
            : values_(values)
        {
        }

        //! Clears all values in the data set.
        void clear()
        {
            std::fill(values_.begin(), values_.end(), ValueType());
        }

        //! Accesses a single value in the data set.
        ValueType &value(int column)
        {
            GMX_ASSERT(column >= 0 && column < static_cast<int>(values_.size()),
                       "Invalid column index");
            return values_[column];
        }

    private:
        ArrayRef<ValueType>  values_;
};

/*! \internal
 * \brief
 * Handle to a single frame data within frame-local data array.
 *
 * Methods in this class do not throw.
 *
 * \see AnalysisDataFrameLocalData
 */
template<typename ValueType>
class AnalysisDataFrameLocalDataHandle
{
    public:
        //! Shorthand for the internal array of values.
        typedef std::vector<ValueType> ValueArray;
        //! Shorthand for a handle to a single data set.
        typedef AnalysisDataFrameLocalDataSetHandle<ValueType> DataSetHandle;

        //! Constructs a handle from specified frame data.
        AnalysisDataFrameLocalDataHandle(const std::vector<int> *dataSetIndices,
                                         ValueArray             *values)
            : dataSetIndices_(dataSetIndices), values_(values)
        {
        }

        //! Returns the number of data sets in the array.
        int dataSetCount() const
        {
            return dataSetIndices_->size() - 1;
        }
        //! Clears all values in the frame.
        void clear()
        {
            std::fill(values_->begin(), values_->end(), ValueType());
        }

        //! Returns a handle for a single data set.
        DataSetHandle dataSet(int dataSet)
        {
            GMX_ASSERT(dataSet >= 0 && dataSet < dataSetCount(),
                       "Invalid data set index");
            const int firstIndex = (*dataSetIndices_)[dataSet];
            const int lastIndex  = (*dataSetIndices_)[dataSet + 1];
            return DataSetHandle(makeArrayRef(*values_).
                                     subArray(firstIndex, lastIndex-firstIndex));
        }
        //! Accesses a single value in the frame.
        ValueType &value(int dataSet, int column)
        {
            GMX_ASSERT(dataSet >= 0 && dataSet < dataSetCount(),
                       "Invalid data set index");
            const int firstIndex = (*dataSetIndices_)[dataSet];
            GMX_ASSERT(column >= 0
                       && column < (*dataSetIndices_)[dataSet+1] - firstIndex,
                       "Invalid column index");
            return (*values_)[firstIndex + column];
        }

    private:
        const std::vector<int> *dataSetIndices_;
        ValueArray             *values_;
};

/*! \internal \brief
 * Container for an array of frame-local values that supports parallel data
 * processing.
 *
 * \tparam ValueType Type of values to store.
 *
 * This class provides a convenient interface to create an array of frame-local
 * data for use in analysis data modules that support parallel processing.
 * The object is initialized by setting the desired dimensionality with
 * setDataSetCount() and setColumnCount(), followed by a call to init(),
 * typically in IAnalysisDataModule::parallelDataStarted(),
 *
 * After initialization, frameData() can be used to access the data for a given
 * frame, independently from other frames.  This works if the assumptions about
 * parallelism hold: if `N` is the parallelization factor given for init() with
 * AnalysisDataParallelOptions::parallelizationFactor(), then frame `i+N` must
 * not be accessed before all processing for frame `i` is finished.
 * Technically, the data for different frames is kept in a ring buffer of size
 * `N`.
 *
 * The data for a frame is not cleared after it is reused for a new frame (but
 * is initially cleared).  This allows using the data for accumulating values
 * over all frames in a lock-free manner.
 *
 * frameDataSet() is provided for convenience when only a single data set
 * needs to be accessed (typically in IAnalysisDataModule::pointsAdded()).
 *
 * Methods in this class do not throw except where indicated.
 *
 * \see AnalysisDataFrameLocalData
 */
template<typename ValueType>
class AnalysisDataFrameLocalData
{
    public:
        //! Shorthand for the internal array of values for a frame.
        typedef std::vector<ValueType> ValueArray;
        //! Shorthand for a handle to a single frame.
        typedef AnalysisDataFrameLocalDataHandle<ValueType> FrameHandle;
        //! Shorthand for a handle to a single data set.
        typedef AnalysisDataFrameLocalDataSetHandle<ValueType> DataSetHandle;

        //! Constructs an empty container with a single data set.
        AnalysisDataFrameLocalData()
        {
            dataSetColumns_.resize(2);
        }

        //! Whether init() has been called.
        bool isInitialized() const { return !values_.empty(); }
        /*! \brief
         * Returns number of independent data frames in this object.
         *
         * This supports looping over all the frame arrays to, e.g., sum them
         * up at the end in accumulation scenarios.
         */
        int frameCount() const { return values_.size(); }

        /*! \brief
         * Sets the number of data sets stored for each frame.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * If not called, there is a single data set in the object.
         * Cannot be called after init().
         */
        void setDataSetCount(int dataSetCount)
        {
            GMX_RELEASE_ASSERT(!isInitialized(),
                               "Cannot change value count after init()");
            GMX_RELEASE_ASSERT(dataSetCount >= 0,
                               "Invalid data set count");
            dataSetColumns_.resize(dataSetCount + 1);
        }
        /*! \brief
         * Sets the number of columns stored for a data set.
         *
         * Must be called for each data set that needs to have values,
         * otherwise there will be zero columns for that data set.
         * Cannot be called after init().
         */
        void setColumnCount(int dataSet, int columnCount)
        {
            GMX_RELEASE_ASSERT(!isInitialized(),
                               "Cannot change value count after init()");
            GMX_RELEASE_ASSERT(dataSet >= 0 && dataSet < static_cast<int>(dataSetColumns_.size()) - 1,
                               "Invalid data set index");
            GMX_RELEASE_ASSERT(columnCount >= 0,
                               "Invalid column count");
            dataSetColumns_[dataSet + 1] = columnCount;
        }

        /*! \brief
         * Initializes the storage to support specified parallelism.
         *
         * \throws std::bad_alloc if out of memory.
         */
        void init(const AnalysisDataParallelOptions &opt)
        {
            GMX_RELEASE_ASSERT(!isInitialized(), "init() called multiple times");
            std::partial_sum(dataSetColumns_.begin(), dataSetColumns_.end(),
                             dataSetColumns_.begin());
            values_.resize(opt.parallelizationFactor());
            typename std::vector<ValueArray>::iterator i;
            for (i = values_.begin(); i != values_.end(); ++i)
            {
                i->resize(dataSetColumns_.back());
            }
        }

        //! Returns a handle to access data for a frame.
        FrameHandle frameData(int frameIndex)
        {
            GMX_ASSERT(frameIndex >= 0, "Invalid frame index");
            GMX_ASSERT(isInitialized(), "Cannot access data before init()");
            return FrameHandle(&dataSetColumns_,
                               &values_[frameIndex % values_.size()]);
        }
        //! Returns a handle to access a single data set within a frame.
        DataSetHandle frameDataSet(int frameIndex, int dataSet)
        {
            return frameData(frameIndex).dataSet(dataSet);
        }

    private:
        /*! \brief
         * Index to find data sets within a per-frame array in `values_`.
         *
         * The first entry is always zero, followed by one entry for each data
         * set.  Before init(), the data set entries hold the numbers set with
         * setColumnCount().  After init(), the data set entries hold the
         * indices of the first column for that data set in the per-frame
         * arrays in `values_`.
         */
        std::vector<int>         dataSetColumns_;
        /*! \brief
         * Data array for each frame.
         *
         * This is a ring buffer whose size is specified by the desired
         * parallelism level.  For each frame, there is a single array of
         * values, where the individual data sets are indexed with
         * `dataSetColumns_`.
         */
        std::vector<ValueArray>  values_;
};

//! \}

} // namespace gmx

#endif
