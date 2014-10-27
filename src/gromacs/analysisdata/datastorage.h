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
/*! \libinternal \file
 * \brief
 * Declares gmx::AnalysisDataStorage.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATASTORAGE_H
#define GMX_ANALYSISDATA_DATASTORAGE_H

#include <vector>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisDataFrameHeader;
class AnalysisDataFrameRef;
class AnalysisDataModuleManager;
class AnalysisDataParallelOptions;

class AnalysisDataStorage;

namespace internal
{
class AnalysisDataStorageImpl;
class AnalysisDataStorageFrameData;
}   // namespace internal

/*! \libinternal \brief
 * Allows assigning values for a data frame in AnalysisDataStorage.
 *
 * This class implements the necessary methods to add new data into the
 * storage.  AnalysisDataStorage::startFrame() returns an object of this type,
 * which can be used to add one or more point sets to that data frame.
 * When all data has been added, finishFrame() needs to be called.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataStorageFrame
{
    public:
        /*! \brief Frees the frame object.
         *
         * Should not be called outside AnalysisDataStorage.
         */
        ~AnalysisDataStorageFrame();

        /*! \brief
         * Select data set that all other methods operate on.
         *
         * \param[in] index  Zero-based data set index to select.
         *
         * With multipoint data, a single point set can only contain values in
         * a single data set.
         * With non-multipoint data, arbitrary sequences of selectDataSet() and
         * setValue() are supported.  The full frame is notified to the modules
         * once it is finished.
         *
         * Does not throw.
         */
        void selectDataSet(int index);

        //! Returns number of columns for the frame.
        int columnCount() const { return columnCount_; }

        /*! \brief
         * Sets value for a column.
         *
         * \param[in] column  Zero-based column index.
         * \param[in] value   Value to set for the column.
         * \param[in] bPresent Present flag to set for the column.
         *
         * If called multiple times for a column (within one point set for
         * multipoint data), old values are overwritten.
         *
         * Does not throw.
         */
        void setValue(int column, real value, bool bPresent = true)
        {
            GMX_ASSERT(column >= 0 && column < columnCount(),
                       "Invalid column index");
            values_[currentOffset_ + column].setValue(value, bPresent);
            bPointSetInProgress_ = true;
        }
        /*! \brief
         * Sets value for a column.
         *
         * \param[in] column  Zero-based column index.
         * \param[in] value   Value to set for the column.
         * \param[in] error   Error estimate to set for the column.
         * \param[in] bPresent Present flag to set for the column.
         *
         * If called multiple times for a column (within one point set for
         * multipoint data), old values are overwritten.
         *
         * Does not throw.
         */
        void setValue(int column, real value, real error, bool bPresent = true)
        {
            GMX_ASSERT(column >= 0 && column < columnCount(),
                       "Invalid column index");
            values_[currentOffset_ + column].setValue(value, error, bPresent);
            bPointSetInProgress_ = true;
        }
        /*! \brief
         * Access value for a column.
         *
         * \param[in] column  Zero-based column index.
         *
         * Should only be called after the column value has been set using
         * setValue(); assigning a value to \c value(i) does not mark the
         * column as set.
         *
         * Does not throw.
         */
        real &value(int column)
        {
            GMX_ASSERT(column >= 0 && column < columnCount(),
                       "Invalid column index");
            return values_[currentOffset_ + column].value();
        }
        /*! \brief
         * Access value for a column.
         *
         * \param[in] column  Zero-based column index.
         *
         * Should only be called after the column value has been set using
         * setValue().
         *
         * Does not throw.
         */
        real value(int column) const
        {
            GMX_ASSERT(column >= 0 && column < columnCount(),
                       "Invalid column index");
            return values_[currentOffset_ + column].value();
        }
        /*! \brief
         * Mark point set as finished for multipoint data.
         *
         * Must be called after each point set for multipoint data, including
         * the last (i.e., no values must be set between the last call to this
         * method and AnalysisDataStorage::finishFrame()).
         * Must not be called for non-multipoint data.
         *
         * After this method has been called, all values appear as not set.
         *
         * May call AnalysisDataModuleManager::notifyPointsAdd() and
         * AnalysisDataModuleManager::notifyParallelPointsAdd(), and may throw
         * any exception these methods throw.
         */
        void finishPointSet();
        /*! \brief
         * Finish storing a frame.
         *
         * Must be called exactly once for each frame returned by startFrame(),
         * after the corresponding call.
         * The frame object must not be accessed after the call.
         *
         * Calls notification methods in AnalysisDataModuleManager, and may
         * throw any exceptions these methods throw.
         */
        void finishFrame();

    private:

        /*! \brief
         * Create a new storage frame.
         *
         * \param[in] data  Data object for which the frame is for
         *      (used for data set and column counts).
         */
        explicit AnalysisDataStorageFrame(const AbstractAnalysisData &data);

        //! Clear all column values from the frame.
        void clearValues();

        //! Implementation data.
        internal::AnalysisDataStorageFrameData *data_;
        //! Values for the currently in-progress point set.
        std::vector<AnalysisDataValue>          values_;

        //! Index of the currently active dataset.
        int                                     currentDataSet_;
        //! Offset of the first value in \a values_ for the current data set.
        int                                     currentOffset_;
        //! Number of columns in the current data set.
        int                                     columnCount_;

        //! Whether any values have been set in the current point set.
        bool                                    bPointSetInProgress_;

        //! Needed for access to the constructor.
        friend class internal::AnalysisDataStorageImpl;
        //! Needed for managing the frame the object points to.
        friend class internal::AnalysisDataStorageFrameData;

        GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisDataStorageFrame);
};

/*! \libinternal \brief
 * Helper class that implements storage of data.
 *
 * This class implements a standard way of storing data to avoid implementing
 * storage in each class derived from AbstractAnalysisData separately.
 * To use this class in a class derived from AbstractAnalysisData, a member
 * variable of this type should be declared and the pure virtual methods
 * forwarded to frameCount(), tryGetDataFrame() and requestStorage().
 * Storage properties should be set up, and then startDataStorage() or
 * startParallelDataStorage() called.
 * New frames can then be added using startFrame(), currentFrame(),
 * finishFrame(), and finishFrameSerial() methods (the last is only necessary
 * if startParallelDataStorage() is used).  When all frames are ready,
 * finishDataStorage() must be called.  These methods (and
 * AnalysisDataStorageFrame::finishPointSet()) take the responsibility of
 * calling all the notification methods in AnalysisDataModuleManager,
 *
 * \todo
 * Proper multi-threaded implementation.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataStorage
{
    public:
        //! Constructs a storage object.
        AnalysisDataStorage();
        ~AnalysisDataStorage();

        /*! \brief
         * Returns the number of ready frames.
         *
         * This method is designed such that calls to
         * AbstractAnalysisData::frameCount() can be directly forwarded to this
         * method.  See that method for more documentation.
         *
         * If this method returns N, this means that the first N frames have
         * all been finished.
         *
         * \see AbstractAnalysisData::frameCount()
         */
        int frameCount() const;
        /*! \brief
         * Implements access to data frames.
         *
         * This method is designed such that calls to
         * AbstractAnalysisData::tryGetDataFrameInternal() can be directly
         * forwarded to this method.  See that method for more documentation.
         *
         * A valid reference for a frame will be returned after finishFrame()
         * has been called for that frame.
         *
         * \see AbstractAnalysisData::tryGetDataFrameInternal()
         */
        AnalysisDataFrameRef tryGetDataFrame(int index) const;
        /*! \brief
         * Implements storage requests.
         *
         * This method is designed such that calls to
         * AbstractAnalysisData::requestStorageInternal() can be directly
         * forwarded to this method.  See that method for more documentation.
         *
         * \see AbstractAnalysisData::requestStorageInternal()
         */
        bool requestStorage(int nframes);

        /*! \brief
         * Start storing data.
         *
         * \param[in] data    AbstractAnalysisData object containing this
         *      storage.
         * \param     modules Module manager for \p data.
         * \exception std::bad_alloc if storage allocation fails.
         *
         * Typically called as \c startDataStorage(this, &moduleManager())
         * from a member of \p data when the data is ready to be started.
         * The storage object will take responsibility of calling all
         * module notification methods in AnalysisDataModuleManager using
         * \p modules.
         *
         * Lifetime of \p data and \p modules must exceed the lifetime of the
         * storage object
         * (typically, the storage object will be a member in \p data).
         *
         * Calls AnalysisDataModuleManager::notifyDataStart(), and throws any
         * exceptions this method throws.
         */
        void startDataStorage(AbstractAnalysisData      *data,
                              AnalysisDataModuleManager *modules);
        /*! \brief
         * Start storing data in parallel.
         *
         * \param[in] data    AbstractAnalysisData object containing this
         *      storage.
         * \param[in] options Parallelization options to use.
         * \param     modules Module manager for \p data.
         * \exception std::bad_alloc if storage allocation fails.
         *
         * Should be called instead of startDataStorage() if the data will be
         * produced in parallel.  Works as startDataStorage(), but additionally
         * initializes the storage and the attached modules to prepare for
         * out-of-order data frames.
         *
         * Calls AnalysisDataModuleManager::notifyParallelDataStart(), and
         * throws any exceptions this method throws.
         */
        void startParallelDataStorage(
            AbstractAnalysisData              *data,
            AnalysisDataModuleManager         *modules,
            const AnalysisDataParallelOptions &options);
        /*! \brief
         * Starts storing a new frame.
         *
         * \param[in] header  Header for the new frame.
         * \retval  Frame object corresponding to the started frame.
         * \exception std::bad_alloc if storage reallocation fails
         *      (only possible if storage of all frames has been requested).
         * \exception APIError if frame is too far in the future.
         *
         * The returned object will be valid until the corresponding
         * finishFrame() call.
         *
         * Must be called exactly once for each frame index.
         *
         * Currently, the implementation only works if the new frame is not too
         * far in the future:
         * If \c i is the index of the last frame such that all frames from
         * 0, ..., \c i have been finished, then \p header().index() should be
         * at most \c parallelizationFactor larger than \c i, where
         * parallelizationFactor is the parallelization factor passed to
         * setParallelOptions().
         * Throws APIError if this constraint is violated.
         *
         * Calls AnalysisDataModuleManager::notifyFrameStart() (in certain
         * cases) and AnalysisDataModuleManager::notifyParallelFrameStart(),
         * and throws any exceptions these methods throw.
         */
        AnalysisDataStorageFrame &startFrame(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Convenience method to start storing a new frame.
         *
         * Identical to \c startFrame(AnalysisDataFrameHeader(index, x, dx));
         */
        AnalysisDataStorageFrame &startFrame(int index, real x, real dx);
        /*! \brief
         * Obtains a frame object for an in-progress frame.
         *
         * \param[in] index  Frame index.
         * \retval  Frame object corresponding to \p index.
         *
         * startFrame() should have been called for the frame with index
         * \p index, and finishFrame() should not yet have been called.
         * Returns the same object as returned by the original startFrame()
         * call for the same index.
         *
         * Does not throw.
         */
        AnalysisDataStorageFrame &currentFrame(int index);
        /*! \brief
         * Convenience method for finishing a data frame.
         *
         * \param[in] index  Frame index.
         *
         * Identical to \c currentFrame(index).finishFrame().
         *
         * \see AnalysisDataStorageFrame::finishFrame()
         */
        void finishFrame(int index);
        /*! \brief
         * Performs in-order sequential processing for a data frame.
         *
         * \param[in] index  Frame index.
         *
         * If startParallelDataStorage() has been called with options that
         * indicate parallelism, this method must be called after
         * `finishFrame(index)` (or the equivalent call in
         * AnalysisDataStorageFrame), such that it is called in the correct
         * order sequentially for each frame.
         *
         * If there is no parallelism, this method does nothing; the equivalent
         * processing is done already during finishFrame().
         */
        void finishFrameSerial(int index);
        /*! \brief
         * Finishes storing data.
         *
         * Calls AnalysisDataModuleManager::notifyDataFinish(), and throws any
         * exceptions this method throws.
         */
        void finishDataStorage();

    private:
        typedef internal::AnalysisDataStorageImpl Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
