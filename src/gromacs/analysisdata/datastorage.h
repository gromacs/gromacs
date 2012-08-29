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
 * Declares gmx::AnalysisDataStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATASTORAGE_H
#define GMX_ANALYSISDATA_DATASTORAGE_H

#include <vector>

#include "../legacyheaders/types/simple.h"

#include "../utility/common.h"
#include "../utility/gmxassert.h"

#include "external/mpp/gmxmpp.h"

#include "dataframe.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisDataFrameHeader;
class AnalysisDataFrameRef;
class AnalysisDataParallelOptions;

class AnalysisDataStorage;

/*! \libinternal \brief
 * Stores a single data frame for AnalysisDataStorage.
 *
 * It also provides access to the frame outside AnalysisDataStorage.
 *
 * It is implemented such that the frame header is always valid, i.e.,
 * header().isValid() returns always true.
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

        //! Returns header for the frame.
        const AnalysisDataFrameHeader &header() const { return header_; }
        //! Returns zero-based index of the frame.
        int frameIndex() const { return header().index(); }
        //! Returns x coordinate for the frame.
        real x() const { return header().x(); }
        //! Returns error in x coordinate for the frame if applicable.
        real dx() const { return header().dx(); }
        //! Returns number of columns for the frame.
        int columnCount() const { return values_.size(); }

        /*! \brief
         * Returns point set reference to currently set values.
         *
         * Does not throw.
         */
        AnalysisDataPointSetRef currentPoints() const;
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
            values_[column].setValue(value, bPresent);
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
            values_[column].setValue(value, error, bPresent);
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
            return values_[column].value();
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
            return values_[column].value();
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
         * Calls AbstractAnalysisData::notifyPointsAdd(), and throws any
         * exception this method throws.
         */
        void finishPointSet();

    private:
        /*! \brief
         * Create a new storage frame.
         *
         * \param     storage      Storage object this frame belongs to.
         * \param[in] columnCount  Number of columns for the frame.
         * \param[in] index        Zero-based index for the frame.
         */
        AnalysisDataStorageFrame(AnalysisDataStorage *storage, int columnCount,
                                 int index);

        //! Clear all column values from the frame.
        void clearValues();

        //! Storage object that contains this frame.
        AnalysisDataStorage    &storage_;
        //! Header for the frame.
        AnalysisDataFrameHeader header_;
        //! Values for the frame.
        std::vector<AnalysisDataValue> values_;

        /*! \brief
         * Needed for full write access to the data and for access to
         * constructor/destructor.
         */
        friend class AnalysisDataStorage;
        friend struct mpi::mpi_type_traits<AnalysisDataStorageFrame>;
        GMX_DISALLOW_COPY_AND_ASSIGN(AnalysisDataStorageFrame);
};

/*! \libinternal \brief
 * Helper class that implements storage of data.
 *
 * This class implements a standard way of storing data to avoid implementing
 * storage in each class derived from AbstractAnalysisData separately.
 * To use this class in a class derived from AbstractAnalysisData, a member
 * variable of this type should be declared and the data storage methods
 * forwarded to tryGetDataFrame() and requestStorage() in that object.
 * Storage properties should be set up, and then startDataStorage() called
 * after calling AbstractAnalysisData::notifyDataStart().
 * New frames can then be added using startFrame(), currentFrame() and
 * finishFrame() methods.  These methods (and
 * AnalysisDataStorageFrame::finishPointSet()) take the responsibility of
 * calling AbstractAnalysisData::notifyFrameStart(),
 * AbstractAnalysisData::notifyPointsAdd() and
 * AbstractAnalysisData::notifyFrameFinish() appropriately.
 *
 * \todo
 * Full support for multipoint data.
 * Currently, multipoint data is only supported in serial pass-through mode
 * without any storage.
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
         * Sets whether the data will be multipoint.
         *
         * \exception APIError if storage has been requested.
         *
         * It is not mandatory to call this method (the multipoint property is
         * automatically detected in startDataStorage()), but currently calling
         * it will produce better diagnostic messages.
         * When full support for multipoint data has been implemented, this
         * method can be removed.
         */
        void setMultipoint(bool bMultipoint);
        /*! \brief
         * Set parallelization options for the storage.
         *
         * \param[in] opt  Parallization options to use.
         *
         * If this method is not called, the storage is set up for serial
         * storage only.
         *
         * Does not throw.
         */
        void setParallelOptions(const AnalysisDataParallelOptions &opt);

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
         * Currently, multipoint data cannot be stored, but all other storage
         * request will be honored.
         *
         * \see AbstractAnalysisData::requestStorageInternal()
         */
        bool requestStorage(int nframes);

        /*! \brief
         * Start storing data.
         *
         * \param  data  AbstractAnalysisData object containing this storage.
         * \exception std::bad_alloc if storage allocation fails.
         * \exception APIError if storage has been requested and \p data is
         *      multipoint.
         *
         * Lifetime of \p data must exceed the lifetime of the storage object
         * (typically, the storage object will be a member in \p data).
         * The storage object will take responsibility of calling
         * AbstractAnalysisData::notifyFrameStart(),
         * AbstractAnalysisData::notifyPointsAdd() and
         * AbstractAnalysisData::notifyFrameFinish() for \p data appropriately.
         *
         * AbstractAnalysisData::notifyDataStart() must have been called for
         * \p data, because that may trigger storage requests from attached
         * modules.
         */
        void startDataStorage(AbstractAnalysisData *data);
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
         * at most \c 2*parallelizationFactor-1 larger than \c i, where
         * parallelizationFactor is the parallelization factor passed to
         * setParallelOptions().
         * Throws APIError if this constraint is violated.
         *
         * Calls AbstractAnalysisData::notifyDataStarted() in certain cases,
         * and throws any exceptions this method throws.
         */
        AnalysisDataStorageFrame &startFrame(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Convenience method to start storing a new frame.
         *
         * Identical to \c startFrame(AnalysisDataFrameHeader(index, x, dx));
         */
        AnalysisDataStorageFrame &startFrame(int index, real x, real dx);
        /*! \brief
         * Obtain a frame object for an in-progress frame.
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
         * Finish storing a frame.
         *
         * \param[in] index  Frame index.
         *
         * Must be called exactly once for each startFrame(), after the
         * corresponding call.
         *
         * Calls notification methods in AbstractAnalysisData, and throws any
         * exceptions these methods throw.
         */
        void finishFrame(int index);
        /*! \brief
         * Convenience method for finishing a data frame.
         *
         * \param[in] frame  Frame object to finish.
         *
         * \p frame should have been previously obtained from startFrame() or
         * currentFrame().  The \p frame reference is no longer valid after the
         * call.
         *
         * Identical to \c finishframe(frame.frameIndex());
         */
        void finishFrame(const AnalysisDataStorageFrame &frame);

        //! Collects all stored frames onto rank 0
        void collectFrames (void);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        /*! \brief
         * Needed because the frame object needs to trigger notifications.
         */
        friend void AnalysisDataStorageFrame::finishPointSet();
        template <class T>
            friend MPI_Datatype mpi::mpi_type_traits<T>::get_type(const T&);
        template <class T>
            friend const size_t mpi::mpi_type_traits<T>::get_size(const T&);
};

} // namespace gmx

namespace mpi
{
template <>
inline MPI_Datatype mpi_type_traits<gmx::AnalysisDataStorageFrame>::get_type(const gmx::AnalysisDataStorageFrame& adsf)
{
	mpi_type_builder builder(adsf);
	builder.add(adsf.header_);
	builder.add(adsf.values_);
	return builder.build();
}
}//mpi

#endif
