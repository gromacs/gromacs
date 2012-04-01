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
/*! \internal \file
 * \brief
 * Declares internal implementation classes for gmx::AnalysisDataStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATASTORAGE_IMPL_H
#define GMX_ANALYSISDATA_DATASTORAGE_IMPL_H

#include <limits>
#include <vector>

#include "gromacs/utility/uniqueptr.h"

#include "datastorage.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisDataPointSetRef;
class AnalysisDataStorageFrame;

/*! \internal \brief
 * Private implementation class for AnalysisDataStorage.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataStorage::Impl
{
    public:
        //! Smart pointer type for managing a stored frame.
        typedef gmx_unique_ptr<AnalysisDataStorageFrame>::type FramePointer;

        /*! \brief
         * Stored information about a single stored frame.
         *
         * Methods in this class do not throw.
         */
        struct StoredFrame
        {
            //! Indicates what operations have been performed on a frame.
            enum Status
            {
                eMissing,  //!< Frame has not yet been started.
                eStarted,  //!< startFrame() has been called.
                eFinished, //!< finishFrame() has been called.
                eNotified  //!< Appropriate notifications have been sent.
            };

            //! Constructs an object that manages a given frame object.
            explicit StoredFrame(AnalysisDataStorageFrame *frame)
                : frame(frame), status(eMissing)
            {
            }
            //! Whether the frame has been started with startFrame().
            bool isStarted() const { return status >= eStarted; }
            //! Whether the frame has been finished with finishFrame().
            bool isFinished() const { return status >= eFinished; }
            //! Whether all notifications have been sent.
            bool isNotified() const { return status >= eNotified; }
            //! Whether the frame is ready to be available outside the storage.
            bool isAvailable() const { return status >= eFinished; }

            /*! \brief
             * Actual frame data.
             *
             * Never NULL.
             */
            FramePointer              frame;
            //! In what state the frame currently is.
            Status                    status;
        };

        //! Shorthand for a list of data frames that are currently stored.
        typedef std::vector<StoredFrame> FrameList;

        Impl();
        ~Impl();

        //! Returns the number of columns in the attached data.
        int columnCount() const;
        //! Returns whether the storage is set to use multipoint data.
        bool isMultipoint() const;
        /*! \brief
         * Whether storage of all frames has been requested.
         *
         * Storage of all frames also works as expected if \a storageLimit_ is
         * used in comparisons directly, but this method should be used to
         * check how to manage \a frames_.
         */
        bool storeAll() const
        {
            return storageLimit_ == std::numeric_limits<int>::max();
        }
        //! Returns the index of the oldest frame that may be currently stored.
        int firstStoredIndex() const;
        /*! \brief
         * Computes index into \a frames_ for accessing frame \p index.
         *
         * \param[in]  index  Zero-based frame index.
         * \retval  -1 if \p index is not available in \a frames_.
         *
         * Does not throw.
         */
        int computeStorageLocation(int index) const;

        /*! \brief
         * Computes an index into \a frames_ that is one past the last frame
         * stored.
         *
         * Does not throw.
         */
        size_t endStorageLocation() const;

        /*! \brief
         * Extends \a frames_ to a new size.
         *
         * \throws std::bad_alloc if out of memory.
         */
        void extendBuffer(AnalysisDataStorage *storage, size_t newSize);
        /*! \brief
         * Remove oldest frame from the storage to make space for a new one.
         *
         * Increments \a firstFrameLocation_ and reinitializes the frame that
         * was made unavailable by this operation.
         *
         * Does not throw.
         *
         * \see frames_
         */
        void rotateBuffer();

        /*! \brief
         * Calls notification method in \a data_.
         *
         * \throws    unspecified  Any exception thrown by
         *      AbstractAnalysisData::notifyPointsAdd().
         */
        void notifyPointSet(const AnalysisDataPointSetRef &points);
        /*! \brief
         * Calls notification methods for new frames.
         *
         * \param[in] firstLocation  First frame to consider.
         * \throws    unspecified  Any exception thrown by frame notification
         *      methods in AbstractAnalysisData.
         *
         * Notifies \a data_ of new frames (from \p firstLocation and after
         * that) if all previous frames have already been notified.
         * Also rotates the \a frames_ buffer as necessary.
         */
        void notifyNextFrames(size_t firstLocation);

        //! Data object to use for notification calls.
        AbstractAnalysisData   *data_;
        /*! \brief
         * Whether the storage has been set to allow multipoint.
         *
         * Should be possible to remove once full support for multipoint data
         * has been implemented;  isMultipoint() can simply return
         * \c data_->isMultipoint() in that case.
         */
        bool                    bMultipoint_;
        /*! \brief
         * Number of past frames that need to be stored.
         *
         * Always non-negative.  If storage of all frames has been requested,
         * this is set to a large number.
         */
        int                     storageLimit_;
        /*! \brief
         * Number of future frames that may need to be started.
         *
         * Should always be at least one.
         *
         * \see AnalysisDataStorage::startFrame()
         */
        int                     pendingLimit_;
        /*! \brief
         * Data frames that are currently stored.
         *
         * If storage of all frames has been requested, this is simply a vector
         * of frames up to the latest frame that has been started.
         * In this case, \a firstFrameLocation_ is always zero.
         *
         * If storage of all frames is not requested, this is a ring buffer of
         * frames of size \c n=storageLimit_+pendingLimit_+1.  If a frame with
         * index \c index is currently stored, its location is
         * \c index%frames_.size().
         * When at most \a storageLimit_ first frames have been finished,
         * this contains storage for the first \c n-1 frames.
         * When more than \a storageLimit_ first frames have been finished,
         * the oldest stored frame is stored in the location
         * \a firstFrameLocation_, and \a storageLimit_ frames starting from
         * this location are the last finished frames.  \a pendingLimit_ frames
         * follow, and some of these may be in progress or finished.
         * There is always one unused frame in the buffer, which is initialized
         * such that when \a firstFrameLocation_ is incremented, it becomes
         * valid.  This makes it easier to rotate the buffer in concurrent
         * access scenarions (which are not yet otherwise implemented).
         */
        FrameList               frames_;
        //! Location of oldest frame in \a frames_.
        size_t                  firstFrameLocation_;
        /*! \brief
         * Index of next frame that will be added to \a frames_.
         *
         * If all frames are not stored, this will be the index of the unused
         * frame (see \a frames_).
         */
        int                     nextIndex_;
};

} // namespace gmx

#endif
