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
 * Implements classes in datastorage.h and paralleloptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "datastorage.h"

#include <limits>
#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/iterator.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

/********************************************************************
 * AnalysisDataParallelOptions
 */

AnalysisDataParallelOptions::AnalysisDataParallelOptions()
    : parallelizationFactor_(1)
{
}


AnalysisDataParallelOptions::AnalysisDataParallelOptions(int parallelizationFactor)
    : parallelizationFactor_(parallelizationFactor)
{
    GMX_RELEASE_ASSERT(parallelizationFactor >= 1,
                       "Invalid parallelization factor");
}

bool AnalysisDataParallelOptions::useMPI() const {
#ifdef GMX_LIB_MPI
    int initialized;
    MPI_Initialized(&initialized);
    return initialized && mpi::comm::world.size() > 1;
#else
    return false;
#endif
}
/********************************************************************
 * AnalysisDataStorage::Impl
 */

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

AnalysisDataStorage::Impl::Impl()
    : data_(NULL), bMultipoint_(false),
      storageLimit_(0), pendingLimit_(1), firstFrameLocation_(0), nextIndex_(0)
{
}


int
AnalysisDataStorage::Impl::columnCount() const
{
    GMX_ASSERT(data_ != NULL, "columnCount() called too early");
    return data_->columnCount();
}


bool
AnalysisDataStorage::Impl::isMultipoint() const
{
    return bMultipoint_;
}


int
AnalysisDataStorage::Impl::firstStoredIndex() const
{
    return frames_[firstFrameLocation_].frame->frameIndex();
}


int
AnalysisDataStorage::Impl::computeStorageLocation(int index) const
{
    if (index < firstStoredIndex() /*|| index >= nextIndex_*/)
    {
        return -1;
    }
    return index % frames_.size();
}


size_t
AnalysisDataStorage::Impl::endStorageLocation() const
{
    if (storeAll())
    {
        return frames_.size();
    }
    if (frames_[0].frame->frameIndex() == 0 || firstFrameLocation_ == 0)
    {
        return frames_.size() - 1;
    }
    return firstFrameLocation_ - 1;
}


void
AnalysisDataStorage::Impl::extendBuffer(AnalysisDataStorage *storage,
                                        size_t newSize)
{
    frames_.reserve(newSize);
    while (frames_.size() < newSize)
    {
        frames_.push_back(StoredFrame(
            new AnalysisDataStorageFrame(storage, columnCount(), nextIndex_)));
        ++nextIndex_;
    }
    // The unused frame should not be included in the count.
    if (!storeAll())
    {
        --nextIndex_;
    }
}


void
AnalysisDataStorage::Impl::rotateBuffer()
{
    GMX_ASSERT(!storeAll(),
               "No need to rotate internal buffer if everything is stored");
    size_t prevFirst = firstFrameLocation_;
    size_t nextFirst = prevFirst + 1;
    if (nextFirst == frames_.size())
    {
        nextFirst = 0;
    }
    firstFrameLocation_ = nextFirst;
    StoredFrame &prevFrame = frames_[prevFirst];
    prevFrame.status = StoredFrame::eMissing;
    prevFrame.frame->header_ = AnalysisDataFrameHeader(nextIndex_ + 1, 0.0, 0.0);
    prevFrame.frame->clearValues();
    ++nextIndex_;
}


void
AnalysisDataStorage::Impl::notifyPointSet(const AnalysisDataPointSetRef &points)
{
    data_->notifyPointsAdd(points);
}


void
AnalysisDataStorage::Impl::notifyNextFrames(size_t firstLocation)
{
    if (!mpi::isMaster())
    {
        return;
    }
    if (firstLocation != firstFrameLocation_)
    {
        // firstLocation can only be zero here if !storeAll() because
        // firstFrameLocation_ is always zero for storeAll()
        int prevIndex =
            (firstLocation == 0 ? frames_.size() - 1 : firstLocation - 1);
        if (!frames_[prevIndex].isNotified())
        {
            return;
        }
    }
    size_t i = firstLocation;
    size_t end = endStorageLocation();
    while (i != end)
    {
        Impl::StoredFrame &storedFrame = frames_[i];
        if (!storedFrame.isFinished())
        {
            break;
        }
        if (storedFrame.status == StoredFrame::eFinished)
        {
            data_->notifyFrameStart(storedFrame.frame->header());
            data_->notifyPointsAdd(storedFrame.frame->currentPoints());
            data_->notifyFrameFinish(storedFrame.frame->header());
            storedFrame.status = StoredFrame::eNotified;
            if (storedFrame.frame->frameIndex() >= storageLimit_)
            {
                rotateBuffer();
            }
        }
        ++i;
        if (!storeAll() && i >= frames_.size())
        {
            i = 0;
        }
    }
}


/********************************************************************
 * AnalysisDataStorageFrame
 */

AnalysisDataStorageFrame::AnalysisDataStorageFrame(AnalysisDataStorage *storage,
                                                   int columnCount, int index)
    : storage_(*storage), header_(index, 0.0, 0.0), values_(columnCount)
{
}


AnalysisDataStorageFrame::~AnalysisDataStorageFrame()
{
}


AnalysisDataPointSetRef
AnalysisDataStorageFrame::currentPoints() const
{
    std::vector<AnalysisDataValue>::const_iterator begin = values_.begin();
    std::vector<AnalysisDataValue>::const_iterator end = values_.end();
    while (begin != end && !begin->isSet())
    {
        ++begin;
    }
    while (end != begin && !(end-1)->isSet())
    {
        --end;
    }
    int firstColumn = (begin != end) ? begin - values_.begin() : 0;
    return AnalysisDataPointSetRef(header_, firstColumn,
                AnalysisDataValuesRef(begin, end));
}


void
AnalysisDataStorageFrame::clearValues()
{
    std::vector<AnalysisDataValue>::iterator i;
    for (i = values_.begin(); i != values_.end(); ++i)
    {
        i->clear();
    }
}


void
AnalysisDataStorageFrame::finishPointSet()
{
    storage_.impl_->notifyPointSet(currentPoints());
    clearValues();
}


/********************************************************************
 * AnalysisDataStorage
 */

AnalysisDataStorage::AnalysisDataStorage()
    : impl_(new Impl())
{
}


AnalysisDataStorage::~AnalysisDataStorage()
{
}


void
AnalysisDataStorage::setMultipoint(bool bMultipoint)
{
    if (bMultipoint && impl_->storageLimit_ > 0)
    {
        GMX_THROW(APIError("Storage of multipoint data not supported"));
    }
    impl_->bMultipoint_ = bMultipoint;
}


void
AnalysisDataStorage::setParallelOptions(const AnalysisDataParallelOptions &opt)
{
    if (opt.useMPI())
    {
        requestStorage(-1);
    }
    impl_->pendingLimit_ = 2 * opt.parallelizationFactor() - 1;
}


AnalysisDataFrameRef
AnalysisDataStorage::tryGetDataFrame(int index) const
{
    if (impl_->isMultipoint())
    {
        return AnalysisDataFrameRef();
    }
    int storageIndex = impl_->computeStorageLocation(index);
    if (storageIndex == -1)
    {
        return AnalysisDataFrameRef();
    }
    const Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    if (!storedFrame.isAvailable())
    {
        return AnalysisDataFrameRef();
    }
    const Impl::FramePointer &frame = storedFrame.frame;
    return AnalysisDataFrameRef(frame->header(), frame->values_);
}


bool
AnalysisDataStorage::requestStorage(int nframes)
{
    if (impl_->isMultipoint())
    {
        return false;
    }

    // Handle the case when everything needs to be stored.
    if (nframes == -1)
    {
        impl_->storageLimit_ = std::numeric_limits<int>::max();
        return true;
    }
    // Check whether an earlier call has requested more storage.
    if (nframes < impl_->storageLimit_)
    {
        return true;
    }
    impl_->storageLimit_ = nframes;
    return true;
}


void
AnalysisDataStorage::startDataStorage(AbstractAnalysisData *data)
{
    //TODO: Multipoint is always notified because it doesn't support storage yet
    //Doesn't work but at least is ignored if anyhow no writer is attached
    //So select can be tested (without e.g. -on option).
    if (mpi::isMaster() || data->isMultipoint())
    {
        data->notifyDataStart();
    }
    // Data needs to be set before calling extendBuffer()
    impl_->data_ = data;
    setMultipoint(data->isMultipoint());
    if (!impl_->storeAll())
    {
        impl_->extendBuffer(this, impl_->storageLimit_ + impl_->pendingLimit_ + 1);
    }
}


AnalysisDataStorageFrame &
AnalysisDataStorage::startFrame(const AnalysisDataFrameHeader &header)
{
    GMX_ASSERT(header.isValid(), "Invalid header");
    Impl::StoredFrame *storedFrame;
    if (impl_->storeAll())
    {
        size_t size = header.index() + 1;
        if (impl_->frames_.size() < size)
        {
            impl_->extendBuffer(this, size);
        }
        storedFrame = &impl_->frames_[header.index()];
    }
    else
    {
        int storageIndex = impl_->computeStorageLocation(header.index());
        if (storageIndex == -1)
        {
            GMX_THROW(APIError("Out of bounds frame index"));
        }
        storedFrame = &impl_->frames_[storageIndex];
    }
    //GMX_RELEASE_ASSERT(!storedFrame->isStarted(),
    //                   "startFrame() called twice for the same frame");
    //GMX_RELEASE_ASSERT(storedFrame->frame->frameIndex() == header.index(),
    //                   "Inconsistent internal frame indexing");
    storedFrame->status = Impl::StoredFrame::eStarted;
    storedFrame->frame->header_ = header;
    if (impl_->isMultipoint())
    {
        impl_->data_->notifyFrameStart(header);
    }
    return *storedFrame->frame;
}


AnalysisDataStorageFrame &
AnalysisDataStorage::startFrame(int index, real x, real dx)
{
    return startFrame(AnalysisDataFrameHeader(index, x, dx));
}


AnalysisDataStorageFrame &
AnalysisDataStorage::currentFrame(int index)
{
    int storageIndex = impl_->computeStorageLocation(index);
    GMX_RELEASE_ASSERT(storageIndex >= 0, "Out of bounds frame index");
    Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    GMX_RELEASE_ASSERT(storedFrame.isStarted(),
                       "currentFrame() called for frame before startFrame()");
    GMX_RELEASE_ASSERT(!storedFrame.isFinished(),
                       "currentFrame() called for frame after finishFrame()");
    GMX_RELEASE_ASSERT(storedFrame.frame->frameIndex() == index,
                       "Inconsistent internal frame indexing");
    return *storedFrame.frame;
}


void
AnalysisDataStorage::finishFrame(int index)
{
    int storageIndex = impl_->computeStorageLocation(index);
    GMX_RELEASE_ASSERT(storageIndex >= 0, "Out of bounds frame index");
    Impl::StoredFrame &storedFrame = impl_->frames_[storageIndex];
    GMX_RELEASE_ASSERT(storedFrame.isStarted(),
                       "finishFrame() called for frame before startFrame()");
    GMX_RELEASE_ASSERT(!storedFrame.isFinished(),
                       "finishFrame() called twice for the same frame");
    GMX_RELEASE_ASSERT(storedFrame.frame->frameIndex() == index,
                       "Inconsistent internal frame indexing");
    storedFrame.status = Impl::StoredFrame::eFinished;
    if (impl_->isMultipoint())
    {
        // TODO: Check that the last point set has been finished
        impl_->data_->notifyFrameFinish(storedFrame.frame->header());
        if (storedFrame.frame->frameIndex() >= impl_->storageLimit_)
        {
            impl_->rotateBuffer();
        }
    }
    else
    {
        impl_->notifyNextFrames(storageIndex);
    }
}

void
AnalysisDataStorage::finishDataStorage()
{
#ifdef GMX_LIB_MPI //otherwise nothing to do
    using mpi::comm;
    Impl timpl;
    std::vector<int> recv_counts;
    std::vector<int> recv_sel;
    if (mpi::isMaster())
    {
        timpl.data_ = impl_->data_;
        recv_counts.resize(comm::world.size());
    }
    // Determine information needed by master
    std::vector<int> selection;
    int recv_size=0;
    // It would be nice to have gather work with in-place. Then
    // this for loop isn't necessary on the master
    // (master wouldn't need to send to itself). Also we could
    // get rid of the whole timpl.
    for (size_t i=0; i<impl_->frames_.size(); i++)
    {
        GMX_ASSERT (mpi::isMaster() || !impl_->frames_[i].isNotified(),
                    "Tried to collect an already notified frame!");
        if (impl_->frames_[i].isFinished())
        {
            recv_size++;
            selection.push_back(i);
        }
    }
    comm::world.gather(&recv_size, (&recv_size)+1, recv_counts.begin(), recv_counts.end(), 0); //gather number of frames on each node

    if (mpi::isMaster())
    {
        size_t totalFrames=0;
        for (int i=0; i<comm::world.size(); i++)
        {
            totalFrames += recv_counts[i];
        }
        recv_sel.resize(totalFrames);
        timpl.extendBuffer(this, std::max(totalFrames, impl_->frames_.size()));
    }
    comm::world.gatherv(selection.begin(), selection.end(),
            recv_sel.begin(), recv_sel.end(),
            recv_counts.begin(), recv_counts.end(), 0); //gather frame numbers of frames
    //gather frames

    comm::world.gatherv(
            make_permutation_iterator(impl_->frames_.begin(),selection.begin()),
            make_permutation_iterator(impl_->frames_.begin(),selection.end()),
            make_permutation_iterator(timpl.frames_.begin(),recv_sel.begin()),
            make_permutation_iterator(timpl.frames_.begin(),recv_sel.end()),
            recv_counts.begin(), recv_counts.end(), 0);

    if (mpi::isMaster())
    {
        impl_->frames_.swap(timpl.frames_);
    }
    impl_->notifyNextFrames(0);
#endif//GMX_LIB_MPI
    if (mpi::isMaster())
    {
        impl_->data_->notifyDataFinish();
    }
}

void
AnalysisDataStorage::finishFrame(const AnalysisDataStorageFrame &frame)
{
    finishFrame(frame.frameIndex());
}

} // namespace gmx

#ifdef GMX_LIB_MPI
namespace mpi
{
//! mpi_type_traits for the enum Status
template <>
inline data_layout
mpi_type_traits<gmx::AnalysisDataStorage::Impl::StoredFrame::Status>::get_layout(
        const gmx::AnalysisDataStorage::Impl::StoredFrame::Status& s)
{
    return
       get_enum_layout<gmx::AnalysisDataStorage::Impl::StoredFrame::Status>(s);
}
//! declaring Status is a named data type
SET_MPI_STATIC(gmx::AnalysisDataStorage::Impl::StoredFrame::Status);

template <>
inline data_layout
mpi_type_traits<gmx::AnalysisDataStorage::Impl::StoredFrame>::get_layout(
        const gmx::AnalysisDataStorage::Impl::StoredFrame& sf)
{
    return struct_layout_builder(sf,2).add(sf.frame).add(sf.status).build();
}
}//mpi

#endif //GMX_LIB_MPI
