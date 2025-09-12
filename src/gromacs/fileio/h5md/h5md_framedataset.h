/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \brief Declarations of H5md frame data set manipulation routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_FRAMEDATASET_H
#define GMX_FILEIO_H5MD_FRAMEDATASET_H

#include <hdf5.h>

#include <memory>

#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_guard.h"

namespace gmx
{

template<typename ValueType>
class ArrayRef;
template<typename ValueType>
class BasicVector;

/*! \brief Class which provides an interface for reading and writing frame data to H5mdDataSet objects.
 *
 * Can be made either by opening an existing data set via its name within a container
 * or from a base H5mdDataSetBase object created by H5mdFrameDataSetBuilder.
 *
 * The data set dimensions and templated type is checked to confirm that it describes framed data.
 * In practice this affirms:
 *
 * 1. That the data set has at least 1 dimension (for frames).
 * 2. For non-primitive types that the remaining dimensions can store the data in row-major order
 *    (e.g. for RVecs the data set dimensions must be [numFrames][...][3] where middle
 *    [...] may be one or more dimensions for the frame).
 *
 * The full data set dimensions of the managed HDF5 data set is for the storage layout
 * of primitive data in the set in row-major order, where the outer dimension (0)
 * is for frame indexing. This function returns the dimensions of a single frame
 * of the data set.
 *
 * For a 1d data set of dimensions [30] the frame dimension is [] (each frame is a single value).
 * For a 3d data set of size [30, 50, 3] the frame dimension is [50, 3].
 *
 * \note For compound types such as BasicVector<float> which consists of 3 floats
 * the frame dimension is slightly more complex. A data set which was created to
 * store 30 frames, each storing 50 BasicVector<float> values has an actual dimension
 * of [30, 50, 3] and the frame dimension returned by this function is [50].
 *
 * I/O methods for this class use 1d ArrayRef buffers of ValueType. The size of the input buffers is
 * checked before the read or write operation. frameDims() is used to get the shape of each frame
 * in terms of ValueType and is used to allocate memory for these buffers.
 *
 * \tparam ValueType Type of values stored in data set (not limited to base primitives).
 */
template<typename ValueType>
class H5mdFrameDataSet : private H5mdDataSetBase<ValueType>
{
public:
    using Base = H5mdDataSetBase<ValueType>;

    //! \brief Constructor to manage a given \p dataSet.
    H5mdFrameDataSet(H5mdDataSetBase<ValueType>&& dataSet);

    //! \brief Constructor to open an existing data set (called \p name, in \p container).
    H5mdFrameDataSet(hid_t container, const char* name);

    //! \brief Destructor.
    ~H5mdFrameDataSet() noexcept;
    //! \brief Move constructor.
    H5mdFrameDataSet(H5mdFrameDataSet<ValueType>&&) noexcept;
    //! \brief Move assignment.
    H5mdFrameDataSet& operator=(H5mdFrameDataSet<ValueType>&&) noexcept;

    /*! \brief Return the dimensions of a single frame of the templated type in data set.
     *
     * The full data set dimensions of the managed HDF5 data set is for the storage layout
     * of primitive data in the set in row-major order, where the outer dimension (0)
     * is for frame indexing. This function returns the dimensions of a single frame
     * of the data set.
     *
     * For a 1d data set of dimensions [30] the frame dimension is []. For a 3d data set of
     * size [30, 50, 3] the frame dimension is [50, 3].
     *
     * \note For compound types such as BasicVector<float> which consists of 3 floats
     * the frame dimension is slightly more complex. A data set which was created to
     * store 30 frames, each storing 50 BasicVector<float> values has an actual dimension
     * of [30, 50, 3] and the frame dimension returned by this function is [50].
     */
    const DataSetDims& frameDims() const;

    /*! \brief Return the number of frames in the data set.
     *
     *  Number of frames refers to the size of the major axis of the data set. For example,
     *  a 1d data set with dimensions [30] obviously has 30 frames. A 3d data set with
     *  dimensions [30, 150, 3] also has 30 frames.
     *
     *  \throws gmx::FileIOError if the data set is 0-dimensional or if there was an error
     *      reading the dimensions.
     */
    hsize_t numFrames() const;

    /*! \brief Read data from frame at \p index into \p values.
     *
     * The output buffer \p values must have a size which is identical to the size
     * of a single frame in the data set as described by frameDims().
     *
     * For a 1d data set of size [30] the frame dimension is scalar, so \p values must
     * store a single value of the base type. For a 3d data set of size [30, 50, 3] the frame
     * dimension is [50, 3] so \p values must be of type ArrayRef<ValueType> with size 150.
     *
     * \note For compound types such as BasicVector<float> which consists of 3 floats
     * the frame dimension is slightly more complex. A data set which was created to
     * store 30 frames, each storing 50 BasicVector<float> values has an actual dimension
     * of [30, 50, 3] and the frame dimension is [50]. Thus \p values must be of type
     * ArrayRef<BasicVector<float>> and store exactly 50 values.
     *
     * \param[in]  index  Frame index to read data from.
     * \param[out] values Container of values to read data into.
     *
     * \throws gmx::FileIOError if the size of \p values does not match the frame dimensions
     *     or if an error occurred when reading the data.
     */
    void readFrame(hsize_t index, ArrayRef<ValueType> values);

    /*! \brief Write data from \p values into the next frame.
     *
     * The input buffer \p values must have a size which is identical to the size
     * of a single frame in the data set as described by frameDims().
     *
     * For a 1d data set of size [30] the frame dimension is scalar, so \p values must
     * store a single value of the base type. For a 3d data set of size [30, 50, 3] the frame
     * dimension is [50, 3] so \p values must be of type ArrayRef<const ValueType> with size 150.
     *
     * \note For compound types such as BasicVector<float> which consists of 3 floats
     * the frame dimension is slightly more complex. A data set which was created to
     * store 30 frames, each storing 50 BasicVector<float> values has an actual dimension
     * of [30, 50, 3] and the frame dimension is [50]. Thus \p values must be of type
     * ArrayRef<BasicVector<float>> and store exactly 50 values.
     *
     * \note HDF5 data sets have a maximum extent which is set when the data set is created.
     * We cannot write a new frame if the maximum extent along the frame axis has been reached.
     * This is by default set to be unlimited for data sets created with \c H5mdFrameDataSetBuilder,
     * in which case there is no limit on the amount of frames we can write (only the available
     * disk space will limit the number of frames).
     *
     * \param[in] values Container of values to write.
     *
     * \throws gmx::FileIOError if the size of \p values does not match the frame dimensions,
     *     if the maximum number of frames set for the data set has already been written,
     *     or if another error occurred when writing the data.
     */
    void writeNextFrame(ArrayRef<const ValueType> values);

private:
    //! \brief Return the data set extents for a given number of frames.
    const DataSetDims& extentForNumFrames(hsize_t numFrames);

    /*! \brief Internal description of a single frame within the data set.
     */
    class FrameDescription
    {
    public:
        FrameDescription(const DataSetDims& dataSetDims);
        ~FrameDescription()                                  = default;
        FrameDescription(const FrameDescription&)            = delete;
        FrameDescription& operator=(const FrameDescription&) = delete;
        FrameDescription(FrameDescription&&)                 = default;
        FrameDescription& operator=(FrameDescription&&)      = default;

        /*! \brief Construct and return a data space for reading a frame at \p frameIndex.
         *
         * This corresponds to `file_space_id` in the hdf5 documentation.
         *
         * \note The returned handle must be closed by the caller to avoid leaking resources.
         */
        hid_t fileDataSpaceForFrame(hsize_t frameIndex, hid_t dataSetHandle) noexcept;

        //! \brief Return the frame dimensions.
        const DataSetDims& dims() const { return dims_; }

        //! \brief Return the dimensions of primitives stored in a single frame.
        const DataSetDims& frameDimsPrimitive() const { return frameDimsPrimitive_; }

        //! \brief Return the memory data space for reading or writing a single frame.
        hid_t memoryDataSpace() const { return memoryDataSpace_; }

        //! \brief Return the total number of templated values in a single frame.
        hsize_t numValues() const { return numValues_; }

    private:
        /*! \brief Dimensions of a single frame of the templated ValueType in the data set.
         *
         * For a data set of floats with internal dimensions [30, 50, 3] this is [50, 3]. For
         * a data set of BasicVector<float> with internal dimensions [30, 50, 3] this is [50].
         */
        DataSetDims dims_;

        /*!< Number of templated values stored in a single frame of the data set.
         */
        hsize_t numValues_;

        /*! \brief Internal dimensions for a single frame of primitives in the data set.
         *
         * Example: For a data set of dimensions [30, 50, 3] this is [1, 50, 3] since
         * frames are the major axis.
         *
         * The dimensions are for the base primitives of the data set, not the templated
         * type. For BasicVector<double> the base type is double, so the frame dimensions
         * describe the storage layout of double[3].
         */
        DataSetDims frameDimsPrimitive_;

        //!< Memory data space for a single frame, called `mem_space_id` in hdf5 documentation.
        hid_t memoryDataSpace_;

        /*! \brief Handle to scope guard which closes the memory data space handle when the destructor is called.
         *
         * Scope guards have no move semantics, so in order for the class to have move semantics
         * we must store this indirectly. For the default move constructor this will be swapped
         * with a nullptr, and thus not result in an extra close of the moved handle. For default
         * move assignment this and the handle will be swapped together, and thus closed by the
         * moved-from object's destructor.
         *
         * \note Must be declared directly below and thus initialized just after the guarded
         * handle. If not the guard may not be active if a subsequent initialization fails,
         * resulting in a memory leak.
         */
        std::unique_ptr<H5mdGuard> memoryDataSpaceGuard_ = std::make_unique<H5mdGuard>(
                sg::make_scope_guard(H5mdCloser(memoryDataSpace_, H5Sclose)));

        /*! \brief Offset for hyperslab used to select frame.
         *
         * Created as a member variable once instead of upon every call to fileDataSpaceForFrame().
         */
        DataSetDims frameOffset_;
    };

    /*! \brief Full extent of the managed data set.
     *
     * At data set creation this is a copy of Base::dims(): the primitive dimensions of the data
     * set. We store this copy to resize the data set when adding new frames without creating a new
     * copy at every call.
     *
     * \note Index = 0 is for the number of frames, do not touch any other value.
     */
    DataSetDims extentDimsPrimitive_;

    //! \brief Description of frames in data set.
    FrameDescription frameDescription_;

    //!< Current number of frames in the data set.
    hsize_t numFrames_;
};

extern template class H5mdFrameDataSet<int32_t>;

extern template class H5mdFrameDataSet<int64_t>;

extern template class H5mdFrameDataSet<float>;

extern template class H5mdFrameDataSet<double>;

extern template class H5mdFrameDataSet<BasicVector<float>>;

extern template class H5mdFrameDataSet<BasicVector<double>>;

} // namespace gmx

#endif // GMX_FILEIO_H5MD_FRAMEDATASET_H
