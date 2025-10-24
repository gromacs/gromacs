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

/*! \brief Declarations of H5md time data block manipulation routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_TIMEDATABLOCK_H
#define GMX_FILEIO_H5MD_TIMEDATABLOCK_H

#include <hdf5.h>

#include <initializer_list>
#include <optional>
#include <string>

#include "external/scope_guard/scope_guard.h"

#include "gromacs/fileio/h5md/h5md_framedataset.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

template<typename ValueType>
class ArrayRef;
template<typename ValueType>
class BasicVector;
template<typename ValueType>
class H5mdTimeDataBlockBuilder;

/*! \brief Class which manages time-dependent data sets of the templated \p ValueType.
 *
 * Time-dependent data in H5md require the coexistance of multiple data sets:
 *
 * 1) A data set for the primary value type that is to be stored,
 * 2) A data set with the simulation step number for each stored value, and
 * 3) (Optionally) A data set with the simulation time for each stored value.
 *
 * This class manages this set of data sets as one, providing routines for writing
 * new frame data as well as reading existing values, steps or times. New data is
 * written to all data sets as one to keep them synchronized.
 *
 * Use the \c H5mdTimeDataBlockBuilder class to construct new objects with
 * certain features. When opening already existing blocks for reading or appending
 * data (i.e. when resuming simulations from checkpoints), use the \c H5mdTimeDataBlock
 * constructor.
 */
template<typename ValueType>
class H5mdTimeDataBlock
{
public:
    // Only make constructors which require the user to pass data sets available
    // to our builder class.
    friend class H5mdTimeDataBlockBuilder<ValueType>;

    /* \brief Open value, step and time data sets in group (called \p name, in \p container).
     *
     * Used when opening an existing H5md trajectory from disk.
     *
     * The data set names within the managed group are opened following the H5md specification:
     *  - value: For the data set storing the values.
     *  - step:  For the data set storing the steps.
     *  - time:  For the data set storing the times. (optional)
     *
     * If the value or step data set cannot be found within the group an error is raised.
     * Since the time data set is optional no error is raised if it is missing.
     *
     * \param[in] container Handle to container in which to open the group.
     * \param[in] name      Name of group to open.
     *
     * \throws FileIOError if the group does not exist, or if the value or step data sets
     * cannot be found within, or if the data set types do not match.
     */
    H5mdTimeDataBlock(hid_t container, const char* name);

    //! \brief Return whether or not this block manages a time data set.
    bool hasTime() const { return timeDataSet_.has_value(); }

    /*! \brief Return the number of frames written to the managed data sets.
     */
    int64_t numFrames() const { return numFrames_; }

    /*! \brief Return the step value of the given \p frameIndex.
     *
     * \param[in] frameIndex Index of frame to get the value for.
     *
     * \returns The step or std::nullopt if \p frameIndex is not within [0, numFrames).
     * \throws  FileIOError if there was an error reading the value.
     */
    std::optional<int64_t> readStepAtIndex(int64_t frameIndex);

    /*! \brief Return the time value of the given \p frameIndex.
     *
     * If the time data set is not managed this always returns std::nullopt.
     *
     * \param[in] frameIndex Index of frame to get the value for.
     *
     * \returns The time if available, or std::nullopt if \p frameIndex is not within
     *          [0, numFrames).
     * \throws  FileIOError if there was an error reading the value.
     */
    std::optional<double> readTimeAtIndex(int64_t frameIndex);

    /*! \brief Read the primary value of the given \p frameIndex into the given buffer if the index exists.
     *
     * \param[in]  frameIndex  Index of frame to get the value for.
     * \param[out] values      Buffer to read the value into.
     *
     * \returns Whether the index is within [0, numframes) and thus if the value was read into the buffer.
     * \throws  FileIOError if there was an error reading the value.
     */
    bool readValueAtIndex(int64_t frameIndex, ArrayRef<ValueType> values);

    /*! \brief Read the frame at \p frameIndex into the given buffers if the index exists.
     *
     * If \p step or \p time are nullptr their values are not read. If a time data set is
     * not managed \p time is never read.
     *
     * \param[in]  frameIndex Index of frame to read.
     * \param[out] values     Buffer to read values into.
     * \param[out] step       Buffer to read step into, or nullptr if not to read.
     * \param[out] time       Buffer to read time into, or nullptr if not to read.
     *
     * \returns Whether the index is within [0,  numFrames) and thus if the frame was read.
     * \throws  FileIOError if there was an error reading any value.
     */
    bool readFrame(int64_t frameIndex, ArrayRef<ValueType> values, int64_t* step, double* time);

    /*! \brief Write given values and step to their data sets as the next frame.
     *
     * \note This overload throws if a time data set is managed, and is designed
     * to be used by callers who know that they are not writing times.
     *
     * \param[in] values Values to write.
     * \param[in] step   Step to write.
     *
     * \throws FileIOError if there was an error writing the values or if a time
     * data set is managed.
     */
    void writeNextFrame(ArrayRef<const ValueType> values, int64_t step);

    /*! \brief Write given values to all managed data sets as the next frame.
     *
     * This function works both when a time data set is managed and when it is not.
     * If a time data set is managed the \p time is written to it. If a time data set
     * is not managed the time is never written, but other values are.
     *
     * \param[in] values Values to write.
     * \param[in] step   Step to write.
     * \param[in] time   Time to write, if managed.
     *
     * \throws FileIOError if there was an error writing the values.
     */
    void writeNextFrame(ArrayRef<const ValueType> values, int64_t step, double time);

private:
    /*! \brief Constructor for the managed data sets.
     *
     * \note A simpler route to create the data sets for a \c H5mdTimeDataBlock group
     * is provided in the \c H5mdTimeDataBlockBuilder class.
     *
     * \param[in] valueDataSet Data set for primary values.
     * \param[in] stepDataSet  Data set for steps.
     * \param[in] timeDataSet  (Optional) Data set for times.
     *
     * \throws FileIOError if the sets do not have the same number of values along
     * their primary axis.
     */
    H5mdTimeDataBlock(H5mdFrameDataSet<ValueType>&&                   valueDataSet,
                      H5mdScalarFrameDataSet<int64_t>&&               stepDataSet,
                      std::optional<H5mdScalarFrameDataSet<double>>&& timeDataSet);

    /*! \brief Private constructor used by the public constructor when opening the data sets.
     *
     * \param[in] blockGroupAndGuard Pair of the group in which to open the internal data sets,
     *                               and a scope guard which closes the group afterwards.
     */
    H5mdTimeDataBlock(decltype(makeH5mdGroupGuard(H5I_INVALID_HID))&& blockGroupAndGuard);

    //! \brief Data set containing the stored value at each frame.
    H5mdFrameDataSet<ValueType> valueDataSet_;

    //! \brief Data set containing the simulation step at each frame.
    H5mdScalarFrameDataSet<int64_t> stepDataSet_;

    //! \brief Data set containing the simulation time at each frame.
    //
    // The H5md specification allows this data set to be optional which we must conform to.
    std::optional<H5mdScalarFrameDataSet<double>> timeDataSet_;

    //! \brief Number of frames written to the data sets.
    int64_t numFrames_;
};

/*! \brief Builder class to create new \c H5mdTimeDataBlock objects.
 *
 * Creates new (empty) value, step and time data sets with the expected names. Also exposes
 * some options for the creation of these data sets.
 *
 * Once all desired options are set, use \c H5mdTimeBlockBuilder::build() to finalize all
 * data sets and return the resulting \c H5mdTimeDataBlock.
 */
template<typename ValueType>
class H5mdTimeDataBlockBuilder
{
public:
    /*! \brief Return a builder to create a \c H5mdTimeDataBlock (called \p name, in \p container).
     *
     * \throws gmx::FileIOError if there already exists an object at the given path.
     */
    H5mdTimeDataBlockBuilder(hid_t container, const std::string& name);

    //! \brief Set dimension for a single frame in the value data set.
    H5mdTimeDataBlockBuilder& withFrameDimension(ArrayRef<const hsize_t> dims);

    //! \copydoc withFrameDimension()
    H5mdTimeDataBlockBuilder& withFrameDimension(std::initializer_list<hsize_t> dims);

    //! \brief Create all data sets and return as a \c H5mdTimeDataBlock.
    H5mdTimeDataBlock<ValueType> build();

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(H5mdTimeDataBlockBuilder);

private:
    //!< Container for the block group in which to create all data sets.
    const hid_t blockGroup_;

    /*! \brief Handle to scope guard which closes the block group handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    H5mdGuard blockGroupGuard_ = sg::make_scope_guard(H5mdCloser(blockGroup_, H5Gclose));

    //!< Builder to create the value data set.
    H5mdFrameDataSetBuilder<ValueType> valueDataSetBuilder_;
};

extern template class H5mdTimeDataBlock<int32_t>;
extern template class H5mdTimeDataBlock<int64_t>;
extern template class H5mdTimeDataBlock<float>;
extern template class H5mdTimeDataBlock<double>;
extern template class H5mdTimeDataBlock<BasicVector<float>>;
extern template class H5mdTimeDataBlock<BasicVector<double>>;

extern template class H5mdTimeDataBlockBuilder<int32_t>;
extern template class H5mdTimeDataBlockBuilder<int64_t>;
extern template class H5mdTimeDataBlockBuilder<float>;
extern template class H5mdTimeDataBlockBuilder<double>;
extern template class H5mdTimeDataBlockBuilder<BasicVector<float>>;
extern template class H5mdTimeDataBlockBuilder<BasicVector<double>>;

} // namespace gmx

#endif // GMX_FILEIO_H5MD_TIMEDATABLOCK_H
