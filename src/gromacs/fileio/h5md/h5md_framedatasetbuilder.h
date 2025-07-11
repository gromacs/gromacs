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

/*! \brief Declarations of H5md frame data set builder routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_FRAMEDATASETBUILDER_H
#define GMX_FILEIO_H5MD_FRAMEDATASETBUILDER_H

#include <hdf5.h>

#include <vector>

#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \brief Builder class for framed H5md data sets.
 *
 * This constructs data sets with the outermost dimension being for frames.
 * The dimension for each frame can be specified using withFrameDimension()
 * and the number of frames in the created data set controlled with withNumFrames()
 * and withMaxNumFrames().
 *
 * By default data sets are created with 0 initial frames and an unlimited maximum
 * number of frames.
 *
 * See also the documentation for H5mdDataSetBuilder.
 *
 * \tparam ValueType Native type to create data set for.
 */
template<typename ValueType>
class H5mdFrameDataSetBuilder : private H5mdDataSetBuilder<ValueType>
{
public:
    using Base = H5mdDataSetBuilder<ValueType>;

    // Inherit constructor with no further changes
    using Base::Base;

    //! \brief Set dimension for a single frame in the data set.
    H5mdFrameDataSetBuilder& withFrameDimension(ArrayRef<const hsize_t> dims)
    {
        frameDims_.assign(dims.begin(), dims.end());
        return *this;
    }

    //! \copydoc withFrameDimension()
    H5mdFrameDataSetBuilder& withFrameDimension(std::initializer_list<hsize_t> dims)
    {
        return withFrameDimension(ArrayRef<const hsize_t>(dims.begin(), dims.end()));
    }

    //! \brief Set number of frames to create the data set with (default is 0).
    H5mdFrameDataSetBuilder& withNumFrames(const hsize_t numFrames)
    {
        numFrames_ = numFrames;
        return *this;
    }

    //! \brief Set maximum number of frames for the created data set (default is unlimited).
    H5mdFrameDataSetBuilder& withMaxNumFrames(const hsize_t maxNumFrames)
    {
        maxNumFrames_ = maxNumFrames;
        return *this;
    }

    //! \brief Create the data set and return a handle to it.
    hid_t build()
    {
        const std::vector<hsize_t> dims = [&]()
        {
            std::vector<hsize_t> dims = { numFrames_ };
            for (const hsize_t d : frameDims_)
            {
                dims.push_back(d);
            }
            return dims;
        }();
        Base::withDimension(dims);

        const std::vector<hsize_t> maxDims = [&]()
        {
            std::vector<hsize_t> maxDims = dims;
            maxDims[0]                   = maxNumFrames_;
            return maxDims;
        }();
        Base::withMaxDimension(maxDims);

        const std::vector<hsize_t> chunkDims = [&]()
        {
            std::vector<hsize_t> chunkDims = dims;
            chunkDims[0]                   = frameChunkSize_;
            return chunkDims;
        }();
        Base::withChunkDimension(chunkDims);

        return Base::build();
    }

private:
    //!< Dimensions of a single frame in the data set.
    std::vector<hsize_t> frameDims_;

    //!< Number of frames to use for each chunk in the data set.
    hsize_t frameChunkSize_ = 1;

    //!< Number of frames to create the data set with.
    hsize_t numFrames_ = 0;

    //!< Maximum number of frames for the data set.
    hsize_t maxNumFrames_ = H5S_UNLIMITED;
};

} // namespace gmx

#endif // GMX_FILEIO_H5MD_FRAMEDATASETBUILDER_H
