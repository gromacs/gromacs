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

/*! \brief Declarations of H5md data set builder routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#ifndef GMX_FILEIO_H5MD_DATASETBUILDER_H
#define GMX_FILEIO_H5MD_DATASETBUILDER_H

#include <hdf5.h>

#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_error.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

//! \brief Data set compression options
enum class H5mdCompression : int
{
    Uncompressed,      //!< No compression
    LosslessNoShuffle, //!< Lossless, only gzip (deflate)
    LosslessShuffle    //!< Lossless, byte shuffle followed by gzip (deflate)
};

// gzip (deflate) compression level: 1 = fastest, 9 = highest compression
constexpr int c_h5mdDeflateCompressionLevel = 1;

/*! \brief Builder class for H5md data sets.
 *
 * This class facilitates setting various options such as dimensions for a data set.
 * When constructing the data set with a call to build() these options are validated
 * and the data set created with a type derived from the templated type.
 *
 * For non-primitive data types such as BasicVector<float> the data set is created
 * with the base primitive as its type and the data type size is added as the innermost
 * dimension to the data storage: This keeps the data layout in row-major order. For
 * example, the storage layout for a data set of BasicVector<float> with dimensions
 * [20][50] is float[20][50][3].
 *
 * \tparam ValueType Native type to create data set for.
 */
template<typename ValueType>
class H5mdDataSetBuilder
{
public:
    //! \brief Construct a builder for a data set with \p name in \p container.
    //
    // \throws gmx::FileIOError if an object with \p name already exists in \p container.
    H5mdDataSetBuilder(const hid_t container, const std::string& name) :
        container_(container), name_(name)
    {
        throwUponH5mdError(
                objectExists(container, name.c_str()),
                "Cannot create data set: object with given name already exists in the container");
        accessPropertyList_   = H5Pcreate(H5P_DATASET_ACCESS);
        creationPropertyList_ = H5Pcreate(H5P_DATASET_CREATE);
    }

    ~H5mdDataSetBuilder() noexcept
    {
        if (handleIsValid(accessPropertyList_))
        {
            H5Pclose(accessPropertyList_);
        }
        if (handleIsValid(creationPropertyList_))
        {
            H5Pclose(creationPropertyList_);
        }
    }

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(H5mdDataSetBuilder);

    //! \brief Set \p compression for data set.
    H5mdDataSetBuilder& withCompression(const H5mdCompression compression)
    {
        switch (compression)
        {
            case H5mdCompression::Uncompressed: break;
            case H5mdCompression::LosslessShuffle:
                H5Pset_shuffle(creationPropertyList_);
                [[fallthrough]];
            case H5mdCompression::LosslessNoShuffle:
                H5Pset_deflate(creationPropertyList_, c_h5mdDeflateCompressionLevel);
                break;
        }

        return *this;
    }

    //! \brief Set data set dimensions (required).
    H5mdDataSetBuilder& withDimension(ArrayRef<const hsize_t> dims)
    {
        dims_.assign(dims.begin(), dims.end());
        return *this;
    }

    //! \copydoc withDimension()
    H5mdDataSetBuilder& withDimension(std::initializer_list<hsize_t> dims)
    {
        return withDimension(ArrayRef<const hsize_t>(dims.begin(), dims.end()));
    }

    //! \brief Set maximum dimension for data set (default: same as data set dimension).
    H5mdDataSetBuilder& withMaxDimension(ArrayRef<const hsize_t> maxDims)
    {
        maxDims_.assign(maxDims.begin(), maxDims.end());
        return *this;
    }

    //! \copydoc withMaxDimension()
    H5mdDataSetBuilder& withMaxDimension(std::initializer_list<hsize_t> maxDims)
    {
        return withMaxDimension(ArrayRef<const hsize_t>(maxDims.begin(), maxDims.end()));
    }

    //! \brief Set chunk dimension for data set (default: same as data set dimension).
    //
    // \note The chunk size must be >=1 along each dimension. Automatic setting of this
    // option adjusts for this but if the chunk dimension is set explicitly through this
    // method a 0 dimension value will result in a throw when calling build().
    H5mdDataSetBuilder& withChunkDimension(ArrayRef<const hsize_t> chunkDims)
    {
        chunkDims_.assign(chunkDims.begin(), chunkDims.end());
        return *this;
    }

    //! \copydoc withChunkDimension()
    H5mdDataSetBuilder& withChunkDimension(std::initializer_list<hsize_t> chunkDims)
    {
        return withChunkDimension(ArrayRef<const hsize_t>(chunkDims.begin(), chunkDims.end()));
    }

    //! \brief Set the data set to use a fixed size string with maximum length \p maxLength.
    //
    // \note The max length must be positive and count the null-terminator character.
    // If neither withMaxStringLength() nor withVariableStringLength() is called,
    // the default is variable-length strings for string data sets.
    H5mdDataSetBuilder& withMaxStringLength(const int maxLength)
    {
        // Use int to prevent the integer overflow if passed a negative value
        throwUponH5mdError(
                maxLength <= 0,
                "Cannot create fixed-size string data set with non-positive maximum length");
        maxStringLength_ = static_cast<size_t>(maxLength);
        return *this;
    }

    //! \brief Set the data set to use variable length strings.
    H5mdDataSetBuilder& withVariableStringLength()
    {
        maxStringLength_ = std::nullopt;
        return *this;
    }

    //! \brief Set \p unit attribute for data set values.
    H5mdDataSetBuilder& withUnit(std::string_view unit)
    {
        unit_ = unit;
        return *this;
    }

    //! \brief Finalize all set options, then build and return the data set.
    //
    // \throws gmx::FileIOError if any set options are incorrect or incompatible with other
    // options, or if an error occurred when creating the data set.
    H5mdDataSetBase<ValueType> build()
    {
        // If we ever need the "SCALAR" data set type (HDF5 data sets for single value storage)
        // this is where we branch off
        throwUponH5mdError(dims_.empty(), "Cannot create data set for 0 dimensions");

        // NOTE: This call finalizes dims_ and other dimension vectors and must be done
        // before they are used to create data set parameters
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(finalizeDataType());

        if (maxDims_.empty())
        {
            maxDims_ = dims_;
        }
        throwUponH5mdError(maxDims_.size() != dims_.size(),
                           "Inconsistent input when creating data set: "
                           "maxDims must be of same dimension as data set dims.");

        if (chunkDims_.empty())
        {
            chunkDims_ = dims_;
            for (hsize_t& d : chunkDims_)
            {
                if (d == 0)
                {
                    d = 1;
                }
            }
        }
        throwUponH5mdError(chunkDims_.size() != dims_.size(),
                           "Inconsistent input when creating data set: "
                           "chunkDims must be of same dimension as data set dims.");

        throwUponH5mdError(H5Pset_chunk(creationPropertyList_, chunkDims_.size(), chunkDims_.data()) < 0,
                           "Cannot set chunk dimensions when creating data set.");
        throwUponH5mdError(
                H5Pset_chunk_opts(creationPropertyList_, H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS) < 0,
                "Cannot set chunk options when creating data set.");

        size_t cacheSize = H5Tget_size(dataType);
        for (const hsize_t d : chunkDims_)
        {
            cacheSize *= d;
        }
        throwUponH5mdError(
                H5Pset_chunk_cache(accessPropertyList_, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, cacheSize, H5D_CHUNK_CACHE_W0_DEFAULT)
                        < 0,
                "Cannot set chunk cache size when creating data set.");

        const auto [dataSpace, dataSpaceGuard] =
                makeH5mdDataSpaceGuard(H5Screate_simple(dims_.size(), dims_.data(), maxDims_.data()));
        throwUponInvalidHid(dataSpace, "Cannot create data space for data set.");

        const hid_t dataSetHandle = H5Dcreate(
                container_, name_.c_str(), dataType, dataSpace, H5P_DEFAULT, creationPropertyList_, accessPropertyList_);
        throwUponInvalidHid(dataSetHandle, "Cannot create data set.");

        if (!unit_.empty())
        {
            setAttribute(dataSetHandle, c_unitAttributeKey, unit_.c_str());
        }

        // Responsibility for closing `dataSetHandle` is taken by the `H5mdDataSetBase<ValueType>`
        return H5mdDataSetBase<ValueType>(dataSetHandle);
    }

private:
    //! \brief Finalize dimensions for and return the data type corresponding to the templated ValueType.
    //
    // \note For some none-primitive types such as BasicVectors the data storage should be in row-major
    // order, with the data type storage being the innermost dimension. We append such dimensions to
    // dimension arrays here.
    hid_t finalizeDataType()
    {
        if constexpr (std::is_same_v<ValueType, gmx::BasicVector<float>>)
        {
            appendDimension(DIM);
            return hdf5DataTypeFor<float>();
        }
        else if constexpr (std::is_same_v<ValueType, gmx::BasicVector<double>>)
        {
            appendDimension(DIM);
            return hdf5DataTypeFor<double>();
        }
        else if constexpr (std::is_same_v<ValueType, std::string> || std::is_same_v<ValueType, char*>)
        {
            if (maxStringLength_.has_value())
            {
                return hdf5DataTypeForFixedSizeString(maxStringLength_.value());
            }
            else
            {
                return hdf5DataTypeForVariableSizeString();
            }
        }
        else
        {
            return hdf5DataTypeFor<ValueType>();
        }
    }

    //! \brief Append \p value to all dimension containers.
    void appendDimension(const hsize_t value)
    {
        dims_.push_back(value);
        if (!maxDims_.empty())
        {
            maxDims_.push_back(value);
        }
        if (!chunkDims_.empty())
        {
            chunkDims_.push_back(value);
        }
    }

    //!< Container in which to create the data set.
    hid_t container_;

    //!< Name of data set.
    std::string name_;

    //!< When set, the maximum length of a string data set (default: variable-length string)
    std::optional<size_t> maxStringLength_ = std::nullopt;

    //!< Dimensions to create data set with.
    std::vector<hsize_t> dims_;

    //!< Maximum dimensions to create data set with.
    std::vector<hsize_t> maxDims_;

    //!< Chunk dimensions of data set storage.
    std::vector<hsize_t> chunkDims_;

    //!< Property list for access options to create data set with.
    hid_t accessPropertyList_;

    //!< Property list for creation options to create data set with.
    hid_t creationPropertyList_;

    //!< Unit for values in data set.
    std::string unit_;
};

} // namespace gmx

#endif // GMX_FILEIO_H5MD_DATASETBUILDER_H
