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

/*! \brief Definitions of H5md data set base class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_datasetbase.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "h5md_error.h"
#include "h5md_guard.h"
#include "h5md_type.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

static hsize_t getNumDims(const hid_t dataSetHandle)
{
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSetHandle));
    return H5Sget_simple_extent_ndims(dataSpace);
}

template<typename ValueType>
class H5mdDataSetBase<ValueType>::Impl
{
public:
    /*! \brief Constructor implementation object from the given \p dataSetHandle.
     *
     * \warning The given \p dataSetHandle is required to be valid. All constructors
     * of the base class must verify this before creating this implementation object.
     */
    Impl(const hid_t dataSetHandle) :
        dataSet_{ dataSetHandle },
        dataType_{ H5Dget_type(dataSet_) },
        nativeDataType_{ H5Tget_native_type(dataType_, H5T_DIR_DEFAULT) },
        numDims_{ getNumDims(dataSet_) }
    {
        // Verify that the primitive type of the data set matches the templated type.
        // For compound types with defined shapes also assert that the data layout
        // is row-major and of correct dimensions.
        if constexpr (std::is_same_v<ValueType, gmx::BasicVector<float>>)
        {
            throwUponH5mdError(numDims_ == 0 || dims().back() != DIM,
                               "Could not open data set: inner dimension of data set must = 3 "
                               "for BasicVector<float>");
            throwUponH5mdError(!valueTypeIsDataType<float>(nativeDataType_),
                               "Could not open data set: compiled type parameter does not match "
                               "the primitive type of the data set");
        }
        else if constexpr (std::is_same_v<ValueType, gmx::BasicVector<double>>)
        {
            throwUponH5mdError(numDims_ == 0 || dims().back() != DIM,
                               "Could not open data set: inner dimension of data set must = 3 "
                               "for BasicVector<double>");
            throwUponH5mdError(!valueTypeIsDataType<double>(nativeDataType_),
                               "Could not open data set: compiled type parameter does not match "
                               "the primitive type of the data set");
        }
        else
        {
            throwUponH5mdError(!valueTypeIsDataType<ValueType>(nativeDataType_),
                               "Could not open data set: compiled type parameter does not match "
                               "the primitive type of the data set");
        }
    }

    Impl(const Impl&)            = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&)                 = delete;
    Impl& operator=(Impl&&)      = delete;

    DataSetDims dims() const
    {
        DataSetDims dataSetDims(numDims_, 0);
        const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet_));
        throwUponH5mdError(H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr) < 0,
                           "Could not read dimensions of data set");
        return dataSetDims;
    }

    //!< Handle to the managed data set.
    const hid_t dataSet_;

    /*! \brief Scope guard which closes the data set handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    const H5mdGuard dataSetGuard_ = sg::make_scope_guard(H5mdCloser(dataSet_, H5Dclose));

    //!< Handle to the data type of the managed data set used when writing the set.
    const hid_t dataType_;

    /*! \brief Scope guard which closes the data type handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    const H5mdGuard dataTypeGuard_ = sg::make_scope_guard(H5mdCloser(dataType_, H5Tclose));

    /*! \brief Handle to the immutable native data type of the managed data set.
     *
     * This may be different from dataType_ when opening a data set from a file which was created
     * on another machine. When reading data from a data set this is used along with dataType_ to
     * convert the data into the correct type of the current compiler target (i.e. native).
     */
    const hid_t nativeDataType_;

    /*! \brief Scope guard which closes the native data type handle when the destructor is called.
     *
     * \note Must be declared directly below and thus initialized just after the guarded
     * handle. If not the guard may not be active if a subsequent initialization fails,
     * resulting in a memory leak.
     */
    const H5mdGuard nativeDataTypeGuard_ = sg::make_scope_guard(H5mdCloser(nativeDataType_, H5Tclose));

    //!< Number of dimensions of data set.
    const hsize_t numDims_;
};

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(const hid_t dataSetHandle)
{
    throwUponInvalidHid(dataSetHandle, "Invalid handle to data set.");
    impl_ = std::make_unique<Impl>(dataSetHandle);
}

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(const hid_t container, const char* name)
{
    const hid_t dataSetHandle = H5Dopen(container, name, H5P_DEFAULT);
    throwUponInvalidHid(dataSetHandle, gmx::formatString("Cannot open data set with name %s.", name));
    impl_ = std::make_unique<Impl>(dataSetHandle);
}

template<typename ValueType>
H5mdDataSetBase<ValueType>::~H5mdDataSetBase() noexcept = default;

template<typename ValueType>
H5mdDataSetBase<ValueType>::H5mdDataSetBase(H5mdDataSetBase<ValueType>&&) noexcept = default;

template<typename ValueType>
H5mdDataSetBase<ValueType>& H5mdDataSetBase<ValueType>::operator=(H5mdDataSetBase<ValueType>&&) noexcept = default;

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::id() const
{
    return impl_->dataSet_;
}

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::dataType() const
{
    return impl_->dataType_;
}

template<typename ValueType>
hid_t H5mdDataSetBase<ValueType>::nativeDataType() const
{
    return impl_->nativeDataType_;
}

template<typename ValueType>
DataSetDims H5mdDataSetBase<ValueType>::dims() const
{
    return impl_->dims();
}

template class H5mdDataSetBase<int32_t>;

template class H5mdDataSetBase<int64_t>;

template class H5mdDataSetBase<float>;

template class H5mdDataSetBase<double>;

template class H5mdDataSetBase<gmx::BasicVector<float>>;

template class H5mdDataSetBase<gmx::BasicVector<double>>;

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
