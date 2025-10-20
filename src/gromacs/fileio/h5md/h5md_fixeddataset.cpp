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

/*! \brief Definitions of H5md fixed data set class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "h5md_fixeddataset.h"

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "h5md_error.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

//! \brief Return the dimensions of the templated ValueType for the primitive type data set dimensions.
template<typename ValueType>
static DataSetDims primitiveDimsToValueTypeDims(const DataSetDims& dims)
{
    if constexpr (std::is_same_v<ValueType, BasicVector<float>>
                  || std::is_same_v<ValueType, BasicVector<double>>)
    {
        // A data set of type and shape BasicVector<T>[50][30] has a primitive type
        // and shape T[50][30][3], so verify that the input primitive dims has this
        // form and return without the inner value.
        throwUponH5mdError(dims.empty(),
                           "Data set dimensions for BasicVector<T> must be at least 1");
        throwUponH5mdError(dims.back() != DIM,
                           "Innermost dimension of data set for BasicVector<T> must be 3");
        return DataSetDims(dims.cbegin(), dims.cend() - 1);
    }
    else if constexpr (std::is_arithmetic_v<ValueType>)
    {
        // For simple primitives the dimensions are unchanged.
        return dims;
    }
    else
    {
        // All of the compiled types should be handled explicitly in this function, so this
        // branch should be unreachable. We verify this via a static_assert check which will
        // not compile if reached. NOTE With some compilers (including GCC 12) all branches
        // are entered and verified even if no templated function instantiation reaches it.
        // To work around this we define a dummy type here to compare our templated type against.
        // Once gcc 13 or higher is required, replace this with
        // static_assert(false,
        class Unreachable;
        static_assert(std::is_same_v<ValueType, Unreachable>,
                      "Unsupported type for primitiveDimsToValueTypeDims()");
    }
}

//! \brief Return the number of values in a data set with templated type dimensions \p dims.
static hsize_t numValuesInDataSet(const DataSetDims& dims)
{
    hsize_t numValues = 1;
    for (const hsize_t d : dims)
    {
        numValues *= d;
    }
    return numValues;
}

template<typename ValueType>
H5mdFixedDataSet<ValueType>::H5mdFixedDataSet(H5mdDataSetBase<ValueType>&& dataSet) :
    Base{ std::move(dataSet) },
    dims_{ primitiveDimsToValueTypeDims<ValueType>(Base::dims()) },
    numValues_{ numValuesInDataSet(dims_) }
{
}

template<typename ValueType>
H5mdFixedDataSet<ValueType>::H5mdFixedDataSet(const hid_t container, const char* name) :
    Base(container, name),
    dims_{ primitiveDimsToValueTypeDims<ValueType>(Base::dims()) },
    numValues_{ numValuesInDataSet(dims_) }
{
}

template<typename ValueType>
const DataSetDims& H5mdFixedDataSet<ValueType>::dims() const
{
    return dims_;
}

template<typename ValueType>
hsize_t H5mdFixedDataSet<ValueType>::numValues() const
{
    return numValues_;
}

template<typename ValueType>
void H5mdFixedDataSet<ValueType>::readData(ArrayRef<ValueType> data) const
{
    throwUponH5mdError(data.size() != numValues_,
                       formatString("Cannot read frame into buffer of incorrect size: "
                                    "size of data set is %llu values but size of buffer is %lu",
                                    static_cast<unsigned long long>(numValues_),
                                    data.size()));
    throwUponH5mdError(
            H5Dread(this->id(), this->nativeDataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0,
            "Error reading data.");
}

template<typename ValueType>
void H5mdFixedDataSet<ValueType>::writeData(ArrayRef<const ValueType> data) const
{
    throwUponH5mdError(data.size() != numValues_,
                       formatString("Cannot write buffer of incorrect size into data set: "
                                    "size of data set is %llu values but size of buffer is %lu",
                                    static_cast<unsigned long long>(numValues_),
                                    data.size()));
    throwUponH5mdError(
            H5Dwrite(this->id(), this->dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0,
            "Error writing data.");
}

template class H5mdFixedDataSet<int32_t>;

template class H5mdFixedDataSet<int64_t>;

template class H5mdFixedDataSet<float>;

template class H5mdFixedDataSet<double>;

template class H5mdFixedDataSet<BasicVector<float>>;

template class H5mdFixedDataSet<BasicVector<double>>;

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
