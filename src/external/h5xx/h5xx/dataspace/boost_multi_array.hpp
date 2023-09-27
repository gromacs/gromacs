/*
 * Copyright © 2014-2015 Klaus Reuter
 * Copyright © 2014      Felix Höfling
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_DATASPACE_MULTI_ARRAY
#define H5XX_DATASPACE_MULTI_ARRAY

#include <algorithm>

#include <h5xx/dataspace.hpp>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

template <typename T>
typename boost::enable_if<is_multi_array<T>, dataspace>::type
create_dataspace(T const& value)
{
    enum { rank = T::dimensionality };
    boost::array<hsize_t, rank> value_dims;
    std::copy(value.shape(), value.shape() + rank, value_dims.begin());
    return dataspace(value_dims);
}

} // namespace h5xx

#endif // ! H5XX_DATASPACE_MULTI_ARRAY
