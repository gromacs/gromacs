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

#ifndef H5XX_DATASPACE_STD_VECTOR
#define H5XX_DATASPACE_STD_VECTOR

#include <vector>

#include <h5xx/dataspace.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, dataspace>::type
create_dataspace(T const& value)
{
    std::vector<hsize_t> value_dims;
    value_dims.push_back(value.size());
    return dataspace(value_dims);
}

} // namespace h5xx

#endif // ! H5XX_DATASPACE_STD_VECTOR
