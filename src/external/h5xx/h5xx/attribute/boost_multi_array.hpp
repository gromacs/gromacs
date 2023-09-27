/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2014 Manuel Dibak
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_ATTRIBUTE_MULTI_ARRAY
#define H5XX_ATTRIBUTE_MULTI_ARRAY

#include <algorithm>

#include <h5xx/attribute/attribute.hpp>
#include <h5xx/ctype.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>
#include <h5xx/attribute/utility.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

/**
* write attribute of multi-dimensional array type
*/
template <typename h5xxObject, typename T>
inline typename boost::enable_if<is_multi_array<T>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    typedef typename T::element value_type;
    enum { rank = T::dimensionality };

    // remove attribute if it exists
    delete_attribute(object, name);

    // create attribute with given dimensions
    boost::array<hsize_t, rank> dims;
    std::copy(value.shape(), value.shape() + rank, dims.begin());
    hid_t type_id = ctype<value_type>::hid();       // this ID must not be closed
    attribute attr(object, name, type_id, dataspace(dims));

    // write attribute
    attr.write(type_id, value.origin());
}

/**
* read attribute of multi-dimensional array type
*/
template <typename T, typename h5xxObject>
inline typename boost::enable_if<is_multi_array<T>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    typedef typename T::element value_type;
    enum { rank = T::dimensionality };

    // open object
    attribute attr(object, name);

    // check if rank of dataspace and rank of array to be returned are matching
    dataspace space(attr);
    if (!(space.rank() == rank)) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has mismatching dataspace");
    }

    // get extents of dataspace
    boost::array<hsize_t, rank> dims = space.extents<rank>();

    // create boost::multi_array of given extents for use as buffer
    boost::array<size_t, rank> shape;
    std::copy(dims.begin(), dims.begin() + rank, shape.begin());
    boost::multi_array<value_type, rank> value(shape);

    // read attribute to buffer
    attr.read(ctype<value_type>::hid(), value.origin());

    return value;
}

} // namespace h5xx

#endif // ! H5XX_ATTRIBUTE_MULTI_ARRAY
