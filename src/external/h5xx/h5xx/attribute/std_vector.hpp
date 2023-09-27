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

#ifndef H5XX_ATTRIBUTE_STD_VECTOR
#define H5XX_ATTRIBUTE_STD_VECTOR

#include <vector>

#include <h5xx/attribute/attribute.hpp>
#include <h5xx/attribute/utility.hpp>
#include <h5xx/ctype.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

/**
 * create and write std::vector attributes of fundamental type
 **/
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_<
    is_vector<T>, boost::is_fundamental<typename T::value_type>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    typedef typename T::value_type value_type;

    // remove attribute if it exists
    delete_attribute(object, name);

    hid_t type_id = ctype<value_type>::hid();       // this ID must not be closed
    boost::array<hsize_t, 1> dims = {{ value.size() }};
    attribute attr(object, name, type_id, dataspace(dims));

    attr.write(type_id, &*value.begin());
}

/**
 * read std::vector attributes of fundamental type
 **/
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_vector<T>, boost::is_fundamental<typename T::value_type>
>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    typedef typename T::value_type value_type;

    // open attribute and check dataspace
    attribute attr(object, name);
    dataspace space(attr);
    if (space.rank() != 1) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has incompatible dataspace");
    }

    // read from opened object
    hsize_t size = space.extents<1>()[0];
    std::vector<value_type> value(size);
    attr.read(ctype<value_type>::hid(), &*value.begin());
    return value;
}

/**
 * create and write std::vector attributes of string type
 * FIXME add string policy
 **/
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_<
    is_vector<T>, boost::is_same<typename T::value_type, std::string>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    size_t size = value.size();
    boost::array<hsize_t, 1> dims = {{ value.size() }};

    // remove attribute if it exists
    delete_attribute(object, name);

    // size of longest string
    size_t str_size = 0;
    for (size_t i = 0; i < size; ++i) {
        str_size = std::max(str_size, value[i].size());
    }

    bool err = false;
    hid_t type_id = H5Tcopy(H5T_C_S1);
    err |= H5Tset_size(type_id, str_size) < 0;
    err |= H5Tset_strpad(type_id, H5T_STR_NULLTERM) < 0;

    attribute attr(object, name, type_id, dataspace(dims));

    std::vector<char> buffer(size * str_size);
    for (size_t i = 0; i < size; ++i){
        value[i].copy(&buffer[i * str_size], str_size);
    }

    attr.write(type_id, &*buffer.begin());
    err |= H5Tclose(type_id) < 0;
    if (err) {
        throw error("writing attribute \"" + name + "\" with ID " + boost::lexical_cast<std::string>(attr.hid()));
    }
}

/**
 * read std::vector attributes of string type
 **/

template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_vector<T>, boost::is_same<typename T::value_type, std::string>
>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    // open object and check dataspace
    attribute attr(object, name);
    dataspace space(attr);
    if (space.rank() != 1) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has incompatible dataspace");
    }

    bool err = false;
    hid_t type_id = attr.get_type();              // attribute's datatype on disk
    hid_t mem_type_id = H5Tcopy(H5T_C_S1);        // create memory datatype of same size
    hsize_t str_size = H5Tget_size(type_id);
    err |= H5Tset_size(mem_type_id, str_size) < 0;

    // read from opened object to buffer
    hsize_t size = space.extents<1>()[0];
    std::vector<char> buffer(str_size * size);
    err |= (H5Aread(attr.hid(), mem_type_id, &*buffer.begin())) < 0;
    if (err) {
        throw error("error while reading attribute \"" + name + "\"");
    }
    // close object
    err |= H5Tclose(mem_type_id) < 0;
    err |= H5Tclose(type_id) < 0;
    if (err) {
        throw error("reading atrribute \"" + name + "\" with ID " + boost::lexical_cast<std::string>(attr.hid()));
    }
    T value;
    value.reserve(size);
    char const* s = &*buffer.begin();
    for (size_t i = 0; i < size; ++i, s += str_size) {
        size_t len = strnlen(s, str_size);       // strings of str_len size are not '\0'-terminated
        value.push_back(std::string(s, len));    // copy len bytes from buffer
    }
    return value;
}

} //namespace h5xx

#endif // ! H5XX_ATTRIBUTE_STD_VECTOR
