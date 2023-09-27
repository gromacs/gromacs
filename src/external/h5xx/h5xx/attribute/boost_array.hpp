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

#ifndef H5XX_ATTRIBUTE_BOOST_ARRAY
#define H5XX_ATTRIBUTE_BOOST_ARRAY

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <h5xx/attribute/attribute.hpp>
#include <h5xx/attribute/utility.hpp>
#include <h5xx/ctype.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/error.hpp>
#include <h5xx/policy/string.hpp>
#include <h5xx/utility.hpp>

#include <vector>

namespace h5xx {

/**
 * create and write fixed-size fundamental array type attribute for h5xx objects
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_fundamental<typename T::value_type>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    typedef typename T::value_type value_type;
    boost::array<hsize_t, 1> dims = {{ T::static_size }};

    delete_attribute(object, name);
    attribute attr(object, name, ctype<value_type>::hid(), dataspace(dims));
    attr.write(ctype<value_type>::hid(), &*value.begin());
}

/**
 * read fixed-size fundamental array type attributes from h5xx objects
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_fundamental<typename T::value_type>
>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    // open attribute
    attribute attr(object, name);
    dataspace space(attr);
    if (space.rank() != 1 || space.extents<1>()[0] != T::static_size) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has mismatching dataspace");
    }

    // read from opened object
    T value;
    attr.read(ctype<typename T::value_type>::hid(), &*value.begin());
    return value;
}

/**
 * create and write std::string arrays
 */
template <typename T, typename h5xxObject, typename StringPolicy>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_same<typename T::value_type, std::string>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value, StringPolicy policy = StringPolicy())
{
    enum { size = T::static_size };
    delete_attribute(object, name);

    // determine the maximum string size
    size_t str_size = 0;
    for (hsize_t i = 0; i < size; ++i) {
        str_size = std::max(str_size, value[i].size());  // include terminating NULL character
    }

    // create type, space and attribute
    boost::array<hsize_t, 1> dims = {{ size }};
    hid_t type_id = policy.make_type(str_size);
    attribute attr(object, name, type_id, dataspace(dims));

    // write attribute
    if (StringPolicy::is_variable_length) {
        std::vector<const char*> buffer(size);
        for (size_t i = 0; i < size; i++) {
            buffer[i] = value[i].c_str();
        }
        attr.write(type_id, &*buffer.begin());
    }
    else {
        std::vector<char> buffer(size * str_size);
        for (size_t i = 0; i < size; ++i) {
            value[i].copy(&buffer[i * str_size], str_size);
        }
        attr.write(type_id, &*buffer.begin());
    }
    H5Tclose(type_id);
}

template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_same<typename T::value_type, std::string>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    write_attribute(object, name, value, policy::string::null_terminated());
}

/**
 * read string array
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_same<typename T::value_type, std::string>
>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    enum { size = T::static_size };
    bool err = false;

    // open object and check dataspace
    attribute attr(object, name);
    dataspace space(attr);
    if (space.rank() != 1 || space.extents<1>()[0] != T::static_size) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has mismatching dataspace");
    }

    hid_t type_id;
    err |= (type_id= attr.get_type()) < 0;     // get copy of the attribute's type

    htri_t is_varlen_str;
    err |= (is_varlen_str = H5Tis_variable_str(type_id)) < 0;
    if (err) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" is not of a valid string type");
    }

    T value;
    if (!is_varlen_str){
        // create memory datatype with matching size and padding,
        // convert space padding to null padding
        hid_t mem_type_id = H5Tcopy(H5T_C_S1);
        hsize_t str_size = H5Tget_size(type_id);
        err |= H5Tset_size(mem_type_id, str_size) < 0;
        if (H5Tget_strpad(type_id) != H5T_STR_NULLTERM) {
            H5Tset_strpad(mem_type_id, H5T_STR_NULLPAD);
        }

        // read from opened object
        std::vector<char> buffer(str_size * size);
        attr.read(mem_type_id, &*buffer.begin());
        err |= H5Tclose(mem_type_id) < 0;

        char const* s = &*buffer.begin();
        for (unsigned int i = 0; i < size; ++i, s += str_size) {
            size_t len = strnlen(s, str_size);
            value[i] = std::string(s, len);
        }
    }
    else {
        hid_t mem_type_id = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
        std::vector<const char*> buffer(size);
        attr.read(mem_type_id, &*buffer.begin());
        for (int i = 0; i < size; i++) {
            value[i] = buffer[i];
        }
        err |= H5Tclose(mem_type_id) < 0;
    }
    err |= H5Tclose(type_id) < 0;

    if (err) {
        throw error("reading atrribute \"" + name);
    }

    return value;
}


/**
* create and write C string array types from h5xx objects
*/
template <typename T, typename h5xxObject, typename StringPolicy>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_same<typename T::value_type, char const*>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value, StringPolicy policy = StringPolicy())
{
    enum { size = T::static_size };
    size_t str_size = 0;
    for (hsize_t i = 0; i < size; ++i) {
        str_size = std::max(str_size, strlen(value[i]));
    }

    // remove attribute if it exists
    delete_attribute(object, name);

    boost::array<hsize_t, 1> dims = {{ size }};
    dataspace space(dims);

    hid_t type_id = policy.make_type(str_size);
    assert(type_id >= 0);
    attribute attr(object, name, type_id, space);

    if (StringPolicy::is_variable_length) {
        attr.write(type_id, &*value.begin());
    }
    else {
        std::vector<char> data(str_size * size);
        for (size_t i = 0; i < size; ++i) {
            strncpy(&data[i * str_size], value[i], str_size);
        }
        attr.write(type_id, &*data.begin());
    }

    if (H5Tclose(type_id) < 0) {
        throw error("closing datatype");
    }
}

template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::mpl::and_<
    is_array<T>, boost::is_same<typename T::value_type, char const*>
>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    write_attribute(object, name, value, policy::string::null_terminated());
}

} //namespace h5xx

#endif // ! H5XX_ATTRIBUTE_BOOST_ARRAY
