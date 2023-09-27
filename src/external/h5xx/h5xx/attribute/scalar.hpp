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

#ifndef H5XX_ATTRIBUTE_SCALAR
#define H5XX_ATTRIBUTE_SCALAR

#include <h5xx/attribute/attribute.hpp>
#include <h5xx/ctype.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>
#include <h5xx/attribute/utility.hpp>
#include <h5xx/policy/string.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

 /**
*create and write fundamental type attribute on h5xx objects
*/
template <typename h5xxObject, typename T>
inline typename boost::enable_if<boost::is_fundamental<T>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    delete_attribute(object, name);
    attribute attr(object, name, ctype<T>::hid(), dataspace(H5S_SCALAR));
    attr.write(ctype<T>::hid(), &value);
}

/**
 * read fundamental type attribute of a given h5xx::object
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_fundamental<T>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    // open attribute and check dataspace
    attribute attr(object, name);
    if (!dataspace(attr).is_scalar()) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has non-scalar dataspace");
    }

    // read attribute
    T value;
    attr.read(ctype<T>::hid(), &value);
    return value;
}

/**
 * create and write std::string attribute for h5xx object
 */
template <typename T, typename h5xxObject, typename StringPolicy> // only in C++11: StringPolicy = policy::string::null_terminated
inline typename boost::enable_if<boost::is_same<T, std::string>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value, StringPolicy policy = StringPolicy())
{
    hid_t type_id = policy.make_type(value.size());
    delete_attribute(object, name);
    attribute attr(object, name, type_id, dataspace(H5S_SCALAR));
    if (StringPolicy::is_variable_length) {
        char const* p = value.c_str();
        attr.write(type_id, &p);
    }
    else {
        attr.write(type_id, value.data());
    }

    if (H5Tclose(type_id) < 0) {
        throw error("closing datatype");
    }
}

template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_same<T, std::string>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T const& value)
{
    write_attribute(object, name, value, policy::string::null_terminated());
}

/**
 * create and write C string attribute for h5xx objects
 */
template <typename T, typename h5xxObject, typename StringPolicy>
inline typename boost::enable_if<boost::is_same<T, char const*>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T value, StringPolicy policy = StringPolicy())
{

    // remove attribute if it exists
    delete_attribute(object, name);
    hid_t type_id = policy.make_type(strlen(value));
    attribute attr(object, name, type_id, dataspace(H5S_SCALAR));

    // write data
    if (StringPolicy::is_variable_length) {
        attr.write(type_id, &value);
    }
    else {
        attr.write(type_id, &*value);
    };

    if (H5Tclose(type_id) < 0) {
        throw error("closing datatype");
    }
}

template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_same<T, char const*>, void>::type
write_attribute(h5xxObject const& object, std::string const& name, T value)
{
    write_attribute(object, name, value, policy::string::null_terminated());
}

/**
 * read std::string attribute for h5xx objects
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_same<T, std::string>, T>::type
read_attribute(h5xxObject const& object, std::string const& name)
{
    // open object
    attribute attr(object, name);
    if (!dataspace(attr).is_scalar()) {
        throw error("attribute \"" + name + "\" of object \"" + get_name(object) + "\" has non-scalar dataspace");
    }

    hid_t type_id = attr.get_type();

    // check if string is of variable or fixed size type
    bool err = false;
    htri_t is_varlen_str;
    err |= (is_varlen_str =  H5Tis_variable_str(type_id)) < 0;
    if (err) {
        throw error("attribute \"" + name + "\" is not a valid string type");
    }

    std::string value;
    if (!is_varlen_str) {
        // read fixed-size string, allocate space in advance and let the HDF5
        // library take care about NULLTERM and NULLPAD strings
        hsize_t size = H5Tget_size(type_id);
        hid_t mem_type_id = H5Tcopy(H5T_C_S1);
        err |= H5Tset_size(mem_type_id, size + 1) < 0;  //one extra character for zero- and space-padded strings
        if (err) {
            throw error("setting size of mem_type_id");
        }
        value.resize(size, std::string::value_type());
        attr.read(mem_type_id, &*value.begin());
        err |= H5Tclose(mem_type_id);

    }  else {
        // read variable-length string, memory will be allocated by HDF5 C
        // library and must be freed by us
        hid_t mem_type_id = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
        char *c_str;
        attr.read(mem_type_id, &c_str);
        if (!err) {
            value = c_str;  // copy '\0'-terminated string
            free(c_str);
        }
    }

    err |= H5Tclose(type_id);
    if (err) {
        throw error("reading atrribute \"" + name + "\" with ID " + boost::lexical_cast<std::string>(attr.hid()));
    }
    return value;
}

} // namespace h5xx

#endif // ! H5XX_ATTRIBUTE_SCALAR
