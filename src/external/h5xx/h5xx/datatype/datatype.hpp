/*
 * Copyright © 2013-2014 Felix Höfling
 * Copyright © 2014-2015 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#include <h5xx/ctype.hpp>

#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#ifndef H5XX_DATATYPE_DATATYPE_HPP
#define H5XX_DATATYPE_DATATYPE_HPP

namespace h5xx {

/**
 * Wrapper class for the HDF5 datatype.
 */
class datatype
{
public:
    /** default constructor */
    datatype() : type_id_(-1) {}

    /** create datatype object from a HDF5 type_id */
    datatype(hid_t type_id) : type_id_(type_id) {}

    /** create datatype object from a Boost multi_array */
    template <class T>
    datatype(T array, typename boost::enable_if<is_multi_array<T> >::type* dummy = 0);

    /** create datatype object from a Boost array */
    template <class T>
    datatype(T array, typename boost::enable_if<is_array<T> >::type* dummy = 0);

    /** create datatype object from a std::vector */
    template <class T>
    datatype(T vector, typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> > >::type* dummy = 0);

    /** return the HDF5 type ID */
    hid_t get_type_id() const;

protected:
    /** HDF5 type ID */
    hid_t type_id_;
};

template <class T>
datatype::datatype(T array, typename boost::enable_if<is_multi_array<T> >::type* dummy)
{
    typedef typename T::element value_type;
    type_id_ = ctype<value_type>::hid();
}

template <class T>
datatype::datatype(T array, typename boost::enable_if<is_array<T> >::type* dummy)
{
    typedef typename T::value_type value_type;
    type_id_ = ctype<value_type>::hid();
}

template <class T>
datatype::datatype(T vector, typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> > >::type* dummy)
{
    typedef typename T::value_type value_type;
    type_id_ = ctype<value_type>::hid();
}

inline hid_t datatype::get_type_id() const
{
    return type_id_;
}

} // namespace h5xx

#endif // ! H5XX_DATATYPE_DATATYPE_HPP
