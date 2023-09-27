/*
 * Copyright © 2014-2015 Klaus Reuter
 * Copyright © 2014      Felix Höfling
 * Copyright © 2014      Manuel Dibak
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_DATASET_SCALAR
#define H5XX_DATASET_SCALAR

#include <h5xx/ctype.hpp>
#include <h5xx/dataset/dataset.hpp>
#include <h5xx/dataset/utility.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/error.hpp>
#include <h5xx/policy/string.hpp>
#include <h5xx/policy/storage.hpp>
#include <h5xx/utility.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {

/**
 * create dataset of fundamental type at h5xx object
 *
 * Use compact storage by default.
 */
template <typename T, typename h5xxObject, typename StoragePolicy>
inline typename boost::enable_if<boost::is_fundamental<T>, dataset>::type
create_dataset(h5xxObject const& object, std::string const& name, StoragePolicy const& storage_policy = StoragePolicy())
{
    if (h5xx::exists_dataset(object, name))
    {
        throw error("dataset \"" + name + "\" of object \"" + get_name(object) + "\" does already exist");
    }
    return dataset(object, name, ctype<T>::hid(), dataspace(H5S_SCALAR), storage_policy);
}

// default value for template parameter StoragePolicy not supported in C++98
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_fundamental<T>, dataset>::type
create_dataset(h5xxObject const& object, std::string const& name)
{
    return create_dataset(object, name, ctype<T>::hid(), dataspace(H5S_SCALAR), policy::storage::compact());
}

/**
 * write fundamental type dataset at h5xx object
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if<boost::is_fundamental<T>, void>::type
write_dataset(h5xxObject const& object, std::string const& name, T const& value)
{
    dataset dset;
    if (h5xx::exists_dataset(object, name))
    {
        dset = dataset(object, name);
        if (!dataspace(dset).is_scalar()) {
            throw error("dataset \"" + name + "\" of object \"" + get_name(object) + "\" has non-scalar dataspace");
        }
    }
    else
    {
        //dset.create(object, name, ctype<T>::hid(), dataspace(H5S_SCALAR));
        throw error("dataset \"" + name + "\" of object \"" + get_name(object) + "\" does not exist");
    }
    dset.write(ctype<T>::hid(), &value);
}

/**
 * read fundamental type dataset from a given h5xx object
 */
template <typename T, typename h5xxObject>
inline typename boost::enable_if<boost::is_fundamental<T>, T>::type
read_dataset(h5xxObject const& object, std::string const& name)
{
    // open dataset and check dataspace
    dataset dset(object, name);
    if (!dataspace(dset).is_scalar()) {
        throw error("dataset \"" + name + "\" of object \"" + get_name(object) + "\" has non-scalar dataspace");
    }
    // read dataset
    T value;
    dset.read(ctype<T>::hid(), &value);
    return value;
}

// TODO : implement string support in analogy to attribute

} // namespace h5xx

#endif // ! H5XX_DATASET_SCALAR
