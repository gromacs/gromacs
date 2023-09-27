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

#ifndef H5XX_DATASET_STD_VECTOR
#define H5XX_DATASET_STD_VECTOR

#include <algorithm>
#include <vector>

#include <h5xx/ctype.hpp>
#include <h5xx/dataset/dataset.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/policy/storage.hpp>
#include <h5xx/utility.hpp>
#include <h5xx/error.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>


namespace h5xx {


/**
 * create dataset from a std::vector of fundamental type
 **/
template <typename h5xxObject, typename T, typename StoragePolicy>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, dataset>::type
create_dataset(h5xxObject const& object, std::string const& name, T const& value,
                   StoragePolicy const& storage_policy = StoragePolicy())
{
    typedef typename T::value_type value_type;
    hid_t type_id = ctype<value_type>::hid(); // this ID must not be closed
    std::vector<hsize_t> dims;
    dims.push_back(value.size());
    return dataset(object, name, type_id, dataspace(dims), storage_policy);
}

/**
 * create dataset from a std::vector of fundamental type, using default storage layout
 **/
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, dataset>::type
create_dataset(h5xxObject const& object, std::string const& name, T const& value)
{
    return create_dataset(object, name, value, h5xx::policy::storage::contiguous());
}



/**
 * write std::vector data to an existing dataset specified by location and name
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(h5xxObject const& object, std::string const& name, T const& value)
{
    dataset dset(object, name);
    write_dataset(dset, value);
}

/**
 * write std::vector data to dataset
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(dataset& dset, T const& value)
{
    typedef typename T::value_type value_type;
    hid_t type_id = ctype<value_type>::hid();
    dset.write(type_id, &*value.begin());
}

/**
 * Write std::vector data to an existing dataset specified by location and name,
 * memory and file locations (hyperslabs) are passed via the dataspace objects.
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(h5xxObject const& object, std::string const& name, T const& value,
                  dataspace const& memspace, dataspace const& filespace)
{
    dataset dset(object, name);
    write_dataset(dset, value, memspace, filespace);
}

/**
 * Write std::vector data, memory and file locations (hyperslabs) are passed via
 * the dataspace objects.
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(dataset& dset, T const& value, dataspace const& memspace, dataspace const& filespace)
{
    typedef typename T::value_type value_type;
    hid_t type_id = ctype<value_type>::hid();
    hid_t mem_space_id = memspace.hid(); //H5S_ALL;
    hid_t file_space_id = filespace.hid();
    hid_t xfer_plist_id = H5P_DEFAULT;
    dset.write(type_id, &*value.begin(), mem_space_id, file_space_id, xfer_plist_id);
}

/**
 * Write std::vector data to an existing dataset specified its location and
 * name, only the file location (hyperslab) is given via a slice object.
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(h5xxObject const& object, std::string const& name, T const& value, slice const& file_slice)
{
    dataset dset(object, name);
    write_dataset(dset, value, file_slice);
}

/**
 * Write std::vector data to an existing dataset, only the file location
 * (hyperslab) is given via a slice object.
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
write_dataset(dataset& dset, T const& value, slice const& file_slice)
{
    // --- create memory dataspace for the complete input array
    h5xx::dataspace memspace = h5xx::create_dataspace(value);
    // --- create file dataspace and select the slice (hyperslab) from it
    h5xx::dataspace filespace(dset);
    filespace.select(file_slice);
    write_dataset(dset, value, memspace, filespace);
}



/**
 * Read std::vector data from an existing dataset specified by location and name.
 * The vector data is resized and overwritten internally.
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
read_dataset(h5xxObject const& object, std::string const& name, T & value)
{
    dataset dset(object, name);
    read_dataset(dset, value);
}

/**
 * Read std::vector data from an existing dataset.
 * The vector data is resized and overwritten internally.
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
read_dataset(dataset & data_set, T & value)
{
    typedef typename T::value_type value_type;
    hid_t type_id = ctype<value_type>::hid();

    dataspace file_space(data_set);
    hssize_t nelem = file_space.get_select_npoints();
    value.clear();
    value.resize(nelem);

    hid_t mem_space_id = H5S_ALL;
    hid_t file_space_id = H5S_ALL;
    hid_t xfer_plist_id = H5P_DEFAULT;

    data_set.read(type_id, &*value.begin(), mem_space_id, file_space_id, xfer_plist_id);
}


/**
 * Read std::vector data from an existing dataset specified by location and name,
 * a slice specifies the data locations to be read in file space.  The vector
 * is not resized internally, the user must resize it in advance to fit the slice.
 */
template <typename h5xxObject, typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
read_dataset(h5xxObject const& object, std::string const& name, T & value, slice const& file_slice)
{
    dataset data_set(object, name);
    read_dataset(data_set, value, file_slice);
}

/**
 * Read std::vector data from an existing dataset, a slice specifies the data
 * locations to be read in file space.  The vector is not resized internally,
 * the user must resize it in advance to fit the slice.
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
read_dataset(dataset & data_set, T & value, slice const& file_slice)
{
    // --- create memory dataspace for the complete input array
    h5xx::dataspace memspace = h5xx::create_dataspace(value);
    // --- create file dataspace and select the slice (hyperslab) from it
    h5xx::dataspace filespace(data_set);
    filespace.select(file_slice);
    // ---
    read_dataset(data_set, value, memspace, filespace);
}

/**
 * Read std::vector data from an existing dataset, dataspace objects for both
 * memory and file allow to specify the locations of the data.  The vector is
 * not resized internally, the user must resize it in advance to fit the
 * dataspace.
 */
template <typename T>
inline typename boost::enable_if< boost::mpl::and_< is_vector<T>, boost::is_fundamental<typename T::value_type> >, void>::type
read_dataset(dataset & data_set, T & value, dataspace const& memspace, dataspace const& filespace)
{
    // Assure that the vector has at least the capacity of the dataspace selection.
    // It is the user's responsibility to allocate memory outside (using value.resize()).
    if (static_cast<hsize_t>(filespace.get_select_npoints()) > value.capacity())
        H5XX_THROW("target vector does not provide enough space to store slice");

    typedef typename T::value_type value_type;
    hid_t type_id = ctype<value_type>::hid();

    hid_t mem_space_id = memspace.hid();
    hid_t file_space_id = filespace.hid();
    hid_t xfer_plist_id = H5P_DEFAULT;

    data_set.read(type_id, &*value.begin(), mem_space_id, file_space_id, xfer_plist_id);
}


} // namespace h5xx

#endif // ! H5XX_DATASET_STD_VECTOR
