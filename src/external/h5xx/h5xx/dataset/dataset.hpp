/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 * Copyright © 2014       Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_DATASET_DATASET_HPP
#define H5XX_DATASET_DATASET_HPP

#include <h5xx/dataset/utility.hpp>
#include <h5xx/dataspace.hpp>
#include <h5xx/datatype.hpp>
#include <h5xx/policy/storage.hpp>

namespace h5xx {

/**
 * Template class to wrap the HDF5 dataset.
 */
class dataset
{
public:
    dataset() : hid_(-1) {}

    /** open existing dataset */
    template <typename h5xxObject>
    dataset(h5xxObject const& object, std::string const& name, hid_t dapl_id = H5P_DEFAULT);

    /** create a new dataset */
    template <typename h5xxObject, typename StoragePolicy>
    dataset(
        h5xxObject const& object, std::string const& name
      , datatype const& dtype, dataspace const& dspace
      , StoragePolicy storage_policy = StoragePolicy()
      , hid_t lcpl_id = H5P_DEFAULT, hid_t dapl_id = H5P_DEFAULT
    );

    /** destructor, implicitly closes the dataset's hid_ */
    ~dataset();

    /**
     * deleted copy constructor
     *
     * Calling the constructor throws an exception. Copying must be elided by
     * return value optimisation. See also "dataset h5xx::move(dataset&)".
     */
    dataset(dataset const& other);

    /**
     * assignment operator
     *
     * Uses the copy-and-swap idiom. Move semantics is achieved in conjunction
     * with "dataset h5xx::move(dataset&)", i.e., the dataset object on the right
     * hand side is empty after move assignment.
     */
    dataset & operator=(dataset other);

    /** deduce dataspace from dataset */
    operator dataspace() const;

    /** return HDF5 object ID */
    hid_t hid() const;

    /** returns true if associated to a valid HDF5 object */
    bool valid() const;

    /** write value to the dataset */
    void write(hid_t type_id, void const* value, hid_t mem_space_id = H5S_ALL, hid_t file_space_id = H5S_ALL, hid_t xfer_plist_id = H5P_DEFAULT);

    /** read from the dataset into the buffer */
    void read(hid_t type_id, void* buffer, hid_t mem_space_id = H5S_ALL, hid_t file_space_id = H5S_ALL, hid_t xfer_plist_id = H5P_DEFAULT);

    /** return copy of dataset's type */
    hid_t get_type() const;

private:
    /** HDF5 handle of the dataset */
    hid_t hid_;

    template <typename h5xxObject>
    friend void swap(h5xxObject& left, h5xxObject& right);
};

template <typename h5xxObject>
dataset::dataset(h5xxObject const& object, std::string const& name, hid_t dapl_id)
  : hid_(-1)
{
    if (h5xx::exists_dataset(object, name))
    {
        hid_ = H5Dopen(object.hid(), name.c_str(), dapl_id);
    }
    if (hid_ < 0)
    {
        throw error("opening dataset \"" + name + "\" at HDF5 object \"" + get_name(object) + "\"");
    }
}

template <typename h5xxObject, typename StoragePolicy>
dataset::dataset(
    h5xxObject const& object, std::string const& name
  , datatype const& dtype, dataspace const& dspace
  , StoragePolicy storage_policy
  , hid_t lcpl_id, hid_t dapl_id
)
  : hid_(-1)
{
    if (h5xx::exists_dataset(object, name))
    {
        throw error("dataset \"" + name + "\" already exists");
    }

    bool err = false;

    // link creation property list
    if (lcpl_id != H5P_DEFAULT) {
        // create missing intermediate groups
        err |= H5Pset_create_intermediate_group(lcpl_id, 1) < 0;
        // TODO derive encoding of dataset name from string 'name'
        // H5Pset_char_encoding(lcpl_id, ...);
    }

    // dataset creation property list: evaluate policies
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    storage_policy.set_storage(dcpl_id);    // set storage and filters

    // create dataset
    err |= (hid_ = H5Dcreate(
        object.hid(),       // hid_t loc_id IN: Location identifier
        name.c_str(),       // const char *name IN: Dataset name
        dtype.get_type_id(), // hid_t dtype_id IN: Datatype identifier
        dspace.hid(),        // hid_t space_id IN: Dataspace identifier
        lcpl_id,            // hid_t lcpl_id IN: Link creation property list
        dcpl_id,            // hid_t dcpl_id IN: Dataset creation property list
        dapl_id             // dataset access property list: size of chunk cache
    )) < 0;
    err |= H5Pclose(dcpl_id) < 0;

    if (err)
    {
        throw error("creating dataset \"" + name + "\"");
    }
}

inline dataset::~dataset()
{
    if (hid_ >= 0) {
        if(H5Dclose(hid_) < 0){
            throw error("closing h5xx::dataset with ID " + boost::lexical_cast<std::string>(hid_));
        }
        hid_ = -1;
    }
}

inline dataset::dataset(dataset const& other)
  : hid_(other.hid_)
{
    // copying would be safe if the exception were disabled.
    throw error("h5xx::dataset can not be copied. Copying must be elided by return value optimisation.");
    H5Iinc_ref(hid_);
}

inline dataset & dataset::operator=(dataset other)
{
    swap(*this, other);
    return *this;
}

inline dataset::operator dataspace() const
{
    if (hid_ < 0) {
        throw error("retrieving dataspace from invalid dataset");
    }

    dataspace dspace;
    if((dspace.hid_ = H5Dget_space(hid_)) < 0) {      // we're friend
        throw error("dataset +\"" + get_name(*this) + "\" has invalid dataspace");
    }
    return dspace;
}

inline void dataset::write(hid_t type_id, void const* value, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id)
{
    if (H5Dwrite(hid_, type_id, mem_space_id, file_space_id, xfer_plist_id, value) < 0)
    {
        throw error("writing dataset");
    }
}

inline void dataset::read(hid_t type_id, void * buffer, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id)
{
    if (H5Dread(hid_, type_id, mem_space_id, file_space_id, xfer_plist_id, buffer) < 0)
    {
        throw error("reading dataset");
    }
}

inline hid_t dataset::get_type() const
{
    hid_t type_id = H5Dget_type(hid_);
    if (type_id < 0)
    {
        throw error("failed to obtain type_id of dataset \"" + get_name(*this) + "\"");
    }
    return type_id;
}

inline hid_t dataset::hid() const
{
    return hid_;
}

inline bool dataset::valid() const
{
    return (hid_ >= 0);
}
// --- END dataset class method implementations ---


/**
 * free function to create datasets
 */
template <typename h5xxObject, typename StoragePolicy>
dataset create_dataset(
    h5xxObject const& object
  , std::string const& name
  , datatype const& dtype
  , dataspace const& dspace
  , StoragePolicy storage_policy = StoragePolicy()
)
{
    return dataset(object, name, dtype, dspace, storage_policy);
}

template <typename h5xxObject>
dataset create_dataset(
    h5xxObject const& object
  , std::string const& name
  , datatype const& dtype
  , dataspace const& dspace
)
{
    return dataset(object, name, dtype, dspace, h5xx::policy::storage::contiguous());
}

} // namespace h5xx

#endif /* ! H5XX_DATASET_DATASET_HPP */
