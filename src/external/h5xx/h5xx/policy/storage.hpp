/*
 * Copyright © 2014-2015 Felix Höfling
 * Copyright © 2014-2015 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_POLICY_STORAGE_HPP
#define H5XX_POLICY_STORAGE_HPP

#include <algorithm>
#include <iterator>
#include <vector>

#include <cstring>

// --- we need a smart pointer to use std::vector as container for the filter pipelines and modifier sets
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <h5xx/error.hpp>
#include <h5xx/h5xx.hpp>
#include <h5xx/policy/filter.hpp>


namespace h5xx {
namespace policy {
namespace storage {



// --- Implementations of HDF5 storage modifier policies (fill_value, track_times) ---

/**
 * Abtract base class for storage modifier policies: Defines the interface and
 * is used to specify a modifier set (in analogy to the filter pipeline).
 */
class storage_modifier_base
{
public:
    virtual void set_storage(hid_t) const = 0;
};

/**
 * storage modifier policy class to define a fill value
 */
class fill_value
  : public storage_modifier_base
{
public:
    template <typename T>
    fill_value(T value_, bool optional_ = false)
    {
        type = ctype<T>::hid();
        ssize_t nbytes = sizeof(T);
        value.resize(nbytes);
        memcpy(&*value.begin(), &value_, nbytes);
        optional = optional_;
    }

    /** set fill_value for given property list */
    virtual void set_storage(hid_t plist) const
    {
        // Preferred a call to H5Pset_fill_value() instead of H5Pset_filter()
        // because of the internal checking it does, cf H5Pdcpl.c
        // TODO : Add H5Pset_fill_time() here, if desired.
        if (H5Pset_fill_value(plist, type, &*value.begin()) < 0)
            if (!optional)
                throw error("setting fill_value failed");
    }
private:
    hid_t type;
    std::vector<char> value;
    bool optional;
};

/**
 * policy class to enable the track times feature
 */
class track_times
  : public storage_modifier_base
{
public:
    track_times(bool optional_ = false)
      : optional(optional_)
    {}

    /** activate track_times for given property list */
    virtual void set_storage(hid_t plist) const
    {
        if (H5Pset_obj_track_times( plist, (hbool_t)true ) < 0)
            if (!optional)
                throw error("setting track_times failed");
    }

private:
    bool optional;
};





// --- Implementations of HDF5 storage layout policies (contiguous, compact, chunked) ---

/**
 * base class for storage policies, provides handling of storage modifier policies
 */
class storage_layout_base
{
public:
    typedef h5xx::policy::storage::storage_modifier_base modifier_base_t;
    typedef std::vector<boost::shared_ptr<modifier_base_t> > modifier_set_t;

    /**
     * activate the modifier pipeline for a given dataset creation property list
     */
    void set_storage_modifiers(hid_t plist) const
    {
        typename modifier_set_t::const_iterator m;
        for (m = modifier_.begin(); m != modifier_.end(); ++m) {
            (*m)->set_storage(plist); // dereference operator required due to shared_ptr usage
        }
    }

    // --- make this class abstract: implementations exist in contiguous(), compact(), chunked()
    virtual void set_storage(hid_t) const = 0;

protected:
    // vector to hold the set of modifiers we want to have applied
    modifier_set_t modifier_;
};


/**
 * policy class to specify a contiguous dataset layout, optionally
 * along with a modifier set (e.g. fill_value, track_times).
 */
class contiguous
        : public storage_layout_base
{
public:
    contiguous() {}

    /**
     * add a modifier to the modifier set
     */
    template <typename ModifierType>
    contiguous set(ModifierType modifier)
    {
        modifier_.push_back( boost::make_shared<ModifierType> (modifier) );
        return *this;
    }


    /** set compact storage for given property list */
    void set_storage(hid_t plist) const
    {
        if (H5Pset_layout(plist, H5D_CONTIGUOUS) < 0) {
            throw error("setting contiguous dataset layout failed");
        }
        // add storage modifiers to property list
        set_storage_modifiers(plist);
    }
};

/**
 * policy class to specify a compact dataset layout, optionally
 * along with a modifier set (e.g. fill_value, track_times).
 */
class compact
        : public storage_layout_base
{
public:
    compact() {}

    /**
     * add a modifier to the modifier set
     */
    template <typename ModifierType>
    compact set(ModifierType modifier)
    {
        modifier_.push_back( boost::make_shared<ModifierType> (modifier) );
        return *this;
    }

    /** set compact storage for given property list */
    void set_storage(hid_t plist) const
    {
        if (H5Pset_layout(plist, H5D_COMPACT) < 0) {
            throw error("setting compact dataset layout failed");
        }
        // add storage modifiers to property list
        set_storage_modifiers(plist);
    }
};

/**
 * policy class to specify a chunked dataset layout of given size, optionally
 * along with a filter pipeline (for, e.g., data compression), optionally
 * along with a modifier set (e.g. fill_value, track_times).
 */
class chunked
        : public storage_layout_base
{
public:
    typedef h5xx::policy::filter::filter_base filter_base_t;
    typedef std::vector<boost::shared_ptr<filter_base_t> > filter_pipeline_t;

    /**
     * Specify the size, in dataset elements, of a chunk in each dimension.
     * The number of dimensions must equal the rank of the dataset.
     */
    chunked(int ndims, hsize_t const* dims)
    {
        std::copy(dims, dims + ndims, std::back_inserter(dims_));
    }

    template <typename ContainerType>
    chunked(ContainerType const& dims)
    {
        std::copy(dims.begin(), dims.end(), std::back_inserter(dims_));
    }

    /** set chunked storage layout for given property list */
    void set_storage(hid_t plist) const
    {
        bool err = false;
        err |= H5Pset_layout(plist, H5D_CHUNKED) < 0;
        err |= H5Pset_chunk(plist, dims_.size(), &*dims_.begin()) < 0;
        if (err) {
            throw error("setting chunked dataset layout failed");
        }

        // add storage modifiers to property list
        set_storage_modifiers(plist);

        // add filter pipeline to property list
        typename filter_pipeline_t::const_iterator f;
        for (f = filter_.begin(); f != filter_.end(); ++f) {
            (*f)->set_filter(plist);  // dereference operator required due to shared_ptr usage
        }
    }

    /**
     * add a filter to the filter pipeline
     */
    template <typename FilterType>
    chunked add(FilterType filter)
    {
        filter_.push_back( boost::make_shared<FilterType> (filter) );
        return *this;
    }

    /**
     * add a modifier to the modifier set
     */
    template <typename ModifierType>
    chunked set(ModifierType modifier)
    {
        modifier_.push_back( boost::make_shared<ModifierType> (modifier) );
        return *this;
    }

private:
    // chunk dimensions
    std::vector<hsize_t> dims_;
    // filter pipeline
    filter_pipeline_t filter_;
};

} //namespace storage
} //namespace policy
} //namespace h5xx

#endif // ! H5XX_POLICY_STORAGE_HPP
