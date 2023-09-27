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

#ifndef H5XX_POLICY_FILTER_HPP
#define H5XX_POLICY_FILTER_HPP

#include <boost/array.hpp>
#include <vector>

#include <h5xx/error.hpp>
#include <h5xx/ctype.hpp>
#include <h5xx/h5xx.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace h5xx {
namespace policy {
namespace filter {

/**
 * abtract base class for filter policies
 *
 * Defines the interface and is used to specify a filter pipeline.
 */
class filter_base
{
public:
    virtual void set_filter(hid_t) const = 0;
};

/**
 * policy class to enable gzip compression of a chunked dataset layout
 *
 * The compression level corresponds to the values of the GNU gzip tool.
 *
 * The filter is optional by default: if the filter result would be larger than
 * the input, then the compression filter returns failure and the uncompressed
 * data is stored in the file.
 */
class deflate
  : public filter_base
{
public:
    deflate(unsigned int level=6, bool optional=true)
      : flags_(optional ? H5Z_FLAG_OPTIONAL : 0), level_(level)
    {}

    /** set deflate filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        if (H5Pset_filter(plist, H5Z_FILTER_DEFLATE, flags_, 1, &level_) < 0) {
            throw error("setting data compression filter (gzip) failed");
        }
    }

private:
    // filter flags as a bit mask
    unsigned int flags_;
    // compression level
    unsigned int level_;
};

/**
 * policy class to enable SZIP compression of a chunked dataset layout
 *
 * The block size must be even and not greater than 32. This parameter affects
 * the compression ratio; the more the data values vary, the smaller this
 * number should be to achieve better performance. For optimal performance, it
 * is recommended that a chunk's fastest-changing dimension be equal to 128
 * times the block size.
 *
 * The filter is optional by default: if the filter result would be larger than
 * the input, then the compression filter returns failure and the uncompressed
 * data is stored in the file.
 */
class szip
  : public filter_base
{
public:
    enum coding_t {
        entropy                 /* entropy coding method: best suited for preprocessed data and small numbers   */
      , nearest_neighbour       /* nearest neighbor coding method: preprocess data, then apply entropy coding */
    };

    szip(unsigned int block_size=16, coding_t coding=nearest_neighbour, bool optional=true)
      : flags_(optional ? H5Z_FLAG_OPTIONAL : 0)
    {
        param_[0] = coding;
        param_[1] = block_size;
        if (block_size > 32 || block_size % 2) {
            throw error("SZIP filter: block size must be even and not greater than 32.");
        }
    }

    /** set szip filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        if (H5Pset_filter(plist, H5Z_FILTER_SZIP, flags_, 2, param_) < 0) {
            throw error("setting data compression filter (SZIP) failed");
        }
    }

private:
    // filter flags as a bit mask
    unsigned int flags_;
    // compression parameters
    unsigned int param_[2];
};

/**
 * policy class to set data shuffling filter for a chunked dataset layout
 */
class shuffle
  : public filter_base
{
public:
    shuffle(bool optional=false)
      : flags_(optional ? H5Z_FLAG_OPTIONAL : 0)
    {}

    /** set data shuffling filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        if (H5Pset_filter(plist, H5Z_FILTER_SHUFFLE, flags_, 0, NULL) < 0) {
            throw error("setting data shuffling filter failed");
        }
    }

private:
    // filter flags as a bit mask
    unsigned int flags_;
};

/**
 * policy class to enable Fletcher32 checksums of a chunked dataset layout
 *
 * The filter can not be made optional.
 */
class fletcher32
  : public filter_base
{
public:
    fletcher32() {}

    /** set fletcher32 filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        if (H5Pset_filter(plist, H5Z_FILTER_FLETCHER32, 0, 0, NULL) < 0) {
            throw error("setting Fletcher32 checksum filter failed");
        }
    }
};

/**
 * policy class to set the scaleoffset filter for a chunked dataset layout
 */
template <typename T>
class scaleoffset
  : public filter_base
{
public:

//    /** enable D-Scaling for floats, passing a scale factor is mandatory */
//    scaleoffset(int scale_factor_,
//                typename boost::enable_if<boost::is_floating_point<T> >::type* dummy = 0)
//        : scale_type(H5Z_SO_FLOAT_DSCALE), scale_factor(scale_factor_)
//    { T dog; }

//    /** enable scaling for integers, scale factor is optional (default: automatic determination) */
//    scaleoffset(int scale_factor_ = H5Z_SO_INT_MINBITS_DEFAULT,
//                typename boost::enable_if<boost::is_integral<T> >::type* dummy = 0)
//        : scale_type(H5Z_SO_INT), scale_factor(scale_factor_)
//    { T dog; }

    /** Temporary workaround, resolution errors pop up when using
     *  enable_if in the constructor implementations above. ??? */
    scaleoffset(int scale_factor_ = 0, bool optional = true)
      : flags_(optional ? H5Z_FLAG_OPTIONAL : 0)
    {
        if ((typeid(T) == typeid(double)) || (typeid(T) == typeid(float)))
            param_[0] = H5Z_SO_FLOAT_DSCALE;
        else if (typeid(T) == typeid(int))
            param_[0] = H5Z_SO_INT;
        else
            throw error("attempting to use an unsupported datatype with the scaleoffset filter");
        // ---
        param_[1] = (unsigned)scale_factor_;
    }

    /** set scaleoffset filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        // see the HDF5 source file H5Pdcpl.c for a reference on how to to the following call
        if ( H5Pset_filter(plist, H5Z_FILTER_SCALEOFFSET, flags_, (size_t)2, param_) < 0) {
            throw error("setting scaleoffset filter failed");
        }
    }
private:
//    H5Z_SO_scale_type_t scale_type;
//    int scale_factor;
    // filter flags as a bit mask
    unsigned int flags_;
    // compression parameters
    unsigned int param_[2];
};

/**
 * policy class to set the nbit filter for a chunked dataset layout
 */
class nbit
  : public filter_base
{
public:
    nbit(bool optional=false)
      : flags_(optional ? H5Z_FLAG_OPTIONAL : 0)
    {}

    /** set nbit filter for given property list */
    virtual void set_filter(hid_t plist) const
    {
        if (H5Pset_filter(plist, H5Z_FILTER_NBIT, flags_, 0, NULL) < 0) {
            throw error("setting nbit filter failed");
        }
    }

private:
    // filter flags as a bit mask
    unsigned int flags_;
};


} //namespace filter
} //namespace policy
} //namespace h5xx

#endif // ! H5XX_POLICY_FILTER_HPP
