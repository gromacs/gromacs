/*
 * Copyright © 2014-2018 Felix Höfling
 * Copyright © 2018      Matthias Werner
 * Copyright © 2014-2015 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_SLICE_HPP
#define H5XX_SLICE_HPP

#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include <sstream>

// --- C regular expression support
#include <sys/types.h>
#include <regex.h>
#include <assert.h>

namespace h5xx {

class slice {
public:
    /**
     * Slice constructor accepting a string with a numpy-like slicing notation.
     */
    slice(std::string const& slice_str);

//   --- constructor with zero offset by default?
//    template <class ArrayType>
//    slice(ArrayType count);

    template <class ArrayType>
    slice(ArrayType offset, ArrayType count);

    template <class ArrayType>
    slice(ArrayType offset, ArrayType count, ArrayType stride);

    template <class ArrayType>
    slice(ArrayType offset, ArrayType count, ArrayType stride, ArrayType block);

    size_t rank() const;

    std::vector<hsize_t> const& get_count() const;
    std::vector<hsize_t> const& get_offset() const;
    std::vector<hsize_t> const& get_stride() const;
    std::vector<hsize_t> const& get_block() const;

    /**
     * Replace -1U in counts with 'extents - offset'. Returns the list of element counts.
     */
    std::vector<hsize_t> get_count(std::vector<hsize_t> const& extents = std::vector<hsize_t>()) const;
    std::vector<hsize_t> get_count_clipped(const std::vector<hsize_t> & extents) const;

private:
    std::vector<hsize_t> offset_, count_, stride_, block_;

    /**
     * Fill offset_, count_, stride_, block_ arrays based on the slice string.
     * When undetermined ranges such as ":" are used, count is set to -1U.
     */
    void parse_string_(std::string const& slice_str);

    // regular expression stuff, used inside the parse_string_() method
    regex_t slicing_fmt_rex_;               // "1:2:3,4,:"
    regex_t colon_rex_;                     // ":"
    regex_t slice_single_rex_;              // "1"
    regex_t slice_range_rex_;               // "1:4"
    regex_t slice_range_stride_rex_;        // "1:4:2"
    regex_t slice_full_range_rex_;          // ":"
    regex_t slice_full_range_ceil_rex_;     // ":2"
    regex_t slice_full_range_floor_rex_;    // "2:"
    regex_t slice_full_range_stride_rex_;   // "::2"

    // auxiliary functions for regular expression handling
    void init_rex_();
    void free_rex_();
    void compile_rex_(regex_t & rex, std::string const& str);
    bool match_rex_(regex_t const& rex, std::string const& str);
};


inline slice::slice(std::string const& slice_str)
{
    // compile regular expressions,
    // TODO share with other instances of the slice class
    // (e.g., using static unique_ptr<regex_t>?)
    init_rex_();

    // parse string, be prepared for exceptions
    try {
        parse_string_(slice_str);
    }
    // release resources in any case
    catch (...) {
        free_rex_();
        throw;
    }
    free_rex_();
}

//template <class ArrayType>
//slice::slice(ArrayType count)
//{
//    std::copy(count.begin(), count.end(), std::back_inserter(count_));
//    offset_.assign(count.size(), 0);
//    stride_.clear();
//    block_.clear();
//}

template <class ArrayType>
slice::slice(ArrayType offset, ArrayType count)
{
    if (offset.size() != count.size()) {
        throw error("slice specification arrays must have identical length");
    }
    std::copy(offset.begin(), offset.end(), std::back_inserter(offset_));
    std::copy(count.begin(),  count.end(),  std::back_inserter(count_));
    stride_.clear();
    block_.clear();
}

template <class ArrayType>
slice::slice(ArrayType offset, ArrayType count, ArrayType stride)
{
    if ((offset.size() != count.size()) || (count.size() != stride.size())) {
        throw error("slice specification arrays must have identical length");
    }
    std::copy(offset.begin(), offset.end(), std::back_inserter(offset_));
    std::copy(count.begin(),  count.end(),  std::back_inserter(count_));
    std::copy(stride.begin(), stride.end(), std::back_inserter(stride_));
    block_.clear();
}

template <class ArrayType>
slice::slice(ArrayType offset, ArrayType count, ArrayType stride, ArrayType block)
{
    if ((offset.size() != count.size()) || (count.size() != stride.size()) || (stride.size() != block.size())) {
        throw error("slice specification arrays must have identical length");
    }
    std::copy(offset.begin(), offset.end(), std::back_inserter(offset_));
    std::copy(count.begin(),  count.end(),  std::back_inserter(count_));
    std::copy(stride.begin(), stride.end(), std::back_inserter(stride_));
    std::copy(block.begin(),  block.end(),  std::back_inserter(block_));
}

inline size_t slice::rank() const
{
    return count_.size();
}

inline std::vector<hsize_t> const& slice::get_offset() const
{
    return offset_;
}

inline std::vector<hsize_t> const& slice::get_count() const
{
    return count_;
}

inline std::vector<hsize_t> const& slice::get_stride() const
{
    return stride_;
}

inline std::vector<hsize_t> const& slice::get_block() const
{
    return block_;
}

inline std::vector<hsize_t> slice::get_count(std::vector<hsize_t> const& extents) const
{
    if (rank() != extents.size()) {
        throw error(std::string("mismatching rank of slice and extents in slice::get_count"));
    }

    // make a copy
    std::vector<hsize_t> count(count_);

    for (size_t i = 0; i < count.size(); i++) {
        if (count[i] == -1U) {
            count[i] = (extents[i] - offset_[i] - 1) / stride_[i] + 1;
        }
    }

    return count;
}

inline void slice::parse_string_(std::string const& slice_str)
{
    // check if slice_str_ contains a valid slicing notation
    if (!match_rex_(slicing_fmt_rex_, slice_str)) {
        throw error( std::string("array slicing format is invalid : ").append(slice_str) );
    }

    // create string vector with slice specification separated for each dimension
    std::vector<std::string> slice_specs = chop(slice_str, ",");

    // decrypt the slice specification, dimension by dimension
    typedef std::vector<std::string>::iterator vsit_t;
    size_t dim = 0;
    for (vsit_t it = slice_specs.begin(); it < slice_specs.end(); ++it)
    {
        std::string spec = *it;

        if (match_rex_(slice_single_rex_, spec)) {
            // "a"
            int i = str2num<int>(spec);
            offset_.push_back(i);
            count_.push_back(1);
            stride_.push_back(1);
        }
        else if (match_rex_(slice_range_rex_, spec)) {
            // "a:b"
            std::vector<std::string> range_nums = chop(spec, ":");
            int lo = str2num<int>( range_nums[0] );
            int hi = str2num<int>( range_nums[1] );
            offset_.push_back(lo);
            count_.push_back(hi - lo);
            stride_.push_back(1);
        }
        else if (match_rex_(slice_range_stride_rex_, spec)) {
            // "a:b:s"
            std::vector<std::string> range_nums = chop(spec, ":");
            int lo = str2num<int>( range_nums[0] );
            int hi = str2num<int>( range_nums[1] );
            int dx = str2num<int>( range_nums[2] );
            offset_.push_back(lo);
            count_.push_back((hi - lo - 1) / dx + 1);
            stride_.push_back(dx);
        }
        else if (match_rex_(slice_full_range_rex_, spec)) {
            // ":"
            offset_.push_back(0);
            count_.push_back(-1U);
            stride_.push_back(1);
        }
        else if (match_rex_(slice_full_range_ceil_rex_, spec)) {
            // ":b"
            std::vector<std::string> range_nums = chop(spec, ":");
            int hi = str2num<int>( range_nums[0] );
            offset_.push_back(0);
            count_.push_back(hi);
            stride_.push_back(1);
        }
        else if (match_rex_(slice_full_range_floor_rex_, spec)) {
            // "a:"
            std::vector<std::string> range_nums = chop(spec, ":");
            int lo = str2num<int>( range_nums[0] );
            offset_.push_back(lo);
            count_.push_back(-1U);
            stride_.push_back(1);
        }
        else if (match_rex_(slice_full_range_stride_rex_, spec)) {
            // "::s"
            std::vector<std::string> range_nums = chop(spec, ":");
            int dx = str2num<int>( range_nums[0] );
            offset_.push_back(0);
            count_.push_back(-1U);
            stride_.push_back(dx);
        }
        else {
            throw error( std::string("invalid slice specification : ").append(spec) );
        }
        ++dim;
    }
}

/**
 * Initialize regular expression patterns.
 */
inline void slice::init_rex_()
{
    compile_rex_(slicing_fmt_rex_,             "^([0-9:]+,)*[0-9:]+$");   // TODO make exact
    compile_rex_(colon_rex_,                   ":");
    compile_rex_(slice_single_rex_,            "^[0-9]+$");
    compile_rex_(slice_range_rex_,             "^[0-9]+:[0-9]+$");
    compile_rex_(slice_range_stride_rex_,      "^[0-9]+:[0-9]+:[0-9]+$");
    compile_rex_(slice_full_range_rex_,        "^:$");
    compile_rex_(slice_full_range_ceil_rex_,   "^:[0-9]+$");
    compile_rex_(slice_full_range_floor_rex_,  "^[0-9]+:$");
    compile_rex_(slice_full_range_stride_rex_, "^::[0-9]+$");
}

/**
 * compile regular expression handle. with error checking
 */
inline void slice::compile_rex_(regex_t & rex, std::string const& str)
{
    int ret = regcomp(&rex, str.c_str(), REG_EXTENDED);
    if (ret != 0) {
        const int errbuf_size = 256;
        char errbuf[errbuf_size];
        regerror(ret, &rex, errbuf, errbuf_size);
        throw error( std::string("regex compilation : ").append(errbuf) );
    }
}

/**
 * free regular expression handles.
 */
inline void slice::free_rex_()
{
    regfree(&slicing_fmt_rex_);
    regfree(&colon_rex_);
    regfree(&slice_single_rex_);
    regfree(&slice_range_rex_);
    regfree(&slice_range_stride_rex_);
    regfree(&slice_full_range_rex_);
    regfree(&slice_full_range_ceil_rex_);
    regfree(&slice_full_range_floor_rex_);
    regfree(&slice_full_range_stride_rex_);
}

/**
 * check if a string matches to a regular expression
 */
inline bool slice::match_rex_(regex_t const& rex, std::string const& str)
{
    int ret = regexec(&rex, str.c_str(), 0, NULL, 0);
    return (ret == 0);
}

inline std::vector<hsize_t> slice::get_count_clipped(const std::vector<hsize_t> & extents) const
{
    std::vector<hsize_t> clipped_count(count_);
    if (clipped_count.size() != extents.size())
        throw error(std::string("mismatch dimensionality of slice and extents in slice::get_count_clipped"));

    for (unsigned int i = 0; i < clipped_count.size(); i++) {
        if (clipped_count[i] == -1U) {
            clipped_count[i] = (extents[i] - offset_[i] - 1)/stride_[i] + 1;
        }
    }

    return clipped_count;
}

} // namespace h5xx


#endif // ! H5XX_SLICE_HPP
