/*
 * Copyright © 2014      Felix Höfling
 * Copyright © 2018      Matthias Werner
 * Copyright © 2014      Manuel Dibak
 * Copyright © 2014-2016 Klaus Reuter
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_DATASPACE_DATASPACE_HPP
#define H5XX_DATASPACE_DATASPACE_HPP

#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/lexical_cast.hpp>

#include <h5xx/hdf5_compat.hpp>
#include <h5xx/error.hpp>
#include <h5xx/utility.hpp>
#include <h5xx/slice.hpp>

namespace h5xx {

// forward declarations
class attribute;
class dataset;

/**
 * Represents an HDF5 dataspace.
 *
 */
class dataspace
{
public:
    /** default constructor */
    dataspace() : hid_(-1) {}

    /** construct dataspace of rank zero */
    dataspace(H5S_class_t type);

    /** construct simple dataspace of given extents, boost::array version */
    template <std::size_t N>
    dataspace(boost::array<std::size_t, N> const& dims);
    template <std::size_t N>
    dataspace(boost::array<hsize_t, N> const& dims);

    /** construct simple dataspace of given extents, std::vector version */
    dataspace(std::vector<std::size_t> const& dims);
    dataspace(std::vector<hsize_t> const& dims);

    /** construct simple dataspace of given extents and maximal extents, boost::array version */
    template <std::size_t N>
    dataspace(boost::array<size_t, N> const& dims, boost::array<size_t, N> const& max_dims);
    template <std::size_t N>
    dataspace(boost::array<hsize_t, N> const& dims, boost::array<hsize_t, N> const& max_dims);

    /** construct simple dataspace of given extents and maximal extents, std::vector version */
    dataspace(std::vector<size_t> const& dims, std::vector<size_t> const& max_dims);
    dataspace(std::vector<hsize_t> const& dims, std::vector<hsize_t> const& max_dims);

// --- generic array versions, cannot be resolved by the compiler currently ---
//    /** construct simple dataspace of given extents */
//    template <class ArrayType>
//    dataspace(ArrayType const& dims);
//
//    /** construct simple dataspace of given extents and maximal extents */
//    template <class ArrayType>
//    dataspace(ArrayType const& dims, ArrayType const& max_dims);

    /**
     * deleted copy constructor
     *
     * Calling the constructor throws an exception. Copying must be elided by
     * return value optimisation. See also "dataspace h5xx::move(dataspace&)".
     */
    dataspace(dataspace const& other);

    /**
     * assignment operator
     *
     * Uses the copy-and-swap idiom. Move semantics is achieved in conjunction
     * with "dataspace h5xx::move(dataspace&)", i.e., the dataspace object on the right
     * hand side is empty after move assignment.
     */
    dataspace & operator=(dataspace other);

    /** default destructor */
    ~dataspace();

    /** return HDF5 object ID */
    hid_t hid() const
    {
        return hid_;
    }

    /** returns true if associated to a valid HDF5 object */
    bool valid() const
    {
        return hid_ >= 0;
    }

    /** returns rank/dimensionality of simple dataspace */
    unsigned int rank() const;

    /** returns extents/dimensions of simple dataspace, optionally the maximum
    *   dimensions are returned in maxdims
    */
//    template <std::size_t N>
//    boost::array<std::size_t, N> extents(hsize_t *maxdims = NULL) const;
    template <std::size_t N>
    boost::array<hsize_t, N> extents(hsize_t *maxdims = NULL) const;

    /** returns extents/dimensions of simple dataspace, optionally the maximum
    *   dimensions are returned in maxdims
    */
//    std::vector<std::size_t> extents(hsize_t *maxdims = NULL) const;
    std::vector<hsize_t> extents(hsize_t *maxdims = NULL) const;

    /** returns true if dataspace is of scalar type */
    bool is_scalar() const;

    /** returns true if dataspace is of simple type */
    bool is_simple() const;

    /**
     * flags for slice (hyperslab) selection
     */
    enum mode
    {
        SET = H5S_SELECT_SET
      , OR = H5S_SELECT_OR
      , AND = H5S_SELECT_AND
      , XOR = H5S_SELECT_XOR
      , NOTB = H5S_SELECT_NOTB
      , NOTA = H5S_SELECT_NOTA
    };

    /**
     * slice/hyperslab selection interface
     */
    void select(slice const& _slice, int mode = SET);

    /**
     * return the number of elements currently selected from the dataspace
     */
    hssize_t get_select_npoints() const;

private:
    /** HDF5 object ID */
    hid_t hid_;

    template <typename h5xxObject>
    friend void swap(h5xxObject& left, h5xxObject& right);

    friend class h5xx::attribute; // method "operator dataspace()"
    friend class h5xx::dataset;   // method "operator dataspace()"
};

inline dataspace::dataspace(dataspace const& other)
  : hid_(other.hid_)
{
    // copying would be safe if the exception were disabled.
    throw error("h5xx::dataspace can not be copied. Copying must be elided by return value optimisation.");
    H5Iinc_ref(hid_);
}

inline dataspace & dataspace::operator=(dataspace other)
{
    swap(*this, other);
    return *this;
}

inline dataspace::dataspace(H5S_class_t type)
{
    if((hid_ = H5Screate(type)) < 0) {
        throw error("creating dataspace");
    }
}

inline dataspace::dataspace(std::vector<std::size_t> const& dims)
{
    std::vector<hsize_t> h_dims = to_hsize_t(dims);
    if ((hid_ = H5Screate_simple(h_dims.size(), &*h_dims.begin(), NULL)) < 0) {
        throw error("creating simple dataspace");
    }
}

inline dataspace::dataspace(std::vector<hsize_t> const& dims)
{
    if ((hid_ = H5Screate_simple(dims.size(), &*dims.begin(), NULL)) < 0) {
        throw error("creating simple dataspace");
    }
}

template <std::size_t N>
dataspace::dataspace(boost::array<std::size_t, N> const& dims)
{
    std::vector<hsize_t> h_dims = to_hsize_t(dims);
    if ((hid_ = H5Screate_simple(N, &*h_dims.begin(), NULL)) < 0) {
        throw error("creating simple dataspace");
    }
}

template <std::size_t N>
dataspace::dataspace(boost::array<hsize_t, N> const& dims)
{
    if ((hid_ = H5Screate_simple(N, &*dims.begin(), NULL)) < 0) {
        throw error("creating simple dataspace");
    }
}

inline dataspace::dataspace(std::vector<size_t> const& dims, std::vector<size_t> const& max_dims)
{
    std::vector<hsize_t> h_dims = to_hsize_t(dims);
    std::vector<hsize_t> h_max_dims = to_hsize_t(max_dims);
    if ((hid_ = H5Screate_simple(h_dims.size(), &*h_dims.begin(), &*h_max_dims.begin())) < 0) {
        throw error("creating simple dataspace");
    }
}

inline dataspace::dataspace(std::vector<hsize_t> const& dims, std::vector<hsize_t> const& max_dims)
{
    if ((hid_ = H5Screate_simple(dims.size(), &*dims.begin(), &*max_dims.begin())) < 0) {
        throw error("creating simple dataspace");
    }
}

template <std::size_t N>
dataspace::dataspace(boost::array<size_t, N> const& dims, boost::array<size_t, N> const& max_dims)
{
    std::vector<hsize_t> h_dims = to_hsize_t(dims);
    std::vector<hsize_t> h_max_dims = to_hsize_t(max_dims);
    if ((hid_ = H5Screate_simple(N, &*h_dims.begin(), &*h_max_dims.begin())) < 0) {
        throw error("creating simple dataspace");
    }
}

template <std::size_t N>
dataspace::dataspace(boost::array<hsize_t, N> const& dims, boost::array<hsize_t, N> const& max_dims)
{
    if ((hid_ = H5Screate_simple(N, &*dims.begin(), &*max_dims.begin())) < 0) {
        throw error("creating simple dataspace");
    }
}

inline dataspace::~dataspace()
{
    if (hid_ >= 0) {
        if(H5Sclose(hid_) < 0)
            throw error("closing h5xx::dataspace with ID " + boost::lexical_cast<std::string>(hid_));
        hid_ = -1;
    }
}

inline unsigned int dataspace::rank() const
{
    if (!valid()) {
        throw error("invalid dataspace");
    }
    int rank = H5Sget_simple_extent_ndims(hid_);
    if (rank < 0) {
        throw error("dataspace has invalid rank");
    }
    return rank;
}

template <std::size_t N>
boost::array<hsize_t, N> dataspace::extents(hsize_t *maxdims) const
{
    boost::array<hsize_t, N> h_dims;
    if (rank() != N) {
        throw error("mismatching dataspace rank");
    }
    if (H5Sget_simple_extent_dims(hid_, &*h_dims.begin(), maxdims) < 0) {
        throw error("determining extents");
    }
    return h_dims;
}
//template <std::size_t N>
//boost::array<size_t, N> dataspace::extents(hsize_t *maxdims) const
//{
//    boost::array<hsize_t, N> h_dims;
//    h_dims = dataspace::extents(maxdims);
//    boost::array<size_t, N> dims = to_size_t(h_dims);
//    return dims;
//}


inline std::vector<hsize_t> dataspace::extents(hsize_t *maxdims) const
{
    std::vector<hsize_t> h_dims(rank());
    if (H5Sget_simple_extent_dims(hid_, &*h_dims.begin(), maxdims) < 0) {
        throw error("determining extents");
    }
    return h_dims;
}
//std::vector<size_t> dataspace::extents(hsize_t *maxdims) const
//{
//    std::vector<hsize_t> h_dims(rank());
//    h_dims = dataspace::extents(maxdims);
//    std::vector<size_t> dims = to_size_t(h_dims);
//    return dims;
//}


inline bool dataspace::is_scalar() const
{
    if (!valid()) {
        return false;
    }
    return H5Sget_simple_extent_type(hid_) == H5S_SCALAR;
}

inline bool dataspace::is_simple() const
{
    if (!valid()) {
        return false;
    }
    return H5Sget_simple_extent_type(hid_) == H5S_SIMPLE;
}

inline void dataspace::select(slice const& s, int mode)
{
    if (!valid()) {
        throw error("invalid dataspace");
    }

    // to interpret the slicing string inside "slice" we need the dataspace's extents
    // slice is kept constant, while we get a local copy of count, which is clipped
    // according to extents
    std::vector<hsize_t> count = s.get_count(extents());

    if (s.rank() != rank()) {
        throw error("dataspace and slice have mismatching rank");
    }
    if (H5Sselect_hyperslab(
            hid_
          , static_cast<H5S_seloper_t>(mode)
          , &*s.get_offset().begin()
          , s.get_stride().size() > 0 ? &*s.get_stride().begin() : NULL
          , &*count.begin()
          , s.get_block().size() > 0 ? &*s.get_block().begin() : NULL
        ) < 0) {
        throw error("H5Sselect_hyperslab");
    }
}

inline hssize_t dataspace::get_select_npoints() const
{
    if (!valid()) {
        throw error("invalid dataspace");
    }
    return H5Sget_select_npoints(hid_);
}

} // namespace h5xx

#endif /* ! H5XX_DATASPACE_DATASPACE_HPP */
