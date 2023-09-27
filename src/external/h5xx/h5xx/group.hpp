/*
 * Copyright © 2008-2019 Felix Höfling
 * Copyright © 2018      Matthias Werner
 * Copyright © 2014-2015 Klaus Reuter
 * Copyright © 2013      Manuel Dibak
 * Copyright © 2008-2010 Peter Colberg
 * All rights reserved.
 *
 * This file is part of h5xx — a C++ wrapper for the HDF5 library.
 *
 * This software may be modified and distributed under the terms of the
 * 3-clause BSD license.  See accompanying file LICENSE for details.
 */

#ifndef H5XX_GROUP_HPP
#define H5XX_GROUP_HPP

#include <h5xx/file.hpp>
#include <h5xx/dataset.hpp>
#include <h5xx/utility.hpp>
#include <h5xx/error.hpp>
#include <h5xx/property.hpp>

namespace h5xx {

// forward declaration
template <typename h5xxObject>
class container;

// this class is meant to replace the H5::Group class
class group
{
public:
    /** default constructor */
    group() : hid_(-1) {}

    /** constructor to open or generate a group in an existing superior group */
    group(group const& other, std::string const& name);

    /** open root group of file */
    group(file const& f);

    /**
     * deleted copy constructor
     *
     * Calling the constructor throws an exception. Copying must be elided by
     * return value optimisation. See also "group h5xx::move(group&)".
     */
    group(group const& other);

    /**
     * assignment operator
     *
     * Uses the copy-and-swap idiom. Move semantics is achieved in conjunction
     * with "group h5xx::move(group&)", i.e., the group object on the right
     * hand side is empty after move assignment.
     */
    group & operator=(group other);

    /** default destructor */
    ~group();

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

    /** open handle to HDF5 group from an object's ID (called by non-default constructors) */
    void open(group const& other, std::string const& name);

    /** close handle to HDF5 group (called by default destructor) */
    void close();

    /** methods to yield container adapters */
    container<h5xx::dataset> datasets() const;
    container<h5xx::group> groups() const;

private:
    /** HDF5 object ID */
    hid_t hid_;

    template <typename h5xxObject>
    friend void swap(h5xxObject& left, h5xxObject& right);
}; // class group

/**
 * iterator template class
 */
template <typename T, bool is_const>
class group_iterator
  : public std::iterator<std::forward_iterator_tag
     , T
     , std::ptrdiff_t
     , typename std::conditional<is_const, T const*, T*>::type
     , typename std::conditional<is_const, T const&, T&>::type
    >
{
public:
    group_iterator() noexcept;
    ~group_iterator() noexcept;

    /** construct iterator over a group */
    group_iterator(group const&) noexcept;

    /** copy constructor, don't copy resource pointer */
    group_iterator(group_iterator const&) noexcept;
//    group_iterator(group_iterator<T, true> const&) noexcept;    // FIXME convert const_iterator to iterator

    /** move constructor */
    group_iterator(group_iterator&&) = default;                 // uses move assignment

    /** copy assignment */
    group_iterator& operator=(group_iterator const& other) = default;   // uses copy constructor

    /** move assignment, special handling of resource pointer */
    group_iterator& operator=(group_iterator&& other);

    /** pre- and post-increment operators */
    group_iterator& operator++();
    group_iterator operator++(int);

    /** returns h5xx-object, reference is to internal copy */
    typename std::conditional<is_const, T const&, T&>::type operator*();
    typename std::conditional<is_const, T const*, T*>::type operator->();

    /** comparison operators
     *
     *  Two iterators are equal if they refer to the same container (parent
     *  group) and if the referred elements have the same distance to the start
     *  of the sequence, i.e., they were reached by the same number of
     *  increments.
     */
    template <bool is_const2>
    bool operator==(group_iterator<T, is_const2> const&) const;

    template <bool is_const2>
    bool operator!=(group_iterator<T, is_const2> const&) const;

    /** return name of current element */
    std::string get_name()
    {
        if(stop_idx_ == 0)
            increment_();

        return name_;
    };

private:
    /** move forward by one step, calls H5Literate */
    bool increment_();

    /** pointer to parent group */
    group const* parent_ = nullptr;

    /**
     * index of element in group as used by H5Literate.
     * If stop_idx_ == -1U, iterator points past the end.
     * stop_idx_ == 0 indicates a freshly constructed iterator, which points at
     * the first element (non-empty group), or past the end (empty group).
     */
    hsize_t stop_idx_ = -1U;

    /** name of HDF5 element pointed to */
    std::string name_;

    /** instance of HDF5 element pointed to */
    T* element_ = nullptr;

    friend class group_iterator<T, !is_const>;

    friend typename container<T>::iterator container<T>::end() const noexcept;
    friend typename container<T>::const_iterator container<T>::cend() const noexcept;
}; // class group_iterator

/**
 * adapter class to convert a group into a container of HDF5 objects of given type
 *
 * Provides an iterator interface and can be used in a range-based loop.
 */
template <typename h5xxObject>
class container
{
public:
    typedef group_iterator<h5xxObject, false> iterator;
    typedef group_iterator<h5xxObject, true> const_iterator;
    container(group const&);

    iterator begin() const noexcept;
    iterator end() const noexcept;

    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

private:
    group const* parent_;
};

// forward declaration
inline bool exists_group(group const& grp, std::string const& name);

/*
 * implementation of h5xx::group
 */
inline group::group(group const& other, std::string const& name)
  : hid_(-1)
{
    open(other, name);
}

inline group::group(file const& f)
{
    hid_ = H5Gopen(f.hid(), "/", H5P_DEFAULT);
    if (hid_ < 0) {
        throw error("opening root group of file \"" + f.name() + "\"");
    }
}

inline group::group(group const& other)
  : hid_(other.hid_)
{
    // copying would be safe if the exception were disabled.
    throw error("h5xx::group can not be copied. Copying must be elided by return value optimisation.");
    H5Iinc_ref(hid_);
}

inline group& group::operator=(group other)
{
    swap(*this, other);
    return *this;
}

inline group::~group()
{
    close();
}

inline void group::open(group const& other, std::string const& name)
{
    if (hid_ >= 0) {
        throw error("h5xx::group object is already in use");
    }

    if (exists_group(other, name)) {
        hid_ = H5Gopen(other.hid(), name.c_str(), H5P_DEFAULT);
    }
    else {
        hid_t lcpl_id = H5Pcreate(H5P_LINK_CREATE);     // create group creation property list
        H5Pset_create_intermediate_group(lcpl_id, 1);   // set intermediate link creation
        hid_ = H5Gcreate(other.hid(), name.c_str(), lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (hid_ < 0){
        throw error("creating or opening group \"" + name + "\"");
    }
}

inline void group::close() {
    if (hid_ >= 0) {
        if(H5Gclose(hid_) < 0){
            throw error("closing h5xx::group with ID " + boost::lexical_cast<std::string>(hid_));
        }
        hid_ = -1;
    }
}

/**
 * return true if group "name" exists in group "grp"
 */
inline bool exists_group(group const& grp, std::string const& name)
{
    hid_t hid = grp.hid();
    H5E_BEGIN_TRY {
        hid = H5Gopen(hid, name.c_str(), H5P_DEFAULT);
        if (hid > 0) {
            H5Gclose(hid);
        }
    } H5E_END_TRY
    return (hid > 0);
}

inline hid_t open_group(hid_t loc_id, std::string const& path)
{
    hid_t group_id;
    H5E_BEGIN_TRY {
        group_id = H5Gopen(loc_id, path.c_str(), H5P_DEFAULT);
    } H5E_END_TRY
    return group_id;
}

inline container<h5xx::dataset> group::datasets() const
{
    return container<h5xx::dataset>(*this);
}

inline container<h5xx::group> group::groups() const
{
    return container<h5xx::group>(*this);
}

/*
 * implementation of adapter classes for group containers
 */
template <typename h5xxObject>
inline container<h5xxObject>::container(group const& grp)
  : parent_(&grp)
{}

template <typename h5xxObject>
inline typename container<h5xxObject>::iterator container<h5xxObject>::begin() const noexcept
{
    return iterator(*parent_);
}

template <typename h5xxObject>
inline typename container<h5xxObject>::const_iterator container<h5xxObject>::cbegin() const noexcept
{
    return const_iterator(*parent_);
}

template <typename h5xxObject>
inline typename container<h5xxObject>::iterator container<h5xxObject>::end() const noexcept
{
    iterator iter(*parent_);
    iter.stop_idx_ = -1U;
    return iter;
}

template <typename h5xxObject>
inline typename container<h5xxObject>::const_iterator container<h5xxObject>::cend() const noexcept
{
    const_iterator iter(*parent_);
    iter.stop_idx_ = -1U;
    return iter;
}

/*
 * free functions begin(), end(), cbegin(), cend() on container<>
 */
template <typename h5xxObject>
inline typename container<h5xxObject>::iterator begin(container<h5xxObject> const& c) noexcept
{
    return c.begin();
}

template <typename h5xxObject>
inline typename container<h5xxObject>::const_iterator cbegin(container<h5xxObject> const& c) noexcept
{
    return c.cbegin();
}

template <typename h5xxObject>
inline typename container<h5xxObject>::iterator end(container<h5xxObject> const& c) noexcept
{
    return c.end();
}

template <typename h5xxObject>
inline typename container<h5xxObject>::const_iterator cend(container<h5xxObject> const& c) noexcept
{
    return c.cend();
}

/*
 * implementation of group_iterator
 */
template <typename T, bool is_const>
inline group_iterator<T, is_const>::group_iterator() noexcept
{}

template <typename T, bool is_const>
inline group_iterator<T, is_const>::group_iterator(group const& group) noexcept
  : parent_(&group)
  , stop_idx_(0)
{}

template <typename T, bool is_const>
inline group_iterator<T, is_const>::group_iterator(group_iterator const& other) noexcept
  : parent_(other.parent_)
  , stop_idx_(other.stop_idx_)
  , name_(other.name_)
  , element_(nullptr)   // don't copy pointer to HDF5 resource
{}

template <typename T, bool is_const>
inline group_iterator<T, is_const>::~group_iterator() noexcept
{
    if (element_) {
        delete element_;
    }
}

template <typename T, bool is_const>
inline group_iterator<T, is_const>& group_iterator<T, is_const>::operator=(group_iterator&& other)
{
    if (this != &other) {
        if (element_) {
            delete element_;                // release HDF5 resource before assignment
        }

        parent_ = other.parent_;
        stop_idx_ = other.stop_idx_;
        name_ = std::move(other.name_);
        element_ = other.element_;

        // leave behind default constructed object,
        // destructor must not release *other.element_
        other.parent_ = nullptr;
        other.stop_idx_ = -1U;
        other.name_ = std::string();
        other.element_ = nullptr;
    }

    return *this;
}

template <typename T, bool is_const> template <bool is_const2>
inline bool group_iterator<T, is_const>::operator==(group_iterator<T, is_const2> const& other) const
{
    if (parent_ != other.parent_) {
        throw h5xx::error("comparing iterators over different groups");
    }

    // ensure that both iterators are fully initialised before the comparison.
    // reason: stop_idx = 0 and 1 (non-empty parent) or -1 (empty parent) are equivalent.
    //
    // FIXME it seems better to have no side effects here, any idea?
    // if the group is non-empty, the following yield true as well: 0 == 1, 1 == 0.
    // for an empty group, 0 == -1 and -1 == 0 hold instead.
    if (stop_idx_ == 0) {
        const_cast<group_iterator*>(this)->increment_();
    }

    if (other.stop_idx_ == 0) {
        const_cast<group_iterator<T, is_const2>*>(&other)->increment_();
    }

    return stop_idx_ == other.stop_idx_;
}

template <typename T, bool is_const> template <bool is_const2>
inline bool group_iterator<T, is_const>::operator!=(group_iterator<T, is_const2> const& other) const
{
    return !(*this == other);
}

template <typename T, bool is_const>
inline group_iterator<T, is_const>& group_iterator<T, is_const>::operator++() // ++it
{
    if (element_) {
        delete element_;
    }
    element_ = nullptr;

    if (!parent_) {
        throw error("cannot increment default constructed h5xx::group_iterator");
    }

    // if needed, fully initialise iterator
    if(stop_idx_ == 0) {
        increment_();
    }

    // perform actual increment
    if(!increment_()) {  // out of range
        stop_idx_ = -1U;
    }

    return *this;
}

template <typename T, bool is_const>
inline group_iterator<T, is_const> group_iterator<T, is_const>::operator++(int) // it++
{
    group_iterator<T, is_const> tmp(*this);
    std::swap(tmp.element_, element_);      // pass HDF5 resource over to 'tmp'
    ++(*this);
    return tmp;
}

template <typename T, bool is_const>
inline typename std::conditional<is_const, T const&, T&>::type group_iterator<T, is_const>::operator*()
{
    if (!parent_) {
        throw error("cannot dereference default constructed h5xx::group_iterator");
    }

    if (stop_idx_ == 0) {       // initialise iterator freshly returned by begin()
        increment_();
    }

    if (stop_idx_ == -1U) {     // iterator is past the end, throw std::out_of_range
        std::string error_msg = "parent group";
        if (parent_ && parent_->valid()) {
            error_msg += " " + h5xx::get_name(*parent_);
        }
        else {
            error_msg = "non-existing " + error_msg;
        }
        throw std::out_of_range(error_msg);
    }

    if (!element_) {
        element_ = new T(*parent_, name_);
    }
    return *element_;
}

template <typename T, bool is_const>
inline typename std::conditional<is_const, T const*, T*>::type group_iterator<T, is_const>::operator->()
{
    return &**this;
}

namespace detail {

// forward declaration of helper function
template <typename T>
herr_t find_name_of_type(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data);

} // namespace detail

template <typename T, bool is_const>
inline bool group_iterator<T, is_const>::increment_()
{
    // if parent_ is not a valid group, set iterator past the end and return false
    if(!parent_ || !parent_->valid()) {
        stop_idx_ = -1U;
        return false;
    }

    herr_t retval = H5Literate(
        parent_->hid(), H5_INDEX_NAME, H5_ITER_INC, &stop_idx_
      , detail::find_name_of_type<T>, &name_
    );

    // evaluate return code
    if(retval == 0) {       // no element of type T was found
        stop_idx_ = -1U;    // set iterator to past-the-end iterator
        name_ = std::string();
    }

    return retval > 0;      // everything went fine
}

namespace detail {

/**
 * determine whether a given HDF5 object has a given type and return its name
 *
 * @return code: success: > 0, wrong type: 0, error: < 0
 */
template <H5O_type_t type>
herr_t find_name_of_type_impl(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data)
{
    H5O_info_t obj_info;

    /** returns non-negative upon success, negative if failed */
    herr_t retval = H5Oget_info_by_name(g_id, name, &obj_info, H5P_DEFAULT);

    /** check retval */
    if(retval >= 0) {
        /** filter for given HDF5 type */
        if(obj_info.type == type) {
            std::string* str_ptr = reinterpret_cast<std::string*>(op_data);
            *str_ptr = name;
            retval++;   // ensure retval is > 0 for short-circuit success
        }
        else {          // if element name was not 'type', continue iteration
            retval = 0;
        }
    }
    else {
        throw error("Cannot get object info of "+std::string(name));
    }

    return retval;
}

template <>
herr_t find_name_of_type<group>(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data)
{
    return find_name_of_type_impl<H5O_TYPE_GROUP>(g_id, name, info, op_data);
}

template <>
herr_t find_name_of_type<group const>(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data)
{
    return find_name_of_type_impl<H5O_TYPE_GROUP>(g_id, name, info, op_data);
}

template <>
herr_t find_name_of_type<dataset>(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data)
{
    return find_name_of_type_impl<H5O_TYPE_DATASET>(g_id, name, info, op_data);
}

template <>
herr_t find_name_of_type<dataset const>(hid_t g_id, char const* name, H5L_info_t const* info, void* op_data)
{
    return find_name_of_type_impl<H5O_TYPE_DATASET>(g_id, name, info, op_data);
}

} // namespace detail
} // namespace h5xx

#endif /* ! H5XX_GROUP_HPP */
