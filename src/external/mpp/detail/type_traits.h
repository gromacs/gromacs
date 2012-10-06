/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
#pragma once

#include <vector>
#include <list>
#ifdef MPP_CXX11_ARRAY
#include <array>
#endif
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#ifdef MPP_CXX11_RVALREF
#include <memory>
#endif
#include <mpi.h>

#ifdef MPP_CXX11_ENUM_TRAIT
#define MPP_CXX11_TYPE_TRAITS
#endif

#ifdef MPP_CXX11_TYPE_TRAITS
#include <type_traits>
#endif

namespace mpi {

#ifdef MPP_CXX11_TYPE_TRAITS
//using stl true/false types
using std::true_type;
using std::false_type;
using std::enable_if;
#else
//define equivalent true/false types
struct false_type { enum { value = false }; };
struct true_type { enum { value = true }; };

//enable_if
template<bool B, class T = void>
struct enable_if {};
template<class T>
struct enable_if<true, T> { typedef T type; };
#endif

class data_layout;
//*****************************************************************************
//                                  MPI Type Traits
//*****************************************************************************
//The enable template parameter is unused for anything but conditional specialization
//(currently only used for enum)
template <class T, class enable=void>
struct mpi_type_traits {
    //return data layout for C++ type. Template specialization and
    //not function overloading is used to make sure no implicit type conversion
    //happens. Should not be called directly but instead the get_type_layout
    //function should be used.
    static inline data_layout get_layout(const T& raw);
    // is_static should be true if get_layout is the same for all objects of a class
    // and the extent of the mpi type returned by get_layout is equal to sizeof(T)
    // also it cannot contain any pointers.
    struct is_static : false_type {};
};

//base class used as shortcut for partial specialization
struct default_mpi_type_trait {
    struct is_static : false_type {};
};

//trait for const T (forwarding to non-const)
template <class T>
struct mpi_type_traits<const T>: mpi_type_traits<T> {};

//*****************************************************************************
//                                  MPI Iterator Traits
//
// is_contiguous: a range of types sequential in memory
//    has to be a sequential container (i.e. C or C++11 array or vector)
//    and the contained type has to be static
//*****************************************************************************

//in general a range is not contiguous
template <class T, class Iterator>
struct is_contiguous_impl : false_type {};

//specialization for vector<T>
template <class T>
struct is_contiguous_impl<T, typename std::vector<T>::iterator> :
    mpi_type_traits<T>::is_static {};

template <class T>
struct is_contiguous_impl<T, typename std::vector<T>::const_iterator> :
    mpi_type_traits<T>::is_static {};

//specialization for array<T,N>
#ifdef MPP_CXX11_ARRAY
//TODO: Ugly hack. Only works if std::array<T,N>::iterator is the same for any N
//The correct solution would require to get N from the iterator and I don't know how.
//(if it fails it doesn't cause wrong results but is slower)
template <class T>
struct is_contiguous_impl<T, typename std::array<T,2>::iterator> :
    mpi_type_traits<T>::is_static {};

template <class T>
struct is_contiguous_impl<T, typename std::array<T,2>::const_iterator> :
    mpi_type_traits<T>::is_static {};
#endif

//forward type to is_contiguous_impl
template <class Iterator>
struct is_contiguous :
    is_contiguous_impl<typename std::iterator_traits<Iterator>::value_type, Iterator> {};

//Specialization for Pointers (used as Iterators)
template <class T>
struct is_contiguous<T*> : mpi_type_traits<T>::is_static {};

template <class T>
struct is_contiguous<T*const> : mpi_type_traits<T>::is_static {};


// Wrapper around MPI_Get_address with the purpose:
//  1) works with MPI_BOTTOM
//  2) remove problem with const correctness (MPI takes non-const pointers even for send)

namespace {

static inline MPI_Aint MPP_Get_address (const void* ptr) {
    MPI_Aint address;
    if (ptr!=MPI_BOTTOM) {
        MPI_Get_address (const_cast<void*>(ptr), &address);
        return address;
    } else {
        return 0;
    }
}
}

//mpp_unique_ptr is std::unique_ptr or boost::shared_ptr if RVALREF is not available
//move is std::move or no-op if shared_ptr is used
#ifdef MPP_CXX11_RVALREF
template<typename T>
struct mpp_unique_ptr { typedef std::unique_ptr<T> type; };
using std::move;
#else
template<typename T>
struct mpp_unique_ptr { typedef boost::shared_ptr<T> type; };
static inline data_layout &move(data_layout &dt) { return dt; } //No-Op
#endif

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  data_layout
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Does RAII for MPI_Datatype. Also contains the size and address.
class data_layout {
    public:
        data_layout() : storage_(new data_layout_storage(MPI_CHAR, false, 0, MPI_BOTTOM)) {}

        data_layout(MPI_Datatype dt, int size, const void* ptr, bool bFree):
            storage_(new data_layout_storage(dt, bFree, size, const_cast<void*>(ptr))) {}

        MPI_Datatype get() const { return storage_->dt_; }
        void* get_ptr() const { return storage_->ptr_; }

        MPI_Aint get_addr() const {
            return MPP_Get_address(storage_->ptr_);
        }
        int get_size() const { return storage_->size_; }
        //change layout into array of previous layout. Only possible for static types.
        //should not be used directly but only by get_layout(Iter begin, Iter end)
        void make_array(int size) { storage_->size_*=size; }
    private:
        struct data_layout_storage {
            data_layout_storage(MPI_Datatype dt, bool bFree, int size, void* ptr) :
                dt_(dt), bFree_(bFree), size_(size), ptr_(ptr) {}
            ~data_layout_storage() {
                if (bFree_) {
                    MPI_Type_free(&dt_);
                }
            }
            MPI_Datatype dt_;
            bool bFree_;
            int size_; //If container length of it. Otherwise 1.
            void* ptr_;
        };
        mpp_unique_ptr<data_layout_storage>::type storage_;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  get_layout:
//  should always be used instead of directly mpi_type_traits::get_layout
//  currently only does caching
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// List of cached layouts.
extern std::vector<data_layout> cached_layouts;
template<class T>
inline data_layout get_type_layout(const T& raw) {
    static MPI_Datatype cached_type = MPI_DATATYPE_NULL;
    static int cached_size = 0;
    //not possible to cache if non-static and not necessary if named.
    if (!mpi_type_traits<T>::is_static::value) {
        return mpi_type_traits<T>::get_layout(raw);
    } else {
        // thread-safety: possible that two threads enter the if statement
        // but that won't cause a problem because cache is correct and both
        // types will be freed.
        if (cached_type == MPI_DATATYPE_NULL) {
            //assigning first to tmp for thread-safety
            data_layout tmp = mpi_type_traits<T>::get_layout(raw);
            cached_type = tmp.get();
            cached_size = tmp.get_size();
            cached_layouts.push_back(move(tmp));
        }
        return data_layout(cached_type, cached_size, &raw, false);
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// struct_layout_builder: For easy creation of layouts for classes/structs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class struct_layout_builder {
public:
    // Creates a builder for the type passed in to the constructor.
    template <class C>
    explicit struct_layout_builder(const C& c, int nelem) :
            base_(mpi_type_traits<C>::is_static::value?&c:NULL), max_ub(0),
            real_extent(mpi_type_traits<C>::is_static::value?sizeof(C):0) {
        layouts_.reserve(nelem);
    }
    // Only needs to be used to specify which members of the class are sent.
    template <class T>
    inline struct_layout_builder& add (const T& t) {
        layouts_.push_back(get_type_layout(t));
        max_ub = std::max(max_ub, (const void*)(&t+1));
        return *this;
    }
    // Returns the build data_layout
    data_layout build() {
        MPI_Datatype dt;
        assert(!layouts_.empty());
        std::vector<int> size;
        std::vector<MPI_Aint> addr;
        std::vector<MPI_Datatype> e_dt;
        MPI_Aint baseAddr = MPP_Get_address(base_);

        bool isStatic = base_!=NULL;
        size.reserve(layouts_.size()+2);
        addr.reserve(layouts_.size()+2);
        e_dt.reserve(layouts_.size()+2);
        assert((size_t)max_ub-(size_t)baseAddr<=real_extent || !isStatic);
        //Add MPI_LB marker if needed (would be correct to always add)
        //needed if smallest relative address is not zero.
        if (isStatic && layouts_[0].get_addr()!=baseAddr) {
            size.push_back(1);
            addr.push_back(0);
            e_dt.push_back(MPI_LB);
        }
        for (size_t i=0;i<layouts_.size();i++) {
            size.push_back(layouts_[i].get_size());
            assert(layouts_[i].get_addr() >= baseAddr || !isStatic);
            addr.push_back(layouts_[i].get_addr() - baseAddr);
            //one could do block compression as in get_range_layout if types are the same
            //or even if types are not the same if running on non-heterogenous hardware (so that MPI_BYTE copy is OK)
            e_dt.push_back(layouts_[i].get());
        }
        //ADD MPI_UB marker if max_ub-base is not equal to the extent
        if (isStatic && real_extent!=(size_t)max_ub-(size_t)baseAddr) {
            size.push_back(1);
            addr.push_back(real_extent);
            e_dt.push_back(MPI_UB);
        }
        MPI_Type_create_struct (addr.size(), &size.front(), &addr.front(), &e_dt.front(), &dt);
        MPI_Type_commit(&dt);
        return data_layout(dt,1,base_,true);
    }
private:
    //Private because it only works correctly for non-static. But needed from get_range_layout,
    //where we don't have access to the container to pass it as c.
    struct_layout_builder() : base_(0), max_ub(0), real_extent(0) {};
    template <class Iter>
    friend data_layout get_range_layout (Iter begin, Iter end, bool bAbsolute=false);
    const void              *base_ , *max_ub;
    std::vector<data_layout> layouts_;
    size_t                   real_extent;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  get layout for range
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//bAbsolute: guarantee absolute addressed MPI datatype. Otherwise it might or might not be absolute.
template <class Iter>
inline data_layout get_range_layout(Iter begin, Iter end, bool bAbsolute) {
    typedef typename std::iterator_traits<Iter>::value_type T;
    Iter beginPlus1 = begin;
    ++beginPlus1;
    if (begin==end) {
        return data_layout();
    } else if (!bAbsolute && (beginPlus1==end || is_contiguous<Iter>::value )) {
        data_layout dl = get_type_layout(*begin);
        int size = std::distance(begin,end);
        dl.make_array(size);
        return dl;
    } else if (mpi_type_traits<T>::is_static::value) {
        std::vector<int> blockLengths;
        std::vector<MPI_Aint> displacements;
        MPI_Aint last_addr = 0;
        data_layout dl = get_type_layout(*begin);
        int esize = dl.get_size();
        typedef typename std::iterator_traits<Iter>::value_type value_type;
        for (Iter i = begin; i != end; ++i) {
            MPI_Aint curr_addr = MPP_Get_address(&*i); //safe because static
            //For static types (last_addr!=0) it is safe to combine adjacent types into
            //one block. By definition their extent is equal to sizeof.
            if (last_addr!=0 &&
                    curr_addr-last_addr == sizeof(value_type)) {
                blockLengths.back() += esize;
            } else {
                blockLengths.push_back(esize);
                displacements.push_back(curr_addr);
            }
            last_addr = curr_addr;
        }
        MPI_Datatype dt;
        MPI_Type_create_hindexed(blockLengths.size(), &blockLengths.front(),
                &displacements.front(), dl.get(), &dt);
        MPI_Type_commit (&dt);
        return data_layout(dt,1,MPI_BOTTOM,true);
    } else {
        struct_layout_builder builder;
        for (Iter i = begin; i != end; ++i) {
            builder.add (*i);
        }
        return builder.build();
    }
}

// This is a shortcut for defining outside the type_trait definition if a
// class is static (look at definition of get_stat for specifically what static means)
#define SET_MPI_STATIC(type) \
    template <> struct mpi_type_traits<type>::is_static : true_type {};

// primitive type traits
template<>
inline data_layout mpi_type_traits<double>::get_layout(const double& t) {
    return data_layout(MPI_DOUBLE, 1, &t, false);
}
SET_MPI_STATIC(double)

template <>
inline data_layout mpi_type_traits<int>::get_layout(const int& raw) {
    return data_layout(MPI_INT, 1, &raw, false);
}
SET_MPI_STATIC(int)

template <>
inline data_layout mpi_type_traits<char>::get_layout(const char& raw) {
    return data_layout(MPI_CHAR, 1, &raw, false);
}
SET_MPI_STATIC(char)

template <>
inline data_layout mpi_type_traits<float>::get_layout(const float& raw) {
    return data_layout(MPI_FLOAT, 1, &raw, false);
}
SET_MPI_STATIC(float)

template <>
inline data_layout mpi_type_traits<long>::get_layout(const long& raw) {
    return data_layout(MPI_LONG, 1, &raw, false);
}
SET_MPI_STATIC(long)

// ... add missing types here ...

// Enum type wrapper that is always available
template<class E>
static inline data_layout get_enum_layout(const E& e) {
    return data_layout(MPI_BYTE, sizeof(E), &e, false);
}

// Enum type trait for C++11.
// Can only be used for code which doesn't have mpi_type_traits for individual enums
// which are required by C++03 compilers (and thus can only be used by code which
// doesn't need to be backward compatible).
#ifdef MPP_CXX11_ENUM_TRAIT
template<class E>
struct mpi_type_traits<E, typename enable_if<std::is_enum<E>::value>::type>
        : default_mpi_type_trait {
    static inline data_layout get_layout(const E &e) {
        return get_enum_layout(e);
    }
    struct is_static : true_type {};
};
#endif

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  smart pointers
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct ptr_mpi_type_trait : default_mpi_type_trait {
    static inline data_layout get_layout(const T& p) {
        data_layout dt = mpi_type_traits<typename T::element_type>::get_layout(*p);
        return dt;
    }
};

#ifdef MPP_CXX11_RVALREF
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::unique_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::unique_ptr<T> > : ptr_mpi_type_trait<std::unique_ptr<T> > {};
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  boost::shared_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<boost::shared_ptr<T> > : ptr_mpi_type_trait<boost::shared_ptr<T> > {};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  boost::scoped_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<boost::scoped_ptr<T> > : ptr_mpi_type_trait<boost::scoped_ptr<T> > {};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  containers
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::vector<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::vector<T> > : default_mpi_type_trait {
    static inline data_layout get_layout(const std::vector<T>& vec) {
        return get_range_layout(vec.begin(), vec.end());
    }
};
#ifdef MPP_CXX11_ARRAY
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::array<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T, size_t N>
struct mpi_type_traits<std::array<T,N> > : default_mpi_type_trait {
    inline static data_layout get_layout(const std::array<T,N>& vec) {
        return get_range_layout(vec.begin(), vec.end());
    }
    struct is_static : mpi_type_traits<T>::is_static {};
};
#endif //MPP_CXX11_ARRAY

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::list<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::list<T> > : default_mpi_type_trait {
    static data_layout get_layout(const std::list<T>& l) {
        return get_range_layout(l.begin(), l.end());
    }
};
}
