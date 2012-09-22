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

#ifdef MPP_CXX11_TYPE_TRAITS
#include <type_traits>
#endif

namespace mpi {

#ifdef MPP_CXX11_TYPE_TRAITS
//using stl true/false types
typedef std::true_type TrueType;
typedef std::false_type FalseType;
#else
//define equivalent true/false types
struct FalseType { enum { value = false }; };
struct TrueType { enum { value = true }; };
#endif

//*****************************************************************************
//                                  MPI Type Traits
//*****************************************************************************
template <class T>
struct mpi_type_traits {
    static inline MPI_Datatype get_type(const T& raw);
    // is_static should be true if get_type is the same for all objects of a class
    // and the extent of the mpi type return by get_type is equal to sizeof(T)
    // also it cannot contain any pointers.
    struct is_static : FalseType {};
    // is_named should be true if get_type returns a MPI named(/primitive) type.
    struct is_named : FalseType {};
};

//base class used as shortcut for partial specialization
struct default_mpi_type_trait {
    struct is_static : FalseType {};
    struct is_named : FalseType {};
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
struct is_contiguous_impl : FalseType {};

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


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  caching wrapper
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// List of cached MPI datatypes.
// Can't use RAII because they need to be freed before MPI_Finalize.
extern std::vector<MPI_Datatype> cached_types;
template<class T>
static inline MPI_Datatype get_type_cached(const T& raw) {
    static MPI_Datatype cache = MPI_DATATYPE_NULL;
    //not possible to cache if non-static and not necessary if named.
    if (!mpi_type_traits<T>::is_static::value || mpi_type_traits<T>::is_named::value) {
        return mpi_type_traits<T>::get_type(raw);
    } else {
        // thread-safety: possible that two threads enter the if statement
        // but that won't cause a problem because cache is correct and both
        // types will be freed.
        if (cache == MPI_DATATYPE_NULL) {
            MPI_Datatype tmp = mpi_type_traits<T>::get_type(raw);
            cached_types.push_back(tmp);
            cache = tmp;
        }
        return cache;
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  datatype - RAII for MPI_Datatype, behavior similar to scoped_ptr
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class datatype_base
{
    public:
        //TODO: replace MPI_CHAR with MPP_TYPE_EMPTY (MPI_Type_contigous(0,MPI_CHAR))
        inline datatype_base() : dt_(MPI_CHAR) { }
        inline const MPI_Datatype& get() const {
            return dt_;
        }
    protected:
        explicit inline datatype_base(MPI_Datatype dt) : dt_(dt) { }
        //swap is used for reset (modeled after scoped_ptr)
        inline void swap(datatype_base & b) {
            std::swap(b.dt_,dt_);
        }
        //no extra data is allowed. Code arrays of datatype to arrays of MPI_Datatype
        MPI_Datatype dt_;
    private:
        datatype_base(datatype_base& ); //disallow copy and assign
        const datatype_base& operator=( const datatype_base& );
};


template <class T>
class datatype : public datatype_base {
    public:
        inline datatype() :  datatype_base() {}
        explicit inline datatype(const T& raw) : datatype_base(get_type_cached<T>(raw)) { }
        inline ~datatype() {
            if (!mpi_type_traits<T>::is_named::value &&
                    !mpi_type_traits<T>::is_static::value &&
                    dt_!=MPI_CHAR) {
                MPI_Type_free(&dt_);
            }
        }
        inline void reset(const T& raw) {
            datatype(raw).swap(*this);
        }
};

// Wrapper around MPI_Get_address with the purpose:
//  1) return pointer based on weather static/non-static as MPI_BOTTOM or real address
//  2) remove problem with const correctness (MPI takes non-const pointers even for send)

template <class C>
inline void* MPP_Get_ptr (const C& location) {
    if (mpi_type_traits<C>::is_static::value) {
        return (void*)(&location);
    } else {
        return MPI_BOTTOM;
    }
}

template <class C>
static inline MPI_Aint MPP_Get_address (const C& location) {
    void* ptr = MPP_Get_ptr(location);
    MPI_Aint address;
    if (ptr!=MPI_BOTTOM) {
        MPI_Get_address (ptr, &address);
        return address;
    } else {
        return 0;
    }
}

// For static classes it would be a good idea to cache those data types
//     in perhaps a static class or map used only by the builder (Similar to
//     how comm::world is).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// mpi_type_builder: For easy creation of custom MPI_Datatypes
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace internal {
template <class InputIterator>
static inline MPI_Datatype get_type_range (InputIterator, InputIterator);
}

class mpi_type_builder {
public:
    // Creates a builder for the type passed in to the constructor.
    template <class C>
    explicit mpi_type_builder(const C& c, int nelem) : base(MPP_Get_address(c)),
            real_extent(sizeof(C)), curr_extent(0) {
        type.reserve(nelem);
        addr.reserve(nelem);
        bFree.reserve(nelem);
    }
    // Only needs to be used to specify which members of the class are sent.
    template <class T>
    inline void add (const T& t) {
        type.push_back(get_type_cached(t));
        MPI_Aint address = MPP_Get_address(t)-base;
        addr.push_back(address);
        bFree.push_back(!mpi_type_traits<T>::is_named::value && !mpi_type_traits<T>::is_static::value);
        curr_extent = address+sizeof(T);
    }
    // Returns an MPI_Datatype for the type that the builder is tied to.
    MPI_Datatype build();
private:
    //Private because it only works correctly for non-static. But needed from get_type_range,
    //where we don't have access to the container to pass it as c.
    mpi_type_builder() : base(0), real_extent(0), curr_extent(0) {};
    template <class InputIterator>
    friend inline MPI_Datatype internal::get_type_range (InputIterator, InputIterator);
    MPI_Aint base;
    std::vector<MPI_Datatype> type;
    std::vector<MPI_Aint>     addr;
    std::vector<bool>         bFree;
    size_t                    real_extent, curr_extent;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  get_type_range
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace internal {

template <class InputIterator>
static inline MPI_Datatype get_type_range (InputIterator begin, InputIterator end) {
    typedef typename std::iterator_traits<InputIterator>::value_type T;
    if (begin==end) {
        MPI_Datatype dt;//TODO better really empty. Should enable in gather to set count always to 1.
        MPI_Type_dup(MPI_CHAR, &dt);
        return dt;
    } else if (is_contiguous<InputIterator>::value) {
        MPI_Datatype dt;
        datatype<T> elem_dt(*begin);
        int length = std::distance(begin,end);
        MPI_Aint addr = MPP_Get_address(*begin);
        MPI_Type_create_hindexed(1, &length, &addr, elem_dt.get(), &dt);
        MPI_Type_commit (&dt);
        return dt;
    } else if (mpi_type_traits<T>::is_static::value) {
        std::vector<int> blockLengths;
        std::vector<MPI_Aint> displacements;
        MPI_Aint last_addr = 0;
        typedef typename std::iterator_traits<InputIterator>::value_type value_type;
        for (InputIterator i = begin; i != end; ++i) {
            MPI_Aint curr_addr = MPP_Get_address(*i);
            //For static types (last_addr!=0) it is safe to combine adjacent types into
            //one block. By definition their extend is equal to sizeof.
            if (last_addr!=0 &&
                    curr_addr-last_addr == sizeof(value_type)) {
                blockLengths.back() += 1;
            } else {
                blockLengths.push_back(1);
                displacements.push_back(curr_addr);
            }
            last_addr = curr_addr;
        }
        MPI_Datatype dt;
        datatype<T> elem_dt(*begin);
        MPI_Type_create_hindexed(blockLengths.size(), &blockLengths.front(),
                &displacements.front(), elem_dt.get(), &dt);
        MPI_Type_commit (&dt);
        return dt;
    } else {
        mpi_type_builder builder;
        for (InputIterator i = begin; i != end; ++i) {
            builder.add (*i);
        }
        return builder.build();
    }
}

} //end internal

//TODO: try to move to other datatypes classes.
class datatype_range : public datatype_base {
    public:
        inline datatype_range() :  datatype_base() {}
        template <class Iter1>
        inline datatype_range(Iter1 begin, Iter1 end) :
            datatype_base(
                internal::get_type_range<Iter1> (begin, end)) { }
        inline ~datatype_range() {
            //doesn't require is_static check because ranges are always dynamic
            //and always create a new datatype which needs to be freed.
            if (dt_!=MPI_CHAR) {
                MPI_Type_free(&dt_);
            }
        }

        template <class Iter1>
        inline void reset(Iter1 begin, Iter1 end) {
            datatype_range(begin, end).swap(*this);
        }
};


// This is a shortcut for defining outside the type_trait definition if a
// class is static (look at definition of get_stat for specifically what static means)
#define SET_MPI_STATIC(type) \
    template <> struct mpi_type_traits<type>::is_static : TrueType {};

// Set type to named and static
#define SET_MPI_PRIMITIVE(type) \
    template <> struct mpi_type_traits<type>::is_named : TrueType {}; \
    SET_MPI_STATIC(type)

// primitive type traits
template<>
inline MPI_Datatype mpi_type_traits<double>::get_type(const double&) {
    return  MPI_DOUBLE;
}
SET_MPI_PRIMITIVE(double)

template <>
inline MPI_Datatype mpi_type_traits<int>::get_type(const int&) {
    return MPI_INT;
}
SET_MPI_PRIMITIVE(int)

template <>
inline MPI_Datatype mpi_type_traits<char>::get_type(const char&) {
    return MPI_CHAR;
}
SET_MPI_PRIMITIVE(char)

template <>
inline MPI_Datatype mpi_type_traits<float>::get_type(const float&) {
    return MPI_FLOAT;
}
SET_MPI_PRIMITIVE(float)

template <>
inline MPI_Datatype mpi_type_traits<long>::get_type(const long&) {
    return MPI_LONG;
}
SET_MPI_PRIMITIVE(long)

// ... add missing types here ...

// Enum type wrapper that is always available
template<class E>
static inline MPI_Datatype get_type_enum(const E& e) {
    MPI_Datatype dt;
    MPI_Type_contiguous(sizeof(E), MPI_BYTE, &dt);
    MPI_Type_commit(&dt);
    return dt;
}

// Enum type trait for C++11 only
#ifdef MPP_CXX11_TYPE_TRAITS
/* FIXME: doesn't compile.
template<class E>
inline typename enable_if<std::is_enum<E>::value, MPI_Datatype>::type mpi_type_traits<E>::get_type(const E &) {
    return get_type_enum<E>();
} */
#endif

template <class T>
struct ptr_mpi_type_trait : default_mpi_type_trait {
    static inline MPI_Datatype get_type(const T& p) {
        if (mpi_type_traits<typename T::element_type>::is_static::value) {
            datatype<typename T::element_type> elem_dt(*p);
            MPI_Datatype dt;
            MPI_Aint addr = MPP_Get_address(*p);
            int length = 1;
            MPI_Type_create_hindexed(1, &length, &addr, elem_dt.get(), &dt);
    //        MPI_Aint lb, extent; //possible alternative.
    //        MPI_Type_get_extent(elem_dt.get(), &lb, &extent);
    //        MPI_Type_create_resized(elem_dt.get(), MPP_Get_address(*p), extent, &dt);
            MPI_Type_commit(&dt);
            return dt;
        } else {
            return mpi_type_traits<typename T::element_type>::get_type(*p);
        }
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
//  std::vector<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::vector<T> > : default_mpi_type_trait {
    static inline MPI_Datatype get_type(const std::vector<T>& vec) {
        return internal::get_type_range (vec.begin(), vec.end());
    }
};
#ifdef MPP_CXX11_ARRAY //currently always wrong - so far not needed
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::array<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T, size_t N>
struct mpi_type_traits<std::array<T,N> > : default_mpi_type_trait {
    inline static MPI_Datatype get_type(const std::array<T,N>& vec) {
        //not using get_type_range because std::array is static.
        MPI_Datatype dt;
        datatype<T> elem_dt(vec[0]);
        MPI_Type_contiguous(vec.size(), elem_dt.get(), &dt);
        MPI_Type_commit(&dt);
        return dt;
    }
    struct is_static : TrueType {};
};
#endif //MPP_CXX11_ARRAY

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::list<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::list<T> > : default_mpi_type_trait {
    static MPI_Datatype get_type(const std::list<T>& l) {
        return internal::get_type_range (l.begin(), l.end());
    }
};
}
