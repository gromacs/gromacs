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

namespace mpi {

//*****************************************************************************
//                                  MPI Type Traits
//*****************************************************************************
template <class T>
struct mpi_type_traits {
    // number of elements of type get_type in the type. Used mainly for vectors of static elements
    // equivalent to blocklength as if used in an mpi call
    static inline size_t get_size(const T& raw) { return 1; }
    static inline const T* get_addr(const T& raw) { return &raw; }
    // true if get_type and get_size is the same for all objects of a class
    // and get_addr is the address of object (no indirection or pointer)
    // and the extend of the mpi type return by get_type is equal to sizeof(T)
    static inline bool is_static() { return false; }
    static inline bool is_named() { return false; } //is a named primitive type
//private: //TODO: currently doesn't compile with private. Is it better to be private?
    static inline MPI_Datatype get_type(const T& raw);

//    template<class C>
//    friend MPI_Datatype get_type(const C& raw);
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  caching wrapper
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// List of cached MPI datatypes.
// Can't use RAII because they need to be freed before MPI_Finalize.
extern std::vector<MPI_Datatype> cached_types;
template<class T>
static inline MPI_Datatype get_type_cached(const T& raw) {
    // TODO: Consider the case where MPI_DATATYPE_NULL is not properly
    // defined when this variable is initialized.
    static MPI_Datatype cache = MPI_DATATYPE_NULL;
    //not possible to cache if non-static and not necessary if named.
    if (!mpi_type_traits<T>::is_static() || mpi_type_traits<T>::is_named()) {
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
//  get_range forward declaration
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace internal
{
template <class T, class InputIterator, class InputIteratorBase>
static inline MPI_Datatype get_type_range_contiguous (InputIterator begin, InputIterator end, InputIteratorBase base);

template <class T, class InputIterator, class InputIteratorBase>
static inline MPI_Datatype get_type_range (InputIterator begin, InputIterator end, InputIteratorBase base);

template <class T, class InputIterator>
static inline MPI_Datatype get_type_range (typename std::vector<T>::const_iterator begin,
                                           typename std::vector<T>::const_iterator end,
                                           InputIterator base) {
    return get_type_range_contiguous<T>(begin, end, base);
}

template <class T, class InputIterator>
static inline MPI_Datatype get_type_range (typename std::vector<T>::iterator begin,
                                           typename std::vector<T>::iterator end,
                                           InputIterator base) {
    return get_type_range_contiguous<T>(begin, end, base);
}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  datatype - RAII for MPI_Datatype, behavior similar to scoped_ptr
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class datatype_base
{
    public:
        inline datatype_base() : dt_(MPI_DATATYPE_NULL) { }
        inline const MPI_Datatype& get() const {
            return dt_;
        }
    protected:
        explicit inline datatype_base(MPI_Datatype dt) : dt_(dt) { }
        inline void swap(datatype_base & b) {
            std::swap(b.dt_,dt_);
        }
        //no extra data is allowed. Code casts datatype* to MPI_Datatype*
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
            if (!mpi_type_traits<T>::is_named() &&
                    !mpi_type_traits<T>::is_static() &&
                    dt_!=MPI_DATATYPE_NULL) {
                MPI_Type_free(&dt_);
            }
        }
        inline void reset(const T& raw) {
            datatype(raw).swap(*this);
        }
};

class datatype_range : public datatype_base {
    public:
        inline datatype_range() :  datatype_base() {}
        template <class Iter1, class Iter2>
        inline datatype_range(Iter1 begin, Iter1 end, Iter2 base) :
            datatype_base(
                    internal::get_type_range<typename std::iterator_traits<Iter1>::value_type> (
                            begin, end, base)) { }

        inline ~datatype_range() {
            if (dt_!=MPI_DATATYPE_NULL) {
                MPI_Type_free(&dt_);
            }
        }

        template <class Iter1, class Iter2>
        inline void reset(Iter1 begin, Iter1 end, Iter2 base) {
            datatype_range(begin, end, base).swap(*this);
        }
};


// This is a shortcut for defining outside the type_trait definition if a
// class is static (look at definition of get_stat for specifically what static means)
#define SET_MPI_STATIC(type) template <> inline bool \
        mpi_type_traits<type>::is_static() \
        {return true;}

// Set type to named and static
#define SET_MPI_PRIMITIVE(type) template <> inline bool \
        mpi_type_traits<type>::is_named() {return true;} \
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


// Wrapper around MPI_Get_address
//  1) remove problem with const correctness
//  2) so that vectors of size 0 don't cause runtime errors.
//TODO: I doubt the usage of address=0 is legal even for size 0.
// Can't this simply be avoided by using MPI_Char for vector of size 0?
// Then part 2 should be removed.
static inline MPI_Aint MPP_Get_address (const void* location) {
    MPI_Aint address;
    if (location != NULL) {
        MPI_Get_address (const_cast<void *>(location), &address);
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
class mpi_type_builder {
    MPI_Aint                  base;
    std::vector<MPI_Datatype> type;
    std::vector<int>          size;
    std::vector<MPI_Aint>     addr;
    std::vector<bool>         bFree;
public:
    // Creates a builder for the type passed in to the constructor.
    template <class C>
    explicit mpi_type_builder(const C& c, int nelem = 0) {
        base = MPP_Get_address(mpi_type_traits<C>::get_addr(c));
        type.reserve(nelem);
        size.reserve(nelem);
        addr.reserve(nelem);
        bFree.reserve(nelem);
    }
    // Only needs to be used to specify which members of the class are sent.
    template <class T>
    inline void add (const T& t) {
        type.push_back(get_type_cached(t));
        size.push_back(mpi_type_traits<T>::get_size(t));
        MPI_Aint address = MPP_Get_address(mpi_type_traits<T>::get_addr(t));
        addr.push_back(address-base);
        bFree.push_back(!mpi_type_traits<T>::is_named() && !mpi_type_traits<T>::is_static());
    }
    // Returns an MPI_Datatype for the type that the builder is tied to.
    MPI_Datatype build();
};

#ifdef MPP_CXX11_RVALREF
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::unique_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::unique_ptr<T> > {
    static inline MPI_Datatype get_type(const std::unique_ptr<T>& p) {
        return mpi_type_traits<T>::get_type(*p);
    }
    static inline size_t get_size(const std::unique_ptr<T>& p) {
        return mpi_type_traits<T>::get_size(*p);
    }
    static inline const T* get_addr(const std::unique_ptr<T>& p) {
        return mpi_type_traits<T>::get_addr(*p);
    }
    static inline bool is_static() { return false; }
    static inline bool is_named() {
        return mpi_type_traits<T>::is_named();
    }
};
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  boost::shared_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<boost::shared_ptr<T> > {
    static inline MPI_Datatype get_type(const boost::shared_ptr<T>& p) {
        return mpi_type_traits<T>::get_type(*p);
    }
    static inline size_t get_size(const boost::shared_ptr<T>& p) {
        return mpi_type_traits<T>::get_size(*p);
    }
    static inline const T* get_addr(const boost::shared_ptr<T>& p) {
        return mpi_type_traits<T>::get_addr(*p);
    }
    static inline bool is_static() { return false; }
    static inline bool is_named() {
        return mpi_type_traits<T>::is_named();
    }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  boost::scoped_ptr<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<boost::scoped_ptr<T> > {
    static inline MPI_Datatype get_type(const boost::scoped_ptr<T>& p) {
        return mpi_type_traits<T>::get_type(*p);
    }
    static inline size_t get_size(const boost::scoped_ptr<T>& p) {
        return mpi_type_traits<T>::get_size(*p);
    }
    static inline const T* get_addr(const boost::scoped_ptr<T>& p) {
        return mpi_type_traits<T>::get_addr(*p);
    }
    static inline bool is_static() { return false; }
    static inline bool is_named() {
        return mpi_type_traits<T>::is_named();
    }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  get_type_range
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//TODO: the whole section should probably be moved somewhere else

//only used inside mpi_type_traits (template specialization is not allowed inside struct)
namespace internal
{
template <class InputIterator, class InputIteratorBase>
static inline MPI_Datatype get_type_range_non_static (InputIterator begin, InputIterator end, InputIteratorBase base) {
    if (begin!=end) {
        mpi_type_builder builder(*base);
        for (InputIterator i = begin; i != end; ++i) {
            builder.add (*i);
        }
        return builder.build();
    } else {
        MPI_Datatype dt; //has to match is_named
        MPI_Type_dup(MPI_CHAR, &dt);
        return dt;
    }
}

//Specialization for continuous iterators like std::vector::iterator.
template <class T, class InputIterator, class InputIteratorBase>
static inline MPI_Datatype get_type_range_contiguous (InputIterator begin, InputIterator end, InputIteratorBase base) {
    //for vector::iterator or array::iterator
    if (begin!=end && mpi_type_traits<T>::is_static()) {
        MPI_Datatype dt;
        datatype<T> elem_dt(*begin);
        int length = end-begin;
        MPI_Aint displ = MPP_Get_address(&*begin) - MPP_Get_address(&*base);
        MPI_Type_create_hindexed(1, &length, &displ, elem_dt.get(), &dt);
        MPI_Type_commit (&dt);
        return dt;
    } else {
        return get_type_range_non_static(begin, end, base);
    }
}

//Create a type for any range of two iterators. The type needs to be used with &*begin as address.
template <class T, class InputIterator, class InputIteratorBase>
static inline MPI_Datatype get_type_range (InputIterator begin, InputIterator end, InputIteratorBase base) {
    if (begin!=end && mpi_type_traits<T>::is_static()) {
        std::vector<int> blockLengths;
        std::vector<MPI_Aint> displacements;
        //mpi_type_traits::get_addr not needed because guaranteed to be same for static types
        MPI_Aint base_addr = MPP_Get_address(&*base);
        MPI_Aint last_addr = base_addr;
        typedef typename std::iterator_traits<InputIterator>::value_type value_type;
        size_t elem_size = mpi_type_traits<value_type>::get_size(*begin);
        for (InputIterator i = begin; i != end; ++i) {
            MPI_Aint curr_addr = MPP_Get_address(&*i);
            if (curr_addr-last_addr == sizeof(value_type)) {
                blockLengths.back() += elem_size;
            } else {
                blockLengths.push_back(elem_size);
                displacements.push_back(curr_addr-base_addr);
            }
            last_addr = curr_addr;
        }
        MPI_Datatype dt;
        datatype<T> elem_dt(*begin);
        MPI_Type_create_hindexed(blockLengths.size(), &blockLengths.front(), &displacements.front(), elem_dt.get(), &dt);
        MPI_Type_commit (&dt);
        return dt;
    } else {
        return get_type_range_non_static(begin, end, base);
    }
}

//simple version of boost::permutation_iterator
template< class ElementIterator, class IndexIterator>
class permutation_iterator
{
private:
        typedef typename std::iterator_traits<ElementIterator>::reference elt_ref;
public:
        permutation_iterator() : m_elt_iter(), m_idx_iter()  {}

        permutation_iterator(ElementIterator x, IndexIterator y)
                : m_elt_iter(x), m_idx_iter(y) {}

        //returns *(m_elt_iter+*m_idx_iter)
        elt_ref operator*() {
            ElementIterator curr = m_elt_iter;
            std::advance(curr, *m_idx_iter);
            return *curr;
        }
        permutation_iterator& operator++() { ++m_idx_iter; return *this; }
        bool operator!=(const permutation_iterator &other) {
            assert(m_elt_iter==other.m_elt_iter);
            return m_idx_iter!=other.m_idx_iter;
        }
private:
        ElementIterator m_elt_iter;
        IndexIterator m_idx_iter;
};
template <class ElementIterator, class IndexIterator>
permutation_iterator<ElementIterator, IndexIterator>
make_permutation_iterator( ElementIterator e, IndexIterator i )
{
    return permutation_iterator<ElementIterator, IndexIterator>( e, i );
}
} //end namespace internal
} //end namespace mpi
namespace std {
template< class ElementIterator, class IndexIterator>
struct iterator_traits<mpi::internal::permutation_iterator<ElementIterator, IndexIterator> > {
        typedef typename iterator_traits<ElementIterator>::value_type value_type;
};
}
namespace mpi {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::vector<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::vector<T> > {

    static inline size_t get_size(const std::vector<T>& vec) {
        if (vec.size()>0) {
            if (mpi_type_traits<T>::is_static()) {
                return vec.size() * mpi_type_traits<T>::get_size(vec.front());
            } else { return 1; }
        } else { return 0; }
    }

    static inline MPI_Datatype get_type(const std::vector<T>& vec) {
        if (mpi_type_traits<T>::is_static()) {
            if (vec.size()>0) {
            return mpi_type_traits<T>::get_type(vec.front());
            } else {
                if (mpi_type_traits<T>::is_named()) {
                    return MPI_CHAR;
                } else {
                    MPI_Datatype dt; //has to match is_named
                    MPI_Type_dup(MPI_CHAR, &dt);
                    return dt;
                }
            }
        } else {
            return get_type_range (vec.begin(), vec.end());
        }
    }

    static inline const T* get_addr(const std::vector<T>& vec) {
        if (vec.size()==0) { return 0; }
        return  mpi_type_traits<T>::get_addr(vec.front());
    }
    template <class InputIterator, class InputIteratorBase>
    static inline MPI_Datatype get_type_range (InputIterator begin, InputIterator end, InputIteratorBase base) {
            return internal::get_type_range<T>(begin, end, base);
    }
    template <class InputIterator>
    static inline MPI_Datatype get_type_range (InputIterator begin, InputIterator end) {
            return internal::get_type_range<T>(begin, end, begin);
    }

    static inline bool is_static() { return false; }
    static inline bool is_named() {
        return mpi_type_traits<T>::is_named() && mpi_type_traits<T>::is_static();
    }
};
#ifdef MPP_CXX11_ARRAY //currently always wrong - so far not needed
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::array<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T, size_t N>
struct mpi_type_traits<std::array<T,N> > {

    inline static size_t get_size(const std::array<T,N>& vec) { return N; }

    inline static MPI_Datatype get_type(const std::array<T,N>& vec) {
        return  mpi_type_traits<T>::get_type( T() );
    }

    static inline const T* get_addr(const std::array<T,N>& vec) {
        return  mpi_type_traits<T>::get_addr( vec.front() );
    }

    static inline bool is_static() {
        //TODO: verify we support is_static with size>1?
        return mpi_type_traits<T>::is_static();
    }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::list<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::list<T> > {

    static inline size_t get_size(const std::list<T>& l) { return 1; }

    static MPI_Datatype get_type(const std::list<T>& l) {
        // we have to get the create an MPI_Datatype containing the offsets
        // of the current object

        // we consider the offsets starting from the first element
        std::vector<MPI_Aint> address( l.size() );
        std::vector<int> dimension( l.size() );
        std::vector<MPI_Datatype> types( l.size() );

        std::vector<int>::iterator dim_it = dimension.begin();
        std::vector<MPI_Aint>::iterator address_it = address.begin();
        std::vector<MPI_Datatype>::iterator type_it = types.begin();

        MPI_Aint base_address = MPP_Get_address(&l.front());

        *(type_it++) = get_type_cached( l.front() ); //TODO use datatype for RAII/free
        *(dim_it++) = mpi_type_traits<T>::get_size( l.front() );
        *(address_it++) = 0;

        typename std::list<T>::const_iterator begin = l.begin();
        ++begin;
        std::for_each(begin, l.cend(), [&](const T& curr) {
            assert( address_it != address.end() &&
                    type_it != types.end() &&
                    dim_it != dimension.end() );

            *address_it = MPP_Get_address(&curr);
            *(address_it++) -= base_address;
            *(type_it++) =  get_type_cached( curr );
            *(dim_it++) = mpi_type_traits<T>::get_size( curr );

        }
        );

        MPI_Datatype list_dt;
        MPI_Type_create_struct(l.size(), &dimension.front(), &address.front(), &types.front(), &list_dt);
        MPI_Type_commit( &list_dt );

        return list_dt;
    }

    static inline const T* get_addr(const std::list<T>& list) {
        return  mpi_type_traits<T>::get_addr( list.front() );
    }

    static inline bool is_static() { return false; }
};
#endif
}
