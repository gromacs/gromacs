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
#ifndef MPP_NO_MPI_INCL
#include <mpi.h>
#endif

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
    MPI_Datatype cache = MPI_DATATYPE_NULL;
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
//  datatype - RAII for MPI_Datatype, behavior similar to scoped_ptr
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
class datatype {
    public:
        explicit inline datatype(const T& raw) : dt_(get_type_cached<T>(raw)) { }
        //TODO: remove - requires alternative solution for get_type_{selection|range}
        explicit inline datatype(MPI_Datatype dt) : dt_(dt) { }
        inline datatype() : dt_(MPI_DATATYPE_NULL) { }
        inline ~datatype() {
            if (!mpi_type_traits<T>::is_named() &&
                    !mpi_type_traits<T>::is_static() &&
                    dt_!=MPI_DATATYPE_NULL) {
                MPI_Type_free(&dt_);
            }
        }
        inline const MPI_Datatype& get() const {
            return dt_;
        }
        inline void reset(const T& raw) {
            datatype(raw).swap(*this);
        }
        inline void reset(MPI_Datatype dt) { //TODO: remove
            datatype(dt).swap(*this);
        }
        inline void swap(datatype & b) {
            std::swap(b.dt_,dt_);
        }
    private:
        //no extra data is allowed. Code casts datatype* to MPI_Datatype*
        MPI_Datatype dt_;
        datatype(datatype& ); //disallow copy and assign
        const datatype& operator=( const datatype& );
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

// Wrapper around MPI_Get_address so that vectors of size 0 don't cause
//     runtime errors.
static inline void MPP_Get_address (void* location, MPI_Aint* address) {
    if (location != NULL) {
        MPI_Get_address (location, address);
    } else {
        *address = 0;
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
        MPP_Get_address(const_cast<void *>(
            static_cast<const void *>(mpi_type_traits<C>::get_addr(c))), &base);
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
        MPI_Aint address;
        MPP_Get_address(const_cast<void *>(static_cast<const void *>(
                mpi_type_traits<T>::get_addr(t))), &address);
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
            return get_type_range (vec, 0, vec.size());
        }
    }

    static inline const T* get_addr(const std::vector<T>& vec) {
        if (vec.size()==0) { return 0; }
        return  mpi_type_traits<T>::get_addr(vec.front());
    }

    static inline MPI_Datatype get_type_range (const std::vector<T>& vec, size_t begin, size_t end) {
        assert(end <= vec.size());
        std::vector<int> sel;
        for (size_t i=0; i<vec.size(); i++) {
            sel.push_back(i);
        }
        return get_type_selection (vec, sel, begin, end);
    }

    //This is for a vector of blocks from each core. Begin should point to the
    //    first element of the block and end is 1 past the end of the block
    static inline MPI_Datatype get_type_selection (const std::vector<T>& vec, const std::vector<int>& sel, const int begin, const int end) {
        if (vec.size()>0 && mpi_type_traits<T>::is_static()) {
            std::vector<int> displacements (1, sel[begin]);
            std::vector<int> blockLengths  (1);
            int count=1;
            //Determine block lengths and displacements
            for (int i=begin; i<end; i++) {
                blockLengths[count-1]++;
                if (i+1 != end && sel[i+1] != sel[i]+1) {
                    count++;
                    blockLengths.push_back(0);
                    displacements.push_back(sel[i+1]);
                }
            }
            MPI_Datatype dt;
            datatype<T> elem_dt(vec.front()); //TODO: check whether works without explicit template
            MPI_Type_indexed (count, &blockLengths.front(), &displacements.front(), elem_dt.get(), &dt);
            MPI_Type_commit (&dt);
            return dt;
        } else {
            mpi_type_builder builder(vec, end-begin);
            for (int i=begin; i<end; i++) {
                builder.add (vec[sel[i]]);
            }
            return builder.build();
        }
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
        //TODO: is this correct?
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

        MPI_Aint base_address;
        MPP_Get_address(const_cast<T*>(&l.front()), &base_address);

        *(type_it++) = get_type_cached( l.front() ); //TODO use datatype for RAII/free
        *(dim_it++) = mpi_type_traits<T>::get_size( l.front() );
        *(address_it++) = 0;

        typename std::list<T>::const_iterator begin = l.begin();
        ++begin;
        std::for_each(begin, l.cend(), [&](const T& curr) {
            assert( address_it != address.end() &&
                    type_it != types.end() &&
                    dim_it != dimension.end() );

            MPP_Get_address(const_cast<T*>(&curr), &*address_it);
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
