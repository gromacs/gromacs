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
    static inline MPI_Datatype get_type(const T& raw);
    // number of elements of type get_type in the type. Used mainly for vectors of static elements
    // equivalent to blocklength as if used in an mpi call
    static inline const size_t get_size(const T& raw) { return 1; }
    static inline const T* get_addr(const T& raw) { return &raw; }
    //true if get_type and get_size is the same for all objects of a class
    //and get_addr is the address of object (no indirection or pointer)
    static inline bool get_static() { return false; }
};

// This is a shortcut for defining outside the type_trait definition if a
// class is static (look at definition of get_stat for specifically what static means)
#define SET_MPI_STATIC(type) template <> inline bool \
        mpi_type_traits<type>::get_static() \
        {return true;}

// primitive type traits
template<>
inline MPI_Datatype mpi_type_traits<double>::get_type(const double&) {
    return  MPI_DOUBLE;
}
SET_MPI_STATIC(double)

template <>
inline MPI_Datatype mpi_type_traits<int>::get_type(const int&) {
    return MPI_INT;
}
SET_MPI_STATIC(int)

template <>
inline MPI_Datatype mpi_type_traits<char>::get_type(const char&) {
    return MPI_CHAR;
}
SET_MPI_STATIC(char)

template <>
inline MPI_Datatype mpi_type_traits<float>::get_type(const float&) {
    return MPI_FLOAT;
}
SET_MPI_STATIC(float)

template <>
inline MPI_Datatype mpi_type_traits<long>::get_type(const long&) {
    return MPI_LONG;
}
SET_MPI_STATIC(long)

// ... add missing types here ...

//Wrapper around MPI_Get_address so that vectors of size 0 don't cause
//    runtime errors.
static inline void MPP_Get_address (void* location, MPI_Aint* address) {
    if (location != NULL) {
        MPI_Get_address (location, address);
    } else {
        *address = 0;
    }
}

//For static classes it would be a good idea to cache those data types
//    in perhaps a static class or map used only by the builder (Similar to
//    how comm::world is).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// mpi_type_builder: For easy creation of custom classes
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class mpi_type_builder {
    MPI_Aint                  base;
    std::vector<MPI_Datatype> type;
    std::vector<int>          size;
    std::vector<MPI_Aint>     addr;
public:
    template <class C>
    explicit mpi_type_builder(const C& c) { MPP_Get_address(const_cast<void *>(
        static_cast<const void *>(mpi_type_traits<C>::get_addr(c))), &base); }
    template <class T>
    inline void add (const T& t) {
        type.push_back(mpi_type_traits<T>::get_type(t));
        size.push_back(mpi_type_traits<T>::get_size(t));
        MPI_Aint address;
        MPP_Get_address(const_cast<void *>(static_cast<const void *>(
                mpi_type_traits<T>::get_addr(t))), &address);
        addr.push_back(address-base);
    }
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
    static inline bool get_static() { return false; }
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
    static inline bool get_static() { return false; }
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
    static inline bool get_static() { return false; }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::vector<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::vector<T> > {

    static inline const size_t get_size(const std::vector<T>& vec) {
        if (vec.size()>0) {
            if (mpi_type_traits<T>::get_static()) {
                return vec.size() * mpi_type_traits<T>::get_size(vec.front());
            } else { return 1; }
        } else { return 0; }
    }

    static inline MPI_Datatype get_type(const std::vector<T>& vec) {
        if (vec.size()>0 && mpi_type_traits<T>::get_static()) {
            return mpi_type_traits<T>::get_type(vec.front());
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
        for (size_t i=begin; i<end; i++) {
            sel.push_back(i);
        }
        return get_type_selection (vec, sel, begin, end);
    }

    static inline MPI_Datatype get_type_selection (const std::vector<T>& vec, const std::vector<int>& sel, size_t begin, size_t end) {
        if (vec.size()>0 && mpi_type_traits<T>::get_static()) {
            std::vector<int> displacements (1);
            std::vector<int> blockLengths  (1);
            MPI_Datatype dt;
            int count=0;
            //Determine block lengths and displacements
            for (unsigned int i=0; i<vec.size()-1; i++) {
                blockLengths[count]++;
                if (sel[i+1] != sel[i]+1) {
                    count++;
                    blockLengths.push_back(0);
                    displacements.push_back(i);
                }
            }
            MPI_Type_indexed (count, &blockLengths.front(), &displacements.front(), mpi_type_traits<T>::get_type(vec.front()), &dt);
            MPI_Type_commit (&dt);
            return dt;
        } else {
            mpi_type_builder builder(vec);
            for (size_t i=begin; i<end; i++) {
                builder.add (vec[sel[i]]);
            }
            return builder.build();
        }
    }

    static inline bool get_static() { return false; }
};
#ifdef MPP_CXX11_ARRAY //currently always wrong - so far not needed
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::array<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T, size_t N>
struct mpi_type_traits<std::array<T,N> > {

    inline static const size_t get_size(const std::array<T,N>& vec) { return N; }

    inline static MPI_Datatype get_type(const std::array<T,N>& vec) {
        return  mpi_type_traits<T>::get_type( T() );
    }

    static inline const T* get_addr(const std::array<T,N>& vec) {
        return  mpi_type_traits<T>::get_addr( vec.front() );
    }

    static inline bool get_static(const std::array<T,N>& vec) {
        //TODO: is this correct?
        return mpi_type_traits<T>::get_static();
    }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  std::list<T> traits
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T>
struct mpi_type_traits<std::list<T> > {

    static inline const size_t get_size(const std::list<T>& l) { return 1; }

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
        MPI_Get_address(const_cast<T*>(&l.front()), &base_address);

        *(type_it++) = mpi_type_traits<T>::get_type( l.front() );
        *(dim_it++) = mpi_type_traits<T>::get_size( l.front() );
        *(address_it++) = 0;

        typename std::list<T>::const_iterator begin = l.begin();
        ++begin;
        std::for_each(begin, l.cend(), [&](const T& curr) {
            assert( address_it != address.end() &&
                    type_it != types.end() &&
                    dim_it != dimension.end() );

            MPI_Get_address(const_cast<T*>(&curr), &*address_it);
            *(address_it++) -= base_address;
            *(type_it++) =  mpi_type_traits<T>::get_type( curr );
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

    static inline bool get_static() { return false; }
};
#endif
}
