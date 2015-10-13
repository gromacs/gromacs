/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
    \ingroup module_data_structures
    \inpublicapi

    \brief AlignedAllocator enables allocator-aware data structures to
           allocate memory with a specific alignment as for example
           needed in SIMD computations.

    \author R. Thomas Ullmann <tullman@gwdg.de>

    \copyright GROMACS license

    \date Feb 2015
 */
#ifndef GMX_IrregArray4D_H
#define GMX_IrregArray4D_H


#include <new>
#include <stdexcept>
#include <cstdlib>
#if __cplusplus >= 201103L
// #include <type_traits>
#endif

/*! \class AlignedAllocator

    \brief allocate aligned memory, e.g., for SIMD computation

    \ingroup module_data_structures
    \inpublicapi

    \tparam   T           data type for which memory will be allocated
    \tparam   Alignment   memory alignment of the data element(s) in byte
 */
#if __cplusplus >= 201103L
template <class T, unsigned alignment = alignof(T)>
class AlignedAllocator;
#else
template <class T, unsigned alignment>
class AlignedAllocator;
#endif

//specialize for void:
template <unsigned alignment> class AlignedAllocator<void, alignment>
{
    public:
        typedef void* pointer;
        typedef const void* const_pointer;
        // reference-to-void members are impossible.
        typedef void value_type;
        template <class U> struct rebind { typedef AlignedAllocator<U, alignment> other; };
};

// ---------------------------------------------------------------------------------------
// __cplusplus >= 201103L does not guarantee availability of std::align
// thus this whitelist
// tested versions, lower may be possible in some cases, add further compilers here?
// GCC >= 5.2
// Clang >= 3.7, Apple Clang does not use the same versioning scheme, excluded for now
// MS VC++ tested 2015 (__MSC_VER == 1900), supported since 2012, 1800 according to MS
// Intel C++ compiler >15
// ---------------------------------------------------------------------------------------
//! \todo replace this macro compiler whitelist by a std::align availability check in the CMake machinery?
#if __cplusplus >= 201103L
 #ifdef __GNUC__
  #if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 50200
   #define GMX_ALIGNED_ALLOCATOR_HAVE_STD_ALIGN
  #endif
 #endif
 #ifdef __MSC_VER
  #if __MSC_VER >= 1900
   #define GMX_ALIGNED_ALLOCATOR_HAVE_STD_ALIGN
  #endif
 #endif
 #ifdef __clang__
  #ifndef __apple_build_version__
   #if (__clang_major__ * 10000 + __clang_minor__ * 100) >= 30700
    #define GMX_ALIGNED_ALLOCATOR_HAVE_STD_ALIGN
   #endif
  #endif
 #endif
 #ifdef __INTEL_COMPILER
  #if __INTEL_COMPILER >= 1500
   #define GMX_ALIGNED_ALLOCATOR_HAVE_STD_ALIGN
  #endif
 #endif
#endif
/* If this macro is defined we should have all C++11 features required for the C++11
   implementation of AlignedAllocator. Otherwise use a fall-back implementation
   employing the GROMACS-specific functions snew_aligned and sfree_aligned           */
#ifdef GMX_ALIGNED_ALLOCATOR_HAVE_STD_ALIGN
// ---------------------------------------------------------------------------------------
#include <memory>

template <class T, unsigned alignment>
class AlignedAllocator : public std::allocator<T>
{
    public:
        //! integer data type used in measuring numbers of elements
        typedef size_t size_type;
        //! integer data type used in measuring differences between pointer values
        typedef std::ptrdiff_t difference_type;

        //! pointer data type corresponding to value_type
        typedef T* pointer;
        //! const pointer data type corresponding to value_type
        typedef const T* const_pointer;

        //! reference data type corresponding to value_type
        typedef T &reference;
        //! const reference data type corresponding to value_type
        typedef const T &const_reference;
        //! data type for which memory will be allocated
        typedef T value_type;

        //! adapt *this allocator to a different data type
        template <class U>
        struct rebind { typedef AlignedAllocator<U, alignment> other; };

        //! default constructor
        AlignedAllocator() throw() { }
        //! copy constructor
        AlignedAllocator(const AlignedAllocator &a) throw()
            : std::allocator<T>(a) { }

        //! copy constructor from different template specializations
        template <class U>
        AlignedAllocator(const AlignedAllocator<U, alignment> &) throw() { }

        //! the destructor
        ~AlignedAllocator() throw() { }

        //! check whether alignment is valid for data type T to avoid undefined behavior of std::align
        //! validity criteria:
        //! a) adopting the C++17 definition, alignment must be a power of two
        //! b) alignment must not be less than the standard alignment requirement of the data type//!
        //!
        //! \param[in]    ali   alignment
        //! \returns      true if the alignement is valid and false otherwise
        bool is_valid_alignment(const unsigned int ali) throw()
        {
            bool result = true;
            if (ali < alignof(value_type))
            {
                result = false;
            }
            else
            {
                unsigned int tmpi      = ali;
                unsigned int remainder = 0;
                while (remainder == 0 && tmpi > 1)
                {
                    remainder = tmpi % 2;
                    tmpi     /= 2;
                }
                if (tmpi != 1)
                {
                    // align is not a power of two
                    result = false;
                }
            }
            return result;
        }

        //! allocate memory on p
        //! \param[in]   n     number of elements for which memory is to be allocated
        pointer allocate(size_type n, AlignedAllocator<void, alignof(void*)>::const_pointer /*hint*/ = 0) // space for n Ts
        {
            unsigned ali = alignment;
            if (ali < alignof(void *))
            {
                ali = alignof(void *);
            }
            if (!is_valid_alignment(ali))
            {
                throw std::invalid_argument("AlignedAllocator::allocate(): invalid alignment");
            }

            // allocate a memory buffer large enough to contain
            // a) the shift of the pointer for alignment needed later when deallocating the memory again
            // b) all n elements of the actual data type to be stored
            // c) the maximum possible amount of the pointer shift for alignment of ali - 1
            const size_t buffer_size = sizeof(difference_type)
                + n * sizeof(value_type)
                + (ali - 1);
            char * p_char = new char[buffer_size];
            // shift the initial pointer position to preserve enough space for saving the alignment shift
            p_char += sizeof(difference_type);
            void * p_void = reinterpret_cast<void*>(p_char);

            // align the pointer by shifting it inside the buffer bounds
            // aligned_buffer_size will be decreased by std::align by the shift in bytes
            size_t aligned_buffer_size = buffer_size;
            p_void = std::align(ali, sizeof(T), p_void, aligned_buffer_size);
            if (p_void == nullptr)
            {
                // try to dispose the unalignable buffer
                try
                {
                    p_char -= sizeof(difference_type);
                    delete [] p_char;
                }
                catch (...)
                {
                    // silence this exception to avoid disguising the the original source of failure
                }
                throw std::bad_alloc();
            }

            // store the alignment shift at the buffer start, directly in front of the first data element
            p_char = reinterpret_cast<char*>(p_void) - sizeof(difference_type);
            difference_type * buffer_start = reinterpret_cast<difference_type*>(p_char);
            *buffer_start = buffer_size - aligned_buffer_size;

            return reinterpret_cast<pointer>(p_void);
        }
        //! deallocate the memory allocated previously on p
        //! \param[in]   p     pointer to the memory location to be deallocated
        void deallocate(pointer p, size_type /*n*/) // deallocate n Ts, don't destroy
        {
            char *p_work = reinterpret_cast<char*>(p);
            // move the pointer back to the unaligned position, i.e., the
            // start of the originally allocated buffer using the amount of the
            // shift previously stored at the buffer end
            p_work -= reinterpret_cast<difference_type &>(*(p_work - sizeof(difference_type)));
            p_work -= sizeof(difference_type);
            delete [] p_work;
        }

        //! construct a new object of type T at the memory location specified by p
        //! \param[in]   p     pointer to the target memory location
        //! \param[in]   val   reference to an object of type T, usually constructed in situ
        void construct(pointer p, const T &val) { new(p) T(val); }      // initialize *p by val
        //! destruct the the object of type T at the memory location specified by p
        //! \param[in]   p     pointer to the target memory location
        void destroy(pointer p) { p->~T(); } // destroy *p but don't deallocate

};

//! comparison operator == for two AlignedAllocator instantiations
//! \param[in]   a    left-hand side operand of the comparison
//! \param[in]   b   right-hand side operand of the comparison
template<class T1, class T2, unsigned alignment1, unsigned alignment2>
bool operator==(const AlignedAllocator<T1, alignment1> & /* a */, const AlignedAllocator<T2, alignment2> & /* b */) noexcept
{
    return (std::is_same<T1, T2>::value && alignment1 == alignment2);
}

//! comparison operator != for two AlignedAllocator instantiations
//! \param[in]   a    left-hand side operand of the comparison
//! \param[in]   b   right-hand side operand of the comparison
template<class T1, class T2, unsigned alignment1, unsigned alignment2>
bool operator!=(const AlignedAllocator<T1, alignment1> &a, const AlignedAllocator<T2, alignment2> &b) noexcept
{
    return !(a == b);
}


// ---------------------------------------------------------------------------------------
#else // end C++11 version, begin pre-C++11 version
// ---------------------------------------------------------------------------------------

 #include "gromacs/utility/smalloc.h"

template <class T, unsigned alignment>
class AlignedAllocator : public std::allocator<T>
{
    public:
        //! integer data type used in measuring numbers of elements
        typedef size_t size_type;
        //! integer data type used in measuring differences between pointer values
        typedef std::ptrdiff_t difference_type;

        //! pointer data type corresponding to value_type
        typedef T* pointer;
        //! const pointer data type corresponding to value_type
        typedef const T* const_pointer;

        //! reference data type corresponding to value_type
        typedef T &reference;
        //! const reference data type corresponding to value_type
        typedef const T &const_reference;
        //! data type for which memory will be allocated
        typedef T value_type;

        //! adapt *this allocator to a different data type
        template <class U>
        struct rebind { typedef AlignedAllocator<U, alignment> other; };

        //! default constructor
        AlignedAllocator() throw() { }
        //! copy constructor
        AlignedAllocator(const AlignedAllocator &a) throw()
            : std::allocator<T>(a) { }

        //! copy constructor from different template specializations
        template <class U>
        AlignedAllocator(const AlignedAllocator<U, alignment> &) throw() { }

        //! the destructor
        ~AlignedAllocator() throw() { }

        //! check whether alignment is valid to avoid undefined behavior of std::align
        //! validity criteria:
        //! a) adopting the C++17 definition, alignment must be a power of two
        //! b) alignment must not be less than the standard alignment requirement of the data type
        bool is_valid_alignment(const unsigned int ali) throw()
        {
            bool result = true;
            if (ali < alignof(value_type))
            {
                result = false;
            }
            else
            {
                unsigned int tmpi      = ali;
                unsigned int remainder = 0;
                while (remainder == 0 && tmpi > 1)
                {
                    remainder = tmpi % 2;
                    tmpi     /= 2;
                }
                if (tmpi != 1)
                {
                    // align is not a power of two
                    result = false;
                }
            }
            return result;
        }

        //! allocate memory on p
        //! \param[in]   n     number of elements for which memory is to be allocated
        pointer allocate(size_type n, AlignedAllocator<void, sizeof(void*)>::const_pointer /*hint*/ = 0) // space for n Ts
        {
            pointer  p;
            unsigned ali = alignment;
            // snew takes p as pointer by reference
            snew_aligned(p, n, ali);
            return p;
        }
        //! deallocate the memory allocated previously on p
        //! \param[in]   p     pointer to the memory location to be deallocated
        void deallocate(pointer p, size_type /*n*/) // deallocate n Ts, don't destroy
        {
            sfree_aligned(p);
        }

        //! construct a new object of type T at the memory location specified by p
        //! \param[in]   p     pointer to the target memory location
        //! \param[in]   val   reference to an object of type T, usually constructed in situ
        void construct(pointer p, const T &val) { new(p) T(val); }      // initialize *p by val
        //! destruct the the object of type T at the memory location specified by p
        //! \param[in]   p     pointer to the target memory location
        void destroy(pointer p) { p->~T(); } // destroy *p but don't deallocate
};

// comparison operator stubs adopted from AnBe's implementation

//! comparison operator == for two AlignedAllocator instantiations
template<class T1, class T2, unsigned alignment>
bool operator==(const AlignedAllocator<T1, alignment> &, const AlignedAllocator<T2, alignment> &) noexcept
{
    return true;
}

//! comparison operator != for two AlignedAllocator instantiations
template<class T1, class T2, unsigned alignment>
bool operator!=(const AlignedAllocator<T1, alignment> &, const AlignedAllocator<T2, alignment> &) noexcept
{
    return false;
}

// ---------------------------------------------------------------------------------------
#endif // end pre-C++11 version
// ---------------------------------------------------------------------------------------

#endif // end GMX_AlignedAlloc_H
