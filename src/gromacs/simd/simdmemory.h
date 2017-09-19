/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief Declares gmx::SimdMemory and supporting traits, for
 * implementing safe SIMD-ized access to vectors of basic types.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 * \inlibraryapi
 */
// TODO eventually this will replace paddedvector.h
#ifndef GMX_SIMD_SIMDMEMORY_H
#define GMX_SIMD_SIMDMEMORY_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

// TODO This is a hack so that GMX_SIMD=None can compile and run this
// code. Find a better approach.
#if GMX_SIMD == 0 && !defined DOXYGEN
using SimdReal = real;
#  if GMX_DOUBLE
#define GMX_SIMD_DOUBLE_WIDTH 1;
#  else
#define GMX_SIMD_FLOAT_WIDTH 1;
#  endif
#endif

/*! \libinternal \brief Type trait to support SimdMemory. */
template <typename T>
struct SimdMemoryTypeTraits {};

/*! \libinternal \brief Type trait to support SimdMemory<real>. */
template <>
struct SimdMemoryTypeTraits<real> {
    //!@{ Trait values
    using ValueType   = real;
    using PointerType = real*;
    using BaseType    = real;
    using SimdType    = SimdReal;
    static constexpr std::size_t s_simdWidth = GMX_SIMD_REAL_WIDTH;
    static constexpr std::size_t s_packSize  = sizeof(ValueType) / sizeof(BaseType);
    //!@}
};

/*! \libinternal \brief Type trait to support SimdMemory<RVec>. */
template <>
struct SimdMemoryTypeTraits<RVec> {
    //!@{ Trait values
    using ValueType   = RVec;
    using PointerType = RVec*;
    using BaseType    = real;
    using SimdType    = BasicVector<SimdReal>;
    static constexpr std::size_t s_simdWidth = GMX_SIMD_REAL_WIDTH;
    static constexpr std::size_t s_packSize  = sizeof(ValueType) / sizeof(BaseType);
    //!@}
};

/*! \libinternal \brief Class to maintain contiguous memory intended for use in
 * SIMD-ized code.
 *
 * The memory is aligned and padded suitably for error-free SIMD-width
 * memory operations. A load of width \c s_simdWidth is safe from any
 * aligned or unaligned address from the storage_ within \c size()
 * elements of T. In particular, such loads are safe from both an
 * aligned load from addresses such as used by storage_.back() (used
 * in PME solve), and possibly DIM successive loads beginning from
 * unaligned addresses such as getUnpaddedArrayRef()[size()-1] (used
 * in bonded, constraint and update code). This is simply ensured by
 * making sure that the padded region has at least as many bytes as \c
 * s_simdWidth * DIM.
 *
 * \tparam T Must be a specialization of SimdMemoryTypeTraits,
 *           typically real or RVec, to provide the necessary
 *           knowledge about the types and prevent erroneous
 *           instantiation of SimdMemory.
 *
 * The external interface, e.g. size(), resizeWithPadding works in
 * terms of T objects, typically real or RVec. Code that needs to
 * access the memory does so by using the get methods that return
 * either an ArrayRef<T> or ArrayRef<SimdType> version.
 *
 * Internally, the memory is a std::vector<T>, so that tooling that
 * understands the behavior of std::vector will work correctly.
 */
template <typename T>
class SimdMemory
{
    private:
        //! Convenience typedef
        using SizeType = std::size_t;
        //! Convenience typedef, typically real or RVec.
        using ValueType = typename SimdMemoryTypeTraits<T>::ValueType;
        //! Convenience typedef, typically real* or RVec*
        using PointerType = typename SimdMemoryTypeTraits<T>::PointerType;
        //! Convenience typedef, typically real
        using BaseType = typename SimdMemoryTypeTraits<T>::BaseType;
        //! Convenience typedef, typically SimdReal
        using SimdType = typename SimdMemoryTypeTraits<T>::SimdType;
        //! Convenience typedef, typically vector<SimdReal, AlignedAllocator<SimdReal> >
        using StorageType = std::vector<T, AlignedAllocator<T> >;

        /*! \brief Helper function for computing the padded size.
         *
         * Internally, the storage is maintained in multiples of T. */
        SizeType computePaddedSize(const SizeType numElements)
        {
            // This is sometimes slightly more padding than strictly
            // required, but is much simpler code.
            SizeType           maxBaseTypeNeeded = (numElements * s_packSize + (DIM+1) * s_simdWidth);
            SizeType           numSimdTypeNeeded = maxBaseTypeNeeded / s_simdWidth;
            return numSimdTypeNeeded * s_simdWidth;
        }

    public:
        //! Convenience value, typically GMX_SIMD_REAL_WIDTH
        static constexpr std::size_t s_simdWidth = SimdMemoryTypeTraits<T>::s_simdWidth;
        //! Convenience value, typically 1 or DIM
        static constexpr std::size_t s_packSize = SimdMemoryTypeTraits<T>::s_packSize;

        //! Constructor.
        SimdMemory() :
            storage_(),
            unpadded_end_()
        {}
        //! Constructor.
        explicit SimdMemory(const SizeType numElements) :
            storage_(computePaddedSize(numElements)),
            unpadded_end_(storage_.data() + numElements)
        {}
        //! Constructor.
        explicit SimdMemory(SimdMemory const &o) :
            storage_(o.storage_),
            unpadded_end_(storage_.data() + o.size())
        {}
        //! Constructor.
        explicit SimdMemory(SimdMemory &&o) :
            storage_(std::move(o.storage_)),
            unpadded_end_(o.unpadded_end_)
        {}
        //! Reserves enough storage for \c numElements Ts, but does not initialize it.
        void reserveWithPadding(const SizeType numElements)
        {
            /* v.reserve(13) should allocate enough memory so that
               v.resize(13) does not reallocate. This means that the
               new extent should be large enough for the padded
               storage for a vector whose size is new_extent. */
            auto newSize        = computePaddedSize(numElements);
            auto oldSize        = storage_.size();
            storage_.reserve(newSize);
            unpadded_end_ = storage_.data() + oldSize;
        }
        //! Resizes the storage to fit \c numElements Ts and default-initializes it.
        void resizeWithPadding(const SizeType numElements)
        {
            auto newSize = computePaddedSize(numElements);
            storage_.resize(newSize);
            unpadded_end_ = storage_.data() + numElements;
        }
        //! Returns the size of the storage in T-sized objects.
        SizeType size() const { return unpadded_end_ - storage_.data(); }
        //! Returns whether the storage is empty.
        bool empty() const { return size() == 0; }
        //! Swaps the contents with those of \c x.
        void swap(SimdMemory &x)
        {
            std::swap(storage_, x.storage_);
            std::swap(unpadded_end_, x.unpadded_end_);
        }
        //! Clears the storage.
        void clear()
        {
            storage_.clear();
            unpadded_end_ = storage_.data();
        }
        //! Returns a view of the storage in Ts, without padding.
        ArrayRef<T> getUnpaddedArrayRef()
        {
            return arrayRefFromPointers<T>(storage_.data(), unpadded_end_);
        }
        //! Returns a const view of the storage in Ts, without padding.
        ConstArrayRef<T> getUnpaddedConstArrayRef() const
        {
            return constArrayRefFromPointers<T>(storage_.data(), unpadded_end_);
        }
        //! Returns a view of the storage in SimdType, with any necessary padding.
        ArrayRef<SimdType> getSimdArrayRef()
        {
            // TODO Check codegen, maybe we care about performance
            return arrayRefFromArray<SimdType>(reinterpret_cast<SimdType *>(storage_.data()), storage_.size() / s_simdWidth / s_packSize);
        }
        //! Returns a const view of the storage in SimdType, with any necessary padding.
        ConstArrayRef<SimdType> getSimdConstArrayRef() const
        {
            return constArrayRefFromArray<SimdType>(reinterpret_cast<const SimdType *>(storage_.data()), storage_.size() / s_simdWidth / s_packSize);
        }
        //! Returns a view of the padding region of the storage in BaseTypes.
        ArrayRef<BaseType> getPaddingArrayRef()
        {
            return arrayRefFromPointers<BaseType>(reinterpret_cast<BaseType*>(unpadded_end_), reinterpret_cast<BaseType*>(storage_.data() + storage_.size()));
        }
        //! Returns a const view of the padding region of the storage in BaseTypes.
        ConstArrayRef<BaseType> getConstPaddingArrayRef() const
        {
            return constArrayRefFromPointers<BaseType>(reinterpret_cast<const BaseType*>(unpadded_end_), reinterpret_cast<const BaseType*>(storage_.data() + storage_.size()));
        }
        /*! \brief Initializes the memory in the padding region to \c value.
         *
         * Avoid calling this routine from inner loops, because the
         * simplistic padding arrangement means that in general the
         * padded region will exceed the range of memory actually
         * accessed. */
        void setPaddingRegionTo(const BaseType value)
        {
            for (auto &element : getPaddingArrayRef())
            {
                element = value;
            }
        }

    private:
        //! Contains the actual memory allocated
        StorageType storage_;
        //! Points to the end of the range of Ts for which storage was requested, ie. before any padding.
        PointerType unpadded_end_;
};

} // namespace gmx

#endif
