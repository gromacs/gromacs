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
/*! \file
 * \brief
 * Declares gmx::Allocator template whose allocation functionality is
 * configured both by type of object allocated and a policy class that
 * configures the necessary matching malloc and free operation.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ALLOCATOR_H
#define GMX_UTILITY_ALLOCATOR_H

#include <cstddef>

#include <memory>
#include <new>

#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

/*! \libinternal \brief Policy-based memory allocator.
 *
 *  \tparam T                 Type of objects to allocate
 *  \tparam AllocationPolicy  Policy of (matching) allocation and deallocation functions.
 *
 * This class can be used for the optional allocator template
 * parameter in standard library containers. It must be configured
 * with both the type of object to allocate, and an AllocationPolicy
 * which effectively wraps a matching pair of malloc and free
 * functions. This permits implementing a family of related allocators
 * e.g. with SIMD alignment, GPU host-side page locking, or perhaps
 * both, in a way that preserves a common programming interface and
 * duplicates minimal code.
 *
 * AllocationPolicy is used as a base class, so that if
 * AllocationPolicy is stateless, then the empty base optimization
 * will ensure that Allocation is also stateless, and objects made
 * with the Allocator will incur no size penalty. (Embedding an
 * AllocationPolicy object incurs a size penalty always, even if the
 * object is empty.) Normally a stateless allocator will be used.
 *
 * However, an AllocationPolicy with state might be desirable for
 * simplifying writing code that needs to allocate suitably for a
 * transfer to a GPU. That code needs to specify an Allocator that can
 * do the right job, which can be stateless. However, if we have code
 * that will not know until run time whether a GPU transfer will
 * occur, then the allocator needs to be aware of the state.  That
 * will increase the size of a container that uses the stateful
 * allocator.
 *
 * \throws std::bad_alloc Instead of a GROMACS exception object, we
 * throw the standard one on allocation failures to make it as
 * compatible as possible with the errors expected by code using the
 * standard library containers.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <class T, typename AllocationPolicy>
class Allocator : public AllocationPolicy
{
    public:
        // The standard library specification for a custom allocator
        // requires these typedefs, with this capitalization/underscoring.
        typedef T              value_type;      //!< Type of allocated elements
        typedef T             &reference;       //!< Reference to allocated elements
        typedef const T       &const_reference; //!< Constant reference to allocated elements
        typedef T *            pointer;         //!< Pointer to allocated elements
        typedef const T *      const_pointer;   //!< Constant pointer to allocated elements
        typedef std::size_t    size_type;       //!< Integer type to use for size of objects
        typedef std::ptrdiff_t difference_type; //!< Type to hold differences between pointers

        // This typedef is required by GROMACS for testing and assertions
        typedef AllocationPolicy allocation_policy; //!< Type of the AllocationPolicy

        /*! \libinternal \brief Standard-required typedef to use allocator with different class.
         *
         *  \tparam U new class
         *
         *  This is used for things like std::list where the size of each link
         *  is larger than the class stored in the link.
         *
         *  Required by the specification for an allocator.
         */
        template <class U>
        struct rebind
        {
            typedef Allocator<U, AllocationPolicy> other; //!< Align class U with our alignment
        };

        /*! \brief Templated copy constructor
         *
         * This template constructor cannot be auto-generated, and is
         * normally unused, except e.g. MSVC2015 standard library uses
         * it in debug mode, presumably to implement some checks.
         */
        template <class U>
        explicit Allocator(const Allocator<U, AllocationPolicy> &) {}

        /*! \brief Constructor
         *
         * No constructor can be auto-generated in the presence of any
         * user-defined constructor, but we want the default constructor.
         */
        Allocator() = default;

        /*! \brief Constructor to accept an AllocationPolicy.
         *
         * This is useful for AllocationPolicies with state.
         */
        Allocator(const AllocationPolicy &p) : AllocationPolicy(p) {}

        /*! \brief Return address of an object
         *
         *  \param r Reference to object of type T
         *  \return Pointer to T memory
         */
        pointer
        address(reference r) const { return &r; }

        /*! \brief Return address of a const object
         *
         *  \param r Const reference to object of type T
         *  \return Pointer to T memory
         */
        const_pointer
        address(const_reference r) const { return &r; }

        /*! \brief Do the actual memory allocation
         *
         *  \param n    Number of elements of type T to allocate. n can be
         *              0 bytes, which will return a non-null properly aligned
         *              and padded pointer that should not be used.
         *  \param hint Optional value returned from previous call to allocate.
         *              For now this is not used.
         *  \return Pointer to allocated memory
         *
         *  \throws std::bad_alloc if the allocation fails.
         */
        pointer
        allocate(std::size_t n, typename std::allocator<void>::const_pointer gmx_unused hint = nullptr)
        {
            void *p = AllocationPolicy::malloc(n*sizeof(T));

            if (p == nullptr)
            {
                throw std::bad_alloc();
            }
            else
            {
                return static_cast<pointer>(p);
            }
        }

        /*! \brief Release memory
         *
         * \param p  Pointer to previously allocated memory returned from allocate()
         * \param n  number of objects previously passed to allocate()
         */
        void
        deallocate(pointer p, std::size_t gmx_unused n)
        {
            AllocationPolicy::free(p);
        }

        //! Return the policy object for this allocator.
        AllocationPolicy getPolicy() const
        {
            return *this;
        }

        /*! \brief Construct an object without allocating memory
         *
         * \tparam Args  Variable-length list of types for constructor args
         * \param p      Adress of memory where to construct object
         * \param args   Variable-length list of arguments to constructor
         */
        template<class ... Args>
        void
        construct(pointer p, Args && ... args) { ::new((void *)p)T(std::forward<Args>(args) ...); }

        /*! \brief Call the destructor of object without releasing memory
         *
         * \param p  Address of memory where to destroy object
         */
        void
        destroy(pointer p) { p->~value_type(); }

        /*! \brief Return largest number of objects that can be allocated
         *
         * This will be set such that the number of objects T multiplied by
         * the size of each object is the largest value that can be represented
         * by size_type.
         */
        std::size_t
        max_size() const { return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T); }

        /*! \brief Return true if two allocators are identical
         *
         * This is a member function of the left-hand-side allocator.
         */
        template<class T2>
        bool
        operator==(const Allocator<T2, AllocationPolicy> &) const { return std::is_same<T, T2>::value; }

        /*! \brief Return true if two allocators are different
         *
         * \param rhs Other allocator.
         *
         * This is a member function of the left-hand-side allocator.
         */
        bool
        operator!=(const Allocator &rhs) const { return !operator==(rhs); }
};

}      // namespace gmx

#endif
