/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
template<class T, typename AllocationPolicy>
class Allocator : public AllocationPolicy
{
public:
    // The standard library specification for a custom allocator
    // requires this typedef, with this capitalization/underscoring.
    typedef T value_type; //!< Type of allocated elements

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
    Allocator(const AllocationPolicy& p) : AllocationPolicy(p) {}

    /*! \brief Do the actual memory allocation
     *
     *  \param n    Number of elements of type T to allocate. n can be
     *              0 bytes, which will return a non-null properly aligned
     *              and padded pointer that should not be used.
     *  \return Pointer to allocated memory
     *
     *  \throws std::bad_alloc if the allocation fails.
     */
    value_type* allocate(std::size_t n)
    {
        void* p = AllocationPolicy::malloc(n * sizeof(T));

        if (p == nullptr)
        {
            throw std::bad_alloc();
        }
        else
        {
            return static_cast<value_type*>(p);
        }
    }

    /*! \brief Release memory
     *
     * \param p  Pointer to previously allocated memory returned from allocate()
     * \param n  number of objects previously passed to allocate()
     */
    void deallocate(value_type* p, std::size_t gmx_unused n) { AllocationPolicy::free(p); }

    /*! \brief Return true if two allocators are identical
     *
     * This is a member function of the left-hand-side allocator.
     * Always true for stateless policies. Has to be defined in the policy for stateful policies.
     * FUTURE: Can be removed with C++17 (is_always_equal)
     *
     * \todo Use std::is_empty_v when CUDA 11 is a requirement.
     */
    template<class T2, class A = AllocationPolicy, typename = std::enable_if_t<std::is_empty<A>::value>>
    bool operator==(const Allocator<T2, AllocationPolicy>& /*unused*/) const
    {
        return true;
    }

    /*! \brief Return true if two allocators are different
     *
     * \param rhs Other allocator.
     *
     * This is a member function of the left-hand-side allocator.
     */
    template<class T2>
    bool operator!=(const Allocator<T2, AllocationPolicy>& rhs) const
    {
        return !(*this == rhs);
    }
};

} // namespace gmx

#endif
