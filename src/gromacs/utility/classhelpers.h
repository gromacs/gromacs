/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * Declares common utility classes and macros.
 *
 * This header contains helpers used to implement classes in the library.
 * It is installed, because the helpers are used in installed headers, but
 * typically users of the library should not need to be aware of these helpers.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_CLASSHELPERS_H
#define GMX_UTILITY_CLASSHELPERS_H

#include <memory>

namespace gmx
{

#ifdef DOXYGEN
/*! \brief
 * Macro to declare a class non-copyable and non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \ingroup module_utility
 */
#    define GMX_DISALLOW_COPY_AND_ASSIGN(ClassName)
#else
#    define GMX_DISALLOW_COPY_AND_ASSIGN(ClassName)      \
        ClassName& operator=(const ClassName&) = delete; \
        ClassName(const ClassName&)            = delete
#endif
/*! \brief
 * Macro to declare a class non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \ingroup module_utility
 */
#define GMX_DISALLOW_ASSIGN(ClassName) ClassName& operator=(const ClassName&) = delete

// clang-format off
#ifdef DOXYGEN
/*! \brief
 * Macro to declare default constructors
 *
 * Intended for copyable interfaces or bases classes which require to create custom
 * destructor (e.g. protected or virtual) but need the default constructors.
 *
 * \ingroup module_utility
 */
#    define GMX_DEFAULT_CONSTRUCTORS(ClassName)
#else
#    define GMX_DEFAULT_CONSTRUCTORS(ClassName)                                                                           \
        ClassName() = default;                                                                                            \
        ClassName& operator=(const ClassName&) = default; /* NOLINT(misc-macro-parentheses,bugprone-macro-parentheses) */ \
        ClassName(const ClassName&) = default;                                                                            \
        ClassName& operator=(ClassName&&) = default; /* NOLINT(misc-macro-parentheses,bugprone-macro-parentheses) */      \
        ClassName(ClassName&&) = default /* NOLINT(misc-macro-parentheses,bugprone-macro-parentheses) */
#endif
//clang-format on

/*! \brief
 * Helper class to manage a pointer to a private implementation class.
 *
 * This helper provides the following benefits (all but the last could also be
 * achieved with std::unique_ptr):
 *  - Automatic memory management: the implementation pointer is freed in
 *    the destructor automatically.  If the destructor is not declared or is
 *    defined inline in the header file, a compilation error occurs instead
 *    of a memory leak or undefined behavior.
 *  - Exception safety in constructors: the implementation pointer is freed
 *    correctly even if the constructor of the containing class throws after
 *    the implementation class is constructed.
 *  - Copy and/or assignment is automatically disallowed if explicit copy
 *    constructor and/or assignment operator is not provided.
 *  - Compiler helps to manage const-correctness: in const methods, it is not
 *    possible to change the implementation class.
 *
 * Move construction and assignment are also disallowed, but can be enabled by
 * providing explicit move constructor and/or assignment.
 *
 * Intended use:
 * \code
   // In exampleclass.h
   class ExampleClass
   {
       public:
           ExampleClass();
           ~ExampleClass(); // Must be defined, must not be defined inline

           // <...>

       private:
           class Impl;

           PrivateImplPointer<Impl> impl_;
   };

   // In exampleclass.cpp

   // <definition of ExampleClass::Impl>

   ExampleClass::ExampleClass()
       : impl_(new Impl)
   {
   }

   ExampleClass::~ExampleClass()
   {
   }
   \endcode
 *
 * Note that ExampleClass::~ExampleClass cannot be declared inline (or
 * generated by the compiler) because the implementation of impl_
 * requires that ExampleClass::Impl be known in size, whereas it has
 * only been forward declared. Only the translation unit where
 * ExampleClass::Impl is declared can define the destructor for
 * ExampleClass (which may be defaulted).
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template<class Impl>
class PrivateImplPointer
{
public:
    //! Allow implicit initialization from nullptr to support comparison.
    PrivateImplPointer(std::nullptr_t) : ptr_(nullptr) {}
    //! Initialize with the given implementation class.
    explicit PrivateImplPointer(Impl* ptr) : ptr_(ptr) {}
    //! \cond
    // Explicitly declared to work around MSVC problems.
    PrivateImplPointer(PrivateImplPointer&& other) noexcept : ptr_(std::move(other.ptr_)) {}
    PrivateImplPointer& operator=(PrivateImplPointer&& other) noexcept
    {
        ptr_ = std::move(other.ptr_);
        return *this;
    }
    //! \endcond

    /*! \brief
     * Sets a new implementation class and destructs the previous one.
     *
     * Needed, e.g., to implement lazily initializable or copy-assignable
     * classes.
     */
    void reset(Impl* ptr) { ptr_.reset(ptr); }
    //! Access the raw pointer.
    Impl* get() { return ptr_.get(); }
    //! Access the implementation class as with a raw pointer.
    Impl* operator->() { return ptr_.get(); }
    //! Access the implementation class as with a raw pointer.
    Impl& operator*() { return *ptr_; }
    //! Access the implementation class as with a raw pointer.
    const Impl* operator->() const { return ptr_.get(); }
    //! Access the implementation class as with a raw pointer.
    const Impl& operator*() const { return *ptr_; }

    //! Allows testing whether the implementation is initialized.
    explicit operator bool() const { return ptr_ != nullptr; }

    //! Tests for equality (mainly useful against nullptr).
    bool operator==(const PrivateImplPointer& other) const { return ptr_ == other.ptr_; }
    //! Tests for inequality (mainly useful against nullptr).
    bool operator!=(const PrivateImplPointer& other) const { return ptr_ != other.ptr_; }

private:
    std::unique_ptr<Impl> ptr_;

    // Copy construction and assignment disabled by the unique_ptr member.
};

} // namespace gmx

#endif
