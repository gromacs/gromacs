/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include <boost/scoped_ptr.hpp>

namespace gmx
{

/*! \brief
 * Macro to declare a class non-copyable and non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \ingroup module_utility
 */
#define GMX_DISALLOW_COPY_AND_ASSIGN(ClassName) \
    private: \
        ClassName &operator=(const ClassName &); \
        ClassName(const ClassName &)
/*! \brief
 * Macro to declare a class non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \ingroup module_utility
 */
#define GMX_DISALLOW_ASSIGN(ClassName) \
    private: \
        ClassName &operator=(const ClassName &)

/*! \brief
 * Helper class to manage a pointer to a private implementation class.
 *
 * This helper provides the following benefits (all but the last could also be
 * achieved with boost::scoped_ptr):
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
 * Intended use:
 * \code
   // In exampleclass.h
   class ExampleClass
   {
       public:
           ExampleClass();
           ~ExampleClass(); // Must not be defined inline

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
 * \inlibraryapi
 * \ingroup module_utility
 */
template <class Impl>
class PrivateImplPointer
{
    public:
        //! Initialize with the given implementation class.
        explicit PrivateImplPointer(Impl *ptr) : ptr_(ptr) {}
        ~PrivateImplPointer() {}

        /*! \brief
         * Sets a new implementation class and destructs the previous one.
         *
         * Needed, e.g., to implement assignable classes.
         */
        void reset(Impl *ptr) { ptr_.reset(ptr); }
        //! Access the raw pointer.
        Impl *get() { return ptr_.get(); }
        //! Access the implementation class as with a raw pointer.
        Impl *operator->() { return ptr_.get(); }
        //! Access the implementation class as with a raw pointer.
        Impl &operator*() { return *ptr_; }
        //! Access the raw pointer.
        const Impl *get() const { return ptr_.get(); }
        //! Access the implementation class as with a raw pointer.
        const Impl *operator->() const { return ptr_.get(); }
        //! Access the implementation class as with a raw pointer.
        const Impl &operator*() const { return *ptr_; }

    private:
        boost::scoped_ptr<Impl> ptr_;

        // Copy and assign disabled by the scoped_ptr member.
};

} // namespace gmx

#endif
