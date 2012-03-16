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
 *
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
/*! \file
 * \brief
 * Declares common utility classes and macros.
 *
 * This header contains helpers used to implement classes in the library.
 * It is installed, because the helpers are used in installed headers, but
 * typically users of the library should not need to be aware of these helpers.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_COMMON_H
#define GMX_UTILITY_COMMON_H

#include <boost/scoped_ptr.hpp>

/*! \libinternal \brief
 * Macro to declare a class non-copyable and non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \inlibraryapi
 */
#define GMX_DISALLOW_COPY_AND_ASSIGN(ClassName) \
    private: \
        ClassName(const ClassName &); \
        ClassName &operator =(const ClassName &)
/*! \libinternal \brief
 * Macro to declare a class non-assignable.
 *
 * For consistency, should appear last in the class declaration.
 *
 * \inlibraryapi
 */
#define GMX_DISALLOW_ASSIGN(ClassName) \
    private: \
        ClassName &operator =(const ClassName &)

namespace gmx
{

/*! \libinternal \brief
 * Helper class to manage a pointer to a private implementation class.
 *
 * This helper provides the following benefits (the first two could also be
 * achieved with boost::scoped_ptr):
 *  - Automatic memory management: the implementation pointer is freed in
 *    the destructor automatically.  If the destructor is not declared or is
 *    defined inline in the header file, a compilation error occurs instead
 *    of a memory leak or undefined behavior.
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
 * \endcode
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

        void reset(Impl *ptr) { ptr_.reset(ptr); }
        Impl *get() { return ptr_.get(); }
        Impl *operator->() { return ptr_.get(); }
        Impl &operator*() { return *ptr_; }
        const Impl *get() const { return ptr_.get(); }
        const Impl *operator->() const { return ptr_.get(); }
        const Impl &operator*() const { return *ptr_; }

    private:
        boost::scoped_ptr<Impl> ptr_;

        // Copy and assign disabled by the scoped_ptr member.
};

} // namespace gmx

#endif
