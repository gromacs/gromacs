/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2019, by the GROMACS development team, led by
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
 * Declares gmx::IOptionValueStore.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_IVALUESTORE_H
#define GMX_OPTIONS_IVALUESTORE_H

namespace gmx
{

template<typename T>
class ArrayRef;

/*! \internal
 * \brief
 * Represents the final storage location of option values.
 *
 * \todo
 * Try to make this more like a write-only interface, getting rid of the need
 * to access the stored values through this interface.  That would simplify
 * things.
 *
 * \ingroup module_options
 */
template<typename T>
class IOptionValueStore
{
public:
    virtual ~IOptionValueStore() {}

    //! Returns the number of values stored so far.
    virtual int valueCount() = 0;
    //! Returns a reference to the actual values.
    virtual ArrayRef<T> values() = 0;
    //! Removes all stored values.
    virtual void clear() = 0;
    //! Reserves memory for additional `count` entries.
    virtual void reserve(size_t count) = 0;
    //! Appends a value to the store.
    virtual void append(const T& value) = 0;
};

} // namespace gmx

#endif
