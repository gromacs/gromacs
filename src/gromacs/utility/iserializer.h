/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Declares a generic serialization interface that supports both directions.
 *
 * \todo Generalize and transfer serialization functionality used in
 *       mrc density file header serialization to here.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ISERIALIZER_H
#define GMX_UTILITY_ISERIALIZER_H

#include <cstddef>

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/*! \libinternal
 * \brief Interface for types that convert standard data types into a
 * form suitable for storage or transfer.
 *
 * Different implementations could suit MPI, file I/O, or in-memory
 * conversion. */
class ISerializer
{
public:
    virtual ~ISerializer() {}
    /*! \brief Returns whether the serializer is reading or
     * writing, because details like memory management vary
     * accordingly. */
    virtual bool reading() const = 0;
    //! \brief Serialize values of different types.
    ///@{
    virtual void doBool(bool* value)                    = 0;
    virtual void doUChar(unsigned char* value)          = 0;
    virtual void doChar(char* value)                    = 0;
    virtual void doUShort(unsigned short* value)        = 0;
    virtual void doInt(int* value)                      = 0;
    virtual void doInt32(int32_t* value)                = 0;
    virtual void doInt64(int64_t* value)                = 0;
    virtual void doFloat(float* value)                  = 0;
    virtual void doDouble(double* value)                = 0;
    virtual void doReal(real* value)                    = 0;
    virtual void doIvec(ivec* value)                    = 0;
    virtual void doRvec(rvec* value)                    = 0;
    virtual void doString(std::string* value)           = 0;
    virtual void doOpaque(char* data, std::size_t size) = 0;
    ///@}

    //! \brief Serialize arrays of values of different types.
    ///@{
    void doBoolArray(bool* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doBool(&(values[i]));
        }
    }
    // Char, UChar and RVec have vector specializations that can be
    // used instead of the default looping.
    virtual void doCharArray(char* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doChar(&(values[i]));
        }
    }
    virtual void doUCharArray(unsigned char* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doUChar(&(values[i]));
        }
    }
    void doUShortArray(unsigned short* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doUShort(&(values[i]));
        }
    }
    void doIntArray(int* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doInt(&(values[i]));
        }
    }
    void doInt32Array(int32_t* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doInt32(&(values[i]));
        }
    }
    void doInt64Array(int64_t* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doInt64(&(values[i]));
        }
    }
    void doFloatArray(float* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doFloat(&(values[i]));
        }
    }
    void doDoubleArray(double* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doDouble(&(values[i]));
        }
    }
    void doRealArray(real* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doReal(&(values[i]));
        }
    }
    void doIvecArray(ivec* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doIvec(&(values[i]));
        }
    }
    virtual void doRvecArray(rvec* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doRvec(&(values[i]));
        }
    }
    ///@}

    //! Serialize enum value with underlying type int
    template<typename EnumType>
    void doEnumAsInt(EnumType* enumValue)
    {
        static_assert(std::is_same<std::underlying_type_t<EnumType>, int>::value,
                      "Only enums with underlying type int are supported.");
        auto castedValue = static_cast<int>(*enumValue);
        doInt(&castedValue);
        *enumValue = static_cast<EnumType>(castedValue);
    }

    //! Serialize array of enum values with underlying type.
    template<typename EnumType>
    void doEnumArrayAsInt(EnumType* values, int elements)
    {
        for (int i = 0; i < elements; i++)
        {
            doEnumAsInt<EnumType>(&(values[i]));
        }
    }
};

} // namespace gmx

#endif
