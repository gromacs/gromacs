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
/*! \internal \file
 * \brief
 * Defines gmx::ISerializer implementation for in-memory serialization.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/inmemoryserializer.h"

#include "config.h"

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

template<typename T>
class CharBuffer
{
public:
    static constexpr size_t ValueSize = sizeof(T);

    explicit CharBuffer(T value) { u.v = value; }
    explicit CharBuffer(const char buffer[]) { std::copy(buffer, buffer + ValueSize, u.c); }

    T value() const { return u.v; }

    void appendTo(std::vector<char>* buffer)
    {
        buffer->insert(buffer->end(), u.c, u.c + ValueSize);
    }

private:
    union
    {
        char c[ValueSize];
        T    v;
    } u;
};

//! Return \c value with the byte order swapped.
template<typename T>
T swapEndian(const T& value)
{
    union
    {
        T                           value_;
        std::array<char, sizeof(T)> valueAsCharArray_;
    } endianessSwappedValue;

    endianessSwappedValue.value_ = value;
    int hiByte                   = sizeof(T) - 1;
    for (int loByte = 0; hiByte > loByte; loByte++, hiByte--)
    {
        std::swap(endianessSwappedValue.valueAsCharArray_[loByte],
                  endianessSwappedValue.valueAsCharArray_[hiByte]);
    }

    return endianessSwappedValue.value_;
}

/*! \brief Change the host-dependent endian settings to either Swap or DoNotSwap.
 *
 * \param endianSwapBehavior input swap behavior, might depend on host.
 *
 * \return Host-independent setting, either Swap or DoNotSwap. */
EndianSwapBehavior setEndianSwapBehaviorFromHost(EndianSwapBehavior endianSwapBehavior)
{
    if (endianSwapBehavior == EndianSwapBehavior::SwapIfHostIsBigEndian)
    {
        return GMX_INTEGER_BIG_ENDIAN ? EndianSwapBehavior::Swap : EndianSwapBehavior::DoNotSwap;
    }
    else if (endianSwapBehavior == EndianSwapBehavior::SwapIfHostIsLittleEndian)
    {
        return GMX_INTEGER_BIG_ENDIAN ? EndianSwapBehavior::DoNotSwap : EndianSwapBehavior::Swap;
    }
    else
    {
        return endianSwapBehavior;
    }
}

} // namespace

/********************************************************************
 * InMemorySerializer
 */

class InMemorySerializer::Impl
{
public:
    Impl(EndianSwapBehavior endianSwapBehavior) :
        endianSwapBehavior_(setEndianSwapBehaviorFromHost(endianSwapBehavior))
    {
    }

    template<typename T>
    void doValue(T value)
    {
        if (endianSwapBehavior_ == EndianSwapBehavior::Swap)
        {
            CharBuffer<T>(swapEndian(value)).appendTo(&buffer_);
        }
        else
        {
            CharBuffer<T>(value).appendTo(&buffer_);
        }
    }
    void doString(const std::string& value)
    {
        doValue<uint64_t>(value.size());
        buffer_.insert(buffer_.end(), value.begin(), value.end());
    }
    void doOpaque(const char* data, std::size_t size)
    {
        buffer_.insert(buffer_.end(), data, data + size);
    }

    std::vector<char>  buffer_;
    EndianSwapBehavior endianSwapBehavior_;
};

InMemorySerializer::InMemorySerializer(EndianSwapBehavior endianSwapBehavior) :
    impl_(new Impl(endianSwapBehavior))
{
}

InMemorySerializer::~InMemorySerializer() = default;

std::vector<char> InMemorySerializer::finishAndGetBuffer()
{
    return std::move(impl_->buffer_);
}

void InMemorySerializer::doBool(bool* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doUChar(unsigned char* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doChar(char* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doUShort(unsigned short* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt(int* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt32(int32_t* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt64(int64_t* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doFloat(float* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doDouble(double* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doReal(real* value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doRvec(rvec* value)
{
    for (real& v : *value)
    {
        doReal(&v);
    }
}

void InMemorySerializer::doIvec(ivec* value)
{
    for (int& v : *value)
    {
        doInt(&v);
    }
}

void InMemorySerializer::doString(std::string* value)
{
    impl_->doString(*value);
}

void InMemorySerializer::doOpaque(char* data, std::size_t size)
{
    impl_->doOpaque(data, size);
}

/********************************************************************
 * InMemoryDeserializer
 */

class InMemoryDeserializer::Impl
{
public:
    explicit Impl(ArrayRef<const char> buffer, bool sourceIsDouble, EndianSwapBehavior endianSwapBehavior) :
        buffer_(buffer),
        sourceIsDouble_(sourceIsDouble),
        pos_(0),
        endianSwapBehavior_(setEndianSwapBehaviorFromHost(endianSwapBehavior))
    {
    }

    template<typename T>
    void doValue(T* value)
    {
        if (endianSwapBehavior_ == EndianSwapBehavior::Swap)
        {
            *value = swapEndian(CharBuffer<T>(&buffer_[pos_]).value());
        }
        else
        {
            *value = CharBuffer<T>(&buffer_[pos_]).value();
        }
        pos_ += CharBuffer<T>::ValueSize;
    }
    void doString(std::string* value)
    {
        uint64_t size = 0;
        doValue<uint64_t>(&size);
        *value = std::string(&buffer_[pos_], size);
        pos_ += size;
    }
    void doOpaque(char* data, std::size_t size)
    {
        std::copy(&buffer_[pos_], &buffer_[pos_ + size], data);
        pos_ += size;
    }

    ArrayRef<const char> buffer_;
    bool                 sourceIsDouble_;
    size_t               pos_;
    EndianSwapBehavior   endianSwapBehavior_;
};

InMemoryDeserializer::InMemoryDeserializer(ArrayRef<const char> buffer,
                                           bool                 sourceIsDouble,
                                           EndianSwapBehavior   endianSwapBehavior) :
    impl_(new Impl(buffer, sourceIsDouble, endianSwapBehavior))
{
}

InMemoryDeserializer::~InMemoryDeserializer() = default;

bool InMemoryDeserializer::sourceIsDouble() const
{
    return impl_->sourceIsDouble_;
}

void InMemoryDeserializer::doBool(bool* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doUChar(unsigned char* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doChar(char* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doUShort(unsigned short* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt(int* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt32(int32_t* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt64(int64_t* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doFloat(float* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doDouble(double* value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doReal(real* value)
{
    if (sourceIsDouble())
    {
        double temp = 0.0;
        doDouble(&temp);
        *value = temp;
    }
    else
    {
        float temp = 0.0;
        doFloat(&temp);
        *value = temp;
    }
}

void InMemoryDeserializer::doRvec(rvec* value)
{
    for (real& v : *value)
    {
        doReal(&v);
    }
}

void InMemoryDeserializer::doIvec(ivec* value)
{
    for (int& v : *value)
    {
        doInt(&v);
    }
}

void InMemoryDeserializer::doString(std::string* value)
{
    impl_->doString(value);
}

void InMemoryDeserializer::doOpaque(char* data, std::size_t size)
{
    impl_->doOpaque(data, size);
}

} // namespace gmx
