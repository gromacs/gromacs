/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "inmemoryserializer.h"

#include <algorithm>
#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

template <typename T>
class CharBuffer
{
    public:
        static const size_t ValueSize = sizeof(T);

        explicit CharBuffer(T value)
        {
            u.v = value;
        }
        explicit CharBuffer(const char buffer[])
        {
            std::copy(buffer, buffer + ValueSize, u.c);
        }

        T value() const { return u.v; }

        void appendTo(std::vector<char> *buffer)
        {
            buffer->insert(buffer->end(), u.c, u.c + ValueSize);
        }

    private:
        union {
            char c[ValueSize];
            T    v;
        } u;
};

}   // namespace

/********************************************************************
 * InMemorySerializer
 */

class InMemorySerializer::Impl
{
    public:
        template <typename T>
        void doValue(T value)
        {
            CharBuffer<T>(value).appendTo(&buffer_);
        }
        void doString(const std::string &value)
        {
            doValue<size_t>(value.size());
            buffer_.insert(buffer_.end(), value.begin(), value.end());
        }

        std::vector<char> buffer_;
};

InMemorySerializer::InMemorySerializer()
    : impl_(new Impl)
{
}

InMemorySerializer::~InMemorySerializer()
{
}

std::vector<char> InMemorySerializer::finishAndGetBuffer()
{
    return std::move(impl_->buffer_);
}

void InMemorySerializer::doBool(bool *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doUChar(unsigned char *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doChar(char *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doUShort(unsigned short *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt(int *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt32(int32_t *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt64(int64_t *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doFloat(float *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doDouble(double *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doReal(real *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doRvec(rvec *value)
{
    for (int d = 0; d < DIM; d++)
    {
        doReal(&(*value)[d]);
    }
}

void InMemorySerializer::doIvec(ivec *value)
{
    for (int d = 0; d < DIM; d++)
    {
        doInt(&(*value)[d]);
    }
}

void InMemorySerializer::doString(std::string *value)
{
    impl_->doString(*value);
}

/********************************************************************
 * InMemoryDeserializer
 */

class InMemoryDeserializer::Impl
{
    public:
        explicit Impl(ArrayRef<const char> buffer, bool sourceIsDouble)
            : buffer_(buffer), sourceIsDouble_(sourceIsDouble), pos_(0)
        {
        }

        template <typename T>
        void doValue(T *value)
        {
            *value = CharBuffer<T>(&buffer_[pos_]).value();
            pos_  += CharBuffer<T>::ValueSize;
        }
        void doString(std::string *value)
        {
            size_t size;
            doValue<size_t>(&size);
            *value = std::string(&buffer_[pos_], size);
            pos_  += size;
        }

        ArrayRef<const char>     buffer_;
        bool                     sourceIsDouble_;
        size_t                   pos_;
};

InMemoryDeserializer::InMemoryDeserializer(ArrayRef<const char> buffer, bool sourceIsDouble)
    : impl_(new Impl(buffer, sourceIsDouble))
{
}

InMemoryDeserializer::~InMemoryDeserializer()
{
}

bool InMemoryDeserializer::sourceIsDouble() const
{
    return impl_->sourceIsDouble_;
}

void InMemoryDeserializer::doBool(bool *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doUChar(unsigned char *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doChar(char *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doUShort(unsigned short *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt(int *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt32(int32_t *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt64(int64_t *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doFloat(float *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doDouble(double *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doReal(real *value)
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

void InMemoryDeserializer::doRvec(rvec *value)
{
    for (int d = 0; d < DIM; d++)
    {
        doReal(&(*value)[d]);
    }
}

void InMemoryDeserializer::doIvec(ivec *value)
{
    for (int d = 0; d < DIM; d++)
    {
        doInt(&(*value)[d]);
    }
}

void InMemoryDeserializer::doString(std::string *value)
{
    impl_->doString(value);
}

} // namespace gmx
