/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

void InMemorySerializer::doInt(int *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt64(gmx_int64_t *value)
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
        explicit Impl(const std::vector<char> &buffer)
            : buffer_(buffer), pos_(0)
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

        const std::vector<char> &buffer_;
        size_t                   pos_;
};

InMemoryDeserializer::InMemoryDeserializer(const std::vector<char> &buffer)
    : impl_(new Impl(buffer))
{
}

InMemoryDeserializer::~InMemoryDeserializer()
{
}

void InMemoryDeserializer::doBool(bool *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doUChar(unsigned char *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt(int *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt64(gmx_int64_t *value)
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

void InMemoryDeserializer::doString(std::string *value)
{
    impl_->doString(value);
}

} // namespace gmx
