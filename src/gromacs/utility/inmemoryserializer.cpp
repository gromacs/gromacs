/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include <memory>
#include <vector>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * InMemorySerializer
 */

class InMemorySerializer::Impl
{
    public:
        Impl() : size_(0)
        {
        }

        template <typename T>
        void doValue(T value)
        {
            // There can be some extra here, but should not matter.
            size_t requiredSpace = sizeof(T) + alignof(T);
            buffer_.resize(size_ + requiredSpace);
            void  *ptr            = &buffer_[size_];
            size_t remainingSpace = requiredSpace;
            if (!std::align(alignof(T), sizeof(T), ptr, remainingSpace))
            {
                GMX_RELEASE_ASSERT(false, "Alignment should always be possible");
            }
            size_ += (requiredSpace - remainingSpace);
            *reinterpret_cast<T *>(ptr) = value;
            size_ += sizeof(T);
        }
        void doString(const std::string &value)
        {
            doValue<size_t>(value.size());
            buffer_.resize(size_ + value.size());
            std::copy(value.begin(), value.end(), &buffer_[size_]);
            size_ += value.size();
        }

        std::vector<char> buffer_;
        size_t            size_;
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
    impl_->buffer_.resize(impl_->size_);
    return std::move(impl_->buffer_);
}

void InMemorySerializer::doUChar(unsigned char *value)
{
    impl_->doValue(*value);
}

void InMemorySerializer::doInt(int *value)
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
            void  *ptr            = const_cast<char *>(&buffer_[pos_]);
            size_t remainingSpace = buffer_.size() - pos_;
            if (!std::align(alignof(T), sizeof(T), ptr, remainingSpace))
            {
                GMX_RELEASE_ASSERT(false, "Alignment should always be possible");
            }
            pos_   = (buffer_.size() - remainingSpace);
            *value = *reinterpret_cast<const T *>(ptr);
            pos_  += sizeof(T);
        }
        void doString(std::string *value)
        {
            size_t size;
            doValue<size_t>(&size);
            *value = std::string(&buffer_[pos_], &buffer_[pos_ + size]);
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

void InMemoryDeserializer::doUChar(unsigned char *value)
{
    impl_->doValue(value);
}

void InMemoryDeserializer::doInt(int *value)
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
