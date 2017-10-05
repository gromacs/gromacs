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
/*! \file
 * \brief
 * Declares implementations for IOptionValueStore.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_VALUESTORE_H
#define GMX_OPTIONS_VALUESTORE_H

#include <vector>

#include "gromacs/options/ivaluestore.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

template <typename T>
class OptionValueStorePlain : public IOptionValueStore<T>
{
    public:
        OptionValueStorePlain(T *store, int *storeCount, int initialCount)
            : count_(initialCount), store_(store), storeCount_(storeCount)
        {
        }

        virtual int valueCount() { return count_; }
        virtual ArrayRef<T> values() { return arrayRefFromArray(store_, count_); }
        virtual void clear()
        {
            count_ = 0;
            if (storeCount_ != nullptr)
            {
                *storeCount_ = count_;
            }
        }
        virtual void reserve(size_t /*count*/)
        {
        }
        virtual void append(const T &value)
        {
            store_[count_] = value;
            ++count_;
            if (storeCount_ != nullptr)
            {
                *storeCount_ = count_;
            }
        }

    private:
        int                          count_;
        T                           *store_;
        int                         *storeCount_;
};

template <typename T>
class OptionValueStoreVector : public IOptionValueStore<T>
{
    public:
        // cppcheck-suppress uninitMemberVar
        explicit OptionValueStoreVector(std::vector<T> *store) : store_(store) {}

        virtual int valueCount() { return static_cast<int>(store_->size()); }
        virtual ArrayRef<T> values() { return *store_; }
        virtual void clear() { store_->clear(); }
        virtual void reserve(size_t count)
        {
            store_->reserve(store_->size() + count);
        }
        virtual void append(const T &value)
        {
            store_->push_back(value);
        }

    private:
        std::vector<T> *store_;
};

// Specialization that works around std::vector<bool> specialities.
template <>
// cppcheck-suppress noConstructor
class OptionValueStoreVector<bool> : public IOptionValueStore<bool>
{
    public:
        explicit OptionValueStoreVector(std::vector<bool> *store) : store_(store)
        {
        }

        virtual int valueCount() { return static_cast<int>(store_->size()); }
        virtual ArrayRef<bool> values()
        {
            return arrayRefFromArray(reinterpret_cast<bool *>(boolStore_.data()),
                                     boolStore_.size());
        }
        virtual void clear()
        {
            boolStore_.clear();
            store_->clear();
        }
        virtual void reserve(size_t count)
        {
            boolStore_.reserve(boolStore_.size() + count);
            store_->reserve(store_->size() + count);
        }
        virtual void append(const bool &value)
        {
            boolStore_.push_back({value});
            store_->push_back(value);
        }

    private:
        struct Bool
        {
            bool value;
        };

        std::vector<Bool>  boolStore_;
        std::vector<bool> *store_;
};

/*! \internal
 * \brief
 * Value storage that does not store anywhere.
 *
 * This is needed because even though the values are not stored anywhere, the
 * code still expects to access them later through valueCount() and values().
 *
 * \ingroup module_options
 */
template <typename T>
class OptionValueStoreNull : public IOptionValueStore<T>
{
    public:
        OptionValueStoreNull() : store_(&vector_) {}

        virtual int valueCount() { return store_.valueCount(); }
        virtual ArrayRef<T> values() { return store_.values(); }
        virtual void clear() { store_.clear(); }
        virtual void reserve(size_t count) { store_.reserve(count); }
        virtual void append(const T &value) { store_.append(value); }

    private:
        std::vector<T>            vector_;
        OptionValueStoreVector<T> store_;
};

} // namespace gmx

#endif
