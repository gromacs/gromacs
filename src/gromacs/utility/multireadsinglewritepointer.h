/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_UTILITY_MULTIREADSINGLEWRITEPOINTER_H
#define GMX_UTILITY_MULTIREADSINGLEWRITEPOINTER_H

#include <memory>

template <class T>
class multiReadSingleWritePointer {
public:
    multiReadSingleWritePointer(T *ptr = nullptr) :
            ptr_(ptr),
            writing_(false)
    {}
    multiReadSingleWritePointer(std::shared_ptr<T> &ptr) :
            ptr_(ptr),
            writing_(false)
    {}
    multiReadSingleWritePointer(std::unique_ptr<T> ptr) :
            ptr_(std::move(ptr)),
            writing_(false)
    {}
    multiReadSingleWritePointer(const multiReadSingleWritePointer<T> &ptr)
    {
        if (ptr.checkWriting()) {
            // no copying while writing
            throw std::runtime_error("Data is being written to.");
        } else {
            this->ptr_ = ptr.ptr_;
            this->writing_ = ptr.writing_;
        }
    }
    multiReadSingleWritePointer& operator = (const multiReadSingleWritePointer<T> &ptr)
    {
        if (ptr.checkWriting()) {
            // no assignment while writing
            throw std::runtime_error("Data is being written to.");
        } else {
            this->ptr_ = ptr.ptr_;
            this->writing_ = ptr.writing_;
        }
        return *this;
    }

    const T* get() const {
        if (checkWriting())
        {
            // no reading while writing
            throw std::runtime_error("Data is being written to.");
        }
        return ptr_.get();
    }

    const T& operator*() const {
        return *this->get();
    }

    const T* operator->() const {
        return this->get();
    }

    explicit operator bool() const {
        return bool(this->get());
    }

    std::shared_ptr<T> write() {
        if (checkWriting())
        {
            // can't have two writing handles out
            throw std::runtime_error("Data is being written to.");
        }
        if (not ptr_.unique()) {
            // Data has reader handles out
            // std::cout << "Making new object." << std::endl;
            ptr_ = std::make_shared<T>(*ptr_);
        }
        writing_ = true;
        return ptr_;
    }
private:
    std::shared_ptr<T> ptr_;
    mutable bool writing_{};
    inline bool checkWriting() const
    {
        if (writing_ and ptr_.unique()) {
            writing_ = false;
        }
        return writing_;
    }
};

template < class T, class U >
bool operator == (const multiReadSingleWritePointer<T>& lhs,
                  const multiReadSingleWritePointer<U>& rhs ) noexcept
{
    return lhs.get() == rhs.get();
}

template < class T, class U >
bool operator != (const multiReadSingleWritePointer<T>& lhs,
                  const multiReadSingleWritePointer<U>& rhs ) noexcept
{
    return !(lhs == rhs);
}

template< class T >
bool operator == (const multiReadSingleWritePointer<T>& lhs,
                  std::nullptr_t rhs ) noexcept
{
    return !bool(lhs);
}

template< class T >
bool operator == (std::nullptr_t lhs,
                  const multiReadSingleWritePointer<T>& rhs ) noexcept
{
    return !bool(rhs);
}

template< class T >
bool operator != (const multiReadSingleWritePointer<T>& lhs,
                  std::nullptr_t rhs ) noexcept
{
    return !bool(lhs);
}

template< class T >
bool operator != (std::nullptr_t lhs,
                  const multiReadSingleWritePointer<T>& rhs ) noexcept
{
    return !bool(rhs);
}

#endif //GMX_UTILITY_MULTIREADSINGLEWRITEPOINTER_H
