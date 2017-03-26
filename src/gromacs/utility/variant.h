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
/*! \libinternal \file
 * \brief
 * Declares gmx::Variant.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_VARIANT_H
#define GMX_UTILITY_VARIANT_H

#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <utility>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \libinternal \brief
 * Represents a dynamically typed value of an arbitrary type.
 *
 * To create a variant, either initialize it as empty, or with the create()
 * method (or the equivalent constructor, if the type parameter can be deduced
 * and is clear to the reader from the context).
 *
 * To query the type of the contents in the variant, use isEmpty(), type(), and
 * isType().
 *
 * To access the value, you need to know the type as a compile-time constant
 * (e.g., through branching based on isType()), and then use cast() or
 * tryCast().
 *
 * Methods in this class do not throw unless otherwise indicated.
 *
 * This provides essentially the same functionality as boost::any.
 *
 * \ingroup module_utility
 */
class Variant
{
    public:
        /*! \brief
         * Creates a variant that holds the given value.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * This method allows explicitly specifying the template argument,
         * contrary to the templated constructor.
         */
        template <typename T>
        static Variant create(const T &value) { return Variant(value); }
        /*! \brief
         * Creates a variant that holds the given value.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * In addition to allowing specifying the template argument, this
         * method avoids copying when move-construction is possible.
         */
        template <typename T>
        static Variant create(T &&value) { return Variant(std::move(value)); }

        //! Creates an empty variant value.
        Variant() {}
        /*! \brief
         * Creates a variant that holds the given value.
         *
         * \throws std::bad_alloc if out of memory.
         */
        template <typename T>
        explicit Variant(const T &value)
            : content_(new Content<typename std::decay<T>::type>(value))
        {
        }
        /*! \brief
         * Creates a variant that holds the given value.
         *
         * \throws std::bad_alloc if out of memory.
         */
        template <typename T>
        explicit Variant(T &&value)
            : content_(new Content<typename std::decay<T>::type>(std::move(value)))
        {
        }
        /*! \brief
         * Creates a deep copy of a variant.
         *
         * \throws std::bad_alloc if out of memory.
         */
        Variant(const Variant &other) : content_(other.cloneContent()) {}
        //! Move-constructs a variant.
        Variant(Variant &&other) noexcept : content_(std::move(other.content_)) {}
        /*! \brief
         * Assigns the variant.
         *
         * \throws std::bad_alloc if out of memory.
         */
        Variant &operator=(const Variant &other)
        {
            content_.reset(other.cloneContent());
            return *this;
        }
        //! Move-assigns the variant.
        Variant &operator=(Variant &&other) noexcept
        {
            content_ = std::move(other.content_);
            return *this;
        }

        //! Whether any value is stored.
        bool isEmpty() const { return content_ == nullptr; }
        //! Returns the dynamic type of the value that is currently stored.
        std::type_index type() const
        {
            const std::type_info &info
                = !isEmpty() ? content_->typeInfo() : typeid(void);
            return std::type_index(info);
        }
        //! Returns whether the type stored matches the template parameter.
        template <typename T>
        bool isType() const
        {
            return !isEmpty() && content_->typeInfo() == typeid(T);
        }

        /*! \brief
         * Tries to get the value as the given type.
         *
         * \tparam T  Type to get.
         * \returns Pointer to the value, or nullptr if the type does not match
         *     the stored value.
         */
        template <typename T>
        const T *tryCast() const
        {
            return isType<T>() ? &static_cast<Content<T> *>(content_.get())->value_ : nullptr;
        }
        /*! \brief
         * Gets the value when the type is known.
         *
         * \tparam T  Type to get (which must match what the variant stores).
         *
         * Asserts if the variant is empty or does not contain the requested type.
         */
        template <typename T>
        const T &cast() const
        {
            const T *value = tryCast<T>();
            GMX_RELEASE_ASSERT(value != nullptr, "Cast to incorrect type");
            return *value;
        }
        /*! \brief
         * Tries to get the value as the given type as a non-const pointer.
         *
         * \tparam T  Type to get.
         * \returns Pointer to the value, or nullptr if the type does not match
         *     the stored value.
         *
         * This method allows modifying the value in-place, which is useful
         * with more complicated data structures.
         */
        template <typename T>
        T *tryCastRef()
        {
            return isType<T>() ? &static_cast<Content<T> *>(content_.get())->value_ : nullptr;
        }
        /*! \brief
         * Gets the value when the type is known as a modifiable reference.
         *
         * \tparam T  Type to get (which must match what the variant stores).
         *
         * Asserts if the variant is empty or does not contain the requested type.
         */
        template <typename T>
        T &castRef()
        {
            T *value = tryCastRef<T>();
            GMX_RELEASE_ASSERT(value != nullptr, "Cast to incorrect type");
            return *value;
        }

    private:
        class IContent
        {
            public:
                virtual ~IContent() {}
                virtual const std::type_info &typeInfo() const = 0;
                virtual IContent *clone() const                = 0;
        };

        template <typename T>
        class Content : public IContent
        {
            public:
                explicit Content(const T &value) : value_(value) {}
                explicit Content(T &&value) : value_(std::move(value)) {}

                virtual const std::type_info &typeInfo() const { return typeid(T); }
                virtual IContent *clone() const { return new Content(value_); }

                T value_;
        };

        //! Creates a deep copy of the content.
        IContent *cloneContent() const
        {
            return content_ != nullptr ? content_->clone() : nullptr;
        }

        std::unique_ptr<IContent> content_;
};

//! \cond libapi
/*! \brief
 * Converts a Variant value to a string.
 *
 * As the name suggests, only some types of "simple" values (such as int) are
 * supported.  Asserts for unsupported types.
 *
 * \ingroup module_utility
 */
std::string simpleValueToString(const Variant &value);
//! \endcond

} // namespace gmx

#endif
