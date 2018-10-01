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
 * Declares classes for building the data structures in keyvaluetree.h.
 *
 * These are separate from the data structures to enforce clear separation of
 * the APIs, and to make the data structure immutable after construction.
 *
 * For the main use case described in \ref page_mdmodules, they are mainly
 * used internally, but currently gmx::KeyValueTreeObjectBuilder (and
 * everything it references) is exposed for more complex transforms through
 * gmx::IKeyValueTreeTransformRules.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_KEYVALUETREEBUILDER_H
#define GMX_UTILITY_KEYVALUETREEBUILDER_H

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/variant.h"

namespace gmx
{

class KeyValueTreeArrayBuilder;
class KeyValueTreeObjectBuilder;

/*! \libinternal \brief
 * Root builder for creating trees that have an object at the root.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeBuilder
{
    public:
        //! Returns a builder for the root object.
        KeyValueTreeObjectBuilder rootObject();

        /*! \brief
         * Builds the final object.
         *
         * The builder should not be accessed after this call.
         */
        KeyValueTreeObject build() { return std::move(root_); }

    private:
        /*! \brief
         * Helper function for other builders to create values of certain type.
         */
        template <typename T>
        static KeyValueTreeValue createValue(const T &value)
        {
            return KeyValueTreeValue(Variant::create<T>(value));
        }
        /*! \brief
         * Helper function for other builders to create default-constructed
         * values.
         */
        template <typename T>
        static KeyValueTreeValue createValue()
        {
            return KeyValueTreeValue(Variant::create<T>(T()));
        }

        KeyValueTreeObject root_;

        //! For access to createValue() methods.
        friend class KeyValueTreeObjectArrayBuilder;
        //! For access to createValue() methods.
        friend class KeyValueTreeObjectBuilder;
        //! For access to createValue() methods.
        template <typename T>
        friend class KeyValueTreeUniformArrayBuilder;
};

/*! \libinternal \brief
 * Builder for KeyValueTreeValue objects.
 *
 * This builder can be constructed directly and can create self-standing
 * KeyValueTreeValue objects.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeValueBuilder
{
    public:
        //! Assigns a scalar value of certain type.
        template <typename T>
        void setValue(const T &value)
        {
            value_ = Variant::create<T>(value);
        }
        //! Assigns a Variant value to the built value.
        void setVariantValue(Variant &&value)
        {
            value_ = std::move(value);
        }
        /*! \brief
         * Returns an object builder for building an object into this value.
         *
         * Any method call in this value builder invalidates the returned
         * builder.
         */
        KeyValueTreeObjectBuilder createObject();
        /*! \brief
         * Returns an array builder for building an array into this value.
         *
         * Any method call in this value builder invalidates the returned
         * builder.
         */
        KeyValueTreeArrayBuilder createArray();

        /*! \brief
         * Builds the final value.
         *
         * The builder should not be accessed after this call.
         */
        KeyValueTreeValue build() { return KeyValueTreeValue(std::move(value_)); }

    private:
        Variant value_;
};

class KeyValueTreeArrayBuilderBase
{
    protected:
        //! Creates an array builder for populating given array object.
        explicit KeyValueTreeArrayBuilderBase(KeyValueTreeArray *array)
            : array_(array)
        {
        }

        //! Appends a raw Variant value to the array.
        KeyValueTreeValue &addRawValue(Variant &&value)
        {
            KeyValueTreeValueBuilder builder;
            builder.setVariantValue(std::move(value));
            array_->values_.push_back(builder.build());
            return array_->values_.back();
        }
        //! Appends a raw KeyValueTreeValue to the array.
        KeyValueTreeValue &addRawValue(KeyValueTreeValue &&value)
        {
            array_->values_.push_back(std::move(value));
            return array_->values_.back();
        }

    private:
        KeyValueTreeArray *array_;
};

class KeyValueTreeArrayBuilder : public KeyValueTreeArrayBuilderBase
{
    public:
        using KeyValueTreeArrayBuilderBase::addRawValue;

    private:
        explicit KeyValueTreeArrayBuilder(KeyValueTreeArray *array)
            : KeyValueTreeArrayBuilderBase(array)
        {
        }

        friend class KeyValueTreeObjectBuilder;
        friend class KeyValueTreeValueBuilder;
};

/*! \libinternal \brief
 * Builder for KeyValueTreeArray objects where all elements are of type `T`.
 *
 * The builder does not own the array being constructed, but instead holds a
 * reference to an object within a tree rooted in KeyValueTreeBuilder or
 * KeyValueTreeValueBuilder.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <typename T>
class KeyValueTreeUniformArrayBuilder : public KeyValueTreeArrayBuilderBase
{
    public:
        //! Appends a value to the array.
        void addValue(const T &value)
        {
            addRawValue(KeyValueTreeBuilder::createValue<T>(value));
        }

    private:
        explicit KeyValueTreeUniformArrayBuilder(KeyValueTreeArray *array)
            : KeyValueTreeArrayBuilderBase(array)
        {
        }

        friend class KeyValueTreeObjectBuilder;
};

/*! \libinternal \brief
 * Builder for KeyValueTreeArray objects where all elements are
 * KeyValueTreeObject objects.
 *
 * The builder does not own the array being constructed, but instead holds a
 * reference to an object within a tree rooted in KeyValueTreeBuilder or
 * KeyValueTreeValueBuilder.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeObjectArrayBuilder : public KeyValueTreeArrayBuilderBase
{
    public:
        /*! \brief
         * Appends an object to the array.
         *
         * The object is created empty and can be built using the returned
         * builder.
         */
        KeyValueTreeObjectBuilder addObject();

    private:
        explicit KeyValueTreeObjectArrayBuilder(KeyValueTreeArray *array)
            : KeyValueTreeArrayBuilderBase(array)
        {
        }

        friend class KeyValueTreeObjectBuilder;
};

/*! \libinternal \brief
 * Builder for KeyValueTreeObject objects.
 *
 * The builder does not own the object being constructed, but instead holds a
 * reference to an object within a tree rooted in KeyValueTreeBuilder or
 * KeyValueTreeValueBuilder.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeObjectBuilder
{
    public:
        //! Adds a property with given key from a KeyValueTreeValue.
        void addRawValue(const std::string &key, KeyValueTreeValue &&value)
        {
            addProperty(key, std::move(value));
        }
        //! Adds a property with given key from a Variant value.
        void addRawValue(const std::string &key, Variant &&value)
        {
            addProperty(key, KeyValueTreeValue(std::move(value)));
        }
        //! Adds a scalar property with given key, type, and value.
        template <typename T>
        void addValue(const std::string &key, const T &value)
        {
            addRawValue(key, KeyValueTreeBuilder::createValue<T>(value));
        }
        /*! \brief
         * Adds an object-valued property with given key.
         *
         * The object is created empty and can be built using the returned
         * builder.
         */
        KeyValueTreeObjectBuilder addObject(const std::string &key)
        {
            auto iter = addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeObject>());
            return KeyValueTreeObjectBuilder(&iter->second);
        }
        /*! \brief
         * Adds a generic array-valued property with given key.
         *
         * The array is created empty and can be built using the returned
         * builder.
         */
        KeyValueTreeArrayBuilder addArray(const std::string &key)
        {
            auto iter = addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeArrayBuilder(&iter->second.asArray());
        }
        /*! \brief
         * Adds an array-valued property with uniform value types with given
         * key.
         *
         * \tparam T  Type for all values in the array.
         *
         * The array is created empty and can be built using the returned
         * builder.
         */
        template <typename T>
        KeyValueTreeUniformArrayBuilder<T> addUniformArray(const std::string &key)
        {
            auto iter = addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeUniformArrayBuilder<T>(&iter->second.asArray());
        }
        /*! \brief
         * Adds an array-valued property with uniform value types with given
         * key and values.
         *
         * \tparam T  Type for all values in the array.
         *
         * The array is created to contain the values from `values`.
         */
        template <typename T>
        void addUniformArray(const std::string &key, std::initializer_list<T> values)
        {
            auto builder = addUniformArray<T>(key);
            for (const auto &value : values)
            {
                builder.addValue(value);
            }
        }
        /*! \brief
         * Adds an array-valued property with objects in the array with given
         * key.
         *
         * The array is created empty and can be built using the returned
         * builder.
         */
        KeyValueTreeObjectArrayBuilder addObjectArray(const std::string &key)
        {
            auto iter = addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeObjectArrayBuilder(&iter->second.asArray());
        }

        //! Whether a property with given key exists.
        bool keyExists(const std::string &key) const { return object_->keyExists(key); }
        //! Returns value for a given key.
        const KeyValueTreeValue &operator[](const std::string &key) const
        {
            return (*object_)[key];
        }
        //! Returns an object builder for an existing object.
        KeyValueTreeObjectBuilder getObjectBuilder(const std::string &key)
        {
            GMX_ASSERT(keyExists(key), "Requested non-existent value");
            GMX_ASSERT((*this)[key].isObject(), "Accessing non-object value as object");
            return KeyValueTreeObjectBuilder(&object_->valueMap_.at(key).asObject());
        }

        /*! \brief
         * Returns whether the given object shares any keys with \p this.
         */
        bool objectHasDistinctProperties(const KeyValueTreeObject &obj) const
        {
            return object_->hasDistinctProperties(obj);
        }
        /*! \brief
         * Merges properties from a given object to `this`.
         *
         * The objects should not share any keys, i.e.,
         * objectHasDistinctProperties() should return `true`.
         */
        void mergeObject(KeyValueTreeObject &&obj)
        {
            GMX_ASSERT(objectHasDistinctProperties(obj),
                       "Trying to merge overlapping object");
            for (auto &prop : obj.valueMap_)
            {
                addRawValue(prop.first, std::move(prop.second));
            }
        }

    private:
        explicit KeyValueTreeObjectBuilder(KeyValueTreeObject *object)
            : object_(object)
        {
        }
        explicit KeyValueTreeObjectBuilder(KeyValueTreeValue *value)
            : object_(&value->asObject())
        {
        }

        std::map<std::string, KeyValueTreeValue>::iterator
        addProperty(const std::string &key, KeyValueTreeValue &&value)
        {
            GMX_RELEASE_ASSERT(!keyExists(key), "Duplicate key value");
            object_->values_.reserve(object_->values_.size() + 1);
            auto iter = object_->valueMap_.insert(std::make_pair(key, std::move(value))).first;
            object_->values_.push_back(KeyValueTreeProperty(iter));
            return iter;
        }

        KeyValueTreeObject *object_;

        friend class KeyValueTreeBuilder;
        friend class KeyValueTreeValueBuilder;
        friend class KeyValueTreeObjectArrayBuilder;
};

/********************************************************************
 * Inline functions that could not be declared within the classes
 */

inline KeyValueTreeObjectBuilder KeyValueTreeBuilder::rootObject()
{
    return KeyValueTreeObjectBuilder(&root_);
}

inline KeyValueTreeObjectBuilder KeyValueTreeValueBuilder::createObject()
{
    value_ = Variant::create<KeyValueTreeObject>(KeyValueTreeObject());
    return KeyValueTreeObjectBuilder(&value_.castRef<KeyValueTreeObject>());
}

inline KeyValueTreeArrayBuilder KeyValueTreeValueBuilder::createArray()
{
    value_ = Variant::create<KeyValueTreeArray>(KeyValueTreeArray());
    return KeyValueTreeArrayBuilder(&value_.castRef<KeyValueTreeArray>());
}

inline KeyValueTreeObjectBuilder KeyValueTreeObjectArrayBuilder::addObject()
{
    auto &value = addRawValue(KeyValueTreeBuilder::createValue<KeyValueTreeObject>());
    return KeyValueTreeObjectBuilder(&value);
}

} // namespace gmx

#endif
