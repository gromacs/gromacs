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
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_KEYVALUETREEBUILDER_H
#define GMX_UTILITY_KEYVALUETREEBUILDER_H

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

class KeyValueTreeBuilder
{
    public:
        KeyValueTreeObjectBuilder rootObject();

        KeyValueTreeObject build() { return std::move(root_); }

    private:
        template <typename T>
        static KeyValueTreeValue createValue(const T &value)
        {
            return KeyValueTreeValue(Variant::create<T>(value));
        }
        template <typename T>
        static KeyValueTreeValue createValue()
        {
            return KeyValueTreeValue(Variant::create<T>(T()));
        }

        KeyValueTreeObject root_;

        friend class KeyValueTreeObjectArrayBuilder;
        friend class KeyValueTreeObjectBuilder;
        template <typename T>
        friend class KeyValueTreeUniformArrayBuilder;
};

class KeyValueTreeValueBuilder
{
    public:
        template <typename T>
        void setValue(const T &value)
        {
            value_ = Variant::create<T>(value);
        }
        void setVariantValue(Variant &&value)
        {
            value_ = std::move(value);
        }
        KeyValueTreeObjectBuilder createObject();
        KeyValueTreeArrayBuilder createArray();

        KeyValueTreeValue build() { return KeyValueTreeValue(std::move(value_)); }

    private:
        Variant value_;
};

class KeyValueTreeArrayBuilderBase
{
    protected:
        explicit KeyValueTreeArrayBuilderBase(KeyValueTreeArray *array)
            : array_(array)
        {
        }

        KeyValueTreeValue &addRawValue(Variant &&value)
        {
            KeyValueTreeValueBuilder builder;
            builder.setVariantValue(std::move(value));
            array_->values_.push_back(builder.build());
            return array_->values_.back();
        }
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

template <typename T>
class KeyValueTreeUniformArrayBuilder : public KeyValueTreeArrayBuilderBase
{
    public:
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

class KeyValueTreeObjectArrayBuilder : public KeyValueTreeArrayBuilderBase
{
    public:
        KeyValueTreeObjectBuilder addObject();

    private:
        explicit KeyValueTreeObjectArrayBuilder(KeyValueTreeArray *array)
            : KeyValueTreeArrayBuilderBase(array)
        {
        }

        friend class KeyValueTreeObjectBuilder;
};

class KeyValueTreeObjectBuilder
{
    public:
        void addRawValue(const std::string &key, KeyValueTreeValue &&value)
        {
            object_->addProperty(key, std::move(value));
        }
        void addRawValue(const std::string &key, Variant &&value)
        {
            object_->addProperty(key, KeyValueTreeValue(std::move(value)));
        }
        template <typename T>
        void addValue(const std::string &key, const T &value)
        {
            addRawValue(key, KeyValueTreeBuilder::createValue<T>(value));
        }
        KeyValueTreeObjectBuilder addObject(const std::string &key)
        {
            auto iter = object_->addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeObject>());
            return KeyValueTreeObjectBuilder(&iter->second);
        }
        KeyValueTreeArrayBuilder addArray(const std::string &key)
        {
            auto iter = object_->addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeArrayBuilder(&iter->second.asArray());
        }
        template <typename T>
        KeyValueTreeUniformArrayBuilder<T> addUniformArray(const std::string &key)
        {
            auto iter = object_->addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeUniformArrayBuilder<T>(&iter->second.asArray());
        }
        KeyValueTreeObjectArrayBuilder addObjectArray(const std::string &key)
        {
            auto iter = object_->addProperty(key, KeyValueTreeBuilder::createValue<KeyValueTreeArray>());
            return KeyValueTreeObjectArrayBuilder(&iter->second.asArray());
        }

        bool objectHasDistinctProperties(const KeyValueTreeObject &obj) const
        {
            return object_->hasDistinctProperties(obj);
        }
        void mergeObject(KeyValueTreeValue &&value)
        {
            mergeObject(std::move(value.asObject()));
        }
        void mergeObject(KeyValueTreeObject &&obj)
        {
            for (auto &prop : obj.valueMap_)
            {
                addRawValue(prop.first, std::move(prop.second));
            }
        }

        bool keyExists(const std::string &key) const { return object_->keyExists(key); }
        const KeyValueTreeValue &getValue(const std::string &key) const
        {
            GMX_ASSERT(keyExists(key), "Requested non-existent value");
            return (*object_)[key];
        }
        KeyValueTreeObjectBuilder getObject(const std::string &key)
        {
            return KeyValueTreeObjectBuilder(&(*object_)[key].asObject());
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
