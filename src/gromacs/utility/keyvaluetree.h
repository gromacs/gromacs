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
 * Declares a data structure for JSON-like structured key-value mapping.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_KEYVALUETREE_H
#define GMX_UTILITY_KEYVALUETREE_H

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/real.h"
#include "gromacs/utility/variant.h"

namespace gmx
{

class KeyValueTreeArray;
class KeyValueTreeObject;
class TextWriter;

/*! \libinternal \brief
 * Identifies an entry in a key-value tree.
 *
 * This class is mainly an internal utility within the key-value tree
 * implementation, but it is exposed on the API level where string-based
 * specification of a location in the tree is necessary.  Constructors are not
 * explicit to allow passing a simple string in contexts where a
 * KeyValueTreePath is expected.
 *
 * The string specifying a location should start with a `/`, followed by the
 * names of the properties separated by `/`.  For example, `/a/b/c` specifies
 * property `c` in an object that is the value of `b` in an object that is the
 * value of `a` at the root of the tree.
 * Currently, there is no support for specifying paths to values within arrays
 * (since none of the places where this is used implement array handling,
 * either).
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreePath
{
    public:
        //! Creates an empty path (corresponds to the root object).
        KeyValueTreePath() = default;
        //! Creates a path from given string representation.
        KeyValueTreePath(const char *path);
        //! Creates a path from given string representation.
        KeyValueTreePath(const std::string &path);

        //! Adds another element to the path, making it a child of the old path.
        void append(const std::string &key) { path_.push_back(key); }
        //! Removes the last element in the path, making it the parent path.
        void pop_back() { return path_.pop_back(); }
        //! Removes and returns the last element in the path.
        std::string pop_last()
        {
            std::string result = std::move(path_.back());
            path_.pop_back();
            return result;
        }

        //! Whether the path is empty (pointing to the root object).
        bool empty() const { return path_.empty(); }
        //! Returns the number of elements (=nesting level) in the path.
        size_t size() const { return path_.size(); }
        //! Returns the i'th path element.
        const std::string &operator[](int i) const { return path_[i]; }
        //! Returns all the path elements.
        const std::vector<std::string> &elements() const { return path_; }

        //! Formats the path as a string for display.
        std::string toString() const;

    private:
        std::vector<std::string> path_;
};

class KeyValueTreeValue
{
    public:
        bool isArray() const;
        bool isObject() const;
        template <typename T>
        bool isType() const { return value_.isType<T>(); }
        std::type_index type() const { return value_.type(); }

        KeyValueTreeArray &asArray();
        KeyValueTreeObject       &asObject();
        const KeyValueTreeArray  &asArray() const;
        const KeyValueTreeObject &asObject() const;
        template <typename T>
        const T                  &cast() const { return value_.cast<T>(); }

        const Variant            &asVariant() const { return value_; }

    private:
        explicit KeyValueTreeValue(Variant &&value) : value_(std::move(value)) {}

        Variant             value_;

        friend class KeyValueTreeBuilder;
        friend class KeyValueTreeObjectBuilder;
        friend class KeyValueTreeValueBuilder;
};

class KeyValueTreeArray
{
    public:
        bool isObjectArray() const
        {
            return std::all_of(values_.begin(), values_.end(),
                               std::mem_fn(&KeyValueTreeValue::isObject));
        }

        const std::vector<KeyValueTreeValue> &values() const { return values_; }

    private:
        std::vector<KeyValueTreeValue> values_;

        friend class KeyValueTreeArrayBuilderBase;
};

class KeyValueTreeProperty
{
    public:
        const std::string &key() const { return value_->first; }
        const KeyValueTreeValue &value() const { return value_->second; }

    private:
        typedef std::map<std::string, KeyValueTreeValue>::const_iterator
            IteratorType;

        explicit KeyValueTreeProperty(IteratorType value) : value_(value) {}

        IteratorType value_;

        friend class KeyValueTreeObject;
};

class KeyValueTreeObject
{
    public:
        KeyValueTreeObject() = default;
        KeyValueTreeObject(const KeyValueTreeObject &other)
        {
            for (const auto &value : other.values_)
            {
                auto iter = valueMap_.insert(std::make_pair(value.key(), value.value())).first;
                values_.push_back(KeyValueTreeProperty(iter));
            }
        }
        KeyValueTreeObject &operator=(KeyValueTreeObject &other)
        {
            KeyValueTreeObject tmp(other);
            std::swap(tmp.valueMap_, valueMap_);
            std::swap(tmp.values_, values_);
            return *this;
        }
        KeyValueTreeObject(KeyValueTreeObject &&)            = default;
        KeyValueTreeObject &operator=(KeyValueTreeObject &&) = default;

        const std::vector<KeyValueTreeProperty> &properties() const { return values_; }

        bool keyExists(const std::string &key) const
        {
            return valueMap_.find(key) != valueMap_.end();
        }
        const KeyValueTreeValue &operator[](const std::string &key) const
        {
            return valueMap_.at(key);
        }

        bool hasDistinctProperties(const KeyValueTreeObject &obj) const;
        void writeUsing(TextWriter *writer) const;

    private:
        KeyValueTreeValue &operator[](const std::string &key)
        {
            return valueMap_.at(key);
        }
        std::map<std::string, KeyValueTreeValue>::iterator
        addProperty(const std::string &key, KeyValueTreeValue &&value)
        {
            GMX_RELEASE_ASSERT(!keyExists(key), "Duplicate key value");
            values_.reserve(values_.size() + 1);
            auto iter = valueMap_.insert(std::make_pair(key, std::move(value))).first;
            values_.push_back(KeyValueTreeProperty(iter));
            return iter;
        }

        std::map<std::string, KeyValueTreeValue> valueMap_;
        std::vector<KeyValueTreeProperty>        values_;

        friend class KeyValueTreeObjectBuilder;
};

/********************************************************************
 * Inline functions that could not be declared within the classes
 */

inline bool KeyValueTreeValue::isArray() const
{
    return value_.isType<KeyValueTreeArray>();
}
inline bool KeyValueTreeValue::isObject() const
{
    return value_.isType<KeyValueTreeObject>();
}
inline const KeyValueTreeArray &KeyValueTreeValue::asArray() const
{
    return value_.cast<KeyValueTreeArray>();
}
inline const KeyValueTreeObject &KeyValueTreeValue::asObject() const
{
    return value_.cast<KeyValueTreeObject>();
}
inline KeyValueTreeArray &KeyValueTreeValue::asArray()
{
    return value_.castRef<KeyValueTreeArray>();
}
inline KeyValueTreeObject &KeyValueTreeValue::asObject()
{
    return value_.castRef<KeyValueTreeObject>();
}

//! \cond libapi
/*! \brief
 * Compares two KeyValueTrees and prints any differences.
 *
 * \ingroup module_utility
 */
void compareKeyValueTrees(TextWriter               *writer,
                          const KeyValueTreeObject &tree1,
                          const KeyValueTreeObject &tree2,
                          real                      ftol,
                          real                      abstol);
//! \endcond

} // namespace gmx

#endif
