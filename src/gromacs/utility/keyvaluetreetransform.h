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
 * Declares utilities for transforming key-value trees.
 *
 * See \ref page_mdmodules for the main use case that these support.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_KEYVALUETREETRANSFORM_H
#define GMX_UTILITY_KEYVALUETREETRANSFORM_H

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/variant.h"

namespace gmx
{

class IKeyValueTreeErrorHandler;
class KeyValueTreeObjectBuilder;

enum class StringCompareType;

class KeyValueTreeTransformResult;
class KeyValueTreeTransformRuleBuilder;
class KeyValueTreeTransformRulesScoped;

namespace internal
{
class KeyValueTreeTransformerImpl;
}

/*! \libinternal \brief
 * Interface to declare rules for transforming key-value trees.
 *
 * This interface is used to add transformation rules for key-value trees.
 * A transformation is a set of rules that is used to map an input key-value
 * tree to an output key-value tree, with possible conversion steps performed
 * in the process.  Currently, each rule maps one item from the source tree to
 * one item in the target tree (it is possible to expand a single value into an
 * object with multiple properties).  See KeyValueTreeTransformRuleBuilder for
 * the kinds of rules currently supported.
 *
 * The main use currently is in converting flat-format mdp files to a
 * structured internal representation.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class IKeyValueTreeTransformRules
{
    public:
        /*! \brief
         * Creates a new rule.
         *
         * Properties of the new rule must be specified using the returned
         * builder.
         */
        virtual KeyValueTreeTransformRuleBuilder addRule() = 0;
        /*! \brief
         * Creates a scoped set of rules, where all rules use a target sub-tree.
         *
         * \param[in] scope Prefix defining the scope in the target tree
         *
         * Any rules added to the returned scope will have `scope` prefixed to
         * their target paths, i.e., it is not possible to produce elements
         * outside the specified subtree.
         */
        virtual KeyValueTreeTransformRulesScoped
        scopedTransform(const KeyValueTreePath &scope) = 0;

    protected:
        ~IKeyValueTreeTransformRules();
};

/*! \libinternal \brief
 * Helper object returned from IKeyValueTreeTransformRules::scopedTransform().
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeTransformRulesScoped
{
    public:
        //! Internal constructor for creating the scope.
        KeyValueTreeTransformRulesScoped(
            internal::KeyValueTreeTransformerImpl *impl,
            const KeyValueTreePath                &prefix);
        //! Supports returning the object from IKeyValueTreeTransformRules::scopedTransform().
        KeyValueTreeTransformRulesScoped(KeyValueTreeTransformRulesScoped &&other);
        //! Supports returning the object from IKeyValueTreeTransformRules::scopedTransform().
        KeyValueTreeTransformRulesScoped &operator=(KeyValueTreeTransformRulesScoped &&other);
        ~KeyValueTreeTransformRulesScoped();

        //! Returns the interface for adding rules to this scope.
        IKeyValueTreeTransformRules *rules();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief
 * Provides methods to specify one transformation rule.
 *
 * \if internal
 * The builder is implemented as a set of nested objects, each of which is
 * provides methods for setting a particular property of the rule.  Setting a
 * property returns another object that has relevant methods for the context.
 * This provides some structure to the methods, and catches at least some types
 * of incorrect rules already at compile time.
 * Additionally, if you use an IDE with completion facilities, it can nicely
 * guide you through which values you need to specify.
 * All values are stored within the main builder object, and the rule is
 * created at the end of the statement.
 * \endif
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class KeyValueTreeTransformRuleBuilder
{
    public:
        /*! \internal \brief
         * Base class used for implementing parameter provider objects.
         */
        class Base
        {
            protected:
                //! Creates a parameter provider object within given builder.
                explicit Base(KeyValueTreeTransformRuleBuilder *builder)
                    : builder_(builder)
                {
                }

                //! The parent builder.
                KeyValueTreeTransformRuleBuilder *builder_;
        };

        /*! \libinternal \brief
         * Properties that can be specified after from().to().
         *
         * \tparam FromType Type specified for from() to map from.
         * \tparam ToType Type specified for to() to map to.
         */
        template <typename FromType, typename ToType>
        class ToValue : public Base
        {
            public:
                //! Creates a parameter provider object within given builder.
                explicit ToValue(KeyValueTreeTransformRuleBuilder *builder)
                    : Base(builder)
                {
                }

                /*! \brief
                 * Specifies the transformation function to convert the value
                 * from FromType to ToType.
                 */
                void transformWith(std::function<ToType(const FromType &)> transform)
                {
                    builder_->addTransformToVariant(
                            [transform] (const Variant &value)
                            {
                                return Variant::create<ToType>(transform(value.cast<FromType>()));
                            });
                }
        };

        /*! \libinternal \brief
         * Properties that can be specified after from().toObject().
         *
         * \tparam FromType Type specified for from() to map from.
         */
        template <typename FromType>
        class ToObject : public Base
        {
            public:
                //! Creates a parameter provider object within given builder.
                explicit ToObject(KeyValueTreeTransformRuleBuilder *builder)
                    : Base(builder)
                {
                }

                /*! \brief
                 * Specifies the transformation function to build the output
                 * object.
                 *
                 * The transform should build the output object with the
                 * provided builder.
                 */
                void transformWith(std::function<void(KeyValueTreeObjectBuilder *, const FromType &)> transform)
                {
                    builder_->addTransformToObject(
                            [transform] (KeyValueTreeObjectBuilder *builder, const Variant &value)
                            {
                                transform(builder, value.cast<FromType>());
                            });
                }
        };

        /*! \libinternal \brief
         * Properties that can be specified after from().
         *
         * \tparam FromType Type specified for from() to map from.
         */
        template <typename FromType>
        class AfterFrom : public Base
        {
            public:
                //! Creates a parameter provider object within given builder.
                explicit AfterFrom(KeyValueTreeTransformRuleBuilder *builder)
                    : Base(builder)
                {
                }

                /*! \brief
                 * Specifies a rule that maps to a value at given path.
                 *
                 * \tparam ToType  Type to map to.
                 * \param[in] path Path to map to.
                 *
                 * It is an error if multiple rules map to the same path, or to
                 * a parent path of the target of an existing rule.
                 * Note that it is possible to have a to() rule map to a child
                 * of a toObject() rule, provided that the path is not created
                 * by the object rule.
                 */
                template <typename ToType>
                ToValue<FromType, ToType> to(const KeyValueTreePath &path)
                {
                    builder_->setToPath(path);
                    return ToValue<FromType, ToType>(builder_);
                }

                /*! \brief
                 * Specifies a rule that maps to an object (collection of named
                 * values) at given path.
                 *
                 * \param[in] path Path to map to.
                 *
                 * It is an error if multiple rules map to the same path, or to
                 * a parent path of the target of an existing rule.
                 * However, it is allowed to have two toObject() rules map to
                 * the same path, provided that the properties they produce are
                 * distinct.
                 */
                ToObject<FromType> toObject(const KeyValueTreePath &path)
                {
                    builder_->setToPath(path);
                    return ToObject<FromType>(builder_);
                }
        };

        //! Internal constructor for creating a builder.
        KeyValueTreeTransformRuleBuilder(internal::KeyValueTreeTransformerImpl *impl,
                                         const KeyValueTreePath                &prefix);
        //! Supports returning the builder from IKeyValueTreeTransformRules::addRule().
        KeyValueTreeTransformRuleBuilder(KeyValueTreeTransformRuleBuilder &&)            = default;
        //! Supports returning the builder from IKeyValueTreeTransformRules::addRule().
        KeyValueTreeTransformRuleBuilder &operator=(KeyValueTreeTransformRuleBuilder &&) = default;
        ~KeyValueTreeTransformRuleBuilder();

        /*! \brief
         * Specifies a rule that maps a value at given path.
         *
         * \tparam FromType Type of value expected at `path`.
         * \param[in] path Path to map in this rule.
         *
         * If the input tree has `path`, but it is not of type `FromType`,
         * the transform will produce an error.
         *
         * It is an error to use the same path in two from() rules.  Similarly,
         * it is an error to use a child path of a path used in a different
         * from() rule.
         */
        template <typename FromType>
        AfterFrom<FromType> from(const KeyValueTreePath &path)
        {
            setFromPath(path);
            setExpectedType(typeid(FromType));
            return AfterFrom<FromType>(this);
        }
        /*! \brief
         * Specifies how strings are matched when matching rules against a path.
         *
         * For properties of the object at `path`, `keyMatchType` is used for
         * string comparison.
         *
         * This rule must be specified first for a path, before any other
         * from() rule specifies the path or a subpath.
         * The rule only applies to immediate properties at the given path, not
         * recursively.
         * It is an error to specify the match type multiple times for a path.
         */
        void keyMatchType(const KeyValueTreePath &path, StringCompareType keyMatchType)
        {
            setFromPath(path);
            setKeyMatchType(keyMatchType);
        }

    private:
        void setFromPath(const KeyValueTreePath &path);
        void setExpectedType(const std::type_index &type);
        void setToPath(const KeyValueTreePath &path);
        void setKeyMatchType(StringCompareType keyMatchType);
        void addTransformToVariant(std::function<Variant(const Variant &)> transform);
        void addTransformToObject(std::function<void(KeyValueTreeObjectBuilder *, const Variant &)> transform);

        class Data;

        internal::KeyValueTreeTransformerImpl *impl_;
        std::unique_ptr<Data>                  data_;
};

class KeyValueTreeTransformer
{
    public:
        KeyValueTreeTransformer();
        ~KeyValueTreeTransformer();

        IKeyValueTreeTransformRules *rules();

        std::vector<KeyValueTreePath> mappedPaths() const;

        KeyValueTreeTransformResult
        transform(const KeyValueTreeObject  &tree,
                  IKeyValueTreeErrorHandler *errorHandler) const;

    private:
        PrivateImplPointer<internal::KeyValueTreeTransformerImpl> impl_;
};

class IKeyValueTreeBackMapping
{
    public:
        virtual ~IKeyValueTreeBackMapping();

        virtual KeyValueTreePath
        originalPath(const KeyValueTreePath &path) const = 0;
};

class KeyValueTreeTransformResult
{
    public:
        KeyValueTreeObject object() { return std::move(object_); }
        const IKeyValueTreeBackMapping &backMapping() const { return *mapping_; }

    private:
        typedef std::unique_ptr<IKeyValueTreeBackMapping> MappingPointer;

        KeyValueTreeTransformResult(KeyValueTreeObject &&object,
                                    MappingPointer     &&mapping)
            : object_(std::move(object)), mapping_(std::move(mapping))
        {
        }

        KeyValueTreeObject  object_;
        MappingPointer      mapping_;

        friend class internal::KeyValueTreeTransformerImpl;
};

} // namespace gmx

#endif
