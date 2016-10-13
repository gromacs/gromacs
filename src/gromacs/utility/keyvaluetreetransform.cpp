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

#include "keyvaluetreetransform.h"

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/********************************************************************
 * IKeyValueTreeTransformRules
 */

IKeyValueTreeTransformRules::~IKeyValueTreeTransformRules()
{
}

/********************************************************************
 * IKeyValueTreeBackMapping
 */

IKeyValueTreeBackMapping::~IKeyValueTreeBackMapping()
{
}

namespace
{

class KeyValueTreeBackMapping : public IKeyValueTreeBackMapping
{
    public:
        class Entry
        {
            public:
                Entry() = default;
                explicit Entry(const std::vector<std::string> &path) : sourcePath_(path) {}

                Entry *getOrCreateChildEntry(const std::string &key)
                {
                    auto iter = childEntries_.find(key);
                    if (iter == childEntries_.end())
                    {
                        iter = childEntries_.insert(std::make_pair(key, Entry())).first;
                    }
                    return &iter->second;
                }
                void setMapping(const std::vector<std::string> &path,
                                const KeyValueTreeValue        &value)
                {
                    if (value.isObject())
                    {
                        const KeyValueTreeObject &object = value.asObject();
                        for (const auto &prop : object.properties())
                        {
                            childEntries_[prop.key()] = Entry(path);
                        }
                    }
                    else
                    {
                        sourcePath_ = path;
                    }
                }

                std::vector<std::string>     sourcePath_;
                std::map<std::string, Entry> childEntries_;
        };

        virtual std::vector<std::string>
        originalPath(const std::vector<std::string> &path) const
        {
            const Entry *entry = &rootEntry_;
            for (const auto &element : path)
            {
                auto iter = entry->childEntries_.find(element);
                if (iter == entry->childEntries_.end())
                {
                    break;
                }
                entry = &iter->second;
            }
            GMX_RELEASE_ASSERT(entry->childEntries_.empty()
                               && !entry->sourcePath_.empty(),
                               "Requested path not uniquely mapped");
            return entry->sourcePath_;
        }

        Entry *rootEntry() { return &rootEntry_; }

    private:
        Entry rootEntry_;
};

}   // namespace

namespace internal
{

/********************************************************************
 * KeyValueTreeTransformerImpl
 */

class KeyValueTreeTransformerImpl : public IKeyValueTreeTransformRules
{
    public:
        class Rule
        {
            public:
                typedef std::function<void(KeyValueTreeValueBuilder *, const KeyValueTreeValue &)>
                    TransformFunction;
                typedef std::map<std::string, Rule, StringCompare> ChildRuleMap;

                explicit Rule(StringCompareType keyMatchType)
                    : childRules_(keyMatchType)
                {
                }

                const Rule *findMatchingChildRule(const std::string &key) const
                {
                    auto iter = childRules_.find(key);
                    if (iter == childRules_.end())
                    {
                        return nullptr;
                    }
                    return &iter->second;
                }
                Rule *getOrCreateChildRule(const std::string &key)
                {
                    auto iter = childRules_.find(key);
                    if (iter == childRules_.end())
                    {
                        return createChildRule(key, StringCompareType::Exact);
                    }
                    return &iter->second;
                }
                Rule *createChildRule(const std::string &key,
                                      StringCompareType  keyMatchType)
                {
                    auto result = childRules_.insert(std::make_pair(key, Rule(keyMatchType)));
                    GMX_RELEASE_ASSERT(result.second,
                                       "Cannot specify key match type after child rules");
                    return &result.first->second;
                }

                void collectMappedPaths(const std::string        &prefix,
                                        std::vector<std::string> *result) const
                {
                    for (const auto &value : childRules_)
                    {
                        std::string path = prefix + "/" + value.first;
                        const Rule &rule = value.second;
                        if (rule.transform_)
                        {
                            result->push_back(path);
                        }
                        else
                        {
                            rule.collectMappedPaths(path, result);
                        }
                    }
                }

                std::vector<std::string>    targetPath_;
                std::string                 targetKey_;
                TransformFunction           transform_;
                ChildRuleMap                childRules_;
        };

        class Transformer
        {
            public:
                Transformer(const KeyValueTreeTransformerImpl &impl)
                    : impl_(impl), backMapping_(new KeyValueTreeBackMapping)
                {
                }

                void transform(const KeyValueTreeObject &tree)
                {
                    doChildTransforms(impl_.rootRule_.get(), tree);
                }

                KeyValueTreeTransformResult result()
                {
                    return KeyValueTreeTransformResult(builder_.build(),
                                                       std::move(backMapping_));
                }

            private:
                void doTransform(const Rule *rule, const KeyValueTreeValue &value);
                void doChildTransforms(const Rule *rule, const KeyValueTreeObject &object);
                void applyTransformedValue(const Rule *rule, KeyValueTreeValue &&value);

                const KeyValueTreeTransformerImpl       &impl_;
                KeyValueTreeBuilder                      builder_;
                std::unique_ptr<KeyValueTreeBackMapping> backMapping_;
                std::vector<std::string>                 context_;
        };

        virtual KeyValueTreeTransformRuleBuilder addRule()
        {
            return KeyValueTreeTransformRuleBuilder(this);
        }

        Rule *getOrCreateRootRule()
        {
            if (rootRule_ == nullptr)
            {
                createRootRule(StringCompareType::Exact);
            }
            return rootRule_.get();
        }
        void createRootRule(StringCompareType keyMatchType)
        {
            GMX_RELEASE_ASSERT(rootRule_ == nullptr,
                               "Cannot specify key match type after child rules");
            rootRule_.reset(new Rule(keyMatchType));
        }

        std::unique_ptr<Rule>  rootRule_;
};

/********************************************************************
 * KeyValueTreeTransformerImpl::Transformer
 */

void KeyValueTreeTransformerImpl::Transformer::doTransform(
        const Rule *rule, const KeyValueTreeValue &value)
{
    if (rule->transform_)
    {
        KeyValueTreeValueBuilder valueBuilder;
        rule->transform_(&valueBuilder, value);
        applyTransformedValue(rule, valueBuilder.build());
        return;
    }
    if (!rule->childRules_.empty())
    {
        doChildTransforms(rule, value.asObject());
    }
}

void KeyValueTreeTransformerImpl::Transformer::doChildTransforms(
        const Rule *rule, const KeyValueTreeObject &object)
{
    for (const auto &prop : object.properties())
    {
        const Rule *childRule = rule->findMatchingChildRule(prop.key());
        if (childRule != nullptr)
        {
            context_.push_back(prop.key());
            doTransform(childRule, prop.value());
            context_.pop_back();
        }
    }
}

void KeyValueTreeTransformerImpl::Transformer::applyTransformedValue(
        const Rule *rule, KeyValueTreeValue &&value)
{
    KeyValueTreeObjectBuilder       objBuilder = builder_.rootObject();
    KeyValueTreeBackMapping::Entry *mapEntry   = backMapping_->rootEntry();
    for (const std::string &key : rule->targetPath_)
    {
        if (objBuilder.keyExists(key))
        {
            objBuilder = objBuilder.getObject(key);
        }
        else
        {
            objBuilder = objBuilder.addObject(key);
        }
        mapEntry = mapEntry->getOrCreateChildEntry(key);
    }
    mapEntry = mapEntry->getOrCreateChildEntry(rule->targetKey_);
    mapEntry->setMapping(context_, value);
    if (objBuilder.keyExists(rule->targetKey_))
    {
        objBuilder.getObject(rule->targetKey_).mergeObject(std::move(value));
    }
    else
    {
        objBuilder.addRawValue(rule->targetKey_, std::move(value));
    }
}

}   // namespace internal

/********************************************************************
 * KeyValueTreeTransformer
 */

KeyValueTreeTransformer::KeyValueTreeTransformer()
    : impl_(new internal::KeyValueTreeTransformerImpl)
{
}

KeyValueTreeTransformer::~KeyValueTreeTransformer()
{
}

IKeyValueTreeTransformRules *KeyValueTreeTransformer::rules()
{
    return impl_.get();
}

std::vector<std::string> KeyValueTreeTransformer::mappedPaths() const
{
    std::vector<std::string> result;
    if (impl_->rootRule_)
    {
        impl_->rootRule_->collectMappedPaths("", &result);
    }
    return result;
}

KeyValueTreeTransformResult
KeyValueTreeTransformer::transform(const KeyValueTreeObject &tree) const
{
    internal::KeyValueTreeTransformerImpl::Transformer transformer(*impl_);
    transformer.transform(tree);
    return transformer.result();
}

/********************************************************************
 * KeyValueTreeTransformRuleBuilder::Data
 */

class KeyValueTreeTransformRuleBuilder::Data
{
    public:
        typedef internal::KeyValueTreeTransformerImpl::Rule Rule;

        Data() : keyMatchType_(StringCompareType::Exact) {}

        void createRule(internal::KeyValueTreeTransformerImpl *impl)
        {
            if (toPath_.empty())
            {
                createRuleWithKeyMatchType(impl);
                return;
            }
            std::vector<std::string>  from = splitDelimitedString(fromPath_.substr(1), '/');
            std::vector<std::string>  to   = splitDelimitedString(toPath_.substr(1), '/');
            Rule                     *rule = impl->getOrCreateRootRule();
            for (const std::string &key : from)
            {
                rule = rule->getOrCreateChildRule(key);
            }
            rule->targetKey_  = to.back();
            to.pop_back();
            rule->targetPath_ = std::move(to);
            rule->transform_  = transform_;
        }

        void createRuleWithKeyMatchType(internal::KeyValueTreeTransformerImpl *impl)
        {
            std::vector<std::string> from = splitDelimitedString(fromPath_.substr(1), '/');
            if (from.empty())
            {
                impl->createRootRule(keyMatchType_);
            }
            else
            {
                std::string lastKey = from.back();
                from.pop_back();
                Rule       *rule = impl->getOrCreateRootRule();
                for (const std::string &key : from)
                {
                    rule = rule->getOrCreateChildRule(key);
                }
                rule->createChildRule(lastKey, keyMatchType_);
            }
        }

        std::string              fromPath_;
        std::string              toPath_;
        Rule::TransformFunction  transform_;
        StringCompareType        keyMatchType_;
};

/********************************************************************
 * KeyValueTreeTransformRuleBuilder
 */

KeyValueTreeTransformRuleBuilder::KeyValueTreeTransformRuleBuilder(
        internal::KeyValueTreeTransformerImpl *impl)
    : impl_(impl), data_(new Data)
{
}

KeyValueTreeTransformRuleBuilder::~KeyValueTreeTransformRuleBuilder()
{
    if (!std::uncaught_exception())
    {
        data_->createRule(impl_);
    }
}

void KeyValueTreeTransformRuleBuilder::setFromPath(const std::string &path)
{
    data_->fromPath_ = path;
}

void KeyValueTreeTransformRuleBuilder::setToPath(const std::string &path)
{
    data_->toPath_ = path;
}

void KeyValueTreeTransformRuleBuilder::setKeyMatchType(StringCompareType keyMatchType)
{
    data_->keyMatchType_ = keyMatchType;
}

void KeyValueTreeTransformRuleBuilder::addTransformToVariant(
        std::function<Variant(const Variant &)> transform)
{
    data_->transform_ =
        [transform] (KeyValueTreeValueBuilder *builder, const KeyValueTreeValue &value)
        {
            builder->setVariantValue(transform(value.asVariant()));
        };
}

void KeyValueTreeTransformRuleBuilder::addTransformToObject(
        std::function<void(KeyValueTreeObjectBuilder *, const Variant &)> transform)
{
    data_->transform_ =
        [transform] (KeyValueTreeValueBuilder *builder, const KeyValueTreeValue &value)
        {
            KeyValueTreeObjectBuilder obj = builder->createObject();
            transform(&obj, value.asVariant());
        };
}

} // namespace gmx
