/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreetransform.h"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>

#include "gromacs/utility/any.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/ikeyvaluetreeerror.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace internal
{
class KeyValueTreeTransformerImpl;
} // namespace internal

/********************************************************************
 * IKeyValueTreeTransformRules
 */

IKeyValueTreeTransformRules::~IKeyValueTreeTransformRules() {}

/********************************************************************
 * KeyValueTreeTransformRulesScoped::Impl
 */

class KeyValueTreeTransformRulesScoped::Impl : public IKeyValueTreeTransformRules
{
public:
    Impl(internal::KeyValueTreeTransformerImpl* impl, const KeyValueTreePath& prefix) :
        impl_(impl), prefix_(prefix)
    {
    }

    KeyValueTreeTransformRuleBuilder addRule() override
    {
        return KeyValueTreeTransformRuleBuilder(impl_, prefix_);
    }

    KeyValueTreeTransformRulesScoped scopedTransform(const KeyValueTreePath& scope) override
    {
        return KeyValueTreeTransformRulesScoped(impl_, prefix_ + scope);
    }

private:
    internal::KeyValueTreeTransformerImpl* impl_;
    KeyValueTreePath                       prefix_;
};

/********************************************************************
 * KeyValueTreeTransformRulesScoped
 */

KeyValueTreeTransformRulesScoped::KeyValueTreeTransformRulesScoped(internal::KeyValueTreeTransformerImpl* impl,
                                                                   const KeyValueTreePath& prefix) :
    impl_(new Impl(impl, prefix))
{
}

KeyValueTreeTransformRulesScoped::KeyValueTreeTransformRulesScoped(KeyValueTreeTransformRulesScoped&&) noexcept = default;

KeyValueTreeTransformRulesScoped&
KeyValueTreeTransformRulesScoped::operator=(KeyValueTreeTransformRulesScoped&&) noexcept = default;

KeyValueTreeTransformRulesScoped::~KeyValueTreeTransformRulesScoped() {}

IKeyValueTreeTransformRules* KeyValueTreeTransformRulesScoped::rules()
{
    return impl_.get();
}

/********************************************************************
 * IKeyValueTreeBackMapping
 */

IKeyValueTreeBackMapping::~IKeyValueTreeBackMapping() {}

namespace
{

class KeyValueTreeBackMapping : public IKeyValueTreeBackMapping
{
public:
    class Entry
    {
    public:
        Entry() = default;
        explicit Entry(const KeyValueTreePath& path) : sourcePath_(path) {}

        Entry* getOrCreateChildEntry(const std::string& key)
        {
            auto iter = childEntries_.find(key);
            if (iter == childEntries_.end())
            {
                iter = childEntries_.insert(std::make_pair(key, Entry())).first;
            }
            return &iter->second;
        }
        void setMapping(const KeyValueTreePath& path, const KeyValueTreeValue& value)
        {
            GMX_RELEASE_ASSERT(sourcePath_.empty(), "Multiple entries map to same path");
            if (value.isObject())
            {
                const KeyValueTreeObject& object = value.asObject();
                for (const auto& prop : object.properties())
                {
                    GMX_RELEASE_ASSERT(!prop.value().isObject(), "Nested objects not implemented");
                    childEntries_[prop.key()] = Entry(path);
                }
            }
            else
            {
                sourcePath_ = path;
            }
        }

        KeyValueTreePath             sourcePath_;
        std::map<std::string, Entry> childEntries_;
    };

    KeyValueTreePath originalPath(const KeyValueTreePath& path) const override
    {
        const Entry* entry = &rootEntry_;
        for (const auto& element : path.elements())
        {
            auto iter = entry->childEntries_.find(element);
            if (iter == entry->childEntries_.end())
            {
                break;
            }
            entry = &iter->second;
        }
        GMX_RELEASE_ASSERT(entry->childEntries_.empty() && !entry->sourcePath_.empty(),
                           "Requested path not uniquely mapped");
        return entry->sourcePath_;
    }

    Entry* rootEntry() { return &rootEntry_; }

private:
    Entry rootEntry_;
};

} // namespace

namespace internal
{

/********************************************************************
 * KeyValueTreeTransformerImpl
 */

class KeyValueTreeTransformerImpl
{
public:
    class Rule
    {
    public:
        typedef std::function<void(KeyValueTreeValueBuilder*, const KeyValueTreeValue&)> TransformFunction;
        typedef std::map<std::string, Rule, StringCompare> ChildRuleMap;

        explicit Rule(StringCompareType keyMatchType) :
            expectedType_(typeid(void)), childRules_(keyMatchType)
        {
        }

        const Rule* findMatchingChildRule(const std::string& key) const
        {
            auto iter = childRules_.find(key);
            if (iter == childRules_.end())
            {
                return nullptr;
            }
            return &iter->second;
        }
        Rule* getOrCreateChildRule(const std::string& key)
        {
            auto iter = childRules_.find(key);
            if (iter == childRules_.end())
            {
                return createChildRule(key, StringCompareType::Exact);
            }
            return &iter->second;
        }
        Rule* createChildRule(const std::string& key, StringCompareType keyMatchType)
        {
            auto result = childRules_.insert(std::make_pair(key, Rule(keyMatchType)));
            GMX_RELEASE_ASSERT(result.second, "Cannot specify key match type after child rules");
            return &result.first->second;
        }

        void collectMappedPaths(const KeyValueTreePath& prefix, std::vector<KeyValueTreePath>* result) const
        {
            for (const auto& value : childRules_)
            {
                KeyValueTreePath path = prefix;
                path.append(value.first);
                const Rule& rule = value.second;
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

        KeyValueTreePath  targetPath_;
        std::string       targetKey_;
        std::type_index   expectedType_;
        TransformFunction transform_;
        ChildRuleMap      childRules_;
    };

    class Transformer
    {
    public:
        explicit Transformer(IKeyValueTreeErrorHandler* errorHandler) :
            errorHandler_(errorHandler), backMapping_(new KeyValueTreeBackMapping)
        {
            if (errorHandler_ == nullptr)
            {
                errorHandler_ = defaultKeyValueTreeErrorHandler();
            }
        }

        void transform(const Rule* rootRule, const KeyValueTreeObject& tree)
        {
            if (rootRule != nullptr)
            {
                doChildTransforms(rootRule, tree);
            }
        }

        KeyValueTreeTransformResult result()
        {
            return KeyValueTreeTransformResult(builder_.build(), std::move(backMapping_));
        }

    private:
        void doTransform(const Rule* rule, const KeyValueTreeValue& value);
        void doChildTransforms(const Rule* rule, const KeyValueTreeObject& object);
        void applyTransformedValue(const Rule* rule, KeyValueTreeValue&& value);

        IKeyValueTreeErrorHandler*               errorHandler_;
        KeyValueTreeBuilder                      builder_;
        std::unique_ptr<KeyValueTreeBackMapping> backMapping_;
        KeyValueTreePath                         context_;
    };

    KeyValueTreeTransformerImpl() : rootScope_(this, KeyValueTreePath()) {}

    Rule* getOrCreateRootRule()
    {
        if (rootRule_ == nullptr)
        {
            createRootRule(StringCompareType::Exact);
        }
        return rootRule_.get();
    }
    void createRootRule(StringCompareType keyMatchType)
    {
        GMX_RELEASE_ASSERT(rootRule_ == nullptr, "Cannot specify key match type after child rules");
        rootRule_ = std::make_unique<Rule>(keyMatchType);
    }

    std::unique_ptr<Rule>            rootRule_;
    KeyValueTreeTransformRulesScoped rootScope_;
};

/********************************************************************
 * KeyValueTreeTransformerImpl::Transformer
 */

void KeyValueTreeTransformerImpl::Transformer::doTransform(const Rule* rule, const KeyValueTreeValue& value)
{
    if (rule->transform_ != nullptr)
    {
        KeyValueTreeValueBuilder valueBuilder;
        try
        {
            if (value.type() != rule->expectedType_)
            {
                // TODO: Better error message.
                GMX_THROW(InvalidInputError("Unexpected type of value"));
            }
            rule->transform_(&valueBuilder, value);
        }
        catch (UserInputError& ex)
        {
            if (!errorHandler_->onError(&ex, context_))
            {
                throw;
            }
            return;
        }
        applyTransformedValue(rule, valueBuilder.build());
        return;
    }
    if (!rule->childRules_.empty())
    {
        doChildTransforms(rule, value.asObject());
    }
}

void KeyValueTreeTransformerImpl::Transformer::doChildTransforms(const Rule*               rule,
                                                                 const KeyValueTreeObject& object)
{
    for (const auto& prop : object.properties())
    {
        const Rule* childRule = rule->findMatchingChildRule(prop.key());
        if (childRule != nullptr)
        {
            context_.append(prop.key());
            doTransform(childRule, prop.value());
            context_.pop_back();
        }
    }
}

void KeyValueTreeTransformerImpl::Transformer::applyTransformedValue(const Rule*         rule,
                                                                     KeyValueTreeValue&& value)
{
    KeyValueTreeObjectBuilder       objBuilder = builder_.rootObject();
    KeyValueTreeBackMapping::Entry* mapEntry   = backMapping_->rootEntry();
    for (const std::string& key : rule->targetPath_.elements())
    {
        if (objBuilder.keyExists(key))
        {
            GMX_RELEASE_ASSERT(objBuilder[key].isObject(),
                               "Inconsistent transform (different items map to same path)");
            objBuilder = objBuilder.getObjectBuilder(key);
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
        GMX_RELEASE_ASSERT(value.isObject(),
                           "Inconsistent transform (different items map to same path)");
        GMX_RELEASE_ASSERT(objBuilder[rule->targetKey_].isObject(),
                           "Inconsistent transform (different items map to same path)");
        objBuilder = objBuilder.getObjectBuilder(rule->targetKey_);
        GMX_RELEASE_ASSERT(objBuilder.objectHasDistinctProperties(value.asObject()),
                           "Inconsistent transform (different items map to same path)");
        objBuilder.mergeObject(std::move(value.asObject()));
    }
    else
    {
        objBuilder.addRawValue(rule->targetKey_, std::move(value));
    }
}

} // namespace internal

/********************************************************************
 * KeyValueTreeTransformer
 */

KeyValueTreeTransformer::KeyValueTreeTransformer() :
    impl_(new internal::KeyValueTreeTransformerImpl)
{
}

KeyValueTreeTransformer::~KeyValueTreeTransformer() {}

IKeyValueTreeTransformRules* KeyValueTreeTransformer::rules()
{
    return impl_->rootScope_.rules();
}

std::vector<KeyValueTreePath> KeyValueTreeTransformer::mappedPaths() const
{
    std::vector<KeyValueTreePath> result;
    if (impl_->rootRule_)
    {
        impl_->rootRule_->collectMappedPaths(KeyValueTreePath(), &result);
    }
    return result;
}

KeyValueTreeTransformResult KeyValueTreeTransformer::transform(const KeyValueTreeObject& tree,
                                                               IKeyValueTreeErrorHandler* errorHandler) const
{
    internal::KeyValueTreeTransformerImpl::Transformer transformer(errorHandler);
    transformer.transform(impl_->rootRule_.get(), tree);
    return transformer.result();
}

/********************************************************************
 * KeyValueTreeTransformRuleBuilder::Data
 */

class KeyValueTreeTransformRuleBuilder::Data
{
public:
    typedef internal::KeyValueTreeTransformerImpl::Rule Rule;

    explicit Data(const KeyValueTreePath& prefix) :
        prefixPath_(prefix),
        expectedType_(typeid(void)),
        keyMatchType_(StringCompareType::Exact),
        keyMatchRule_(false)
    {
    }

    void createRule(internal::KeyValueTreeTransformerImpl* impl)
    {
        if (keyMatchRule_)
        {
            createRuleWithKeyMatchType(impl);
            return;
        }
        GMX_RELEASE_ASSERT(transform_ != nullptr, "Transform has not been specified");
        Rule* rule = impl->getOrCreateRootRule();
        for (const std::string& key : fromPath_.elements())
        {
            GMX_RELEASE_ASSERT(rule->targetKey_.empty(),
                               "Cannot specify multiple rules from a single path");
            rule = rule->getOrCreateChildRule(key);
        }
        GMX_RELEASE_ASSERT(rule->targetKey_.empty(),
                           "Cannot specify multiple rules from a single path");
        rule->targetKey_    = toPath_.pop_last();
        rule->targetPath_   = std::move(toPath_);
        rule->expectedType_ = expectedType_;
        rule->transform_    = transform_;
    }

    void createRuleWithKeyMatchType(internal::KeyValueTreeTransformerImpl* impl)
    {
        if (fromPath_.empty())
        {
            impl->createRootRule(keyMatchType_);
        }
        else
        {
            std::string lastKey = fromPath_.pop_last();
            Rule*       rule    = impl->getOrCreateRootRule();
            for (const std::string& key : fromPath_.elements())
            {
                rule = rule->getOrCreateChildRule(key);
            }
            rule->createChildRule(lastKey, keyMatchType_);
        }
    }

    const KeyValueTreePath  prefixPath_;
    KeyValueTreePath        fromPath_;
    KeyValueTreePath        toPath_;
    std::type_index         expectedType_;
    Rule::TransformFunction transform_;
    StringCompareType       keyMatchType_;
    bool                    keyMatchRule_;
};

/********************************************************************
 * KeyValueTreeTransformRuleBuilder
 */

KeyValueTreeTransformRuleBuilder::KeyValueTreeTransformRuleBuilder(internal::KeyValueTreeTransformerImpl* impl,
                                                                   const KeyValueTreePath& prefix) :
    impl_(impl), data_(new Data(prefix))
{
}

// TODO If the function called here would ever throw
// (e.g. std::bad_alloc) then std::terminate will be called.
// Alternatively, createRule could catch all exceptions and terminate
// but that's the same for an end user. So suppressing the clang-tidy
// warning is about as bad as any alternative.
KeyValueTreeTransformRuleBuilder::~KeyValueTreeTransformRuleBuilder() // NOLINT(bugprone-exception-escape)
{
    if (!std::uncaught_exceptions())
    {
        data_->createRule(impl_);
    }
}

void KeyValueTreeTransformRuleBuilder::setFromPath(const KeyValueTreePath& path)
{
    data_->fromPath_ = path;
}

void KeyValueTreeTransformRuleBuilder::setExpectedType(const std::type_index& type)
{
    data_->expectedType_ = type;
}

void KeyValueTreeTransformRuleBuilder::setToPath(const KeyValueTreePath& path)
{
    data_->toPath_ = data_->prefixPath_ + path;
}

void KeyValueTreeTransformRuleBuilder::setKeyMatchType(StringCompareType keyMatchType)
{
    data_->keyMatchType_ = keyMatchType;
    data_->keyMatchRule_ = true;
}

void KeyValueTreeTransformRuleBuilder::addTransformToAny(const std::function<Any(const Any&)>& transform)
{
    data_->transform_ = [transform](KeyValueTreeValueBuilder* builder, const KeyValueTreeValue& value) {
        builder->setAnyValue(transform(value.asAny()));
    };
}

void KeyValueTreeTransformRuleBuilder::addTransformToObject(
        const std::function<void(KeyValueTreeObjectBuilder*, const Any&)>& transform)
{
    data_->transform_ = [transform](KeyValueTreeValueBuilder* builder, const KeyValueTreeValue& value) {
        KeyValueTreeObjectBuilder obj = builder->createObject();
        transform(&obj, value.asAny());
    };
}

} // namespace gmx
