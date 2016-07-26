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
/*! \libinternal \file
 * \brief
 * Declares a data structure for JSON-like structured key-value mapping.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_KEYVALUETREETRANSFORM_H
#define GMX_UTILITY_KEYVALUETREETRANSFORM_H

#include <functional>
#include <string>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/variant.h"

namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;

class KeyValueTreeTransformRuleBuilder;

namespace internal
{
class KeyValueTreeTransformerImpl;
}

class IKeyValueTreeTransformRules
{
    public:
        virtual KeyValueTreeTransformRuleBuilder addRule() = 0;

    protected:
        ~IKeyValueTreeTransformRules();
};

class KeyValueTreeTransformer
{
    public:
        KeyValueTreeTransformer();
        ~KeyValueTreeTransformer();

        IKeyValueTreeTransformRules *rules();

        KeyValueTreeObject transform(const KeyValueTreeObject &tree);

    private:
        PrivateImplPointer<internal::KeyValueTreeTransformerImpl> impl_;
};

class KeyValueTreeTransformRuleBuilder
{
    public:
        class Base
        {
            protected:
                explicit Base(KeyValueTreeTransformRuleBuilder *builder)
                    : builder_(builder)
                {
                }

                KeyValueTreeTransformRuleBuilder *builder_;
        };

        template <typename T>
        class ToObject : public Base
        {
            public:
                explicit ToObject(KeyValueTreeTransformRuleBuilder *builder)
                    : Base(builder)
                {
                }

                void transformWith(std::function<void(KeyValueTreeObjectBuilder *, const T &)> transform)
                {
                    builder_->addTransform(
                            [transform] (KeyValueTreeObjectBuilder *builder, const Variant &value)
                            {
                                transform(builder, value.cast<T>());
                            });
                }
        };

        template <typename T>
        class AfterFrom : public Base
        {
            public:
                explicit AfterFrom(KeyValueTreeTransformRuleBuilder *builder)
                    : Base(builder)
                {
                }

                ToObject<T> toObject(const std::string &path)
                {
                    builder_->setToPath(path);
                    return ToObject<T>(builder_);
                }
        };

        explicit KeyValueTreeTransformRuleBuilder(internal::KeyValueTreeTransformerImpl *impl);
        KeyValueTreeTransformRuleBuilder(KeyValueTreeTransformRuleBuilder &&)            = default;
        KeyValueTreeTransformRuleBuilder &operator=(KeyValueTreeTransformRuleBuilder &&) = default;
        ~KeyValueTreeTransformRuleBuilder();

        template <typename T>
        AfterFrom<T> from(const std::string &path)
        {
            setFromPath(path);
            return AfterFrom<T>(this);
        }

    private:
        void setFromPath(const std::string &path);
        void setToPath(const std::string &path);
        void addTransform(std::function<void(KeyValueTreeObjectBuilder *, const Variant &)> transform);

        class Data;

        internal::KeyValueTreeTransformerImpl *impl_;
        std::unique_ptr<Data>                  data_;
};

} // namespace gmx

#endif
