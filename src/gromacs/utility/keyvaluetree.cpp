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

#include "gromacs/utility/keyvaluetree.h"

#include <cstdint>

#include <algorithm>
#include <map>
#include <string>
#include <typeindex>
#include <vector>

#include "gromacs/utility/compare.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

namespace
{

//! Helper function to split a KeyValueTreePath to its components
std::vector<std::string> splitPathElements(const std::string& path)
{
    GMX_ASSERT(!path.empty() && path[0] == '/', "Paths to KeyValueTree should start with '/'");
    return splitDelimitedString(path.substr(1), '/');
}

} // namespace

/********************************************************************
 * KeyValueTreePath
 */

KeyValueTreePath::KeyValueTreePath(const char* path) : path_(splitPathElements(path)) {}

KeyValueTreePath::KeyValueTreePath(const std::string& path) : path_(splitPathElements(path)) {}

std::string KeyValueTreePath::toString() const
{
    return std::string("/") + joinStrings(path_, "/");
}

/********************************************************************
 * KeyValueTreeObject
 */

bool KeyValueTreeObject::hasDistinctProperties(const KeyValueTreeObject& obj) const
{
    for (const auto& prop : obj.values_)
    {
        if (keyExists(prop.key()))
        {
            GMX_RELEASE_ASSERT(!prop.value().isArray(), "Comparison of arrays not implemented");
            if (prop.value().isObject() && valueMap_.at(prop.key()).isObject())
            {
                return valueMap_.at(prop.key()).asObject().hasDistinctProperties(prop.value().asObject());
            }
            return false;
        }
    }
    return true;
}

/********************************************************************
 * Key value tree dump
 */

//! \cond libapi
void dumpKeyValueTree(TextWriter* writer, const KeyValueTreeObject& tree)
{
    for (const auto& prop : tree.properties())
    {
        const auto& value = prop.value();
        if (value.isObject())
        {
            writer->writeString(prop.key());
            writer->writeLine(":");
            int oldIndent = writer->wrapperSettings().indent();
            writer->wrapperSettings().setIndent(oldIndent + 2);
            dumpKeyValueTree(writer, value.asObject());
            writer->wrapperSettings().setIndent(oldIndent);
        }
        else if (value.isArray()
                 && std::all_of(value.asArray().values().begin(),
                                value.asArray().values().end(),
                                [](const auto& elem) { return elem.isObject(); }))
        {
            // Array containing only objects
            writer->writeString(prop.key());
            writer->writeLine(":");
            int oldIndent = writer->wrapperSettings().indent();
            writer->wrapperSettings().setIndent(oldIndent + 2);
            for (const auto& elem : value.asArray().values())
            {
                dumpKeyValueTree(writer, elem.asObject());
            }
            writer->wrapperSettings().setIndent(oldIndent);
        }
        else
        {
            int indent = writer->wrapperSettings().indent();
            writer->writeString(formatString("%*s", -(33 - indent), prop.key().c_str()));
            writer->writeString(" = ");
            if (value.isArray())
            {
                writer->writeString("[");
                for (const auto& elem : value.asArray().values())
                {
                    GMX_RELEASE_ASSERT(
                            !elem.isObject() && !elem.isArray(),
                            "Only arrays of simple types and array of objects are implemented. "
                            "Arrays of arrays and mixed arrays are not supported.");
                    writer->writeString(" ");
                    writer->writeString(simpleValueToString(elem));
                }
                writer->writeString(" ]");
            }
            else
            {
                writer->writeString(simpleValueToString(value));
            }
            writer->writeLine();
        }
    }
}
//! \endcond

/********************************************************************
 * Key value tree comparison
 */

namespace
{

class CompareHelper
{
public:
    CompareHelper(TextWriter* writer, real ftol, real abstol) :
        writer_(writer), ftol_(ftol), abstol_(abstol)
    {
    }
    void compareObjects(const KeyValueTreeObject& obj1, const KeyValueTreeObject& obj2)
    {
        for (const auto& prop1 : obj1.properties())
        {
            currentPath_.append(prop1.key());
            if (obj2.keyExists(prop1.key()))
            {
                compareValues(prop1.value(), obj2[prop1.key()]);
            }
            else
            {
                handleMissingKeyInSecondObject(prop1.value());
            }
            currentPath_.pop_back();
        }
        for (const auto& prop2 : obj2.properties())
        {
            currentPath_.append(prop2.key());
            if (!obj1.keyExists(prop2.key()))
            {
                handleMissingKeyInFirstObject(prop2.value());
            }
            currentPath_.pop_back();
        }
    }

private:
    void compareValues(const KeyValueTreeValue& value1, const KeyValueTreeValue& value2)
    {
        if (value1.type() == value2.type())
        {
            if (value1.isObject())
            {
                compareObjects(value1.asObject(), value2.asObject());
            }
            else if (value1.isArray())
            {
                GMX_RELEASE_ASSERT(false, "Array comparison not implemented");
            }
            else if (!areSimpleValuesOfSameTypeEqual(value1, value2))
            {
                writer_->writeString(currentPath_.toString());
                writer_->writeLine(formatString(" (%s - %s)",
                                                simpleValueToString(value1).c_str(),
                                                simpleValueToString(value2).c_str()));
            }
        }
        else if ((value1.isType<double>() && value2.isType<float>())
                 || (value1.isType<float>() && value2.isType<double>()))
        {
            const bool  firstIsDouble = value1.isType<double>();
            const float v1 = firstIsDouble ? value1.cast<double>() : value1.cast<float>();
            const float v2 = firstIsDouble ? value2.cast<float>() : value2.cast<double>();
            if (!equal_float(v1, v2, ftol_, abstol_))
            {
                writer_->writeString(currentPath_.toString());
                writer_->writeLine(formatString(" (%e - %e)", v1, v2));
            }
        }
        else
        {
            handleMismatchingTypes(value1, value2);
        }
    }

    bool areSimpleValuesOfSameTypeEqual(const KeyValueTreeValue& value1, const KeyValueTreeValue& value2) const
    {
        GMX_ASSERT(value1.type() == value2.type(), "Caller should ensure that types are equal");
        if (value1.isType<bool>())
        {
            return value1.cast<bool>() == value2.cast<bool>();
        }
        else if (value1.isType<int>())
        {
            return value1.cast<int>() == value2.cast<int>();
        }
        else if (value1.isType<int64_t>())
        {
            return value1.cast<int64_t>() == value2.cast<int64_t>();
        }
        else if (value1.isType<double>())
        {
            return equal_double(value1.cast<double>(), value2.cast<double>(), ftol_, abstol_);
        }
        else if (value1.isType<float>())
        {
            return equal_float(value1.cast<float>(), value2.cast<float>(), ftol_, abstol_);
        }
        else if (value1.isType<std::string>())
        {
            return value1.cast<std::string>() == value2.cast<std::string>();
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Unknown value type");
            return false;
        }
    }

    void handleMismatchingTypes(const KeyValueTreeValue& /* value1 */,
                                const KeyValueTreeValue& /* value2 */)
    {
        writer_->writeString(currentPath_.toString());
        writer_->writeString(" type mismatch");
    }

    void handleMissingKeyInFirstObject(const KeyValueTreeValue& value)
    {
        const std::string message = formatString("%s (missing - %s)",
                                                 currentPath_.toString().c_str(),
                                                 formatValueForMissingMessage(value).c_str());
        writer_->writeLine(message);
    }
    void handleMissingKeyInSecondObject(const KeyValueTreeValue& value)
    {
        const std::string message = formatString("%s (%s - missing)",
                                                 currentPath_.toString().c_str(),
                                                 formatValueForMissingMessage(value).c_str());
        writer_->writeLine(message);
    }

    static std::string formatValueForMissingMessage(const KeyValueTreeValue& value)
    {
        if (value.isObject() || value.isArray())
        {
            return "present";
        }
        return simpleValueToString(value);
    }

    KeyValueTreePath currentPath_;
    TextWriter*      writer_;
    real             ftol_;
    real             abstol_;
};

} // namespace

//! \cond libapi
void compareKeyValueTrees(TextWriter*               writer,
                          const KeyValueTreeObject& tree1,
                          const KeyValueTreeObject& tree2,
                          real                      ftol,
                          real                      abstol)
{
    CompareHelper helper(writer, ftol, abstol);
    helper.compareObjects(tree1, tree2);
}
//! \endcond

} // namespace gmx
