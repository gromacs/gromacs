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

#include "keyvaluetree.h"

#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

namespace
{

//! Helper function to split a KeyValueTreePath to its components
std::vector<std::string> splitPathElements(const std::string &path)
{
    GMX_ASSERT(!path.empty() && path[0] == '/',
               "Paths to KeyValueTree should start with '/'");
    return splitDelimitedString(path.substr(1), '/');
}

//! Helper function to format a simple KeyValueTreeValue.
std::string formatSingleValue(const KeyValueTreeValue &value)
{
    if (value.isType<float>())
    {
        return formatString("%g", value.cast<float>());
    }
    else if (value.isType<double>())
    {
        return formatString("%g", value.cast<double>());
    }
    GMX_RELEASE_ASSERT(false, "Unknown value type");
    return std::string();
}

}   // namespace

/********************************************************************
 * KeyValueTreePath
 */

KeyValueTreePath::KeyValueTreePath(const std::string &path)
    : path_(splitPathElements(path))
{
}

std::string KeyValueTreePath::toString() const
{
    return "/" + joinStrings(path_, "/");
}

/********************************************************************
 * KeyValueTreeObject
 */

void KeyValueTreeObject::writeUsing(TextWriter *writer) const
{
    for (const auto &prop : properties())
    {
        const auto &value = prop.value();
        if (value.isObject())
        {
            writer->writeString(prop.key());
            writer->writeLine(":");
            int oldIndent = writer->wrapperSettings().indent();
            writer->wrapperSettings().setIndent(oldIndent + 2);
            value.asObject().writeUsing(writer);
            writer->wrapperSettings().setIndent(oldIndent);
        }
        else
        {
            int indent = writer->wrapperSettings().indent();
            writer->writeString(formatString("%*s", -(33-indent), prop.key().c_str()));
            writer->writeString(" = ");
            if (value.isArray())
            {
                writer->writeString("[");
                for (const auto &elem : value.asArray().values())
                {
                    GMX_RELEASE_ASSERT(!elem.isObject() && !elem.isArray(),
                                       "Arrays of objects not currently implemented");
                    writer->writeString(" ");
                    writer->writeString(formatSingleValue(elem));
                }
                writer->writeString(" ]");
            }
            else
            {
                writer->writeString(formatSingleValue(value));
            }
            writer->writeLine();
        }
    }
}

} // namespace gmx
