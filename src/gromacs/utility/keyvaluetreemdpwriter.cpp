/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Defines a function to write a flat key-value tree to look like
 * old-style mdp output.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreemdpwriter.h"

#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

void writeKeyValueTreeAsMdp(TextWriter* writer, const KeyValueTreeObject& tree)
{
    for (const auto& prop : tree.properties())
    {
        const auto& value = prop.value();
        GMX_RELEASE_ASSERT(!value.isObject(), "Only flat key-value trees can be written as mdp");

        // Recognize a special key prefix that identifies comment
        // lines. This mechanism is not pretty, but our plan is to
        // write key-value trees, rather than old-style mdp, and
        // comments will need different handling then.
        if (prop.key().compare(0, 7, "comment") == 0)
        {
            GMX_RELEASE_ASSERT(prop.value().isType<std::string>(),
                               "Comments must have string-typed values");
            auto comment = prop.value().cast<std::string>();
            // TODO Consider implementing an MdpTextWriter that can
            // format an array of strings suitably, e.g. by prefixing
            // each line like "; %s". We'd need to implement arrays of
            // objects, such a comment key to be able to refer to an
            // array of strings. Also, we'd want to have a plan to
            // write such comments in whatever replaces the mdp
            // format.
            writer->writeLine(comment);
        }
        else
        {
            writer->writeString(formatString("%-24s = ", prop.key().c_str()));
            if (value.isArray())
            {
                bool first = true;
                for (const auto& elem : value.asArray().values())
                {
                    GMX_RELEASE_ASSERT(!elem.isObject() && !elem.isArray(),
                                       "Arrays of objects not currently implemented");
                    if (!first)
                    {
                        writer->writeString(" ");
                    }
                    writer->writeString(simpleValueToString(elem));
                    first = false;
                }
            }
            else
            {
                writer->writeString(simpleValueToString(value));
            }
            writer->writeLine();
        }
    }
}

} // namespace gmx
