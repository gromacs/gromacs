/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief
 * Implements classes in argsbuilder.h
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "argsbuilder.h"

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/string2.h"

namespace gmx
{

namespace test
{

/********************************************************************
 * ArgsBuilder
 */

ArgsBuilder::ArgsBuilder(const char* _programName) : arguments()
{
    std::string programName(_programName);
    std::string::size_type first_space = programName.find(' ');
    if (std::string::npos != first_space)
    {
        // Cater for two-part tool names like "gmx energy"
        arguments.push_back(programName.substr(0, first_space));
        arguments.push_back(programName.substr(first_space + 1));

        // TODO Cater for more than a single space between the two
        // parts of the name.
    }
    else
    {
        arguments.push_back(programName);
    }
}
    
ArgsBuilder::~ArgsBuilder()
{
    std::vector<char *>::iterator thisArg = theArgvHolder.begin();
    for(; theArgvHolder.end() != thisArg; ++thisArg)
    {
        sfree(*thisArg);
    }
}

void
ArgsBuilder::visitSubSection(const Options &section)
{
    gmx::OptionsIterator iterator(section);
    iterator.acceptSubSections(this);
    iterator.acceptOptions(this);
}

void
ArgsBuilder::visitOption(const OptionInfo &option)
{
    /* Options might be defined that have not been set. Only those
     * that have been set should be included in the arguments
     * list. */
    if(option.isSet())
    {
        arguments.push_back("-" + option.name());
        arguments.push_back(option.formatValue(0));
    }
}

int
ArgsBuilder::getArgc()
{
    return arguments.size();
}
    
char **
ArgsBuilder::getArgv()
{
    std::vector<std::string>::iterator thisArg = arguments.begin();
    for(int i = 0; arguments.end() != thisArg; ++thisArg, ++i)
    {
        // TODO find out how to avoid using gmx_strdup() and sfree().
        theArgvHolder.push_back(gmx_strdup(thisArg->c_str()));
        theArgv.push_back(theArgvHolder.back());
    }
    return &(theArgv[0]);
}
 
} // namespace test
} // namespace gmx
