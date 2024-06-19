/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements classes from stringstream.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/stringstream.h"

#include <cstddef>

#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

void StringOutputStream::write(const char* str)
{
    str_.append(str);
}

void StringOutputStream::close() {}

StringInputStream::StringInputStream(const std::string& input) : input_(input), pos_(0) {}

StringInputStream::StringInputStream(const std::vector<std::string>& input) :
    input_(joinStrings(input.begin(), input.end(), "\n")), pos_(0)
{
    input_.append("\n");
}

StringInputStream::StringInputStream(ArrayRef<const char* const> const& input) :
    input_(joinStrings(input.begin(), input.end(), "\n")), pos_(0)
{
    input_.append("\n");
}

bool StringInputStream::readLine(std::string* line)
{
    if (pos_ == input_.size())
    {
        line->clear();
        return false;
    }
    else
    {
        size_t newpos = input_.find('\n', pos_);
        if (newpos == std::string::npos)
        {
            newpos = input_.size();
        }
        else
        {
            // To include the newline as well!
            newpos += 1;
        }
        line->assign(input_.substr(pos_, newpos - pos_));
        pos_ = newpos;
        return true;
    }
}

} // namespace gmx
