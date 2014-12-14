/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements function in select_energy.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "select_energy.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"

#include "analysismodule.h"

namespace gmx
{

namespace energyanalysis
{

void select_energies(ConstArrayRef<EnergyNameUnit> eNU,
                     bool                          bVerbose,
                     TextInputStream              *input,
                     std::vector<int>             &set)
{
    std::vector<std::string> newnm;

    for (auto enu : eNU)
    {
        /* Insert dashes in all the names */
        std::string            buf = enu.energyName;
        std::string::size_type index;
        do
        {
            index = buf.find(' ');
            if (index != buf.npos)
            {
                buf[index] = '-';
            }
        }
        while (index != buf.npos);
        newnm.push_back(buf);
    }

    //bVerbose = bVerbose && (terms.size() == 0);
    if (bVerbose)
    {
        bool         bLong = false;
        unsigned int j     = 0;
        unsigned int k     = 0;

        printf("\n");
        printf("Select the terms you want from the following list by\n");
        printf("selecting either (part of) the name or the number or a combination.\n");
        printf("End your selection with 0, an empty line or Ctrl-D.\n");
        printf("-------------------------------------------------------------------\n");

        for (auto enu : eNU)
        {
            if (j == 0)
            {
                bLong = false;
                for (unsigned int kk = k; kk < k+4; kk++)
                {
                    if (kk < eNU.size() && (enu.energyName.size() > 14))
                    {
                        bLong = true;
                    }
                }
            }
            else
            {
                printf(" ");
            }
            if (!bLong)
            {
                printf("%3u  %-14s", k+1, newnm[k].c_str());
                j++;
                if (j == 4)
                {
                    j = 0;
                }
            }
            else
            {
                printf("%3u  %-34s", k+1, newnm[k].c_str());
                j++;
                if (j == 2)
                {
                    j = 0;
                }
            }
        }
        printf("\n\n");
    }

    set.clear();
    bool         done = false;
    std::string  line;
    while (!done && input->readLine(&line))
    {
        std::vector<std::string> subs = splitString(line);
        if (subs.size() == 0)
        {
            done = true;
        }
        else
        {
            for (std::vector<std::string>::iterator sub = subs.begin(); !done && (sub < subs.end()); ++sub)
            {
                char *endptr;
                // First check whether the input is an integer
                errno = 0;
                unsigned long kk = std::strtoul(sub->c_str(), &endptr, 0);
                if ((errno == ERANGE) || (errno == EINVAL))
                {
                    // Not an integer, now check strings
                    kk = 0;
                    for (std::vector<std::string>::iterator nn = newnm.begin();
                         (nn < newnm.end()); ++nn, ++kk)
                    {
                        if (0 == gmx_strcasecmp(nn->c_str(), sub->c_str()))
                        {
                            break;
                        }
                    }
                }
                if (0 == kk)
                {
                    // Time to finish up
                    done = true;
                }
                else if (kk <= newnm.size())
                {
                    set.push_back(static_cast<int>(kk-1));
                }
                else
                {
                    GMX_THROW(InvalidInputError("Invalid energy selection " + *sub));
                }
            }
        }
    }
    input->close();
}

} // namespace energyanalysis

} // namespace gmx
