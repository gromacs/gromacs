/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>
#include "gromacs/energyanalysis/select.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/cstringutil.h"

static void bool2set(const std::vector<bool> bE, std::vector<int> &set)
{
    for (unsigned int i = 0; (i < bE.size()); i++)
    {
        if (bE[i])
        {
            set.push_back(i);
        }
    }
}

void select_by_index(std::vector<std::string> nm,
                     std::vector<int>        &set)
{
    // TODO: getenv should be removed!
    bool bVerbose = ((getenv("VERBOSE")) != NULL);
    if (bVerbose)
    {
        int kk = 1;

        fprintf(stderr, "Select the terms you want from the following list\n");
        fprintf(stderr, "End your selection with 0\n");
        for (std::vector<std::string>::iterator k = nm.begin(); (k < nm.end()); )
        {
            for (int j = 0; (j < 4) && (k < nm.end()); j++, ++k)
            {
                fprintf(stderr, " %3d=%14s", kk++, k->c_str());
            }
            fprintf(stderr, "\n");
        }
    }

    // Initiate the bE array with falsehood
    std::vector<bool> bE;
    for (unsigned int kk = 0; (kk < nm.size()); kk++)
    {
        bE.push_back(false);
    }

    unsigned int n;
    do
    {
        if (1 != scanf("%d", &n))
        {
            // TODO: use c++ tools!
            gmx_fatal(FARGS, "Error reading user input");
        }
        if ((n > 0) && (n <= nm.size()))
        {
            bE[n-1] = true;
        }
    }
    while (n != 0);

    bool2set(bE, set);
}

static void chomp(char *buf)
{
    int len = strlen(buf);

    while ((len > 0) && (buf[len-1] == '\n'))
    {
        buf[len-1] = '\0';
        len--;
    }
}

void select_by_name(std::vector<std::string> nm, std::vector<int> &set)
{
    int                      j, i, nmatch, nss;
    bool                     bVerbose;
    const char              *fm4   = "%3d  %-14s";
    const char              *fm2   = "%3d  %-34s";
    std::vector<std::string> newnm;

    // TODO: remove getenv
    bVerbose = ((getenv("VERBOSE")) == NULL);

    if (bVerbose)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Select the terms you want from the following list by\n");
        fprintf(stderr, "selecting either (part of) the name or the number or a combination.\n");
        fprintf(stderr, "End your selection with an empty line or a zero.\n");
        fprintf(stderr, "-------------------------------------------------------------------\n");
    }

    j = 0;
    for (unsigned int k = 0; k < nm.size(); k++)
    {
        /* Insert dashes in all the names */
        std::string            buf = nm[k];
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

        bool bLong = false;
        if (bVerbose)
        {
            if (j == 0)
            {
                if (k > 0)
                {
                    fprintf(stderr, "\n");
                }
                bLong = false;
                for (unsigned int kk = k; kk < k+4; kk++)
                {
                    if (kk < nm.size() && (nm[kk].size() > 14))
                    {
                        bLong = true;
                    }
                }
            }
            else
            {
                fprintf(stderr, " ");
            }
            if (!bLong)
            {
                fprintf(stderr, fm4, k+1, newnm[k].c_str());
                j++;
                if (j == 4)
                {
                    j = 0;
                }
            }
            else
            {
                fprintf(stderr, fm2, k+1, newnm[k].c_str());
                j++;
                if (j == 2)
                {
                    j = 0;
                }
            }
        }
    }
    if (bVerbose)
    {
        fprintf(stderr, "\n\n");
    }

    // Initiate the bE array with falsehood
    std::vector<bool> bE;
    for (unsigned int kk = 0; (kk < nm.size()); kk++)
    {
        bE.push_back(false);
    }

    bool bEOF = false;
    char buf[320];
    while (!bEOF && (fgets2(buf, 320-1, stdin)))
    {
        /* Remove newlines */
        chomp(buf);

        /* Remove leading and trailing spaces */
        trim(buf);

        /* Empty line means end of input */
        bEOF = (strlen(buf) == 0);
        if (!bEOF)
        {
            char *ptr = buf;
            do
            {
                if (!bEOF)
                {
                    /* First try to read an integer */
                    unsigned int nind;
                    nss   = sscanf(ptr, "%10u", &nind);
                    if (nss == 1)
                    {
                        /* Zero means end of input */
                        if (nind == 0)
                        {
                            bEOF = true;
                        }
                        else if ((1 <= nind) && (nind <= nm.size()))
                        {
                            bE[nind-1] = true;
                        }
                        else
                        {
                            fprintf(stderr, "number %d is out of range\n", nind);
                        }
                    }
                    else
                    {
                        /* Now try to read a string */
                        i      = strlen(ptr);
                        nmatch = 0;
                        for (unsigned int nind = 0; nind < nm.size(); nind++)
                        {
                            if (gmx_strcasecmp(newnm[nind].c_str(), ptr) == 0)
                            {
                                bE[nind] = true;
                                nmatch++;
                            }
                        }
                        if (nmatch == 0)
                        {
                            i      = strlen(ptr);
                            nmatch = 0;
                            for (unsigned int nind = 0; nind < nm.size(); nind++)
                            {
                                if (gmx_strncasecmp(newnm[nind].c_str(),
                                                    ptr, i) == 0)
                                {
                                    bE[nind] = true;
                                    nmatch++;
                                }
                            }
                            if (nmatch == 0)
                            {
                                fprintf(stderr, "String '%s' does not match anything\n", ptr);
                            }
                        }
                    }
                }
                /* Look for the first space, and remove spaces from there */
                if ((ptr = strchr(ptr, ' ')) != NULL)
                {
                    trim(ptr);
                }
            }
            while (!bEOF && (ptr && (strlen(ptr) > 0)));
        }
    }

    bool2set(bE, set);
}
