/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "gmxpre.h"

#include "plistwrapper.h"

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype)
{
    return std::find_if(plist.begin(), plist.end(),
                        [ftype](PlistWrapper &p) 
                        { return p.getFtype() == ftype; });
}

PlistWrapperIterator SearchPlist(std::vector<PlistWrapper> &plist, InteractionType itype)
{
    return std::find_if(plist.begin(), plist.end(),
                        [itype](PlistWrapper &p) 
                        { return p.interactionType() == itype; });
}

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype)
{
    return std::count_if(plist.begin(), plist.end(),
                         [ftype](PlistWrapper &p) 
                         { return p.getFtype() == ftype; });
}

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int                  ftype,
                   const int                  alist[])
{
    int nra;
    std::vector<PlistWrapper>::iterator pwi = SearchPlist(plist_, ftype);

    if (plist_.end() != pwi)
    {
        nra = interaction_function[ftype].nratoms;
        switch (nra)
        {
            case 2:
                /* Remove bonds, if present */
                for (auto j = pwi->beginParam(); (j < pwi->endParam()); ++j)
                {
                    if (((j->a[0] == alist[0]) &&
                         (j->a[1] == alist[1])) ||
                        ((j->a[1] == alist[0]) &&
                         (j->a[0] == alist[1])))
                    {
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Removing bond beteen atoms %d %d\n",
                                    alist[0], alist[1]);
                        }
                        j = pwi->eraseParam(j);
                        break;
                    }
                }
           break;
            case 3:
                /* Remove angle, if present */
                for (auto j = pwi->beginParam(); (j < pwi->endParam()); ++j)
                {
                    if (j->a[1] == alist[1])
                    {
                        if (((j->a[0] == alist[0]) &&
                             (j->a[2] == alist[2])) ||
                            ((j->a[2] == alist[0]) &&
                             (j->a[0] == alist[2])))
                        {
                            if (nullptr != debug)
                            {
                                fprintf(debug, "Removing angle beteen atoms %d %d %d\n",
                                        alist[0], alist[1], alist[2]);
                            }
                            j = pwi->eraseParam(j);
                            break;
                        }
                    }
                }
	       break;
            case 4:
                /* Remove dihedral, if present. Allow wildcard in alist[3] (specified as -1) */
                for (auto j = pwi->beginParam(); (j < pwi->endParam()); ++j)
                {
                    if (((j->a[0] == alist[0]) &&
                         (j->a[1] == alist[1]) &&
                         (j->a[2] == alist[2]) &&
                         ((alist[3] == -1) || (j->a[3] == alist[3]))) ||
                        ((j->a[3] == alist[0]) &&
                         (j->a[2] == alist[1]) &&
                         (j->a[1] == alist[2]) &&
                         ((alist[3] == -1) || (j->a[0] == alist[3]))) ||
                        ((j->a[1] == alist[0]) &&
                         (j->a[2] == alist[1]) &&
                         (j->a[3] == alist[2]) &&
                         (alist[3] == -1)))
                    {
                        if (nullptr != debug)
                        {
                            fprintf(debug, "Removing dihedral beteen atoms %d %d %d %d\n",
                                    alist[0], alist[1], alist[2], alist[3]);
                        }
                        j = pwi->eraseParam(j);
                        break;
                    }
                }
	       break;
            default:
                fprintf(stderr, "Don't know how to remove params from type %s\n",
                        interaction_function[ftype].name);
        }
    }
}
void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        int                        ftype,
                        InteractionType            itype,
                        const t_param             &p)
{
    std::vector<PlistWrapper>::iterator pwi = SearchPlist(plist, ftype);

    if (plist.end() == pwi)
    {
        PlistWrapper pw(itype, ftype);
        plist.push_back(pw);
        pwi = plist.end() - 1;
    }
    pwi->addParam(p);
}

}
