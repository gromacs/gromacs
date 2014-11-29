/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <vector>
#include "gmxpre.h"
#include <stdio.h>
#include <string.h>
#include "gromacs/utility/fatalerror.h"
#include "plistwrapper.h"

namespace alexandria {

std::vector<PlistWrapper>::iterator SearchPlist(std::vector<PlistWrapper> &plist, int ftype)
{
    std::vector<PlistWrapper>::iterator pw;
    for(pw = plist.begin(); (pw < plist.end() && pw->getFtype() != ftype); ++pw)
    {
        ;
    }
    return pw;
}

unsigned int CountPlist(std::vector<PlistWrapper> &plist, int ftype)
{
    std::vector<PlistWrapper>::iterator pw = SearchPlist(plist, ftype);
    if (plist.end() != pw)
    {
        return pw->nParam();
    }
    return 0;
}

void delete_params(std::vector<PlistWrapper> &plist_,
                   const int etype,
                   const int alist[])
{
    int nra;
    std::vector<PlistWrapper>::iterator pwi = SearchPlist(plist_, etype);

    if (plist_.end() == pwi)
    {
        fprintf(stderr, "Cannot find function type %s in plist\n", 
                interaction_function[etype].name);
        return;
    }
    nra = interaction_function[etype].nratoms;
    switch (nra)
    {
        case 2:
            /* Remove bonds, if present */
            for (ParamIterator j = pwi->beginParam();
                 (j < pwi->endParam()); ++j)
            {
                if (((j->a[0] == alist[0]) &&
                     (j->a[1] == alist[1])) ||
                    ((j->a[1] == alist[0]) &&
                     (j->a[0] == alist[1])))
                {
                    if (NULL != debug)
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
            for (ParamIterator j = pwi->beginParam();
                 (j < pwi->endParam()); ++j)
            {
                if (j->a[1] == alist[1])
                {
                    if (((j->a[0] == alist[0]) &&
                         (j->a[2] == alist[2])) ||
                        ((j->a[2] == alist[0]) &&
                         (j->a[0] == alist[2])))
                    {
                        if (NULL != debug)
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
            for (ParamIterator j = pwi->beginParam();
                 (j < pwi->endParam()); ++j)
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
                    if (NULL != debug)
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
                    interaction_function[etype].name);
    }
}

void add_param_to_plist(std::vector<PlistWrapper> &plist,
                        const int ftype,
                        const t_param &p)
{
    std::vector<PlistWrapper>::iterator pwi = SearchPlist(plist, ftype);
    
    if (plist.end() == pwi)
    {
        PlistWrapper pw(ftype);
        plist.push_back(pw);
        pwi = plist.end() - 1;
    }
    pwi->addParam(p);
}

}
