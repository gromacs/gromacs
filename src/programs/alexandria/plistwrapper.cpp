/*! \internal \brief
 * Implements part of the alexandria program.
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

    if (plist_.end() == pwi)
    {
        if (NULL != debug)
        {
            fprintf(debug, "Cannot find function type %s in plist\n",
                    interaction_function[ftype].name);
        }
        return;
    }
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
            for (auto j = pwi->beginParam(); (j < pwi->endParam()); ++j)
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
                    interaction_function[ftype].name);
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
