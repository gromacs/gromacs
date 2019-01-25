/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "topshake.h"

#include <cctype>
#include <cmath>

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int count_hydrogens (char ***atomname, gmx::ArrayRef<const int> a)
{
    if (!atomname)
    {
        gmx_fatal(FARGS, "Cannot call count_hydrogens with no atomname (%s %d)",
                  __FILE__, __LINE__);
    }

    int nh = 0;
    for (const auto &i : a)
    {
        if (toupper(**(atomname[i])) == 'H')
        {
            nh++;
        }
    }

    return nh;
}

void make_shake (gmx::ArrayRef<t_params> plist, t_atoms *atoms, int nshake)
{
    char          ***info = atoms->atomname;
    t_param          p;
    real             b_ij, b_jk;

    if (nshake != eshNONE)
    {
        switch (nshake)
        {
            case eshHBONDS:
                printf("turning H bonds into constraints...\n");
                break;
            case eshALLBONDS:
                printf("turning all bonds into constraints...\n");
                break;
            case eshHANGLES:
                printf("turning all bonds and H angles into constraints...\n");
                break;
            case eshALLANGLES:
                printf("turning all bonds and angles into constraints...\n");
                break;
            default:
                gmx_fatal(FARGS, "Invalid option for make_shake (%d)", nshake);
        }

        if ((nshake == eshHANGLES) || (nshake == eshALLANGLES))
        {
            /* Add all the angles with hydrogens to the shake list
             * and remove them from the bond list
             */
            for (int ftype = 0; (ftype < F_NRE); ftype++)
            {
                if (interaction_function[ftype].flags & IF_BTYPE)
                {
                    t_params *bonds = &(plist[ftype]);

                    for (int ftype_a = 0; (bonds->nr() > 0 && ftype_a < F_NRE); ftype_a++)
                    {
                        if (interaction_function[ftype_a].flags & IF_ATYPE)
                        {
                            t_params *pr = &(plist[ftype_a]);

                            for (int i = 0; (i < pr->nr()); )
                            {
                                t_param *ang = &(pr->param[i]);
#ifdef DEBUG
                                printf("Angle: %d-%d-%d\n", ang->ai(), ang->aj(), ang->ak());
#endif
                                int numhydrogens = count_hydrogens(info, ang->a);
                                if ((nshake == eshALLANGLES) ||
                                    (numhydrogens > 1) ||
                                    (numhydrogens == 1 && toupper(**(info[ang->a[1]])) == 'O'))
                                {
                                    /* Can only add hydrogen angle shake, if the two bonds
                                     * are constrained.
                                     * append this angle to the shake list
                                     */
                                    p.ai() = ang->ai();
                                    p.aj() = ang->ak();

                                    /* Calculate length of constraint */
                                    bool bFound = false;
                                    b_ij   = b_jk = 0.0;
                                    for (int j = 0; !bFound && (j < bonds->nr()); j++)
                                    {
                                        t_param *bond = &(bonds->param[j]);
                                        if (((bond->ai() == ang->ai()) &&
                                             (bond->aj() == ang->aj())) ||
                                            ((bond->ai() == ang->aj()) &&
                                             (bond->aj() == ang->ai())))
                                        {
                                            b_ij = bond->c0();
                                        }
                                        if (((bond->ai() == ang->ak()) &&
                                             (bond->aj() == ang->aj())) ||
                                            ((bond->ai() == ang->aj()) &&
                                             (bond->aj() == ang->ak())))
                                        {
                                            b_jk = bond->c0();
                                        }
                                        bFound = (b_ij != 0.0) && (b_jk != 0.0);
                                    }
                                    if (bFound)
                                    {
                                        /* apply law of cosines */
                                        p.c0() = std::sqrt( b_ij*b_ij + b_jk*b_jk -
                                                            2.0*b_ij*b_jk*cos(DEG2RAD*ang->c0()) );
                                        p.c1() = p.c0();
#ifdef DEBUG
                                        printf("p: %d, q: %d, dist: %12.5e\n", p.ai(), p.aj(), p.c0());
#endif
                                        add_param_to_list (&(plist[F_CONSTR]), &p);
                                        /* move the last bond to this position */
                                        pr->param[i] = pr->param.back();
                                        pr->param.erase(pr->param.end() - 1);
                                    }
                                }
                                else
                                {
                                    i++;
                                }
                            }
                        } /* if IF_ATYPE */
                    }     /* for ftype_A */
                }         /* if IF_BTYPE */
            }             /* for ftype */
        }                 /* if shake angles */

        /* Add all the bonds with hydrogens to the shake list
         * and remove them from the bond list
         */
        for (int ftype = 0; (ftype < F_NRE); ftype++)
        {
            if (interaction_function[ftype].flags & IF_BTYPE)
            {
                t_params            *pr = &(plist[ftype]);
                std::vector<t_param> cleanList;
                for (const auto &entry : pr->param)
                {
                    if ( (nshake != eshHBONDS) ||
                         (count_hydrogens (info, entry.a) > 0) )
                    {
                        /* append this bond to the shake list */
                        p.ai() = entry.ai();
                        p.aj() = entry.aj();
                        p.c0() = entry.c0();
                        p.c1() = entry.c2();
                        add_param_to_list (&(plist[F_CONSTR]), &p);
                        cleanList.push_back(entry);
                    }
                }
                pr->param = cleanList;
            }
        }
    }
}
