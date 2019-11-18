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

#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int count_hydrogens(char*** atomname, int nra, gmx::ArrayRef<const int> a)
{
    int nh;

    if (!atomname)
    {
        gmx_fatal(FARGS, "Cannot call count_hydrogens with no atomname (%s %d)", __FILE__, __LINE__);
    }

    nh = 0;
    for (int i = 0; (i < nra); i++)
    {
        if (toupper(**(atomname[a[i]])) == 'H')
        {
            nh++;
        }
    }

    return nh;
}

void make_shake(gmx::ArrayRef<InteractionsOfType> plist, t_atoms* atoms, int nshake)
{
    char*** info = atoms->atomname;
    real    b_ij, b_jk;
    if (nshake != eshNONE)
    {
        switch (nshake)
        {
            case eshHBONDS: printf("turning H bonds into constraints...\n"); break;
            case eshALLBONDS: printf("turning all bonds into constraints...\n"); break;
            case eshHANGLES: printf("turning all bonds and H angles into constraints...\n"); break;
            case eshALLANGLES: printf("turning all bonds and angles into constraints...\n"); break;
            default: gmx_fatal(FARGS, "Invalid option for make_shake (%d)", nshake);
        }

        if ((nshake == eshHANGLES) || (nshake == eshALLANGLES))
        {
            /* Add all the angles with hydrogens to the shake list
             * and remove them from the bond list
             */
            for (int ftype = 0; (ftype < F_NRE); ftype++)
            {
                if (interaction_function[ftype].flags & IF_CHEMBOND)
                {
                    InteractionsOfType* bonds = &(plist[ftype]);

                    for (int ftype_a = 0; (gmx::ssize(*bonds) > 0 && ftype_a < F_NRE); ftype_a++)
                    {
                        if (interaction_function[ftype_a].flags & IF_ATYPE)
                        {
                            InteractionsOfType* pr = &(plist[ftype_a]);

                            for (auto parm = pr->interactionTypes.begin();
                                 parm != pr->interactionTypes.end();)
                            {
                                const InteractionOfType* ang = &(*parm);
#ifdef DEBUG
                                printf("Angle: %d-%d-%d\n", ang->ai(), ang->aj(), ang->ak());
#endif
                                int numhydrogens = count_hydrogens(info, 3, ang->atoms());
                                if ((nshake == eshALLANGLES) || (numhydrogens > 1)
                                    || (numhydrogens == 1 && toupper(**(info[ang->aj()])) == 'O'))
                                {
                                    /* Can only add hydrogen angle shake, if the two bonds
                                     * are constrained.
                                     * append this angle to the shake list
                                     */
                                    std::vector<int> atomNumbers = { ang->ai(), ang->ak() };

                                    /* Calculate length of constraint */
                                    bool bFound = false;
                                    b_ij = b_jk = 0.0;
                                    for (const auto& bond : bonds->interactionTypes)
                                    {
                                        if (((bond.ai() == ang->ai()) && (bond.aj() == ang->aj()))
                                            || ((bond.ai() == ang->aj()) && (bond.aj() == ang->ai())))
                                        {
                                            b_ij = bond.c0();
                                        }
                                        if (((bond.ai() == ang->ak()) && (bond.aj() == ang->aj()))
                                            || ((bond.ai() == ang->aj()) && (bond.aj() == ang->ak())))
                                        {
                                            b_jk = bond.c0();
                                        }
                                        bFound = (b_ij != 0.0) && (b_jk != 0.0);
                                    }
                                    if (bFound)
                                    {
                                        real param = std::sqrt(b_ij * b_ij + b_jk * b_jk
                                                               - 2.0 * b_ij * b_jk
                                                                         * cos(DEG2RAD * ang->c0()));
                                        std::vector<real> forceParm = { param, param };
                                        if (ftype == F_CONNBONDS || ftype_a == F_CONNBONDS)
                                        {
                                            gmx_fatal(FARGS,
                                                      "Can not constrain all angles when they "
                                                      "involved bonds of type %s",
                                                      interaction_function[F_CONNBONDS].longname);
                                        }
                                        /* apply law of cosines */
#ifdef DEBUG
                                        printf("p: %d, q: %d, dist: %12.5e\n", atomNumbers[0],
                                               atomNumbers[1], forceParm[0]);
#endif
                                        add_param_to_list(&(plist[F_CONSTR]),
                                                          InteractionOfType(atomNumbers, forceParm));
                                        /* move the last bond to this position */
                                        *parm = *(pr->interactionTypes.end() - 1);
                                        pr->interactionTypes.erase(pr->interactionTypes.end() - 1);
                                    }
                                }
                                else
                                {
                                    ++parm;
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
                InteractionsOfType* pr = &(plist[ftype]);
                for (auto parm = pr->interactionTypes.begin(); parm != pr->interactionTypes.end();)
                {
                    if ((nshake != eshHBONDS) || (count_hydrogens(info, 2, parm->atoms()) > 0))
                    {
                        /* append this bond to the shake list */
                        std::vector<int>  atomNumbers = { parm->ai(), parm->aj() };
                        std::vector<real> forceParm   = { parm->c0(), parm->c2() };
                        add_param_to_list(&(plist[F_CONSTR]), InteractionOfType(atomNumbers, forceParm));
                        parm = pr->interactionTypes.erase(parm);
                    }
                    else
                    {
                        ++parm;
                    }
                }
            }
        }
    }
}
