/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "coolstuff.h"

#include "config.h"

#include <cstdlib>
#include <cstring>

/* This file is completely threadsafe - keep it that way! */

#include "gromacs/fileio/strdb.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool be_cool(void)
{
    /* Yes, it is bad to check the environment variable every call,
     * but we dont call this routine often, and it avoids using
     * a mutex for locking the variable...
     */
#if GMX_COOL_QUOTES
    return (getenv("GMX_NO_QUOTES") == NULL);
#else
    /*be uncool*/
    return FALSE;
#endif
}

static void pukeit(const char *db, const char *defstring, char *retstring,
                   int retsize, int *cqnum)
{
    FILE     *fp;
    char    **help;
    int       i, nhlp;
    gmx_rng_t rng;

    if (be_cool() && ((fp = low_libopen(db, FALSE)) != NULL))
    {
        nhlp = fget_lines(fp, &help);
        /* for libraries we can use the low-level close routines */
        gmx_ffclose(fp);
        rng    = gmx_rng_init(gmx_rng_make_seed());
        *cqnum = static_cast<int>(nhlp*gmx_rng_uniform_real(rng));
        gmx_rng_destroy(rng);
        if (std::strlen(help[*cqnum]) >= STRLEN)
        {
            help[*cqnum][STRLEN-1] = '\0';
        }
        std::strncpy(retstring, help[*cqnum], retsize);
        for (i = 0; (i < nhlp); i++)
        {
            sfree(help[i]);
        }
        sfree(help);
    }
    else
    {
        *cqnum = -1;
        std::strncpy(retstring, defstring, retsize);
    }
}

void bromacs(char *retstring, int retsize)
{
    int dum;

    pukeit("bromacs.dat",
           "Groningen Machine for Chemical Simulation",
           retstring, retsize, &dum);
}

void cool_quote(char *retstring, int retsize, int *cqnum)
{
    char *tmpstr;
    char *ptr;
    int   tmpcq, *p;

    if (cqnum != NULL)
    {
        p = cqnum;
    }
    else
    {
        p = &tmpcq;
    }

    /* protect audience from explicit lyrics */
    snew(tmpstr, retsize+1);
    pukeit("gurgle.dat", "Thanx for Using GROMACS - Have a Nice Day",
           tmpstr, retsize-2, p);

    if ((ptr = std::strchr(tmpstr, '_')) != NULL)
    {
        *ptr = '\0';
        ptr++;
        sprintf(retstring, "\"%s\" %s", tmpstr, ptr);
    }
    else
    {
        std::strcpy(retstring, tmpstr);
    }
    sfree(tmpstr);
}
