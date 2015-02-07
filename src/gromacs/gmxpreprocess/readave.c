/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec.h"
#include "gromacs/utility/smalloc.h"
#include "readir.h"
#include "names.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/trnio.h"
#include "txtdump.h"

#include "gmx_fatal.h"

extern char **read_aveparams(int *ninp_p, t_inpfile **inp_p, t_ave *ave,
                             warninp_t wi)
{
    int         ninp, g, m;
    t_inpfile  *inp;
    const char *tmp;
    char      **grpbuf;
    char        buf[STRLEN];
    char        warn_buf[STRLEN];
    dvec        vec;
    t_avegrp   *aveg;

    ninp   = *ninp_p;
    inp    = *inp_p;

    /* read averaging parameters */
    CTYPE("Number of averaging groups");
    ITYPE("ave-ngroups",     ave->nave, 1);

    if (ave->nave < 1)
    {
        gmx_fatal(FARGS, "ave-ngroups should be >= 1");
    }

    snew(ave->grp, ave->nave);

    /* Read the averaging groups */
    snew(grpbuf, ave->nave);
    for (g = 0; g < ave->nave; g++)
    {
        aveg = &ave->grp[g];
        snew(grpbuf[g], STRLEN);
        CTYPE("Averaging group name");
        sprintf(buf, "ave-group%d", g);
        STYPE(buf, grpbuf[g], "");
    }

    *ninp_p   = ninp;
    *inp_p    = inp;

    return grpbuf;
}

extern void make_averaging_groups(t_ave *ave, char **avegnames, t_blocka *grps, char **gnames)
{
    int       g, ig = -1, i;
    t_avegrp *aveg;


    for (g = 0; g < ave->nave; g++)
    {
        aveg      = &ave->grp[g];
        ig        = search_string(avegnames[g], grps->nr, gnames);
        aveg->nat = grps->index[ig+1] - grps->index[ig];

        if (aveg->nat > 0)
        {
            fprintf(stderr, "Averaging group %d '%s' has %d atoms\n", g, avegnames[g], aveg->nat);
            snew(aveg->ind, aveg->nat);
            for (i = 0; i < aveg->nat; i++)
            {
                aveg->ind[i] = grps->a[grps->index[ig]+i];
            }
        }
        else
        {
            gmx_fatal(FARGS, "Averaging group %d '%s' is empty", g, avegnames[g]);
        }
    }
}
