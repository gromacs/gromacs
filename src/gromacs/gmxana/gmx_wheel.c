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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/writeps.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool *bPhobics(int nres, char *resnm[])
{
    int       i, nb;
    char    **cb;
    gmx_bool *bb;

    nb = get_strings("phbres.dat", &cb);
    snew(bb, nres);

    for (i = 0; (i < nres); i++)
    {
        if (search_str(nb, cb, resnm[i]) != -1)
        {
            bb[i] = TRUE;
        }
    }
    return bb;
}

void wheel(const char *fn, int nres, char *resnm[], int r0, real rot0, char *title)
{
    const real fontsize  = 16;
    const real gray      = 0.9;
    const real fontasp   = 0.6;
    const real fontwidth = fontsize*fontasp;

    t_psdata   out;
    int        i, sl, slen;
    real       ring, inner, outer;
    real       xc, yc, box;
    gmx_bool  *bPh;
    char     **rnms;
    char       sign;

    inner = 75.0;
    slen  = 0;
    snew(rnms, nres);
    for (i = 0; (i < nres); i++)
    {
        snew(rnms[i], 256);
        sl   = strlen(resnm[i]);
        sign = resnm[i][sl-1];
        if ((sign == '+') || (sign == '-'))
        {
            resnm[i][sl-1] = '\0';
        }
        sprintf(rnms[i], "%s-%d", resnm[i], i+r0);
        if ((sign == '+') || (sign == '-'))
        {
            sl            = strlen(rnms[i]);
            rnms[i][sl]   = sign;
            rnms[i][sl+1] = '\0';
        }

        slen = max(slen, (int)strlen(rnms[i]));
    }
    ring  = (2+slen)*fontwidth;
    outer = inner+ring;
    box   = inner*1.5+(1+(nres / 18))*ring;

    bPh = bPhobics(nres, resnm);

    out = ps_open(fn, 0, 0, 2.0*box, 2.0*box);
    xc  = box;
    yc  = box;

    ps_font(out, efontHELV, 1.5*fontsize);
    ps_translate(out, xc, yc);
    if (title)
    {
        ps_ctext(out, 0, -fontsize*1.5/2.0, title, eXCenter);
    }
    ps_font(out, efontHELV, fontsize);
    ps_rotate(out, rot0);
    for (i = 0; (i < nres); )
    {
        if (bPh[i])
        {
            ps_color(out, gray, gray, gray);
            ps_fillarcslice(out, 0, 0, inner, outer, -10, 10);
            ps_color(out, 0, 0, 0);
        }
        ps_arcslice(out, 0, 0, inner, outer, -10, 10);

        ps_ctext(out, inner+fontwidth, -fontsize/2.0, rnms[i], eXLeft);
        ps_rotate(out, -100);
        i++;

        if ((i % 18) == 0)
        {
            inner  = outer;
            outer += ring;
        }
    }
    ps_close(out);
}

void wheel2(const char *fn, int nres, char *resnm[], real rot0, char *title)
{
    const real fontsize  = 14;
    const real gray      = 0.9;
    const real fontasp   = 0.45;
    const int  angle     = 9;
    const real fontwidth = fontsize*fontasp;

    t_psdata   out;
    int        i, slen;
    real       ring, inner, outer;
    real       xc, yc, box;

    inner = 60.0;
    slen  = 0;
    for (i = 0; (i < nres); i++)
    {
        slen = max(slen, (int)strlen(resnm[i]));
    }
    fprintf(stderr, "slen = %d\n", slen);
    ring  = (slen)*fontwidth;
    outer = inner+ring;
    box   = (1+(nres / (2*angle)))*outer;

    out = ps_open(fn, 0, 0, 2.0*box, 2.0*box);
    xc  = box;
    yc  = box;

    ps_font(out, efontHELV, 1.5*fontsize);
    ps_translate(out, xc, yc);
    ps_color(out, 0, 0, 0);
    if (title)
    {
        ps_ctext(out, 0, -fontsize*1.5/2.0, title, eXCenter);
    }
    ps_font(out, efontHELV, fontsize);

    ps_rotate(out, rot0);
    for (i = 0; (i < nres); )
    {
        if ((i % 5) == 4)
        {
            ps_color(out, gray, gray, 1.0);
            ps_fillarcslice(out, 0, 0, inner, outer, -angle, angle);
            ps_color(out, 0, 0, 0);
        }
        ps_arcslice(out, 0, 0, inner, outer, -angle, angle);

        ps_ctext(out, inner+fontwidth, -fontsize/2.0, resnm[i], eXLeft);
        ps_rotate(out, -2*angle);
        i++;

        if ((i % (2*angle)) == 0)
        {
            inner  = outer;
            outer += ring;
        }
    }
    ps_close(out);
}

int gmx_wheel(int argc, char *argv[])
{
    const char     *desc[] = {
        "[THISMODULE] plots a helical wheel representation of your sequence.",
        "The input sequence is in the [REF].dat[ref] file where the first line contains",
        "the number of residues and each consecutive line contains a residue "
        "name."
    };
    output_env_t    oenv;
    static real     rot0  = 0;
    static gmx_bool bNum  = TRUE;
    static char    *title = NULL;
    static int      r0    = 1;
    t_pargs         pa [] = {
        { "-r0",  FALSE, etINT, {&r0},
          "The first residue number in the sequence" },
        { "-rot0", FALSE, etREAL, {&rot0},
          "Rotate around an angle initially (90 degrees makes sense)" },
        { "-T",   FALSE, etSTR, {&title},
          "Plot a title in the center of the wheel (must be shorter than 10 characters, or it will overwrite the wheel)" },
        { "-nn",  FALSE, etBOOL, {&bNum},
          "Toggle numbers" }
    };
    t_filenm        fnm[] = {
        { efDAT, "-f", NULL,  ffREAD  },
        { efEPS, "-o", NULL,  ffWRITE }
    };
#define NFILE asize(fnm)

    int    i, nres;
    char **resnm;

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    for (i = 1; (i < argc); i++)
    {
        if (strcmp(argv[i], "-r0") == 0)
        {
            r0 = strtol(argv[++i], NULL, 10);
            fprintf(stderr, "First residue is %d\n", r0);
        }
        else if (strcmp(argv[i], "-rot0") == 0)
        {
            rot0 = strtod(argv[++i], NULL);
            fprintf(stderr, "Initial rotation is %g\n", rot0);
        }
        else if (strcmp(argv[i], "-T") == 0)
        {
            title = gmx_strdup(argv[++i]);
            fprintf(stderr, "Title will be '%s'\n", title);
        }
        else if (strcmp(argv[i], "-nn") == 0)
        {
            bNum = FALSE;
            fprintf(stderr, "No residue numbers\n");
        }
        else
        {
            gmx_fatal(FARGS, "Incorrect usage of option %s", argv[i]);
        }
    }

    nres = get_lines(ftp2fn(efDAT, NFILE, fnm), &resnm);
    if (bNum)
    {
        wheel(ftp2fn(efEPS, NFILE, fnm), nres, resnm, r0, rot0, title);
    }
    else
    {
        wheel2(ftp2fn(efEPS, NFILE, fnm), nres, resnm, rot0, title);
    }

    return 0;
}
