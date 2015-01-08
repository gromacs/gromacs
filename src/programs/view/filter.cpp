/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2013, The GROMACS development team.
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

#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "dialogs.h"
#include "xdlghi.h"

t_filter *init_filter(t_atoms *atoms, const char *fn, int natom_trx)
{
    t_filter *f;
    int       g, i;

    snew(f, 1);
    if (fn != NULL)
    {
        f->grps = init_index(fn, &f->grpnames);
    }
    else
    {
        snew(f->grps, 1);
        snew(f->grps->index, 1);
        analyse(atoms, f->grps, &f->grpnames, false, false);
    }
    snew(f->bDisable, f->grps->nr);
    for (g = 0; g < f->grps->nr; g++)
    {
        for (i = f->grps->index[g]; i < f->grps->index[g+1] && !f->bDisable[g]; i++)
        {
            f->bDisable[g] = (f->grps->a[i] >= natom_trx);
        }
    }

    snew(f->bShow, f->grps->nr);

    return f;
}

static void FilterCB(t_x11 *x11, int dlg_mess, int /*item_id*/,
                     char *set, void *data)
{
    int       nset;
    t_filter *f;
    t_gmx    *gmx;
    t_dlg    *dlg;

    gmx = (t_gmx *)data;
    dlg = gmx->dlgs[edFilter];
    f   = gmx->filter;

#ifdef DEBUG
    printf("item_id: %d, set: %s\n", item_id, set);
#endif
    switch (dlg_mess)
    {
        case DLG_SET:
            if (set)
            {
                if (sscanf(set, "%10d", &nset) == 1)
                {
                    f->bShow[nset] = !f->bShow[nset];
                }
            }
            break;
        case DLG_EXIT:
            HideDlg(dlg);
            write_gmx(x11, gmx, IDDOFILTER);
            break;
    }
}

t_dlg *select_filter(t_x11 *x11, t_gmx *gmx)
{
    static const char *title = "Group";
    static const char *dummy = "\"FALSE\"";
    static const char *ok    = "\"Ok\"";
    FILE              *tmp;
    t_dlg             *dlg;
    char               tmpfile[STRLEN];
    int                i, j, k, len, tlen, ht, ncol, nrow, x0;

    len = strlen(title);
    for (i = 0; (i < (int)gmx->filter->grps->nr); i++)
    {
        len = std::max(len, (int)strlen(gmx->filter->grpnames[i]));
    }
    len += 2;

    ncol = 1+(gmx->filter->grps->nr / 15);
    nrow = gmx->filter->grps->nr/ncol;
    if (nrow*ncol < gmx->filter->grps->nr)
    {
        nrow++;
    }
    if (ncol > 1)
    {
        ht = 1+(nrow+1)*2+3;
    }
    else
    {
        ht = 1+(gmx->filter->grps->nr+1)*2+3;
    }
    strcpy(tmpfile, "filterXXXXXX");
    gmx_tmpnam(tmpfile);
#ifdef DEBUG
    fprintf(stderr, "file: %s\n", tmpfile);
#endif
    if ((tmp = fopen(tmpfile, "w")) == NULL)
    {
        sprintf(tmpfile, "%ctmp%cfilterXXXXXX", DIR_SEPARATOR, DIR_SEPARATOR);
        gmx_tmpnam(tmpfile);
        if ((tmp = fopen(tmpfile, "w")) == NULL)
        {
            gmx_fatal(FARGS, "Can not open tmp file %s", tmpfile);
        }
    }
    tlen = 1+ncol*(1+len);
    fprintf(tmp, "grid %d %d {\n\n", tlen, ht);

    for (k = j = 0, x0 = 1; (j < ncol); j++, x0 += len+1)
    {
        fprintf(tmp, "group \"%s-%d\" %d 1 %d %d {\n", title, j+1, x0, len, ht-5);
        for (i = 0; (i < nrow) && (k < gmx->filter->grps->nr); i++, k++)
        {
            if (!gmx->filter->bDisable[k])
            {
                fprintf(tmp, "checkbox \"%s\" \"%d\" %s %s %s\n",
                        gmx->filter->grpnames[k], k, dummy, dummy, dummy);
            }
            else
            {
                fprintf(tmp, "statictext { \"  %s\" } \"%d\" %s %s %s\n",
                        gmx->filter->grpnames[k], k, dummy, dummy, dummy);
            }
        }
        fprintf(tmp, "}\n\n");
    }
    fprintf(tmp, "simple 1 %d %d 2 {\n", ht-3, tlen-2);
    fprintf(tmp, "defbutton %s %s %s %s %s\n", ok, ok, dummy, dummy, dummy);
    fprintf(tmp, "}\n\n}\n");
    fclose(tmp);

    dlg = ReadDlg(x11, gmx->wd->self, title, tmpfile,
                  0, 0, true, false, FilterCB, gmx);

    remove(tmpfile);

    return dlg;
}
