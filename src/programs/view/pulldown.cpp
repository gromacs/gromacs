/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2013, The GROMACS development team.
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
#include "gmxpre.h"

#include "pulldown.h"

#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/smalloc.h"

#include "popup.h"
#include "x11.h"

static bool PDCallBack(t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_pulldown *pd;
    int         i, x, x1, y, nsel;

    pd = (t_pulldown *)data;
    y  = pd->wd.height;
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, x11->fg);
            XDrawLine(x11->disp, w, x11->gc, 0, y-1, pd->wd.width, y-1);
            for (i = 0; (i < pd->nmenu); i++)
            {
                XDrawString(x11->disp, pd->wd.self, x11->gc, pd->xpos[i], x11->font->ascent,
                            pd->title[i], strlen(pd->title[i]));
            }
            break;
        case ButtonPress:
            if (pd->nsel == -1)
            {
                x = event->xbutton.x;
                for (nsel = 0; (pd->xpos[nsel+1] < x) && (nsel < pd->nmenu-1); nsel++)
                {
                    ;
                }
                pd->nsel = nsel;
                x1       = std::max(0, std::min(pd_width(pd)-menu_width(pd->m[nsel]), pd->xpos[nsel]));
                show_menu(x11, pd->m[nsel], x1, y+1, false);
            }
            break;
        case ButtonRelease:
            hide_pd(x11, pd);
            break;
        default:
            break;
    }
    return false;
}

t_pulldown *init_pd(t_x11 *x11, Window Parent, int width,
                    unsigned long fg, unsigned long bg,
                    int nmenu, int *nsub, t_mentry *ent[], const char **title)
{
    t_pulldown *pd;
    int         i;

    snew(pd, 1);
    pd->title = title;
    pd->nmenu = nmenu;
    pd->nsel  = -1;
    snew(pd->m, nmenu);
    snew(pd->xpos, nmenu+1);
    pd->xpos[0] = 5;
    for (i = 1; (i <= nmenu); i++)
    {
        pd->xpos[i] = 20+pd->xpos[i-1]+
            XTextWidth(x11->font, title[i-1], strlen(title[i-1]));
    }
    if (pd->xpos[nmenu] > width)
    {
        printf("Menu too wide\n");
    }

    InitWin(&(pd->wd), 0, 0, width, XTextHeight(x11->font)+2, 0, "PullDown");
    pd->wd.self = XCreateSimpleWindow(x11->disp, Parent,
                                      pd->wd.x, pd->wd.y,
                                      pd->wd.width, pd->wd.height,
                                      pd->wd.bwidth, fg, bg);
    x11->RegisterCallback(x11, pd->wd.self, Parent, PDCallBack, pd);
    x11->SetInputMask(x11, pd->wd.self, ExposureMask | ButtonPressMask |
                      OwnerGrabButtonMask | ButtonReleaseMask);
    XMapWindow(x11->disp, pd->wd.self);

    for (i = 0; (i < nmenu); i++)
    {
        pd->m[i] = init_menu(x11, Parent, fg, bg, nsub[i], ent[i], 1);
    }

    return pd;
}

void hide_pd(t_x11 *x11, t_pulldown *pd)
{
    if (pd->nsel != -1)
    {
        hide_menu(x11, pd->m[pd->nsel]);
    }
    pd->nsel = -1;
}

void check_pd_item(t_pulldown *pd, int nreturn, bool bStatus)
{
    int i;

    for (i = 0; (i < pd->nmenu); i++)
    {
        check_menu_item(pd->m[i], nreturn, bStatus);
    }
}

void done_pd(t_x11 *x11, t_pulldown *pd)
{
    int i;

    for (i = 0; (i < pd->nmenu); i++)
    {
        done_menu(x11, pd->m[i]);
    }
    x11->UnRegisterCallback(x11, pd->wd.self);
    sfree(pd->m);
    sfree(pd->xpos);
}

int pd_width(t_pulldown *pd)
{
    int i, w;

    w = 0;
    for (i = 0; (i < pd->nmenu); i++)
    {
        w = std::max(w, menu_width(pd->m[i]));
    }
    w = std::max(w, pd->xpos[pd->nmenu]);
    return w;
}

int pd_height(t_pulldown *pd)
{
    int i, h;

    h = 0;
    for (i = 0; (i < pd->nmenu); i++)
    {
        h = std::max(h, menu_height(pd->m[i]));
    }

    return h;
}
