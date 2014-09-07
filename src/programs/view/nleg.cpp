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

#include "nleg.h"

#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/types/rgb.h"
#include "gromacs/utility/smalloc.h"

#include "buttons.h"

typedef struct {
    const char    *tp;
    unsigned long *col;
    t_rgb          rgb;
} t_atomcolor;

static t_atomcolor ac[] = {
    { "O",  &LIGHTRED,     { 1,  0,  0   } },
    { "N",  &LIGHTCYAN,    { 0,  0,  1   } },
    { "NA", &LIGHTGREY,    { 0.6, 0.6, 0.6 } },
    { "S",  &YELLOW,       { 1,  1,  0   } },
    { "C",  &LIGHTGREEN,   { 0,  1,  0   } },
    { "CL", &VIOLET,       { 1,  0,  1   } },
    { "F",  &LIGHTGREY,    { 0.6, 0.6, 0.6 } },
    { "Z",  &LIGHTGREY,    { 0.6, 0.6, 0.6 } },
    { "P",  &LIGHTBLUE,    { 0.4, 0.4, 1.0 } },
    { "H",  &WHITE,        { 0.8, 0.8, 0.8 } }
};
#define NAC asize(ac)

int search_ac(const char *type)
{
    unsigned int i, nb, mij, best, besti;

    best  = 0;
    besti = 0;
    if (NULL != type)
    {
        for (i = 0; (i < NAC); i++)
        {
            mij = std::min((int)strlen(type), (int)strlen(ac[i].tp));
            for (nb = 0; (nb < mij); nb++)
            {
                if (type[nb] != ac[i].tp[nb])
                {
                    break;
                }
            }
            if (nb > best)
            {
                best  = nb;
                besti = i;
            }
        }
    }
    return besti;
}

unsigned long Type2Color(const char *type)
{
    int i;

    i = search_ac(type);

    return *(ac[i].col);
}

t_rgb *Type2RGB(const char *type)
{
    int i;

    i = search_ac(type);

    return &(ac[i].rgb);
}

void DrawLegend(t_x11 *x11, t_windata *Win)
{
#define NLAB 6
#define COLS 3
    static const char *lab[NLAB] = { "C", "O", "H", "S", "N", "P" };
    int                i, i0, dh, dw, w, y, x1, x0;
    unsigned long      cind;
    real               h_2;

    XClearWindow(x11->disp, Win->self);
    w   = Win->width;
    h_2 = Win->height/(2.0*NLAB/COLS);
    dh  = h_2-2;
    dw  = dh;

    for (i = 0; (i < NLAB); i++)
    {
        i0   = i % (NLAB/COLS);
        x0   = (i / (NLAB/COLS))*(Win->width/COLS)+AIR;
        x1   = x0+2*dw+AIR;
        cind = Type2Color(lab[i]);
        XSetForeground(x11->disp, x11->gc, cind);
        y = ((2*i0+1)*h_2);
        XFillRectangle (x11->disp, Win->self, x11->gc, x0, y-dh, 2*dw, 2*dh);
        XSetForeground(x11->disp, x11->gc, WHITE);
        TextInRect(x11, Win->self, lab[i], x1, y-dh, w-x1, 2*dh,
                   eXLeft, eYCenter);
    }
    XSetForeground(x11->disp, x11->gc, x11->fg);
}

static bool LegWCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_legendwin *lw;

    lw = (t_legendwin *)data;
    switch (event->type)
    {
        case Expose:
            DrawLegend(x11, &lw->wd);
            break;
        default:
            break;
    }
    return false;
}

t_legendwin *init_legw(t_x11 *x11, Window Parent,
                       int x, int y, int width, int height,
                       unsigned long fg, unsigned long bg)
{
    t_legendwin *lw;

    snew(lw, 1);
    InitWin(&lw->wd, x, y, width, height, 1, "Legend Window");
    lw->wd.self = XCreateSimpleWindow(x11->disp, Parent, x, y, 1, 1, 1, fg, bg);
    x11->RegisterCallback(x11, lw->wd.self, Parent, LegWCallBack, lw);
    x11->SetInputMask(x11, lw->wd.self, ExposureMask);

    return lw;
}

void map_legw(t_x11 *x11, t_legendwin *lw)
{
    XMapWindow(x11->disp, lw->wd.self);
}


void done_legw(t_x11 *x11, t_legendwin *lw)
{
    x11->UnRegisterCallback(x11, lw->wd.self);
    sfree(lw);
}
