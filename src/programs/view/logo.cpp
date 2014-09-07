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

#include "logo.h"

#include <stdlib.h>

#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "Xstuff.h"
#include "xutil.h"

typedef struct {
    int            x, y, rad;
    unsigned long *col;
} t_circle;

typedef struct {
    const char  *text;
    int          y, h;
    XFontStruct *fnt;
} t_mess;

void show_logo(t_x11 *x11, t_logo *logo)
{
    XMapWindow(x11->disp, logo->wd.self);
    XMapSubwindows(x11->disp, logo->wd.self);
}

void hide_logo(t_x11 *x11, t_logo *logo)
{
    XUnmapWindow(x11->disp, logo->wd.self);
}

static bool LogoCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    /* Assume window is 100x110 */
    static bool bFirst = true;
#define CSIZE 9
#define NSIZE 8
#define OSIZE 9
#define HSIZE 7
#define YOFFS 30
    static t_circle c[] = {
        { 10, YOFFS+12, CSIZE, &LIGHTGREEN },
        { 20, YOFFS+22, CSIZE, &LIGHTGREEN },
        { 20, YOFFS+34, OSIZE, &LIGHTRED   },
        { 30, YOFFS+12, NSIZE, &LIGHTCYAN  },
        { 30, YOFFS+ 2, HSIZE, &WHITE     },
        { 40, YOFFS+22, CSIZE, &LIGHTGREEN },
        { 40, YOFFS+34, CSIZE, &LIGHTGREEN },
        { 50, YOFFS+12, CSIZE, &LIGHTGREEN },
        { 50, YOFFS,    OSIZE, &LIGHTRED   },
        { 60, YOFFS+22, NSIZE, &LIGHTCYAN  },
        { 60, YOFFS+32, HSIZE, &WHITE     },
        { 70, YOFFS+12, CSIZE, &LIGHTGREEN },
        { 80, YOFFS+22, CSIZE, &LIGHTGREEN },
        { 80, YOFFS+34, OSIZE, &LIGHTRED   },
        { 90, YOFFS+12, NSIZE, &LIGHTCYAN  },
        { 90, YOFFS+ 2, HSIZE, &WHITE      },
        {100, YOFFS+22, CSIZE, &LIGHTGREEN }
    };
    static int      lines[] = {
        0, 1, 1, 2, 1, 3, 3, 4, 3, 5, 5, 6, 5, 7, 7, 8, 7, 9,
        9, 10, 9, 11, 11, 12, 12, 13, 12, 14, 14, 15, 14, 16
    };
#define COFFS 70
    static t_mess   Mess[] = {
        { "GROMACS",                         0,       20, NULL },
        { NULL,                             16,        9, NULL },
        { "Copyright (c) 1991-2013",        COFFS+ 2,  9, NULL },
        { "D.v.d.Spoel, E.Lindahl, B.Hess", COFFS+11,  9, NULL },
        { "& Groningen University ",        COFFS+20,  9, NULL },
        { "click to dismiss",               COFFS+31,  8, NULL }
    };
#define NMESS asize(Mess)
    int             i;
    t_logo         *logo;
    t_windata      *wd;

    logo = (t_logo *)data;
    wd   = &(logo->wd);
    if (bFirst)
    {
        const real wfac = wd->width/110.0;
        const real hfac = wd->height/110.0;
        for (i = 0; (i < asize(c)); i++)
        {
            c[i].x *= wfac;
            c[i].y *= hfac;
        }
        Mess[1].text = GromacsVersion();
        for (i = 0; (i < NMESS); i++)
        {
            Mess[i].y  *= hfac;
            Mess[i].h  *= hfac;
            Mess[i].fnt = (i == 0) ? logo->bigfont : (i == NMESS-1) ? x11->font :
                logo->smallfont;
        }
        bFirst = false;
    }
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, WHITE);
            XSetLineAttributes(x11->disp, x11->gc, 3, LineSolid, CapNotLast, JoinRound);
            for (i = 0; (i < asize(lines)); i += 2)
            {
                XDrawLine(x11->disp, wd->self, x11->gc,
                          c[lines[i]].x, c[lines[i]].y, c[lines[i+1]].x, c[lines[i+1]].y);
            }
            XSetLineAttributes(x11->disp, x11->gc, 1, LineSolid, CapNotLast, JoinRound);
            for (i = 0; (i < asize(c)); i++)
            {
                XSetForeground(x11->disp, x11->gc, *(c[i].col));
                XFillCircle(x11->disp, wd->self, x11->gc, c[i].x, c[i].y, c[i].rad);
            }
            XSetForeground(x11->disp, x11->gc, BLACK);
            XDrawRectangle(x11->disp, wd->self, x11->gc, 2, 2, wd->width-5, wd->height-5);
            for (i = 0; (i < NMESS); i++)
            {
                SpecialTextInRect(x11, Mess[i].fnt, wd->self, Mess[i].text,
                                  0, Mess[i].y, wd->width, Mess[i].h,
                                  eXCenter, eYCenter);
            }
            XSetForeground(x11->disp, x11->gc, x11->fg);
            break;
        case ButtonPress:
            hide_logo(x11, logo);
            return logo->bQuitOnClick;
            break;
        default:
            break;
    }

    return false;
}

t_logo *init_logo(t_x11 *x11, Window parent, bool bQuitOnClick)
{
    static const char *bfname[] = {
        "-b&h-lucida-bold-i-normal-sans-34-240-100-100-p-215-iso8859-1",
        "-b&h-lucida-bold-i-normal-sans-26-190-100-100-p-166-iso8859-1",
        "lucidasans-bolditalic-24",
        "lucidasans-italic-24",
        "10x20",
        "fixed"
    };
#define NBF asize(bfname)
    static const char *sfname[] = {
        "lucidasans-bold-18",
        "10x20",
        "fixed"
    };
#define NSF asize(sfname)
    unsigned int       i;
    unsigned long      bg;
    char              *newcol;
    t_logo            *logo;

    snew(logo, 1);
    logo->bQuitOnClick = bQuitOnClick;
    InitWin(&logo->wd, 0, 0, 360, 270, 1, "GROMACS");
    bg = LIGHTGREY;
    if ((newcol = getenv("GMX_LOGO_COLOR")) != NULL)
    {
        GetNamedColor(x11, newcol, &bg);
    }
    logo->wd.self = XCreateSimpleWindow(x11->disp, parent,
                                        logo->wd.x, logo->wd.y,
                                        logo->wd.width, logo->wd.height,
                                        logo->wd.bwidth, WHITE, bg);
    for (i = 0, logo->bigfont = NULL; (i < NBF); i++)
    {
        if ((logo->bigfont = XLoadQueryFont(x11->disp, bfname[i])) != NULL)
        {
            break;
        }
    }
    if (i == NBF)
    {
        perror(bfname[i-1]);
        exit(1);
    }
#ifdef DEBUG
    fprintf(stderr, "Big Logofont: %s\n", bfname[i]);
#endif
    for (i = 0, logo->smallfont = NULL; (i < NSF); i++)
    {
        if ((logo->smallfont = XLoadQueryFont(x11->disp, sfname[i])) != NULL)
        {
            break;
        }
    }
    if (i == NSF)
    {
        perror(sfname[i-1]);
        exit(1);
    }
#ifdef DEBUG
    fprintf(stderr, "Small Logofont: %s\n", sfname[i]);
#endif
    x11->RegisterCallback(x11, logo->wd.self, parent, LogoCallBack, logo);
    x11->SetInputMask(x11, logo->wd.self, ButtonPressMask | ExposureMask);

    return logo;
}

void done_logo(t_x11 *x11, t_logo *logo)
{
    x11->UnRegisterCallback(x11, logo->wd.self);
}
