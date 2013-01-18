/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include <math.h>
#include <smalloc.h>
#include <sysstuff.h>
#include <string2.h>
#include <Xstuff.h>
#include "xutil.h"
#include "futil.h"

typedef struct {
    XColor    xc;
    t_windata wd;
} t_col;

typedef struct {
    t_windata  wd;
    t_windata  but;
    int        ncol;
    t_col     *col;
} t_sc;

static gmx_bool ColCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_col     *col;
    XColor    *xc;
    int        r, g, b;
    t_windata *wd;

    col = (t_col *)data;
    xc  = &(col->xc);
    wd  = &(col->wd);
    switch (event->type)
    {
        case ButtonPress:
            r = xc->red >> 8;
            g = xc->green >> 8;
            b = xc->blue >> 8;
            printf("%6d%6d%6d\t\t%-20s\n", r, g, b, wd->text);
            break;
    }
    return FALSE;
}

static gmx_bool BCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_sc  *sc;

    sc = (t_sc *)data;
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, sc->col[sc->ncol-1].xc.pixel);
            TextInWin(x11, &sc->but, sc->but.text, eXCenter, eYCenter);
            break;
        case ConfigureNotify:
            sc->but.width  = event->xconfigure.width;
            sc->but.height = event->xconfigure.height;
            break;
        case EnterNotify:
            LightBorder(x11->disp, sc->but.self, sc->col[sc->ncol-1].xc.pixel);
            break;
        case LeaveNotify:
            LightBorder(x11->disp, sc->but.self, sc->col[0].xc.pixel);
            break;
        case ButtonPress:
            x11->UnRegisterCallback(x11, sc->wd.self);
            x11->UnRegisterCallback(x11, sc->but.self);
            return TRUE;
    }
    return FALSE;
}

static gmx_bool scCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_sc  *sc;
    int    i;
    int    nx, ny, X, Y, Dx, Dy, th;
    float  x, y, dx, dy;

    th = XTextHeight(x11->font)+6;
    sc = (t_sc *)data;
    switch (event->type)
    {
        case ConfigureNotify:
            x = sc->wd.width = event->xconfigure.width;
            y = sc->wd.height = event->xconfigure.height-th;
            if ((x >= 100) && (y >= 100))
            {
#ifdef DEBUG
                printf("Disp: %x, But: %x, Font: %x x: %d, th: %d\n",
                       x11->disp, sc->but.self, x11->font, (int)x, th);
#endif
                XResizeWindow(x11->disp, sc->but.self, (int)x-4, th-4);
                nx = sqrt(sc->ncol);
                ny = sc->ncol/nx;
                while (nx*ny < sc->ncol)
                {
                    ny++;
                }
                dx = ((float)x) / nx;
                dy = ((float)y) / ny;
                for (i = 0; (i < sc->ncol); i++)
                {
                    X  = x = (i%nx) * dx;
                    Y  = y = (i/nx) * dy;
                    Dx = (int)(x+dx)-X;
                    Dy = (int)(y+dy)-Y;
                    XMoveWindow(x11->disp, sc->col[i].wd.self, X, th+Y);
                    XResizeWindow(x11->disp, sc->col[i].wd.self, Dx, Dy);
                }
            }
            break;
    }
    return FALSE;
}

static int col_comp(const void *c1, const void *c2)
{
    XColor *x1, *x2;
    int     dr, dg;
    int     diff;

    x1 = &(((t_col *)c1)->xc);
    x2 = &(((t_col *)c2)->xc);

    /* sort on intensity */
    diff = (x1->red+x1->green+x1->blue)-(x2->red+x2->green+x2->blue);

    if (diff == 0)
    {
        if ((dr = (x1->red-x2->red)) != 0)
        {
            return dr;
        }
        if ((dg = (x1->green-x2->green)) != 0)
        {
            return dg;
        }
        return x1->blue-x2->blue;
    }
    else
    {
        return diff;
    }
}

static void read_col(t_x11 *x11, t_sc *sc, char *rgb)
{
    FILE   *fp;
    char    buf[STRLEN], name[STRLEN], *dummy;
    int     i, j;
    int     r, g, b;
    XColor  xc, exact;
    t_col  *cols;

    if ((fp = fopen(rgb, "r")) == NULL)
    {
        perror(rgb);
        exit(1);
    }
    do
    {
        dummy = NULL;
        if (fscanf(fp, "%d%d%d", &r, &g, &b) == 3)
        {
            dummy = fgets2(buf, STRLEN-1, fp);
            if (dummy)
            {
                trim(buf);
                /* Filter out colours with names of two words */
                sscanf(buf, "%s", name);
                /* Filter out duplicate colours (grey and gray) */
                if (strstr(name, "gray") == 0)
                {
                    if (XAllocNamedColor(x11->disp, x11->cmap, name, &xc, &exact))
                    {
                        srenew(sc->col, ++sc->ncol);
                        sc->col[sc->ncol-1].xc = xc;
#ifdef DEBUG
                        printf("color %3d: %s\n", sc->ncol, buf);
#endif
                        InitWin(&(sc->col[sc->ncol-1].wd), 0, 0, 1, 1, 0, name);
                    }
                }
            }
        }
    }
    while (dummy);
    fclose(fp);
    if (sc->ncol)
    {
        qsort(sc->col, sc->ncol, sizeof(sc->col[0]), col_comp);
    }
    /* Now filter out doubles */
    cols = &(sc->col[0]);
    for (i = 1, j = 0; (i < sc->ncol); i++)
    {
        if ((cols[i].xc.red   != cols[j].xc.red) ||
            (cols[i].xc.green != cols[j].xc.green) ||
            (cols[i].xc.blue  != cols[j].xc.blue))
        {
            j++;
            cols[j] = cols[i];
        }
    }
    sc->ncol = j;
}

static t_sc *init_sc(t_x11 *x11, Window Parent, char *rgb)
{
    t_sc   *sc;
    Window  w;
    int     i;

    snew(sc, 1);
    InitWin(&sc->wd, 0, 0, 400, 300, 1, "Show Colours");
    sc->wd.self = XCreateSimpleWindow(x11->disp, Parent, sc->wd.x, sc->wd.y,
                                      sc->wd.width, sc->wd.height,
                                      sc->wd.bwidth, WHITE, BLACK);
    x11->RegisterCallback(x11, sc->wd.self, Parent, scCallBack, sc);
    x11->SetInputMask(x11, sc->wd.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask);
    InitWin(&sc->but, 0, 0, sc->wd.width-2, XTextHeight(x11->font)+2, 1, "Quit");
    sc->but.self = XCreateSimpleWindow(x11->disp, sc->wd.self, sc->but.x, sc->but.y,
                                       sc->but.width, sc->but.height,
                                       sc->but.bwidth, BLACK, BLACK);
    x11->RegisterCallback(x11, sc->but.self, sc->but.self, BCallBack, sc);
    x11->SetInputMask(x11, sc->but.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask | EnterWindowMask |
                      LeaveWindowMask);

    read_col(x11, sc, rgb);
    fprintf(stderr, "%d colors found\n", sc->ncol);
    fprintf(stderr, "%6s%6s%6s\t\t%-20s\n", "Red", "Green", "Blue", "name");
    for (i = 0; (i < sc->ncol); i++)
    {
        sc->col[i].wd.self = w = XCreateSimpleWindow(x11->disp, sc->wd.self, 0, 0, 1, 1, 0,
                                                     BLACK, sc->col[i].xc.pixel);
        x11->RegisterCallback(x11, w, sc->wd.self, ColCallBack, &(sc->col[i]));
        x11->SetInputMask(x11, w, ButtonPressMask);
    }

    return sc;
}

int
main(int argc, char *argv[])
{
    t_x11 *x11;
    t_sc  *sc;
    char  *fn;

    x11 = GetX11(&argc, argv);
    if (argc > 1)
    {
        fn = argv[1];
    }
    else
    {
        fn = "/usr/lib/X11/rgb.txt";
    }
    if (!gmx_fexist(fn))
    {
        fprintf(stderr, "Usage: %s rgb.txt\n", argv[0]);
        fprintf(stderr, "rgb.txt is usually somewhere in your X windows directories.\n");
        exit(1);
    }
    sc = init_sc(x11, x11->root, fn);
    XMapWindow(x11->disp, sc->wd.self);
    XMapSubwindows(x11->disp, sc->wd.self);
    x11->MainLoop(x11);
    x11->CleanUp(x11);

    return 0;
}
