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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "macros.h"
#include "Xstuff.h"
#include "copyrite.h"
#include "xutil.h"
#include "futil.h"
#include "x11.h"
#include "smalloc.h"
#include "statutil.h"
#include "rama.bm"
#include "nrama.h"

#define MAXDEG 360

enum {
    ebQuit, ebStart, ebStop, ebRewind, ebGly, ebutNR
};

static real scx, scy;
#define SX(x) ((int)(((x)+M_PI)*scx))
#define SY(y) ((int)((M_PI-(y))*scy))

enum {
    esStop, esGo, esNR
};

typedef struct {
    int            status;
    gmx_bool       bShowGly;
    gmx_bool      *bIsGly;
    t_windata      wd;
    t_windata      xrwd;
    t_xrama       *xr;
    t_windata      but[ebutNR];
} t_app;

static void plot_pp(t_x11 *x11, Window w, t_phipsi *pp, t_dih dih[])
{
    int x0, y0;
    int th = (XTextHeight(x11->font)+6)/2;

    x0 = SX(dih[pp->iphi].ang);
    y0 = SY(dih[pp->ipsi].ang);
    XFillRectangle(x11->disp, w, x11->gc, x0-1, y0-1, 4, 4);
    /* Draw Label ? */
    if (pp->bShow)
    {
        TextInRect(x11, w, pp->label, x0+6, y0-th, 30, 2*th, eXLeft, eYCenter);
    }
}

static gmx_bool label_pp(t_x11 *x11, Window w, int npp, t_phipsi pp[],
                         t_dih dih[], int mx, int my)
{
    int d, md, x0, y0;
    int i, imin;

    imin = -1;
    md   = 16;
    for (i = 0; (i < npp); i++)
    {
        x0 = SX(dih[pp[i].iphi].ang);
        y0 = SY(dih[pp[i].ipsi].ang);
        d  = (mx-x0)*(mx-x0)+(my-y0)*(my-y0);
        if (d < md)
        {
            md   = d;
            imin = i;
        }
    }
    if (imin != -1)
    {
        pp[imin].bShow = !pp[imin].bShow;
        return TRUE;
    }
    return FALSE;
}

static gmx_bool xrCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_app   *app;
    t_xrama *xr;
    char     buf[256];
    int      i;

    (void)XTextHeight(x11->font);
    app = (t_app *)data;
    xr  = app->xr;
    scx = app->xrwd.width/(2.0*M_PI);
    scy = app->xrwd.height/(2.0*M_PI);
    switch (event->type)
    {
        case Expose:
            XClearWindow(x11->disp, app->xrwd.self);
            XDrawLine(x11->disp, app->xrwd.self, x11->gc,
                      SX(0), SY(-M_PI)+1, SX(0), SY(M_PI)-1);
            XDrawLine(x11->disp, app->xrwd.self, x11->gc,
                      SX(-M_PI)+1, SY(0), SX(M_PI)-1, SY(0));
            TextInRect(x11, app->xrwd.self, "Phi", SX(M_PI)-50, SY(0)+4, 46, 20, eXRight, eYTop);
            TextInRect(x11, app->xrwd.self, "Psi", SX(0)+4, 4, 46, 20, eXLeft, eYTop);
            for (i = 0; (i < xr->npp); i++)
            {
                if (app->bShowGly || !app->bIsGly[i])
                {
                    plot_pp(x11, app->xrwd.self, &(xr->pp[i]), xr->dih);
                }
            }
            break;
        case ButtonPress:
            if (label_pp(x11, app->xrwd.self, xr->npp, xr->pp, xr->dih,
                         event->xbutton.x, event->xbutton.y))
            {
                ExposeWin(x11->disp, app->xrwd.self);
            }
            break;
        case ConfigureNotify:
            app->xrwd.width  = event->xconfigure.width;
            app->xrwd.height = event->xconfigure.height;
            break;
    }
    if (app->status == esGo)
    {
        if (!new_data(app->xr))
        {
            app->status = ebStop;
        }
        else
        {
            ExposeWin(x11->disp, app->xrwd.self);
            sprintf(buf, "Rama: t=%.2f", app->xr->t);
            XSetStandardProperties(x11->disp, app->wd.self, buf,
                                   "Rama", 0, NULL, 0, NULL);

        }
    }
    return FALSE;
}

static gmx_bool appCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_app  *app;
    int     win;

    app = (t_app *)data;
    for (win = 0; (win < ebutNR); win++)
    {
        if (app->but[win].self == w)
        {
            break;
        }
    }
    if (win == ebutNR)
    {
        return FALSE;
    }

    switch (event->type)
    {
        case Expose:
            TextInWin(x11, &(app->but[win]), app->but[win].text, eXCenter, eYCenter);
            break;
        case ButtonPress:
            switch (win)
            {
                case ebQuit:
                    exit(1);
                case ebStart:
                    app->status = esGo;
                    ExposeWin(x11->disp, app->xrwd.self);
                    break;
                case ebStop:
                    app->status = esStop;
                    break;
                case ebRewind:
                    rewind_trj(app->xr->traj);
                    break;
                case ebGly:
                    app->bShowGly = !app->bShowGly;
                    ExposeWin(x11->disp, app->xrwd.self);
                    break;
                default:
                    XBell(x11->disp, 50);
                    break;
            }
            break;
    }
    return FALSE;
}

static void size_app(t_x11 *x11, t_app *app)
{
    int i, dx, th;

    th = XTextHeight(x11->font)+4;
    dx = app->wd.width/ebutNR;
    for (i = 0; (i < ebutNR); i++)
    {
        app->but[i].width  = dx-4;
        app->but[i].height = th+4;
        XMoveResizeWindow(x11->disp, app->but[i].self, i*dx+2, 2, dx-4, th+4);
    }
    XMoveResizeWindow(x11->disp, app->xrwd.self, 2, th+10,
                      app->wd.width-6, app->wd.height-th-10-4);
}

static gmx_bool mainCallBack(struct t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_app *app;
    int    wt, ht;

    app = (t_app *)data;
    switch (event->type)
    {
        case ConfigureNotify:
            wt = event->xconfigure.width;
            ht = event->xconfigure.height;
            if ((app->wd.width != wt) || (app->wd.height != ht))
            {
                fprintf(stderr, "New wxh = %dx%d\n", wt, ht);
                app->wd.width  = wt;
                app->wd.height = ht;
                size_app(x11, app);
            }
            break;
    }
    return FALSE;
}

static t_xrama *init_xrama(t_x11 *x11, Window Parent, int y0, t_app *app)
{
    t_xrama *xr;

    snew(xr, 1);

    InitWin(&(app->xrwd), 2, y0, MAXDEG+1, MAXDEG+1, 1, "Ramachandran Movie");
    app->xrwd.self = XCreateSimpleWindow(x11->disp, Parent, app->xrwd.x, app->xrwd.y,
                                         app->xrwd.width, app->xrwd.height,
                                         app->xrwd.bwidth, x11->fg, x11->bg);
    x11->RegisterCallback(x11, app->xrwd.self, Parent, xrCallBack, app);
    x11->SetInputMask(x11, app->xrwd.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask);

    return xr;
}

static t_app *init_app(t_x11 *x11, int argc, char *argv[])
{
    static const char *but_nm[ebutNR] = { "Quit", "Start", "Stop", "Rewind", "Toggle Gly" };
    XSizeHints         hints;
    Pixmap             pm;
    int                th;
    t_app             *app;
    t_windata         *wd;
    int                i, dx;

    snew(app, 1);
    th = XTextHeight(x11->font)+4;
    InitWin(&(app->wd), 0, 0, MAXDEG+6, MAXDEG+6+th+6, 0, "Ramachandran Movie");
    dx           = app->wd.width/ebutNR;
    app->wd.self = XCreateSimpleWindow(x11->disp, x11->root, app->wd.x, app->wd.y,
                                       app->wd.width, app->wd.height,
                                       app->wd.bwidth, x11->fg, x11->bg);
    x11->RegisterCallback(x11, app->wd.self, x11->root, mainCallBack, app);
    x11->SetInputMask(x11, app->wd.self, StructureNotifyMask);
    hints.flags = 0;
    pm          = XCreatePixmapFromBitmapData(x11->disp, x11->root, (char *)rama_bits,
                                              rama_width, rama_height, WHITE, BLACK, 1);
    XSetStandardProperties(x11->disp, app->wd.self, app->wd.text,
                           "Rama", pm, argv, argc, &hints);
    x11->RegisterCallback(x11, app->wd.self, x11->root, appCallBack, app);
    x11->SetInputMask(x11, app->wd.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask);

    app->xr = init_xrama(x11, app->wd.self, th+6, app);
    for (i = 0; (i < ebutNR); i++)
    {
        wd = &(app->but[i]);
        InitWin(wd, i*dx+2, 2, dx-4, th, 1, but_nm[i]);
        wd->self = XCreateSimpleWindow(x11->disp, app->wd.self,
                                       wd->x, wd->y, wd->width, wd->height,
                                       wd->bwidth, x11->fg, x11->bg);
        x11->RegisterCallback(x11, wd->self, app->wd.self, appCallBack, app);
        x11->SetInputMask(x11, wd->self, ButtonPressMask | ExposureMask);
    }
    return app;
}

static void mk_gly(t_app *app)
{
    int i;

    snew(app->bIsGly, app->xr->npp);
    for (i = 0; (i < app->xr->npp); i++)
    {
        if (strstr(app->xr->pp[i].label, "GLY") != NULL)
        {
            app->bIsGly[i] = TRUE;
        }
    }
}

int main(int argc, char *argv[])
{
    const char  *desc[] = {
        "[TT]g_xrama[tt] shows a Ramachandran movie, that is, it shows",
        "the Phi/Psi angles as a function of time in an X-Window.[PAR]"
        "Static Phi/Psi plots for printing can be made with [TT]g_rama[tt].[PAR]",
        "Some of the more common X command line options can be used:[BR]",
        "[TT]-bg[tt], [TT]-fg[tt] change colors, [TT]-font fontname[tt], changes the font."
    };

    output_env_t oenv;
    t_x11       *x11;
    t_topology  *ramatop;
    t_app       *app;
    t_filenm     fnm[] = {
        { efTRX, "-f", NULL, ffREAD },
        { efTPX, NULL, NULL, ffREAD }
    };
#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm, 0, NULL,
                      asize(desc), desc, 0, NULL, &oenv);


    if ((x11 = GetX11(&argc, argv)) == NULL)
    {
        fprintf(stderr, "Can't open display, set your DISPLAY environment variable\n");
        exit(1);
    }
    XSetForeground(x11->disp, x11->gc, x11->fg);
    app = init_app(x11, argc, argv);

    ramatop = init_rama(oenv, ftp2fn(efTRX, NFILE, fnm), ftp2fn(efTPX, NFILE, fnm),
                        app->xr, 3);
    mk_gly(app);

    XMapWindow(x11->disp, app->wd.self);
    XMapSubwindows(x11->disp, app->wd.self);
    x11->MainLoop(x11);
    x11->CleanUp(x11);

    thanx(stderr);

    return 0;
}
