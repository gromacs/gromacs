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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "futil.h"
#include "macros.h"
#include "smalloc.h"
#include "xutil.h"
#include "copyrite.h"
#include "statutil.h"
#include "gmx_fatal.h"

/* Units are meter and second	*/

typedef struct {
    int           id;            /* Identification      */
    float         x, xold;
    float         v;             /* Position and velocity	*/
    float         vwanted;       /* Wants to drive at this speed	*/
    float         acc;           /* Acceleration			*/
    float         brake;         /* Break			*/
    int           lane, oldlane; /* Currently in lane		*/
    gmx_bool      bBrake;        /* Currently on the brakes	*/
    unsigned long col;           /* Colour			*/
    unsigned long roof;          /* Roof Colour			*/
} t_car;

typedef struct {
    int   nlane;    /* Number of lanes on highway	*/
    int   metres;   /* Road length			*/
    float dt;       /* Time step			*/
    float min_dist; /* Min distance cars can come	*/
    int   sleep;    /* How long to sleep in between updates */
} t_input;

static const char *Driving[]      = { "Start", "Stop"   };
static const char *Fogs[]         = { "Fog",  "No Fog" };
enum buttons                        {
    Quit,   StartStop,   Fog, NBUT
};
static const char *but_name[NBUT] = { "Quit", "Start", "Fog"  };

typedef struct {
    int           ncars;
    t_car        *cars;
    t_input       ir;
    int           step;
    gmx_bool      bDriving; /* Are we driving ?		*/
    gmx_bool      bFog;     /* Is it foggy ?		*/
    t_windata     main;
    t_windata     win;
    t_windata     but[NBUT];
} t_xhighway;

int read_input(t_x11 *x11, const char *fn, t_car **cars, t_input *ir)
{
    FILE   *in;
    int     i, n;
    char    buf[100], b2[100];
    t_car  *c;

    in = ffopen(fn, "r");
    if (fscanf(in, "%d %d %f %f  %d",
               &ir->nlane, &ir->metres, &ir->dt, &ir->min_dist, &ir->sleep) != 5)
    {
        gmx_fatal(FARGS, "Not enough parameters in %s line 1", fn);
    }
    if (fscanf(in, "%d", &n) != 1)
    {
        gmx_fatal(FARGS, "Not enough parameters in %s line 2", fn);
    }
    snew(*cars, n);

    for (i = 0; (i < n); i++)
    {
        c       = &((*cars)[i]);
        c->id   = i;
        c->lane = 0;
        if (fscanf(in, "%f %f %f %f %f %s %s", &(c->x), &(c->v), &(c->vwanted),
                   &(c->acc), &(c->brake), buf, b2) != 7)
        {
            gmx_fatal(FARGS, "Not enough parameters in %s line %d", fn, 3+i);
        }
        x11->GetNamedColor(x11, buf, &(c->col));
        x11->GetNamedColor(x11, b2, &(c->roof));
    }
    ffclose(in);

    return n;
}

static float get_dist(int ncars, t_car cars[], int which, gmx_bool bFog,
                      int dir, int lane, int metres, int *nearest)
{
    int   i, near;
    float dist, nd;

    if (dir < 0)
    {
        dist = -metres;
    }
    else
    {
        dist = metres;
    }
    near = -1;

    for (i = 0; (i < ncars); i++)
    {
        if ((i != which) && (cars[i].oldlane == lane))
        {
            nd = cars[i].xold-cars[which].xold;
            if ((nd < 0) && (dir > 0))
            {
                nd += metres;
            }
            else if ((nd > 0) && (dir < 0))
            {
                nd -= metres;
            }

            if (!bFog || (fabs(nd) < 50))
            {
                if (dir < 0)
                {
                    if (nd > dist)
                    {
                        dist = nd;
                        near = i;
                    }
                }
                else if (dir > 0)
                {
                    if (nd < dist)
                    {
                        dist = nd;
                        near = i;
                    }
                }
            }
        }
    }
    *nearest = near;
    return fabs(dist);
}

void simulate(t_x11 *x11, t_xhighway *xhw,
              int ncars, t_car cars[], t_input *ir)
{
    int   i, n_bef, n_bef1, n_beh;
    float dist, distf, distb;

    for (i = 0; (i < ncars); i++)
    {
        cars[i].xold    = cars[i].x;
        cars[i].oldlane = cars[i].lane;
    }
    for (i = 0; (i < ncars); i++)
    {
        cars[i].bBrake = FALSE;
        dist           = get_dist(ncars, cars, i, xhw->bFog,
                                  1, cars[i].lane, ir->metres, &n_bef);
        if (dist < ir->min_dist)
        {
            distf = get_dist(ncars, cars, i, xhw->bFog,
                             1, cars[i].lane+1, ir->metres, &n_bef1);
            distb = get_dist(ncars, cars, i, xhw->bFog,
                             -1, cars[i].lane+1, ir->metres, &n_beh);
            if ((cars[i].lane < ir->nlane-1) && (distb >= ir->min_dist) &&
                (distf >= ir->min_dist))
            {
                cars[i].lane += 1;
            }
            else
            {
                /* Use brakes */
                cars[i].v -= cars[i].brake*ir->dt;
                if (cars[i].v < 0)
                {
                    cars[i].v = 0;
                }
                if (n_bef != -1)
                {
                    if ((cars[i].v < cars[n_bef].v) && (dist > ir->min_dist/2))
                    {
                        cars[i].v = cars[n_bef].v;
                    }
                }
                cars[i].bBrake = TRUE;
            }
        }
        else if ((cars[i].lane > 0) && (cars[i].v == cars[i].vwanted))
        {
            /* Check if I can go right again */
            dist = get_dist(ncars, cars, i, xhw->bFog,
                            1, cars[i].lane-1, ir->metres, &n_bef);
            distb = get_dist(ncars, cars, i, xhw->bFog,
                             -1, cars[i].lane-1, ir->metres, &n_beh);
            if ((dist >= ir->min_dist) && (distb >= ir->min_dist))
            {
                cars[i].lane -= 1;
            }
        }

        cars[i].x += cars[i].v*ir->dt;
        if (cars[i].x > ir->metres)
        {
            cars[i].x -= ir->metres;
        }
        if (!cars[i].bBrake && (cars[i].v < cars[i].vwanted))
        {
            cars[i].v += cars[i].acc*ir->dt;
            if (cars[i].v > cars[i].vwanted)
            {
                cars[i].v = cars[i].vwanted;
            }
        }
    }
    /* Detect Crashes */
    /* Plot */
    usleep(xhw->ir.sleep);
    ExposeWin(x11->disp, xhw->win.self);
}

static void Configure(t_xhighway *xhw)
{
    Window self;
    int    i, h, w, dh;
    float  dw;

    dh = 20;
    h  = xhw->main.height;
    w  = xhw->main.width;
    dw = ((float)(w-2))/NBUT-4;
    for (i = 0; (i < NBUT); i++)
    {
        t_windata *wd = &(xhw->but[i]);

        self = wd->self;
        InitWin(wd, 2+i*(dw+4), 2, dw, dh, 1, but_name[i]);
        wd->self = self;
    }
    self = xhw->win.self;
    InitWin(&xhw->win, 2, dh+6, w-6, h-dh-10, 1, xhw->main.text);
    xhw->win.self = self;
}

static void draw_car(Display *disp, Window wd, GC gc,
                     t_car *car, int w0, int h0)
{
    const int w  = 30;
    const int h  = 14;
    const int wr = 10;
    const int hr = 8;
    int       j, w1, h1;
    int       jmax, hmax;

    w1 = w0-w / 2;
    h1 = h0-h / 2;

    /* Carosserie */
    XSetForeground(disp, gc, car->col);
    XFillRectangle(disp, wd, gc, w1, h1, w, h);

    /* Dak */
    XSetForeground(disp, gc, car->roof);
    XFillRectangle(disp, wd, gc, w0-wr/2, h0-hr/2, wr, hr);

    /* Achterlicht */
    if (car->bBrake)
    {
        XSetForeground(disp, gc, YELLOW);
        jmax = 5;
        hmax = 5;
    }
    else
    {
        XSetForeground(disp, gc, LIGHTRED);
        jmax = 3;
        hmax = 3;
    }
    for (j = 1; (j < jmax); j++)
    {
        int w11 = w1-1-j;
        int h11 = h1-1;
        XDrawLine(disp, wd, gc, w11, h11,       w11, h11+hmax);
        XDrawLine(disp, wd, gc, w11, h11+h-hmax, w11, h11+h);
    }

    /* Voorlicht */
    XSetForeground(disp, gc, WHITE);
    for (j = 1; (j < 3); j++)
    {
        int w11 = w1+w+j;
        int h11 = h1-1;
        XDrawLine(disp, wd, gc, w11, h11,    w11, h11+3);
        XDrawLine(disp, wd, gc, w11, h11+h-3, w11, h11+h);
    }
}

static gmx_bool xhwCallBack(struct t_x11 *x11, XEvent *event, Window wd, void *data)
{
    t_xhighway     *xhw;
    t_windata      *win;
    float           sx;
    int             i;
    static     int  nyy = 0;
    static     int *yy;

    xhw = (t_xhighway *)data;
    win = &(xhw->win);

    if (nyy == 0)
    {
        nyy = 2*xhw->ir.nlane+1;
        snew(yy, nyy);
    }
    for (i = 0; (i < nyy); i++)
    {
        yy[i] = ((float) i*win->height)/(nyy-1);
    }

    switch (event->type)
    {
        case Expose:
        {
            if (wd == win->self)
            {
                sx = (float)win->width  / xhw->ir.metres;

                XClearWindow(x11->disp, win->self);
                XSetForeground(x11->disp, x11->gc, WHITE);

                for (i = 2; (i < nyy-1); i += 2)
                {
                    XDrawLine(x11->disp, win->self, x11->gc, 0, yy[i], win->width-1, yy[i]);
                }

                for (i = 0; (i < xhw->ncars); i++)
                {
                    t_car *car = &(xhw->cars[i]);
                    int    w1  = car->x*sx;
                    int    h1  = yy[1+2*(xhw->ir.nlane-1-car->lane)];

                    draw_car(x11->disp, win->self, x11->gc, car, w1, h1);
                }
                if (xhw->bDriving)
                {
                    simulate(x11, xhw, xhw->ncars, xhw->cars, &xhw->ir);
                }
            }
            break;
        }
        case ConfigureNotify:
            if (wd == xhw->main.self)
            {
                xhw->main.width  = event->xconfigure.width;
                xhw->main.height = event->xconfigure.height;
                Configure(xhw);
                for (i = 0; (i < NBUT); i++)
                {
                    XMoveResizeWindow(x11->disp, xhw->but[i].self,
                                      xhw->but[i].x, xhw->but[i].y,
                                      xhw->but[i].width, xhw->but[i].height);
                }
                XMoveResizeWindow(x11->disp, win->self,
                                  win->x, win->y, win->width, win->height);
            }
            else if (wd == win->self)
            {
                win->width  = event->xconfigure.width;
                win->height = event->xconfigure.height;
            }
            break;
        case ButtonPress:
            return TRUE;
        default:
            break;
    }
    return FALSE;
}

static gmx_bool butCallBack(struct t_x11 *x11, XEvent *event, Window wd, void *data)
{
    XSetWindowAttributes attr;
    t_xhighway          *xhw;
    t_windata           *win;
    const char          *label;
    int                  i;

    xhw = (t_xhighway *)data;
    for (i = 0; (i < NBUT); i++)
    {
        if (xhw->but[i].self == wd)
        {
            break;
        }
    }
    if (i == NBUT)
    {
        fprintf(stderr, "Incorrect window: %x in butcallback\n", (unsigned)wd);
        return FALSE;
    }
    win = &(xhw->but[i]);

    switch (event->type)
    {
        case Expose:
            XClearWindow(x11->disp, win->self);
            switch (i)
            {
                case StartStop:
                    label = Driving[xhw->bDriving];
                    break;
                case Fog:
                    label = Fogs[xhw->bFog];
                    break;
                default:
                    label = win->text;
            }
            XSetForeground(x11->disp, x11->gc, WHITE);
            TextInWin(x11, win, label, eXCenter, eYCenter);
            break;

        case ConfigureNotify:
            win->width  = event->xconfigure.width;
            win->height = event->xconfigure.height;
            break;
        case ButtonPress:
            switch (i)
            {
                case Quit:
                    return TRUE;
                case StartStop:
                    xhw->bDriving = 1-xhw->bDriving;
                    ExposeWin(x11->disp, win->self);
                    if (xhw->bDriving)
                    {
                        ExposeWin(x11->disp, xhw->win.self);
                    }
                    break;
                case Fog:
                    xhw->bFog = 1-xhw->bFog;
                    if (xhw->bFog)
                    {
                        attr.background_pixel = DARKGREY;
                    }
                    else
                    {
                        attr.background_pixel = BLACK;
                    }
                    XChangeWindowAttributes(x11->disp, xhw->win.self, CWBackPixel, &attr);
                    /*ExposeWin(x11->disp,win->self);*/
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
    }
    return FALSE;
}

t_xhighway *GetXHW(t_x11 *x11, const char *infile)
{
    t_xhighway *xhw;
    int         i, h, dh, w;
    char        progname[STRLEN];

    snew(xhw, 1);
    xhw->ncars = read_input(x11, infile, &(xhw->cars), &xhw->ir);

    h  = xhw->ir.nlane*40;
    dh = 20;
    w  = 752;
    strncpy(progname, Program(), STRLEN-1);
    InitWin(&xhw->main, 0, 0, w, h+dh+7, 1, progname);
    xhw->main.self = XCreateSimpleWindow(x11->disp, x11->root,
                                         xhw->main.x, xhw->main.y,
                                         xhw->main.width, xhw->main.height,
                                         xhw->main.bwidth, WHITE, BLACK);
    x11->RegisterCallback(x11, xhw->main.self, 0, xhwCallBack, xhw);
    x11->SetInputMask(x11, xhw->main.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask);

    Configure(xhw);

    for (i = 0; (i < NBUT); i++)
    {
        t_windata *wd = &(xhw->but[i]);

        wd->self = XCreateSimpleWindow(x11->disp, xhw->main.self,
                                       wd->x, wd->y,
                                       wd->width, wd->height,
                                       wd->bwidth, WHITE, BLACK);
        x11->RegisterCallback(x11, wd->self, xhw->main.self,
                              butCallBack, xhw);
        x11->SetInputMask(x11, wd->self, ButtonPressMask | ExposureMask |
                          StructureNotifyMask);

    }
    xhw->win.self = XCreateSimpleWindow(x11->disp, xhw->main.self,
                                        xhw->win.x, xhw->win.y,
                                        xhw->win.width, xhw->win.height,
                                        xhw->win.bwidth, WHITE, BLACK);
    x11->RegisterCallback(x11, xhw->win.self, 0, xhwCallBack, xhw);
    x11->SetInputMask(x11, xhw->win.self, ButtonPressMask | ExposureMask |
                      StructureNotifyMask);

    return xhw;
}

int main(int argc, char *argv[])
{
    const char  *desc[] = {
        "[TT]highway[tt] is the GROMCAS highway simulator. It is an X-windows",
        "gadget that shows a (periodic) Autobahn with a user defined",
        "number of cars. Fog can be turned on or off to increase the",
        "number of crashes. Nice for a background CPU-eater. A sample",
        "input file is in [TT]$GMXDATA/top/highway.dat[tt]"
    };
    output_env_t oenv;
    t_x11       *x11;
    t_xhighway  *xhw;
    t_filenm     fnm[] = {
        { efDAT, "-f", "highway", ffREAD }
    };
#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, 0, NFILE, fnm,
                      0, NULL, asize(desc), desc, 0, NULL, &oenv);

    if ((x11 = GetX11(&argc, argv)) == NULL)
    {
        fprintf(stderr, "Can't connect to X Server.\n"
                "Check your DISPLAY environment variable\n");
        exit(1);
    }
    xhw = GetXHW(x11, opt2fn("-f", NFILE, fnm));

    XMapWindow(x11->disp, xhw->main.self);
    XMapSubwindows(x11->disp, xhw->main.self);
    x11->MainLoop(x11);
    x11->CleanUp(x11);

    thanx(stderr);

    return 0;
}
