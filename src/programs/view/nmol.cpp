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

#include "nmol.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "3dview.h"
#include "buttons.h"
#include "manager.h"
#include "xutil.h"

#define MSIZE 4

static bool MWCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_molwin *mw;
    Window    To;
    XEvent    letter;

    mw                          = (t_molwin *)data;
    To                          = mw->wd.Parent;
    letter.type                 = ClientMessage;
    letter.xclient.display      = x11->disp;
    letter.xclient.window       = To;
    letter.xclient.message_type = 0;
    letter.xclient.format       = 32;
    switch (event->type)
    {
        case Expose:
            /* Do Not draw anything, but signal parent instead, he will
             * coordinate drawing.
             */
            letter.xclient.data.l[0] = IDDRAWMOL;
            letter.xclient.data.l[1] = Button1;
            XSendEvent(x11->disp, To, True, 0, &letter);
            break;
        case ButtonPress:
#ifdef DEBUG
            printf("Molwindow: Buttonpress\n");
#endif
            letter.xclient.data.l[0] = IDLABEL;
            letter.xclient.data.l[1] = (long)event->xbutton.button;
            letter.xclient.data.l[2] = event->xbutton.x;
            letter.xclient.data.l[3] = event->xbutton.y;
            XSendEvent(x11->disp, To, True, 0, &letter);
            break;
        case ConfigureNotify:
            mw->wd.width  = event->xconfigure.width;
            mw->wd.height = event->xconfigure.height;
            break;
        default:
            break;
    }
    return false;
}

void set_def (t_molwin *mw, int ePBC, matrix box)
{
    mw->bShowHydrogen = true;
    mw->bond_type     = eBFat;
    mw->ePBC          = ePBC;
    mw->boxtype       = esbRect;
    mw->realbox       = TRICLINIC(box) ? esbTri : esbRect;
}

t_molwin *init_mw(t_x11 *x11, Window Parent,
                  int x, int y, int width, int height,
                  unsigned long fg, unsigned long bg,
                  int ePBC, matrix box)
{
    t_molwin *mw;

    snew(mw, 1);
    set_def(mw, ePBC, box);

    InitWin(&mw->wd, x, y, width, height, 1, "Mol Window");

    mw->wd.Parent = Parent;
    mw->wd.self   = XCreateSimpleWindow(x11->disp, Parent, x, y, width, height, 1, fg, bg);
    x11->RegisterCallback(x11, mw->wd.self, Parent, MWCallBack, mw);
    x11->SetInputMask(x11, mw->wd.self,
                      ExposureMask | StructureNotifyMask |
                      ButtonPressMask);
    return mw;
}

void map_mw(t_x11 *x11, t_molwin *mw)
{
    XMapWindow(x11->disp, mw->wd.self);
}

bool toggle_hydrogen(t_x11 *x11, t_molwin *mw)
{
    mw->bShowHydrogen = !mw->bShowHydrogen;
    ExposeWin(x11->disp, mw->wd.self);

    return mw->bShowHydrogen;
}

void set_bond_type(t_x11 *x11, t_molwin *mw, int bt)
{
    if (bt != mw->bond_type)
    {
        mw->bond_type = bt;
        ExposeWin(x11->disp, mw->wd.self);
    }
}

void set_box_type (t_x11 *x11, t_molwin *mw, int bt)
{
#ifdef DEBUG
    fprintf(stderr, "mw->boxtype = %d, bt = %d\n", mw->boxtype, bt);
#endif
    if (bt != mw->boxtype)
    {
        if ((bt == esbTrunc && mw->realbox == esbTri) || bt == esbTri || bt == esbNone)
        {
            mw->boxtype = bt;
            ExposeWin(x11->disp, mw->wd.self);
        }
        else
        {
            fprintf(stderr, "Can not change rectangular box to truncated octahedron\n");
        }
    }
}

void done_mw(t_x11 *x11, t_molwin *mw)
{
    x11->UnRegisterCallback(x11, mw->wd.self);
    sfree(mw);
}

static void draw_atom(Display *disp, Window w, GC gc,
                      atom_id ai, iv2 vec2[], unsigned long col[], int size[],
                      bool bBall, bool bPlus)
{
    int xi, yi;

    xi = vec2[ai][XX];
    yi = vec2[ai][YY];
    XSetForeground(disp, gc, col[ai]);
    if (bBall)
    {
        XFillCircle(disp, w, gc, xi, yi, size[ai]-1);
        XSetForeground(disp, gc, BLACK);
        XDrawCircle(disp, w, gc, xi, yi, size[ai]);
        /*    XSetForeground(disp,gc,WHITE);
           XFillCircle(disp,w,gc,xi+4,yi-4,4); */
    }
    else if (bPlus)
    {
        XDrawLine(disp, w, gc, xi-MSIZE, yi, xi+MSIZE+1, yi);
        XDrawLine(disp, w, gc, xi, yi-MSIZE, xi, yi+MSIZE+1);
    }
    else
    {
        XDrawLine(disp, w, gc, xi-1, yi, xi+1, yi);
    }

}

/* Global variables */
static rvec gl_fbox, gl_hbox, gl_mhbox;

static void my_init_pbc(matrix box)
{
    int i;

    for (i = 0; (i < DIM); i++)
    {
        gl_fbox[i]  =  box[i][i];
        gl_hbox[i]  =  gl_fbox[i]*0.5;
        gl_mhbox[i] = -gl_hbox[i];
    }
}

static bool local_pbc_dx(rvec x1, rvec x2)
{
    int  i;
    real dx;

    for (i = 0; (i < DIM); i++)
    {
        dx = x1[i]-x2[i];
        if (dx > gl_hbox[i])
        {
            return false;
        }
        else if (dx <= gl_mhbox[i])
        {
            return false;
        }
    }
    return true;
}

static void draw_bond(Display *disp, Window w, GC gc,
                      atom_id ai, atom_id aj, iv2 vec2[],
                      rvec x[], unsigned long col[], int size[], bool bBalls)
{
    unsigned long   ic, jc;
    int             xi, yi, xj, yj;
    int             xm, ym;

    if (bBalls)
    {
        draw_atom(disp, w, gc, ai, vec2, col, size, true, false);
        draw_atom(disp, w, gc, aj, vec2, col, size, true, false);
    }
    else
    {
        if (local_pbc_dx(x[ai], x[aj]))
        {
            ic = col[ai];
            jc = col[aj];
            xi = vec2[ai][XX];
            yi = vec2[ai][YY];
            xj = vec2[aj][XX];
            yj = vec2[aj][YY];

            if (ic != jc)
            {
                xm = (xi+xj) >> 1;
                ym = (yi+yj) >> 1;

                XSetForeground(disp, gc, ic);
                XDrawLine(disp, w, gc, xi, yi, xm, ym);
                XSetForeground(disp, gc, jc);
                XDrawLine(disp, w, gc, xm, ym, xj, yj);
            }
            else
            {
                XSetForeground(disp, gc, ic);
                XDrawLine(disp, w, gc, xi, yi, xj, yj);
            }
        }
    }
}

int compare_obj(const void *a, const void *b)
{
    t_object *oa, *ob;
    real      z;

    oa = (t_object *)a;
    ob = (t_object *)b;

    z = oa->z-ob->z;

    if (z < 0)
    {
        return 1;
    }
    else if (z > 0)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

void z_fill(t_manager *man, real *zz)
{
    t_object *obj;
    int       i;

    for (i = 0, obj = man->obj; (i < man->nobj); i++, obj++)
    {
        switch (obj->eO)
        {
            case eOSingle:
                obj->z = zz[obj->ai];
                break;
            case eOBond:
            case eOHBond:
                obj->z = (zz[obj->ai] + zz[obj->aj]) * 0.5;
                break;
            default:
                break;
        }
    }
}

int filter_vis(t_manager *man)
{
    int          i, nobj, nvis, nhide;
    atom_id      ai;
    bool         bAdd, *bVis;
    t_object    *obj;
    t_object    *newobj;

    nobj = man->nobj;
    snew(newobj, nobj);
    obj   = man->obj;
    bVis  = man->bVis;
    nvis  = 0;
    nhide = nobj-1;
    for (i = 0; (i < nobj); i++, obj++)
    {
        ai   = obj->ai;
        bAdd = bVis[ai];
        if (bAdd)
        {
            if (obj->eO != eOSingle)
            {
                bAdd = bVis[obj->aj];
            }
        }
        if (bAdd)
        {
            newobj[nvis++] = *obj;
        }
        else
        {
            newobj[nhide--] = *obj;
        }
    }
    sfree(man->obj);
    man->obj = newobj;

    return nvis;
}

void draw_objects(Display *disp, Window w, GC gc, int nobj,
                  t_object objs[], iv2 vec2[], rvec x[],
                  unsigned long col[], int size[], bool bShowHydro, int bond_type,
                  bool bPlus)
{
    bool         bBalls;
    int          i;
    t_object    *obj;

    bBalls = false;
    switch (bond_type)
    {
        case eBThin:
            XSetLineAttributes(disp, gc, 1, LineSolid, CapNotLast, JoinRound);
            break;
        case eBFat:
            XSetLineAttributes(disp, gc, 3, LineSolid, CapNotLast, JoinRound);
            break;
        case eBVeryFat:
            XSetLineAttributes(disp, gc, 5, LineSolid, CapNotLast, JoinRound);
            break;
        case eBSpheres:
            bBalls = true;
            bPlus  = false;
            break;
        default:
            gmx_fatal(FARGS, "Invalid bond_type selected: %d\n", bond_type);
            break;
    }
    for (i = 0; (i < nobj); i++)
    {
        obj = &(objs[i]);
        switch (obj->eO)
        {
            case eOSingle:
                draw_atom(disp, w, gc, obj->ai, vec2, col, size, bBalls, bPlus);
                break;
            case eOBond:
                draw_bond(disp, w, gc, obj->ai, obj->aj, vec2, x, col, size, bBalls);
                break;
            case eOHBond:
                if (bShowHydro)
                {
                    draw_bond(disp, w, gc, obj->ai, obj->aj, vec2, x, col, size, bBalls);
                }
                break;
            default:
                break;
        }
    }
    XSetLineAttributes(disp, gc, 1, LineSolid, CapNotLast, JoinRound);
}

static void v4_to_iv2(vec4 x4, iv2 v2, int x0, int y0, real sx, real sy)
{
    real inv_z;

    inv_z  = 1.0/x4[ZZ];
    v2[XX] = x0+sx*x4[XX]*inv_z;
    v2[YY] = y0-sy*x4[YY]*inv_z;
}

static void draw_box(t_x11 *x11, Window w, t_3dview *view, matrix box,
                     int x0, int y0, real sx, real sy, int boxtype)
{
    rvec        rect_tri[8] =  {
        { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
        { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 }
    };
    int         tr_bonds[12][2] = {
        { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
        { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
        { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
    };
    static int *edge = NULL;
    int         i, j, k, i0, i1;
    rvec        corner[NCUCEDGE], box_center;
    vec4        x4;
    iv2         vec2[NCUCEDGE];

    calc_box_center(view->ecenter, box, box_center);
    if (boxtype == esbTrunc)
    {
        calc_compact_unitcell_vertices(view->ecenter, box, corner);
        if (edge == NULL)
        {
            edge = compact_unitcell_edges();
        }

        for (i = 0; (i < NCUCEDGE); i++)
        {
            gmx_mat4_transform_point(view->proj, corner[i], x4);
            v4_to_iv2(x4, vec2[i], x0, y0, sx, sy);
        }
        XSetForeground(x11->disp, x11->gc, YELLOW);
        for (i = 0; i < NCUCEDGE; i++)
        {
            i0 = edge[2*i];
            i1 = edge[2*i+1];
            XDrawLine(x11->disp, w, x11->gc,
                      vec2[i0][XX], vec2[i0][YY], vec2[i1][XX], vec2[i1][YY]);
        }
    }
    else
    {
        if (boxtype == esbRect)
        {
            for (j = 0; (j < DIM); j++)
            {
                box_center[j] -= 0.5*box[j][j];
            }
        }
        else
        {
            for (i = 0; (i < DIM); i++)
            {
                for (j = 0; (j < DIM); j++)
                {
                    box_center[j] -= 0.5*box[i][j];
                }
            }
        }
        for (i = 0; (i < 8); i++)
        {
            clear_rvec(corner[i]);
            for (j = 0; (j < DIM); j++)
            {
                if (boxtype == esbTri)
                {
                    for (k = 0; (k < DIM); k++)
                    {
                        corner[i][k] += rect_tri[i][j]*box[j][k];
                    }
                }
                else
                {
                    corner[i][j] = rect_tri[i][j]*box[j][j];
                }
            }
            rvec_inc(corner[i], box_center);
            gmx_mat4_transform_point(view->proj, corner[i], x4);
            v4_to_iv2(x4, vec2[i], x0, y0, sx, sy);
        }
        if (debug)
        {
            pr_rvecs(debug, 0, "box", box, DIM);
            pr_rvecs(debug, 0, "corner", corner, 8);
        }
        XSetForeground(x11->disp, x11->gc, YELLOW);
        for (i = 0; (i < 12); i++)
        {
            i0 = tr_bonds[i][0];
            i1 = tr_bonds[i][1];
            XDrawLine(x11->disp, w, x11->gc,
                      vec2[i0][XX], vec2[i0][YY], vec2[i1][XX], vec2[i1][YY]);
        }
    }
}

void set_sizes(t_manager *man)
{
    for (int i = 0; i < man->natom; i++)
    {
        if (man->bVis[i])
        {
            man->size[i] = 180*man->vdw[i];
        }
    }
}

void draw_mol(t_x11 *x11, t_manager *man)
{
    static char tstr[2][20];
    static int  ntime = 0;
    t_windata  *win;
    t_3dview   *view;
    t_molwin   *mw;
    int         i, x0, y0, nvis;
    iv2        *vec2;
    real        sx, sy;
    vec4        x4;

    if (!man->status)
    {
        return;
    }

    view = man->view;
    mw   = man->molw;

    win = &(mw->wd);

    vec2 = man->ix;
    x0   = win->width/2;
    y0   = win->height/2;
    sx   = win->width/2*view->sc_x;
    sy   = win->height/2*view->sc_y;

    my_init_pbc(man->box);

    for (i = 0; (i < man->natom); i++)
    {
        if (man->bVis[i])
        {
            gmx_mat4_transform_point(view->proj, man->x[i], x4);
            man->zz[i] = x4[ZZ];
            v4_to_iv2(x4, vec2[i], x0, y0, sx, sy);
        }
    }
    set_sizes(man);

    z_fill (man, man->zz);

    /* Start drawing */
    XClearWindow(x11->disp, win->self);

    /* Draw Time */
    sprintf(tstr[ntime], "Time: %.3f ps", man->time);
    if (strcmp(tstr[ntime], tstr[1-ntime]) != 0)
    {
        set_vbtime(x11, man->vbox, tstr[ntime]);
        ntime = 1-ntime;
    }

    if (mw->boxtype != esbNone)
    {
        draw_box(x11, win->self, view, man->box, x0, y0, sx, sy, mw->boxtype);
    }

    /* Should sort on Z-Coordinates here! */
    nvis = filter_vis(man);
    if (nvis && man->bSort)
    {
        qsort(man->obj, nvis, sizeof(man->obj[0]), compare_obj);
    }

    /* Draw the objects */
    draw_objects(x11->disp, win->self, x11->gc,
                 nvis, man->obj, man->ix, man->x, man->col, man->size,
                 mw->bShowHydrogen, mw->bond_type, man->bPlus);

    /* Draw the labels */
    XSetForeground(x11->disp, x11->gc, WHITE);
    for (i = 0; (i < man->natom); i++)
    {
        if (man->bLabel[i] && man->bVis[i])
        {
            XDrawString(x11->disp, win->self, x11->gc, vec2[i][XX]+2, vec2[i][YY]-2,
                        man->szLab[i], strlen(man->szLab[i]));
        }
    }

    XSetForeground(x11->disp, x11->gc, x11->fg);
}
