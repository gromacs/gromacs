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

#include "popup.h"

#include <math.h>
#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/smalloc.h"

#include "x11.h"
#include "xutil.h"

bool ChildCallBack(t_x11 *x11, XEvent *event, Window w, void *data)
{
    t_child   *child;
    t_mentry  *m;
    t_windata *wd;
    XEvent     letter;

    child = (t_child *)data;
    m     = child->m;
    wd    = &(child->wd);
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, x11->fg);
            TextInRect(x11, w, m->str, 16, 0, wd->width-16-2, wd->height-2,
                       eXLeft, eYCenter);
            if (m->bChecked)
            {
                int y = x11->font->ascent;
                XDrawLine(x11->disp, w, x11->gc, 2, (y*2)/3, 6, y);
                XDrawLine(x11->disp, w, x11->gc, 3, (y*2)/3, 7, y);
                XDrawLine(x11->disp, w, x11->gc, 7, y, 12, 2);
            }
            break;
        case EnterNotify:
            LightBorder(x11->disp, w, x11->fg);
            break;
        case LeaveNotify:
            LightBorder(x11->disp, w, x11->bg);
            break;
        case ButtonRelease:
            letter.type                 = ClientMessage;
            letter.xclient.display      = x11->disp;
            letter.xclient.window       = m->send_to ? m->send_to : child->Parent;
            letter.xclient.message_type = 0;
            letter.xclient.format       = 32;
            letter.xclient.data.l[0]    = m->nreturn;
            XSendEvent(x11->disp, letter.xclient.window, True, 0, &letter);
            break;
        default:
            break;
    }
    return false;
}

bool MenuCallBack(t_x11 *x11, XEvent *event, Window /*w*/, void *data)
{
    t_menu *m;

    m = (t_menu *)data;
    switch (event->type)
    {
        case Expose:
            /* Nothing to be done */
            if (m->bGrabbed)
            {
                m->bGrabbed =
                    GrabOK(stderr, XGrabPointer(x11->disp, m->wd.self, True,
                                                ButtonReleaseMask, GrabModeAsync,
                                                GrabModeAsync, m->wd.self, None, CurrentTime));
            }
            break;
        case ButtonRelease:
            hide_menu(x11, m);
            break;
        case ClientMessage:
            event->xclient.window = m->Parent;
            XSendEvent(x11->disp, m->Parent, True, 0, event);
            break;
        default:
            break;
    }
    return false;
}

t_menu *init_menu(t_x11 *x11, Window Parent, unsigned long fg, unsigned long bg,
                  int nent, t_mentry ent[], int ncol)
{
    int        i, mlen, mht, area, ht;
    int        j, k, l;
    int        frows, fcol;
    t_menu    *m;
    t_child   *kid;
    t_windata *w;

    snew(m, 1);
    m->nitem  = nent;
    m->Parent = Parent;

    /* Calculate dimensions of the menu */
    mlen = 0;
    for (i = 0; (i < nent); i++)
    {
        mlen = std::max(mlen, XTextWidth(x11->font, ent[i].str, strlen(ent[i].str)));
    }
    mht = XTextHeight(x11->font);
    /* Now we have the biggest single box, add a border of 2 pixels */
    mlen += 20; /* We need extra space at the left for checkmarks */
    mht  += 4;
    /* Calculate the area of the menu */
    area = mlen*mht;
    ht   = sqrt(area);
    /* No the number of rows per column, only beyond 8 rows */
    if (ncol == 0)
    {
        if (nent > 8)
        {
            frows = (1+ht/mht);
        }
        else
        {
            frows = nent;
        }
        fcol = nent/frows;
    }
    else
    {
        fcol  = ncol;
        frows = nent/ncol;
        if (nent % ncol)
        {
            frows++;
        }
    }
    InitWin(&(m->wd), 10, 10, fcol*mlen, frows*mht, 1, "Menu");
    snew(m->item, nent);
    m->wd.self = XCreateSimpleWindow(x11->disp, Parent,
                                     m->wd.x, m->wd.y,
                                     m->wd.width, m->wd.height,
                                     m->wd.bwidth, fg, bg);
    x11->RegisterCallback(x11, m->wd.self, Parent, MenuCallBack, m);
    x11->SetInputMask(x11, m->wd.self, ExposureMask |
                      OwnerGrabButtonMask | ButtonReleaseMask);

    for (j = l = 0; (j < fcol); j++)
    {
        for (k = 0; (k < frows) && (l < nent); k++, l++)
        {
            kid         = &(m->item[l]);
            kid->m      = &(ent[l]);
            kid->Parent = Parent;
            w           = &(kid->wd);
            InitWin(w, j*mlen, k*mht, mlen-2, mht-2, 1, NULL);
            w->self = XCreateSimpleWindow(x11->disp, m->wd.self,
                                          w->x, w->y, w->width, w->height,
                                          w->bwidth, bg, bg);
            x11->RegisterCallback(x11, w->self, m->wd.self,
                                  ChildCallBack, kid);
            x11->SetInputMask(x11, w->self,
                              ButtonPressMask | ButtonReleaseMask |
                              OwnerGrabButtonMask | ExposureMask |
                              EnterWindowMask | LeaveWindowMask);
        }
    }

    return m;
}

void show_menu(t_x11 *x11, t_menu *m, int x, int y, bool bGrab)
{
    XMoveWindow(x11->disp, m->wd.self, x, y);
    m->bGrabbed = bGrab;
    XMapWindow(x11->disp, m->wd.self);
    XMapSubwindows(x11->disp, m->wd.self);
}

void hide_menu(t_x11 *x11, t_menu *m)
{
    if (m->bGrabbed)
    {
        XUngrabPointer(x11->disp, CurrentTime);
    }
    XUnmapWindow(x11->disp, m->wd.self);
}

void check_menu_item(t_menu *m, int nreturn, bool bStatus)
{
    int i;

    for (i = 0; (i < m->nitem); i++)
    {
        if (m->item[i].m->nreturn == nreturn)
        {
            m->item[i].m->bChecked = bStatus;
        }
    }
}

void done_menu(t_x11 *x11, t_menu *m)
{
    int i;

    for (i = 0; (i < m->nitem); i++)
    {
        x11->UnRegisterCallback(x11, m->item[i].wd.self);
    }
    sfree(m->item);
    x11->UnRegisterCallback(x11, m->wd.self);
    sfree(m);
}

int menu_width(t_menu *m)
{
    return m->wd.width;
}

int menu_height(t_menu *m)
{
    return m->wd.height;
}
