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

#include "xdlgitem.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "Xstuff.h"

#define BUFSIZE 16

static t_dlgitem *newitem(void)
{
    t_dlgitem *item;

    snew(item, 1);

    return item;
}

/*****************************
 *
 * Window Procedures and helpful functions
 *
 ****************************/
static void ShowCaret(t_x11 *x11, t_dlgitem *dlgitem)
{
    t_edittext *et;

    if (dlgitem->type == edlgET)
    {
        int x, y1, y2;

        et = &(dlgitem->u.edittext);
        x  = XTextWidth(x11->font, dlgitem->win.text, strlen(dlgitem->win.text))+XCARET+
            XTextWidth(x11->font, (char*) &(et->buf[et->strbegin]), et->pos);
        y1 = (dlgitem->win.height-XTextHeight(x11->font))/2;
        y2 = (dlgitem->win.height-y1);
        y1--, y2++;
        XDrawLine(x11->disp, dlgitem->win.self, x11->gc, x-XCARET, y1, x+XCARET, y1);
        XDrawLine(x11->disp, dlgitem->win.self, x11->gc, x, y1, x, y2);
        XDrawLine(x11->disp, dlgitem->win.self, x11->gc, x-XCARET, y2, x+XCARET, y2);
    }
}

static void HideCaret(t_x11 *x11, t_dlgitem *dlgitem)
{
    XSetForeground(x11->disp, x11->gc, x11->bg);
    ShowCaret(x11, dlgitem);
    XSetForeground(x11->disp, x11->gc, x11->fg);
}

static int DefWndProc(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    XComposeStatus status;
    KeySym         keysym;
    char           c[BUFSIZE+1];

#ifdef DEBUG
    printf("DefWndProc\n");
#endif
    switch (event->type)
    {
        case Expose:
        case ButtonPress:
        case KeyPress:
            if (HelpPressed(event))
            {
                return HELPPRESSED;
            }
            else
            {
                XLookupString(&(event->xkey), c, BUFSIZE, &keysym, &status);
                if ((keysym == XK_Return) || (keysym == XK_KP_Enter))
                {
                    return ENTERPRESSED;
                }
            }
            break;
        case EnterNotify:
            dlgitem->win.bFocus = true;
            ShowCaret(x11, dlgitem);
            /*    LightBorder(x11->disp,dlgitem->win.self,x11->fg); */
            break;
        case LeaveNotify:
            dlgitem->win.bFocus = false;
            HideCaret(x11, dlgitem);
            /*    LightBorder(x11->disp,dlgitem->win.self,x11->bg); */
            break;
        default:
            XBell(x11->disp, 50);
    }
    return ITEMOK;
}

static int WndProcBN(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_windata *win;
    int        x, w, th;

    if (dlgitem->type != edlgBN)
    {
        gmx_incons("button processing");
    }
    win = &(dlgitem->win);
    w   = XTextWidth(x11->font, win->text, strlen(win->text));
    x   = (win->width-w)/2;
    th  = XTextHeight(x11->font)+OFFS_Y;
    switch (event->type)
    {
        case Expose:
            RectWin(x11->disp, x11->gc, win, x11->fg);
            TextInRect(x11, win->self, win->text, 0, 0, win->width, th, eXCenter, eYCenter);
            break;
        case ButtonPress:
            return BNPRESSED;
        case EnterNotify:
            XDrawLine(x11->disp, win->self, x11->gc, x-1, th, x+w, th);
            break;
        case LeaveNotify:
            XSetForeground(x11->disp, x11->gc, x11->bg);
            XDrawLine(x11->disp, win->self, x11->gc, x-1, th, x+w, th);
            XSetForeground(x11->disp, x11->gc, x11->fg);
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

static int WndProcRB(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_radiobutton *rb;
    t_windata     *win;
    int            x, y, rad;

    if (dlgitem->type != edlgRB)
    {
        gmx_incons("radiobutton processing");
    }
    rb  = &(dlgitem->u.radiobutton);
    win = &(dlgitem->win);

    rad = win->height/3;
    x   = rad;
    y   = win->height/2;
    switch (event->type)
    {
        case Expose:
            XClearArea(x11->disp, win->self, x-rad, y-rad, x+rad, y+rad, False);
            if (rb->bSelect)
            {
                /* Filled */
                XFillCircle(x11->disp, win->self, x11->gc, x, y, rad);
            }
            XDrawCircle(x11->disp, win->self, x11->gc, x, y, rad);
            x += rad+OFFS_X;
            TextInRect(x11, win->self, win->text, x, 0, win->width-x, win->height,
                       eXLeft, eYCenter);
            break;
        case ButtonPress:
            if (!rb->bSelect)
            {
                return RBPRESSED;
            }
            XBell(x11->disp, 50);
            break;
        case EnterNotify:
        case LeaveNotify:
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

static int WndProcGB(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_windata *win;
    int        x, y;

    if (dlgitem->type != edlgGB)
    {
        gmx_incons("gb processing");
    }
    win = &(dlgitem->win);

    x = XTextWidth(x11->font, win->text, strlen(win->text));
    y = XTextHeight(x11->font);
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, x11->fg);
            XDrawRoundRect(x11->disp, win->self, x11->gc, 0, y/2,
                           win->width-1, win->height-y/2-1);
            XClearArea(x11->disp, win->self, OFFS_X, 0, x+OFFS_X, y, False);
            TextInRect(x11, win->self, win->text, 2*OFFS_X, 0, x, y, eXCenter, eYCenter);
            break;
        case EnterNotify:
        case LeaveNotify:
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

static int WndProcCB(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_checkbox *cb;
    t_windata  *win;
    int         x, y, w, h;

    if (dlgitem->type != edlgCB)
    {
        gmx_incons("check box processing");
    }
    cb  = &(dlgitem->u.checkbox);
    win = &(dlgitem->win);

    x = 0;
    y = win->height/7;
    w = 5*y;
    h = 5*y;
    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, x11->fg);
            XClearArea(x11->disp, win->self, x, y, w, h, False);
            XDrawRectangle(x11->disp, win->self, x11->gc, x, y, w, h);
            if (cb->bChecked)
            {
                XDrawLine(x11->disp, win->self, x11->gc, x, y, x+w, y+h);
                XDrawLine(x11->disp, win->self, x11->gc, x+w, y, x, y+h);
            }
            x = w+OFFS_X;
            TextInRect(x11, win->self, win->text, x, 0, win->width-x, win->height,
                       eXLeft, eYCenter);
            break;
        case ButtonPress:
            cb->bChecked = !cb->bChecked;
            return CBPRESSED;
        case EnterNotify:
        case LeaveNotify:
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

static int WndProcST(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_statictext *st;
    t_windata    *win;
    int           i, dy;

    if (dlgitem->type != edlgST)
    {
        gmx_incons("st processing");
    }
    st  = &(dlgitem->u.statictext);
    win = &(dlgitem->win);

    switch (event->type)
    {
        case Expose:
            dy = XTextHeight(x11->font)+OFFS_Y;
            for (i = 0; (i < st->nlines); i++)
            {
                TextInRect(x11, win->self, st->lines[i],
                           0, OFFS_Y+i*dy, win->width, dy, eXLeft, eYCenter);
            }
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

static bool insert(char *s, char c, int *pos)
{
    int i, sl;

    if (isprint(c))
    {
        sl = strlen(s);
        /* +1 for zero termination */
        for (i = sl+1; (i > *pos); i--)
        {
            s[i+1] = s[i];
        }
        s[*pos] = c;
        (*pos)++;
        return true;
    }
    return false;
}

static bool my_backspace(char *s, int *pos)
{
    int i, sl;

    sl = strlen(s);
    if ((sl > 0) && ((*pos) > 0))
    {
        for (i = *pos-1; (i < sl); i++)
        {
            s[i] = s[i+1];
        }
        (*pos) = std::max(0, (*pos)-1);
        return true;
    }
    return false;
}

static bool my_delete(char *s, int *pos)
{
    int i, sl;

    sl = strlen(s);
    if ((sl > 0) && ((*pos) < sl))
    {
        for (i = *pos; (i < sl); i++)
        {
            s[i] = s[i+1];
        }
        return true;
    }
    return false;
}

static int WndProcET(t_x11 *x11, t_dlgitem *dlgitem, XEvent *event)
{
    t_edittext  *et;
    t_windata   *win;
    KeySym       keysym;
    char         c[BUFSIZE+1], *bp;
    char         scrbuf[STRLEN];
    int          i;
    int          xp, xtitle, ewidth;

    if (dlgitem->type != edlgET)
    {
        gmx_incons("st processing");
    }
    et  = &(dlgitem->u.edittext);
    win = &(dlgitem->win);

    /* Copy string part that is visible into screen buffer */
    for (i = 0; (i < et->buflen); i++)
    {
        scrbuf[i] = et->buf[i+et->strbegin];
    }
    scrbuf[i] = '\0';

    switch (event->type)
    {
        case Expose:
            XSetForeground(x11->disp, x11->gc, x11->fg);
            xtitle = XTextWidth(x11->font, win->text, strlen(win->text));
            ewidth = win->width-xtitle;
            TextInRect(x11, win->self, win->text,
                       0, 0, xtitle-1, win->height, eXLeft, eYCenter);
            XClearArea(x11->disp, win->self, xtitle, 0, ewidth+XCARET, win->height, False);
            TextInRect(x11, win->self, scrbuf,
                       xtitle+XCARET, 0, ewidth, win->height, eXLeft, eYCenter);
#ifdef DEBUG
            printf("Expose\n");
#endif
            if (win->bFocus)
            {
                ShowCaret(x11, dlgitem);
            }
            break;
        case ButtonPress:
            /* Calculate new position for caret */
            et->pos = strlen(et->buf);
            bp      = gmx_strdup(et->buf);
            xp      = event->xbutton.x-XTextWidth(x11->font, win->text, strlen(win->text))-
                XCARET;
            while ((et->pos > 0) && (XTextWidth(x11->font, bp, strlen(bp)) > xp))
            {
                et->pos--;
                bp[et->pos] = '\0';
            }
            sfree(bp);
            et->bChanged = true;
            return ETCHANGED;
        case KeyPress:
            /* Check for HelpKey */
            if (HelpPressed(event))
            {
                return DefWndProc(x11, dlgitem, event);
            }
            XLookupString(&(event->xkey), c, BUFSIZE, &keysym, NULL);
#ifdef DEBUG
            printf("Keysym: %x\n", keysym);
#endif
            switch (keysym)
            {
                case XK_Delete:
                    if (my_delete(et->buf, &(et->pos)))
                    {
                        et->bChanged = true;
                        return ETCHANGED;
                    }
                    else
                    {
                        XBell(x11->disp, 50);
                    }
                    break;
                case XK_BackSpace:
                    if (my_backspace(et->buf, &(et->pos)))
                    {
                        et->bChanged = true;
                        return ETCHANGED;
                    }
                    else
                    {
                        XBell(x11->disp, 50);
                    }
                    break;
                case XK_KP_Enter:
                case XK_Return:
                    return ENTERPRESSED;
                case XK_Home:
                    et->pos      = 0;
                    et->strbegin = 0;
                    et->bChanged = true;
                    return ETCHANGED;
                case XK_End:
                    if (strlen(et->buf) <= (unsigned int)et->buflen)
                    {
                        et->pos = strlen(et->buf);
                    }
                    else
                    {
                        et->pos      = et->buflen;
                        et->strbegin = strlen(et->buf)-et->buflen;
                    }
                    et->bChanged = true;
                    return ETCHANGED;
                case XK_Left:
                    et->pos      = std::max(0, et->pos-1);
                    et->strbegin = std::min(et->strbegin, et->pos);
                    et->bChanged = true;
                    return ETCHANGED;
                case XK_Right:
                    if ((et->pos < et->buflen) &&
                        (et->strbegin+et->buflen > (int)strlen(et->buf)))
                    {
                        et->pos++;
                    }
                    else if ((et->buflen   < (int)strlen(et->buf)) &&
                             (et->strbegin < (int)strlen(et->buf)-et->buflen))
                    {
                        et->strbegin++;
                    }
                    else
                    {
                        break;
                    }
                    et->bChanged = true;
                    return ETCHANGED;
                default:
                    if (keysym < 256)
                    {
                        if (insert(et->buf, c[0], &(et->pos)))
                        {
                            et->bChanged = true;
                            return ETCHANGED;
                        }
                    }
                    XBell(x11->disp, 50);
                    break;
            }
            break;
        case LeaveNotify:
            win->bFocus = false;
            HideCaret(x11, dlgitem);
            if (et->bChanged)
            {
                et->bChanged = false;
            }
            break;
        default:
            return DefWndProc(x11, dlgitem, event);
    }
    return ITEMOK;
}

/*****************************
 *
 * Routines to create dialog items, all items have an id
 * which you can use to extract info. It is possible to have
 * multiple items with the same id but it may then not be possible
 * to extract information.
 * All routines take the position relative to the parent dlg
 * and the size and border width.
 * If the width and height are set to zero initially, they will
 * be calculated and set by the routine. With the dlgitem manipulation
 * routines listed below, the application can then move the items around
 * on the dlg box, and if wished resize them.
 *
 ****************************/
t_dlgitem *CreateButton(t_x11 *x11,
                        const char *szLab, bool bDef, t_id id, t_id groupid,
                        int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;
    char      *lab;

    dlgitem = newitem();
    if (h == 0)
    {
        h = XTextHeight(x11->font)+2*OFFS_Y;
    }
    if (w == 0)
    {
        w = XTextWidth(x11->font, szLab, strlen(szLab))+2*OFFS_X;
    }
    if (bDef)
    {
        snew(lab, strlen(szLab)+7); /* 6 for >> << and 1 for \0 */
        sprintf(lab, ">> %s <<", szLab);
    }
    else
    {
        lab = gmx_strdup(szLab);
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, szLab);
    sfree(lab);
    dlgitem->ID                = id;
    dlgitem->GroupID           = groupid;
    dlgitem->type              = edlgBN;
    dlgitem->u.button.bDefault = bDef;
    dlgitem->WndProc           = WndProcBN;

    return dlgitem;
}

t_dlgitem *CreateRadioButton(t_x11 *x11,
                             const char *szLab, bool bSet, t_id id,
                             t_id groupid,
                             int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;

    dlgitem = newitem();
    if (h == 0)
    {
        h = XTextHeight(x11->font)+OFFS_Y;
    }
    if (w == 0)
    {
        w = XTextWidth(x11->font, szLab, strlen(szLab))+OFFS_X+h;
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, szLab);
    dlgitem->ID                    = id;
    dlgitem->GroupID               = groupid;
    dlgitem->type                  = edlgRB;
    dlgitem->u.radiobutton.bSelect = bSet;
    dlgitem->WndProc               = WndProcRB;

    return dlgitem;
}

t_dlgitem *CreateGroupBox(t_x11 *x11,
                          const char *szLab, t_id id,
                          int nitems, t_id items[],
                          int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;

    dlgitem = newitem();
    if (h == 0)
    {
        h = XTextHeight(x11->font)+OFFS_Y;
    }
    if (w == 0)
    {
        w = XTextWidth(x11->font, szLab, strlen(szLab))+2*OFFS_X;
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, szLab);
    dlgitem->GroupID           = id;
    dlgitem->ID                = id;
    dlgitem->type              = edlgGB;
    dlgitem->u.groupbox.nitems = nitems;
    snew(dlgitem->u.groupbox.item, nitems);
    memcpy((char *)dlgitem->u.groupbox.item, (char *)items,
           nitems*sizeof(items[0]));
    dlgitem->WndProc = WndProcGB;

    return dlgitem;
}

t_dlgitem *CreateCheckBox(t_x11 *x11,
                          const char *szLab, bool bCheckedInitial, t_id id,
                          t_id groupid,
                          int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;

    dlgitem = newitem();
    if (h == 0)
    {
        h = XTextHeight(x11->font)+OFFS_Y;
    }
    if (w == 0)
    {
        w = XTextWidth(x11->font, szLab, strlen(szLab))+OFFS_X+h;
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, szLab);
    dlgitem->ID                  = id;
    dlgitem->GroupID             = groupid;
    dlgitem->type                = edlgCB;
    dlgitem->u.checkbox.bChecked = bCheckedInitial;
    dlgitem->WndProc             = WndProcCB;

    return dlgitem;
}

t_dlgitem *CreatePixmap(Pixmap pm, t_id id,
                        t_id /*groupid*/, int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;

    dlgitem = newitem();
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, NULL);
    dlgitem->ID          = id;
    dlgitem->type        = edlgPM;
    dlgitem->u.pixmap.pm = pm;
    dlgitem->WndProc     = DefWndProc;

    return dlgitem;
}

t_dlgitem *CreateStaticText(t_x11 *x11,
                            int nlines, const char * const *lines, t_id id,
                            t_id groupid,
                            int x0, int y0, int w, int h, int bw)
{
    t_dlgitem *dlgitem;
    int        i;

    dlgitem = newitem();
    if (h == 0)
    {
        h = (XTextHeight(x11->font)+OFFS_Y)*nlines+OFFS_Y;
    }
    if (w == 0)
    {
        for (i = 0; (i < nlines); i++)
        {
            w = std::max(w, XTextWidth(x11->font, lines[i], strlen(lines[i])));
        }
        w += 2*OFFS_X;
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, NULL);
    dlgitem->ID                  = id;
    dlgitem->GroupID             = groupid;
    dlgitem->type                = edlgST;
    dlgitem->u.statictext.nlines = nlines;
    snew(dlgitem->u.statictext.lines, nlines);
    for (i = 0; (i < nlines); i++)
    {
        dlgitem->u.statictext.lines[i] = gmx_strdup(lines[i]);
    }
    dlgitem->WndProc = WndProcST;

    return dlgitem;
}

t_dlgitem *CreateEditText(t_x11 *x11,
                          const char *title,
                          int screenbuf, char *buf, t_id id, t_id groupid,
                          int x0, int y0, int w, int h, int bw)
{
    t_dlgitem  *dlgitem;
    t_edittext *et;

    dlgitem = newitem();
    if (h == 0)
    {
        h = XTextHeight(x11->font)+OFFS_Y;
    }
    if (w == 0)
    {
        char *test;

        snew(test, screenbuf);
        memset(test, 'w', screenbuf);
        w = XTextWidth(x11->font, test, screenbuf)+
            XTextWidth(x11->font, title, strlen(title))+
            2*XCARET+2*OFFS_X;
        sfree(test);
    }
    InitWin(&(dlgitem->win), x0, y0, w, h, bw, title);
    dlgitem->ID      = id;
    dlgitem->GroupID = groupid;
    dlgitem->type    = edlgET;
    et               = &(dlgitem->u.edittext);
    snew(et->buf, STRLEN);
    strcpy(et->buf, buf);
    et->buflen       = screenbuf;
    et->strbegin     = 0;
    et->bChanged     = false;
    dlgitem->WndProc = WndProcET;

    return dlgitem;
}

#define SC(src) (strlen(src) ? gmx_strdup(src) : NULL)

void SetDlgitemOpts(t_dlgitem *dlgitem, bool bUseMon,
                    char *set, char *get, char *help)
{
    dlgitem->bUseMon = bUseMon;
    dlgitem->set     = SC(set);
    dlgitem->get     = SC(get);
    dlgitem->help    = SC(help);
#ifdef DEBUG
    printf("Help is: '%s'\n", dlgitem->help);
#endif
}
