/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013, by the GROMACS development team, led by
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

#ifndef _xdlg_h
#define _xdlg_h

#include "xdlgitem.h"

#define DLG_SHOW         (1<<0)
#define DLG_HIDE         (1<<1)
#define DLG_SHOWANDHIDE  (DLG_SHOW | DLG_HIDE)
#define DLG_SYSTEMMODAL  (1<<2)
#define DLG_APPLMODAL    (1<<3)
#define DLG_HIDEONBUTTON (1<<4)
#define DLG_FREEONBUTTON (1<<5)

enum {
    DLG_SET, DLG_EXIT
};

typedef void DlgCallback (t_x11 *x11, int dlg_mess, int item_id,
                          char *set, void *data);
/* User function that can be called by the dialog box. All setting of
 * check-boxes and radio-buttons etc., is done by the dialog manager,
 * the user can let himself be informed about mouse activity also.
 */

typedef struct {
    t_x11                *x11;        /* All about X                */
    t_windata             win;        /* The position and size of the window    */
    char                 *title;      /* Window name                */
    Window                wDad;       /* The parent window          */
    unsigned int          xmax, ymax; /* Dimensions of parent window        */
    unsigned long         flags;      /* Flags for display          */
    unsigned long         fg, bg;     /* The colours                */
    bool                  bPop;       /* Should we pop the mouse back   */
    bool                  bGrab;      /* Have we grabbed the mouse ?        */
    int                   nitem;      /* The number of items            */
    t_dlgitem           **dlgitem;    /* The array of item pointers         */
    DlgCallback          *cb;         /* User call back function		*/
    void                 *data;       /* User data				*/
} t_dlg;

/*****************************
 *
 * Routine to create the DLG structure, returns NULL on failure
 * cb and data may be NULL.
 *
 ****************************/
t_dlg *CreateDlg(t_x11 *x11, Window Parent, const char *title,
                 int x0, int y0, int w, int h, int bw,
                 DlgCallback *cb, void *data);

/*****************************
 *
 * Routine to add an item to the dialog box
 * The pointer to the item is copied to the dlg struct,
 * the item itself may not be freed until the dlg is done with
 *
 ****************************/
void AddDlgItem(t_dlg *dlg, t_dlgitem *item);

void AddDlgItems(t_dlg *dlg, int nitem, t_dlgitem *item[]);

/*****************************
 *
 * Routines to manipulate items on a dialog box
 * They return true on succes, false otherwise
 * false will mean most of the time, that item id was not found
 *
 ****************************/
bool QueryDlgItemSize(t_dlg *dlg, t_id id, int *w, int *h);

bool QueryDlgItemPos(t_dlg *dlg, t_id id, int *x0, int *y0);

int QueryDlgItemX(t_dlg *dlg, t_id id);

int QueryDlgItemY(t_dlg *dlg, t_id id);

int QueryDlgItemW(t_dlg *dlg, t_id id);

int QueryDlgItemH(t_dlg *dlg, t_id id);

bool SetDlgItemSize(t_dlg *dlg, t_id id, int w, int h);

bool SetDlgItemPos(t_dlg *dlg, t_id id, int x0, int y0);

void SetDlgSize(t_dlg *dlg, int w, int h, bool bAutoPosition);

/*****************************
 *
 * Routines to extract information from the dlg proc
 * after dlg is exec'ed
 *
 ****************************/
bool IsCBChecked(t_dlg *dlg, t_id id);

t_id RBSelected(t_dlg *dlg, int gid);

int  EditTextLen(t_dlg *dlg, t_id id);

char *EditText(t_dlg *dlg, t_id id);

/*****************************
 *
 * Routines to do internal things
 *
 ****************************/
t_dlgitem *FindWin(t_dlg *dlg, Window win);

t_dlgitem *FindItem(t_dlg *dlg, t_id id);

void HelpDlg(t_dlg *dlg);

void HelpNow(t_dlg *dlg, t_dlgitem *dlgitem);

void NoHelp(t_dlg *dlg);

/*****************************
 *
 * Exececute the dialog box procedure
 * Returns when a button is pushed.
 * return value is the ID of the button
 *
 ****************************/
void ShowDlg(t_dlg *dlg);

void HideDlg(t_dlg *dlg);

void FreeDlgItem(t_dlg *dlg, t_id id);

void FreeDlg(t_dlg *dlg);

#endif  /* _xdlg_h */
