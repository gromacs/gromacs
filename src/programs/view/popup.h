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

#ifndef _popup_h
#define _popup_h

#include "x11.h"
#include "xutil.h"

typedef struct {
    Window      send_to;  /* Window to send messages to			*/
    int         nreturn;  /* Value returned when entry is selected  */
    bool        bChecked; /* Indicate whether menu item is check-marked	*/
    const char *str;      /* Text for menu entry				*/
} t_mentry;

typedef struct {
    t_windata  wd;     /* The window struct				*/
    t_mentry  *m;      /* The menu entry				*/
    Window     Parent; /* Parent window id				*/
} t_child;

typedef struct {
    t_windata     wd;       /* The window struct                */
    Window        Parent;   /* The parent of the menu                   */
    int           nitem;    /* The number of menu items			*/
    t_child      *item;     /* Array of child windows                   */
    bool          bGrabbed; /* Did this menu grab the pointer?              */
} t_menu;

extern t_menu *init_menu(t_x11 *x11, Window Parent, unsigned long fg, unsigned long bg,
                         int nent, t_mentry ent[], int ncol);
/* This routine will create a popup menu. It will create a
 * a base window, and child windows for all the items.
 * If ncol != 0 then ncol columns of items will be created;
 * otherwise the routine will try to evenly space the menu, eg. if there
 * are 20 items then the menu will be 2x10 entries, depending on the
 * string lengths.
 * !!!
 * !!! Do not destroy the ent structure while using this menu
 * !!!
 * The routine will create the windows but not map them. That is, this
 * routine can be called once at the beginning of a program. When a menu
 * has to be shown, call show_menu.
 */

extern void show_menu(t_x11 *x11, t_menu *m, int x, int y, bool bGrab);
/* Show the menu in m at (x,y)
 * This will popup the menu, and when a button is released in the
 * menu send a ClientMessage to the Parent window of the menu
 * specifying the selected menu item in xclient.data.l[0].
 * bGrab specifies whether or not to grab the pointer.
 */

extern void hide_menu(t_x11 *x11, t_menu *m);
/* Unmaps the window for m, hides the window */

extern void check_menu_item(t_menu *m, int nreturn, bool bStatus);
/* Set the bChecked field in the menu item with return code
 * nreturn to bStatus. This function must always be called when
 * the bChecked flag has to changed.
 */

extern void done_menu(t_x11 *x11, t_menu *m);
/* This routine destroys the menu m, and unregisters it with x11 */

extern int menu_width(t_menu *m);
/* Return the width of the window */

extern int menu_height(t_menu *m);
/* Return the height of the window */

#endif  /* _popup_h */
