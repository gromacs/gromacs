/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef _pulldown_h
#define _pulldown_h

#include "popup.h"

typedef struct
{
    t_windata    wd;
    int          nmenu;
    int          nsel;
    int*         xpos;
    t_menu**     m;
    const char** title;
} t_pulldown;

extern t_pulldown* init_pd(t_x11*        x11,
                           Window        Parent,
                           int           width,
                           unsigned long fg,
                           unsigned long bg,
                           int           nmenu,
                           int*          nsub,
                           t_mentry*     ent[],
                           const char**  title);
/* nmenu is the number of submenus, title are the titles of
 * the submenus, nsub are the numbers of entries in each submenu
 * ent are the entries in the pulldown menu, analogous to these in the
 * popup menu.
 * The Parent is the parent window, the width is the width of the parent
 * window. The Menu is constructed as a bar on the topside of the window
 * (as usual). It calculates it own height by the font size.
 * !!!
 * !!! Do not destroy the ent structure, or the titles, while using
 * !!! the menu.
 * !!!
 * When the menu is selected, a ClientMessage will be sent to the Parent
 * specifying the selected item in xclient.data.l[0].
 */

extern void hide_pd(t_x11* x11, t_pulldown* pd);
/* Hides any menu that is still on the screen when it shouldn't */

extern void check_pd_item(t_pulldown* pd, int nreturn, bool bStatus);
/* Set the bChecked field in the pd item with return code
 * nreturn to bStatus. This function must always be called when
 * the bChecked flag has to changed.
 */

extern void done_pd(t_x11* x11, t_pulldown* pd);
/* This routine destroys the menu pd, and unregisters it with x11 */

extern int pd_width(t_pulldown* pd);
/* Return the width of the window */

extern int pd_height(t_pulldown* pd);
/* Return the height of the window */

#endif /* _pulldown_h */
