/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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

#ifndef _x11_h
#define _x11_h

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"

#include "Xstuff.h"

/* These colours will be mapped to black on a monochrome screen */
extern unsigned long BLACK, BLUE, GREEN, CYAN, RED, BROWN, GREY, DARKGREY;

/* These colours will be mapped to white on a monochrome screen */
extern unsigned long LIGHTBLUE, LIGHTGREY, LIGHTGREEN, LIGHTCYAN,
                     LIGHTRED, VIOLET, YELLOW, WHITE;

#define CBARGS (struct t_x11 *x11, XEvent *event, Window w, void *data)
/* Callback function. Return false to continue, true to exit */

typedef struct t_x11 {
    Display            *disp;
    XFontStruct        *font;
    GC                  gc;
    Window              root;
    char               *dispname;
    FILE               *console;
    int                 screen, depth;
    Colormap            cmap;
    unsigned long       fg, bg;
    char               *title;
    struct t_wlist     *wlist;
    void        (*GetNamedColor)(struct t_x11 *x11, const char *name, unsigned long *col);
    void        (*MainLoop)(struct t_x11 *x11);
    void        (*RegisterCallback)(struct t_x11 *x11, Window w, Window Parent,
                                    bool cb CBARGS, void *data);
    void        (*UnRegisterCallback)(struct t_x11 *x11, Window w);
    void        (*SetInputMask)(struct t_x11 *x11, Window w, unsigned long mask);
    unsigned long       (*GetInputMask)(struct t_x11 *x11, Window w);
    void        (*CleanUp)(struct t_x11 *x11);
    void        (*Flush)(struct t_x11 *x11);
} t_x11;

typedef bool CallBack CBARGS;

typedef struct t_wlist {
    Window                 w;      /* The window itself			*/
    Window                 Parent; /* It's parent window			*/
    CallBack              *cb;     /* Call back function			*/
    unsigned long          mask;   /* Input mask				*/
    void                  *data;   /* User data struct			*/
    struct t_wlist        *next;
} t_wlist;

t_x11 *GetX11(int *argc, char *argv[]);
/* x11 is a struct / function-set that manages a number of windows.
 * more or (presumably) less like Xt does, but since x11 uses only
 * Xlib calls, it is *PORTABLE* software.
 *
 * The x11 struct is in principle Object Oriented, in that the functions
 * are member of the struct. This makes the software a little more
 * managable. Because of portability I decided not to use C++, even
 * though it would be much nicer to work with in the X-Bizz.
 *
 * Here's the description of how to use the x11 struct
 * 1. Call the GetX11 routine, with the argc and argv from your main.
 *    This will sort out the X-arguments on the command line and remove
 *    them from the command line. When the routine returns, only the
 *    application specific arguments should be left. Thi opens the
 *    display, selects a font, creates a Graphics Context and also sets
 *    the colours listed above in the global variables.
 * 2. Call x11->RegisterCallback for each window you want to have
 *    managed by x11. You have to create a Callback routine for your
 *    application that handles *ONE* event at a time. The idea is that
 *    each window has it's own Callback which is not polluted by code
 *    for other windows, but it is of course entirely possible to have
 *    one Callback routine for a number of windows (eg. when you need
 *    to know something about your children).
 * 3. Call x11->SetInputMask. This comes in place of the normal
 *    XSelectInput, because it enables x11 to manually decide which
 *    events are passed to the windows. With the x11->GetInputMask,
 *    x11->SetInputMask combination, a child window can temporarily
 *    disable mouse and keyboard input for it's parent, while allowing
 *    redraw events to pass through for instance. Hereby a simple way
 *    for creating application modal child windows is implemented.
 * 4. Call x11->MainLoop. This will call every callback function as
 *    appropriate. When a window receives a message, that makes it decide
 *    to terminate it should call x11->UnRegisterCallback, in order to
 *    tell the x11 Manager that it does not want to receive any more
 *    events. It is up to the window to destroy itself. The MainLoop
 *    routine exits when there are no more windows to manage, i.e. when
 *    all routines have called UnRegisterCallback, OR when one Callback
 *    routine returns non-zero (true).
 * 5. Call x11->CleanUp. This closes the display, and frees all
 *    memory allocated by x11 before.
 */

extern void GetNamedColor(t_x11 *x11, const char *name, unsigned long *col);

#endif  /* _x11_h */
