/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _x11_h
#define _x11_h

static char *SRCID_x11_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) x11.h 1.6 12/16/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "Xstuff.h"

/* These colours will be mapped to black on a monochrome screen */
extern ulong BLACK,BLUE,GREEN,CYAN,RED,BROWN,GREY,DARKGREY;

/* These colours will be mapped to white on a monochrome screen */
extern ulong LIGHTBLUE,LIGHTGREY,LIGHTGREEN,LIGHTCYAN,
             LIGHTRED,VIOLET,YELLOW,WHITE;

typedef enum { ecbOK } ecbReturn;

#define CBARGS (struct t_x11 *x11,XEvent *event, Window w, void *data)
/* Callback function. Return FALSE to continue, TRUE to exit */

typedef struct t_x11 {
  Display     *disp;
  XFontStruct *font;
  GC          gc;
  Window      root;
  char        *dispname;
  FILE        *console;
  int         screen,depth;
  Colormap    cmap;
  ulong       fg,bg;
  char        *title;
  struct t_wlist *wlist;
  void        (*GetNamedColor)(struct t_x11 *x11,char *name,ulong *col);
  void        (*MainLoop)(struct t_x11 *x11);
  void        (*RegisterCallback)(struct t_x11 *x11,Window w,Window Parent,
				  bool cb CBARGS, void *data);
  void        (*UnRegisterCallback)(struct t_x11 *x11, Window w);
  void        (*SetInputMask)(struct t_x11 *x11, Window w, ulong mask);
  ulong       (*GetInputMask)(struct t_x11 *x11, Window w);
  void        (*CleanUp)(struct t_x11 *x11);
  void        (*Flush)(struct t_x11 *x11);
} t_x11;

typedef bool CallBack CBARGS;

typedef struct t_wlist {
  Window         w;		/* The window itself			*/
  Window         Parent;	/* It's parent window			*/
  CallBack       *cb;		/* Call back function			*/
  ulong          mask;		/* Input mask				*/
  void           *data;		/* User data struct			*/
  struct t_wlist *next;
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
 *    routine returns non-zero (TRUE).
 * 5. Call x11->CleanUp. This closes the display, and frees all 
 *    memory allocated by x11 before.
 */

extern void GetNamedColor(t_x11 *x11,char *name,ulong *col);

#endif	/* _x11_h */
