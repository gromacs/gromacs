/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef	_xutil_h
#define	_xutil_h

#ifdef HAVE_IDENT
#ident	"@(#) xutil.h 1.5 11/11/92"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "writeps.h"
#include "Xstuff.h"
#include "x11.h"

#define OFFS_X 	        4
#define OFFS_Y 	        4

typedef struct {
  Window self,Parent;
  ulong  color;
  char   *text;
  bool   bFocus;
  int    x,y,width,height,bwidth;
  Cursor cursor;
} t_windata;

extern int CheckWin(Window win,char *file, int line);

#define CheckWindow(win) CheckWin(win,__FILE__,__LINE__)

extern void LightBorder(Display *disp, Window win, ulong color);

extern void SpecialTextInRect(t_x11 *x11,XFontStruct *font,Drawable win,
			      char *s,int x,int y,int width,int height,
			      eXPos eX,eYPos eY);

extern void TextInRect(t_x11 *x11, Drawable win,
		       char *s, int x, int y, int width, int height,
		       eXPos eX, eYPos eY);

extern void TextInWin(t_x11 *x11, t_windata *win, char *s, eXPos eX, eYPos eY);

extern void InitWin(t_windata *win, int x0,int y0, int w, int h, int bw, char *text);

extern void FreeWin(Display *disp, t_windata *win);

extern void ExposeWin(Display *disp,Window win);

extern void RectWin(Display *disp, GC gc, t_windata *win, ulong color);

extern void XDrawRoundRect(Display *disp, Window win, GC gc,
			   int x, int y, int w, int h);

extern void RoundRectWin(Display *disp, GC gc, t_windata *win,
			 int offsx, int offsy,ulong color);

extern void PushMouse(Display *disp, Window dest, int x, int y);

extern void PopMouse(Display *disp);

extern bool HelpPressed(XEvent *event);

extern bool GrabOK(FILE *out, int err);
/* Return TRUE if grab succeeded, prints a message to out 
 * and returns FALSE otherwise.
 */

#endif	/* _xutil_h */
