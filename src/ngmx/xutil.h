/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _xutil_h
#define _xutil_h

#include "typedefs.h"
#include "writeps.h"
#include "Xstuff.h"
#include "x11.h"

#define OFFS_X 	        4
#define OFFS_Y 	        4

typedef struct {
  Window self,Parent;
  unsigned long  color;
  char   *text;
  gmx_bool   bFocus;
  int    x,y,width,height,bwidth;
  Cursor cursor;
} t_windata;

extern int CheckWin(Window win,const char *file, int line);

#define CheckWindow(win) CheckWin(win,__FILE__,__LINE__)

extern void LightBorder(Display *disp, Window win, unsigned long color);

extern void SpecialTextInRect(t_x11 *x11,XFontStruct *font,Drawable win,
			      const char *s,int x,int y,int width,int height,
			      eXPos eX,eYPos eY);

extern void TextInRect(t_x11 *x11, Drawable win,
		       const char *s, int x, int y, int width, int height,
		       eXPos eX, eYPos eY);

extern void TextInWin(t_x11 *x11, t_windata *win, const char *s, eXPos eX, eYPos eY);

extern void InitWin(t_windata *win, int x0,int y0, int w, int h, int bw, const char *text);

extern void FreeWin(Display *disp, t_windata *win);

extern void ExposeWin(Display *disp,Window win);

extern void RectWin(Display *disp, GC gc, t_windata *win, unsigned long color);

extern void XDrawRoundRect(Display *disp, Window win, GC gc,
			   int x, int y, int w, int h);

extern void RoundRectWin(Display *disp, GC gc, t_windata *win,
			 int offsx, int offsy,unsigned long color);

extern void PushMouse(Display *disp, Window dest, int x, int y);

extern void PopMouse(Display *disp);

extern gmx_bool HelpPressed(XEvent *event);

extern gmx_bool GrabOK(FILE *out, int err);
/* Return TRUE if grab succeeded, prints a message to out 
 * and returns FALSE otherwise.
 */

#endif	/* _xutil_h */
