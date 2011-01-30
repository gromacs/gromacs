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

#ifndef _logo_h
#define _logo_h

#include <x11.h>
#include <xutil.h>

typedef struct {
  XFontStruct *bigfont;
  XFontStruct *smallfont;
  t_windata   wd;
  gmx_bool bQuitOnClick;
} t_logo;

extern void show_logo(t_x11 *x11,t_logo *logo);

extern void hide_logo(t_x11 *x11,t_logo *logo);

extern t_logo *init_logo(t_x11 *x11,Window parent,gmx_bool bQuitOnClick);

extern void done_logo(t_x11 *x11,t_logo *logo);

#endif	/* _logo_h */
