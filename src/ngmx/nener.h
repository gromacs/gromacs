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

#ifndef _nener_h
#define _nener_h

#include "x11.h"
#include "xutil.h"
#include "popup.h"

typedef struct {
  t_windata   wd;		/* Window struct 		*/
  t_menu      *selener;		/* The Select energy menu	*/
  int         nre,nwidth;	/* The number of terms 		*/
  int	      nlast;		/* The last frame added		*/	
  int         etype;		/* The term selected		*/
  real        **e;		/* The energy array		*/
} t_enerwin;

extern t_enerwin *init_ew(t_x11 *x11,Window Parent,
			  int x,int y,int width,int height,
			  unsigned long fg,unsigned long bg);

extern void map_ewin(t_x11 *x11,t_enerwin *ew);

extern void add_ener(t_x11 *x11,t_enerwin *ew,t_energy e[]);

extern void rewind_ener(t_x11 *x11,t_enerwin *ew);

extern void done_ew(t_x11 *x11,t_enerwin *ew);

#endif	/* _nener_h */
