/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _xmb_h
#define _xmb_h

static char *SRCID_xmb_h = "$Id$";
#ifdef HAVE_IDENT
#ident	"@(#) xmb.h 1.3 9/29/92"
#endif /* HAVE_IDENT */

#include <x11.h>

#define MB_OK              1
#define MB_CANCEL          (1<<1)
#define MB_OKCANCEL        (MB_OK | MB_CANCEL)
#define MB_YES             (1<<2)
#define MB_NO              (1<<3)
#define MB_YESNO           (MB_YES | MB_NO)
#define MB_ICONSTOP        (1<<16)
#define MB_ICONINFORMATION (1<<17)
#define MB_ICONEXCLAMATION (1<<18)
#define MB_ICONGMX         (1<<19)
#define MB_SYSTEMMODAL     (1<<20)
#define MB_APPLMODAL       (1<<21)
#define MB_DONTSHOW        (1<<22)

t_dlg *MessageBox(t_x11 *x11, Window Parent, char *title,
		  int nlines, char *lines[], unsigned long Flags,
		  DlgCallback *cb, void *data);

#endif	/* _xmb_h */
