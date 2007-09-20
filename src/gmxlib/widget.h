/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _widget_h
#define _widget_h

#include "typedefs.h"
#include <Xm/Xm.h>

typedef int windex;

extern int      nwidget(void);
extern windex   add_widget(Widget new_widget,char *desc);
extern Widget   get_widget(windex win);
extern windex   get_windex(Widget www);
extern Widget   get_parent(windex win);
extern void     set_parent(windex win,Widget parent);
extern void     set_windex_orignm(windex win,char *orignm);
extern char     *get_windex_orignm(windex win);
extern XmString get_widget_desc(Widget www);
extern bool     have_windex_desc(windex www);
extern bool     get_windex_popup(windex win);
extern void     set_windex_popup(windex win,bool bPopup);

extern int      get_widget_ftp(Widget www);
extern void     set_widget_ftp(windex win,int ftp);

extern char     *get_widget_dir(windex win);
extern void     set_widget_dir(Widget www,XmString label);

extern Widget   get_widget_other(windex win,bool bFail);
extern void     set_widget_other(windex win,Widget www);

extern void     mk_desc_callbacks(void);
extern XmString char2xms(const char *ptr);
extern char     *xms2char(XmString xms);

#endif
