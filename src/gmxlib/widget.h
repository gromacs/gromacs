/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _widget_h
#define _widget_h

static char *SRCID_widget_h = "$Id$";

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
extern XmString char2xms(char *ptr);
extern char     *xms2char(XmString xms);

#endif
