/*
 * $Id$
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "widget.h"
#include "smalloc.h"
#include "fatal.h"

typedef struct {
  Widget   w,other,parent;
  bool     bDesc,bPopup;
  XmString desc;
  char     *orignm;
  char     *directory;
  int      ftp;
} t_widget;

static   t_widget *w=NULL;
static   int      nwindex=0,maxwindex=0;
XmString empty_str;

int nwidget(void)
{
  return nwindex;
}

windex add_widget(Widget new_widget,char *desc)
{
  int i;
  
  if (nwindex == maxwindex) {
    maxwindex += 8;
    srenew(w,maxwindex);
    for(i=nwindex; (i<maxwindex); i++) {
      w[i].ftp   = -1;
      w[i].w     = 0;
      w[i].other = 0;
      w[i].parent= 0;
      w[i].bDesc = FALSE;
      w[i].bPopup= FALSE;
      w[i].orignm    = NULL;
      w[i].directory = NULL;
    }
  }
  w[nwindex].w    = new_widget;
  if (desc) {
    w[nwindex].desc  = char2xms(desc);
    w[nwindex].bDesc = TRUE;
  }
  if (debug)
    fprintf(debug,"%s,%d: Successfully added widget %d (%s)\n",__FILE__,__LINE__,
	    nwindex,desc ? desc : "");
  nwindex++;
  
  return nwindex-1;
}

Widget get_parent(windex win)
{
  range_check(win,0,nwindex);
  if (w[win].parent == 0)
    gmx_fatal(FARGS,"No parent widget known for widget %d. I'm an orphan!",win);

  return w[win].parent;
}

void set_parent(windex win,Widget parent)
{
  range_check(win,0,nwindex);
  if (w[win].parent != 0)
    gmx_fatal(FARGS,"Parent widget already set for widget %d",win);
  
  w[win].parent = parent;
}

void set_windex_orignm(windex win,char *orignm)
{
  range_check(win,0,nwindex);

  if (orignm)
    w[win].orignm = strdup(orignm);
  else
    w[win].orignm = NULL;
}

char *get_windex_orignm(windex win)
{
  range_check(win,0,nwindex);
  
  return w[win].orignm;
}

void set_windex_popup(windex win,bool bPopup)
{
  range_check(win,0,nwindex);
  
  w[win].bPopup = bPopup;
}

bool get_windex_popup(windex win)
{
  range_check(win,0,nwindex);
  
  return w[win].bPopup;
}

Widget get_widget(windex win)
{
  range_check(win,0,nwindex);
  
  return w[win].w;
}

windex get_windex(Widget www)
{
  int i;
  
  for(i=0; (i<nwindex); i++)
    if (w[i].w == www)
      return i;
  gmx_fatal(FARGS,"No such widget %x\n",www);
  
  return -1;
}

XmString get_widget_desc(Widget www)
{
  int i;
  
  for(i=0; (i<nwindex); i++)
    if ((w[i].w == www) && w[i].bDesc)
      return w[i].desc;
  
  return empty_str;
}

bool have_windex_desc(windex www)
{
  range_check(www,0,nwindex);
  
  return w[www].bDesc;
}

int get_widget_ftp(Widget www)
{
  int i;
  
  for(i=0; (i<nwindex); i++)
    if (w[i].w == www)
      return w[i].ftp;
  
  return -1;
}

char *get_widget_dir(windex win)
{
  range_check(win,0,nwindex);

  return w[win].directory ? w[win].directory : "";
}

void set_widget_ftp(windex win,int ftp)
{
  range_check(win,0,nwindex);

  w[win].ftp = ftp;
}

void set_widget_dir(Widget www,XmString label)
{
  Arg  args[4];
  int  i,narg;
  char *ptr,*clab,tmp;

  i = get_windex(www);
  
  if (i < nwindex) {
    clab = xms2char(label);
    if (w[i].directory)
      sfree(w[i].directory);
    /* Check for last directory slash */
    if ((ptr = strrchr(clab,'/')) != NULL) {
      /* check whether there is more than the directory */
      if (ptr[1] != '\0') {
	tmp = ptr[1];
	ptr[1] = '\0';
	w[i].directory = strdup(clab);
	ptr[1] = tmp;
      }
      /* Increase the pointer beyond the slash */
      ptr++;
    }
    else {
      w[i].directory = NULL;
      ptr = clab;
    }
    if (strlen(ptr) > 0) {
      narg = 0;
      XtSetArg(args[narg],XmNvalue, ptr); narg++;
      XtSetValues(www,args,narg); 
      /* Set the toggle button if we have an optional one */
      if (w[i].other != 0) {
	narg = 0;
	XtSetArg(args[narg], XmNset, True); narg++;
	XtSetValues(w[i].other,args,narg); 
      }
    }
  }
}

Widget get_widget_other(windex win,bool bFail)
{
  range_check(win,0,nwindex);
  if ((w[win].other == 0) && bFail)
    gmx_fatal(FARGS,"Calling the wrong window: %d, no other widget\n",win);
  
  return w[win].other;
}
   
void set_widget_other(windex win,Widget www)
{
  range_check(win,0,nwindex);
  w[win].other = www;
}
   
XmString char2xms(const char *ptr)
{
  char tmp[4096];

  strncpy(tmp,ptr,4095);

  return XmStringCreate(tmp,XmSTRING_DEFAULT_CHARSET);
}

char *xms2char(XmString xms)
{
  /* This is NOT fool proof */
  char *str;

  /* set str to point to the text */
  XmStringGetLtoR(xms,XmSTRING_DEFAULT_CHARSET,&str);

  return str;
}

