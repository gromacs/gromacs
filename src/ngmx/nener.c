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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <smalloc.h>
#include <macros.h>
#include <names.h>
#include "nener.h"
#include "buttons.h"

static void DrawEGraph(t_x11 *x11,t_enerwin *ew)
{
  t_windata *wd;
  int       i,EHeight,EZero;
  real      epr,scale,MaxE,MinE;
  char      maxstr[80];
  int       y;

  wd=&(ew->wd);
  /* Clear */
  XClearWindow(x11->disp,wd->self);

  /* Calculate boundaries */
  MaxE=MinE=ew->e[ew->etype][0];
  for (i=1; (i<ew->nlast); i++) {
    MaxE=max(ew->e[ew->etype][i],MaxE);
    MinE=min(ew->e[ew->etype][i],MinE);
  }

  /* Print title */
  epr=max(fabs(MaxE),fabs(MinE));
  sprintf(maxstr,"%.0f",epr);
  EHeight=XTextHeight(x11->font)+AIR;
  TextInRect(x11,wd->self,EType[ew->etype],AIR,0,
	     wd->width-2*AIR,EHeight,eXLeft,eYCenter);
  TextInRect(x11,wd->self,maxstr,AIR,0,
	     wd->width-2*AIR,EHeight,eXRight,eYCenter);
  XDrawLine(x11->disp, wd->self,x11->gc,0,EHeight,wd->width,EHeight);
  
  if (ew->nlast==0)
    return;

  if (fabs(MaxE-MinE) < 1e-5)
    return;
  
  EZero=(wd->height-EHeight)/2;
  scale=EZero/(real) epr;
  EZero+=EHeight;
  XDrawLine(x11->disp,wd->self,x11->gc,0,EZero,wd->width,EZero);
  
  for(i=0; (i<ew->nlast); i++) {
    y=ew->e[ew->etype][i]*scale;
    if (y)
      XDrawLine(x11->disp,wd->self,x11->gc,i,EZero,i,EZero-y);
  }
}

static gmx_bool EWCallBack(t_x11 *x11,XEvent *event, Window w, void *data)
{
  t_enerwin *ew;
  int       i,x,y,width;

  return FALSE;
  ew=(t_enerwin *)data;
  switch(event->type) {
  case Expose:
    XSetForeground(x11->disp,x11->gc,WHITE);
    DrawEGraph(x11,ew);
    XSetForeground(x11->disp,x11->gc,x11->fg);
    break;
  case ConfigureNotify:
    ew->wd.x=event->xconfigure.x;
    ew->wd.y=event->xconfigure.y;
    ew->wd.width=event->xconfigure.width;
    ew->wd.height=event->xconfigure.height;
    if (ew->wd.width > ew->nwidth) {
      ew->nwidth=ew->wd.width;
      for (i=0; (i<ew->nre); i++)
	srenew(ew->e[i],ew->nwidth);
    }
    break;
  case ButtonPress:
    x=event->xbutton.x;
    y=ew->wd.y+event->xbutton.y;
    width=menu_width(ew->selener);
    x=min(x+ew->wd.x,ew->wd.x+ew->wd.width-width);
    printf("Showing at %d,%d, width %d\n",x,y,width);
    show_menu(x11,ew->selener,x,y,TRUE);
    break;
  case ClientMessage:
    ew->etype=event->xclient.data.l[0];
    ExposeWin(x11->disp,ew->wd.self);
    /* no break */
  case ButtonRelease:
    hide_menu(x11,ew->selener);
    break;
  default:
    break;
  }
  return FALSE;
}

static void create_selener(t_x11 *x11,t_enerwin *ew,Window Parent)
{
  static t_mentry *se;
  int    i;

  snew(se,ew->nre);
  for(i=0; (i<ew->nre); i++) {
    se[i].send_to=ew->wd.self;
    se[i].nreturn=i;
    se[i].bChecked=FALSE;
    se[i].str=EType[i];
  }
  ew->selener=init_menu(x11,Parent,x11->fg,x11->bg,ew->nre,se,1);
}

t_enerwin *init_ew(t_x11 *x11,Window Parent,
		   int x,int y,int width,int height,
		   unsigned long fg,unsigned long bg)
{
  t_enerwin *ew;
  int       i;
  
  snew(ew,1);
  ew->etype=0;
  ew->nlast=0;
  ew->nwidth=width;
  ew->nre=F_NRE;
  snew(ew->e,ew->nre);
  for(i=0; (i<ew->nre); i++)
    snew(ew->e[i],width);
  InitWin(&ew->wd,x,y,width,height,1,"Ener Window");
  ew->wd.self=XCreateSimpleWindow(x11->disp,Parent,x,y,1,1,1,fg,bg);
  x11->RegisterCallback(x11,ew->wd.self,Parent,EWCallBack,ew);
  x11->SetInputMask(x11,ew->wd.self,ExposureMask | ButtonPressMask |
		    ButtonReleaseMask |  StructureNotifyMask |
		    OwnerGrabButtonMask);
  create_selener(x11,ew,Parent);

  return ew;
}

void map_ewin(t_x11 *x11,t_enerwin *ew)
{
  XMapWindow(x11->disp,ew->wd.self);
}

void add_ener(t_x11 *x11,t_enerwin *ew,t_energy e[])
{
  int i,j,w;
  
  w=ew->nwidth/2;
  if (ew->nlast >= ew->nwidth) {
    for(j=0; (j<ew->nre); j++)
      for(i=0; (i<w); i++)
	ew->e[j][i]=ew->e[j][i+w];
    ew->nlast=w;
  }

  for(j=0; (j<ew->nre); j++) {
    ew->e[j][ew->nlast]=e[j].e;
  }
  ew->nlast++;
  ExposeWin(x11->disp,ew->wd.self);
}

void rewind_ener(t_x11 *x11,t_enerwin *ew)
{
  ew->nlast=0;
  ExposeWin(x11->disp,ew->wd.self);
}

void done_ew(t_x11 *x11,t_enerwin *ew)
{
  done_menu(x11,ew->selener);
  x11->UnRegisterCallback(x11,ew->wd.self);
  sfree(ew);
}


