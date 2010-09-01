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
#include <typedefs.h>
#include <macros.h>
#include <smalloc.h>
#include <string.h>
#include "nload.h"
#include "buttons.h"

void DrawLoad(t_x11 *x11,t_windata *Win,int nloads,int *loadinfo)
{
  static char *Strings[] = { "Unbalance","Single Node","Your Ad Here ?"};
  int  i,y0,bwidth,boff,bar,bmax,bmin,ym,yh;
  int  *lb;
  real bav,bscale;
  char maxstr[6];

  return;
  
  XClearWindow(x11->disp, Win->self);
  y0=XTextHeight(x11->font)+AIR;
  yh=(Win->height-y0)/2;
  ym=y0+yh;
  XSetForeground(x11->disp,x11->gc,WHITE);
  XDrawLine(x11->disp,Win->self,x11->gc,0,y0,Win->width,y0);
    
  if (nloads >= 2) {
    TextInRect(x11,Win->self,Strings[0],AIR,0,Win->width-2*AIR,y0,
	       eXLeft,eYCenter);
    if (loadinfo[0] == 0) {
      nloads--;
      lb=&loadinfo[1];
    }
    else {
      lb=loadinfo;
      if (loadinfo[nloads-1] == 0) 
	nloads--;
    }
    bwidth = (Win->width) / nloads;
    boff   = (Win->width % nloads)/2;
    bav    = 0.0; 
    
    bmax=bmin=lb[0];
    
    for (i=1; (i<nloads); i++) {
      bmax = max (bmax,lb[i]);
      bmin = min (bmin,lb[i]);
      bav += lb[i];
    }
    bav/=nloads;
    bscale = (yh-2)/max(fabs(bmax-bav),fabs(bav-bmin));
    sprintf(maxstr,"(%d%%)",(int)(100.0*(bmax-bav)/bav));
    TextInRect(x11,Win->self,maxstr,AIR,0,Win->width-2*AIR,y0,
	       eXRight,eYCenter);

    XDrawLine(x11->disp,Win->self,x11->gc,0,ym,Win->width,ym);
    if (bmax-bmin) {
      for(i=0; i<nloads; i++) {
	bar=(lb[i]-bav)*bscale;
	if (bar != 0) {
	  if (bar > 0)
	    XFillRectangle(x11->disp,Win->self,x11->gc,
			   (i*bwidth)+boff+1,ym-bar+1,bwidth-2,bar);
	  else
	    XFillRectangle(x11->disp,Win->self,x11->gc,
			   (i*bwidth)+boff+1,ym,bwidth-2,-bar);
	}
      }
      
    }
  }
  else {
    TextInRect(x11,Win->self,Strings[1],AIR,0,Win->width,y0,eXLeft,eYCenter);
    TextInRect(x11,Win->self,Strings[2],AIR,y0,Win->width,
	       Win->height-y0,eXLeft,eYCenter);
  }
  XSetForeground(x11->disp,x11->gc,x11->fg);
}

static gmx_bool LWCallBack(t_x11 *x11,XEvent *event, Window w, void *data)
{
  t_loadwin *lw;

  lw=(t_loadwin *)data;
  switch(event->type) {
  case Expose:
    DrawLoad(x11,&lw->wd,lw->nnodes,lw->load);
    break;
  default:
    break;
  }
  return FALSE;
}

t_loadwin *init_lw(t_x11 *x11,Window Parent,
		   int x,int y,int width,int height,
		   unsigned long fg,unsigned long bg)
{
  t_loadwin *lw;
  
  snew(lw,1);
  snew(lw->load,MAXNODES);
  lw->nnodes=1;
  InitWin(&lw->wd,x,y,width,height,1,"Load Window");
  lw->wd.self=XCreateSimpleWindow(x11->disp,Parent,x,y,1,1,1,fg,bg);
  x11->RegisterCallback(x11,lw->wd.self,Parent,LWCallBack,lw);
  x11->SetInputMask(x11,lw->wd.self,ExposureMask);

  return lw;
}

void map_lw(t_x11 *x11,t_loadwin *lw)
{
  XMapWindow(x11->disp,lw->wd.self);
}

void set_load(t_x11 *x11,t_loadwin *lw,int nnodes,int load[])
{
  int  i;
  gmx_bool bChange=FALSE;

  lw->nnodes=nnodes;
  for(i=0; (i<nnodes); i++)
    if (lw->load[i] != load[i]) {
      bChange=TRUE;
      lw->load[i]=load[i];
    }
  if (bChange)
    ExposeWin(x11->disp,lw->wd.self);
}

void done_lw(t_x11 *x11,t_loadwin *lw)
{
  x11->UnRegisterCallback(x11,lw->wd.self);
  sfree(lw);
}

