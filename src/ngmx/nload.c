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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_nload_c = "$Id$";

#include <math.h>
#include <typedefs.h>
#include <macros.h>
#include <smalloc.h>
#include <string.h>
#include "nload.h"
#include "buttons.h"

void DrawLoad(t_x11 *x11,t_windata *Win,int nloads,int *loadinfo)
{
  static char *Strings[] = { "Unbalance","Single Processor","Your Ad Here ?"};
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

static bool LWCallBack(t_x11 *x11,XEvent *event, Window w, void *data)
{
  t_loadwin *lw;

  lw=(t_loadwin *)data;
  switch(event->type) {
  case Expose:
    DrawLoad(x11,&lw->wd,lw->nprocs,lw->load);
    break;
  default:
    break;
  }
  return FALSE;
}

t_loadwin *init_lw(t_x11 *x11,Window Parent,
		   int x,int y,int width,int height,
		   ulong fg,ulong bg)
{
  t_loadwin *lw;
  
  snew(lw,1);
  snew(lw->load,MAXPROC);
  lw->nprocs=1;
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

void set_load(t_x11 *x11,t_loadwin *lw,int nprocs,int load[])
{
  int  i;
  bool bChange=FALSE;

  lw->nprocs=nprocs;
  for(i=0; (i<nprocs); i++)
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

