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
static char *SRCID_pulldown_c = "$Id$";

#include <string.h>
#include <smalloc.h>
#include <x11.h>
#include <macros.h>
#include "popup.h"
#include "pulldown.h"

static bool PDCallBack(t_x11 *x11,XEvent *event,Window w,void *data)
{
  t_pulldown *pd;
  int        i,x,x1,y,nsel;

  pd=(t_pulldown *)data;
  y=pd->wd.height;
  switch(event->type) {
  case Expose:
    XSetForeground(x11->disp,x11->gc,x11->fg);
    XDrawLine(x11->disp,w,x11->gc,0,y-1,pd->wd.width,y-1);
    for(i=0; (i<pd->nmenu); i++)
      XDrawString(x11->disp,pd->wd.self,x11->gc,pd->xpos[i],x11->font->ascent,
		  pd->title[i],strlen(pd->title[i]));
    break;
  case ButtonPress:
    if (pd->nsel==-1) {
      x=event->xbutton.x;
      for(nsel=0; (pd->xpos[nsel+1] < x) && (nsel < pd->nmenu-1); nsel++);
      pd->nsel=nsel;
      x1=max(0,min(pd_width(pd)-menu_width(pd->m[nsel]),pd->xpos[nsel]));
      show_menu(x11,pd->m[nsel],x1,y+1,FALSE);
    }
    break;
  case ButtonRelease:
    hide_pd(x11,pd);
    break;
  default:
    break;
  }
  return FALSE;
}

t_pulldown *init_pd(t_x11 *x11,Window Parent,int width,int height,
		    ulong fg,ulong bg,
		    int nmenu,int *nsub,t_mentry *ent[],char **title)
{
  t_pulldown *pd;
  int        i;

  snew(pd,1);
  pd->title=title;
  pd->nmenu=nmenu;
  pd->nsel=-1;
  snew(pd->m,nmenu);
  snew(pd->xpos,nmenu+1);
  pd->xpos[0]=5;
  for(i=1; (i<=nmenu); i++)
    pd->xpos[i]=20+pd->xpos[i-1]+
      XTextWidth(x11->font,title[i-1],strlen(title[i-1]));
  if (pd->xpos[nmenu] > width) 
    printf("Menu too wide\n");

  InitWin(&(pd->wd),0,0,width,XTextHeight(x11->font)+2,0,"PullDown");
  pd->wd.self=XCreateSimpleWindow(x11->disp,Parent,
				  pd->wd.x, pd->wd.y,
				  pd->wd.width,pd->wd.height,
				  pd->wd.bwidth,fg,bg);
  x11->RegisterCallback(x11,pd->wd.self,Parent,PDCallBack,pd);
  x11->SetInputMask(x11,pd->wd.self,ExposureMask | ButtonPressMask | 
		    OwnerGrabButtonMask | ButtonReleaseMask);
  XMapWindow(x11->disp,pd->wd.self);

  for(i=0; (i<nmenu); i++) 
    pd->m[i]=init_menu(x11,Parent,fg,bg,nsub[i],ent[i],1);

  return pd;
}

void hide_pd(t_x11 *x11,t_pulldown *pd)
{
  if (pd->nsel != -1)
    hide_menu(x11,pd->m[pd->nsel]);
  pd->nsel=-1;
}

void check_pd_item(t_pulldown *pd,int nreturn,bool bStatus)
{
  int i;

  for(i=0; (i<pd->nmenu); i++)
    check_menu_item(pd->m[i],nreturn,bStatus);
}

void done_pd(t_x11 *x11,t_pulldown *pd)
{
  int i;

  for(i=0; (i<pd->nmenu); i++)
    done_menu(x11,pd->m[i]);
  x11->UnRegisterCallback(x11,pd->wd.self);
  sfree(pd->m);
  sfree(pd->xpos);
}

int pd_width(t_pulldown *pd)
{
  int i,w;
  
  w=0;
  for(i=0; (i<pd->nmenu); i++)
    w=max(w,menu_width(pd->m[i]));
  w=max(w,pd->xpos[pd->nmenu]);
  return w;
}

int pd_height(t_pulldown *pd)
{
  int i,h;
  
  h=0;
  for(i=0; (i<pd->nmenu); i++)
    h=max(h,menu_height(pd->m[i]));

  return h;
}

