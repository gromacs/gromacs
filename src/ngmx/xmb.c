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
static char *SRCID_xmb_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "typedefs.h"
#include "macros.h"
#include "Xstuff.h"
#include "x11.h"
#include "xdlg.h"
#include "xmb.h"
#include "fatal.h"
#include "gromacs.bm"
#include "stop.bm"
#include "info.bm"
#include "alert.bm"

#define ID_BOX     -3
#define ID_ICON    -2
#define ID_TEXT    -1

static bmchar *icon_bits=NULL;
static int     icon_width=0;
static int     icon_height=0;
static ulong   icon_fg=0;
static ulong   icon_bg=0;

void SetIcon(unsigned char *bits, int w, int h, ulong fg, ulong bg)
{
  icon_bits=(bmchar *)bits;
  icon_width=w;
  icon_height=h;
  icon_fg=fg;
  icon_bg=bg;
}

t_dlg *MessageBox(t_x11 *x11, Window Parent, char *title,
		  int nlines, char *lines[], ulong Flags,
		  DlgCallback *cb, void *data)
{
  t_dlg         *dlg;
  int		width,nicon;
  int           x,y,x0;
  ulong         nFlag;
  ulong         bg;
  
  /* Check flags for inconsistencies */
  if (((Flags & MB_OK) && (Flags & MB_YES)) 	||
      ((Flags & MB_NO) && (Flags & MB_CANCEL))  ||
      (!(Flags & MB_OK) && !(Flags & MB_YES))) {
    fprintf(stderr,"Invalid button selection in MessageBox\n");
    exit(1);
  }
  nicon=0;
  if (Flags & MB_ICONSTOP) nicon++;
  if (Flags & MB_ICONINFORMATION) nicon++;
  if (Flags & MB_ICONEXCLAMATION) nicon++;
  if (Flags & MB_ICONGMX) nicon++;
  if (nicon > 1) 
    fatal_error(0,"More than one (%d) icon selected in MessageBox",nicon);
  /* Input seems ok */
  bg=x11->bg;
  if (nicon > 0) {
    if (Flags & MB_ICONSTOP)
      SetIcon(stop_bits,stop_width,stop_height,RED,bg);
    if (Flags & MB_ICONINFORMATION)
      SetIcon(info_bits,info_width,info_height,BLUE,bg);
    if (Flags & MB_ICONEXCLAMATION)
      SetIcon(alert_bits,alert_width,alert_height,GREEN,bg);
    if (Flags & MB_ICONGMX)
      SetIcon(gromacs_bits,gromacs_width,gromacs_height,BLUE,bg);
  }
  
  dlg=CreateDlg(x11,Parent,title,0,0,0,0,3,x11->fg,bg,cb,data);
  x=2*OFFS_X;
  if (nicon > 0) {
    AddDlgItem(dlg,CreatePixmap
	       (x11,XCreatePixmapFromBitmapData
		(x11->disp,dlg->win.self,icon_bits,icon_width,icon_height,
		 icon_fg,icon_bg,x11->depth),
		ID_ICON,ID_BOX,2*OFFS_X,2*OFFS_Y,icon_width,icon_height,0));
    x+=QueryDlgItemW(dlg,ID_ICON)+2*OFFS_X;
  }
  
  AddDlgItem(dlg,CreateStaticText(x11,nlines,lines,ID_TEXT,ID_BOX,
				  x,2*OFFS_Y,0,0,0));

  y=QueryDlgItemY(dlg,ID_TEXT)+QueryDlgItemH(dlg,ID_TEXT);
  if (nicon > 0) {
    int yi;
    yi=QueryDlgItemY(dlg,ID_ICON)+QueryDlgItemH(dlg,ID_ICON);
    if (yi > y)
      SetDlgItemPos(dlg,ID_TEXT,x,2*OFFS_Y+(yi-y)/2);
    else
      SetDlgItemPos(dlg,ID_ICON,2*OFFS_X,2*OFFS_Y+(y-yi)/2);
    y=max(y,yi);
  }
  x+=QueryDlgItemW(dlg,ID_TEXT)+2*OFFS_X;
  y+=2*OFFS_Y;
  width=(x-8*OFFS_X)/2;
  
  if (((Flags & MB_OKCANCEL) == MB_OKCANCEL) ||
      ((Flags & MB_YESNO) == MB_YESNO))
    x0=2*OFFS_X;
  else
    x0=(x-width)/2;

#define CB(name,butx,id) AddDlgItem(dlg,CreateButton(x11,name,\
						     TRUE,id,ID_BOX,\
						     butx,y,width,0,0))
  if (Flags & MB_OK) CB("OK",x0,MB_OK);
  if (Flags & MB_CANCEL) CB("Cancel",x/2+2*OFFS_X,MB_CANCEL);
  if (Flags & MB_YES) CB("Yes",x0,MB_YES);
  if (Flags & MB_NO) CB("No",x/2+2*OFFS_X,MB_NO);

  SetDlgSize(dlg,x,y+2*OFFS_Y+
	     QueryDlgItemH(dlg,(Flags & MB_OK) ? MB_OK : MB_YES),TRUE);

  if (Flags & MB_SYSTEMMODAL)
    nFlag=DLG_SYSTEMMODAL;
  else if (Flags & MB_APPLMODAL)
    nFlag=DLG_APPLMODAL;
  else
    nFlag=0;
  nFlag=nFlag | DLG_FREEONBUTTON;
  dlg->flags=nFlag;

  if (!(Flags & MB_DONTSHOW))
    ShowDlg(dlg);

  return dlg;
}
