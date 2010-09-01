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

#include <string.h>

#include "gmx_fatal.h"
#include "string2.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "xdlghi.h"
#include "fgrid.h"

t_dlgitem **CreateRadioButtonGroup(t_x11 *x11, char *szTitle, 
				   t_id GroupID, int nrb, t_id rb[],
				   int nSelect,
				   char *szRB[], int x0,int y0)
     /* This routine creates a radio button group at the
      * specified position. The return values is a pointer to an
      * array of dlgitems, the array has length (nrb+1) with the +1
      * because of the groupbox.
      * nSelect is the ordinal of the selected button.
      */
{
  t_dlgitem **dlgitem;
  int x,y,w;
  int i;
  
  snew(dlgitem,nrb+1);
  dlgitem[0]=CreateGroupBox(x11,szTitle,GroupID,nrb,rb,x0,y0,0,0,0);
  x=x0+2*OFFS_X;
  y=dlgitem[0]->win.y+dlgitem[0]->win.height;
  w=0;
  for (i=0; (i<nrb); i++) {
    dlgitem[i+1]=CreateRadioButton(x11,szRB[i],(i==nSelect),
				   rb[i],GroupID,x,y,0,0,0);
    y+=dlgitem[i+1]->win.height+OFFS_Y;
    w=max(w,dlgitem[i+1]->win.width);
  }
  for (i=0; (i<nrb); i++)
    dlgitem[i+1]->win.width=w;
  dlgitem[0]->win.width=w+4*OFFS_X;
  dlgitem[0]->win.height=y-y0;
  
  return dlgitem;
}

t_dlgitem **CreateDlgitemGroup(t_x11 *x11, const char *szTitle, 
			       t_id GroupID, int x0, int y0,
			       int nitem, ...)
     /* This routine creates a dlgitem group at the
      * specified position. The return values is a pointer to an
      * array of dlgitems, the array has length (nitem+1) with the +1
      * because of the groupbox.
      */
{
  va_list   ap;
  
  t_dlgitem **dlgitem;
  t_id      *ids;
  edlgitem  edlg;
  char      *name;
  gmx_bool      bBool;
  Pixmap    pm;
  int       nlines,buflen;
  char      *buf,**lines;
  int       x,y,w,i;
  
  va_start(ap,nitem);
  
  snew(dlgitem,nitem+1);
  snew(ids,nitem);
  x=x0+2*OFFS_X;
  dlgitem[0]=CreateGroupBox(x11,szTitle,GroupID,nitem,ids,x0,y0,0,0,0);
  y=dlgitem[0]->win.y+dlgitem[0]->win.height;
  w=0;
  for (i=0; (i<nitem); i++) {
    edlg=(edlgitem)va_arg(ap,int);
    ids[i]=va_arg(ap,int);
    switch (edlg) {
    case edlgBN:
      name=va_arg(ap,char *);
      bBool=va_arg(ap,int);
      dlgitem[i+1]=CreateButton(x11,name,bBool,ids[i],GroupID,x,y,0,0,0);
      break;
    case edlgRB:
      name=va_arg(ap,char *);
      bBool=va_arg(ap,int);
      dlgitem[i+1]=CreateRadioButton(x11,name,bBool,ids[i],GroupID,x,y,0,0,0);
      break;
    case edlgCB:
      name=va_arg(ap,char *);
      bBool=va_arg(ap,int);
      dlgitem[i+1]=CreateCheckBox(x11,name,bBool,ids[i],GroupID,x,y,0,0,0);
      break;
    case edlgPM:
      pm=va_arg(ap,Pixmap);
      dlgitem[i+1]=CreatePixmap(x11,pm,ids[i],GroupID,x,y,0,0,0);
      break;
    case edlgST:
      nlines=va_arg(ap,int);
      lines=va_arg(ap,char **);
      dlgitem[i+1]=CreateStaticText(x11,nlines,lines,ids[i],GroupID,
				    x,y,0,0,0);
      break;
    case edlgET:
      name=va_arg(ap,char *);
      buflen=va_arg(ap,int);
      buf=va_arg(ap,char *);
      dlgitem[i+1]=CreateEditText(x11,name,buflen,buf,ids[i],
				  GroupID,x,y,0,0,0);
      break;
    case edlgGB:
    default:
      gmx_fatal(FARGS,"Invalid dlgitem type: %d\n",edlg);
    }
    y+=dlgitem[i+1]->win.height+OFFS_Y;
    w=max(w,dlgitem[i+1]->win.width);
  }
  va_end(ap);
  sfree(dlgitem[0]->u.groupbox.item);
  sfree(dlgitem[0]->win.text);
  dlgitem[0]=CreateGroupBox(x11,szTitle,GroupID,nitem,ids,x0,y0,0,0,0);
  for (i=0; (i<nitem); i++)
    dlgitem[i+1]->win.width=w;
  dlgitem[0]->win.width=w+4*OFFS_X;
  dlgitem[0]->win.height=y-y0;
  return dlgitem;
}

static void AddDlgItemGroups(t_dlg *dlg, int gridx, int gridy, 
			     t_dlgitemlist **grid, gmx_bool bAutoPosition)
{
  t_dlgitemlist *item;
  int x1,y1,w1,h1;
  int x,y,dw,dh;
  float w,h;
  
  w=h=0;
  for(x=0; (x<gridx); x++)
    for(y=0; (y<gridy); y++) {
      item=&(grid[x][y]);
      if (item->nitem) {
	if (!item->list) {
	  printf("Error: empty list with non-empty nitem (%d)\n",item->nitem);
	  printf("       at grid point: %d,%d\n",x,y);
	  printf("       with size: %dx%d\n",item->w,item->h);
	  exit(1);
	}
	else {
	  AddDlgItems(dlg,item->nitem,item->list);
	  dw=item->w;
	  dh=item->h;
	  w=max(w,((float) QueryDlgItemW(dlg,item->list[0]->ID))/dw);
	  h=max(h,((float) QueryDlgItemH(dlg,item->list[0]->ID))/dh);
	}
      }
    }
  w1=gridx*w;
  h1=gridy*h;
  SetDlgSize(dlg,w1,h1,bAutoPosition);
#ifdef DEBUG
  printf("Dimensions of grid cell: %8.3f x %8.3f\n",w,h);
  printf("Dimensions of window:    %d x %d\n",w1,h1);
#endif
  
  for(x=0; (x<gridx); x++) 
    for(y=0; (y<gridy); y++) {
      item=&(grid[x][y]);
      if (item->nitem) {
	x1=x*w;
	y1=y*h;
	w1=item->w*w;
	h1=item->h*h;
#ifdef DEBUG
	printf("New size: %d x %d at %d, %d\n",w1,h1,x1,y1);
#endif
	SetDlgItemSize(dlg,item->list[0]->ID,w1,h1);
	SetDlgItemPos(dlg,item->list[0]->ID,x1,y1);
      }
    }
}

static t_dlgitemlist **NewDlgitemList(int w, int h)
{
  int i,j;
  t_dlgitemlist **grid;
  
  snew(grid,w);
  for(i=0; (i<w); i++) {
    snew(grid[i],h);
    for (j=0; (j<h); j++) {
      grid[i][j].nitem=0;
      grid[i][j].list=NULL;
    }
  }
  return grid;
}

static void AddListItem(t_dlgitemlist *list, t_dlgitem *item)
{
  srenew(list->list,++list->nitem);
  list->list[list->nitem-1]=item;
}

static void AddListFItem(t_x11 *x11, t_dlgitemlist *list, 
			 t_fitem *fitem, t_id GroupID, t_id *ID,
			 int x, int *y, int *w,gmx_bool bUseMon)
{
  int i,iSel,slen;
  char buf[STRLEN];
  
  switch(fitem->edlg) {
  case edlgBN:
    AddListItem
      (list,CreateButton(x11,fitem->name[0],fitem->bDef,(*ID)++,GroupID,
			 x,(*y),0,0,0));
    break;
  case edlgRB:
    strcpy(buf,fitem->def);
    iSel=-1;
    for (i=0; (i<fitem->nname); i++) {
      char buf2[100];

      strcpy(buf2,fitem->name[i]);
      buf2[strlen(buf)]='\0'; /* truncate itemname */
      if (gmx_strcasecmp(buf2,buf)==0)
	iSel=i;
    }

    for(i=0; (i<fitem->nname); i++) {
      AddListItem(list,
		  CreateRadioButton(x11,fitem->name[i],(iSel==i),
				    (*ID)++,GroupID,x,(*y),0,0,0));
      (*y)+=list->list[list->nitem-1]->win.height+OFFS_Y;
      (*w)=max((*w),list->list[list->nitem-1]->win.width);
      SetDlgitemOpts(list->list[list->nitem-1],bUseMon,
		     fitem->set,fitem->get,fitem->help);
    }
    break;
  case edlgCB: {
    gmx_bool bCheck;

    bCheck=gmx_strcasecmp(fitem->def,"TRUE")==0;
    AddListItem(list,CreateCheckBox(x11,fitem->name[0],bCheck,
				    (*ID)++,GroupID,x,(*y),0,0,0));
    break;
  }
  case edlgST:
    AddListItem(list,
		CreateStaticText(x11,fitem->nname,fitem->name,(*ID)++,
				 GroupID,x,(*y),0,0,0));
    break;
  case edlgET:
    slen=strlen(fitem->name[0])+strlen(fitem->def);
    AddListItem(list,CreateEditText(x11,fitem->name[0],slen,fitem->def,
				    (*ID)++,GroupID,x,(*y),0,0,0));
    break;
  case edlgPM:
  case edlgGB:
  default:
    gmx_fatal(FARGS,"Invalid list->list type: %d\n",fitem->edlg);
  }
  SetDlgitemOpts(list->list[list->nitem-1],bUseMon,
		 fitem->set,fitem->get,fitem->help);
  
  if (fitem->edlg != edlgRB) {
    (*y)+=list->list[list->nitem-1]->win.height+OFFS_Y;
    (*w)=max((*w),list->list[list->nitem-1]->win.width);
  }
}

static void AddListFGroup(t_x11 *x11, t_dlgitemlist **grid,
			  t_fgroup *fgroup, t_id *ID,gmx_bool bUseMon)
{
  int i;
  t_id GroupID,*ids;
  t_dlgitemlist *item;
  int x,y,w;
  
  GroupID=(*ID)++;
  item=&(grid[fgroup->x][fgroup->y]);
  AddListItem(item,CreateGroupBox(x11,fgroup->name,GroupID,
				  0,NULL,0,0,0,0,0));
  x=2*OFFS_X;
  y=item->list[0]->win.y+item->list[0]->win.height;
  w=0;
  for(i=0; (i<fgroup->nfitem); i++)
    AddListFItem(x11,item,fgroup->fitem[i],GroupID,ID,x,&y,&w,bUseMon);
  
  w=max(w,item->list[0]->win.width+4*OFFS_X);
  sfree(item->list[0]->u.groupbox.item);
  sfree(item->list[0]->win.text);
  snew(ids,item->nitem);
  for(i=0; (i<item->nitem-1); i++)
    ids[i]=GroupID+i+1;
  item->list[0]=
    CreateGroupBox(x11,fgroup->name,GroupID,item->nitem-1,ids,
		   2*OFFS_X,2*OFFS_Y,w+2*OFFS_X,y,0);
  sfree(ids);
  item->w=fgroup->w;
  item->h=fgroup->h;
}

static void AddListFSimple(t_x11 *x11, t_dlgitemlist **grid,
			   t_fsimple *fsimple, t_id *ID,gmx_bool bUseMon)
{
  t_dlgitemlist *item;
  int x,y,w;
  
  item=&(grid[fsimple->x][fsimple->y]);
  x=0;
  y=0;
  w=0;
  AddListFItem(x11,item,fsimple->fitem,*ID,ID,x,&y,&w,bUseMon);
  item->w=fsimple->w;
  item->h=fsimple->h;
}

t_dlg *ReadDlg(t_x11 *x11,Window Parent, const char *title,
	       unsigned long fg, unsigned long bg, const char *infile, 
	       int x0, int y0, gmx_bool bAutoPosition,gmx_bool bUseMon,
	       DlgCallback *cb,void *data)
{
  t_fgrid       *fgrid;
  t_dlgitemlist **grid;
  t_dlg         *dlg;
  int           i;
  t_id          ID;
  
  fgrid=FGridFromFile(infile);
  dlg=CreateDlg(x11,Parent,title,x0,y0,0,0,0,fg,bg,cb,data);
  grid=NewDlgitemList(fgrid->w,fgrid->h);
  ID=0;
  
  for(i=0; (i<fgrid->nfgroup); i++)
    AddListFGroup(x11,grid,fgrid->fgroup[i],&ID,bUseMon);
  for(i=0; (i<fgrid->nfsimple); i++)
    AddListFSimple(x11,grid,fgrid->fsimple[i],&ID,bUseMon);
  AddDlgItemGroups(dlg,fgrid->w,fgrid->h,grid,bAutoPosition);

  DoneFGrid(fgrid);

  return dlg;
}



