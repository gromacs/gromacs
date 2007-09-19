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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <Xstuff.h>
#include <xutil.h>
#include <smalloc.h>
#include <macros.h>
#include <futil.h>
#include <string2.h>

#define YSPACE 2

typedef struct {
  t_windata   wd;		/* Window structure			*/
  int         nlines,top;	/* Number of lines, current top line	*/
  char        **lines;		/* The strings				*/
  int         wheight,wwidth;	/* The size of the window in chars	*/
  XFontStruct *font;		/* Font					*/
  unsigned long       fg,bg;		/* Colours				*/
} t_scrollw;

static void calc_scrollw(t_scrollw *sw,int w,int h)
{
  sw->wd.width=w;
  sw->wd.height=h;
  sw->wheight=h/(YSPACE+XTextHeight(sw->font));
  sw->wwidth=w/XTextWidth(sw->font,"W",1);
}

static bool SWCallback(t_x11 *x11,XEvent *event,Window w,void *data)
{
  t_scrollw *sw;
  int       i,y,nl,barw,btop,bheight;
  real      h,frac;

  sw=(t_scrollw *)data;

  /* Calc some bar data */
  barw=20;
  h=XTextHeight(sw->font)+YSPACE;
  frac=min(1.0,((real)sw->wheight)/((real)sw->nlines));
  btop=((((real)sw->top)/((real)sw->nlines)))*(sw->wd.height);
  bheight=frac*sw->wd.height;
  bheight-=bheight/h;

  switch(event->type) {
  case Expose:
    nl=min(sw->nlines,sw->top+sw->wheight);
    y=0;

    XClearWindow(x11->disp,w);
#ifdef DEBUG
    printf("btop: %d, bheight: %d, frac: %e, h: %e\n",btop,bheight,frac,h);
#endif
    /* Draw the bar */
    XSetForeground(x11->disp,x11->gc,LIGHTGREY);
    XFillRectangle(x11->disp,w,x11->gc,2,btop+2,barw-4,bheight-4);
    XDrawLine(x11->disp,w,x11->gc,barw,0,barw,sw->wd.height);

    /* Draw the text */
    XSetForeground(x11->disp,x11->gc,sw->fg);
    for(i=sw->top; (i<nl); i++) {
      SpecialTextInRect(x11,sw->font,w,
			sw->lines[i],barw+2,y,sw->wd.width-barw-4,(int)h,
			eXLeft,eYCenter);
      y+=h;
    }
    XSetForeground(x11->disp,x11->gc,x11->fg);
    break;
  case ConfigureNotify:
    calc_scrollw(sw,event->xconfigure.width,event->xconfigure.height);
    break;
  case ButtonPress:
    if (event->xbutton.x < barw) {
      int y=event->xbutton.y;

      if (sw->nlines > sw->wheight) {
	if (y<btop) 
	  sw->top=max(0,sw->top-1);
	else if (y>btop+bheight) 
	  sw->top=min(sw->nlines-sw->wheight,sw->top+1);
	else
	  break;
	ExposeWin(x11->disp,sw->wd.self);
      }
    }
    break;
  default:
    break;
  }

  return FALSE;
}

t_scrollw *init_scrollw(t_x11 *x11,Window parent,int x,int y,int w,int h,
			unsigned long fg,unsigned long bg)
{
  t_scrollw *sw;

  snew(sw,1);

  InitWin(&sw->wd,x,y,w,h,1,"Scroll Window");
  sw->fg=fg;
  sw->bg=bg;
  sw->font=x11->font;
  sw->wd.self=XCreateSimpleWindow(x11->disp,parent,x,y,w,h,
				  sw->wd.bwidth,fg,bg);
  x11->RegisterCallback(x11,sw->wd.self,parent,SWCallback,sw);
  x11->SetInputMask(x11,sw->wd.self,ExposureMask | ButtonPressMask |
		    StructureNotifyMask);
  calc_scrollw(sw,w,h);

  return sw;
}

void show_scrollw(t_x11 *x11,t_scrollw *sw)
{
  XMapWindow(x11->disp,sw->wd.self);
}

char *tab2spc(char *buf)
{
  char *buf2;
  char *new;
  int  i,j;

  snew(buf2,8*strlen(buf)+1);
  for(i=j=0; (buf[i]!='\0'); i++)
    if (buf[i]=='\t') 
      do {
	buf2[j++]=' ';
      } while ((j % 8)!=0);
    else
      buf2[j++]=buf[i];
  buf2[j]='\0';
  new=strdup(buf2);
  sfree(buf2);
  return new;
}

void read_lines(FILE *in,t_scrollw *sw)
{
  char buf[1024];

  while (fgets2(buf,1023,in)) {
    sw->nlines++;
    srenew(sw->lines,sw->nlines);
    sw->lines[sw->nlines-1]=tab2spc(buf);
  }
}

void main(int argc, char *argv[])
{
  t_x11     *x11;
  t_scrollw *sw;

  if ((x11=GetX11(&argc,argv))==NULL) {
    fprintf(stderr,"No X!\n");
    exit(1);
  }
  sw=init_scrollw(x11,x11->root,0,0,600,200,WHITE,BLACK);
  read_lines(stdin,sw);
  show_scrollw(x11,sw);
  x11->MainLoop(x11);

  x11->CleanUp(x11);
}
