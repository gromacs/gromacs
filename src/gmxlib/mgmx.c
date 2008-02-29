/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

#include <string.h>
#include <ctype.h>

#include <Xm/Xm.h>
#include <Xm/ScrolledW.h>
#include <Xm/ScrollBar.h>
#include <Xm/PushB.h>
#include <Xm/PushBP.h>
#include <Xm/ToggleB.h>
#include <Xm/ArrowB.h>
#include <Xm/CascadeB.h>
#include <Xm/Separator.h>
#include <Xm/DrawnB.h>
#include <Xm/DialogS.h>
#include <Xm/Scale.h>
#include <Xm/Frame.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/RowColumn.h>
#include <Xm/Label.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include <Xm/List.h>
#include <Xm/DrawingA.h>
#include <Xm/MenuShell.h>
#include <Xm/MessageB.h>
#include <X11/Shell.h>
#include <X11/StringDefs.h>

#include "typedefs.h"
#include "copyrite.h"
#include "smalloc.h"
#include "string2.h"
#include "statutil.h"
#include "wman.h"
#include "macros.h"
#include "widget.h"
#include "names.h"

typedef struct {
  char *label;
  char *desc;
  void (*cbfn)();
} t_button;

/* global variables */
#define  NARGS  20
static Arg      args[NARGS];
static Widget   FdlgCaller;
static int      helpw=-1,aboutw=-1,descw=-1,fdlgw=-1,gmxDialog;
static int      *pa_index,*pa_set_index,*fnm_index;
static bool     bFdlgUp=FALSE,bDone=FALSE,bDescSet=FALSE;
extern XmString empty_str;
static XmString desc_str;

#define GMXDLG    "gmxdlg"
#define GMXTOGGLE "gmxtoggle"
#define GMXEDIT   "gmxedit"
#define GMXBUTTON "gmxbutton"
#define GMXPD     "gmxpd"
#define GMXDESC   "gmxdesc"
#define GMXSEP    "gmxsep"
#define GMXHELP   "gmxhelp"
#define GMXSCROLL "gmxscroll"
#define GMXFILE   "gmxfile"

/* Height handling routines */
static windex TopW=0,NewTopW=0;
static void set_top_windex(windex t)
{
  NewTopW = t;
}

static windex top_windex(void)
{
  if (TopW == 0)
    gmx_fatal(FARGS,"No top widget");
  return TopW;
}

static void new_windex_row(void)
{
  if ((NewTopW == 0) && (TopW != 0)) {
    if (debug)
      fprintf(debug,"Empty windex row!\n");
  }
  else {
    TopW    = NewTopW;
    NewTopW = 0;
  }
}

/******************************************************************
 *
 *      C A L L B A C K       R O U T I N E S
 *
 ******************************************************************/

static void cancel_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  printf("Maybe next time...\n");
  thanx(stderr);
  exit(0);
}

static void ok_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  bDone = TRUE;
}

static void help_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  XtManageChild(get_widget(helpw));
}

static void help_ok_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  XtUnmanageChild(get_widget(helpw));
}

static void about_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  XtManageChild(get_widget(aboutw));
}

static void about_ok_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  XtUnmanageChild(get_widget(aboutw));
}

static void file_callback(Widget www,caddr_t client_data,caddr_t call_data)
{
  int      ftp,narg;
  Widget   fdlg;
  XmString xms;

  if (fdlgw != -1) {
    if (bFdlgUp)
      fprintf(stderr,"Only one file selector box at a time!\n");
    else {
      fdlg = get_widget(fdlgw);
      if ((ftp = get_widget_ftp(www)) != -1) {
	xms  = char2xms(ftp2filter(ftp));
	narg = 0;
	XtSetArg(args[narg],XmNdirMask,xms);        narg++;
	XtSetArg(args[narg],XmNpattern,xms);        narg++;
	XtSetArg(args[narg],XmNlistUpdated, False); narg++;
	
	XtSetValues(fdlg,args,narg);
	XmStringFree(xms);
      }
      XtManageChild(fdlg);
      bFdlgUp    = TRUE;
      FdlgCaller = get_widget_other(get_windex(www),TRUE);
    }
  }
}

static void file_ok_callback(Widget www,int *which,
			     XmFileSelectionBoxCallbackStruct *xmf)
{
  if (bFdlgUp) {
    if ((xmf->reason == XmCR_OK) && (xmf->length > 0))
      set_widget_dir(FdlgCaller,xmf->value);
  
    XtUnmanageChild(get_widget(fdlgw));
    bFdlgUp = FALSE;
  }
}

static void file_cancel_callback(Widget www,int *which,
				 XmFileSelectionBoxCallbackStruct *xmf)
{
  if (bFdlgUp) {
    XtUnmanageChild(get_widget(fdlgw));
    bFdlgUp = FALSE;
  }
}

static void enter_callback(Widget www,int *which,caddr_t call_data)
{
  int narg;

  if (descw != -1) {  
    if (have_windex_desc(get_windex(www))) {
      narg=0;
      XtSetArg(args[narg],XmNlabelString,get_widget_desc(www));  narg++;
      XtSetValues(get_widget(descw),args,narg);
      bDescSet = TRUE;
    }
  }
}

static void leave_callback(Widget www,int *which,caddr_t call_data)
{
  int narg;

  if (descw != -1) {  
    if (bDescSet) {
      narg=0;
      XtSetArg(args[narg],XmNlabelString,desc_str); narg++;
      XtSetValues(get_widget(descw),args,narg);
      bDescSet = FALSE;
    }
  }
}

/************************************************************************
 *
 *    S E T T I N G     U P    T H E   D I A L O G B O X   
 *
 ************************************************************************/
 
static void mk_desc_handlers(void)
{
  Widget www;
  int    i;
  
  for(i=0; (i<nwidget()); i++) {
    if (have_windex_desc(i)) {
      www  = get_widget(i);
      XtAddEventHandler(www,EnterWindowMask,True,
			(XtEventHandler) enter_callback,&i);
      XtAddEventHandler(www,LeaveWindowMask,True,
			(XtEventHandler) leave_callback,&i);
    }
  }
  empty_str = char2xms("");
}

static int mk_toggle(int parent,char *title,int top,int left,
		     bool bStatus,char *desc)
{
  int    narg;
  
  narg = 0;
  if (top == parent) {
    XtSetArg(args[narg],XmNtopAttachment, XmATTACH_FORM);   narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment, XmATTACH_WIDGET);  narg++;
    XtSetArg(args[narg],XmNtopWidget,     get_widget(top));  narg++;
  }
  XtSetArg(args[narg],XmNleftAttachment, XmATTACH_POSITION);  narg++;
  XtSetArg(args[narg],XmNleftPosition,   left+1);             narg++;
  XtSetArg(args[narg],XmNrightAttachment,XmATTACH_NONE);      narg++;
  /*XtSetArg(args[narg],XmNindicatorType,  XmONE_OF_MANY);      narg++;*/
  /*XtSetArg(args[narg],XmNvisibleWhenOff, False);              narg++;*/
  /*XtSetArg(args[narg],XmNselectPixmap*/
  XtSetArg(args[narg],XmNindicatorOn,  True);      narg++;
  if (bStatus) {
    XtSetArg(args[narg],XmNset,             True);            narg++;
  }
  XtSetArg(args[narg],XmNlabelString, char2xms(title)); narg++;
  
  return add_widget(XtCreateWidget(GMXTOGGLE,xmToggleButtonWidgetClass,
				   get_widget(parent),args,narg),desc);
}

static void mk_editor(int paindex,int parent,int top,int left,
		      char *label,char *initial_value,char *desc)
{
  enum { nwcTOGGLE, nwcTEXT, NWC };
  WidgetClass wc[NWC];
  char   *wlab[NWC];
  int    rleft[NWC] = { 1, 22 };
  int    rright[NWC]= { 21, 49 };
  int    ww[NWC];
  int    j,narg;
  
  /* Create & Position the label */
  wc[nwcTOGGLE]   = xmToggleButtonWidgetClass;
  wc[nwcTEXT]     = xmTextFieldWidgetClass;
  wlab[nwcTOGGLE] = GMXTOGGLE;
  wlab[nwcTEXT]   = GMXEDIT;
  
  for(j=0; (j<NWC); j++) {  
    narg = 0;
    if (top == parent) {
      XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
    }
    else {
      XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
      XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
    }
    XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
    XtSetArg(args[narg],XmNleftPosition,      left+rleft[j]);     narg++;
    if (j == nwcTOGGLE) {
      XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_NONE);   narg++;
    }
    else {
      XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_POSITION); narg++;
      XtSetArg(args[narg],XmNrightPosition,     left+rright[j]);    narg++;
    }
    XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
    if (j == nwcTEXT) {
      XtSetArg(args[narg],XmNvalue,           initial_value);     narg++;
    }
    else {
      XtSetArg(args[narg],XmNlabelString,     char2xms(label)); narg++;
    }
    if (debug)
      fprintf(debug,"There are %d args for %s\n",narg,wlab[j]);
    ww[j] = add_widget(XtCreateWidget(wlab[j],wc[j],
				      get_widget(parent),args,narg),desc);
  }
  pa_set_index[paindex] = ww[nwcTOGGLE];
  pa_index[paindex]     = ww[nwcTEXT];
}

static void mk_enumerated(int parent,int top,int left,
			  int *wlabel,int *wtextf,
			  char *label,char **but,char *desc)
{
  enum { nwcTOGGLE, nwcBUTTON, NWC };
  WidgetClass wc[NWC];
  int    pd;
  char   *wlab[NWC],buf[256];
  int    rleft[NWC] = { 1,  22 };
  int    rright[NWC]= { 21, 49 };
  int    ww[NWC];
  int    narg,i,j,wi;

  /* Create & Position the label */
  wc[nwcTOGGLE]   = xmToggleButtonWidgetClass;
  wc[nwcBUTTON]   = 0; 
  wlab[nwcTOGGLE] = GMXTOGGLE;
  wlab[nwcBUTTON] = GMXDLG;

  for(j=0; (j<NWC); j++) {  
    narg = 0;
    XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
    XtSetArg(args[narg],XmNleftPosition,      left+rleft[j]);     narg++;
    if (j != nwcTOGGLE) {
      XtSetArg(args[narg],XmNrightAttachment, XmATTACH_POSITION); narg++;
      XtSetArg(args[narg],XmNrightPosition,   left+rright[j]);    narg++;
    }
    if (top == parent) {
      XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
    }
    else {
      XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
      XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
    }
    XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
    
    if (j == nwcBUTTON) {
      ww[j] = add_widget(XmCreateOptionMenu(get_widget(parent),wlab[j],
					    args,narg),desc);
    }
    else {
      XtSetArg(args[narg], XmNlabelString, char2xms(label)); narg++;
      ww[j] = add_widget(XtCreateWidget(wlab[j],wc[j],
					get_widget(parent),args,narg),desc);
    }
  }
  /* Create the popup menu father */
  pd = add_widget(XmCreatePulldownMenu(get_widget(parent),GMXPD,NULL,0),
		  desc);
  
  /* Now create the popup menu children */
  set_windex_popup(pd,TRUE);
  set_parent(pd,get_widget(ww[nwcBUTTON]));
  set_widget_other(pd,get_widget(ww[nwcTOGGLE]));
  
  for(i=1; (but[i] != NULL); i++) {
    sprintf(buf,"%s = %s",desc,but[i]);
    narg = 0;
    XtSetArg(args[narg],XmNlabelString,char2xms(but[i])); narg++;
    wi = add_widget(XtCreateWidget(GMXBUTTON,xmPushButtonWidgetClass,
				   get_widget(pd),args,narg),buf);
    set_parent(wi,get_widget(pd));
  }

  /* Tell the option menu what to do */  
  XtVaSetValues(get_widget(ww[nwcBUTTON]),
		/*XmNlabelString, str,*/
		XmNsubMenuId, get_widget(pd),
		XmNentryBorder, 2,
		XmNwhichButton, 1,
		NULL);
  *wtextf = ww[nwcBUTTON];
}

static void mk_buttons(int parent,int top,int nb,t_button bbb[])
{
  int  i,narg,nw,left,right,bw;
  real dx;

  nw      = 1;
  dx      = (100-(nb+1)*nw)/(real)nb;
  right   = nw;
  for(i=0; (i<nb); i++) {
    left  = nw*(i+1)+i*dx;
    right = left+dx;
    narg  = 0;
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);      narg++;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(top));      narg++;
    if (i == 0) {
      XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);      narg++;
    }
    else {
      XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_POSITION);  narg++;
      XtSetArg(args[narg],XmNleftPosition,    left);               narg++;
    }
    if (i == nb-1) {
      XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);      narg++;
    } 
    else {
      XtSetArg(args[narg],XmNrightAttachment, XmATTACH_POSITION);  narg++;
      XtSetArg(args[narg],XmNrightPosition,   right);              narg++;
    }
    XtSetArg(args[narg],XmNlabelString,char2xms(bbb[i].label));    narg++;
    bw = add_widget(XtCreateWidget(GMXBUTTON,xmPushButtonWidgetClass,
				   get_widget(parent),args,narg),
		    bbb[i].desc);
    XtAddCallback(get_widget(bw),XmNactivateCallback,
		  (XtCallbackProc) bbb[i].cbfn,NULL);
  }
}

static XmString xs_str_array_to_xmstr(char *header,int ndesc,char *desc[])
{
  int      i;
  XmString xmstr;
  char     *cptr,*nlptr;
  const char *ptr;

  if (ndesc <= 0)
    return empty_str;
  else {
    xmstr = char2xms(header);
    for(i=0; (i<ndesc); i++) {
      xmstr = XmStringConcat(xmstr,XmStringSeparatorCreate());
      ptr   = check_tty(desc[i]);
      cptr  = wrap_lines(ptr,70,0,FALSE);
      ptr   = cptr;
      while ((nlptr = strchr(ptr,'\n')) != NULL) {
	*nlptr='\0';
	xmstr = XmStringConcat(xmstr,char2xms(ptr));
	xmstr = XmStringConcat(xmstr,XmStringSeparatorCreate());
	ptr   = nlptr+1;
	while (*ptr == '\n')
	  ptr++;
      }
      xmstr = XmStringConcat(xmstr,char2xms(ptr));
      sfree(cptr);
    }
    return xmstr;
  }
}

static windex mk_separator(windex parent,windex topw)
{
  int narg,sep;
  
  narg = 0;
  if (nwidget() > 0) {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
    if (topw < 0) 
      topw = nwidget()-1;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(topw));  narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);   narg++;
  }
  XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);   narg++;
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);   narg++;
  
  sep = add_widget(XtCreateWidget(GMXSEP,xmSeparatorWidgetClass,
				  get_widget(parent),args,narg),NULL);
  
  set_top_windex(sep);
				  
  return sep;
}

static int mk_helplabel(int parent,int top)
{
  int  narg;
  char buf[] = "Place the mouse over an item to get information";
  
  desc_str = char2xms(buf);
  narg = 0;
  XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET); narg++;
  XtSetArg(args[narg],XmNtopWidget,       get_widget(top)); narg++;
  XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);   narg++;
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);   narg++;
  XtSetArg(args[narg],XmNalignment,       XmALIGNMENT_BEGINNING);  narg++;
  XtSetArg(args[narg],XmNlabelString,     desc_str);        narg++;

  return add_widget(XmCreateLabel(get_widget(parent),GMXDESC,args,narg),NULL);
}

static windex mk_filedlgs(int parent,int top,int nfile,
			  t_filenm fnm[],int fnm_index[])
{
  enum { nwcLABEL, nwcTEXT, nwcFDLG, NWC };
  WidgetClass wc[NWC];
  int    left[NWC]  = { 1,  12, 39 };
  int    right[NWC] = { 11, 38, 49 };
  char   **wname;
  int    i,j,ftp,narg,dx,www[NWC];
  Widget topw;
  char   *fn,dbuf[256],buf[255];

  wc[nwcTEXT]   = xmTextFieldWidgetClass;
  wc[nwcFDLG]   = xmPushButtonWidgetClass;
  
  /* Make the compiler happy */
  topw = (Widget) 0;
  
  snew(wname,NWC);
  wname[nwcLABEL]   = GMXTOGGLE;
  wname[nwcTEXT]    = GMXEDIT;
  wname[nwcFDLG]    = GMXBUTTON;
  for(i=0; (i<nfile); i++) {
    ftp = fnm[i].ftp;
    dx  = (i % 2)*50;
    sprintf(dbuf,"%s [%s]",ftp2desc(ftp),fileopt(fnm[i].flag,buf,255));
    if ((fn = strrchr(fnm[i].fns[0],'/')) == NULL)
      fn = fnm[i].fns[0];
    else
      fn++;
    
    if (is_optional(&(fnm[i])))
      wc[nwcLABEL] = xmToggleButtonWidgetClass;
    else
      wc[nwcLABEL]  = xmLabelWidgetClass;
      
    for(j=0; (j<NWC); j++) {
      narg    = 0;
      if (i < 2) {
	if (top == parent) {
	  XtSetArg(args[narg],XmNtopAttachment, XmATTACH_FORM); narg++;
	}
	else {
	  XtSetArg(args[narg],XmNtopAttachment, XmATTACH_WIDGET); narg++;
	  XtSetArg(args[narg],XmNtopWidget,     get_widget(top)); narg++;
	}
      }
      else {
	XtSetArg(args[narg],XmNtopAttachment, XmATTACH_WIDGET);    narg++;
	XtSetArg(args[narg],XmNtopWidget,     topw);               narg++;
      }
      
      XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_POSITION);  narg++;
      XtSetArg(args[narg],XmNleftPosition,    dx+left[j]);         narg++;
      
      
      if (j != nwcLABEL) {
	XtSetArg(args[narg],XmNrightAttachment, XmATTACH_POSITION);  narg++;
	XtSetArg(args[narg],XmNrightPosition,   dx+right[j]);         narg++;
      }
      switch (j) {
      case nwcLABEL:
	XtSetArg(args[narg],XmNlabelString,   char2xms(fnm[i].opt)); narg++;
	break;
      case nwcFDLG:
	XtSetArg(args[narg],XmNlabelString,   char2xms("Browse")); narg++;
	break;
      case nwcTEXT:
	XtSetArg(args[narg],XmNvalue,         fnm[i].fns[0]);   narg++;
	break;
      }
      www[j] = add_widget(XtCreateWidget(wname[j],wc[j],
					 get_widget(parent),args,narg),dbuf);
    }
    fnm_index[i] = www[nwcTEXT];
    set_windex_orignm(www[nwcTEXT],fnm[i].fns[0]);
    set_widget_ftp(www[nwcFDLG],ftp);
    set_widget_other(www[nwcFDLG],get_widget(www[nwcTEXT]));
    
    if (is_optional(&(fnm[i])))
      set_widget_other(www[nwcTEXT],get_widget(www[nwcLABEL]));
      
    XtAddCallback(get_widget(www[nwcFDLG]),XmNactivateCallback,
		  (XtCallbackProc) file_callback,NULL);
    if ((i % 2) == 1)
      topw = get_widget(www[nwcTEXT]);
  }
  sfree(wname);
  
  if (nfile > 0) 
    return mk_separator(parent,nwidget()-2);
  else 
    return 0;
}

static bool motif_hidden(t_pargs *pa)
{
  return (is_hidden(pa) || 
	  (strcmp(pa->option,"-X") == 0) ||
	  (strcmp(pa->option,"-h") == 0));
}

static windex mk_pargs(int parent,int npargs,t_pargs pa[],int pa_index[],
		       windex topw)
{
  /**************************************************************
   *
   * Create all the children in two passes:
   * 1st the checkboxes (etBOOL variables)
   * 2nd the editable fields
   * Return the index of the separator (or the last created widget)
   *
   **************************************************************/
  int    icb,i,npa,dummy,separator,nbool,nedit,nenum,nelem,left;
  real   bwidth;
  char   buf[132],descbuf[132];

  /* Make the compiler happy */
  separator = 0;
  
  /* First round just count booleans and editors */
  nbool = nedit = nenum = 0;
  for(i=0; (i<npargs); i++) {
    if (!motif_hidden(&(pa[i]))) {
      if (pa[i].type == etBOOL)
	nbool++;
      else if (pa[i].type == etENUM)
	nenum++;
      else
	nedit++;
    }
  }
  if (nbool > 6)
    nbool = 6;
  bwidth = (100.0/nbool);

  if (nenum > 2)
    nenum = 2;
    
  set_top_windex(topw);
      
  for(icb=0; (icb<3); icb++) {
    new_windex_row();
    npa  = 0;
    for(i=0; (i<npargs); i++) {
      topw  = top_windex();
      nelem = 0;
      if (!motif_hidden(&(pa[i]))) {
	sprintf(descbuf,"%s (%s)",pa[i].desc,argtp[pa[i].type]);
	switch (pa[i].type) {
	case etBOOL:
	  if (icb == 0) {
	    nelem = nbool;
	    left  = ((npa % nbool)*bwidth)+0.5;
	    if (debug)
	      fprintf(debug,"%s,%d: nbool %d, bwidth %g, left %d, topw %d\n",
		      __FILE__,__LINE__,nbool,bwidth,left,topw);
	    pa_index[i] = mk_toggle(parent,pa[i].option,topw,
				    left,*(pa[i].u.b),descbuf);
	  }
	  break;
	case etENUM:
	  if (icb == 1) {
	    nelem = nenum;
	    left  = (npa % nenum)*50;
	    mk_enumerated(parent,topw,left,&dummy,&(pa_index[i]),
			  pa[i].option,pa[i].u.c,descbuf); 
	  }
	  break;
	default:
	  if (icb == 2) {
	    nelem = 2;
	    switch (pa[i].type) {
	    case etREAL:
	      sprintf(buf,"%g",*(pa[i].u.r));
	      break;
	    case etRVEC:
	      sprintf(buf,"%g %g %g",(*pa[i].u.rv)[XX],
		      (*pa[i].u.rv)[YY],(*pa[i].u.rv)[ZZ]);
	      break;
	    case etINT:
	      sprintf(buf,"%d",*(pa[i].u.i));
	      break;
	    case etSTR:
	      if (*(pa[i].u.c) != NULL)
		strcpy(buf,*(pa[i].u.c));
	      else
		buf[0] = '\0';
	      break;
	    }
	    if (debug)
	      fprintf(debug,"%s,%d: buf = %s\n",__FILE__,__LINE__,buf);
	    left = (npa % nelem)*50;
	    mk_editor(i,parent,topw,left,pa[i].option,buf,descbuf); 
	  }
	  break;
	}
      }
      if (nelem) {
	npa ++;
	set_top_windex(pa_index[i]);
	if ((npa % nelem) == 0)
	  new_windex_row();
      }
    }
    
    if (npa > 0) {
      new_windex_row();
      separator = mk_separator(parent,top_windex());
    }
  }
  return separator;
}

static void append_str(char **buf,int *blen,int *maxlen,char *str,
		       int indent)
{
#define DELTA 256
  int  i,slen,width=80;
  char *nptr;
  const char *ptr;

  ptr = check_tty(str);
  if (indent > 0) {
    slen=strlen(ptr);
    snew(nptr,slen+8);
    sprintf(nptr,"* %s",ptr);
    str = wrap_lines(nptr,width,indent,FALSE);
  }
  else
    str = wrap_lines(ptr,width,indent,FALSE);
  
  /*while ((ptr = strstr(str,"\n\n")) != 0)
    *ptr = ' ';
    */
  slen = strlen(str);
    
  while (*blen+slen+1 > *maxlen) {
    srenew((*buf),*maxlen+DELTA);
    for(i=(*maxlen); (i<(*maxlen)+DELTA); i++)
      (*buf)[i] = '\0';
    (*maxlen) += DELTA;
  }
  strcat(*buf,str);
  strcat(*buf,"\n");
  *blen+=slen+1;
#undef DELTA
}

static char *concat_str(char *dtitle,int ndesc,char *desc[],
			char *btitle,int nbugs,char *bugs[])
{
  char *descer;
  char *ptr  = NULL;
  char *buf;
  int  i,blen=0,maxlen=0,dlen;

  append_str(&ptr,&blen,&maxlen,dtitle,0);
  buf=strdup(dtitle);
  for(i=0; (buf[i] != '\0'); i++)
    buf[i] = '-';
  append_str(&ptr,&blen,&maxlen,buf,0);
  sfree(buf);
  if (ndesc == 0) 
    append_str(&ptr,&blen,&maxlen,"none?",0);

  /* Count the length of the strings together */    
  dlen   = 0;
  for(i=0; (i<ndesc); i++) {
    dlen += strlen(desc[i])+1;
  }
  snew(descer,dlen+1);
  for(i=0; (i<ndesc); i++) {
    strcat(descer,desc[i]);
    if (i < ndesc-1)
      strcat(descer," ");
  }
  append_str(&ptr,&blen,&maxlen,descer,0);
  sfree(descer);
  if (nbugs > 0) {
    append_str(&ptr,&blen,&maxlen," ",0);
    append_str(&ptr,&blen,&maxlen,btitle,0);
    buf=strdup(btitle);
    for(i=0; (buf[i] != '\0'); i++)
      buf[i] = '-';
    append_str(&ptr,&blen,&maxlen,buf,0);
    sfree(buf);
  }
  for(i=0; (i<nbugs); i++) {
    append_str(&ptr,&blen,&maxlen,bugs[i],2);
  }  
  return ptr;
}

static int low_mk_help(Widget parent,
		       char *dtitle,int ndesc,char *desc[],
		       char *btitle,int nbugs,char *bugs[],
		       XtCallbackProc ok_callback)
{
  windex   text,sep,ok,sw;
  char     buf[256],*ptr;
  int      narg,awin;

  /* Create the mother of all help windows */
  sprintf(buf,"Gromacs Help - %s",ShortProgram());
  narg = 0;
  XtSetArg(args[narg],XmNdialogTitle,char2xms(buf)); narg++;
  awin  = add_widget(XmCreateFormDialog(parent,GMXHELP,args,narg),buf);
  
  ptr = concat_str(dtitle,ndesc,desc,btitle,nbugs,bugs);
  
  /* Now create the contents */
  narg   = 0; 
  XtSetArg(args[narg],XmNheight,          480); narg++;
  XtSetArg(args[narg],XmNwidth,           570); narg++;
  XtSetArg(args[narg],XmNeditMode,        XmMULTI_LINE_EDIT); narg++;
  XtSetArg(args[narg],XmNeditable,        FALSE);             narg++;
  XtSetArg(args[narg],XmNvalue,           ptr);               narg++;
  XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
  XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);     narg++;
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);     narg++;
  XtSetArg(args[narg],XmNscrollingPolicy, XmAUTOMATIC);       narg++;
  XtSetArg(args[narg],XmNscrollBarDisplayPolicy,XmAUTOMATIC); narg++;
  text   = add_widget(XmCreateScrolledText(get_widget(awin),GMXHELP,
					   args,narg),
		      "There is supposed to be useful information in the help window");
  sw     = add_widget(XtParent(get_widget(text)),NULL);
  sep    = mk_separator(awin,sw);
  
  narg   = 0;
  XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);    narg++;
  XtSetArg(args[narg],XmNtopWidget,       get_widget(sep));    narg++;
  XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);      narg++;
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);      narg++;
  XtSetArg(args[narg],XmNalignment,       XmALIGNMENT_CENTER); narg++;
  XtSetArg(args[narg],XmNbottomAttachment,XmATTACH_FORM);      narg++;
  XtSetArg(args[narg],XmNlabelString,     char2xms("OK"));     narg++;
  ok     = add_widget(XmCreatePushButton(get_widget(awin),
					 GMXBUTTON,args,narg),
		      "Press OK to close the helpwindow");
  XtAddCallback(get_widget(ok),XmNactivateCallback,ok_callback,NULL);
		
  /* Now manage all except the mother of all help windows */
  XtManageChild(get_widget(ok));
  XtManageChild(get_widget(sw));
  XtManageChild(get_widget(text));
  XtManageChild(get_widget(sep));
  
  return awin;
}

static int mk_help(Widget base,int ndesc,char *desc[],int nbugs,char *bugs[])
{
  return low_mk_help(base,"DESCRIPTION:",ndesc,desc,
		     "DIAGNOSTICS:",nbugs,bugs,(XtCallbackProc) help_ok_callback);
}

static int mk_about(Widget base)
{
  char *about[] = {
    "Starting with version 2.0 GROMACS has an X/Motif graphical user interface (GUI) to all",
    "programs. The command line interface is translated automatically into",
    "a dialog box which should make life a bit easier for those new to the",
    "programs. A drawback of this approach is that the options are usually",
    "short, and therefore not very desriptive. In the GUI there is a line with",
    "extra information at the bottom, which informs you about the option under",
    "the mouse cursor.[PAR]",
    "Note that in the following description all possible elements of the GUI",
    "are described, but your program need not have all these option types.[PAR]",
    "In the upper pane of the dialog box you find the file options. These can",
    "be manually edited or using a File selector box under the [BB]Browse[bb]",
    "button. Note that some files are optional, only if you tick the box before the",
    "option the program will actually use the file. If you use the File selector",
    "box to select a file the option will be used automatically. If you change your",
    "mind about the option you can turn it off again using the tick box.[PAR]",
    "In the second pane of the dialog box you will find a number of Toggle buttons",
    "with which you can turn flags on and off.[PAR]",
    "In the third pane of the dialog box there are options which can take",
    "a limited set of values (enumerated in programmers jargon). This is implemented",
    "in the GUI using a popup menu.[PAR]",
    "In the fourth pane you find the options which require you to type (yuckie!)",
    "a numeric or textual argument (note that the description window indicates",
    "which kind of argument). The validity of what you type is *not* checked",
    "by the GUI as this does not know what the options mean. Instead the program",
    "using the GUI will (hopefully) verify that your input is meaningful.[PAR]",
    "In the fifth pane you find the usual buttons, [BB]OK[bb], [BB]Cancel[bb],",
    "[BB]Help[bb], and [BB]About[bb] which presumably don't need any explanation.[PAR]",
    "The GUI was written by David van der Spoel (comments to spoel@xray.bmc.uu.se)[PAR]",
    "And hey:[BR]"
  };
#define NABOUT asize(about)

  static char *mbugs[] = {
    "The file selector box generates a core dump under Linux/lesstif0.86",
    "The file selector box does not work with multiple file selections yet",
    "It took about 1500 lines of pretty ugly C code to get this dialog box working"
  };

  char tmpstr[256];
  int dum;
  
  return low_mk_help(base,
		     "ABOUT THE GROMACS MOTIF USER INTERFACE",asize(about),about,
		     "PROBLEMS IN THE GUI",asize(mbugs),mbugs,
		     (XtCallbackProc) about_ok_callback);
}

static void mk_fdlg(Widget base)
{
  int narg;
  
  narg  = 0;
  XtSetArg(args[narg],XmNdialogTitle,char2xms("GMX File selector")); narg++;
  fdlgw = add_widget(XmCreateFileSelectionDialog(base,GMXFILE,args,narg),NULL);
  XtAddCallback(get_widget(fdlgw),XmNokCallback,
		(XtCallbackProc) file_ok_callback,NULL);
  XtAddCallback(XmFileSelectionBoxGetChild(get_widget(fdlgw),
					   XmDIALOG_CANCEL_BUTTON),
		XmNactivateCallback,
		(XtCallbackProc) file_cancel_callback,NULL);
  XtUnmanageChild(XmFileSelectionBoxGetChild(get_widget(fdlgw),
					     XmDIALOG_HELP_BUTTON));

  /* Make the filter fixed... */
  narg = 0;
  XtSetArg(args[narg],XmNeditable,False); narg++;
  XtSetValues(XmFileSelectionBoxGetChild(get_widget(fdlgw),
					 XmDIALOG_FILTER_TEXT),args,narg);
}

static void mk_gui(Widget gmxBase,
		   int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
		   int ndesc,char *desc[],int nbugs,char *bugs[])
{
  /****************************************************************
   *
   *  M A K E   G U I   F O R   G R O M A C S
   * 
   ****************************************************************/
  t_button bbb[] = {
    { "OK",     "Press OK to accept current settings and continue",
      &ok_callback },
    { "Cancel", "Press Cancel to quit the program",
      &cancel_callback },
    { "Help",   "Press Help for more information about the program",
      &help_callback },
    { "About",  "Press About for information about using this dialog box",
      &about_callback },
  };
#define NBUT asize(bbb)
  int  sep1,sep2,widg0;
  int  i;

  /* Make the compiler happy */
  sep1   = 0;
    
  /* Create the help window */
  helpw  = mk_help(gmxBase,ndesc,desc,nbugs,bugs); 

  /* Create the about window */
  aboutw = mk_about(gmxBase);
		   
  /* Create the file selector box */
  mk_fdlg(gmxBase);
  
  /* Create the dialog box! */
  gmxDialog = add_widget(XmCreateForm(gmxBase,GMXDLG,NULL,0),NULL);
  XtManageChild(get_widget(gmxDialog));

  widg0 = nwidget();
  
  /* Create buttons for file dialogboxes */
  if (nfile > 0) {
    snew(fnm_index,nfile);
    sep1 = mk_filedlgs(gmxDialog,gmxDialog,nfile,fnm,fnm_index);
  }
  
  /* Create the checkboxes and editable fields */
  if (npargs > 0) {
    snew(pa_index,npargs);
    snew(pa_set_index,npargs);
    sep1 = mk_pargs(gmxDialog,npargs,pa,pa_index,sep1);
  }
  else
    sep1 = nwidget()-1;
  
  /* Create & count the buttons */  
  mk_buttons(gmxDialog,sep1,NBUT,bbb);
  
  /* Create & Position the separator */
  sep2 = mk_separator(gmxDialog,-1);
  
  /* Create help label */
  descw = mk_helplabel(gmxDialog,sep2);

  /* Make eventhandlers for the help line */  
  mk_desc_handlers();
  
  /* Give the children a father or mother */
  for(i=widg0; (i<nwidget()); i++)
    if (!get_windex_popup(i))
      XtManageChild(get_widget(i));
}

/***************************************************************************
 *
 *                         M A I N      L O O P
 *
 ***************************************************************************/

static void MyMainLoop(XtAppContext appcontext,Widget gmxBase,
		       int nfile,t_filenm fnm[],int npargs,t_pargs pa[])
{
  int        i,narg;
  Widget     www,wsub;
  XEvent     event;
  XmString   xms;
  Boolean    xmrb,bbb;
  char       *fn,buf[256],*ptr;
  double     ddd,dx,dy,dz;
  int        iii;
  
  while (!bDone) {
    XtAppNextEvent(appcontext,&event);
    XtDispatchEvent(&event);
  }
  /* Extract all the information from the X widgets */
  for(i=0; (i<nfile); i++) {
    www  = get_widget(fnm_index[i]);
    narg = 0;
    XtSetArg(args[narg], XmNvalue, &fn); narg++;
    XtGetValues(www,args,narg);
    sprintf(buf,"%s%s",get_widget_dir(fnm_index[i]),fn);
    XtFree(fn);
    sfree(fnm[i].fns[0]);
    fnm[i].fns[0] = strdup(buf);
    if (is_optional(&(fnm[i]))) {
      www = get_widget_other(fnm_index[i],FALSE);
      if (www != 0) {
	narg = 0;
	XtSetArg(args[narg],XmNset,&xmrb); narg++;
	XtGetValues(www,args,narg);
	if (xmrb)
	  fnm[i].flag = fnm[i].flag | ffSET;
      }
      else
	gmx_fatal(FARGS,"No toggle button for optional file (option %s)",
		    fnm[i].opt);
      if (strcmp(fnm[i].fns[0],get_windex_orignm(fnm_index[i])) != 0) {
	if (debug) {
	  fprintf(debug,"File corr. to option %s has been modified from\n"
		  "'%s' to '%s'\n",fnm[i].opt,
		  get_windex_orignm(fnm_index[i]),fnm[i].fns[0]);
	}
	fnm[i].flag = fnm[i].flag | ffSET;
      }
    }
    if (debug)
      fprintf(debug,"%s,%d: File is now %s\n",__FILE__,__LINE__,buf);
    
  }
  for(i=0; (i<npargs); i++) {
    if (!is_hidden(&pa[i])) {
      if (pa[i].type == etBOOL) {
	narg = 0;
	XtSetArg(args[narg],XmNset, &bbb);            narg++;
	XtGetValues(get_widget(pa_index[i]),args,narg);
	*(pa[i].u.b) = (bbb == True);
	if (debug)
	  fprintf(debug,"%s,%d: Boolean (windex %d) %s = %s\n",__FILE__,__LINE__,
		  pa_index[i],pa[i].option,bool_names[*(pa[i].u.b)]);
      }
      else {
	/* Check whether it is set */
	narg = 0;
	XtSetArg(args[narg],XmNset, &bbb);            narg++;
	XtGetValues(get_widget(pa_set_index[i]),args,narg);
	pa[i].bSet = (bbb == True);
	
	/* Now extract the value */
	if (pa[i].type == etENUM) {
	  /* First get the selected widget */
	  narg = 0;
	  XtSetArg(args[narg],XmNmenuHistory, &wsub);            narg++;
	  XtGetValues(get_widget(pa_index[i]),args,narg);
	  /* Now get it's label! */
	  narg = 0;
	  XtSetArg(args[narg],XmNlabelString, &xms);            narg++;
	  XtGetValues(wsub,args,narg);
	  ptr = xms2char(xms);
	}
	else {
	  narg = 0;
	  XtSetArg(args[narg],XmNvalue, &ptr);            narg++;
	  XtGetValues(get_widget(pa_index[i]),args,narg);
	}
	if (debug)
	  fprintf(debug,"%s,%d: I found option %s value %s\n",
		  __FILE__,__LINE__,pa[i].option,ptr);
	switch (pa[i].type) {
	case etREAL:
	  if (sscanf(ptr,"%lf",&ddd) != 1)
	    fprintf(stderr,"Warning: invalid entry (%s) for real value %s, using default %g\n",
		    ptr,pa[i].option,*(pa[i].u.r));
	  else
	    *pa[i].u.r = ddd;
	  break;
	case etRVEC:
	  if (sscanf(ptr,"%lf%lf%lf",&dx,&dy,&dz) != 3)
	    fprintf(stderr,"Warning: invalid entry (%s) for rvec value %s, using default (%g,%g,%g)\n",
		    ptr,pa[i].option,(*pa[i].u.rv)[XX],(*pa[i].u.rv)[YY],
		    (*pa[i].u.rv)[ZZ]);
	  else {
	    (*pa[i].u.rv)[XX] = dx;
	    (*pa[i].u.rv)[YY] = dy;
	    (*pa[i].u.rv)[ZZ] = dz;
	  }
	  break;
	case etINT:
	  if (sscanf(ptr,"%d",&iii) != 1)
	    fprintf(stderr,"Warning: invalid entry (%s) for integer value %s, using default %d\n",
		    ptr,pa[i].option,*(pa[i].u.i));
	  else
	    *pa[i].u.i = iii;
	  break;
	case etSTR:
	  *(pa[i].u.c) = strdup(ptr);
	  break;
	case etENUM:
	  pa[i].u.c[0] = strdup(ptr);
	  break;
	}
	sfree(ptr);
      }
    }
  }
  /* Clean up windows */
  XtUnmanageChild(get_widget(gmxDialog));
  XtUnmanageChild(get_widget(fdlgw));
  XtUnmanageChild(get_widget(helpw));
  XtUnrealizeWidget(gmxBase);
  XtDestroyApplicationContext(appcontext);
}

void gmx_gui(int *argc,char *argv[],
	     int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
	     int ndesc,char *desc[],int nbugs,char *bugs[])
{
  Widget           gmxBase;
  XtAppContext     appcontext;
  String           Fallbacks[] = {
    "*gmx*background:           lightgrey",
    /*"*gmxdlg.background:        lightgrey",
    "*gmxbutton.background:     lightskyblue",
    "*gmxtoggle.background:     lightgrey",
    "*gmxedit.background:       lightgoldenrod1",
    "*gmxdesc.background:       lightsalmon1",
    "*gmxsep.background:        darkgrey",
    "*gmxhelp.background:       lightgoldenrod1",
    "*gmxhelpSW.background:     lightgrey",
    "*gmxhelpSW.*.background:   lightgrey",
    "*gmxfile.*.background:     lightgrey",
    "*gmxfile.Text.background:  lightgoldenrod1",*/
    NULL
  };
  
  /* Initialize toolkit and parse command line options. */
  gmxBase = XtOpenApplication(&appcontext,"gmx",NULL,0,argc,argv,Fallbacks,
			      applicationShellWidgetClass,NULL,0);
			      
  mk_gui(gmxBase,nfile,fnm,npargs,pa,ndesc,desc,nbugs,bugs);
  XtRealizeWidget(gmxBase);
  MyMainLoop(appcontext,gmxBase,nfile,fnm,npargs,pa);
}
