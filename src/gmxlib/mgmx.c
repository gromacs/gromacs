#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <Xm/Xm.h>
#include <Xm/ScrolledW.h>
#include <Xm/ScrollBar.h>
#include <Xm/PushB.h>
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
static int      helpw=-1,descw=-1,fdlgw=-1,gmxDialog;
static int      *pa_index,*fnm_index;
static bool     bFdlgUp=FALSE,bDone=FALSE;
extern XmString empty_str;

static void help_callback();
static void help_ok_callback();
static void cancel_callback();
static void ok_callback();
static void enter_callback();
static void leave_callback();
static void file_callback();
static void file_ok_callback();
static void file_cancel_callback();

void mk_desc_callbacks()
{
  Widget www;
  int    i,narg;
  
  for(i=0; (i<nwidget()); i++) {
    if (have_windex_desc(i)) {
      www  = get_widget(i);
      XtAddEventHandler(www,EnterWindowMask,True,enter_callback,&i);
      XtAddEventHandler(www,LeaveWindowMask,True,leave_callback,&i);
    }
  }
  empty_str = char2xms("");
}

static int mk_toggle(int parent,char *title,int top,int left,int right,
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
  XtSetArg(args[narg],XmNleftPosition,   left);               narg++;
  XtSetArg(args[narg],XmNrightAttachment,XmATTACH_POSITION);  narg++;
  XtSetArg(args[narg],XmNrightPosition,  right);              narg++;
  
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_NONE);     narg++;
  if (bStatus) {
    XtSetArg(args[narg],XmNset,             True);            narg++;
  }
  return add_widget(XtCreateWidget(title,xmToggleButtonWidgetClass,
				   get_widget(parent),args,narg),desc);
}

static void mk_editor(int parent,int top,int left,int right,
		      int *wlabel,int *wtextf,
		      char *label,char *initial_value,char *desc)
{
  int    narg,lab_right,but_left;
  
  /* Create & Position the label */
  lab_right = left + 0.3*(right-left);
  but_left  = lab_right + 2;
  narg = 0;
  XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNleftPosition,      left);              narg++;
  if (top == parent) {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
  }
  XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
  *wlabel = add_widget(XtCreateWidget(label,xmLabelWidgetClass,
				      get_widget(parent),args,narg),desc);
  
  /* Create & Position the textfield */
  narg = 0;
  XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNleftPosition,      but_left);          narg++;
  XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNrightPosition,     right);             narg++;
  if (top == parent) {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
  }
  XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
  XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_NONE);     narg++;
  XtSetArg(args[narg],XmNvalue,             initial_value);     narg++;
  *wtextf = add_widget(XtCreateWidget(initial_value,xmTextFieldWidgetClass,
				      get_widget(parent),args,narg),desc);
}

char **mk_but_array(char *val,int *n)
{
  int  i,j,nn;
  char **but=NULL;
  char buf[256];
  bool bCp;
  
  nn  = 0;
  j   = 0;
  bCp = FALSE;
  for(i=0; (val[i] != '\0'); i++) {
    if (!isspace(val[i])) {
      if (!bCp) {
	/* New string */
	bCp = TRUE;
	j   = 0;
      }
      buf[j++] = val[i];
    }
    else if (bCp) {
      bCp = FALSE;
      buf[j++] = '\0';
    }
    srenew(but,nn+1);
    but[nn] = strdup(buf);
    nn++;
  }
  if (debug) {
    fprintf(debug,"Found %d elements in '%s'\n",nn,val);
    for(i=0; (i<nn); i++)
      fprintf(debug,"Elem[%5d] = %s\n",i,but[i]);
  }
  *n = nn;
  return but;
}

static void mk_enumerated(int parent,int top,int left,int right,
			  int *wlabel,int *wtextf,
			  char *label,char *initial_value,char *desc)
{
  int  narg,lab_right,but_left,nbut,i;
  char **but;
  
  /* Create & Position the label */
  lab_right = left + 0.3*(right-left);
  but_left  = lab_right + 2;
  narg = 0;
  XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNleftPosition,      left);              narg++;
  if (top == parent) {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
  }
  XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
  *wlabel = add_widget(XtCreateWidget(label,xmLabelWidgetClass,
				      get_widget(parent),args,narg),desc);
  
  /* Create & Position the radiobutton */
  narg = 0;
  XtSetArg(args[narg],XmNleftAttachment,    XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNleftPosition,      but_left);          narg++;
  XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_POSITION); narg++;
  XtSetArg(args[narg],XmNrightPosition,     right);             narg++;
  if (top == parent) {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_FORM);     narg++;
  }
  else {
    XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET);   narg++;
    XtSetArg(args[narg],XmNtopWidget,       get_widget(top));   narg++;
  }
  XtSetArg(args[narg],XmNbottomAttachment,  XmATTACH_NONE);     narg++;
  XtSetArg(args[narg],XmNrightAttachment,   XmATTACH_NONE);     narg++;
  /*XtSetArg(args[narg],XmNvalue,             initial_value);     narg++;*/
  XtSetArg(args[narg],XmNradioBehavior,     True);              narg++;
  *wtextf = add_widget(XtCreateWidget(initial_value,xmRowColumnWidgetClass,
				      get_widget(parent),args,narg),desc);
  but = mk_but_array(initial_value,&nbut);
  for(i=0; (i<nbut); i++) {
    (void) add_widget(XtCreateWidget(but[i],xmPushButtonWidgetClass,
				     get_widget(*wtextf),NULL,0),but[i]);
  }
}

static void mk_buttons(int parent,int top,int nb,t_button bbb[])
{
  int   i,narg,nw,left,right,bw;
  float dx;

  nw      = 1;
  dx      = (100-(nb+1)*nw)/nb;
  for(i=0; (i<nb); i++) {
    left  = nw*(i+1)+i*dx;
    right = i*nw+(i+1)*(dx+nw);
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
    
    bw = add_widget(XtCreateWidget(bbb[i].label,xmPushButtonWidgetClass,
				   get_widget(parent),args,narg),
		    bbb[i].desc);
    XtAddCallback(get_widget(bw),XmNactivateCallback,bbb[i].cbfn,NULL);
  }
}

XmString xs_str_array_to_xmstr(char *header,int ndesc,char *desc[])
{
  int      i;
  XmString xmstr;
  char     *ptr,*cptr,*nlptr;
  
  if (ndesc <= 0)
    return empty_str;
  else {
    xmstr = char2xms(header);
    for(i=0; (i<ndesc); i++) {
      xmstr = XmStringConcat(xmstr,XmStringSeparatorCreate());
      ptr   = check_tty(desc[i]);
      cptr  = wrap_lines(ptr,70,0);
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

int mk_separator(int parent,int topw)
{
  int narg;
  
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
  
  return add_widget(XtCreateWidget("separator",xmSeparatorWidgetClass,
				   get_widget(parent),args,narg),NULL);
}

int mk_helplabel(int parent,int top)
{
  int narg;
  
  narg = 0;
  XtSetArg(args[narg],XmNtopAttachment,   XmATTACH_WIDGET); narg++;
  XtSetArg(args[narg],XmNtopWidget,       get_widget(top)); narg++;
  XtSetArg(args[narg],XmNleftAttachment,  XmATTACH_FORM);   narg++;
  XtSetArg(args[narg],XmNrightAttachment, XmATTACH_FORM);   narg++;
  XtSetArg(args[narg],XmNalignment,XmALIGNMENT_BEGINNING);  narg++;
  
  return add_widget(XtCreateWidget("  ",xmLabelWidgetClass,
				   get_widget(parent),args,narg),NULL);
}

void mk_filedlgs(int parent,int top,int nfile,t_filenm fnm[],int fnm_index[])
{
#define NWC 3
  WidgetClass wc[NWC];
  int    left[NWC]  = { 1,  8, 40 };
  int    right[NWC] = { 8, 40, 49 };
  char   **wname;
  int    i,j,ftp,narg,dx,www[NWC];
  Widget topw;
  char   *fn,dbuf[256];
  
  wc[0] = xmLabelWidgetClass;
  wc[1] = xmTextFieldWidgetClass;
  wc[2] = xmPushButtonWidgetClass;
  snew(wname,NWC);
  for(i=0; (i<nfile); i++) {
    ftp = fnm[i].ftp;
    dx  = (i % 2)*50;
    sprintf(dbuf,"%s [%s]",ftp2desc(ftp),fileopt(fnm[i].flag));
    if ((fn = strrchr(fnm[i].fn,'/')) == NULL)
      fn = fnm[i].fn;
    else
      fn++;
    wname[0]  = fnm[i].opt;
    wname[1]  = fn;
    wname[2]  = "Browse";
    
    for(j=0; (j<3); j++) {
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
      /*XtSetArg(args[narg],XmNrightAttachment, XmATTACH_POSITION);  narg++;
	XtSetArg(args[narg],XmNrightPosition,   dx+right[j]);        narg++;*/
      if (j == 1) {
	XtSetArg(args[narg],XmNvalue,         wname[j]);           narg++;
      }
      www[j] = add_widget(XtCreateWidget(wname[j],wc[j],
					 get_widget(parent),args,narg),dbuf);
    }
    fnm_index[i] = www[1];
    set_widget_ftp(www[2],ftp);
    set_widget_other(www[2],get_widget(www[1]));
    XtAddCallback(get_widget(www[2]),XmNactivateCallback,file_callback,NULL);
    if ((i % 2) == 1)
      topw = get_widget(www[1]);
  }
  
  if (nfile > 0) 
    (void) mk_separator(parent,nwidget()-2);
  sfree(wname);
}

int mk_pargs(int parent,int npargs,t_pargs pa[],int pa_index[])
{
  /**************************************************************
   *
   * Create all the children in two passes:
   * 1st the checkboxes (etBOOL variables)
   * 2nd the editable fields
   * Return the index of the separator (or the last created widget)
   *
   **************************************************************/
  int    icb,i,npa,dummy,separator,nbool,nedit,left,topw;
  real   bwidth;
  char   buf[132],descbuf[132];

  separator = nwidget()-1;  

  /* First round just count booleans and editors */
  nbool = nedit = 0;
  for(i=0; (i<npargs); i++) {
    if (!is_hidden(&(pa[i]))) {
      if (pa[i].type == etBOOL)
	nbool++;
      else
	nedit++;
    }
  }
  if (nbool > 6)
    nbool = 5;
  bwidth = (100.0/nbool);
    
  for(icb=0; (icb<2); icb++) {
    npa  = 0;
    topw = nwidget()-1;
    for(i=0; (i<npargs); i++) {
      if (!is_hidden(&(pa[i]))) {
	if ((icb == 0) && (pa[i].type == etBOOL)) {
	  left = 1+(npa % nbool)*bwidth;
	  if (debug)
	    fprintf(debug,"%s,%d: nbool %d, bwidth %g, left %d, topw %d\n",
		    __FILE__,__LINE__,nbool,bwidth,left,topw+nbool*(npa/nbool));
	  pa_index[i] = mk_toggle(parent,pa[i].option,
				  topw + nbool*(npa/nbool),
				  left,left+bwidth-2,
				  *(pa[i].u.b),pa[i].desc);
	  npa ++;
	}
	else if ((icb == 1) && (pa[i].type != etBOOL)) {
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
	  sprintf(descbuf,"%s (%s)",pa[i].desc,argtp[pa[i].type]);
	  if (debug)
	    fprintf(debug,"%s,%d: buf = %s\n",__FILE__,__LINE__,buf);
	  left = 1+(npa % 2)*50;
	  if (pa[i].type == etSTR) {
	    mk_enumerated(parent,topw+4*(npa/2),
			  left,left+47,
			  &dummy,&(pa_index[i]),
			  pa[i].option,buf,descbuf); 
	  }
	  else {
	    mk_editor(parent,topw+4*(npa/2),
		      left,left+47,
		      &dummy,&(pa_index[i]),
		      pa[i].option,buf,descbuf); 
	  }
	  npa ++;
	}
      }
    }
    if (npa > 0)
      separator = mk_separator(parent,-1);
  }
  return separator;
}

void mk_help(Widget parent,int ndesc,char *desc[],int nbugs,char *bugs[])
{
  Widget   label,dialog,ok;
  XmString xmstr;
  int      narg;

  /* Convert the strings */  
  xmstr = xs_str_array_to_xmstr("DESCRIPTION:",ndesc,desc);
  if (nbugs > 0) {
    xmstr = XmStringConcat(xmstr,XmStringSeparatorCreate());
    xmstr = XmStringConcat(xmstr,XmStringSeparatorCreate());
    xmstr = XmStringConcat(xmstr,xs_str_array_to_xmstr("DIAGNOSTICS:",
						       nbugs,bugs));
  }
  if (debug)
    fprintf(debug,"%s,%d: Created xmstr\n",__FILE__,__LINE__);
  narg  = 0;
  XtSetArg(args[narg],XmNautoUnmanage,FALSE);  narg++;
  XtSetArg(args[narg],XmNmessageString,xmstr); narg++;
  /*XtSetArg(args[narg],XmNdialogType,XmDIALOG_INFORMATION); narg++;*/
  /*XtSetArg(args[narg],XmNdialogTitle,char2xms("Gromacs Help")); narg++;*/
  dialog = XmCreateMessageDialog(parent,"Help",args,narg);
  helpw  = add_widget(dialog,NULL);
  XtUnmanageChild(XmMessageBoxGetChild(dialog,XmDIALOG_CANCEL_BUTTON));
  XtUnmanageChild(XmMessageBoxGetChild(dialog,XmDIALOG_HELP_BUTTON));
  if (debug)
    fprintf(debug,"%s,%d: Unmanaged children\n",__FILE__,__LINE__);
  ok = XmMessageBoxGetChild(dialog,XmDIALOG_OK_BUTTON);
  XtAddCallback(ok,XmNactivateCallback,help_ok_callback,NULL);
  
  label = XmMessageBoxGetChild(dialog,XmDIALOG_MESSAGE_LABEL);
  narg  = 0;
  XtSetArg(args[narg],XmNalignment,XmALIGNMENT_BEGINNING); narg++;
  XtSetValues(label,args,narg);
  if (debug)
    fprintf(debug,"%s,%d: Aligned help\n",__FILE__,__LINE__);
}

void mk_fdlg(Widget base)
{
  int narg;
  
  narg  = 0;
  fdlgw = add_widget(XmCreateFileSelectionDialog(base,"GMX File Selector",
						 args,narg),NULL);
  XtAddCallback(get_widget(fdlgw),XmNokCallback,file_ok_callback,NULL);
  XtAddCallback(XmFileSelectionBoxGetChild(get_widget(fdlgw),
					   XmDIALOG_CANCEL_BUTTON),
		XmNactivateCallback,file_cancel_callback,NULL);
  XtUnmanageChild(XmFileSelectionBoxGetChild(get_widget(fdlgw),
					     XmDIALOG_HELP_BUTTON));
  /*XtUnmanageChild(XmFileSelectionBoxGetChild(get_widget(fdlgw),
    XmDIALOG_APPLY_BUTTON));*/
}

static void mk_gui(Widget gmxBase,
		   int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
		   int ndesc,char *desc[],int nbugs,char *bugs[])
{
  t_button bbb[] = {
    { "OK",     "Press OK to accept current settings and continue",
      &ok_callback },
    { "Cancel", "Press Cancel to quit the program",
      &cancel_callback },
    { "Help",   "Press Help for more information about the program",
      &help_callback }
  };
#define NBUT asize(bbb)
  int  sep1,sep2,widg0;
  int  i,narg;
  
  /* Create the help window */
  mk_help(gmxBase,ndesc,desc,nbugs,bugs);

  /* Create the file selector box */
  mk_fdlg(gmxBase);
  
  /* Create the dialog box! */
  gmxDialog = add_widget(XmCreateForm(gmxBase,"gmxDialog",NULL,0),NULL);
  XtManageChild(get_widget(gmxDialog));

  widg0 = nwidget();
  
  /* Create buttons for file dialogboxes */
  if (nfile > 0) {
    snew(fnm_index,nfile);
    mk_filedlgs(gmxDialog,gmxDialog,nfile,fnm,fnm_index);
  }
  
  /* Create the checkboxes and editable fields */
  if (npargs > 0) {
    snew(pa_index,npargs);
    sep1 = mk_pargs(gmxDialog,npargs,pa,pa_index);
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
  mk_desc_callbacks();
  
  /* Give the children a father or mother */
  for(i=widg0; (i<nwidget()); i++)
    XtManageChild(get_widget(i));
}

/******************************************************************
 *
 *      C A L L B A C K       R O U T I N E S
 *
 ******************************************************************/

static void cancel_callback(Widget w,caddr_t client_data,caddr_t call_data)
{
  printf("Maybe next time...\n");
  thanx(stdout);
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

static void enter_callback(Widget www,int *which,caddr_t call_data)
{
  int narg;

  if (descw != -1) {  
    narg=0;
    XtSetArg(args[narg],XmNlabelString,get_widget_desc(www));  narg++;
    XtSetValues(get_widget(descw),args,narg);
  }
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
	/*XtSetArg(args[narg],XmNpattern,xms);        narg++;*/
	XtSetArg(args[narg],XmNlistUpdated, False); narg++;
	
	XtSetValues(fdlg,args,narg);
	XmStringFree(xms);
      }
      XtManageChild(fdlg);
      bFdlgUp    = TRUE;
      FdlgCaller = get_widget_other(get_windex(www));
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

static void leave_callback(Widget www,int *which,caddr_t call_data)
{
  int narg;

  if (descw != -1) {  
    narg=0;
    XtSetArg(args[narg],XmNlabelString,empty_str); narg++;
    XtSetValues(get_widget(descw),args,narg);
  }
}

void MyMainLoop(XtAppContext appcontext,Widget gmxBase,
		int nfile,t_filenm fnm[],int npargs,t_pargs pa[])
{
  int      i,narg;
  Widget   www;
  XEvent   event;
  XmString xms;
  Boolean  bbb;
  char     *fn,buf[256],*ptr;
  double   ddd,dx,dy,dz;
  int      iii;
  
  while (!bDone) {
    /*XtAppNextEvent(appcontext,&event);*/
    XtNextEvent(&event);
    XtDispatchEvent(&event);
  }
  /* Extract all the information from the X widgets */
  for(i=0; (i<nfile); i++) {
    www  = get_widget(fnm_index[i]);
    narg = 0;
    /*xms  = char2xms("");
      XtSetArg(args[narg], XmNlabelString, &xms); narg++;
      XtGetValues(www,args,narg);
      fn = xms2char(xms);
    */
    XtSetArg(args[narg], XmNvalue, &fn); narg++;
    XtGetValues(www,args,narg);
    sprintf(buf,"%s%s",get_widget_dir(fnm_index[i]),fn);
    XtFree(fn);
    sfree(fnm[i].fn);
    fnm[i].fn = strdup(buf);
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
	narg = 0;
	XtSetArg(args[narg],XmNvalue, &ptr);            narg++;
	XtGetValues(get_widget(pa_index[i]),args,narg);
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
  /*XtCloseDisplay();*/
}

void gmx_gui(int *argc,char *argv[],
	     int nfile,t_filenm fnm[],int npargs,t_pargs pa[],
	     int ndesc,char *desc[],int nbugs,char *bugs[])
{
  Widget       gmxBase;
  XtAppContext appcontext;
  
  /* Initialize toolkit and parse command line options. */
  gmxBase = XtInitialize("GROMACS","GROMACS",NULL,0,argc,argv);
  /*gmxBase = XtVaAppInitialize(&appcontext,"GROMACS",NULL,0,argc,argv,NULL);*/
  mk_gui(gmxBase,nfile,fnm,npargs,pa,ndesc,desc,nbugs,bugs);
  XtRealizeWidget(gmxBase);
  MyMainLoop(appcontext,gmxBase,nfile,fnm,npargs,pa);
}


