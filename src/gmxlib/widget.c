#include <stdio.h>
#include "widget.h"
#include "smalloc.h"
#include "fatal.h"

typedef struct {
  Widget   w,other;
  bool     bDesc;
  XmString desc;
  char     *directory;
  int      ftp;
} t_widget;

static   t_widget *w=NULL;
static   int      nwindex=0,maxwindex=0;
XmString empty_str;

int nwidget()
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
      w[i].bDesc = FALSE;
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

static void widget_range_check(windex win)
{
  if (!((win>=0) && (win<nwindex)))
    fatal_error(0,"Widget index %d out of range, nwindex = %d",win,nwindex);
}

Widget get_widget(windex win)
{
  widget_range_check(win);
  
  return w[win].w;
}

windex get_windex(Widget www)
{
  int i;
  
  for(i=0; (i<nwindex); i++)
    if (w[i].w == www)
      return i;
  fatal_error(0,"No such widget %x\n",www);
  
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
  widget_range_check(www);
  
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
  widget_range_check(win);

  return w[win].directory ? w[win].directory : "";
}

void set_widget_ftp(windex win,int ftp)
{
  widget_range_check(win);

  w[win].ftp = ftp;
}

void set_widget_dir(Widget www,XmString label)
{
  Arg  args[4];
  int  i,narg;
  char *ptr,*clab,tmp;
  XmString xms;

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
    }
  }
}

Widget get_widget_other(windex win)
{
  widget_range_check(win);
  if (w[win].other == 0)
    fatal_error(0,"Calling the wrong window: %d, no other widget\n",win);
  
  return w[win].other;
}
   
void set_widget_other(windex win,Widget www)
{
  widget_range_check(win);
  w[win].other = www;
}
   
XmString char2xms(char *ptr)
{
  return XmStringCreate(ptr,XmSTRING_DEFAULT_CHARSET);
}

char *xms2char(XmString xms)
{
  /* This is NOT fool proof */
  char *str;

  /* set str to point to the text */
  XmStringGetLtoR(xms,XmSTRING_DEFAULT_CHARSET,&str);

  return str;
}

