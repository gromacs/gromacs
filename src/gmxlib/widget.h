#ifndef _widget_h
#define _widget_h

#include "typedefs.h"
#include <Xm/Xm.h>

typedef int windex;

extern int      nwidget(void);
extern windex   add_widget(Widget new_widget,char *desc);
extern Widget   get_widget(windex win);
extern windex   get_windex(Widget www);
extern XmString get_widget_desc(Widget www);
extern bool     have_windex_desc(windex www);
extern char     *get_widget_directory(windex win);
extern int      get_widget_ftp(Widget www);
extern void     set_widget_ftp(windex win,int ftp);
extern void     set_widget_dir(Widget www,XmString label);
extern void     mk_desc_callbacks(void);
extern XmString char2xms(char *ptr);
extern char     *xms2char(XmString xms);

#endif
