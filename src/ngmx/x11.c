/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
#include "sysstuff.h"
#include "string.h"
#include "assert.h"
#include "smalloc.h"
#include "Xstuff.h"
#include "typedefs.h"
#include "x11.h"

/* These colours will be mapped to black on a monochrome screen */
ulong BLACK,BLUE,GREEN,CYAN,RED,BROWN,GREY,DARKGREY;

/* These colours will be mapped to white on a monochrome screen */
ulong LIGHTBLUE,LIGHTGREEN,LIGHTGREY,LIGHTCYAN,LIGHTRED,VIOLET,YELLOW,WHITE;

static XFontStruct *XLQF(FILE *err, Display *disp, char *name)
{
  XFontStruct *font=XLoadQueryFont(disp,name);
#ifdef DEBUG
  if (font != NULL) 
    fprintf(err, "Loaded font %s\n",name);
#endif
  return font;
}
  
static XFontStruct *GetFont(FILE *err, Display *disp, char *name)
{
  static char *fontnames[] = { 
    "sansb12","8x13bold","8x13",
    "9x15","6x13","fixed" 
    };
#define MAXNAMES (sizeof(fontnames)/sizeof(fontnames[0]))
  int i;
  XFontStruct *font;
  int count;
  char **fontlist;
  bool bFont = FALSE;

  if (name)
    bFont=(bool) (font=XLQF(err,disp,name));
  else
    font=NULL;
    
  for (i=0; (!bFont && (i<MAXNAMES)); i++) 
    bFont=(bool) (font=XLQF(err,disp,fontnames[i]));

  if (!bFont) {
    fontlist=XListFonts(disp,"?",1,&count);
    if (count!=0) 
      bFont=(bool) (font=XLQF(err,disp,fontlist[0]));
  }
  if (!bFont) 
    fprintf (err, "Cannot load any suitable font\n");
  return font;
}

static GC GetGC(Display *disp, XFontStruct *font)
{
  XGCValues     values;

  values.font = font->fid;
  values.foreground = WhitePixel(disp,DefaultScreen(disp));

  return XCreateGC(disp,DefaultRootWindow(disp),GCForeground|GCFont,&values);
}

void GetNamedColor(t_x11 *x11,char *name,ulong *col)
{
  /* If name is found than col set to that colour else col is unchanged */
  XColor exact,clr;

  if (XAllocNamedColor(x11->disp,x11->cmap,name,&clr,&exact))
    *col=clr.pixel;
  else
    fprintf(x11->console,"No colour %s\n",name);
}

static t_wlist *GetWList(t_x11 *x11, Window w)
{
  t_wlist *curs;

  curs=x11->wlist;
  while (curs && (curs->w != w))
    curs=curs->next;

  return curs;
}

typedef struct {
  Window w;
  bool   b;
} t_peek;

static Bool TestEvent(Display *disp,XEvent *event,char *arg)
{
  t_peek *tp;

  fprintf(stderr,"TestEvent\n");
  tp=(t_peek *)arg;
  if ((event->xany.window==tp->w) && (event->type==ConfigureNotify)) {
    tp->b=TRUE;
    return True;
  }
  return False;
}

static void MainLoop(t_x11 *x11)
{
  bool    bReturn;
  XEvent  event;
  t_peek  tp;
  t_wlist *curs;
  Window  w;

  for (bReturn=FALSE; (!bReturn); ) {
    if (x11->wlist) {
      XNextEvent(x11->disp,&event);
      w=event.xany.window;
      tp.w=w;
      tp.b=FALSE;
      curs=GetWList(x11,w);
      if (!curs)
	bReturn=TRUE;
      if (!bReturn) {
	switch (event.type) {
	case Expose:
	  /* Filter out expose events with non-zero count field */
	  if (event.xexpose.count != 0)
	    curs=NULL;
	  break;
	case ConfigureNotify:
	  /* Check if more are coming... 
	  if (XCheckTypedWindowEvent(x11->disp,w,ConfigureNotify,&config))
	    curs=NULL; */
	  break;
	default:
	  break;
	}
	if (curs)
	  bReturn=(*curs->cb)(x11,&event,w,curs->data);
      }
    }
  }
}

static void RegisterCallback(t_x11 *x11,Window w,Window Parent,
			     CallBack cb, void *data)
{
  t_wlist *curs,*new;

  snew(new,1);
  new->w=w;
  new->Parent=Parent;
  new->cb=cb;
  new->mask=0;
  new->data=data;
  new->next=NULL;

  if (x11->wlist) {
    curs=x11->wlist;
    while(curs->next)
      curs=curs->next;
    curs->next=new;
  }
  else
    x11->wlist=new;
}

static void UnRegisterCallback(t_x11 *x11, Window w)
{
  t_wlist *curs;

  curs=x11->wlist;
  if (curs) {
    if (curs->w==w) {
      x11->wlist=curs->next;
      sfree(curs);
    }
    else {
      while (curs->next && (curs->next->w != w))
	curs=curs->next;
      if (curs->next) {
	t_wlist *tmp=curs->next;

	curs->next=curs->next->next;
	sfree(tmp);
      }
    }
  }
}

static void SetInputMask(t_x11 *x11, Window w, ulong mask)
{
  t_wlist *curs;

  curs=GetWList(x11,w);
  if (curs) {
    curs->mask=mask;
    XSelectInput(x11->disp,w,mask);
  }
  else 
    fprintf(x11->console,"No such window (%d)\n",(int)w);
}

static ulong GetInputMask(t_x11 *x11, Window w)
{
  t_wlist *curs;

  curs=GetWList(x11,w);
  if (curs)
    return curs->mask;
  else
    return 0;
}

static void CleanUp(t_x11 *x11) 
{
  t_wlist *curs;
  
  curs=x11->wlist;
  while (curs) {
    x11->wlist=curs->next;
    XDestroyWindow(x11->disp,curs->w);
    sfree(curs);
    curs=x11->wlist;
  }
  XCloseDisplay(x11->disp);
}

static void Xrm(int *argc, char *argv[])
{
  /*
  static XrmOptionDescRec opTable[] = {
    {"-background",   "*background",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-bd",           "*borderColor",   
       XrmoptionSepArg, (caddr_t) NULL},
    {"-bg",           "*background",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-borderwidth",  "*TopLevelShell.borderwidth",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-bordercolor",   "*borderColor",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-bw",            "*TopLevelShell.borderColor",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-display",       ".display",    
       XrmoptionSepArg, (caddr_t) NULL},
    {"-fg",            "*foreground",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-fn",            "*font",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-font",          "*font",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-foreground",    "*foreground",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-geometry",      ".TopLevelShell.geometry",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-iconic",        ".TopLevelShell.iconic",
       XrmoptionNoArg,  (caddr_t) "on"},
    {"-name",          ".name",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-reverse",       "*reverseVideo",
       XrmoptionNoArg,  (caddr_t) "on"},
    {"-rv",            "*reverseVideo",
       XrmoptionNoArg,  (caddr_t) "on"},
    {"-synchronous",   ".synchronous",
       XrmoptionNoArg,  (caddr_t) "on"},
    {"-title",         ".TopLevelShell.title",
       XrmoptionSepArg, (caddr_t) NULL},
    {"-xrm",            NULL,
       XrmoptionSepArg, (caddr_t) NULL},
  };
#define TABLELENGTH (sizeof(opTable)/sizeof(opTable[0]))
  XrmInitialize();*/
}

static void Flush(t_x11 *x11)
{
  fflush(x11->console);
}

t_x11 *GetX11(int *argc, char *argv[])
{
  static char *v_name[] = {
    "DirectColor","TrueColor", "PseudoColor",
    "StaticColor","GrayScale", "StaticGray"
    };
  static int v_class[] = {
    DirectColor,TrueColor, PseudoColor,
    StaticColor,GrayScale, StaticGray
    };
#define NCLASS (sizeof(v_class)/sizeof(v_class[0]))

  XVisualInfo v_info;
  t_x11       *x11;
  int         ARGC;
  char        **ARGV;
  char        *display;
  char        *fontname;
  char        *title,*FG=NULL,*BG=NULL;
  bool        bVerbose=FALSE;
  int         i;

  title=strdup(argv[0]);

  /* First check environment */
  fontname=getenv("GMXFONT");
  display=getenv("DISPLAY");

  snew(ARGV,*argc);
  ARGC=1;
  for(i=1; (i < *argc); i++) {
    if (argv[i][0]=='-') {
      if (strlen(argv[i]) > 1) {
	if ((*argc)>i+1)
	  switch(argv[i][1]) {
	  case 'b':
	    BG=argv[++i];
	    break;
	  case 'd':
	    display=argv[++i];
	    break;
	  case 'f':
	    switch(argv[i][2]) {
	    case 'o':
	      fontname=argv[++i];
	      break;
	    case 'g':
	      FG=argv[++i];
	      break;
	    }
	    break;
	  case 't':
	    sfree(title);
	    title=strdup(argv[++i]);
	    break;
	  case 'v':
	    bVerbose=TRUE;
	    break;
	  default:
	    ARGV[ARGC++]=argv[i];
	    break;
	  }
      }
    }
    else 
      ARGV[ARGC++]=argv[i];
  }
  for (i=1; (i<ARGC); i++)
    argv[i]=ARGV[i];
  *argc=ARGC;
  argv[ARGC]=NULL;

  snew(x11,1);
  x11->dispname=display;
  if (bVerbose)
    x11->console=stderr;
  else
    if ((x11->console=fopen("/dev/null","w"))== NULL)
      x11->console=stderr;

  if ((x11->disp=XOpenDisplay(display))==NULL) {
    if (bVerbose)
      fprintf(x11->console,"Display %s invalid\n",display);
    return NULL;
  }
  
  if ((x11->font=GetFont(x11->console,x11->disp,fontname))==NULL)
    return NULL;
  if ((x11->gc=GetGC(x11->disp,x11->font))==NULL)
    return NULL;

  x11->root=DefaultRootWindow(x11->disp);
  x11->screen=DefaultScreen(x11->disp);
  x11->depth=DefaultDepth(x11->disp,x11->screen);
  x11->cmap=DefaultColormap(x11->disp,x11->screen);

  /* These colours will be mapped to black on a monochrome screen */
  x11->fg=BLACK=BLUE=GREEN=CYAN=RED=BROWN=GREY=DARKGREY=
    BlackPixel(x11->disp,x11->screen);

  /* These colours will be mapped to white on a monochrome screen */
  x11->bg=
    LIGHTBLUE=LIGHTGREY=LIGHTGREEN=LIGHTCYAN=LIGHTRED=VIOLET=YELLOW=WHITE=
      WhitePixel(x11->disp,x11->screen);

  if (x11->depth > 1) {
    /* Not B & W, Look what kind of screen we've got... */
    for (i=0; (i < NCLASS); i++)
      if (!XMatchVisualInfo(x11->disp,x11->screen,x11->depth,
			    v_class[i],&v_info))
	break;
    if ((i==4) || (i==5)) 
      fprintf(x11->console,"Greyscale screen, using B & W only\n");
    else {
      /* We have real color! */
      fprintf(x11->console,"%s screen with depth %d.\n",
	     (i==NCLASS)?"Unknown":v_name[i],x11->depth);
      GetNamedColor(x11,"midnight blue",&BLUE);
      GetNamedColor(x11,"DarkGreen",&GREEN);
      GetNamedColor(x11,"SeaGreen",&CYAN);
      GetNamedColor(x11,"red4",&RED);
      GetNamedColor(x11,"Gray",&GREY);
      GetNamedColor(x11,"Gray",&DARKGREY);
      GetNamedColor(x11,"LightGray",&LIGHTGREY);
      GetNamedColor(x11,"green",&LIGHTGREEN);
      GetNamedColor(x11,"cyan",&LIGHTCYAN);
      GetNamedColor(x11,"tomato1",&LIGHTRED);
      GetNamedColor(x11,"violet",&VIOLET);
      GetNamedColor(x11,"yellow",&YELLOW);
      GetNamedColor(x11,"brown",&BROWN);
      GetNamedColor(x11,"CornFlowerBlue",&LIGHTBLUE);
    }
  }
  else
    fprintf(x11->console,"Monochrome screen.\n");

  /* We should use Xrm here... */
  if (FG)
    GetNamedColor(x11,FG,&(x11->fg));
  else
    x11->fg=BLACK;
  if (BG)
    GetNamedColor(x11,BG,&(x11->bg));
  else
    x11->bg=LIGHTGREY;
  x11->title=strdup(title);
  sfree(title);
  x11->wlist=NULL;
  x11->GetNamedColor=&GetNamedColor;
  x11->MainLoop=&MainLoop;
  x11->RegisterCallback=&RegisterCallback;
  x11->UnRegisterCallback=&UnRegisterCallback;
  x11->SetInputMask=&SetInputMask;
  x11->GetInputMask=&GetInputMask;
  x11->CleanUp=&CleanUp;
  x11->Flush=&Flush;

  x11->Flush(x11);

  return x11;
}

