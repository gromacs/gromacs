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
 * GROup of MAchos and Cynical Suckers
 */

#ifndef	_xdlg_h
#define	_xdlg_h

#ifdef HAVE_IDENT
#ident	"@(#) xdlg.h 1.3 9/29/92"
#endif /* HAVE_IDENT */

#include <xdlgitem.h>

#define DLG_SHOW         (1<<0)
#define DLG_HIDE         (1<<1)
#define DLG_SHOWANDHIDE  (DLG_SHOW | DLG_HIDE)
#define DLG_SYSTEMMODAL  (1<<2)
#define DLG_APPLMODAL    (1<<3)
#define DLG_HIDEONBUTTON (1<<4)
#define DLG_FREEONBUTTON (1<<5)

enum { DLG_SET, DLG_EXIT };

typedef void DlgCallback(t_x11 *x11,int dlg_mess,int item_id,
			 char *set,void *data);
/* User function that can be called by the dialog box. All setting of
 * check-boxes and radio-buttons etc., is done by the dialog manager,
 * the user can let himself be informed about mouse activity also.
 */

typedef struct {
  t_x11         *x11;           /* All about X 				*/
  t_windata    	win;     	/* The position and size of the window 	*/
  char          *title;         /* Window name 				*/
  Window	wDad;		/* The parent window 			*/
  uint          xmax,ymax;      /* Dimensions of parent window 		*/
  ulong         flags;          /* Flags for display 			*/
  ulong         fg,bg;          /* The colours 				*/
  bool          bPop;           /* Should we pop the mouse back 	*/
  bool          bGrab;          /* Have we grabbed the mouse ? 		*/
  int         	nitem;		/* The number of items 			*/
  t_dlgitem   	**dlgitem;	/* The array of item pointers 		*/
  DlgCallback   *cb;		/* User call back function		*/
  void          *data;		/* User data				*/
} t_dlg;

/*****************************
 *
 * Routine to create the DLG structure, returns NULL on failure
 * cb and data may be NULL.
 *
 ****************************/
t_dlg *CreateDlg(t_x11 *x11,Window Parent,char *title,
		 int x0,int y0,int w,int h,int bw,ulong fg,ulong bg,
		 DlgCallback *cb,void *data);

/*****************************
 *
 * Routine to add an item to the dialog box
 * The pointer to the item is copied to the dlg struct,
 * the item itself may not be freed until the dlg is done with
 *
 ****************************/
void AddDlgItem(t_dlg *dlg,t_dlgitem *new);

void AddDlgItems(t_dlg *dlg,int nitem,t_dlgitem *new[]);

/*****************************
 *
 * Routines to manipulate items on a dialog box
 * They return TRUE on succes, FALSE otherwise
 * FALSE will mean most of the time, that item id was not found
 *
 ****************************/
bool QueryDlgItemSize(t_dlg *dlg,t_id id,int *w,int *h);

bool QueryDlgItemPos(t_dlg *dlg,t_id id,int *x0,int *y0);

int QueryDlgItemX(t_dlg *dlg, t_id id);

int QueryDlgItemY(t_dlg *dlg, t_id id);

int QueryDlgItemW(t_dlg *dlg, t_id id);

int QueryDlgItemH(t_dlg *dlg, t_id id);

bool SetDlgItemSize(t_dlg *dlg,t_id id,int w,int h);

bool SetDlgItemPos(t_dlg *dlg,t_id id,int x0,int y0);

void SetDlgSize(t_dlg *dlg,int w,int h, bool bAutoPosition);

/*****************************
 *
 * Routines to extract information from the dlg proc
 * after dlg is exec'ed
 *
 ****************************/
bool IsCBChecked(t_dlg *dlg,t_id id);

t_id RBSelected(t_dlg *dlg,int gid);

int  EditTextLen(t_dlg *dlg,t_id id);

char *EditText(t_dlg *dlg,t_id id);

/*****************************
 *
 * Routines to do internal things
 *
 ****************************/
t_dlgitem *FindWin(t_dlg *dlg, Window win);

t_dlgitem *FindItem(t_dlg *dlg, t_id id);

void HelpDlg(t_dlg *dlg);

void HelpNow(t_dlg *dlg, t_dlgitem *dlgitem);

void NoHelp(t_dlg *dlg);

/*****************************
 *
 * Exececute the dialog box procedure
 * Returns when a button is pushed.
 * return value is the ID of the button
 *
 ****************************/
void ShowDlg(t_dlg *dlg);

void HideDlg(t_dlg *dlg);

void FreeDlgItem(t_dlg *dlg, t_id id);

void FreeDlg(t_dlg *dlg);

#endif	/* _xdlg_h */
