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

#ifndef _dialogs_h
#define _dialogs_h

static char *SRCID_dialogs_h = "$Id$";

#include "xdlg.h"
#include "pulldown.h"
#include "manager.h"
#include "logo.h"

typedef struct {
  bool  bMapped;
  t_dlg *dlg;
} t_dialogs;

typedef enum { edExport, edBonds, edFilter, edNR } eDialogs;
	 
typedef enum {
  emQuit, emHelp, emAbout, emNotImplemented, emNR
} eMBoxes;

typedef enum {
  eExpGromos, eExpPDB, eExpNR
} eExport;

typedef struct {
  char         confout[256];	/* Export file			*/
  int          ExpMode;		/* Export mode			*/
  t_dlg        **dlgs;	        /* Temporary storage for dlgs	*/
  int          which_mb;        /* Which mb is visible          */
  t_dlg        **mboxes;        /* id for message boxes         */
  t_filter     *filter; 	/* Filter for visibility etc.	*/
  t_windata    *wd;		/* The main window		*/
  t_pulldown   *pd;		/* The pull-down menu		*/
  t_manager    *man;		/* The manager			*/
  /*t_statrec    *sr;*/		/* The statistics dlg		*/
  t_logo       *logo;           /* The gromacs logo             */
} t_gmx;

enum { 
  IDNEW,IDOPEN,IDOPENED,IDCLOSE,IDIMPORT,IDEXPORT,IDDOEXPORT,IDQUIT,IDTERM,
  IDEDITTOP,IDEDITCOORDS,IDEDITPARAMS,
  IDGROMPP,IDRUNMD,IDDOGROMPP,IDGSTAT,IDDOGSTAT,IDDORUNMD,
  IDFILTER,IDDOFILTER,
  IDANIMATE,IDSHOWBOX,IDRMPBC,IDHYDROGEN,IDLABELSOFF,IDRESETVIEW,IDPHOTO,
  IDDUMPWIN,IDDODUMP,
  IDBONDOPTS,IDTHIN,IDFAT,IDVERYFAT,IDBALLS,
  IDNOBOX,IDRECTBOX,IDTRIBOX,IDTOBOX,
  IDBOND,IDANGLE,IDDIH,IDRMS,IDRDF,IDENERGIES,IDCORR,
  IDHELP,IDABOUT,
  
  /* Last line specifies how many IDs there are */
  IDMENUNR
  };

extern void run_grompp(t_gmx *gmx);

extern void run_mdrun(t_gmx *gmx);

extern void write_gmx(t_x11 *x11,t_gmx *gmx,int mess);

/*extern void run_sr(t_statrec *sr);

extern t_statrec *init_sr();*/

extern void init_dlgs(t_x11 *x11,t_gmx *gmx);

extern void show_mb(t_gmx *gmx,int mb);

extern void done_dlgs(t_gmx *gmx);

extern void edit_file(char *fn);

extern t_filter *init_filter(t_atoms *atoms, char *fn);

extern t_dlg *select_filter(t_x11 *x11,t_gmx *gmx);

#endif
