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

#ifndef _xdlghi_h
#define _xdlghi_h

static char *SRCID_xdlghi_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) xdlghi.h 1.2 9/29/92"
#endif /* HAVE_IDENT */

#include <stdarg.h>
#include "Xstuff.h"
#include "x11.h"
#include "xdlg.h"

typedef struct {
  int       nitem;
  int       w,h;
  t_dlgitem **list;
} t_dlgitemlist;

extern t_dlgitem **CreateRadioButtonGroup(t_x11 *x11, char *szTitle, 
					  t_id GroupID, int nrb, t_id rb[],
					  int nSelect,
					  char *szRB[], int x0,int y0);
/* This routine creates a radio button group at the
 * specified position. The return values is a pointer to an
 * array of dlgitems, the array has length (nrb+1) with the +1
 * because of the groupbox.
 * nSelect is the ordinal of the selected button.
 */

extern t_dlgitem **CreateDlgitemGroup(t_x11 *x11, char *szTitle, 
				      t_id GroupID, int x0, int y0,
				      int nitem, ...);
/* This routine creates a dlgitem group at the
 * specified position. The return values is a pointer to an
 * array of dlgitems, the array has length (nitem+1) with the +1
 * because of the groupbox.
 */

extern t_dlg *ReadDlg(t_x11 *x11,Window Parent, char *title,
		      unsigned long fg, unsigned long bg, char *infile, 
		      int x0, int y0, bool bAutoPosition,bool bUseMon,
		      DlgCallback *cb,void *data);
/* Read a dialog box from a template file */

#endif	/* _xdlghi_h */
