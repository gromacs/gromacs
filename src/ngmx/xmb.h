/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

#ifndef _xmb_h
#define _xmb_h

static char *SRCID_xmb_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) xmb.h 1.3 9/29/92"
#endif /* HAVE_IDENT */

#include <x11.h>

#define MB_OK              1
#define MB_CANCEL          (1<<1)
#define MB_OKCANCEL        (MB_OK | MB_CANCEL)
#define MB_YES             (1<<2)
#define MB_NO              (1<<3)
#define MB_YESNO           (MB_YES | MB_NO)
#define MB_ICONSTOP        (1<<16)
#define MB_ICONINFORMATION (1<<17)
#define MB_ICONEXCLAMATION (1<<18)
#define MB_ICONGMX         (1<<19)
#define MB_SYSTEMMODAL     (1<<20)
#define MB_APPLMODAL       (1<<21)
#define MB_DONTSHOW        (1<<22)

t_dlg *MessageBox(t_x11 *x11, Window Parent, char *title,
		  int nlines, char *lines[], ulong Flags,
		  DlgCallback *cb, void *data);

#endif	/* _xmb_h */
