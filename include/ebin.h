/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _ebin_h
#define _ebin_h

static char *SRCID_ebin_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) ebin.h 1.10 5/2/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"
	
typedef struct {
  int      nener;
  char     **enm;
  t_energy *e;
} t_ebin;

enum { eprNORMAL, eprAVER, eprRMS, eprNR };

extern t_ebin *mk_ebin(void);
/* Create an energy bin */

extern int get_ebin_space(t_ebin *eb,int nener,char *enm[]);
/* Create space in the energy bin and register names.
 * The enm array must be static, because the contents are not copied,
 * but only the pointers.
 * Function returns an index number that must be used in subsequent
 * calls to add_ebin.
 */

extern void add_ebin(t_ebin *eb,int index,int nener,real ener[],int step);
/* Add nener reals (eg. energies, box-lengths, pressures) to the
 * energy bin at position index. 
 */
 
extern void pr_ebin(FILE *fp,t_ebin *eb,int index,int nener,int nperline,
		    int prmode,int tsteps,bool bPrHead);
/* Print the contents of the energy bin. If nener = -1 ALL energies from
 * index to the end will be printed. We will print nperline entries on a text
 * line (advisory <= 5). prmode may be any of the above listed enum values.
 * tsteps is used only when eprAVER or eprRMS is set.
 * If bPrHead than the header is printed.
 */



#endif	/* _ebin_h */
