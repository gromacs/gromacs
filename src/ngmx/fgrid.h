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
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _fgrid_h
#define _fgrid_h

static char *SRCID_fgrid_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) fgrid.h 1.3 9/29/92"
#endif /* HAVE_IDENT */

#include <typedefs.h>
#include <xdlg.h>

typedef struct {
  edlgitem edlg;
  bool bDef;
  int  nname;
  char **name;
  char *set,*get,*def,*help;
} t_fitem;

typedef struct {
  char *name;
  int x,y,w,h;
  int nfitem;
  t_fitem **fitem;
} t_fgroup;

typedef struct {
  int x,y,w,h;
  t_fitem *fitem;
} t_fsimple;

typedef struct {
  int       w,h;
  int       nfgroup;
  t_fgroup  **fgroup;
  int       nfsimple;
  t_fsimple **fsimple;
} t_fgrid;

typedef enum {
  eGRIDEXP, eACCOEXP, eACCCEXP, eGRPEXP, eITEMEXP, eSAMEPOINT, 
  eTOOWIDE, eTOOHIGH, eQUOTE,   eNOVALS 
  } eDLGERR;

extern void ReadDlgErr(char *infile, eDLGERR err, char *s);

extern t_fgrid *FGridFromFile(char *infile);

extern void DoneFGrid(t_fgrid *fgrid);

extern void DumpFGrid(t_fgrid *fgrid);

extern void ReadQuoteString(char *infile, FILE *in, char *buf);

#endif	/* _fgrid_h */
