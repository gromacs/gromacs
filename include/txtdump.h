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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _txtdump_h
#define _txtdump_h

static char *SRCID_txtdump_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) txtdump.h 1.20 12/16/92"
#endif /* HAVE_IDENT */

#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"

#define LINE_WIDTH	80
#define RMARGIN		10
#define	USE_WIDTH	((LINE_WIDTH)-(RMARGIN))
#define INDENT		3

extern char *atomname(t_atoms *a,int i);
/* Return pointer to a buffer which holds the atomname in the
 * form resname resnr atomname. Pointer can be freed afterwards.
 */
extern void pr_shownumbers(bool bShow);
extern int available(FILE *fp,void *p,char *title);
extern int pr_indent(FILE *fp,int n);
extern int pr_title(FILE *fp,int indent,char *title);
extern int pr_title_n(FILE *fp,int indent,char *title,int n);
extern int pr_title_nxn(FILE *fp,int indent,char *title,int n1,int n2);
extern void pr_ivec(FILE *fp,int indent,char *title,int vec[],int n);
extern void pr_ivecs(FILE *fp,int indent,char *title,ivec vec[],int n);
extern void pr_rvec(FILE *fp,int indent,char *title,real vec[],int n);
extern void pr_rvecs(FILE *fp,int indent,char *title,rvec vec[],int n);
extern void pr_rvecs_len(FILE *fp,int indent,char *title,rvec vec[],int n);
extern void pr_block(FILE *fp,int indent,char *title,t_block *block);
extern void pr_iparams(FILE *fp,t_functype ftype,t_iparams *iparams);
extern void pr_idef(FILE *fp,int indent,char *title,t_idef *idef);
extern void pr_inputrec(FILE *fp,int indent,char *title,t_inputrec *ir);
extern void pr_top(FILE *fp,int indent,char *title,t_topology *top);
/*
 * This routine prints out a (human) readable representation of 
 * the topology to the file fp. Ident specifies the number of 
 * spaces the text should be indented. Title is used to print a 
 * header text.
 */
extern void pr_header(FILE *fp,int indent,char *title,t_tpxheader *sh);
     /*
      * This routine prints out a (human) readable representation of
      * a header to the file fp. Ident specifies the number of spaces
      * the text should be indented. Title is used to print a header text.
      */



#endif	/* _txtdump_h */
