/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Getting the Right Output Means no Artefacts in Calculating Stuff
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
