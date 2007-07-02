/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _txtdump_h
#define _txtdump_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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
extern int pr_indent(FILE *fp,int n);
extern int available(FILE *fp,void *p,int indent,const char *title);
extern int pr_title(FILE *fp,int indent,const char *title);
extern int pr_title_n(FILE *fp,int indent,const char *title,int n);
extern int pr_title_nxn(FILE *fp,int indent,const char *title,int n1,int n2);
extern void pr_ivec(FILE *fp,int indent,const char *title,int vec[],int n, bool bShowNumbers);
extern void pr_ivecs(FILE *fp,int indent,const char *title,ivec vec[],int n, bool bShowNumbers);
extern void pr_rvec(FILE *fp,int indent,const char *title,real vec[],int n, bool bShowNumbers);
extern void pr_rvecs(FILE *fp,int indent,const char *title,rvec vec[],int n);
extern void pr_rvecs_len(FILE *fp,int indent,const char *title,rvec vec[],int n);
extern void pr_reals(FILE *fp,int indent,const char *title,real vec[],int n);
extern void pr_block(FILE *fp,int indent,const char *title,t_block *block,bool bShowNumbers);
extern void pr_ilist(FILE *fp,int indent,const char *title,
		     t_idef *idef,t_ilist *ilist, bool bShowNumbers);
extern void pr_iparams(FILE *fp,t_functype ftype,t_iparams *iparams);
extern void pr_idef(FILE *fp,int indent,const char *title,t_idef *idef, bool bShowNumbers);
extern void pr_inputrec(FILE *fp,int indent,const char *title,t_inputrec *ir);
extern void pr_atoms(FILE *fp,int indent,const char *title,t_atoms *atoms, 
		     bool bShownumbers);
extern void pr_atomtypes(FILE *fp,int indent,const char *title,
			 t_atomtypes *atomtypes,bool bShowNumbers);
extern void pr_top(FILE *fp,int indent,const char *title,t_topology *top, bool bShowNumbers);
/*
 * This routine prints out a (human) readable representation of 
 * the topology to the file fp. Ident specifies the number of 
 * spaces the text should be indented. Title is used to print a 
 * header text.
 */
extern void pr_header(FILE *fp,int indent,const char *title,t_tpxheader *sh);
      /*
      * This routine prints out a (human) readable representation of
      * a header to the file fp. Ident specifies the number of spaces
      * the text should be indented. Title is used to print a header text.
      */

extern void pr_commrec(FILE *fp,int indent,t_commrec *cr);

#endif	/* _txtdump_h */
