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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _readir_h
#define _readir_h

#include "typedefs.h"
#include "string2.h"
#include "readinp.h"
#include "grompp.h"

enum { eshNONE, eshHBONDS, eshALLBONDS, eshHANGLES, eshALLANGLES, eshNR };

static const char *constraints[eshNR+1]    = { 
    "none", "h-bonds", "all-bonds", "h-angles", "all-angles", NULL 
  };

typedef struct {
  int warnings;
  int nshake;
  real fourierspacing;
  int nprocs;
  int splitalg;
  char *title;
  char *cpp;
  char *include;
  char *define;
  char *SolventOpt;
  bool bGenVel;
  bool bGenPairs;
  real tempi;
  int  seed;
  int  eDisre;
  bool bOrire;
  int  eDihre;
  bool bMorse;
  char *wall_atomtype[2];
  bool pull_start;
} t_gromppopts;


extern void init_ir(t_inputrec *ir, t_gromppopts *opts);
/* Initiate stuff */

extern void check_ir(t_inputrec *ir, t_gromppopts *opts,int *nerror);
/* Validate inputrec data.
 * Fatal errors will be added to nerror.
 */
extern int search_string(char *s,int ng,char *gn[]);
/* Returns the index of string s in the index groups */

extern void double_check(t_inputrec *ir,matrix box,t_molinfo *mol,int *nerror);
/* Do more checks */

extern void triple_check(char *mdparin,t_inputrec *ir,t_topology *sys,
			 int *nerror);
/* Do even more checks */

extern void get_ir(char *mdparin,char *mdparout,
		   t_inputrec *ir,t_gromppopts *opts,int *nerror);
/* Read the input file, and retrieve data for inputrec.
 * More data are read, but the are only evaluated when the next
 * function is called. Also prints the input file back to mdparout.
 * Add errors no nerror.
 */
 
extern void do_index(char *ndx,
		     t_symtab   *symtab,
		     t_atoms    *atoms,bool bVerbose,
		     t_inputrec *ir,t_idef *idef,
		     rvec *v);
/* Read the index file and assign grp numbers to atoms.
 * If v is not NULL, the velocities will be scaled to the correct number
 * of degrees of freedom.
 */

/* Routines In readpull.c */

extern char **read_pullparams(int *ninp_p,t_inpfile **inp,
			      t_pull *pull,bool *bStart,int *nerror_p);
/* Reads the pull parameters, returns a list of the pull group names */

extern void make_pull_groups(t_pull *pull,char **pgnames,
			     t_block *grps,char **gnames);
/* Process the pull parameters after reading the index groups */

extern void set_pull_init(t_inputrec *ir,t_atoms *atoms,rvec *x,matrix box,
			  bool bStart);
/* Prints the initial pull group distances in x.
 * If bStart adds the distance to the initial reference location.
 */

#endif	/* _readir_h */
