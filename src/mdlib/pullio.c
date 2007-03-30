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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "princ.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pull.h"
#include "pull_internal.h"
#include "string.h"


/*  None of the routines in this file check to see if we are
    doing parallel pulling. It is the responsibility of the 
    calling function to ensure that these routines only get
    called by the master node.
*/

#define MAX_PULL_GROUPS 21

static void pull_print_x(FILE *out,t_pull *pull, real t) 
{
  int g,j;
  
  /* Print time */
  fprintf(out, "%f", t);

  if (pull->grp[0].nat && pull->cyl_r0 == 0) {
    /* Print COM of reference group */
    for(j=0; j<DIM; ++j) {
      if (pull->dim[j] != 0.0)
	fprintf(out, "\t%f", pull->grp[0].x[j]);
    }
  }
  
  /* Print out position of pulled groups */
  for (g=1; g<1+pull->ngrp; g++) {
    for (j=0;j<3;++j) {
      if (pull->dim[j] != 0)
	fprintf(out,"\t%f",pull->grp[g].x[j]);
    }
  }
  fprintf(out,"\n");
}

static void pull_print_f(FILE *out,t_pull *pull, real t) 
{
  int g,d;

  fprintf(out, "%f\t", t);

  for(g=1; g<1+pull->ngrp; g++) {
    if (pull->eGeom == epullgPOS) {
      for(d=0; d<DIM; d++)
	if (pull->dim[d])
	  fprintf(out,"\t%f",pull->grp[g].f[d]);
    } else {
      fprintf(out,"\t%f",pull->grp[g].f_scal);
    }
  }
  fprintf(out,"\n");
}

void pull_print_output(t_pull *pull, int step, real time)
{
  if ((pull->nstxout != 0) && (step % pull->nstxout == 0))
    pull_print_x(pull->out_x,pull,time);

  if ((pull->nstfout != 0) && (step % pull->nstfout == 0))
    pull_print_f(pull->out_f,pull,time);
}

static void string2dvec(char buf[], dvec nums)
{
  if (sscanf(buf,"%lf%lf%lf",&nums[0],&nums[1],&nums[2]) != 3)
    gmx_fatal(FARGS,"Expected three numbers at input line %s",buf);
}

static void init_pullgrp(t_pullgrp *pg,char *name,char *wbuf,
			 bool bRef,int eGeom,char *s_vec)
{
  double d;
  int n;

  pg->name = strdup(name);
  pg->nweight = 0;
  while (sscanf(wbuf,"%lf %n",&d,&n) == 1) {
    if (pg->nweight == 0) {
      snew(pg->weight,1);
    } else {
      srenew(pg->weight,pg->nweight+1);
    }
    pg->weight[pg->nweight++] = d;
    wbuf += n;
  }
  if (!bRef) {
    if (eGeom != epullgDIST) {
      string2dvec(s_vec,pg->vec);
      if (eGeom == epullgDIR || eGeom == epullgCYL)
	/* Normalize the direction vector */
	dsvmul(1/dnorm(pg->vec),pg->vec,pg->vec);
    }
  }
}

void read_pullparams(t_pull *pull, char *infile, char *outfile) 
{
  t_inpfile *inp;
  int ninp,i,nchar;
  char *tmp,*ptr;
  char dummy[STRLEN],buf[STRLEN];
  char grpbuf[MAX_PULL_GROUPS][STRLEN],
    wbuf[MAX_PULL_GROUPS][STRLEN],
    bf[MAX_PULL_GROUPS][STRLEN],
    pos[MAX_PULL_GROUPS][STRLEN],
    pulldim[STRLEN],pulldim1[STRLEN];
  char VecTemp[MAX_PULL_GROUPS][STRLEN], InitTemp[MAX_PULL_GROUPS][STRLEN];

  t_pullgrp *pgrp;

  int nerror = 0;

  /* read input parameter file */
  fprintf(stderr,"Reading parameter file %s\n",infile);
  inp=read_inpfile(infile,&ninp);

  CTYPE("GENERAL");
  CTYPE("Pull type: no, umbrella, constraint");
  EETYPE("pull",            pull->ePull, epull_names, &nerror, TRUE);
  CTYPE("Pull geometry: distance, direction, cylinder or position");
  EETYPE("pull_geometry",   pull->eGeom, epullg_names, &nerror, TRUE);
  CTYPE("Select components for the pull vector. default: Y Y Y");
  STYPE("pull_dim",         pulldim, "Y Y Y");
  CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
  RTYPE("pull_r1",          pull->cyl_r1, 0.0);
  CTYPE("Switch from r1 to r0 in case of dynamic reaction force");
  RTYPE("pull_r0",          pull->cyl_r0, 0.0);
  RTYPE("constraint_tolerance",            pull->constr_tol, 1E-6);
  EETYPE("pull_start",      pull->bStart, yesno_names, &nerror, TRUE);
  ITYPE("pull_nstxout",     pull->nstxout, 10);
  ITYPE("pull_nstfout",     pull->nstfout,  1);

  CCTYPE("PULL GROUP OPTIONS");
  CTYPE("Number of pull groups");
  ITYPE("ngroups",          pull->ngrp,1);
  if (pull->ngrp < 1 || pull->ngrp > MAX_PULL_GROUPS-1)
    gmx_fatal(FARGS,"ngroups should be between 1 and %d",MAX_PULL_GROUPS-1);

  snew(pull->grp,pull->ngrp+1);

  /* pull group options */
  CTYPE("Group name, weight (default all 1), vector, init, rate (nm/ps), kJ/(mol*nm^2)");
  for(i=0; i<pull->ngrp+1; i++) {
    pgrp = &pull->grp[i];
    sprintf(buf,"group%d",i);
    STYPE(buf,              grpbuf[i], "");
    sprintf(buf,"weights%d",i);
    STYPE(buf,              wbuf[i], "");
    sprintf(buf,"pbcatom%d",i);
    ITYPE(buf,              pgrp->pbcatom, 0);
    if (i > 0) {
      sprintf(buf,"vec%d",i);
      STYPE(buf,              VecTemp[i], "0.0 0.0 0.0");
      sprintf(buf,"init%d",i);
      RTYPE(buf,              pgrp->init, 0.0);
      sprintf(buf,"rate%d",i);
      RTYPE(buf,              pgrp->rate, 0.0);
      sprintf(buf,"k%d",i);
      RTYPE(buf,              pgrp->k, 0.0);
    }
  }

  write_inpfile(outfile,ninp,inp,TRUE);
  for(i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  if (PULL_CYL(pull)) {
    pull->dim[0] = 0;
    pull->dim[1] = 0;
    pull->dim[2] = 1;
  } else {
    ptr = pulldim;
    for(i=0; i<DIM; i++) {
      if (sscanf(ptr,"%s%n",pulldim1,&nchar) != 1)
	gmx_fatal(FARGS,"Less than 3 pull dimensions given in pulldim: '%s'",
		  pulldim);
      
      if (strncasecmp(pulldim1,"N",1) == 0) {
	pull->dim[i] = 0;
      } else if (strncasecmp(pulldim1,"Y",1) == 0) {
	pull->dim[i] = 1;
      } else {
	gmx_fatal(FARGS,"Please use Y(ES) or N(O) for pulldim only (not %s)",
		  pulldim1);
      }
      ptr += nchar;
    }
  }

  fprintf(stderr,"Using distance components %d %d %d\n",
	  pull->dim[0],pull->dim[1],pull->dim[2]);

  if (!strcmp(grpbuf[0], "")) {
    fprintf(stderr,"Pull code using absolute reference.\n");
  } else {
    fprintf(stderr,"Reference Group:  %s\n",grpbuf[0]);
  }
  
  /* sort out the groups */
  fprintf(stderr,"Using %d pull groups\n",pull->ngrp);
  for(i=0; i<pull->ngrp+1; i++)
    fprintf(stderr," %s",grpbuf[i]);
  fprintf(stderr,"\n");
  
  /* Initialize the pull groups */
  for(i=0; i<pull->ngrp+1; i++) {
    init_pullgrp(&pull->grp[i],grpbuf[i],wbuf[i],i==0,pull->eGeom,VecTemp[i]);
  }
}

void print_pull_header(FILE *out,t_pull * pull)
{
  /* print header information */
  int i,j;

  if (pull->ePull == epullUMBRELLA)
    fprintf(out,"# UMBRELLA\t3.0\n");
  else if (pull->ePull == epullCONSTRAINT)
    fprintf(out,"# CONSTRAINT\t3.0\n");

  fprintf(out,"# Component selection:");
  for(i=0;i<3;++i) {
    fprintf(out," %d",pull->dim[i]);
  }
  fprintf(out,"\n");

  fprintf(out,"# nstxcout %d\n",pull->nstxout);
  fprintf(out,"# nstfcout %d\n",pull->nstfout);

  fprintf(out,"# Ref. Group '%s'\n",pull->grp[0].name);

  fprintf(out,"# Nr. of pull groups %d\n",pull->ngrp);
  for(i=1;i<pull->ngrp+1;i++) {
    fprintf(out,"# Group %d '%s'",i,pull->grp[i].name);

    fprintf(out, "  dir %f %f %f  init %f  rate %f",
	    pull->grp[i].vec[XX],
	    pull->grp[i].vec[YY],
	    pull->grp[i].vec[ZZ],
	    pull->grp[i].init,
	    pull->grp[i].rate);

    if (pull->ePull == epullUMBRELLA) {
      fprintf(out, "  k %f",pull->grp[i].k);
    } else if (pull->ePull == epullCONSTRAINT) {
      fprintf(out,"  Pos.");
      for (j=0;j<3;++j) {
	if (pull->dim[j] != 0)
	  fprintf(out," %f",pull->grp[i].x[j]);
      }
    }
    fprintf(out,"\n");
  }
  fprintf(out,"#####\n");
}
