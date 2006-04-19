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

#define MAX_PULL_GROUPS 20

void print_afm(t_pull *pull, int step, real t) 
{
  int i,j;
  
  /* Do we need to do output? */
  if (step % pull->nSkip) return;
  
  /* Print time */
  fprintf(pull->out, "%f", t);
  
  /* Print COM of reference group */
  for(j=0; j<DIM; ++j) {
    if(pull->dims[j] != 0.0) {
      fprintf(pull->out, "\t%f", pull->ref.x_unc[j]);
    }
  }
  
  /* Print out position of pulled groups */
  for (i=0;i<pull->ngrp;i++) {
    for (j=0;j<3;++j) {
      if (pull->dims[j] != 0.0) {
	fprintf(pull->out,"\t%f\t%f",
		pull->grp[i].x_unc[j],pull->grp[i].spring[j]);
      }
    }
  }
  fprintf(pull->out,"\n");
}

void print_constraint(t_pull *pull, int step, real t) 
{
  int i;

  /* step < 0 means shakefirst, yuk, yuk */
  if (step < 0 || step % pull->nSkip) return;

  fprintf(pull->out, "%f\t", t);

  for(i=0;i<pull->ngrp;i++) {
    /* Print the constraint force */
    fprintf(pull->out,"\t%f",pull->grp[i].f[ZZ]);
  }
  fprintf(pull->out,"\n");
}

void print_umbrella(t_pull *pull, int step, real t)
{
  int i,m;

  /* Do we need to do any pulling ? */
  if (step % pull->nSkip) return;

  fprintf(pull->out, "%f", t);
  
  /* Print deviation of pulled group from desired position */
  for (i=0; i<pull->ngrp; ++i) {    /* Loop over pulled groups */
    for (m=0; m<DIM; ++m) {             /* Loop over dimensions */
      if (pull->dims[m] != 0.0) {
	fprintf(pull->out,"\t%f",-pull->grp[i].spring[m]);
      }
    }
  }
  fprintf(pull->out,"\n");
}

static void init_pullgrp(t_pullgrp *pg,char *name,char *wbuf,real UmbCons)
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
  pg->UmbCons = UmbCons;
}

static void string2dvec(char buf[], dvec nums)
{
  if (sscanf(buf,"%lf%lf%lf",&nums[0],&nums[1],&nums[2]) != 3)
    gmx_fatal(FARGS,"Expected three numbers at input line %s",buf);
}

void read_pullparams(t_pull *pull, char *infile, char *outfile) 
{
  t_inpfile *inp;
  int ninp,i,nchar;
  char *tmp,*ptr;
  char dummy[STRLEN],buf[STRLEN];
  char refbuf[STRLEN],grpbuf[MAX_PULL_GROUPS][STRLEN],
    refwbuf[STRLEN],wbuf[MAX_PULL_GROUPS][STRLEN],
    bf[MAX_PULL_GROUPS][STRLEN],
    pos[MAX_PULL_GROUPS][STRLEN],
    pulldim[STRLEN],pulldim1[STRLEN],condir[STRLEN];
  char DirTemp[MAX_PULL_GROUPS][STRLEN], InitTemp[MAX_PULL_GROUPS][STRLEN];
  real constr_d0[MAX_PULL_GROUPS],constr_rate[MAX_PULL_GROUPS];
  real AfmRate[MAX_PULL_GROUPS],AfmK[MAX_PULL_GROUPS],UmbCons[MAX_PULL_GROUPS];

  int bReverse; int tmpref; int tmprun; 

  static const char *runtypes[ePullruntypeNR+1] = { 
    "afm", "constraint", "umbrella", NULL
  };
  enum {
    erefCom, erefComT0, erefDyn, erefDynT0, erefNR
  };
  static const char *reftypes[erefNR+1] = {
    "com", "com_t0", "dynamic", "dynamic_t0", NULL
  };
  enum {
    everboseYES, everboseNO, everboseNR
  };
  static const char *verbosetypes[erefNR+1] = {
    "no", "yes", NULL
  };
  int nerror = 0;

  /* read input parameter file */
  fprintf(stderr,"Reading parameter file %s\n",infile);
  inp=read_inpfile(infile,&ninp);

  /* general options */
  CTYPE("GENERAL");
  EETYPE("verbose",         pull->bVerbose, verbosetypes, &nerror, TRUE);
  ITYPE("Skip steps",       pull->nSkip,1);
  CTYPE("Runtype: afm, constraint, umbrella");
  EETYPE("runtype",         tmprun, runtypes, &nerror, TRUE);
  CTYPE("Number of pull groups");
  ITYPE("ngroups",          pull->ngrp,1);
  if (pull->ngrp < 1 || pull->ngrp > MAX_PULL_GROUPS)
    gmx_fatal(FARGS,"ngroups should be between 1 and %d",MAX_PULL_GROUPS);
  CTYPE("Groups to be pulled");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"group_%d",i+1);
    STYPE(buf,              grpbuf[i], "");
  }
  CTYPE("The group for the reaction force.");
  STYPE("reference_group",  refbuf, "");
  CTYPE("Weights for all atoms in each group (default all 1)");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"weights_%d",i+1);
    STYPE(buf,              wbuf[i], "");
  }
  STYPE("reference_weights", refwbuf, "");
  CTYPE("Ref. type: com, com_t0, dynamic, dynamic_t0");
  EETYPE("reftype",         tmpref, reftypes, &nerror, TRUE);
  CTYPE("Use running average for reflag steps for com calculation");
  ITYPE("reflag",           pull->reflag, 1);
  CTYPE("Select components for the pull vector. default: Y Y Y");
  STYPE("pulldim",          pulldim, "Y Y Y");

  /* options for dynamic reference groups */
  CTYPE("DYNAMIC REFERENCE GROUP OPTIONS");
  CTYPE("Cylinder radius for dynamic reaction force groups (nm)");
  RTYPE("r",                pull->r, 0.0);
  CTYPE("Switch from r to rc in case of dynamic reaction force");
  RTYPE("rc",   pull->rc,   0.0);
  CTYPE("Update frequency for dynamic reference groups (steps)");
  ITYPE("update",           pull->update, 1);
    

  /* constraint run options */
  CCTYPE("CONSTRAINT RUN OPTIONS");
  CTYPE("Direction, default: 0 0 0, no direction");
  STYPE("constraint_direction",        condir, "0.0 0.0 0.0");
  CTYPE("Constraint distance (nm), default: 0, use starting distance");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"constraint_distance%d",i+1);
    RTYPE(buf,              constr_d0[i], 0.0);
  }
  CTYPE("Rate of chance of the constraint length, in nm/ps");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"constraint_rate%d",i+1);
    RTYPE(buf,              constr_rate[i], 0);
  }
  CTYPE("Tolerance of constraints, in nm");
  RTYPE("constraint_tolerance",            pull->constr_tol, 1E-6);

   /* options for AFM type pulling simulations */
  CCTYPE("AFM OPTIONS");
  CTYPE("Pull rates in nm/ps");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"afm_rate%d",i+1);
    RTYPE(buf,              AfmRate[i], 0.0);
  }
  CTYPE("Force constants in kJ/(mol*nm^2)");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"afm_k%d",i+1);
    RTYPE(buf,              AfmK[i], 0.0);
  }
  CTYPE("Directions");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"afm_dir%d",i+1);
    STYPE(buf,              DirTemp[i], "0.0 0.0 1.0");
  }
  CTYPE("Initial spring positions");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"afm_init%d",i+1);
    STYPE(buf,              InitTemp[i], "0.0 0.0 0.0");
  }
  
  /* umbrella sampling options */
  CCTYPE("UMBRELLA SAMPLING OPTIONS");
  CTYPE("Force constants for umbrella sampling in kJ/(mol*nm^2)");
  CTYPE("Centers of umbrella potentials with respect to reference:");
  CTYPE("Ref - Pull.");
  for(i=0; i<pull->ngrp; i++) {
    sprintf(buf,"K%d",i+1);
    RTYPE(buf,              UmbCons[i], 0.0);
    sprintf(buf,"Pos%d",i+1);
    STYPE(buf,              pos[i], "0.0 0.0 0.0");
  }

  write_inpfile(outfile,ninp,inp,TRUE);
  for(i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  pull->runtype = (t_pullruntype)tmprun;
  pull->reftype = (t_pullreftype)tmpref;

  /* sort out the groups */
  fprintf(stderr,"Groups:");
  for(i=0; i<pull->ngrp; i++)
    fprintf(stderr," %s",grpbuf[i]);
  fprintf(stderr,"\n");

  if (!strcmp(refbuf, "")) {
    if (pull->runtype == eConstraint)
      gmx_fatal(FARGS, "Constraint forces require a reference group to be specified.\n");
    pull->AbsoluteRef = TRUE;
    fprintf(stderr, "Pull code using absolute reference.\n");
  }
  else if (pull->AbsoluteRef == FALSE) {
    fprintf(stderr,"Reference Group:  %s\n",refbuf);
  }
  
  fprintf(stderr,"Using %d pull groups\n",pull->ngrp);
  
  /* Make the pull groups */
  init_pullgrp(&pull->ref,refbuf,refwbuf,0);
  snew(pull->grp,pull->ngrp);
  for(i=0; i<pull->ngrp; i++)
    init_pullgrp(&pull->grp[i],grpbuf[i],wbuf[i],UmbCons[i]);
  
  if (pull->runtype == eConstraint) {
    for (i=0; i<pull->ngrp; ++i) {
      pull->grp[i].constr_d0 = constr_d0[i];
      pull->grp[i].constr_rate = constr_rate[i];
    }
  }

  if (pull->runtype == eAfm) {
    for (i=0; i<pull->ngrp; ++i) {
      pull->grp[i].AfmRate = AfmRate[i];
      pull->grp[i].AfmRate = AfmRate[i];
      pull->grp[i].AfmK = AfmK[i];
      string2dvec(DirTemp[i], pull->grp[i].AfmVec);
      string2dvec(InitTemp[i], pull->grp[i].AfmInit);
    }
  }

  if(pull->runtype == eUmbrella) {
    if(! (pull->reftype==erefCom || pull->reftype==erefComT0))
      gmx_fatal(FARGS,"Umbrella sampling currently only works with COM and COM_T0 reftypes");
    for(i=0;i<pull->ngrp;++i) {
      string2dvec(pos[i],pull->grp[i].UmbPos);
    }
  }
  /* End Modification */

  ptr = pulldim;
  for(i=0; i<DIM; i++) {
    if (sscanf(ptr,"%s%n",pulldim1,&nchar) != 1)
      gmx_fatal(FARGS,"Less than 3 pull dimensions given in pulldim: '%s'",
		  pulldim);
    
    if (strncasecmp(pulldim1,"N",1) == 0) {
      pull->dims[i] = 0;
    } else if (strncasecmp(pulldim1,"Y",1) == 0) {
      pull->dims[i] = 1;
    } else {
      gmx_fatal(FARGS,"Please use Y(ES) or N(O) for pulldim only (not %s)",
		  pulldim1);
    }
    ptr += nchar;
  }

  fprintf(stderr,"Using distance components %d %d %d\n",
	  pull->dims[0],pull->dims[1],pull->dims[2]);

  string2dvec(condir,pull->dir);
  if (pull->runtype != eConstraint ||
      (pull->dir[XX] == 0 && pull->dir[YY] == 0 && pull->dir[ZZ] == 0)) {
    pull->bDir = FALSE;
    if (pull->runtype == eConstraint)
      fprintf(stderr,"Not using directional constraining");
  } else {
    pull->bDir = TRUE;
    dsvmul(1/dnorm(pull->dir),pull->dir,pull->dir);
    fprintf(stderr,"Using directional constraining %5.2f %5.2f %5.2f\n",
    pull->dir[XX],pull->dir[YY],pull->dir[ZZ]);
  }

  if(pull->r > 0.001)
    pull->bCyl = TRUE;
  else
    pull->bCyl = FALSE;
}

void print_pull_header(t_pull * pull)
{
  /* print header information */
  int i,j;

  if(pull->runtype==eUmbrella)
    fprintf(pull->out,"# UMBRELLA\t3.0\n");
  else if(pull->runtype==eAfm)
    fprintf(pull->out,"# AFM\t3.0\n");
  else if(pull->runtype==eConstraint)
    fprintf(pull->out,"# CONSTRAINT\t3.0\n");

  fprintf(pull->out,"# Component selection:");
  for(i=0;i<3;++i) {
    fprintf(pull->out," %d",pull->dims[i]);
  }
  fprintf(pull->out,"\n");

  fprintf(pull->out,"# nSkip %d\n",pull->nSkip);

  fprintf(pull->out,"# Ref. Group '%s'\n",pull->ref.name);

  fprintf(pull->out,"# Nr. of pull groups %d\n",pull->ngrp);
  for(i=0;i<pull->ngrp;i++) {
    fprintf(pull->out,"# Group %d '%s'",i+1,pull->grp[i].name);

    if (pull->runtype == eAfm) {
      fprintf(pull->out, "  afmVec %f %f %f  AfmRate %f  AfmK %f",
	      pull->grp[i].AfmVec[XX],
	      pull->grp[i].AfmVec[YY],
	      pull->grp[i].AfmVec[ZZ],
	      pull->grp[i].AfmRate,
	      pull->grp[i].AfmK);
    } else if (pull->runtype == eUmbrella) {
      fprintf(pull->out,"  Umb. Pos.");
      for (j=0;j<3;++j) {
	if (pull->dims[j] != 0)
	  fprintf(pull->out," %f",pull->grp[i].UmbPos[j]);
      }
      fprintf(pull->out,"  Umb. Cons.");
      for (j=0;j<3;++j) {
	if (pull->dims[j] != 0)
	  fprintf(pull->out," %f",pull->grp[i].UmbCons);
      }
    } else if (pull->runtype == eConstraint) {
      fprintf(pull->out,"  Pos.");
      for (j=0;j<3;++j) {
	if (pull->dims[j] != 0)
	  fprintf(pull->out," %f",pull->grp[i].x_unc[j]);
      }
    }
    fprintf(pull->out,"\n");
  }
  fprintf(pull->out,"#####\n");

}
