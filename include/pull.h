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

#ifndef _pull_h
#define _pull_h

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

#include "vec.h"
#include "typedefs.h"

/* This file contains datatypes and function declarations necessary 
   for mdrun to interface with the pull code */

typedef enum {
  eAfm, eConstraint, eUmbrella, ePullruntypeNR
} t_pullruntype;
typedef enum {
  eCom, eComT0, eDyn, eDynT0
} t_pullreftype;

typedef struct {
  atom_id    *idx;     /* indices of pull atoms in full coordinate array */
  int        nweight;  /* The number of weights read from the param file */
  real       *weight;  /* weights (use all 1 when nweight==0) */
  int        ngx;      /* pull group size */
  char       *name;    /* pull group name */
  real       wscale;   /* scaling factor for the weights: sum w m/sum w w m */
  real       invtm;    /* inverse total mass of the group: 1/wscale sum w m */
  rvec       *x0;      /* pull group coordinates at t=0 */
  rvec       *xp;      /* pull group coordinates at previous step */
  dvec       x_ref;    /* reference position */
  dvec       x_unc;    /* center of mass before constraining */
  dvec       x_con;    /* center of mass, obeying constraints */
  dvec       xprev;    /* position of coms in last written structure */
  dvec       f;        /* force due to the pulling/constraining */
  dvec       spring;   /* coordinates of the spring (eAfm) */
  dvec       dir;      /* direction of constraint */
  dvec       xtarget;  /* target coordinates for structure generation */
  dvec       *comhist; /* com over the last nhist steps (for running aver) */
  dvec       AfmVec;   /* Vector to pull along for AFM */
  real       AfmK;     /* Force constant to use for AFM */
  real       AfmRate;  /* Pull rate in nm/ps */
  dvec       AfmInit;  /* Initial sprint posistions for AFM */
  dvec       UmbPos;   /* center of umbrella potential */
  real       UmbCons;  /* force constant of umbrella potential */
} t_pullgrp; 

typedef struct {
  t_pullgrp  ref;         /* reference group, reaction force grps */
  int        ngrp;        /* number of groups */
  t_pullgrp  *grp;        /* groups to pull/restrain/etc/ */
  t_pullgrp  *dyna;       /* dynamic groups for use with local constraints */
  t_pullruntype  runtype; /* start, afm, constraint, umbrella, test */
  t_pullreftype  reftype; /* com, com_t0, dynamic, dynamic_t0 */
  ivec       dims;        /* used to select components for constraint */
  int        bDir;        /* use only the direction dir */
  dvec       dir;         /* direction */
  real       r;           /* radius of cylinder for dynamic COM */
  real       rc;          /* radius of cylinder including switch length */
  real       constr_rate; /* rate of change of the constraint length (nm/ps) */
  real       constr_tol;  /* absolute tolerance for constraints in (nm) */
  bool       bPull;       /* true if we're doing any pulling */
  bool       bCyl;        /* true if we're using dynamic ref. groups */
  FILE       *out;        /* output file for pull data */
  int        update;      /* update frequency for dynamic grps */
  int        reflag;      /* running average over reflag steps for com */
  bool       AbsoluteRef; /* Reference is in absolute coordinates */
  bool       bVerbose;    /* be loud and noise */
  int        nSkip;       /* only write output every nSkip steps */
} t_pull;

/* main pull routine that controls all the action */
extern void pull(t_pull *pull,    /* all pull data */
                 rvec *x,         /* coordinates, changed by constraint run */ 
                 rvec *f,         /* forces, changed by Afm run */
                 matrix box,               
                 t_topology *top, /* needed to write out coordinate files */   
                 real dt,         /* time step */
                 int step,        /* step number in simulation */
                 real time,       /* time, only used for printing output */
                 t_mdatoms *md,   /* masses and charges of all atoms */
                 int start,       /* number of first atom belonging to this node */
                 int homenr,      /* number of atoms that belong to this node */
                 t_commrec * cr   /* Stuff for communication */
                );


/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
extern void init_pull(FILE *log,  
                      int nfile,       
                      t_filenm fnm[], /* standard filename struct */
                      t_pull *pull,   /* all pull data */
                      rvec *x,        /* all coordinates */
                      t_mdatoms *md,  /* masses and charges of all atoms */
		      ivec nFreeze[], /* the freeze dimensions */
                      matrix box,     
                      int start,      /* startinig index of this node */
                      int homenr,     /* number of atoms on this node */
                      t_commrec * cr  /* struct for communication info */
                      );
#endif
