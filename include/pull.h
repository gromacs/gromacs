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
  eStart, eAfm, eConstraint, eUmbrella, eTest
} t_runtype;
typedef enum {
  eCom, eComT0, eDyn, eDynT0
} t_reftype;

typedef struct {
  int        n;         /* number of groups */
  atom_id    **idx;     /* indices of pull atoms in full coordinate array */
  real       **weights; /* position depended weight (switch function) */
  int        *ngx;      /* pull group sizes */
  char       **grps;    /* pull group names */
  real       *tmass;    /* total mass of the groups */
  rvec       **x0;      /* pull group coordinates at t=0 */
  rvec       **xp;      /* pull group coordinates at previous step */
  rvec       *x_ref;    /* reference positions */
  rvec       *x_unc;    /* center of mass before constraining */
  rvec       *x_con;    /* center of mass, obeying constraints */
  rvec       *xprev;    /* position of coms in last written structure */
  rvec       *f;        /* forces due to the pulling/constraining */
  rvec       *spring;   /* coordinates of the springs (eAfm) */
  rvec       *dir;      /* direction of constraint */
  real       *d_ref;    /* reference distance  */
  rvec       *xtarget;  /* target coordinates for structure generation */
  rvec       **comhist; /* com over the last nhist steps (for running aver) */
} t_pullgrps; 

typedef struct {
  t_pullgrps dyna;      /* dynamic groups for use with local constraints */
  t_pullgrps pull;      /* groups to pull/restrain/etc/ */
  t_pullgrps ref;       /* reference group, reaction force grps */
  t_runtype  runtype;   /* start, afm, constraint, umbrella, test */
  t_reftype  reftype;   /* com, com_t0, dynamic, dynamic_t0 */
  rvec       dims;      /* used to select components for constraint */
  rvec       coor;      /* reaction coordinate */
  real       r;         /* radius of cylinder for dynamic COM */
  real       rc;        /* radius of cylinder including switch length */
  real       start_k0;  /* starting force constant */
  real       start_k1;  /* ending force constant */
  real       k_rate;    /* switch between k0 and k1 over this many steps */
  int        k_step;    /* current step */
  real       xlt_incr;  /* write out structure every xlt_incr nm */
  real       tolerance; /* tolerance for reaching desired coordinates (nm) */
  real       constr_tol;/* absolute tolerance for constraints in (nm) */
  bool       bPull;     /* true if we're doing any pulling */
  bool       bCyl;      /* true if we're using dynamic ref. groups */
  bool       bReverse;  /* reverse reference direction */
  FILE       *out;      /* output file for pull data */
  real       k;         /* force constant for atoms */
  real       rate;      /* pull rate, in nm/timestep */
  int        update;    /* update frequency for dynamic grps */
  int        reflag;    /* running average over reflag steps for com */
  bool       bVerbose;  /* be loud and noise */
  rvec       UmbPos[4]; /* center of umbrella potentials */
  real       UmbCons[4];/* force constant of umbrella potential */
  int        nSkip;     /* only write output every nSkip steps */
  bool       bCompress; /* compress output */
  int        start_nout;/* should we output starting structure? */
  bool       bFirst;    /* is this the first step for dynamic ref group? */
} t_pull;

/* main pull routine that controls all the action */
extern void pull(t_pull *pull,    /* all pull data */
                 rvec *x,         /* coordinates, changed by constraint run */ 
                 rvec *f,         /* forces, changed by Afm run */
                 matrix box,               
                 t_topology *top, /* needed to write out coordinate files */   
                 real dt,         /* time step */
                 int step,        /* step number in simulation */
                 int natoms,      /* total number of atoms on this processor */
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
                      matrix box,     
                      int start,      /* startinig index of this node */
                      int homenr,     /* number of atoms on this node */
                      t_commrec * cr  /* struct for communication info */
                      );
#endif
