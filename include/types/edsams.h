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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _TYPES_EDSAMS_H_
#define _TYPES_EDSAMS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


typedef struct {
  rvec *x;
  rvec *transvec;
  rvec *forces_cartesian;
} t_edlocals;

typedef struct {
  int 		neig;		/* nr of eigenvectors            	*/
  int           *ieig;          /* index nrs of eigenvectors            */
  real          *stpsz;         /* stepsizes (per eigenvector)          */
  rvec          **vec;          /* eigenvector components               */
  real          *xproj;         /* instantaneous x projections          */
  real          *vproj;         /* instantaneous v projections          */
  real          *fproj;         /* instantaneous f projections          */
  real          *refproj;       /* starting or target projecions        */
  real          radius;         /* instantaneous radius                 */
} t_eigvec;

typedef struct {
  t_eigvec 	mon;		/* only monitored, no constraints       */
  t_eigvec 	linfix;		/* fixed linear constraints             */
  t_eigvec 	linacc;		/* acceptance linear constraints        */
  t_eigvec 	radfix;		/* fixed radial constraints (exp)       */
  t_eigvec 	radacc;		/* acceptance radial constraints (exp)  */  
  t_eigvec 	radcon;		/* acceptance rad. contraction constr.  */
} t_edvecs;

typedef struct {
  int 		nr;		/* Nr. of atoms        			*/
  int           *anrs;          /* Index numbers                        */
  rvec          *x;             /* Positions                            */
  matrix        box;            /* Box lenghts                          */
  real 	        *sqrtm;	        /* masses used for mass-weighting of analysis, only used in sav*/
} t_edx;

typedef struct { 
  real deltaF0;
  bool bHarmonic;
  real tau;
  real deltaF;
  real Efl;
  real kT; 
  real Vfl;
  real dt;
  real constEfl;
  real alpha2; 
  int flood_id;
  t_edlocals loc;
  t_eigvec      vecs;          /* use flooding for these               */
} t_edflood;

typedef struct t_ed_local* p_ed_local; 
/* handle for structure of local variables cannot be accessed outside edsam.c */

struct _t_edpar {
  int 		nini;		/* Total Nr of atoms    		*/
  int           ned;            /* Nr of atoms in essdyn                */
  bool          fitmas;         /* true if trans fit with cm            */
  bool          pcamas;
  int           presteps;       /* number of steps to run without any perturbations ... just monitoring */
  int           outfrq;         /* freq (in steps) of writing output    */
  int           logfrq;         /* freq (in steps) of writing to log    */
  int           maxedsteps;     /* max nr of steps per cycle            */
  t_edx         sref;           /* reference positions                  */
  t_edx         sav;            /* average positions                    */
  t_edvecs      vecs;           /* eigenvectors                         */
  real          slope;          /* minimal slope in acceptance radexp   */
  t_edx         star;           /* target positions                     */
  t_edx         sori;           /* origin positions                     */
  int           nmass;          /* Nr of masses                         */
  int           *masnrs;        /* index nrs for atoms with masses      */
  real          *mass;          /* atomic masses                        */
  real          tmass;          /* total mass                           */
  int           nfit;           /* Number of atoms to use for rot fit   */
  int           *fitnrs;        /* index nrs of atoms to use for rot fit*/
  bool          bNeedDoEdsam;    /* if any of the options mon,linfix,..,radcon is used (i.e. anything apart from flood) */
  t_edflood     flood;          /* parameters especially for flooding */
  FILE          *edo;           /* output file                          */
  p_ed_local  local;          /* handle to local buffer */
  rvec        *x_unc;
  struct _t_edpar       *next_edi;
};

typedef struct _t_edpar t_edpar;

typedef struct {
  bool 		bEdsam;		/* Do ED sampling?			*/
  char          *edinam; 	/* name of ED sampling input file       */
  char          *edonam;        /*                     output           */
  t_edpar       *edpar;
} t_edsamyn;


#endif


