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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  bool 		bEdsam;		/* Do ED sampling?			*/
  char          *edinam; 	/* name of ED sampling input file       */
  char          *edonam;        /*                     output           */
} t_edsamyn;

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
} t_edx;


typedef struct {
  int 		nini;		/* Total Nr of atoms    		*/
  int           npro;           /* Nr of protein atoms                  */
  int           ned;            /* Nr of atoms in essdyn                */
  bool          selmas;         /* true if trans fit with cm            */
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
  FILE          *edo;           /* output file                          */
} t_edpar;

