/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef real rvec5[5];

/* Distance restraining stuff */
typedef struct {
  int  dr_weighting;  /* Weighting of pairs in one restraint        */
  bool dr_bMixed;     /* Use sqrt of the instantaneous times        *
		       * the time averaged violation                */
  real dr_fc;	      /* Force constant for disres,                 *
		       * which is multiplied by a (possibly)        *
		       * different factor for each restraint        */
  real dr_tau;	      /* Time constant for disres		    */
  real ETerm;         /* multiplication factor for time averaging   */
  real ETerm1;        /* 1 - ETerm1                                 */
  real exp_min_t_tau; /* Factor for slowly switching on the force   */
  int  npr;           /* The number of distance restraint pairs     */
  real sumviol;       /* The sum of violations                      */
  real *rav;          /* The calculated time averaged distance      */
  real *rt;           /* The calculated instantaneous distance      */
} t_disresdata;

/* Orientation restraining stuff */
typedef struct {
  real   fc;          /* Force constant for the restraints                 */
  real   edt;         /* Multiplication factor for time averaging          */
  real   edt1;        /* 1 - edt                                           */
  real   exp_min_t_tau; /* Factor for slowly switching on the force        */
  int    nr;          /* The number of orientation restraints              */
  int    nex;         /* The number of experiments                         */
  tensor *S;          /* Array of order tensors for each expiriment (nexp) */
  rvec5  *Dins;       /* The order matrix for all restraints (nr x 5)      */
  rvec5  *Dtav;       /* The order matrix for all restraints (nr x 5) tav. */
  real   *oinsl;      /* The calculated instantaneous orientations         */
  real   *oins;       /* The calculated emsemble averaged orientations     */
  real   *otav;       /* The calculated time and ensemble averaged orient. */
  real   rmsviol;     /* The weighted (using kfac) RMS violation           */
  rvec5  *tmp;        /* An array of temporary 5-vectors (nex);            */ 
  real   ***TMP;      /* An array of temporary 5x5 matrices (nex);         */
} t_oriresdata;

/* 
 * Data struct used in the force calculation routines
 * for storing information which is needed in following steps
 * (for instance for time averaging in distance retraints)
 * or for storing output, since force routines only return the potential.
 */
typedef struct {
  t_disresdata disres;
  t_oriresdata orires;
} t_fcdata;
