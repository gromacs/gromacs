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

#ifndef _superb_h
#define _superb_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <typedefs.h>

typedef struct {
  t_block *grps;	/* The group members			*/
  char    **name;	/* The group names			*/
  int     *nrdf;	/* Nr of degrees of freedom in a group	*/
  real    *temp;	/* Coupling temperature	per group	*/
  rvec    *acc;		/* Acceleration per group		*/
  tensor  *ekin;	/* Array of energy tensors...		*/
  rvec	  *u;           /* Mean velocities of home particles    */
  atom_id *invgrp;      /* Group number for each atom           */
} t_superblock;

extern t_superblock *init_grps(FILE *log,int left,int right,int nodeid,int nnodes,
			       char *gfile,bool bMaster);
/* Read a superblock structure from gfile. Do communication if
 * necessary.
 */

#endif	/* _superb_h */
