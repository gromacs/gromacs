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

#ifndef _init_h
#define _init_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "mvdata.h"
#include "nsb.h"
#include "vsite.h"

typedef enum
{
  LIST_SCALARS	=0001,
  LIST_PARM	=0002,
  LIST_TOP	=0004,
  LIST_X	=0010,
  LIST_V	=0020,
  LIST_F	=0040,
  LIST_LOAD	=0100
} t_listitem;

extern void check_nnodes_top(char *fn,t_topology *top,int nnodes);
/* Verify whether this tpr file is for nnodes nodes, and quit if not */

extern void init_single(FILE *log,
                        t_parm *parm, char *tpbfile, t_topology *top,
			rvec **x,rvec **v,t_mdatoms **mdatoms,
			t_nsborder *nsb);
     /*
      * Allocates space for the topology (top), the coordinates x, the
      * velocities v, masses mass. Reads the parameters, topology,
      * coordinates and velocities from the file specified in tpbfile
      */

extern void distribute_parts(int left,int right,int pid,int nprocs,
                             t_parm *parm,char *tpbfile,int nstDlb);
     /*
      * Reads the parameters, topology, coordinates and velocities for the
      * multi processor version of the program from the file specified in
      * parm->files[STATUS_NM]. This file should also contain a so called
      * split descriptor which describes how to distribute particles over
      * the system. It then selects for all subsystems the appropriate data
      * and sends this to the processor using the left and right channels.
      * At last it sends its own subsystem down the ring where it is buffered.
      * Its own buffers for reading the data from the file are freed, and it
      * is now possible to reload this processor from the ring by using the
      * init_parts() routine.
      * The routine also creates a renum array which can be used for writing
      * out the x,v and f for analysis purpose.
      */

extern void init_parts(FILE *log,t_commrec *cr,
		       t_parm *parm,t_topology *top,
		       rvec **x,rvec **v,t_mdatoms **mdatoms,
		       t_nsborder *nsb,int list,
		       bool *bParallelVsites,
		       t_comm_vsites *vsitecomm);
     /*
      * Loads the data for a simulation from the ring. Parameters, topology
      * coordinates, velocities, and masses are initialised equal to using
      * init_single() in the single processor version. The extra argument
      * f_add is allocated to use for the update of the forces, the load
      * array specifies in which part of the x and f array the subsystems
      * of the other processors are located. Homenr0, homenr1, nparts0 and
      * nparts1 are necessary to calculate the non bonded interaction using
      * the symmetry and thus calculating every force only once. List is a facility
      * for logging (and debugging). One can decide to print none or a set of
      * selected parameters to the file specified by log. Parameters are
      * printed by or-ing the corresponding items from t_listitem. A 0 (zero)
      * specifies that nothing is to be printed on the file. The function
      * returns the number of shifts over the ring to perform to calculate
      * all interactions.
      */

extern void write_parm(FILE *log,char *title,int pid,t_parm *parm);
/* Write parm for debugging */

#endif	/* _init_h */
