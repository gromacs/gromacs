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
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _pull_h
#define _pull_h

static char *SRCID_pull_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec.h"
#include "typedefs.h"

#define NINT(x) (x<0?((int)((x)-0.5)):((int)((x)+0.5)))
#define DEBUG_START fprintf(stderr,"\n\nDEBUG\n");
#define DEBUG_END   fprintf(stderr,"\nEND DEBUG\n\n");

/* print to output file for the various types of runs */
extern void print_umbrella(t_pull *pull, int step);
extern void print_afm(t_pull *pull, int step);
extern void print_constraint(t_pull *pull,rvec *force,int step,matrix box,
			     int niter);
extern void print_start(t_pull *pull, int step);


/* calculate center of mass of index group, making sure it's inside the box.
   function returns total mass.  
*/
extern real calc_com(rvec x[],       /* coordinates of all atoms in system */ 
		     int gnx,        /* size of index group */
		     atom_id *index, /* indices of atoms to be used for com */
		     t_mdatoms *md,  /* all atoms */
		     rvec com,       /* calculated center of mass */
		     matrix box);    


/* calculate center of mass of all atoms x[], index needed to get the right
   masses from the atom array. function returns total mass.*/
extern real calc_com2(rvec x[],       /* coordinates to calc. com from */
		      int gnx,        /* nr. of atom in group  */
		      atom_id *index, /* indices of x[] in all atoms */
		      t_mdatoms *md,  /* all atoms */
		      rvec com,       /* calculated center of mass */
		      matrix box); 


/* calculate a running average for center of mass */
extern void calc_running_com(t_pull *pull);


/* calculate the center of mass from the true coordinates, without
   corrections for pbc */
extern void correct_t0_pbc(t_pull *pull, 
			   rvec x[], 
			   t_mdatoms *md,
			   matrix box);


/* parse a string for 3 numbers and put them in rvec */
extern void string2rvec(char *buf,
			rvec x);

/* read the parameter file .ppa and write out what was read in */
extern void read_pullparams(t_pull *pull,
			    char *infile,
			    char *outfile); 

/* find all atoms in group pull->idx[pull->n] that are inside a cylinder
   with as origin com[i][x],com[i][y] with radius pull->r and possibly
   a switch function pull->rc. Remember their weight. Now each group i
   has its own reference group (HOW?) with com defined as 
   Sum(wi*mi*ri)/(Sum(wi*mi). Basically, build an index structure with
   the reference groups for the groups i, plus an array with the 
   weight factors for each of the atoms in those index groups? 
   */
extern void make_refgrps(t_pull *pull,
			 matrix box,
			 t_mdatoms *md);


/* write a numbered .gro file in procedure to make starting structures */
extern void dump_conf(t_pull *pull,
		      rvec x[],        /* all coordinates */
		      matrix box,      /* box */
		      t_topology *top, /* names and residue info */
		      int nout,        /* sequence number of this file */
		      real time);      /* time in simulation */


/* main pull routine that controls all the action */
extern void pull(t_pull *pull,    /* all pull data */
		 rvec *x,         /* coordinates, changed by constraint run */ 
		 rvec *f,         /* forces, changed by Afm run */
		 matrix box,               
		 t_topology *top, /* needed to write out coordinate files */   
		 real dt,         /* time step */
		 int step,        /* step number in simulation */
		 int natoms,      /* total number of atoms on this processor */
		 t_mdatoms *md);  /* masses and charges of all atoms */


/* get memory and initialize the fields of pull that still need it, and
   do runtype specific initialization */
extern void init_pull(FILE *log,  
		      int nfile,       
		      t_filenm fnm[], /* standard filename struct */
		      t_pull *pull,   /* all pull data */
		      rvec *x,        /* all coordinates */
		      t_mdatoms *md,  /* masses and charges of all atoms */
		      matrix box);

#endif
