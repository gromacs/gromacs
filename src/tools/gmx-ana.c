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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "copyrite.h"
#include "gmx-ana.h"

typedef int (*ana_func)(int argc,char *argv[]);
typedef struct {
  char *name;
  ana_func f;
} t_ana_func;

int main(int argc,char *argv[]) 
{
  int i,j,k; 

  t_ana_func af[] = {
    { "wham", gmx_wham },
    { "analyze", gmx_analyze },
    { "anaeig", gmx_anaeig },
    { "angle", gmx_angle },
    { "bond", gmx_bond },
    { "bundle", gmx_bundle },
    { "chi", gmx_chi },
    { "cluster", gmx_cluster },
    { "confrms", gmx_confrms },
    { "covar", gmx_covar },
    { "density", gmx_density },
    { "dielectric", gmx_dielectric },
    { "dih", gmx_dih },
    { "dipoles", gmx_dipoles },
    { "disre", gmx_disre },
    { "dist", gmx_dist },
    { "dyndom", gmx_dyndom },
    { "enemat", gmx_enemat },
    { "energy", gmx_energy },
    { "lie", gmx_lie },
    { "filter", gmx_filter },
    { "gyrate", gmx_gyrate },
    { "h2order", gmx_h2order },
    { "hbond", gmx_hbond },
    { "helix", gmx_helix },
    { "mindist", gmx_mindist },
    { "msd", gmx_msd },
    { "morph", gmx_morph },
    { "nmeig", gmx_nmeig },
    { "nmens", gmx_nmens },
    { "order", gmx_order },
    { "potential", gmx_potential },
    { "rama", gmx_rama },
    { "rdf", gmx_rdf },
    { "rms", gmx_rms },
    { "rmsdist", gmx_rmsdist },
    { "rmsf", gmx_rmsf },
    { "rotacf", gmx_rotacf },
    { "saltbr", gmx_saltbr },
    { "sas", gmx_sas },
    { "sgangle", gmx_sgangle },
    { "sorient", gmx_sorient },
    { "tcaf", gmx_tcaf },
    { "traj", gmx_traj },
    { "velacc", gmx_velacc },
    { "clustsize", gmx_clustsize },
    { "mdmat", gmx_mdmat }
  };

#define NAF asize(af)
  for(i=1; (i<argc-1); i++) {
      if (strcmp(argv[i],"-type") == 0) {
	  for(j=0; (j<NAF); j++) {
	      if (strcmp(af[j].name,argv[i+1]) == 0) {
		  for(k=i+2; (k<argc); k++)
		      argv[k-2] = argv[k];
		  return af[j].f(argc-2,argv);
	      }
	  }
	  if (j == NAF) {
	      fatal_error(0,"Invalid argument for -type; %s",argv[i+1]);
	  }
      }
  }
  CopyRight(stderr,argv[0]);

  fprintf(stderr,"Usage: %s -type [",argv[0]);
  for(i=0; (i<NAF-1); i++) 
    fprintf(stderr," %s |",af[i].name);
  fprintf(stderr," %s ]\n",af[i].name);
  
  return -1;
}
