/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "gmx_fatal.h"
#include "macros.h"
#include "gmx_qhop_db.h"

int main(int argc,char *argv[])
{
  gmx_qhop_db db;
  char *donors[] = { "H3O+", "ACE" };
  char *acceptors[] = { "H2O", "GLU" };
  t_qhop_parameters qp;
  int i,j;
  
  if ((db = gmx_qhop_db_read("ffoplsaa")) == NULL) 
    gmx_fatal(FARGS,"Can not read qhop database information");
  if (gmx_qhop_db_write("koe.dat",db) != 1)
    gmx_fatal(FARGS,"Can not write qhop database information");
  
  for(i=0; (i<asize(donors)); i++) {
    for(j=0; (j<asize(acceptors)); j++) {
      if (gmx_qhop_db_get_parameters(db,donors[i],acceptors[j],&qp) == 1) {
	printf("Found qhop parameters for donor %s and acceptor %s\n",
	       donors[i],acceptors[j]);
      }
      else {
	printf("Could not find qhop parameters for donor %s and acceptor %s\n",
	       donors[i],acceptors[j]);
      }
    }
  }
  
  if (gmx_qhop_db_done(db) != 1)
    gmx_fatal(FARGS,"Error destroying qhop data");
    
  return 0;
}
