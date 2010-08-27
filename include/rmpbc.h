/*
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

#ifndef _rmpbc_h
#define _rmpbc_h

#include "typedefs.h"
	
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct gmx_rmpbc *gmx_rmpbc_t;
  
  gmx_rmpbc_t gmx_rmpbc_init(t_idef *idef,int ePBC,int natoms,
				    matrix box);
  
  void gmx_rmpbc_done(gmx_rmpbc_t gpbc);

  void gmx_rmpbc(gmx_rmpbc_t gpbc,int natoms,matrix box,rvec x[]);
  /* Correct coordinates x for atoms within every molecule for the periodic
   * boundary conditions such that every molecule is whole.
   * natoms is the size x and can be smaller than the number 
   * of atoms in idef, but should only contain complete molecules.
   * When ePBC=-1, the type of pbc is guessed from the box matrix.
   */

  void gmx_rmpbc_copy(gmx_rmpbc_t gpbc,int natoms,matrix box,rvec x[],
			     rvec x_s[]);
  /* As gmx_rmpbc, but outputs in x_s and does not modify x. */

  void gmx_rmpbc_trxfr(gmx_rmpbc_t gpbc,t_trxframe *fr);
  /* As gmx_rmpbc but operates on a t_trxframe data structure. */

  /*void rm_pbc(t_idef *idef,int ePBC,int natoms,
    matrix box,rvec x[],rvec x_s[]);*/
  /* Convenience function that still holds a static variable. */
  
  void rm_gropbc(t_atoms *atoms,rvec x[],matrix box);
  /* Simple routine for use in analysis tools that just have a pdb or 
   * similar file.
   */
  
#ifdef __cplusplus
}
#endif
 
#endif	/* _rmpbc_h */
