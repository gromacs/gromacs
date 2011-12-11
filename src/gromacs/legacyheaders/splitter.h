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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
 
#ifndef _splitter_h
#define _splitter_h

#include "typedefs.h"
#include "types/inputrec.h"

#ifdef __cplusplus
extern "C" {
#endif

void split_top(FILE *fp,int nnodes,gmx_localtop_t *top,
		      t_inputrec *ir,t_block *mols,
		      real *capacity,int *mulitnr_cgs,int **multinr_nre,
		      int *left_range, int *right_range);
/* Split the topology (blocks and forces, based on charge groups 
 * and shake blocks.
 * The capacity is releated to the capacity of each node. If all numbers are 
 * equal, load will be distributed equally. If not some (the higher ones)
 * will get more than others. The sum of capacities should be 1.
 * Info is written to the file pointer fp.
 */

void gen_sblocks(FILE *fp,int at_start,int at_end,
			t_idef *idef,t_blocka *sblock,
			gmx_bool bSettle);
/* Generate shake blocks from the constraint list. Set bSettle to yes for shake
 * blocks including settles. You normally do not want this.
 */

#ifdef __cplusplus
}
#endif

#endif
