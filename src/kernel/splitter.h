/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_splitter_h = "$Id$";

extern void split_top(bool bVerbose,int nprocs,t_topology *top,real *capacity);
/* Split the topology (blocks and forces, based on charge groups 
 * and shake blocks.
 * The capacity is releated to the capacity of each node. If all numbers are 
 * equal, load will be distributed equally. If not some (the higher ones)
 * will get more than others. The sum of capacities should be 1.
 */

extern void gen_sblocks(bool bVerbose,int natoms,t_idef *idef,t_block *sblock,
			bool bSettle);
/* Generate shake blocks from the constraint list. Set bSettle to yes for shake
 * blocks including settles. You normally do not want this.
 */
