/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_block_tx_h = "$Id$";

        
extern void _blocktx(int dest,int nelem,int size,void *data);
extern void _blockrx(int src,int nelem,int size,void *data);

#define blocktx(dest,dta)	_blocktx((dest),1,sizeof(dta),&(dta))
#define blockrx(src,dta)	_blockrx((src),1,sizeof(dta),&(dta))
#define nblocktx(dest,n,dta)	_blocktx((dest),1,(n)*(sizeof(*dta)),dta)
#define nblockrx(src,n,dta)	_blockrx((src),1,(n)*(sizeof(*dta)),(dta))

extern void mv_block(int dest,t_block *block);
extern void ld_block(int src,t_block *block);
/* Send and receive an t_block type */
