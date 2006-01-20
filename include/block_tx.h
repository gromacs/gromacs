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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern void _blocktx(const t_commrec *cr,int dest,
		     int nelem,int size,void *data);
extern void _blockrx(const t_commrec *cr,int src,
		     int nelem,int size,void *data);

#define blocktx(cr,dest,dta)	_blocktx((cr),(dest),1,sizeof(dta),&(dta))
#define blockrx(cr,src,dta)	_blockrx((cr),(src),1,sizeof(dta),&(dta))
#define nblocktx(cr,dest,n,dta)	_blocktx((cr),(dest),1,(n)*(sizeof(*dta)),dta)
#define nblockrx(cr,src,n,dta)	_blockrx((cr),(src),1,(n)*(sizeof(*dta)),(dta))

extern void mv_block(const t_commrec *cr,int dest,t_block *block);
extern void ld_block(const t_commrec *cr,int src,t_block *block);
/* Send and receive an t_block type */
