/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
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
/* This file is completely threadsafe - keep it that way! */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "network.h"
#include "block_tx.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "main.h"
        
void _blocktx(int dest,int nelem,int size,void *data)
{
  int i;
  char *buf=data;
  
  if ((data==NULL) && (size > 0))
    gmx_fatal(FARGS,"TX: Null pointer (size=%d)!\n",size);

  for (i=0; i<nelem; i++) {
    gmx_txs(dest,&size,sizeof(size));
    if (size > 0)
      gmx_txs(dest,buf,size);
    buf+=size;      
  }
}

void _blockrx(int src,int nelem,int size,void *data)
{
  int i,len;
  char *buf=data;
  
  if ((data==NULL) && (size > 0))
    gmx_fatal(FARGS,"RX: Null pointer (size=%d)!\n",size);
  for (i=0; i<nelem; i++) {
    gmx_rxs(src,&len,sizeof(len));
    if (size!=len)
      gmx_fatal(FARGS,"%d: size=%d, len=%d, rx_count=%d\n",
		  0,size,len,0);
    if (len > 0)
      gmx_rxs(src,buf,len);
    buf+=len;      
  }
}

void mv_block(int dest,t_block *block)
{
  nblocktx(dest,MAXNODES,block->multinr);
#ifdef DEBUG
  fprintf(stdlog,"mv multinr\n");
#endif
  blocktx(dest,block->nr);
#ifdef DEBUG
  fprintf(stdlog,"mv block->nr (%d)\n",block->nr);
#endif
  nblocktx(dest,block->nr+1,block->index);
#ifdef DEBUG
  fprintf(stdlog,"mv block->index\n");
#endif
  blocktx(dest,block->nra);
#ifdef DEBUG
  fprintf(stdlog,"mv block->nra (%d)\n",block->nra);
#endif
  nblocktx(dest,block->nra,block->a);
#ifdef DEBUG
  fprintf(stdlog,"mv block->a\n");
#endif
}

void ld_block(int src,t_block *block)
{
  nblockrx(src,MAXNODES,block->multinr);
#ifdef DEBUG
  fprintf(stdlog,"ld multinr\n");
#endif
  blockrx(src,block->nr);
#ifdef DEBUG
  fprintf(stdlog,"ld block->nr (%d)\n",block->nr);
#endif
  snew(block->index,block->nr+1);
#ifdef DEBUG
  fprintf(stdlog,"block->index=%x\n",block->index);
#endif
  nblockrx(src,block->nr+1,block->index);
#ifdef DEBUG
  fprintf(stdlog,"ld block->index\n");
#endif
  blockrx(src,block->nra);
#ifdef DEBUG
  fprintf(stdlog,"ld block->nra (%d)\n",block->nra);
#endif
  snew(block->a,block->nra);
#ifdef DEBUG
  fprintf(stdlog,"block->a=%x\n",block->a);
#endif
  nblockrx(src,block->nra,block->a);
#ifdef DEBUG
  fprintf(stdlog,"ld block->a\n");
#endif
}

