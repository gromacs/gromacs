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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_block_tx_c = "$Id$";
#include "network.h"
#include "block_tx.h"
#include "fatal.h"
#include "smalloc.h"
#include "main.h"
#include "assert.h"
        
void _blocktx(int dest,int nelem,int size,void *data)
{
  int i;
  char *buf=data;
  
  if ((data==NULL) && (size > 0))
    fatal_error(0,"TX: Null pointer (size=%d)!\n",size);
#ifdef DEBUG
#ifdef _amb_
  if (size > 256000) 
    fprintf(stdlog,"Trying to send %d bytes!\n",size);
#endif
#endif

  for (i=0; i<nelem; i++)
    {
      gmx_txs(dest,record(size));
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
    fatal_error(0,"RX: Null pointer (size=%d)!\n",size);
#ifdef DEBUG
#ifdef _amb_
  if (size > 256000) 
    fprintf(stdlog,"Trying to receive %d bytes!\n",size);
#endif
#endif
  for (i=0; i<nelem; i++)
    {
      gmx_rxs(src,record(len));
      if (size!=len)
        fatal_error(0,"%d: size=%d, len=%d, rx_count=%d\n",
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

