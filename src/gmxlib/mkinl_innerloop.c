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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is NOT threadsafe, but it is only used to create
 * the innerloops during the build process, so it will never be
 * executed by multiple threads.
 */

#include "mkinl.h"

void unpack_inner_data(bool calcdist,bool calcforce)
{
    int nflop = 0;
    
    assign("jnr", "%s%s", ARRAY(jjnr,k), bC ? "" : "+1");
    if(DO_FORCE || calcdist)
      assign("j3",bC ? "3*jnr" : "3*jnr-2");
    
    if(arch.vectorcpu && calcforce)    /* prefetch x is turned off in this case */
      assign("kk","%d*(k-nj0)%s",(loop.sol==SOL_WATERWATER) ? 9 : 3 ,bC ? "" : "+1");
}



void fetch_coord(void)
{
  int j,offset;

  comment("Fetching coordinates");

  for(j=1;j<=loop.nj;j++) {
    offset=3*(j-1);
    assign("jx%d",_array("pos", "j3+%d",offset),j);
    assign("jy%d",_array("pos", "j3+%d",offset+1),j);
    assign("jz%d",_array("pos", "j3+%d",offset+2),j);
  }
}
	

void prefetch_forces(void)
{

  int j;
  if(DO_FORCE) {  
    comment("Prefetching forces");
  
    for(j=1;j<=loop.nj;j++) {
      assign("fjx%d", _array("faction","j3+%d",3*(j-1)),j);
      assign("fjy%d", _array("faction","j3+%d",3*(j-1)+1),j);
      assign("fjz%d", _array("faction","j3+%d",3*(j-1)+2),j);
    }
  }
  /* no flops */
}



void unpack_vector_machine_forces(bool calcdist,bool calcforce)
{
  /* this shouldnt be used with prefetching */
  if(DO_FORCE) {
    comment("Unpack forces for vectorization on a vector machine");
    vector_pragma();
    start_loop("k","nj0","nj1");
    
    unpack_inner_data(calcdist,calcforce);  
    
    assign(ARRAY(fbuf,kk), ARRAY(faction,j3));
    assign(ARRAY(fbuf,kk+1), ARRAY(faction,j3+1));
    assign(ARRAY(fbuf,kk+2), ARRAY(faction,j3+2));
    
    if(loop.sol==SOL_WATERWATER) {
      assign(ARRAY(fbuf,kk+3), ARRAY(faction,j3+3));
      assign(ARRAY(fbuf,kk+4), ARRAY(faction,j3+4));
      assign(ARRAY(fbuf,kk+5), ARRAY(faction,j3+5));
      assign(ARRAY(fbuf,kk+6), ARRAY(faction,j3+6));
      assign(ARRAY(fbuf,kk+7), ARRAY(faction,j3+7));
      assign(ARRAY(fbuf,kk+8), ARRAY(faction,j3+8));
    }
    end_loop();
  }
}



void inner_loop(bool calcdist, bool calcforce)
{
  int nflop=0;
  char buf[50],buf2[50];
  
  if(arch.vectorcpu && calcforce)
    unpack_vector_machine_forces(calcdist,calcforce);
  
  comment("Inner loop (over j-particles) starts right here");
    
  vector_pragma();
  
  start_loop("k", "nj0" , "nj1");
  
  unpack_inner_data(calcdist,calcforce);

  fetch_coord();
  
  /* non-vectorized coordinate prefetching is issued from 
   * calc_interactions() in mkinl_interactions.c due to 
   * timing optimization.
   */

  if(calcdist)
    nflop += calc_dist(); 
  
  if(calcforce) { 
    
    nflop += calc_rinv_and_rinvsq();
    
    if(DO_PREFETCH) 
      prefetch_forces();

    nflop += calc_interactions();
    
    /* Only print vanilla loop count */
    if(!DO_VECTORIZE && !arch.vectorcpu) {
      sprintf(buf,"Innerloop of %s costs %d flops",loopname,nflop);
      comment(buf);
    }
  }
  end_loop();  
}
