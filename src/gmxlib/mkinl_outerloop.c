/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This file is NOT threadsafe, but it is only used to create
 * the innerloops during the build process, so it will never be
 * executed by multiple threads.
 */

#include <string.h>
#include "mkinl.h"

/* Routines to start/end various types of loops */

void start_thread_loop(void)
{
  int i;
  char buf1[100];
    comment("Thread loop start. Get limits for the outer loop");

    if(bC) {
      code("do {");
      IND += 2;
      code("pthread_mutex_lock(mtx);");
      assign("nn0","*count");
      assign("nn1","nn0+(nri-nn0)/(2*nthreads)+3"); 
      /* take sucessively smaller chunks */
      assign("*count","nn1");
      code("pthread_mutex_unlock(mtx);");
      code("if(nn1>nri) nn1=nri;");
    } else {
      char space[25];
      for(i=0;i<IND-6;i++)
	space[i]=' ';
      space[i]=0;
      
      sprintf(buf1,"   10 %scall inlsync(nri,nthreads,count,nn0,nn1,mtx)\n",space);
      strcat(codebuffer,buf1);
      /* since f77 use call-by-reference we can send the pointer of the
       * count variable and the mutex to c without fortran knowing about it!
       */
      IND += 2;
      code("if(nn1.gt.nri) nn1=nri");
    }
}


void end_thread_loop(void)
{
  IND -= 2;
  if(bC) 
    code("} while (nn1<nri);"); 
  else
    code("if(nn1.lt.nri) goto 10");
}



void start_vectorized_calcdist_loop(void)
{
  switch(arch.threads) {
  case THREAD_MUTEX:
    newline();
    code( bC ? "while(nn0<nn1) {" : "do while(nn0.le.nn1)");
   IND += 2;
    assign("m", bC ? "0" : "1");
    assign("m3", bC ? "0" : "1");
    start_loop("n","nn0","nn1");
    break;
    
  case THREAD_STRIDE:
    assign("nn0", bC ? "mythread" : "1+mythread");
    newline();
    code( bC ? "while(nn0<nri) {" : "do while(nn0.le.nri)");
    IND += 2;
    assign("m", bC ? "0" : "1");
    assign("m3", bC ? "0" : "1");
    start_stride_loop("n","nn0","nri","nthreads");
    break;
    
  default: /* no threads */
    assign("nn0", bC ? "0" : "1");
    code( bC ? "while(nn0<nri) {" : "do while(nn0.le.nri)");
    IND += 2;    
    assign("m", bC ? "0" : "1");
    assign("m3", bC ? "0" : "1");
    start_loop("n","nn0","nri");
    break;
  }
}



void end_vectorized_calcdist_loop(void)
{
  char buf1[100];
  end_loop(); /* close the for loop */
  
  if(bC) 
    assign("nn2","n");
  else {
    char space[25];
    int i;
    
    for(i=0;i<25;i++)
      space[i]=' ';
    space[IND-6]=0;
    sprintf(buf1,"   20 %snn2  =  n - 1\n",space);   
    strcat(codebuffer,buf1);
  }
}



void start_vectorized_calcforce_loop(void)
{
  if(arch.threads==THREAD_STRIDE)
    start_stride_loop("n","nn0","nn2","nthreads");
  else
    start_loop("n","nn0","nn2");
}


void end_vectorized_calcforce_loop(void)
{
  end_loop(); /* close the 2nd nlist for loop */

  assign("nn0", bC ? "nn2" : "nn2+1");

  end_loop(); /* close the do-while before the first loop */
}


/* Data unpacking/updating routines */
void unpack_outer_indices(bool calcdist)
{
  /* This routine calculates the length of the neighbourlist
   * and unpacks the shift vector. It is only called once for
   * each outerloop iteration. (Even if we have several I
   * particles)
   */

  if(calcdist) {
    comment("Unpack shift vector");
    assign("is3","3*%s%s",ARRAY(shift,n), bC ? "" : "+1");  

    assign("shX",ARRAY(shiftvec,is3));
    assign("shY",ARRAY(shiftvec,is3+1));
    assign("shZ",ARRAY(shiftvec,is3+2));
  }

  if(calcdist || DO_VECTORIZE) {
    assign("ii", bC ? "%s" : "%s+1",ARRAY(iinr,n));
    if(DO_FORCE || calcdist)
      assign("ii3", bC ? "3*ii" : "3*ii-2");
  }
  
  if(!(calcdist && DO_VECTORIZE)) {
    comment("Local variables for energy");
    if(loop.coul)
      assign("vctot","nul");
    if (loop.vdw) 
      assign("vnbtot","nul");
  }

  if(calcdist) {
    comment("Bounds for the innerloop");
    assign("nj0","%s%s",ARRAY(jindex,n), bC ? "" : "+1");
    assign("nj1",ARRAY(jindex,n+1));
  }
/* no flops */
}

  
int unpack_outer_data(bool calcdist, bool calcforce)
{
  /* This routine calculates coordinates and other data
   * for an i particle. It is called once before each
   * innerloop iteration. For water loops we assign three
   * atoms in parallel. For simplewater or general solvent
   * loops we loop over this and the innerloop.
   */
  int nflop = 0;
  int i;

  comment("Unpack I particle");  
  
  if(calcdist)
    for(i=1;i<=loop.ni;i++) {
      add("ix%d","shX", _array("pos","ii3+%d",3*(i-1)),i);
      add("iy%d","shY", _array("pos","ii3+%d",3*(i-1)+1),i);
      add("iz%d","shZ", _array("pos","ii3+%d",3*(i-1)+2),i);
      nflop += 3;
    }

  if(calcforce) {
    /* Need charges */
    if(loop.coul && (!DO_WATER || arch.simplewater)) {
      comment("Charge of i particle divided by 4 pi eps0");
      assign("iqA","facel*%s",  ARRAY(charge,ii));
      nflop++;
      if (loop.free) {
	assign("iqB","facel*%s",ARRAY(chargeB,ii));
	nflop++;
      }
    }
    /* Need VDW parameters */
    if(!DO_WATER && (loop.vdw || DO_SOFTCORE)) {
      comment("Atom types for Van der Waals interactions");
     
      assign("ntiA","%d*ntype*%s",N_VDWPARAM, ARRAY(type,ii));
      if (loop.free) 
	assign("ntiB","%d*ntype*%s",N_VDWPARAM, ARRAY(typeB,ii));
    }    

    /* zero local i forces */
    if(DO_FORCE) 
      for(i=1;i<=loop.ni;i++) {
	assign("fix%d","nul",i);
	assign("fiy%d","nul",i);
	assign("fiz%d","nul",i);
      }  
  }
  return nflop;
}


int update_outer_data(void)
{
  int nflop = 0;
  int i;
  char buf[3][256];

  if(DO_FORCE) {
    comment("Update forces on i particles");
    buf[XX][0]=buf[YY][0]=buf[ZZ][0]=0;

    for(i=1;i<=loop.ni;i++) {
      increment(_array("faction","ii3+%d",3*(i-1)),"fix%d",i);
      increment(_array("faction","ii3+%d",3*(i-1)+1),"fiy%d",i);
      increment(_array("faction","ii3+%d",3*(i-1)+2),"fiz%d",i);
      
      sprintf(buf[XX]+strlen(buf[XX]),"%sfix%d", (i==1) ? "" : "+",i);
      sprintf(buf[YY]+strlen(buf[YY]),"%sfiy%d", (i==1) ? "" : "+",i);
      sprintf(buf[ZZ]+strlen(buf[ZZ]),"%sfiz%d", (i==1) ? "" : "+",i);
      nflop += 6;
    }
    increment(ARRAY(fshift,is3),  buf[XX]);
    increment(ARRAY(fshift,is3+1),buf[YY]);
    increment(ARRAY(fshift,is3+2),buf[ZZ]);  
    /* we add 3 flops too many above, so dont do it here */
  }
  return nflop;    
}


void outer_loop()
{
  vdw_t vdwsave;
  coul_t coulsave;
  bool dodist,doforce;
  char buf[100];
  int nflop = 0;
  
  /* When we are doing explicit syncing with mutexes we have an
   * outermost loop where we allocate chunks of neighbourlists
   * to the threads 
   */
  if(arch.threads==THREAD_MUTEX)
    start_thread_loop();

  /* start the outer loop */
  if(arch.threads==THREAD_STRIDE)
    start_stride_loop("n", bC ? "mythread" : "mythread+1","nri","nthreads");
  else
    start_loop("n",bC ? "0" : "1","nri");
  
  if(loop.sol==SOL_MNO) {
    assign("nstot", _array("nsatoms" , "3*n%s", bC ? "" : "-2"));
    assign("nsvdwc", _array("nsatoms" , "3*n%s", bC ? "+1" : "-1"));
    assign("nscoul", _array("nsatoms" , "3*n%s", bC ? "+2" : ""));
  }
  
  if(DO_VECTORIZE) {
    /* When using vectorized routines we first do an innerloop and
     * calc the distance, then repeat it and do the force.
     */
    dodist=TRUE;
    doforce=FALSE;     /* only distance calculation */
    
    unpack_outer_indices(dodist);
    
    /* For simplewater or general solvent we loop over the
     * i particles in the molecule (all atoms for calcdist)
     */
    assign("m", bC ? "0" : "1");
    if(DO_FORCE)
      assign("m3", bC ? "0" : "1");

    if(loop.sol==SOL_MNO) {
      if(loop.vdw && loop.coul)
	start_loop("s", bC ? "0" : "1", "nstot");
      else if(loop.coul && !loop.vdw)
	start_loop("s", bC ? "0" : "1", "nscoul");
      else if(!loop.coul && loop.vdw)
	start_loop("s",bC ? "0" : "1", "nsvdwc");
    } else if(DO_WATER && arch.simplewater)
      start_loop("s",bC ? "0" : "1", "3");
    
    nflop += unpack_outer_data(dodist,doforce);
    
    inner_loop(dodist,doforce);
    
    if(loop.sol==SOL_MNO || (DO_WATER && arch.simplewater)) {
      increment("ii","1");
      increment("ii3","3");
      end_loop();
    }
    
    /* and maybe a second iteration over vdw-only atoms - 
     * if we are skipping the coul-only ones 
     */
    if(loop.sol==SOL_MNO && !loop.coul && loop.vdw) {
      /* skip some indices */
      increment("ii","nscoul-nsvdwc");
      assign("ii3", bC ? "3*ii" : "3*ii-2");
      start_loop("s", bC ? "nscoul" : "nscoul+1" , "nstot" );
      nflop += unpack_outer_data(dodist,doforce);
      
      inner_loop(dodist,doforce);
      increment("ii","1");
      increment("ii3","3");
      end_loop();
    }
    
    comment("Perform vectorized calculations");
    call_vectorized_routines();
    
    assign("m", bC ? "0" : "1");
    assign("m3", bC ? "0" : "1");
    
    comment("Starting outer loop to calc forces");
  } else 
    comment("Outer loop starts here");
    
  
  /* Now we calc the force. Unless we used the vectorization above
   * we also need to calculate the distance
   */
  
  doforce=TRUE;
  dodist=!DO_VECTORIZE;
    
  unpack_outer_indices(dodist);
  
  /* First general solvent loop - vdw and coulomb , or coulomb/vdw
   * only over the first atoms with both interactions if the other one is disabled */
  if(loop.sol==SOL_MNO) {
    if(loop.vdw) /* maybe also coulomb */
      start_loop("s", bC ? "0" : "1" , "nsvdwc");
    else if(loop.coul) /* coul only */
      start_loop("s", bC ? "0" : "1" , "nscoul");
  } else if(DO_WATER && arch.simplewater) {
    /* coul and vdw part of simplewater loop */
    if(!loop.vdw)
      start_loop("s",bC ? "0" : "1", "3");
  }
  
  /* This is also the core for the non-water/solvent innerloop */
  nflop += unpack_outer_data(dodist,doforce);
  inner_loop(dodist,doforce);
  nflop += update_outer_data();
  
  if(loop.sol==SOL_MNO || (DO_WATER && arch.simplewater && !loop.vdw)) {
    /* only updated when we are doing a for loop */
    increment("ii","1");
    increment("ii3","3");
    end_loop();
  }

  /* Second part of general solvent or simplewater - loop over
   * coul-only atoms. However, if we didnt do vdw above we looped
   * over all coulomb atoms already there.
   */
  if(loop.vdw && (loop.sol==SOL_MNO || (DO_WATER && arch.simplewater))) {
    if(loop.coul) {
      /* dont do VDW on coul-only atoms */
      vdwsave=loop.vdw;
      loop.vdw=VDW_NO;
      if(loop.sol==SOL_MNO) 
      	start_loop("s", bC ? "nsvdwc" : "nsvdwc+1" , "nscoul" );
      else {/* simplewater */
	/* increment indices here, since we didnt loop over the single vdwc atom */
	increment("ii","1");
	increment("ii3","3");
	start_loop("s", bC ? "1" : "2" , "3" );
      }
      nflop += unpack_outer_data(dodist,doforce);
      inner_loop(dodist,doforce);
      nflop += update_outer_data();
      loop.vdw=vdwsave;
      increment("ii","1");
      increment("ii3","3");
      end_loop();
    }
  }
  
  /* And, finally, for general solvent loops - do a possible third
   * iteration over vdw only atoms
   */
  if(loop.sol==SOL_MNO && loop.vdw) {
    /* if we are not doing coul we might have to skip some coul-only
     * atoms. Update the indices accordingly!
     */
    if(!loop.coul) {
      increment("ii","nscoul-nsvdwc");
      assign("ii3", bC ? "3*ii" : "3*ii-2");
    }
    start_loop("s", bC ? "nscoul" : "nscoul+1" , "nstot" );
    
    coulsave=loop.coul;
    loop.coul=COUL_NO;
    loop.invsqrt=loop.vectorize_invsqrt ||
      ((loop.coul && (loop.coul_needs_rinv || loop.coul_needs_r)) ||
       (loop.vdw && (loop.vdw_needs_rinv || loop.vdw_needs_r)));
    loop.recip=!loop.invsqrt;

    nflop += unpack_outer_data(dodist,doforce);
    inner_loop(dodist,doforce);
    nflop += update_outer_data();
    
    loop.coul=coulsave;
    loop.invsqrt=
      (loop.coul && (loop.coul_needs_rinv || loop.coul_needs_r)) ||
      (loop.vdw && (loop.vdw_needs_rinv || loop.vdw_needs_r));
    loop.recip=!loop.invsqrt;
    
    increment("ii","1");
    increment("ii3","3");
    end_loop();

  }
  
  if(arch.threads==THREAD_MUTEX)
    end_thread_loop();
  
  comment("Update energies");  
  assign("ggid","%s%s", ARRAY(gid,n), bC ? "" : "+1");
  
  if(loop.coul) 
    increment(ARRAY(Vc,ggid),"vctot");  
  if (loop.vdw) 
    increment(ARRAY(Vnb,ggid),"vnbtot");
  
  end_loop();  
  
  if (loop.free) {
    increment("%sdvdlambda","dvdl", bC ? "*" : "");
    nflop ++;
  }
  /* Only print vanilla loop count */
  if(!DO_VECTORIZE && !arch.vectorcpu) {
    sprintf(buf,"Outerloop of %s costs %d flops",loopname,nflop);
    comment(buf);    
  }
}
