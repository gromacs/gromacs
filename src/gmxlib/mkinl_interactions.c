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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "mkinl.h"
#include <string.h>

char forcebuf[255];
char tabforcebuf[255];


int calc_scalar_force(int i,int j)
{
  int ij = 10*i+j;
  int nflop = 0;

  /* calculate scalar force */
  if(!DO_SOFTCORE) {
    /* For softcore loops fs is calculated directly in the routine,
     * due to the special construction with rinva and rivb. We also
     * fix the extra flop added for each buffer. */
    comment("Calculate the scalar force");
    if(strlen(tabforcebuf) && strlen(forcebuf)) {
      /* both table and non-table interactions */
      assign("fs%d","((%s)*rinv%d-(%s))*rinv%d",ij,forcebuf,ij,tabforcebuf,ij);
      nflop ++; /* 3-2 flops */
    }
    else if(strlen(tabforcebuf))
      /* only table interactions */
      assign("fs%d","-(%s)*rinv%d",ij,tabforcebuf,ij); /* 1-1 flops */
    else if(strlen(forcebuf))
      /* only non-table interactions */
      assign("fs%d","(%s)*rinvsq%d",ij,forcebuf,ij); /* 1-1 flops */
  }
  return nflop;
}


int update_inner_forces(int i,int j)
{
  int ij = 10*i+j;
  int nflop = 0;
  int m,jidxnr,offset;
  char src[25],dest[25];


  comment("Convert scalar force to cartesian coords");
  if(DO_VECTORIZE) {
    assign("tx%d","%s*fs%d",ij,ARRAY(drbuf,m3),ij);
    increment("m3","1");
    assign("ty%d","%s*fs%d",ij,ARRAY(drbuf,m3),ij);
    increment("m3","1");
    assign("tz%d","%s*fs%d",ij,ARRAY(drbuf,m3),ij);
    increment("m3","1");
  } else {
    assign("tx%d","dx%d*fs%d",ij,ij,ij);
    assign("ty%d","dy%d*fs%d",ij,ij,ij);
    assign("tz%d","dz%d*fs%d",ij,ij,ij);
  }

  /* i forces */
  increment("fix%d","tx%d",i,ij);
  increment("fiy%d","ty%d",i,ij);
  increment("fiz%d","tz%d",i,ij);
  nflop += 6;

  /* and finally, j forces */
  for(m=0;m<DIM;m++) {
    
    offset = 3*(j-1)+m;
    
    /* determine source array or variable */
    if(!DO_PREFETCH && i==1) /* source is faction or fbuf */      
      sprintf(src, arch.vectorcpu ? _array("fbuf","kk+%d",offset) : _array("faction","j3+%d",offset));
    else  /* source is fj */
      sprintf(src,"fj%c%d",m+'x',j);
    
    /* determine destination array or variable */
    if(i==loop.ni) /* dest is faction or fbuf */
      sprintf(dest, _array("faction", "j3+%d",offset));
    else  /* dest is fj */
      sprintf(dest,"fj%c%d",m+'x',j);
    
    assign(dest,"%s-t%c%d%d",src,m+'x',i,j); /* dest=src-tx */
    nflop++;
  }
  return nflop;
}


int table_index(char *rin)
{
    assign("rt",rin);
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined DOUBLE && defined USE_X86TRUNC)
    if(bC) 
      code("x86trunc(rt,n0);");
    else
#endif
      assign("n0","rt"); /* also used for g77 */
    
    assign("eps","rt-n0");
    assign("eps2","eps*eps");  
    return 2;
}


int extract_table(char *pos)
{
    char buf[25];
    int nflop=0;

    assign("nnn",pos);
    assign("Y",ARRAY(VFtab,nnn));
    assign("F",ARRAY(VFtab,nnn+1)); 
    sprintf(buf,"eps*%s",ARRAY(VFtab,nnn+2));
    assign("Geps",buf);
    sprintf(buf,"eps2*%s",ARRAY(VFtab,nnn+3));
    assign("Heps2",buf);
    assign("Fp","F+Geps+Heps2"); 
    assign("VV","Y+eps*Fp"); 

    nflop = 6;

    if(DO_FORCE || DO_SOFTCORE) {
      assign("FF","Fp+Geps+two*Heps2");
      nflop += 3;
    }
    return nflop;
}



static int assign_chargeproduct(int i,int j,char *qq)
{
  /* Assign and return a string with the name of 
   * the variable containing the i charge * j charge. */
  int nflop = 0;

  if(loop.free==FREE_LAMBDA) {
    assign("qqA","iqA*%s",ARRAY(charge,jnr));
    assign("qqB","iqB*%s",ARRAY(chargeB,jnr));
    sprintf(qq,"qqq");
    nflop+=2;
  } else if(loop.free==FREE_SOFTCORE) {
    sprintf(qq,"iq*jq");
    nflop++;
  } else {
    switch(loop.sol) {
    case SOL_WATER:
      assign("qq","%s*%s", arch.simplewater ? "iqA" : ((i==1) ? "qO" : "qH"),ARRAY(charge,jnr));
      nflop += 1;
      sprintf(qq,"qq");
      break;
    case SOL_WATERWATER:
      if(arch.simplewater) {
	  sprintf(qq,"iqA*%s", (j==1) ? "qO" : "qH");
      } else {
      if(i==1 && j==1)
	sprintf(qq,"qqOO");
      else if(i==1 || j==1)
	sprintf(qq,"qqOH");
      else
	sprintf(qq,"qqHH");
      }
      break;
    default:
      assign("qq","iqA*%s",ARRAY(charge,jnr)); 
      nflop += 1;
      sprintf(qq,"qq");
    }
  }
  return nflop;
}


int do_softcore(int i, int j)
{
    int nflop = 0;
    char fabuf[50],fbbuf[50],dvdlbuf[100];
    char ifbuf[100];
    int ij=10*i+j;
    
    fabuf[0]=fbbuf[0]=dvdlbuf[0]=0;
    
    comment("Tabulated Softcore free energy interactions");    

    assign("tjA","ntiA+%d*%s%s", N_VDWPARAM, ARRAY(type,jnr) , bC ? "" : "+1");
    assign("tjB","ntiB+%d*%s%s", N_VDWPARAM, ARRAY(typeB,jnr) , bC ? "" : "+1");
    assign("c6a",ARRAY(nbfp,tjA));
    assign("c6b",ARRAY(nbfp,tjB));
    assign( DO_BHAM ? "cexp1a" : "c12a",ARRAY(nbfp,tjA+1));
    assign( DO_BHAM ? "cexp1b" : "c12b",ARRAY(nbfp,tjB+1));
    if(DO_BHAM) {
      assign("cexp2a",ARRAY(nbfp,tjA+2));
      assign("cexp2b",ARRAY(nbfp,tjB+2));
      assign("sigma6a","defsigma6");
      assign("sigma6b","defsigma6");       
    } 

    if(!DO_BHAM) {
      start_if( bC ? "(c6a > 0) && (c12a > 0)" : 
		"(c6a.gt.0).and.(c12a.gt.0)");
      /* check if we can compute sigma_a from c6/c12 */
      assign("sigma6a","c12a/c6a");
      do_else(); /* use default value for sigma_a */
    }
    assign("sigma6a","defsigma6");
    if(!DO_BHAM) {
      end_if(); /* end of check for sigma_a */
    
      start_if( bC ? "(c6b > 0) && (c12b > 0)" : 
		"(c6b.gt.0).and.(c12b.gt.0)");
      /* corresponding check for sigma_b */
      assign("sigma6b","c12b/c6b");
      do_else();
    }
    assign("sigma6b","defsigma6");
    if(!DO_BHAM)
      end_if(); /* end of sigma_b check */
    
    assign("rfour","rsq%d*rsq%d",ij,ij);
    assign("rsix","rfour*rsq%d",ij);

    ifbuf[0]=0;
    
    if(DO_COULTAB) {
      assign("qqA","iqA*%s",ARRAY(charge,jnr));
      sprintf(ifbuf, bC ? "(qqA != 0)" : "(qqA.ne.0)");
    }
    if(DO_VDWTAB) {
      if(strlen(ifbuf)>0)
	strcat(ifbuf, bC ? " || " : ".or.");
      strcat(ifbuf, bC ? "(c6a > 0) || " : "(c6a.gt.0).or.");
      if(DO_BHAM)
	strcat(ifbuf, bC ? "(cexp1a > 0)" : "(cexp1a.gt.0)");
      else
	strcat(ifbuf, bC ? "(c12a > 0)" : "(c12a.gt.0)");
    }
    
    start_if(ifbuf);
    /* do state A */
    assign("rA", bC ? "pow(Alpha*sigma6a*lam2+rsix,onesixth)" :
	   "cpow(Alpha*sigma6a*lam2+rsix,onesixth)");
    assign("rinva","1.0/rA");
    assign("rinv5a","rinva*rinva");
    assign("rinv5a","rinv5a*rinv5a*rinva");
     
    /* Do the table lookup on r_A */
    comment("Lookup on rA");
    nflop += 22;
    nflop += table_index("rA*tabscale");    
    assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
    
    if(DO_COULTAB) {
      comment("Coulomb table");
      nflop += extract_table("n1");
      assign("VVCa","qqA*VV");
      assign("FFCa","qqA*tabscale*FF");
      nflop += 4;
    }
    if(DO_VDWTAB) {
      comment("Dispersion");
      nflop += extract_table( DO_COULTAB ? "n1+4" :"n1");
      assign("VVDa","c6a*VV");
      assign("FFDa","c6a*tabscale*FF");
      nflop += 3;
      if(DO_BHAM) {
	comment("Buckingham repulsion");
	/* Make a new lookup index for the exponential table */
	nflop += table_index("cexp2a*rA*tabscale");
	assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
	extract_table( DO_COULTAB ? "n1+8" : "n1+4");
	assign("VVRa","cexp1a*VV");
	assign("FFRa","cexp1a*cexp2a*exptabscale*FF");
	nflop += 6;
      } else {
	comment("Repulsion");
	extract_table( DO_COULTAB ? "n1+8" : "n1+4");	
	assign("VVRa","c12a*VV");
	assign("FFRa","c12a*tabscale*FF");
	nflop += 3;
      }
    }
    do_else(); /* If all A interaction parameters are 0 */
    if(DO_COULTAB) {
      assign("VVCa","0");
      assign("FFCa","0");
    }
    if(DO_VDWTAB) {
      assign("VVDa","0");
      assign("FFDa","0");
      assign("VVRa","0");
      assign("FFRa","0");
    }
    assign("rinv5a","0");    
    end_if(); /* finished A */

    ifbuf[0]=0;
 
    /* now do B */
    if(DO_COULTAB) {
      assign("qqB","iqB*%s",ARRAY(chargeB,jnr));
      sprintf(ifbuf, bC ? "(qqB != 0)" : "(qqB.ne.0)");
    }
    if(DO_VDWTAB) {
      if(strlen(ifbuf)>0)
	strcat(ifbuf, bC ? " || " : ".or.");
      strcat(ifbuf, bC ? "(c6b > 0) || " : "(c6b.gt.0).or.");
      if(DO_BHAM)
	strcat(ifbuf, bC ? "(cexp1b > 0)" : "(cexp1b.gt.0)");
      else
	strcat(ifbuf, bC ? "(c12b > 0)" : "(c12b.gt.0)");
    }
    
    start_if(ifbuf);
    
    assign("rB", bC ? "pow(Alpha*sigma6b*L12+rsix,onesixth)" :
	   "cpow(Alpha*sigma6b*L12+rsix,onesixth)");
    assign("rinvb","1.0/rB");
    assign("rinv5b","rinvb*rinvb");
    assign("rinv5b","rinv5b*rinv5b*rinvb");
    /* do the B table */
    comment("Lookup on rB");
    nflop += 1 + table_index("rB*tabscale");    
    
    assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
    if(DO_COULTAB) {
      comment("Coulomb table");
      nflop += extract_table("n1");
      assign("VVCb","qqB*VV");
      assign("FFCb","qqB*tabscale*FF");
      nflop += 4;
    }
    if(DO_VDWTAB) {
      comment("Dispersion");
      nflop += extract_table( DO_COULTAB ? "n1+4" : "n1");
      assign("VVDb","c6b*VV");
      assign("FFDb","c6b*tabscale*FF");
      nflop += 3;
      if(DO_BHAM) {
	comment("Buckingham repulsion");
	/* Make a new lookup index for the exponential table */
	nflop += 2 + table_index("cexp2b*rB*tabscale");
	assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
	nflop += extract_table( DO_COULTAB ? "n1+8" : "n1+4");
	assign("VVRb","cexp1b*VV");
	assign("FFRb","cexp1b*cexp2b*exptabscale*FF");
	nflop += 4;
      } else {
	comment("Repulsion");	
	nflop += extract_table( DO_COULTAB ? "n1+8" : "n1+4");
	assign("VVRb","c12b*VV");
	assign("FFRb","c12b*tabscale*FF");
	nflop += 3;
      }
    }
    do_else(); /* If all B interaction parameters are 0 */
    if(DO_COULTAB) {
      assign("VVCb","0");
      assign("FFCb","0");
    }
    if(DO_VDWTAB) {
      assign("VVDb","0");
      assign("FFDb","0");
      assign("VVRb","0");
      assign("FFRb","0");
    }
    assign("rinv5b","0");
    end_if(); /* finished B */
    
    /* OK. Now we have all potential and force lookup values,
     * all that is left is to calculate the actual force, potential
     * and dv/dlambda.
     */    
    if(DO_COULTAB) {
      add_to_buffer(fabuf,"FFCa");
      add_to_buffer(fbbuf,"FFCb");
      assign("vcoul","lambda*VVCb+L1*VVCa");
      add_to_buffer(dvdlbuf,"VVCb-VVCa");
      nflop += 6;
    }
    if(DO_VDWTAB) {		
      increment("vnbtot","lambda*(VVDb+VVRb)+L1*(VVDa+VVRa)");
      add_to_buffer(dvdlbuf,"VVDb+VVRb-VVDa-VVRa");
      add_to_buffer(fabuf,"FFDa+FFRa");
      add_to_buffer(fbbuf,"FFDb+FFRb");
      nflop += 13;
    }    
    assign("Fa","-(%s)",fabuf);  
    assign("Fb","-(%s)",fbbuf);
    assign("fs%d","(L1*Fa*rinv5a + lambda*Fb*rinv5b)*rfour",ij);    
    add_to_buffer(dvdlbuf,"onethird*Alpha*lambda*L1*(Fb*sigma6b*rinv5b-Fa*sigma6a*rinv5a)");
    increment("dvdl",dvdlbuf);
    nflop += 16;
    return nflop;
}



int do_table(int i,int j)
{
  int ij = 10*i+j;
  char qq[100],buf[100];
  char f1buf[100],f2buf[100];  
  int nflop = 0;

  f1buf[0]=f2buf[0]=0;
  
  comment("Tabulated interactions");
  sprintf(buf,"r%d*tabscale",ij);
  nflop += table_index(buf) + 1;

  assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
  if(loop.sol==SOL_MNO && !loop.coul && table_element_size==12)
    assign("n1","n1+4");
    
  if(DO_COULTAB) {
    /* determine the name of and set the charge parameter q_i*q_j */
    nflop += assign_chargeproduct(i,j,qq);
    if(loop.free) {
      assign("qqq","L1*qqA + lambda*qqB");
      nflop += 3;
    }
    nflop += extract_table("n1");
    assign("vcoul","%s*VV",qq);
    nflop ++;
    if(DO_FORCE) {
      assign("fijC","%s*FF",qq);
      nflop ++;
    }
    if(loop.free) {
      increment("dvdl","(qqB - qqA)*VV",j,j);
      nflop += 7;
    }
    if(DO_FORCE) {
      add_to_buffer(f1buf,"fijC");
      nflop ++;
    }
  }
  /* only do VDW for the first atom in water interactions */
  if(DO_VDWTAB && 
     (!DO_WATER || (loop.sol==SOL_WATER && i==1) || (i==1 && j==1))) {
    /* first determine nonbonded parameters */
    if(loop.free) {
      assign("tjA","ntiA+%d*%s%s", N_VDWPARAM, ARRAY(type,jnr) , bC ? "" : "+1");
      assign("tjB","ntiB+%d*%s%s", N_VDWPARAM, ARRAY(typeB,jnr) , bC ? "" : "+1");
      assign("c6a",ARRAY(nbfp,tjA));
      assign("c6b",ARRAY(nbfp,tjB));
      assign("c6","L1*c6a + lambda*c6b");
      nflop += 3;
      assign( DO_BHAM ? "cexp1a" : "c12a" ,ARRAY(nbfp,tjA+1));
      assign( DO_BHAM ? "cexp1b" : "c12b" ,ARRAY(nbfp,tjB+1));
      if(DO_BHAM) {
	assign("cexp2a",ARRAY(nbfp,tjA+2));
	assign("cexp2b",ARRAY(nbfp,tjB+2));
      } else {
	assign("c12","L1*c12a + lambda*c12b");
	nflop += 3;
      }
    } else if(loop.sol!=SOL_WATERWATER) {
      assign("tjA","ntiA+%d*%s%s", N_VDWPARAM, ARRAY(type,jnr) , bC ? "" : "+1");
      assign("c6",ARRAY(nbfp,tjA));
      assign( DO_BHAM ? "cexp1" : "c12",ARRAY(nbfp,tjA+1));
      if(DO_BHAM) 
	assign("cexp2",ARRAY(nbfp,tjA+2));
    }    
    comment("Dispersion");
    nflop += extract_table( DO_COULTAB ? "n1+4" : "n1");
    assign("vnb6","c6*VV");
    nflop++;
    if(DO_FORCE) {
      assign("fijD","c6*FF");
      nflop++;
    }
    if(loop.free) {
      increment("dvdl","(c6b - c6a)*VV");
      nflop += 3;
    }
    if(DO_FORCE) {
     add_to_buffer(f1buf,"fijD");
     nflop++;
    }
    if(DO_BHAM) {
      comment("Buckingham repulsion");
      /* Make a new lookup index for the exponential table */
      sprintf(buf, loop.free ? "cexp2a*r%d*exptabscale" :
	      "cexp2*r%d*exptabscale",ij);
      nflop += 2 + table_index(buf);
      assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
      nflop += extract_table( DO_COULTAB ? "n1+8" : "n1+4");

      if(loop.free) {
	assign("vnbexpa","cexp1a*VV");
        if(DO_FORCE)
	  assign("fijRa","cexp1a*cexp2a*FF");	
	sprintf(buf,"cexp2b*r%d*exptabscale",ij);
	nflop += table_index(buf);
	assign("n1","%d*n0%s",table_element_size, bC ? "" : "+1");
	nflop += extract_table( DO_COULTAB ? "n1+8" : "n1+4");
	
	assign("vnbexpb","cexp1b*VV");
	if(DO_FORCE) {
	  assign("fijRb","cexp1b*cexp2b*FF");
	  assign("fijR","L1*fijRa + lambda*fijRb");
	}
	increment("vnbtot","vnb6 + L1*vnbexpa + lambda*vnbexpb");
	increment("dvdl","vnbexpb - vnbexpa");
	nflop += 18;
      } else {
	assign("vnbexp","cexp1*VV");
	nflop++;
	if(DO_FORCE) {
	  assign("fijR","cexp1*cexp2*FF");
	  nflop += 2;
	}
	increment("vnbtot","vnb6 + vnbexp");
      }
      if(DO_FORCE) {
	add_to_buffer(f2buf,"fijR");
	nflop++;
      }
    } else {
      comment("Lennard-Jones repulsion");   
      nflop += extract_table( DO_COULTAB ? "n1+8" : "n1+4");
      assign("vnb12","c12*VV");
      if(DO_FORCE) {
	assign("fijR","c12*FF");
	nflop++;
      }
      increment("vnbtot","vnb6 + vnb12");
      if(loop.free) {
	increment("dvdl","(c12b - c12a)*VV");
	nflop += 3;
      }
      add_to_buffer(f1buf,"fijR");
      nflop += 4;
    }
  }

  if(DO_FORCE) {
    sprintf(buf,"(%s)*tabscale",f1buf);
    if(strlen(f1buf)) {
      add_to_buffer(tabforcebuf,buf);
    }   
    sprintf(buf,"(%s)*exptabscale",f2buf);  
    if(strlen(f2buf)) {
      add_to_buffer(tabforcebuf,buf);
    }
    nflop+=2;
  }
  return nflop;
}


int do_coul(int i, int j)
{
    int nflop = 0;
    char qq[16];
    int ij = 10*i+j;
    
    comment("Coulomb");
  
    /* determine the name of and set the charge product iq*jq,
     * depending on wheter we do water loops or not. */
    nflop += assign_chargeproduct(i,j,qq);
    
    /* this short code is the actual interaction! */
    assign("vcoul","%s*rinv%d",qq,ij);      
    /* concatenate the potential to the force buffer string.
     * This is later multiplied by 1/r to get the force once
     * all interactions for this pair have been calculated
     */
    nflop += 1;
    if(DO_FORCE) {
      nflop ++;
      add_to_buffer(forcebuf,"vcoul");
    }
    return nflop;
}


int do_rf(int i, int j)
{
    int nflop = 0;
    char qq[16],buf[50];
    int ij = 10*i+j;
    
    comment("Reaction field coulomb");
  
    nflop += assign_chargeproduct(i,j,qq);

    /* this short code is the actual interaction! */
    assign("krsq","krf*rsq%d",ij);
    assign("vcoul","%s*(rinv%d+krsq-crf)",qq,ij);

    /* concatenate the potential to the force buffer string.
     * This is later multiplied by 1/r to get the force once
     * all interactions for this pair have been calculated
     */
    nflop += 4;
    if(DO_FORCE) {
      sprintf(buf,"%s*(rinv%d-two*krsq)",qq,ij);
      add_to_buffer(forcebuf,buf);
      nflop += 4;
    }
    return nflop;
}


int do_lj(int i, int j)
{
    int ij = 10*i+j;
    int nflop=0;

    comment("LJ");
    /* calculate r^-6 from r^-2 */
    assign("rinvsix","rinvsq%d*rinvsq%d*rinvsq%d",ij,ij,ij);
    /* for water-water loops this is only called for O-O interactions,
     * and we have already extracted those c6/c12
     */
    if(loop.sol==SOL_WATERWATER) {
      assign("vnb6","c6*rinvsix");
      assign("vnb12","c12*rinvsix*rinvsix");
    } else {
      assign("tjA","ntiA+2*%s%s", ARRAY(type,jnr) , bC ? "" : "+1");
      assign("vnb6","rinvsix*%s",ARRAY(nbfp,tjA));
      assign("vnb12","rinvsix*rinvsix*%s",ARRAY(nbfp,tjA+1));
    }

    nflop = 7;
    if(DO_FORCE) {
      add_to_buffer(forcebuf,"twelve*vnb12-six*vnb6");
      nflop += 4;
    }
    increment("vnbtot","vnb12-vnb6");
    return nflop;
}


int do_bham(int i, int j)
{
  int ij = 10*i+j;
  int nflop=0;

  comment("Buckingham");
  
  /* calculate r^-6 from r^-2 */
  assign("rinvsix","rinvsq%d*rinvsq%d*rinvsq%d",ij,ij,ij);
  if(loop.sol==SOL_WATERWATER) {
    assign("vnb6","c6*rinvsix");
    assign("br","cexp2*r%d",ij);
    assign("vnbexp","cexp1*exp(-br)"); /* V=C_a*exp(-C_b*r) */
  } else {
    /* fortran indices start at one, instead of 0 for C.
     * And there are three nb parameters for buckingham
     */
    assign("tjA","ntiA+3*%s%s", ARRAY(type,jnr) , bC ? "" : "+1");
    assign("vnb6","rinvsix*%s",ARRAY(nbfp,tjA));
    assign("br","r%d*%s",ij,ARRAY(nbfp,tjA+2));
    assign("vnbexp","exp(-br)*%s",ARRAY(nbfp,tjA+1));
  }
  nflop=8;
  if(DO_FORCE) {
    add_to_buffer(forcebuf,"br*vnbexp-six*vnb6");    
    nflop += 4;
  }
  increment("vnbtot","vnbexp-vnb6");
  return nflop;
}


int calc_interactions(void)
{
  /* Remember to set need_invsqrt, need_reciprocal and need_exp 
   * to the correct needs for your new interaction already 
   * in mkinl.c. Then add new entries in the switch statements
   * and call your interaction implementation in a separate routine.
   * Have a look at e.g. do_coul to see how it is done!
   */
  int nflop = 0;
  int i,j,ij;
  bool do_vdw;
  /* once we get here the distance measures 
   * we requested in mkinl.c are calculated, but we might
   * have to extract them from the buffers
   * used for vectorization.
   */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      ij=10*i+j;
      /* Dont do vdw for coul-only water atoms. For the general
       * solvent we change the value of loop.vdw outside this loop */
      do_vdw = loop.vdw && 
	(!DO_WATER || (loop.sol==SOL_WATER && i==1) ||
	 (loop.sol==SOL_WATERWATER && i==1 && j==1));
      
      forcebuf[0]=tabforcebuf[0]=0;
      
      if(DO_VECTORIZE) {
	/* get rinv or rinvsq from temp buffer */
	assign( loop.invsqrt ? "rinv%d" : "rinvsq%d",
		OVERWRITE_RSQ ? ARRAY(buf1,m) : ARRAY(buf2,m), ij);
	/* Do we need to copy rsq too? */
	if((loop.coul_needs_rsq) || (loop.vdw_needs_rsq && do_vdw))
	  assign( "rsq%d" , ARRAY(buf1,m), ij);
	/* And what about r? */
	if((loop.coul_needs_r) || (loop.vdw_needs_r && do_vdw))
	  assign( "r%d" ,"%s*%s",ij,ARRAY(buf1,m),ARRAY(buf2,m));
	increment("m","1");
      } else { /* no vectorization */
	/* If we aren't trying to hide the latency of invsqrt 
	 * we calculate it here just before doing the
	 * interaction for each pair of atoms!
	 */
	if(opt.delay_invsqrt) {
	  assign( loop.invsqrt ? "rinv%d" : "rinvsq%d",
		  loop.invsqrt ? "1.0/sqrt(rsq%d)" :
		  "1.0/rsq%d",ij,ij);
	  nflop++;
	}
	if((loop.coul_needs_r) || (loop.vdw_needs_r && do_vdw)) {
	  assign( "r%d","rsq%d*rinv%d", ij,ij,ij);
	  nflop++;
	}
      }
      /* Do we need rinvsq when calculating rinv? */
      if(((loop.coul_needs_rinvsq) || 
	  (loop.vdw_needs_rinvsq && do_vdw)) && !loop.recip) {
	assign("rinvsq%d","rinv%d*rinv%d",ij,ij,ij);
	nflop++;
      }
      /* which interactions should we call ? */
      if(do_vdw) 
	switch(loop.vdw) {
	case VDW_LJ:
	  nflop += do_lj(i,j);
	  break;
	case VDW_BHAM:
	  nflop += do_bham(i,j);
	  break;
	case VDW_TAB:
	case VDW_BHAMTAB:
	  nflop += DO_SOFTCORE ? do_softcore(i,j) : do_table(i,j);
	  break;
	default:
	  comment("No nonbonded interaction");
	  break;
	}
      
      switch(loop.coul) {
      case COUL_NORMAL:
	nflop += do_coul(i,j);
	break;
      case COUL_RF:
	nflop += do_rf(i,j);
	break;
      case COUL_TAB:
	if(!(do_vdw && DO_VDWTAB))
	  nflop += DO_SOFTCORE ? do_softcore(i,j) : do_table(i,j);
	/* coul table was done simultaneous with
	 * vdw if the latter was tabulated.
	 */
	break;
      default:
	comment("No coulomb interaction");	  
	break;
      }
      /* The total coulomb energy is incremented after calculating
       * the scalar force to hide a flop.
       */

      if(DO_FORCE)
	nflop += calc_scalar_force(i,j);

      if(loop.coul) {
	increment("vctot","vcoul");
	nflop ++;
      }
      if(DO_FORCE)
	nflop += update_inner_forces(i,j);
    }       
  return nflop;
}



