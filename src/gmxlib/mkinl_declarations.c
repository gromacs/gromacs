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
#include <mkinl_fortrandata.h>

void init_block_data(void) 
{
  /* block to initialize the table data - write directly to file */
  fprintf(output,"%sblock data\n",indent());
  fprintf(output,"%sinteger*4 I\n",indent());

  if(arch.gmx_invsqrt) 
    fprintf(output,"%scommon /finvsqrtdata/ finvsqrtexptab(256),finvsqrtfracttab(4096)\n",indent());
  
  if(arch.gmx_invsqrt)
    fprintf(output,"%sinteger*4 finvsqrtexptab,finvsqrtfracttab\n",indent());

  if(arch.gmx_invsqrt)
    fprintf(output,finvsqrtdata);   
  
  fprintf(output,"%send\n\n\n",indent());
}


void file_header(void)
{ 
  edit_warning("mkinl");
  
  if (bC) {
    if(arch.threads)
      fprintf(output,"#include <pthread.h>\n"); /* must come first */
    fprintf(output,"#include <stdio.h>\n"
	    "#include <math.h>\n"
	    "#include \"typedefs.h\"\n"
	    "#include \"inner.h\"\n"
#ifdef USE_AXP_ASM 
	    "#include \"axp_asm.h\"\n"
#endif
	    "#include <vec.h>\n");
  } else {
    /* Write some full fortran functions, performing the
     * software invsqrt code.
     * Corresponding c table routines are found in vec.h
     */

    fortran_invsqrt();
    if(arch.gmx_invsqrt)
      init_block_data();
  } 

#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined DOUBLE && defined USE_X86TRUNC)
  /* This is a very fast inline assembly truncation routine
   * The slow step an most x86 chips is to store and restore
   * the control word when we want to do truncation. Here it
   * is instead stored first in each innerloop, then we
   * dont care about it until the end of the loop when it is restored
   * But beware, dont try to do any roundings in the loop!
   */
  if(bC)
    fprintf(output,"\n\n#define x86trunc(a,b) asm(\"fld %%1 ; fistpl %%0\" : \"=m\" (*&b) : \"f\" (a));\n\n");
#endif
#ifdef HAVE_LIBMASSV_ANY
  if(bC)
    fprintf(output,
	    "void vsrec(float *, float *, int *);\n"
	    "void vsrsqrt(float *, float *, int *);\n"
	    "void vrec(double *, double *, int *);\n"
	    "void vrsqrt(double *, double *, int *);\n");
#endif
  
}


void function_info(void)
{
  char s0[60],s1[60],s2[60],s3[60],s4[60],buf[512];

  /* Add a nice header before each routine telling what it does */
  

  switch(loop.do_force) {
  case TRUE:
    sprintf(s0,"Calculated");
    break;
  case FALSE:
    sprintf(s0,"Not calculated");
    break;
  }

  switch(loop.coul) {
  case COUL_NO:
    sprintf(s1,"Not calculated");
    break;
  case COUL_NORMAL:
    sprintf(s1,"Normal");
    break;
  case COUL_RF:
    sprintf(s1,"Reaction field");
    break;
  case COUL_TAB:
    sprintf(s1,"Tabulated");
    break;
  default:
    printf("Error: Bad value of loop.coul encountered: %d\n",loop.coul);
    exit(-1);
    break;
  }
  
  switch(loop.vdw) {
  case VDW_NO:
    sprintf(s2,"Not calculated");
    break;
  case VDW_LJ:
    sprintf(s2,"Lennard-Jones");
    break;
  case VDW_BHAM:
    sprintf(s2,"Buckingham");
    break;
  case VDW_TAB:
    sprintf(s2,"Tabulated");
    break;
  case VDW_BHAMTAB:
    sprintf(s2,"Tabulated Buckingham");
    break;
  default:
    printf("Error: Bad value of loop.vdw encountered\n");
    exit(-1);
    break;
  }

  switch(loop.sol) {
  case SOL_NO:
    sprintf(s3,"No");
    break;
  case SOL_MNO:
    sprintf(s3,"general M:N:O solvent - other atom");
    break;
  case SOL_WATER:
    sprintf(s3,"Water (%d atoms) - other atom",WATERATOMS);
    break;
  case SOL_WATERWATER:
    sprintf(s3,"Water (%d atoms) - water (%d atoms)",WATERATOMS,WATERATOMS);
    break;
  default:
    printf("Error: Bad value of loop.sol encountered\n");
    exit(-1);
    break;
  }

    switch(loop.free) {
  case FREE_NO:
    sprintf(s4,"No");
    break;
  case FREE_LAMBDA:
    sprintf(s4,"Lambda (alpha=0)");
    break;
  case FREE_SOFTCORE:
    sprintf(s4,"Softcore");
    break;
  default:
    printf("Error: Bad value of loop.free encountered\n");
    exit(-1);
    break;
  }

  if(bC) 
    sprintf(buf,
	    "  /**********************************************************\n" 
	    "   * This is gromacs innerloop %s\n"
            "   * Forces:      %s\n"
	    "   * Coulomb:     %s\n"
	    "   * Nonbonded:   %s\n"
	    "   * Solvent opt: %s\n"
	    "   * Free energy: %s\n"
	    "   **********************************************************/"
	    "\n",loopname,s0,s1,s2,s3,s4);
  else
    sprintf(buf,
	    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
	    "C    This is gromacs innerloop %1s\n" 
	    "C    Forces:      %s\n"
	    "C    Coulomb:     %s\n"
	    "C    Nonbonded:   %s\n"
	    "C    Solvent opt: %s\n"
	    "C    Free energy: %s\n"
	    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	    "\n",loopname,s0,s1,s2,s3,s4);
  strcat(header,buf);
}

void func_header(char *appendname)
{
  char buf[512];

  sprintf(buf,"\n\n\n\n");
  strcat(header,buf);
  sprintf(loopname,"%sinl%1x%1x%1x%1x%s",
          DO_FORCE ? "" : "mc",
	  loop.coul,loop.vdw,loop.sol,loop.free,appendname);
  
  function_info();
  if(bC)
    sprintf(buf,"void %s(",loopname);
  else
    sprintf(buf,"%ssubroutine %s(",indent(),loopname);
  strcat(header,buf);
}


void func_args()
{
  /* We declare EVERYTHING here, even redundant stuff.
   * (Only exception is we cannot declare a variable twice)
   * Unreferenced things will be removed automatically 
   * by the buffer flushing routine!
   * The order is still significant, though. The remaining
   * stuff will be declared in that order in the argument
   * list, influencing the calling format.
   */
  
  /* 1: COMMON_ARGS */
  /* number of i particles, their indices and neighbourlists
   * with j particles
   */
  declare_int("nri");
  declare_const_int_vector("iinr");
  declare_const_int_vector("jindex");
  declare_const_int_vector("jjnr");
  /* shift indices, vectors, force shifts and group indices */
  declare_const_int_vector("shift");
  declare_const_real_vector("shiftvec");
  declare_real_vector("fshift");
  declare_const_int_vector("gid");


  /* coordinates, forces, 
   * charges and VDW types (also for B state with free energy)
   */
  declare_const_real_vector("pos");
  declare_real_vector("faction");
  /* buffer to remove dependencies on vector machines */
  declare_real_vector("fbuf");
  /* Thread parallelization stuff */
  declare_int("mythread");
  declare_int("nthreads");
  /* Count is a global variable, but to make it easier
   * for the compiler to inline the mutex routines we pass
   * it as a pointer. 
   */
  declare_int( bC ? "*count" : "count");
  /* mask the mutex pointer as an integer passed by
   * reference when using fortran. 
   */
  declare_other( bC ? "pthread_mutex_t" : "integer",
		 bC ? "*mtx" : "mtx");

  /* Temporary work arrays used for vectorization
   */
  declare_real_vector("drbuf");
  declare_real_vector("buf1");
  declare_real_vector("buf2");

  /* 2: COUL_ARGS */
  declare_const_real_vector("charge");
  declare_real("facel");
  declare_real_vector("Vc");
  /* RF_ARGS - Reaction field constant */
  declare_real("krf");
  declare_real("crf");

  /* 3: LJ_ARGS */
  declare_const_int_vector("type");
  declare_int("ntype");
  declare_const_real_vector("nbfp");
  declare_real_vector("Vnb");

  /* 4: TABLE ARGUMENTS */
  declare_real("tabscale");
  declare_const_real_vector("VFtab");
  declare_real("exptabscale");
  
  /* 5: FREE ENERGY */
  declare_real("lambda");
  declare_real( bC ? "*dvdlambda" : "dvdlambda");
  /* Free energy charge stuff */
  declare_const_real_vector("chargeB");
  /* Free energy lj stuff */
  declare_const_int_vector("typeB");
  /* softcore thingies */
  declare_real("Alpha");
  declare_real("defsigma6");

  /* Vector with molecule info for General M:N solvent stuff */
  declare_const_int_vector("nsatoms");
  
  nargs=ndecl;
  /* set the number of arguments to the current number of declarations.
   * Further declarations will be local variables.
   */
}

    
void func_localvars()
{ 
  int i,j;
  char buf[100];
  /* Declare EVERYTHING here, even redundant stuff.
   * (But you may only declare a variable once).
   * Unreferenced things will be removed automatically 
   * by the buffer flushing routine.
   * The order is still significant, though. The remaining
   * stuff will be declared in that order in the variable list.
   */

  /* Fortran functions */
  if(!bC) {
    declare_real("invsqrt");
    declare_real("rec");
  }
  
  /* variables */
  declare_int("ii");  /* i particle nr */
  declare_int("k");   
  declare_int("kk");  /* vector machine stuff */ 
  declare_int("n");   /* index in i particle list */
  declare_int("nn0"); /* limits for i particle loop */
  declare_int("nn1");
  declare_int("nn2");
  declare_int("nj0"); /* limits for inner loop */
  declare_int("nj1");
  declare_int("is3"); /* shift vector index */
  declare_int("ggid");
  declare_int("ii3"); 
  declare_int("m");
  declare_int("m3");
  declare_int("s");
  declare_int("n0");
  declare_int("n1");
  declare_int("nnn");
  declare_real("krsq");
  declare_real("qqA");
  declare_real("qqB");
  declare_real("qqq");
  declare_real("qO");
  declare_real("qH");
  declare_real("qqOO");
  declare_real("qqOH");
  declare_real("qqHH");
  declare_real("qq");
  declare_real("iq");
  declare_real("jq");
  declare_real("rinvsix");
  declare_real("cexp1");
  declare_real("cexp2");
  declare_real("cexp1a");
  declare_real("cexp2a");
  declare_real("cexp1b");
  declare_real("cexp2b");
  declare_real("c6");
  declare_real("c12");
  declare_real("c6a");
  declare_real("c12a");
  declare_real("c6b");
  declare_real("c12b");
  declare_real("vctot");
  declare_real("vcoul");
  declare_real("rt");
  declare_real("eps");
  declare_real("eps2");
  declare_real("Y");
  declare_real("F");
  declare_real("Geps");
  declare_real("Heps2");
  declare_real("Fp");
  declare_real("VV");
  declare_real("FF");
  declare_real("fijC");
  declare_real("fijD");
  declare_real("fijR");
  declare_real("dvdl");
  declare_real("L1");
  declare_real("L12");
  declare_real("lam2");
  declare_real("sigma6a");
  declare_real("sigma6b");
  declare_real("rA");
  declare_real("rB");
  declare_real("rfour");
  declare_real("rsix");
  declare_real("rinva");
  declare_real("rinvb");
  declare_real("rinv5a");
  declare_real("rinv5b");
  declare_real("Fa");
  declare_real("Fb");
  declare_real("vnbexpa");
  declare_real("vnbexpb");
  declare_real("fijRa");
  declare_real("fijRb");
  declare_real("vnb6");
  declare_real("vnb12");
  declare_real("vnbexp");
  declare_real("vnbtot");
  declare_real("br");
  declare_real("VVCa");
  declare_real("VVCb");
  declare_real("FFCa");
  declare_real("FFCb");
  declare_real("VVDa");
  declare_real("VVDb");
  declare_real("FFDa");
  declare_real("FFDb");
  declare_real("VVRa");
  declare_real("VVRb");
  declare_real("FFRa");
  declare_real("FFRb");
  declare_int("x86_cw");
  declare_int("x86_cwsave");
  
  /* neighbourlist indices */
 
  declare_int("j3");
  declare_int("jnr");
  
  /* shift vectors */
  declare_real("shX");
  declare_real("shY");
  declare_real("shZ");

  /* constants */
  declare_const_real("nul",0.0);
  declare_const_real("one",1.0);
  declare_const_real("two",2.0);
  declare_const_real("three",3.0);
  declare_const_real("six",6.0);
  declare_const_real("twelve",12.0);
  declare_const_real("onesixth",1.0/6.0);
  declare_const_real("onethird",1.0/3.0);
  
  /* i coordinates */
  for(i=1;i<=loop.ni;i++) {
    sprintf(buf,"ix%d",i);
    declare_real(buf);
    sprintf(buf,"iy%d",i);
    declare_real(buf);
    sprintf(buf,"iz%d",i);
    declare_real(buf);
  }
  /* i charges */
  declare_real("iqA");
  declare_real("iqB");

  declare_int("ntiA");
  declare_int("ntiB");
  declare_int("tjA");
  declare_int("tjB");

  /* i forces */
  for(i=1;i<=loop.ni;i++) {
    sprintf(buf,"fix%d",i);
    declare_real(buf);
    sprintf(buf,"fiy%d",i);
    declare_real(buf);
    sprintf(buf,"fiz%d",i);
    declare_real(buf);
  }

  /* j coordinates (also prefetched) */
  for(j=1;j<=loop.nj;j++) {
    sprintf(buf,"jx%d",j);
    declare_real(buf);
    sprintf(buf,"jy%d",j);
    declare_real(buf);
    sprintf(buf,"jz%d",j);
    declare_real(buf);
  }
  /* j charges and types */
    declare_real("jqA");
    declare_real("jqB");
  
  /* i-j vectorial distance */  
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"dx%d",10*i+j);
      declare_real(buf);
      sprintf(buf,"dy%d",10*i+j);
      declare_real(buf);
      sprintf(buf,"dz%d",10*i+j);
      declare_real(buf);
    }

  /* square distance */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"rsq%d",10*i+j);
      declare_real(buf);
    }
  /* inverse dist */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"rinv%d",10*i+j);
      declare_real(buf);
    }
  /* inverse square dist */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"rinvsq%d",10*i+j);
      declare_real(buf);
    }
  /* distance */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"r%d",10*i+j);
      declare_real(buf);
    }

  /* scalar force */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"fs%d",10*i+j);
      declare_real(buf);
    }

  /* i-j temp. force */
  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      sprintf(buf,"tx%d",10*i+j);
      declare_real(buf);
      sprintf(buf,"ty%d",10*i+j);
      declare_real(buf);
      sprintf(buf,"tz%d",10*i+j);
      declare_real(buf);
    }
  /* j forces */
  for(j=1;j<=loop.nj;j++) {
    sprintf(buf,"fjx%d",j);
    declare_real(buf);
    sprintf(buf,"fjy%d",j);
    declare_real(buf);
    sprintf(buf,"fjz%d",j);
    declare_real(buf);
  }  

  /* M:N solvent stuff */
  declare_int("nstot");
  declare_int("nsvdwc");
  declare_int("nscoul");

  /* Wrapper routine to call C pow() from f77 */
  declare_real("cpow");
  
  if(DO_INLINE_INVSQRT)
    invsqrt_vars();
}


void init_water_data(void)
{ 
  comment("initiation of water data");
  assign("ii","%s%s", bC ? ARRAY(iinr,0) : ARRAY(iinr,1) , bC ? "" : "+1");
  
  /* assign the charge combinations for OO,OH and HH */
  /* we're always doing coulomb */
  
  if(!arch.simplewater  || loop.sol==SOL_WATERWATER) {
    assign("qO", ARRAY(charge,ii));
    assign("qH", ARRAY(charge,ii+1));   
  }
  
  if(!arch.simplewater && loop.sol==SOL_WATERWATER) {
    assign("qqOO","facel*qO*qO");
    assign("qqOH","facel*qO*qH");
    assign("qqHH","facel*qH*qH");
  } else if(!arch.simplewater) {
    assign("qO","qO*facel");
    assign("qH","qH*facel");
  }
  if(arch.simplewater)
    comment("For simplewater, qO and qH are j charges (thus no facel)");
  
  /* assign the nonbonded combination for the Oxygen-oxygen interactions */
  if(loop.vdw) {
   
    assign("ntiA","%d*ntype*%s",N_VDWPARAM, ARRAY(type,ii));
    
    if(loop.sol==SOL_WATERWATER) {
      assign("tjA","ntiA+%d*%s%s",N_VDWPARAM,  ARRAY(type,ii), bC ? "" : "+1");
      assign("c6",ARRAY(nbfp,tjA));
      if(DO_BHAM) {
	assign("cexp1",ARRAY(nbfp,tjA+1));
	assign("cexp2",ARRAY(nbfp,tjA+2));
      } else
	assign("c12",ARRAY(nbfp,tjA+1));
    }
  }
}


void func_init_vars()
{
  newline();
  
  if(loop.free) {
    assign("dvdl","nul");
    assign("L1","one - lambda");
    if(loop.free==FREE_SOFTCORE) {
      assign("lam2","lambda*lambda");
      assign("L12","L1*L1");
    }
  } 
  if(DO_WATER)
    init_water_data();
    
  /* It should be faster to divide once by the number of
   * elements in each loop iteration than to multiply once
   * for each nlist when using vectorization buffers.
   * (But for MNO solvent we have to do it)
   */

#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined DOUBLE && defined USE_X86TRUNC)
  /* store cw */
  if(bC && DO_TAB) {
    strcat(codebuffer,"asm(\"fnstcw %%0\" : \"=m\" (*&x86_cwsave));\n");
    strcat(codebuffer,"x86_cw = x86_cwsave | 3072;\n");
    strcat(codebuffer,"asm(\"fldcw %%0\" : : \"m\" (*&x86_cw));\n");
  }
#endif
  newline(); 
}   


