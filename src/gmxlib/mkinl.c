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
static char *SRCID_mkinl_c = "$Id$";
#include <string.h>
#include "mkinl.h"
#include <types/simple.h>

/* This program generates inner-loop source code for the GROMACS
 * Molecular Dynamics package. The innerloops are nearly optimal, i.e.
 * there are no superfluous declarations or statements in the code.
 * The loops are generated in either C or Fortran 77, according to
 * the first argument when calling the routine ("c" or "fortran")
 *
 * The following defines, usually entered in Makefile.$GMXCPU, affect
 * the code:
 *
 * ARCHITECTURAL AND COMPILER OPTIONS:
 *
 * -DUSE_VECTOR        Generates code more suitable for vector
 *                     computers.
 *
 * -DCRAY_PRAGMA       Insert a Cray compatible pragma that instructs
 *                     the compiler to ignore dependencies in loops
 *                     Contact us if you use it, we are interested in
 *                     the effect of these options on real vector machines!
 *
 * -DSIMPLEWATER       Expand the water loop to have 3 innerloops. This is
 *                     a workaround for SGI compilers.
 *
 * -DUSE_IBM_MASS      Use the vectorized invsqrt and reciprocal routines
 *                     from IBM's mathematical accelerated subsystem
 *                     library.
 *
 * -DUSE_AXP_ASM       Use gromacs assembly routines to perform
 *                     vectorized invsqrt on dec/compaq alpha chips.
 *
 * automatic on gcc:   Change the x86 control word once in each
 *                     routine to speed up truncation. This will
 *                     increase table routine performance significantly.
 *
 * -DTHREADS_MUTEX     Make code suitable for shared-memory multithreading
 * -DTHREADS_STRIDE    parallelization and determines the way to 
 * (NOT SUPPORTED YET) synchronize threads. Striding loops over N processors 
 *                     means they treat every N:th i particle. Using mutexes 
 *                     is more flexible, but might be slower (or faster ;-).
 *                     Mutexes also enables unequal load balancing between
 *                     the processors due to e.g. other jobs on a machine. 
 *                     Only choose one of them!
 *
 * -DSOFTWARE_INVSQRT  Use the gromacs implementation of the inverse
 *                     square root.
 *                     This is probably a good thing to do on most 
 *                     computers except SGI and ibm power2/power3 which
 *                     have special hardware invsqrt instructions.
 *                     (but the PPC604 available in some ibm smp nodes doesnt)
 *
 * INNERLOOP OPTIMIZATION OPTIONS:  
 * 
 * -DDONT_INLINE_INVSQRT
 *                     Turns off inlining of the gromacs inverse square root 
 *                     on architectures which use
 *                     it, which comes handy for debugging. It is
 *                     probably slightly slower, though.
 *
 * -DPREFETCH        Prefetching of forces in normal, solvent, water 
 * -DPREFETCH_S      and water-water loops. 
 * -DPREFETCH_W
 * -DPREFETCH_WW     
 *
 * -DVECTORIZE_INVSQRT    Vectorize the corresponding part of the force
 * -DVECTORIZE_INVSQRT_S  calculation. This can avoid cache trashing, and
 * -DVECTORIZE_INVSQRT_W  in some cases it makes it possible to use vector
 * -DVECTORIZE_INVSQRT_WW intrinsic functions. (also on scalar machines)
 * -DVECTORIZE_RECIP      Water loops are only used when we include coulomb
 *                        interactions (requiring sqrt), so the reciprocal  
 *                        vectorization is only meaningful on normal loops.
 *
 * 
 * -DDECREASE_LOOKUP_LATENCY  Hide 1/sqrt lookup table latency by doing them
 *                            in parallel for water & water-water loops.
 *
 * The special water-water loops can be turned off in ns.c by defining
 *
 * -DDISABLE_WATERWATER_LOOPS
 *
 * Note that all these options should be given in the main makefile, not here,
 * since they affect code in mkinl.h and other places.
 *
 * This code-generating source is necessarily somewhat complex and
 * large since it in essence is a simple pre-compiler,
 * but it is not very complicated, just a lot of if statements.
 * If you want to hack into it, the first thing to do is to turn comments
 * on by setting -DKEEP_COMMENTS, and remove all fancy optimizations. 
 *
 * For the vanilla loops without special optimizations like vectorization
 * we also generate loop counts for the inner and outer loops. 
 *
 * The program flow goes something like this:
 * 
 * 1. main() checks the defines and fills global variable structures 
 *    with optimizations and architectural settings.
 *    The routine loops over the different innerloops to be created.
 *    For each loop we set the loop-specific options, create the 
 *    header and all local variables/initializations.
 *    
 * 2. The outer loop (over neighbourlists) is started by calling
 *    outer_loop(). 
 * 
 * 3. The outer loop routine calls the inner_loop routine which 
 *    constructs the innermost code.
 * 
 * 4. That's it - the loop is done.
 *
 * Erik Lindahl, David van der Spoel 1999-2000
 */

char loopname[STRINGSIZE];

arch_t arch;        /* architectural options */
opt_t  opt;         /* optimization options */
loop_t loop;        /* options for the loop currently being built */

int table_element_size;

void set_loop_options(void) 
{
  /* Determine what we need to calculate, based on the interaction
   * and optimzation features. When implementing a new interaction
   * you should ideally only have to add it to the enumerated lists in
   * mkinl.h, change the lines below so you get the powers of the
   * distance you want, and finally implement the interaction in
   * mkinl_interactions.c
   */
  loop.coul_needs_rinv=TRUE;
  /* all current coulomb stuff needs 1/r */
  loop.coul_needs_rinvsq=(loop.coul && !DO_COULTAB) && DO_FORCE;
  /* non-table coul also need r^-2 , except the MC loops */
  loop.coul_needs_rsq=(DO_RF || (loop.coul && loop.free));
  /* reaction field and softcore need r^2 */
  loop.coul_needs_r=DO_COULTAB;
  /* tabulated coulomb needs r */ 

  /* table nb needs 1/r */
  loop.vdw_needs_rinv=(loop.vdw && DO_VDWTAB);
  /* all other need r^-2 */
  loop.vdw_needs_rinvsq=(loop.vdw && !DO_VDWTAB);
  loop.vdw_needs_rsq=(loop.vdw && loop.free);
  /* softcore needs r^2 */
  loop.vdw_needs_r=(DO_BHAM || DO_VDWTAB || (loop.vdw && loop.free));
  /* softcore and buckingham need r */

  /* Now, what shall we calculate? */ 
  loop.invsqrt=
    (loop.coul && (loop.coul_needs_rinv || loop.coul_needs_r)) ||
    (loop.vdw && (loop.vdw_needs_rinv || loop.vdw_needs_r));

  loop.vectorize_invsqrt=loop.invsqrt && !DO_SOFTCORE &&
    ((loop.sol==SOL_NO && (opt.vectorize_invsqrt & TYPE_NORMAL)) ||
     (loop.sol==SOL_MNO && (opt.vectorize_invsqrt & TYPE_SOLVENT)) ||
     (loop.sol==SOL_WATER && (opt.vectorize_invsqrt & TYPE_WATER)) ||
     (loop.sol==SOL_WATERWATER && (opt.vectorize_invsqrt & TYPE_WATERWATER)));

  loop.recip=!loop.invsqrt;

  loop.vectorize_recip=loop.recip && opt.vectorize_recip && !DO_SOFTCORE;
     
  /* If we are vectorizing the invsqrt it is probably faster to
   * include the vdw-only atoms in that loop and square the
   * result to get rinvsq. Just tell the generator we want rinv for
   * those atoms too:
   */
  loop.vdw_needs_rinv=DO_VDWTAB || (loop.coul && loop.vectorize_invsqrt);
  
  table_element_size=0;
  if(DO_COULTAB)
    table_element_size+=4;
  if(DO_VDWTAB)
    table_element_size+=8;

  /* set the number of i and j atoms to treat in each iteration */
  loop.ni = (DO_WATER && !arch.simplewater) ? WATERATOMS : 1;
  loop.nj = (loop.sol==SOL_WATERWATER) ? WATERATOMS : 1;
}


void make_func(char *appendname)
{
  set_loop_options();
  
  /* make name and information */  
  func_header(appendname);
  /* function call arguments */
  func_args();
  /* our local variables */
  func_localvars();
  /* initiation of local variables */
  func_init_vars();
  /* start the loop over i particles (neighborlists) */    
  outer_loop();
  /* The innerloop creation is called from the outerloop creation */
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined DOUBLE && defined USE_X86TRUNC)
  if(bC && DO_TAB && DO_FORCE)
    strcat(codebuffer,"asm(\"fldcw %%0\" : : \"m\" (*&x86_cwsave));\n");
#endif
  close_func();
}


int main(int argc,char *argv[])
{
 
  int i,j,nfunc,nlines,maxfunc;
  bool bSep14;
  char fn[32];
  
  if((argc==2) && !strcmp(argv[1],"c"))
    bC=TRUE;
  else if((argc==2) && !strcmp(argv[1],"fortran"))
    bC=FALSE;
  else {
    fprintf(stderr,"Usage: %s language\n"
                   "Currently supported languages:  c  fortran\n",argv[0]); 
    exit(-1);
  }

  sprintf(fn,"%s",bC ? "innerc.c" : "innerf.f");
  
  fprintf(stderr,">>> This is the GROMACS code generator for MD & MC inner loops\n"
	  ">>> It will generate %s precision %s code in file %s\n",
	  (prec == 8) ? "double" : "single",
	  bC ? "C" : "Fortran 77",fn);

#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined DOUBLE && defined USE_X86TRUNC)
  fprintf(stderr,">>> Using fast inline assembly gcc/x86 truncation. Since we are changing\n"
	         "    the control word this might affect the numerical result slightly.\n");
#endif  
  
#ifdef SOFTWARE_INVSQRT
  arch.gmx_invsqrt = TRUE;
  fprintf(stderr,">>> Using gromacs invsqrt code\n");
#else
  arch.gmx_invsqrt = FALSE;
#endif
    
#if ((defined SOFTWARE_INVSQRT) && !defined DONT_INLINE)
  opt.inline_invsqrt = TRUE;
  fprintf(stderr,">>> Inlining gromacs invsqrt code\n");
#else
  opt.inline_invsqrt = FALSE;  
#endif

#ifdef SIMPLEWATER
  arch.simplewater = TRUE;
  fprintf(stderr,">>> Will expand the solvent optimized loops to have 3 inner loops\n");
#else
  fprintf(stderr,">>> Using normal solvent optimized loops\n");
  arch.simplewater = FALSE;
#endif

  fprintf(stderr,">>> Prefetching forces in loops: ");
#if (!defined PREFETCH) && (!defined PREFETCH_S) && \
    (!defined PREFETCH_W) && (!defined PREFETCH_WW)
  fprintf(stderr,"none");
#endif
  opt.prefetch=0;
#ifdef PREFETCH
  fprintf(stderr,"normal ");
  opt.prefetch |= TYPE_NORMAL;
#endif
#ifdef PREFETCH_S
  fprintf(stderr,"solvent ");
  opt.prefetch |= TYPE_SOLVENT;
#endif
#ifdef PREFETCH_W
  fprintf(stderr,"water ");
  opt.prefetch |= TYPE_WATER;
#endif
#ifdef PREFETCH_WW
  fprintf(stderr,"water-water ");
  opt.prefetch |= TYPE_WATERWATER;
#endif
  fprintf(stderr,"\n");  

#ifdef DECREASE_LOOKUP_LATENCY
  opt.decrease_lookup_latency = TRUE;
  fprintf(stderr,">>> Will try to decrease table lookup latencies\n");
#else
  opt.decrease_lookup_latency = FALSE;
#endif

#ifdef THREADS_MUTEX
  arch.threads=THREAD_MUTEX;
  fprintf(stderr,">>> Creating thread parallel inner loops\n");
  fprintf(stderr,">>> Synchronizing threads with mutex\n");
#elif defined THREADS_STRIDE
  arch.threads=THREAD_STRIDE;
  fprintf(stderr,">>> Creating thread parallel inner loops\n");
  fprintf(stderr,">>> Striding loops over threads\n");
#else  /* no threads */
  arch.threads = THREAD_NO;
  fprintf(stderr,">>> Nonthreaded inner loops\n");
#endif

#ifdef __sgi
  opt.delay_invsqrt = TRUE;
  fprintf(stderr,">>> Delaying invsqrt and reciprocal\n");
#else
  opt.delay_invsqrt = FALSE;
#endif
  

  fprintf(stderr,">>> Vectorizing invsqrt in loops:");
#if (!defined VECTORIZE_INVSQRT) && (!defined VECTORIZE_INVSQRT_S) && \
    (!defined VECTORIZE_INVSQRT_W) && (!defined VECTORIZE_INVSQRT_WW)
  fprintf(stderr,"none");
#endif
  opt.vectorize_invsqrt=0;
#ifdef VECTORIZE_INVSQRT
  fprintf(stderr,"normal ");
  opt.vectorize_invsqrt |= TYPE_NORMAL;
#endif
#ifdef VECTORIZE_INVSQRT_S
  fprintf(stderr,"solvent ");
  opt.vectorize_invsqrt |= TYPE_SOLVENT;
#endif
#ifdef VECTORIZE_INVSQRT_W
  fprintf(stderr,"water ");
  opt.vectorize_invsqrt |= TYPE_WATER;
#endif
#ifdef VECTORIZE_INVSQRT_WW
  fprintf(stderr,"water-water ");
  opt.vectorize_invsqrt |= TYPE_WATERWATER;
#endif
  fprintf(stderr,"\n");  


  
#ifdef VECTORIZE_RECIP
  opt.vectorize_recip = TRUE;
  fprintf(stderr,">>> Vectorizing the reciprocal calculation\n");
#else
  opt.vectorize_recip = FALSE;
#endif

#ifdef USE_VECTOR
  arch.vectorcpu = TRUE;
  fprintf(stderr,">>> Will generate better vectorizable code (hopefully)\n");
#else
  arch.vectorcpu = FALSE;
#endif 
  
#ifdef CRAY_PRAGMA
  arch.cray_pragma = TRUE;
  fprintf(stderr,">>> Using cray-compatible pragma\n");
#else
  arch.cray_pragma = FALSE;
#endif
 
  if(arch.vectorcpu && arch.threads) {
    fprintf(stderr,"Error: Can't use threads on a vector architecture\n");
    exit(-1);
  }
  if(arch.vectorcpu || opt.prefetch) {
    fprintf(stderr,"Error: Prefetching on a vector architecture is bad\n");
    exit(-1);
  }
  if(opt.delay_invsqrt && arch.gmx_invsqrt) {
    fprintf(stderr,
	    "Warning: Can't delay invsqrt when using gromacs code.\n"
	    "         Turning off the delayed invsqrt\n");
    opt.delay_invsqrt=FALSE;
  }
 
  if ((output = fopen(fn,"w")) == NULL) {
    file_error(fn);
    exit(1);
  }

  if(opt.vectorize_invsqrt || arch.threads) {
    bSep14=TRUE;
    maxfunc=162;
  } else {
    bSep14=FALSE;
    maxfunc=156;
    /* there are both MC and MD loops now */
  }
  
  /* Allocate memory for the temporary code and variable buffers,
     and set indentation */
  init_metacode();
  /* Add include files, edit warning and fortran invsqrt routines */
  file_header();
  nfunc=0;

  /* Loop over all combinations to construct innerloops */
  loop.do_force=TRUE;
  for(i=0;i<2;i++) { /* force and non-force loops */
    for(loop.coul=COUL_NO;loop.coul<COUL_NR;loop.coul++)
      for(loop.vdw=VDW_NO;loop.vdw<VDW_NR;loop.vdw++)
	for(loop.sol=SOL_NO;loop.sol<SOL_NR;loop.sol++)
	  for(loop.free=FREE_NO;loop.free<FREE_NR;loop.free++) {
	    /* Exclude some impossible or unnecessary
	       combinations. */

	    if(!loop.coul && !loop.vdw)
	      continue;   /* no interactions at all to calculate */
	    if(loop.free &&
	       ((loop.coul && !DO_COULTAB) ||
		(loop.vdw && !DO_VDWTAB) ||
		DO_SOL || DO_WATER))
	      continue;  /* Always tabulate free energy, w/o solvent
			    opt. */
	    if(DO_WATER && !loop.coul)
	      continue;    /* No point with LJ only water loops */	      
	    /* (but we have LJ only general solvent loops) */

	    /* Write metacode to buffer */
	    make_func("");
	    flush_buffers();
	    nfunc++;
	    fprintf(stderr,"\rProgress: %2d%%",100*nfunc/maxfunc);
	  }
    loop.do_force=FALSE; /* do the non-force loops next round */
  }
  /* And maybe the extra, special 1-4 loops without fancy optimizations */

  if(bSep14) {
    opt.vectorize_invsqrt=FALSE; 
    arch.threads=THREAD_NO;
    loop.coul=COUL_TAB;
    loop.vdw=VDW_TAB;
    loop.sol=SOL_NO;
    loop.do_force=TRUE;
    for(i=0;i<2;i++) {
      for(loop.free=FREE_NO;loop.free<FREE_NR;loop.free++) {
	make_func("n");
	flush_buffers();
	nfunc++;
	fprintf(stderr,"\rProgress: %2d%%",100*nfunc/maxfunc);
      }
      loop.do_force=FALSE;
    }
  }
 
  fclose(output);
  nlines = count_lines(fn);
  fprintf(stderr,"\r>>> Generated %d lines of code in %d functions\n",nlines,nfunc);
  fprintf(stderr,"(%s may take a while to compile)\n\n",  bC ? "innerc.c" : "innerf.f");
  return 0;
}

