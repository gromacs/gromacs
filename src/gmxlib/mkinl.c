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
 * Green Red Orange Magenta Azure Cyan Skyblue
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
 * -DFAST_X86TRUNC     Change the x86 control word once in each
 *                     routine to speed up truncation. This will
 *                     increase table routine performance significantly,
 *                     but it might also affect the floating point
 *                     precision. Check your output to see if you
 *                     can live with it.
 *
 * -DUSE_SSE_AND_3DNOW Include assembly loops which can utilize fast 
 *                     single precision optimizations in newer Intel and
 *                     Amd processors. The code will still run on any
 *                     x86 machine, but to make use of the special loops you
 *                     will need e.g. redhat>=6.2 and a pentium III/amd k6 
 *                     or later running a kernel with support for
 *                     these instructions. It only affects single precision.
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
 * -DGMX_INVSQRT       Use the gromacs implementation of the inverse
 *                     square root.
 *                     This is probably a good thing to do on most 
 *                     computers except SGI and ibm power2/power3 which
 *                     have special hardware invsqrt instructions.
 *                     (but the PPC604 available in some ibm smp nodes doesnt)
 *
 * -DGMX_RECIP         Same as invsqrt, but this is for 1/x. This
 *                     should NOT be used without checking - most 
 *                     computers have a more efficient hardware
 *                     routine, but it may help on SOME x86 boxes.
 *
 * INNERLOOP OPTIMIZATION OPTIONS:  
 * 
 * -DINLINE_GMXCODE    This inlines the gromacs inverse square root 
 *                     and/or reciprocal on architectures which use
 *                     it, saving some calling time in the inner loops.  
 *
 * -DPREFETCH_F        Prefetching of forces. The two extra options 
 * -DPREFETCH_F_W      additionally does it for single water and 
 * -DPREFETCH_F_WW     water-water loops.
 * 
 * -DPREFETCH_X        Prefetches coordinates.
 * -DPREFETCH_X_W      In general, more prefetching is better, but on
 * -DPREFETCH_X_WW     processors with few registers (x86) it can decrease
 *                     performance when we run out of them.
 *                     This is even more probable for the water-water loops,
 *                     so they have special options.
 *
 * -DVECTORIZE_INVSQRT    Vectorize the corresponding part of the force
 * -DVECTORIZE_INVSQRT_W  calculation. This can avoid cache trashing, and
 * -DVECTORIZE_INVSQRT_WW in some cases it makes it possible to use vector
 * -DVECTORIZE_RECIP      intrinsic functions. (also on scalar machines)
 *                        The two extra invsqrt options control water
 *                        loops in the same way as the prefetching.           
 *
 * -DDECREASE_SQUARE_LATENCY  Try to hide latencies by doing some operations
 * -DDECREASE_LOOKUP_LATENCY  in parallel rather than serial in the innermost loop.
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
  loop.coul_needs_rinvsq=(loop.coul && !DO_COULTAB);
  /* non-table coul also need r^-2 */
  loop.coul_needs_rsq=(DO_RF || (loop.coul && loop.free));
  /* reaction field and softcore need r^2 */
  loop.coul_needs_r=DO_COULTAB;
  /* tabulated coulomb needs r */ 

  /* If we are vectorizing the invsqrt it is probably faster to
   * include the vdw-only atoms in that loop and square the
   * result to get rinvsq. Just tell the generator we want rinv for
   * those atoms too:
   */
  loop.vdw_needs_rinv=DO_VDWTAB || (loop.coul && loop.vectorize_invsqrt);
  
  /* table nb needs 1/r */
  loop.vdw_needs_rinvsq=(loop.vdw && !DO_VDWTAB);
  /* all other need r^-2 */
  loop.vdw_needs_rsq=(loop.vdw && loop.free);
  /* softcore needs r^2 */
  loop.vdw_needs_r=(DO_BHAM || DO_VDWTAB || (loop.vdw && loop.free));
  /* softcore and buckingham need r */

  /* Now, what shall we calculate? */ 
  loop.invsqrt=
    (loop.coul && (loop.coul_needs_rinv || loop.coul_needs_r)) ||
    (loop.vdw && (loop.vdw_needs_rinv || loop.vdw_needs_r));
  
  loop.recip=!loop.invsqrt;

  loop.vectorize_invsqrt=loop.invsqrt && 
    (opt.vectorize_invsqrt==YES_WW || 
     (opt.vectorize_invsqrt==YES_W && loop.sol!=SOL_WATERWATER) ||
     (opt.vectorize_invsqrt==YES && !DO_WATER));

  loop.vectorize_recip=loop.recip && opt.vectorize_recip;
     
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
#if (defined __GNUC__ && defined _lnx_ && defined FAST_X86TRUNC)
  if(bC && DO_TAB)
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
  
  fprintf(stderr,">>> This is the GROMACS innerloop code generator\n"
	  ">>> It will generate %s precision %s code in file %s\n",
	  (prec == 8) ? "double" : "single",
	  bC ? "C" : "Fortran 77",fn);

#ifdef USE_SSE_AND_3DNOW
  arch.sse_and_3dnow=FALSE; /* no support - yet */
  fprintf(stderr,">>> Including x86 assembly loops with SSE and 3DNOW instructions\n");
#else
  arch.sse_and_3dnow=FALSE;
#endif

#if (defined __GNUC__ && defined _lnx_)
#ifdef FAST_X86TRUNC
  fprintf(stderr,">>> Using fast inline assembly gcc/x86 truncation. Since we are changing\n"
	         "    the control word this might affect the numerical result slightly.\n");
#else
  fprintf(stderr,">>> Using normal x86 truncation\n");
#endif  
#endif  
  
#ifdef GMX_INVSQRT
  arch.gmx_invsqrt = TRUE;
  fprintf(stderr,">>> Using gromacs invsqrt code\n");
#else
  arch.gmx_invsqrt = FALSE;
#endif
    
#ifdef GMX_RECIP
  arch.gmx_recip = TRUE;
  fprintf(stderr,">>> Using gromacs reciprocal code\n");
#else
  arch.gmx_recip = FALSE;
#endif
    
#if ((defined GMX_INVSQRT || defined GMX_RECIP) && defined INLINE_GMXCODE)
  opt.inline_gmxcode = TRUE;
  fprintf(stderr,">>> Inlining gromacs invsqrt and/or reciprocal code\n");
#else
  opt.inline_gmxcode = FALSE;  
#endif

#ifdef SIMPLEWATER
  arch.simplewater = TRUE;
  fprintf(stderr,">>> Will expand the solvent optimized loops to have 3 inner loops\n");
#else
  fprintf(stderr,">>> Using normal solvent optimized loops\n");
  arch.simplewater = FALSE;
#endif

#ifdef PREFETCH_X 
  opt.prefetch_x = YES;
  fprintf(stderr,">>> Prefetching coordinates\n");
#elif defined PREFETCH_X_W
  opt.prefetch_x = YES_W;
  fprintf(stderr,">>> Prefetching coordinates (also single water loops)\n");
#elif defined PREFETCH_X_WW
  opt.prefetch_x = YES_WW;
  fprintf(stderr,">>> Prefetching coordinates (all loops)\n");
#else
  opt.prefetch_x = NO;
#endif

#ifdef PREFETCH_F
  opt.prefetch_f = YES;
  fprintf(stderr,">>> Prefetching forces\n");
#elif defined PREFETCH_F_W
  opt.prefetch_f = YES_W;
  fprintf(stderr,">>> Prefetching forces (also single water loops)\n");
#elif defined PREFETCH_F_WW
  opt.prefetch_f = YES_WW;
  fprintf(stderr,">>> Prefetching forces (all loops)\n");
#else
  opt.prefetch_f = NO;
#endif

#ifdef DECREASE_SQUARE_LATENCY
  opt.decrease_square_latency = TRUE;
  fprintf(stderr,">>> Will try to decrease square distance calculation latency\n");
#else
  opt.decrease_square_latency = FALSE;
#endif

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

#ifdef _SGI_
  opt.delay_invsqrt = TRUE;
  fprintf(stderr,">>> Delaying invsqrt and reciprocal\n");
#else
  opt.delay_invsqrt = FALSE;
#endif 

#ifdef VECTORIZE_INVSQRT
  opt.vectorize_invsqrt = YES;
  fprintf(stderr,">>> Vectorizing the invsqrt calculation\n");
#elif defined VECTORIZE_INVSQRT_W
  opt.vectorize_invsqrt = YES_W;
  fprintf(stderr,">>> Vectorizing the invsqrt calculation (also single water loops)\n");
#elif defined VECTORIZE_INVSQRT_WW
  opt.vectorize_invsqrt = YES_WW;
  fprintf(stderr,">>> Vectorizing the invsqrt calculation (all loops)\n");
#else
  opt.vectorize_invsqrt = NO;
#endif

#ifdef VECTORIZE_RECIP
  opt.vectorize_recip = TRUE;
  fprintf(stderr,">>> Vectorizing the reciprocal calculation\n");
#else
  opt.vectorize_recip = FALSE;
#endif

#ifdef USE_VECTOR
  arch.vector = TRUE;
  fprintf(stderr,">>> Will generate better vectorizable code (hopefully)\n");
#else
  arch.vector = FALSE;
#endif 
  
#ifdef CRAY_PRAGMA
  arch.cray_pragma = TRUE;
  fprintf(stderr,">>> Using cray-compatible pragma\n");
#else
  arch.cray_pragma = FALSE;
#endif
 
  if(arch.vector && arch.threads) {
    fprintf(stderr,"Error: Can't use threads on a vector architecture\n");
    exit(-1);
  }
  if(arch.vector && (opt.prefetch_x || opt.prefetch_f)) {
    fprintf(stderr,"Error: Prefetching on a vector architecture is bad\n");
    exit(-1);
  }
  if(opt.delay_invsqrt && arch.gmx_invsqrt) {
    fprintf(stderr,"Error: Can't delay invsqrt when using gromacs code.\n");
    exit(-1);
  }
  if(prec!=4 && arch.sse_and_3dnow) {
    fprintf(stderr,"Error: SSE/3DNOW loops can only be used in single precision\n");
    fprintf(stderr,"       (Not our choice - blame Intel and AMD.)\n");    
    exit(-1);
  }
 
  if ((output = fopen(fn,"w")) == NULL) {
    file_error(fn);
    exit(1);
  }

  if(opt.vectorize_invsqrt || arch.threads) {
    bSep14=TRUE;
    maxfunc=81;
  } else {
    bSep14=FALSE;
    maxfunc=78;
    /* 74 different loops present right now... */
  }
  
  /* Allocate memory for the temporary code and variable buffers,
     and set indentation */
  init_metacode();
  /* Add include files, edit warning and fortran invsqrt routines */
  file_header();
  nfunc=0;

  /* Loop over all combinations to construct innerloops */
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

  /* And maybe some extra, special 1-4 loops without some fancy optimizations */

  if(bSep14) {
    opt.vectorize_invsqrt=FALSE; 
    arch.threads=THREAD_NO;
    loop.coul=COUL_TAB;
    loop.vdw=VDW_TAB;
    loop.sol=SOL_NO;
    for(loop.free=FREE_NO;loop.free<FREE_NR;loop.free++) {
      make_func("n");
      flush_buffers();
      nfunc++;
      fprintf(stderr,"\rProgress: %2d%%",100*nfunc/maxfunc);
    }
  }
 
  fclose(output);
  nlines = count_lines(fn);
  fprintf(stderr,"\r>>> Generated %d lines of code in %d functions\n\n",nlines,nfunc);
  return 0;
}

