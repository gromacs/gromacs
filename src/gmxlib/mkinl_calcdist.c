static char *SRCID_mkinl_calcdist_c = "";
#include "mkinl.h"
#include <string.h>


void call_vectorized_routines()
{
  /* m was updated a last time after the assignments. Reduce
   * it if we are using fortran, so we avoid referencing m-1
   * repeatedly (and some functions need an adress to a variable
   * rather than a value)
   */
  if(!bC)
    decrement("m","1");
  
  /* some machines perform better if the can vectorize the
   * square distance calculation. */
  if(opt.decrease_square_latency) {
    assign("m3", bC ? "0" : "1");
    start_loop("k", bC ? "0" : "1" ,"m");
    assign("dx11",ARRAY(drbuf,m3));
    assign("dy11",ARRAY(drbuf,m3+1));
    assign("dz11",ARRAY(drbuf,m3+2));
    assign(ARRAY(buf1,k),"dx11*dx11+dy11*dy11+dz11*dz11");
    increment("m3","3");
    end_loop();
  }
  
  /* always start with vectorized invsqrt if we should do it
   */
  if(loop.vectorize_invsqrt) {
    /* if we need the sqrt (non-reciprocal) later
     * we have to save the square values in buf1 and write to buf2,
     * otherwise we can gain some speed by writing back to the
     * same array buf1
     */
#ifdef USE_IBM_MASS
    /* IBM provide special fast single and double routines in their MASS libraries */
#ifdef DOUBLE
    code( bC ? "vrsqrt(%s,buf1,&m);" : "call vrsqrt(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );	  
#else
    code( bC ? "vsrsqrt(%s,buf1,&m);" : "call vsrsqrt(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif /* End of IBM routines */
#elif defined(USE_AXP_ASM)
    /* Compaq provide a math library too, but not an inverse square root or
     * a fast reciprocal. Our own assembly routines are much faster! */
#ifdef DOUBLE
    code( bC ? "sqrtiv_(buf1,%s,&m);" : "call sqrtiv(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#else
    code( bC ? "ssqrtiv_(buf1,%s,&m);" : "call ssqrtiv(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif /* End of code for compaq alpha */
#else 
    /* this is a vanilla GMX routine, suitable for most chips w/o hardware sqrt.
     * It is quite fast though - about 18 clocks on intel! */
    if(arch.gmx_invsqrt)
      code( bC ? "vecinvsqrt(buf1,%s,m);" : "call vecinvsqrt(buf1,%s,m)",
	    OVERWRITE_RSQ ? "buf1" : "buf2" );
    else {
      /* on some architectures, though, the compiler can 
       * link with appropriate libs if only the code is easily recognized
       * as vectorizable. Do a simple loop over the vector:
       */
      start_loop("n", bC ? "0" : "1" ,"m");
      assign(OVERWRITE_RSQ ? ARRAY(buf1,n) : ARRAY(buf2,n) , "1.0/sqrt(%s)",ARRAY(buf1,n));
      end_loop();
    }
#endif
  } else if(loop.vectorize_recip) {
      /* Only do this if we didnt do invsqrt. If we did,
       * it is probably faster to square the value when
       * we use it.
       */
#if defined(USE_IBM_MASS)
    /* IBM stuff again. They provide a nice reciprocal */
#ifdef DOUBLE
    code( bC ? "vrec(%s,buf1,&m);" : "call vrec(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#else
    code( bC ? "vsrec(%s,buf1,&m);" : "call vsrec(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif
#else
	/* on some architectures, though, the compiler can 
	 * link with appropriate libs if the code is easily recognized as
	 * vectorizable.
	 */
    vector_pragma();
    start_loop("n", bC ? "0" : "1" ,"m");
    assign(OVERWRITE_RSQ ? ARRAY(buf1,n) : ARRAY(buf2,n),
	   "1.0/%s",ARRAY(buf1,n));
    end_loop();
#endif
  }
}


int calc_rsq(int i, int j)
{
  char lm[12],lm2[32],sqbuf[64];
  int  m;

  sqbuf[0]=0;
  for(m='x'; (m<='z'); m++) {
    sprintf(lm,"d%c%d%d",m,i,j);
    sprintf(lm2,"%s*%s",lm,lm);
    add_to_buffer(sqbuf,lm2);
  }
  if(DO_VECTORIZE) {
    assign(ARRAY(buf1,m),sqbuf);
    increment("m","1");
  } else 
    assign("rsq%d%d",sqbuf,i,j);
  return 5;
}


int calc_dist()
{
  int nflop=0;
  int i,j,m;
  char buf[10];
  /* First calculate dr, and possibly store it */

  comment("Calculate distance vector"); 

  for(i=1;i<=loop.ni;i++)
    for(j=1;j<=loop.nj;j++) {
      /* first calculate dr */
      for(m='x'; (m<='z'); m++) { 
	subtract("d%c%d%d","i%c%d","j%c%d",m,i,j,m,i,m,j);
	nflop ++;
	if(DO_VECTORIZE) {
	  assign(ARRAY(drbuf,m3),"d%c%d%d",m,i,j);
	  increment("m3","1");
	}
      }
      if(!opt.decrease_square_latency)
	nflop += calc_rsq(i,j);
      else if(DO_VECTORIZE) /* calc square separately later, but we  */
	increment("m","1"); /* need to know the number of items      */
    }
  
  if(opt.decrease_square_latency && !DO_VECTORIZE)
    for(i=1;i<=loop.ni;i++)
      for(j=1;j<=loop.nj;j++)
	nflop += calc_rsq(i,j);

  return nflop;
}


int calc_rinv_and_rinvsq()
{
  int nflop=0;
  int i,j;
  
  if(!DO_VECTORIZE && !opt.delay_invsqrt) {
    if(loop.invsqrt)
      nflop += calc_invsqrt();
    else if(loop.recip) {
      for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++)
	  assign("rinvsq%d", "1.0/rsq%d",i*10+j,i*10+j);
      if(loop.sol==SOL_WATER)
	nflop=3;
      else if(loop.sol==SOL_WATERWATER)
	nflop=9;
      else
	nflop=1;
    }
  }
  return nflop;
}
