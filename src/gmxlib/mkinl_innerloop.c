static char *SRCID_mkinl_innerloop_c = "";

#include "mkinl.h"

void unpack_inner_data(bool calcdist,bool calcforce)
{
    int nflop = 0;
    
    comment("Update indices");

    if(DO_PREFETCH_X && !DO_VECTORIZE) {
      assign("jnr","nextjnr");
	
      if(calcdist) {
	assign("j3","nextj3");
      } else if(calcforce)  /* only force, no need for j prefetch? */
	assign("j3", bC ? "3*jnr" : "3*jnr-2");

    } else { /* no prefetch, or vectorized */
      assign("jnr", "%s%s", ARRAY(jjnr,k), bC ? "" : "+1");
      assign("j3",bC ? "3*jnr" : "3*jnr-2");
    }

    if(arch.vector && calcforce)    /* prefetch x is turned off in this case */
      assign("kk","%d*(k-nj0)%s",(loop.sol==SOL_WATERWATER) ? 9 : 3 ,bC ? "" : "+1");
    
}



void move_coords(void)
{
  int j;

  /* move prefetched coords into real ones */
  for(j=1;j<=loop.nj;j++) {
    assign("jx%d","nextjx%d",j,j);
    assign("jy%d","nextjy%d",j,j);
    assign("jz%d","nextjz%d",j,j);
  }
}


void fetch_coord(bool bPrefetch)
{
  int j,offset;

  comment("Fetching coordinates");
  
  if(DO_VECTORIZE && bPrefetch) {
    /* update indices */
    assign("nextjnr","%s%s", ARRAY(jjnr,nj0), bC ? "" : "+1");
      assign("nextj3", bC ? "3*nextjnr" : "3*nextjnr-2");  
    
      move_coords();
      /* get new coords */
      for(j=1;j<=loop.nj;j++) {
	offset=3*(j-1);
	assign("nextjx%d",_array("pos", "j3+%d",offset),j);
	assign("nextjy%d",_array("pos", "j3+%d",offset+1),j);
	assign("nextjz%d",_array("pos", "j3+%d",offset+2),j);
      }
      
  } else {
    if(bPrefetch)
      assign("nextj3", bC ? "3*nextjnr" : "3*nextjnr-2");  

    for(j=1;j<=loop.nj;j++) {
      offset=3*(j-1);
      assign("jx%d",_array("pos", bPrefetch ? "nextj3+%d" : "j3+%d",offset),j);
      assign("jy%d",_array("pos", bPrefetch ? "nextj3+%d" : "j3+%d",offset+1),j);
      assign("jz%d",_array("pos", bPrefetch ? "nextj3+%d" : "j3+%d",offset+2),j);
    }
  }
}
	

void prefetch_forces()
{

  int j;
  
  comment("Prefetching forces");
  
  for(j=1;j<=loop.nj;j++) {
    assign("fjx%d", _array("faction","j3+%d",3*(j-1)),j);
    assign("fjy%d", _array("faction","j3+%d",3*(j-1)+1),j);
    assign("fjz%d", _array("faction","j3+%d",3*(j-1)+2),j);
  }
  /* no flops */
}



void unpack_vector_machine_forces(bool calcdist,bool calcforce)
{
  /* this shouldnt be used with prefetching */
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


void prefetch_first_coord(void) 
{
  int j;
  char buf[50];  
 
  comment("Prefetch first round of x data");
  
  /* When using prefetching of coordinates in the vectorized loop
   * we need the dummy storage variables nextjx, etc.
   */
  if(DO_VECTORIZE) {
    assign("jnr","%s%s", ARRAY(jjnr,nj0), bC ? "" : "+1");
    assign("j3", bC ? "3*jnr" : "3*jnr-2");  
    for(j=1;j<=loop.nj;j++) {
      assign("nextjx%d", _array("pos","j3+%d",3*(j-1)),j);
      assign("nextjy%d", _array("pos","j3+%d",3*(j-1)+1),j);
      assign("nextjz%d", _array("pos","j3+%d",3*(j-1)+2),j);
    }
  } else {
    assign("nextjnr","%s%s", ARRAY(jjnr,nj0), bC ? "" : "+1");
    assign("nextj3", bC ? "3*nextjnr" : "3*nextjnr-2");  
    for(j=1;j<=loop.nj;j++) {
      assign("jx%d",_array("pos","nextj3+%d",3*(j-1)),j);
      assign("jy%d",_array("pos","nextj3+%d",3*(j-1)+1),j);
      assign("jz%d",_array("pos","nextj3+%d",3*(j-1)+2),j);
    }
  }
}


void last_inner_iteration(bool calcdist,bool calcforce)
{
    int i;

    if(DO_PREFETCH_X) {
      comment("Treat the last round of prefetched data outside the loop");
      if(!DO_VECTORIZE) {	
	assign("jnr","nextjnr");
	assign("j3", calcdist ? "nextj3" : "3*jnr");
      } 

      if(calcdist) {
	if(DO_VECTORIZE)
	  move_coords();
	calc_dist();
      }

      if(calcforce) {
	calc_rinv_and_rinvsq();

	if(DO_PREFETCH_F) 
	  prefetch_forces();

	calc_interactions(FALSE);
      }
    }
}


void inner_loop(bool calcdist, bool calcforce)
{
  int nflop=0;
  char buf[50],buf2[50];
  
  if(arch.vector && calcforce)
    unpack_vector_machine_forces(calcdist,calcforce);
  
  comment("Inner loop (over j-particles) starts right here");
  
  if(DO_PREFETCH_X && calcdist)
    prefetch_first_coord();
  
  vector_pragma();
  
  start_loop("k", (DO_PREFETCH_X && calcdist) ? "nj0+1" : "nj0" , "nj1");
  
  unpack_inner_data(calcdist,calcforce);

  if((!DO_PREFETCH_X || DO_VECTORIZE) && calcdist)
    fetch_coord(DO_PREFETCH_X);
  
  /* non-vectorized coordinate prefetching is issued from 
   * calc_interactions() in mkinl_interactions.c due to 
   * timing optimization.
   */

  if(calcdist)
    nflop += calc_dist(); 
  
  if(calcforce) { 
    
    nflop += calc_rinv_and_rinvsq();
    
    if(DO_PREFETCH_F) 
      prefetch_forces();

    if(DO_PREFETCH_X && calcdist)
      assign("nextjnr","%s%s", ARRAY(jjnr,k), bC ? "" : "+1");

    nflop += calc_interactions(calcdist && DO_PREFETCH_X);
    
    /* Only print vanilla loop count */
    if(!DO_VECTORIZE && !DO_PREFETCH_X && !arch.sse_and_3dnow && !arch.vector) {
      sprintf(buf,"Innerloop of %s costs %d flops",loopname,nflop);
      comment(buf);
    }
  }
  end_loop();

  if(DO_PREFETCH_X && calcdist)
    last_inner_iteration(calcdist,calcforce);
}
