static char *SRCID_mkinl_innerloop_c = "";

#include "mkinl.h"

void unpack_inner_data(bool calcdist,bool calcforce)
{
    int nflop = 0;
    
    comment("Update indices");

    if(DO_PREFETCH_X) {
      assign("jnr","nextjnr");
      assign("nextjnr","%s%s", ARRAY(jjnr,k), bC ? "" : "+1");
	
      if(calcdist) {
	assign("j3","nextj3");
	assign("nextj3", bC ? "3*nextjnr" : "3*nextjnr-2");
      } else if(calcforce)  /* only force, no need for j prefetch? */
	assign("j3", bC ? "3*jnr" : "3*jnr-2");

    } else { /* no prefetch */
      assign("jnr", "%s%s", ARRAY(jjnr,k), bC ? "" : "+1");
      assign("j3",bC ? "3*jnr" : "3*jnr-2");
    }

    if(arch.vector && calcforce)    /* prefetch x is turned off in this case */
      assign("kk","%d*(k-nj0)%s",(loop.sol==SOL_WATERWATER) ? 9 : 3 ,bC ? "" : "+1");
    
}


void move_coord(bool calcdist, bool calcforce)
{
  int j;
  
  comment("Moving prefetched coords to present ones");

  if(calcdist) 
    for(j=1;j<=loop.nj;j++) {
      assign("jx%d","nextjx%d",j,j);
      assign("jy%d","nextjy%d",j,j);
      assign("jz%d","nextjz%d",j,j);
    }

  if(loop.sol!=SOL_WATERWATER) {
    if(calcforce) {
      if(loop.coul) {
	assign("jqA","nextjqA");
	if(loop.free)
	  assign("jqB","nextjqB");
      }
      if(loop.vdw || DO_SOFTCORE) {
	assign("tpA","nexttpA");
	if(loop.free)
	  assign("tpB","nexttpB");
      }
    }
  }
  
}


void fetch_coord(bool calcdist, bool calcforce)
{
    int j,offset;
    char jjname[16]; /* (next)j3 */
    char jname[16];  /* (next)jnr */
    
    if(DO_PREFETCH_X) 
      move_coord(calcdist,calcforce);    
   
    comment("Fetching coordinates, charge and type");

    for(j=1;j<=loop.nj;j++) {
      
      offset=3*(j-1);
      
      sprintf(jjname,"%sj3", DO_PREFETCH_X ? "next" : "");
      sprintf(jname,"%sjnr", DO_PREFETCH_X ? "next" : ""); 

      if(calcdist) {
	assign( DO_PREFETCH_X ? "nextjx%d" : "jx%d",_array("pos","%s+%d",jjname,offset),j);
	assign( DO_PREFETCH_X ? "nextjy%d" : "jy%d",_array("pos","%s+%d",jjname,offset+1),j);
	assign( DO_PREFETCH_X ? "nextjz%d" : "jz%d",_array("pos","%s+%d",jjname,offset+2),j);
      }
      
    }
    
    if(calcforce && loop.sol!=SOL_WATERWATER) {
      if(loop.coul) {	
	assign( DO_PREFETCH_X ? "nextjqA" : "jqA", _array("charge",jname),j);
	if(loop.free)
	  assign( DO_PREFETCH_X ? "nextjqB" : "jqB", _array("chargeB",jname), j);
      }
      if(loop.vdw || DO_SOFTCORE) {
	assign( DO_PREFETCH_X ? "nexttpA" : "tpA", _array("type",jname), j);
	if(loop.free)
	  assign( DO_PREFETCH_X ? "nexttpB" : "tpB", _array("typeB",jname), j);
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


void prefetch_first_data(bool calcdist,bool calcforce) 
{
  int j;
  char buf[50];  
 
  comment("Prefetch first round of x data");
  
  assign("nextjnr","%s%s", ARRAY(jjnr,nj0), bC ? "" : "+1");

  if(calcdist) {
      assign("nextj3", bC ? "3*nextjnr" : "3*nextjnr-2");  
      for(j=1;j<=loop.nj;j++) {
	assign("nextjx%d",_array("pos","nextj3+%d",3*(j-1)),j);
	assign("nextjy%d",_array("pos","nextj3+%d",3*(j-1)+1),j);
	assign("nextjz%d",_array("pos","nextj3+%d",3*(j-1)+2),j);
      }
  }
  
  if(calcforce && loop.sol!=SOL_WATERWATER) {
    if(loop.coul) {
      assign("nextjqA",ARRAY(charge,nextjnr));
      if(loop.free)
	assign("nextjqB",ARRAY(chargeB,nextjnr));	    
    }
    if(loop.vdw || DO_SOFTCORE) {
      assign("nexttpA", ARRAY(type,nextjnr));
      if(loop.free)
	assign("nexttpB", ARRAY(typeB,nextjnr));	      
    }
  }
}


void last_inner_iteration(bool calcdist, bool calcforce)
{
    int i;

    if(DO_PREFETCH_X) {
      comment("Treat the last round of prefetched data outside the loop");
      if(calcforce) {	
	assign("jnr","nextjnr");
	assign("j3", calcdist ? "nextj3" : "3*jnr");
      }
      
      move_coord(calcdist,calcforce);
      if(calcdist)
	calc_dist();
      
      if(calcforce) {
	if(DO_PREFETCH_F) 
	  prefetch_forces();
	
	calc_rinv_and_rinvsq();
	calc_interactions();
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
  
  if(DO_PREFETCH_X)
    prefetch_first_data(calcdist,calcforce);
  
  vector_pragma();
  
  start_loop("k","nj0","nj1");
  
  unpack_inner_data(calcdist,calcforce);

  /* this also fetches charge and LJ parameters */
  fetch_coord(calcdist,calcforce);
  
  if(calcdist)
    nflop += calc_dist(); 
  
  if(calcforce) { 
    if(DO_PREFETCH_F) 
      prefetch_forces();
    
    nflop += calc_rinv_and_rinvsq();
    
    nflop += calc_interactions();
    
    /* Only print vanilla loop count */
    if(!DO_VECTORIZE && !DO_PREFETCH_X && !arch.sse_and_3dnow && !arch.vector) {
      sprintf(buf,"Innerloop of %s costs %d flops",loopname,nflop);
      comment(buf);
    }
  }
  end_loop();

  last_inner_iteration(calcdist,calcforce);
}
