static char *SRCID_mkinl_invsqrt_c = "";
#include "mkinl.h"
#include <string.h>

void fortran_invsqrt_seed(void)
{
  comment("Invsqrt variables");

  declare_int4( (prec==8) ? "finvsqrtdoubletab" : "finvsqrtsingletab"); 

  declare_intreal("bval,tt");
  declare_real("fval");

  code("equivalence(bval,fval)");

  code( (prec==8) ? "common /finvsqrtdoubledata/ finvsqrtdoubletab(4096)" : "common /finvsqrtsingledata/ finvsqrtsingletab(4096)");
}



/* perform the first stage table lookup step of the invsqrt */
void invsqrt_table_lookup(char *text)
{
  assign( bC ? "conv%s.fval" : "fval%s","rsq%s",text,text);

  assign("tA%s", bC ? "(conv%s.bval >> 1 )" : "ishft(bval%s,-1)",text,text);
  assign("tB%s","tt - tA%s",text,text);

#ifdef DOUBLE
  assign("tC%s", bC ? "(tB%s >> 40)" : "ishft(tB%s,-40)",text,text);
#else
  assign("tC%s", bC ? "(tB%s >> 11)" : "ishft(tB%s,-11)",text,text);
#endif

  if(bC) {
      code("tC%s &= ((1 << 12)-1);",text);
      assign("tD%s", (prec==8) ? "cinvsqrtdoubletab[tC%s]" : "cinvsqrtsingletab[tC%s]",text,text);
  } else {
      assign("tC%s","iand(tC%s,ishft(1,12)-1)",text,text);
      assign("tD%s", (prec==8) ? "finvsqrtdoubletab(tC%s+1)" : "finvsqrtsingletab(tC%s+1)",text,text);
  }

#ifdef DOUBLE
  assign("tD%s", bC ? "(tD%s << 26)" : "ishft(tD%s,26)",text,text);
#endif
  assign( bC ? "conv%s.bval" : "bval%s","tB%s - tD%s",text,text,text);

  assign("lu%s", bC ? "conv%s.fval" : "fval%s",text,text);
}


int invsqrt_iterate(bool bFirst)
{
  int nflop=0;
  int i,j;
  bool repeat=FALSE;

#ifdef DOUBLE
  if(bFirst)
    repeat=TRUE;
#endif

  comment("N-R iteration to improve accuracy");
  newline();

  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("tmp%d","rsq%d","lu%d",i*10+j,i*10+j,i*10+j);

  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("tmp%d","tmp%d","lu%d",i*10+j,i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("tmp%d","tmp%d","half",i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      subtract("tmp%d","onehalf","tmp%d",i*10+j,i*10+j);
   
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply(repeat ? "lu%d" : "rinv%d","tmp%d","lu%d",i*10+j,i*10+j,i*10+j);

  if(repeat) 
    invsqrt_iterate(FALSE);
  
  return nflop;
}


void whole_invsqrt(char *name)
{
  invsqrt_table_lookup(name);
  multiply("tmp%s","rsq%s","lu%s",name,name,name);
  multiply("tmp%s","tmp%s","lu%s",name,name,name);
  multiply("tmp%s","tmp%s","half",name,name);
  subtract("tmp%s","onehalf","tmp%s",name,name);
#ifdef DOUBLE
  multiply("lu%s","tmp%s","lu%s",name,name,name);
  multiply("tmp%s","rsq%s","lu%s",name,name,name);
  multiply("tmp%s","tmp%s","lu%s",name,name,name);
  multiply("tmp%s","tmp%s","half",name,name);
  subtract("tmp%s","onehalf","tmp%s",name,name);
#endif
  multiply("rinv%s","tmp%s","lu%s",name,name,name);
}


int inlined_invsqrt(void)
{
  int nflop=0;
  int i,j,k;
  char buf[100];
  comment("Invsqrt table lookup");

  if(opt.decrease_lookup_latency) {
    for(j=1;j<=loop.nj;j++)
      for(i=1;i<=loop.ni;i++) {
	sprintf(buf,"%d",i*10+j);
	invsqrt_table_lookup(buf);
      }
    nflop+=invsqrt_iterate(TRUE);  
  } else {
    for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++) {
	  sprintf(buf,"%d",i*10+j);
	  whole_invsqrt(buf);
	}
  }

  return nflop;
}




int calc_invsqrt()
{
    int nflop=0;
    int i,j,tmpflop;

    /* check if it should be inlined */
    if(opt.inline_invsqrt)
      inlined_invsqrt();
    else {    /* otherwise do the ordinary one */
      for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++)
	  assign("rinv%d", arch.gmx_invsqrt ? 
		 "invsqrt(rsq%d)" : "1.0/sqrt(rsq%d)",i*10+j,i*10+j);
    }
#ifdef DOUBLE
    tmpflop=9;
#else
    tmpflop=5;
#endif
    
    if(!arch.gmx_invsqrt)
      tmpflop=1;
    if(loop.sol==SOL_NO)
      nflop=tmpflop;
    else if(loop.sol!=SOL_WATERWATER)
      nflop=tmpflop*3;
    else
      nflop=tmpflop*9;
    return nflop;
}



void invsqrt_vars()
{
  int i,j,ij;
  char buf[100];

  if(!bC) 
    declare_int4( (prec==8) ? "finvsqrtdoubletab" : "finvsqrtsingletab");
  
  declare_intreal("tt");  

  declare_const_real("onehalf",1.5);
  declare_const_real("half",0.5);  

  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) {
      ij=10*i+j;
      sprintf(buf,"lu%d",ij);
      declare_real(buf);
      sprintf(buf,"tmp%d",ij);
      declare_real(buf);
      sprintf(buf,"tA%d",ij);
      declare_intreal(buf);
      sprintf(buf,"tB%d",ij);
      declare_intreal(buf);
      sprintf(buf,"tC%d",ij);
      declare_intreal(buf);
      sprintf(buf,"tD%d",ij);
      declare_intreal(buf);

      if(bC) {
	sprintf(buf,"conv%d",ij);
	declare_other("t_convert",buf);
      } else {
	sprintf(buf,"bval%d",ij);
	declare_intreal(buf);
	sprintf(buf,"fval%d",ij);
	declare_real(buf);
	sprintf(buf,"equivalence(bval%d,fval%d)",ij,ij);
	code(buf);
      }
    }

  if(!bC)
      code( (prec==8) ? "common /finvsqrtdoubledata/ finvsqrtdoubletab(4096)" : 
	   "common /finvsqrtsingledata/ finvsqrtsingletab(4096)");

  if(bC) {
#ifdef DOUBLE
    assign("tt","0x5fe80000");
    assign("tt","tt<<32");
#else
    assign("tt","0x1fc00000");
#endif
  } else { /* f77 */
#ifdef DOUBLE
    assign("tt","X'5fe80000'");
    assign("tt","ishft(tt,32)");
#else
    assign("tt","X'1fc00000'");
#endif
  }

}


void fortran_invsqrt()
{
  newline();
  comment("fortran invsqrt routine");
  strcat(header,"      function invsqrt(");
  declare_real("x");
  nargs=ndecl; 
  declare_real("invsqrt");
  declare_int4( (prec==8) ? "finvsqrtdoubletab" : "finvsqrtsingletab");
  
  declare_intreal("bval");
  declare_real("fval");
  
  code("equivalence(bval,fval)");
  
  declare_intreal("tt");
  declare_intreal("tA");
  declare_intreal("tB");
  declare_intreal("tC");
  declare_intreal("tD");
  
  declare_real("f1");
  declare_real("y");
  declare_const_real("half",0.5);
  declare_const_real("onehalf",1.5);
  
#ifdef DOUBLE
  declare_real("y2");
#endif
  
  code( (prec==8) ? "common /finvsqrtdoubledata/ finvsqrtdoubletab(4096)" : "common /finvsqrtsingledata/ finvsqrtsingletab(4096)");
  
#ifdef DOUBLE
  assign("tt","X'5fe80000'");
  assign("tt","ishft(tt,32)");
#else
  assign("tt","X'1fc00000'");
#endif
  
  assign("fval","x");
  assign("tA","ishft(bval,-1)");
  assign("tB","tt - tA");
#ifdef DOUBLE
  assign("tC","ishft(tB,-40)");
#else
  assign("tC","ishft(tB,-11)");
#endif
  assign("tC","iand(tC,ishft(1,12)-1)");
  
  assign("tD", (prec==8) ? "finvsqrtdoubletab(tC+1)" : "finvsqrtsingletab(tC+1)");
  assign("f1","fval * half");
#ifdef DOUBLE
  assign("tD","ishft(tD,26)");
#endif
  assign("bval","tB - tD");
  assign("y","fval");
  
#ifdef DOUBLE
  assign("y2","y * (onehalf - f1 * y * y )");
  assign("invsqrt","y2 * (onehalf - f1 * y2 * y2 )");
#else
  assign("invsqrt","y * (onehalf - f1 * y * y )");
#endif
  code("return");
  code("end");
  flush_buffers();
  
  /* block to initialize the table data - write directly to file */
  fprintf(output,"%sblock data\n",indent());
  if(prec==8) 
    fprintf(output,"%scommon /finvsqrtdoubledata/ finvsqrtdoubletab(4096)\n",indent());
  else
    fprintf(output,"%scommon /finvsqrtsingledata/ finvsqrtsingletab(4096)\n",indent());
  
  fprintf(output,"%sinteger*4 I\n",indent());
  
  
  if(prec==8) 
    fprintf(output,"%sinteger*4 finvsqrtdoubletab\n"
	    "%sinclude 'finvsqrtdata_double.inc'\n",indent(),indent());
  else
    fprintf(output,"%sinteger*4 finvsqrtsingletab\n"
	    "%sinclude 'finvsqrtdata_single.inc'\n",indent(),indent());   
  
  fprintf(output,"%send\n\n\n",indent());
}




void fortran_vecinvsqrt()
{
  newline();
  comment("fortran vectorized invsqrt routine");

  strcat(header,"      subroutine vecinvsqrt(");
  declare_real_vector("in");
  declare_real_vector("out");
  declare_int("n");
  nargs=ndecl;

  declare_int4( (prec==8) ? "finvsqrtdoubletab" : "finvsqrtsingletab");

  declare_int("i");
  declare_intreal("bval");
  declare_real("fval");
  code("equivalence(bval,fval)");

  declare_intreal("tt");
  declare_intreal("tA");
  declare_intreal("tB");
  declare_intreal("tC");
  declare_intreal("tD");

  declare_real("f1");
  declare_real("y");
  declare_real("x");
  declare_const_real("half",0.5);
  declare_const_real("onehalf",1.5);
  
#ifdef DOUBLE
  declare_real("y2");
#else
  
#endif
  code( (prec==8) ? "common /finvsqrtdoubledata/ finvsqrtdoubletab(4096)" : "common /finvsqrtsingledata/ finvsqrtsingletab(4096)");
  
#ifdef DOUBLE
  assign("tt","X'5fe80000'");
  assign("tt","ishft(tt,32)");
#else
  assign("tt","X'1fc00000'");
#endif

  start_loop("i","1","n");
  
  assign("x","in(i)");
  assign("fval","x");
  assign("tA","ishft(bval,-1)");
  assign("tB","tt - tA");
#ifdef DOUBLE
  assign("tC","ishft(tB,-40)");
#else
  assign("tC","ishft(tB,-11)");
#endif
  assign("tC","iand(tC,ishft(1,12)-1)");
 
  assign("tD", (prec==8) ? "finvsqrtdoubletab(tC+1)" : "finvsqrtsingletab(tC+1)");
  assign("f1","fval * half");
#ifdef DOUBLE
  assign("tD","ishft(tD,26)");
#endif
  assign("bval","tB - tD");
  assign("y","fval");

#ifdef DOUBLE
  assign("y2","y * (onehalf - f1 * y * y )");
  assign("out(i)","y2 * (onehalf - f1 * y2 * y2 )");
#else
  assign("out(i)","y * (onehalf - f1 * y * y )");
#endif
  end_loop();
  code("return");
  code("end");
  flush_buffers();
}
