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

#include "mkinl.h"
#include <string.h>



/* perform the first stage table lookup step of the invsqrt */
void invsqrt_table_lookup(char *text)
{
  assign( bC ? "bitpattern%s.fval" : "fval%s","rsq%s",text,text);

  if(bC) {
    assign("iexp%s","EXP_ADDR(bitpattern%s.bval)",text,text);
    assign("addr%s","FRACT_ADDR(bitpattern%s.bval)",text,text);
    assign("result%s.bval","cinvsqrtexptab[iexp%s] | cinvsqrtfracttab[addr%s]",text,text,text);
  } else {
    assign("iexp%s","rshift(and(bval%s,expmask),expshift)",text,text);
    assign("addr%s","rshift(and(bval%s,or(fractmask,explsb)),fractshift)",text,text);
    assign("result%s","or(finvsqrtexptab(iexp%s+1),finvsqrtfracttab(addr%s+1))",text,text,text);
  }
}


int invsqrt_iterate(void)
{
  int nflop=0;
  int i,j;

  comment("N-R iteration to improve accuracy");
  newline();

  if(bC)
    for(j=1;j<=loop.nj;j++)
      for(i=1;i<=loop.ni;i++) 
	assign("lu%d", "result%d.fval",i*10+j,i*10+j);
	
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("rsqlu%d","rsq%d","lu%d",i*10+j,i*10+j,i*10+j);
  
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("tmp%d","rsqlu%d","lu%d",i*10+j,i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++)
      assign("rinv%d","half*lu%d*(three-tmp%d)",i*10+j,i*10+j,i*10+j);
  
#ifdef DOUBLE
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("rsqlu%d","rsq%d","rinv%d",i*10+j,i*10+j,i*10+j);

  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("tmp%d","rsqlu%d","rinv%d",i*10+j,i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++)
      assign("rinv%d","half*rinv%d*(three-tmp%d)",i*10+j,i*10+j,i*10+j);  
#endif
  return nflop;
}


void whole_invsqrt(char *name)
{
  invsqrt_table_lookup(name);
  if(bC)
    assign("lu%s", "result%s.fval",name,name);
  assign("rinv%s","(half*lu%s*(three-((rsq%s*lu%s)*lu%s)))",
	 name,name,name,name,name);
#ifdef DOUBLE
  assign("rinv%s","(half*rinv%s*(three-((rsq%s*rinv%s)*rinv%s)))",
	 name,name,name,name,name); 
#endif
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
    nflop+=invsqrt_iterate();  
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
    if(opt.inline_invsqrt && arch.gmx_invsqrt)
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

  if(!bC) {
    declare_int4("finvsqrtexptab");
    declare_int4("finvsqrtfracttab");
  }

  declare_const_real("half",0.5);
  declare_const_int("nexp",256);
  declare_const_int("nfract",4096);
  declare_const_int("maxfract",8388607);
  declare_const_int("fractshift",12);
  declare_const_int("fractmask",8388607);
  declare_const_int("fractf",1065353216);
  declare_const_int("expshift",23);
  declare_const_int("expmask",2139095040);
  declare_const_int("explsb",8388608);
  
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) {
      ij=10*i+j;
      sprintf(buf,"lu%d",ij);
      declare_real4(buf);
      sprintf(buf,"tmp%d",ij);
      declare_real(buf);
      sprintf(buf,"rsqlu%d",ij);
      declare_real(buf);
      sprintf(buf,"iexp%d",ij);
      declare_int4(buf);
      sprintf(buf,"addr%d",ij);
      declare_int4(buf);

      if(bC) {
	sprintf(buf,"bitpattern%d",ij);
	declare_other("t_convert",buf);
	sprintf(buf,"result%d",ij);
	declare_other("t_convert",buf);
      } else {
	sprintf(buf,"bval%d",ij);
	declare_int4(buf);
	sprintf(buf,"fval%d",ij);
	declare_real4(buf);
	sprintf(buf,"result%d",ij);
	declare_int4(buf);
	sprintf(buf,"equivalence(bval%d,fval%d)",ij,ij);
	code(buf);
	sprintf(buf,"equivalence(result%d,lu%d)",ij,ij);
	code(buf);
      }
    }

  if(!bC)
      code("common /finvsqrtdata/ finvsqrtexptab(256),finvsqrtfracttab(4096)");

}


void fortran_invsqrt()
{
  /* First the nonvectorized version */
#ifdef SOFTWARE_INVSQRT
  newline();
  comment("fortran invsqrt routine");
  strcat(header,"      function invsqrt(");
  declare_real4("x");
  nargs=ndecl; 
  declare_real("invsqrt");
  
  declare_int4("finvsqrtexptab"); 
  declare_int4("finvsqrtfracttab"); 

  declare_int4("bval");
  declare_real4("fval");
  declare_int4("result");
  declare_real4("lu");
  declare_int("iexp");
  declare_int("addr");
    
  code("equivalence(bval,fval)");
  code("equivalence(result,lu)");
  
#ifdef DOUBLE
  declare_real("y");
#endif
  declare_const_real("half",0.5);
  declare_const_real("three",3.0);
  declare_const_int("nexp",256);
  declare_const_int("nfract",4096);
  declare_const_int("maxfract",8388607);
  declare_const_int("fractshift",12);
  declare_const_int("fractmask",8388607);
  declare_const_int("fractf",1065353216);
  declare_const_int("expshift",23);
  declare_const_int("expmask",2139095040);
  declare_const_int("explsb",8388608);
  
  
  code("common /finvsqrtdata/ finvsqrtexptab(nexp),finvsqrtfracttab(nfract)");
  
  assign("fval","x");
  assign("iexp","rshift(and(bval,expmask),expshift)");

  assign("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
  assign("result","or(finvsqrtexptab(iexp+1),finvsqrtfracttab(addr+1))");
#ifdef DOUBLE
  assign("y","(half*lu*(three-((x*lu)*lu)))");
  assign("invsqrt","(half*y*(three-((x*y)*y)))");
#else
  assign("invsqrt","(half*lu*(three-((x*lu)*lu)))");
#endif
  
  code("return");
  code("end");
  flush_buffers();

#endif
  
  /* And now the vectorized version */
  newline();
  comment("fortran vectorized invsqrt routine");
  strcat(header,"      subroutine f77vecinvsqrt(");
  declare_real_vector("indata");
  declare_real_vector("utdata");
  declare_int("n");
  nargs=ndecl;
  
  declare_int("i");  
#ifdef SOFTWARE_SQRT
  declare_int4("finvsqrtexptab"); 
  declare_int4("finvsqrtfracttab"); 

  declare_real4("x");
  declare_int4("bval");
  declare_real4("fval");
  declare_int4("result");
  declare_real4("lu");
  declare_int("iexp");
  declare_int("addr");

  
  code("equivalence(bval,fval)");
  code("equivalence(result,lu)");
  
#ifdef DOUBLE
  declare_real("y");
#endif
  declare_const_real("half",0.5);
  declare_const_real("three",3.0);
  declare_const_int("nexp",256);
  declare_const_int("nfract",4096);
  declare_const_int("maxfract",8388607);
  declare_const_int("fractshift",12);
  declare_const_int("fractmask",8388607);
  declare_const_int("fractf",1065353216);
  declare_const_int("expshift",23);
  declare_const_int("expmask",2139095040);
  declare_const_int("explsb",8388608);  
  
  code("common /finvsqrtdata/ finvsqrtexptab(nexp),finvsqrtfracttab(nfract)");

  start_loop("i","1","n");
  assign("x","indata(i)");
  assign("fval","x");
  assign("iexp","rshift(and(bval,expmask),expshift)");

  assign("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
  assign("result","or(finvsqrtexptab(iexp+1),finvsqrtfracttab(addr+1))");
#ifdef DOUBLE
  assign("y","(half*lu*(three-((x*lu)*lu)))");
  assign("utdata(i)","(half*y*(three-((x*y)*y)))");
#else
  assign("utdata(i)","(half*lu*(three-((x*lu)*lu)))");
#endif
  end_loop();

#else /* no SOFTWARE_SQRT */
  start_loop("i","1","n");
  assign("utdata(i)","1.0/sqrt(indata(i))");
  end_loop();
#endif
  code("end");
  flush_buffers();
}
