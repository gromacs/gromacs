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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "mkinl.h"
#include <string.h>


/* perform the first stage table lookup step of the recip */
void recip_table_lookup(char *text)
{
  assign( bC ? "bitpattern%s.fval" : "fval%s","rsq%s",text,text);

  if(bC) {
    assign("iexp%s","EXP_ADDR(bitpattern%s.bval)",text,text);
    assign("addr%s","FRACT_ADDR(bitpattern%s.bval)",text,text);
    assign("result%s.bval","crecipexptab[iexp%s] | crecipfracttab[addr%s]",text,text,text);
  } else {
    assign("iexp%s","rshift(and(bval%s,expmask),expshift)",text,text);
    assign("addr%s","rshift(and(bval%s,or(fractmask,explsb)),fractshift)",text,text);
    assign("result%s","or(frecipexptab(iexp%s+1),frecipfracttab(addr%s+1))",text,text,text);
  }
}


int recip_iterate(void)
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
      subtract("tmp%d","two","rsqlu%d",i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++)
      multiply("rinvsq%d","lu%d","tmp%d",i*10+j,i*10+j,i*10+j);
  
#ifdef DOUBLE
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      multiply("rsqlu%d","rsq%d","rinvsq%d",i*10+j,i*10+j,i*10+j);

  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++) 
      subtract("tmp%d","two","rsqlu%d",i*10+j,i*10+j);
      
  for(j=1;j<=loop.nj;j++)
    for(i=1;i<=loop.ni;i++)
      multiply("rinvsq%d","rinvsq%d","tmp%d",i*10+j,i*10+j,i*10+j);  
#endif
  return nflop;
}


void whole_recip(char *name)
{
  recip_table_lookup(name);
  if(bC)
    assign("lu%s", "result%s.fval",name,name);
  assign("rinvsq%s","lu%s*(two-rsq%s*lu%s)",name,name,name,name);

#ifdef DOUBLE
  assign("rinvsq%s","rinvsq%s*(two-rsq%s*rinvsq%s)",name,name,name,name); 
#endif
}


int inlined_recip(void)
{
  int nflop=0;
  int i,j,k;
  char buf[100];
  comment("Recip table lookup");

  if(opt.decrease_lookup_latency) {
    for(j=1;j<=loop.nj;j++)
      for(i=1;i<=loop.ni;i++) {
	sprintf(buf,"%d",i*10+j);
	recip_table_lookup(buf);
      }
    nflop+=recip_iterate();  
  } else {
    for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++) {
	  sprintf(buf,"%d",i*10+j);
	  whole_recip(buf);
	}
  }

  return nflop;
}




int calc_recip()
{
    int nflop=0;
    int i,j,tmpflop;
    /* check if it should be inlined */
    if(opt.inline_gmxcode && arch.gmx_recip)
      inlined_recip();
    else {    /* otherwise do the ordinary one */
      for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++)
	  assign("rinvsq%d", arch.gmx_recip ? 
		 "recip(rsq%d)" : "1.0/(rsq%d)",i*10+j,i*10+j);
    }
#ifdef DOUBLE
    tmpflop=6;
#else
    tmpflop=3;
#endif
    
    if(!arch.gmx_recip)
      tmpflop=1;
    if(loop.sol==SOL_NO)
      nflop=tmpflop;
    else if(loop.sol!=SOL_WATERWATER)
      nflop=tmpflop*3;
    else
      nflop=tmpflop*9;
    return nflop;
}



void recip_vars()
{
  int i,j,ij;
  char buf[100];

  if(!bC) {
    declare_int4("frecipexptab");
    declare_int4("frecipfracttab");
  }

  if(!DO_INLINE_INVSQRT) {
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
	declare_real(buf);
	sprintf(buf,"tmp%d",ij);
	declare_real(buf);
	sprintf(buf,"rsqlu%d",ij);
	declare_real(buf);
	sprintf(buf,"iexp%d",ij);
	declare_intreal(buf);
	sprintf(buf,"addr%d",ij);
	declare_intreal(buf);
	
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
  }

  if(!bC)
      code("common /frecipdata/ frecipexptab(256),frecipfracttab(4096)");

}

void fortran_recip()
{
  /* First the nonvectorized version */
#ifdef SOFTWARE_RECIP
  newline();
  comment("fortran recip routine");
  strcat(header,"      function recip(");
  declare_real4("x");
  nargs=ndecl; 
  declare_real("recip");
  
  declare_int4("frecipexptab"); 
  declare_int4("frecipfracttab"); 

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
  declare_const_real("two",2.0);
  declare_const_int("nexp",256);
  declare_const_int("nfract",4096);
  declare_const_int("maxfract",8388607);
  declare_const_int("fractshift",12);
  declare_const_int("fractmask",8388607);
  declare_const_int("fractf",1065353216);
  declare_const_int("expshift",23);
  declare_const_int("expmask",2139095040);
  declare_const_int("explsb",8388608);
  
  
  code("common /frecipdata/ frecipexptab(nexp),frecipfracttab(nfract)");
  
  assign("fval","x");
  assign("iexp","rshift(and(bval,expmask),expshift)");

  assign("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
  assign("result","or(frecipexptab(iexp+1),frecipfracttab(addr+1))");
#ifdef DOUBLE
  assign("y","lu*(two-x*lu)");
  assign("recip","y*(two-x*y)");
#else
  assign("recip","lu*(two-x*lu)");
#endif
  
  code("return");
  code("end");
  flush_buffers();
#endif
  
  /* And now the vectorized version */
  newline();
  comment("fortran vectorized recip routine");
  strcat(header,"      subroutine f77vecrecip(");
  declare_real_vector("indata");
  declare_real_vector("utdata");
  declare_int("n");
  
  nargs=ndecl; 
  declare_int("i");
#ifdef SOFTWARE_RECIP
  declare_real4("x");
  
  declare_int4("frecipexptab"); 
  declare_int4("frecipfracttab"); 

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
  declare_const_real("two",2.0);
  declare_const_int("nexp",256);
  declare_const_int("nfract",4096);
  declare_const_int("maxfract",8388607);
  declare_const_int("fractshift",12);
  declare_const_int("fractmask",8388607);
  declare_const_int("fractf",1065353216);
  declare_const_int("expshift",23);
  declare_const_int("expmask",2139095040);
  declare_const_int("explsb",8388608);
  
  
  code("common /frecipdata/ frecipexptab(nexp),frecipfracttab(nfract)");
  
  start_loop("i","1","n");
  assign("x","indata(i)");
  assign("fval","x");
  assign("iexp","rshift(and(bval,expmask),expshift)");

  assign("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
  assign("result","or(frecipexptab(iexp+1),frecipfracttab(addr+1))");
#ifdef DOUBLE
  assign("y","lu*(two-x*lu)");
  assign("utdata(i)","y*(two-x*y)");
#else
  assign("utdata(i)","lu*(two-x*lu)");
#endif
  end_loop();
  

#else /* no SOFTWARE_RECIP */
  start_loop("i","1","n");
  assign("utdata(i)","1.0/indata(i)");
  end_loop();
#endif
  code("end");
  flush_buffers();
}




