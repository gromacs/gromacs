/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_inner_altivec_c = "$Id$";

#include <ppc_altivec.h>

#include<stdio.h>
void inl0100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float fs,nul;
  vector float dx,dy,dz;
  vector float vnbtot,c6,c12;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvsq,rsq,rinvsix,vnb6,vnb12;

  int n,k,k0,ii,is3,ii3,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int ntiA,tja,tjb,tjc,tjd;

  nul=vec_zero();
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void inl0300_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float fs,nul,tsc;
  vector float dx,dy,dz;
  vector float vnbtot,c6,c12;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,r,rsq;
  vector float VVd,FFd,VVr,FFr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_4_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
      fs              = vec_nmsub(c6,FFd,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_madd(vec_madd(fs,tsc,nul),rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_2_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
      fs              = vec_nmsub(c6,FFd,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_madd(vec_madd(fs,tsc,nul),rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      do_1_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
      fs              = vec_nmsub(c6,FFd,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_madd(vec_madd(fs,tsc,nul),rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void inl1000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vcoul,fs,nul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,rinvsq,rsq;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      vcoul           = vec_madd(qq,rinv,nul);
      fs              = vec_madd(vcoul,rinvsq,nul);
      vctot           = vec_add(vctot,vcoul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      vcoul           = vec_madd(qq,rinv,nul);
      fs              = vec_madd(vcoul,rinvsq,nul);
      vctot           = vec_add(vctot,vcoul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);     
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      vcoul           = vec_madd(qq,rinv,nul);
      fs              = vec_madd(vcoul,rinvsq,nul);
      vctot           = vec_add(vctot,vcoul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl1100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vcoul,fs,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,rinvsq,rsq,rinvsix,vnb6,vnb12;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoul           = vec_madd(qq,rinv,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_add(vctot,vcoul);
      fs              = vec_madd(vec_twelve(),vnb12,vcoul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoul           = vec_madd(qq,rinv,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_add(vctot,vcoul);
      fs              = vec_madd(vec_twelve(),vnb12,vcoul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoul           = vec_madd(qq,rinv,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_add(vctot,vcoul);
      fs              = vec_madd(vec_twelve(),vnb12,vcoul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    

    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}




void inl2000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vkrf,vcrf,krsq,vcoul,fs,nul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,rinvsq,rsq;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);

      fs              = vec_nmsub(vec_two(),krsq,rinv);
      vctot           = vec_madd(qq,vcoul,vctot);
      fs              = vec_madd(fs,qq,nul);
      fs              = vec_madd(fs,rinvsq,nul);

      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      fs              = vec_nmsub(vec_two(),krsq,rinv);
      vctot           = vec_madd(qq,vcoul,vctot);
      fs              = vec_madd(fs,qq,nul);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);     
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      fs              = vec_nmsub(vec_two(),krsq,rinv);
      vctot           = vec_madd(qq,vcoul,vctot);
      fs              = vec_madd(fs,qq,nul);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl2100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vkrf,vcrf,krsq,vcoul,fs,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,rinvsq,rsq,rinvsix,vnb6,vnb12;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_nmsub(vec_two(),krsq,rinv);  /* rinv-2*krsq */
      fs              = vec_madd(qq,fs,nul);          /* qq*(rinv-2*krsq) */
      fs              = vec_madd(vec_twelve(),vnb12,fs); 
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_nmsub(vec_two(),krsq,rinv);  /* rinv-2*krsq */
      fs              = vec_madd(qq,fs,nul);          /* qq*(rinv-2*krsq) */
      fs              = vec_madd(vec_twelve(),vnb12,fs); 
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fs              = vec_nmsub(vec_two(),krsq,rinv);  /* rinv-2*krsq */
      fs              = vec_madd(qq,fs,nul);          /* qq*(rinv-2*krsq) */
      fs              = vec_madd(vec_twelve(),vnb12,fs); 
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinvsq,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    

    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}






void inl3000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,tsc,vcoul,fs,nul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,r,rsq,VVc,FFc;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      do_4_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs              = vec_nmsub(qq,FFc,nul);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      do_2_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs              = vec_nmsub(qq,FFc,nul);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);     
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      do_1_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs              = vec_nmsub(qq,FFc,nul);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl3100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vcoul,tsc,fs,fs2,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12,VVc,FFc;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,r,rinvsq,rsq,rinvsix,vnb6,vnb12;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_4_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs2             = vec_madd(qq,FFc,nul);   /* fijC */
      vctot           = vec_madd(qq,VVc,vctot);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinv,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fs              = vec_nmsub(fs2,tsc,fs);
      fs              = vec_madd(fs,rinv,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_2_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs2             = vec_madd(qq,FFc,nul);   /* fijC */
      vctot           = vec_madd(qq,VVc,vctot);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      vnbtot          = vec_add(vnbtot,vnb12);
      fs              = vec_madd(fs,rinv,nul);
      fs              = vec_nmsub(fs2,tsc,fs);
      fs              = vec_madd(fs,rinv,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_1_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc);
      fs2             = vec_madd(qq,FFc,nul);   /* fijC */
      vctot           = vec_madd(qq,VVc,vctot);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fs              = vec_madd(vec_twelve(),vnb12,nul);
      fs              = vec_nmsub(vec_six(),vnb6,fs);
      fs              = vec_madd(fs,rinv,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fs              = vec_nmsub(fs2,tsc,fs);
      fs              = vec_madd(fs,rinv,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    

    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void inl3300_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float fs,nul,tsc;
  vector float dx,dy,dz,vfacel,vcoul,vctot;
  vector float vnbtot,c6,c12,iq,qq;
  vector float fix,fiy,fiz;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinv,r,rsq;
  vector float VVc,FFc,VVd,FFd,VVr,FFr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    fix        = nul;
    fiy        = nul;
    fiz        = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_4_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc,&VVd,&FFd,&VVr,&FFr);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_nmsub(qq,FFc,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c6,FFd,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_2_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc,&VVd,&FFd,&VVr,&FFr);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_nmsub(qq,FFc,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c6,FFd,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      do_1_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&FFc,&VVd,&FFd,&VVr,&FFr);
      vctot           = vec_madd(qq,VVc,vctot);
      fs              = vec_nmsub(qq,FFc,nul);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fs              = vec_nmsub(c6,FFd,fs);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fs              = vec_nmsub(c12,FFr,fs);
      fs              = vec_madd(fs,tsc,nul);
      fs              = vec_madd(fs,rinv,nul);
      fix             = vec_madd(fs,dx,fix); /* +=fx */
      fiy             = vec_madd(fs,dy,fiy); /* +=fy */
      fiz             = vec_madd(fs,dz,fiz); /* +=fz */
      dx              = vec_nmsub(dx,fs,nul); /* -fx */
      dy              = vec_nmsub(dy,fs,nul); /* -fy */
      dz              = vec_nmsub(dz,fs,nul); /* -fz */
      transpose_3_to_1(dx,dy,dz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
    tmp1 = vec_add(tmp1,tmp3);
    tmp2 = vec_add(tmp2,tmp4);
    tmp1 = vec_add(tmp1,tmp2);    
    add_xyz_to_mem(faction+ii3,tmp1);
    add_xyz_to_mem(fshift+is3,tmp1);

    add_vector_to_float(Vnb+gid[n],vnbtot);
    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void inl1020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float fsO,fsH1,fsH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rinvsqH1,rinvsqH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      fsO             = vec_madd(vcoulO,rinvsqO,nul);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);

      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
 
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      fsO             = vec_madd(vcoulO,rinvsqO,nul);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);
      
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      fsO             = vec_madd(vcoulO,rinvsqO,nul);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      fsO             = vec_madd(vcoulO,rinvsqO,nul);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void inl1120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float vnbtot,c6,c12,rinvsix,vnb6,vnb12;
  vector float fsO,fsH1,fsH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rinvsqH1,rinvsqH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(vec_twelve(),vnb12,vcoulO);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      /* load 3 j charges and multiply by iq */
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(vec_twelve(),vnb12,vcoulO);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(vec_twelve(),vnb12,vcoulO);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vcoulO          = vec_madd(qqO,rinvO,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vcoulH1         = vec_madd(qqH,rinvH1,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(vec_twelve(),vnb12,vcoulO);
      vcoulH2         = vec_madd(qqH,rinvH2,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_madd(vcoulH1,rinvsqH1,nul);
      fsH2            = vec_madd(vcoulH2,rinvsqH2,nul);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      vctot           = vec_add(vctot,vcoulO);
      vcoulH1         = vec_add(vcoulH1,vcoulH2);
      vctot           = vec_add(vctot,vcoulH1);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void inl2020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float vkrf,vcrf;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float fsO,fsH1,fsH2,krsqO,krsqH1,krsqH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rinvsqH1,rinvsqH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);

  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,qqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 3 j charges and multiply by iq */
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,qqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,qqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,qqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl2120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float vkrf,vcrf;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float vnbtot,c6,c12,rinvsix,vnb6,vnb12;
  vector float fsO,fsH1,fsH2,krsqO,krsqH1,krsqH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rinvsqH1,rinvsqH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(qqO,fsO,nul);          
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(vec_twelve(),vnb12,fsO); 
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(qqO,fsO,nul);          
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(vec_twelve(),vnb12,fsO); 
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(qqO,fsO,nul);          
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(vec_twelve(),vnb12,fsO); 
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsqH1        = vec_madd(rinvH1,rinvH1,nul);
      rinvsqH2        = vec_madd(rinvH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fsO             = vec_nmsub(vec_two(),krsqO,rinvO);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_madd(qqO,fsO,nul);          
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(vec_twelve(),vnb12,fsO); 
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      fsH1            = vec_nmsub(vec_two(),krsqH1,rinvH1);
      fsH2            = vec_nmsub(vec_two(),krsqH2,rinvH2);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      fsO             = vec_madd(fsO,rinvsqO,nul);
      fsH1            = vec_madd(fsH1,qqH,nul);
      fsH2            = vec_madd(fsH2,qqH,nul);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
      fsH1            = vec_madd(fsH1,rinvsqH1,nul);
      fsH2            = vec_madd(fsH2,rinvsqH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void inl3020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float fsO,fsH1,fsH2,tsc,VVcO,FFcO,VVcH1,FFcH1,VVcH2,FFcH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rO,rH1,rH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
    
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      do_4_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_nmsub(qqO,FFcO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsO             = vec_madd(fsO,tsc,nul);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);
      
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);

      /* load 3 j charges and multiply by iq */
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_3_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_nmsub(qqO,FFcO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsO             = vec_madd(fsO,tsc,nul);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);
      
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
    
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      do_2_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_nmsub(qqO,FFcO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsO             = vec_madd(fsO,tsc,nul);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);
 
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);

      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      do_1_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_nmsub(qqO,FFcO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsO             = vec_madd(fsO,tsc,nul);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);
      
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl3120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float vnbtot,c6,c12,rinvsix,rinvsqO,vnb6,vnb12;
  vector float fsO,fsH1,fsH2,tsc,VVcO,FFcO,VVcH1,FFcH1,VVcH2,FFcH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rO,rH1,rH2,rsqO,rsqH1,rsqH2;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_4_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      vnb6            = vec_madd(c6,rinvsix,nul);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_madd(vec_twelve(),vnb12,nul);
      tmp1            = vec_madd(qqO,FFcO,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      fsO             = vec_nmsub(tmp1,tsc,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      /* load 3 j charges and multiply by iq */
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_3_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      vnb6            = vec_madd(c6,rinvsix,nul);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_madd(vec_twelve(),vnb12,nul);
      tmp1            = vec_madd(qqO,FFcO,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      fsO             = vec_nmsub(tmp1,tsc,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_2_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      vnb6            = vec_madd(c6,rinvsix,nul);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_madd(vec_twelve(),vnb12,nul);
      tmp1            = vec_madd(qqO,FFcO,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      fsO             = vec_nmsub(tmp1,tsc,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_1_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO);
      do_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      vnb6            = vec_madd(c6,rinvsix,nul);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      fsO             = vec_madd(vec_twelve(),vnb12,nul);
      tmp1            = vec_madd(qqO,FFcO,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      fsO             = vec_nmsub(vec_six(),vnb6,fsO);
      vnbtot          = vec_sub(vnbtot,vnb6);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_nmsub(qqH,FFcH1,nul);
      fsH2            = vec_nmsub(qqH,FFcH2,nul);
      fsO             = vec_nmsub(tmp1,tsc,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      fsH1            = vec_madd(fsH1,tsc,nul);
      fsH2            = vec_madd(fsH2,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void inl3320_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float tsc;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,vcoulO,vcoulH1,vcoulH2,nul;
  vector float vnbtot,c6,c12;
  vector float fsO,fsH1,fsH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float rinvO,rinvH1,rinvH2,rsqO,rsqH1,rsqH2;
  vector float rO,rH1,rH2,VVcO,FFcO,VVcH1,FFcH1,VVcH2,FFcH2,VVd,FFd,VVr,FFr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  vfacel=load_float_and_splat(&facel);

  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    fiOx       = nul;
    fiOy       = nul;
    fiOz       = nul;
    fiH1x      = nul;
    fiH1y      = nul;
    fiH1z      = nul;
    fiH2x      = nul;
    fiH2y      = nul;
    fiH2z      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_4_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO,&VVd,&FFd,&VVr,&FFr);
      do_4_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_4_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      fsO             = vec_madd(qqO,FFcO,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fsO             = vec_madd(c6,FFd,fsO);
      fsH1            = vec_madd(qqH,FFcH1,nul);
      fsH2            = vec_madd(qqH,FFcH2,nul);
      fsO             = vec_madd(c12,FFr,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fsO             = vec_nmsub(fsO,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsH1            = vec_nmsub(fsH1,tsc,nul);
      fsH2            = vec_nmsub(fsH2,tsc,nul);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
      add_xyz_to_mem(faction+j3d,tmp4);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
  
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_3_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO,&VVd,&FFd,&VVr,&FFr);
      do_3_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_3_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      fsO             = vec_madd(qqO,FFcO,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fsO             = vec_madd(c6,FFd,fsO);
      fsH1            = vec_madd(qqH,FFcH1,nul);
      fsH2            = vec_madd(qqH,FFcH2,nul);
      fsO             = vec_madd(c12,FFr,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fsO             = vec_nmsub(fsO,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsH1            = vec_nmsub(fsH1,tsc,nul);
      fsH2            = vec_nmsub(fsH2,tsc,nul);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
      add_xyz_to_mem(faction+j3c,tmp3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_2_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO,&VVd,&FFd,&VVr,&FFr);
      do_2_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_2_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      fsO             = vec_madd(qqO,FFcO,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fsO             = vec_madd(c6,FFd,fsO);
      fsH1            = vec_madd(qqH,FFcH1,nul);
      fsH2            = vec_madd(qqH,FFcH2,nul);
      fsO             = vec_madd(c12,FFr,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fsO             = vec_nmsub(fsO,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsH1            = vec_nmsub(fsH1,tsc,nul);
      fsH2            = vec_nmsub(fsH2,tsc,nul);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);

      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
      add_xyz_to_mem(faction+j3a,tmp1);
      add_xyz_to_mem(faction+j3b,tmp2);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_1_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&FFcO,&VVd,&FFd,&VVr,&FFr);
      do_1_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
      do_1_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      fsO             = vec_madd(qqO,FFcO,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      fsO             = vec_madd(c6,FFd,fsO);
      fsH1            = vec_madd(qqH,FFcH1,nul);
      fsH2            = vec_madd(qqH,FFcH2,nul);
      fsO             = vec_madd(c12,FFr,fsO);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      fsO             = vec_nmsub(fsO,tsc,nul);
      vctot           = vec_madd(qqH,VVcH2,vctot);
      fsH1            = vec_nmsub(fsH1,tsc,nul);
      fsH2            = vec_nmsub(fsH2,tsc,nul);
      fsO             = vec_madd(fsO,rinvO,nul);
      fsH1            = vec_madd(fsH1,rinvH1,nul);
      fsH2            = vec_madd(fsH2,rinvH2,nul);
 
      fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
      dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
      fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
      dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
      fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
      dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
      fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
      dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
      fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
      dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
      fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
      dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
      fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
      dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
      fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
      dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
      fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
      dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */

      transpose_3_to_1(dOx,dOy,dOz,&tmp1);
      add_xyz_to_mem(faction+j3a,tmp1);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}




void inl1030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11,rinvsq12,rinvsq13;
  vector float rinvsq21,rinvsq22,rinvsq23;
  vector float rinvsq31,rinvsq32,rinvsq33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,qqOOt,qqOHt,qqHHt;

 

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_madd(rinv11,qqOO,nul);
      vc12            = vec_madd(rinv12,qqOH,nul);
      vc13            = vec_madd(rinv13,qqOH,nul);
      vc21            = vec_madd(rinv21,qqOH,nul);
      vc22            = vec_madd(rinv22,qqHH,nul);
      vc23            = vec_madd(rinv23,qqHH,nul);
      vc31            = vec_madd(rinv31,qqOH,nul);
      vc32            = vec_madd(rinv32,qqHH,nul);
      vc33            = vec_madd(rinv33,qqHH,nul);

      fs11            = vec_madd(vc11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);
      
      fs11            = vec_madd(vc11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);
      
      fs11            = vec_madd(vc11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);
      
      fs11            = vec_madd(vc11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void inl1130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11,rinvsq12,rinvsq13;
  vector float rinvsq21,rinvsq22,rinvsq23;
  vector float rinvsq31,rinvsq32,rinvsq33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnb6,vnb12,vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) { 
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);
      
      vnb6            = vec_madd(c6,rinvsix,nul);
      rinvsix         = vec_madd(rinvsix,rinvsix,nul);
      vc11            = vec_madd(rinv11,qqOO,nul);
      vc12            = vec_madd(rinv12,qqOH,nul);
      vc13            = vec_madd(rinv13,qqOH,nul);
      vnb12           = vec_madd(c12,rinvsix,nul);
      vc21            = vec_madd(rinv21,qqOH,nul);
      vc22            = vec_madd(rinv22,qqHH,nul);
      vc23            = vec_madd(rinv23,qqHH,nul);
      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vc31            = vec_madd(rinv31,qqOH,nul);
      vc32            = vec_madd(rinv32,qqHH,nul);
      vc33            = vec_madd(rinv33,qqHH,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs11            = vec_add(fs11,vc11);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      vctot           = vec_add(vctot,vc32);

      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);
      vctot           = vec_add(vctot,vc21);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);
      
      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);
      
      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);

      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs11            = vec_add(fs11,vc11);
      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);
      
      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);

      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs11            = vec_add(fs11,vc11);
      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);
      
      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vc11            = vec_madd(rinv11,qqOOt,nul);
      vc12            = vec_madd(rinv12,qqOHt,nul);
      vc13            = vec_madd(rinv13,qqOHt,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vc21            = vec_madd(rinv21,qqOHt,nul);
      vc22            = vec_madd(rinv22,qqHHt,nul);
      vc23            = vec_madd(rinv23,qqHHt,nul);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vc31            = vec_madd(rinv31,qqOHt,nul);
      vc32            = vec_madd(rinv32,qqHHt,nul);
      vc33            = vec_madd(rinv33,qqHHt,nul);

      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs11            = vec_add(fs11,vc11);
      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(vc12,rinvsq12,nul);
      fs13            = vec_madd(vc13,rinvsq13,nul);
      fs21            = vec_madd(vc21,rinvsq21,nul);
      fs22            = vec_madd(vc22,rinvsq22,nul);
      fs23            = vec_madd(vc23,rinvsq23,nul);
      fs31            = vec_madd(vc31,rinvsq31,nul);
      fs32            = vec_madd(vc32,rinvsq32,nul);
      fs33            = vec_madd(vc33,rinvsq33,nul);

      vctot           = vec_add(vctot,vc11);
      vc12            = vec_add(vc12,vc13);
      vc21            = vec_add(vc21,vc22);
      vc23            = vec_add(vc23,vc31);
      vc32            = vec_add(vc32,vc33);
      vctot           = vec_add(vctot,vc12);
      vc21            = vec_add(vc21,vc23);
      vctot           = vec_add(vctot,vc32);
      vctot           = vec_add(vctot,vc21); 

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
     
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);
    
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
 
}




void inl2030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11,rinvsq12,rinvsq13;
  vector float rinvsq21,rinvsq22,rinvsq23;
  vector float rinvsq31,rinvsq32,rinvsq33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,vkrf,vcrf;
  vector float krsq11,krsq12,krsq13,krsq21,krsq22,krsq23,krsq31,krsq32,krsq33;
  vector float qqOOt,qqOHt,qqHHt;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(fs11,qqOO,nul);
      fs12            = vec_madd(fs12,qqOH,nul);
      fs13            = vec_madd(fs13,qqOH,nul);
      fs21            = vec_madd(fs21,qqOH,nul);
      fs22            = vec_madd(fs22,qqHH,nul);
      fs23            = vec_madd(fs23,qqHH,nul);
      fs31            = vec_madd(fs31,qqOH,nul);
      fs32            = vec_madd(fs32,qqHH,nul);
      fs33            = vec_madd(fs33,qqHH,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOO,vc11,vctot);
      vctot           = vec_madd(qqOH,vc12,vctot);
      vctot           = vec_madd(qqOH,vc13,vctot);
      vctot           = vec_madd(qqOH,vc21,vctot);
      vctot           = vec_madd(qqHH,vc22,vctot);
      vctot           = vec_madd(qqHH,vc23,vctot);
      vctot           = vec_madd(qqOH,vc31,vctot);
      vctot           = vec_madd(qqHH,vc32,vctot);
      vctot           = vec_madd(qqHH,vc33,vctot);

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(fs11,qqOOt,nul);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);

      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);
      
      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(fs11,qqOOt,nul);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);
      
      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(fs11,qqOOt,nul);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl2130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11,rinvsq12,rinvsq13;
  vector float rinvsq21,rinvsq22,rinvsq23;
  vector float rinvsq31,rinvsq32,rinvsq33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33,vkrf,vcrf;
  vector float krsq11,krsq12,krsq13,krsq21,krsq22,krsq23,krsq31,krsq32,krsq33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnb6,vnb12,vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      fs11            = vec_madd(qqOO,fs11,nul);
      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(vec_twelve(),vnb12,fs11);
      fs12            = vec_madd(fs12,qqOH,nul);
      fs13            = vec_madd(fs13,qqOH,nul);
      fs21            = vec_madd(fs21,qqOH,nul);
      fs22            = vec_madd(fs22,qqHH,nul);
      fs23            = vec_madd(fs23,qqHH,nul);
      fs31            = vec_madd(fs31,qqOH,nul);
      fs32            = vec_madd(fs32,qqHH,nul);
      fs33            = vec_madd(fs33,qqHH,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOO,vc11,vctot);
      vctot           = vec_madd(qqOH,vc12,vctot);
      vctot           = vec_madd(qqOH,vc13,vctot);
      vctot           = vec_madd(qqOH,vc21,vctot);
      vctot           = vec_madd(qqHH,vc22,vctot);
      vctot           = vec_madd(qqHH,vc23,vctot);
      vctot           = vec_madd(qqOH,vc31,vctot);
      vctot           = vec_madd(qqHH,vc32,vctot);
      vctot           = vec_madd(qqHH,vc33,vctot);

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      fs11            = vec_madd(qqOOt,fs11,nul);
      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(vec_twelve(),vnb12,fs11);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);
      
      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      fs11            = vec_madd(qqOOt,fs11,nul);
      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(vec_twelve(),vnb12,fs11);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      rinvsq12        = vec_madd(rinv12,rinv12,nul);
      rinvsq13        = vec_madd(rinv13,rinv13,nul);
      rinvsq21        = vec_madd(rinv21,rinv21,nul);
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsq22        = vec_madd(rinv22,rinv22,nul);
      rinvsq23        = vec_madd(rinv23,rinv23,nul);
      rinvsq31        = vec_madd(rinv31,rinv31,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      rinvsq32        = vec_madd(rinv32,rinv32,nul);
      rinvsq33        = vec_madd(rinv33,rinv33,nul);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);

      fs11            = vec_nmsub(vec_two(),krsq11,rinv11);
      vc11            = vec_add(rinv11,krsq11);
      vc12            = vec_add(rinv12,krsq12);
      vc13            = vec_add(rinv13,krsq13);
      vc21            = vec_add(rinv21,krsq21);
      vc22            = vec_add(rinv22,krsq22);
      vc23            = vec_add(rinv23,krsq23);
      vc31            = vec_add(rinv31,krsq31);
      vc32            = vec_add(rinv32,krsq32);
      vc33            = vec_add(rinv33,krsq33);

      fs11            = vec_madd(qqOOt,fs11,nul);
      vc11            = vec_sub(vc11,vcrf);
      vc12            = vec_sub(vc12,vcrf);
      vc13            = vec_sub(vc13,vcrf);
      vc21            = vec_sub(vc21,vcrf);
      vc22            = vec_sub(vc22,vcrf);
      vc23            = vec_sub(vc23,vcrf);
      vc31            = vec_sub(vc31,vcrf);
      vc32            = vec_sub(vc32,vcrf);
      vc33            = vec_sub(vc33,vcrf);

      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs12            = vec_nmsub(vec_two(),krsq12,rinv12);
      fs13            = vec_nmsub(vec_two(),krsq13,rinv13);
      fs21            = vec_nmsub(vec_two(),krsq21,rinv21);
      fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
      fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
      fs31            = vec_nmsub(vec_two(),krsq31,rinv31);
      fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
      fs33            = vec_nmsub(vec_two(),krsq33,rinv33);

      fs11            = vec_madd(vec_twelve(),vnb12,fs11);
      fs12            = vec_madd(fs12,qqOHt,nul);
      fs13            = vec_madd(fs13,qqOHt,nul);
      fs21            = vec_madd(fs21,qqOHt,nul);
      fs22            = vec_madd(fs22,qqHHt,nul);
      fs23            = vec_madd(fs23,qqHHt,nul);
      fs31            = vec_madd(fs31,qqOHt,nul);
      fs32            = vec_madd(fs32,qqHHt,nul);
      fs33            = vec_madd(fs33,qqHHt,nul);

      fs11            = vec_madd(fs11,rinvsq11,nul);
      fs12            = vec_madd(fs12,rinvsq12,nul);
      fs13            = vec_madd(fs13,rinvsq13,nul);
      fs21            = vec_madd(fs21,rinvsq21,nul);
      fs22            = vec_madd(fs22,rinvsq22,nul);
      fs23            = vec_madd(fs23,rinvsq23,nul);
      fs31            = vec_madd(fs31,rinvsq31,nul);
      fs32            = vec_madd(fs32,rinvsq32,nul);
      fs33            = vec_madd(fs33,rinvsq33,nul);
      
      vctot           = vec_madd(qqOOt,vc11,vctot);
      vctot           = vec_madd(qqOHt,vc12,vctot);
      vctot           = vec_madd(qqOHt,vc13,vctot);
      vctot           = vec_madd(qqOHt,vc21,vctot);
      vctot           = vec_madd(qqHHt,vc22,vctot);
      vctot           = vec_madd(qqHHt,vc23,vctot);
      vctot           = vec_madd(qqOHt,vc31,vctot);
      vctot           = vec_madd(qqHHt,vc32,vctot);
      vctot           = vec_madd(qqHHt,vc33,vctot);
      
      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void inl3030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33;
  
  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,tsc;
  vector float VV11c,FF11c,VV12c,FF12c,VV13c,FF13c;
  vector float VV21c,FF21c,VV22c,FF22c,VV23c,FF23c;
  vector float VV31c,FF31c,VV32c,FF32c,VV33c,FF33c;
  vector float qqOOt,qqOHt,qqHHt;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_4_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_4_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_4_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_4_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_4_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_4_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_4_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_4_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_4_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      fs11            = vec_nmsub(qqOO,FF11c,nul);
      fs12            = vec_nmsub(qqOH,FF12c,nul);
      fs13            = vec_nmsub(qqOH,FF13c,nul);
      fs21            = vec_nmsub(qqOH,FF21c,nul);
      fs22            = vec_nmsub(qqHH,FF22c,nul);
      fs23            = vec_nmsub(qqHH,FF23c,nul);
      fs31            = vec_nmsub(qqOH,FF31c,nul);
      fs32            = vec_nmsub(qqHH,FF32c,nul);
      fs33            = vec_nmsub(qqHH,FF33c,nul);
      
      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_3_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_3_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_3_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_3_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_3_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_3_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_3_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_3_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_3_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_2_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_2_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_2_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_2_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_2_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_2_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_2_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_2_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_2_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_1_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_1_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_1_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_1_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_1_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_1_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_1_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_1_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_1_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);

    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void inl3130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[]) 
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33,tsc,VVc,FFc;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33,fs11c;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnb6,vnb12,vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;
  vector float VV11c,FF11c,VV12c,FF12c,VV13c,FF13c;
  vector float VV21c,FF21c,VV22c,FF22c,VV23c,FF23c;
  vector float VV31c,FF31c,VV32c,FF32c,VV33c,FF33c;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
   
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_4_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_4_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_4_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_4_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_4_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_4_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_4_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_4_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_4_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnb6            = vec_madd(c6,rinvsix,nul);
      vnb12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),nul);
      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11c           = vec_nmsub(qqOO,FF11c,nul);
      fs12            = vec_nmsub(qqOH,FF12c,nul);
      fs13            = vec_nmsub(qqOH,FF13c,nul);
      fs21            = vec_nmsub(qqOH,FF21c,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs22            = vec_nmsub(qqHH,FF22c,nul);
      fs23            = vec_nmsub(qqHH,FF23c,nul);
      fs31            = vec_nmsub(qqOH,FF31c,nul);
      fs32            = vec_nmsub(qqHH,FF32c,nul);
      fs11            = vec_madd(fs11,rinv11,nul);
      fs33            = vec_nmsub(qqHH,FF33c,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);

      fs11            = vec_madd(fs11c,tsc,fs11);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_3_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_3_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_3_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_3_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_3_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_3_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_3_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_3_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_3_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11c           = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs11            = vec_madd(fs11,rinv11,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11c,tsc,fs11);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      
      
      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_2_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_2_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_2_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_2_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_2_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_2_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_2_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_2_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_2_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11c           = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs11            = vec_madd(fs11,rinv11,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11c,tsc,fs11);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_1_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c,&FF11c);
      do_1_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_1_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_1_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_1_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_1_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_1_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_1_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_1_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnb6            = vec_madd(c6t,rinvsix,nul);
      vnb12           = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),nul);
      fs11            = vec_madd(vec_twelve(),vnb12,nul);
      fs11c           = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(vec_six(),vnb6,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs11            = vec_madd(fs11,rinv11,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      vnbtot          = vec_add(vnbtot,vnb12);
      vnbtot          = vec_sub(vnbtot,vnb6);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11c,tsc,fs11);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);
    
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void inl3330_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     float fshift[],
	     int gid[],
	     float pos[],
	     float faction[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
  vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
  vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12;
  vector float vnb6,vnb12,vnbtot,tsc,qqOOt,qqOHt,qqHHt,c6t,c12t;
  vector float VV11c,FF11c,VV12c,FF12c,VV13c,FF13c;
  vector float VV21c,FF21c,VV22c,FF22c,VV23c,FF23c;
  vector float VV31c,FF31c,VV32c,FF32c,VV33c,FF33c;
  vector float VVd,FFd,VVr,FFr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    fix1       = nul;
    fiy1       = nul;
    fiz1       = nul;
    fix2       = nul;
    fiy2       = nul;
    fiz2       = nul;
    fix3       = nul;
    fiy3       = nul;
    fiz3       = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_4_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
      do_4_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_4_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_4_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_4_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_4_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_4_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_4_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_4_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);

      fs11            = vec_nmsub(qqOO,FF11c,nul);
      fs12            = vec_nmsub(qqOH,FF12c,nul);
      fs13            = vec_nmsub(qqOH,FF13c,nul);
      fs21            = vec_nmsub(qqOH,FF21c,nul);
      fs11            = vec_nmsub(c6,FFd,fs11);
      fs22            = vec_nmsub(qqHH,FF22c,nul);
      fs23            = vec_nmsub(qqHH,FF23c,nul);
      fs31            = vec_nmsub(qqOH,FF31c,nul);
      fs32            = vec_nmsub(qqHH,FF32c,nul);
      fs33            = vec_nmsub(qqHH,FF33c,nul);
      fs11            = vec_nmsub(c12,FFr,fs11);

      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
      add_force_to_4_water(faction+j3a,faction+j3b,faction+j3c,faction+j3d,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_3_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
      do_3_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_3_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_3_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_3_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_3_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_3_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_3_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_3_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(c6t,FFd,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      fs11            = vec_nmsub(c12t,FFr,fs11);

      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);            

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);

      add_force_to_3_water(faction+j3a,faction+j3b,faction+j3c,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_2_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
      do_2_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_2_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_2_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_2_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_2_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_2_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_2_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_2_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(c6t,FFd,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      fs11            = vec_nmsub(c12t,FFr,fs11);

      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
      add_force_to_2_water(faction+j3a,faction+j3b,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_1_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
      do_1_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
      do_1_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
      do_1_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
      do_1_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
      do_1_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
      do_1_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
      do_1_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
      do_1_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);

      fs11            = vec_nmsub(qqOOt,FF11c,nul);
      fs12            = vec_nmsub(qqOHt,FF12c,nul);
      fs13            = vec_nmsub(qqOHt,FF13c,nul);
      fs21            = vec_nmsub(qqOHt,FF21c,nul);
      fs11            = vec_nmsub(c6t,FFd,fs11);
      fs22            = vec_nmsub(qqHHt,FF22c,nul);
      fs23            = vec_nmsub(qqHHt,FF23c,nul);
      fs31            = vec_nmsub(qqOHt,FF31c,nul);
      fs32            = vec_nmsub(qqHHt,FF32c,nul);
      fs33            = vec_nmsub(qqHHt,FF33c,nul);
      fs11            = vec_nmsub(c12t,FFr,fs11);

      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);

      fs11            = vec_madd(fs11,tsc,nul);
      fs12            = vec_madd(fs12,tsc,nul);
      fs13            = vec_madd(fs13,tsc,nul);
      fs21            = vec_madd(fs21,tsc,nul);
      fs22            = vec_madd(fs22,tsc,nul);
      fs23            = vec_madd(fs23,tsc,nul);
      fs31            = vec_madd(fs31,tsc,nul);
      fs32            = vec_madd(fs32,tsc,nul);
      fs33            = vec_madd(fs33,tsc,nul);
      
      fs11            = vec_madd(fs11,rinv11,nul);
      fs12            = vec_madd(fs12,rinv12,nul);
      fs13            = vec_madd(fs13,rinv13,nul);
      fs21            = vec_madd(fs21,rinv21,nul);
      fs22            = vec_madd(fs22,rinv22,nul);
      fs23            = vec_madd(fs23,rinv23,nul);
      fs31            = vec_madd(fs31,rinv31,nul);
      fs32            = vec_madd(fs32,rinv32,nul);
      fs33            = vec_madd(fs33,rinv33,nul);      

      fix1            = vec_madd(fs11,dx11,fix1);
      fiy1            = vec_madd(fs11,dy11,fiy1);
      fiz1            = vec_madd(fs11,dz11,fiz1);
      fix2            = vec_madd(fs21,dx21,fix2);
      fiy2            = vec_madd(fs21,dy21,fiy2);
      fiz2            = vec_madd(fs21,dz21,fiz2);
      fix3            = vec_madd(fs31,dx31,fix3);
      fiy3            = vec_madd(fs31,dy31,fiy3);
      fiz3            = vec_madd(fs31,dz31,fiz3);

      fix1            = vec_madd(fs12,dx12,fix1);
      fiy1            = vec_madd(fs12,dy12,fiy1);
      fiz1            = vec_madd(fs12,dz12,fiz1);
      fix2            = vec_madd(fs22,dx22,fix2);
      fiy2            = vec_madd(fs22,dy22,fiy2);
      fiz2            = vec_madd(fs22,dz22,fiz2);
      fix3            = vec_madd(fs32,dx32,fix3);
      fiy3            = vec_madd(fs32,dy32,fiy3);
      fiz3            = vec_madd(fs32,dz32,fiz3);

      fix1            = vec_madd(fs13,dx13,fix1);
      fiy1            = vec_madd(fs13,dy13,fiy1);
      fiz1            = vec_madd(fs13,dz13,fiz1);
      fix2            = vec_madd(fs23,dx23,fix2);
      fiy2            = vec_madd(fs23,dy23,fiy2);
      fiz2            = vec_madd(fs23,dz23,fiz2);
      fix3            = vec_madd(fs33,dx33,fix3);
      fiy3            = vec_madd(fs33,dy33,fiy3);
      fiz3            = vec_madd(fs33,dz33,fiz3);

      fjx1            = vec_nmsub(fs11,dx11,nul);
      fjy1            = vec_nmsub(fs11,dy11,nul);
      fjz1            = vec_nmsub(fs11,dz11,nul);
      fjx2            = vec_nmsub(fs12,dx12,nul);
      fjy2            = vec_nmsub(fs12,dy12,nul);
      fjz2            = vec_nmsub(fs12,dz12,nul);
      fjx3            = vec_nmsub(fs13,dx13,nul);
      fjy3            = vec_nmsub(fs13,dy13,nul);
      fjz3            = vec_nmsub(fs13,dz13,nul);

      fjx1            = vec_nmsub(fs21,dx21,fjx1);
      fjy1            = vec_nmsub(fs21,dy21,fjy1);
      fjz1            = vec_nmsub(fs21,dz21,fjz1);
      fjx2            = vec_nmsub(fs22,dx22,fjx2);
      fjy2            = vec_nmsub(fs22,dy22,fjy2);
      fjz2            = vec_nmsub(fs22,dz22,fjz2);
      fjx3            = vec_nmsub(fs23,dx23,fjx3);
      fjy3            = vec_nmsub(fs23,dy23,fjy3);
      fjz3            = vec_nmsub(fs23,dz23,fjz3);

      fjx1            = vec_nmsub(fs31,dx31,fjx1);
      fjy1            = vec_nmsub(fs31,dy31,fjy1);
      fjz1            = vec_nmsub(fs31,dz31,fjz1);
      fjx2            = vec_nmsub(fs32,dx32,fjx2);
      fjy2            = vec_nmsub(fs32,dy32,fjy2);
      fjz2            = vec_nmsub(fs32,dz32,fjz2);
      fjx3            = vec_nmsub(fs33,dx33,fjx3);
      fjy3            = vec_nmsub(fs33,dy33,fjy3);
      fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
      add_force_to_1_water(faction+j3a,
			   fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3);
    }
    /* update outer data */
    update_i_water_forces(faction+ii3,fshift+is3,
			  fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3);

    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void mcinl0100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float nul;
  vector float dx,dy,dz;
  vector float vnbtot,c6,c12;
  vector float rinvsq,rsq,rinvsix;

  int n,k,k0,ii,is3,ii3,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int ntiA,tja,tjb,tjc,tjd;

  nul=vec_zero();
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinvsq          = do_recip(rsq);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
    }
    /* update outer data */
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void mcinl0300_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float nul,tsc;
  vector float dx,dy,dz;
  vector float vnbtot,c6,c12;
  vector float rinv,r,rsq;
  vector float VVd,VVr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_vonly_4_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_vonly_2_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      do_vonly_1_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
    }
    /* update outer data */
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void mcinl1000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,nul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float rinv,rsq;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      vctot           = vec_madd(qq,rinv,vctot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      vctot           = vec_madd(qq,rinv,vctot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      vctot           = vec_madd(qq,rinv,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl1100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12;
  vector float rinv,rinvsq,rsq,rinvsix;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qq,rinv,vctot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qq,rinv,vctot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qq,rinv,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}




void mcinl2000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vkrf,vcrf,krsq,nul,vcoul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float rinv,rsq;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl2100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,vkrf,vcrf,krsq,vcoul,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12;
  vector float rinv,rinvsq,rsq,rinvsix;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      krsq            = vec_madd(vkrf,rsq,nul);
      vcoul           = vec_add(rinv,krsq);
      vcoul           = vec_sub(vcoul,vcrf);
      vctot           = vec_madd(qq,vcoul,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}






void mcinl3000_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,tsc,nul;
  vector float dx,dy,dz;
  vector float vctot,qq,iq;
  vector float rinv,r,rsq,VVc;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    iq         = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 4 j charges and multiply by iq */
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				 charge+jnrc,charge+jnrd),iq,nul);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 2 j charges and multiply by iq */
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      /* load 1 j charge and multiply by iq */
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl3100_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float vfacel,tsc,nul;
  vector float dx,dy,dz;
  vector float vnbtot,vctot,qq,iq,c6,c12,VVc;
  vector float rinv,r,rinvsq,rsq,rinvsix;
  
  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      rinvsq          = vec_madd(rinv,rinv,nul);     
      r               = vec_madd(rinv,rsq,nul);
      rinvsix         = vec_madd(rinvsq,rinvsq,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq,nul);
      tja             = ntiA+2*type[jnra];
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r,tsc,nul),&VVc);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void mcinl3300_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix,iy,iz,shvec;
  vector float fs,nul,tsc;
  vector float dx,dy,dz,vfacel,vctot;
  vector float vnbtot,c6,c12,iq,qq;
  vector float rinv,r,rsq;
  vector float VVc,VVd,VVr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;

  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  vfacel=load_float_and_splat(&facel);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    shvec      = load_xyz(shiftvec+is3);
    ii         = iinr[n];
    ii3        = 3*ii;
    ix         = load_xyz(pos+ii3);
    vnbtot     = nul;
    vctot      = nul;
    ix         = vec_add(ix,shvec);    
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    splat_xyz_to_vectors(ix,&ix,&iy,&iz);
    ntiA       = 2*ntype*type[ii];
    iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
				  charge+jnrc,charge+jnrd),iq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_vonly_4_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&VVd,&VVr);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
    }
    if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_vonly_2_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&VVd,&VVr);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      k              += 2;
    }
    if((nj1-nj0)%2) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
      dx              = vec_sub(ix,dx);
      dy              = vec_sub(iy,dy);
      dz              = vec_sub(iz,dz);
      rsq             = vec_madd(dx,dx,nul);
      rsq             = vec_madd(dy,dy,rsq);
      rsq             = vec_madd(dz,dz,rsq);
      rinv            = do_invsqrt(rsq);
      r               = vec_madd(rinv,rsq,nul);
      qq = vec_madd(load_1_float(charge+jnra),iq,nul);
      tja             = ntiA+2*type[jnra];
      load_1_pair(nbfp+tja,&c6,&c12);
      do_vonly_1_ljctable_coul_and_lj(VFtab,vec_madd(r,tsc,nul),&VVc,&VVd,&VVr);
      vctot           = vec_madd(qq,VVc,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
    }
    /* update outer data */
    add_vector_to_float(Vnb+gid[n],vnbtot);
    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void mcinl1020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rsqO,rsqH1,rsqH2;
  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);

      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void mcinl1120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float vnbtot,c6,c12,rinvsix;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rsqO,rsqH1,rsqH2;  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      /* load 3 j charges and multiply by iq */
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vctot           = vec_madd(qqO,rinvO,vctot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqH,rinvH1,vctot);
      vctot           = vec_madd(qqH,rinvH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void mcinl2020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float vkrf,vcrf;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float krsqO,krsqH1,krsqH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq,vcoulO,vcoulH1,vcoulH2;
  vector float rinvO,rinvH1,rinvH2,rsqO,rsqH1,rsqH2;  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);

  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      /* load 3 j charges and multiply by iq */
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl2120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float vkrf,vcrf;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul,vcoulO,vcoulH1,vcoulH2;
  vector float vnbtot,c6,c12,rinvsix;
  vector float krsqO,krsqH1,krsqH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rinvsqO,rsqO,rsqH1,rsqH2;  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      krsqO           = vec_madd(vkrf,rsqO,nul);
      krsqH1          = vec_madd(vkrf,rsqH1,nul);
      krsqH2          = vec_madd(vkrf,rsqH2,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vcoulO          = vec_add(rinvO,krsqO);
      vcoulH1         = vec_add(rinvH1,krsqH1);
      vcoulH2         = vec_add(rinvH2,krsqH2);
      vcoulO          = vec_sub(vcoulO,vcrf);
      vcoulH1         = vec_sub(vcoulH1,vcrf);
      vcoulH2         = vec_sub(vcoulH2,vcrf);
      vctot           = vec_madd(qqO,vcoulO,vctot);
      vctot           = vec_madd(qqH,vcoulH1,vctot);
      vctot           = vec_madd(qqH,vcoulH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void mcinl3020_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float tsc,VVcO,VVcH1,VVcH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rO,rH1,rH2,rsqO,rsqH1,rsqH2;  

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  iqO        = vec_madd(load_float_and_splat(charge+iinr[0]),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
    
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);

      /* load 3 j charges and multiply by iq */
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
    
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);

      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl3120_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float vnbtot,c6,c12,rinvsix,rinvsqO;
  vector float tsc,VVcO,VVcH1,VVcH2;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rO,rH1,rH2,rsqO,rsqH1,rsqH2;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      /* load 3 j charges and multiply by iq */
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rinvsqO         = vec_madd(rinvO,rinvO,nul);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);
      rinvsix         = vec_madd(rinvsqO,rinvsqO,nul);
      rinvsix         = vec_madd(rinvsix,rinvsqO,nul);
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rO,tsc,nul),&VVcO);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void mcinl3320_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float tsc;
  vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
  vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
  vector float vfacel,nul;
  vector float vnbtot,c6,c12;
  vector float vctot,qqO,qqH,iqO,iqH,jq;
  vector float rinvO,rinvH1,rinvH2,rsqO,rsqH1,rsqH2;
  vector float rO,rH1,rH2,VVcO,VVcH1,VVcH2,VVd,VVr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;
  int tja,tjb,tjc,tjd;
  
  nul=vec_zero();
  tsc=load_float_and_splat(&tabscale);
  vfacel=load_float_and_splat(&facel);

  ii         = iinr[0];
  iqO        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
  iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
  ntiA       = 2*ntype*type[ii];
  
  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
				 &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),
		       load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
      tjd             = ntiA+2*type[jnrd];
      /* load 4 j charges and multiply by iq */
      jq=load_4_float(charge+jnra,charge+jnrb,charge+jnrc,charge+jnrd);
      load_4_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,nbfp+tjd,&c6,&c12);
      do_vonly_4_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&VVd,&VVr);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      transpose_4_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),
		       load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      tjc             = ntiA+2*type[jnrc];
  
      load_3_pair(nbfp+tja,nbfp+tjb,nbfp+tjc,&c6,&c12);
      jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
      do_vonly_3_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&VVd,&VVr);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      transpose_2_to_3(load_xyz(pos+j3a),
		       load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      tjb             = ntiA+2*type[jnrb];
      /* load 2 j charges and multiply by iq */
      jq=load_2_float(charge+jnra,charge+jnrb);
      load_2_pair(nbfp+tja,nbfp+tjb,&c6,&c12);
      do_vonly_2_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&VVd,&VVr);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
      dOx             = vec_sub(iOx,dH2x);
      dOy             = vec_sub(iOy,dH2y);
      dOz             = vec_sub(iOz,dH2z);
      dH1x            = vec_sub(iH1x,dH2x);
      dH1y            = vec_sub(iH1y,dH2y);
      dH1z            = vec_sub(iH1z,dH2z);
      dH2x            = vec_sub(iH2x,dH2x);
      dH2y            = vec_sub(iH2y,dH2y);
      dH2z            = vec_sub(iH2z,dH2z);
      
      rsqO            = vec_madd(dOx,dOx,nul);
      rsqH1           = vec_madd(dH1x,dH1x,nul);
      rsqH2           = vec_madd(dH2x,dH2x,nul);
      rsqO            = vec_madd(dOy,dOy,rsqO);
      rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
      rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
      rsqO            = vec_madd(dOz,dOz,rsqO);
      rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
      rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
      do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
      rO              = vec_madd(rsqO,rinvO,nul);
      rH1             = vec_madd(rsqH1,rinvH1,nul);
      rH2             = vec_madd(rsqH2,rinvH2,nul);      
      tja             = ntiA+2*type[jnra];
      /* load 1 j charges and multiply by iq */
      jq=load_1_float(charge+jnra);
      load_1_pair(nbfp+tja,&c6,&c12);
      do_vonly_1_ljctable_coul_and_lj(VFtab,vec_madd(rO,tsc,nul),&VVcO,&VVd,&VVr);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
      qqO             = vec_madd(iqO,jq,nul);
      qqH             = vec_madd(iqH,jq,nul);
      vctot           = vec_madd(qqO,VVcO,vctot);
      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vctot           = vec_madd(qqH,VVcH1,vctot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      vctot           = vec_madd(qqH,VVcH2,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}




void mcinl1030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;

  vector float vfacel,vcoul1,vcoul2,vcoul3,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,qqOOt,qqOHt,qqHHt;

 

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];

    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      vctot           = vec_madd(rinv11,qqOO,vctot);
      vctot           = vec_madd(rinv12,qqOH,vctot);
      vctot           = vec_madd(rinv13,qqOH,vctot);
      vctot           = vec_madd(rinv21,qqOH,vctot);
      vctot           = vec_madd(rinv22,qqHH,vctot);
      vctot           = vec_madd(rinv23,qqHH,vctot);
      vctot           = vec_madd(rinv31,qqOH,vctot);
      vctot           = vec_madd(rinv32,qqHH,vctot);
      vctot           = vec_madd(rinv33,qqHH,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}


void mcinl1130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11;

  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnb6,vnb12,vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) { 
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
  
      rinvsq11        = vec_madd(rinv11,rinv11,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      rinvsix         = vec_madd(rinvsix,rinvsix,nul);
      vctot           = vec_madd(rinv11,qqOO,vctot);
      vctot           = vec_madd(rinv12,qqOH,vctot);
      vctot           = vec_madd(rinv13,qqOH,vctot);
      vnbtot          = vec_madd(c12,rinvsix,vnbtot);
      vctot           = vec_madd(rinv21,qqOH,vctot);
      vctot           = vec_madd(rinv22,qqHH,vctot);
      vctot           = vec_madd(rinv23,qqHH,vctot);
      vctot           = vec_madd(rinv31,qqOH,vctot);
      vctot           = vec_madd(rinv32,qqHH,vctot);
      vctot           = vec_madd(rinv33,qqHH,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      rinvsix         = vec_madd(rinvsix,rinvsix,nul);
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vnbtot          = vec_madd(c12t,rinvsix,vnbtot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      rinvsix         = vec_madd(rinvsix,rinvsix,nul);
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vnbtot          = vec_madd(c12t,rinvsix,vnbtot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      rinvsix         = vec_madd(rinvsix,rinvsix,nul);
      vctot           = vec_madd(rinv11,qqOOt,vctot);
      vctot           = vec_madd(rinv12,qqOHt,vctot);
      vctot           = vec_madd(rinv13,qqOHt,vctot);
      vnbtot          = vec_madd(c12t,rinvsix,vnbtot);
      vctot           = vec_madd(rinv21,qqOHt,vctot);
      vctot           = vec_madd(rinv22,qqHHt,vctot);
      vctot           = vec_madd(rinv23,qqHHt,vctot);
      vctot           = vec_madd(rinv31,qqOHt,vctot);
      vctot           = vec_madd(rinv32,qqHHt,vctot);
      vctot           = vec_madd(rinv33,qqHHt,vctot);
    }
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
 
}




void mcinl2030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf)
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;

  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,vkrf,vcrf;
  vector float krsq11,krsq12,krsq13,krsq21,krsq22,krsq23,krsq31,krsq32,krsq33;
  vector float qqOOt,qqOHt,qqHHt;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);      

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);
      
      vctot           = vec_madd(qqOO,rinv11,vctot);
      vctot           = vec_madd(qqOH,rinv12,vctot);
      vctot           = vec_madd(qqOH,rinv13,vctot);
      vctot           = vec_madd(qqOH,rinv21,vctot);
      vctot           = vec_madd(qqHH,rinv22,vctot);
      vctot           = vec_madd(qqHH,rinv23,vctot);
      vctot           = vec_madd(qqOH,rinv31,vctot);
      vctot           = vec_madd(qqHH,rinv32,vctot);
      vctot           = vec_madd(qqHH,rinv33,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);      

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);
      
      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);      

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);
      
      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);      

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);
      
      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl2130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float krf,
	     float crf,
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11,vkrf,vcrf;
  vector float krsq11,krsq12,krsq13,krsq21,krsq22,krsq23,krsq31,krsq32,krsq33;

  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  vkrf=load_float_and_splat(&krf);
  vcrf=load_float_and_splat(&crf);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);

      vctot           = vec_madd(qqOO,rinv11,vctot);
      vctot           = vec_madd(qqOH,rinv12,vctot);
      vctot           = vec_madd(qqOH,rinv13,vctot);
      vctot           = vec_madd(qqOH,rinv21,vctot);
      vctot           = vec_madd(qqHH,rinv22,vctot);
      vctot           = vec_madd(qqHH,rinv23,vctot);
      vctot           = vec_madd(qqOH,rinv31,vctot);
      vctot           = vec_madd(qqHH,rinv32,vctot);
      vctot           = vec_madd(qqHH,rinv33,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);

      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);

      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      krsq11          = vec_madd(vkrf,rsq11,nul);
      krsq12          = vec_madd(vkrf,rsq12,nul);
      krsq13          = vec_madd(vkrf,rsq13,nul);
      krsq21          = vec_madd(vkrf,rsq21,nul);
      krsq22          = vec_madd(vkrf,rsq22,nul);
      krsq23          = vec_madd(vkrf,rsq23,nul);
      krsq31          = vec_madd(vkrf,rsq31,nul);
      krsq32          = vec_madd(vkrf,rsq32,nul);
      krsq33          = vec_madd(vkrf,rsq33,nul);
      
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);

      rinv11          = vec_add(rinv11,krsq11);
      rinv12          = vec_add(rinv12,krsq12);
      rinv13          = vec_add(rinv13,krsq13);
      rinv21          = vec_add(rinv21,krsq21);
      rinv22          = vec_add(rinv22,krsq22);
      rinv23          = vec_add(rinv23,krsq23);
      rinv31          = vec_add(rinv31,krsq31);
      rinv32          = vec_add(rinv32,krsq32);
      rinv33          = vec_add(rinv33,krsq33);

      rinv11          = vec_sub(rinv11,vcrf);
      rinv12          = vec_sub(rinv12,vcrf);
      rinv13          = vec_sub(rinv13,vcrf);
      rinv21          = vec_sub(rinv21,vcrf);
      rinv22          = vec_sub(rinv22,vcrf);
      rinv23          = vec_sub(rinv23,vcrf);
      rinv31          = vec_sub(rinv31,vcrf);
      rinv32          = vec_sub(rinv32,vcrf);
      rinv33          = vec_sub(rinv33,vcrf);

      vctot           = vec_madd(qqOOt,rinv11,vctot);
      vctot           = vec_madd(qqOHt,rinv12,vctot);
      vctot           = vec_madd(qqOHt,rinv13,vctot);
      vctot           = vec_madd(qqOHt,rinv21,vctot);
      vctot           = vec_madd(qqHHt,rinv22,vctot);
      vctot           = vec_madd(qqHHt,rinv23,vctot);
      vctot           = vec_madd(qqOHt,rinv31,vctot);
      vctot           = vec_madd(qqHHt,rinv32,vctot);
      vctot           = vec_madd(qqHHt,rinv33,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



void mcinl3030_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  
  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,tsc;
  vector float VV11c,VV12c,VV13c;
  vector float VV21c,VV22c,VV23c;
  vector float VV31c,VV32c,VV33c;
  vector float qqOOt,qqOHt,qqHHt;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  qO        = load_float_and_splat(charge+iinr[0]);
  qH        = load_float_and_splat(charge+iinr[0]+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_4_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_3_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);
     
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_2_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);
      
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_1_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);
      
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
  }
}



void mcinl3130_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[]) 
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
  vector float rinvsq11;
  vector float vc11,vc12,vc13,vc21,vc22,vc23,vc31,vc32,vc33,tsc,VVc;

  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
  vector float vnbtot,qqOOt,qqOHt,qqHHt,c6t,c12t;
  vector float VV11c,VV12c,VV13c;
  vector float VV21c,VV22c,VV23c;
  vector float VV31c,VV32c,VV33c;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
   
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_vonly_4_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_4_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_nmsub(c6,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_vonly_3_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_3_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
   } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_vonly_2_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_2_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      rinvsq11        = vec_madd(rinv11,rinv11,nul);
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 
      rinvsix         = vec_madd(rinvsix,rinvsq11,nul);

      do_vonly_1_ctable_coul(VFtab,vec_madd(r11,tsc,nul),&VV11c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_1_ctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_nmsub(c6t,rinvsix,vnbtot);
      vnbtot          = vec_madd(c12t,vec_madd(rinvsix,rinvsix,nul),vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}


void mcinl3330_altivec(
	     int nri,
	     int iinr[],
	     int jindex[],
	     int jjnr[],
	     int shift[],
	     float shiftvec[],
	     int gid[],
	     float pos[],
	     float charge[],
	     float facel,
	     float Vc[],
	     int type[],
	     int ntype,
	     float nbfp[],
	     float Vnb[],
	     float tabscale,
	     float VFtab[])
{
  vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
  vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

  vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
  vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
  vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

  vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
  vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
  vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;

  vector float vfacel,nul;
  vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12;
  vector float vnbtot,tsc,qqOOt,qqOHt,qqHHt,c6t,c12t;
  vector float VV11c,VV12c,VV13c;
  vector float VV21c,VV22c,VV23c;
  vector float VV31c,VV32c,VV33c;
  vector float VVd,VVr;

  int n,k,k0,ii,is3,ii3,ntiA,nj0,nj1;
  int jnra,jnrb,jnrc,jnrd,tp,tj;
  int j3a,j3b,j3c,j3d;

  nul=vec_zero();
  vfacel=load_float_and_splat(&facel);
  tsc=load_float_and_splat(&tabscale);
  ii        = iinr[0];
  qO        = load_float_and_splat(charge+ii);
  qH        = load_float_and_splat(charge+ii+1);
  qqOO      = vec_madd(qO,qO,nul);
  qqOH      = vec_madd(qO,qH,nul);
  qqHH      = vec_madd(qH,qH,nul);
  qqOO      = vec_madd(qqOO,vfacel,nul);
  qqOH      = vec_madd(qqOH,vfacel,nul);
  qqHH      = vec_madd(qqHH,vfacel,nul);
  tp        = 2*type[ii];
  tj        = (ntype+1)*tp;
  load_1_pair(nbfp+tj,&c6,&c12);
  c6        = vec_splat(c6,0);
  c12       = vec_splat(c12,0);

  for(n=0;n<nri;n++) {
    is3        = 3*shift[n];
    ii         = iinr[n];
    ii3        = 3*ii;
    load_1_water_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
				 &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
    vctot      = nul;
    vnbtot     = nul;
    nj0        = jindex[n];
    nj1        = jindex[n+1];
    
    for(k=nj0; k<(nj1-3); k+=4) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      jnrd            = jjnr[k+3];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      j3d             = 3*jnrd;
      load_4_water(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_4_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&VVd,&VVr);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_4_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_madd(c6,VVd,vnbtot);
      vnbtot          = vec_madd(c12,VVr,vnbtot);
      vctot           = vec_madd(qqOO,VV11c,vctot);
      vctot           = vec_madd(qqOH,VV12c,vctot);
      vctot           = vec_madd(qqOH,VV13c,vctot);
      vctot           = vec_madd(qqOH,VV21c,vctot);
      vctot           = vec_madd(qqHH,VV22c,vctot);
      vctot           = vec_madd(qqHH,VV23c,vctot);
      vctot           = vec_madd(qqOH,VV31c,vctot);
      vctot           = vec_madd(qqHH,VV32c,vctot);
      vctot           = vec_madd(qqHH,VV33c,vctot);
    } 
    if(k<(nj1-2)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      jnrc            = jjnr[k+2];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      j3c             = 3*jnrc;
      load_3_water(pos+j3a,pos+j3b,pos+j3c,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,4);
      qqOHt           = vec_sld(qqOH,nul,4);
      qqHHt           = vec_sld(qqHH,nul,4);
      c6t             = vec_sld(c6,nul,4);
      c12t            = vec_sld(c12,nul,4);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);

      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_3_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&VVd,&VVr);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_3_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    } else if(k<(nj1-1)) {
      jnra            = jjnr[k];
      jnrb            = jjnr[k+1];
      j3a             = 3*jnra;
      j3b             = 3*jnrb;
      load_2_water(pos+j3a,pos+j3b,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,8);
      qqOHt           = vec_sld(qqOH,nul,8);
      qqHHt           = vec_sld(qqHH,nul,8);
      c6t             = vec_sld(c6,nul,8);
      c12t            = vec_sld(c12,nul,8);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_2_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&VVd,&VVr);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_2_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    } else if(k<nj1) {
      jnra            = jjnr[k];
      j3a             = 3*jnra;
      load_1_water(pos+j3a,
		   &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
      qqOOt           = vec_sld(qqOO,nul,12);
      qqOHt           = vec_sld(qqOH,nul,12);
      qqHHt           = vec_sld(qqHH,nul,12);
      c6t             = vec_sld(c6,nul,12);
      c12t            = vec_sld(c12,nul,12);

      dx11            = vec_sub(ix1,jx1);
      dx12            = vec_sub(ix1,jx2);
      dx13            = vec_sub(ix1,jx3);
      dy11            = vec_sub(iy1,jy1);
      dy12            = vec_sub(iy1,jy2);
      dy13            = vec_sub(iy1,jy3);
      dz11            = vec_sub(iz1,jz1);
      dz12            = vec_sub(iz1,jz2);
      dz13            = vec_sub(iz1,jz3);
      dx21            = vec_sub(ix2,jx1);
      dx22            = vec_sub(ix2,jx2);
      dx23            = vec_sub(ix2,jx3);
      dy21            = vec_sub(iy2,jy1);
      dy22            = vec_sub(iy2,jy2);
      dy23            = vec_sub(iy2,jy3);
      dz21            = vec_sub(iz2,jz1);
      dz22            = vec_sub(iz2,jz2);
      dz23            = vec_sub(iz2,jz3);
      dx31            = vec_sub(ix3,jx1);
      dx32            = vec_sub(ix3,jx2);
      dx33            = vec_sub(ix3,jx3);
      dy31            = vec_sub(iy3,jy1);
      dy32            = vec_sub(iy3,jy2);
      dy33            = vec_sub(iy3,jy3);
      dz31            = vec_sub(iz3,jz1);
      dz32            = vec_sub(iz3,jz2);
      dz33            = vec_sub(iz3,jz3);

      rsq11           = vec_madd(dx11,dx11,nul);
      rsq12           = vec_madd(dx12,dx12,nul);
      rsq13           = vec_madd(dx13,dx13,nul);
      rsq21           = vec_madd(dx21,dx21,nul);
      rsq22           = vec_madd(dx22,dx22,nul);
      rsq23           = vec_madd(dx23,dx23,nul);
      rsq31           = vec_madd(dx31,dx31,nul);
      rsq32           = vec_madd(dx32,dx32,nul);
      rsq33           = vec_madd(dx33,dx33,nul);
      rsq11           = vec_madd(dy11,dy11,rsq11);
      rsq12           = vec_madd(dy12,dy12,rsq12);
      rsq13           = vec_madd(dy13,dy13,rsq13);
      rsq21           = vec_madd(dy21,dy21,rsq21);
      rsq22           = vec_madd(dy22,dy22,rsq22);
      rsq23           = vec_madd(dy23,dy23,rsq23);
      rsq31           = vec_madd(dy31,dy31,rsq31);
      rsq32           = vec_madd(dy32,dy32,rsq32);
      rsq33           = vec_madd(dy33,dy33,rsq33);
      rsq11           = vec_madd(dz11,dz11,rsq11);
      rsq12           = vec_madd(dz12,dz12,rsq12);
      rsq13           = vec_madd(dz13,dz13,rsq13);
      rsq21           = vec_madd(dz21,dz21,rsq21);
      rsq22           = vec_madd(dz22,dz22,rsq22);
      rsq23           = vec_madd(dz23,dz23,rsq23);
      rsq31           = vec_madd(dz31,dz31,rsq31);
      rsq32           = vec_madd(dz32,dz32,rsq32);
      rsq33           = vec_madd(dz33,dz33,rsq33);

      do_9_invsqrt(rsq11,rsq12,rsq13,
		   rsq21,rsq22,rsq23,
		   rsq31,rsq32,rsq33,
		   &rinv11,&rinv12,&rinv13,
		   &rinv21,&rinv22,&rinv23,
		   &rinv31,&rinv32,&rinv33);
      
      r11             = vec_madd(rsq11,rinv11,nul); 
      r12             = vec_madd(rsq12,rinv12,nul); 
      r13             = vec_madd(rsq13,rinv13,nul); 
      r21             = vec_madd(rsq21,rinv21,nul); 
      r22             = vec_madd(rsq22,rinv22,nul); 
      r23             = vec_madd(rsq23,rinv23,nul); 
      r31             = vec_madd(rsq31,rinv31,nul); 
      r32             = vec_madd(rsq32,rinv32,nul); 
      r33             = vec_madd(rsq33,rinv33,nul); 

      do_vonly_1_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
				 &VV11c,&VVd,&VVr);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
      do_vonly_1_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

      vnbtot          = vec_madd(c6t,VVd,vnbtot);
      vnbtot          = vec_madd(c12t,VVr,vnbtot);
      vctot           = vec_madd(qqOOt,VV11c,vctot);
      vctot           = vec_madd(qqOHt,VV12c,vctot);
      vctot           = vec_madd(qqOHt,VV13c,vctot);
      vctot           = vec_madd(qqOHt,VV21c,vctot);
      vctot           = vec_madd(qqHHt,VV22c,vctot);
      vctot           = vec_madd(qqHHt,VV23c,vctot);
      vctot           = vec_madd(qqOHt,VV31c,vctot);
      vctot           = vec_madd(qqHHt,VV32c,vctot);
      vctot           = vec_madd(qqHHt,VV33c,vctot);
    }
    /* update outer data */
    add_vector_to_float(Vc+gid[n],vctot);
    add_vector_to_float(Vnb+gid[n],vnbtot);
  }
}



