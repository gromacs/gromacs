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
 * GRoups of Organic Molecules in ACtion for Science
 */


#ifndef _ppc_altivec_h
#define _ppc_altivec_h

#include<stdio.h>

static char *SRCID_ppc_altivec_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


static void printvec(vector float v)
{
  int i;
  printf(" ");

  for(i=0;i<4;i++)
    printf("%8.5f ",*(((float *)&v)+i));
  printf("\n");
}


/* Altivec utility routines */
static void set_non_java_mode(void)
{
  vector unsigned short vsr1,vsr2;
  vsr1=vec_mfvscr();
  vsr2=(vector unsigned short)vec_sl(vec_splat_u32(1),vec_splat_u32(16));
  vsr1=vec_or(vsr1,vsr2);
  vec_mtvscr(vsr1);
}  

/* Simple numerical constants */
static inline vector float vec_zero(void)
{
  return vec_ctf(vec_splat_u32(0),0);
}

static inline vector float vec_half(void)
{  /* 0.5 */
  return vec_ctf(vec_splat_u32(1),1);
}

static inline vector float vec_one(void)
{
  return vec_ctf(vec_splat_u32(1),0);
}

static inline vector float vec_two(void)
{
  return vec_ctf(vec_splat_u32(2),0);
}

static inline vector float vec_three(void)
{
  return vec_ctf(vec_splat_u32(3),0);
}

static inline vector float vec_six(void)
{
  return vec_ctf(vec_splat_u32(6),0);
}

static inline vector float vec_twelve(void)
{
  return vec_ctf(vec_splat_u32(12),0);
}


/* load three consecutive floats. We never access the fourth,
 * so this is safe even at the end of an array. 
 */
static inline vector float load_xyz(float *address)
{
  vector float c1,c2,c3;
  vector unsigned char perm;
  
  perm              = vec_lvsl( 0, address ); 
  c1                = vec_lde( 0, address );
  c2                = vec_lde( 4, address );
  c3                = vec_lde( 8, address );
  c1                = vec_perm(c1,c1,perm);
  c2                = vec_perm(c2,c2,perm);
  c3                = vec_perm(c3,c3,perm);
  c2                = vec_sld(c2,c2,4);
  c3                = vec_sld(c3,c3,8);
  c1                = vec_mergeh(c1,c3);
  
  return vec_mergeh(c1,c2);
}


/* load four consecutive floats (from unaligned address) */
static inline vector float load_vector_unaligned(float *address)
{
  vector unsigned char perm;
  vector float low,high;
  
  perm              = vec_lvsl( 0, (int *) address ); 
  low               = vec_ld(  0, address ); 
  high              = vec_ld( 16, address ); 
  
  return vec_perm(low,high,perm);
}

/* load a single float and spread it to all vector elements */
static inline vector float load_float_and_splat(float *address)
{
  vector unsigned char perm;
  vector float tmp;

  tmp               = vec_lde(0,address);
  perm              = vec_lvsl(0,address);
  tmp               = vec_perm(tmp,tmp,perm);

  return vec_splat(tmp,0);
}


/* load four non-consecutive floats into vector [mem1 mem2 mem3 mem4] */
static inline vector float load_4_float(float *float1,
					float *float2,
					float *float3,
					float *float4)
{
  vector unsigned char xshift = vec_lvsl( 12, float1 ); 
  vector unsigned char yshift = vec_lvsl( 12, float2 ); 
  vector unsigned char zshift = vec_lvsl( 0, float3 ); 
  vector unsigned char wshift = vec_lvsl( 0, float4 ); 
  
  vector float X = vec_lde( 0, float1 ); 
  vector float Y = vec_lde( 0, float2 ); 
  vector float Z = vec_lde( 0, float3 ); 
  vector float W = vec_lde( 0, float4 ); 
  
  X = vec_perm( X, X, xshift); 
  Y = vec_perm( Y, Y, yshift); 
  Z = vec_perm( Z, Z, zshift); 
  W = vec_perm( W, W, wshift); 
  
  X = vec_mergeh( X, Y ); 
  Z = vec_mergeh( Z, W ); 

  return vec_sld( X, Z, 8 );
}

/* load four non-consecutive floats into vector [mem1 mem2 mem3 mem4] */
static inline vector float load_3_float(float *float1,
					float *float2,
					float *float3)
{
  vector unsigned char xshift = vec_lvsl( 12, float1 ); 
  vector unsigned char yshift = vec_lvsl( 12, float2 ); 
  vector unsigned char zshift = vec_lvsl( 0, float3 ); 

  vector float X = vec_lde( 0, float1 ); 
  vector float Y = vec_lde( 0, float2 ); 
  vector float Z = vec_lde( 0, float3 ); 
  
  X = vec_perm( X, X, xshift); 
  Y = vec_perm( Y, Y, yshift); 
  Z = vec_perm( Z, Z, zshift); 
  
  X = vec_mergeh( X, Y ); 

  return vec_sld( X, Z, 8 );
}


/* load two non-consecutive floats into vector [mem1 mem2  0 0] */
static inline vector float load_2_float(float *float1,
					float *float2)
{
  vector unsigned char xshift = vec_lvsl( 8, float1 ); 
  vector unsigned char yshift = vec_lvsl( 8, float2 ); 
  
  vector float X = vec_lde( 0, float1 ); 
  vector float Y = vec_lde( 0, float2 ); 
  
  X = vec_perm( X, X, xshift); 
  Y = vec_perm( Y, Y, yshift); 
  
  return vec_mergel( X, Y ); 
}

/* load one float into vector [mem1 0 0 0] */
static inline vector float load_1_float(float *float1)
{
  vector unsigned char xshift = vec_lvsl( 4, float1 ); 
  
  vector float X = vec_lde( 0, float1 ); 
  X = vec_perm( X, X, xshift); 
    
  return vec_sld(X,vec_zero(),12);
}




/* load lj parameters for four atoms and store in two vectors. 
 * The nbfp memory is aligned on an 8-byte boundary and consists
 * of pairs, so the two parameters are either in the upper or 
 * lower half of a 16-byte structure.
 */
static inline void load_4_pair(float *pair1,
			       float *pair2,
			       float *pair3,
			       float *pair4,
			       vector float *c6,
			       vector float *c12)
{
  vector float X = vec_ld( 0,pair1); /* c6a c12a  */
  vector float Y = vec_ld( 0,pair2); /* c6b c12b  */
  vector float Z = vec_ld( 0,pair3); /* c6c c12c  */
  vector float W = vec_ld( 0,pair4); /* c6d c12d  */
  vector unsigned char perm1 = vec_lvsl(0,pair1);
  vector unsigned char perm2 = vec_lvsl(0,pair2);
  vector unsigned char perm3 = vec_lvsl(0,pair3);
  vector unsigned char perm4 = vec_lvsl(0,pair4);
  X = vec_perm(X,X,perm1);
  Y = vec_perm(Y,Y,perm2);
  Z = vec_perm(Z,Z,perm3);
  W = vec_perm(W,W,perm4);

  X = vec_mergeh(X,Z); /* c6a c6c c12a c12c */
  Y = vec_mergeh(Y,W); /* c6b c6d c12b c12d */
  
  *c6  = vec_mergeh(X,Y);
  *c12 = vec_mergel(X,Y);  
}


/* load lj parameters for three atoms and store in two vectors 
 * See load_4_pair for alignment requirements.
 */
static inline void load_3_pair(float *pair1,
			       float *pair2,
			       float *pair3,
			       vector float *c6,
			       vector float *c12)
{
  vector float X = vec_ld( 0,pair1); /* c6a c12a  */
  vector float Y = vec_ld( 0,pair2); /* c6b c12b  */
  vector float Z = vec_ld( 0,pair3); /* c6c c12c  */
  vector unsigned char perm1 = vec_lvsl(0,pair1);
  vector unsigned char perm2 = vec_lvsl(0,pair2);
  vector unsigned char perm3 = vec_lvsl(0,pair3);
  X = vec_perm(X,X,perm1);
  Y = vec_perm(Y,Y,perm2);
  Z = vec_perm(Z,Z,perm3);

  X = vec_mergeh(X,Z); /* c6a c6c c12a c12c */
  Y = vec_mergeh(Y,vec_zero()); /* c6b  0  c12b  0 */
  
  *c6  = vec_mergeh(X,Y);
  *c12 = vec_mergel(X,Y);  
}


/* load lj parameters for two atoms and store in two vectors 
 * See load_4_pair for alignment requirements.
 */
static inline void load_2_pair(float *pair1,
				float *pair2,
				vector float *c6,
				vector float *c12)
{
  vector float X = vec_ld( 0,pair1); /* c6a c12a  */
  vector float Y = vec_ld( 0,pair2); /* c6b c12b  */
  vector unsigned char perm1 = vec_lvsl(0,pair1);
  vector unsigned char perm2 = vec_lvsl(0,pair2);
  X = vec_perm(X,X,perm1);
  Y = vec_perm(Y,Y,perm2);

  X = vec_mergeh(X,vec_zero()); /* c6a 0 c12a 0 */
  Y = vec_mergeh(Y,vec_zero()); /* c6b 0 c12b 0 */
  
  *c6  = vec_mergeh(X,Y);
  *c12 = vec_mergel(X,Y);  
}

/* load lj parameters for one atom and store in two vectors 
 * See load_4_pair for alignment requirements.
 */
static inline void load_1_pair(float *pair1,
			       vector float *c6,
			       vector float *c12)
{
  vector float X = vec_ld( 0,pair1); /* c6a c12a  */
  vector unsigned char perm1 = vec_lvsl(0,pair1);
  X = vec_perm(X,X,perm1);
  X = vec_mergeh(X,vec_zero()); /* c6a 0 c12a 0 */
  
  *c6  = vec_mergeh(X,vec_zero());
  *c12 = vec_mergel(X,vec_zero());  
}					   


/*
 * permute [x y z ?] into [x x x x] [y y y y] [z z z z] 
 * xreg can be the same as xyzreg, but not yreg/zreg!
 */
static inline void splat_xyz_to_vectors(vector float xyzreg,
					vector float *xreg,
					vector float *yreg,
					vector float *zreg)
{
  *zreg                = vec_splat(xyzreg,2);
  *yreg                = vec_splat(xyzreg,1);
  *xreg                = vec_splat(xyzreg,0);
}


/* The following transpose routines do not fill with zero padding:
 *
 * 1_to_3
 * 2_to_3
 * 1_to_4
 * 2_to_4
 *
 * ...but the rest do. 
 */


/* move 4*[x y z ?] into [x1 x2 x3 x4] [y1 y2 y3 y4] [z1 z2 z3 z4]  
 */
static inline void transpose_4_to_3(vector float xyz1,
				    vector float xyz2,
				    vector float xyz3,
				    vector float xyz4,
				    vector float *xvector,
				    vector float *yvector,
				    vector float *zvector)
{
  *yvector    = vec_mergeh(xyz1,xyz3);       /* [x1 x3 y1 y3] */
  *zvector    = vec_mergeh(xyz2,xyz4);       /* [x2 x4 y2 y4] */
  xyz1       = vec_mergel(xyz1,xyz3);       /* [z1 z3  ?  ?] */
  xyz2       = vec_mergel(xyz2,xyz4);       /* [z2 z4  ?  ?] */
  *xvector    = vec_mergeh(*yvector,*zvector); /* [x1 x2 x3 x4] */
  *yvector    = vec_mergel(*yvector,*zvector); /* [y1 y2 y3 y4] */
  *zvector    = vec_mergeh(xyz1,xyz2);       /* [z1 z2 z3 z4] */
}


/* move 2*[x y z ?] into [x1 x2 ? ?] [y1 y2 ? ?] [z1 z2 ? ?]  
 */
static inline void transpose_2_to_3(vector float xyz1,
				    vector float xyz2,
				    vector float *xvector,
				    vector float *yvector,
				    vector float *zvector)
{
  *xvector = vec_mergeh(xyz1,xyz2);            /* x1 x2 y1 y2 */
  *zvector = vec_mergel(xyz1,xyz2);            /* z1 z2  ?  ? */
  *yvector = vec_sld(*xvector,*xvector,8); /* y1 y2 0 0 */
}


/* move [x y z ?] into [x1 ? ? ?] [y1 ? ? ?] [z1 ? ? ?]  
 */
static inline void transpose_1_to_3(vector float xyz1,
				    vector float *xvector,
				    vector float *yvector,
				    vector float *zvector)
{
  /* simply use splat, since elem 2,3,4 dont matter. */
  *xvector           = vec_splat(xyz1,0);   /* x1 x1 x1 x1 */
  *yvector           = vec_splat(xyz1,1);   /* y1 y1 y1 y1 */
  *zvector           = vec_splat(xyz1,2);   /* z1 z1 z1 z1 */
}




/* move [x1 x2 x3 x4] [y1 y2 y3 y4] [z1 z2 z3 z4] to 4* [x y z 0]
 */
static inline void transpose_3_to_4(vector float xvector,
				    vector float yvector,
				    vector float zvector,
				    vector float *xyz1,
				    vector float *xyz2,
				    vector float *xyz3,
				    vector float *xyz4)
{
  vector float tmp1,tmp2;
  vector float nul=vec_zero();

  *xyz2       = vec_mergeh(xvector,zvector); /* [x1 z1 x2 z2 ] */
  tmp1       = vec_mergeh(yvector,nul);     /* [y1  0 y2  0 ] */
  *xyz4       = vec_mergel(xvector,zvector); /* [x3 z3 x4 z4 ] */
  tmp2       = vec_mergel(yvector,nul);     /* [y3  0 y4  0 ] */
  
  *xyz1       = vec_mergeh(*xyz2,tmp1);       /* [x1 y1 z1 0 ] */
  *xyz2       = vec_mergel(*xyz2,tmp1);       /* [x2 y2 z2 0 ] */
  *xyz3       = vec_mergeh(*xyz4,tmp2);       /* [x3 y3 z3 0 ] */
  *xyz4       = vec_mergel(*xyz4,tmp2);       /* [x4 y4 z4 0 ] */
}



/* move [x1 x2 ? ?] [y1 y2 ? ?] [z1 z2 ? ?] to 2* [x y z 0]
 */
static inline void transpose_3_to_2(vector float xvector,
				    vector float yvector,
				    vector float zvector,
				    vector float *xyz1,
				    vector float *xyz2)
{
  vector float tmp1;

  *xyz2       = vec_mergeh(xvector,zvector); /* [x1 z1 x2 z2 ] */
  tmp1       = vec_mergeh(yvector,vec_zero());     /* [y1  0 y2  0 ] */
  
  *xyz1       = vec_mergeh(*xyz2,tmp1);       /* [x1 y1 z1 0 ] */
  *xyz2       = vec_mergel(*xyz2,tmp1);       /* [x2 y2 z2 0 ] */
}


/* move [x1 ? ? ?] [y1 ? ? ?] [z1 ? ? ?] to 1* [x y z 0]
 */
static inline void transpose_3_to_1(vector float xvector,
				    vector float yvector,
				    vector float zvector,
				    vector float *xyz1)
{
  *xyz1       = vec_mergeh(xvector,zvector); /* [x1 z1 ? ? ] */
  yvector     = vec_mergeh(yvector,vec_zero()); /* [y1 0 ? 0] */
  *xyz1       = vec_mergeh(*xyz1,yvector);   /* [x1 y1 z1 0] */
}



static inline void transpose_4_to_4(vector float in1,
				    vector float in2,
				    vector float in3,
				    vector float in4,
				    vector float *out1,
				    vector float *out2,
				    vector float *out3,
				    vector float *out4)
{
  *out2    = vec_mergeh(in1,in3);       /* [x1 x3 y1 y3] */
  *out3    = vec_mergeh(in2,in4);       /* [x2 x4 y2 y4] */
  in1       = vec_mergel(in1,in3);      /* [z1 z3 w1 w3] */
  in2       = vec_mergel(in2,in4);      /* [z2 z4 w2 w4] */

  *out1    = vec_mergeh(*out2,*out3);   /* [x1 x2 x3 x4] */
  *out2    = vec_mergel(*out2,*out3);   /* [y1 y2 y3 y4] */
  *out3    = vec_mergeh(in1,in2);       /* [z1 z2 z3 z4] */
  *out4    = vec_mergel(in1,in2);       /* [w1 w2 w3 w4] */   
}


static inline void transpose_4_to_2(vector float in1,
				    vector float in2,
				    vector float in3,
				    vector float in4,
				    vector float *out1,
				    vector float *out2)
{
  vector float tmp;

  tmp      = vec_mergeh(in1,in3);   /* [x1 x3 y1 y3] */ 
  *out2    = vec_mergeh(in2,in4);   /* [x2 x4 y2 y4] */ 
  *out1    = vec_mergeh(tmp,*out2); /* [x1 x2 x3 x4] */ 
  *out2    = vec_mergel(tmp,*out2); /* [y1 y2 y3 y4] */ 
}

static inline void transpose_4_to_1(vector float in1,
				    vector float in2,
				    vector float in3,
				    vector float in4,
				    vector float *out1)
{
  vector float tmp;

  tmp      = vec_mergeh(in1,in3);   /* [x1 x3  ?  ?] */ 
  *out1    = vec_mergeh(in2,in4);   /* [x2 x4  ?  ?] */ 
  *out1    = vec_mergeh(tmp,*out1); /* [x1 x2 x3 x4] */ 
}



static inline void transpose_2_to_4(vector float in1,
				    vector float in2,
				    vector float *out1,
				    vector float *out2,
				    vector float *out3,
				    vector float *out4)
{
  *out1    = vec_mergeh(in1,in2);       /* [x1 x2 y1 y2] */
  *out3    = vec_mergel(in1,in2);       /* [z1 z2 w1 w2] */
  *out2    = vec_sld(*out1,*out1,8);/* [y1 y2  0 0] */
  *out4    = vec_sld(*out3,*out3,8);/* [w1 w2  0 0] */  
}


static inline void transpose_1_to_4(vector float in1,
				    vector float *out1,
				    vector float *out2,
				    vector float *out3,
				    vector float *out4)
{
  *out1    = vec_splat(in1,0); /* [x1 x1 x1 x1] */
  *out2    = vec_splat(in1,1); /* [y1 y1 y1 y1] */ 
  *out3    = vec_splat(in1,2); /* [z1 z1 z1 z1] */
  *out4    = vec_splat(in1,3); /* [w1 w1 w1 w1] */
}


/* Add the contents in xyz to an unaligned address location.
 */
static inline void add_xyz_to_mem(float *address,vector float vdata)
{
  vector float c1,c2,c3;

  c1                = vec_lde( 0, address);
  c2                = vec_lde( 4, address);
  c3                = vec_lde( 8, address);
  c1                = vec_add(c1,vec_splat(vdata,0));
  c2                = vec_add(c2,vec_splat(vdata,1));
  c3                = vec_add(c3,vec_splat(vdata,2));  
  vec_ste( c1,  0, address);
  vec_ste( c2,  4, address);
  vec_ste( c3,  8, address);
}


static inline void add_vector_to_float(float *address,vector float vdata)
{
  vector float tmp;
  
  tmp              = vec_sld(vdata,vdata,8);
  vdata            = vec_add(vdata,tmp);
  tmp              = vec_sld(vdata,vdata,4);
  vdata            = vec_add(vdata,tmp);
  tmp              = vec_lde(0,address);
  /* all four positions in vdata contain the sum */
  tmp              = vec_add(tmp,vdata);
  vec_ste(tmp,0,address);
}


/*
 * load entire water molecule to vectors. Water is special; 
 * the allocated memory is 8-byte aligned, so with nine atoms
 * we can always read 3 16-byte chunks (rounded down by vec_ld)
 * without going outside the current page. The extra elements 
 * before and after the water dont matter as long as we dont
 * try to write to them.
 */
static inline void load_1_water_shift_and_splat(float *address,
						float *shiftvec,
						vector float *Ox,
						vector float *Oy,
						vector float *Oz,
						vector float *H1x,
						vector float *H1y,
						vector float *H1z,
						vector float *H2x,
						vector float *H2y,
						vector float *H2z)
{
  vector unsigned char perm;
  vector float c1,c2,c3,sha,shb,sh1,sh2,sh3;

  /* load shift */
  shb               = load_xyz(shiftvec);  /* [ shX shY shZ -] */
  /* load the coordinates */
  perm              = vec_lvsl( 0, (int *) address ); 
  c1                = vec_ld(  0, address ); 
  c2                = vec_ld( 16, address ); 
  c3                = vec_ld( 32, address ); 
  
  sha               = vec_sld(shb,shb,12);    /* [ - shX shY shZ ] */

  sh1               = vec_sld(sha,shb,4);  /* [ shX shY shZ shX ] */
  sh2               = vec_sld(sha,shb,8);  /* [ shY shZ shX shY ] */
  sh3               = vec_sld(sha,shb,12);  /* [ shZ shX shY shZ ] */

  c1                = vec_perm(c1,c2,perm);      /*  Ox  Oy  Oz H1x */
  c2                = vec_perm(c2,c3,perm);      /* H1y H1z H2x H2y */
  c3                = vec_perm(c3,c3,perm);      /* H2z   -   -   - */
  c1                = vec_add(c1,sh1);
  c2                = vec_add(c2,sh2);
  c3                = vec_add(c3,sh3);

  *Ox               = vec_splat(c1,0);
  *Oy               = vec_splat(c1,1);
  *Oz               = vec_splat(c1,2);
  *H1x              = vec_splat(c1,3);
  *H1y              = vec_splat(c2,0);
  *H1z              = vec_splat(c2,1);
  *H2x              = vec_splat(c2,2);
  *H2y              = vec_splat(c2,3);
  *H2z              = vec_splat(c3,0);
}

static inline void load_4_water(float *address1,
				float *address2,
				float *address3,
				float *address4,
				vector float *Ox,
				vector float *Oy,
				vector float *Oz,
				vector float *H1x,
				vector float *H1y,
				vector float *H1z,
				vector float *H2x,
				vector float *H2y,
				vector float *H2z)
{
  vector unsigned char perm;
  vector float tmp1,tmp2,tmp3;
  
  vector float c1a,c2a,c3a;
  vector float c1b,c2b,c3b;
  vector float c1c,c2c,c3c;
  vector float c1d,c2d,c3d;
  
  /* load the coordinates */
  perm              = vec_lvsl( 0, address1 ); 
  tmp1              = vec_ld(  0, address1 ); 
  tmp2              = vec_ld( 16, address1 ); 
  tmp3              = vec_ld( 32, address1 ); 
  c1a               = vec_perm(tmp1,tmp2,perm);     /*  Oxa  Oya  Oza H1xa */
  c2a               = vec_perm(tmp2,tmp3,perm);     /* H1ya H1za H2xa H2ya */
  c3a               = vec_perm(tmp3,tmp3,perm);     /* H2za   -   -   - */


  perm              = vec_lvsl( 0, address3 );
  tmp1              = vec_ld(  0, address3 );
  tmp2              = vec_ld( 16, address3 );
  tmp3              = vec_ld( 32, address3 );
  c1c               = vec_perm(tmp1,tmp2,perm);     /*  Oxc  Oyc  Ozc H1xc */
  c2c               = vec_perm(tmp2,tmp3,perm);     /* H1yc H1zc H2xc H2yc */
  c3c               = vec_perm(tmp3,tmp3,perm);     /* H2zc   -   -   - */


  perm              = vec_lvsl( 0, address2 );
  tmp1              = vec_ld(  0, address2 );
  tmp2              = vec_ld( 16, address2 );
  tmp3              = vec_ld( 32, address2 );
  c1b               = vec_perm(tmp1,tmp2,perm);     /*  Oxb  Oyb  Ozb H1xb */
  c2b               = vec_perm(tmp2,tmp3,perm);     /* H1yb H1zb H2xb H2yb */
  c3b               = vec_perm(tmp3,tmp3,perm);     /* H2zb   -   -   - */


  perm              = vec_lvsl( 0, address4 );
  tmp1              = vec_ld(  0, address4 );
  tmp2              = vec_ld( 16, address4 );
  tmp3              = vec_ld( 32, address4 );
  c1d               = vec_perm(tmp1,tmp2,perm);     /*  Oxd  Oyd  Ozd H1xd */
  c2d               = vec_perm(tmp2,tmp3,perm);     /* H1yd H1zd H2xd H2yd */
  c3d               = vec_perm(tmp3,tmp3,perm);     /* H2zd   -   -   - */

  /* permute things */
  tmp1              = vec_mergeh(c1a,c1c);          /*  Oxa  Oxc  Oya  Oyc */
  c1a               = vec_mergel(c1a,c1c);          /*  Oza  Ozc H1xa H1xc */
  c1c               = vec_mergeh(c1b,c1d);          /*  Oxb  Oxd  Oyb  Oyd */
  c1b               = vec_mergel(c1b,c1d);          /*  Ozb  Ozd H1xb H1xd */

  c1d               = vec_mergeh(c2a,c2c);          /* H1ya H1yc H1za H1zc */
  c2a               = vec_mergel(c2a,c2c);          /* H2xa H2xc H2ya H2yc */
  c2c               = vec_mergeh(c2b,c2d);          /* H1yb H1yd H1zb H1zd */
  c2b               = vec_mergel(c2b,c2d);          /* H2xb H2xd H2yb H2yd */

  c3a               = vec_mergeh(c3a,c3c);          /* H2za H2zc   -    -  */
  c3b               = vec_mergeh(c3b,c3d);          /* H2zb H2zd   -    -  */
  
  *Ox               = vec_mergeh(tmp1,c1c);         /*  Oxa  Oxb  Oxc  Oxd */
  *Oy               = vec_mergel(tmp1,c1c);         /*  Oya  Oyb  Oyc  Oyd */
  *Oz               = vec_mergeh(c1a,c1b);          /*  Oza  Ozb  Ozc  Ozd */
  *H1x              = vec_mergel(c1a,c1b);          /* H1xa H1xb H1xc H1xd */
  *H1y              = vec_mergeh(c1d,c2c);          /* H1ya H1yb H1yc H1yd */
  *H1z              = vec_mergel(c1d,c2c);          /* H1za H1zb H1zc H1zd */
  *H2x              = vec_mergeh(c2a,c2b);          /* H2xa H2xb H2xc H2xd */
  *H2y              = vec_mergel(c2a,c2b);          /* H2ya H2yb H2yc H2yd */
  *H2z              = vec_mergeh(c3a,c3b);          /* H2za H2zb H2zc H2zd */
}



static inline void load_3_water(float *address1,
				float *address2,
				float *address3,
				vector float *Ox,
				vector float *Oy,
				vector float *Oz,
				vector float *H1x,
				vector float *H1y,
				vector float *H1z,
				vector float *H2x,
				vector float *H2y,
				vector float *H2z)
{
  vector unsigned char perm;
  vector float tmp1,tmp2,tmp3;
  vector float c1a,c2a,c3a;
  vector float c1b,c2b,c3b;
  vector float c1c,c2c,c3c;
  
  /* load the coordinates */
  perm              = vec_lvsl( 0, (int *) address1 ); 
  tmp1              = vec_ld(  0, address1 ); 
  tmp2              = vec_ld( 16, address1 ); 
  tmp3              = vec_ld( 32, address1 ); 
  c1a               = vec_perm(tmp1,tmp2,perm);     /*  Oxa  Oya  Oza H1xa */
  c2a               = vec_perm(tmp2,tmp3,perm);     /* H1ya H1za H2xa H2ya */
  c3a               = vec_perm(tmp3,tmp3,perm);     /* H2za   -   -   - */

  perm              = vec_lvsl( 0, (int *) address3 ); 
  tmp1              = vec_ld(  0, address3 ); 
  tmp2              = vec_ld( 16, address3 ); 
  tmp3              = vec_ld( 32, address3 ); 
  c1c               = vec_perm(tmp1,tmp2,perm);     /*  Oxc  Oyc  Ozc H1xc */
  c2c               = vec_perm(tmp2,tmp3,perm);     /* H1yc H1zc H2xc H2yc */
  c3c               = vec_perm(tmp3,tmp3,perm);     /* H2zc   -   -   - */

  perm              = vec_lvsl( 0, (int *) address2 ); 
  tmp1              = vec_ld(  0, address2 ); 
  tmp2              = vec_ld( 16, address2 ); 
  tmp3              = vec_ld( 32, address2 );  
  c1b               = vec_perm(tmp1,tmp2,perm);     /*  Oxb  Oyb  Ozb H1xb */
  c2b               = vec_perm(tmp2,tmp3,perm);     /* H1yb H1zb H2xb H2yb */
  c3b               = vec_perm(tmp3,tmp3,perm);     /* H2zb   -   -   - */

  /* permute things */
  tmp1              = vec_mergeh(c1a,c1c);          /*  Oxa  Oxc  Oya  Oyc */
  c1a               = vec_mergel(c1a,c1c);          /*  Oza  Ozc H1xa H1xc */
  c1c               = vec_mergeh(c1b,c1b);          /*  Oxb   -   Oyb   -  */
  c1b               = vec_mergel(c1b,c1b);          /*  Ozb   -  H1xb   -  */

  tmp2              = vec_mergeh(c2a,c2c);          /* H1ya H1yc H1za H1zc */
  c2a               = vec_mergel(c2a,c2c);          /* H2xa H2xc H2ya H2yc */
  c2c               = vec_mergeh(c2b,c2b);          /* H1yb   -  H1zb   -  */
  c2b               = vec_mergel(c2b,c2b);          /* H2xb   -  H2yb   -  */

  c3a               = vec_mergeh(c3a,c3c);          /* H2za H2zc   -    -  */
  
  *Ox               = vec_mergeh(tmp1,c1c);         /*  Oxa  Oxb  Oxc   -  */
  *Oy               = vec_mergel(tmp1,c1c);         /*  Oya  Oyb  Oyc   -  */
  *Oz               = vec_mergeh(c1a,c1b);          /*  Oza  Ozb  Ozc   -  */
  *H1x              = vec_mergel(c1a,c1b);          /* H1xa H1xb H1xc   -  */
  *H1y              = vec_mergeh(tmp2,c2c);         /* H1ya H1yb H1yc   -  */
  *H1z              = vec_mergel(tmp2,c2c);         /* H1za H1zb H1zc   -  */
  *H2x              = vec_mergeh(c2a,c2b);          /* H2xa H2xb H2xc   -  */
  *H2y              = vec_mergel(c2a,c2b);          /* H2ya H2yb H2yc   -  */
  *H2z              = vec_mergeh(c3a,c3b);          /* H2za H2zb H2zc   -  */
}



static inline void load_2_water(float *address1,
				float *address2,
				vector float *Ox,
				vector float *Oy,
				vector float *Oz,
				vector float *H1x,
				vector float *H1y,
				vector float *H1z,
				vector float *H2x,
				vector float *H2y,
				vector float *H2z)
{
  vector unsigned char perm;
  vector float tmp1,tmp2,tmp3;
  vector float c1a,c2a,c3a;
  vector float c1b,c2b,c3b;
  
  /* load the coordinates */
  perm              = vec_lvsl( 0, (int *) address1 ); 
  tmp1              = vec_ld(  0, address1 ); 
  tmp2              = vec_ld( 16, address1 ); 
  tmp3              = vec_ld( 32, address1 ); 
  c1a               = vec_perm(tmp1,tmp2,perm);     /*  Oxa  Oya  Oza H1xa */
  c2a               = vec_perm(tmp2,tmp3,perm);     /* H1ya H1za H2xa H2ya */
  c3a               = vec_perm(tmp3,tmp3,perm);     /* H2za   -   -   - */

  perm              = vec_lvsl( 0, (int *) address2 ); 
  tmp1              = vec_ld(  0, address2 ); 
  tmp2              = vec_ld( 16, address2 ); 
  tmp3              = vec_ld( 32, address2 );  
  c1b               = vec_perm(tmp1,tmp2,perm);     /*  Oxb  Oyb  Ozb H1xb */
  c2b               = vec_perm(tmp2,tmp3,perm);     /* H1yb H1zb H2xb H2yb */
  c3b               = vec_perm(tmp3,tmp3,perm);     /* H2zb   -   -   - */

  /* never mind what we get in the two remaining elements */
  *Ox               = vec_mergeh(c1a,c1b);          /*  Oxa  Oxb  Oya  Oyb */ 
  *H1y              = vec_mergeh(c2a,c2b);          /* H1ya H1yb H1za H1zb */
  *H2z              = vec_mergeh(c3a,c3b);          /* H2za H2zb   -    -  */  
  *Oz               = vec_mergel(c1a,c1b);          /*  Oza  Ozb H1xa H1xb */ 
  *H2x              = vec_mergel(c2a,c2b);          /* H2xa H2xb H2ya H2yb */
  *Oy               = vec_sld(*Ox,*Ox,8);
  *H1z              = vec_sld(*H1y,*H1y,8); 
  *H1x              = vec_sld(*Oz,*Oz,8);
  *H2y              = vec_sld(*H2x,*H2x,8);
}


static inline void load_1_water(float *address1,
				vector float *Ox,
				vector float *Oy,
				vector float *Oz,
				vector float *H1x,
				vector float *H1y,
				vector float *H1z,
				vector float *H2x,
				vector float *H2y,
				vector float *H2z)
{
  vector unsigned char perm;
  vector float tmp1,tmp2,tmp3;

  /* load the coordinates */
  perm              = vec_lvsl( 0, (int *) address1 ); 
  tmp1              = vec_ld(  0, address1 ); 
  tmp2              = vec_ld( 16, address1 ); 
  tmp3              = vec_ld( 32, address1 ); 
  *Ox               = vec_perm(tmp1,tmp2,perm);     /*  Ox  Oy  Oz H1x */
  *H1y              = vec_perm(tmp2,tmp3,perm);     /* H1y H1z H2x H2y */
  *H2z              = vec_perm(tmp3,tmp3,perm);     /* H2z   -   -   - */

  /* just splat things... never mind that we fill all cells :-) */
  *Oy               = vec_splat(*Ox,1);
  *Oz               = vec_splat(*Ox,2);
  *H1x              = vec_splat(*Ox,3);
  *H1z              = vec_splat(*H1y,1);
  *H2x              = vec_splat(*H1y,2);
  *H2y              = vec_splat(*H1y,3);
}



static inline void update_i_water_forces(float *water,
					 float *shiftvec,
					 vector float Ox,
					 vector float Oy,
					 vector float Oz,
					 vector float H1x,
					 vector float H1y,
					 vector float H1z,
					 vector float H2x,
					 vector float H2y,
					 vector float H2z)
{
  vector float l1,l2,l3;
  vector unsigned char perm;
  vector unsigned int mask,ox00;
  vector float nul=vec_zero();
  
  ox00=(vector unsigned int)vec_splat_s32(0);
  /* load */
  perm              = vec_lvsr( 0, water ); 
  l1                = vec_ld(  0, water);
  l2                = vec_ld( 16, water);
  l3                = vec_ld( 32, water);
  mask              = vec_perm(ox00,
			       (vector unsigned int)vec_splat_s32(-1),perm);
  
  /* accumulate the forces */
  Ox               = vec_add(Ox,vec_sld(Ox,Ox,8));   /*  Ox  Ox' - - */
  Oy               = vec_add(Oy,vec_sld(Oy,Oy,8));   /*  Oy  Oy' - - */
  Oz               = vec_add(Oz,vec_sld(Oz,Oz,8));   /*  Oz  Oz' - - */
  H1x              = vec_add(H1x,vec_sld(H1x,H1x,8));  /* H1x H1x' - - */
  H1y              = vec_add(H1y,vec_sld(H1y,H1y,8));  /* H1y H1y' - - */
  H1z              = vec_add(H1z,vec_sld(H1z,H1z,8));  /* H1z H1z' - - */
  H2x              = vec_add(H2x,vec_sld(H2x,H2x,8));  /* H2x H2x' - - */
  H2y              = vec_add(H2y,vec_sld(H2y,H2y,8));  /* H2y H2y' - - */
  H2z              = vec_add(H2z,vec_sld(H2z,H2z,8));  /* H2z H2z' - - */

  Ox               = vec_mergeh(Ox,Oz);        /*  Ox  Oz  Ox'  Oz' */
  Oy               = vec_mergeh(Oy,H1x);       /*  Oy H1x  Oy' H1x' */
  H1y              = vec_mergeh(H1y,H2x);      /* H1y H2x H1y' H2x' */
  H1z              = vec_mergeh(H1z,H2y);      /* H1z H2y H1z' H2y' */
  H2z              = vec_mergeh(H2z,nul);       /* H2z  0  H2z'  0   */
  
  Ox          = vec_add(Ox,vec_sld(Ox,Ox,8));   /* Ox Oz - - */
  Oy          = vec_add(Oy,vec_sld(Oy,Oy,8));   /* Oy H1x - - */
  H1y         = vec_add(H1y,vec_sld(H1y,H1y,8));   /* H1y H2x - - */
  H1z         = vec_add(H1z,vec_sld(H1z,H1z,8));   /* H1z H2y - - */
  H2z         = vec_add(H2z,vec_sld(H2z,H2z,8));   /* H2z 0 ? 0 */

  Ox          = vec_mergeh(Ox,Oy);   /*  Ox  Oy  Oz H1x */
  H1y         = vec_mergeh(H1y,H1z); /* H1y H1z H2x H2y */
  H2z         = vec_mergeh(H2z,nul);  /* H2z  0   0   0  */

  Oy          = vec_sld(nul,Ox,12);  /* 0   Ox  Oy  Oz  */
  Oz          = vec_sld(Ox,H1y,8);   /* -  H1x H1y H1z  */
  H1x         = vec_sld(H1y,H2z,4);  /* -  H2x H2y H2z  */

  H2x         = vec_perm(nul,Ox,perm);   /* The part to add to l1 */
  H2y         = vec_perm(Ox,H1y,perm);  /* The part to add to l2 */
  H2z         = vec_perm(H1y,H2z,perm); /* The part to add to l3 */

  H2x         = vec_add(l1,H2x);
  H2y         = vec_add(l2,H2y);
  H2z         = vec_add(l3,H2z);  
   
  H2x         = vec_sel(l1,H2x,mask);
  mask        = vec_sld(ox00,mask,12);
  H2z         = vec_sel(H2z,l3,mask);
  
  /* store */
  vec_st(H2x,  0, water);
  vec_st(H2y, 16, water);
  vec_st(H2z, 32, water);
  
  Oy          = vec_add(Oy,Oz);
  Oy          = vec_add(Oy,H1x);     
  Oy          = vec_sld(Oy,nul,4);   /* x   y  z  0 */

  add_xyz_to_mem(shiftvec,Oy);
}

static inline void add_force_to_4_water(float *address1,
					float *address2,
					float *address3,
					float *address4,
					vector float Ox,
					vector float Oy,
					vector float Oz,
					vector float H1x,
					vector float H1y,
					vector float H1z,
					vector float H2x,
					vector float H2y,
					vector float H2z)
{
  vector float low,medium,high;
  vector unsigned char perm;
  vector unsigned int mask,oxFF,ox00;
  vector float tmp1,tmp2,tmp3,tmp4;
  vector float nul=vec_zero();

  oxFF=(vector unsigned int)vec_splat_s32(-1);
  ox00=(vector unsigned int)vec_splat_s32(0);

  tmp1              = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
  Ox               = vec_mergel(Ox,Oz);    /*  Oxc  Ozc  Oxd  Ozd */
  Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
  Oy               = vec_mergel(Oy,H1x);   /*  Oyc H1xc  Oyd H1xd */
  H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
  H1y              = vec_mergel(H1y,H2x);  /* H1yc H2xc H1yd H2xd */
  H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
  H1z              = vec_mergel(H1z,H2y);  /* H1zc H2yc H1zd H2yd */
  H2y              = vec_mergeh(H2z,nul);   /* H2za   0  H2zb   0  */
  H2z              = vec_mergel(H2z,nul);   /* H2zc   0  H2zd   0  */

  tmp2             = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
  Oz               = vec_mergel(tmp1,Oz);   /*  Oxb  Oyb  Ozb H1xb */
  tmp1             = vec_mergeh(Ox,Oy);    /*  Oxc  Oyc  Ozc H1xc */
  Ox               = vec_mergel(Ox,Oy);    /*  Oxd  Oyd  Ozd H1xd */
  Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
  H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
  H2x              = vec_mergeh(H1y,H1z);  /* H1yc H1zc H2xc H2yc */
  H1y              = vec_mergel(H1y,H1z);  /* H1yd H1zd H2xd H2yd */
  H1z              = vec_mergeh(H2y,nul);   /* H2za   0    0    0  */
  H2y              = vec_mergel(H2y,nul);   /* H2zb   0    0    0  */
  tmp3             = vec_mergeh(H2z,nul);   /* H2zc   0    0    0  */
  H2z              = vec_mergel(H2z,nul);   /* H2zd   0    0    0  */
  
  /* move into position, load and add */  
  perm             = vec_lvsr( 0, (int *) address1 ); 
  low              = vec_ld(  0, address1);
  medium           = vec_ld( 16, address1);
  high             = vec_ld( 32, address1);
  mask             = vec_perm(ox00,oxFF,perm);
  tmp4             = vec_add(vec_perm(nul,tmp2,perm),low);
  tmp2             = vec_add(vec_perm(tmp2,Oy,perm),medium);
  Oy               = vec_add(vec_perm(Oy,H1z,perm),high);
  vec_st(vec_sel(low,tmp4,mask),  0, address1);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp2, 16, address1);
  vec_st(vec_sel(Oy,high,mask), 32, address1);

  perm             = vec_lvsr( 0, (int *) address2 ); 
  low              = vec_ld(  0, address2);
  medium           = vec_ld( 16, address2);
  high             = vec_ld( 32, address2);
  mask             = vec_perm(ox00,oxFF,perm);
  H1z              = vec_add(vec_perm(nul,Oz,perm),low);
  Oz               = vec_add(vec_perm(Oz,H1x,perm),medium);
  H1x              = vec_add(vec_perm(H1x,H2y,perm),high);
  vec_st(vec_sel(low,H1z,mask),  0, address2);
  mask        = vec_sld(ox00,mask,12);
  vec_st(Oz, 16, address2);
  vec_st(vec_sel(H1x,high,mask), 32, address2);

  perm             = vec_lvsr( 0, (int *) address3 ); 
  low              = vec_ld(  0, address3);
  medium           = vec_ld( 16, address3);
  high             = vec_ld( 32, address3);
  mask             = vec_perm(ox00,oxFF,perm);
  H2y              = vec_add(vec_perm(nul,tmp1,perm),low);
  tmp1             = vec_add(vec_perm(tmp1,H2x,perm),medium);
  H2x              = vec_add(vec_perm(H2x,tmp3,perm),high);
  vec_st(vec_sel(low,H2y,mask),  0, address3);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp1, 16, address3);
  vec_st(vec_sel(H2x,high,mask), 32, address3);

  perm             = vec_lvsr( 0, (int *) address4 ); 
  low              = vec_ld(  0, address4);
  medium           = vec_ld( 16, address4);
  high             = vec_ld( 32, address4);
  mask             = vec_perm(ox00,oxFF,perm);
  H2x              = vec_add(vec_perm(nul,Ox,perm),low);
  Ox               = vec_add(vec_perm(Ox,H1y,perm),medium);
  H1y              = vec_add(vec_perm(H1y,H2z,perm),high);
  vec_st(vec_sel(low,H2x,mask),  0, address4);
  mask        = vec_sld(ox00,mask,12);
  vec_st(Ox, 16, address4);
  vec_st(vec_sel(H1y,high,mask), 32, address4);
}



static inline void add_force_to_3_water(float *address1,
					float *address2,
					float *address3,
					vector float Ox,
					vector float Oy,
					vector float Oz,
					vector float H1x,
					vector float H1y,
					vector float H1z,
					vector float H2x,
					vector float H2y,
					vector float H2z)
{
  vector float low,medium,high;
  vector unsigned char perm;
  vector unsigned int mask,oxFF,ox00;
  vector float tmp1,tmp2,tmp3;

  vector float nul=vec_zero();
  oxFF=(vector unsigned int)vec_splat_s32(-1);
  ox00=(vector unsigned int)vec_splat_s32(0);
  
  tmp1              = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
  Ox               = vec_mergel(Ox,Oz);    /*  Oxc  Ozc   ?    ?  */
  Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
  Oy               = vec_mergel(Oy,H1x);   /*  Oyc H1xc   ?    ?  */
  H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
  H1y              = vec_mergel(H1y,H2x);  /* H1yc H2xc   ?    ?  */
  H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
  H1z              = vec_mergel(H1z,H2y);  /* H1zc H2yc   ?    ?  */
  H2y              = vec_mergeh(H2z,nul);   /* H2za   0  H2zb   0  */
  H2z              = vec_mergel(H2z,nul);   /* H2zc   0    ?    0  */

  tmp2              = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
  Oz               = vec_mergel(tmp1,Oz);   /*  Oxb  Oyb  Ozb H1xb */
  tmp1              = vec_mergeh(Ox,Oy);    /*  Oxc  Oyc  Ozc H1xc */
  Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
  H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
  H2x              = vec_mergeh(H1y,H1z);  /* H1yc H1zc H2xc H2yc */
  H1z              = vec_mergeh(H2y,nul);   /* H2za   0    0    0  */
  H2y              = vec_mergel(H2y,nul);   /* H2zb   0    0    0  */
  tmp3              = vec_mergeh(H2z,nul);   /* H2zc   0    0    0  */
 
  /* move into position, load and add */  
  perm             = vec_lvsr( 0, (int *) address1 ); 
  low              = vec_ld(  0, address1);
  medium           = vec_ld( 16, address1);
  high             = vec_ld( 32, address1);
  mask             = vec_perm(ox00,oxFF,perm);
  H1y              = vec_add(vec_perm(nul,tmp2,perm),low);
  tmp2             = vec_add(vec_perm(tmp2,Oy,perm),medium);
  Oy               = vec_add(vec_perm(Oy,H1z,perm),high);
  vec_st(vec_sel(low,H1y,mask),  0, address1);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp2, 16, address1);
  vec_st(vec_sel(Oy,high,mask), 32, address1);

  perm             = vec_lvsr( 0, (int *) address2 ); 
  low              = vec_ld(  0, address2);
  medium           = vec_ld( 16, address2);
  high             = vec_ld( 32, address2);
  mask             = vec_perm(ox00,oxFF,perm);
  H1z              = vec_add(vec_perm(nul,Oz,perm),low);
  Oz               = vec_add(vec_perm(Oz,H1x,perm),medium);
  H1x              = vec_add(vec_perm(H1x,H2y,perm),high);
  vec_st(vec_sel(low,H1z,mask),  0, address2);
  mask        = vec_sld(ox00,mask,12);
  vec_st(Oz, 16, address2);
  vec_st(vec_sel(H1x,high,mask), 32, address2);

  perm             = vec_lvsr( 0, (int *) address3 ); 
  low              = vec_ld(  0, address3);
  medium           = vec_ld( 16, address3);
  high             = vec_ld( 32, address3);
  mask             = vec_perm(ox00,oxFF,perm);
  H2y              = vec_add(vec_perm(nul,tmp1,perm),low);
  tmp1             = vec_add(vec_perm(tmp1,H2x,perm),medium);
  H2x              = vec_add(vec_perm(H2x,tmp3,perm),high);
  vec_st(vec_sel(low,H2y,mask),  0, address3);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp1, 16, address3);
  vec_st(vec_sel(H2x,high,mask), 32, address3);
}




static inline void add_force_to_2_water(float *address1,
					float *address2,
					vector float Ox,
					vector float Oy,
					vector float Oz,
					vector float H1x,
					vector float H1y,
					vector float H1z,
					vector float H2x,
					vector float H2y,
					vector float H2z)
{
  vector float low,medium,high;
  vector unsigned char perm;
  vector unsigned int mask,oxFF,ox00;
  vector float tmp1,tmp2,tmp3;

  vector float nul=vec_zero();
  oxFF=(vector unsigned int)vec_splat_s32(-1);
  ox00=(vector unsigned int)vec_splat_s32(0);

  tmp1              = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
  Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
  H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
  H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
  H2y              = vec_mergeh(H2z,nul);   /* H2za   0  H2zb   0  */

  tmp2              = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
  Oz               = vec_mergel(tmp1,Oz);   /*  Oxb  Oyb  Ozb H1xb */
  Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
  H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
  H1z              = vec_mergeh(H2y,nul);   /* H2za   0    0    0  */
  H2y              = vec_mergel(H2y,nul);   /* H2zb   0    0    0  */
  
  /* move into position and add */
  perm             = vec_lvsr( 0, (int *) address1 ); 
  low              = vec_ld(  0, address1);
  medium           = vec_ld( 16, address1);
  high             = vec_ld( 32, address1);
  mask             = vec_perm(ox00,oxFF,perm);
  H2x              = vec_add(vec_perm(nul,tmp2,perm),low);
  tmp2             = vec_add(vec_perm(tmp2,Oy,perm),medium);
  Oy               = vec_add(vec_perm(Oy,H1z,perm),high);
  vec_st(vec_sel(low,H2x,mask),  0, address1);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp2, 16, address1);
  vec_st(vec_sel(Oy,high,mask), 32, address1);

  perm             = vec_lvsr( 0, (int *) address2 ); 
  low              = vec_ld(  0, address2);
  medium           = vec_ld( 16, address2);
  high             = vec_ld( 32, address2);
  mask             = vec_perm(ox00,oxFF,perm);
  H1z              = vec_add(vec_perm(nul,Oz,perm),low);
  Oz               = vec_add(vec_perm(Oz,H1x,perm),medium);
  H1x              = vec_add(vec_perm(H1x,H2y,perm),high);
  vec_st(vec_sel(low,H1z,mask),  0, address2);
  mask        = vec_sld(ox00,mask,12);
  vec_st(Oz, 16, address2);
  vec_st(vec_sel(H1x,high,mask), 32, address2);
}



static inline void add_force_to_1_water(float *address1,
					vector float Ox,
					vector float Oy,
					vector float Oz,
					vector float H1x,
					vector float H1y,
					vector float H1z,
					vector float H2x,
					vector float H2y,
					vector float H2z)
{
  vector float low,medium,high;
   vector unsigned char perm;
  vector unsigned int mask,oxFF,ox00;
  vector float tmp1,tmp2,tmp3;

  vector float nul=vec_zero();
  oxFF=(vector unsigned int)vec_splat_s32(-1);
  ox00=(vector unsigned int)vec_splat_s32(0);

  /* load */
  perm              = vec_lvsr( 0, (int *) address1 ); 
  low               = vec_ld(  0, address1);
  medium            = vec_ld( 16, address1);
  high              = vec_ld( 32, address1);

  tmp1              = vec_mergeh(Ox,Oz);    /*  Oxa  Oza   ?    ?  */
  Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa   ?    ?  */
  H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa   ?    ?  */
  H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya   ?    ?  */
  H2y              = vec_mergeh(H2z,nul);   /* H2za   0    ?    0  */

  tmp2              = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
  Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
  H1z              = vec_mergeh(H2y,nul);   /* H2za   0    0    0  */
  
  /* move into position and add */
  perm             = vec_lvsr( 0, (int *) address1 ); 
  low              = vec_ld(  0, address1);
  medium           = vec_ld( 16, address1);
  high             = vec_ld( 32, address1);
  mask             = vec_perm(ox00,oxFF,perm);
  H2x              = vec_add(vec_perm(nul,tmp2,perm),low);
  tmp2             = vec_add(vec_perm(tmp2,Oy,perm),medium);
  Oy               = vec_add(vec_perm(Oy,H1z,perm),high);
  vec_st(vec_sel(low,H2x,mask),  0, address1);
  mask        = vec_sld(ox00,mask,12);
  vec_st(tmp2, 16, address1);
  vec_st(vec_sel(Oy,high,mask), 32, address1);
}



static inline void do_4_ctable_coul(float *VFtab,
				   vector float rtab,
				   vector float *VV,
				   vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
  int idx1,idx2,idx3,idx4;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);
  idx4     = *(((int *)&vidx)+3);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
    tab3   = vec_ld( 0, VFtab+idx3);
    G      = vec_ld(16, VFtab+idx3);
    tab3   = vec_sld(tab3,G,8);
    tab4   = vec_ld( 0, VFtab+idx4);
    H      = vec_ld(16, VFtab+idx4);
    tab4   = vec_sld(tab4,H,8);
  } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
    tab3=vec_ld(0, VFtab+idx3);
    tab4=vec_ld(0, VFtab+idx4);
  }  
  
  /* table data is aligned */
  transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);
  
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/* do only coulomb, but on a table with both coulomb and lj data */
static inline void do_4_ljctable_coul(float *VFtab,
				      vector float rtab,
				      vector float *VV,
				      vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
  int idx1,idx2,idx3,idx4;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);
  idx4     = *(((int *)&vidx)+3);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
    tab3   = vec_ld( 0, VFtab+idx3);
    G      = vec_ld(16, VFtab+idx3);
    tab3   = vec_sld(tab3,G,8);
    tab4   = vec_ld( 0, VFtab+idx4);
    H      = vec_ld(16, VFtab+idx4);
    tab4   = vec_sld(tab4,H,8);
  } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
    tab3=vec_ld(0, VFtab+idx3);
    tab4=vec_ld(0, VFtab+idx4);
  }  

  transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_4_ljtable_lj(float *VFtab,
				   vector float rtab,
				   vector float *VVdisp,
				   vector float *FFdisp,
				   vector float *VVrep,
				   vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  vector float tabd1,tabd2,tabd3,tabd4;
  vector float tabr1,tabr2,tabr3,tabr4;

  int    idx1,idx2,idx3,idx4;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(3));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);
  idx4     = *(((int *)&vidx)+3);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { 
    /* not 16 byte aligned, i.e. must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabd1  =  vec_ld( 0, VFtab+idx1);
    tabr1  =  vec_ld(16, VFtab+idx1);
    Y      =  vec_ld(32, VFtab+idx1);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabd2  =  vec_ld( 0, VFtab+idx2);
    tabr2  =  vec_ld(16, VFtab+idx2);
    F      =  vec_ld(32, VFtab+idx2);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);

    tabd3  =  vec_ld( 0, VFtab+idx3);
    tabr3  =  vec_ld(16, VFtab+idx3);
    G      =  vec_ld(32, VFtab+idx3);
    tabd3  =  vec_sld(tabd3,tabr3,8);
    tabr3  =  vec_sld(tabr3,G,8);

    tabd4  =  vec_ld( 0, VFtab+idx4);
    tabr4  =  vec_ld(16, VFtab+idx4);
    H      =  vec_ld(32, VFtab+idx4);
    tabd4  =  vec_sld(tabd4,tabr4,8);
    tabr4  =  vec_sld(tabr4,H,8);
  } else { /* 16 byte aligned */
    tabd1  = vec_ld( 0, VFtab+idx1);
    tabr1  = vec_ld(16, VFtab+idx1);
    tabd2  = vec_ld( 0, VFtab+idx2);
    tabr2  = vec_ld(16, VFtab+idx2);
    tabd3  = vec_ld( 0, VFtab+idx3);
    tabr3  = vec_ld(16, VFtab+idx3);
    tabd4  = vec_ld( 0, VFtab+idx4);
    tabr4  = vec_ld(16, VFtab+idx4);
  }
    
  transpose_4_to_4(tabd1,tabd2,tabd3,tabd4,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_4_to_4(tabr1,tabr2,tabr3,tabr4,&Y,&F,&G,&H);
 
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_4_ljctable_coul_and_lj(float *VFtab,
					     vector float rtab,
					     vector float *VVcoul,
					     vector float *FFcoul,
					     vector float *VVdisp,
					     vector float *FFdisp,
					     vector float *VVrep,
					     vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
  vector float tabd3,tabr3,tabc3,tabd4,tabr4,tabc4;
  int idx1,idx2,idx3,idx4;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);
  idx4     = *(((int *)&vidx)+3);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) {
    /* not 16-byte aligned, but must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabc1  =  vec_ld( 0, VFtab+idx1);
    tabd1  =  vec_ld(16, VFtab+idx1);
    tabr1  =  vec_ld(32, VFtab+idx1);
    Y      =  vec_ld(48, VFtab+idx1);
    tabc1  =  vec_sld(tabc1,tabd1,8);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabc2  =  vec_ld( 0, VFtab+idx2);
    tabd2  =  vec_ld(16, VFtab+idx2);
    tabr2  =  vec_ld(32, VFtab+idx2);
    F      =  vec_ld(48, VFtab+idx2);
    tabc2  =  vec_sld(tabc2,tabd2,8);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);

    tabc3  =  vec_ld( 0, VFtab+idx3);
    tabd3  =  vec_ld(16, VFtab+idx3);
    tabr3  =  vec_ld(32, VFtab+idx3);
    G      =  vec_ld(48, VFtab+idx3);
    tabc3  =  vec_sld(tabc3,tabd3,8);
    tabd3  =  vec_sld(tabd3,tabr3,8);
    tabr3  =  vec_sld(tabr3,G,8);

    tabc4  =  vec_ld( 0, VFtab+idx4);
    tabd4  =  vec_ld(16, VFtab+idx4);
    tabr4  =  vec_ld(32, VFtab+idx4);
    H      =  vec_ld(48, VFtab+idx4);
    tabc4  =  vec_sld(tabc4,tabd4,8);
    tabd4  =  vec_sld(tabd4,tabr4,8);
    tabr4  =  vec_sld(tabr4,H,8);
  } else { /* 16 byte aligned */
    tabc1  = vec_ld( 0, VFtab+idx1);
    tabd1  = vec_ld(16, VFtab+idx1);
    tabr1  = vec_ld(32, VFtab+idx1);
    tabc2  = vec_ld( 0, VFtab+idx2);
    tabd2  = vec_ld(16, VFtab+idx2);
    tabr2  = vec_ld(32, VFtab+idx2);
    tabc3  = vec_ld( 0, VFtab+idx3);
    tabd3  = vec_ld(16, VFtab+idx3);
    tabr3  = vec_ld(32, VFtab+idx3);
    tabc4  = vec_ld( 0, VFtab+idx4);
    tabd4  = vec_ld(16, VFtab+idx4);
    tabr4  = vec_ld(32, VFtab+idx4);
  }
  transpose_4_to_4(tabc1,tabc2,tabc3,tabc4,&Y,&F,&G,&H);
   
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFcoul   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

 transpose_4_to_4(tabd1,tabd2,tabd3,tabd4,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

 transpose_4_to_4(tabr1,tabr2,tabr3,tabr4,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_3_ctable_coul(float *VFtab,
				    vector float rtab,
				    vector float *VV,
				    vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
  int idx1,idx2,idx3;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
    tab3   = vec_ld( 0, VFtab+idx3);
    G      = vec_ld(16, VFtab+idx3);
    tab3   = vec_sld(tab3,G,8);
   } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
    tab3=vec_ld(0, VFtab+idx3);
   }  

  /* table data is aligned */
  transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);
  

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/* do only coulomb, but on a table with both coulomb and lj data */
static inline void do_3_ljctable_coul(float *VFtab,
				      vector float rtab,
				      vector float *VV,
				      vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
  int idx1,idx2,idx3;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());


  if(((unsigned int)VFtab)%16) { 
    /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
    tab3   = vec_ld( 0, VFtab+idx3);
    G      = vec_ld(16, VFtab+idx3);
    tab3   = vec_sld(tab3,G,8);
  } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
    tab3=vec_ld(0, VFtab+idx3);
  }  
  
  /* table data is aligned */
  transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_3_ljtable_lj(float *VFtab,
				   vector float rtab,
				   vector float *VVdisp,
				   vector float *FFdisp,
				   vector float *VVrep,
				   vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  int    idx1,idx2,idx3;
  vector float tabd1,tabd2,tabd3;
  vector float tabr1,tabr2,tabr3;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(3));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { 
    /* not 16 byte aligned, i.e. must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabd1  =  vec_ld( 0, VFtab+idx1);
    tabr1  =  vec_ld(16, VFtab+idx1);
    Y      =  vec_ld(32, VFtab+idx1);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabd2  =  vec_ld( 0, VFtab+idx2);
    tabr2  =  vec_ld(16, VFtab+idx2);
    F      =  vec_ld(32, VFtab+idx2);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);

    tabd3  =  vec_ld( 0, VFtab+idx3);
    tabr3  =  vec_ld(16, VFtab+idx3);
    G      =  vec_ld(32, VFtab+idx3);
    tabd3  =  vec_sld(tabd3,tabr3,8);
    tabr3  =  vec_sld(tabr3,G,8);
  } else { /* 16 byte aligned */
    tabd1  = vec_ld( 0, VFtab+idx1);
    tabr1  = vec_ld(16, VFtab+idx1);
    tabd2  = vec_ld( 0, VFtab+idx2);
    tabr2  = vec_ld(16, VFtab+idx2);
    tabd3  = vec_ld( 0, VFtab+idx3);
    tabr3  = vec_ld(16, VFtab+idx3);
  }
    
  transpose_3_to_4(tabd1,tabd2,tabd3,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_3_to_4(tabr1,tabr2,tabr3,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_3_ljctable_coul_and_lj(float *VFtab,
					     vector float rtab,
					     vector float *VVcoul,
					     vector float *FFcoul,
					     vector float *VVdisp,
					     vector float *FFdisp,
					     vector float *VVrep,
					     vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
  vector float tabd3,tabr3,tabc3;
  int idx1,idx2,idx3;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  idx3     = *(((int *)&vidx)+2);
  
  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) {
    /* not 16-byte aligned, but must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabc1  =  vec_ld( 0, VFtab+idx1);
    tabd1  =  vec_ld(16, VFtab+idx1);
    tabr1  =  vec_ld(32, VFtab+idx1);
    Y      =  vec_ld(48, VFtab+idx1);
    tabc1  =  vec_sld(tabc1,tabd1,8);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabc2  =  vec_ld( 0, VFtab+idx2);
    tabd2  =  vec_ld(16, VFtab+idx2);
    tabr2  =  vec_ld(32, VFtab+idx2);
    F      =  vec_ld(48, VFtab+idx2);
    tabc2  =  vec_sld(tabc2,tabd2,8);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);

    tabc3  =  vec_ld( 0, VFtab+idx3);
    tabd3  =  vec_ld(16, VFtab+idx3);
    tabr3  =  vec_ld(32, VFtab+idx3);
    G      =  vec_ld(48, VFtab+idx3);
    tabc3  =  vec_sld(tabc3,tabd3,8);
    tabd3  =  vec_sld(tabd3,tabr3,8);
    tabr3  =  vec_sld(tabr3,G,8);  
  } else { /* 16 byte aligned */
    tabc1  = vec_ld( 0, VFtab+idx1);
    tabd1  = vec_ld(16, VFtab+idx1);
    tabr1  = vec_ld(32, VFtab+idx1);
    tabc2  = vec_ld( 0, VFtab+idx2);
    tabd2  = vec_ld(16, VFtab+idx2);
    tabr2  = vec_ld(32, VFtab+idx2);
    tabc3  = vec_ld( 0, VFtab+idx3);
    tabd3  = vec_ld(16, VFtab+idx3);
    tabr3  = vec_ld(32, VFtab+idx3);
  }
  transpose_3_to_4(tabc1,tabc2,tabc3,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFcoul   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_3_to_4(tabd1,tabd2,tabd3,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_3_to_4(tabr1,tabr2,tabr3,&Y,&F,&G,&H);
 
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}




static inline void do_2_ctable_coul(float *VFtab,
				    vector float rtab,
				    vector float *VV,
				    vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2;
  int idx1,idx2;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
   } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
   }  

  /* table data is aligned */
  transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/* do only coulomb, but on a table with both coulomb and lj data */
static inline void do_2_ljctable_coul(float *VFtab,
				      vector float rtab,
				      vector float *VV,
				      vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1,tab2;
  int idx1,idx2;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());


  if(((unsigned int)VFtab)%16) { 
    /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
    tab2   = vec_ld( 0, VFtab+idx2);
    F      = vec_ld(16, VFtab+idx2);
    tab2   = vec_sld(tab2,F,8);
  } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
    tab2=vec_ld(0, VFtab+idx2);
  }  
  
  /* table data is aligned */
  transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_2_ljtable_lj(float *VFtab,
				   vector float rtab,
				   vector float *VVdisp,
				   vector float *FFdisp,
				   vector float *VVrep,
				   vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  int    idx1,idx2;
  vector float tabd1,tabd2;
  vector float tabr1,tabr2;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(3));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { 
    /* not 16 byte aligned, i.e. must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabd1  =  vec_ld( 0, VFtab+idx1);
    tabr1  =  vec_ld(16, VFtab+idx1);
    Y      =  vec_ld(32, VFtab+idx1);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabd2  =  vec_ld( 0, VFtab+idx2);
    tabr2  =  vec_ld(16, VFtab+idx2);
    F      =  vec_ld(32, VFtab+idx2);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);
  } else { /* 16 byte aligned */
    tabd1  = vec_ld( 0, VFtab+idx1);
    tabr1  = vec_ld(16, VFtab+idx1);
    tabd2  = vec_ld( 0, VFtab+idx2);
    tabr2  = vec_ld(16, VFtab+idx2);
  }
    
  transpose_2_to_4(tabd1,tabd2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_2_to_4(tabr1,tabr2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_2_ljctable_coul_and_lj(float *VFtab,
					     vector float rtab,
					     vector float *VVcoul,
					     vector float *FFcoul,
					     vector float *VVdisp,
					     vector float *FFdisp,
					     vector float *VVrep,
					     vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
  int idx1,idx2;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  idx2     = *(((int *)&vidx)+1);
  
  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) {
    /* not 16-byte aligned, but must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabc1  =  vec_ld( 0, VFtab+idx1);
    tabd1  =  vec_ld(16, VFtab+idx1);
    tabr1  =  vec_ld(32, VFtab+idx1);
    Y      =  vec_ld(48, VFtab+idx1);
    tabc1  =  vec_sld(tabc1,tabd1,8);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);

    tabc2  =  vec_ld( 0, VFtab+idx2);
    tabd2  =  vec_ld(16, VFtab+idx2);
    tabr2  =  vec_ld(32, VFtab+idx2);
    F      =  vec_ld(48, VFtab+idx2);
    tabc2  =  vec_sld(tabc2,tabd2,8);
    tabd2  =  vec_sld(tabd2,tabr2,8);
    tabr2  =  vec_sld(tabr2,F,8);
  } else { /* 16 byte aligned */
    tabc1  = vec_ld( 0, VFtab+idx1);
    tabd1  = vec_ld(16, VFtab+idx1);
    tabr1  = vec_ld(32, VFtab+idx1);
    tabc2  = vec_ld( 0, VFtab+idx2);
    tabd2  = vec_ld(16, VFtab+idx2);
    tabr2  = vec_ld(32, VFtab+idx2);
  }
  transpose_2_to_4(tabc1,tabc2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFcoul   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_2_to_4(tabd1,tabd2,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_2_to_4(tabr1,tabr2,&Y,&F,&G,&H);
 
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}



static inline void do_1_ctable_coul(float *VFtab,
				    vector float rtab,
				    vector float *VV,
				    vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1;
  int idx1;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
   } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
   }  

  /* table data is aligned */
  transpose_1_to_4(tab1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/* do only coulomb, but on a table with both coulomb and lj data */
static inline void do_1_ljctable_coul(float *VFtab,
				      vector float rtab,
				      vector float *VV,
				      vector float *FF)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2,tab1;
  int idx1;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());


  if(((unsigned int)VFtab)%16) { 
    /* not 16-byte aligned, but must be 8 byte. */
    tab1   = vec_ld( 0, VFtab+idx1);
    Y      = vec_ld(16, VFtab+idx1);
    tab1   = vec_sld(tab1,Y,8);
  } else { /* aligned */
    tab1=vec_ld(0, VFtab+idx1);
  }  
  
  /* table data is aligned */
  transpose_1_to_4(tab1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_1_ljtable_lj(float *VFtab,
				   vector float rtab,
				   vector float *VVdisp,
				   vector float *FFdisp,
				   vector float *VVrep,
				   vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  int    idx1;
  vector float tabd1;
  vector float tabr1;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_sl(vidx,vec_splat_u32(3));

  idx1     = *((int *)&vidx);

  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) { 
    /* not 16 byte aligned, i.e. must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabd1  =  vec_ld( 0, VFtab+idx1);
    tabr1  =  vec_ld(16, VFtab+idx1);
    Y      =  vec_ld(32, VFtab+idx1);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);
  } else { /* 16 byte aligned */
    tabd1  = vec_ld( 0, VFtab+idx1);
    tabr1  = vec_ld(16, VFtab+idx1);
  }
    
  transpose_1_to_4(tabd1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_1_to_4(tabr1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


static inline void do_1_ljctable_coul_and_lj(float *VFtab,
					     vector float rtab,
					     vector float *VVcoul,
					     vector float *FFcoul,
					     vector float *VVdisp,
					     vector float *FFdisp,
					     vector float *VVrep,
					     vector float *FFrep)
{
  vector signed int vidx;
  vector float Y,F,G,H,eps,eps2;
  vector float tabd1,tabr1,tabc1;
  int idx1;

  vidx     = vec_cts(rtab,0); 
  vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
  vidx     = vec_sl(vidx,vec_splat_u32(2));

  idx1     = *((int *)&vidx);
  
  eps      = vec_sub(rtab,vec_floor(rtab));
  eps2     = vec_madd(eps,eps,vec_zero());

  if(((unsigned int)VFtab)%16) {
    /* not 16-byte aligned, but must be 8 byte. */
    /* use Y,F,G,H as temp storage */
    tabc1  =  vec_ld( 0, VFtab+idx1);
    tabd1  =  vec_ld(16, VFtab+idx1);
    tabr1  =  vec_ld(32, VFtab+idx1);
    Y      =  vec_ld(48, VFtab+idx1);
    tabc1  =  vec_sld(tabc1,tabd1,8);
    tabd1  =  vec_sld(tabd1,tabr1,8);
    tabr1  =  vec_sld(tabr1,Y,8);
  } else { /* 16 byte aligned */
    tabc1  = vec_ld( 0, VFtab+idx1);
    tabd1  = vec_ld(16, VFtab+idx1);
    tabr1  = vec_ld(32, VFtab+idx1);
  }
  transpose_1_to_4(tabc1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFcoul   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_1_to_4(tabd1,&Y,&F,&G,&H);

  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFdisp   = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */

  transpose_1_to_4(tabr1,&Y,&F,&G,&H);
 
  F         = vec_madd(G,eps,F);           /* F + Geps   */
  H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
  F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
  *VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
  F         = vec_madd(G,eps,F);           /* Fp + Geps */
  *FFrep    = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}







static inline vector float do_recip(vector float rsq)
{
  vector float tmp,lu;

  lu        = vec_re(rsq);
  tmp       = vec_nmsub(rsq,lu,vec_two());
  return vec_madd(lu,tmp,vec_zero());
}

static inline vector float do_invsqrt(vector float rsq)
{
  vector float lu,tmpA,tmpB;

  lu        = vec_rsqrte(rsq);
  tmpA      = vec_madd(lu,lu,vec_zero());
  tmpB      = vec_madd(lu,vec_half(),vec_zero());
  tmpA      = vec_nmsub(rsq,tmpA,vec_one());
  return vec_madd(tmpA,tmpB,lu);
}

static inline void do_3_invsqrt(vector float rsq1,
				vector float rsq2,
				vector float rsq3,
				vector float *rinvsq1,
				vector float *rinvsq2,
				vector float *rinvsq3)
{
  vector float lu1,lu2,lu3;
  vector float tmpA1,tmpA2,tmpA3,tmpB1,tmpB2,tmpB3;
  vector float nul=vec_zero();
  vector float half=vec_half();
  vector float one=vec_one();

  lu1       = vec_rsqrte(rsq1);
  lu2       = vec_rsqrte(rsq2);
  lu3       = vec_rsqrte(rsq3);
  tmpA1     = vec_madd(lu1,lu1,nul);
  tmpA2     = vec_madd(lu2,lu2,nul);
  tmpA3     = vec_madd(lu3,lu3,nul);
  tmpB1     = vec_madd(lu1,half,nul);
  tmpB2     = vec_madd(lu2,half,nul);
  tmpB3     = vec_madd(lu3,half,nul);
  tmpA1     = vec_nmsub(rsq1,tmpA1,one);
  tmpA2     = vec_nmsub(rsq2,tmpA2,one);
  tmpA3     = vec_nmsub(rsq3,tmpA3,one);
  *rinvsq1  = vec_madd(tmpA1,tmpB1,lu1);
  *rinvsq2  = vec_madd(tmpA2,tmpB2,lu2);
  *rinvsq3  = vec_madd(tmpA3,tmpB3,lu3);
}

static inline void do_9_invsqrt(vector float rsq1,
				vector float rsq2,
				vector float rsq3,
				vector float rsq4,
				vector float rsq5,
				vector float rsq6,
				vector float rsq7,
				vector float rsq8,
				vector float rsq9,
				vector float *rinv1,
				vector float *rinv2,
				vector float *rinv3,
				vector float *rinv4,
				vector float *rinv5,
				vector float *rinv6,
				vector float *rinv7,
				vector float *rinv8,
				vector float *rinv9)		
{
  vector float lu1,lu2,lu3,lu4,lu5,lu6,lu7,lu8,lu9;
  vector float tmpA1,tmpA2,tmpA3,tmpA4,tmpA5,tmpA6,tmpA7,tmpA8,tmpA9;
  vector float tmpB1,tmpB2,tmpB3,tmpB4,tmpB5,tmpB6,tmpB7,tmpB8,tmpB9;
  vector float nul=vec_zero();
  vector float half=vec_half();
  vector float one=vec_one();

  lu1       = vec_rsqrte(rsq1);
  lu2       = vec_rsqrte(rsq2);
  lu3       = vec_rsqrte(rsq3);
  lu4       = vec_rsqrte(rsq4);
  lu5       = vec_rsqrte(rsq5);
  lu6       = vec_rsqrte(rsq6);
  lu7       = vec_rsqrte(rsq7);
  lu8       = vec_rsqrte(rsq8);
  lu9       = vec_rsqrte(rsq9);
  tmpA1     = vec_madd(lu1,lu1,nul);
  tmpA2     = vec_madd(lu2,lu2,nul);
  tmpA3     = vec_madd(lu3,lu3,nul);
  tmpA4     = vec_madd(lu4,lu4,nul);
  tmpA5     = vec_madd(lu5,lu5,nul);
  tmpA6     = vec_madd(lu6,lu6,nul);
  tmpA7     = vec_madd(lu7,lu7,nul);
  tmpA8     = vec_madd(lu8,lu8,nul);
  tmpA9     = vec_madd(lu9,lu9,nul);
  tmpB1     = vec_madd(lu1,half,nul);
  tmpB2     = vec_madd(lu2,half,nul);
  tmpB3     = vec_madd(lu3,half,nul);
  tmpB4     = vec_madd(lu4,half,nul);
  tmpB5     = vec_madd(lu5,half,nul);
  tmpB6     = vec_madd(lu6,half,nul);
  tmpB7     = vec_madd(lu7,half,nul);
  tmpB8     = vec_madd(lu8,half,nul);
  tmpB9     = vec_madd(lu9,half,nul);
  tmpA1     = vec_nmsub(rsq1,tmpA1,one);
  tmpA2     = vec_nmsub(rsq2,tmpA2,one);
  tmpA3     = vec_nmsub(rsq3,tmpA3,one);
  tmpA4     = vec_nmsub(rsq4,tmpA4,one);
  tmpA5     = vec_nmsub(rsq5,tmpA5,one);
  tmpA6     = vec_nmsub(rsq6,tmpA6,one);
  tmpA7     = vec_nmsub(rsq7,tmpA7,one);
  tmpA8     = vec_nmsub(rsq8,tmpA8,one);
  tmpA9     = vec_nmsub(rsq9,tmpA9,one);
  *rinv1  = vec_madd(tmpA1,tmpB1,lu1);
  *rinv2  = vec_madd(tmpA2,tmpB2,lu2);
  *rinv3  = vec_madd(tmpA3,tmpB3,lu3);
  *rinv4  = vec_madd(tmpA4,tmpB4,lu4);
  *rinv5  = vec_madd(tmpA5,tmpB5,lu5);
  *rinv6  = vec_madd(tmpA6,tmpB6,lu6);
  *rinv7  = vec_madd(tmpA7,tmpB7,lu7);
  *rinv8  = vec_madd(tmpA8,tmpB8,lu8);
  *rinv9  = vec_madd(tmpA9,tmpB9,lu9);
  /* 36 invsqrt in about 48 cycles due to pipelining ... pretty fast :-) */
}

void inl0100_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],int type[],int ntype,float nbfp[],
		     float Vnb[]);
void inl0300_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],int type[],int ntype,float nbfp[],
		     float Vnb[],float tabscale,float VFtab[]);
void inl1000_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[]);
void inl1020_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[]);
void inl1030_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[]);
void inl1100_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[]);
void inl1120_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[]);
void inl1130_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[]);
void inl2000_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf);
void inl2020_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf);
void inl2030_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf);
void inl2100_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf, int type[],int ntype,
		     float nbfp[],float Vnb[]);
void inl2120_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf, int type[],int ntype,
		     float nbfp[],float Vnb[]);
void inl2130_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float krf, float crf, int type[],int ntype,
		     float nbfp[],float Vnb[]);
void inl3000_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float tabscale,float VFtab[]); 
void inl3020_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float tabscale,float VFtab[]);
void inl3030_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     float tabscale,float VFtab[]);
void inl3100_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale, float VFtab[]);
void inl3120_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale, float VFtab[]);
void inl3130_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale, float VFtab[]);
void inl3300_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale,float VFtab[]);
void inl3320_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale,float VFtab[]);
void inl3330_altivec(int nri,int iinr[],int jindex[],int jjnr[],int shift[],
		     float shiftvec[],float fshift[],int gid[],float pos[],
		     float faction[],float charge[],float facel,float Vc[],
		     int type[],int ntype,float nbfp[],float Vnb[],
		     float tabscale,float VFtab[]);

#endif /* ppc_altivec.h */
