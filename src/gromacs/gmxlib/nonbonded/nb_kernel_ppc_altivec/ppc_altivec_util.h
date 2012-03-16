/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _ALTIVEC_UTIL_H_
#define _ALTIVEC_UTIL_H_

/** @file ppc_altivec_util.h
 *
 *  @brief Altivec utility functions for optimized kernels.
 *
 * This file contains static inline utility functions that accomplish
 * tasks like loading/storing coordinates and forces, generating constaints,
 * and loading table data.
 *
 * Due to all the static functions it might take a while to compile files
 * that include this header, but from a performance point of view it makes
 * a tremendous difference.
 */
#include<stdio.h>

/* altivec.h must be included on vanilla gcc-4.0,
 * but not on Apple gcc or the IBM compilers.
 */
#ifdef HAVE_ALTIVEC_H
#include <altivec.h>
#endif

/** Write contents of a SIMD FP variable on standard out.
 *
 * @internal
 * 
 * @param v    SIMD floating-point variable to print.
 */
static void 
printvec(vector float v)
{
	int i;
	printf(" ");

	for(i=0;i<4;i++)
		printf("%8.5f ",*(((float *)&v)+i));
	printf("\n");
}


/** Set SIMD unit to use non-java rounding mode.
 *
 * @internal
 * 
 * On most PowerPC processors, FP operations take an additional clock 
 * cycle when the default java rounding mode is used. We couldn't care less,
 * so we can save a couple of percent of runtime by using classical IEEE mode.
 */
static void 
set_non_java_mode(void)
{
	vector unsigned short vsr1,vsr2;
	vector unsigned int tmp;

	vsr1=vec_mfvscr();
	tmp=vec_sl(vec_splat_u32(1),vec_splat_u32(8));
	vsr2=(vector unsigned short)vec_sl(tmp,vec_splat_u32(8));
	vsr1=vec_or(vsr1,vsr2);
	vec_mtvscr(vsr1);
}  

/** Create the SIMD FP constant 0.0
 *
 *  This routine returns a SIMD variable filled with 0.0 in all four elements,
 *  without loading any data from memory.
 *
 *  @return SIMD FP 0.0
 */
static inline vector float 
vec_zero(void)
{
	return vec_ctf(vec_splat_u32(0),0);
}

/** Create the SIMD FP constant 0.5
*
*  This routine returns a SIMD variable filled with 0.5 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 0.5
*/
static inline vector float 
vec_half(void)
{
	return vec_ctf(vec_splat_u32(1),1);
}


/** Create the SIMD FP constant 1.0
*
*  This routine returns a SIMD variable filled with 1.0 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 1.0
*/
static inline vector float 
vec_one(void)
{
	return vec_ctf(vec_splat_u32(1),0);
}


/** Create the SIMD FP constant 2.0
*
*  This routine returns a SIMD variable filled with 2.0 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 2.0
*/
static inline vector float 
vec_two(void)
{
	return vec_ctf(vec_splat_u32(2),0);
}


/** Create the SIMD FP constant 3.0
*
*  This routine returns a SIMD variable filled with 3.0 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 3.0
*/
static inline vector float 
vec_three(void)
{
	return vec_ctf(vec_splat_u32(3),0);
}


/** Create the SIMD FP constant 6.0
*
*  This routine returns a SIMD variable filled with 6.0 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 6.0
*/
static inline vector float 
vec_six(void)
{
	return vec_ctf(vec_splat_u32(6),0);
}

/** Create the SIMD FP constant 12.0
*
*  This routine returns a SIMD variable filled with 12.0 in all four elements,
*  without loading any data from memory.
*
*  @return SIMD FP 12.0
*/
static inline vector float 
vec_twelve(void)
{
	return vec_ctf(vec_splat_u32(12),0);
}


/** Load 3 floats from memory into elements 0-2 of a SIMD variable.
 *
 *  This routine loads 3 floating-point values from memory, which does not
 *  have to be aligned, and returns a SIMD variable with the values in the
 *  lower three elements.
 *
 *  @param address Pointer to values in memory.
 *
 *  @return SIMD FP variable with values in lower three elements.
 */
static inline vector float 
load_xyz(float *address)
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


/** Load four floats from unaligned memory 
*
*  This routine loads 4 floating-point values from memory, which does not
*  have to be aligned, and returns a SIMD variable with the values.
*
*  @param address Pointer to values in memory.
*
*  @return SIMD FP variable with values.
*/
static inline vector float 
load_vector_unaligned(float *address)
{
	vector unsigned char perm;
	vector float low,high;
  
	perm              = vec_lvsl( 0, (int *) address ); 
	low               = vec_ld(  0, address ); 
	high              = vec_ld( 16, address ); 
  
	return vec_perm(low,high,perm);
}

/** Load FP variable from memory, spread into all elements of SIMD variable.
*
*  This routine loads a single floating-point values from memory, which does not
*  have to be aligned, and returns a SIMD variable with this values in all
*  four elements.
*
*  @param address Pointer to value in memory.
*
*  @return SIMD FP variable with values in all elements.
*/
static inline vector float 
load_float_and_splat(float *address)
{
	vector unsigned char perm;
	vector float tmp;

	tmp               = vec_lde(0,address);
	perm              = vec_lvsl(0,address);
	tmp               = vec_perm(tmp,tmp,perm);

	return vec_splat(tmp,0);
}


/** Load 4 non-consecutive floats into a single SIMD variable.
*
*  This routine loads four floating-point values from different memory 
*  locations, which do not have to be aligned, and returns a SIMD variable 
*  with the four elements.
*
*  @param float1  Pointer to first value.
*  @param float2  Pointer to first value.
*  @param float3  Pointer to first value.
*  @param float4  Pointer to first value.
*
*  @return SIMD FP variable with the four values.
*/
static inline vector float 
load_4_float(float *float1,
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

/** Load 3 non-consecutive floats into a single SIMD variable.
*
*  This routine loads three floating-point values from different memory 
*  locations, which do not have to be aligned, and returns a SIMD variable 
*  with the values.
*
*  @param float1  Pointer to first value.
*  @param float2  Pointer to first value.
*  @param float3  Pointer to first value.
*
*  @return SIMD FP variable with values in lower three elements. The fourth
*          element is undefined.
*/
static inline vector float 
load_3_float(float *float1,
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


/** Load 2 non-consecutive floats into a single SIMD variable.
*
*  This routine loads two floating-point values from different memory 
*  locations, which do not have to be aligned, and returns a SIMD variable 
*  with the values.
*
*  @param float1  Pointer to first value.
*  @param float2  Pointer to first value.
*
*  @return SIMD FP variable with values in lower two elements. 
*          Elements 3 and 4 are undefined.
*/
static inline vector float 
load_2_float(float *float1,
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

/** Load a single float into a single SIMD variable.
*
*  This routine loads a floating-point value from memory,
*  which do not have to be aligned, and returns a SIMD variable 
*  with the value.
*
*  @param float1  Pointer to first value.
*
*  @return SIMD FP variable with value in lowest element. 
*          Elements 2, 3 and 4 are zeroed.
*/
static inline vector float 
load_1_float(float *float1)
{
	vector unsigned char xshift = vec_lvsl( 4, float1 ); 
  
	vector float X = vec_lde( 0, float1 ); 
	X = vec_perm( X, X, xshift); 
    
	return vec_sld(X,vec_zero(),12);
}



/** Store 4 floats from a SIMD variable to 4 different memory locations.
*
*  This routine stores four floating-point values to different memory 
*  locations, which do not have to be aligned.
*
*  @param v       SIMD variable with the values.
*  @param float1  Pointer to memory where first element will be stored.
*  @param float2  Pointer to memory where second element will be stored.
*  @param float3  Pointer to memory where third element will be stored.
*  @param float4  Pointer to memory where fourth element will be stored.
*
*/
static inline void 
store_4_float(vector float v,
              float *float1,
              float *float2,
              float *float3,
              float *float4)
{
	vector float e0 = vec_splat(v,0); 
	vector float e1 = vec_splat(v,1);
	vector float e2 = vec_splat(v,2);
  
	vector float f0 = vec_lde( 0, float1 ); 
	vector float f1 = vec_lde( 0, float2 ); 
	vector float f2 = vec_lde( 0, float3 ); 
	vector float f3 = vec_lde( 0, float4 ); 
  
	v               = vec_splat(v,3);
  
	e0              = vec_add(e0,f0);
	e1              = vec_add(e1,f1);
	e2              = vec_add(e2,f2);
	v               = vec_add(v,f3);
  
	vec_ste(e0, 0, float1);
	vec_ste(e1, 0, float2);
	vec_ste(e2, 0, float3);
	vec_ste(v, 0, float4);
}


/** Store 3 floats from a SIMD variable to 3 different memory locations.
*
*  This routine stores three floating-point values to different memory 
*  locations, which do not have to be aligned.
*
*  @param v       SIMD variable with the values.
*  @param float1  Pointer to memory where first element will be stored.
*  @param float2  Pointer to memory where second element will be stored.
*  @param float3  Pointer to memory where third element will be stored.
*
*/
static inline void 
store_3_float(vector float v,
              float *float1,
              float *float2,
              float *float3)
{
	vector float e0 = vec_splat(v,0); 
	vector float e1 = vec_splat(v,1);
  
	vector float f0 = vec_lde( 0, float1 ); 
	vector float f1 = vec_lde( 0, float2 ); 
	vector float f2 = vec_lde( 0, float3 ); 
  
	v               = vec_splat(v,2);
  
	e0              = vec_add(e0,f0);
	e1              = vec_add(e1,f1);
	v               = vec_add(v,f2);
  
	vec_ste(e0, 0, float1);
	vec_ste(e1, 0, float2);
	vec_ste(v, 0, float3);
}

/** Store 2 floats from a SIMD variable to 2 different memory locations.
*
*  This routine stores two floating-point values to different memory 
*  locations, which do not have to be aligned.
*
*  @param v       SIMD variable with the values.
*  @param float1  Pointer to memory where first element will be stored.
*  @param float2  Pointer to memory where second element will be stored.
*/
static inline void 
store_2_float(vector float    v,
              float *         float1,
              float *         float2)
{
	vector float e0 = vec_splat(v,0); 
  
	vector float f0 = vec_lde( 0, float1 ); 
	vector float f1 = vec_lde( 0, float2 ); 

	v               = vec_splat(v,1);

	e0              = vec_add(e0,f0);
	v               = vec_add(v,f1);
  
	vec_ste(e0, 0, float1);
	vec_ste(v, 0, float2);
}


/** Store first element of SIMD variable to memory.
*
*  This routine stores the first floating-point entry of a SIMD variable
*  to memory, which do not have to be aligned.
*
*  @param v       SIMD variable with the value.
*  @param float1  Pointer to memory where first element will be stored.
*/
static inline void 
store_1_float(vector float     v,
              float *          float1)
{
	vector float f0 = vec_lde( 0, float1 ); 
  
	v               = vec_splat(v,0);
	v               = vec_add(v,f0);
  
	vec_ste(v, 0, float1);
}


/** Load four pairs of floats and store into two SIMD variables.
 *
 * This routine is designed to load LJ parameters for four atoms and store in 
 * two vectors, one for c6 and one for c12.
 *
 * The vdwparam memory is aligned on an 8-byte boundary and consists
 * of pairs, so the two parameters are either in the upper or 
 * lower half of a 16-byte structure. The first value in the pair is the c6
 * parameter, and the second c12.
 *
 * @param  pair1   Pointer to LJ parameters for first atom.
 * @param  pair2   Pointer to LJ parameters for second atom.
 * @param  pair3   Pointer to LJ parameters for third atom.
 * @param  pair4   Pointer to LJ parameters for fourth atom.
 * @param  c6      Output: SIMD value with c6 for the four atoms.
 * @param  c12     Output: SIMD value with c12 for the four atoms.
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


/** Load three pairs of floats and store into two SIMD variables.
*
* This routine is designed to load LJ parameters for three atoms and store in 
* two vectors, one for c6 and one for c12.
*
* The vdwparam memory is aligned on an 8-byte boundary and consists
* of pairs, so the two parameters are either in the upper or 
* lower half of a 16-byte structure. The first value in the pair is the c6
* parameter, and the second c12.
*
* @param  pair1   Pointer to LJ parameters for first atom.
* @param  pair2   Pointer to LJ parameters for second atom.
* @param  pair3   Pointer to LJ parameters for third atom.
* @param  c6      Output: SIMD value with c6 for the atoms in elements 0-2.
* @param  c12     Output: SIMD value with c12 for the atoms in elements 0-2.
*/
static inline void 
load_3_pair(float *pair1,
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


/** Load two pairs of floats and store into two SIMD variables.
*
* This routine is designed to load LJ parameters for two atoms and store in 
* two vectors, one for c6 and one for c12.
*
* The vdwparam memory is aligned on an 8-byte boundary and consists
* of pairs, so the two parameters are either in the upper or 
* lower half of a 16-byte structure. The first value in the pair is the c6
* parameter, and the second c12.
*
* @param  pair1   Pointer to LJ parameters for first atom.
* @param  pair2   Pointer to LJ parameters for second atom.
* @param  c6      Output: SIMD value with c6 for the atoms in elements 0,1.
* @param  c12     Output: SIMD value with c12 for the atoms in elements 0,1.
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

/** Load a pair of floats and store into two SIMD variables.
*
* This routine is designed to load LJ parameters for a single atom and store in 
* two vectors, one for c6 and one for c12.
*
* The vdwparam memory is aligned on an 8-byte boundary and consists
* of pairs, so the two parameters are either in the upper or 
* lower half of a 16-byte structure. The first value in the pair is the c6
* parameter, and the second c12.
*
* @param  pair1   Pointer to LJ parameters for atom.
* @param  c6      Output: SIMD value with c6 for the atom in lowest element.
* @param  c12     Output: SIMD value with c12 for the atom in lowest element.
*/
static inline void 
load_1_pair(float *pair1,
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


/** Spread a coordinate triplet to three separate vectors.
 *
 * This routine permutes a SIMD variable [x y z ?] into three separate 
 * variables: [x x x x] [y y y y] [z z z z] .
 * *
 * @param xyzreg  SIMD variable containing [x y z ?]
 * @param xreg    Output: SIMD variable [x x x x]
 * @param yreg    Output: SIMD variable [y y y y]
 * @param zreg    Output: SIMD variable [z z z z]
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


/** Transpose four coordinate triplets into common X/Y/Z variables.
*
* This routine permutes four SIMD variables [x y z ?] into 
* [x1 x2 x3 x4] [y1 y2 y3 y4] [z1 z2 z3 z4]  
*
* @param xyz1     SIMD variable containing first coordinate triplet.
* @param xyz2     SIMD variable containing second coordinate triplet.
* @param xyz3     SIMD variable containing third coordinate triplet.
* @param xyz4     SIMD variable containing fourth coordinate triplet.
* @param xvector  Output: SIMD variable [x1 x2 x3 x4]
* @param yvector  Output: SIMD variable [y1 y2 y3 y4]
* @param zvector  Output: SIMD variable [z1 z2 z3 z4]
*/
static inline void 
transpose_4_to_3(vector float xyz1,
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


/** Transpose two coordinate triplets into common X/Y/Z variables.
*
* This routine permutes two SIMD variables [x y z ?] into 
* [x1 x2 ? ?] [y1 y2 ? ?] [z1 z2 ? ?]  
*
*
* @param xyz1     SIMD variable containing first coordinate triplet.
* @param xyz2     SIMD variable containing second coordinate triplet.
* @param xvector  Output: SIMD variable [x1 x2 ? ?]
* @param yvector  Output: SIMD variable [y1 y2 ? ?]
* @param zvector  Output: SIMD variable [z1 z2 ? ?]
*/
static inline void 
transpose_2_to_3(vector float xyz1,
                 vector float xyz2,
                 vector float *xvector,
                 vector float *yvector,
                 vector float *zvector)
{
	*xvector = vec_mergeh(xyz1,xyz2);            /* x1 x2 y1 y2 */
	*zvector = vec_mergel(xyz1,xyz2);            /* z1 z2  ?  ? */
	*yvector = vec_sld(*xvector,*xvector,8); /* y1 y2 0 0 */
}


/** Transpose a single coordinate triplets into separate X/Y/Z variables.
*
* This routine permutes a SIMD variablea [x y z ?] into 
* [x1 ? ? ?] [y1 ? ? ?] [z1 ? ? ?]  
*
*
* @param xyz1     SIMD variable containing coordinate triplet.
* @param xvector  Output: SIMD variable [x1 ? ? ?]
* @param yvector  Output: SIMD variable [y1 ? ? ?]
* @param zvector  Output: SIMD variable [z1 ? ? ?]
*/
static inline void 
transpose_1_to_3(vector float xyz1,
                 vector float *xvector,
                 vector float *yvector,
                 vector float *zvector)
{
	/* simply use splat, since elem 2,3,4 dont matter. */
	*xvector           = vec_splat(xyz1,0);   /* x1 x1 x1 x1 */
	*yvector           = vec_splat(xyz1,1);   /* y1 y1 y1 y1 */
	*zvector           = vec_splat(xyz1,2);   /* z1 z1 z1 z1 */
}




/** Transpose three separate X/Y/Z variables into 4 coordinate triplets.
*
* This routine permutes three separate SIMD coordinate vectors 
* [x1 x2 x3 x4] [y1 y2 y3 y4] [z1 z2 z3 z4] into 4 triplets [x y z 0].
*
* The fourth position of the output triplets will be set to zero.
*
* @param xvector  SIMD variable [x1 x2 x3 x4]
* @param yvector  SIMD variable [y1 y2 y3 y4]
* @param zvector  SIMD variable [z1 z2 z3 z4]
* @param xyz1     Output: SIMD variable containing first coordinate triplet.
* @param xyz2     Output: SIMD variable containing second coordinate triplet.
* @param xyz3     Output: SIMD variable containing third coordinate triplet.
* @param xyz4     Output: SIMD variable containing fourth coordinate triplet.
*/
static inline void 
transpose_3_to_4(vector float xvector,
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




/** Transpose three separate X/Y/Z variables into 2 coordinate triplets.
*
* This routine permutes three separate SIMD coordinate vectors 
* [x1 x2 ? ?] [y1 y2 ? ?] [z1 z2 ? ?] into 2 triplets [x y z 0].
*
* The fourth position of the output triplets will be set to zero.
*
* @param xvector  SIMD variable [x1 x2 ? ?]
* @param yvector  SIMD variable [y1 y2 ? ?]
* @param zvector  SIMD variable [z1 z2 ? ?]
* @param xyz1     Output: SIMD variable containing first coordinate triplet.
* @param xyz2     Output: SIMD variable containing second coordinate triplet.
*/
static inline void 
transpose_3_to_2(vector float xvector,
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


/** Transpose three separate X/Y/Z variables into a single coordinate triplet.
*
* This routine permutes three separate SIMD coordinate vectors 
* [x1 ? ? ?] [y1 ? ? ?] [z1 ? ? ?] into a single triplet [x y z 0].
*
* The fourth position of the output triplet will be set to zero.
*
* @param xvector  SIMD variable [x1 ? ? ?]
* @param yvector  SIMD variable [y1 ? ? ?]
* @param zvector  SIMD variable [z1 ? ? ?]
* @param xyz1     Output: SIMD variable containing coordinate triplet.
*/
static inline void 
transpose_3_to_1(vector float xvector,
                 vector float yvector,
                 vector float zvector,
                 vector float *xyz1)
{
	*xyz1       = vec_mergeh(xvector,zvector); /* [x1 z1 ? ? ] */
	yvector     = vec_mergeh(yvector,vec_zero()); /* [y1 0 ? 0] */
	*xyz1       = vec_mergeh(*xyz1,yvector);   /* [x1 y1 z1 0] */
}



/** Full 4-to-4 transpose of SIMD variables.
*
* If the input SIMD vectors are pictured as rows in a 4x4 matrix, the output
* result will be the column vectors of the same matrix.
*
* @param in1      First row vector.
* @param in2      Second row vector.
* @param in3      Third row vector.
* @param in4      Fourth row vector.
* @param out1     Output: First column vector.
* @param out2     Output: Second column vector.
* @param out3     Output: Third column vector.
* @param out4     Output: Fourth column vector.
*/
static inline void 
transpose_4_to_4(vector float in1,
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


/** Transpose the first two elements of 4 SIMD variables to new variables.
*
* If the input SIMD vectors are pictured as rows in a 4x4 matrix, the output
* result will be the first two column vectors of the same matrix.
*
* @param in1      First row vector.
* @param in2      Second row vector.
* @param in3      Third row vector.
* @param in4      Fourth row vector.
* @param out1     Output: First column vector.
* @param out2     Output: Second column vector.
*/
static inline void 
transpose_4_to_2(vector float in1,
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



/** Transpose the first element of 4 SIMD variables to new variables.
*
* If the input SIMD vectors are pictured as rows in a 4x4 matrix, the output
* result will be the first column vector of the same matrix.
*
* @param in1      First row vector.
* @param in2      Second row vector.
* @param in3      Third row vector.
* @param in4      Fourth row vector.
* @param out1     Output: First column vector.
*/
static inline void 
transpose_4_to_1(vector float in1,
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



/** Transpose two SIMD variables into the two first elements of 4 new ones.
*
* If the input SIMD vectors are pictured as the two first rows in a 4x4 matrix, 
* the output result will be the column column vectors of the same matrix.
*
* @param in1      First row vector.
* @param in2      Second row vector.
* @param out1     Output: First column vector.
* @param out2     Output: Second column vector.
* @param out3     Output: Third column vector.
* @param out4     Output: Fourth column vector.
*/
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


/** Transpose a SIMD variable into the first element of 4 new ones.
*
* If the input SIMD vectors are pictured as the first row in a 4x4 matrix, 
* the output result will be the column column vectors of the same matrix.
*
* @param in1      First row vector.
* @param out1     Output: First column vector.
* @param out2     Output: Second column vector.
* @param out3     Output: Third column vector.
* @param out4     Output: Fourth column vector.
*/
static inline void
transpose_1_to_4(vector float in1,
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


/** Add first three elements of a SIMD variable to three memory values.
 *
 * This routine adds each of the three lowest elements in a FP SIMD variable
 * to a corresponding location in (unaligned). Note that we do not accumulate
 * the SIMD entries - three separate values in the SIMD variable are added to 
 * three separate (but adjacent) memory locations.
 *
 * @param address Pointer to the float triplet in memory to add to.
 * @param vdata   SIMD variable whose three lowest elements should be added
 *                to corresponding positions in memory.
 */
static inline void 
add_xyz_to_mem(float *         address,
               vector float    vdata)
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


/** Add the elements of a SIMD variable to four memory values.
*
* This routine adds each of the elements in a FP SIMD variable
* to a corresponding location in (unaligned). Note that we do not accumulate
* the SIMD entries - the separate values in the SIMD variable are added to 
* four separate (but adjacent) memory locations.
*
* @param address Pointer to the 4 float values in memory to add to.
* @param vdata   SIMD variable whose elements should be added
*                to corresponding positions in memory.
*/
static inline void 
add_vector_to_float(float *         address,
                    vector float    vdata)
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


/** Load coords of a 3-atom water (9 floats) to all elements in SIMD vectors.
 *
 * We use our knowledge about water to improve the loading performance.
 * The allocated memory is always 8-byte aligned, so with nine atoms
 * we can always read 3 16-byte chunks (rounded down by vec_ld)
 * without going outside the current page (which might incur a page fault). 
 * The extra elements before and after the water dont matter as long as we dont
 * try to write to them.
 *
 * @param address    Address of first coord in water (oxygen x).
 * @param shiftvec   Address of shift vector x/y/z triplet to add to
 *                   all atom coordinates.
 * @param Ox         Output: Oxygen x coordinate in all elements.
 * @param Oy         Output: Oxygen y coordinate in all elements.
 * @param Oz         Output: Oxygen z coordinate in all elements.
 * @param H1x        Output: First hydrogen x coordinate in all elements.
 * @param H1y        Output: First hydrogen y coordinate in all elements.
 * @param H1z        Output: First hydrogen z coordinate in all elements.
 * @param H2x        Output: Second hydrogen x coordinate in all elements.
 * @param H2y        Output: Second hydrogen y coordinate in all elements.
 * @param H2z        Output: Second hydrogen z coordinate in all elements.
 */
static inline void 
load_1_3atoms_shift_and_splat(float *address,
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


/** Load coords of a 4-atom water (12 floats) to all elements in SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with twelve atoms
* we can always read 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address    Address of first coord in 4-atom water (oxygen x).
* @param shiftvec   Address of shift vector x/y/z triplet to add to
*                   all atom coordinates.
* @param Ox         Output: Oxygen x coordinate in all elements.
* @param Oy         Output: Oxygen y coordinate in all elements.
* @param Oz         Output: Oxygen z coordinate in all elements.
* @param H1x        Output: First hydrogen x coordinate in all elements.
* @param H1y        Output: First hydrogen y coordinate in all elements.
* @param H1z        Output: First hydrogen z coordinate in all elements.
* @param H2x        Output: Second hydrogen x coordinate in all elements.
* @param H2y        Output: Second hydrogen y coordinate in all elements.
* @param H2z        Output: Second hydrogen z coordinate in all elements.
* @param Mx         Output: Virtual site x coordinate in all elements.
* @param My         Output: Virtual site y coordinate in all elements.
* @param Mz         Output: Virtual site z coordinate in all elements.
*/
static inline void load_1_4atoms_shift_and_splat(float *address,
                                                 float *shiftvec,
                                                 vector float *Ox,
                                                 vector float *Oy,
                                                 vector float *Oz,
                                                 vector float *H1x,
                                                 vector float *H1y,
                                                 vector float *H1z,
                                                 vector float *H2x,
                                                 vector float *H2y,
                                                 vector float *H2z,
                                                 vector float *Mx,
                                                 vector float *My,
                                                 vector float *Mz)
{
	vector unsigned char perm;
	vector float c1,c2,c3,c4,sha,shb,sh1,sh2,sh3;

	/* load shift data */
	shb               = load_xyz(shiftvec);  /* [ shX shY shZ -] */
	/* load the coordinates. */
	c1                = vec_ld(  0, address ); 
	c2                = vec_ld( 16, address ); 
	c3                = vec_ld( 32, address ); 
	if(((unsigned int)address) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		c4                = vec_ld( 48, address ); 
		perm              = vec_lvsl( 0, (int *) address ); 
		c1                = vec_perm(c1,c2,perm);      /*  Ox  Oy  Oz H1x */
		c2                = vec_perm(c2,c3,perm);      /* H1y H1z H2x H2y */
		c3                = vec_perm(c3,c4,perm);      /* H2z  Mx  My  Mz */
	}

	sha               = vec_sld(shb,shb,12);       /* [  -  shX shY shZ ] */
	sh1               = vec_sld(sha,shb,4);        /* [ shX shY shZ shX ] */
	sh2               = vec_sld(sha,shb,8);        /* [ shY shZ shX shY ] */
	sh3               = vec_sld(sha,shb,12);       /* [ shZ shX shY shZ ] */

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
	*Mx              = vec_splat(c3,1);
	*My              = vec_splat(c3,2);
	*Mz              = vec_splat(c3,3);
}


/** Load coords of 4 3-atom waters (4*9 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with nine atoms
* we can always read 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param address3   Address of first coord in third water (oxygen x).
* @param address4   Address of first coord in fourth water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
*/
static inline void
load_4_3atoms(float *address1,
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

	perm              = vec_lvsl( 0, address2 );
	tmp1              = vec_ld(  0, address2 );
	tmp2              = vec_ld( 16, address2 );
	tmp3              = vec_ld( 32, address2 );
	c1b               = vec_perm(tmp1,tmp2,perm);     /*  Oxb  Oyb  Ozb H1xb */
	c2b               = vec_perm(tmp2,tmp3,perm);     /* H1yb H1zb H2xb H2yb */
	c3b               = vec_perm(tmp3,tmp3,perm);     /* H2zb   -   -   - */

	perm              = vec_lvsl( 0, address3 );
	tmp1              = vec_ld(  0, address3 );
	tmp2              = vec_ld( 16, address3 );
	tmp3              = vec_ld( 32, address3 );
	c1c               = vec_perm(tmp1,tmp2,perm);     /*  Oxc  Oyc  Ozc H1xc */
	c2c               = vec_perm(tmp2,tmp3,perm);     /* H1yc H1zc H2xc H2yc */
	c3c               = vec_perm(tmp3,tmp3,perm);     /* H2zc   -   -   - */

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



/** Load coords of 4 4-atom waters (4*12 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param address3   Address of first coord in third water (oxygen x).
* @param address4   Address of first coord in fourth water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
* @param Mx         Output: Virtual site x coordinates.
* @param My         Output: Virtual site y coordinates.
* @param Mz         Output: Virtual site z coordinates.
*/
static inline void 
load_4_4atoms(float *address1,
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
              vector float *H2z,
              vector float *Mx,
              vector float *My,
              vector float *Mz)
{
	vector unsigned char perm;
	vector float tmp1,tmp2,tmp4;
  
	vector float c1a,c2a,c3a;
	vector float c1b,c2b,c3b;
	vector float c1c,c2c,c3c;
	vector float c1d,c2d,c3d;
 
	/* load the coordinates */
	c1a                 = vec_ld(  0, address1 ); 
	c2a                 = vec_ld( 16, address1 ); 
	c3a                 = vec_ld( 32, address1 ); 
	if(((unsigned int)address1) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address1 );
		perm              = vec_lvsl( 0, address1 ); 
		c1a               = vec_perm(c1a,c2a,perm);      /*  Oxa  Oya  Oza H1xa */
		c2a               = vec_perm(c2a,c3a,perm);      /* H1ya H1za H2xa H2ya */
		c3a               = vec_perm(c3a,tmp4,perm);     /* H2za  Mxa  Mya  Mza */
	}

	c1b                 = vec_ld(  0, address2 ); 
	c2b                 = vec_ld( 16, address2 ); 
	c3b                 = vec_ld( 32, address2 ); 
	if(((unsigned int)address2) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address2 );
		perm              = vec_lvsl( 0, address2 ); 
		c1b               = vec_perm(c1b,c2b,perm);      /*  Oxb  Oyb  Ozb H1xb */
		c2b               = vec_perm(c2b,c3b,perm);      /* H1yb H1zb H2xb H2yb */
		c3b               = vec_perm(c3b,tmp4,perm);     /* H2zb  Mxb  Myb  Mzb */
	}

	c1c                 = vec_ld(  0, address3 ); 
	c2c                 = vec_ld( 16, address3 ); 
	c3c                 = vec_ld( 32, address3 ); 
	if(((unsigned int)address3) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address3 );
		perm              = vec_lvsl( 0, address3 ); 
		c1c               = vec_perm(c1c,c2c,perm);     /*  Oxc  Oyc  Ozc H1xc */
		c2c               = vec_perm(c2c,c3c,perm);     /* H1yc H1zc H2xc H2yc */
		c3c               = vec_perm(c3c,tmp4,perm);    /* H2zc  Mxc  Myc  Mzc */
	}

	c1d                 = vec_ld(  0, address4 ); 
	c2d                 = vec_ld( 16, address4 ); 
	c3d                 = vec_ld( 32, address4 ); 
	if(((unsigned int)address4) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address4 );
		perm              = vec_lvsl( 0, address4 ); 
		c1d               = vec_perm(c1d,c2d,perm);      /*  Oxd  Oyd  Ozd H1xd */
		c2d               = vec_perm(c2d,c3d,perm);      /* H1yd H1zd H2xd H2yd */
		c3d               = vec_perm(c3d,tmp4,perm);     /* H2zd  Mxd  Myd  Mzd */
	}

	/* permute things */
	tmp1              = vec_mergeh(c1a,c1c);          /*  Oxa  Oxc  Oya  Oyc */
	c1a               = vec_mergel(c1a,c1c);          /*  Oza  Ozc H1xa H1xc */
	c1c               = vec_mergeh(c1b,c1d);          /*  Oxb  Oxd  Oyb  Oyd */
	c1b               = vec_mergel(c1b,c1d);          /*  Ozb  Ozd H1xb H1xd */

	c1d               = vec_mergeh(c2a,c2c);          /* H1ya H1yc H1za H1zc */
	c2a               = vec_mergel(c2a,c2c);          /* H2xa H2xc H2ya H2yc */
	c2c               = vec_mergeh(c2b,c2d);          /* H1yb H1yd H1zb H1zd */
	c2b               = vec_mergel(c2b,c2d);          /* H2xb H2xd H2yb H2yd */

	tmp2              = vec_mergeh(c3a,c3c);          /* H2za H2zc  Mxa  Mxc */
	c3a               = vec_mergel(c3a,c3c);          /*  Mya  Myc  Mza  Mzc */
	c3c               = vec_mergeh(c3b,c3d);          /* H2zb H2zd  Mxb  Mxd */
	c3b               = vec_mergel(c3b,c3d);          /*  Myb  Myd  Mzb  Mzd */

  
	*Ox               = vec_mergeh(tmp1,c1c);         /*  Oxa  Oxb  Oxc  Oxd */
	*Oy               = vec_mergel(tmp1,c1c);         /*  Oya  Oyb  Oyc  Oyd */
	*Oz               = vec_mergeh(c1a,c1b);          /*  Oza  Ozb  Ozc  Ozd */
	*H1x              = vec_mergel(c1a,c1b);          /* H1xa H1xb H1xc H1xd */
	*H1y              = vec_mergeh(c1d,c2c);          /* H1ya H1yb H1yc H1yd */
	*H1z              = vec_mergel(c1d,c2c);          /* H1za H1zb H1zc H1zd */
	*H2x              = vec_mergeh(c2a,c2b);          /* H2xa H2xb H2xc H2xd */
	*H2y              = vec_mergel(c2a,c2b);          /* H2ya H2yb H2yc H2yd */

	*H2z              = vec_mergeh(tmp2,c3c);          /* H2za H2zb H2zc H2zd */
	*Mx               = vec_mergel(tmp2,c3c);          /*  Mxa  Mxb  Mxc  Mxd */
	*My               = vec_mergeh(c3a,c3b);           /*  Mya  Myb  Myc  Myd */
	*Mz               = vec_mergel(c3a,c3b);           /*  Mza  Mzb  Mzc  Mzd */
}



/** Load coords of 3 3-atom waters (3*9 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with nine atoms
* we can always read 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param address3   Address of first coord in third water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
*/
static inline void 
load_3_3atoms(float *address1,
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


/** Load coords of 3 4-atom waters (3*12 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param address3   Address of first coord in third water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
* @param Mx         Output: Virtual site x coordinates.
* @param My         Output: Virtual site y coordinates.
* @param Mz         Output: Virtual site z coordinates.
*/
static inline void 
load_3_4atoms(float *address1,
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
              vector float *H2z,
              vector float *Mx,
              vector float *My,
              vector float *Mz)
{
	vector unsigned char perm;
	vector float tmp1,tmp2,tmp3,tmp4;
	vector float c1a,c2a,c3a;
	vector float c1b,c2b,c3b;
	vector float c1c,c2c,c3c;
  
	/* load the coordinates */
	c1a                 = vec_ld(  0, address1 ); 
	c2a                 = vec_ld( 16, address1 ); 
	c3a                 = vec_ld( 32, address1 ); 
	if(((unsigned int)address1) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address1 );
		perm              = vec_lvsl( 0, address1 ); 
		c1a               = vec_perm(c1a,c2a,perm);      /*  Oxa  Oya  Oza H1xa */
		c2a               = vec_perm(c2a,c3a,perm);      /* H1ya H1za H2xa H2ya */
		c3a               = vec_perm(c3a,tmp4,perm);     /* H2za  Mxa  Mya  Mza */
	}

	c1b                 = vec_ld(  0, address2 ); 
	c2b                 = vec_ld( 16, address2 ); 
	c3b                 = vec_ld( 32, address2 ); 
	if(((unsigned int)address2) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address2 );
		perm              = vec_lvsl( 0, address2 ); 
		c1b               = vec_perm(c1b,c2b,perm);      /*  Oxb  Oyb  Ozb H1xb */
		c2b               = vec_perm(c2b,c3b,perm);      /* H1yb H1zb H2xb H2yb */
		c3b               = vec_perm(c3b,tmp4,perm);     /* H2zb  Mxb  Myb  Mzb */
	}

	c1c                 = vec_ld(  0, address3 ); 
	c2c                 = vec_ld( 16, address3 ); 
	c3c                 = vec_ld( 32, address3 ); 
	if(((unsigned int)address3) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address3 );
		perm              = vec_lvsl( 0, address3 ); 
		c1c               = vec_perm(c1c,c2c,perm);     /*  Oxc  Oyc  Ozc H1xc */
		c2c               = vec_perm(c2c,c3c,perm);     /* H1yc H1zc H2xc H2yc */
		c3c               = vec_perm(c3c,tmp4,perm);    /* H2zc  Mxc  Myc  Mzc */
	}

	/* permute things */
	tmp1              = vec_mergeh(c1a,c1c);          /*  Oxa  Oxc  Oya  Oyc */
	c1a               = vec_mergel(c1a,c1c);          /*  Oza  Ozc H1xa H1xc */
	c1c               = vec_mergeh(c1b,c1b);          /*  Oxb   -   Oyb   -  */
	c1b               = vec_mergel(c1b,c1b);          /*  Ozb   -  H1xb   -  */

	tmp2              = vec_mergeh(c2a,c2c);          /* H1ya H1yc H1za H1zc */
	c2a               = vec_mergel(c2a,c2c);          /* H2xa H2xc H2ya H2yc */
	c2c               = vec_mergeh(c2b,c2b);          /* H1yb   -  H1zb   -  */
	c2b               = vec_mergel(c2b,c2b);          /* H2xb   -  H2yb   -  */

	tmp3              = vec_mergeh(c3a,c3c);          /* H2za H2zc  Mxa  Mxc */
	c3a               = vec_mergel(c3a,c3c);          /*  Mya  Myc  Mza  Mzc */
	c3c               = vec_mergeh(c3b,c3b);          /* H2zb   -   Mxb   -  */
	c3b               = vec_mergel(c3b,c3b);          /*  Myb   -   Mzb   -  */

	*Ox               = vec_mergeh(tmp1,c1c);         /*  Oxa  Oxb  Oxc   -  */
	*Oy               = vec_mergel(tmp1,c1c);         /*  Oya  Oyb  Oyc   -  */
	*Oz               = vec_mergeh(c1a,c1b);          /*  Oza  Ozb  Ozc   -  */
	*H1x              = vec_mergel(c1a,c1b);          /* H1xa H1xb H1xc   -  */
	*H1y              = vec_mergeh(tmp2,c2c);         /* H1ya H1yb H1yc   -  */
	*H1z              = vec_mergel(tmp2,c2c);         /* H1za H1zb H1zc   -  */
	*H2x              = vec_mergeh(c2a,c2b);          /* H2xa H2xb H2xc   -  */
	*H2y              = vec_mergel(c2a,c2b);          /* H2ya H2yb H2yc   -  */
  
	*H2z              = vec_mergeh(tmp3,c3c);         /* H2za H2zb H2zc   -  */
	*Mx               = vec_mergel(tmp3,c3c);         /*  Mxa  Mxb  Mxc   -  */
	*My               = vec_mergeh(c3a,c3b);          /*  Mya  Myb  Myc   -  */
	*Mz               = vec_mergel(c3a,c3b);          /*  Mza  Mzb  Mzc   -  */
}



/** Load coords of 2 3-atom waters (2*9 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with nine atoms
* we can always read 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
*/
static inline void 
load_2_3atoms(float *address1,
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


/** Load coords of 2 4-atom waters (2*12 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param address2   Address of first coord in second water (oxygen x).
* @param Ox         Output: Oxygen x coordinates.
* @param Oy         Output: Oxygen y coordinates.
* @param Oz         Output: Oxygen z coordinates.
* @param H1x        Output: First hydrogen x coordinates.
* @param H1y        Output: First hydrogen y coordinates.
* @param H1z        Output: First hydrogen z coordinates.
* @param H2x        Output: Second hydrogen x coordinates.
* @param H2y        Output: Second hydrogen y coordinates.
* @param H2z        Output: Second hydrogen z coordinates.
* @param Mx         Output: Virtual site x coordinates.
* @param My         Output: Virtual site y coordinates.
* @param Mz         Output: Virtual site z coordinates.
*/
static inline void
load_2_4atoms(float *address1,
              float *address2,
              vector float *Ox,
              vector float *Oy,
              vector float *Oz,
              vector float *H1x,
              vector float *H1y,
              vector float *H1z,
              vector float *H2x,
              vector float *H2y,
              vector float *H2z,
              vector float *Mx,
              vector float *My,
              vector float *Mz)
{
	vector unsigned char perm;
	vector float tmp4;
	vector float c1a,c2a,c3a;
	vector float c1b,c2b,c3b;
  
	/* load the coordinates */
  
	/* load the coordinates */
	c1a                 = vec_ld(  0, address1 ); 
	c2a                 = vec_ld( 16, address1 ); 
	c3a                 = vec_ld( 32, address1 ); 
	if(((unsigned int)address1) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address1 );
		perm              = vec_lvsl( 0, address1 ); 
		c1a               = vec_perm(c1a,c2a,perm);      /*  Oxa  Oya  Oza H1xa */
		c2a               = vec_perm(c2a,c3a,perm);      /* H1ya H1za H2xa H2ya */
		c3a               = vec_perm(c3a,tmp4,perm);     /* H2za  Mxa  Mya  Mza */
	}

	c1b                 = vec_ld(  0, address2 ); 
	c2b                 = vec_ld( 16, address2 ); 
	c3b                 = vec_ld( 32, address2 ); 
	if(((unsigned int)address2) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address2 );
		perm              = vec_lvsl( 0, address2 ); 
		c1b               = vec_perm(c1b,c2b,perm);      /*  Oxb  Oyb  Ozb H1xb */
		c2b               = vec_perm(c2b,c3b,perm);      /* H1yb H1zb H2xb H2yb */
		c3b               = vec_perm(c3b,tmp4,perm);     /* H2zb  Mxb  Myb  Mzb */
	}

	/* never mind what we get in the two remaining elements */
	*Ox               = vec_mergeh(c1a,c1b);          /*  Oxa  Oxb  Oya  Oyb */ 
	*H1y              = vec_mergeh(c2a,c2b);          /* H1ya H1yb H1za H1zb */
	*H2z              = vec_mergeh(c3a,c3b);          /* H2za H2zb  Mxa  Mxb */  
	*My               = vec_mergel(c3a,c3b);          /*  Mya  Myb  Mza  Mzb */
	*Oz               = vec_mergel(c1a,c1b);          /*  Oza  Ozb H1xa H1xb */ 
	*H2x              = vec_mergel(c2a,c2b);          /* H2xa H2xb H2ya H2yb */
	*Oy               = vec_sld(*Ox,*Ox,8);
	*H1z              = vec_sld(*H1y,*H1y,8); 
	*H1x              = vec_sld(*Oz,*Oz,8);
	*H2y              = vec_sld(*H2x,*H2x,8);
	*Mx               = vec_sld(*H2z,*H2z,8);
	*Mz               = vec_sld(*My,*My,8);
}


/** Load coords of a 3-atom water (9 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with nine atoms
* we can always read 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param Ox         Output: Oxygen x coordinate in lowest element.
* @param Oy         Output: Oxygen y coordinate in lowest element.
* @param Oz         Output: Oxygen z coordinate in lowest element.
* @param H1x        Output: First hydrogen x coordinate in lowest element.
* @param H1y        Output: First hydrogen y coordinate in lowest element.
* @param H1z        Output: First hydrogen z coordinate in lowest element.
* @param H2x        Output: Second hydrogen x coordinate in lowest element.
* @param H2y        Output: Second hydrogen y coordinate in lowest element.
* @param H2z        Output: Second hydrogen z coordinate in lowest element.
*/
static inline void 
load_1_3atoms(float *address1,
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


/** Load coords of a 4-atom water (12 floats) to SIMD vectors.
*
* We use our knowledge about water to improve the loading performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to write to them.
*
* @param address1   Address of first coord in first water (oxygen x).
* @param Ox         Output: Oxygen x coordinate in lowest element.
* @param Oy         Output: Oxygen y coordinate in lowest element.
* @param Oz         Output: Oxygen z coordinate in lowest element.
* @param H1x        Output: First hydrogen x coordinate in lowest element.
* @param H1y        Output: First hydrogen y coordinate in lowest element.
* @param H1z        Output: First hydrogen z coordinate in lowest element.
* @param H2x        Output: Second hydrogen x coordinate in lowest element.
* @param H2y        Output: Second hydrogen y coordinate in lowest element.
* @param H2z        Output: Second hydrogen z coordinate in lowest element.
* @param Mx         Output: Virtual site x coordinate in lowest element.
* @param My         Output: Virtual site y coordinate in lowest element.
* @param Mz         Output: Virtual site z coordinate in lowest element.
*/
static inline void 
load_1_4atoms(float *address1,
              vector float *Ox,
              vector float *Oy,
              vector float *Oz,
              vector float *H1x,
              vector float *H1y,
              vector float *H1z,
              vector float *H2x,
              vector float *H2y,
              vector float *H2z,
              vector float *Mx,
              vector float *My,
              vector float *Mz)
{
	vector unsigned char perm;
	vector float tmp4;

	/* load the coordinates */
	*Ox                 = vec_ld(  0, address1 ); 
	*H1y                = vec_ld( 16, address1 ); 
	*H2z                = vec_ld( 32, address1 ); 
	if(((unsigned int)address1) & 0xf) { 
		/*  Load extra data and rotate if not aligned */
		tmp4              = vec_ld( 48, address1 );
		perm              = vec_lvsl( 0, address1 ); 
		*Ox               = vec_perm(*Ox,*H1y,perm);     /*  Oxa  Oya  Oza H1xa */
		*H1y              = vec_perm(*H1y,*H2z,perm);    /* H1ya H1za H2xa H2ya */
		*H2z              = vec_perm(*H2z,tmp4,perm);    /* H2za  Mxa  Mya  Mza */
	}

	/* just splat things... never mind that we fill all cells */
	*Oy               = vec_splat(*Ox,1);
	*Oz               = vec_splat(*Ox,2);
	*H1x              = vec_splat(*Ox,3);
	*H1z              = vec_splat(*H1y,1);
	*H2x              = vec_splat(*H1y,2);
	*H2y              = vec_splat(*H1y,3);
	*Mx               = vec_splat(*H2z,1);
	*My               = vec_splat(*H2z,2);
	*Mz               = vec_splat(*H2z,3);
}



/** Accumulate forces from SIMD to single 3-atom water forces in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 9 atoms
* we can always read/write 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* All four elements of each SIMD vector will be accumulated before adding to
* memory, and the X/Y/Z forces on all atoms are also added (accumulated) to
* the shifted force in memory.
*
* @param water      Address of first force in the water (oxygen x).
* @param fshift     Address of shifted force triplet to increment in memory.
* @param Ox         Output: Oxygen x forces (elements will be accumulated). 
* @param Oy         Output: Oxygen y forces (elements will be accumulated). 
* @param Oz         Output: Oxygen z forces (elements will be accumulated). 
* @param H1x        Output: 1st hydrogen x forces (elements will be accumulated). 
* @param H1y        Output: 1st hydrogen y forces (elements will be accumulated). 
* @param H1z        Output: 1st hydrogen z forces (elements will be accumulated). 
* @param H2x        Output: 2nd hydrogen x forces (elements will be accumulated). 
* @param H2y        Output: 2nd hydrogen y forces (elements will be accumulated). 
* @param H2z        Output: 2nd hydrogen z forces (elements will be accumulated). 
*/
static inline void 
update_i_3atoms_forces(float *water,
                       float *fshift,
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

	add_xyz_to_mem(fshift,Oy);
}


/** Accumulate forces from SIMD to single 4-atom water forces in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read/write 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* All four elements of each SIMD vector will be accumulated before adding to
* memory, and the X/Y/Z forces on all atoms are also added (accumulated) to
* the shifted force in memory.
*
* @param water      Address of first force in the water (oxygen x).
* @param fshift     Address of shifted force triplet to increment in memory.
* @param Ox         Output: Oxygen x forces (elements will be accumulated). 
* @param Oy         Output: Oxygen y forces (elements will be accumulated). 
* @param Oz         Output: Oxygen z forces (elements will be accumulated). 
* @param H1x        Output: 1st hydrogen x forces (elements will be accumulated). 
* @param H1y        Output: 1st hydrogen y forces (elements will be accumulated). 
* @param H1z        Output: 1st hydrogen z forces (elements will be accumulated). 
* @param H2x        Output: 2nd hydrogen x forces (elements will be accumulated). 
* @param H2y        Output: 2nd hydrogen y forces (elements will be accumulated). 
* @param H2z        Output: 2nd hydrogen z forces (elements will be accumulated). 
* @param Mx         Output: Virtual site x forces (elements will be accumulated). 
* @param My         Output: Virtual site y forces (elements will be accumulated). 
* @param Mz         Output: Virtual site z forces (elements will be accumulated). 
*/
static inline void 
update_i_4atoms_forces(float *water,
                       float *fshift,
                       vector float Ox,
                       vector float Oy,
                       vector float Oz,
                       vector float H1x,
                       vector float H1y,
                       vector float H1z,
                       vector float H2x,
                       vector float H2y,
                       vector float H2z,
                       vector float Mx,
                       vector float My,
                       vector float Mz)
{
	vector float l1,l2,l3,l4;
	vector unsigned char perm;
	vector unsigned int mask,ox00;
	vector float nul=vec_zero();

	ox00=(vector unsigned int)vec_splat_s32(0);
	/* load */
	perm              = vec_lvsr( 0, water ); 
	l1                = vec_ld(  0, water);
	l2                = vec_ld( 16, water);
	l3                = vec_ld( 32, water);
	if(((unsigned int)water) & 0xf) {
		l4            = vec_ld( 48, water );
	} else {
		l4            = nul;
	}
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
	Mx               = vec_add(Mx,vec_sld(Mx,Mx,8));     /*  Mx  Mx' - - */
	My               = vec_add(My,vec_sld(My,My,8));     /*  My  My' - - */
	Mz               = vec_add(Mz,vec_sld(Mz,Mz,8));     /*  Mz  Mz' - - */

	Ox               = vec_mergeh(Ox,Oz);        /*  Ox  Oz  Ox'  Oz' */
	Oy               = vec_mergeh(Oy,H1x);       /*  Oy H1x  Oy' H1x' */
	H1y              = vec_mergeh(H1y,H2x);      /* H1y H2x H1y' H2x' */
	H1z              = vec_mergeh(H1z,H2y);      /* H1z H2y H1z' H2y' */
	H2z              = vec_mergeh(H2z,My);       /* H2z  My H2z'  My' */
	Mx               = vec_mergeh(Mx,Mz);        /*  Mx  Mz  Mx'  Mz' */
  
	Ox          = vec_add(Ox,vec_sld(Ox,Ox,8));      /* Ox Oz - - */
	Oy          = vec_add(Oy,vec_sld(Oy,Oy,8));      /* Oy H1x - - */
	H1y         = vec_add(H1y,vec_sld(H1y,H1y,8));   /* H1y H2x - - */
	H1z         = vec_add(H1z,vec_sld(H1z,H1z,8));   /* H1z H2y - - */
	H2z         = vec_add(H2z,vec_sld(H2z,H2z,8));   /* H2z My - - */
	Mx          = vec_add(Mx,vec_sld(Mx,Mx,8));      /* Mx Mz - - */

	Ox          = vec_mergeh(Ox,Oy);   /*  Ox  Oy  Oz H1x */
	H1y         = vec_mergeh(H1y,H1z); /* H1y H1z H2x H2y */
	H2z         = vec_mergeh(H2z,Mx);  /* H2z  Mx  My  Mz  */

	Oy          = vec_sld(nul,Ox,12);  /* 0   Ox  Oy  Oz  */
	Oz          = vec_sld(Ox,H1y,8);   /* -  H1x H1y H1z  */
	H1x         = vec_sld(H1y,H2z,4);  /* -  H2x H2y H2z  */
	/* H2z is already  H2z  Mx  My  Mz  */

	H2x         = vec_perm(nul,Ox,perm);   /* The part to add to l1 */
	H2y         = vec_perm(Ox,H1y,perm);  /* The part to add to l2 */
	Mx          = vec_perm(H1y,H2z,perm); /* The part to add to l3 */
	My          = vec_perm(H2z,nul,perm); /* The part to (maybe) add to l4 */

	H2x         = vec_add(l1,H2x);
	H2y         = vec_add(l2,H2y);
	Mx          = vec_add(l3,Mx);  
   
	H2x         = vec_sel(l1,H2x,mask);

	/* store */
	vec_st(H2x,  0, water);
	vec_st(H2y, 16, water);
	vec_st(Mx, 32, water);

	if(((unsigned int)water) & 0xf) {
		My          = vec_add(l4,My);
		My          = vec_sel(My,l4,mask);
		vec_st(My, 48, water);
	}
  
	Oy          = vec_add(Oy,Oz);
	H1x         = vec_add(H1x,H2z);
	Oy          = vec_add(Oy,H1x);
	Oy          = vec_sld(Oy,nul,4);   /* x   y  z  0 */

	add_xyz_to_mem(fshift,Oy);
}


/** Add forces from SIMD elements to 4 separate 3-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 9 atoms
* we can always read/write 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param address3   Address of first force in the third water (oxygen x).
* @param address4   Address of first force in the fourth water (oxygen x).
* @param Ox         Output: Oxygen x forces for the four waters. 
* @param Oy         Output: Oxygen y forces for the four waters.
* @param Oz         Output: Oxygen z forces for the four waters.
* @param H1x        Output: 1st hydrogen x forces for the four waters.
* @param H1y        Output: 1st hydrogen y forces for the four waters.
* @param H1z        Output: 1st hydrogen z forces for the four waters.
* @param H2x        Output: 2nd hydrogen x forces for the four waters.
* @param H2y        Output: 2nd hydrogen y forces for the four waters.
* @param H2z        Output: 2nd hydrogen z forces for the four waters.
*/
static inline void 
add_force_to_4_3atoms(float *address1,
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
  
	/*  load add and save water 1 */  
	perm             = vec_lvsr( 0, (int *) address1 ); 
	low              = vec_ld(  0, address1);
	medium           = vec_ld( 16, address1);
	high             = vec_ld( 32, address1);
	mask             = vec_perm(ox00,oxFF,perm);
	tmp4              = vec_add(vec_perm(nul,tmp2,perm),low);
	tmp2             = vec_add(vec_perm(tmp2,Oy,perm),medium);
	Oy               = vec_add(vec_perm(Oy,H1z,perm),high);
	vec_st(vec_sel(low,tmp4,mask),  0, address1);
	mask        = vec_sld(ox00,mask,12);
	vec_st(tmp2, 16, address1);
	vec_st(vec_sel(Oy,high,mask), 32, address1);

	/*  load add and save water 2 */  
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

	/*  load add and save water 3 */  
	perm             = vec_lvsr( 0, (int *) address3 );
	low               = vec_ld(  0, address3);
	medium            = vec_ld( 16, address3);
	high              = vec_ld( 32, address3);
	mask              = vec_perm(ox00,oxFF,perm);
	tmp4              = vec_add(vec_perm(nul,tmp1,perm),low);
	tmp1             = vec_add(vec_perm(tmp1,H2x,perm),medium);
	H2x              = vec_add(vec_perm(H2x,tmp3,perm),high);
	vec_st(vec_sel(low,tmp4,mask),  0, address3);
	mask         = vec_sld(ox00,mask,12);
	vec_st(tmp1, 16, address3);
	vec_st(vec_sel(H2x,high,mask), 32, address3);

	/*  load add and save water 4 */  
	perm              = vec_lvsr( 0, (int *) address4 ); 
	low               = vec_ld(  0, address4);
	medium            = vec_ld( 16, address4);
	high              = vec_ld( 32, address4);
	mask              = vec_perm(ox00,oxFF,perm);
	H2x              = vec_add(vec_perm(nul,Ox,perm),low);
	Ox               = vec_add(vec_perm(Ox,H1y,perm),medium);
	H1y              = vec_add(vec_perm(H1y,H2z,perm),high);
	vec_st(vec_sel(low,H2x,mask),  0, address4);
	mask         = vec_sld(ox00,mask,12);
	vec_st(Ox, 16, address4);
	vec_st(vec_sel(H1y,high,mask), 32, address4);
}



/** Add forces from SIMD elements to 4 separate 4-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read/write 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param address3   Address of first force in the third water (oxygen x).
* @param address4   Address of first force in the fourth water (oxygen x).
* @param Ox         Output: Oxygen x forces for the four waters. 
* @param Oy         Output: Oxygen y forces for the four waters.
* @param Oz         Output: Oxygen z forces for the four waters.
* @param H1x        Output: 1st hydrogen x forces for the four waters.
* @param H1y        Output: 1st hydrogen y forces for the four waters.
* @param H1z        Output: 1st hydrogen z forces for the four waters.
* @param H2x        Output: 2nd hydrogen x forces for the four waters.
* @param H2y        Output: 2nd hydrogen y forces for the four waters.
* @param H2z        Output: 2nd hydrogen z forces for the four waters.
* @param Mx         Output: Virtual site x forces for the four waters.
* @param My         Output: Virtual site y forces for the four waters.
* @param Mz         Output: Virtual site z forces for the four waters.
*/
static inline void 
add_force_to_4_4atoms(float *address1,
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
                      vector float H2z,
                      vector float Mx,
                      vector float My,
                      vector float Mz)
{
	vector float l1,l2,l3,l4;
	vector unsigned char perm;
	vector unsigned int mask,oxFF,ox00;
	vector float tmp1,tmp2,tmp3,tmp4;
	vector float nul=vec_zero();

	oxFF=(vector unsigned int)vec_splat_s32(-1);
	ox00=(vector unsigned int)vec_splat_s32(0);
  
	tmp1             = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
	Ox               = vec_mergel(Ox,Oz);    /*  Oxc  Ozc  Oxd  Ozd */
	Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
	Oy               = vec_mergel(Oy,H1x);   /*  Oyc H1xc  Oyd H1xd */
	H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
	H1y              = vec_mergel(H1y,H2x);  /* H1yc H2xc H1yd H2xd */

	H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
	H1z              = vec_mergel(H1z,H2y);  /* H1zc H2yc H1zd H2yd */
	H2y              = vec_mergeh(H2z,My);   /* H2za  Mya H2zb  Myb */
	H2z              = vec_mergel(H2z,My);   /* H2zc  Myc H2zd  Myd */
	My               = vec_mergeh(Mx,Mz);    /*  Mxa  Mza  Mxb  Mzb */
	Mx               = vec_mergel(Mx,Mz);    /*  Mxc  Mzc  Mxd  Mzd */

	tmp2             = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
	Oz               = vec_mergel(tmp1,Oz);   /*  Oxb  Oyb  Ozb H1xb */
	tmp1             = vec_mergeh(Ox,Oy);    /*  Oxc  Oyc  Ozc H1xc */
	Ox               = vec_mergel(Ox,Oy);    /*  Oxd  Oyd  Ozd H1xd */
	Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
	H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
	H2x              = vec_mergeh(H1y,H1z);  /* H1yc H1zc H2xc H2yc */
	H1y              = vec_mergel(H1y,H1z);  /* H1yd H1zd H2xd H2yd */

	H1z              = vec_mergeh(H2y,My);   /* H2za  Mxa  Mya  Mza  */
	H2y              = vec_mergel(H2y,My);   /* H2zb  Mxb  Myb  Mzb  */
	tmp3             = vec_mergeh(H2z,Mx);   /* H2zc  Mxc  Myc  Mzc  */
	H2z              = vec_mergel(H2z,Mx);   /* H2zd  Mxd  Myd  Mzd  */
  
	/*  load add and save water 1 */  
	l1               = vec_ld(  0, address1);
	l2               = vec_ld( 16, address1);
	l3               = vec_ld( 32, address1);
	if((unsigned int)address1 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address1);
		perm             = vec_lvsr( 0, (int *) address1 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp2,perm),l1);
		tmp2             = vec_add(vec_perm(tmp2,Oy,perm),l2);
		Oy               = vec_add(vec_perm(Oy,H1z,perm),l3);
		H1z              = vec_add(vec_perm(H1z,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address1);
		vec_st(tmp2, 16, address1);
		vec_st(Oy,   32, address1);
		vec_st(vec_sel(H1z,l4,mask), 48, address1);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp2),  0, address1);
		vec_st(vec_add(l2,Oy),   16, address1);
		vec_st(vec_add(l3,H1z),  32, address1);
	}


	/*  load add and save water 2 */  
	l1               = vec_ld(  0, address2);
	l2               = vec_ld( 16, address2);
	l3               = vec_ld( 32, address2);
	if((unsigned int)address2 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address2);
		perm             = vec_lvsr( 0, (int *) address2 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,Oz,perm),l1);
		H1z              = vec_add(vec_perm(Oz,H1x,perm),l2);
		Oz               = vec_add(vec_perm(H1x,H2y,perm),l3);
		H1x              = vec_add(vec_perm(H2y,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address2);
		vec_st(H1z, 16, address2);
		vec_st(Oz,  32, address2);
		vec_st(vec_sel(H1x,l4,mask), 48, address2);
	} else {
		/* aligned */
		vec_st(vec_add(l1,Oz),   0, address2);
		vec_st(vec_add(l2,H1x), 16, address2);
		vec_st(vec_add(l3,H2y), 32, address2);
	}

	/*  load add and save water 3 */  
	l1               = vec_ld(  0, address3);
	l2               = vec_ld( 16, address3);
	l3               = vec_ld( 32, address3);
	if((unsigned int)address3 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address3);
		perm             = vec_lvsr( 0, (int *) address3 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp1,perm),l1);
		tmp1             = vec_add(vec_perm(tmp1,H2x,perm),l2);
		H2x              = vec_add(vec_perm(H2x,tmp3,perm),l3);
		tmp3             = vec_add(vec_perm(tmp3,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address3);
		vec_st(tmp1, 16, address3);
		vec_st(H2x,  32, address3);
		vec_st(vec_sel(tmp3,l4,mask), 48, address3);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp1),   0, address3);
		vec_st(vec_add(l2,H2x), 16, address3);
		vec_st(vec_add(l3,tmp3), 32, address3);
	}

	/*  load add and save water 4 */  
	l1               = vec_ld(  0, address4);
	l2               = vec_ld( 16, address4);
	l3               = vec_ld( 32, address4);
	if((unsigned int)address4 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address4);
		perm             = vec_lvsr( 0, (int *) address4 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,Ox,perm),l1);
		Ox               = vec_add(vec_perm(Ox,H1y,perm),l2);
		H1y              = vec_add(vec_perm(H1y,H2z,perm),l3);
		H2z              = vec_add(vec_perm(H2z,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address4);
		vec_st(Ox,  16, address4);
		vec_st(H1y,  32, address4);
		vec_st(vec_sel(H2z,l4,mask), 48, address4);
	} else {
		/* aligned */
		vec_st(vec_add(l1,Ox),   0, address4);
		vec_st(vec_add(l2,H1y), 16, address4);
		vec_st(vec_add(l3,H2z), 32, address4);
	}
}


/** Add forces from SIMD elements to 3 separate 3-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 9 atoms
* we can always read/write 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc. The contents of the fourth element is ignored in
* this routine.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param address3   Address of first force in the third water (oxygen x).
* @param Ox         Output: Oxygen x forces for the three waters. 
* @param Oy         Output: Oxygen y forces for the three waters.
* @param Oz         Output: Oxygen z forces for the three waters.
* @param H1x        Output: 1st hydrogen x forces for the three waters.
* @param H1y        Output: 1st hydrogen y forces for the three waters.
* @param H1z        Output: 1st hydrogen z forces for the three waters.
* @param H2x        Output: 2nd hydrogen x forces for the three waters.
* @param H2y        Output: 2nd hydrogen y forces for the three waters.
* @param H2z        Output: 2nd hydrogen z forces for the three waters.
*/
static inline void 
add_force_to_3_3atoms(float *address1,
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





/** Add forces from SIMD elements to 3 separate 4-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read/write 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param address3   Address of first force in the third water (oxygen x).
* @param Ox         Output: Oxygen x forces for the three waters. 
* @param Oy         Output: Oxygen y forces for the three waters.
* @param Oz         Output: Oxygen z forces for the three waters.
* @param H1x        Output: 1st hydrogen x forces for the three waters.
* @param H1y        Output: 1st hydrogen y forces for the three waters.
* @param H1z        Output: 1st hydrogen z forces for the three waters.
* @param H2x        Output: 2nd hydrogen x forces for the three waters.
* @param H2y        Output: 2nd hydrogen y forces for the three waters.
* @param H2z        Output: 2nd hydrogen z forces for the three waters.
* @param Mx         Output: Virtual site x forces for the three waters.
* @param My         Output: Virtual site y forces for the three waters.
* @param Mz         Output: Virtual site z forces for the three waters.
*/
static inline void
add_force_to_3_4atoms(float *address1,
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
                      vector float H2z,
                      vector float Mx,
                      vector float My,
                      vector float Mz)
{
	vector float l1,l2,l3,l4;
	vector unsigned char perm;
	vector unsigned int mask,oxFF,ox00;
	vector float tmp1,tmp2,tmp3,tmp4;
	vector float nul=vec_zero();

	oxFF=(vector unsigned int)vec_splat_s32(-1);
	ox00=(vector unsigned int)vec_splat_s32(0);
  
	tmp1             = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
	Ox               = vec_mergel(Ox,Oz);    /*  Oxc  Ozc   -    -  */
	Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
	Oy               = vec_mergel(Oy,H1x);   /*  Oyc H1xc   -    -  */
	H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
	H1y              = vec_mergel(H1y,H2x);  /* H1yc H2xc   -    -  */

	H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
	H1z              = vec_mergel(H1z,H2y);  /* H1zc H2yc   -    -  */
	H2y              = vec_mergeh(H2z,My);   /* H2za  Mya H2zb  Myb */
	H2z              = vec_mergel(H2z,My);   /* H2zc  Myc   -    -  */
	My               = vec_mergeh(Mx,Mz);    /*  Mxa  Mza  Mxb  Mzb */
	Mx               = vec_mergel(Mx,Mz);    /*  Mxc  Mzc   -    -  */

	tmp2             = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
	Oz               = vec_mergel(tmp1,Oz);   /*  Oxb  Oyb  Ozb H1xb */
	tmp1             = vec_mergeh(Ox,Oy);    /*  Oxc  Oyc  Ozc H1xc */
	Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
	H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
	H2x              = vec_mergeh(H1y,H1z);  /* H1yc H1zc H2xc H2yc */
	H1z              = vec_mergeh(H2y,My);   /* H2za  Mxa  Mya  Mza  */
	H2y              = vec_mergel(H2y,My);   /* H2zb  Mxb  Myb  Mzb  */
	tmp3             = vec_mergeh(H2z,Mx);   /* H2zc  Mxc  Myc  Mzc  */

  
	/*  load add and save water 1 */  
	l1               = vec_ld(  0, address1);
	l2               = vec_ld( 16, address1);
	l3               = vec_ld( 32, address1);
	if((unsigned int)address1 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address1);
		perm             = vec_lvsr( 0, (int *) address1 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp2,perm),l1);
		tmp2             = vec_add(vec_perm(tmp2,Oy,perm),l2);
		Oy               = vec_add(vec_perm(Oy,H1z,perm),l3);
		H1z              = vec_add(vec_perm(H1z,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address1);
		vec_st(tmp2, 16, address1);
		vec_st(Oy,   32, address1);
		vec_st(vec_sel(H1z,l4,mask), 48, address1);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp2),  0, address1);
		vec_st(vec_add(l2,Oy),   16, address1);
		vec_st(vec_add(l3,H1z),  32, address1);
	}


	/*  load add and save water 2 */  
	l1               = vec_ld(  0, address2);
	l2               = vec_ld( 16, address2);
	l3               = vec_ld( 32, address2);
	if((unsigned int)address2 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address2);
		perm             = vec_lvsr( 0, (int *) address2 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,Oz,perm),l1);
		H1z              = vec_add(vec_perm(Oz,H1x,perm),l2);
		Oz               = vec_add(vec_perm(H1x,H2y,perm),l3);
		H1x              = vec_add(vec_perm(H2y,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address2);
		vec_st(H1z, 16, address2);
		vec_st(Oz,  32, address2);
		vec_st(vec_sel(H1x,l4,mask), 48, address2);
	} else {
		/* aligned */
		vec_st(vec_add(l1,Oz),   0, address2);
		vec_st(vec_add(l2,H1x), 16, address2);
		vec_st(vec_add(l3,H2y), 32, address2);
	}

	/*  load add and save water 3 */  
	l1               = vec_ld(  0, address3);
	l2               = vec_ld( 16, address3);
	l3               = vec_ld( 32, address3);
	if((unsigned int)address3 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address3);
		perm             = vec_lvsr( 0, (int *) address3 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp1,perm),l1);
		tmp1             = vec_add(vec_perm(tmp1,H2x,perm),l2);
		H2x              = vec_add(vec_perm(H2x,tmp3,perm),l3);
		tmp3             = vec_add(vec_perm(tmp3,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address3);
		vec_st(tmp1, 16, address3);
		vec_st(H2x,  32, address3);
		vec_st(vec_sel(tmp3,l4,mask), 48, address3);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp1),   0, address3);
		vec_st(vec_add(l2,H2x), 16, address3);
		vec_st(vec_add(l3,tmp3), 32, address3);
	}
}




/** Add forces from SIMD elements to 2 separate 3-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 9 atoms
* we can always read/write 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc. The contents of the third and fourth element 
* are ignored in this routine.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param Ox         Output: Oxygen x forces for the two waters. 
* @param Oy         Output: Oxygen y forces for the two waters.
* @param Oz         Output: Oxygen z forces for the two waters.
* @param H1x        Output: 1st hydrogen x forces for the two waters.
* @param H1y        Output: 1st hydrogen y forces for the two waters.
* @param H1z        Output: 1st hydrogen z forces for the two waters.
* @param H2x        Output: 2nd hydrogen x forces for the two waters.
* @param H2y        Output: 2nd hydrogen y forces for the two waters.
* @param H2z        Output: 2nd hydrogen z forces for the two waters.
*/
static inline void
add_force_to_2_3atoms(float *address1,
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
	vector float tmp1,tmp2;

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


/** Add forces from SIMD elements to 2 separate 4-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read/write 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water, etc.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param address2   Address of first force in the second water (oxygen x).
* @param Ox         Output: Oxygen x forces for the two waters. 
* @param Oy         Output: Oxygen y forces for the two waters.
* @param Oz         Output: Oxygen z forces for the two waters.
* @param H1x        Output: 1st hydrogen x forces for the two waters.
* @param H1y        Output: 1st hydrogen y forces for the two waters.
* @param H1z        Output: 1st hydrogen z forces for the two waters.
* @param H2x        Output: 2nd hydrogen x forces for the two waters.
* @param H2y        Output: 2nd hydrogen y forces for the two waters.
* @param H2z        Output: 2nd hydrogen z forces for the two waters.
* @param Mx         Output: Virtual site x forces for the two waters.
* @param My         Output: Virtual site y forces for the two waters.
* @param Mz         Output: Virtual site z forces for the two waters.
*/
static inline void 
add_force_to_2_4atoms(float *address1,
                      float *address2,
                      vector float Ox,
                      vector float Oy,
                      vector float Oz,
                      vector float H1x,
                      vector float H1y,
                      vector float H1z,
                      vector float H2x,
                      vector float H2y,
                      vector float H2z,
                      vector float Mx,
                      vector float My,
                      vector float Mz)
{
	vector float l1,l2,l3,l4;
	vector unsigned char perm;
	vector unsigned int mask,oxFF,ox00;
	vector float tmp1,tmp2,tmp4;

	vector float nul=vec_zero();
	oxFF=(vector unsigned int)vec_splat_s32(-1);
	ox00=(vector unsigned int)vec_splat_s32(0);

	tmp1             = vec_mergeh(Ox,Oz);    /*  Oxa  Oza  Oxb  Ozb */
	Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa  Oyb H1xb */
	H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa H1yb H2xb */
	H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya H1zb H2yb */
	H2y              = vec_mergeh(H2z,My);   /* H2za  Mya  H2zb Myb */
	My               = vec_mergeh(Mx,Mz);    /*  Mxa  Mza  Mxb  Mzb */

	tmp2             = vec_mergeh(tmp1,Oz);  /*  Oxa  Oya  Oza H1xa */
	Oz               = vec_mergel(tmp1,Oz);  /*  Oxb  Oyb  Ozb H1xb */
	Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
	H1x              = vec_mergel(H1x,H2x);  /* H1yb H1zb H2xb H2yb */
	H1z              = vec_mergeh(H2y,My);   /* H2za  Mxa  Mya  Mza */
	H2y              = vec_mergel(H2y,My);   /* H2zb  Mxb  Myb  Mzb */

 
	/*  load add and save water 1 */  
	l1               = vec_ld(  0, address1);
	l2               = vec_ld( 16, address1);
	l3               = vec_ld( 32, address1);
	if((unsigned int)address1 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address1);
		perm             = vec_lvsr( 0, (int *) address1 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp2,perm),l1);
		tmp2             = vec_add(vec_perm(tmp2,Oy,perm),l2);
		Oy               = vec_add(vec_perm(Oy,H1z,perm),l3);
		H1z              = vec_add(vec_perm(H1z,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address1);
		vec_st(tmp2, 16, address1);
		vec_st(Oy,   32, address1);
		vec_st(vec_sel(H1z,l4,mask), 48, address1);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp2),  0, address1);
		vec_st(vec_add(l2,Oy),   16, address1);
		vec_st(vec_add(l3,H1z),  32, address1);
	}


	/*  load add and save water 2 */  
	l1               = vec_ld(  0, address2);
	l2               = vec_ld( 16, address2);
	l3               = vec_ld( 32, address2);
	if((unsigned int)address2 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address2);
		perm             = vec_lvsr( 0, (int *) address2 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,Oz,perm),l1);
		H1z              = vec_add(vec_perm(Oz,H1x,perm),l2);
		Oz               = vec_add(vec_perm(H1x,H2y,perm),l3);
		H1x              = vec_add(vec_perm(H2y,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address2);
		vec_st(H1z, 16, address2);
		vec_st(Oz,  32, address2);
		vec_st(vec_sel(H1x,l4,mask), 48, address2);
	} else {
		/* aligned */
		vec_st(vec_add(l1,Oz),   0, address2);
		vec_st(vec_add(l2,H1x), 16, address2);
		vec_st(vec_add(l3,H2y), 32, address2);
	}
}



/** Add forces from first SIMD elements to a 3-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 9 atoms
* we can always read/write 3 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water. The remaining elements are ignored in this routine.
*
* @param address1   Address of first force in the first water (oxygen x).
* @param Ox         Output: Oxygen x forces for the water. 
* @param Oy         Output: Oxygen y forces for the water.
* @param Oz         Output: Oxygen z forces for the water.
* @param H1x        Output: 1st hydrogen x forces for the water.
* @param H1y        Output: 1st hydrogen y forces for the water.
* @param H1z        Output: 1st hydrogen z forces for the water.
* @param H2x        Output: 2nd hydrogen x forces for the water.
* @param H2y        Output: 2nd hydrogen y forces for the water.
* @param H2z        Output: 2nd hydrogen z forces for the water.
*/
static inline void
add_force_to_1_3atoms(float *address1,
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
	vector float tmp1,tmp2;

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


/** Add forces from first SIMD elements to a 4-atom waters in memory.
*
* We use our knowledge about water to improve the load/store performance.
* The allocated memory is always 8-byte aligned, so with 12 atoms
* we can always read/write 4 16-byte chunks (rounded down by vec_ld)
* without going outside the current page (which might incur a page fault). 
* The extra elements before and after the water dont matter as long as we dont
* try to modify them.
*
* The SIMD vectors will be transposed, so that the first element is added
* to the first water. The remaining elements are ignored in this routine.
*
* @param address1   Address of first force in the water (oxygen x).
* @param Ox         Output: Oxygen x forces for the water. 
* @param Oy         Output: Oxygen y forces for the water.
* @param Oz         Output: Oxygen z forces for the water.
* @param H1x        Output: 1st hydrogen x forces for the water.
* @param H1y        Output: 1st hydrogen y forces for the water.
* @param H1z        Output: 1st hydrogen z forces for the water.
* @param H2x        Output: 2nd hydrogen x forces for the water.
* @param H2y        Output: 2nd hydrogen y forces for the water.
* @param H2z        Output: 2nd hydrogen z forces for the water.
* @param Mx         Output: Virtual site x forces for the water.
* @param My         Output: Virtual site y forces for the water.
* @param Mz         Output: Virtual site z forces for the water.
*/
static inline void
add_force_to_1_4atoms(float *address1,
                                         vector float Ox,
                                         vector float Oy,
                                         vector float Oz,
                                         vector float H1x,
                                         vector float H1y,
                                         vector float H1z,
                                         vector float H2x,
                                         vector float H2y,
                                         vector float H2z,
                                         vector float Mx,
                                         vector float My,
                                         vector float Mz)
{
	vector float l1,l2,l3,l4;
	vector unsigned char perm;
	vector unsigned int mask,oxFF,ox00;
	vector float tmp1,tmp2,tmp4;

	vector float nul=vec_zero();
	oxFF=(vector unsigned int)vec_splat_s32(-1);
	ox00=(vector unsigned int)vec_splat_s32(0);

	tmp1              = vec_mergeh(Ox,Oz);    /*  Oxa  Oza   ?    ?  */
	Oz               = vec_mergeh(Oy,H1x);   /*  Oya H1xa   ?    ?  */
	H1x              = vec_mergeh(H1y,H2x);  /* H1ya H2xa   ?    ?  */
	H2x              = vec_mergeh(H1z,H2y);  /* H1za H2ya   ?    ?  */
	H2y              = vec_mergeh(H2z,My);   /* H2za  Mya   ?    0  */
	My               = vec_mergeh(Mx,Mz);    /*  Mxa  Mza    ?  ? */

	tmp2              = vec_mergeh(tmp1,Oz);   /*  Oxa  Oya  Oza H1xa */
	Oy               = vec_mergeh(H1x,H2x);  /* H1ya H1za H2xa H2ya */
	H1z              = vec_mergeh(H2y,My);   /* H2za Mxa Mya Mza */
  
	/*  load, add and save water 1 */  
	l1               = vec_ld(  0, address1);
	l2               = vec_ld( 16, address1);
	l3               = vec_ld( 32, address1);
	if((unsigned int)address1 & 0xf) {
		/* unaligned */
		l4               = vec_ld( 48, address1);
		perm             = vec_lvsr( 0, (int *) address1 ); 
		mask             = vec_perm(ox00,oxFF,perm);
		tmp4             = vec_add(vec_perm(nul,tmp2,perm),l1);
		tmp2             = vec_add(vec_perm(tmp2,Oy,perm),l2);
		Oy               = vec_add(vec_perm(Oy,H1z,perm),l3);
		H1z              = vec_add(vec_perm(H1z,nul,perm),l4);  
		vec_st(vec_sel(l1,tmp4,mask),  0, address1);
		vec_st(tmp2, 16, address1);
		vec_st(Oy,   32, address1);
		vec_st(vec_sel(H1z,l4,mask), 48, address1);
	} else {
		/* aligned */
		vec_st(vec_add(l1,tmp2),  0, address1);
		vec_st(vec_add(l2,Oy),   16, address1);
		vec_st(vec_add(l3,H1z),  32, address1);
	}
}



/** Load Coulomb-only potential and force from Coulomb-only table for 4 values.
 *
 *  This routine returns the potential and force values from the table
 *  provided, for the table distance rtab.
 *
 *  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
 *                   format (see manual) with four floats per table point.
 *  @param rtab      Table lookup value in floating point. Note that this is 
 *                   NOT the interaction distance, but the distance multiplied
 *                   with the table scale, e.g. r*500 if you have 500 table
 *                   points per nm.
 *  @param VV        Output: potential for the four lookup values.
 *  @param FF        Output: force (negative potential derivative) for the
 *                   four lookup values.
 */
static inline void 
do_4_ctable_coul(float *VFtab,
                 vector float rtab,
                 vector float *VV,
                 vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
	int idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];

	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
	tab4=vec_ld(0, VFtab+idx4);

  
	/* table data is aligned */
	transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);
  
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


/** Load Coulomb-only potential and force from Coul+LJ table for 4 values.
*
*  This routine returns the coulomb potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the four lookup values.
*  @param FF        Output: Coulomb force (negative potential derivative) for 
*                   the four lookup values.
*/
static inline void do_4_ljctable_coul(float *VFtab,
									  vector float rtab,
									  vector float *VV,
									  vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
	int idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
	tab4=vec_ld(0, VFtab+idx4);
  
	transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


/** Load LJ-only potential and force from LJ-only table for 4 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the four lookup values.
*/
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
	tabd3  = vec_ld( 0, VFtab+idx3);
	tabr3  = vec_ld(16, VFtab+idx3);
	tabd4  = vec_ld( 0, VFtab+idx4);
	tabr4  = vec_ld(16, VFtab+idx4);
      
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


/** Load LJ-only potential and force from Coul+LJ table for 4 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the four lookup values.
*/
static inline void do_4_ljctable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];

	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
	tabd4  = vec_ld(16, VFtab+idx4);
	tabr4  = vec_ld(32, VFtab+idx4);
    
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


/** Load Coul+LJ potential and force from Coul+LJ table for 4 values.
*
*  This routine returns both Coulomb and LJ potential and force values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the four lookup values.
*  @param FFcoul    Output: Coulomb force (negative potential derivative) for 
*                   the four lookup values.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the four lookup values.
*/
static inline void 
do_4_ljctable_coul_and_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
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



/** Load Coulomb-only potential and force from Coulomb-only table for 3 values.
*
*  This routine returns the potential and force values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the three lookup values
*                   (highest element will not be used).
*  @param FF        Output: force (negative potential derivative) for the
*                   three lookup values (highest element will not be used).
*/
static inline void
do_3_ctable_coul(float *VFtab,
                 vector float rtab,
                 vector float *VV,
                 vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
	int idx1,idx2,idx3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
  
	/* table data is aligned */
	transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);
  

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/** Load Coulomb-only potential and force from Coul+LJ table for 3 values.
*
*  This routine returns the coulomb potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest element will not be used.
*  @param FF        Output: Coulomb force (negative potential derivative) for 
*                   the lookup values. Highest element will not be used.
*/
static inline void 
do_3_ljctable_coul(float *VFtab,
                   vector float rtab,
                   vector float *VV,
                   vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
	int idx1,idx2,idx3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
    
	/* table data is aligned */
	transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


/** Load LJ-only potential and force from LJ-only table for 3 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the three lookup values (highest element unused).
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the three lookup values (highest element unused).
*/
static inline void
do_3_ljtable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
	tabd3  = vec_ld( 0, VFtab+idx3);
	tabr3  = vec_ld(16, VFtab+idx3);
      
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



/** Load LJ-only potential and force from Coul+LJ table for 3 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the three lookup values in the lowest three elements.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the three lookup values in the lowest three elements.
*/
static inline void
do_3_ljctable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
      
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


/** Load Coul+LJ potential and force from Coul+LJ table for 3 values.
*
*  This routine returns both Coulomb and LJ potential and force values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest three elements.
*  @param FFcoul    Output: Coulomb force (negative potential derivative) for 
*                   the lookup values in the lowest three elements.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest three elements.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup values in the lowest three elements.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest three elements.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the lookup values in the lowest three elements.
*/
static inline void 
do_3_ljctable_coul_and_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabc2  = vec_ld( 0, VFtab+idx2);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabc3  = vec_ld( 0, VFtab+idx3);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
  
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




/** Load Coulomb-only potential and force from Coulomb-only table for 2 values.
*
*  This routine returns the potential and force values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the lookup values in elements 0,1.
*  @param FF        Output: force (negative potential derivative) for the
*                   lookup values (highest 2 elements will not be used).
*/
static inline void 
do_2_ctable_coul(float *VFtab,
                 vector float rtab,
                 vector float *VV,
                 vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2;
	int idx1,idx2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
    vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
  
	/* table data is aligned */
	transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/** Load Coulomb-only potential and force from Coul+LJ table for 2 values.
*
*  This routine returns the coulomb potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest two elements will not be used.
*  @param FF        Output: Coulomb force (negative potential derivative) for 
*                   the lookup values. Highest two elements will not be used.
*/
static inline void 
do_2_ljctable_coul(float *VFtab,
                   vector float rtab,
                   vector float *VV,
                   vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2;
	int idx1,idx2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
    
	/* table data is aligned */
	transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


/** Load LJ-only potential and force from LJ-only table for 2 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values
*                   in elements 0,1.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup values in elements 0,1.
*  @param VVrep     Output: LJ repulsive potential for the lookup values
*                   in elements 0,1.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the two lookup values (highest two elements unused).
*/
static inline void
do_2_ljtable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
      
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


/** Load LJ-only potential and force from Coul+LJ table for 2 values.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values
*                   in elements 0,1.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the two lookup values in the lowest two elements.
*  @param VVrep     Output: LJ repulsive potential for the lookup values
*                   in elements 0,1.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the two lookup values in the lowest two elements.
*/
static inline void 
do_2_ljctable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
      
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


/** Load Coul+LJ potential and force from Coul+LJ table for 2 values.
*
*  This routine returns both Coulomb and LJ potential and force values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest two elements.
*  @param FFcoul    Output: Coulomb force (negative potential derivative) for 
*                   the lookup values in the lowest two elements.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest two elements.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup values in the lowest two elements.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest two elements.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the lookup values in the lowest two elements.
*/
static inline void 
do_2_ljctable_coul_and_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabc2  = vec_ld( 0, VFtab+idx2);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
  
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



/** Load Coulomb-only potential and force from Coulomb-only table for 1 value.
*
*  This routine returns the potential and force values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the lookup value 
*                   (only lowest element will be used).
*  @param FF        Output: force (negative potential derivative) for the
*                   lookup value (only lowest element will be used).
*/
static inline void
do_1_ctable_coul(float *VFtab,
                 vector float rtab,
                 vector float *VV,
                 vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1;
	int idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);

	/* table data is aligned */
	transpose_1_to_4(tab1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}

/** Load Coulomb-only potential and force from Coul+LJ table for 1 value.
*
*  This routine returns the coulomb potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest three elements will not be used.
*  @param FF        Output: Coulomb force (negative potential derivative) for 
*                   the lookup values. Highest three elements will not be used.
*/
static inline void do_1_ljctable_coul(float *VFtab,
									  vector float rtab,
									  vector float *VV,
									  vector float *FF)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1;
	int idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
    
	/* table data is aligned */
	transpose_1_to_4(tab1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
	F         = vec_madd(G,eps,F);           /* Fp + Geps */
	*FF       = vec_madd(vec_two(),H,F);     /* Fp + Geps + 2.0*Heps2 */
}


/** Load LJ-only potential and force from LJ-only table for 1 value.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup value. (only uses lowest element).
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the lookup value (only uses lowest element).
*/
static inline void
do_1_ljtable_lj(float *VFtab,
                vector float rtab,
                vector float *VVdisp,
                vector float *FFdisp,
                vector float *VVrep,
                vector float *FFrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1;
	vector float tabr1;
	int    idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
      
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



/** Load LJ-only potential and force from Coul+LJ table for 1 value.
*
*  This routine returns the LJ potential and force values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup value in the lowest element.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the lookup value in the lowest element.
*/
static inline void 
do_1_ljctable_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
      
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


/** Load Coul+LJ potential and force from Coul+LJ table for 1 value.
*
*  This routine returns both Coulomb and LJ potential and force values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest element.
*  @param FFcoul    Output: Coulomb force (negative potential derivative) for 
*                   the lookup values in the lowest element.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest element.
*  @param FFdisp    Output: LJ dispersion force (negative potential derivative) 
*                   for the lookup values in the lowest element.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest element.
*  @param FFrep     Output: LJ repulsive force (negative potential derivative) 
*                   for the lookup values in the lowest element.
*/
static inline void 
do_1_ljctable_coul_and_lj(float *VFtab,
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

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);

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




/** Load Coulomb-only potential from Coulomb-only table for 4 values.
*
*  This routine returns the potential values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the four lookup values.
*/
static inline void 
do_vonly_4_ctable_coul(float *VFtab,
                       vector float rtab,
                       vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
	int idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
	tab4=vec_ld(0, VFtab+idx4);
  
	/* table data is aligned */
	transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);
  
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}



/** Load Coulomb-only potential from Coul+LJ table for 4 values.
*
*  This routine returns the coulomb potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the four lookup values.
*/
static inline void 
do_vonly_4_ljctable_coul(float *VFtab,
                         vector float rtab,
                         vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3,tab4;
	int idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
	tab4=vec_ld(0, VFtab+idx4);
  
	transpose_4_to_4(tab1,tab2,tab3,tab4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from LJ-only table for 4 values.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void
do_vonly_4_ljtable_lj(float *VFtab,
                      vector float rtab,
                      vector float *VVdisp,
                      vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabd2,tabd3,tabd4;
	vector float tabr1,tabr2,tabr3,tabr4;

	int    idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
	tabd3  = vec_ld( 0, VFtab+idx3);
	tabr3  = vec_ld(16, VFtab+idx3);
	tabd4  = vec_ld( 0, VFtab+idx4);
	tabr4  = vec_ld(16, VFtab+idx4);
    
	transpose_4_to_4(tabd1,tabd2,tabd3,tabd4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_4_to_4(tabr1,tabr2,tabr3,tabr4,&Y,&F,&G,&H);
 
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from Coul+LJ table for 4 values.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void 
do_vonly_4_ljctable_lj(float *VFtab,
                       vector float rtab,
                       vector float *VVdisp,
                       vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabd2,tabd3,tabd4;
	vector float tabr1,tabr2,tabr3,tabr4;

	int    idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Tab must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
	tabd4  = vec_ld(16, VFtab+idx4);
	tabr4  = vec_ld(32, VFtab+idx4);
      
	transpose_4_to_4(tabd1,tabd2,tabd3,tabd4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_4_to_4(tabr1,tabr2,tabr3,tabr4,&Y,&F,&G,&H);
 
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coul+LJ potential from Coul+LJ table for 4 values.
*
*  This routine returns both Coulomb and LJ potential values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the four lookup values.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void 
do_vonly_4_ljctable_coul_and_lj(float *VFtab,
                                vector float rtab,
                                vector float *VVcoul,
                                vector float *VVdisp,
                                vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
	vector float tabd3,tabr3,tabc3,tabd4,tabr4,tabc4;
	int idx1,idx2,idx3,idx4;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
	idx4     = conv.i[3];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
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
  
	transpose_4_to_4(tabc1,tabc2,tabc3,tabc4,&Y,&F,&G,&H);
   
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_4_to_4(tabd1,tabd2,tabd3,tabd4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_4_to_4(tabr1,tabr2,tabr3,tabr4,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coulomb-only potential from Coulomb-only table for 3 values.
*
*  This routine returns the potential values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the three lookup values
*                   (highest element will not be used).
*/
static inline void 
do_vonly_3_ctable_coul(float *VFtab,
                       vector float rtab,
                       vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
	int idx1,idx2,idx3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
  
	/* table data is aligned */
	transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);
  

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coulomb-only potential from Coul+LJ table for 3 values.
*
*  This routine returns the coulomb potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest element will not be used.
*/
static inline void do_vonly_3_ljctable_coul(float *VFtab,
											vector float rtab,
											vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2,tab3;
	int idx1,idx2,idx3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
	tab3=vec_ld(0, VFtab+idx3);
    
	/* table data is aligned */
	transpose_3_to_4(tab1,tab2,tab3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}



/** Load LJ-only potential from LJ-only table for 3 values.
 *
 *  This routine returns the LJ potential values from the table
 *  provided, for the table distance rtab. 
 *
 *  @param VFtab     Table data, should be a LJ-only table in Gromacs
 *                   format (see manual) with 8 floats per table point.
 *  @param rtab      Table lookup value in floating point. Note that this is 
 *                   NOT the interaction distance, but the distance multiplied
 *                   with the table scale, e.g. r*500 if you have 500 table
 *                   points per nm.
 *  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
 *  @param VVrep     Output: LJ repulsive potential for the four lookup values.
 */
static inline void 
do_vonly_3_ljtable_lj(float *VFtab,
                      vector float rtab,
                      vector float *VVdisp,
                      vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1,idx2,idx3;
	vector float tabd1,tabd2,tabd3;
	vector float tabr1,tabr2,tabr3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    }  conv;
        
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];

	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
	tabd3  = vec_ld( 0, VFtab+idx3);
	tabr3  = vec_ld(16, VFtab+idx3);
      
	transpose_3_to_4(tabd1,tabd2,tabd3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_3_to_4(tabr1,tabr2,tabr3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from Coul+LJ table for 3 values.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void 
do_vonly_3_ljctable_lj(float *VFtab,
                       vector float rtab,
                       vector float *VVdisp,
                       vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1,idx2,idx3;
	vector float tabd1,tabd2,tabd3;
	vector float tabr1,tabr2,tabr3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
      
	transpose_3_to_4(tabd1,tabd2,tabd3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_3_to_4(tabr1,tabr2,tabr3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coul+LJ potential from Coul+LJ table for 3 values.
*
*  This routine returns both Coulomb and LJ potential values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest three elements.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest three elements.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest three elements.
*/
static inline void 
do_vonly_3_ljctable_coul_and_lj(float *VFtab,
                                vector float rtab,
                                vector float *VVcoul,
                                vector float *VVdisp,
                                vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
	vector float tabd3,tabr3,tabc3;
	int idx1,idx2,idx3;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
	idx3     = conv.i[2];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabc2  = vec_ld( 0, VFtab+idx2);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
	tabc3  = vec_ld( 0, VFtab+idx3);
	tabd3  = vec_ld(16, VFtab+idx3);
	tabr3  = vec_ld(32, VFtab+idx3);
  
	transpose_3_to_4(tabc1,tabc2,tabc3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_3_to_4(tabd1,tabd2,tabd3,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_3_to_4(tabr1,tabr2,tabr3,&Y,&F,&G,&H);
 
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}




/** Load Coulomb-only potential from Coulomb-only table for 2 values.
*
*  This routine returns the potential values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the lookup values in elements 0,1.
*/
static inline void 
do_vonly_2_ctable_coul(float *VFtab,
                       vector float rtab,
                       vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2;
	int idx1,idx2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
  
	/* table data is aligned */
	transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}

/** Load Coulomb-only potential from Coul+LJ table for 2 values.
*
*  This routine returns the coulomb potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest two elements will not be used.
*/
static inline void 
do_vonly_2_ljctable_coul(float *VFtab,
                         vector float rtab,
                         vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1,tab2;
	int idx1,idx2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
	tab2=vec_ld(0, VFtab+idx2);
    
	/* table data is aligned */
	transpose_2_to_4(tab1,tab2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from LJ-only table for 2 values.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values
*                   in elements 0,1.
*  @param VVrep     Output: LJ repulsive potential for the lookup values
*                   in elements 0,1.
*/
static inline void 
do_vonly_2_ljtable_lj(float *VFtab,
                      vector float rtab,
                      vector float *VVdisp,
                      vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1,idx2;
	vector float tabd1,tabd2;
	vector float tabr1,tabr2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
	tabd2  = vec_ld( 0, VFtab+idx2);
	tabr2  = vec_ld(16, VFtab+idx2);
      
	transpose_2_to_4(tabd1,tabd2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_2_to_4(tabr1,tabr2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from Coul+LJ table for 2 values.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values
*                   in elements 0,1.
*  @param VVrep     Output: LJ repulsive potential for the lookup values
*                   in elements 0,1.
*/
static inline void 
do_vonly_2_ljctable_lj(float *VFtab,
                       vector float rtab,
                       vector float *VVdisp,
                       vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1,idx2;
	vector float tabd1,tabd2;
	vector float tabr1,tabr2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
      
	transpose_2_to_4(tabd1,tabd2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_2_to_4(tabr1,tabr2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coul+LJ potential from Coul+LJ table for 2 values.
*
*  This routine returns both Coulomb and LJ potential values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest two elements.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest two elements.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest two elements.
*/
static inline void 
do_vonly_2_ljctable_coul_and_lj(float *VFtab,
                                vector float rtab,
                                vector float *VVcoul,
                                vector float *VVdisp,
                                vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabr1,tabc1,tabd2,tabr2,tabc2;
	int idx1,idx2;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
	idx2     = conv.i[1];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
	tabc2  = vec_ld( 0, VFtab+idx2);
	tabd2  = vec_ld(16, VFtab+idx2);
	tabr2  = vec_ld(32, VFtab+idx2);
  
	transpose_2_to_4(tabc1,tabc2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_2_to_4(tabd1,tabd2,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_2_to_4(tabr1,tabr2,&Y,&F,&G,&H);
 
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}



/** Load Coulomb-only potential from Coulomb-only table for 1 value.
*
*  This routine returns the potential values from the table
*  provided, for the table distance rtab.
*
*  @param VFtab     Table data, should be a Coulomb-only table in Gromacs
*                   format (see manual) with four floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: potential for the lookup value 
*                   (only lowest element will be used).
*/
static inline void 
do_vonly_1_ctable_coul(float *VFtab,
                       vector float rtab,
                       vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1;
	int idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);

	/* table data is aligned */
	transpose_1_to_4(tab1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}



/** Load Coulomb-only potential from Coul+LJ table for 1 value.
*
*  This routine returns the coulomb potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in water-water loops, where the table might be a combination of 
*  Coulomb and LJ, but most interactions (8 of 9) are coulomb-only.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VV        Output: Coulomb potential for the lookup values.
*                   The highest three elements will not be used.
*/
static inline void 
do_vonly_1_ljctable_coul(float *VFtab,
                         vector float rtab,
                         vector float *VV)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2,tab1;
	int idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tab1=vec_ld(0, VFtab+idx1);
    
	/* table data is aligned */
	transpose_1_to_4(tab1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VV       = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from LJ-only table for 1 value.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a LJ-only table in Gromacs
*                   format (see manual) with 8 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void 
do_vonly_1_ljtable_lj(float *VFtab,
                      vector float rtab,
                      vector float *VVdisp,
                      vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1;
	vector float tabd1;
	vector float tabr1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_sl(vidx,vec_splat_u32(3));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld( 0, VFtab+idx1);
	tabr1  = vec_ld(16, VFtab+idx1);
      
	transpose_1_to_4(tabd1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_1_to_4(tabr1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load LJ-only potential from Coul+LJ table for 1 value.
*
*  This routine returns the LJ potential values from the table
*  provided, for the table distance rtab. This routine is intended to be
*  used in Tip4p water-water loops, where the virtual site only has LJ.
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVdisp    Output: LJ dispersion potential for the four lookup values.
*  @param VVrep     Output: LJ repulsive potential for the four lookup values.
*/
static inline void 
do_vonly_1_ljctable_lj(float *VFtab,
                       vector float rtab,
                       vector float *VVdisp,
                       vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	int    idx1;
	vector float tabd1;
	vector float tabr1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
      
	transpose_1_to_4(tabd1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_1_to_4(tabr1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}


/** Load Coul+LJ potential from Coul+LJ table for 1 value.
*
*  This routine returns both Coulomb and LJ potential values from the 
*  table provided, for the table distance rtab. 
*
*  @param VFtab     Table data, should be a Coul+LJ table in Gromacs
*                   format (see manual) with 12 floats per table point.
*  @param rtab      Table lookup value in floating point. Note that this is 
*                   NOT the interaction distance, but the distance multiplied
*                   with the table scale, e.g. r*500 if you have 500 table
*                   points per nm.
*  @param VVcoul    Output: Coulomb potential for the lookup values
*                   in the lowest element.
*  @param VVdisp    Output: LJ dispersion potential for the lookup values in
*                   the lowest element.
*  @param VVrep     Output: LJ repulsive potential for the lookup values in
*                   the lowest element.
*/
static inline void 
do_vonly_1_ljctable_coul_and_lj(float *VFtab,
                                vector float rtab,
                                vector float *VVcoul,
                                vector float *VVdisp,
                                vector float *VVrep)
{
	vector signed int vidx;
	vector float Y,F,G,H,eps,eps2;
	vector float tabd1,tabr1,tabc1;
	int idx1;

    /* necessary to avoid aliasing optimization problems */
    union
    {
        vector signed int   v;
        int                 i[4];
    } conv;
	
	vidx     = vec_cts(rtab,0); 
	vidx     = vec_add(vidx,vec_sl(vidx,vec_splat_u32(1))); /* multiply by 3 */
	vidx     = vec_sl(vidx,vec_splat_u32(2));

    conv.v   = vidx;
    
	idx1     = conv.i[0];
    
	eps      = vec_sub(rtab,vec_floor(rtab));
	eps2     = vec_madd(eps,eps,vec_zero());

	/* Table must be aligned */
	tabc1  = vec_ld( 0, VFtab+idx1);
	tabd1  = vec_ld(16, VFtab+idx1);
	tabr1  = vec_ld(32, VFtab+idx1);
  
	transpose_1_to_4(tabc1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVcoul   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_1_to_4(tabd1,&Y,&F,&G,&H);

	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVdisp   = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */

	transpose_1_to_4(tabr1,&Y,&F,&G,&H);
 
	F         = vec_madd(G,eps,F);           /* F + Geps   */
	H         = vec_madd(H,eps2,vec_zero()); /* Heps2 */
	F         = vec_add(F,H);                /* F + Geps + Heps2 (=Fp) */
	*VVrep    = vec_madd(eps,F,Y);           /* VV = Y + eps*Fp */
}




/** Calculate 1/x.
*
* @param rsq   Input argument.
*
* @return      1.0/rsq.
*/
static inline vector float do_recip(vector float rsq)
{
	vector float tmp,lu;

	lu        = vec_re(rsq);
	tmp       = vec_nmsub(rsq,lu,vec_two());
	return vec_madd(lu,tmp,vec_zero());
}


/** Calculate an inverse square root.
 *
 * @param rsq   Input argument.
 *
 * @return      1.0/sqrt(rsq)
 */
static inline vector float 
do_invsqrt(vector float rsq)
{
	vector float lu,tmpA,tmpB;

	lu        = vec_rsqrte(rsq);
	tmpA      = vec_madd(lu,lu,vec_zero());
	tmpB      = vec_madd(lu,vec_half(),vec_zero());
	tmpA      = vec_nmsub(rsq,tmpA,vec_one());
	return vec_madd(tmpA,tmpB,lu);
}


/** Calculate 3 inverse square roots.
*
* @param rsq1   Input argument.
* @param rsq2   Input argument.
* @param rsq3   Input argument.
* @param rinv1  Output: 1.0/sqrt(rsq1)
* @param rinv2  Output: 1.0/sqrt(rsq2)
* @param rinv3  Output: 1.0/sqrt(rsq3)
*/
static inline void 
do_3_invsqrt(vector float rsq1,
             vector float rsq2,
             vector float rsq3,
             vector float *rinv1,
             vector float *rinv2,
             vector float *rinv3)
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
	*rinv1    = vec_madd(tmpA1,tmpB1,lu1);
	*rinv2    = vec_madd(tmpA2,tmpB2,lu2);
	*rinv3    = vec_madd(tmpA3,tmpB3,lu3);
}

/** Calculate 9 inverse square roots (unrolled).
 *
 * @param rsq1   Input argument.
 * @param rsq2   Input argument.
 * @param rsq3   Input argument.
 * @param rsq4   Input argument.
 * @param rsq5   Input argument.
 * @param rsq6   Input argument.
 * @param rsq7   Input argument.
 * @param rsq8   Input argument.
 * @param rsq9   Input argument.
 * @param rinv1  Output: 1.0/sqrt(rsq1)
 * @param rinv2  Output: 1.0/sqrt(rsq2)
 * @param rinv3  Output: 1.0/sqrt(rsq3)
 * @param rinv4  Output: 1.0/sqrt(rsq4)
 * @param rinv5  Output: 1.0/sqrt(rsq5)
 * @param rinv6  Output: 1.0/sqrt(rsq6)
 * @param rinv7  Output: 1.0/sqrt(rsq7)
 * @param rinv8  Output: 1.0/sqrt(rsq8)
 * @param rinv9  Output: 1.0/sqrt(rsq9)
 */
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




/** Clear the highest element in a SIMD variable.
*
*  @param v  SIMD variable to clear element 3 of.
*/
static inline void 
zero_highest_element_in_vector(vector float *v)
{
	vector signed int zero = (vector signed int) vec_zero();
  
	*v = (vector float)vec_sel((vector signed int)*v,zero,
							   (vector unsigned int)vec_sld(zero,
															vec_splat_s32(-1),
															4));
}

/** Clear the upper half in a SIMD variable.
*
*  @param v  SIMD variable to clear elements 2-3 of.
*/
static inline void 
zero_highest_2_elements_in_vector(vector float *v)
{
	vector signed int zero = (vector signed int) vec_zero();
  
	*v = (vector float)vec_sel((vector signed int)*v,zero,
							   (vector unsigned int)vec_sld(zero,
															vec_splat_s32(-1),
															8));
}


/** Clear the highest three elements in a SIMD variable.
*
*  @param v  SIMD variable to clear elements 1-3 of.
*/
static inline void
zero_highest_3_elements_in_vector(vector float *v)
{
	vector signed int zero = (vector signed int) vec_zero();
  
	*v = (vector float)vec_sel((vector signed int)*v,zero,
							   (vector unsigned int)vec_sld(zero,
															vec_splat_s32(-1),
															12));
}


/** Clear the highest element in 3 SIMD variables.
*
*  @param v1  SIMD variable to clear element 3 of.
*  @param v2  SIMD variable to clear element 3 of.
*  @param v3  SIMD variable to clear element 3 of.
*/
static inline void 
zero_highest_element_in_3_vectors(vector float *v1,
                                  vector float *v2,
                                  vector float *v3)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),4);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
}


/** Clear the highest element in 4 SIMD variables.
*
*  @param v1  SIMD variable to clear element 3 of.
*  @param v2  SIMD variable to clear element 3 of.
*  @param v3  SIMD variable to clear element 3 of.
*  @param v4  SIMD variable to clear element 3 of.
*/
static inline void 
zero_highest_element_in_4_vectors(vector float *v1,
                                  vector float *v2,
                                  vector float *v3,
								  vector float *v4)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),4);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
}



/** Clear the upper half in 3 SIMD variables.
*
*  @param v1  SIMD variable to clear elements 2-3 of.
*  @param v2  SIMD variable to clear elements 2-3 of.
*  @param v3  SIMD variable to clear elements 2-3 of.
*/
static inline void 
zero_highest_2_elements_in_3_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),8);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
}

/** Clear the upper half in 4 SIMD variables.
*
*  @param v1  SIMD variable to clear elements 2-3 of.
*  @param v2  SIMD variable to clear elements 2-3 of.
*  @param v3  SIMD variable to clear elements 2-3 of.
*  @param v4  SIMD variable to clear elements 2-3 of.
*/
static inline void 
zero_highest_2_elements_in_4_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3,
                                     vector float *v4)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),8);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
}




/** Clear the highest three element in 3 SIMD variables.
*
*  @param v1  SIMD variable to clear elements 1-3 of.
*  @param v2  SIMD variable to clear elements 1-3 of.
*  @param v3  SIMD variable to clear elements 1-3 of.
*/
static inline void 
zero_highest_3_elements_in_3_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),
															 12);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
}


/** Clear the highest three element in 4 SIMD variables.
*
*  @param v1  SIMD variable to clear elements 1-3 of.
*  @param v2  SIMD variable to clear elements 1-3 of.
*  @param v3  SIMD variable to clear elements 1-3 of.
*  @param v4  SIMD variable to clear elements 1-3 of.
*/
static inline void 
zero_highest_3_elements_in_4_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3,
                                     vector float *v4)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),
															 12);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
}


/** Clear the highest element in 9 SIMD variables.
*
*  @param v1  SIMD variable to clear last element of.
*  @param v2  SIMD variable to clear last element of.
*  @param v3  SIMD variable to clear last element of.
*  @param v4  SIMD variable to clear last element of.
*  @param v5  SIMD variable to clear last element of.
*  @param v6  SIMD variable to clear last element of.
*  @param v7  SIMD variable to clear last element of.
*  @param v8  SIMD variable to clear last element of.
*  @param v9  SIMD variable to clear last element of.
*/
static inline void
zero_highest_element_in_9_vectors(vector float *v1,
                                  vector float *v2,
                                  vector float *v3,
                                  vector float *v4,
                                  vector float *v5,
                                  vector float *v6,
                                  vector float *v7,
                                  vector float *v8,
                                  vector float *v9)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask = (vector unsigned int)vec_sld(zero,
															vec_splat_s32(-1),
															4);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
	*v5 = (vector float)vec_sel((vector signed int)*v5,zero,mask);
	*v6 = (vector float)vec_sel((vector signed int)*v6,zero,mask);
	*v7 = (vector float)vec_sel((vector signed int)*v7,zero,mask);
	*v8 = (vector float)vec_sel((vector signed int)*v8,zero,mask);
	*v9 = (vector float)vec_sel((vector signed int)*v9,zero,mask);
}


/** Clear the highest two elements in 9 SIMD variables.
 *
 *  @param v1  SIMD variable to clear upper half of.
 *  @param v2  SIMD variable to clear upper half of.
 *  @param v3  SIMD variable to clear upper half of.
 *  @param v4  SIMD variable to clear upper half of.
 *  @param v5  SIMD variable to clear upper half of.
 *  @param v6  SIMD variable to clear upper half of.
 *  @param v7  SIMD variable to clear upper half of.
 *  @param v8  SIMD variable to clear upper half of.
 *  @param v9  SIMD variable to clear upper half of.
 */
static inline void 
zero_highest_2_elements_in_9_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3,
                                     vector float *v4,
                                     vector float *v5,
                                     vector float *v6,
                                     vector float *v7,
                                     vector float *v8,
                                     vector float *v9)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),
															 8);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
	*v5 = (vector float)vec_sel((vector signed int)*v5,zero,mask);
	*v6 = (vector float)vec_sel((vector signed int)*v6,zero,mask);
	*v7 = (vector float)vec_sel((vector signed int)*v7,zero,mask);
	*v8 = (vector float)vec_sel((vector signed int)*v8,zero,mask);
	*v9 = (vector float)vec_sel((vector signed int)*v9,zero,mask);
}


/** Clear the highest three elements in 9 SIMD variables.
*
*  @param v1  SIMD variable to clear elements 1-3 of.
*  @param v2  SIMD variable to clear elements 1-3 of.
*  @param v3  SIMD variable to clear elements 1-3 of.
*  @param v4  SIMD variable to clear elements 1-3 of.
*  @param v5  SIMD variable to clear elements 1-3 of.
*  @param v6  SIMD variable to clear elements 1-3 of.
*  @param v7  SIMD variable to clear elements 1-3 of.
*  @param v8  SIMD variable to clear elements 1-3 of.
*  @param v9  SIMD variable to clear elements 1-3 of.
*/
static inline void 
zero_highest_3_elements_in_9_vectors(vector float *v1,
                                     vector float *v2,
                                     vector float *v3,
                                     vector float *v4,
                                     vector float *v5,
                                     vector float *v6,
                                     vector float *v7,
                                     vector float *v8,
                                     vector float *v9)
{
	vector signed int zero = (vector signed int) vec_zero();
	vector unsigned int mask  = (vector unsigned int)vec_sld(zero,
															 vec_splat_s32(-1),
															 12);
				     
	*v1 = (vector float)vec_sel((vector signed int)*v1,zero,mask);
	*v2 = (vector float)vec_sel((vector signed int)*v2,zero,mask);
	*v3 = (vector float)vec_sel((vector signed int)*v3,zero,mask);
	*v4 = (vector float)vec_sel((vector signed int)*v4,zero,mask);
	*v5 = (vector float)vec_sel((vector signed int)*v5,zero,mask);
	*v6 = (vector float)vec_sel((vector signed int)*v6,zero,mask);
	*v7 = (vector float)vec_sel((vector signed int)*v7,zero,mask);
	*v8 = (vector float)vec_sel((vector signed int)*v8,zero,mask);
	*v9 = (vector float)vec_sel((vector signed int)*v9,zero,mask);
}

#endif /* _ALTIVEC_UTIL_ */



