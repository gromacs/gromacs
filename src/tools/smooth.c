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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cdist.h"

#define NMRLEN  99.9
#define OVERLAP(a_lb,a_ub,b_lb,b_ub) !((a_ub)<(b_lb)) || ((b_ub)<(a_lb))
/* Bound smoothing for cdist*/
/* Adam Kirrander 990119 */

real _mysqrt(real x,char *fn,int line)
{
  if (x < 0) {
    if (debug)
      fprintf(debug,"MYSQRT: negative argument %g called from %s, line %d\n",
	      x,fn,line);
    return 0;
  }
  else
    return sqrt(x);
}

#define mysqrt(x) _mysqrt(x,__FILE__,__LINE__)

void check_bounds(int ai,int aj,real lb,real ub)
{
  if (lb > ub)
    gmx_fatal(FARGS,"Strange bounds for %d-%d: lb=%f, ub=%f",
		ai+1,aj+1,lb,ub);
}

void check_len_(int ai,int aj,real lb,real ub,real len,char *file,int line)
{
  if (len != NMRLEN) 
    if (((len > ub) || (len < lb)) && len > 0.0)
      gmx_fatal(FARGS,"%s, line %d: Ideal len. outside bounds for %d-%d:"
		  " lb=%f, ub=%f, len=%f",
		  file,line,ai+1,aj+1,lb,ub,len);
}

#define check_len(ai,aj,lb,ub,len) check_len_(ai,aj,lb,ub,len,__FILE__,__LINE__)

void check_tri(t_dist *d,int natoms,int i,int j,int k)
{
  real a_lb,a_ub,a_len,b_lb,b_ub,b_len,c_lb,c_ub,c_len;
  
  a_lb  = d_lb(d,natoms,i,j);
  a_ub  = d_ub(d,natoms,i,j);
  a_len = d_len(d,natoms,i,j);
  
  b_lb  = d_lb(d,natoms,i,k);
  b_ub  = d_ub(d,natoms,i,k);
  b_len = d_len(d,natoms,i,k);
  
  c_lb  = d_lb(d,natoms,j,k);
  c_ub  = d_ub(d,natoms,j,k);
  c_len = d_len(d,natoms,j,k);

  check_len(i,j,a_lb,a_ub,a_len);
  check_len(i,k,b_lb,b_ub,b_len);
  check_len(j,k,c_lb,c_ub,c_len);
  
  check_bounds(i,j,a_lb,a_ub);
  check_bounds(i,k,b_lb,b_ub);
  check_bounds(j,k,c_lb,c_ub);
  
  /* Check if triangle inequality is fulfilled */
  if ((a_lb > (b_ub+c_ub)) || 
      (b_lb > (a_ub+c_ub)) || 
      (c_lb > (a_ub+b_ub))) {
    fprintf(stderr,"Impossible triangle (%d,%d,%d):\n"\
	    "a_lb=%10.5f, a_ub=%10.5f, a_len=%10.5f\n"\
	    "b_lb=%10.5f, b_ub=%10.5f, b_len=%10.5f\n"
	    "c_lb=%10.5f, c_ub=%10.5f, c_len=%10.5f\n",\
	    i,j,k,a_lb,a_ub,a_len,b_lb,b_ub,b_len,c_lb,c_ub,c_len);
    exit(1);
  }
}

double triangle_upper_bound (t_dist *d,int natoms,real tol)
{
  double possible,ntotal=0;
  int    i,j,k,nsmooth,innerloop,count=0,nloops=0;
  real   a_lb,a_ub,a_len,b_lb,b_ub,b_len,c_lb,c_ub,c_len;

  if (debug)
    fprintf(debug,"Just entered triangle_upper_bound!\n");

  /* Loop over all triangles until no smoothing occurs */
  do {
    nsmooth=0;
    nloops=nloops+1;
    for (i=0; (i < natoms); i++) {
      for (j=i+1; (j < natoms); j++) {
	if (!dist_set(d,natoms,i,j)) continue;
	for (k=j+1; (k < natoms); k++) {

	  count=count+1; 

	  /* Should make use of distances between zero-weighted if that can 
	     help us to define better distances for e.g. loops. But make sure
	     we don't waste time smoothing zero-weighted tripples...*/

	  /* Are all distances defined? */
	  if (!dist_set(d,natoms,i,k) || !dist_set(d,natoms,j,k)) continue;
	  
	  /* smooth the triangle */
	  do {
	    /* Reset inner counter */
	    innerloop=0;

	    check_tri(d,natoms,i,j,k);
	    
	    /* Read bounds */
	    a_lb  = d_lb(d,natoms,i,j);
	    a_ub  = d_ub(d,natoms,i,j);
	    a_len = d_len(d,natoms,i,j);
	    
	    b_lb  = d_lb(d,natoms,i,k);
	    b_ub  = d_ub(d,natoms,i,k);
	    b_len = d_len(d,natoms,i,k);

	    c_lb  = d_lb(d,natoms,j,k);
	    c_ub  = d_ub(d,natoms,j,k);
	    c_len = d_len(d,natoms,j,k);

	    /* Check if triangle inequality is fulfilled */
	    if (a_ub > (b_ub+c_ub+tol)) {
	      set_dist(d,natoms,i,j,a_lb,b_ub+c_ub,a_len);
	      a_ub = d_ub(d,natoms,i,j);
	      innerloop++;
	    }
	    if (b_ub > (a_ub+c_ub+tol)) {
	      set_dist(d,natoms,i,k,b_lb,a_ub+c_ub,b_len);
	      b_ub = d_ub(d,natoms,i,k);
	      innerloop++;
	    }
	    if (c_ub > (a_ub+b_ub+tol)) {
	      set_dist(d,natoms,j,k,c_lb,a_ub+b_ub,c_len);
	      c_ub = d_ub(d,natoms,j,k); 
	      innerloop++;
	    }

	    nsmooth += innerloop;
	  } 
	  while (innerloop > 0);
	}
      }
    }
    ntotal += nsmooth;
  }
  while (nsmooth>0);
  
  possible = natoms*(natoms-1.0)*(natoms-2.0)/6;
  fprintf(stderr,"Checked %d triples (of %.0f) in %d rounds of ub"
	  " triangle smoothing.\n",
	  count/nloops,possible,nloops);
  fprintf(stderr,"Smoothed %g upper bounds with triagonal ineq.\n",ntotal);
  
  return ntotal;
}


double triangle_lower_bound (t_dist *d,int natoms,real tol)
{
  double ntotal = 0;
  int    i,j,k,nsmooth,innerloop;
  real   a_lb,a_ub,a_len,b_lb,b_ub,b_len,c_lb,c_ub,c_len,new_lb;

  /*fprintf(stderr,"Just entered triangle_lower_bound!\n");*/

  /* Loop over all triangles until no smoothing occurs */
  do {
    nsmooth=0;
    for (i=0; (i < natoms); i++) {
      for (j=i+1; (j < natoms); j++) {
	if (!dist_set(d,natoms,i,j)) continue;
	for (k=j+1; (k < natoms); k++) {

	  /* Are all distances defined? If not, continue */
	  if (!dist_set(d,natoms,i,k) || !dist_set(d,natoms,j,k)) continue;


	  /* smooth the triangle */
	  do {
	    /* Reset inner counter */
	    innerloop=0;

	    /* Read bounds */
	    a_lb  = d_lb(d,natoms,i,j);
	    a_ub  = d_ub(d,natoms,i,j);
	    a_len = d_len(d,natoms,i,j);
	    
	    b_lb = d_lb(d,natoms,i,k);
	    b_ub = d_ub(d,natoms,i,k);
	    b_len = d_len(d,natoms,i,k);

	    c_lb  = d_lb(d,natoms,j,k);
	    c_ub  = d_ub(d,natoms,j,k);
	    c_len = d_len(d,natoms,j,k);

	    /* Smooth lower bound */
	    if (!OVERLAP(b_lb,b_ub,c_lb,c_ub)) {
	      new_lb = min(fabs(b_lb-c_ub),fabs(c_lb-b_ub));
	      if ( tol < new_lb - a_lb) {
		set_dist(d,natoms,i,j,new_lb,a_ub,a_len);
		nsmooth++;
		innerloop++;
		/*fprintf(stderr,"Smoothed %d-%d",i,j);*/
	      }
	    }
	    if (!OVERLAP(a_lb,a_ub,c_lb,c_ub)) {
	      new_lb = min(fabs(a_lb-c_ub),fabs(c_lb-a_ub));
	      if ( tol < new_lb - b_lb) {
		set_dist(d,natoms,i,k,new_lb,b_ub,b_len);
		nsmooth++;
		innerloop++;
		/*fprintf(stderr,"Smoothed %d-%d",i,k);*/
	      }
	    }
	    if (!OVERLAP(a_lb,a_ub,b_lb,b_ub)) {
	      new_lb = min(fabs(a_lb-b_ub),fabs(b_lb-a_ub));
	      /*Programming error? */
	      /* if ( (a_ub+b_ub) < c_ub ) { */
	      /* New if-statement 990609 Adam */
	      if ( tol < new_lb - c_lb) {
		set_dist(d,natoms,j,k,new_lb,c_ub,c_len);
		nsmooth++;
		innerloop++;
		/*fprintf(stderr,"Smoothed %d-%d",j,k);*/
	      }
	    }
	      
	    /* Check that bounds make sense */
	    check_bounds(i,j,a_lb,a_ub);
	    check_bounds(i,k,b_lb,b_ub);
	    check_bounds(j,k,c_lb,c_ub);
	  } 
	  while (innerloop > 0);
	}
      }
    }
    ntotal += nsmooth;
  }
  while (nsmooth>0);
  fprintf(stderr,"Smoothed %g lower bounds with triagonal ineq.\n",ntotal);
  return ntotal;
}

double do_triangle (t_dist *d,t_atoms *atoms,real tol)
{
  double ntot=0;
  /* Do triangular smoothing only */
  int natoms = atoms->nr;
  
  fprintf(stderr,"Starting triangle smoothing\n");
  ntot += triangle_upper_bound (d,natoms,tol);
  ntot += triangle_lower_bound (d,natoms,tol);
  
  return ntot;
}

/* Tetrangle-bound smoothing for cdist.
   Adam Kirrander 990209 */

/* According Havel,Kuntz and Crippen 1983 Bull.Math.Biol */




/* FUNCTIONS ASSOCIATED WITH TRIANGLE_LIMITED */
/* Routines to check if triangle inequalty limits are attainable between
   points 1 and 4, given that the 5 other distances are defined.
   Sets logical TUL=true if the upper limit of distance 1-4 is defined by
   triangle inequalty, and TLL=true if corresponding lower limit is set by 
   the triangle ineq. Written to be used in conjunction with tetragonal
   bound smoothing in cdist.

   The main routine is trianglelimited, and it is the one that should be 
   called. The rest are supporting routines.
   Ref.: Bull.Math.Biol.Vol45,Nr5,665-720,Havel,Crippen,Kuntz,1983 */

/* Adam Kirrander    990209  */

/*#define OVERLAP(a_lb,a_ub,b_lb,b_ub) !((a_ub<b_lb) || (b_ub<a_lb))*/
/* ((a_ub >= b_lb) && (b_ub >= a_lb)) */
#define apex(AC,AD,BC,BD) mysqrt((BC*(AD*AD-BD*BD)+BD*(AC*AC-BC*BC))/(BC+BD))


bool _triangle_bound_attainable(real BC,real BD,real AClb,real ACub,
				real ADlb,real ADub,real ABlb,real ABub,
				char *fn,int line)
{
  /* Assume triangle with A,C and D at corners and B at C-D. 
     Given the distances AC,AD,BC,BD it calculates the distance AB 
     at its min and max, which corresponds to AC and AD being min and 
     max respectively. 
     If this intervall overlaps with the set bounds for AB, it follows that 
     the points C,B and D can be made colinear and thus are restricted by the
     triangle inequality. */

  real ABmin,ABmax;

  if (debug) 
    fprintf(debug,"%s, line %d: ABlb: %g, ABub: %g, AClb: %g, ACub: %g, ADlb: %g, ADub: %g, BC: %g, BD: %g\n",
	    fn,line,ABlb,ABub,AClb,ACub,ADlb,ADub,BC,BD);
  ABmax = apex(ACub,ADub,BC,BD);
    /* Triangle can't be formed, i.e. there is a risk all 4 points can
     become colinear */
  if ( (AClb+ADlb) < (BC+BD) )
    return ( ABlb <= ABmax );
  
  /* I am not sure the above truly is the optimal solution, but it should
   * be the safe solution 
   */

  ABmin = apex(AClb,ADlb,BC,BD);

  return (OVERLAP(ABlb,ABub,ABmin,ABmax));
}

#define triangle_bound_attainable(BC,BD,AClb,ACub,ADlb,ADub,ABlb,ABub) \
       _triangle_bound_attainable(BC,BD,AClb,ACub,ADlb,ADub,ABlb,ABub, \
                                  __FILE__,__LINE__)


void triangle_limited(int *ATMS,int natoms,t_dist *d,bool *TUL, bool *TLL)
{
  /* Given 4 atoms it checks if the triangle inequality lower or upper bounds 
     for the distance 1-4 are attainable. */
  /* Situation looks something like the following: */
  /* Can:                                          */
  /*   1                      4                    */
  /*    \\\               ////                     */
  /*        2------------3                         */
  /* become:                                       */
  /*   1---------------3------4                    */
  /*    \\\\\\\ ///////                            */
  /*           2                                   */


  int N1=ATMS[0],N2=ATMS[1],N3=ATMS[2],N4=ATMS[3];
  real d12[2],d13[2],d14[2],d23[2],d24[2],d34[2];


  /* Initiate and read bounds */
  *TUL = FALSE;
  *TLL = FALSE;
  
  d12[0] = d_lb(d,natoms,N1,N2);
  d12[1] = d_ub(d,natoms,N1,N2);
  
  d13[0] = d_lb(d,natoms,N1,N3);
  d13[1] = d_ub(d,natoms,N1,N3);
  
  d14[0] = d_lb(d,natoms,N1,N4);
  d14[1] = d_ub(d,natoms,N1,N4);
  
  d23[0] = d_lb(d,natoms,N2,N3);
  d23[1] = d_ub(d,natoms,N2,N3);

  d24[0] = d_lb(d,natoms,N2,N4);
  d24[1] = d_ub(d,natoms,N2,N4);
  
  d34[0] = d_lb(d,natoms,N3,N4);
  d34[1] = d_ub(d,natoms,N3,N4);
  
  /* Check if UPPER triangle inequality limit is attainable. */
  if ( d12[1] + d24[1] < d13[1] + d34[1] ) {
    /* Corresponds to N1,N2 and N4 colinear.
       N1=D,N2=B,N3=A,N4=C ;in notation of subroutines */
    *TUL = triangle_bound_attainable(d24[1],d12[1],d34[0],d34[1],
				     d13[0],d13[1],d23[0],d23[1]);
  }
  else if ( d12[1] + d24[1] > d13[1] + d34[1] ) {
    /* Corresponds to N1,N3 and N4 colinear.
       N1=D,N2=A,N3=B,N4=C ;in notation of subroutines */
    *TUL = triangle_bound_attainable(d34[1],d13[1],d24[0],d24[1],
				     d12[0],d12[1],d23[0],d23[1]);
  }
  /* if N2 and N3 can superimpose TUL is true by necessity (AB=0) */
  else if (d23[0] == 0) *TUL = TRUE;
  
  /* Check if LOWER triangle inequality limit is attainable */ 
  if ( (d13[0] - d34[1] <= 0) && (d34[0] - d13[1] <= 0) &&
       (d24[0] - d12[1] <= 0) && (d12[0] - d24[1] <= 0) ) {
    /* if all inverse triangle ineq. limits are zero then */
    *TLL = TRUE;
  }
  else if (d12[0] - d24[1] > 0) {
    /* Corresponds to N1,N4 and N2 colinear.
       A=N3,B=N4,C=N1,D=N2 */
    *TLL = triangle_bound_attainable(d24[1],d12[0]-d24[1],d23[0],d23[1],
				     d13[0],d13[1],d34[0],d34[1]);
  }
  else if (d24[0] - d12[1] > 0) {
    /* Corresponds to N2,N1 and N4 colinear.
       A=N3,B=N1,C=N2,D=N4 */
    *TLL = triangle_bound_attainable(d12[1],d24[0]-d12[1],d23[0],d23[1],
				     d34[0],d34[1],d13[0],d13[1]);
  }
  else if (d13[0] - d34[1] > 0) {
    /* Corresponds to N1,N4 and N3 colinear.
       A=N2,B=N4,C=N3,D=N1 */
    *TLL = triangle_bound_attainable(d34[1],d13[0]-d34[1],d23[0],d23[1],
				     d12[0],d12[1],d24[0],d24[1]);
  }
  else if (d34[0] - d13[1] > 0) {
    /* Corresponds to N3,N1 and N4 colinear.
       A=N2,B=N1,C=N3,D=N4 */
    *TLL = triangle_bound_attainable(d13[1],d34[0]-d13[1],d23[0],d23[1],
				     d24[0],d24[1],d12[0],d12[1]);
  }

}
/* END OF FUNCTIONS >TRIANGLE_LIMITED< */



/* FUNCTIONS ASSOCIATED WITH TETRAGONAL-SMOOTHING */ 

void minmax(real array[],int SIZE,int *min,int *max)
{
  /* Finds min and max (indices thereof) in array.
     Early version that only accepts real. 
     Fix it!
     Adam Kirrander 990211                         */

  int i;  
  *min = 0;
  *max = 0;

  for (i=1; (i<SIZE); i++) {
    if (array[i] < array[*min]) 
      *min = i;
    else if (array[i] > array[*max]) 
      *max = i;
  }
}

void swap (int *x,int *y)
{
  /* Swaps x and y.
     Early version that only can swap integers.
     Fix it!
     Adam Kirrander 990211                         */

  int temp;
  temp = *x;
  *x = *y;
  *y = temp;
}


bool alldistset (int i,int j,int k,int l,int natoms,t_dist *d)
{
  /* Returns FALSE if all distances are not set */
  if (!dist_set(d,natoms,i,j)) return FALSE;
  if (!dist_set(d,natoms,i,k)) return FALSE;
  if (!dist_set(d,natoms,i,l)) return FALSE;
  if (!dist_set(d,natoms,j,k)) return FALSE;  
  if (!dist_set(d,natoms,j,l)) return FALSE;
  if (!dist_set(d,natoms,k,l)) return FALSE;

  return TRUE;
}

int nr_nb (int i,int j,int k,int l,int natoms,t_dist *d)
{
  /* Counts the number of nb's */
  int nnb=0;

  if (d_len(d,natoms,i,j) == 0.0) nnb++;
  if (d_len(d,natoms,i,k) == 0.0) nnb++;
  if (d_len(d,natoms,i,l) == 0.0) nnb++;
  if (d_len(d,natoms,j,k) == 0.0) nnb++;
  if (d_len(d,natoms,j,l) == 0.0) nnb++;
  if (d_len(d,natoms,k,l) == 0.0) nnb++;

  return nnb;
}

void sort_atoms(int natoms,t_dist *d,int *ATMS)
{
  /* Sorts atoms. The nb-atoms end up in ATMS[0] and ATMS[3]. */
  
  if (d_len(d,natoms,ATMS[0],ATMS[1]) == 0.0)
    swap(&ATMS[1],&ATMS[3]);
  else if (d_len(d,natoms,ATMS[0],ATMS[2]) == 0.0) 
    swap(&ATMS[2],&ATMS[3]);
  else if (d_len(d,natoms,ATMS[1],ATMS[2]) == 0.0) {
    swap(&ATMS[0],&ATMS[1]);
    swap(&ATMS[2],&ATMS[3]);
  }
  else if (d_len(d,natoms,ATMS[1],ATMS[3]) == 0.0) 
    swap(&ATMS[1],&ATMS[0]);
  else if (d_len(d,natoms,ATMS[2],ATMS[3]) == 0.0) 
    swap(&ATMS[2],&ATMS[0]);
  
  /* Put the two middle ones in order (ambivalent procedure, necessary?) */
  if (d_len(d,natoms,ATMS[0],ATMS[1]) > d_len(d,natoms,ATMS[0],ATMS[2])) 
    swap(&ATMS[1],&ATMS[2]);
}

real solve_tetrangle(real dAB,real dAC,real dBC,real dBD,real dCD,int cosphi)
{
  /* Solves tetrangle eq. for distance AD.
     cosphi=cos(phi), where phi is the dihedral angle
     cosphi=1 corresponds to min, cosphi=-1 to max
     eq. solved and optimized by Maple, watch out for bugs */

  real t1,t2,t3,t4,t5,t7,t8,t9,t11,t12,t13,t14,t15,t17,t20;
  
  t1 = dAB*dAB;
  t2 = dBC*dBC;
  t3 = t1*t2;
  t4 = t1*t1;
  t5 = dAC*dAC;
  t7 = t2*t2;
  t8 = t5*t2;
  t9 = t5*t5;
  t11 = dCD*dCD;
  t12 = t2*t11;
  t13 = dBD*dBD;
  t14 = t2*t13;
  t15 = t11*t11;
  t17 = t13*t13;
  t20 = mysqrt((-2.0*t3+t4-2.0*t1*t5+t7-2.0*t8+t9)*(-2.0*t12+t7-2.0*t14+t15
						    -2.0*t11*t13+t17));
  return mysqrt(-(cosphi*t20-t8-t14+t7-t3-t1*t11+t1*t13-t12+t5*t11-t5*t13)/(2*t2));
  /* t20 = mysqrt((-2.0*(t3+t1*t5+t8)     + t4 + t7  + t9)*
     (-2.0*(t12+t14+t11*t13) + t7 + t15 + t17));
     return mysqrt(-0.5*(cosphi*t20-t8-t14+t7-t3-t1*t11+t1*t13-t12+t5*t11-t5*t13))/dBC;*/
}


int tetrangle_bound(int *ATMS,int cosphi,int natoms,t_dist *d,real tol)
{
  int sol=0,D1,D2,D3,D4,D5,nsmooth=0,dmin,dmax;
  real dAD[32],dAB,dAC,dBC,dBD,dCD,present_lb,present_ub,present_len;
  
  /*fprintf(stderr,"entering tetrangle_bound\n");*/

  /* Try all 32 combinations of bounds */
  for (D1=0; (D1<2); D1++) {
    for (D2=0; (D2<2); D2++) {
      for (D3=0; (D3<2); D3++) {
	for (D4=0; (D4<2); D4++) {
	  for (D5=0; (D5<2); D5++) {
	    dAB = (D1) ? d_lb(d,natoms,ATMS[0],ATMS[1]):
	      d_ub(d,natoms,ATMS[0],ATMS[1]);
	    dAC = (D2) ? d_lb(d,natoms,ATMS[0],ATMS[2]):
	      d_ub(d,natoms,ATMS[0],ATMS[2]);
	    dBC = (D3) ? d_lb(d,natoms,ATMS[1],ATMS[2]):
	      d_ub(d,natoms,ATMS[1],ATMS[2]);
	    dBD = (D4) ? d_lb(d,natoms,ATMS[1],ATMS[3]):
	      d_ub(d,natoms,ATMS[1],ATMS[3]);
	    dCD = (D5) ? d_lb(d,natoms,ATMS[2],ATMS[3]):
	      d_ub(d,natoms,ATMS[2],ATMS[3]);

	    dAD[sol] = solve_tetrangle(dAB,dAC,dBC,dBD,dCD,cosphi);
	    if (debug)
	      fprintf(debug,"dAD[%d]=%10.5f\n",sol,dAD[sol]);
	    sol++;
	  }
	}
      }
    }
  }

  /* we need these to re-set one of the bounds for the distance */
  present_len = d_len(d,natoms,ATMS[0],ATMS[3]);
  present_lb  = d_lb(d,natoms,ATMS[0],ATMS[3]);
  present_ub  = d_ub(d,natoms,ATMS[0],ATMS[3]);

  /* Set the new bound(s) */
  minmax(dAD,32,&dmin,&dmax);
  if (debug)
    fprintf(debug,"Max=dAD[%d] (%12g), Min=dAD[%d] (%12g)\n",
	    dmax,dAD[dmax],dmin,dAD[dmin]);
  /* Set new upper limit  */
  if ((cosphi == -1)  &&  ((present_ub - dAD[dmax]) > tol) ) {
    set_dist(d,natoms,ATMS[0],ATMS[3],present_lb,dAD[dmax],present_len);
    if (debug)
      fprintf(debug,"Corrected %d-%d-%d-%d ub to %10.5f\n",ATMS[0]+1,
	      ATMS[1]+1,ATMS[2]+1,ATMS[3]+1,dAD[dmax]);
    nsmooth++;
  }
  /* Set new lower limit  */
  else if ((cosphi == 1)  &&  ( tol < (dAD[dmin] - present_lb))) { 
    set_dist(d,natoms,ATMS[0],ATMS[3],dAD[dmin],present_ub,present_len);
    if (debug) 
      fprintf(debug,"Corrected %d-%d-%d-%d lb to %10.5f\n",ATMS[0]+1,
	      ATMS[1]+1,ATMS[2]+1,ATMS[3]+1,dAD[dmin]);
    nsmooth++;
  }

  /* Check the new bounds */
  if (d_ub(d,natoms,ATMS[0],ATMS[3]) < d_lb(d,natoms,ATMS[0],ATMS[3]))
    gmx_fatal(FARGS,"Fatal error for tetr. smooth of (%d,%d) ub=%f, lb=%f",
		ATMS[0],ATMS[3],d_ub(d,natoms,ATMS[0],ATMS[3]),
		d_lb(d,natoms,ATMS[0],ATMS[3]));
  
  return nsmooth;
}
  


int tetrangle (t_dist *d,int natoms,real tol)
{

  int i,j,k,l,AT[4],nsmooth=0,nnonbonds1=0,nnonbonds2=0;
  bool TLL=FALSE,TUL=FALSE;

  fprintf(stderr,"Starting tetrangle smoothing\n");
  /* Loops over all sets of four atoms */
  for (i=0;i<(natoms-3);i++) {
    for (j=i+1;j<(natoms-2);j++) {
      if (!dist_set(d,natoms,i,j)) continue;
      for (k=j+1;k<(natoms-1);k++) {
	if (!dist_set(d,natoms,i,k)) continue;
	if (!dist_set(d,natoms,j,k)) continue;
	nnonbonds1 = 0;
	if (d_len(d,natoms,i,j) == 0.0) nnonbonds1++;
	if (d_len(d,natoms,i,k) == 0.0) nnonbonds1++;
	if (d_len(d,natoms,j,k) == 0.0) nnonbonds1++;	
	if (nnonbonds1 > 1) continue;
	for (l=k+1;l<natoms;l++) {
	  if (!dist_set(d,natoms,i,l)) continue;
	  if (!dist_set(d,natoms,j,l)) continue;
	  if (!dist_set(d,natoms,k,l)) continue;
	  nnonbonds2 = 0;
	  if (d_len(d,natoms,i,l) == 0.0) nnonbonds2++;
	  if (d_len(d,natoms,j,l) == 0.0) nnonbonds2++;
	  if (d_len(d,natoms,k,l) == 0.0) nnonbonds2++;
	  if ( (nnonbonds1+nnonbonds2) != 1) continue;
	 
	  /*fprintf(stderr,"Trying %d-%d-%d-%d\n",i+1,j+1,k+1,l+1);*/

	  /* The two following subroutines have been nested into the loops
	     above. Doesn't look as good, but should be substantially 
	     faster. */
	  /*	  if (!alldistset(i,j,k,l,natoms,d)) continue; */
	  /*       if (nr_nb(i,j,k,l,natoms,d) != 1)  continue; */

	  /* Copy indices to array AT for easier handling & sort them */
	  AT[0] = i;
	  AT[1] = j;
	  AT[2] = k;
	  AT[3] = l;
	  /*fprintf(stderr,"OK'd atoms    (%d,%d,%d,%d)\n",i+1,j+1,k+1,l+1);*/
	  sort_atoms(natoms,d,AT);
	  /*fprintf(stderr,"After sorting (%d,%d,%d,%d)\n",AT[0]+1,AT[1]+1,AT[2]+1,AT[3]+1);*/

	  /* Are the distances limited by the triangle-ineq.? */
	  triangle_limited(AT,natoms,d,&TUL,&TLL);
	  /*if (TUL) fprintf(stderr,"TUL=true\n");
	    else fprintf(stderr,"TUL=false\n");
	    if (TLL) fprintf(stderr,"TLL=true\n");
	    else fprintf(stderr,"TLL=false\n");	  */

	  /* If not, correct bounds according to tetrangle-ineq. */
	  /* upper limit */
	  if (!TUL) nsmooth += tetrangle_bound(AT,-1,natoms,d,tol);
	  /* lower limit */
	  if (!TLL) nsmooth += tetrangle_bound(AT,1,natoms,d,tol);
	  
	  /* Comment: If we were to calculate upper and lower bound 
	     simultaneously we could make the code a bit faster     */
	}
      }
    }
  }
  return nsmooth;
}

/* END OF FUNCTIONS >TETRANGLE-SMOOTH< */



void do_smooth (t_dist *d,t_atoms *atoms,real tol)
{
  /* Triangular and tetragonal bound smoothing.
     Do we need to nestle triangular and tetragonal 
     smoothing or can we do them separately? THINK! */

  /* In this case I think we can get away with this, because
     we are not using distances that we correct with tetragonal
     ineq. to smooth other distances. The triangle smoothing
     should then spread the tighter bounds to other distances,
     not corrected by the teragonal ineq.
     (Of course, if you correct one distance by tetragonal ineq.
     it will become shorter, and other distances that are involved
     in triagonal ineq. with it will be affected.)
  */

  /* It should only do 2 loops, and in the second no tetragonal
     correction will be made. This is because (in the current version)
     we only correct nb's that are connected by bonded distances, and 
     these are unlikely to change during triagonal smoothing. */

  /* If this turns out to be true => SKIP DO-LOOP AND JUST HAVE
     TRIAGONAL-TETRAGONAL-TRIAGONAL */
  
  int ntri;
  int nloops = 0,nsmooth,natoms = atoms->nr;

  if (debug)
    fprintf(debug,"Starting tetragonal smoothing routine\n");
    
  /* ntri = do_triangle(d,atoms,tol);
     do { 
     nloops += 1; 
     
     nsmooth = tetrangle (d,natoms,tol);
     fprintf(stderr,"Smoothed %d distances with tetr. ineq.\n",nsmooth);
     if (nsmooth > 0)
     ntri = do_triangle(d,atoms,tol);
     else
     ntri = 0;
     
     fprintf(stderr,"Finished %d loop(s) of smoothing.\n",nloops);
     } while (ntri > 0); 
  */
  ntri    = do_triangle(d,atoms,tol);
  nsmooth = tetrangle (d,natoms,tol);
  fprintf(stderr,"Smoothed %d distances with tetr. ineq.\n",nsmooth);
  ntri    = do_triangle(d,atoms,tol);
  
  fprintf(stderr,"Finished smoothing.\n");
}
