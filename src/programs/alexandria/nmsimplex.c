/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * Program: nmsimplex.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * 11/3/97
 * $Id: nmsimplex.c,v 1.2 2007/07/23 13:46:38 mike Exp $
 *
 * gcc -Wall -o nmsimplex nmsimplex.c -lm
 * 
 * An implementation of the Nelder-Mead simplex method applied to
 * Rosenbrock's function.
 *
 * Copyright (c) 1997-2007 <Michael F. Hutt>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 * Jan. 6, 1999 
 * Modified to conform to the algorithm presented
 * in Margaret H. Wright's paper on Direct Search Methods.
 *
 * Jul. 23, 2007
 * Fixed memory leak.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nmsimplex.h"

#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */

static void print_progress(FILE *fp,int iter,int n,double **v,double *f)
{
    int i,j;
    
	if (NULL != fp)
    {
        if (iter == 0)
            fprintf(fp,"Initial values\n");
        else
            fprintf(fp,"Iteration %d\n",iter);
        
        for (j=0;j<=n;j++) 
        {
            fprintf(fp,"%2d",j+1);
            for (i=0;i<n;i++) 
            {
                fprintf(fp,"  %10g",v[j][i]);
            }
            fprintf(fp,"  chi2: %10g\n",f[j]);
	    }
    }
}

int nmsimplex(FILE *fp,
              void *data,
              nm_target_func func,
              double start[],
              int n,
              double toler,
              double scale,
              int    MAX_IT,
              double *chi2)
{
	
	int vs;         /* vertex with smallest value */
	int vh;         /* vertex with next smallest value */
	int vg;         /* vertex with largest value */
	int bConv;      /* Are we converged yet? */
	
	int i,j,m,row;
	int k;   	/* track the number of function evaluations */
	int itr;	/* track the number of iterations */
	
	double **v;     /* holds vertices of simplex */
	double pn,qn;   /* values used to create initial simplex */
	double *f;      /* value of function at each vertex */
	double fr;      /* value of function at reflection point */
	double fe;      /* value of function at expansion point */
	double fc;      /* value of function at contraction point */
	double *vr;     /* reflection - coordinates */
	double *ve;     /* expansion - coordinates */
	double *vc;     /* contraction - coordinates */
	double *vm;     /* centroid - coordinates */
	double min;
	
	double fsum,favg,s,cent;
	
	/* dynamically allocate arrays */
	bConv = 0;
		
	/* allocate the rows of the arrays */
	v =  (double **) malloc ((n+1) * sizeof(double *));
	f =  (double *) malloc ((n+1) * sizeof(double));
	vr = (double *) malloc (n * sizeof(double));
	ve = (double *) malloc (n * sizeof(double));  
	vc = (double *) malloc (n * sizeof(double));  
	vm = (double *) malloc (n * sizeof(double));  
	
	/* allocate the columns of the arrays */
	for (i=0;i<=n;i++) {
		v[i] = (double *) malloc (n * sizeof(double));
	}
	
	
	/* create the initial simplex */
	/* assume one of the vertices is 0,0 */
	
	pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));
	qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));
	
	for (i=0;i<n;i++) {
		v[0][i] = start[i];
	}
	
	for (i=1;i<=n;i++) {
		for (j=0;j<n;j++) {
			if (i-1 == j) {
				v[i][j] = pn + start[j];
			}
			else {
				v[i][j] = qn + start[j];
			}
		}
	}
	
	/* find the initial function values */
	for (j=0;j<=n;j++) {
        f[j] = func(data,v[j]);
	}
	
	k = n+1;
	
	/* print out the initial values */
    print_progress(fp,0,n,v,f);
	
	/* begin the main loop of the minimization */
	for (itr=1; (itr<=MAX_IT) && !bConv; itr++) {     
		/* find the index of the largest value */
		vg=0;
		for (j=0;j<=n;j++) {
			if (f[j] > f[vg]) {
				vg = j;
			}
		}
		
		/* find the index of the smallest value */
		vs=0;
		for (j=0;j<=n;j++) {
			if (f[j] < f[vs]) {
				vs = j;
			}
		}
		
		/* find the index of the second largest value */
		vh=vs;
		for (j=0;j<=n;j++) {
			if (f[j] > f[vh] && f[j] < f[vg]) {
				vh = j;
			}
		}
		
		/* calculate the centroid */
		for (j=0;j<=n-1;j++) {
			cent=0.0;
			for (m=0;m<=n;m++) {
				if (m!=vg) {
					cent += v[m][j];
				}
			}
			vm[j] = cent/n;
		}
		
		/* reflect vg to new vertex vr */
		for (j=0;j<=n-1;j++) {
			/*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
			vr[j] = vm[j]+ALPHA*(vm[j]-v[vg][j]);
		}
		fr = func(data,vr);
		k++;
		
		if (fr < f[vh] && fr >= f[vs]) {
			for (j=0;j<=n-1;j++) {
				v[vg][j] = vr[j];
			}
			f[vg] = fr;
		}
		
		/* investigate a step further in this direction */
		if ( fr <  f[vs]) {
			for (j=0;j<=n-1;j++) {
				/*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
				ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
			}
			fe = func(data,ve);
			k++;
			
			/* by making fe < fr as opposed to fe < f[vs], 			   
			   Rosenbrocks function takes 63 iterations as opposed 
			   to 64 when using double variables. */
			
			if (fe < fr) {
				for (j=0;j<=n-1;j++) {
					v[vg][j] = ve[j];
				}
				f[vg] = fe;
			}
			else {
				for (j=0;j<=n-1;j++) {
					v[vg][j] = vr[j];
				}
				f[vg] = fr;
			}
		}
		
		/* check to see if a contraction is necessary */
		if (fr >= f[vh]) {
			if (fr < f[vg] && fr >= f[vh]) {
				/* perform outside contraction */
				for (j=0;j<=n-1;j++) {
					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
					vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
				}
				fc = func(data,vc);
				k++;
			}
			else {
				/* perform inside contraction */
				for (j=0;j<=n-1;j++) {
					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
					vc[j] = vm[j]-BETA*(vm[j]-v[vg][j]);
				}
				fc = func(data,vc);
				k++;
			}
			
			
			if (fc < f[vg]) {
				for (j=0;j<=n-1;j++) {
					v[vg][j] = vc[j];
				}
				f[vg] = fc;
			}
			/* at this point the contraction is not successful,
			   we must halve the distance from vs to all the 
			   vertices of the simplex and then continue.
			   10/31/97 - modified to account for ALL vertices. 
			*/
			else {
				for (row=0;row<=n;row++) {
					if (row != vs) {
						for (j=0;j<=n-1;j++) {
							v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
						}
					}
				}
				f[vg] = func(data,v[vg]);
				k++;
				f[vh] = func(data,v[vh]);
				k++;
				
				
			}
		}
		
		/* print out the value at each iteration */
        print_progress(fp,itr,n,v,f);
		
		/* test for convergence */
		fsum = 0.0;
		for (j=0;j<=n;j++) {
			fsum += f[j];
		}
		favg = fsum/(n+1);
		s = 0.0;
		for (j=0;j<=n;j++) {
			s += pow((f[j]-favg),2.0)/(n);
		}
		s = sqrt(s);
		bConv = (s < toler);
	}
	/* end main loop of the minimization */
	
	/* find the index of the smallest value */
    if (bConv)
    {
        vs=0;
        for (j=0;j<=n;j++) {
            if (f[j] < f[vs]) {
                vs = j;
            }
        }
        
        for (j=0;j<n;j++) 
        {
            start[j] = v[vs][j];
        }
        min=func(data,v[vs]);
        k++;
        
        if (NULL != fp)
        {
            fprintf(fp,"\nA local minimum was found at chi2 = %g\n",min);
            fprintf(fp,"after %d fuction calls and %d iterations.\n",k,itr); 
            fprintf(fp,"Parameters:");
            for (j=0;j<n;j++) {
                fprintf(fp," %10g",v[vs][j]);
            }
            fprintf(fp,"\n\n");
        }
    }
    else if (NULL != fp)
    {
        fprintf(fp,"\nNM simplex did not converge. Size is %g, tolerance %g.\n",
                s,toler); 
    }
	free(f);
	free(vr);
	free(ve);
	free(vc);
	free(vm);
	for (i=0;i<=n;i++) {
		free (v[i]);
	}
	free(v);
	
	*chi2 = min;
	
	return bConv;
}

#ifdef TEST_NMSIMPLEX
double rosen(void *data,double x[])
{
	return (100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(1.0-x[0])*(1.0-x[0]));
}


int main()
{
	double start[] = {-1.2,1.0};
	double min;
	int i,bConv;
	
	bConv = nmsimplex(NULL,stdout,rosen,start,2,1.0e-4,1,1000,&min);
	
	for (i=0;i<2;i++) 
    {
		printf("%f\n",start[i]);
	}
	return 0;
}

#endif
