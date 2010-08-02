/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

#include <math.h>

#include "typedefs.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "calcgrid.h"

#define facNR 4
const int factor[facNR] = {2,3,5,7};

static void make_list(int start_fac,int *ng,int ng_max,int *n_list,int **list)
{
    int i,fac;
  
    if (*ng < ng_max)
    {
        if (*n_list % 100 == 0)
        {
            srenew(*list,*n_list+100);
        }
        (*list)[*n_list] = *ng;
        (*n_list)++;
        
        for(i=start_fac; i<facNR; i++)
        {
            fac = factor[i];
            /* The choice of grid size is based on benchmarks of fftw
             * and the need for a lot of factors for nice DD decomposition.
             * The base criterion is that a grid size is not included
             * when there is a larger grid size that produces a faster 3D FFT.
             * Allow any power for 2, two for 3 and 5, but only one for 7.
             * Three for 3 are ok when there is also a factor of 2.
             * Two factors of 5 are not allowed with a factor of 3 or 7.
             * A factor of 7 does not go with a factor of 5, 7 or 9.
             */
            if ((fac == 2) ||
                (fac == 3 && (*ng % 9 != 0 ||
                              (*ng % 2 == 0 && *ng % 27 != 0))) ||
                (fac == 5 && *ng % 15 != 0 && *ng % 25 != 0) ||
                (fac == 7 && *ng % 5 != 0 && *ng % 7 != 0 && *ng % 9 != 0))
            {
                *ng *= fac;
                make_list(i,ng,ng_max,n_list,list);
                *ng /= fac;
            }
        }
    }
}

static int list_comp(const void *a,const void *b)
{
  return (*((int *)a) - *((int *)b));
}

real calc_grid(FILE *fp,matrix box,real gr_sp,
               int *nx,int *ny,int *nz)
{
    int  d,n[DIM];
    int  i,j,nmin[DIM];
    rvec box_size,spacing;
    real max_spacing;
    int  ng_max,ng;
    int  n_list,*list;

    if (gr_sp <= 0)
    {
        gmx_fatal(FARGS,"invalid fourier grid spacing: %g",gr_sp);
    }

    /* New grid calculation setup:
     *
     * To maintain similar accuracy for triclinic PME grids as for rectangular
     * ones, the max grid spacing should set along the box vectors rather than
     * cartesian X/Y/Z directions. This will lead to slightly larger grids, but
     * it is much better than having to go to pme_order=6.
     *
     * Thus, instead of just extracting the diagonal elements to box_size[d], we
     * now calculate the cartesian length of the vectors.
     *
     * /Erik Lindahl, 20060402.
     */
    for(d=0; d<DIM; d++)
    {
        box_size[d] = 0;
        for(i=0;i<DIM;i++)
        {
            box_size[d] += box[d][i]*box[d][i];
        }
        box_size[d] = sqrt(box_size[d]);
    }
    
    n[XX] = *nx;
    n[YY] = *ny;
    n[ZZ] = *nz;
    
    ng = 1;
    ng_max = 1;
    for(d=0; d<DIM; d++)
    {
        nmin[d] = (int)(box_size[d]/gr_sp + 0.999);
        if (2*nmin[d] > ng_max)
        {
            ng_max = 2*nmin[d];
        }
    }
    n_list=0;
    list=NULL;
    make_list(0,&ng,ng_max,&n_list,&list);
    
    if ((*nx<=0) || (*ny<=0) || (*nz<=0))
    {
        if (NULL != fp)
        {
            fprintf(fp,"Calculating fourier grid dimensions for%s%s%s\n",
                    *nx > 0 ? "":" X",*ny > 0 ? "":" Y",*nz > 0 ? "":" Z");
        }
    }
    
    qsort(list,n_list,sizeof(list[0]),list_comp);
    if (debug)
    {
        for(i=0; i<n_list; i++)
            fprintf(debug,"grid: %d\n",list[i]);
    }
        
    for(d=0; d<DIM; d++)
    {
        for(i=0; (i<n_list) && (n[d]<=0); i++)
        {
            if (list[i] >= nmin[d])
            {
                n[d] = list[i];
            }
        }
    }
    
    sfree(list);
    
    max_spacing = 0;
    for(d=0; d<DIM; d++)
    {
        spacing[d] = box_size[d]/n[d];
        if (spacing[d] > max_spacing)
        {
            max_spacing = spacing[d];
        }
    }
    *nx = n[XX];
    *ny = n[YY];
    *nz = n[ZZ];
    if (NULL != fp)
    {
        fprintf(fp,"Using a fourier grid of %dx%dx%d, spacing %.3f %.3f %.3f\n",
                *nx,*ny,*nz,spacing[XX],spacing[YY],spacing[ZZ]);
    }

    return max_spacing;
}


