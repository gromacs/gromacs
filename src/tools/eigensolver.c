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

#ifdef GMX_MPI
#include <mpi.h>
#endif

#include "types/simple.h"
#include "smalloc.h"
#include "gmx_fatal.h"

#include "sparsematrix.h"
#include "eigensolver.h"

#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

#include "gmx_lapack.h"
#include "gmx_arpack.h"

void
eigensolver(real *   a,
            int      n,
            int      index_lower,
            int      index_upper,
            real *   eigenvalues,
            real *   eigenvectors)
{
    int *   isuppz;
    int     lwork,liwork;
    int     il,iu,m,iw0,info;
    real    w0,abstol;
    int *   iwork;
    real *  work;
    real    vl,vu;
    char *  jobz;
    
    if(index_lower<0)
        index_lower = 0;
    
    if(index_upper>=n)
        index_upper = n-1;
    
    /* Make jobz point to the character "V" if eigenvectors
     * should be calculated, otherwise "N" (only eigenvalues).
     */   
    jobz = (eigenvectors != NULL) ? "V" : "N";

    /* allocate lapack stuff */
    snew(isuppz,2*n);
    vl = vu = 0;
    
    /* First time we ask the routine how much workspace it needs */
    lwork  = -1;
    liwork = -1;
    abstol = 0;
    
    /* Convert indices to fortran standard */
    index_lower++;
    index_upper++;
    
    /* Call LAPACK routine using fortran interface. Note that we use upper storage,
     * but this corresponds to lower storage ("L") in Fortran.
     */    
#ifdef GMX_DOUBLE
    F77_FUNC(dsyevr,DSYEVR)(jobz,"I","L",&n,a,&n,&vl,&vu,&index_lower,&index_upper,
                            &abstol,&m,eigenvalues,eigenvectors,&n,
                            isuppz,&w0,&lwork,&iw0,&liwork,&info);
#else
    F77_FUNC(ssyevr,SSYEVR)(jobz,"I","L",&n,a,&n,&vl,&vu,&index_lower,&index_upper,
                            &abstol,&m,eigenvalues,eigenvectors,&n,
                            isuppz,&w0,&lwork,&iw0,&liwork,&info);
#endif

    if(info != 0)
    {
        sfree(isuppz);
        gmx_fatal(FARGS,"Internal errror in LAPACK diagonalization.");        
    }
    
    lwork = w0;
    liwork = iw0;
    
    snew(work,lwork);
    snew(iwork,liwork);
    
    abstol = 0;
    
#ifdef GMX_DOUBLE
    F77_FUNC(dsyevr,DSYEVR)(jobz,"I","L",&n,a,&n,&vl,&vu,&index_lower,&index_upper,
                            &abstol,&m,eigenvalues,eigenvectors,&n,
                            isuppz,work,&lwork,iwork,&liwork,&info);
#else
    F77_FUNC(ssyevr,SSYEVR)(jobz,"I","L",&n,a,&n,&vl,&vu,&index_lower,&index_upper,
                            &abstol,&m,eigenvalues,eigenvectors,&n,
                            isuppz,work,&lwork,iwork,&liwork,&info);
#endif
    
    sfree(isuppz);
    sfree(work);
    sfree(iwork);
    
    if(info != 0)
    {
        gmx_fatal(FARGS,"Internal errror in LAPACK diagonalization.");
    }
    
}


void 
sparse_eigensolver(gmx_sparsematrix_t *    A,
                   int                     neig,
                   real *                  eigenvalues,
                   real *                  eigenvectors,
                   int                     maxiter)
{
    int      iwork[80];
    int      iparam[11];
    int      ipntr[11];
    real *   resid;
    real *   workd;
    real *   pIN;
    real *   pOUT;
    real *   workl;
    real *   v;
    real *   buf1;
    real *   buf2;
    int      n;
    int      ido,info,lworkl,i,ncv,dovec;
    real     abstol;
    int *    select;
    int      iter;
    int      nglobal;
#ifdef GMX_MPI
    int      *nlocal;
    int     *nstart;
    int      nnodes;
    int      rank,size;
    int      scratch;
    int      comm = 0;
#endif
    
    if(eigenvectors != NULL)
        dovec = 1;
    else
        dovec = 0;
    
    nglobal  = A->nrow;

#ifdef GMX_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    scratch = nglobal % size;

    for(i=0;i<rank;i++)
    {
        nlocal[i] = nglobal / size;
        nstart[i] = i*nlocal[i];
        if(rank<scratch)
        {
            nlocal[i]++;
            nstart[i] += i;
        }
    }
    n = nlocal[rank];
    snew(buf1,nglobal);
    snew(buf2,nglobal);
#else
    n = nglobal;
#endif
    
    ncv = 2*neig;
    
    if(ncv>nglobal)
        ncv=nglobal;
    
    for(i=0;i<11;i++)
        iparam[i]=ipntr[i]=0;
	
	iparam[0] = 1;       /* Don't use explicit shifts */
	iparam[2] = maxiter; /* Max number of iterations */
	iparam[6] = 1;       /* Standard symmetric eigenproblem */
    
	lworkl = ncv*(8+ncv);
    snew(resid,n);
    snew(workd,(3*n+4));
    snew(workl,lworkl);
    snew(select,ncv);
    snew(v,n*ncv);

    /* Use machine tolerance - roughly 1e-16 in double precision */
    abstol = 0;
    
 	ido = info = 0;
    fprintf(stderr,"Calculating Ritz values and Lanczos vectors, max %d iterations...\n",maxiter);
    
    iter = 1;
	do {
#ifdef GMX_MPI
        /* MPI */
#ifdef GMX_DOUBLE
        F77_FUNC(pdsaupd,PDSAUPD)(&comm,&ido, "I", &n, "SA", &neig, &abstol, 
                                resid, &ncv, v, &n, iparam, ipntr, 
                                workd, iwork, workl, &lworkl, &info);
#else
        F77_FUNC(pssaupd,PSSAUPD)(&comm,&ido, "I", &n, "SA", &neig, &abstol, 
                                resid, &ncv, v, &n, iparam, ipntr, 
                                workd, iwork, workl, &lworkl, &info);
#endif

        if(ido==-1 || ido==1)
        {
            pIN  = workd+ipntr[0]-1;
            pOUT = workd+ipntr[1]-1;
        
            for(i=0;i<nlocal[rank];i++)
                buf1[nstart[rank]+i]=pIN[i];
        
            /* COMMUNICATE STUFF */
            for(i=0;i<size;i++)
            {
                /* overwrite stuff */
                /* Broadcast from process i to all other processes */
                MPI_Bcast(buf1+nstart[i],nlocal[i],MPI_DOUBLE,i,MPI_COMM_WORLD);
            }
        
            gmx_sparsematrix_vector_multiply_partial(A,buf1,buf2,nstart[rank],nlocal[rank]);
            
            /* COMMUNICATE BACK */
            for(i=0;i<size;i++)
            {
                /* Add values */
                /* Reduce (sum) to process i from all other processes */
                MPI_Reduce(buf2+nstart[i],pOUT,nlocal[i],MPI_DOUBLE,MPI_SUM,i,MPI_COMM_WORLD);
            }            
        }
        
#else
        /* NOT MPI */
#ifdef GMX_DOUBLE
            F77_FUNC(dsaupd,DSAUPD)(&ido, "I", &n, "SA", &neig, &abstol, 
                                    resid, &ncv, v, &n, iparam, ipntr, 
                                    workd, iwork, workl, &lworkl, &info);
#else
            F77_FUNC(ssaupd,SSAUPD)(&ido, "I", &n, "SA", &neig, &abstol, 
                                    resid, &ncv, v, &n, iparam, ipntr, 
                                    workd, iwork, workl, &lworkl, &info);
#endif
            if(ido==-1 || ido==1)
                gmx_sparsematrix_vector_multiply(A,workd+ipntr[0]-1, workd+ipntr[1]-1);
#endif

#ifdef GMX_MPI
            if(rank==0)
                fprintf(stderr,"\rIteration %4d: %3d out of %3d Ritz values converged.",iter++,iparam[4],neig);
#else
            fprintf(stderr,"\rIteration %4d: %3d out of %3d Ritz values converged.",iter++,iparam[4],neig);
#endif            
	} while(info==0 && (ido==-1 || ido==1));

	
#ifdef GMX_MPI
    if(rank==0)
    {    
#endif
    fprintf(stderr,"\n");
	
    if(info==1)
    {
	    gmx_fatal(FARGS,
                  "Maximum number of iterations (%d) reached in Arnoldi\n"
                  "diagonalization, but only %d of %d eigenvectors converged.\n",
                  maxiter,iparam[4],neig);
    }
	else if(info!=0)
    {
        gmx_fatal(FARGS,"Unspecified error from Arnoldi diagonalization:%d\n",info);
    }
	
	info = 0;
	/* Extract eigenvalues and vectors from data */
    fprintf(stderr,"Calculating eigenvalues and eigenvectors...\n");
    
    
#ifdef GMX_MPI
    }
#ifdef GMX_DOUBLE
    F77_FUNC(pdseupd,PDSEUPD)(&comm,&dovec, "A", select, eigenvalues, eigenvectors, 
                            &n, NULL, "I", &n, "SA", &neig, &abstol, 
                            resid, &ncv, v, &n, iparam, ipntr, 
                            workd, workl, &lworkl, &info);
#else
    F77_FUNC(psseupd,PSSEUPD)(&comm,&dovec, "A", select, eigenvalues, eigenvectors, 
                            &n, NULL, "I", &n, "SA", &neig, &abstol, 
                            resid, &ncv, v, &n, iparam, ipntr, 
                            workd, workl, &lworkl, &info);
#endif
    sfree(buf1);
    sfree(buf2);
#else
#ifdef GMX_DOUBLE
    F77_FUNC(dseupd,DSEUPD)(&dovec, "A", select, eigenvalues, eigenvectors, 
			    &n, NULL, "I", &n, "SA", &neig, &abstol, 
			    resid, &ncv, v, &n, iparam, ipntr, 
			    workd, workl, &lworkl, &info);
#else
    F77_FUNC(sseupd,SSEUPD)(&dovec, "A", select, eigenvalues, eigenvectors, 
			    &n, NULL, "I", &n, "SA", &neig, &abstol, 
			    resid, &ncv, v, &n, iparam, ipntr, 
			    workd, workl, &lworkl, &info);
#endif
#endif	
    sfree(v);
    sfree(resid);
    sfree(workd);
    sfree(workl);  
    sfree(select);    
}


