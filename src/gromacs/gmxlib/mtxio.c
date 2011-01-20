/*
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

/* This module provides routines to read/write sparse or full storage
 * matrices from/to files. It is normally used for the Hessian matrix
 * in normal mode analysis.
 */

#include "xdrf.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "mtxio.h"
#include "gmxfio.h"


/* Just a number to identify our file type */
#define GMX_MTXIO_MAGIC_NUMBER  0x34ce8fd2

#define GMX_MTXIO_FULL_MATRIX     0
#define GMX_MTXIO_SPARSE_MATRIX   1



/* Matrix file format definition:
 *
 * All entries are stored in XDR format.
 *
 * 1. Magic number integer, should be GMX_MTXIO_MAGIC_NUMBER
 * 2. An XDR string specifying the Gromacs version used to generate the file.
 * 3. Integer to denote precision. 1 if double, 0 if single precision.
 * 4. Two integers specifying number of rows and columns.
 * 5. Integer to denote storage type:
*    GMX_MTXIO_FULL_MATRIX or GMX_MTXIO_SPARSE_MATRIX
 *
 * 6. Matrix data.
 *    a) In case of full matrix, this is nrow*ncol floating-point values.
 *    b) In case of sparse matrix the data is:
 *       - Integer specifying compressed_symmetric format (1=yes, 0=no)
 *       - Integer specifying number of rows (again)
 *       - nrow integers specifying the number of data entries on each row ("ndata")
 *       - All the actual entries for each row, stored contiguous.
 *         Each entry consists of an integer column index and floating-point data value.
 */

void gmx_mtxio_write(const char *             filename,
                     int                      nrow,
                     int                      ncol,
                     real *                   full_matrix,
                     gmx_sparsematrix_t *     sparse_matrix)
{
    t_fileio *fio;
    XDR *   xd;
    int     i,j,prec;
    gmx_bool    bDum = TRUE;
    gmx_bool    bRead = FALSE;
    size_t  sz;
    
    if(full_matrix!=NULL && sparse_matrix!=NULL)
    {
        gmx_fatal(FARGS,"Both full AND sparse matrix specified to gmx_mtxio_write().\n");
    }
    
    fio = gmx_fio_open(filename,"w");
    gmx_fio_checktype(fio);
    xd = gmx_fio_getxdr(fio);
    
    /* Write magic number */
    i = GMX_MTXIO_MAGIC_NUMBER;
    gmx_fio_do_int(fio, i);
    
    /* Write generating Gromacs version */
    gmx_fio_write_string(fio, GromacsVersion());
    
    /* Write 1 for double, 0 for single precision */
    if(sizeof(real)==sizeof(double))
        prec = 1;
    else
        prec = 0;
    gmx_fio_do_int(fio, prec);
    
    gmx_fio_do_int(fio, nrow);
    gmx_fio_do_int(fio, ncol);
    
    if(full_matrix!=NULL)
    {
        /* Full matrix storage format */
        i = GMX_MTXIO_FULL_MATRIX;
        gmx_fio_do_int(fio, i);
        sz = nrow*ncol;        
        bDum=gmx_fio_ndo_real(fio, full_matrix,sz);
    }
    else
    {
        /* Sparse storage */
        i = GMX_MTXIO_SPARSE_MATRIX;
        gmx_fio_do_int(fio, i);

        gmx_fio_do_gmx_bool(fio, sparse_matrix->compressed_symmetric);
        gmx_fio_do_int(fio, sparse_matrix->nrow);
        if(sparse_matrix->nrow != nrow)
        {
            gmx_fatal(FARGS,"Internal inconsistency in sparse matrix.\n");
        }
        bDum=gmx_fio_ndo_int(fio, sparse_matrix->ndata,sparse_matrix->nrow);
        for(i=0;i<sparse_matrix->nrow;i++)
        {
            for(j=0;j<sparse_matrix->ndata[i];j++)
            {
                gmx_fio_do_int(fio, sparse_matrix->data[i][j].col);
                gmx_fio_do_real(fio, sparse_matrix->data[i][j].value);
            }
        }
    }
    gmx_fio_close(fio);
}


void
gmx_mtxio_read (const char *            filename,
                int *                   nrow,
                int *                   ncol,
                real **                 full_matrix,
                gmx_sparsematrix_t **   sparse_matrix)
{
    t_fileio  *fio;
    XDR *   xd;
    int     i,j,prec;
    gmx_bool    bDum = TRUE;
    gmx_bool    bRead = TRUE;
    char    gmxver[256];
    size_t  sz;
    
    fio = gmx_fio_open(filename,"r");
    gmx_fio_checktype(fio);
    xd = gmx_fio_getxdr(fio);
    
    /* Read and check magic number */
    i = GMX_MTXIO_MAGIC_NUMBER;
    gmx_fio_do_int(fio, i);

    if(i!=GMX_MTXIO_MAGIC_NUMBER)
    {
        gmx_fatal(FARGS,
                  "No matrix data found in file. Note that the Hessian matrix format changed\n"
                  "in Gromacs 3.3 to enable portable files and sparse matrix storage.\n");
    }
    
    /* Read generating Gromacs version */
    gmx_fio_do_string(fio, gmxver);
    
    /* Write 1 for double, 0 for single precision */
    if(sizeof(real)==sizeof(double))
        prec = 1;
    else
        prec = 0;
    gmx_fio_do_int(fio, prec);

    fprintf(stderr,"Reading %s precision matrix generated by Gromacs %s\n",
            (prec == 1) ? "double" : "single",gmxver);
    
    gmx_fio_do_int(fio, i);
    *nrow=i;
    gmx_fio_do_int(fio, i);
    *ncol=i;

    gmx_fio_do_int(fio, i);
    
    if(i==GMX_MTXIO_FULL_MATRIX)
    {
        printf("Full matrix storage format, nrow=%d, ncols=%d\n",*nrow,*ncol);

        sz = (*nrow) * (*ncol);
        snew((*full_matrix),sz);
        bDum=gmx_fio_ndo_real(fio, (*full_matrix),sz);
    }
    else
    {
        /* Sparse storage */
        printf("Sparse matrix storage format, nrow=%d, ncols=%d\n",*nrow,*ncol);

        snew((*sparse_matrix),1);
        gmx_fio_do_gmx_bool(fio, (*sparse_matrix)->compressed_symmetric);
        gmx_fio_do_int(fio, (*sparse_matrix)->nrow);        
        if((*sparse_matrix)->nrow != *nrow)
        {
            gmx_fatal(FARGS,"Internal inconsistency in sparse matrix.\n");
        }
        snew((*sparse_matrix)->ndata,(*sparse_matrix)->nrow);
        snew((*sparse_matrix)->nalloc,(*sparse_matrix)->nrow);
        snew((*sparse_matrix)->data,(*sparse_matrix)->nrow);
        bDum=gmx_fio_ndo_int(fio, (*sparse_matrix)->ndata,
                             (*sparse_matrix)->nrow);

        for(i=0;i<(*sparse_matrix)->nrow;i++)
        {
            (*sparse_matrix)->nalloc[i] = (*sparse_matrix)->ndata[i] + 10;
            snew(((*sparse_matrix)->data[i]),(*sparse_matrix)->nalloc[i]);
            
            for(j=0;j<(*sparse_matrix)->ndata[i];j++)
            {
                gmx_fio_do_int(fio, (*sparse_matrix)->data[i][j].col);
                gmx_fio_do_real(fio, (*sparse_matrix)->data[i][j].value);
            }
        }
    }
    gmx_fio_close(fio);
}



