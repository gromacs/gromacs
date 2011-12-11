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



/* This is the new improved and thread safe version of gmxfio.  */


/* WARNING WARNING WARNING WARNING
   The data types used here are PRIVATE to gmxfio routines. DO NOT use them
   directly in your own code, but use the external functions provided in 
   include/gmxfio.h 

   If you don't heed this warning, your code will suddenly stop working 
   at some point in the not-so-distant future. 

   WARNING WARNING WARNING WARNING */


/* XDR should be available on all platforms now, 
 * but we keep the possibility of turning it off...
 */
#define USE_XDR



/* the reader/writer functions  for t_iotype */
typedef gmx_bool read_func(t_fileio *fio, void *item, int nitem, int eio,
                       const char *desc,const char *srcfile,int line);
typedef gmx_bool write_func(t_fileio *fio, const void *item, int nitem, int eio,
                        const char *desc,const char *srcfile,int line);


/* these are pointers to the actual reading & writing functions */
typedef struct
{
    read_func *nread;
    write_func *nwrite;
} t_iotype;



struct t_fileio
{
    FILE *fp; /* the file pointer */
    const t_iotype *iotp;  /* file type */
    gmx_bool bOpen,  /* the file is open */
         bRead,  /* the file is open for reading */
         bDouble, /* write doubles instead of floats */
         bDebug, /* the file ops should come with debug info */
         bStdio, /* the file is actually stdin or stdout */
         bLargerThan_off_t,  /* if the file position is largen than off_t 
                                could hold */
         bReadWrite; /* the file is open for reading and writing */
    char *fn; /* the file name */
    XDR *xdr; /* the xdr data pointer */
    enum xdr_op xdrmode; /* the xdr mode */
    int iFTP; /* the file type identifier */

    const char *comment; /* a comment string for debugging */

    t_fileio *next, *prev; /* next and previous file pointers in the
                              linked list */
#ifdef GMX_THREADS
    tMPI_Lock_t  mtx;  /* content locking mutex. This is a fast lock
                          for performance reasons: in some cases every
                          single byte that gets read/written requires
                          a lock */
#endif
}; 



extern const t_iotype asc_iotype;
extern const t_iotype bin_iotype;
extern const t_iotype xdr_iotype;
extern const t_iotype dummy_iotype;

extern const char *eioNames[eioNR];



#define GMX_FIO_BUFLEN 256

/* make a debug string if that is requested in the fio */
const char *gmx_fio_dbgstr(t_fileio *fio, const char *desc, char *buf);
/* check the number of items against the allowed number of items */
void gmx_fio_check_nitem(t_fileio *fio, int eio, int nitem, const char *file, 
                         int line);
/* check the output type against allowed values */
void gmx_fio_fe(t_fileio *fio, int eio, const char *desc, const char *srcfile, 
                int line);

/* lock/unlock the mutex associated with a fio  */
void gmx_fio_lock(t_fileio *fio);
void gmx_fio_unlock(t_fileio *fio);



