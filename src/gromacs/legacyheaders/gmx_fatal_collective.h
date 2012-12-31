/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _fatal_collective_h
#define _fatal_collective_h

#include "types/simple.h"
#include "types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif


    void
    gmx_fatal_collective(int f_errno,const char *file,int line,
                         const t_commrec *cr,gmx_domdec_t *dd,
                         const char *fmt,...);
    /* As gmx_fatal declared in gmx_fatal.h,
     * but only the master process prints the error message.
     * This should only be called one of the following two situations:
     * 1) On all nodes in cr->mpi_comm_mysim, with cr!=NULL,dd==NULL.
     * 2) On all nodes in dd->mpi_comm_all,   with cr==NULL,dd!=NULL.
     * This will call MPI_Finalize instead of MPI_Abort when possible,
     * This is useful for handling errors in code that is executed identically
     * for all processes.
     */


#ifdef __cplusplus
}
#endif

#endif  /* _fatal_collective_h */
