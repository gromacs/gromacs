/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif


/*! \brief Lock innerloop mutex and read lists indices
 *
 *  \internal
 *
 *  This routine is only used when both Fortran innerloops
 *  and threads are enabled.
 *
 *  Since the Fortran77 standard does not include support
 *  for POSIX threads, we call this routine instead which
 *
 *    -  Locks the provided mutex
 *    -  Reads the counter from memory
 *    -  Advances the counter in successively smaller chunks
 *    -  Releases the mutex
 *
 *  In other words, it performs exactly the same action as
 *  we do natively in the nonbonded kernel outer loop 
 *  when using C language for the kernels.
 *
 * Fortran does not know anything about a mutex, but since
 * arguments are passed by reference we mask it as a pointer
 * to an integer in the Fortran code.
 *
 * \param mtx       Pointer to the mutex to use, masked as int
 * \param count     Pointer to the outer loop counter
 * \param nri       Total umber of (outer loop) neighborlists
 * \param nthreads  Number of working threads
 * \param nn0       Returned value: Low index to use for outerloop 
 * \param nn1       Returned value: High index to use for outerloop
 *
 * \warning There is one possible cause of problems. Some
 *          fortran compilers make all variables static by
 *          default, and that will obviously screw up
 *          multithreading in a major way. If your Fortran
 *          compiler does this you only have two alternatives:
 *          Either find the flag to turn it off, or compile
 *          Gromacs without any Fortran nonbonded kernels.
 */
void 
F77_FUNC_(nb_kernel_f77sync,NB_KERNEL_F77SYNC)
     (int *       mtx,
      int *       count,
      int *       nri,
      int *       nthreads,
      int *       nn0,
      int *       nn1)
{
  int n0,n1;
#ifdef HAVE_PTHREADS_H  
  pthread_mutex_lock((pthread_mutex_t *)mtx);
#endif
  n0              = *count;
  /* Take successively smaller chunks */
  n1              = n0+((*nri)-n0)/(2*(*nthreads))+3;
  *count          = n1;
#ifdef HAVE_PTHREADS_H
  pthread_mutex_unlock((pthread_mutex_t *)mtx);
#endif
  n0              = *count;
  *nn0            = n0;
  *nn1            = n1;
}
