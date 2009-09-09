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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _enxio_h
#define _enxio_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef __cplusplus
external "C" {
#endif

  /**************************************************************
   *
   * The routines in the corresponding c-file enxio.c
   * are based on the lower level routines in gmxfio.c
   * The integer file pointer returned from open_enx
   * can also be used with the routines in gmxfio.h
   *
   **************************************************************/

#include "sysstuff.h"
#include "typedefs.h"
#include "pbc.h"
  
  typedef struct {
    char *name;
    char *unit;
  } gmx_enxnm_t;
  
  /* 
   * Index for the additional blocks in the energy file.
   * Blocks can be added without sacrificing backward and forward
   * compatibility of the energy files.
   */
  enum {
    enxOR,   /* Time and ensemble averaged data for orientation restraints */
    enxORI,  /* Instantaneous data for orientation restraints              */
    enxORT,  /* Order tensor(s) for orientation restraints                 */
    enxNR    /* Total number of extra blocks in the current code,
              * note that the enxio code can read files written by
	      * future code which contain more blocks.
	      */
  };

  typedef struct {
    double   t;	            /* Timestamp of this frame	                     */
    gmx_large_int_t step;        /* MD step	   		                     */
    gmx_large_int_t nsteps;      /* The number of steps between frames            */
    gmx_large_int_t nsum;        /* The number of terms for the sums in ener      */
    int      nre;           /* Number of energies			     */
    int      ndisre;        /* Number of distance restraints	             */
    int      nblock;        /* Number of following energy blocks              */
    int      *nr;           /* Number of things in additional blocks (nblock) */
    int      e_size;        /* Size (in bytes) of energies		     */
    int      d_size;        /* Size (in bytes) of disre blocks              */
    int      nr_alloc;      /* Allocated size of nr and block                 */
    int      e_alloc;       /* Allocated size (in elements) of ener           */
    int      d_alloc;       /* Allocated size (in elements) of rav and rt     */
    int      *b_alloc;      /* Allocated size (in elements) of each block     */
    t_energy *ener;         /* The energies                                   */
    real     *disre_rm3tav; /* Time averaged data for distance restraints     */
    real     *disre_rt;     /* Instantaneous data for distance restraints     */
    real     **block;       /* Additional energy blocks (nblock x b_alloc[b]) */
  } t_enxframe;

  /* file handle */
  typedef struct ener_file *ener_file_t;

  /* 
   * An energy file is read like this:
   *
   * ener_file_t fp;
   * t_enxframe *fr;
   *
   * fp = open_enx(...);
   * do_enxnms(fp,...);
   * snew(fr,1);
   * while (do_enx(fp,fr)) {
   * ...
   * }
   * free_enxframe(fr);
   * sfree(fr);
   */
  
  /* New energy reading and writing interface */
  extern void free_enxframe(t_enxframe *fr);
  /* Frees all allocated memory in fr */

  extern ener_file_t open_enx(const char *fn,const char *mode);

  extern int enx_file_pointer(const ener_file_t ef);

  extern void close_enx(ener_file_t ef);
  
  extern void do_enxnms(ener_file_t ef,int *nre,gmx_enxnm_t **enms);
  
  extern void free_enxnms(int n,gmx_enxnm_t *nms);
  /* Frees nms and all strings in it */

  extern bool do_enx(ener_file_t ef,t_enxframe *fr);
  /* Reads enx_frames, memory in fr is (re)allocated if necessary */

  extern void get_enx_state(const char *fn, real t,
			    gmx_groups_t *groups, t_inputrec *ir,
			    t_state *state);
  /*
   * Reads state variables from enx file fn at time t.
   * atoms and ir are required for determining which things must be read.
   * Currently pcoupl and tcoupl state are read from enx.
   */
  
#ifdef __cplusplus
}
#endif

#endif	/* _enerio_h */
