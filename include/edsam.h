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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _edsam_h
#define _edsam_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern void ed_open(int nfile,t_filenm fnm[],t_edsamyn *edyn, t_commrec *cr);
extern void init_edsam(FILE *log,t_topology *top,
		       t_mdatoms *md,int start,int homenr, t_commrec *cr,
		       rvec x[],matrix box, 
		       t_edsamyn *edyn,t_edpar *edi);
extern void read_edi(t_edsamyn *edyn,t_edpar *edi,int nr_mdatoms);
extern int  read_edint(FILE *file);
extern int  read_edint2(FILE *file);
extern real read_edreal(FILE *file);
extern void read_edx(FILE *file,int number,int *anrs,rvec *x);
extern void read_edvecs(FILE *in,int nr,t_edvecs *vecs);
extern void read_edvec(FILE *in,int nr,t_eigvec *tvec);
extern void scan_edvec(FILE *in,int nr,rvec *vec);
extern real fitit(int nr, rvec *x,t_edpar *edi,rvec *transvec,
		       matrix rmat);
extern void do_edfit(int natoms,rvec *xp,rvec *x,matrix R); 
extern void put_in_origin(int nr,rvec *x,int nmass,int *masnrs,
			       real *mass,real tmass);
extern void project(rvec *x,t_edpar *edi,char *mode);
extern void do_project(rvec *x, t_eigvec *vec, t_edpar *edi,char *mode);
extern void projectx(t_edpar *edi,rvec *x,t_eigvec *vec);
extern real do_projectx(t_edpar *edi,rvec *x,rvec *vec);
extern real calc_radius(t_eigvec *vec);
extern void do_edsam(FILE *log,t_topology *top,t_inputrec *ir,int step,
		     t_mdatoms *md,int start,int homenr, t_commrec *cr,
                     rvec x[],rvec xold[],rvec x_unc[],rvec f[],matrix box,
                     t_edsamyn *edyn,t_edpar *edi,bool bHave_force);
extern void rmfit(int ned,rvec *x,rvec *transvec,matrix rotmat);
extern void rotate_vec(int nr,rvec *x,matrix rotmat);
extern void ed_cons(rvec *x,t_edpar *edi,int step);
extern void do_linfix(rvec *x,t_edpar *edi,int step);
extern void do_linacc(rvec *x,t_edpar *edi);
extern void do_radfix(rvec *x,t_edpar *edi,int step);
extern void do_radacc(rvec *x,t_edpar *edi);
extern void do_radcon(rvec *x,t_edpar *edi);
extern void write_edo(t_edpar *edi,int step,real rmsd);
extern void write_proj(FILE *out,t_edpar *edi,char *mode);
extern void do_write_proj(FILE *out,t_eigvec *vec,char *mode);
extern void write_edidx(FILE *out,t_edpar *edi);
#endif	/* _edsam_h */






