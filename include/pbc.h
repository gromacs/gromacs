/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _pbc_h
#define _pbc_h

static char *SRCID_pbc_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) pbc.h 1.17 2/2/97"
#endif /* HAVE_IDENT */
#include "sysstuff.h"
#include "typedefs.h"

#ifdef CPLUSPLUS
extern "C" { 
#endif

#define BOX_MARGIN 1.001
/* margin factor for checking if the box is too skewed */

#define TRICLINIC(box) (box[YY][XX]!=0 || box[ZZ][XX]!=0 || box[ZZ][YY]!=0)

#define NTRICIMG 14
#define NCUCVERT 24
#define NCUCEDGE 36

extern char *check_box(matrix box);
/* Returns NULL if the box is supported by Gromacs.
 * Otherwise is returns a string with the problem.
 */

extern void init_pbc(matrix box,bool bTruncOct);
/* Initiate the periodic boundary conditions. Set bTruncOct to
 * TRUE when using a truncated octahedron box.
 */

extern void pbc_dx(rvec x1, rvec x2, rvec dx);
/* Calculate the correct distance vector from x1 and x2 and put it in
 * dx. init_pbc must be called before ever calling this routine
 * (this is done by put_charge_groups_in_box).
 */

extern bool image_rect(ivec xi,ivec xj,ivec box_size,
		       real rlong2,int *shift,real *r2);
/* Calculate the distance between xi and xj for a rectangular box.
 * When the distance is SMALLER than rlong2 return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

extern bool image_tri(ivec xi,ivec xj,imatrix box,
		      real rlong2,int *shift,real *r2);
/* Calculate the distance between xi and xj for a triclinic box.
 * When the distance is SMALLER than rlong2 return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

extern bool image_cylindric(ivec xi,ivec xj,ivec box_size,real rlong2,
			    int *shift,real *r2);
/* Calculate the distance between xi and xj for a rectangular box
 * using a cylindric cutoff for long-range only.
 * When the distance is SMALLER than rlong2 (in X and Y dir.)
 * return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

extern void calc_shifts(matrix box,rvec box_size,rvec shift_vec[],
			bool bTruncOct);
/* This routine calculates ths shift vectors necessary to use the
 * ns routine. Note that for the truncated octahedron case too many
 * shift vectors can be calculated: The ones for which exactly
 * 2 of the k,l,m indexes are not 0 (12 vectors lying along the box
 * edges. This can be compensated for by removing all the shift_vecs with
 * (k+l+m) even. This is a feature of the way in which the counting is 
 * done. It implies that, when using truncated octahedron,
 * the shift codes 1,3,5,7,9,11,15,17,19,21,23,25 should never occur,
 * that is, every second entry, EXCEPT the central box.
 */

extern void calc_cgcm(FILE *log,int cg0,int cg1,t_block *cgs,
		      rvec pos[],rvec cg_cm[]);
/* Routine to compute centers of geometry of charge groups. No periodicity
 * is used.
 */
   
extern void put_charge_groups_in_box (FILE *log,int cg0,int cg1,bool bTruncOct,
				      matrix box,rvec box_size,t_block *cgs,
				      rvec pos[],
				      rvec cg_cm[]);
			    
/* This routine puts charge groups in the periodic box, keeping them
 * together. When bTruncOct==TRUE a truncated octahedron
 * box is used. There are no checks: the first element of the box matrix
 * is taken to be the box edge.
 */

extern void calc_box_center(matrix box,rvec box_center);
/* Calculates the center of the box */

extern void calc_triclinic_images(matrix box,rvec img[]);
/* Calculates the NTRICIMG box images */

extern void calc_compact_unitcell_vertices(matrix box,rvec vert[]);
/* Calculates the NCUCVERT vertices of a compact unitcell */

extern int *compact_unitcell_edges(void);
/* Return an array of unitcell edges of length NCUCEDGE*2,
 * this is an index in vert[], which is calculated by calc_unitcell_vertices.
 * The index consists of NCUCEDGE pairs of vertex indices.
 * The index does not change, so it needs to be retrieved only once.
 */

extern void put_atoms_in_box(matrix box,int natoms,rvec x[]);
/* This puts ALL atoms in the box, not caring about charge groups!
 * Also works for triclinic cells.
 */

extern void put_atoms_in_triclinic_unitcell(matrix box,int natoms,rvec x[]);
/* This puts ALL atoms in the triclinic unit cell, centered around the
 * box center as calculated by calc_box_center.
 */

extern void put_atoms_in_compact_unitcell(matrix box,int natoms,rvec x[]);
/* This puts ALL atoms at the closest distance for the center of the box
 * as calculated by calc_box_center.
 */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _pbc_h */
