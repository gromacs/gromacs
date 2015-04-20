/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_PBCUTIL_PBC_H
#define GMX_PBCUTIL_PBC_H

#include <stdio.h>

#include "gromacs/legacyheaders/types/commrec_fwd.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of combinations of single triclinic box vectors
 * required to shift atoms that are within a brick of the size of
 * the diagonal of the box to within the maximum cut-off distance.
 */
#define MAX_NTRICVEC 12

typedef struct t_pbc {
    int        ePBC;
    int        ndim_ePBC;
    int        ePBCDX;
    int        dim;
    matrix     box;
    rvec       fbox_diag;
    rvec       hbox_diag;
    rvec       mhbox_diag;
    real       max_cutoff2;
    gmx_bool   bLimitDistance;
    real       limit_distance2;
    int        ntric_vec;
    ivec       tric_shift[MAX_NTRICVEC];
    rvec       tric_vec[MAX_NTRICVEC];
} t_pbc;

#define TRICLINIC(box) (box[YY][XX] != 0 || box[ZZ][XX] != 0 || box[ZZ][YY] != 0)

#define NTRICIMG 14
#define NCUCVERT 24
#define NCUCEDGE 36

enum {
    ecenterTRIC, /* 0.5*(a+b+c)                  */
    ecenterRECT, /* (0.5*a[x],0.5*b[y],0.5*c[z]) */
    ecenterZERO, /* (0,0,0)                      */
    ecenterDEF = ecenterTRIC
};

struct t_graph;

int ePBC2npbcdim(int ePBC);
/* Returns the number of dimensions that use pbc, starting at X */

int inputrec2nboundeddim(t_inputrec *ir);
/* Returns the number of dimensions in which
 * the coordinates of the particles are bounded, starting at X.
 */

void dump_pbc(FILE *fp, t_pbc *pbc);
/* Dump the contents of the pbc structure to the file */

const char *check_box(int ePBC, matrix box);
/* Returns NULL if the box is supported by Gromacs.
 * Otherwise is returns a string with the problem.
 * When ePBC=-1, the type of pbc is guessed from the box matrix.
 */

real max_cutoff2(int ePBC, matrix box);
/* Returns the square of the maximum cut-off allowed for the box,
 * taking into account that the grid neighborsearch code and pbc_dx
 * only check combinations of single box-vector shifts.
 */

int guess_ePBC(matrix box);
/* Guesses the type of periodic boundary conditions using the box */

gmx_bool correct_box(FILE *fplog, int step, tensor box, struct t_graph *graph);
/* Checks for un-allowed box angles and corrects the box
 * and the integer shift vectors in the graph (if graph!=NULL) if necessary.
 * Returns TRUE when the box was corrected.
 */

int ndof_com(t_inputrec *ir);
/* Returns the number of degrees of freedom of the center of mass */

void set_pbc(t_pbc *pbc, int ePBC, matrix box);
/* Initiate the periodic boundary conditions.
 * pbc_dx will not use pbc and return the normal difference vector
 * when one or more of the diagonal elements of box are zero.
 * When ePBC=-1, the type of pbc is guessed from the box matrix.
 */

t_pbc *set_pbc_dd(t_pbc *pbc, int ePBC,
                  gmx_domdec_t *dd, gmx_bool bSingleDir, matrix box);
/* As set_pbc, but additionally sets that correct distances can
 * be obtained using (combinations of) single box-vector shifts.
 * Should be used with pbc_dx_aiuc.
 * If dd!=NULL pbc is not used for directions
 * with dd->nc[i]==1 with bSingleDir==TRUE or
 * with dd->nc[i]<=2 with bSingleDir==FALSE.
 * Returns pbc when pbc operations are required, NULL otherwise.
 */

void pbc_dx(const t_pbc *pbc, const rvec x1, const rvec x2, rvec dx);
/* Calculate the correct distance vector from x2 to x1 and put it in dx.
 * set_pbc must be called before ever calling this routine.
 *
 * For triclinic boxes pbc_dx does not necessarily return the shortest
 * distance vector. If pbc->bLimitDistance=TRUE an atom pair with
 * distance vector dx with norm2(dx) > pbc->limit_distance2 could
 * have a shorter distance, but not shorter than sqrt(pbc->limit_distance2).
 * pbc->limit_distance2 is always larger than max_cutoff2(box).
 * For the standard rhombic dodecahedron and truncated octahedron
 * pbc->bLimitDistance=FALSE and thus all distances are correct.
 */

int pbc_dx_aiuc(const t_pbc *pbc, const rvec x1, const rvec x2, rvec dx);
/* Calculate the correct distance vector from x2 to x1 and put it in dx,
 * This function can only be used when all atoms are in the rectangular
 * or triclinic unit-cell.
 * Returns the ishift required to shift x1 at closest distance to x2;
 * i.e. if 0<=ishift<SHIFTS then x1 - x2 + shift_vec[ishift] = dx
 * (see calc_shifts below on how to obtain shift_vec)
 * set_pbc_dd or set_pbc must be called before ever calling this routine.
 */
void pbc_dx_d(const t_pbc *pbc, const dvec x1, const dvec x2, dvec dx);
/* As pbc_dx, but for double precision vectors.
 * set_pbc must be called before ever calling this routine.
 */

gmx_bool image_rect(ivec xi, ivec xj, ivec box_size,
                    real rlong2, int *shift, real *r2);
/* Calculate the distance between xi and xj for a rectangular box.
 * When the distance is SMALLER than rlong2 return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

gmx_bool image_tri(ivec xi, ivec xj, imatrix box,
                   real rlong2, int *shift, real *r2);
/* Calculate the distance between xi and xj for a triclinic box.
 * When the distance is SMALLER than rlong2 return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

gmx_bool image_cylindric(ivec xi, ivec xj, ivec box_size, real rlong2,
                         int *shift, real *r2);
/* Calculate the distance between xi and xj for a rectangular box
 * using a cylindric cutoff for long-range only.
 * When the distance is SMALLER than rlong2 (in X and Y dir.)
 * return TRUE, return
 * the shift code in shift and the distance in r2. When the distance is
 * >= rlong2 return FALSE;
 * It is assumed that rlong2 is scaled the same way as the ivecs xi and xj.
 */

void calc_shifts(matrix box, rvec shift_vec[]);
/* This routine calculates ths shift vectors necessary to use the
 * ns routine.
 */

void calc_box_center(int ecenter, matrix box, rvec box_center);
/* Calculates the center of the box.
 * See the description for the enum ecenter above.
 */

void calc_triclinic_images(matrix box, rvec img[]);
/* Calculates the NTRICIMG box images */

void calc_compact_unitcell_vertices(int ecenter, matrix box,
                                    rvec vert[]);
/* Calculates the NCUCVERT vertices of a compact unitcell */

int *compact_unitcell_edges(void);
/* Return an array of unitcell edges of length NCUCEDGE*2,
 * this is an index in vert[], which is calculated by calc_unitcell_vertices.
 * The index consists of NCUCEDGE pairs of vertex indices.
 * The index does not change, so it needs to be retrieved only once.
 */

void put_atoms_in_box_omp(int ePBC, matrix box, int natoms, rvec x[]);
/* This wrapper function around put_atoms_in_box() with the ugly manual
 * workload splitting is needed toavoid silently introducing multithreading
 * in tools.
 * */


void put_atoms_in_box(int ePBC, matrix box, int natoms, rvec x[]);
/* These routines puts ONE or ALL atoms in the box, not caring
 * about charge groups!
 * Also works for triclinic cells.
 */

void put_atoms_in_triclinic_unitcell(int ecenter, matrix box,
                                     int natoms, rvec x[]);
/* This puts ALL atoms in the triclinic unit cell, centered around the
 * box center as calculated by calc_box_center.
 */

const char *put_atoms_in_compact_unitcell(int ePBC, int ecenter,
                                          matrix box,
                                          int natoms, rvec x[]);
/* This puts ALL atoms at the closest distance for the center of the box
 * as calculated by calc_box_center.
 * Will return NULL is everything went ok and a warning string if not
 * all atoms could be placed in the unitcell. This can happen for some
 * triclinic unitcells, see the comment at pbc_dx above.
 * When ePBC=-1, the type of pbc is guessed from the box matrix.
 */

#ifdef __cplusplus
}
#endif

#endif
