/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef GMX_FILEIO_TPXIO_H
#define GMX_FILEIO_TPXIO_H


/**************************************************************
 *
 * The routines in the corresponding c-file tpxio.c
 * are based on the lower level routines in gmxfio.c
 * The integer file pointer returned from open_tpx
 * can also be used with the routines in gmxfio.h
 *
 **************************************************************/
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/state.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_mtop_t;
struct t_atoms;
struct t_block;
struct t_topology;

typedef struct
{
    int   bIr;       /* Non zero if input_rec is present		*/
    int   bBox;      /* Non zero if a box is present			*/
    int   bTop;      /* Non zero if a topology is present		*/
    int   bX;        /* Non zero if coordinates are present		*/
    int   bV;        /* Non zero if velocities are present		*/
    int   bF;        /* Non zero if forces are present		*/

    int   natoms;    /* The total number of atoms			*/
    int   ngtc;      /* The number of temperature coupling groups    */
    real  lambda;    /* Current value of lambda			*/
    int   fep_state; /* Current value of the alchemical state --
                      * not yet printed out.  */
    /*a better decision will eventually (5.0 or later) need to be made
       on how to treat the alchemical state of the system, which can now
       vary through a simulation, and cannot be completely described
       though a single lambda variable, or even a single state
       index. Eventually, should probably be a vector. MRS*/
} t_tpxheader;

/*
 * These routines handle reading and writing of preprocessed
 * topology files in any of the following formats:
 * TPR : topology in XDR format, portable accross platforms
 * TRR : trajectory in XDR format (non compressed)
 *
 * Files are written in the precision with which the source are compiled,
 * but double and single precision can be read by either.
 */

t_fileio *open_tpx(const char *fn, const char *mode);
/* Return an file pointer corresponding to the file you have just opened */

void close_tpx(t_fileio *fio);
/*  Close the file corresponding to fio */

void read_tpxheader(const char *fn, t_tpxheader *tpx, gmx_bool TopOnlyOK,
                    int *version, int *generation);
/* Read the header from a tpx file and then close it again.
 * By setting TopOnlyOK to true, it is possible to read future
 * versions too (we skip the changed inputrec), provided we havent
 * changed the topology description. If it is possible to read
 * the inputrec it will still be done even if TopOnlyOK is TRUE.
 *
 * The version and generation if the topology (see top of tpxio.c)
 * are returned in the two last arguments.
 */

void write_tpx_state(const char *fn,
                     t_inputrec *ir, t_state *state, struct gmx_mtop_t *mtop);
/* Write a file, and close it again.
 */

void read_tpx_state(const char *fn,
                    t_inputrec *ir, t_state *state, rvec *f,
                    struct gmx_mtop_t *mtop);
int read_tpx(const char *fn,
             t_inputrec *ir, matrix box, int *natoms,
             rvec *x, rvec *v, rvec *f, struct gmx_mtop_t *mtop);
/* Read a file, and close it again.
 * When step, t or lambda are NULL they will not be stored.
 * Returns ir->ePBC, if it could be read from the file.
 */

int read_tpx_top(const char *fn,
                 t_inputrec *ir, matrix box, int *natoms,
                 rvec *x, rvec *v, rvec *f, struct t_topology *top);
/* As read_tpx, but for the old t_topology struct */

gmx_bool fn2bTPX(const char *file);
/* return if *file is one of the TPX file types */

gmx_bool read_tps_conf(const char *infile, char *title, struct t_topology *top,
                       int *ePBC, rvec **x, rvec **v, matrix box, gmx_bool bMass);
/* Read title, top.atoms, x, v (if not NULL) and box from an STX file,
 * memory for atoms, x and v will be allocated.
 * Return TRUE if a complete topology was read.
 * If infile is a TPX file read the whole top,
 * else if bMass=TRUE, read the masses into top.atoms from the mass database.
 */

void tpx_make_chain_identifiers(struct t_atoms *atoms, struct t_block *mols);

#ifdef __cplusplus
}
#endif

#endif
