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

#ifndef GMX_FILEIO_TRXIO_H
#define GMX_FILEIO_TRXIO_H

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/readinp.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0 /* avoid screwing up indentation */
}
#endif

struct gmx_mtop_t;
struct t_atoms;
struct t_topology;
struct t_trxframe;

/* a dedicated status type contains fp, etc. */
typedef struct t_trxstatus t_trxstatus;

/* I/O function types */

/************************************************
 *             Trajectory functions
 ************************************************/

int prec2ndec(real prec);
/* Convert precision in 1/(nm) to number of decimal places */

/*! \brief Convert number of decimal places to trajectory precision in
 * 1/(nm) */
real ndec2prec(int ndec);

void clear_trxframe(struct t_trxframe *fr, gmx_bool bFirst);
/* Set all content gmx_booleans to FALSE.
 * When bFirst = TRUE, set natoms=-1, all pointers to NULL
 *                     and all data to zero.
 */

void set_trxframe_ePBC(struct t_trxframe *fr, int ePBC);
/* Set the type of periodic boundary conditions, ePBC=-1 is not set */

int nframes_read(t_trxstatus *status);
/* Returns the number of frames read from the trajectory */

int write_trxframe_indexed(t_trxstatus *status, struct t_trxframe *fr, int nind,
                           const atom_id *ind, gmx_conect gc);
/* Write an indexed frame to a TRX file, see write_trxframe. gc may be NULL */

int write_trxframe(t_trxstatus *status, struct t_trxframe *fr, gmx_conect gc);
/* Write a frame to a TRX file.
 * Only entries for which the gmx_boolean is TRUE will be written,
 * except for step, time, lambda and/or box, which may not be
 * omitted for certain trajectory formats.
 * The precision for .xtc and .gro is fr->prec, when fr->bPrec=FALSE,
 * the precision is set to 1000.
 * gc is important for pdb file writing only and may be NULL.
 */

int write_trx(t_trxstatus *status, int nind, const atom_id *ind, struct t_atoms *atoms,
              int step, real time, matrix box, rvec x[], rvec *v,
              gmx_conect gc);
/* Write an indexed frame to a TRX file.
 * v can be NULL.
 * atoms can be NULL for file types which don't need atom names.
 */

void trjtools_gmx_prepare_tng_writing(const char               *filename,
                                      char                      filemode,
                                      t_trxstatus              *in,
                                      t_trxstatus             **out,
                                      const char               *infile,
                                      const int                 natoms,
                                      const struct gmx_mtop_t  *mtop,
                                      const atom_id            *index,
                                      const char               *index_group_name);
/* Sets up *out for writing TNG. If *in != NULL and contains a TNG trajectory
 * some data, e.g. molecule system, will be copied over from *in to *out.
 * If *in == NULL a file name (infile) of a TNG file can be provided instead
 * and used for copying data to *out.
 * If there is no TNG input natoms is used to create "implicit atoms" (no atom
 * or molecular data present). If natoms == -1 the number of atoms are
 * not known (or there is already a TNG molecule system to copy, in which case
 * natoms is not required anyhow). If an group of indexed atoms are written
 * natoms must be the length of index. index_group_name is the name of the
 * index group.
 */

/*! \brief Write a trxframe to the TNG file in status.
 *
 * This function is needed because both t_trxstatus and
 * tng_trajectory_t are encapsulated, so client trajectory-writing
 * code with a t_trxstatus can't just call the TNG writing
 * function. */
void write_tng_frame(t_trxstatus        *status,
                     struct t_trxframe  *fr);

void close_trx(t_trxstatus *status);
/* Close trajectory file as opened with read_first_x, read_frist_frame
 * or open_trx. Identical to close_trj.
 */

t_trxstatus *open_trx(const char *outfile, const char *filemode);
/* Open a TRX file and return an allocated status pointer */

t_fileio *trx_get_fileio(t_trxstatus *status);
/* get a fileio from a trxstatus */

float trx_get_time_of_final_frame(t_trxstatus *status);
/* get time of final frame. Only supported for TNG and XTC */

gmx_bool bRmod_fd(double a, double b, double c, gmx_bool bDouble);
/* Returns TRUE when (a - b) MOD c = 0, using a margin which is slightly
 * larger than the float/double precision.
 */

#ifdef GMX_DOUBLE
#define bRmod(a, b, c) bRmod_fd(a, b, c, TRUE)
#else
#define bRmod(a, b, c) bRmod_fd(a, b, c, FALSE)
#endif

int check_times2(real t, real t0, gmx_bool bDouble);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin or t MOD dt = t0,
 *          0 if tbegin <= t <=tend+margin,
 *          1 if t>tend
 * where margin is 0.1*min(t-tp,tp-tpp), if this positive, 0 otherwise.
 * tp and tpp should be the time of the previous frame and the one before.
 * The mod is done with single or double precision accuracy depending
 * on the value of bDouble.
 */

int check_times(real t);
/* This routine checkes if the read-in time is correct or not;
 * returns -1 if t<tbegin,
 *          0 if tbegin <= t <=tend,
 *          1 if t>tend
 */





/* For trxframe.flags, used in trxframe read routines.
 * When a READ flag is set, the field will be read when present,
 * but a frame might be returned which does not contain the field.
 * When a NEED flag is set, frames not containing the field will be skipped.
 */
#define TRX_READ_X    (1<<0)
#define TRX_NEED_X    (1<<1)
#define TRX_READ_V    (1<<2)
#define TRX_NEED_V    (1<<3)
#define TRX_READ_F    (1<<4)
#define TRX_NEED_F    (1<<5)
/* Useful for reading natoms from a trajectory without skipping */
#define TRX_DONT_SKIP (1<<6)

/* For trxframe.not_ok */
#define HEADER_NOT_OK (1<<0)
#define DATA_NOT_OK   (1<<1)
#define FRAME_NOT_OK  (HEADER_NOT_OK | DATA_NOT_OK)

int read_first_frame(const output_env_t oenv, t_trxstatus **status,
                     const char *fn, struct t_trxframe *fr, int flags);
/* Read the first frame which is in accordance with flags, which are
 * defined further up in this file.
 * Returns natoms when succeeded, 0 otherwise.
 * Memory will be allocated for flagged entries.
 * The flags are copied to fr for subsequent calls to read_next_frame.
 * Returns TRUE when succeeded, FALSE otherwise.
 */

gmx_bool read_next_frame(const output_env_t oenv, t_trxstatus *status,
                         struct t_trxframe *fr);
/* Reads the next frame which is in accordance with fr->flags.
 * Returns TRUE when succeeded, FALSE otherwise.
 */

int read_first_x(const output_env_t oenv, t_trxstatus **status,
                 const char *fn, real *t, rvec **x, matrix box);
/* These routines read first coordinates and box, and allocates
 * memory for the coordinates, for a trajectory file.
 * The routine returns the number of atoms, or 0 when something is wrong.
 * The integer in status should be passed to calls of read_next_x
 */

gmx_bool read_next_x(const output_env_t oenv, t_trxstatus *status, real *t, rvec x[], matrix box);
/* Read coordinates and box from a trajectory file. Return TRUE when all well,
 * or FALSE when end of file (or last frame requested by user).
 * status is the integer set in read_first_x.
 */

void close_trj(t_trxstatus *status);
/* Close trajectory file as opened with read_first_x, read_first_frame
 * or open_trx. Identical to close_trx.
 */

void rewind_trj(t_trxstatus *status);
/* Rewind trajectory file as opened with read_first_x */

struct t_topology *read_top(const char *fn, int *ePBC);
/* Extract a topology data structure from a topology file.
 * If ePBC!=NULL *ePBC gives the pbc type.
 */

#ifdef __cplusplus
}
#endif

#endif
