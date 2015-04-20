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

#ifndef GMX_FILEIO_GMXFIO_H
#define GMX_FILEIO_GMXFIO_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
/* types */


/* Enumerated for different items in files */
enum {
    eitemHEADER, eitemIR, eitemBOX,
    eitemTOP, eitemX, eitemV, eitemF, eitemNR
};

/* Enumerated for data types in files */
enum {
    eioREAL, eioFLOAT, eioDOUBLE, eioINT, eioINT64,
    eioUCHAR, eioNUCHAR, eioUSHORT,
    eioRVEC, eioNRVEC, eioIVEC, eioSTRING, eioNR
};

typedef struct t_fileio t_fileio;

/* NOTE ABOUT THREAD SAFETY:

   The functions are all thread-safe, provided that two threads don't
   do something silly like closing the same file, or one thread
   accesses a file that has been closed by another.
 */

/********************************************************
 * Open and Close
 ********************************************************/

t_fileio *gmx_fio_open(const char *fn, const char *mode);
/* Open a new file for reading or writing.
 * The file type will be deduced from the file name.
 */

int gmx_fio_close(t_fileio *fp);
/* Close the file corresponding to fp (if not stdio)
 * The routine will exit when an invalid fio is handled.
 * Returns 0 on success.
 */

int gmx_fio_fp_close(t_fileio *fp);
/* Close the file corresponding to fp without closing the FIO entry
 * Needed e.g. for trxio because the FIO entries are used to store
 * additional data.
 * NOTE that the fp still needs to be properly closed with gmx_fio_close().
 * The routine will exit when an invalid fio is handled.
 * Returns 0 on success.
 */


/* Open a file, return a stream, record the entry in internal FIO object */
FILE* gmx_fio_fopen(const char *fn, const char *mode);

/* Close a file previously opened with gmx_fio_fopen.
 * Do not mix these calls with standard fopen/fclose ones!
 * Returns 0 on success.  */
int gmx_fio_fclose(FILE *fp);



/********************************************************
 * Change properties of the open file
 ********************************************************/

void gmx_fio_setprecision(t_fileio *fio, gmx_bool bDouble);
/* Select the floating point precision for reading and writing files */

char *gmx_fio_getname(t_fileio *fio);
/* Return the filename corresponding to the fio index */

int gmx_fio_getftp(t_fileio *fio);
/* Return the filetype corresponding to the fio index.
    There is as of now no corresponding setftp function because the file
    was opened as a specific file type and changing that midway is most
    likely an evil hack. */

void gmx_fio_setdebug(t_fileio *fio, gmx_bool bDebug);
/* Set the debug mode */

gmx_bool gmx_fio_getdebug(t_fileio *fio);
/* Return  whether debug mode is on in fio  */

gmx_bool gmx_fio_getread(t_fileio *fio);
/* Return  whether read mode is on in fio  */


void gmx_fio_checktype(t_fileio *fio);
/* Check whether the fio is of a sane type */

/***************************************************
 * FILE Operations
 ***************************************************/

void gmx_fio_rewind(t_fileio *fio);
/* Rewind the file in fio */

int gmx_fio_flush(t_fileio *fio);
/* Flush the fio, returns 0 on success */

int gmx_fio_fsync(t_fileio *fio);
/* fsync the fio, returns 0 on success.
   NOTE: don't use fsync function unless you're absolutely sure you need it
   because it deliberately interferes with the OS's caching mechanisms and
   can cause dramatically slowed down IO performance. Some OSes (Linux,
   for example), may implement fsync as a full sync() point. */

gmx_off_t gmx_fio_ftell(t_fileio *fio);
/* Return file position if possible */

int gmx_fio_seek(t_fileio *fio, gmx_off_t fpos);
/* Set file position if possible, quit otherwise */

FILE *gmx_fio_getfp(t_fileio *fio);
/* Return the file pointer itself */


/* Element with information about position in a currently open file.
 * gmx_off_t should be defined by autoconf if your system does not have it.
 * If you do not have it on some other platform you do not have largefile
 * support at all, and you can define it to int (or better, find out how to
 * enable large files).  */
typedef struct
{
    char          filename[STRLEN];
    gmx_off_t     offset;
    unsigned char chksum[16];
    int           chksum_size;
}
gmx_file_position_t;

int gmx_fio_get_output_file_positions(gmx_file_position_t ** outputfiles,
                                      int                   *nfiles );
/* Return the name and file pointer positions for all currently open
 * output files. This is used for saving in the checkpoint files, so we
 * can truncate output files upon restart-with-appending.
 *
 * For the first argument you should use a pointer, which will be set to
 * point to a list of open files.
 */

t_fileio *gmx_fio_all_output_fsync(void);
/* fsync all open output files. This is used for checkpointing, where
   we need to ensure that all output is actually written out to
   disk.
   This is most important in the case of some networked file systems,
   where data is not synced with the file server until close() or
   fsync(), so data could remain in cache for days.
   Note the caveats reported with gmx_fio_fsync().

    returns: NULL if no error occurred, or a pointer to the first file that
             failed if an error occurred
 */


int gmx_fio_get_file_md5(t_fileio *fio, gmx_off_t offset,
                         unsigned char digest[]);


int xtc_seek_frame(t_fileio *fio, int frame, int natoms);

int xtc_seek_time(t_fileio *fio, real time, int natoms, gmx_bool bSeekForwardOnly);


/* Add this to the comment string for debugging */
void gmx_fio_set_comment(t_fileio *fio, const char *comment);

/* Remove previously set comment */
void gmx_fio_unset_comment(t_fileio *fio);




/********************************************************
 * Read and write
 ********************************************************/


/* basic reading & writing */
gmx_bool gmx_fio_reade_real(t_fileio *fio, real *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_float(t_fileio *fio, float *item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_double(t_fileio *fio, double *item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_int(t_fileio *fio, int *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_int64(t_fileio *fio, gmx_int64_t *item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_uchar(t_fileio *fio, unsigned char *item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_ushort(t_fileio *fio, unsigned short *item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_rvec(t_fileio *fio, rvec *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_ivec(t_fileio *fio, ivec *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_reade_string(t_fileio *fio, char *item,
                              const char *desc, const char *srcfile, int line);

gmx_bool gmx_fio_writee_real(t_fileio *fio, real item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_float(t_fileio *fio, float item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_double(t_fileio *fio, double item,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_int(t_fileio *fio, int item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_int64(t_fileio *fio, gmx_int64_t item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_uchar(t_fileio *fio, unsigned char item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_ushort(t_fileio *fio, unsigned short item,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_rvec(t_fileio *fio, rvec *item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_ivec(t_fileio *fio, ivec *item,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_writee_string(t_fileio *fio, const char *item,
                               const char *desc, const char *srcfile, int line);

/* reading or writing, depending on the file's opening mode string */
gmx_bool gmx_fio_doe_real(t_fileio *fio, real *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_float(t_fileio *fio, float *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_double(t_fileio *fio, double *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_gmx_bool(t_fileio *fio, gmx_bool *item,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_int(t_fileio *fio, int *item,
                         const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_int64(t_fileio *fio, gmx_int64_t *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_uchar(t_fileio *fio, unsigned char *item,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_ushort(t_fileio *fio, unsigned short *item,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_rvec(t_fileio *fio, rvec *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_ivec(t_fileio *fio, ivec *item,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_doe_string(t_fileio *fio, char *item,
                            const char *desc, const char *srcfile, int line);




/* array reading & writing */
gmx_bool gmx_fio_nreade_real(t_fileio *fio, real *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_float(t_fileio *fio, float *item, int n,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_double(t_fileio *fio, double *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_int(t_fileio *fio, int *item, int n,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_int64(t_fileio *fio, gmx_int64_t *item, int n,
                              const char *desc, const char *srcfile,
                              int line);
gmx_bool gmx_fio_nreade_uchar(t_fileio *fio, unsigned char *item, int n,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_ushort(t_fileio *fio, unsigned short *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_rvec(t_fileio *fio, rvec *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_ivec(t_fileio *fio, ivec *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nreade_string(t_fileio *fio, char *item[], int n,
                               const char *desc, const char *srcfile, int line);

gmx_bool gmx_fio_nwritee_real(t_fileio *fio, const real *item, int n,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_float(t_fileio *fio, const float *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_double(t_fileio *fio, const double *item, int n,
                                const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_int(t_fileio *fio, const int *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_int64(t_fileio *fio,
                               const gmx_int64_t *item, int n,
                               const char *desc, const char *srcfile,
                               int line);
gmx_bool gmx_fio_nwritee_uchar(t_fileio *fio, const unsigned char *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_ushort(t_fileio *fio, const unsigned short *item, int n,
                                const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_rvec(t_fileio *fio, const rvec *item, int n,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_ivec(t_fileio *fio, const ivec *item, int n,
                              const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_nwritee_string(t_fileio *fio, const char *item[], int n,
                                const char *desc, const char *srcfile, int line);

gmx_bool gmx_fio_ndoe_real(t_fileio *fio, real *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_float(t_fileio *fio, float *item, int n,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_double(t_fileio *fio, double *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_gmx_bool(t_fileio *fio, gmx_bool *item, int n,
                               const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_int(t_fileio *fio, int *item, int n,
                          const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_int64(t_fileio *fio, gmx_int64_t *item, int n,
                            const char *desc, const char *srcfile,
                            int line);
gmx_bool gmx_fio_ndoe_uchar(t_fileio *fio, unsigned char *item, int n,
                            const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_ushort(t_fileio *fio, unsigned short *item, int n,
                             const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_rvec(t_fileio *fio, rvec *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_ivec(t_fileio *fio, ivec *item, int n,
                           const char *desc, const char *srcfile, int line);
gmx_bool gmx_fio_ndoe_string(t_fileio *fio, char *item[], int n,
                             const char *desc, const char *srcfile, int line);



/* convenience macros */
#define gmx_fio_read_real(fio, item)           gmx_fio_reade_real(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_float(fio, item)          gmx_fio_reade_float(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_double(fio, item)         gmx_fio_reade_double(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_int(fio, item)            gmx_fio_reade_int(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_int64(fio, item)          gmx_fio_reade_int64(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_uchar(fio, item)          gmx_fio_reade_uchar(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_ushort(fio, item)         gmx_fio_reade_ushort(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_rvec(fio, item)           gmx_fio_reade_rvec(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_ivec(fio, item)           gmx_fio_reade_ivec(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_read_string(fio, item)         gmx_fio_reade_string(fio, item, (#item), __FILE__, __LINE__)

#define gmx_fio_write_real(fio, item)           gmx_fio_writee_real(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_float(fio, item)          gmx_fio_writee_float(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_double(fio, item)         gmx_fio_writee_double(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_int(fio, item)            gmx_fio_writee_int(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_int64(fio, item)          gmx_fio_writee_int64(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_uchar(fio, item)          gmx_fio_writee_uchar(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_ushort(fio, item)         gmx_fio_writee_ushort(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_rvec(fio, item)           gmx_fio_writee_rvec(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_ivec(fio, item)           gmx_fio_writee_ivec(fio, item, (#item), __FILE__, __LINE__)
#define gmx_fio_write_string(fio, item)         gmx_fio_writee_string(fio, item, (#item), __FILE__, __LINE__)

#define gmx_fio_do_real(fio, item)              gmx_fio_doe_real(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_float(fio, item)             gmx_fio_doe_float(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_double(fio, item)            gmx_fio_doe_double(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_gmx_bool(fio, item)          gmx_fio_doe_gmx_bool(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_int(fio, item)               gmx_fio_doe_int(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_int64(fio, item)             gmx_fio_doe_int64(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_uchar(fio, item)             gmx_fio_doe_uchar(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_ushort(fio, item)            gmx_fio_doe_ushort(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_rvec(fio, item)              gmx_fio_doe_rvec(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_ivec(fio, item)              gmx_fio_doe_ivec(fio, &item, (#item), __FILE__, __LINE__)
#define gmx_fio_do_string(fio, item)            gmx_fio_doe_string(fio, item, (#item), __FILE__, __LINE__)




#define gmx_fio_nread_real(fio, item, n)            gmx_fio_nreade_real(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_float(fio, item, n)           gmx_fio_nreade_float(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_double(fio, item, n)          gmx_fio_nreade_double(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_int(fio, item, n)             gmx_fio_nreade_int(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_int64(fio, item, n)           gmx_fio_nreade_int64(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_uchar(fio, item, n)           gmx_fio_nreade_uchar(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_ushort(fio, item, n)          gmx_fio_nreade_ushort(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_rvec(fio, item, n)            gmx_fio_nreade_rvec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_ivec(fio, item, n)            gmx_fio_nreade_ivec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nread_string(fio, item, n)          gmx_fio_nreade_string(fio, item, n, (#item), __FILE__, __LINE__)

#define gmx_fio_nwrite_real(fio, item, n)           gmx_fio_nwritee_real(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_float(fio, item, n)          gmx_fio_nwritee_float(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_double(fio, item, n)         gmx_fio_nwritee_double(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_int(fio, item, n)            gmx_fio_nwritee_int(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_int64(fio, item, n)          gmx_fio_nwritee_int64(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_uchar(fio, item, n)          gmx_fio_nwritee_uchar(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_ushort(fio, item, n)         gmx_fio_nwritee_ushort(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_rvec(fio, item, n)           gmx_fio_nwritee_rvec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_ivec(fio, item, n)           gmx_fio_nwritee_ivec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_nwrite_string(fio, item, n)         gmx_fio_nwritee_string(fio, item, n, (#item), __FILE__, __LINE__)

#define gmx_fio_ndo_real(fio, item, n)              gmx_fio_ndoe_real(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_float(fio, item, n)             gmx_fio_ndoe_float(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_double(fio, item, n)            gmx_fio_ndoe_double(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_gmx_bool(fio, item, n)          gmx_fio_ndoe_gmx_bool(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int(fio, item, n)               gmx_fio_ndoe_int(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_int64(fio, item, n)             gmx_fio_ndoe_int64(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_uchar(fio, item, n)             gmx_fio_ndoe_uchar(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ushort(fio, item, n)            gmx_fio_ndoe_ushort(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_rvec(fio, item, n)              gmx_fio_ndoe_rvec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_ivec(fio, item, n)              gmx_fio_ndoe_ivec(fio, item, n, (#item), __FILE__, __LINE__)
#define gmx_fio_ndo_string(fio, item, n)            gmx_fio_ndoe_string(fio, item, n, (#item), __FILE__, __LINE__)

#ifdef __cplusplus
}
#endif

#endif
