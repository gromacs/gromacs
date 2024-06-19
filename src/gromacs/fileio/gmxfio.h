/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_FILEIO_GMXFIO_H
#define GMX_FILEIO_GMXFIO_H

#include <cstdio>

#include <array>
#include <filesystem>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"

typedef struct t_fileio t_fileio;

/* NOTE ABOUT THREAD SAFETY:

   The functions are all thread-safe, provided that two threads don't
   do something silly like closing the same file, or one thread
   accesses a file that has been closed by another.
 */

/********************************************************
 * Open and Close
 ********************************************************/

t_fileio* gmx_fio_open(const std::filesystem::path& fn, const char* mode);
/* Open a new file for reading or writing.
 * The file type will be deduced from the file name.
 */

int gmx_fio_close(t_fileio* fp);
/* Close the file corresponding to fp (if not stdio)
 * The routine will exit when an invalid fio is handled.
 * Returns 0 on success.
 */

int gmx_fio_fp_close(t_fileio* fp);
/* Close the file corresponding to fp without closing the FIO entry
 * Needed e.g. for trxio because the FIO entries are used to store
 * additional data.
 * NOTE that the fp still needs to be properly closed with gmx_fio_close().
 * The routine will exit when an invalid fio is handled.
 * Returns 0 on success.
 */


/* Open a file, return a stream, record the entry in internal FIO object */
FILE* gmx_fio_fopen(const std::filesystem::path& fn, const char* mode);

/* Close a file previously opened with gmx_fio_fopen.
 * Do not mix these calls with standard fopen/fclose ones!
 * Returns 0 on success.  */
int gmx_fio_fclose(FILE* fp);


/********************************************************
 * Change properties of the open file
 ********************************************************/

std::filesystem::path gmx_fio_getname(t_fileio* fio);
/* Return the filename corresponding to the fio index */

int gmx_fio_getftp(t_fileio* fio);
/* Return the filetype corresponding to the fio index.
    There is as of now no corresponding setftp function because the file
    was opened as a specific file type and changing that midway is most
    likely an evil hack. */

gmx_bool gmx_fio_getread(t_fileio* fio);
/* Return  whether read mode is on in fio  */

/***************************************************
 * FILE Operations
 ***************************************************/

void gmx_fio_rewind(t_fileio* fio);
/* Rewind the file in fio */

int gmx_fio_flush(t_fileio* fio);
/* Flush the fio, returns 0 on success */

int gmx_fio_fsync(t_fileio* fio);
/* fsync the fio, returns 0 on success.
   NOTE: don't use fsync function unless you're absolutely sure you need it
   because it deliberately interferes with the OS's caching mechanisms and
   can cause dramatically slowed down IO performance. Some OSes (Linux,
   for example), may implement fsync as a full sync() point. */

gmx_off_t gmx_fio_ftell(t_fileio* fio);
/* Return file position if possible */

int gmx_fio_seek(t_fileio* fio, gmx_off_t fpos);
/* Set file position if possible, quit otherwise */

FILE* gmx_fio_getfp(t_fileio* fio);
/* Return the file pointer itself */


/* Element with information about position in a currently open file.
 * gmx_off_t should be defined by autoconf if your system does not have it.
 * If you do not have it on some other platform you do not have largefile
 * support at all, and you can define it to int (or better, find out how to
 * enable large files).  */
struct gmx_file_position_t
{
    char                          filename[STRLEN] = { 0 };
    gmx_off_t                     offset           = 0;
    std::array<unsigned char, 16> checksum         = { { 0 } };
    int                           checksumSize     = 0;
};

/*! \brief Return data about output files.
 *
 * This is used for handling data stored in the checkpoint files, so
 * we can truncate output files upon restart-with-appending. */
std::vector<gmx_file_position_t> gmx_fio_get_output_file_positions();

t_fileio* gmx_fio_all_output_fsync();
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


int gmx_fio_get_file_md5(t_fileio* fio, gmx_off_t offset, std::array<unsigned char, 16>* checksum);


int xtc_seek_time(t_fileio* fio, real time, int natoms, gmx_bool bSeekForwardOnly);


#endif
