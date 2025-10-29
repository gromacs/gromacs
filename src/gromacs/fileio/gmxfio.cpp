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
#include "gmxpre.h"

#include "gmxfio.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <mutex>
#include <string>
#include <vector>

#include "gromacs/fileio/xdrf.h"

#if HAVE_IO_H
#    include <io.h>
#endif
#ifdef HAVE_UNISTD_H
#    include <unistd.h>
#endif

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/md5.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "gmxfio_impl.h"

/* This is the new improved and thread safe version of gmxfio. */


/* the list of open files is a linked list, with a dummy element at its head;
       it is initialized when the first file is opened. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static t_fileio* open_files = nullptr;


/* this mutex locks the open_files structure so that list modifications
   are not concurrent with other list operations

   For now, we use this as a coarse grained lock on all file
   insertion/deletion operations and traversing the list,
   because it makes avoiding deadlocks
   easier, and adds almost no overhead: the only overhead is during
   opening and closing of files, or during global operations like
   iterating along all open files. All these cases should be rare
   during the simulation. */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex open_file_mutex;

using Lock = std::lock_guard<std::mutex>;

/******************************************************************
 *
 * Internal functions:
 *
 ******************************************************************/

static int gmx_fio_int_flush(t_fileio* fio)
{
    int rc = 0;

    if (fio->fp)
    {
        rc = std::fflush(fio->fp);
    }

    return rc;
}

/* make a dummy head element, assuming we locked open_file_mutex. */
static void gmx_fio_make_dummy()
{
    if (!open_files)
    {
        open_files     = new t_fileio{};
        open_files->fp = nullptr;
        open_files->fn.clear();
        open_files->next = open_files;
        open_files->prev = open_files;
    }
}


/***********************************************************************
 *
 * FILE LIST OPERATIONS
 *
 ***********************************************************************/


/* insert a new t_fileio into the list */
static void gmx_fio_insert(t_fileio* fio)
{
    t_fileio* prev;
    Lock      openFilesLock(open_file_mutex);
    gmx_fio_make_dummy();

    prev = open_files->prev;

    /* now do the actual insertion: */
    fio->next        = open_files;
    open_files->prev = fio;
    prev->next       = fio;
    fio->prev        = prev;
}

/* remove a t_fileio from the list.
   NOTE: We assume that the open_file_mutex has been locked */
static void gmx_fio_remove(t_fileio* fio)
{
    fio->prev->next = fio->next;
    fio->next->prev = fio->prev;
    /* and make sure we point nowhere in particular */
    fio->next = fio->prev = fio;
}


/* get the first open file, or NULL if there is none.
   Assumes open_file_mutex is locked. */
static t_fileio* gmx_fio_get_first()
{
    t_fileio* ret;

    gmx_fio_make_dummy();

    ret = open_files->next;


    /* check whether there were any to begin with */
    if (ret == open_files)
    {
        /* after this, the open_file pointer should never change */
        ret = nullptr;
    }


    return ret;
}

/* get the next open file, or NULL if there is none.
   Assumes open_file_mutex is locked. */
static t_fileio* gmx_fio_get_next(t_fileio* fio)
{
    t_fileio* ret;

    ret = fio->next;
    /* check if that was the last one */
    if (fio->next == open_files)
    {
        ret = nullptr;
    }

    return ret;
}


/*****************************************************************
 *
 *                     EXPORTED SECTION
 *
 *****************************************************************/
t_fileio* gmx_fio_open(const std::filesystem::path& fn, const char* mode)
{
    t_fileio* fio = nullptr;
    char      newmode[5];
    gmx_bool  bRead, bReadWrite;

    /* sanitize the mode string */
    if (std::strncmp(mode, "r+", 2) == 0)
    {
        std::strcpy(newmode, "r+");
    }
    else if (mode[0] == 'r')
    {
        std::strcpy(newmode, "r");
    }
    else if (std::strncmp(mode, "w+", 2) == 0)
    {
        std::strcpy(newmode, "w+");
    }
    else if (mode[0] == 'w')
    {
        std::strcpy(newmode, "w");
    }
    else if (std::strncmp(mode, "a+", 2) == 0)
    {
        std::strcpy(newmode, "a+");
    }
    else if (mode[0] == 'a')
    {
        std::strcpy(newmode, "a");
    }
    else
    {
        gmx_fatal(FARGS, "DEATH HORROR in gmx_fio_open, mode is '%s'", mode);
    }

    /* Check if it should be opened as a binary file */
    if (!ftp_is_text(fn2ftp(fn)))
    {
        std::strcat(newmode, "b");
    }

    fio        = new t_fileio{};
    bRead      = (newmode[0] == 'r' && newmode[1] != '+');
    bReadWrite = (newmode[1] == '+');
    fio->fp    = nullptr;
    fio->xdr   = nullptr;
    if (!fn.empty())
    {
        if (fn2ftp(fn) == efTNG || fn2ftp(fn) == efH5MD)
        {
            gmx_incons("gmx_fio_open may not be used to open TNG or H5MD files");
        }
        fio->iFTP = fn2ftp(fn);
        fio->fn   = fn;

        fio->fp = gmx_ffopen(fn, newmode);
        /* If this file type is in the list of XDR files, open it like that */
        if (ftp_is_xdr(fio->iFTP))
        {
            /* determine the XDR direction */
            if (newmode[0] == 'w' || newmode[0] == 'a')
            {
                fio->xdrmode = XDR_ENCODE;
            }
            else
            {
                fio->xdrmode = XDR_DECODE;
            }
            snew(fio->xdr, 1);
            xdrstdio_create(fio->xdr, fio->fp, fio->xdrmode);
        }

        /* for appending seek to end of file to make sure ftell gives correct position
         * important for checkpointing */
        if (newmode[0] == 'a')
        {
            gmx_fseek(fio->fp, 0, SEEK_END);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Cannot open file with empty filename");
    }

    fio->bRead      = bRead;
    fio->bReadWrite = bReadWrite;
    fio->bDouble    = (sizeof(real) == sizeof(double));

    /* and now insert this file into the list of open files. */
    gmx_fio_insert(fio);
    return fio;
}

static int gmx_fio_close_inner(t_fileio* fio)
{
    int rc = 0;

    if (fio->xdr != nullptr)
    {
        xdr_destroy(fio->xdr);
        sfree(fio->xdr);
    }

    if (fio->fp != nullptr)
    {
        rc = gmx_ffclose(fio->fp); /* fclose returns 0 if happy */
    }

    return rc;
}

int gmx_fio_close(t_fileio* fio)
{
    int rc = 0;

    Lock openFilesLock(open_file_mutex);
    /* first remove it from the list */
    gmx_fio_remove(fio);
    rc = gmx_fio_close_inner(fio);
    delete fio;

    return rc;
}

/* close only fp but keep FIO entry. */
int gmx_fio_fp_close(t_fileio* fio)
{
    int rc = 0;
    if (fio->xdr == nullptr)
    {
        rc      = gmx_ffclose(fio->fp); /* fclose returns 0 if happy */
        fio->fp = nullptr;
    }

    return rc;
}

FILE* gmx_fio_fopen(const std::filesystem::path& fn, const char* mode)
{
    FILE*     ret;
    t_fileio* fio;

    fio = gmx_fio_open(fn, mode);
    ret = fio->fp;

    return ret;
}

int gmx_fio_fclose(FILE* fp)
{
    t_fileio* cur;
    int       rc = -1;

    Lock openFilesLock(open_file_mutex);
    cur = gmx_fio_get_first();
    while (cur)
    {
        if (cur->fp == fp)
        {
            rc = gmx_fio_close_inner(cur);
            gmx_fio_remove(cur);
            delete cur;
            break;
        }
        cur = gmx_fio_get_next(cur);
    }

    return rc;
}

/*! \brief Internal variant of get_file_md5
 *
 * \return -1 any time a checksum cannot be computed, otherwise the
 *            length of the data from which the checksum was computed. */
static int gmx_fio_int_get_file_md5(t_fileio* fio, gmx_off_t offset, std::array<unsigned char, 16>* checksum)
{
    /*1MB: large size important to catch almost identical files */
    constexpr size_t maximumChecksumInputSize = 1048576;
    md5_state_t      state;
    gmx_off_t        readLength;
    gmx_off_t        seekOffset;

    seekOffset = offset - maximumChecksumInputSize;
    if (seekOffset < 0)
    {
        seekOffset = 0;
    }
    readLength = offset - seekOffset;

    if (!fio->fp)
    {
        // It's not an error if the file isn't open.
        return -1;
    }
    if (!fio->bReadWrite)
    {
        // It's not an error if the file is open in the wrong mode.
        //
        // TODO It is unclear why this check exists. The bReadWrite
        // flag is true when the file-opening mode included "+" but we
        // only need read and seek to be able to compute the
        // md5sum. Other requirements (e.g. that we can truncate when
        // doing an appending restart) should be expressed in a
        // different way, but it is unclear whether that is part of
        // the logic here.
        return -1;
    }

    if (gmx_fseek(fio->fp, seekOffset, SEEK_SET))
    {
        // It's not an error if file seeking fails. (But it could be
        // an issue when moving a checkpoint from one platform to
        // another, when they differ in their support for seeking, and
        // so can't agree on a checksum for appending).
        gmx_fseek(fio->fp, 0, SEEK_END);
        return -1;
    }

    std::vector<unsigned char> buf(maximumChecksumInputSize);
    // The fread puts the file position back to offset.
    if (static_cast<gmx_off_t>(std::fread(buf.data(), 1, readLength, fio->fp)) != readLength)
    {
        // Read an unexpected length. This is not a fatal error; the
        // md5sum check to prevent overwriting files is not vital.
        if (std::ferror(fio->fp))
        {
            fprintf(stderr, "\nTrying to get md5sum: %s: %s\n", fio->fn.string().c_str(), std::strerror(errno));
        }
        else if (!std::feof(fio->fp))
        {
            fprintf(stderr,
                    "\nTrying to get md5sum: Unknown reason for short read: %s\n",
                    fio->fn.string().c_str());
        }

        gmx_fseek(fio->fp, 0, SEEK_END);
        return -1;
    }
    // Return the file position to the end of the file.
    gmx_fseek(fio->fp, 0, SEEK_END);

    if (debug)
    {
        fprintf(debug, "chksum %s readlen %ld\n", fio->fn.string().c_str(), static_cast<long int>(readLength));
    }

    gmx_md5_init(&state);
    gmx_md5_append(&state, buf.data(), readLength);
    *checksum = gmx_md5_finish(&state);
    return readLength;
}


/*
 * fio: file to compute md5 for
 * offset: starting pointer of region to use for md5
 * digest: return array of md5 sum
 */
int gmx_fio_get_file_md5(t_fileio* fio, gmx_off_t offset, std::array<unsigned char, 16>* checksum)
{
    int ret;

    ret = gmx_fio_int_get_file_md5(fio, offset, checksum);

    return ret;
}

static int gmx_fio_int_get_file_position(t_fileio* fio, gmx_off_t* offset)
{
    /* Flush the file, so we are sure it is written */
    if (gmx_fio_int_flush(fio))
    {
        char buf[STRLEN];
        sprintf(buf, "Cannot write file '%s'; maybe you are out of disk space?", fio->fn.string().c_str());
        gmx_file(buf);
    }

    /* We cannot count on XDR being able to write 64-bit integers,
       so separate into high/low 32-bit values.
       In case the filesystem has 128-bit offsets we only care
       about the first 64 bits - we'll have to fix
       this when exabyte-size output files are common...
     */
    *offset = gmx_ftell(fio->fp);

    return 0;
}

std::vector<gmx_file_position_t> gmx_fio_get_output_file_positions()
{
    std::vector<gmx_file_position_t> outputfiles;
    t_fileio*                        cur;

    Lock openFilesLock(open_file_mutex);
    cur = gmx_fio_get_first();
    while (cur)
    {
        /* Skip the checkpoint files themselves, since they could be open when
           we call this routine... */
        if (!cur->bRead && cur->iFTP != efCPT)
        {
            outputfiles.emplace_back();
            std::strncpy(outputfiles.back().filename, cur->fn.string().data(), STRLEN - 1);

            /* Get the file position */
            gmx_fio_int_get_file_position(cur, &outputfiles.back().offset);
            if (!GMX_FAHCORE)
            {
                outputfiles.back().checksumSize = gmx_fio_int_get_file_md5(
                        cur, outputfiles.back().offset, &outputfiles.back().checksum);
            }
        }

        cur = gmx_fio_get_next(cur);
    }

    return outputfiles;
}


std::filesystem::path gmx_fio_getname(t_fileio* fio)
{
    std::filesystem::path ret;
    ret = fio->fn;

    return ret;
}

int gmx_fio_getftp(t_fileio* fio)
{
    int ret;

    ret = fio->iFTP;

    return ret;
}

void gmx_fio_rewind(t_fileio* fio)
{
    if (fio->xdr)
    {
        xdr_destroy(fio->xdr);
        frewind(fio->fp);
        xdrstdio_create(fio->xdr, fio->fp, fio->xdrmode);
    }
    else
    {
        frewind(fio->fp);
    }
}


int gmx_fio_flush(t_fileio* fio)
{
    int ret;

    ret = gmx_fio_int_flush(fio);

    return ret;
}


/* fsync the fio, returns 0 on success.
   NOTE: don't use fsync function unless you're absolutely sure you need it
   because it deliberately interferes with the OS's caching mechanisms and
   can cause dramatically slowed down IO performance. Some OSes (Linux,
   for example), may implement fsync as a full sync() point. */
static int gmx_fio_int_fsync(t_fileio* fio)
{
    int rc = 0;

    if (fio->fp)
    {
        rc = gmx_fsync(fio->fp);
    }
    return rc;
}

t_fileio* gmx_fio_all_output_fsync()
{
    t_fileio* ret = nullptr;
    t_fileio* cur;

    Lock openFilesLock(open_file_mutex);
    cur = gmx_fio_get_first();
    while (cur)
    {
        if (!cur->bRead)
        {
            /* if any of them fails, return failure code */
            int rc = gmx_fio_int_fsync(cur);
            if (rc != 0 && !ret)
            {
                ret = cur;
            }
        }
        cur = gmx_fio_get_next(cur);
    }

    /* in addition, we force these to be written out too, if they're being
       redirected. We don't check for errors because errors most likely mean
       that they're not redirected. */
    std::fflush(stdout);
    std::fflush(stderr);
#if HAVE_FSYNC
    /* again, fahcore defines HAVE_FSYNC and fsync() */
    fsync(STDOUT_FILENO);
    fsync(STDERR_FILENO);
#endif

    return ret;
}


gmx_off_t gmx_fio_ftell(t_fileio* fio)
{
    gmx_off_t ret = 0;

    if (fio->fp)
    {
        ret = gmx_ftell(fio->fp);
    }
    return ret;
}

int gmx_fio_seek(t_fileio* fio, gmx_off_t fpos)
{
    int rc;

    if (fio->fp)
    {
        rc = gmx_fseek(fio->fp, fpos, SEEK_SET);
    }
    else
    {
        gmx_file(fio->fn.string());
    }
    return rc;
}

FILE* gmx_fio_getfp(t_fileio* fio)
{
    FILE* ret = nullptr;

    if (fio->fp)
    {
        ret = fio->fp;
    }
    return ret;
}

gmx_bool gmx_fio_getread(t_fileio* fio)
{
    gmx_bool ret;

    ret = fio->bRead;

    return ret;
}

int xtc_seek_time(t_fileio* fio, real time, int natoms, gmx_bool bSeekForwardOnly)
{
    int ret;

    ret = xdr_xtc_seek_time(time, fio->fp, fio->xdr, natoms, bSeekForwardOnly);

    return ret;
}
