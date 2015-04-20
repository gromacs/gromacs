/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Low-level wrappers for OS-specific file handling with some \Gromacs
 * customizations.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FUTIL_H
#define GMX_UTILITY_FUTIL_H

#include <limits.h>
#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/*! \def GMX_PATH_MAX
 * \brief
 * Maximum path length, if the OS provides one, otherwise a fixed constant.
 */
#ifdef PATH_MAX
#  define GMX_PATH_MAX PATH_MAX
#elif defined MAX_PATH
#  define GMX_PATH_MAX MAX_PATH
#else
#  define GMX_PATH_MAX 4096
#endif

/** \Gromacs definition to use instead of `off_t`. */
typedef gmx_int64_t    gmx_off_t;

/*! \brief
 * Turn off buffering for output files (which is default) for debugging
 * purposes.
 *
 * This only has effect on files opened with gmx_ffopen().
 */
void gmx_disable_file_buffering(void);

/*! \brief
 * Enables backups with the specified number of maximum backups.
 *
 * If \p count == 0, disables backups.  If not called, this is the default.
 * If \p count == -1, reads the count from an environment variable.
 *
 * This is not currently thread-safe, as it is only called during
 * initialization code.
 */
void gmx_set_max_backup_count(int count);

/*! \brief
 * Check whether a path exists.
 *
 * \returns `TRUE` when \p fname exists.
 *
 * Note that this returns `TRUE` even if \p fname is a directory instead of a
 * file.
 */
gmx_bool gmx_fexist(const char *fname);

/*! \brief
 * Makes a backup of file if the file exists.
 */
void make_backup(const char *file);

/*! \brief
 * Opens a file, with \Gromacs-specific additions.
 *
 * If the file is in compressed format, opens a pipe which uncompresses the
 * file on the fly.  For this to work, gmx_ffclose() and frewind() should
 * always be used for files opened with gmx_ffopen() instead of fclose() and
 * rewind().  For compressed files, the \p file parameter should be passed
 * without the compressed extension; if that file is not found, then a few
 * compression extensions are tried.
 * Creates a backup if a file opened for writing already exists before
 * overwriting it.
 * A fatal error results if the file cannot be opened, for whatever reason.
 */
FILE *gmx_ffopen(const char *file, const char *mode);

/** Closes a file opened with gmx_ffopen(). */
int gmx_ffclose(FILE *fp);

/*! \brief
 * Wraps rewind() for files opened with gmx_ffopen().
 *
 * A fatal error results if this function is called for a pipe (a compressed
 * input file).
 */
void frewind(FILE *fp);

/** OS-independent 64-bit fseek(). */
int gmx_fseek(FILE *stream, gmx_off_t offset, int whence);

/** OS-independent 64-bit ftell(). */
gmx_off_t gmx_ftell(FILE *stream);

/** OS-independent truncate(). */
int gmx_truncate(const char *filename, gmx_off_t length);

/*! \brief
 * Finds full path for a library file.
 *
 * Searches first in the current directory, and then in the configured library
 * directories.
 * Fatal error results if the file is not found in any location.
 * The caller is responsible of freeing the returned string.
 */
char *gmxlibfn(const char *file);

/*! \brief
 * Opens a library file for reading.
 *
 * Works as gmxlibfn(), except that it opens the file and returns a file
 * handle.
 */
FILE *libopen(const char *file);

/*! \brief
 * More flexible gmxlibfn().
 *
 * Works as gmxlibfn(), but provides control whether the current working
 * directory is searched or not, and whether a missing file is a fatal error or
 * not.
 */
char *low_gmxlibfn(const char *file, gmx_bool bAddCWD, gmx_bool bFatal);

/*! \brief
 * Alternative for libopen() that optionally does not exit.
 *
 * Works as libopen(), but provides control whether a missing file is a fatal
 * error or not.
 */
FILE *low_libopen(const char *file, gmx_bool bFatal);


/*! \brief
 * Creates unique name for temp file (wrapper around mkstemp).
 *
 * \p buf should be at least 7 bytes long
 */
void gmx_tmpnam(char *buf);

/*! \brief
 * OS-independent rename().
 *
 * Renames/moves a file atomically, if the OS makes that available.
 */
int gmx_file_rename(const char *oldname, const char *newname);

/*! \brief
 * Copies a file (data only) oldname to newname.
 *
 * If \p copy_if_empty is `FALSE`, the file won't be copied if it's empty.
 */
int gmx_file_copy(const char *oldname, const char *newname, gmx_bool copy_if_empty);

/*! \brief
 * OS-independent fsync().
 *
 * Only use this during checkpointing!
 */
int gmx_fsync(FILE *fp);

/*! \brief
 * OS-independent chdir().
 *
 * Exits with a fatal error if changing the directory fails.
 */
void gmx_chdir(const char *directory);
/*! \brief
 * OS-independent getcwd().
 *
 * Exits with a fatal error if the call fails.
 */
void gmx_getcwd(char *buffer, size_t size);

#ifdef __cplusplus
}

namespace gmx
{

class DataFileFinder;

/*! \brief
 * Gets a finder for locating data files from share/top/.
 *
 * \returns Finder set with setLibraryFileFinder(), or a default finder.
 *
 * If setLibraryFileFinder() has not been called (or a `NULL` finder has been
 * set), a default finder is returned.
 * The default finder searches data files from the directory identified by the
 * global program context; it does not respect GMXLIB environment variable.
 * Calling initForCommandLine() sets a finder that respects GMXLIB.
 *
 * Does not throw.
 *
 * See setLibraryFileFinder() for thread safety.
 *
 * \ingroup module_utility
 */
const DataFileFinder &getLibraryFileFinder();
/*! \brief
 * Sets a finder for location data files from share/top/.
 *
 * \param[in] finder  finder to set
 *     (can be NULL to restore the default finder).
 *
 * The library does not take ownership of \p finder.
 * The provided object must remain valid until the global instance is changed
 * by another call to setLibraryFileFinder().
 *
 * The global instance is used by gmxlibfn() and libopen().
 *
 * This method is not thread-safe.  See setProgramContext(); the same
 * constraints apply here as well.
 *
 * Does not throw.
 */
void setLibraryFileFinder(const DataFileFinder *finder);

} // namespace gmx
#endif

#endif
