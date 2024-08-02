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

#include "gromacs/utility/futil.h"

#include "config.h"

#include <fcntl.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <mutex>
#include <string>
#include <system_error>
#include <tuple>

#include <sys/stat.h>
#include <sys/types.h>

#include "gromacs/utility/fileptr.h"
#include "gromacs/utility/unique_cptr.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h>
#endif
#if GMX_NATIVE_WINDOWS
#    include <direct.h> // For _chdir() and _getcwd()
#    include <io.h>
#    include <windows.h>
#endif

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/datafilefinder.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* we keep a linked list of all files opened through pipes (i.e.
   compressed or .gzipped files. This way we can distinguish between them
   without having to change the semantics of reading from/writing to files)
 */
typedef struct t_pstack
{
    FILE*            fp;
    struct t_pstack* prev;
} t_pstack;

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static t_pstack* pstack = nullptr;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static bool bUnbuffered = false;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static int s_maxBackupCount = 0;

/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex pstack_mutex;

using Lock = std::lock_guard<std::mutex>;

namespace gmx
{
namespace
{
//! Global library file finder; stores the object set with setLibraryFileFinder().
const DataFileFinder* g_libFileFinder; //NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
//! Default library file finder if nothing is set.
const DataFileFinder g_defaultLibFileFinder;
} // namespace

const DataFileFinder& getLibraryFileFinder()
{
    if (g_libFileFinder != nullptr)
    {
        return *g_libFileFinder;
    }
    return g_defaultLibFileFinder;
}

void setLibraryFileFinder(const DataFileFinder* finder)
{
    g_libFileFinder = finder;
}

} // namespace gmx

void gmx_disable_file_buffering()
{
    bUnbuffered = true;
}

void gmx_set_max_backup_count(int count)
{
    if (count < 0)
    {
        const char* env = getenv("GMX_MAXBACKUP");
        if (env != nullptr)
        {
            // TODO: Check that the value is converted properly.
            count = strtol(env, nullptr, 10);
            if (count < 0)
            {
                count = 0;
            }
        }
        else
        {
            // Use a reasonably low value for countmax; we might
            // generate 4-5 files in each round, and we don't
            // want to hit directory limits of 1024 or 2048 files.
            count = 99;
        }
    }
    s_maxBackupCount = count;
}

static void push_ps(FILE* fp)
{
    t_pstack* ps = nullptr;

    Lock pstackLock(pstack_mutex);

    snew(ps, 1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;
}

#if GMX_FAHCORE
#    ifdef gmx_ffclose
#        undef gmx_ffclose
#    endif
#endif
#if (!HAVE_PIPES && !defined(__native_client__))
static FILE* popen(const char* /* nm */, const char* /* mode */)
{
    gmx_impl("Sorry no pipes...");

    return NULL;
}

static int pclose(FILE* /* fp */)
{
    gmx_impl("Sorry no pipes...");

    return 0;
}
#endif /* !HAVE_PIPES && !defined(__native_client__) */

int gmx_ffclose(FILE* fp)
{
    int ret = 0;

    Lock pstackLock(pstack_mutex);

    t_pstack* ps = pstack;
    if (ps == nullptr)
    {
        if (fp != nullptr)
        {
            ret = fclose(fp);
        }
    }
    else if (ps->fp == fp)
    {
        if (fp != nullptr)
        {
            ret = pclose(fp);
        }
        pstack = pstack->prev;
        sfree(ps);
    }
    else
    {
        while ((ps->prev != nullptr) && (ps->prev->fp != fp))
        {
            ps = ps->prev;
        }
        if ((ps->prev != nullptr) && ps->prev->fp == fp)
        {
            if (ps->prev->fp != nullptr)
            {
                ret = pclose(ps->prev->fp);
            }
            t_pstack* tmp = ps->prev;
            ps->prev      = ps->prev->prev;
            sfree(tmp);
        }
        else
        {
            if (fp != nullptr)
            {
                ret = fclose(fp);
            }
        }
    }

    return ret;
}


void frewind(FILE* fp)
{
    Lock pstackLock(pstack_mutex);

    t_pstack* ps = pstack;
    while (ps != nullptr)
    {
        if (ps->fp == fp)
        {
            fprintf(stderr, "Cannot rewind compressed file!\n");
            return;
        }
        ps = ps->prev;
    }
    rewind(fp);
}

int gmx_fseek(FILE* stream, gmx_off_t offset, int whence)
{
#if HAVE_FSEEKO
    return fseeko(stream, offset, whence);
#else
#    if HAVE__FSEEKI64
    return _fseeki64(stream, offset, whence);
#    else
    return fseek(stream, offset, whence);
#    endif
#endif
}

gmx_off_t gmx_ftell(FILE* stream)
{
#if HAVE_FSEEKO
    return ftello(stream);
#else
#    if HAVE__FSEEKI64
#        ifndef __MINGW32__
    return _ftelli64(stream);
#        else
    return ftello64(stream);
#        endif
#    else
    return ftell(stream);
#    endif
#endif
}

int gmx_truncate(const std::filesystem::path& filename, gmx_off_t length)
{
    std::error_code errorCode;
    std::filesystem::resize_file(filename, length, errorCode);
    return errorCode.value();
}

static FILE* uncompress(const std::filesystem::path& fn, const char* mode)
{
    FILE*       fp  = nullptr;
    std::string buf = "uncompress -c < " + fn.string();
    fprintf(stderr, "Going to execute '%s'\n", buf.c_str());
    if ((fp = popen(buf.c_str(), mode)) == nullptr)
    {
        gmx_open(fn.string());
    }
    push_ps(fp);

    return fp;
}

static FILE* gunzip(const std::filesystem::path& fn, const char* mode)
{
    FILE*       fp  = nullptr;
    std::string buf = "gunzip -c < ";
    buf += fn.string();
    fprintf(stderr, "Going to execute '%s'\n", buf.c_str());
    if ((fp = popen(buf.c_str(), mode)) == nullptr)
    {
        gmx_open(fn.string());
    }
    push_ps(fp);

    return fp;
}

bool gmx_fexist(const std::filesystem::path& fname)
{
    if (fname.empty())
    {
        return false;
    }
    std::error_code errorCode;
    return std::filesystem::exists(fname, errorCode);
}

static std::filesystem::path backup_fn(const std::filesystem::path& file)
{
    int count = 1;

    auto        directory = file.parent_path();
    auto        fn        = file.filename();
    std::string buf;
    if (directory.empty())
    {
        directory = ".";
    }
    do
    {
        buf = gmx::formatString("%s/#%s.%d#", directory.string().c_str(), fn.string().c_str(), count);
        count++;
    } while ((count <= s_maxBackupCount) && gmx_fexist(buf));

    /* Arbitrarily bail out */
    if (count > s_maxBackupCount)
    {
        /* TODO: The error message is only accurate for code that starts with
         * Gromacs command-line interface. */
        gmx_fatal(FARGS,
                  "Won't make more than %d backups of %s for you.\n"
                  "The env.var. GMX_MAXBACKUP controls this maximum, -1 disables backups.",
                  s_maxBackupCount,
                  fn.string().c_str());
    }

    return buf;
}

void make_backup(const std::filesystem::path& name)
{
    if (s_maxBackupCount <= 0)
    {
        return;
    }
    if (gmx_fexist(name))
    {
        auto            backup = backup_fn(name);
        std::error_code errorCode;
        std::filesystem::rename(name, backup, errorCode);
        if (errorCode.value() == 0)
        {
            fprintf(stderr,
                    "\nBack Off! I just backed up %s to %s\n",
                    name.string().c_str(),
                    backup.string().c_str());
        }
        else
        {
            fprintf(stderr,
                    "\nSorry couldn't backup %s to %s\n",
                    name.string().c_str(),
                    backup.string().c_str());
        }
    }
}

FILE* gmx_ffopen(const std::filesystem::path& file, const char* mode)
{
    FILE* ff = nullptr;

    if (file.empty())
    {
        return nullptr;
    }

    if (mode[0] == 'w')
    {
        make_backup(file);
    }

    bool bRead = (mode[0] == 'r' && mode[1] != '+');
    if (!bRead || gmx_fexist(file))
    {
        if ((ff = fopen(file.string().c_str(), mode)) == nullptr)
        {
            gmx_file(file.string());
        }
        /* Check whether we should be using buffering (default) or not
         * (for debugging)
         */
        const char* bufsize = nullptr;
        if (bUnbuffered || ((bufsize = getenv("GMX_LOG_BUFFER")) != nullptr))
        {
            /* Check whether to use completely unbuffered */
            const int bs = bUnbuffered ? 0 : strtol(bufsize, nullptr, 10);
            if (bs <= 0)
            {
                setbuf(ff, nullptr);
            }
            else
            {
                // Note: this leaks memory, because one has to free ptr after closing the file.
                char* ptr = nullptr;
                snew(ptr, bs + 8);
                if (setvbuf(ff, ptr, _IOFBF, bs) != 0)
                {
                    gmx_file("Buffering File");
                }
            }
        }
    }
    else
    {
        auto compressedFileName = file;
        compressedFileName.concat(".Z");
        if (gmx_fexist(compressedFileName))
        {
            ff = uncompress(compressedFileName, mode);
        }
        else
        {
            compressedFileName = file;
            compressedFileName.concat(".gz");
            if (gmx_fexist(compressedFileName))
            {
                ff = gunzip(compressedFileName, mode);
            }
            else
            {
                gmx_file(file.string());
            }
        }
    }
    return ff;
}

namespace gmx
{

std::filesystem::path findLibraryFile(const std::filesystem::path& filename, bool bAddCWD, bool bFatal)
{
    std::filesystem::path result;
    try
    {
        const DataFileFinder& finder = getLibraryFileFinder();
        result                       = finder.findFile(
                DataFileOptions(filename).includeCurrentDir(bAddCWD).throwIfNotFound(bFatal));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    return result;
}

FilePtr openLibraryFile(const std::filesystem::path& filename, bool bAddCWD, bool bFatal)
{
    FilePtr fp;
    try
    {
        const DataFileFinder& finder = getLibraryFileFinder();
        fp = finder.openFile(DataFileOptions(filename).includeCurrentDir(bAddCWD).throwIfNotFound(bFatal));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    return fp;
}

} // namespace gmx

/*! \brief Use mkstemp (or similar function to make a new temporary
 * file and (on non-Windows systems) return a file descriptor to it.
 *
 * Note: not thread-safe on non-Windows systems
 *
 * \todo Use std::string and std::vector<char>. */
static int makeTemporaryFilename(char* buf)
{
    int len = 0;

    if ((len = strlen(buf)) < 7)
    {
        gmx_fatal(FARGS, "Buf passed to gmx_tmpnam must be at least 7 bytes long");
    }
    for (int i = len - 6; (i < len); i++)
    {
        buf[i] = 'X';
    }
    /* mktemp is dangerous and we should use mkstemp instead, but
     * since windows doesnt support it we have to separate the cases.
     * 20090307: mktemp deprecated, use iso c++ _mktemp instead.
     */
#if GMX_NATIVE_WINDOWS
    _mktemp(buf);
    if (buf == NULL)
    {
        gmx_fatal(FARGS, "Error creating temporary file %s: %s", buf, strerror(errno));
    }
    int fd = 0;
#else
    int fd = mkstemp(buf);

    /* mkstemp creates 0600 files - respect umask instead */
    mode_t currUmask = umask(0);
    umask(currUmask);
    fchmod(fd, 0666 & ~currUmask);

    if (fd < 0)
    {
        gmx_fatal(FARGS, "Error creating temporary file %s: %s", buf, strerror(errno));
    }
#endif
    return fd;
}
// TODO use std::string
void gmx_tmpnam(char* buf)
{
    int fd = makeTemporaryFilename(buf);
#if !GMX_NATIVE_WINDOWS
    close(fd);
#endif
}

// TODO use std::string
FILE* gmx_fopen_temporary(char* buf)
{
    FILE* fpout = nullptr;
    int   fd    = makeTemporaryFilename(buf);

#if GMX_NATIVE_WINDOWS
    if ((fpout = fopen(buf, "w")) == NULL)
    {
        gmx_fatal(FARGS, "Cannot open temporary file %s", buf);
    }
#else
    if ((fpout = fdopen(fd, "w")) == nullptr)
    {
        gmx_fatal(FARGS, "Cannot open temporary file %s", buf);
    }
#endif

    return fpout;
}

void gmx_file_rename(const std::filesystem::path& oldname, const std::filesystem::path& newname)
{
    std::error_code errorCode;
    std::filesystem::rename(oldname, newname, errorCode);
#if GMX_FAHCORE
    /* This just lets the F@H checksumming system know about the rename */
    if (errorCode.value() == 0)
    {
        fcRename(oldname.string().c_str(), newname.string().c_str());
    }
#endif
    if (errorCode.value() != 0)
    {
        auto errorMsg = gmx::formatString(
                "Failed to rename %s to %s.", oldname.string().c_str(), newname.string().c_str());
        GMX_THROW(gmx::FileIOError(errorMsg));
    }
}

int gmx_file_copy(const std::filesystem::path& oldname, const std::filesystem::path& newname, bool copy_if_empty)
{
    if (!std::filesystem::exists(oldname))
    {
        return 1;
    }

    std::error_code errorCode;
    if (!std::filesystem::is_empty(oldname) || copy_if_empty)
    {
        std::filesystem::copy_file(
                oldname, newname, std::filesystem::copy_options::overwrite_existing, errorCode);
        return errorCode.value();
    }
    return 0;
}


int gmx_fsync(FILE* fp)
{
    int rc = 0;

    {
        /* get the file number */
#if HAVE_FILENO
        int fn = fileno(fp);
#elif HAVE__FILENO
        int fn = _fileno(fp);
#else
        GMX_UNUSED_VALUE(fp);
        int fn = -1;
#endif

        /* do the actual fsync */
        if (fn >= 0)
        {
#if HAVE_FSYNC
            rc = fsync(fn);
#elif HAVE__COMMIT
            rc = _commit(fn);
#endif
        }
    }

    /* We check for these error codes this way because POSIX requires them
       to be defined, and using anything other than macros is unlikely: */
#ifdef EINTR
    /* we don't want to report an error just because fsync() caught a signal.
       For our purposes, we can just ignore this. */
    if (rc && errno == EINTR)
    {
        rc = 0;
    }
#endif
#ifdef EINVAL
    /* we don't want to report an error just because we tried to fsync()
       stdout, a socket or a pipe. */
    if (rc && errno == EINVAL)
    {
        rc = 0;
    }
#endif
    return rc;
}

void gmx_chdir(const std::filesystem::path& directory)
{
    std::error_code errorCode;
    std::filesystem::current_path(directory, errorCode);
    if (errorCode.value() != 0)
    {
        auto message = gmx::formatString("Cannot change directory to '%s'. Reason: %s",
                                         directory.string().c_str(),
                                         errorCode.message().c_str());
        GMX_THROW(gmx::FileIOError(message));
    }
}

std::filesystem::path gmx_getcwd()
{
    return std::filesystem::current_path();
}
