/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "futil.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <tuple>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

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
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/mutex.h"
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

static t_pstack* pstack           = nullptr;
static bool      bUnbuffered      = false;
static int       s_maxBackupCount = 0;

/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
static gmx::Mutex pstack_mutex;

using Lock = gmx::lock_guard<gmx::Mutex>;

namespace gmx
{
namespace
{
//! Global library file finder; stores the object set with setLibraryFileFinder().
const DataFileFinder* g_libFileFinder;
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
    t_pstack* ps;

    Lock pstackLock(pstack_mutex);

    snew(ps, 1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;
}

#if GMX_FAHCORE
/* don't use pipes!*/
#    define popen fah_fopen
#    define pclose fah_fclose
#    define SKIP_FFOPS 1
#else
#    ifdef gmx_ffclose
#        undef gmx_ffclose
#    endif
#    if (!HAVE_PIPES && !defined(__native_client__))
static FILE* popen(const char* nm, const char* mode)
{
    gmx_impl("Sorry no pipes...");

    return NULL;
}

static int pclose(FILE* fp)
{
    gmx_impl("Sorry no pipes...");

    return 0;
}
#    endif /* !HAVE_PIPES && !defined(__native_client__) */
#endif     /* GMX_FAHCORE */

int gmx_ffclose(FILE* fp)
{
#ifdef SKIP_FFOPS
    return fclose(fp);
#else
    t_pstack *ps, *tmp;
    int       ret = 0;

    Lock pstackLock(pstack_mutex);

    ps = pstack;
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
            tmp      = ps->prev;
            ps->prev = ps->prev->prev;
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
#endif
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

int gmx_truncate(const std::string& filename, gmx_off_t length)
{
#if GMX_NATIVE_WINDOWS
    FILE* fp = fopen(filename.c_str(), "rb+");
    if (fp == NULL)
    {
        return -1;
    }
#    ifdef _MSC_VER
    int rc = _chsize_s(fileno(fp), length);
#    else
    int rc = _chsize(fileno(fp), length);
#    endif
    fclose(fp);
    return rc;
#else
    return truncate(filename.c_str(), length);
#endif
}

static FILE* uncompress(const std::string& fn, const char* mode)
{
    FILE*       fp;
    std::string buf = "uncompress -c < " + fn;
    fprintf(stderr, "Going to execute '%s'\n", buf.c_str());
    if ((fp = popen(buf.c_str(), mode)) == nullptr)
    {
        gmx_open(fn);
    }
    push_ps(fp);

    return fp;
}

static FILE* gunzip(const std::string& fn, const char* mode)
{
    FILE*       fp;
    std::string buf = "gunzip -c < ";
    buf += fn;
    fprintf(stderr, "Going to execute '%s'\n", buf.c_str());
    if ((fp = popen(buf.c_str(), mode)) == nullptr)
    {
        gmx_open(fn);
    }
    push_ps(fp);

    return fp;
}

gmx_bool gmx_fexist(const std::string& fname)
{
    FILE* test;

    if (fname.empty())
    {
        return FALSE;
    }
    test = fopen(fname.c_str(), "r");
    if (test == nullptr)
    {
/*Windows doesn't allow fopen of directory - so we need to check this seperately */
#if GMX_NATIVE_WINDOWS
        DWORD attr = GetFileAttributes(fname.c_str());
        return (attr != INVALID_FILE_ATTRIBUTES) && (attr & FILE_ATTRIBUTE_DIRECTORY);
#else
        return FALSE;
#endif
    }
    else
    {
        fclose(test);
        return TRUE;
    }
}

static std::string backup_fn(const std::string& file)
{
    int count = 1;

    std::string directory = gmx::Path::getParentPath(file);
    std::string fn        = gmx::Path::getFilename(file);
    std::string buf;
    if (directory.empty())
    {
        directory = ".";
    }
    do
    {
        buf = gmx::formatString("%s/#%s.%d#", directory.c_str(), fn.c_str(), count);
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
                  s_maxBackupCount, fn.c_str());
    }

    return buf;
}

void make_backup(const std::string& name)
{
    if (s_maxBackupCount <= 0)
    {
        return;
    }
    if (gmx_fexist(name))
    {
        auto backup = backup_fn(name);
        if (rename(name.c_str(), backup.c_str()) == 0)
        {
            fprintf(stderr, "\nBack Off! I just backed up %s to %s\n", name.c_str(), backup.c_str());
        }
        else
        {
            fprintf(stderr, "\nSorry couldn't backup %s to %s\n", name.c_str(), backup.c_str());
        }
    }
}

FILE* gmx_ffopen(const std::string& file, const char* mode)
{
#ifdef SKIP_FFOPS
    return fopen(file, mode);
#else
    FILE*    ff = nullptr;
    gmx_bool bRead;
    int      bs;

    if (file.empty())
    {
        return nullptr;
    }

    if (mode[0] == 'w')
    {
        make_backup(file);
    }

    bRead = (mode[0] == 'r' && mode[1] != '+');
    if (!bRead || gmx_fexist(file))
    {
        if ((ff = fopen(file.c_str(), mode)) == nullptr)
        {
            gmx_file(file);
        }
        /* Check whether we should be using buffering (default) or not
         * (for debugging)
         */
        const char* bufsize = nullptr;
        if (bUnbuffered || ((bufsize = getenv("GMX_LOG_BUFFER")) != nullptr))
        {
            /* Check whether to use completely unbuffered */
            if (bUnbuffered)
            {
                bs = 0;
            }
            else
            {
                bs = strtol(bufsize, nullptr, 10);
            }
            if (bs <= 0)
            {
                setbuf(ff, nullptr);
            }
            else
            {
                char* ptr;
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
        std::string compressedFileName = file;
        compressedFileName += ".Z";
        if (gmx_fexist(compressedFileName))
        {
            ff = uncompress(compressedFileName, mode);
        }
        else
        {
            compressedFileName = file;
            compressedFileName += ".gz";
            if (gmx_fexist(compressedFileName))
            {
                ff = gunzip(compressedFileName, mode);
            }
            else
            {
                gmx_file(file);
            }
        }
    }
    return ff;
#endif
}

namespace gmx
{

std::string findLibraryFile(const std::string& filename, bool bAddCWD, bool bFatal)
{
    std::string result;
    try
    {
        const DataFileFinder& finder = getLibraryFileFinder();
        result                       = finder.findFile(
                DataFileOptions(filename).includeCurrentDir(bAddCWD).throwIfNotFound(bFatal));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    return result;
}

std::string findLibraryFile(const char* filename, bool bAddCWD, bool bFatal)
{
    return findLibraryFile(std::string(filename), bAddCWD, bFatal);
}

FilePtr openLibraryFile(const std::string& filename, bool bAddCWD, bool bFatal)
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

FilePtr openLibraryFile(const char* filename, bool bAddCWD, bool bFatal)
{
    return openLibraryFile(std::string(filename), bAddCWD, bFatal);
}

} // namespace gmx

/*! \brief Use mkstemp (or similar function to make a new temporary
 * file and (on non-Windows systems) return a file descriptor to it.
 *
 * \todo Use std::string and std::vector<char>. */
static int makeTemporaryFilename(char* buf)
{
    int len;

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
    int fd;
#if GMX_NATIVE_WINDOWS
    _mktemp(buf);
    if (buf == NULL)
    {
        gmx_fatal(FARGS, "Error creating temporary file %s: %s", buf, strerror(errno));
    }
    fd = 0;
#else
    fd = mkstemp(buf);

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

int gmx_file_rename(const char* oldname, const char* newname)
{
#if !GMX_NATIVE_WINDOWS
    /* under unix, rename() is atomic (at least, it should be). */
    return rename(oldname, newname);
#else
    if (MoveFileEx(oldname, newname, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
    {
        return 0;
    }
    else
    {
        return 1;
    }
#endif
}

int gmx_file_copy(const char* oldname, const char* newname, gmx_bool copy_if_empty)
{
    gmx::FilePtr in(fopen(oldname, "rb"));
    if (!in)
    {
        return 1;
    }

    /* If we don't copy when empty, we postpone opening the file
       until we're actually ready to write. */
    gmx::FilePtr out;
    if (copy_if_empty)
    {
        out.reset(fopen(newname, "wb"));
        if (!out)
        {
            return 1;
        }
    }

    /* the full copy buffer size: */
    constexpr int     FILECOPY_BUFSIZE = 1 << 16;
    std::vector<char> buf(FILECOPY_BUFSIZE);

    while (!feof(in.get()))
    {
        size_t nread;

        nread = fread(buf.data(), sizeof(char), FILECOPY_BUFSIZE, in.get());
        if (nread > 0)
        {
            size_t ret;
            if (!out)
            {
                /* so this is where we open when copy_if_empty is false:
                   here we know we read something. */
                out.reset(fopen(newname, "wb"));
                if (!out)
                {
                    return 1;
                }
            }
            ret = fwrite(buf.data(), sizeof(char), nread, out.get());
            if (ret != nread)
            {
                return 1;
            }
        }
        if (ferror(in.get()))
        {
            return 1;
        }
    }
    return 0;
}


int gmx_fsync(FILE* fp)
{
    int rc = 0;

#if GMX_FAHCORE
    /* the fahcore defines its own os-independent fsync */
    rc = fah_fsync(fp);
#else /* GMX_FAHCORE */
    {
        int fn;

        /* get the file number */
#    if HAVE_FILENO
        fn = fileno(fp);
#    elif HAVE__FILENO
        fn = _fileno(fp);
#    else
        fn = -1;
#    endif

        /* do the actual fsync */
        if (fn >= 0)
        {
#    if HAVE_FSYNC
            rc = fsync(fn);
#    elif HAVE__COMMIT
            rc = _commit(fn);
#    endif
        }
    }
#endif /* GMX_FAHCORE */

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

void gmx_chdir(const char* directory)
{
#if GMX_NATIVE_WINDOWS
    int rc = _chdir(directory);
#else
    int   rc   = chdir(directory);
#endif
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Cannot change directory to '%s'. Reason: %s", directory, strerror(errno));
    }
}

void gmx_getcwd(char* buffer, size_t size)
{
#if GMX_NATIVE_WINDOWS
    char* pdum = _getcwd(buffer, size);
#else
    char* pdum = getcwd(buffer, size);
#endif
    if (pdum == nullptr)
    {
        gmx_fatal(FARGS, "Cannot get working directory. Reason: %s", strerror(errno));
    }
}
