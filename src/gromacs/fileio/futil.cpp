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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef HAVE_DIRENT_H
/* POSIX */
#include <dirent.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef GMX_NATIVE_WINDOWS
#include <windows.h>
#include <direct.h>
#include <io.h>
#endif

/* Windows file stuff, only necessary for visual studio */
#ifdef _MSC_VER
#include <windows.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/network.h"

#include "gromacs/fileio/futil.h"
#include "gromacs/fileio/path.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/* we keep a linked list of all files opened through pipes (i.e.
   compressed or .gzipped files. This way we can distinguish between them
   without having to change the semantics of reading from/writing to files)
 */
typedef struct t_pstack {
    FILE            *fp;
    struct t_pstack *prev;
} t_pstack;

static t_pstack    *pstack      = NULL;
static gmx_bool     bUnbuffered = FALSE;

/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
static tMPI_Thread_mutex_t pstack_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

void no_buffers(void)
{
    bUnbuffered = TRUE;
}

void push_ps(FILE *fp)
{
    t_pstack *ps;

    tMPI_Thread_mutex_lock(&pstack_mutex);

    snew(ps, 1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;

    tMPI_Thread_mutex_unlock(&pstack_mutex);
}

#ifdef GMX_FAHCORE
/* don't use pipes!*/
#define popen fah_fopen
#define pclose fah_fclose
#define SKIP_FFOPS 1
#else
#ifdef gmx_ffclose
#undef gmx_ffclose
#endif
#if (!defined(HAVE_PIPES) && !defined(__native_client__))
static FILE *popen(const char *nm, const char *mode)
{
    gmx_impl("Sorry no pipes...");

    return NULL;
}

static int pclose(FILE *fp)
{
    gmx_impl("Sorry no pipes...");

    return 0;
}
#endif /* !defined(HAVE_PIPES) && !defined(__native_client__) */
#endif /* GMX_FAHCORE */

int gmx_ffclose(FILE *fp)
{
#ifdef SKIP_FFOPS
    return fclose(fp);
#else
    t_pstack *ps, *tmp;
    int       ret = 0;

    tMPI_Thread_mutex_lock(&pstack_mutex);

    ps = pstack;
    if (ps == NULL)
    {
        if (fp != NULL)
        {
            ret = fclose(fp);
        }
    }
    else if (ps->fp == fp)
    {
        if (fp != NULL)
        {
            ret = pclose(fp);
        }
        pstack = pstack->prev;
        sfree(ps);
    }
    else
    {
        while ((ps->prev != NULL) && (ps->prev->fp != fp))
        {
            ps = ps->prev;
        }
        if ((ps->prev != NULL) && ps->prev->fp == fp)
        {
            if (ps->prev->fp != NULL)
            {
                ret = pclose(ps->prev->fp);
            }
            tmp      = ps->prev;
            ps->prev = ps->prev->prev;
            sfree(tmp);
        }
        else
        {
            if (fp != NULL)
            {
                ret = fclose(fp);
            }
        }
    }

    tMPI_Thread_mutex_unlock(&pstack_mutex);
    return ret;
#endif
}


#ifdef rewind
#undef rewind
#endif

void frewind(FILE *fp)
{
    tMPI_Thread_mutex_lock(&pstack_mutex);

    t_pstack *ps = pstack;
    while (ps != NULL)
    {
        if (ps->fp == fp)
        {
            fprintf(stderr, "Cannot rewind compressed file!\n");
            tMPI_Thread_mutex_unlock(&pstack_mutex);
            return;
        }
        ps = ps->prev;
    }
    rewind(fp);
    tMPI_Thread_mutex_unlock(&pstack_mutex);
}

int gmx_fseek(FILE *stream, gmx_off_t offset, int whence)
{
#ifdef HAVE_FSEEKO
    return fseeko(stream, offset, whence);
#else
#ifdef HAVE__FSEEKI64
    return _fseeki64(stream, offset, whence);
#else
    return fseek(stream, offset, whence);
#endif
#endif
}

gmx_off_t gmx_ftell(FILE *stream)
{
#ifdef HAVE_FSEEKO
    return ftello(stream);
#else
#ifdef HAVE__FSEEKI64
#ifndef __MINGW32__
    return _ftelli64(stream);
#else
    return ftello64(stream);
#endif
#else
    return ftell(stream);
#endif
#endif
}


gmx_bool is_pipe(FILE *fp)
{
    tMPI_Thread_mutex_lock(&pstack_mutex);

    t_pstack *ps = pstack;
    while (ps != NULL)
    {
        if (ps->fp == fp)
        {
            tMPI_Thread_mutex_unlock(&pstack_mutex);
            return TRUE;
        }
        ps = ps->prev;
    }
    tMPI_Thread_mutex_unlock(&pstack_mutex);
    return FALSE;
}


static FILE *uncompress(const char *fn, const char *mode)
{
    FILE *fp;
    char  buf[256];

    sprintf(buf, "uncompress -c < %s", fn);
    fprintf(stderr, "Going to execute '%s'\n", buf);
    if ((fp = popen(buf, mode)) == NULL)
    {
        gmx_open(fn);
    }
    push_ps(fp);

    return fp;
}

static FILE *gunzip(const char *fn, const char *mode)
{
    FILE *fp;
    char  buf[256];

    sprintf(buf, "gunzip -c < %s", fn);
    fprintf(stderr, "Going to execute '%s'\n", buf);
    if ((fp = popen(buf, mode)) == NULL)
    {
        gmx_open(fn);
    }
    push_ps(fp);

    return fp;
}

gmx_bool gmx_fexist(const char *fname)
{
    FILE *test;

    if (fname == NULL)
    {
        return FALSE;
    }
    test = fopen(fname, "r");
    if (test == NULL)
    {
        /*Windows doesn't allow fopen of directory - so we need to check this seperately */
        #ifdef GMX_NATIVE_WINDOWS
        DWORD attr = GetFileAttributes(fname);
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


gmx_bool gmx_fexist_master(const char *fname, t_commrec *cr)
{
    gmx_bool bExist;

    if (SIMMASTER(cr))
    {
        bExist = gmx_fexist(fname);
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(bExist), &bExist, cr);
    }
    return bExist;
}

gmx_bool gmx_eof(FILE *fp)
{
    char     data[4];
    gmx_bool beof;

    if (is_pipe(fp))
    {
        return feof(fp);
    }
    else
    {
        if ((beof = fread(data, 1, 1, fp)) == 1)
        {
            gmx_fseek(fp, -1, SEEK_CUR);
        }
        return !beof;
    }
}

static char *backup_fn(const char *file, int count_max)
{
    /* Use a reasonably low value for countmax; we might
     * generate 4-5 files in each round, and we dont
     * want to hit directory limits of 1024 or 2048 files.
     */
#define COUNTMAX 99
    int          i, count = 1;
    char        *directory, *fn;
    char        *buf;

    if (count_max == -1)
    {
        count_max = COUNTMAX;
    }

    smalloc(buf, GMX_PATH_MAX);

    for (i = strlen(file)-1; ((i > 0) && (file[i] != DIR_SEPARATOR)); i--)
    {
        ;
    }
    /* Must check whether i > 0, i.e. whether there is a directory
     * in the file name. In that case we overwrite the / sign with
     * a '\0' to end the directory string .
     */
    if (i > 0)
    {
        directory    = gmx_strdup(file);
        directory[i] = '\0';
        fn           = gmx_strdup(file+i+1);
    }
    else
    {
        directory    = gmx_strdup(".");
        fn           = gmx_strdup(file);
    }
    do
    {
        sprintf(buf, "%s/#%s.%d#", directory, fn, count);
        count++;
    }
    while ((count <= count_max) && gmx_fexist(buf));

    /* Arbitrarily bail out */
    if (count > count_max)
    {
        gmx_fatal(FARGS, "Won't make more than %d backups of %s for you.\n"
                  "The env.var. GMX_MAXBACKUP controls this maximum, -1 disables backups.",
                  count_max, fn);
    }

    sfree(directory);
    sfree(fn);

    return buf;
}

gmx_bool make_backup(const char * name)
{
    char * env;
    int    count_max;
    char * backup;

#ifdef GMX_FAHCORE
    return FALSE; /* skip making backups */
#else

    if (gmx_fexist(name))
    {
        env = getenv("GMX_MAXBACKUP");
        if (env != NULL)
        {
            count_max = strtol(env, NULL, 10);
            if (count_max == -1)
            {
                /* Do not make backups and possibly overwrite old files */
                return TRUE;
            }
        }
        else
        {
            /* Use the default maximum */
            count_max = -1;
        }
        backup = backup_fn(name, count_max);
        if (rename(name, backup) == 0)
        {
            fprintf(stderr, "\nBack Off! I just backed up %s to %s\n",
                    name, backup);
        }
        else
        {
            fprintf(stderr, "Sorry couldn't backup %s to %s\n", name, backup);
            return FALSE;
        }
        sfree(backup);
    }
    return TRUE;
#endif
}

FILE *gmx_ffopen(const char *file, const char *mode)
{
#ifdef SKIP_FFOPS
    return fopen(file, mode);
#else
    FILE    *ff = NULL;
    char     buf[256], *bufsize = 0, *ptr;
    gmx_bool bRead;
    int      bs;

    if (file == NULL)
    {
        return NULL;
    }

    if (mode[0] == 'w')
    {
        make_backup(file);
    }
    where();

    bRead = (mode[0] == 'r' && mode[1] != '+');
    strcpy(buf, file);
    if (!bRead || gmx_fexist(buf))
    {
        if ((ff = fopen(buf, mode)) == NULL)
        {
            gmx_file(buf);
        }
        where();
        /* Check whether we should be using buffering (default) or not
         * (for debugging)
         */
        if (bUnbuffered || ((bufsize = getenv("GMX_LOG_BUFFER")) != NULL))
        {
            /* Check whether to use completely unbuffered */
            if (bUnbuffered)
            {
                bs = 0;
            }
            else
            {
                bs = strtol(bufsize, NULL, 10);
            }
            if (bs <= 0)
            {
                setbuf(ff, NULL);
            }
            else
            {
                snew(ptr, bs+8);
                if (setvbuf(ff, ptr, _IOFBF, bs) != 0)
                {
                    gmx_file("Buffering File");
                }
            }
        }
        where();
    }
    else
    {
        sprintf(buf, "%s.Z", file);
        if (gmx_fexist(buf))
        {
            ff = uncompress(buf, mode);
        }
        else
        {
            sprintf(buf, "%s.gz", file);
            if (gmx_fexist(buf))
            {
                ff = gunzip(buf, mode);
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

/* Our own implementation of dirent-like functionality to scan directories. */
struct gmx_directory
{
#if defined(GMX_NATIVE_WINDOWS)
    intptr_t             windows_handle;
    struct _finddata_t   finddata;
    int                  first;
#elif defined(HAVE_DIRENT_H)
    DIR  *               dirent_handle;
#else
    int                  dummy;
#endif
};


int
gmx_directory_open(gmx_directory_t *p_gmxdir, const char *dirname)
{
    struct gmx_directory *  gmxdir;
    int                     rc;

    snew(gmxdir, 1);

    *p_gmxdir = gmxdir;

#if defined(GMX_NATIVE_WINDOWS)
    if (dirname != NULL && strlen(dirname) > 0)
    {
        char *     tmpname;
        int        len;

        len = strlen(dirname);
        snew(tmpname, len+3);

        strncpy(tmpname, dirname, len+1);

        /* Remove possible trailing directory separator */
        if (tmpname[len] == '/' || tmpname[len] == '\\')
        {
            tmpname[len] = '\0';
        }

        /* Add wildcard */
        strcat(tmpname, "/*");

        gmxdir->first = 1;
        if ( (gmxdir->windows_handle = _findfirst(tmpname, &gmxdir->finddata)) > 0L)
        {
            rc = 0;
        }
        else
        {
            if (errno == EINVAL)
            {
                sfree(gmxdir);
                *p_gmxdir = NULL;
                rc        = EINVAL;
            }
            else
            {
                rc        = 0;
            }
        }
    }
    else
    {
        rc = EINVAL;
    }
#elif defined(HAVE_DIRENT_H)
    if ( (gmxdir->dirent_handle = opendir(dirname)) != NULL)
    {
        rc = 0;
    }
    else
    {
        sfree(gmxdir);
        *p_gmxdir = NULL;
        rc        = EINVAL;
    }
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n"
              "In the very unlikely event this is not a compile-time mistake you could consider\n"
              "implementing support for your platform in futil.c, but contact the developers\n"
              "to make sure it's really necessary!\n");
    rc = -1;
#endif
    return rc;
}


int
gmx_directory_nextfile(gmx_directory_t gmxdir, char *name, int maxlength_name)
{
    int                     rc;

#if defined(GMX_NATIVE_WINDOWS)
    if (gmxdir != NULL)
    {
        if (gmxdir->windows_handle <= 0)
        {

            name[0] = '\0';
            rc      = ENOENT;
        }
        else if (gmxdir->first == 1)
        {
            strncpy(name, gmxdir->finddata.name, maxlength_name);
            rc            = 0;
            gmxdir->first = 0;
        }
        else
        {
            if (_findnext(gmxdir->windows_handle, &gmxdir->finddata) == 0)
            {
                strncpy(name, gmxdir->finddata.name, maxlength_name);
                rc      = 0;
            }
            else
            {
                name[0] = '\0';
                rc      = ENOENT;
            }
        }
    }
    else
    {
        name[0] = '\0';
        rc      = EINVAL;
    }
#elif defined(HAVE_DIRENT_H)
    struct dirent *         direntp_large;
    struct dirent *         p;


    if (gmxdir != NULL && gmxdir->dirent_handle != NULL)
    {
        /* On some platforms no space is present for d_name in dirent.
         * Since d_name is guaranteed to be the last entry, allocating
         * extra space for dirent will allow more size for d_name.
         * GMX_MAX_PATH should always be >= the max possible d_name.
         */
        smalloc(direntp_large, sizeof(*direntp_large) + GMX_PATH_MAX);
        rc = readdir_r(gmxdir->dirent_handle, direntp_large, &p);

        if (p != NULL && rc == 0)
        {
            strncpy(name, direntp_large->d_name, maxlength_name);
        }
        else
        {
            name[0] = '\0';
            rc      = ENOENT;
        }
        sfree(direntp_large);
    }
    else
    {
        name[0] = '\0';
        rc      = EINVAL;
    }
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n");
    rc = -1;
#endif
    return rc;
}


int
gmx_directory_close(gmx_directory_t gmxdir)
{
    int                     rc;
#if defined(GMX_NATIVE_WINDOWS)
    rc = (gmxdir != NULL) ? _findclose(gmxdir->windows_handle) : EINVAL;
#elif defined(HAVE_DIRENT_H)
    rc = (gmxdir != NULL) ? closedir(gmxdir->dirent_handle) : EINVAL;
#else
    gmx_fatal(FARGS,
              "Source compiled without POSIX dirent or windows support - cannot scan directories.\n");
    rc = -1;
#endif

    sfree(gmxdir);
    return rc;
}


char *low_gmxlibfn(const char *file, gmx_bool bAddCWD, gmx_bool bFatal)
{
    bool bEnvIsSet = false;
    try
    {
        if (bAddCWD && gmx_fexist(file))
        {
            return gmx_strdup(file);
        }
        else
        {
            std::string  libpath;
            // GMXLIB can be a path.
            const char  *lib = getenv("GMXLIB");
            if (lib != NULL)
            {
                bEnvIsSet = true;
                libpath   = lib;
            }
            else
            {
                libpath = gmx::getProgramContext().defaultLibraryDataPath();
            }

            std::vector<std::string>                 pathEntries;
            gmx::Path::splitPathEnvironment(libpath, &pathEntries);
            std::vector<std::string>::const_iterator i;
            for (i = pathEntries.begin(); i != pathEntries.end(); ++i)
            {
                std::string testPath = gmx::Path::join(*i, file);
                if (gmx::Path::exists(testPath))
                {
                    return gmx_strdup(testPath.c_str());
                }
            }
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    if (bFatal)
    {
        if (bEnvIsSet)
        {
            gmx_fatal(FARGS,
                      "Library file %s not found %sin your GMXLIB path.",
                      file, bAddCWD ? "in current dir nor " : "");
        }
        else
        {
            gmx_fatal(FARGS,
                      "Library file %s not found %sin default directories.\n"
                      "(You can set the directories to search with the GMXLIB path variable)",
                      file, bAddCWD ? "in current dir nor " : "");
        }
    }
    return NULL;
}

FILE *low_libopen(const char *file, gmx_bool bFatal)
{
    FILE *ff;
    char *fn;

    fn = low_gmxlibfn(file, TRUE, bFatal);

    if (fn == NULL)
    {
        ff = NULL;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "Opening library file %s\n", fn);
        }
        ff = fopen(fn, "r");
    }
    sfree(fn);

    return ff;
}

char *gmxlibfn(const char *file)
{
    return low_gmxlibfn(file, TRUE, TRUE);
}

FILE *libopen(const char *file)
{
    return low_libopen(file, TRUE);
}

void gmx_tmpnam(char *buf)
{
    int i, len;

    if ((len = strlen(buf)) < 7)
    {
        gmx_fatal(FARGS, "Buf passed to gmx_tmpnam must be at least 7 bytes long");
    }
    for (i = len-6; (i < len); i++)
    {
        buf[i] = 'X';
    }
    /* mktemp is dangerous and we should use mkstemp instead, but
     * since windows doesnt support it we have to separate the cases.
     * 20090307: mktemp deprecated, use iso c++ _mktemp instead.
     */
#ifdef GMX_NATIVE_WINDOWS
    _mktemp(buf);
#else
    int fd = mkstemp(buf);

    switch (fd)
    {
        case EINVAL:
            gmx_fatal(FARGS, "Invalid template %s for mkstemp", buf);
            break;
        case EEXIST:
            gmx_fatal(FARGS, "mkstemp created existing file", buf);
            break;
        case EACCES:
            gmx_fatal(FARGS, "Permission denied for opening %s", buf);
            break;
        default:
            break;
    }
    close(fd);
#endif
    /* name in Buf should now be OK */
}

int gmx_truncatefile(char *path, gmx_off_t length)
{
#ifdef GMX_NATIVE_WINDOWS
    /* Microsoft visual studio does not have "truncate" */
    HANDLE        fh;
    LARGE_INTEGER win_length;

    win_length.QuadPart = length;

    fh = CreateFile(path, GENERIC_READ | GENERIC_WRITE, 0, NULL,
                    OPEN_EXISTING, 0, NULL);
    SetFilePointerEx(fh, win_length, NULL, FILE_BEGIN);
    SetEndOfFile(fh);
    CloseHandle(fh);

    return 0;
#else
    return truncate(path, length);
#endif
}


int gmx_file_rename(const char *oldname, const char *newname)
{
#ifndef GMX_NATIVE_WINDOWS
    /* under unix, rename() is atomic (at least, it should be). */
    return rename(oldname, newname);
#else
    if (MoveFileEx(oldname, newname,
                   MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
    {
        return 0;
    }
    else
    {
        return 1;
    }
#endif
}

int gmx_file_copy(const char *oldname, const char *newname, gmx_bool copy_if_empty)
{
/* the full copy buffer size: */
#define FILECOPY_BUFSIZE (1<<16)
    FILE *in  = NULL;
    FILE *out = NULL;
    char *buf;

    snew(buf, FILECOPY_BUFSIZE);

    in = fopen(oldname, "rb");
    if (!in)
    {
        goto error;
    }

    /* If we don't copy when empty, we postpone opening the file
       until we're actually ready to write. */
    if (copy_if_empty)
    {
        out = fopen(newname, "wb");
        if (!out)
        {
            goto error;
        }
    }

    while (!feof(in))
    {
        size_t nread;

        nread = fread(buf, sizeof(char), FILECOPY_BUFSIZE, in);
        if (nread > 0)
        {
            size_t ret;
            if (!out)
            {
                /* so this is where we open when copy_if_empty is false:
                   here we know we read something. */
                out = fopen(newname, "wb");
                if (!out)
                {
                    goto error;
                }
            }
            ret = fwrite(buf, sizeof(char), nread, out);
            if (ret != nread)
            {
                goto error;
            }
        }
        if (ferror(in))
        {
            goto error;
        }
    }
    sfree(buf);
    fclose(in);
    fclose(out);
    return 0;
error:
    sfree(buf);
    if (in)
    {
        fclose(in);
    }
    if (out)
    {
        fclose(out);
    }
    return 1;
#undef FILECOPY_BUFSIZE
}


int gmx_fsync(FILE *fp)
{
    int rc = 0;

#ifdef GMX_FAHCORE
    /* the fahcore defines its own os-independent fsync */
    rc = fah_fsync(fp);
#else /* GMX_FAHCORE */
    {
        int fn = -1;

        /* get the file number */
#if defined(HAVE_FILENO)
        fn = fileno(fp);
#elif defined(HAVE__FILENO)
        fn = _fileno(fp);
#endif

        /* do the actual fsync */
        if (fn >= 0)
        {
#if (defined(HAVE_FSYNC))
            rc = fsync(fn);
#elif (defined(HAVE__COMMIT))
            rc = _commit(fn);
#endif
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

void gmx_chdir(const char *directory)
{
#ifdef GMX_NATIVE_WINDOWS
    int rc = _chdir(directory);
#else
    int rc = chdir(directory);
#endif
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Cannot change directory to '%s'. Reason: %s",
                  directory, strerror(errno));
    }
}

void gmx_getcwd(char *buffer, size_t size)
{
#ifdef GMX_NATIVE_WINDOWS
    char *pdum = _getcwd(buffer, size);
#else
    char *pdum = getcwd(buffer, size);
#endif
    if (pdum == NULL)
    {
        gmx_fatal(FARGS, "Cannot get working directory. Reason: %s",
                  strerror(errno));
    }
}
