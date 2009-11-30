/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <direct.h>
#include <io.h>
#endif

#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "statutil.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* Windows file stuff, only necessary for visual studio */
#ifdef _MSC_VER
#include "windows.h"
#endif

/* we keep a linked list of all files opened through pipes (i.e. 
   compressed or .gzipped files. This way we can distinguish between them
   without having to change the semantics of reading from/writing to files) 
   */
typedef struct t_pstack {
    FILE   *fp;
    struct t_pstack *prev;
} t_pstack;

static t_pstack *pstack=NULL;
static bool     bUnbuffered=FALSE;

#ifdef GMX_THREADS
/* this linked list is an intrinsically globally shared object, so we have
   to protect it with mutexes */
static tMPI_Thread_mutex_t pstack_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif

void no_buffers(void)
{
    bUnbuffered=TRUE;
}

void push_ps(FILE *fp)
{
    t_pstack *ps;

#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    snew(ps,1);
    ps->fp   = fp;
    ps->prev = pstack;
    pstack   = ps;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
}

#ifdef GMX_FAHCORE
/* redefine ffclose */
#define ffclose fah_fclose
#else
#ifdef ffclose
#undef ffclose
#endif
#endif


#ifndef HAVE_PIPES
static FILE *popen(const char *nm,const char *mode)
{
    gmx_impl("Sorry no pipes...");

    return NULL;
}

static int pclose(FILE *fp)
{
    gmx_impl("Sorry no pipes...");

    return 0;
}
#endif



int ffclose(FILE *fp)
{
    t_pstack *ps,*tmp;
    int ret=0;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    if (ps == NULL) {
        if (fp != NULL) 
            ret = fclose(fp);
    }
    else if (ps->fp == fp) {
        if (fp != NULL)
            ret = pclose(fp);
        pstack=pstack->prev;
        sfree(ps);
    }
    else {
        while ((ps->prev != NULL) && (ps->prev->fp != fp))
            ps=ps->prev;
        if (ps->prev->fp == fp) {
            if (ps->prev->fp != NULL)
                ret = pclose(ps->prev->fp);
            tmp=ps->prev;
            ps->prev=ps->prev->prev;
            sfree(tmp);
        }
        else {
            if (fp != NULL)
                ret = fclose(fp);
        }
    }
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
    return ret;
}

#ifdef rewind
#undef rewind
#endif

void frewind(FILE *fp)
{
    t_pstack *ps;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    while (ps != NULL) {
        if (ps->fp == fp) {
            fprintf(stderr,"Cannot rewind compressed file!\n");
#ifdef GMX_THREADS
            tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
            return;
        }
        ps=ps->prev;
    }
    rewind(fp);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
}

bool is_pipe(FILE *fp)
{
    t_pstack *ps;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&pstack_mutex);
#endif

    ps=pstack;
    while (ps != NULL) {
        if (ps->fp == fp) {
#ifdef GMX_THREADS
            tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
            return TRUE;
        }
        ps=ps->prev;
    }
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&pstack_mutex);
#endif
    return FALSE;
}


static FILE *uncompress(const char *fn,const char *mode)
{
    FILE *fp;
    char buf[256];

    sprintf(buf,"uncompress -c < %s",fn);
    fprintf(stderr,"Going to execute '%s'\n",buf);
    if ((fp=popen(buf,mode)) == NULL)
        gmx_open(fn);
    push_ps(fp);

    return fp;
}

static FILE *gunzip(const char *fn,const char *mode)
{
    FILE *fp;
    char buf[256];

    sprintf(buf,"gunzip -c < %s",fn);
    fprintf(stderr,"Going to execute '%s'\n",buf);
    if ((fp=popen(buf,mode)) == NULL)
        gmx_open(fn);
    push_ps(fp);

    return fp;
}

bool gmx_fexist(const char *fname)
{
    FILE *test;

    if (fname == NULL)
        return FALSE;
    test=fopen(fname,"r");
    if (test == NULL) 
        return FALSE;
    else {
        ffclose(test);
        return TRUE;
    }
}


bool gmx_fexist_master(const char *fname, t_commrec *cr)
{
  bool bExist;
  
  if (SIMMASTER(cr)) 
  {
      bExist = gmx_fexist(fname);
  }
  if (PAR(cr)) 
  {
      gmx_bcast(sizeof(bExist),&bExist,cr);
  }
  return bExist;
}

bool gmx_eof(FILE *fp)
{
    char data[4];
    bool beof;

    if (is_pipe(fp))
        return feof(fp);
    else {
        if ((beof=fread(data,1,1,fp))==1)
            fseek(fp,-1,SEEK_CUR);
        return !beof;
    }
}

char *backup_fn(const char *file)
{
    /* Use a reasonably low value for countmax; we might
     * generate 4-5 files in each round, and we dont
     * want to hit directory limits of 1024 or 2048 files.
     */
#define COUNTMAX 128
    int         i,count=1;
    char        *directory,*fn;
    char        *buf;

    smalloc(buf, GMX_PATH_MAX);

    for(i=strlen(file)-1; ((i > 0) && (file[i] != '/')); i--)
        ;
    /* Must check whether i > 0, i.e. whether there is a directory
     * in the file name. In that case we overwrite the / sign with
     * a '\0' to end the directory string .
     */
    if (i > 0) {
        directory    = strdup(file);
        directory[i] = '\0';
        fn           = strdup(file+i+1);
    }
    else {
        directory    = strdup(".");
        fn           = strdup(file);
    }
    do {
        sprintf(buf,"%s/#%s.%d#",directory,fn,count);
        count++;
    } while ((count < COUNTMAX) && gmx_fexist(buf));

    /* Arbitrarily bail out */
    if (count == COUNTMAX) 
        gmx_fatal(FARGS,"Won't make more than %d backups of %s for you",
                COUNTMAX,fn);

    sfree(directory);
    sfree(fn);

    return buf;
}

bool make_backup(const char * name)
{
    char * backup;

#ifdef GMX_FAHCORE
    return FALSE; /* skip making backups */
#else

    if(gmx_fexist(name)) {
        backup = backup_fn(name);
        if(rename(name, backup) == 0) {
            fprintf(stderr, "\nBack Off! I just backed up %s to %s\n",
                    name, backup);
        } else {
            fprintf(stderr, "Sorry couldn't backup %s to %s\n", name, backup);
            return FALSE;
        }
        sfree(backup);
    }
    return TRUE;
#endif
}

FILE *ffopen(const char *file,const char *mode)
{
    FILE *ff=NULL;
    char buf[256],*bf,*bufsize=0,*ptr;
    bool bRead;
    int  bs;

    if (mode[0]=='w') {
        make_backup(file);
    }
    where();

    bRead= (mode[0]=='r'&&mode[1]!='+');
    strcpy(buf,file);
    if (gmx_fexist(buf) || !bRead) {
        if ((ff=fopen(buf,mode))==NULL)
            gmx_file(buf);
        where();
        /* Check whether we should be using buffering (default) or not
         * (for debugging)
         */
        if (bUnbuffered || ((bufsize=getenv("LOG_BUFS")) != NULL)) {
            /* Check whether to use completely unbuffered */
            if (bUnbuffered)
                bs = 0;
            else
                bs=strtol(bufsize, NULL, 10); 
            if (bs <= 0)
                setbuf(ff,NULL); 
            else {
                snew(ptr,bs+8);
                if (setvbuf(ff,ptr,_IOFBF,bs) != 0)
                    gmx_file("Buffering File");
            }
        }
        where();
    }
    else {
        sprintf(buf,"%s.Z",file);
        if (gmx_fexist(buf)) {
            ff=uncompress(buf,mode);
        }
        else {
            sprintf(buf,"%s.gz",file);
            if (gmx_fexist(buf)) {
                ff=gunzip(buf,mode);
            }
            else 
                gmx_file(file);
        }
    }
    return ff;
}



bool search_subdirs(const char *parent, char *libdir)
{
    char *ptr;
    bool found;

    /* Search a few common subdirectory names for the gromacs library dir */
    sprintf(libdir,"%s%cshare%ctop%cgurgle.dat",parent,
            DIR_SEPARATOR,DIR_SEPARATOR,DIR_SEPARATOR);
    found=gmx_fexist(libdir);
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,
                DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs-%s%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,VERSION,
                DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    
    if(!found) {
        sprintf(libdir,"%s%cshare%cgromacs%cgromacs-%s%ctop%cgurgle.dat",parent,
                DIR_SEPARATOR,DIR_SEPARATOR,DIR_SEPARATOR,
                VERSION,DIR_SEPARATOR,DIR_SEPARATOR);
        found=gmx_fexist(libdir);
    }    

    /* Remove the gurgle.dat part from libdir if we found something */
    if(found) {
        ptr=strrchr(libdir,DIR_SEPARATOR); /* slash or backslash always present, no check necessary */
        *ptr='\0';
    }
    return found;
}


/* Check if the program name begins with "/" on unix/cygwin, or
 * with "\" or "X:\" on windows. If not, the program name
 * is relative to the current directory.
 */
static bool filename_is_absolute(char *name)
{
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    return ((name[0] == DIR_SEPARATOR) || ((strlen(name)>3) && strncmp(name+1,":\\",2)));
#else
    return (name[0] == DIR_SEPARATOR);
#endif
}

bool get_libdir(char *libdir)
{
    char bin_name[512];
    char buf[512];
    char full_path[GMX_PATH_MAX];
    char test_file[GMX_PATH_MAX];
    char system_path[GMX_PATH_MAX];
    char *dir,*ptr,*s,*pdum;
    bool found=FALSE;
    int i;

    /* First - detect binary name */
    strncpy(bin_name,Program(),512);

    /* On windows & cygwin we need to add the .exe extension
     * too, or we wont be able to detect that the file exists
     */
#if (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64 || defined __CYGWIN__ || defined __CYGWIN32__)
    if(strlen(bin_name)<3 || strncasecmp(bin_name+strlen(bin_name)-4,".exe",4))
        strcat(bin_name,".exe");
#endif

    /* Only do the smart search part if we got a real name */
    if (NULL!=bin_name && strncmp(bin_name,"GROMACS",512)) {

        if (!strchr(bin_name,DIR_SEPARATOR)) {
            /* No slash or backslash in name means it must be in the path - search it! */
            s=getenv("PATH");

            /* Add the local dir since it is not in the path on windows */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
            pdum=_getcwd(system_path,sizeof(system_path)-1);
#else
            pdum=getcwd(system_path,sizeof(system_path)-1);
#endif
            strcat(system_path,PATH_SEPARATOR);
            if (s != NULL)
                strcat(system_path,s);
            s=system_path;
            found=FALSE;
            while(!found && (dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL)
            {
                sprintf(full_path,"%s%c%s",dir,DIR_SEPARATOR,bin_name);
                found=gmx_fexist(full_path);
            }
            if (!found)
                return FALSE;
        } else if (!filename_is_absolute(bin_name)) {
            /* name contains directory separators, but 
             * it does not start at the root, i.e.
             * name is relative to the current dir 
             */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
            pdum=_getcwd(buf,sizeof(buf)-1);
#else
            pdum=getcwd(buf,sizeof(buf)-1);
#endif
            strncpy(full_path,buf,GMX_PATH_MAX);
            strcat(full_path,"/");
            strcat(full_path,bin_name);
        } else {
            strncpy(full_path,bin_name,GMX_PATH_MAX);
        }

        /* Now we should have a full path and name in full_path,
         * but on unix it might be a link, or a link to a link to a link..
         */
#if (!defined WIN32 && !defined _WIN32 && !defined WIN64 && !defined _WIN64)
        while( (i=readlink(full_path,buf,sizeof(buf)-1)) > 0 ) {
            buf[i]='\0';
            /* If it doesn't start with "/" it is relative */
            if (buf[0]!=DIR_SEPARATOR) {
                strncpy(strrchr(full_path,DIR_SEPARATOR)+1,buf,GMX_PATH_MAX);
            } else
                strncpy(full_path,buf,GMX_PATH_MAX);
        }
#endif

        /* Remove the executable name - it always contains at least one slash */
        *(strrchr(full_path,DIR_SEPARATOR)+1)='\0';
        /* Now we have the full path to the gromacs executable.
         * Use it to find the library dir. 
         */
        found=FALSE;
        while(!found && ( (ptr=strrchr(full_path,DIR_SEPARATOR)) != NULL ) ) {
            *ptr='\0';
            found=search_subdirs(full_path,libdir);
        }
    }
    /* End of smart searching. If we didn't find it in our parent tree,
     * or if the program name wasn't set, at least try some standard 
     * locations before giving up, in case we are running from e.g. 
     * a users home directory. This only works on unix or cygwin...
     */
#if ((!defined WIN32 && !defined _WIN32 && !defined WIN64 && !defined _WIN64) || defined __CYGWIN__ || defined __CYGWIN32__)
    if(!found) 
        found=search_subdirs("/usr/local",libdir);
    if(!found) 
        found=search_subdirs("/usr",libdir);
    if(!found) 
        found=search_subdirs("/opt",libdir);
#endif
    return found;
}


char *low_libfn(const char *file, bool bFatal)
{
    char *ret=NULL;
    char *lib,*dir;
    char buf[1024];
    char libpath[GMX_PATH_MAX];
    bool env_is_set=FALSE;
    char   *s,tmppath[GMX_PATH_MAX];
    bool found;

    /* GMXLIB can be a path now */
    lib=getenv("GMXLIB");
    if (lib != NULL) {
        env_is_set=TRUE;
        strncpy(libpath,lib,GMX_PATH_MAX);
    } 
    else if (!get_libdir(libpath))
        strncpy(libpath,GMXLIBDIR,GMX_PATH_MAX);

    if (gmx_fexist(file))
    {
        ret=strdup(file);
    }
    else 
    {
        found=FALSE;
        strncpy(tmppath,libpath,GMX_PATH_MAX);
        s=tmppath;
        while(!found && (dir=gmx_strsep(&s, PATH_SEPARATOR)) != NULL )
        {
            sprintf(buf,"%s%c%s",dir,DIR_SEPARATOR,file);
            found=gmx_fexist(buf);
        }
        if (bFatal && !found) 
        {
            if (env_is_set) 
                gmx_fatal(FARGS,"Library file %s not found in current dir nor in your GMXLIB path.\n",file);
            else
                gmx_fatal(FARGS,"Library file %s not found in current dir nor in default directories.\n"
                        "(You can set the directories to search with the GMXLIB path variable)",file);
        }
        ret=strdup(buf);
    }

    return ret;
}




FILE *low_libopen(const char *file,bool bFatal)
{
    FILE *ff;
    char *fn;

    fn=low_libfn(file,bFatal);

    if (fn==NULL) {
        ff=NULL;
    } else {
        if (bFatal)
            fprintf(stderr,"Opening library file %s\n",fn);
        ff=fopen(fn,"r");
    }
    sfree(fn);

    return ff;
}

char *libfn(const char *file)
{
    return low_libfn(file,TRUE);
}

FILE *libopen(const char *file)
{
    return low_libopen(file,TRUE);
}

void gmx_tmpnam(char *buf)
{
    int i,len,fd;

    if ((len = strlen(buf)) < 7)
        gmx_fatal(FARGS,"Buf passed to gmx_tmpnam must be at least 7 bytes long");
    for(i=len-6; (i<len); i++) {
        buf[i] = 'X';
    }
    /* mktemp is dangerous and we should use mkstemp instead, but 
     * since windows doesnt support it we have to separate the cases.
     * 20090307: mktemp deprecated, use iso c++ _mktemp instead.
     */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
    _mktemp(buf);
#else
    fd = mkstemp(buf);

    switch (fd) {
        case EINVAL:
            gmx_fatal(FARGS,"Invalid template %s for mkstemp",buf);
            break;
        case EEXIST:
            gmx_fatal(FARGS,"mkstemp created existing file",buf);
            break;
        case EACCES: 
            gmx_fatal(FARGS,"Permission denied for opening %s",buf);
            break;
        default:
            break;
    }   
    close(fd);
#endif
    /* name in Buf should now be OK */
}

    int
gmx_truncatefile(char *path, off_t length)
{
#ifdef _MSC_VER
    /* Microsoft visual studio does not have "truncate" */
    HANDLE fh;
    LARGE_INTEGER win_length;

    win_length.QuadPart = length;

    fh = CreateFile(path,GENERIC_READ | GENERIC_WRITE,0,NULL,
            OPEN_EXISTING,0,NULL);
    SetFilePointerEx(fh,win_length,NULL,FILE_BEGIN);
    SetEndOfFile(fh);
    CloseHandle(fh);

    return 0;
#else
    return truncate(path,length);
#endif
}
