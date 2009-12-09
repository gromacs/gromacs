/*
 * history.c
 *
 *  Created on: Nov 27, 2009
 *      Author: rschulz
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include <string.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "smalloc.h"
#include "filenm.h"
#include "futil.h"
#include "md5.h"
#include "copyrite.h"

int _argc;
char **_argv;

/*
 * data:
 * 
 * Test:
 * 
 * 
 * replace all fgets,fgetc,scanf,fscanf calls with gmx_scanf
 * use vsprintf to add it to cmdline
*/


typedef struct {
    FILE* file;
    char* mode;
    char* fn;
    md5_byte_t md5sum[16];
} t_hist_file;

#define MAX_FILES 200
#define MAX_STDINPUT 200
int _nfile=0;
t_hist_file _files[MAX_FILES];
#define MAX_LINES 200
char* histbuf[MAX_LINES];
int nhistbuf=0;
FILE* histfile;

void histopenfile(FILE* file, const char* fn, const char* mode) {
    _files[_nfile].file = file;
    _files[_nfile].fn = strdup(fn);
    _files[_nfile].mode = strdup(mode);
    _nfile++;
}

char *
make_message(const char *fmt, ...)
{
    /* Guess we need no more than 100 bytes. */
    int n, size = 100;
    char *p, *np;
    va_list ap;

    if ((p = malloc(size)) == NULL)
        return NULL;

    while (1) {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf(p, size, fmt, ap);
        va_end(ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            return p;
        /* Else try again with more space. */
        if (n > -1)    /* glibc 2.1 */
            size = n+1; /* precisely what is needed */
        else           /* glibc 2.0 */
            size *= 2;  /* twice the old size */
        if ((np = realloc (p, size)) == NULL) {
            free(p);
            return NULL;
        } else {
            p = np;
        }
    }
}


void histaddinput(char* str)
{
    char* npos = rindex(str,'\n');
    while(npos && strlen(str)>0)
    {
        histbuf[nhistbuf++]=make_message("INP: %.*s\n",npos-str,str);
        str=npos+1;
        npos = rindex(str,'\n');
    }
    if (strlen(str)>0)
    {
        histbuf[nhistbuf++]=make_message("INP: %s\n",str);   
    }
}

static FILE* print_file(int i) 
{
    const char* fnm = _files[i].fn;
    const char* mode = _files[i].mode;
    struct stat sbuf;
    FILE* file = _files[i].file;
#define CPT_CHK_LEN  1048576 
    md5_state_t state;
    md5_byte_t digest[16];
    md5_byte_t buf[CPT_CHK_LEN+4];
    const char* cmd;
    char digest_str[33];
    bool bOpened=FALSE;
    int read_len;
    int j,ret;
    
    if (file==NULL)
    {
        file = fopen(fnm,"r");
        if (file==NULL)
        {
            gmx_fatal(FARGS,"Failed to open file: %s.",fnm);
        }
        bOpened = TRUE;
    } else {
        if (mode[0]!='r' || mode[1]!='+' || 
            !fseek(file, 0, SEEK_SET))
        {
            ret=fclose(file);  /* if can't write or seek does not work reopen */
            if(ret!=0) 
            {
                gmx_fatal(FARGS,"Failed to close file: %s.",fnm);
            }
            file = fopen(fnm,"r");
            if(file==NULL)
            {
                gmx_fatal(FARGS,"Failed to open file: %s.",fnm);
            }
        }
    }
    
    ret=fstat(fileno(file),&sbuf);
    if (ret!=0)
    {
        gmx_fatal(FARGS,"Failed to get file stat for file: %s",fnm);
    }
    
    if (sbuf.st_size>CPT_CHK_LEN)  /*file bigger -> read more*/
    {
        read_len=fread(buf,1,CPT_CHK_LEN/2,file);
        fseek(file,-CPT_CHK_LEN/2,SEEK_END);
        read_len+=fread(buf+CPT_CHK_LEN/2,1,CPT_CHK_LEN/2,file);
        read_len+=sizeof(off_t);
        *((gmx_large_int_t*)(buf+CPT_CHK_LEN))=sbuf.st_size;
    }
    else 
    {
        read_len=fread(buf,1,CPT_CHK_LEN,file);
    }
    if (bOpened)
    {
        ret=fclose(file);
        if (ret!=0) 
        {
            gmx_fatal(FARGS,"Failed to close file: %s.",fnm);
        }
    }
    
    md5_init(&state);
    md5_append(&state, buf, read_len);
    md5_finish(&state, digest);
    
    /* checking for already printed identical file 
     * md5sum is only set if the file has already been printed to log*/
    for (j=0;j<_nfile;j++) 
    {
        if (strcmp(_files[j].fn,fnm)==0 && _files[j].mode[0]==mode[0] 
            && memcmp(_files[j].md5sum,digest,sizeof(digest))==0)
        {
            return file;
        }
    }
    
    memcpy(_files[i].md5sum,digest,sizeof(digest));
        
    if (mode[0]=='r')  /*TODO: what about r+? is this reading or writing (as in win_truncate, checkpoint chksum_file)*/
    {
        cmd="IN : ";
    }
    else 
    {
        cmd="OUT: ";
    }
    for (j=0; j<16; j++)
    {
        sprintf(digest_str+j*2,"%02x",digest[j]);
    }
    histbuf[nhistbuf++]=make_message("%s%s %s\n",cmd,fnm,digest_str);
    
    return file;
}

int  histclosefile(FILE** file) {
    int i;
    if (file)
    {
        for (i=0;i<_nfile;i++)
        {
            if (*file==_files[i].file)
            {
                *file = print_file(i);
                _files[i].file=NULL;
            }
        }
    }
    return 0;
}

/* return: has to be freed */
char* getuser() {
    const char* name = NULL;
#ifdef __MINGW32__
    name = getenv("USERNAME");
#endif
    if (!name || !*name)
        name = getenv("LOGNAME");

    if (!name || !*name)
        name = getenv("USER");

#ifdef HAVE_UNISTD_H
#ifndef __MINGW32__
    if (!name || !*name) {
        struct passwd *p = getpwuid(getuid());
        if (p && p->pw_name && *p->pw_name)
            name = p->pw_name;
    }
#endif
#endif
    if (name==NULL) 
    {
        return NULL;
    }
    else 
    {
        return strdup(name);
    }
}    

void init_history(int argc, char** argv) {
    int i;

    _argc=argc;
    snew(_argv,argc);
    for (i=0;i<argc;i++) 
    {
        _argv[i]=strdup(argv[i]);
    }
    
}


void print_history() {
    int i,j,ret;
    char* lfn;
    const char* histfnm;
    time_t t;
    char host[256];
    char *user;
    char pwd[GMX_PATH_MAX];

    for (i=0;i<_nfile;i++) {
        if (_files[i].file!=NULL) {
            fprintf(stderr,"BUG: %s was not closed correctly\n",_files[i].fn);
            print_file(i);
            /*if (_files[i].file!=NULL) //don't close because it might already be
                                       //closed without calling histclosefile   
            {   
                ret = fclose(_files[i].file); 
                if (ret!=0)  { gmx_fatal(FARGS,"Failed to close");}
            }*/
        }
        sfree(_files[i].fn);
        sfree(_files[i].mode);
    }
    histfnm = getenv("GMX_HISTFILE");
    if (histfnm==NULL)
    {
        histfnm=".gmx_history";
    }
    if (strcmp(histfnm,"")==0)  /*don't write history*/
    {
        return;
    }
    histfile = fopen(histfnm,"a");  
    if (histfile==NULL)
    {
        histfile = fopen(histfnm,"w");
        if (histfile==NULL)
        {
            gmx_fatal(errno,__FILE__,__LINE__,"Could not open history file: %s.",histfnm);
        }
    }
    gmx_lock(histfile);

    time(&t);
    user = getuser();
#ifdef HAVE_UNISTD_H
    if( gethostname(host,255) != 0)
    {
        sprintf(host,"unknown");
    }
    if (getcwd(pwd,sizeof(pwd))==NULL)
    {
        strcpy(pwd,"unknown");
    }
#else
    sprintf(host,"unknown");
    strcpy(pwd,"unknown");
#endif  

    fprintf(histfile,"\nCMD: ");
    for (i=0;i<_argc;i++)
    {
        fprintf(histfile,"%s ", _argv[i]);
    }
    
    fprintf(histfile,"\nPWD: %s\nBY : %s@%s %s",pwd,user,host,ctime(&t));
    sfree(user);
    
#if (defined BUILD_TIME && defined BUILD_USER) 
    fprintf(histfile,"VER: %s %s %s\n",VERSION,BUILD_USER,BUILD_TIME);
#endif

    for (i=0;i<nhistbuf;i++)
    {
        fwrite(histbuf[i],1,strlen(histbuf[i]),histfile);
        sfree(histbuf[i]);
    }
    ret=fclose(histfile);
    if (ret!=0)
    {
        gmx_fatal(FARGS,"Failed to close file: %s.",histfnm);
    }
}
