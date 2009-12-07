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

//#define USING_ARGS 1  /*define for using command line arguments, undef for using ffopen, fflclose calls*/

/* currently does not work without USING_ARGS: segfaults 
 * mainly because at many places fclose is directly called instead of ffclose*/
/* while it has the advantage that it doesn't show unused optional files
 * it currently shows some input files double - would have to be filtered
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
int _ninput=0;
char* _stdinput[MAX_STDINPUT];

FILE* histfile;

void histopenfile(FILE* file, const char* fn, const char* mode) {
    _files[_nfile].file = file;
    _files[_nfile].fn = strdup(fn);
    _files[_nfile].mode = strdup(mode);
    _nfile++;
}

void histaddinput(char* str)
{
    char* npos = rindex(str,'\n');
    if (npos)
    {
        _stdinput[_ninput] = strndup(str,npos-str);
    }
    else
    {
        _stdinput[_ninput] = strdup(str);
    }
    _ninput++;
}
/*
static FILE* print_file(const char* fnm, const char* mode, FILE* file) {
*/
static FILE* print_file(int i) 
{
    const char* fnm = _files[i].fn;
    const char* mode = _files[i].mode;
    FILE* file = _files[i].file;
#define CPT_CHK_LEN  1048576 
    md5_state_t state;
    md5_byte_t digest[16];
    md5_byte_t buf[CPT_CHK_LEN];
    bool bOpened=FALSE;
    int read_len;
    int j;
    
    if (file==NULL)
    {
        file = fopen(fnm,"r");
        bOpened = TRUE;
    } else {
        if (mode[0]!='r' || mode[1]!='+' || 
            !fseek(file, 0, SEEK_SET))
        {
            fclose(file);  /* if can't write or seek does not work reopen */
            file = fopen(fnm,"r");
        }
    }
    read_len = fread(buf,1,CPT_CHK_LEN,file);
    if (bOpened)
    {
        fclose(file);
    }
    
    md5_init(&state);
    md5_append(&state, buf, read_len);
    md5_finish(&state, digest);
    
    /* looking for earlier identical files 
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
        fprintf(histfile,"IN : ");
    }
    else 
    {
        fprintf(histfile,"OUT: ");
    }
    fprintf(histfile,"%s ", fnm);
    for (j=0; j<16; j++)
    {
        fprintf(histfile,"%02x",digest[j]);
    }
    fprintf(histfile,"\n");
    
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

/* can we pass just the copy without making a copy? 
 * we need to call init_history from parse_common_args before argc/argv is changed
 * at this point fnm is not yet filled with the correct fnms
 * thus we try it just with the pointers first
 */ 
void init_history(int argc, char** argv, int nfile, t_filenm *fnm) {
    int i,j;
    char* lfn;
    time_t t;
    char host[256];
    char *user;
    char pwd[GMX_PATH_MAX]="unknown";

    _argc=argc;
    snew(_argv,argc);
    for (i=0;i<argc;i++) 
    {
        _argv[i]=strdup(argv[i]);
    }
    histfile = fopen(".gmx_history","a");

    time(&t);
    user = getuser();
#ifdef HAVE_UNISTD_H
    if( gethostname(host,255) != 0)
    {
        sprintf(host,"unknown");
    }
    getwd(pwd);
#else
    sprintf(host,"unknown");
#endif  

    fprintf(histfile,"\nCMD: ");
    if (_ninput>0)
    {
        fprintf(histfile,"echo ");
        for (i=0;i<_ninput;i++)
        {
            fprintf(histfile,"%s ",_stdinput[i]);
            sfree(_stdinput[i]);
        }
        fprintf(histfile,"| ");
    }
    for (i=0;i<_argc;i++)
    {
        fprintf(histfile,"%s ", _argv[i]);
    }
    
    fprintf(histfile,"\nPWD: %s\nBY : %s@%s %s",pwd,user,host,ctime(&t));
    sfree(user);
    
#if (defined BUILD_TIME && defined BUILD_USER) 
    fprintf(histfile,"VER: %s %s %s\n",VERSION,BUILD_USER,BUILD_TIME);
#endif

}
/*
 * data:
 * md5 each file (not 1st 1MB but 1/2 from beginning and 1/2 from end + filesize)
 * 
 * Functionality
 * file locking on history file
 * writing to history file
 * location from ENV_VARIABLE (not writing if this is "")
 * history read tool
 * 
 * write stdin to cmdline: echo ... | 
 * replace all fgets,fgetc,scanf,fscanf calls with gmx_scanf
 * use vsprintf to add it to cmdline
*/


void print_history() {
    int i;
    for (i=0;i<_nfile;i++) {
        if (_files[i].file!=NULL) {
            fprintf(stderr,"BUG: %s was not closed correctly\n",_files[i].fn);
            print_file(i);
            if (_files[i].file!=NULL)
            {   
                fclose(_files[i].file);
            }
        }
        sfree(_files[i].fn);
        sfree(_files[i].mode);
    }
    fclose(histfile);
}
