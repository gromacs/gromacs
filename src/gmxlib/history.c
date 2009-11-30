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

#ifdef USING_ARGS
int _nfile;
t_filenm *_fnm;
#else
typedef struct {
    FILE* file;
    char* mode;
    char* fn;
} t_hist_file;

#define MAX_FILES 200
int _nfile=0;
t_hist_file _files[MAX_FILES];
#endif

void histopenfile(FILE* file, const char* fn, const char* mode) {
#ifndef USING_ARGS
    _files[_nfile].file = file;
    _files[_nfile].fn = strdup(fn);
    _files[_nfile].mode = strdup(mode);
    _nfile++;
#endif
}

static FILE* print_file(const char* fnm, const char* mode, FILE* file) {
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
    
    if (mode[0]=='r')  /*TODO: what about r+? is this reading or writing (as in win_truncate, checkpoint chksum_file)*/
    {
        printf("IN : ");
    }
    else 
    {
        printf("OUT: ");
    }
    printf("%s ", fnm);
    for (j=0; j<16; j++)
    {
        printf("%02x",digest[j]);
    }
    printf("\n");
    
    return file;
}

int  histclosefile(FILE** file) {
#ifndef USING_ARGS
    int i;
    if (file)
    {
        for (i=0;i<_nfile;i++)
        {
            if (*file==_files[i].file)
            {
                *file = print_file(_files[i].fn,_files[i].mode,_files[i].file);
                _files[i].file=NULL;
                sfree(_files[i].fn);
                sfree(_files[i].mode);
            }
        }
    }
#endif
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
    int i;
    _argc=argc;
    snew(_argv,argc);
    for (i=0;i<argc;i++) 
    {
        _argv[i]=strdup(argv[i]);
    }
#ifdef USING_ARGS
    _nfile=nfile;
    _fnm=fnm;
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

#ifdef USING_ARGS
void print_files() {
    int i,j;
    char* lfn;
    for (i=0; i<_nfile;i++) 
      {
          /*if ((!(_fnm[i].flag & ffOPT)) || (_fnm[i].flag & ffSET))*/
          /* if a file is optional it is currently unknown whether it will be read or written
           * (usually) a ndx is only read if it is also set
           * but e.g. in trjconv the tpx might also be read if it is not set
           * for writing it is anyhow not known
           * testing whether the file exists is just the best we can do at the moment without
           * the applications notifiying the history lib which files were read or written
           * This probably could be best done by adding a flag ffUSED which is set if an 
           * optional file was read/written.
           * Even for writing test file is not optimal because file might be there thus files
           * which were not generated in this step are added as output file.
           * 
           * workaround: use access and modify time. But that does not work if files are changed
           * also by other processes. This is also a problem for the md5sum
           * to make better for md5sum: keep filehandle and seek instead of reopening. 
           * But not clear how to do that. 
           * Could be done by putting it in ffopen and ffclose
           * this would also automatically get a list of really all files read/written
           */
          for (j=0;j<_fnm[i].nfiles;j++)
          {
              if (_fnm[i].flag & ffLIB)
              {
                  lfn = low_libfn(_fnm[i].fns[j],FALSE);
                  if (lfn!=NULL && gmx_fexist(lfn))
                  {
                      print_file(lfn,(_fnm[i].flag&ffREAD)?"r":"w",NULL);
                      sfree(lfn);
                  }
              }
              else
              {
                  if (gmx_fexist(_fnm[i].fns[j]))
                  {
                      print_file(_fnm[i].fns[j],(_fnm[i].flag&ffREAD)?"r":"w",NULL);
                  }
              }
          }
      }
}
#endif

void print_history() {
    int i,j;
    char* lfn;
    time_t t;
    char host[256];
    char *user;
    char pwd[GMX_PATH_MAX]="unknown";

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

    printf("CMD: ");
    for (i=0;i<_argc;i++)
    {
        printf("%s ", _argv[i]);
    }
    
    printf("\nPWD: %s\nBY : %s@%s %s",pwd,user,host,ctime(&t));
    sfree(user);
    
#ifdef USING_ARGS
    print_files();
#else
    for (i=0;i<_nfile;i++) {
        if (_files[i].file!=NULL) {
            fprintf(stderr,"BUG: %s was not closed correctly\n",_files[i].fn);
            print_file(_files[i].fn,_files[i].mode,_files[i].file);
            if (_files[i].file!=NULL)
            {   
                fclose(_files[i].file);
            }
            sfree(_files[i].fn);
            sfree(_files[i].mode);
        }
    }
#endif
#if (defined BUILD_TIME && defined BUILD_USER) 
    printf("VER: %s %s %s\n",VERSION,BUILD_USER,BUILD_TIME);
#endif
    
}
