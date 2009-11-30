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


#include <ctype.h>
#include "sysstuff.h"
#include "macros.h"
#include "string2.h"
#include "smalloc.h"
#include "pbc.h"
#include "statutil.h"
#include "names.h"
#include "vec.h"
#include "futil.h"
#include "wman.h"
#include "tpxio.h"
#include "gmx_fatal.h"
#include "network.h"
#include "vec.h"
#include "mtop_util.h"
#include "gmxfio.h"
#include "history.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* used for npri */
#ifdef __sgi
#include <sys/schedctl.h>
#include <sys/sysmp.h>
#endif

/* The source code in this file should be thread-safe. 
      Please keep it that way. */

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* read only time names */
static const real timefactors[] =   { 0,  1e3,  1, 1e-3, 1e-6, 1e-9, 1e-12, 0 };
static const real timeinvfactors[] ={ 0, 1e-3,  1,  1e3,  1e6,  1e9,  1e12, 0 };
static const char *time_units_str[] = { NULL, "fs", "ps", "ns", "us", 
                                        "ms", "s", NULL };
static const char *time_units_xvgr[] = { NULL, "fs", "ps", "ns", 
                                         "\\8m\\4s", "ms", "s" };



struct output_env
{
    int time_unit; /*  the time unit as index for the above-defined arrays */
    bool view;  /* view of file requested */
    bool xvgr_codes; /* xmgrace-style legends etc. in file output */
#if 0
    const char *program; /* the program name */
    char *cmdline; /* the re-assembled command line */
#endif
};

/* inherently globally shared names: */
static const char *program_name=NULL;
static char *cmd_line=NULL;

#ifdef GMX_THREADS
/* For now, some things here are simply not re-entrant, so
   we have to actively lock them. */
static tMPI_Thread_mutex_t init_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


/****************************************************************
 *
 *            E X P O R T E D   F U N C T I O N S
 *
 ****************************************************************/

void init_output_env(output_env_t oenv,  int argc, char *argv[],
                     bool view, bool xvgr_codes, const char *timenm)
{
    set_output_env(oenv, view, xvgr_codes, timenm);
}

void set_output_env(output_env_t oenv,  bool view, bool xvgr_codes, 
                    const char *timenm)
{
    int i;
    int time_unit=2; /* the default is ps */
    
    if (timenm)
    {
        i=1;
        while(time_units_str[i])
        {
            if (strcmp(timenm, time_units_str[i])==0)
                break;
            i++;
        }
        if (time_units_str[i])
            time_unit=i;
    }
    oenv->time_unit=time_unit;
    oenv->view=view;
    oenv->xvgr_codes=xvgr_codes;
}


const char *ShortProgram(void)
{
    const char *pr,*ret;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&init_mutex);
#endif
    pr=ret=program_name; 
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&init_mutex);
#endif
    if ((pr=strrchr(ret,'/')) != NULL)
        ret=pr+1;
    /*else
        ret=ret;*/
    return ret;
}

const char *Program(void)
{
    const char *ret;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&init_mutex);
#endif
    ret=program_name; 
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&init_mutex);
#endif
    return ret;
}

const char *command_line(void)
{
    const char *ret;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&init_mutex);
#endif
    ret=cmd_line; 
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&init_mutex);
#endif
    return ret;
}

void set_program_name(const char *argvzero)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&init_mutex);
#endif
    /* When you run a dynamically linked program before installing
     * it, libtool uses wrapper scripts and prefixes the name with "lt-".
     * Until libtool is fixed to set argv[0] right, rip away the prefix:
     */
    if (program_name == NULL)
    {
        if(strlen(argvzero)>3 && !strncmp(argvzero,"lt-",3))
            program_name=strdup(argvzero+3);
        else
            program_name=strdup(argvzero);
    }
    if (program_name == NULL)
        program_name="GROMACS";
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&init_mutex);
#endif
}


void set_command_line(int argc, char *argv[])
{
    int i;
    size_t cmdlength;

#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&init_mutex);
#endif
    if (cmd_line==NULL)
    {
        cmdlength = strlen(argv[0]);
        for (i=1; i<argc; i++) 
        {
            cmdlength += strlen(argv[i]);
        }
        
        /* Fill the cmdline string */
        snew(cmd_line,cmdlength+argc+1);
        for (i=0; i<argc; i++) 
        {
            strcat(cmd_line,argv[i]);
            strcat(cmd_line," ");
        }
    }
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&init_mutex);
#endif

}


bool bRmod_fd(double a, double b, double c, bool bDouble)
{
    int iq;
    double tol;
    
    tol = 2*(bDouble ? GMX_DOUBLE_EPS : GMX_FLOAT_EPS);
    
    iq = (a - b + tol*a)/c;
    
    if (fabs(a - b - c*iq) <= tol*fabs(a))
        return TRUE;
    else
        return FALSE;
}

int check_times2(real t,real t0,real tp, real tpp, bool bDouble)
{
    int  r;
    real margin;
    
#ifndef GMX_DOUBLE
    /* since t is float, we can not use double precision for bRmod */
    bDouble = FALSE;
#endif
    
    if (t-tp>0 && tp-tpp>0)
        margin = 0.1*min(t-tp,tp-tpp);
    else
        margin = 0;
    
    r=-1;
    if ((!bTimeSet(TBEGIN) || (t >= rTimeValue(TBEGIN)))  &&
        (!bTimeSet(TEND)   || (t <= rTimeValue(TEND)))) {
        if (bTimeSet(TDELTA) && !bRmod_fd(t,t0,rTimeValue(TDELTA),bDouble))
            r = -1;
        else
            r = 0;
    }
    else if (bTimeSet(TEND) && (t >= rTimeValue(TEND)))
        r = 1;
    if (debug) 
        fprintf(debug,"t=%g, t0=%g, b=%g, e=%g, dt=%g: r=%d\n",
                t,t0,rTimeValue(TBEGIN),rTimeValue(TEND),rTimeValue(TDELTA),r);
    return r;
}

int check_times(real t)
{
    return check_times2(t,t,t,t,FALSE);
}



const char *get_time_unit(const output_env_t oenv)
{
    return time_units_str[oenv->time_unit];
}

const char *get_time_label(const output_env_t oenv)
{
    char *label;
    snew(label, 20);
    
    sprintf(label,"Time (%s)",time_units_str[oenv->time_unit] ? 
            time_units_str[oenv->time_unit]: "ps");
    
    return label;
}

const char *get_xvgr_tlabel(const output_env_t oenv)
{
    char *label;
    snew(label, 20);
    
    sprintf(label,"Time (%s)", time_units_xvgr[oenv->time_unit] ?
            time_units_xvgr[oenv->time_unit] : "ps");
    
    return label;
}


real get_time_factor(const output_env_t oenv)
{
    return timefactors[oenv->time_unit];
}

real get_time_invfactor(const output_env_t oenv)
{
    return timeinvfactors[oenv->time_unit];
}

real conv_time(const output_env_t oenv, real time)
{
    return time*timefactors[oenv->time_unit];
}


void conv_times(const output_env_t oenv, int n, real *time)
{
    int i;
    double fact=timefactors[oenv->time_unit];
    
    if (fact!=1.)
        for(i=0; i<n; i++)
            time[i] *= fact;
}

bool get_view(const output_env_t oenv)
{
    return oenv->view;
}

bool get_print_xvgr_codes(const output_env_t oenv)
{
    return oenv->xvgr_codes;
}



static void set_default_time_unit(const char *time_list[])
{
    int i,j;
    const char *select=getenv("GMXTIMEUNIT");
    if (!select)
        select="ps";
    
    i=1;
    while(time_list[i] && strcmp(time_list[i], select)!=0)
        i++;
    if (strcmp(time_list[i], select)==0) {
        /* swap the values and set time_list[0]*/
        time_list[0]=time_list[i];
        time_list[i]=time_list[1];
        time_list[1]=time_list[0];
    }
}

/***** T O P O L O G Y   S T U F F ******/

t_topology *read_top(const char *fn,int *ePBC)
{
    int        epbc,natoms;
    t_topology *top;
    
    snew(top,1);
    epbc = read_tpx_top(fn,NULL,NULL,&natoms,NULL,NULL,NULL,top);
    if (ePBC)
        *ePBC = epbc;
    
    return top;
}

/*************************************************************
 *
 *           P A R S I N G   S T U F F
 *
 *************************************************************/

static void usage(const char *type,const char *arg)
{
    if (arg != NULL)
        gmx_fatal(FARGS,"Expected %s argument for option %s\n",type,arg);
}

int iscan(int argc,char *argv[],int *i)
{
    int var;
    
    if (argc > (*i)+1) {
        if (!sscanf(argv[++(*i)],"%d",&var))
            usage("an integer",argv[(*i)-1]);
    } else
        usage("an integer",argv[*i]);
    
    return var;
}

gmx_large_int_t istepscan(int argc,char *argv[],int *i)
{
    gmx_large_int_t var;
    
    if (argc > (*i)+1) {
        if (!sscanf(argv[++(*i)],gmx_large_int_pfmt,&var))
            usage("an integer",argv[(*i)-1]);
    } else
        usage("an integer",argv[*i]);
    
    return var;
}

double dscan(int argc,char *argv[],int *i)
{
    double var;
    
    if (argc > (*i)+1) {
        if (!sscanf(argv[++(*i)],"%lf",&var))
            usage("a real",argv[(*i)-1]);
    } else
        usage("a real",argv[*i]);
    
    return var;
}

char *sscan(int argc,char *argv[],int *i)
{
    if (argc > (*i)+1) 
    {
        if ( (argv[(*i)+1][0]=='-') && (argc > (*i)+2) && 
           (argv[(*i)+2][0]!='-') )
        {
            fprintf(stderr,"Possible missing string argument for option %s\n\n",
                    argv[*i]);
        }
    } 
    else
        usage("a string",argv[*i]);
    
    return argv[++(*i)];
}

int nenum(const char *const enumc[])
{
    int i;
    
    i=1;
    /* we *can* compare pointers directly here! */
    while(enumc[i] && enumc[0]!=enumc[i])
        i++;
    
    return i;
}

static void pdesc(char *desc)
{
    char *ptr,*nptr;
    
    ptr=desc;
    if ((int)strlen(ptr) < 70)
        fprintf(stderr,"\t%s\n",ptr);
    else {
        for(nptr=ptr+70; (nptr != ptr) && (!isspace(*nptr)); nptr--)
            ;
        if (nptr == ptr)
            fprintf(stderr,"\t%s\n",ptr);
        else {
            *nptr='\0';
            nptr++;
            fprintf(stderr,"\t%s\n",ptr);
            pdesc(nptr);
        }
    }
}

static FILE *man_file(const output_env_t oenv,const char *mantp)
{
    FILE   *fp;
    char   buf[256];
    const char *pr;
    const char *program_name=Program();
    
    if ( (pr=strrchr(program_name,'/')) == NULL)
        pr=program_name;
    else 
        pr+=1;
    
    if (strcmp(mantp,"ascii") != 0)
        sprintf(buf,"%s.%s",pr,mantp);
    else
        sprintf(buf,"%s.txt",pr);
    fp = gmx_fio_fopen(buf,"w");
    
    return fp;
}

static int add_parg(int npargs,t_pargs *pa,t_pargs *pa_add)
{
    memcpy(&(pa[npargs]),pa_add,sizeof(*pa_add));
    
    return npargs+1;
}

static char *mk_desc(t_pargs *pa, const char *time_unit_str)
{
    char *newdesc=NULL,*ndesc=NULL,*nptr=NULL;
    const char*ptr=NULL;
    int  len,k;
    
    /* First compute length for description */
    len = strlen(pa->desc)+1;
    if ((ptr = strstr(pa->desc,"HIDDEN")) != NULL)
        len += 4;
    if (pa->type == etENUM) {
        len += 10;
        for(k=1; (pa->u.c[k] != NULL); k++) {
            len += strlen(pa->u.c[k])+12;
        }
    }
    snew(newdesc,len);
    
    /* add label for hidden options */
    if (is_hidden(pa)) 
        sprintf(newdesc,"[hidden] %s",ptr+6);
    else
        strcpy(newdesc,pa->desc);
    
    /* change '%t' into time_unit */
#define TUNITLABEL "%t"
#define NTUNIT strlen(TUNITLABEL)
    if (pa->type == etTIME)
        while( (nptr=strstr(newdesc,TUNITLABEL)) != NULL ) {
            nptr[0]='\0';
            nptr+=NTUNIT;
            len+=strlen(time_unit_str)-NTUNIT;
            snew(ndesc,len);
            strcpy(ndesc,newdesc);
            strcat(ndesc,time_unit_str);
            strcat(ndesc,nptr);
            sfree(newdesc);
            newdesc=ndesc;
            ndesc=NULL;
        }
#undef TUNITLABEL
#undef NTUNIT
    
    /* Add extra comment for enumerateds */
    if (pa->type == etENUM) {
        strcat(newdesc,": ");
        for(k=1; (pa->u.c[k] != NULL); k++) {
            strcat(newdesc,"[TT]");
            strcat(newdesc,pa->u.c[k]);
            strcat(newdesc,"[tt]");
            /* Print a comma everywhere but at the last one */
            if (pa->u.c[k+1] != NULL) {
                if (pa->u.c[k+2] == NULL)
                    strcat(newdesc," or ");
                else
                    strcat(newdesc,", ");
            }
        }
    }
    return newdesc;
}


void parse_common_args(int *argc,char *argv[],unsigned long Flags,
		       int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
		       int ndesc,const char **desc,
		       int nbugs,const char **bugs,
                       output_env_t *oenv)
{
    bool bHelp=FALSE,bHidden=FALSE,bQuiet=FALSE;
    const char *manstr[] = { NULL, "no", "html", "tex", "nroff", "ascii", 
                            "completion", "py", "xml", "wiki", NULL };
    const char *time_units[] = { NULL, "ps", "fs", "ns", "us", "ms", "s", 
                                NULL };
    int  nicelevel=0,mantp=0,npri=0,debug_level=0;
    char *deffnm=NULL;
    real tbegin=0,tend=0,tdelta=0;
    bool bView=FALSE, bXvgrCodes=TRUE;
    
    t_pargs *all_pa=NULL;
    
    t_pargs npri_pa   = { "-npri", FALSE, etINT,   {&npri},
    "HIDDEN Set non blocking priority (try 128)" };
    t_pargs nice_pa   = { "-nice", FALSE, etINT,   {&nicelevel}, 
    "Set the nicelevel" };
    t_pargs deffnm_pa = { "-deffnm", FALSE, etSTR, {&deffnm}, 
    "Set the default filename for all file options" };
    t_pargs begin_pa  = { "-b",    FALSE, etTIME,  {&tbegin},        
    "First frame (%t) to read from trajectory" };
    t_pargs end_pa    = { "-e",    FALSE, etTIME,  {&tend},        
    "Last frame (%t) to read from trajectory" };
    t_pargs dt_pa     = { "-dt",   FALSE, etTIME,  {&tdelta},        
    "Only use frame when t MOD dt = first time (%t)" };
    t_pargs view_pa   = { "-w",    FALSE, etBOOL,  {&bView},
    "View output xvg, xpm, eps and pdb files" };
    t_pargs code_pa   = { "-xvgr", FALSE, etBOOL,  {&bXvgrCodes},
    "Add specific codes (legends etc.) in the output xvg files for the xmgrace program" };
    t_pargs time_pa   = { "-tu",   FALSE, etENUM,  {time_units},
    "Time unit" };
    /* Maximum number of extra arguments */
#define EXTRA_PA 16
    
    t_pargs pca_pa[] = {
        { "-h",    FALSE, etBOOL, {&bHelp},     
        "Print help info and quit" }, 
        { "-hidden", FALSE, etBOOL, {&bHidden},
        "HIDDENPrint hidden options" },
        { "-quiet",FALSE, etBOOL, {&bQuiet},
        "HIDDENDo not print help info" },
        { "-man",  FALSE, etENUM,  {manstr},
        "HIDDENWrite manual and quit" },
        { "-debug",FALSE, etINT, {&debug_level},
        "HIDDENWrite file with debug information, 1: short, 2: also x and f" },
    };
#define NPCA_PA asize(pca_pa)
    FILE *fp;  
    bool bPrint,bExit,bXvgr;
    int  i,j,k,npall,max_pa,cmdlength;
    char *ptr,*newdesc;
    const char *envstr;
    
#define FF(arg) ((Flags & arg)==arg)

    init_history(*argc,argv,nfile,fnm);
    
    snew(*oenv, 1);
    
    cmdlength = strlen(argv[0]);
    /* Check for double arguments */
    for (i=1; (i<*argc); i++) 
    {
        cmdlength += strlen(argv[i]);
        if (argv[i] && (strlen(argv[i]) > 1) && (!isdigit(argv[i][1]))) 
        {
            for (j=i+1; (j<*argc); j++) 
            {
                if ( (argv[i][0]=='-') && (argv[j][0]=='-') && 
                    (strcmp(argv[i],argv[j])==0) ) 
                {
                    if (FF(PCA_NOEXIT_ON_ARGS))
                        fprintf(stderr,"Double command line argument %s\n",
                                argv[i]);
                    else
                        gmx_fatal(FARGS,"Double command line argument %s\n",
                                  argv[i]);
                }
            }
        }
    }
    debug_gmx();
    set_program_name(argv[0]);
    set_command_line(*argc, argv);
    
    /* Handle the flags argument, which is a bit field 
     * The FF macro returns whether or not the bit is set
     */
    bPrint        = !FF(PCA_SILENT);
    
    /* Check ALL the flags ... */
    max_pa = NPCA_PA + EXTRA_PA + npargs+1;
    snew(all_pa,max_pa);
    
    for(i=npall=0; (i<NPCA_PA); i++)
        npall = add_parg(npall,all_pa,&(pca_pa[i]));
    
#ifdef __sgi
    envstr = getenv("GMXNPRIALL");
    if (envstr)
        npri=strtol(envstr,NULL,10);
    if (FF(PCA_BE_NICE)) {
        envstr = getenv("GMXNPRI");
        if (envstr)
            npri=strtol(envstr,NULL,10);
    }
    npall = add_parg(npall,all_pa,&npri_pa);
#endif
    
    if (FF(PCA_BE_NICE)) 
        nicelevel=19;
    npall = add_parg(npall,all_pa,&nice_pa);
    
    if (FF(PCA_CAN_SET_DEFFNM)) 
        npall = add_parg(npall,all_pa,&deffnm_pa);   
    if (FF(PCA_CAN_BEGIN)) 
        npall = add_parg(npall,all_pa,&begin_pa);
    if (FF(PCA_CAN_END))
        npall = add_parg(npall,all_pa,&end_pa);
    if (FF(PCA_CAN_DT))
        npall = add_parg(npall,all_pa,&dt_pa);
    if (FF(PCA_TIME_UNIT)) {
        set_default_time_unit(time_units);
        npall = add_parg(npall,all_pa,&time_pa);
    } 
    if (FF(PCA_CAN_VIEW)) 
        npall = add_parg(npall,all_pa,&view_pa);
    
    bXvgr = FALSE;
    for(i=0; (i<nfile); i++)
        bXvgr = bXvgr ||  (fnm[i].ftp == efXVG);
    if (bXvgr)
        npall = add_parg(npall,all_pa,&code_pa);
    
    /* Now append the program specific arguments */
    for(i=0; (i<npargs); i++)
        npall = add_parg(npall,all_pa,&(pa[i]));
    
    /* set etENUM options to default */
    for(i=0; (i<npall); i++)
        if (all_pa[i].type==etENUM)
            all_pa[i].u.c[0]=all_pa[i].u.c[1];
    
    
    /* set program name, command line, and default values for output options */
    init_output_env(*oenv, *argc, argv, bView, bXvgrCodes, time_units[1]);
   
    /* Now parse all the command-line options */
    get_pargs(argc,argv,npall,all_pa,FF(PCA_KEEP_ARGS));
    
    if (FF(PCA_CAN_SET_DEFFNM) && (deffnm!=NULL))
        set_default_file_name(deffnm);
    
    /* Parse the file args */
    parse_file_args(argc,argv,nfile,fnm,FF(PCA_KEEP_ARGS),!FF(PCA_NOT_READ_NODE));
    
    /* Open the debug file */
    if (debug_level > 0) {
        char buf[256];
        
        if (gmx_mpi_initialized())
            sprintf(buf,"%s%d.debug",ShortProgram(),gmx_node_rank());
        else
            sprintf(buf,"%s.debug",ShortProgram());
        
        init_debug(debug_level,buf);
        fprintf(stderr,"Opening debug file %s (src code file %s, line %d)\n",
                buf,__FILE__,__LINE__);
    }
    
    /* Now copy the results back... */
    for(i=0,k=npall-npargs; (i<npargs); i++,k++) 
        memcpy(&(pa[i]),&(all_pa[k]),(size_t)sizeof(pa[i]));


    for(i=0; (i<npall); i++)
        all_pa[i].desc = mk_desc(&(all_pa[i]), get_time_unit(*oenv) );
   
    bExit = bHelp || (strcmp(manstr[0],"no") != 0);
    
#if (defined __sgi && USE_SGI_FPE)
    doexceptions();
#endif
    
    /* Set the nice level */
#ifdef __sgi
    if (npri != 0 && !bExit) {
        schedctl(MPTS_RTPRI,0,npri);
    }
#endif 
    
#ifdef HAVE_UNISTD_H
    
#ifndef GMX_NO_NICE
    /* The some system, e.g. the catamount kernel on cray xt3 do not have nice(2). */
    if (nicelevel != 0 && !bExit)
    {
#ifdef GMX_THREADS
        static bool nice_set=FALSE; /* only set it once */
        tMPI_Thread_mutex_lock(&init_mutex);
        if (!nice_set)
        {
#endif
            i=nice(nicelevel); /* assign ret value to avoid warnings */
#ifdef GMX_THREADS
            nice_set=TRUE;
        }
        tMPI_Thread_mutex_unlock(&init_mutex);
#endif
    }
#endif
#endif
    
    set_output_env(*oenv, bView, bXvgrCodes, time_units[0]);
    /*init_time_factor(tms);*/
   
    
    if (!(FF(PCA_QUIET) || bQuiet )) {
        if (bHelp)
            write_man(stderr,"help",Program(),ndesc,desc,nfile,
                      fnm,npall,all_pa, nbugs,bugs,bHidden);
        else if (bPrint) {
            pr_fns(stderr,nfile,fnm);
            print_pargs(stderr,npall,all_pa,FALSE);
        }
    }
    
    if (strcmp(manstr[0],"no") != 0) {
        if(!strcmp(manstr[0],"completion")) {
            /* one file each for csh, bash and zsh if we do completions */
            fp=man_file(*oenv,"completion-zsh");
            write_man(fp,"completion-zsh",Program(),ndesc,desc,nfile,
                      fnm, npall,all_pa,nbugs,bugs,bHidden);
            gmx_fio_fclose(fp);
            fp=man_file(*oenv,"completion-bash");
            write_man(fp,"completion-bash",Program(),ndesc,desc,nfile,
                      fnm, npall,all_pa,nbugs,bugs,bHidden);
            gmx_fio_fclose(fp);
            fp=man_file(*oenv,"completion-csh");
            write_man(fp,"completion-csh",Program(),ndesc,desc,nfile,
                      fnm, npall,all_pa,nbugs,bugs,bHidden);
            gmx_fio_fclose(fp);
        } else {
            fp=man_file(*oenv,manstr[0]);
            write_man(fp,manstr[0],Program(),ndesc,desc,nfile,fnm,
                      npall, all_pa,nbugs,bugs,bHidden);
            gmx_fio_fclose(fp);
        }
    }
    
    /* convert time options, must be done after printing! */
    
    for(i=0; i<npall; i++) {
        if ((all_pa[i].type == etTIME) && (*all_pa[i].u.r >= 0)) {
            *all_pa[i].u.r *= get_time_invfactor(*oenv);
        }
    }
    
    /* Extract Time info from arguments */
    if (FF(PCA_CAN_BEGIN) && opt2parg_bSet("-b",npall,all_pa))
        setTimeValue(TBEGIN,opt2parg_real("-b",npall,all_pa));
    
    if (FF(PCA_CAN_END) && opt2parg_bSet("-e",npall,all_pa))
        setTimeValue(TEND,opt2parg_real("-e",npall,all_pa));
    
    if (FF(PCA_CAN_DT) && opt2parg_bSet("-dt",npall,all_pa))
        setTimeValue(TDELTA,opt2parg_real("-dt",npall,all_pa));
    
    /* clear memory */
    sfree(all_pa);
    
    if (!FF(PCA_NOEXIT_ON_ARGS)) {
        if (*argc > 1) {
            gmx_cmd(argv[1]);
        }
    } 
    if (bExit) {
        if (gmx_parallel_env_initialized())
            /*gmx_abort(gmx_node_rank(),gmx_node_num(),0);*/
            gmx_finalize();
        exit(0);
    }
#undef FF
}

