/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_statutil_c = "$Id$";

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
#include "assert.h"
#include "fatal.h"

/* used for npri */
#ifdef _SGI_
#include <sys/schedctl.h>
#include <sys/sysmp.h>
#endif

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* Globals for trajectory input */
real         tbegin=-1.0,tend=-1.0;
static bool  bView=FALSE;
static ulong uFlags=0;
static char  *program=NULL;

#define FF(arg) ((uFlags & arg)==arg)

char *ShortProgram(void)
{
  char *pr;
  
  if (program) {
    if ((pr=strrchr(program,'/')) != NULL)
      return pr+1;
    else
      return program;
  }
  else
    return "GROMACS";
}

char *Program(void)
{
  if (program)
    return program;
  else
    return "GROMACS";
}

/****************************************************************
 *
 *            E X P O R T E D   F U N C T I O N S
 *
 ****************************************************************/

int check_times(real t) 
{
  if ((((tbegin >= 0.0) && (t >= tbegin)) || (tbegin == -1.0)) &&
      (((tend   >= 0.0) && (t <= tend))   || (tend   == -1.0))) {
    return 0;
  }
  else if ((tend != -1.0) && (t>=tend))
    return 1;
  return -1;
}

#ifdef NOTOBSOLETE

/**** C O O R D I N A T E S   A N D   V  E L O C  I T I E S ****/
int read_first_x_v(FILE *status,real *t,rvec **x,rvec **v,matrix box)
{
  t_statheader sh;
  real lambda;
  int  step,nre;

  fprintf(stderr,"\nReading statusfile, version: %s\n",rd_header(status,&sh));
  snew(*v,sh.natoms);
  snew(*x,sh.natoms);
  rd_hstatus(status,&sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
	     &sh.natoms,*x,*v,NULL,
	     &nre,NULL,NULL);

  return sh.natoms;
}

bool read_next_x_v(FILE *status,real *t,int natoms,rvec x[],rvec v[],matrix box)
{
  t_statheader sh;
  real lambda;
  int  step,nre;
  bool bV;
  bool bX;

  while (!eof(status)) {
    rd_header(status,&sh);
    bX=sh.x_size;
    bV=sh.v_size;
    rd_hstatus(status,&sh,&step,t,&lambda,NULL,NULL,NULL,NULL,
	       &sh.natoms,bX ? x : NULL,bV ? v : NULL,NULL,
	       &nre,NULL,NULL);
    if ((check_times(*t)==0) && (bV) && (bX))
      return TRUE;
    if (check_times(*t) > 0)
      return FALSE;
  }
  return FALSE;
}

#endif

/***** T O P O L O G Y   S T U F F ******/

t_topology *read_top(char *fn)
{
  int        step,natoms;
  real       t,lambda;
  t_topology *top;

  snew(top,1);
  read_tpx(fn,&step,&t,&lambda,NULL,NULL,
	   &natoms,NULL,NULL,NULL,top);
  
  return top;
}

void mk_single_top(t_topology *top)
{
  int i;

  for(i=0; (i<ebNR); i++)
    top->blocks[i].multinr[0]=top->blocks[i].multinr[MAXPROC-1];
  for(i=0; (i<F_NRE); i++)
    top->idef.il[i].multinr[0]=top->idef.il[i].multinr[MAXPROC-1];
}

/*************************************************************
 *
 *           P A R S I N G   S T U F F
 *
 *************************************************************/
 
int iscan(int argc,char *argv[],int *i)
{
  int var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%d",&var))
      usage(argv[0],argv[*i]);
  }
  else
    usage(argv[0],argv[*i]);

  return var;
}

double dscan(int argc,char *argv[],int *i)
{
  double var;

  if (argc > (*i)+1) {
    if (!sscanf(argv[++(*i)],"%lf",&var))
      usage(argv[0],argv[*i]);
  }
  else
    usage(argv[0],argv[*i]);

  return var;
}

char *sscan(int argc,char *argv[],int *i)
{
  /* 
     if (strlen(argv[*i]) > 2) 
     return &(argv[*i][2]);
     else
     */
  if (argc > (*i)+1) 
    return argv[++(*i)];
  else
    usage(argv[0],argv[*i]);

  return NULL; /* For the compiler! */
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

void usage(char *prog,char *arg)
{
  if (arg != NULL)
    fprintf(stderr,"Conflicting argument for program %s: %s\n",prog,arg);
  exit(1);
}

bool bDoView(void)
{
  return bView;
}

static FILE *man_file(char *program,char *mantp)
{
  FILE   *fp;
  char   buf[256],*pr;

  if ((pr=strrchr(program,'/')) == NULL)
    pr=program;
  else 
    pr+=1;
    
  if (strcmp(mantp,"ascii") != 0)
    sprintf(buf,"%s.%s",pr,mantp);
  else
    sprintf(buf,"%s.txt",pr);
  fp = ffopen(buf,"w");
  
  return fp;
}

static int add_parg(int npargs,t_pargs **pa,t_pargs *pa_add)
{
  srenew((*pa),npargs+1);
  memcpy(&((*pa)[npargs]),pa_add,sizeof(*pa_add));
  
  return npargs+1;
}

static char *mk_desc(t_pargs *pa)
{
  char *newdesc,*ptr;
  int  len,k;
  
  /* First compute length for description */
  len = strlen(pa->desc)+1;
  if ((ptr = strstr(pa->desc,"HIDDEN")) != NULL)
    len += 4;
  if (pa->type == etENUM) {
    len += 10;
    for(k=1; (pa->u.c[k] != NULL); k++) {
      len += strlen(pa->u.c[k])+4;
    }
  }
  snew(newdesc,len);
  
  if (is_hidden(pa)) 
    sprintf(newdesc,"[hidden] %s",ptr+6);
  else
    strcpy(newdesc,pa->desc);
    
  /* Add extra comment for enumerateds */
  if (pa->type == etENUM) {
    strcat(newdesc,": ");
    for(k=1; (pa->u.c[k] != NULL); k++) {
      strcat(newdesc,pa->u.c[k]);
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

void parse_common_args(int *argc,char *argv[],ulong Flags,bool bNice,
		       int nfile,t_filenm fnm[],int npargs,t_pargs *pa,
		       int ndesc,char **desc,int nbugs,char **bugs)
{
  static bool bHelp=FALSE,bHidden=FALSE,bQuiet=FALSE;
  static char *manstr[]      = { NULL, "no", "html", "tex", "java", "ascii", NULL };
  static char *not_nicestr[] = { NULL, "0", "4", "10", "19", NULL };
  static char *nicestr[]     = { NULL, "19", "10", "4", "0", NULL };
  static char *not_npristr[] = { NULL, "0", "100", "200", "250", NULL };
  static char *npristr[]     = { NULL, "250", "200", "100", "0", NULL };
  static int  nicelevel=0,mantp=0,npri=0;
  static bool bExcept=FALSE,bGUI=FALSE,bDebug=FALSE;
  
  FILE *fp;  
  bool bPrint;
  int  i,j,k,npall;
  char *ptr,*newdesc;
  
  t_pargs *all_pa=NULL;
  
  t_pargs motif_pa = { "-X",    FALSE, etBOOL, &bGUI,
		       "Use dialog box GUI to edit command line options" };
  t_pargs fpe_pa   = { "-exception", FALSE, etBOOL, &bExcept,
		       "HIDDENTurn on exception handling" };
  t_pargs npri_paX = { "-npri", FALSE, etENUM,  not_npristr,
		       "Set non blocking priority" };
  t_pargs npri_pa  = { "-npri", FALSE, etINT,  &npri,
		       "Set non blocking priority (try 250)" };
  t_pargs nice_paX = { "-nice", FALSE, etENUM, not_nicestr, 
		       "Set the nicelevel" };
  t_pargs nice_pa  = { "-nice", FALSE, etINT,  &nicelevel, 
		       "Set the nicelevel" };
  t_pargs begin_pa = { "-b",    FALSE, etREAL, &tbegin,        
		       "First frame (ps) to read from trajectory" };
  t_pargs end_pa   = { "-e",    FALSE, etREAL, &tend,        
		       "Last frame (ps) to read from trajectory" };
  t_pargs view_pa  = { "-w",    FALSE, etBOOL, &bView,     
		       "View output using xvgr or ghostview" };
  
  t_pargs pca_pa[] = {
    { "-h",    FALSE, etBOOL, &bHelp,     
      "Print help info and quit" }, 
    { "-hidden", FALSE, etBOOL, &bHidden,
      "HIDDENPrint hidden options" },
    { "-quiet",FALSE, etBOOL, &bQuiet,
      "HIDDENDo not print help info" },
    { "-man",  FALSE, etENUM,  manstr,
      "HIDDENWrite manual and quit" },
    { "-debug",FALSE, etBOOL, &bDebug,
      "HIDDENWrite file with debug information" },
  };
#define NPCA_PA asize(pca_pa)

  /* Check for double arguments */
  for (i=1; (i<*argc); i++)
    if ( !isdigit(argv[i][1]) )
      for (j=i+1; (j<*argc); j++)
	if ( (argv[i][0]=='-') && (argv[j][0]=='-') && 
	(strcmp(argv[i],argv[j])==0) )
	  fatal_error(0,"Double command line argument %s\n",argv[i]);
  
  /* Handle the flags argument, which is a bit field 
   * The FF macro returns whether or not the bit is set
   */
  uFlags        = Flags;
  bPrint        = !FF(PCA_SILENT);
  
  /* Check whether we should have GUI or not */
#ifdef HAVE_MOTIF
  bGUI = (getenv("GMXMOTIF") != NULL);
  for(i=1; (i<*argc); i++) {
    if (strcmp(argv[i],"-X") == 0)
      bGUI = TRUE;
    else if (strcmp(argv[i],"-noX") == 0)
      bGUI = FALSE;
  }
  if (bGUI)
    bQuiet = TRUE;
#endif
  
  /* First do file stuff */
  if (!program)
    program  = strdup(argv[0]);

  /* Parse the file args */
  parse_file_args(argc,argv,nfile,fnm,FF(PCA_KEEP_ARGS));
  
  /* Check ALL the flags ... */
  snew(all_pa,NPCA_PA+npargs);
  for(i=npall=0; (i<NPCA_PA); i++)
    npall = add_parg(npall,&(all_pa),&(pca_pa[i]));

  /* Motif options */
#ifdef HAVE_MOTIF
  npall = add_parg(npall,&(all_pa),&motif_pa);
#endif

#ifdef _SGI_
#ifdef USE_SGI_FPE
  npall = add_parg(npall,&(all_pa),&fpe_pa);
#endif
#ifndef NO_NICE
  if (FF(PCA_SET_NPRI))
    npri=atoi(getenv("GMXNPRI"));
  if (bGUI) {
    if (npri)
      npri_paX.u.c = npristr;
    npall = add_parg(npall,&(all_pa),&npri_paX);
  }
  else
    npall = add_parg(npall,&(all_pa),&npri_pa);
#endif
#endif

#ifndef NO_NICE
  if (bGUI) {
    /* Automatic nice or scheduling options */
    if (bNice) 
      nice_paX.u.c = nicestr;
    npall = add_parg(npall,&(all_pa),&nice_paX);
  }
  else {
    if (bNice) 
      nicelevel=19;
    npall = add_parg(npall,&(all_pa),&nice_pa);
  }
#endif

  if (FF(PCA_CAN_BEGIN)) 
    npall = add_parg(npall,&(all_pa),&begin_pa);
  if (FF(PCA_CAN_END))
    npall = add_parg(npall,&(all_pa),&end_pa);
  if (FF(PCA_CAN_VIEW))
    npall = add_parg(npall,&(all_pa),&view_pa);

  /* Now append the program specific arguments */
  for(i=0; (i<npargs); i++)
    npall = add_parg(npall,&(all_pa),&(pa[i]));

  /* set etENUM options to default */
  for(i=0; (i<npall); i++)
    if (all_pa[i].type==etENUM)
      all_pa[i].u.c[0]=all_pa[i].u.c[1];
  
  /* Now parse all the command-line options */
  get_pargs(argc,argv,npall,all_pa,FF(PCA_KEEP_ARGS));

  /* Open the debug file */
  if (bDebug) {
    char buf[256];
    
    sprintf(buf,"%s.log",ShortProgram());
    init_debug(buf);
    fprintf(debug,"%s (this file) opened in file %s, line %d\n",
	    buf,__FILE__,__LINE__);
  }

#ifdef HAVE_MOTIF
  /* Now we have parsed the command line arguments. If the user wants it
   * we can now plop up a GUI dialog box to edit options.
   */
  if (bGUI) {
    gmx_gui(argc,argv,nfile,fnm,npall,all_pa,ndesc,desc,nbugs,bugs);
  }
#endif

  /* Now copy the results back... */
  for(i=0,k=npall-npargs; (i<npargs); i++,k++) 
    memcpy(&(pa[i]),&(all_pa[k]),(size_t)sizeof(pa[i]));
  
  for(i=0; (i<npall); i++) {
    all_pa[i].desc = mk_desc(&(all_pa[i]));
  }
    
#ifdef _SGI_
#ifdef USE_SGI_FPE
  /* Install exception handler if necessary */
  if (bExcept)
    doexceptions();
#endif
#endif
  
  /* Set the nice level */
#ifdef _SGI_
  if (bGUI)
    if (npri)
      sscanf(npristr[0],"%d",&npri);
    else
      sscanf(not_npristr[0],"%d",&npri);
  if (npri != 0) {
    (void) schedctl(MPTS_RTPRI,0,npri);
  }
  else
#endif
#ifndef NO_NICE
    if (bGUI) {
      if (bNice)
	sscanf(nicestr[0],"%d",&nicelevel);
      else
	sscanf(not_nicestr[0],"%d",&nicelevel);
    }
  if (nicelevel != 0)
    nice(nicelevel);
  
#endif

  if (!(FF(PCA_QUIET) || bQuiet )) {
    if (bHelp)
      write_man(stdout,"help",program,ndesc,desc,nfile,fnm,npall,all_pa,
		nbugs,bugs,bHidden);
    else if (bPrint) {
      pr_fns(stdout,nfile,fnm);
      print_pargs(stdout,npall,all_pa);
    }
  }

  if (strcmp(manstr[0],"no") != 0) {
    fp=man_file(program,manstr[0]);
    write_man(fp,manstr[0],program,ndesc,desc,nfile,fnm,npall,all_pa,
	      nbugs,bugs,bHidden);
    fclose(fp);
    exit(0);
  }
  for(i=0; i<npall; i++)
    sfree(all_pa[i].desc);
  sfree(all_pa);

  if (!FF(PCA_NOEXIT_ON_ARGS)) {
    if (*argc > 1) {
      for(i=1; (i<*argc); i++) 
	fprintf(stderr,"Unknown argument: %s\n",argv[i]);
      fprintf(stderr,"\nProgram %s halted\n",Program());
      exit(1);
    }
  } 
  if (bHelp)
    exit(0);
}
