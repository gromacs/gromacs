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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_filenm_c = "$Id$";

#include <string.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "string2.h"
#include "fatal.h"
#include "filenm.h"
#include "futil.h"
#include "wman.h"
#include "macros.h"

/* Use bitflag ... */
#define is_set(fn) ((fn.flag & ffSET) != 0)
#define is_opt(fn) ((fn.flag & ffOPT) != 0)
#define un_set(fn) (fn.flag = (fn.flag & ~ffSET))
#define do_set(fn) (fn.flag = (fn.flag |  ffSET))

enum { eftASC, eftBIN, eftXDR, eftGEN, eftNR };

/* To support multiple file types with one general (eg TRX) we have 
 * these arrays.
 */
static    int trxs[]={
#ifdef USE_XDR 
  efXTC, efTRR, 
#endif
  efTRJ, efGRO, efPDB, efG87 };
#define NTRXS asize(trxs)

static    int trns[]={ 
#ifdef USE_XDR
  efTRR, 
#endif
  efTRJ };
#define NTRNS asize(trns)

static    int stos[]={ efGRO, efPDB, efBRK, efENT};
#define NSTOS asize(stos)

static    int stxs[]={ efGRO, efPDB, efBRK, efENT,
#ifdef USE_XDR 
		       efTPR, 
#endif 
		       efTPB, efTPA };
#define NSTXS asize(stxs)

static    int enxs[]={ 
#ifdef USE_XDR
  efEDR, 
#endif
  efENE };
#define NENXS asize(enxs)

static    int tpxs[]={ 
#ifdef USE_XDR
  efTPR, 
#endif
  efTPB, efTPA };
#define NTPXS asize(tpxs)

typedef struct {
  int  ftype;
  char *ext;
  char *defnm;
  char *defopt;
  char *descr;
} t_deffile;

/* this array should correspond to the enum in include/types/filenm.h */
static t_deffile deffile[efNR] = {
  { eftASC, ".mdp", "grompp", "-f", "grompp input file with MD parameters"   },
  { eftASC, ".gip", "genion", "-f", "genion input file with parameters etc." },
  { eftASC, ".gct", "gct",    "-f", "general coupling stuff"                 },
  { eftASC, ".gpp", "group",  "-g", "group related input for grompp"         },
  { eftGEN, ".???", "traj",   "-f", "Generic trajectory: xtc trr trj gro pdb g87"},
  { eftGEN, ".???", "traj",   NULL, "Full precision trajectory: trr trj"     },
  { eftXDR, ".trr", "traj",   NULL, "Trajectory in portable xdr format"      },
  { eftBIN, ".trj", "traj",   NULL, "Trajectory file (cpu specific)"         },
  { eftXDR, ".xtc", "traj",   NULL, "Compressed trajectory (portable xdr format)"},
  { eftASC, ".g87", "gtraj",  NULL, "Gromos-87 ASCII trajectory format"      },
  { eftGEN, ".???", "ener",   NULL, "Generic energy: edr ene"                },
  { eftXDR, ".edr", "ener",   NULL, "Energy file in portable XDR format"     },
  { eftBIN, ".ene", "ener",   NULL, "Energy file"                            },
  { eftGEN, ".???", "conf",   "-c", "Generic structure: gro pdb tpr tpb tpa" },
  { eftGEN, ".???", "out",    "-o", "Generic structure: gro pdb"             },
  { eftASC, ".gro", "conf",   "-c", "Coordinate file in Gromos-87 format"    },
  { eftASC, ".pdb", "eiwit",  "-f", "Protein data bank file"                 },
  { eftASC, ".brk", "eiwit",  "-f", "Brookhaven data bank file"              },
  { eftASC, ".ent", "eiwit",  "-f", "Entry in the protein date bank"         },
  { eftASC, ".log", "run",    "-l", "Log file from MD/LD/EM/NM run"          },
  { eftASC, ".xvg", "graph",  "-o", "xvgr file as produced by analysis tools"},
  { eftASC, ".out", "hello",  "-o", "Generic output file"                    },
  { eftASC, ".ndx", "index",  "-n", "Index file",                            },
  { eftASC, ".top", "topol",  "-p", "Topology file"                          },
  { eftASC, ".itp", "topinc", NULL, "Include file for topology"              },
  { eftGEN, ".???", "topol",  "-s", "Generic run input: tpr tpb tpa"         },
  { eftXDR, ".tpr", "topol",  "-s", "portable xdr run input file"            },
  { eftASC, ".tpa", "topol",  "-s", "Ascii run input file"                   },
  { eftBIN, ".tpb", "topol",  "-s", "Binary run input file"                  },
  { eftASC, ".tex", "doc",    "-o", "LaTeX file"                             },
  { eftASC, ".rtp", "residue",NULL, "Residue Type file used by pdb2gmx"      },
  { eftASC, ".atp", "atomtp", NULL, "Atomtype file used by pdb2gmx"          },
  { eftASC, ".hdb", "polar",  NULL, "Hydrogen data base"                     },
  { eftASC, ".dat", "nnnice", NULL, "Generic data file"                      },
  { eftASC, ".dlg", "user",   NULL, "Dialog Box data for ngmx"               },
  { eftASC, ".map", "ss",     NULL, "File that maps matrix data to colors"   },
  { eftASC, ".eps", "plot",   NULL, "Encapsulated PostScript (tm) file"      },
  { eftASC, ".mat", "ss",     NULL, "Matrix Data file"			     },
  { eftASC, ".m2p", "ps",     NULL, "Input file for mat2ps"                  },
  { eftASC, ".vdw", "radii",  NULL, "Database containing Van der Waals radii"},
  { eftBIN, ".mtx", "hessian","-m", "Hessian matrix"                         },
  { eftASC, ".vec", "eigvec", NULL, "(Eigen)Vectors"                         },
  { eftASC, ".edi", "sam",    NULL, "ED sampling input"                      },
  { eftASC, ".edo", "sam",    NULL, "ED sampling output"                     },
  { eftASC, ".hat", "gk",     NULL, "Fourier transform of spread function"   },
  { eftASC, ".xpm", "root",   NULL, "X PixMap compatible matrix file"        }
};

char *ftp2ext(int ftp)
{
  if ((0 <= ftp) && (ftp < efNR))
    return deffile[ftp].ext+1;
  else
    return NULL;
}

char *ftp2desc(int ftp)
{
  if ((0 <= ftp) && (ftp < efNR))
    return deffile[ftp].descr;
  else
    return NULL;
}

char *ftp2ftype(int ftp)
{
  if ((ftp >= 0) && (ftp < efNR)) {
    switch (deffile[ftp].ftype) {
    case eftASC: return "ASCII";
    case eftBIN: return "Binary";
    case eftXDR: return "XDR portable";
    case eftGEN: return "";
    default: fatal_error(0,"Some error");
      break;
    }
  }
  return NULL;
}

char *ftp2defnm(int ftp)
{
  static char buf[256];
  
  if ((0 <= ftp) && (ftp < efNR)) {
    sprintf(buf,"%s",deffile[ftp].defnm);
    return buf;
  }
  else
    return NULL;
}

void pr_def(FILE *fp,int ftp)
{
  t_deffile *df;
  char c;
  
  df=&(deffile[ftp]);
  switch (df->ftype) {
  case eftASC: c='A';
    break;
  case eftBIN: c='B';
    break;
  case eftXDR: c='P';
    break;
  case eftGEN: c=' ';
    break;
  default: fatal_error(0,"Some error, some ints %d %d",ftp,efNR);
    break;
  }
  fprintf(fp,"%8s & %5s & %c & %5s & ",
	  df->defnm,df->ext,c,df->defopt ? df->defopt : "");
  fprintf(fp,"%s \\\\\n",check_tex(df->descr));
}

void pr_fn(FILE *fp,t_filenm *tfn)
{
}

void pr_fns(FILE *fp,int nf,t_filenm tfn[])
{
  int  i,j;
  char buf[256],*wbuf;
#define OPTLEN 4
#define NAMELEN 14
  fprintf(fp,"%6s %12s  %-12s  %s\n",
	  "Option","Filename","Type","Description");
  fprintf(fp,"------------------------------------------------------------\n");
  for(i=0; (i<nf); i++) {
    sprintf(buf, "%4s %14s  %-12s  %s\n", tfn[i].opt,tfn[i].fn,
	    fileopt(tfn[i].flag),deffile[tfn[i].ftp].descr);
    if ( (strlen(tfn[i].opt)>OPTLEN) && 
	 (strlen(tfn[i].opt)<=((OPTLEN+NAMELEN)-strlen(tfn[i].fn))) ) {
      for(j=strlen(tfn[i].opt); 
	  j<strlen(buf)-(strlen(tfn[i].opt)-OPTLEN)+1; j++)
	buf[j]=buf[j+strlen(tfn[i].opt)-OPTLEN];
    }
    wbuf=wrap_lines(buf,80,35);
    fprintf(fp,wbuf);
    sfree(wbuf);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

static void check_opts(int nf,t_filenm fnm[])
{
  int       i;
  t_deffile *df;
  
  for(i=0; (i<nf); i++) {
    df=&(deffile[fnm[i].ftp]);
    if (fnm[i].opt == NULL) {
      if (df->defopt == NULL)
	fatal_error(0,"No default cmd-line option for %s (type %d)\n",
		    deffile[fnm[i].ftp].ext,fnm[i].ftp);
      else
	fnm[i].opt=df->defopt;
    }
  }
}

int fn2ftp(char *fn)
{
  int  i,len;
  char *feptr,*eptr;
  
  if (!fn)
    return efNR;

  len=strlen(fn);
  if ((len >= 4) && (fn[len-4] == '.'))
    feptr=&(fn[len-4]);
  else
    return efNR;
  
  for(i=0; (i<efNR); i++)
    if ((eptr=deffile[i].ext) != NULL)
      if (strcasecmp(feptr,eptr)==0)
	break;
      
  return i;
}

static void set_extension(char *buf,int ftp)
{
  int len;
  t_deffile *df;
  
  df=&(deffile[ftp]);
  len=strlen(buf);
  
  /* check if extension is already at end of filename */
  if (len >= 4) 
    if (strcasecmp(&(buf[len-4]),df->ext) == 0) 
      return;
  
  strcat(buf,df->ext);
}

void set_grpfnm(t_filenm *fnm,char *name,int nopts,int ftps[])
{
  char      buf[256];
  int       i,ti;
  bool      bSet;
  
  /* First check whether we have a valid filename already */
  /* ti is the extension type */
  ti=fn2ftp(name); 
  
  for(i=0; (i<nopts); i++) {
    /* check if it is a trajectory */
    if (ti == ftps[i]) 
      break;
  }
  if (i != nopts) {
    /* If it is a trajectory file, return */
    fnm->fn=strdup(name);
    
    return;
  }
  bSet=FALSE;
  if (fnm->flag & ffREAD) { 
    /* for input-files only: search for filenames in the directory */ 
    for(i=0; (i<nopts) && !bSet; i++) {
      ti=ftps[i];
      if (name) {
	strcpy(buf,name);
	set_extension(buf,ti);
      }
      else {
	strcpy(buf,ftp2defnm(ti));
	set_extension(buf,ti);
      }
      if (fexist(buf))
	bSet=TRUE;
    }
  }
  if (!bSet) {
    ti=ftps[0];
    if (name)
      strcpy(buf,name);
    else
      strcpy(buf,ftp2defnm(fnm->ftp));
    set_extension(buf,ti);
  }
  fnm->fn=strdup(buf);
}

static void set_trxnm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NTRXS,trxs);
}

static void set_trnnm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NTRNS,trns);
}

static void set_stonm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NSTOS,stos);
}

static void set_stxnm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NSTXS,stxs);
}

static void set_enxnm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NENXS,enxs);
}

static void set_tpxnm(t_filenm *fnm,char *name)
{
  set_grpfnm(fnm,name,NTPXS,tpxs);
}


static void set_filenm(t_filenm *fnm,char *name)
{
  /* Set the default filename, extension and option for those fields that 
   * are not already set. An extension is added if not present, if fn = NULL
   * or empty, the default filename is given.
   */
  char      buf[256];

  if ((fnm->ftp < 0) || (fnm->ftp >= efNR))
    fatal_error(0,"file type out of range (%d)",fnm->ftp);

  if (fnm->ftp == efTRX) {
    set_trxnm(fnm,name);
  }
  else if (fnm->ftp == efTRN) {
    set_trnnm(fnm,name);
  } 
  else if (fnm->ftp == efSTX) {
    set_stxnm(fnm,name);
  }
  else if (fnm->ftp == efSTO) {
    set_stonm(fnm,name);
  }
  else if (fnm->ftp == efTPX) {
    set_tpxnm(fnm,name);
  }
  else if (fnm->ftp == efENX) {
    set_enxnm(fnm,name);
  }
  else {
    if (name != NULL) {
      strcpy(buf,name);
    }
    else {
      strcpy(buf,deffile[fnm->ftp].defnm);
    }
    set_extension(buf,fnm->ftp);
    
    fnm->fn=strdup(buf);
  }
}

static void set_filenms(int nf,t_filenm fnm[])
{
  int i;

  for(i=0; (i<nf); i++)
    if (!is_set(fnm[i]))
      set_filenm(&(fnm[i]),fnm[i].fn);
}

void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
		     bool bKeep)
{
  int  i,j;
  bool *bRemove;

  check_opts(nf,fnm);
  
  for(i=0; (i<nf); i++)
    un_set(fnm[i]);
    
  if (*argc > 1) {
    snew(bRemove,(*argc)+1);
    i=1;
    do {
      for(j=0; (j<nf); j++) {
	if (strcmp(argv[i],fnm[j].opt) == 0) {
	  do_set(fnm[j]);
	  bRemove[i]=TRUE;
	  i++;
	  if ((i < *argc) && (argv[i][0] != '-')) {
	    set_filenm(&fnm[j],argv[i]);
	    bRemove[i]=TRUE;
	    i++;
	  }
	  else
	    set_filenm(&fnm[j],fnm[j].fn);

	  break;
	}
      }
      /* No file found corresponding to option argv[i] */
      if (j == nf)
	i++;
    } while (i < *argc);
    
    if (!bKeep) {
      /* Remove used entries */
      for(i=j=0; (i<=*argc); i++) {
	if (!bRemove[i])
	  argv[j++]=argv[i];
      }
      (*argc)=j-1;
    }
    sfree(bRemove);
  }
  
  set_filenms(nf,fnm);
}

char *opt2fn(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      return fnm[i].fn;

  fprintf(stderr,"No option %s\n",opt);
  return NULL;
}

char *ftp2fn(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return fnm[i].fn;
  
  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}

bool ftp2bSet(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      return (bool) is_set(fnm[i]);
      
  fprintf(stderr,"ftp2bSet: No filetype %s\n",deffile[ftp].ext);
  
  return FALSE;
}

bool opt2bSet(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      return (bool) is_set(fnm[i]);

  fprintf(stderr,"No option %s\n",opt);
  
  return FALSE;
}

char *opt2fn_null(char *opt,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (strcmp(opt,fnm[i].opt)==0)
      if (is_opt(fnm[i]) && !is_set(fnm[i]))
	return NULL;
      else
	return fnm[i].fn;

  fprintf(stderr,"No option %s\n",opt);
  return NULL;
}

char *ftp2fn_null(int ftp,int nfile,t_filenm fnm[])
{
  int i;
  
  for(i=0; (i<nfile); i++)
    if (ftp == fnm[i].ftp)
      if (is_opt(fnm[i]) && !is_set(fnm[i]))
	return NULL;
      else
	return fnm[i].fn;
  
  fprintf(stderr,"ftp2fn: No filetype %s\n",deffile[ftp].ext);
  return NULL;
}

static void add_filters(char *filter,int *n,int nf,int ftp[])
{
  char buf[8];
  int  i;
  
  for(i=0; (i<nf); i++) {
    sprintf(buf,"*%s",ftp2ext(ftp[i]));
    if (*n > 0)
      strcat(filter,",");
    strcat(filter,buf);
    (*n) ++;
  }
}

char *ftp2filter(int ftp)
{
  int    n;
  static char filter[128];

  filter[0] = '\0';  
  n         = 0;
  switch (ftp) {
  case efENX:
    add_filters(filter,&n,NENXS,enxs);
    break;
  case efTRX:
    add_filters(filter,&n,NTRXS,trxs);
    break;
  case efTRN:
    add_filters(filter,&n,NTRNS,trns);
    break;
  case efSTO:
    add_filters(filter,&n,NSTOS,stos);
    break;
  case efSTX:
    add_filters(filter,&n,NSTXS,stxs);
    break;
  case efTPX:
    add_filters(filter,&n,NTPXS,tpxs);
    break;
  default:
    sprintf(filter,"*%s",ftp2ext(ftp));
    break;
  }
  return filter;
}
