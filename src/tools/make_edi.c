/*
 * $Id$
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 *
 * The make_edi program was generously contributed by Oliver Lange, based
 * on the code from g_anaeig. You can reach him as olange@gwdg.de.
 *
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
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */

#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "readinp.h"
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "pdbio.h"
#include "confio.h"
#include "tpxio.h"
#include "matio.h"
#include "mshift.h"
#include "xvgr.h"
#include "do_fit.h"
#include "rmpbc.h"
#include "txtdump.h"
#include "eigio.h"
#include "edsam.h"
#include "index.h" 

void make_t_edx(t_edx *edx, int natoms, rvec *pos, atom_id index[]) {
  edx->nr=natoms;
  edx->anrs=index;
  edx->x=pos;
}

void write_t_edx(FILE *fp, t_edx edx, char *comment) {
 /*here we copy only the pointers into the t_edx struct
  no data is copied and edx.box is ignored  */
 int i;
  fprintf(fp,"#%s \n %d \n",comment,edx.nr);
  for (i=0;i<edx.nr;i++) {
    fprintf(fp,"%d  %f  %f  %f\n",(edx.anrs)[i]+1,(edx.x)[i][XX],(edx.x)[i][YY],(edx.x)[i][ZZ]);
  }
}

int sscan_list(int *list[], char *str, char *listname) {
 /*this routine scans a string of the form 1,3-6,9 and returns the
    selected numbers (in this case 1 3 4 5 6 9) in NULL-terminated array of integers.
    memory for this list will be allocated  in this routine -- sscan_list expects *list to
    be a NULL-Pointer

    listname is a string used in the errormessage*/


   int i,istep;
   char c;
   char *pos,*startpos,*step;
   int n=strlen(str);

   /*enums to define the different lexical stati */
   enum { sBefore, sNumber, sMinus, sRange, sZero, sSmaller, sError, sSteppedRange };
  
   int status=sBefore; /*status of the deterministic automat to scan str   */
   int number=0;
   int end_number=0;
    
   char *start=NULL; /*holds the string of the number behind a ','*/
   char *end=NULL; /*holds the string of the number behind a '-' */
  
   int nvecs=0; /* counts the number of vectors in the list*/

   step=NULL;
   snew(pos,n+4);
   startpos=pos;
   strcpy(pos,str);
   pos[n]=',';  
   pos[n+1]='1';
   pos[n+2]='\0';

   if (*list!=NULL) {
      gmx_fatal(FARGS,"list argument of sscan_list should be NULL"); }
     
   while ((c=*pos)!=0) {
     switch(status) {
        /* expect a number */
        case sBefore: if (isdigit(c)) {
              start=pos;
              status=sNumber;
              break;
           } else status=sError; break;

        /* have read a number, expect ',' or '-' */
        case sNumber: if (c==',') {
             /*store number*/
             srenew(*list,nvecs+1);
             (*list)[nvecs++]=number=atoi(start);
             status=sBefore;
             if (number==0)
                 status=sZero;
              break;
             } else if (c=='-') { status=sMinus; break; }
             else if (isdigit(c)) break;
             else status=sError; break;

        /* have read a '-' -> expect a number */
     case sMinus: 
       if (isdigit(c)) {
	 end=pos;
	 status=sRange; break;
       } else status=sError; break;
       
     case sSteppedRange:
       if (isdigit(c)) {
	 if (step) {
	   status=sError; break; 
	 } else
	   step=pos;
	 status=sRange;
	 break;
       } else status=sError; break;
       
        /* have read the number after a minus, expect ',' or ':' */
        case sRange:
            if (c==',') {
               /*store numbers*/
               end_number=atoi(end);
               number=atoi(start);
               status=sBefore;
               if (number==0) {
                  status=sZero; break;
               }
               if (end_number<=number) {
                  status=sSmaller; break;
               }
               srenew(*list,nvecs+end_number-number+1);
	       if (step) {
		 istep=atoi(step);
		 step=NULL;
	       } else istep=1;
               for (i=number;i<=end_number;i+=istep)
                    (*list)[nvecs++]=i;
               break;
            } else if (c==':') {
	      status = sSteppedRange;
	      break;
	    } else if (isdigit(c)) break; else status=sError;  break;

       /* format error occured */
       case sError:
       gmx_fatal(FARGS,"Error in the list of eigenvectors for %s at pos %d with char %c",listname,pos-startpos,*(pos-1)); 

       /* logical error occured */
       case sZero:
               gmx_fatal(FARGS,"Error in the list of eigenvectors for %s at pos %d: eigenvector 0 is not valid",listname,pos-startpos);
       case sSmaller:
               gmx_fatal(FARGS,"Error in the list of eigenvectors for %s at pos %d: second index %d is not bigger than %d",listname,pos-startpos,end_number,number);

     }
   ++pos; /* read next character */
   } /*scanner has finished */

   /* append zero to list of eigenvectors */
   srenew(*list,nvecs+1);
   (*list)[nvecs]=0;
   sfree(startpos);
   return nvecs;
} /*sscan_list*/
  
void write_eigvec(FILE* fp, int natoms, int eig_list[], rvec** eigvecs,int nvec, char *grouptitle,real steps[]) {
/* eig_list is a zero-terminated list of indices into the eigvecs array.
   eigvecs are coordinates of eigenvectors
   grouptitle to write in the comment line
   steps  -- array with stepsizes for evLINFIX, evLINACC and evRADACC
*/

  int n=0,i; rvec x;
  real sum;
  while (eig_list[n++]);  /*count selected eigenvecs*/
  
  fprintf(fp,"# NUMBER OF EIGENVECTORS + %s\n %d\n",grouptitle,n-1);
  
  /* write list of eigenvector indicess */
  for(n=0;eig_list[n];n++) {
    if (steps)
      fprintf(fp,"%8d   %f\n",eig_list[n],steps[n]); 
    else 
      fprintf(fp,"%8d   %f\n",eig_list[n],1.0);
  }
  n=0;
    
  /* dump coordinates of the selected eigenvectors */
  while (eig_list[n]) {
    sum=0;
    for (i=0; i<natoms; i++) {
      if (eig_list[n]>nvec)
	gmx_fatal(FARGS,"Selected eigenvector %d is higher than maximum number %d of available eigenvectors",eig_list[n],nvec);
      copy_rvec(eigvecs[eig_list[n]-1][i],x);
      sum+=norm2(x);
      fprintf(fp,"%8.5f %8.5f %8.5f\n",x[XX],x[YY],x[ZZ]);      
    };    
    n++;
  };
};


/*enum referring to the different lists of eigenvectors*/
enum { evLINFIX, evLINACC, evRADFIX, evRADACC, evRADCON , evMON, evEND };
#define MAGIC 666


void write_the_whole_thing(FILE* fp, t_edpar *edpars, rvec** eigvecs, int nvec, int *eig_listen[], real* evStepList[]) {
/* write edi-file */

    /*Header*/
    fprintf(fp,"#MAGIC\n %d \n#NINI\n %d\n#NPRO\n %d\n#SELMAS\n %d\n",
        MAGIC,edpars->nini,edpars->npro,edpars->selmas);
    fprintf(fp,"#OUTFRQ\n %d\n#LOGFRQ\n %d\n#MAXLEN\n %d\n#SLOPECRIT\n %f\n",
        edpars->outfrq,edpars->logfrq,edpars->maxedsteps,edpars->slope);

    /* Average and reference positions */
    write_t_edx(fp,edpars->sref,"NREF, XREF");
    write_t_edx(fp,edpars->sav,"NAV, XAV");

    fprintf(fp,"#NED\n %d\n",edpars->ned);
    /*Eigenvectors */
 
    write_eigvec(fp, edpars->ned, eig_listen[evMON],eigvecs,nvec,"COMPONENTS GROUP 1",NULL);
    write_eigvec(fp, edpars->ned, eig_listen[evLINFIX],eigvecs,nvec,"COMPONENTS GROUP 2",evStepList[evLINFIX]);
    write_eigvec(fp, edpars->ned, eig_listen[evLINACC],eigvecs,nvec,"COMPONENTS GROUP 3",evStepList[evLINACC]);
    write_eigvec(fp, edpars->ned, eig_listen[evRADFIX],eigvecs,nvec,"COMPONENTS GROUP 4",evStepList[evRADFIX]);
    write_eigvec(fp, edpars->ned, eig_listen[evRADACC],eigvecs,nvec,"COMPONENTS GROUP 5",NULL);
    write_eigvec(fp, edpars->ned, eig_listen[evRADCON],eigvecs,nvec,"COMPONENTS GROUP 6",NULL);

    /*Target and Origin positions */
    write_t_edx(fp,edpars->star,"NTARGET, XTARGET");
    write_t_edx(fp,edpars->sori,"NORIGIN, XORIGIN");
}; 

int read_conffile(char *confin,char *title,rvec *x[]) {
/* read coordinates out of STX file  */
  int natoms;
  t_atoms  confat;
  matrix box;
  printf("read coordnumber from file %s\n",confin);
  get_stx_coordnum(confin,&natoms);
  printf("number of coordinates in file %d\n",natoms);
/*  if (natoms != ncoords)
     gmx_fatal(FARGS,"number of coordinates in coordinate file (%s, %d)\n"
           "             does not match topology (= %d)",
           confin,natoms,ncoords);
  else {*/
    /* make space for coordinates and velocities */
    init_t_atoms(&confat,natoms,FALSE);
    printf("init_t\n");
    snew(*x,natoms);
    read_stx_conf(confin,title,&confat,*x,NULL,box);
    return natoms;
};   


static real *scan_vecparams(char *str,char * par, int nvecs)
{
  char   f0[256],f1[256];             /*format strings adapted every pass of the loop*/
  double d,tcap=0;
  int    i;
  real   *vec_params;

  snew(vec_params,nvecs);
  if (str) {
    f0[0] = '\0';
    for(i=0; (i<nvecs); i++) {
      strcpy(f1,f0);  /*f0 is the format string for the "to-be-ignored" numbers*/
      strcat(f1,"%lf"); /*and f1 to read the actual number in this pass of the loop*/
      if (sscanf(str,f1,&d) != 1)
	gmx_fatal(FARGS,"Not enough elements for %s parameter (I need %d)",par,nvecs);
      vec_params[i] = d;
      tcap += d;
      strcat(f0,"%*s");
    }
  }  
  return vec_params;
}    


void init_edx(t_edx *edx) {
  edx->nr=0;
  snew(edx->x,1);
  snew(edx->anrs,1);
};

void filter2edx(t_edx *edx,int nindex, atom_id index[],int ngro, atom_id igro[],rvec *x,char* structure) {
/* filter2edx copies coordinates from x to edx which are given in index
*/
  
   int pos,i;
   int ix=edx->nr;
   edx->nr+=nindex;
   srenew(edx->x,edx->nr);
   srenew(edx->anrs,edx->nr);
   for (i=0;i<nindex;i++,ix++) {
         for (pos=0; pos<ngro-1 && igro[pos]!=index[i] ; ++pos) {};  /*search element in igro*/
         if (igro[pos]!=index[i])
              gmx_fatal(FARGS,"Couldn't find atom with index %d in structure %s",index[i],structure);
         edx->anrs[ix]=index[i];
         copy_rvec(x[pos],edx->x[ix]);
   };
};

void get_structure(t_atoms *atoms,char *IndexFile,char *StructureFile,t_edx *edx,int nfit,
                    atom_id ifit[],int natoms, atom_id index[]) {


  atom_id *igro;  /*index corresponding to target or origin structure*/
  int ngro;
  int ntar;
  rvec *xtar;
  char  title[STRLEN];
  char* grpname;
  

  ntar=read_conffile(StructureFile,title,&xtar);
  printf("Select an index group of %d elements that corresponds to the atoms in the structure file %s\n",
            ntar,StructureFile);
  get_index(atoms,IndexFile,1,&ngro,&igro,&grpname);
  if (ngro!=ntar)
     gmx_fatal(FARGS,"You selected an index group with %d elements instead of %d",ngro,ntar);
  init_edx(edx);
  filter2edx(edx,nfit,ifit,ngro,igro,xtar,StructureFile);
  if (ifit!=index) /*if fit structure is different append these coordinates, too -- don't mind duplicates*/
     filter2edx(edx,natoms,index,ngro,igro,xtar,StructureFile);
};

int main(int argc,char *argv[])
{

  static char *desc[] = {
    "[TT]make_edi[tt] generates an ED-sampling input file to be used with mdrun",
    "based on eigenvectors of a covariance matrix ([TT]g_covar[tt]) or from a", 
    "Normal Modes anaysis ([TT]g_nmeig[tt]).",
    "ED-sampling can be used to manipulate the position along collective coordinates",
    "(eigenvectors) of (biological) macromolecules during a simulation. Particularly,",
    "it may be used to enhance the sampling efficiency of MD simulations by stimulating",
    "the system to explore new regions along these collective coordinates. A number",
    "of different algorithms are implemented to drive the system along the eigenvectors",
    "([TT]-linfix[tt], [TT]-linacc[tt], [TT]-radfix[tt], [TT]-radacc[tt], [TT]-radcon[tt]),",
    "to keep the position along a certain (set of) coordinate(s) fixed ([TT]-linfix[tt]),",
    "or to only monitor the projections of the positions, velocities and forces onto",
    "these coordinates([TT]-mon[tt]).[PAR]"
    "References:[BR]",
    "A. Amadei, A.B.M. Linssen, B.L. de Groot, D.M.F. van Aalten and ",
    "H.J.C. Berendsen; An efficient method for sampling the essential subspace ",
    "of proteins., J. Biomol. Struct. Dyn. 13:615-626 (1996)[BR]",
    "B.L. de Groot, A. Amadei, D.M.F. van Aalten and H.J.C. Berendsen; ",
    "Towards an exhaustive sampling of the configurational spaces of the ",
    "two forms of the peptide hormone guanylin,"
     "J. Biomol. Struct. Dyn. 13 : 741-751 (1996)[BR]",
    "B.L. de Groot, A.Amadei, R.M. Scheek, N.A.J. van Nuland and H.J.C. Berendsen; ",
    "An extended sampling of the configurational space of HPr from E. coli",
    "PROTEINS: Struct. Funct. Gen. 26: 314-322 (1996)",
    "[PAR]You will be prompted for one or more index groups that correspond to the eigenvectors,",
    "reference structure, target positions, etc.[PAR]",

    "[TT]-mon[tt]: monitor projections of x, v and f onto selected eigenvectors.[PAR]",
    "[TT]-linfix[tt]: perform fixed-step linear expansion along selected eigenvectors.[PAR]",
    "[TT]-linacc[tt]: perform acceptance linear expansion along selected eigenvectors.",
    "(steps in the desired directions will be accepted, others will be rejected).[PAR]",
    "[TT]-radfix[tt]: perform fixed-step radius expansion along selected eigenvectors.[PAR]",
    "[TT]-radacc[tt]: perform acceptance radius expansion along selected eigenvectors.",
    "(steps in the desired direction will be accepted, others will be rejected).",
    "Note: by default the starting MD structure will be taken as origin of the first",
    "expansion cycle for radius expansion. If [TT]-ori[tt] is specified, you will be able",
    "to read in a structure file that defines an external origin.[PAR]"
    "[TT]-radcon[tt]: perform acceptance radius contraction along selected eigenvectors",
    "towards a target structure specified with [TT]-tar[tt]."
    "NOTE: each eigenvector can be selected only once. [PAR]"
    "[TT]-outfrq[tt]: frequency (in steps) of writing out projections etc.[PAR]",
    "[TT]-logfrq[tt]: frequency (in steps) of writing out statistics to log file.[PAR]",
    "[TT]-slope[tt]: minimal slope in acceptance radius expansion. A new expansion",
    "cycle will be started if the spontaneous increase of the radius (in nm/step)",
    "is less than the value specified.[PAR]" 
    "[TT]-maxedsteps[tt]: maximum number of steps per cycle in radius expansion",
    "before a new cycle is started.[PAR]"
    "Note on the parallel implementation: since ED sampling is a 'global' thing",
    "(collective coordinates etc), at least on the 'protein' side, ED sampling",
    "is not very parallel-friendly from an implentation point of view (it would",
    "require much extra communication to fully parallelize the algorithms).",
    "Fortunately, however, a typical parallel protein simulation in gromacs has",
    "most or all protein coordinates on one processor (the master) and has only",
    "other atoms (solvent, lipid, ions etc) on the other processors. With such a",
    "setup, ED sampling will still work. If the atoms over which ED sampling should ",
    "be performed are spread over multiple processors, a fatal error will result.[PAR]"
    "All output of mdrun (specify with -eo) is written to a .edo file (some extra",
    "information is written to the log file of mdrun too, actually). The .edo format",
    "is a simple ASCII file that should be easy to parse with standard unix tools",
    "like awk. A script (parse_edo) can be downloaded from contribution section at",
    " www.gromacs.org to extract information from the",
    ".edo files for your convinience. In short, the header defines which information",
    "can be expected in the rest of the .edo file. After the header, per step the",
    "following information is present: [PAR]",
    "* the step number[BR]",
    "* RMSD (for atoms in fitting prior to calculating ED constr.)[BR]",
    "* projections of the positions onto selected eigenvectors[BR]",
    "* projections of the velocities onto selected eigenvectors[BR]",
    "* projections of the forces onto selected eigenvectors",
    "[PAR]",
    "All projections are in the same order as in the header, so if you have e.g.",
    "2 groups (say one group over which acceptance radius expansion is performed,",
    "and another for which the projections are merely monitored) then you first",
    "get the position projections for each of the 2 groups, then the velocities",
    "and then the forces. Radii are not explicitly written to the .edo file, as",
    "they can be readily projected back from the positions. Alternatively, radii",
    "may be 'grepped from the log file."
  };
  
  t_edpar edi_params;  /*save all the params in this struct and then save it in an edi file.
                         ignoring fields nmass,massnrs,mass,tmass,nfit,fitnrs,edo*/


  static int  first=1,last=8,skip=1,nextr=2;
  static real max=0.0;
  static bool bSplit=FALSE;
   enum  { evSTEPEND = evRADFIX + 1}; 
  static char* evSelections[evEND]= {NULL,NULL,NULL,NULL,NULL};
  static char* evOptions[evEND] = {"-linfix","-linacc","-radfix","-radacc","-radcon","-mon"};
  static char* evParams[evSTEPEND] ={NULL,NULL};
  static char* evStepOptions[evSTEPEND] = {"-linstep","-accdir","-radstep"};
  static real* evStepList[evSTEPEND];
  static real  radfix=0.0;
  static int* listen[evEND];
  t_pargs pa[] = {
    { evOptions[evMON], FALSE, etSTR, {&evSelections[evMON]},
      "Indices of eigenvectors  for projections of x, v and f (e.g. 1,2-5,9) or 1-100:10 means 1 11 21 31 ... 91 " },
    { evOptions[evLINFIX], FALSE, etSTR, {&evSelections[evLINFIX]},
      "Indices of eigenvectors for fixed increment linear sampling" },
    { evOptions[evLINACC], FALSE, etSTR, {&evSelections[evLINACC]},
      "Indices of eigenvectors for acceptance linear sampling" },
    { evOptions[evRADFIX], FALSE, etSTR, {&evSelections[evRADFIX]},
      "Indices of eigenvectors for fixed increment radius expansion" },
    { evOptions[evRADACC], FALSE, etSTR, {&evSelections[evRADACC]},
      "Indices of eigenvectors for acceptance radius expansion" },
    { evOptions[evRADCON], FALSE, etSTR, {&evSelections[evRADCON]},
      "Indices of eigenvectors for acceptance radius contraction" },
    { "-outfrq", FALSE, etINT, {&edi_params.outfrq},
      "freqency (in steps) of writing output in .edo file" },
    { "-logfrq", FALSE, etINT, {&edi_params.logfrq},
      "frequency (in steps) of writing to log" },
    { "-slope", FALSE, etREAL, { &edi_params.slope},
      "minimal slope in acceptance radius expamsion"},
    { "-maxedsteps", FALSE, etINT, {&edi_params.maxedsteps},
      "max nr of steps per cycle" },
    { evStepOptions[evLINFIX], FALSE, etSTR, {&evParams[evLINFIX]},
      "Stepsizes (nm/step) for fixed increment linear sampling (put in quotes! \"1.0 2.3 5.1 -3.1\")"},
    { evStepOptions[evLINACC], FALSE, etSTR, {&evParams[evLINACC]},
      "Directions for acceptance linear sampling - only sign counts! (put in quotes! \"-1 +1 -1.1\")"},
    { evStepOptions[evRADFIX], FALSE, etREAL, {&radfix},
      "Stepsize (nm/step) for fixed increment radius expansion"}
  };
#define NPA asize(pa)

  rvec       *xref1;
  bool       bDMA1;
  int        nvec1,*eignr1=NULL;
  rvec       *xav1,**eigvec1=NULL;
  rvec       *xtar,*xori;
  t_atoms    *atoms=NULL;
  int natoms;
  char       *grpname,*indexfile;
  int        i;
  atom_id    *index,*ifit;
  int        nfit;
  int ev_class; /* parameter _class i.e. evMON, evRADFIX etc. */
  int nvecs;


  char       *EdiFile;
  char       *TargetFile;
  char       *OriginFile;


/*to read topology file*/
   t_topology top;
   t_topology tartop,oritop;
   matrix     tarbox,oribox;
   char       title[STRLEN];
   matrix     topbox;
   rvec       *xtop;
   bool bTop, bM, bFit1;

  t_filenm fnm[] = {
    { efTRN, "-v",    "eigenvec",    ffREAD  },
    { efTPS, NULL,    NULL,          ffREAD },
    { efNDX, NULL,    NULL,  ffOPTRD },
    { efSTX, "-tar", "target", ffOPTRD},
    { efSTX, "-ori", "origin", ffOPTRD},
    { efEDI, "-o", "sam", ffWRITE }
  };
#define NFILE asize(fnm)
 edi_params.outfrq=100; edi_params.logfrq=100; edi_params.slope=0.0; edi_params.maxedsteps=0;
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv, 0 ,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  indexfile=ftp2fn_null(efNDX,NFILE,fnm);
  EdiFile=ftp2fn(efEDI,NFILE,fnm);
  TargetFile      = opt2fn_null("-tar",NFILE,fnm);
  OriginFile      = opt2fn_null("-ori",NFILE,fnm);
 

  for (ev_class=0; ev_class<evEND; ++ev_class)
     if (opt2parg_bSet(evOptions[ev_class],NPA,pa)) {
        /*get list of eigenvectors*/
        nvecs=sscan_list(&listen[ev_class],opt2parg_str(evOptions[ev_class],NPA,pa),evOptions[ev_class]);
        if (ev_class<evSTEPEND-1)
          /*if apropriate get list of stepsizes for these eigenvectors*/
          if (opt2parg_bSet(evStepOptions[ev_class],NPA,pa)) 
            evStepList[ev_class]=
              scan_vecparams(opt2parg_str(evStepOptions[ev_class],NPA,pa),evStepOptions[ev_class],nvecs);
          else { /*if list is not given fill with zeros */
            snew(evStepList[ev_class],nvecs);
            for (i=0; i<nvecs; i++) 
               evStepList[ev_class][i]=0.0;
          }
        else if (ev_class == evRADFIX && opt2parg_bSet(evStepOptions[ev_class],NPA,pa)) {
              snew(evStepList[ev_class],nvecs);
              for (i=0; i<nvecs; i++)
                 evStepList[ev_class][i]=radfix;

             } else {}; /*to avoid ambiguity   */
      } else { /* if there are no eigenvectors for this option set list to zero */
        snew(listen[ev_class],1);
        listen[ev_class][0]=0;
  };

  /* print the interpreted list of eigenvectors - to give some feedback*/
  for (ev_class=0; ev_class<evEND; ++ev_class) {
     printf("list %s consist of the indices:",evOptions[ev_class]);
     i=0;
     while(listen[ev_class][i])
        printf("%d ",listen[ev_class][i++]);
     printf("\n");
  }

 
  /*read eigenvectors from eigvec.trr*/
  read_eigenvectors(opt2fn("-v",NFILE,fnm),&natoms,&bFit1,
          &xref1,&edi_params.selmas,&xav1,&bDMA1,&nvec1,&eignr1,&eigvec1);

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),
               title,&top,&xtop,NULL,topbox,0);
  atoms=&top.atoms;
  edi_params.npro=0;


  printf("\nSelect an index group of %d elements that corresponds to the eigenvectors\n",natoms);
  get_index(atoms,indexfile,1,&i,&index,&grpname); /*if indexfile != NULL parameter 'atoms' is ignored */
  if (i!=natoms) {
      gmx_fatal(FARGS,"you selected a group with %d elements instead of %d",
       		  i,natoms);
  }
  printf("\n");


  if (xref1==NULL && bFit1) {   /* if g_covar used different coordinate groups to fit and to do the PCA */
	  printf("\nNote: the structure in %s should be the same\n"
		 "      as the one used for the fit in g_covar\n",ftp2fn(efTPS,NFILE,fnm));
	  printf("\nSelect the index group that was used for the least squares fit in g_covar\n");
	  get_index(atoms,indexfile,1,&nfit,&ifit,&grpname);
    snew(xref1,nfit);
    for (i=0;i<nfit;i++)
      copy_rvec(xtop[ifit[i]],xref1[i]);
  } else {
     nfit=natoms;
     ifit=index;
  };

 
  /*number of eigenvectors*/
  edi_params.ned=natoms;

  

  /*number of system atoms  */
  edi_params.nini=atoms->nr;

 
  /*store reference and average structure in edi_params*/
  make_t_edx(&edi_params.sref,nfit,xref1,ifit);
  make_t_edx(&edi_params.sav,natoms,xav1,index);

                                                         
  /*store target positions in edi_params*/
  if (opt2bSet("-tar",NFILE,fnm)) {
    get_structure(atoms,indexfile,TargetFile,&edi_params.star,nfit,
                   ifit,natoms,index);
 } else
      make_t_edx(&edi_params.star,0,NULL,index);
     /*store origin positions*/
 if (opt2bSet("-ori",NFILE,fnm)) {
     get_structure(atoms,indexfile,OriginFile,&edi_params.sori,nfit,
                   ifit,natoms,index);
 } else
      make_t_edx(&edi_params.sori,0,NULL,index);
  
  /*write edi-file*/
  write_the_whole_thing(ffopen(EdiFile,"w"), &edi_params, eigvec1,nvec1, listen, evStepList);
  thanx(stderr);
  return 0;
}

