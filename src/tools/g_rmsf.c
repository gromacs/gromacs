/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_rmsf_c = "$Id$";

#include "smalloc.h"
#include "math.h"
#include "macros.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pdbio.h"
#include "tpxio.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "do_fit.h"
#include "princ.h"
#include "rmpbc.h"
#include "confio.h"

static int calc_xav(bool bAverX,
		    char *fn,rvec xref[],t_topology *top,matrix box,
		    real w_rls[],int isize,atom_id index[])
{
  rvec *x,*xav;
  rvec xcm;
  real t,rmsd,xxx,tfac;
  int  i,m,natoms,status,teller;
  
  /* remove pbc */
  rm_pbc(&(top->idef),top->atoms.nr,box,xref,xref);

  /* set center of mass to zero */
  sub_xcm(xref,isize,index,top->atoms.atom,xcm,FALSE);

  /* read first frame  */
  natoms = read_first_x(&status,fn,&t,&x,box);
  if (natoms != top->atoms.nr) 
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top->atoms.nr,natoms);
		
  if (bAverX) {
    snew(xav,natoms);
    teller = 0;
    do {
      /* remove periodic boundary */
      rm_pbc(&(top->idef),top->atoms.nr,box,x,x);
      
      /* set center of mass to zero */
      sub_xcm(x,isize,index,top->atoms.atom,xcm,FALSE);
      
      /* fit to reference structure */
      do_fit(top->atoms.nr,w_rls,xref,x);
      
      for(i=0; (i<natoms); i++)
	rvec_inc(xav[i],x[i]);
      
      teller++;
    } while (read_next_x(status,&t,natoms,x,box));

    tfac = 1.0/teller;    
    rmsd = 0;
    for(i=0; (i<natoms); i++) {
      for(m=0; (m<DIM); m++) {
	xxx        = xav[i][m] * tfac;
	rmsd      += sqr(xref[i][m] - xxx);
	xref[i][m] = xxx;
      }
    }
    rmsd = sqrt(rmsd/natoms);    
    sfree(xav);
    fprintf(stderr,"Computed average structure. RMSD with reference is %g nm\n",
	    rmsd);
  }
  sfree(x);

  rewind_trj(status);
  
  return status;
}

static real find_pdb_bfac(t_atoms *atoms,char *resnm,int resnr,char *atomnm)
{
  char rresnm[8];
  int i;
  
  strcpy(rresnm,resnm);
  rresnm[3]='\0';
  for(i=0; (i<atoms->nr); i++) {
    if ((resnr == atoms->atom[i].resnr) &&
	(strcmp(*atoms->resname[resnr],rresnm) == 0) &&
	(strstr(*atoms->atomname[i],atomnm) != NULL))
      break;
  }
  if (i == atoms->nr) {
    fprintf(stderr,"\rCan not find %s%d-%s in pdbfile\n",
	    rresnm,resnr,atomnm);
    return 0.0;
  }
    
  return atoms->pdbinfo[i].bfac;
}

void correlate_aniso(char *fn,t_atoms *ref,t_atoms *calc)
{
  FILE *fp;
  int  i,j;
  
  fp = xvgropen(fn,"Correlation between X-Ray and Computed Uij","X-Ray","Computed");
  for(i=0; (i<ref->nr); i++) {
    if (ref->pdbinfo[i].bAnisotropic) {
      for(j=U11; (j<=U23); j++)
	fprintf(fp,"%10d  %10d\n",ref->pdbinfo[i].uij[j],calc->pdbinfo[i].uij[j]);
    }
  }
  fclose(fp);
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_rmsf computes the root mean square fluctuation (RMSF, i.e. standard ",
    "deviation) of atomic positions ",
    "after first fitting to a reference frame.[PAR]",
    "When the (optional) pdb file is given, the RMSF values are converted",
    "to B-factor values and plotted with the experimental data.[PAR]",
    "With option -aver the average coordinates will be calculated and used",
    "as reference for fitting (which is useless usually). ",
    "They are also saved to a gro file (which may be usefull).[PAR]",
    "With the option -aniso g_rmsf will compute anisotropic temperature factors",
    "and then it will also output average coordinates and a pdb file with ANISOU",
    "records (corresonding to the -oq option). Please note that the U values",
    "are orientation dependent, so before comparison with experimental data",
    "you should verify that you fit to the experimental coordinates.[PAR]",
    "When a pdb input file is passed to the program and the [TT]-aniso[tt]",
    "flag is set",
    "a correlation plot of the Uij will be created, if any anisotropic",
    "temperature factors are present in the pdb file."
  };
  static bool bAverX=FALSE,bAniso=FALSE;
  t_pargs pargs[] = { 
    { "-aver", FALSE, etBOOL, {&bAverX},
      "Calculate average coordinates first. Requires reading the coordinates twice" },
    { "-aniso",FALSE, etBOOL, {&bAniso},
      "Compute anisotropic termperature factors" }
  };
  int          step,nre,natom,natoms,i,g,m,teller=0;
  real         t,lambda,*w_rls,*w_rms;
  
  t_tpxheader  header;
  t_inputrec   ir;
  t_topology   top;
  t_atoms      *pdbatoms,*refatoms;
  bool         bCont;

  matrix       box,pdbbox;
  rvec         *x,*pdbx,*xref;
  int          status,npdbatoms,res0;
  char         buf[256];
  char         title[STRLEN];
  
  FILE         *fp;               /* the graphics file */
  int          resnr;

  bool         bReadPDB;  
  atom_id      *index;
  int          isize;
  char         *grpnames;

  real         bfac,pdb_bfac;
  matrix       *U=NULL;
  atom_id      aid;
  rvec         *rmsf_xx,*rmsf_x,tmp;
  real         *rmsf;
  int          d;
  real         count=0;
  rvec         xcm;

  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,     ffREAD  },
    { efTRX, "-f",  NULL,     ffREAD  },
    { efPDB, "-q",  NULL,     ffOPTRD },
    { efPDB, "-oq", "anisou", ffOPTWR },
    { efNDX, NULL,  NULL,     ffOPTRD },
    { efXVG, "-o",  "rmsf",   ffWRITE },
    { efXVG, "-oc", "correl", ffOPTWR },
    { efSTO, "-ox", "xaver",  ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pargs),pargs,asize(desc),desc,0,NULL);

  if (bAniso) {
    if (!bAverX) {
      fprintf(stderr,
	      "Anisotropic temperature factors require calculation of average\n"
	      "coordinates before computing RMSF.");
      bAverX = TRUE;
    }
  }
  bReadPDB = ftp2bSet(efPDB,NFILE,fnm);
  
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header);
  snew(x,header.natoms);
  snew(xref,header.natoms);
  snew(w_rls,header.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	   box,&natom,xref,NULL,NULL,&top);

  /* Set box type*/
  init_pbc(box,FALSE);
  
  fprintf(stderr,"Select group(s) for root mean square calculation\n");
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,&grpnames);

  /* Set the weight */
  for(i=0;(i<isize);i++) 
    w_rls[index[i]]=top.atoms.atom[index[i]].m;

  /* Malloc the rmsf arrays */
  snew(rmsf_xx,isize);
  snew(rmsf_x, isize);
  snew(rmsf,isize);
  
  if (bAniso)
    snew(U,isize);
  
  if (bReadPDB) {
    get_stx_coordnum(opt2fn("-q",NFILE,fnm),&npdbatoms);
    snew(pdbatoms,1);
    snew(refatoms,1);
    init_t_atoms(pdbatoms,npdbatoms,TRUE);
    init_t_atoms(refatoms,npdbatoms,TRUE);
    snew(pdbx,npdbatoms);
    /* Read coordinates twice */
    read_stx_conf(opt2fn("-q",NFILE,fnm),title,pdbatoms,pdbx,NULL,pdbbox);    
    read_stx_conf(opt2fn("-q",NFILE,fnm),title,refatoms,pdbx,NULL,pdbbox);    
  }
  else {
    pdbatoms  = &top.atoms;
    refatoms  = &top.atoms;
    pdbx      = xref;
    npdbatoms = pdbatoms->nr;
    copy_mat(box,pdbbox);
  }
  
  /* This routine computes average coordinates. Input in xref is the
   * reference structure, output the average.
   * Also input is the file name, returns the file number.
   */
  status = calc_xav(bAverX,
		    ftp2fn(efTRX,NFILE,fnm),xref,&top,box,w_rls,isize,index);
    
  if (bAverX)
    write_sto_conf_indexed(opt2fn("-ox",NFILE,fnm),
			   "Average coords generated by g_rmsf",
			   &top.atoms,xref,NULL,box,isize,index);
  
  /* Now read the trj again to compute fluctuations */
  teller = 0;
  while (read_next_x(status,&t,natom,x,box)) {
    /* Remove periodic boundary */
    rm_pbc(&(top.idef),top.atoms.nr,box,x,x);
    
    /* Set center of mass to zero */
    sub_xcm(x,isize,index,top.atoms.atom,xcm,FALSE);
    
    /* Fit to reference structure */
    do_fit(top.atoms.nr,w_rls,xref,x);
 
    /* Print time of frame*/
    if ((teller % 10) == 0)
      fprintf(stderr,"\r %5.2f",t);
    
    /* Calculate Anisotropic U Tensor */  
    if (bAniso) {
      for(i=0;(i<isize);i++) {
	aid = index[i];
	rvec_sub(x[aid],xref[aid],tmp);
	for(d=0; (d<DIM); d++)
	  for(m=0; (m<DIM); m++)
	    U[i][d][m] += tmp[d]*tmp[m];
      }
    }
    /* Calculate RMS Fluctuations */
    for(i=0;(i<isize);i++) {
      aid = index[i];
      for(d=0;(d<DIM);d++) {
	rmsf_xx[i][d]+=x[aid][d]*x[aid][d];
	rmsf_x[i][d] +=x[aid][d];
      }
    }
    count += 1.0;
    teller++;
  } 
  close_trj(status);

  if (bAniso) {
    for(i=0;(i<isize);i++) {
      aid = index[i];
      for(d=0; (d<DIM); d++)
	for(m=0; (m<DIM); m++)
	  U[i][d][m] /= count;
      pdbatoms->pdbinfo[aid].bAnisotropic = TRUE;
      pdbatoms->pdbinfo[aid].uij[U11] = 1e6*U[i][XX][XX];
      pdbatoms->pdbinfo[aid].uij[U22] = 1e6*U[i][YY][YY];
      pdbatoms->pdbinfo[aid].uij[U33] = 1e6*U[i][ZZ][ZZ];
      pdbatoms->pdbinfo[aid].uij[U12] = 1e6*U[i][XX][YY];
      pdbatoms->pdbinfo[aid].uij[U13] = 1e6*U[i][XX][ZZ];
      pdbatoms->pdbinfo[aid].uij[U23] = 1e6*U[i][YY][ZZ];
    }
    sfree(U);
  }
  for(i=0;(i<isize);i++) {
    rmsf[i] = (rmsf_xx[i][XX]/count - sqr(rmsf_x[i][XX]/count)+ 
	       rmsf_xx[i][YY]/count - sqr(rmsf_x[i][YY]/count)+ 
	       rmsf_xx[i][ZZ]/count - sqr(rmsf_x[i][ZZ]/count)); 
    if (bAniso)
      pdbatoms->pdbinfo[index[i]].bfac = 800*M_PI/3.0*rmsf[i];
  }
  
  /* Write RMSF output */
  if (bReadPDB) {
    bfac = 8.0*M_PI*M_PI/3.0*100;
    fp   = xvgropen(ftp2fn(efXVG,NFILE,fnm),"B-Factors",
		    "Atom","A\\b\\S\\So\\N\\S 2");
    for(i=0;(i<isize);i++) {
      resnr    = top.atoms.atom[index[i]].resnr;
      pdb_bfac = find_pdb_bfac(pdbatoms,*(top.atoms.resname[resnr]),resnr,
			       *(top.atoms.atomname[index[i]]));
      
      fprintf(fp,"%5d  %10.5f  %10.5f\n",i,rmsf[i]*bfac,pdb_bfac);
    }
  }
  else {
    fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),"RMS fluctuation","Atom","nm");
    for(i=0;(i<isize);i++) 
      fprintf(fp,"%5d %8.4f\n",i,sqrt(rmsf[i]));
  }
  fclose(fp);
  
  /* Write a pdb file with anisou records */
  if (bAniso) {
    write_sto_conf(opt2fn("-oq",NFILE,fnm),title,pdbatoms,pdbx,NULL,pdbbox);
    correlate_aniso(opt2fn("-oc",NFILE,fnm),refatoms,pdbatoms);
    xvgr_file(opt2fn("-oc",NFILE,fnm),"-nxy");
  }
  xvgr_file(opt2fn("-o",NFILE,fnm),"-nxy");
    
  thanx(stdout);
  
  return 0;
}
