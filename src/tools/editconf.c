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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_editconf_c = "$Id$";

#include <math.h>
#include <string.h>
#include <ctype.h>
#include "pdbio.h"
#include "confio.h"
#include "symtab.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "strdb.h"
#include "rdgroup.h"
#include "vec.h"
#include "typedefs.h"
#include "gstat.h"
#include "strdb.h"
#include "rdgroup.h"
#include "physics.h"
#include "mass.h"
#include "3dview.h"

typedef struct {
  char   sanm[12];
  int    natm;
  int    nw;
  char   anm[6][12];
}  t_simat;

typedef struct {
  char     reso[12];
  char     resn[12];
  int      nsatm;
  t_simat sat[3];
} t_simlist;
static char *pdbtp[epdbNR]={"ATOM  ","HETATM"};
static char *pdbformat=
"%6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";

real calc_mass(t_atoms *atoms)
{
  real tmass;
  int i,j,k;

  tmass = 0;
  for(i=0; (i<atoms->nr); i++) {
    if ((atoms->atom[i].m = get_mass(*atoms->resname[atoms->atom[i].resnr], 
				     *atoms->atomname[i])) == 0.0) {
      if ( ((*atoms->atomname[i])[0]=='H') ||
	   (isdigit((*atoms->atomname[i])[0]) && 
	    ((*atoms->atomname[i])[1]=='H')) ) {
	atoms->atom[i].m=1.008; /* proton mass */
      } else {
	atoms->atom[i].m=12.0110; /* carbon mass */
      }
    }
    tmass += atoms->atom[i].m;
  }

  return tmass;
}

void calc_geom(int natom, rvec *x, rvec geom_center, rvec min, rvec max)
{
  int i,j;
  
  clear_rvec(geom_center);
  for (j=0; (j<DIM); j++)
    min[j]=max[j]=x[0][j];
  for (i=0; (i<natom); i++) {
    rvec_inc(geom_center,x[i]);
    for (j=0; (j<DIM); j++) {
      if (x[i][j] < min[j]) min[j]=x[i][j];
      if (x[i][j] > max[j]) max[j]=x[i][j];
    }
  }
  svmul(1./natom,geom_center,geom_center);
}

void center_conf(int natom, rvec *x, rvec center, rvec geom_cent)
{
  int       i;
  rvec shift;
  
  rvec_sub(center,geom_cent,shift);

  printf("shift     : %6.3f %6.3f %6.3f\n",
	 shift[XX],shift[YY],shift[ZZ]);

  for (i=0; (i<natom); i++) 
    rvec_inc(x[i], shift);
}

void scale_conf(int natom,rvec x[],matrix box,rvec scale)
{
  int i,j;
  
  for(i=0; (i<natom); i++) {
    for (j=0; (j<DIM); j++)
      x[i][j] *= scale[j];
  }
  for (j=0; (j<DIM); j++)
    box[j][j] *= scale[j];
}

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,d;
  
  /* check periodic boundary */
  for(d=0;(d<DIM);d++) {
    for(n=1;(n<atoms->nr);n++) {
      dist = x[n][d]-x[n-1][d];
      if ( fabs(dist) > 0.9 * box[d][d]  ) {
	if ( dist >  0 )
	  x[n][d]-=box[d][d];
	else
	  x[n][d]+=box[d][d];
      } 	
    }
  }
}

void read_bfac(char *fn, int *n_bfac, double **bfac_val, int **bfac_nr)
{
  int  i;
  char **bfac_lines;

  *n_bfac = get_lines(fn, &bfac_lines);
  snew(*bfac_val, *n_bfac);
  snew(*bfac_nr, *n_bfac);
  fprintf(stderr, "Reading %d B-factors from %s\n",*n_bfac,fn);
  for(i=0; (i<*n_bfac); i++) {
    /*fprintf(stderr, "Line %d: %s",i,bfac_lines[i]);*/
    sscanf(bfac_lines[i],"%d %lf",&(*bfac_nr)[i],&(*bfac_val)[i]);
    /*fprintf(stderr," nr %d val %g\n",(*bfac_nr)[i],(*bfac_val)[i]);*/
  }
  
}

void set_pdb_conf_bfac(int natoms,int nres,t_atoms *atoms,rvec x[],
		       matrix box,int n_bfac,double *bfac,int *bfac_nr,
		       bool peratom, bool bLegend)
{
  FILE *out;
  real bfac_min,bfac_max;
  int  i,n;
  bool found;
  char buf[120];

  bfac_max=-1e10;
  bfac_min=1e10;
  for(i=0; (i<n_bfac); i++) {
    if (bfac_nr[i]-1>=atoms->nres)
      peratom=TRUE;
    if ((bfac_nr[i]-1<0) || (bfac_nr[i]-1>=atoms->nr))
      fatal_error(0,"Index of B-Factor %d is out of range: %d (%g)",
		  i+1,bfac_nr[i],bfac[i]);
    if (bfac[i] > bfac_max) 
      bfac_max = bfac[i];
    if (bfac[i] < bfac_min) 
      bfac_min = bfac[i];
  }
  while ( (bfac_max > 99.99) || (bfac_min < -99.99) ) {
    fprintf(stderr,"Range of values for B-factors too large (min %g, max %g) "
	    "will scale down a factor 10\n",bfac_min,bfac_max);
    for(i=0; (i<n_bfac); i++)
      bfac[i] /= 10;
    bfac_max /= 10;
    bfac_min /= 10;
  }
  while ( (abs(bfac_max) < 0.5) && (abs(bfac_min) < 0.5) ) {
    fprintf(stderr,"Range of values for B-factors too small (min %g, max %g) "
	    "will scale up a factor 10\n",bfac_min,bfac_max);
    for(i=0; (i<n_bfac); i++)
      bfac[i] *= 10;
    bfac_max *= 10;
    bfac_min *= 10;
  }
  
  for(i=0; (i<natoms); i++)
    atoms->pdbinfo[i].bfac=0;
  
  if (!peratom) {
    fprintf(stderr,"Will attach %d B-factors to %d residues\n",
	    n_bfac,nres);
    for(i=0; (i<n_bfac); i++) {
      found=FALSE;
      for(n=0; (n<natoms); n++)
	if ( bfac_nr[i] == (atoms->atom[n].resnr+1) ) {
	  atoms->pdbinfo[n].bfac=bfac[i];
	  found=TRUE;
	}
      if (!found) {
	sprintf(buf,"Residue nr %d not found\n",bfac_nr[i]);
	warning(buf);
      }
    }
  } else {
    fprintf(stderr,"Will attach %d B-factors to %d atoms\n",n_bfac,natoms);
    for(i=0; (i<n_bfac); i++) {
      atoms->pdbinfo[bfac_nr[i]-1].bfac=bfac[i];
    }
  }
}

void pdb_legend(FILE *out,int natoms,int nres,t_atoms *atoms,rvec x[])
{
  real bfac_min,bfac_max,xmin,ymin,zmin;
  int  i;
  char buf[256];
  
  bfac_max=-1e10;
  bfac_min=1e10;
  xmin = 1e10;
  ymin = 1e10;
  zmin = 1e10;
  for (i=0; (i<natoms); i++) {
    xmin     = min(xmin,x[i][XX]);
    ymin     = min(ymin,x[i][YY]);
    zmin     = min(zmin,x[i][ZZ]);
    bfac_min = min(bfac_min,atoms->pdbinfo[i].bfac);
    bfac_max = max(bfac_max,atoms->pdbinfo[i].bfac);
  }
  fprintf(stderr,"B-factors range from %g to %g\n",bfac_min,bfac_max);
  sprintf(buf,"%s","LEG");
  buf[3]='\0';
  for (i=1; (i<12); i++) {
    fprintf(out,pdbformat,
	    "ATOM  ",natoms+1+i,"CA",buf,' ',nres+1,
	    (xmin+(i*0.12))*10,ymin*10,zmin*10,1.0,
	    bfac_min+ ((i-1.0)*(bfac_max-bfac_min)/10) );
  }
}

void simple_rotate_conf(int natom,rvec x[],rvec v[],matrix box,char *dir)
{
  real dum;
  
  switch (dir[0]) {
  case 'X':
    rotate_conf(natom,x,v,90,0,0);
    dum = box[ZZ][ZZ];
    box[ZZ][ZZ] = -box[YY][YY];
    box[YY][YY] = dum;
    break;
  case 'Y':
    rotate_conf(natom,x,v,0,90,0);
    dum = box[XX][XX];
    box[XX][XX] = -box[ZZ][ZZ];
    box[ZZ][ZZ] = dum;
    break;
  case 'Z':
    rotate_conf(natom,x,v,0,0,90);
    dum = box[YY][YY];
    box[YY][YY] = -box[XX][XX];
    box[XX][XX] = dum;
    break;
  }
  fprintf(stderr,"WARNING Rotating not debugged CHECK OUTPUT...\n");
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "editconf converts generic structure format to [TT].gro[tt] or",
    "[TT].pdb[tt].[PAR]",
    "A number of options is present to modify the coordinates",
    "and box. [TT]-d[tt], [TT]-dc[tt] and [TT]-box[tt] modify the box and",
    "center the coordinates relative to the new box.",
    "[TT]-dc[tt] takes precedent over [TT]-d[tt]. [TT]-box[tt]",
    "takes precedent over [TT]-dc[tt] and [TT]-d[tt].[PAR]",
    "Scaling is applied before any of the other operations are",
    "performed. Boxes can be scaled to give a certain density (option",
    "[TT]-density[tt]).[PAR]",
    "Groups are selected after all operations have been applied.[PAR]",
    "Periodicity can be removed in a crude manner.",
    "It is important that the box sizes at the bottom of your input file",
    "are correct when the periodicity is to be removed.[PAR]",
    "When writing [TT].pdb[tt] files, B-factors can be",
    "added with the [TT]-bf[tt] option. B-factors are read",
    "from a file with with following format: first line states number of",
    "entries in the file, next lines state an index",
    "followed by a B-factor. The B-factors will be attached per residue",
    "unless an index is larger than the number of residues or unless the",
    "[TT]-atom[tt] option is set. Obviously, any type of numeric data can",
    "be added instead of B-factors. [TT]-legend[tt] will produce",
    "a row of CA atoms with B-factors ranging from the minimum to the",
    "maximum value found, effectively making a legend for viewing.[PAR]",
    "Finally with option [TT]-label[tt] editconf can add a chain identifier",
    "to a pdb file, which can be useful for analysis with e.g. rasmol."
  };
  static char *bugs[] = {
    "For complex molecules, the periodicity removal routine may break down, "
    "in that case you can use trjconv"
  };
  static real dist   = 0.0,rbox=0.0;
  static bool bNDEF=FALSE,bRMPBC=FALSE,bCenter=FALSE;
  static bool peratom=FALSE,bLegend=FALSE;
  static rvec scale={1.0,1.0,1.0},newbox={0.0,0.0,0.0};
  static real rho=1000.0;
  static rvec center={0.0,0.0,0.0};
  static char *label="A";
  static char  *cRotate[] = { NULL, "X", "Y", "Z", NULL };  
  t_pargs pa[] = {
    { "-ndef", FALSE, etBOOL, &bNDEF, 
      "Choose output from default index groups" },    
    { "-d", FALSE, etREAL, &dist, 
      "Distance between the solute and the rectangular box" }, 
    { "-dc", FALSE, etREAL, &dist,
      "Distance between the solute and the cubic box" },
    { "-box", FALSE, etRVEC, &newbox, "Size of box" },
    { "-c", FALSE, etBOOL, &bCenter,
      "Center molecule in box (implied by -d -dc -box)" },
    { "-center", FALSE, etRVEC, &center, "Coordinates of geometrical center"},
    { "-scale", FALSE, etRVEC, &scale, "Scaling factor" },
    { "-density",FALSE, etREAL, &rho, 
      "Density (g/l) of the output box achieved by scaling" },
    { "-pbc",  FALSE, etBOOL, &bRMPBC, 
      "Remove the periodicity (make molecule whole again)" },
    { "-atom", FALSE, etBOOL, &peratom, "Force B-factor attachment per atom" },
    { "-legend",FALSE,etBOOL, &bLegend, "Make B-factor legend" },
    { "-label", FALSE, etSTR, &label,   "Add chain label for all residues" },
    { "-rotate", FALSE, etENUM, cRotate,
      "Perform a rotation of the box around the given coordinate" }
  };
#define NPA asize(pa)

  FILE      *out;
  char      *infile,*outfile,title[STRLEN];
  int       outftp,natom,i,j,n_bfac;
  double    *bfac=NULL;
  int       *bfac_nr=NULL;
  t_atoms   atoms;
  char      *groupnames;
  int       isize;
  atom_id   *index;
  rvec      *x,*v,gc,min,max,size;
  matrix    box;
  bool      bSetSize,bCubic,bDist,bSetCenter;
  bool      bHaveV,bScale,bRho,bRotate;
  real      xs,ys,zs,xcent,ycent,zcent,d;
  t_filenm fnm[] = {
    { efSTX, "-f", NULL, ffREAD },
    { efNDX, "-n", NULL, ffOPTRD },
    { efSTO, NULL, NULL, ffWRITE },
    { efDAT, "-bf", "bfact", ffOPTRD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,NPA,pa,
		    asize(desc),desc,asize(bugs),bugs);

  bSetSize  = opt2parg_bSet("-box" ,NPA,pa);
  bSetCenter= opt2parg_bSet("-center" ,NPA,pa);
  bCubic    = opt2parg_bSet("-dc",NPA,pa);
  bDist     = opt2parg_bSet("-d" ,NPA,pa) || bCubic;
  bCenter   = bCenter || bDist || bSetCenter || bSetSize;
  bScale    = opt2parg_bSet("-scale" ,NPA,pa);
  bRho      = opt2parg_bSet("-density",NPA,pa);
  bRotate   = opt2parg_bSet("-rotate",NPA,pa);
  if (bScale && bRho)
    fprintf(stderr,"WARNING: setting -density overrides -scale");
  bScale  = bScale || bRho;
  
  infile=ftp2fn(efSTX,NFILE,fnm);
  outfile=ftp2fn(efSTO,NFILE,fnm);
  outftp=fn2ftp(outfile);
  
  get_stx_coordnum(infile,&natom);
  init_t_atoms(&atoms,natom,TRUE);
  snew(x,natom);
  snew(v,natom);
  read_stx_conf(infile,title,&atoms,x,v,box);
  printf("Read %d atoms\n",atoms.nr); 

  bHaveV=FALSE;
  for (i=0; (i<natom) && !bHaveV; i++)
    for (j=0; (j<DIM) && !bHaveV; j++)
      bHaveV=bHaveV || (v[i][j]!=0);
  printf("%selocities found\n",bHaveV?"V":"No v");
  
  /* remove pbc */
  if (bRMPBC) 
    rm_gropbc(&atoms,x,box);
    
  /* calc geometrical center and max and min coordinates and size */
  calc_geom(natom, x, gc, min, max);
  rvec_sub(max, min, size);
  printf("size      : %6.3f %6.3f %6.3f\n", size[XX], size[YY], size[ZZ]);
  printf("center    : %6.3f %6.3f %6.3f\n", gc[XX], gc[YY], gc[ZZ]);
  printf("box       : %6.3f %6.3f %6.3f\n", 
	 box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
    
  if ( bScale ) {
    /* scale coordinates and box */
    if (bRho) {
      /* Compute scaling constant */
      real vol,mass,dens;
      
      vol = det(box);
      mass = calc_mass(&atoms);
      dens = (mass*AMU)/(vol*NANO*NANO*NANO);
      fprintf(stderr,"Volume  of input %g (nm3)\n",vol);
      fprintf(stderr,"Mass    of input %g (a.m.u.)\n",mass);
      fprintf(stderr,"Density of input %g (g/l)\n",dens);
      if (vol==0.0)
	fatal_error(0,"Cannot scale density with zero box\n");
      if (mass==0.0)
	fatal_error(0,"Cannot scale density with zero mass\n");

      scale[XX] = scale[YY] = scale[ZZ] = pow(dens/rho,1.0/3.0);
      fprintf(stderr,"Scaling all box edges by %g\n",scale[XX]);
    }
    scale_conf(atoms.nr,x,box,scale);
  
    /* recalc geometrical center and max and min coordinates and size */
    calc_geom(natom, x, gc, min, max);
    rvec_sub(max, min, size);
    printf("new size  : %6.3f %6.3f %6.3f\n", size[XX], size[YY], size[ZZ]);
  }
  
  /* calculate new boxsize */
  if (bDist) 
    for (i=0; (i<DIM); i++)
      box[i][i]=size[i]+2*dist;
  if (bSetSize)
    for (i=0; (i<DIM); i++)
      box[i][i]=newbox[i];
  if (bCubic) {
    d=box[XX][XX];
    for (i=1; (i<DIM); i++)
      if (box[i][i]>d) d=box[i][i];
    for (i=0; (i<DIM); i++)
      box[i][i]=d;
  }
  /* calculate new coords for geometrical center */
  if (!bSetCenter) 
    for (i=0; (i<DIM); i++)
      center[i]=box[i][i]/2;

  /* center molecule on 'center' */
  if (bCenter)
    center_conf(natom, x,center,gc);
    
  /* rotate it, if necessary */
  if (bRotate)
    simple_rotate_conf(natom,x,v,box,cRotate[0]);
    
  /* print some */
  if (bCenter || bScale)
    printf("new center: %6.3f %6.3f %6.3f\n", 
	   center[XX],center[YY],center[ZZ]);
  if ( bScale || bDist || bSetSize )
    printf("new box   : %6.3f %6.3f %6.3f\n", 
	   box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
  
  if (opt2bSet("-n",NFILE,fnm) || bNDEF) {
    get_index(&atoms,opt2fn_null("-n",NFILE,fnm),
	      1,&isize,&index,&groupnames);
    if (opt2bSet("-bf",NFILE,fnm))
      fatal_error(0,"combination not implemented: -bf -n  or -bf -ndef");
    else
      write_sto_conf_indexed(outfile,title,&atoms,x,bHaveV?v:NULL,box,
			     isize,index); 
  }
  else {
    if (outftp == efGRO) {
      out=ffopen(outfile,"w");
      write_hconf(out,title,&atoms,x,bHaveV?v:NULL,box); 
      fclose(out);
    } else {
      out=ffopen(outfile,"w");
      if (opt2bSet("-bf",NFILE,fnm)) {
	read_bfac(opt2fn("-bf",NFILE,fnm),&n_bfac,&bfac,&bfac_nr);
	set_pdb_conf_bfac(atoms.nr,atoms.nres,&atoms,x,box,
			  n_bfac,bfac,bfac_nr,peratom,bLegend);
      }
      if (opt2parg_bSet("-label",NPA,pa)) {
	for(i=0; (i<atoms.nr); i++) 
	  atoms.atom[i].chain=label[0];
      }
      write_pdbfile(out,title,&atoms,x,box,0,!bLegend);
      if (bLegend)
	pdb_legend(out,atoms.nr,atoms.nres,&atoms,x);
      fclose(out);
    }  
  }
  
  thanx(stdout);
  
  return 0;
}

