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

void set_pdb_conf_bfac(int natoms,int nres,t_pdbatom pdba[],rvec x[],
		       matrix box,int n_bfac,double *bfac,int *bfac_nr,
		       bool peratom, bool perres, bool bLegend)
{
  FILE *out;
  real bfac_min,bfac_max;
  int  i,n;
  bool found;
  char buf[120];

  bfac_max=-1e10;
  bfac_min=1e10;
  for(i=0; (i<n_bfac); i++) {
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
    pdba[i].bfac=0;
  
  if ( (n_bfac == nres) || (perres) ) {
    fprintf(stderr,"Will attach %d B-factors to %d residues\n",
	    n_bfac,nres);
    for(i=0; (i<n_bfac); i++) {
      found=FALSE;
      for(n=0; (n<natoms); n++)
	if ( bfac_nr[i] == (pdba[n].resnr+1) ) {
	  pdba[n].bfac=bfac[i];
	  found=TRUE;
	}
      if (!found) {
	sprintf(buf,"Residue nr %d not found\n",bfac_nr[i]);
	warning(buf);
      }
    }
  } else if ( (n_bfac == natoms) || (peratom) ){
    fprintf(stderr,"Will attach %d B-factors to %d atoms\n",n_bfac,natoms);
    for(i=0; (i<n_bfac); i++) {
      found=FALSE;
      for(n=0; (n<natoms); n++)
	if ( bfac_nr[i] == pdba[n].atomnr ) {
	  pdba[n].bfac=bfac[i];
	  found=TRUE;
	}
      if (!found) {
	sprintf(buf,"Atom nr %d not found\n",bfac_nr[i]);
	warning(buf);
      }
    }
  } else
    fatal_error(0,"Number of B-factors (%d) does not match number of atoms "
		"(%d) or residues (%d) and no attachment type (atom or "
		"residue) specified.",n_bfac,natoms,nres);
}

void pdb_legend(FILE *out,int natoms,int nres,t_pdbatom pdba[])
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
    xmin     = min(xmin,pdba[i].x[XX]);
    ymin     = min(ymin,pdba[i].x[YY]);
    zmin     = min(zmin,pdba[i].x[ZZ]);
    bfac_min = min(bfac_min,pdba[i].bfac);
    bfac_max = max(bfac_max,pdba[i].bfac);
  }
  fprintf(stderr,"B-factors range from %g to %g\n",bfac_min,bfac_max);
  sprintf(buf,"%s","LEG");
  buf[3]='\0';
  for (i=1; (i<12); i++) {
    fprintf(out,pdbformat,
	    "ATOM  ",natoms+1+i,"CA",buf,' ',nres+1,
	    xmin+(i*1.2),ymin,zmin,1.0,
	    bfac_min+ ((i-1.0)*(bfac_max-bfac_min)/10) );
  }
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "editconf converts generic structure format to [TT].gro[tt] or",
    "[TT].pdb[tt].[PAR]",
    "A number of options is present to modify the coordinates",
    "and box. [TT]-d[tt], [TT]-dc[tt] and [TT]-b[tt] modify the box and",
    "center the coordinates relative to the new box.",
    "[TT]-dc[tt] takes precedent over [TT]-d[tt]. [TT]-b[tt]",
    "takes precedent over [TT]-dc[tt] and [TT]-d[tt].",
    "[TT]-bx[tt], [TT]-by[tt], [TT]-bz[tt], each override [TT]-b[tt]",
    "for one coordinate.",
    "[TT]-cx[tt], [TT]-cy[tt], [TT]-cz[tt] override the center for",
    "one coordinate.[PAR]",
    "Scaling is applied before any of the other operations are",
    "performed. Boxes can be scaled to give a certain density (option",
    "[TT]-rho[tt][PAR]",
    "Groups are selected after all operations have been applied.[PAR]",
    "Periodicity can be removed in a crude manner.",
    "It is important that the box sizes at the bottom of your input file",
    "are correct when the periodicity is to be removed.[PAR]",
    "When writing [TT].pdb[tt] files, B-factors can be",
    "added per atom or per residue. B-factors are read",
    "from a file with with following format: first line states number of",
    "entries in the file, next lines state either residue or atom",
    "number followed by B-factor. Obiously, any type of numeric data can",
    "be displayed in stead of B-factors. [TT]-legend[tt] will produce",
    "a row of CA atoms with B-factors ranging from the minimum to the",
    "maximum value found, effectively making a legend for viewing.[PAR]",
    "Finally with option [TT]-label[tt] editconf can add a chain identifier",
    "to a pdb file, which can be useful for analysis using e.g. rasmol."
  };
  static char *bugs[] = {
    "For complex molecules, the periodicity removal routine may break down,",
    "in that case you can use trjconv"
  };
  static real dist   = 0.0,rbox=0.0;
  static bool bNDEF=FALSE,bRMPBC=FALSE,bCenter=FALSE;
  static bool peratom=FALSE,perres=FALSE,bLegend=FALSE;
  static rvec scale={1.0,1.0,1.0},newbox={0.0,0.0,0.0};
  static real rho=1000.0;
  static rvec center={0.0,0.0,0.0};
  static char *label="A";
  
  t_pargs pa[] = {
    { "-ndef", FALSE, etBOOL, &bNDEF, 
      "Choose output from default index groups" },    
    { "-d", FALSE, etREAL, &dist, 
      "Distance between the solute and the rectangular box "
      "(default don't change box)" }, 
    { "-dc", FALSE, etREAL, &dist,
      "Distance between the solute and the cubic box "
      "(default don't change box)" },
    { "-b", FALSE, etREAL,  &rbox,
      "size of the cubic box (default don't change box)" },
    { "-center", FALSE, etBOOL, &bCenter,
      "Center molecule in box (implied by -d -dc -b)" },
    { "-cx", FALSE, etREAL, &center[XX], "x coordinate of geometrical center"},
    { "-cy", FALSE, etREAL, &center[YY], "y coordinate of geometrical center"},
    { "-cz", FALSE, etREAL, &center[ZZ], "z coordinate of geometrical center"},
    { "-bx", FALSE, etREAL, &newbox[XX], "x size of box" },
    { "-by", FALSE, etREAL, &newbox[YY], "y size of box" },
    { "-bz", FALSE, etREAL, &newbox[ZZ], "z size of box" },
    { "-sx", FALSE, etREAL, &scale[XX], "Scale factor for x coordinate" },
    { "-sy", FALSE, etREAL, &scale[YY], "Scale factor for y coordinate" },
    { "-sz", FALSE, etREAL, &scale[ZZ], "Scale factor for z coordinate" },
    { "-rho",FALSE, etREAL, &rho, 
      "Density (g/l) of the output box achieved by scaling" },
    { "-pbc",  FALSE, etBOOL, &bRMPBC, 
      "Remove the periodicity (make molecule whole again)" },
    { "-atom", FALSE, etBOOL, &peratom, "Attach B-factors per atom" },
    { "-res",  FALSE, etBOOL, &perres,  "Attach B-factors per residue" },
    { "-legend",FALSE,etBOOL, &bLegend, "Make B-factor legend" },
    { "-label", FALSE, etSTR, &label,   "Add chain label for all residues" }
  };
#define NPA asize(pa)

  FILE      *out;
  char      title[STRLEN];
  int       natom,i,j,n_bfac;
  double    *bfac=NULL;
  int       *bfac_nr=NULL;
  t_atoms   atoms;
  char      *groupnames;
  int       isize;
  atom_id   *index;
  rvec      *x,*v,gc,min,max,size;
  matrix    box;
  bool      bSetSizeAll,bSetSize[DIM],bCubic,bDist,bSetCenter[DIM];
  bool      bHaveV,bScale,bRho;
  real      xs,ys,zs,xcent,ycent,zcent,d;
  t_filenm fnm[] = {
    { efSTX, "-f", NULL, ffREAD },
    { efNDX, "-n", NULL, ffOPTRD },
    { efGRO, "-o", "out", ffOPTWR },
    { efPDB, "-op", NULL, ffOPTWR },
    { efDAT, "-bf", "bfact", ffOPTRD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,NPA,pa,
		    asize(desc),desc,asize(bugs),bugs);

#define ANY(s) (s[XX] || s[YY] || s[ZZ])
#define ALL(s) (s[XX] && s[YY] && s[ZZ])
  bSetSizeAll =opt2parg_bSet("-b", NPA,pa);
  bSetSize[XX]=opt2parg_bSet("-bx",NPA,pa);
  bSetSize[YY]=opt2parg_bSet("-by",NPA,pa);
  bSetSize[ZZ]=opt2parg_bSet("-bz",NPA,pa);
  bSetCenter[XX]=opt2parg_bSet("-cx",NPA,pa);
  bSetCenter[YY]=opt2parg_bSet("-cy",NPA,pa);
  bSetCenter[ZZ]=opt2parg_bSet("-cz",NPA,pa);
  bCubic  = opt2parg_bSet("-dc",NPA,pa);
  bDist   = ( bCubic || opt2parg_bSet("-d",NPA,pa) );
  bCenter = ( bCenter || bDist || bSetSizeAll || ANY(bSetCenter) );
  bScale  = ( opt2parg_bSet("-sx",NPA,pa) || 
	     opt2parg_bSet("-sy",NPA,pa) || 
	     opt2parg_bSet("-sz",NPA,pa) );
  bRho    =  opt2parg_bSet("-rho",NPA,pa);
  if (bScale && bRho)
    fprintf(stderr,"WARNING: setting -rho overrides -sx, -sy and -sz");
  bScale  = bScale || bRho;
  
  /* set newbox size */
  if (bSetSizeAll) {
    for (i=0; (i<DIM); i++)
      if (!bSetSize[i]) {
	newbox[i]=rbox;
	bSetSize[i]=TRUE;
      }
  }
  
  atoms.nr=0;
  atoms.nres=0;
  get_stx_coordnum(fnm[0].fn,&natom);
  snew(atoms.atomname,natom);
  snew(atoms.resname,natom);
  snew(atoms.atom,natom);
  snew(x,natom);
  snew(v,natom);
  read_stx_conf(fnm[0].fn,title,&atoms,x,v,box);

  bHaveV=FALSE;
  for (i=0; (i<natom) && !bHaveV; i++)
    for (j=0; (j<DIM) && !bHaveV; j++)
      bHaveV=bHaveV || (v[i][j]!=0);
  printf("%svelocities found\n",bHaveV?"":"No ");
  
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
  for (i=0; (i<DIM); i++)
    if (bSetSize[i])
      box[i][i]=newbox[i];
  if (bCubic) {
    d=box[XX][XX];
    for (i=1; (i<DIM); i++)
      if (box[i][i]>d) d=box[i][i];
    for (i=0; (i<DIM); i++)
      box[i][i]=d;
  }
  /* calculate new coords for geometrical center */
  for (i=0; (i<DIM); i++)
    if (!bSetCenter[i]) 
      center[i]=box[i][i]/2;
  

  /* center molecule on 'center' */
  if (bCenter)
    center_conf(natom, x,center,gc);
  /* print some */
  if (bCenter || bScale)
    printf("new center: %6.3f %6.3f %6.3f\n", 
	   center[XX],center[YY],center[ZZ]);
  if ( bScale || bDist || ANY(bSetSize) )
    printf("new box   : %6.3f %6.3f %6.3f\n", 
	   box[XX][XX], box[YY][YY], box[ZZ][ZZ]);

  if (opt2bSet("-n",NFILE,fnm) || bNDEF) {
    get_index(&atoms,opt2fn_null("-n",NFILE,fnm),
	    1,&isize,&index,&groupnames);
    if (opt2bSet("-o",NFILE,fnm) || !opt2bSet("-op",NFILE,fnm)) {
      out=opt2FILE("-o",NFILE,fnm,"w");
      write_hconf_indexed(out,title,&atoms,isize,index,x,bHaveV?v:NULL,box); 
      fclose(out);
    }
    if (opt2bSet("-op",NFILE,fnm)) {
      if (opt2bSet("-bf",NFILE,fnm)) {
	fatal_error(0,"combination not implemented: -bf -n  or -bf -ndef");
      } else {
	write_pdb_conf_indexed(opt2fn("-op",NFILE,fnm),title,
			       &atoms,x,box,isize,index); 
      }
    }
  }
  else {
    if (opt2bSet("-o",NFILE,fnm) || !opt2bSet("-op",NFILE,fnm)) {
      out=opt2FILE("-o",NFILE,fnm,"w");
      write_hconf(out,title,&atoms,x,bHaveV?v:NULL,box); 
      fclose(out);
    }
    if (opt2bSet("-op",NFILE,fnm)) {
      FILE *out;
      t_pdbatom *pdba = atoms2pdba(&atoms,x);
      
      if (opt2bSet("-bf",NFILE,fnm)) {
	read_bfac(opt2fn("-bf",NFILE,fnm),&n_bfac,&bfac,&bfac_nr);
	set_pdb_conf_bfac(atoms.nr,atoms.nres,pdba,x,box,
			  n_bfac,bfac,bfac_nr,peratom,perres,bLegend);
      }
      if (opt2parg_bSet("-label",NPA,pa)) {
	for(i=0; (i<atoms.nr); i++) 
	  pdba[i].chain=label[0];
      }
      
      out=opt2FILE("-op",NFILE,fnm,"w");
      print_pdbatoms(out,title,atoms.nr,pdba,box);
      if (bLegend)
	pdb_legend(out,atoms.nr,atoms.nres,pdba);
      fclose(out);
      sfree(pdba);
    }  
  }
  
  thanx(stdout);
  
  return 0;
}

