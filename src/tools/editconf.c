/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include <math.h>
#include <string.h>
#include "pdbio.h"
#include "confio.h"
#include "symtab.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "typedefs.h"
#include "gstat.h"

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

void calc_geom(int natom, rvec *pdba, rvec geom_center, rvec min, rvec max)
{
  int i,j;
  
  clear_rvec(geom_center);
  for (j=0; (j<DIM); j++)
    min[j]=max[j]=pdba[0][j];
  for (i=0; (i<natom); i++) {
    rvec_inc(geom_center,pdba[i]);
    for (j=0; (j<DIM); j++) {
      if (pdba[i][j] < min[j]) min[j]=pdba[i][j];
      if (pdba[i][j] > max[j]) max[j]=pdba[i][j];
    }
  }
  svmul(1./natom,geom_center,geom_center);
}

void center_conf(int natom, rvec *pdba, rvec center, rvec geom_cent)
{
  int       i;
  rvec shift;
  
  rvec_sub(center,geom_cent,shift);

  printf("shift     : %6.3f %6.3f %6.3f\n",
	 shift[XX],shift[YY],shift[ZZ]);

  for (i=0; (i<natom); i++) 
    rvec_inc(pdba[i], shift);
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

void print_pdbatoms_x(FILE *out,int natom,t_pdbatom pdba[])
{
  int i;
  char buf[12];
  
  for(i=0; (i<natom); i++) {
    sprintf(buf,"%s",pdba[i].resnm);
    buf[3]='\0';
    fprintf(out,pdbformat,
            pdbtp[pdba[i].pdbtp],pdba[i].atomnr + 1,pdba[i].atomnm,
            buf,pdba[i].chain,pdba[i].resnr + 1,
            10*pdba[i].x[XX],10*pdba[i].x[YY],10*pdba[i].x[ZZ],
            pdba[i].dummy,pdba[i].bfac);
  }
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

void write_pdb_conf_bfac(char *fn, t_atoms *atoms, rvec x[], matrix box,
			 int n_bfac,double *bfac,int *bfac_nr,
			 bool peratom, bool perres, bool bLegend)
{
  FILE *out;
  real bfac_min,bfac_max,xmin,ymin,zmin;
  int  i,n,natoms;
  bool found;
  t_pdbatom *pdba;
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
  
  pdba=atoms2pdba(atoms,x);  
  natoms=atoms->nr;
  for(i=0; (i<natoms); i++)
    pdba[i].bfac=0;
  
  if ( (n_bfac == atoms->nres) || (perres) ) {
    fprintf(stderr,"Will attach %d B-factors to %d residues\n",
	    n_bfac,atoms->nres);
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
		"residue) specified.",n_bfac,natoms,atoms->nres);
  
  out=ffopen(fn,"w");
  print_pdbatoms(out,natoms,pdba,box);
  if (bLegend) {
    fprintf(stderr,"B-factors range from %g to %g\n",bfac_min,bfac_max);
    xmin = 1e10;
    ymin = 1e10;
    zmin = 1e10;
    for (i=0; (i<atoms->nr); i++) {
      if (pdba[i].x[XX] < xmin) xmin = pdba[i].x[XX];
      if (pdba[i].x[YY] < ymin) ymin = pdba[i].x[YY];
      if (pdba[i].x[ZZ] < zmin) zmin = pdba[i].x[ZZ];
    }
    sprintf(buf,"%s","LEG");
    fprintf(out,"REMARK    NOW IT'S LEGEND TIME\n");
    buf[3]='\0';
    for (i=1; (i<12); i++) {
      fprintf(out,pdbformat,
	      "ATOM  ",atoms->nr+1+i,"CA",buf,' ',atoms->nres+1,
	      xmin+(i*1.2),ymin,zmin,1.0,
	      bfac_min+ ((i-1.0)*(bfac_max-bfac_min)/10) );
    }
  }
  fclose(out);
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
    "performed.[PAR]",
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
    "maximum value found, effectively making a legend for viewing."
  };
  static char *bugs[] = {
    "For complex molecules, the periodicity removal routine may break down,",
    "in that case you can use trjconv"
  };
  static real dist   = 0.0,rbox=0.0;
  static bool bNDEF=FALSE,bRMPBC=FALSE,bCenter=FALSE;
  static bool peratom=FALSE,perres=FALSE,bLegend=FALSE;
  static rvec scale={1.0,1.0,1.0},newbox={0.0,0.0,0.0};
  static rvec center={0.0,0.0,0.0};
  t_pargs pa[] = {
    { "-ndef", FALSE, etBOOL, &bNDEF, "Choose output from default index groups" },    
    { "-d", FALSE, etREAL, &dist, 
	"Distance between the solute and the rectangular box (default don't change box)" }, 
    { "-dc", FALSE, etREAL, &dist,
      "Distance between the solute and the cubic box (default don't change box)" },
    { "-b", FALSE, etREAL,  &rbox,
	  "size of the cubic box (default don't change box)" },
    { "-center", FALSE, etBOOL, &bCenter,
	"Center molecule in box (implied by -d -dc -b)" },
    { "-cx", FALSE, etREAL, &center[XX], "x coordinate of geometrical center" },
    { "-cy", FALSE, etREAL, &center[YY], "y coordinate of geometrical center" },
    { "-cz", FALSE, etREAL, &center[ZZ], "z coordinate of geometrical center" },
    { "-bx", FALSE, etREAL, &newbox[XX], "x size of box" },
    { "-by", FALSE, etREAL, &newbox[YY], "y size of box" },
    { "-bz", FALSE, etREAL, &newbox[ZZ], "z size of box" },
    { "-sx", FALSE, etREAL, &scale[XX], "Scale factor for x coordinate" },
    { "-sy", FALSE, etREAL, &scale[YY], "Scale factor for y coordinate" },
    { "-sz", FALSE, etREAL, &scale[ZZ], "Scale factor for z coordinate" },
    { "-pbc",  FALSE, etBOOL, &bRMPBC, "Remove the periodicity (make molecule whole again)" },
    { "-atom", FALSE, etBOOL, &peratom, "Attach B-factors per atom" },
    { "-res",  FALSE, etBOOL, &perres,  "Attach B-factors per residue" },
    { "-legend",FALSE,etBOOL, &bLegend, "Make B-factor legend"}
  };
#define NPA asize(pa)

  FILE      *out;
  char      title[STRLEN];
  int       natom,i,n_bfac;
  double    *bfac=NULL;
  int       *bfac_nr=NULL;
  t_atoms   atoms;
  char      *groupnames;
  int       isize;
  atom_id   *index;
  rvec      *x,*v,gc,min,max,size;
  matrix    box;
  bool      bSetSizeAll,bSetSize[DIM],bCubic,bDist,bSetCenter[DIM],bScale;
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
  bScale = ( opt2parg_bSet("-sx",NPA,pa) || 
	     opt2parg_bSet("-sy",NPA,pa) || 
	     opt2parg_bSet("-sz",NPA,pa) );
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
      write_hconf_indexed(out,cool_quote(),&atoms,isize,index,x,v,box); 
      fclose(out);
    }
    if (opt2bSet("-op",NFILE,fnm)) {
      if (opt2bSet("-bf",NFILE,fnm)) {
	fatal_error(0,"combination not implemented: -bf -n  or -bf -ndef");
	/*
	read_bfac(opt2fn("-bf",NFILE,fnm),&n_bfac,&bfac,&bfac_nr);
	write_pdb_conf_bfac(opt2fn("-op",NFILE,fnm),
			    &atoms,x,box,n_bfac,bfac,bfac_nr,
			    peratom,perres,bLegend);
	*/
      } else {
	write_pdb_conf_indexed(opt2fn("-op",NFILE,fnm),
			       &atoms,x,box,isize,index); 
      }
    }
  }
  else {
    if (opt2bSet("-o",NFILE,fnm) || !opt2bSet("-op",NFILE,fnm)) {
      out=opt2FILE("-o",NFILE,fnm,"w");
      write_hconf(out,cool_quote(),&atoms,x,v,box); 
      fclose(out);
    }
    if (opt2bSet("-op",NFILE,fnm)) {
      if (opt2bSet("-bf",NFILE,fnm)) {
	read_bfac(opt2fn("-bf",NFILE,fnm),&n_bfac,&bfac,&bfac_nr);
	write_pdb_conf_bfac(opt2fn("-op",NFILE,fnm),
			    &atoms,x,box,n_bfac,bfac,bfac_nr,
			    peratom,perres,bLegend);
      } else {
	write_pdb_conf(opt2fn("-op",NFILE,fnm),&atoms,x,box,FALSE); 
      }
    }  
  }
  
  thanx(stdout);
  
  return 0;
}

