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
#include "gbutil.h"
#include "strdb.h"
#include "rdgroup.h"
#include "physics.h"
#include "atomprop.h"
#include "addconf.h"
#include "tpxio.h"
#include "pbc.h"

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

real calc_mass(t_atoms *atoms,bool bGetMass)
{
  real tmass;
  int i;

  tmass = 0;
  for(i=0; (i<atoms->nr); i++) {
    if (bGetMass)
      atoms->atom[i].m = get_mass(*atoms->resname[atoms->atom[i].resnr], 
				  *atoms->atomname[i]);
    tmass += atoms->atom[i].m;
  }

  return tmass;
}

real calc_geom(int isize,atom_id *index,
	       rvec *x, rvec geom_center, rvec min, rvec max,bool bDiam)
{
  real    diam2,d;
  char    *grpnames;
  int     ii,i,j;
  
  clear_rvec(geom_center);
  diam2 = 0;
  if (isize == 0) {
    clear_rvec(min);
    clear_rvec(max);
  } else {
    if (index)
      ii=index[0];
    else
      ii=0;
    for (j=0; j<DIM; j++)
      min[j]=max[j]=x[ii][j];
    for (i=0; i<isize; i++) {
      if (index)
      ii = index[i];
      else
	ii = i;
      rvec_inc(geom_center,x[ii]);
      for (j=0; j<DIM; j++) {
	if (x[ii][j] < min[j]) min[j]=x[ii][j];
	if (x[ii][j] > max[j]) max[j]=x[ii][j];
      }
      if (bDiam) {
	if (index) 
	for (j=i+1; j<isize; j++) {
	  d = distance2(x[ii],x[index[j]]);
	  diam2 = max(d,diam2);
	}
	else
	  for (j=i+1; j<isize; j++) {
	  d = distance2(x[i],x[j]);
	  diam2 = max(d,diam2);
	  }
      }
    }
    svmul(1.0/isize,geom_center,geom_center);
  }
  
  return sqrt(diam2);
}

void center_conf(int natom, rvec *x, rvec center, rvec geom_cent)
{
  int       i;
  rvec shift;
  
  rvec_sub(center,geom_cent,shift);

  printf("    shift       :%7.3f%7.3f%7.3f (nm)\n",
	 shift[XX],shift[YY],shift[ZZ]);

  for (i=0; (i<natom); i++) 
    rvec_inc(x[i], shift);
}

void scale_conf(int natom,rvec x[],matrix box,rvec scale)
{
  int i,j;
  
  for(i=0; i<natom; i++) {
    for (j=0; j<DIM; j++)
      x[i][j] *= scale[j];
  }
  for (i=0; i<DIM; i++)
    for (j=0; j<DIM; j++)
      box[i][j] *= scale[j];
}

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,m,d;
  
  /* check periodic boundary */
  for(n=1;(n<atoms->nr);n++) {
    for(m=DIM-1; m>=0; m--) {
      dist = x[n][m]-x[n-1][m];
      if (fabs(dist) > 0.9*box[m][m]) { 
	if ( dist >  0 )
	  for(d=0; d<=m; d++)
	    x[n][d] -= box[m][d];
	else
	  for(d=0; d<=m; d++)
	    x[n][d] += box[m][d];
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

void set_pdb_conf_bfac(int natoms,int nres,t_atoms *atoms,
		       int n_bfac,double *bfac,int *bfac_nr,
		       bool peratom)
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
  while ( (fabs(bfac_max) < 0.5) && (fabs(bfac_min) < 0.5) ) {
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
void visualize_images(char *fn,matrix box)
{
  t_atoms atoms;
  rvec    *img;
  char    *c,*ala;
  int     nat,i;

  nat = NTRICIMG+1;
  init_t_atoms(&atoms,nat,FALSE);
  atoms.nr = nat;
  snew(img,nat);
  c = "C";
  ala = "ALA";
  for(i=0; i<nat; i++) {
    atoms.atomname[i] = &c;
    atoms.atom[i].resnr = i;
    atoms.resname[i] = &ala;
    atoms.atom[i].chain = 'A'+i/NCUCVERT;
  }
  calc_triclinic_images(box,img+1);

  write_sto_conf(fn,"Images",&atoms,img,NULL,box); 

  free_t_atoms(&atoms);
  sfree(img);
}


void visualize_box(FILE *out,int a0,int r0,matrix box,rvec gridsize)
{
  int     *edge;
  rvec    *vert,shift;
  int     nx,ny,nz,nbox,nat;
  int     i,j,x,y,z;
  int     rectedge[24] = { 0,1, 1,3, 3,2, 0,2, 0,4, 1,5, 3,7, 2,6, 4,5, 5,7, 7,6, 6,4 };

  a0++;
  r0++;
  
  nx = (int)(gridsize[XX]+0.5);
  ny = (int)(gridsize[YY]+0.5);
  nz = (int)(gridsize[ZZ]+0.5);
  nbox = nx*ny*nz;
  if (TRICLINIC(box)) {
    nat = nbox*NCUCVERT;
    snew(vert,nat);
    calc_compact_unitcell_vertices(box,vert);
    j = 0;
    for(z=0; z<nz; z++)
      for(y=0; y<ny; y++)
	for(x=0; x<nx; x++) {
	  for(i=0; i<DIM; i++)
	    shift[i] = x*box[0][i]+y*box[1][i]+z*box[2][i];
	  for(i=0; i<NCUCVERT; i++) {
	    rvec_add(vert[i],shift,vert[j]);
	    j++;
	  }
	}
    
    for(i=0; i<nat; i++)
      fprintf(out,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f\n",
	      "ATOM",a0+i,"C","BOX",'K'+i/NCUCVERT,r0+i,
	      10*vert[i][XX],10*vert[i][YY],10*vert[i][ZZ]);
    
    edge = compact_unitcell_edges();
    for(j=0; j<nbox; j++)
      for(i=0; i<NCUCEDGE; i++)
	fprintf(out,"CONECT%5d%5d\n",
		a0 + j*NCUCVERT + edge[2*i],
		a0 + j*NCUCVERT + edge[2*i+1]);
    
    sfree(vert);
  } else {
    i=0;
    for(z=0; z<=1; z++)
      for(y=0; y<=1; y++)
	for(x=0; x<=1; x++) {
	  fprintf(out,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f\n",
	      "ATOM",a0+i,"C","BOX",'K'+i/8,r0+i,
	      x*10*box[XX][XX],y*10*box[YY][YY],z*10*box[ZZ][ZZ]);
	  i++;
	}
    for(i=0; i<24; i+=2)
      fprintf(out,"CONECT%5d%5d\n",a0+rectedge[i],a0+rectedge[i+1]);
  }
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "editconf converts generic structure format to [TT].gro[tt], [TT].g96[tt]",
    "or [TT].pdb[tt].[PAR]",
    "The box can be modified with options [TT]-box[tt] and [TT]-d[tt], both",
    "will center the system in the box.",
    "Option [TT]-bt[tt] determines the box type: [TT]rect[tt] is a",
    "rectangular box, [TT]cubic[tt] is a cubic box, [TT]dodecahedron[tt] is",
    "a rhombic dodecahedron and [TT]truncoct[tt] is a truncated octahedron.",
    "The last two are special cases of a triclinic box.",
    "The length of the three box vectors of the truncated octahedron is the",
    "shortest distance between two opposite hexagons.",
    "The volume of a dodecahedron is 0.71 and that of a truncated octahedron",
    "is 0.77 of that of a cubic box with the same periodic image distance.",
    "Option [TT]-box[tt] requires only",
    "one value for a cubic box, dodecahedron and a truncated octahedron.",
    "With [TT]-d[tt] and [TT]rect[tt] the size of the system in the x, y",
    "and z directions is used. With [TT]-d[tt] and [TT]cubic[tt],",
    "[TT]dodecahedron[tt] or [TT]truncoct[tt] the diameter of the system",
    "is used, which is the largest distance between two atoms.[PAR]",
    "When [TT]-n[tt] or [TT]-ndef[tt] is set, a group",
    "can be selected for calculating the size and the geometric center,",
    "otherwise the whole system is used.[PAR]",
    "[TT]-rotate[tt] rotates the coordinates and velocities.",
    "[TT]-princ[tt] aligns the principal axes of the system along the",
    "coordinate axes, this may allow you to decrease the box volume,",
    "but beware that molecules can rotate significantly in a nanosecond.[PAR]",
    "Scaling is applied before any of the other operations are",
    "performed. Boxes can be scaled to give a certain density (option",
    "[TT]-density[tt]). A special feature of the scaling option, when the",
    "factor -1 is given in one dimension, one obtains a mirror image,",
    "mirrored in one of the plains, when one uses -1 in three dimensions",
    "a point-mirror image is obtained.[PAR]",
    "Groups are selected after all operations have been applied.[PAR]",
    "Periodicity can be removed in a crude manner.",
    "It is important that the box sizes at the bottom of your input file",
    "are correct when the periodicity is to be removed.[PAR]",
    "The program can optionally rotate the solute molecule to align the",
    "molecule along its principal axes ([TT]-rotate[tt])[PAR]",
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
    "With the option -mead a special pdb file for the MEAD electrostatics",
    "program (Poisson-Boltzmann solver) can be made. A further prerequisite",
    "is that the input file is a run input file.",
    "The B-factor field is then filled with the Van der Waals radius",
    "of the atoms while the occupancy field will hold the charge.[PAR]",
    "Finally with option [TT]-label[tt] editconf can add a chain identifier",
    "to a pdb file, which can be useful for analysis with e.g. rasmol."
  };
  static char *bugs[] = {
    "For complex molecules, the periodicity removal routine may break down, "
    "in that case you can use trjconv"
  };
  static real dist=0.0,rbox=0.0,to_diam=0.0;
  static bool bNDEF=FALSE,bRMPBC=FALSE,bCenter=FALSE;
  static bool peratom=FALSE,bLegend=FALSE,bOrient=FALSE,bMead=FALSE;
  static rvec scale={1.0,1.0,1.0},newbox={0.0,0.0,0.0};
  static real rho=1000.0,rvdw=0.12;
  static rvec center={0.0,0.0,0.0},rotangles={0.0,0.0,0.0};
  static char *btype[]={ NULL, "rect", "cubic", "dodecahedron", "truncoct", NULL },*label="A";
  static rvec visbox={0,0,0};
  t_pargs pa[] = {
    { "-ndef",   FALSE, etBOOL, {&bNDEF}, 
      "Choose output from default index groups" },
    { "-visbox",    FALSE, etRVEC, {visbox}, 
      "HIDDENVisualize a grid of boxes, -1 visualizes the 14 box images" },
    { "-bt",   FALSE, etENUM, {btype}, 
      "Box type for -box and -d" },
    { "-box",    FALSE, etRVEC, {&newbox}, "Box vector lengths" },
    { "-d",      FALSE, etREAL, {&dist}, 
      "Distance between the solute and the box" },
    { "-c",      FALSE, etBOOL, {&bCenter},
      "Center molecule in box (implied by -box and -d)" },
    { "-center", FALSE, etRVEC, {&center}, "Coordinates of geometrical center"},
    { "-rotate", FALSE, etRVEC, {rotangles},
      "Rotation around the X, Y and Z axes in degrees" },
    { "-princ",  FALSE, etBOOL, {&bOrient}, "Orient molecule(s) along their principal axes" },
    { "-scale",  FALSE, etRVEC, {&scale}, "Scaling factor" },
    { "-density",FALSE, etREAL, {&rho}, 
      "Density (g/l) of the output box achieved by scaling" },
    { "-pbc",    FALSE, etBOOL, {&bRMPBC}, 
      "Remove the periodicity (make molecule whole again)" },
    { "-mead",   FALSE, etBOOL, {&bMead},
      "Store the charge of the atom in the occupancy field and the radius of the atom in the B-factor field" },
    { "-rvdw",   FALSE, etREAL, {&rvdw},
      "Default Van der Waals radius if one can not be found in the database" },
    { "-atom",   FALSE, etBOOL, {&peratom}, "Force B-factor attachment per atom" },
    { "-legend", FALSE, etBOOL, {&bLegend}, "Make B-factor legend" },
    { "-label",  FALSE, etSTR,  {&label},   "Add chain label for all residues" }
  };
#define NPA asize(pa)

  FILE       *out;
  char       *infile,*outfile,title[STRLEN];
  int        outftp,natom,i,j,n_bfac;
  double     *bfac=NULL;
  int        *bfac_nr=NULL;
  t_topology *top;
  t_atoms    atoms;
  char       *grpname,*sgrpname;
  int        isize,ssize;
  atom_id    *index,*sindex;
  rvec       *x,*v,gc,min,max,size;
  matrix     box;
  bool       bIndex,bSetSize,bCubic,bDist,bSetCenter;
  bool       bHaveV,bScale,bRho,bRotate,bCalcGeom,bCalcDiam;
  real       xs,ys,zs,xcent,ycent,zcent,diam,d;
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

  bIndex    = opt2bSet("-n",NFILE,fnm) || bNDEF;
  bSetSize  = opt2parg_bSet("-box" ,NPA,pa);
  bSetCenter= opt2parg_bSet("-center" ,NPA,pa);
  bDist     = opt2parg_bSet("-d" ,NPA,pa);
  bCenter   = bCenter || bDist || bSetCenter || bSetSize;
  bScale    = opt2parg_bSet("-scale" ,NPA,pa);
  bRho      = opt2parg_bSet("-density",NPA,pa);
  bRotate   = opt2parg_bSet("-rotate",NPA,pa);
  if (bScale && bRho)
    fprintf(stderr,"WARNING: setting -density overrides -scale");
  bScale    = bScale || bRho;
  bCalcGeom = bCenter || bRotate || bOrient || bScale;
  bCalcDiam = btype[0][0]=='c' || btype[0][0]=='d' || btype[0][0]=='t';
  
  infile  = ftp2fn(efSTX,NFILE,fnm);
  outfile = ftp2fn(efSTO,NFILE,fnm);
  outftp  = fn2ftp(outfile);
  if (bMead && (outftp != efPDB))
    fatal_error(0,"Output file should be a .pdb file"
		" when using the -mead option\n");
  if (bMead && !((fn2ftp(infile) == efTPR) || 
		 (fn2ftp(infile) == efTPA) ||
		 (fn2ftp(infile) == efTPB)))
    fatal_error(0,"Input file should be a .tp[abr] file"
		" when using the -mead option\n");
  
  get_stx_coordnum(infile,&natom);
  init_t_atoms(&atoms,natom,TRUE);
  snew(x,natom);
  snew(v,natom);
  read_stx_conf(infile,title,&atoms,x,v,box);
  printf("Read %d atoms\n",atoms.nr); 

  if (bMead) {
    top = read_top(infile);
    if (atoms.nr != top->atoms.nr)
      fatal_error(0,"Atom numbers don't match (%d vs. %d)",
		  atoms.nr,top->atoms.nr);
    for(i=0; (i<atoms.nr); i++) {
      atoms.pdbinfo[i].occup = top->atoms.atom[i].q;
      /* Factor of 10 for Angstroms */
      atoms.pdbinfo[i].bfac  = 
	10*get_vdw(*top->atoms.resname[top->atoms.atom[i].resnr],
		   *top->atoms.atomname[i],rvdw);
    }
  }
  bHaveV=FALSE;
  for (i=0; (i<natom) && !bHaveV; i++)
    for (j=0; (j<DIM) && !bHaveV; j++)
      bHaveV=bHaveV || (v[i][j]!=0);
  printf("%selocities found\n",bHaveV?"V":"No v");

  if (visbox[0] > 0) {
    if (bIndex)
      fatal_error(0,"Sorry, can not visualize box with index groups");
    if (outftp != efPDB)
      fatal_error(0,"Sorry, can only visualize box with a pdb file");
  } else if (visbox[0] == -1)
    visualize_images("images.pdb",box);

  /* remove pbc */
  if (bRMPBC) 
    rm_gropbc(&atoms,x,box);

  if (bCalcGeom) {
    if (bIndex) {
      fprintf(stderr,"\nSelect a group for determining the system size:\n");
      get_index(&atoms,ftp2fn_null(efNDX,NFILE,fnm),
		1,&ssize,&sindex,&sgrpname);
    } else {
      ssize = atoms.nr;
      sindex = NULL;
    }
    diam=calc_geom(ssize,sindex,x,gc,min,max,bCalcDiam);
    rvec_sub(max, min, size);
    printf("    system size :%7.3f%7.3f%7.3f (nm)\n",
	   size[XX], size[YY], size[ZZ]);
    if (bCalcDiam)
      printf("    diameter    :%7.3f               (nm)\n",diam);
    printf("    center      :%7.3f%7.3f%7.3f (nm)\n", gc[XX], gc[YY], gc[ZZ]);
    printf("    box vectors :%7.3f%7.3f%7.3f (nm)\n", 
	   norm(box[XX]), norm(box[YY]), norm(box[ZZ]));
    printf("    box angles  :%7.2f%7.2f%7.2f (degrees)\n",
	   RAD2DEG*acos(cos_angle_no_table(box[YY],box[ZZ])),
	   RAD2DEG*acos(cos_angle_no_table(box[XX],box[ZZ])),
	   RAD2DEG*acos(cos_angle_no_table(box[XX],box[YY])));
    printf("    box volume  :%7.2f               (nm^3)\n",det(box));
  }
  
  if (bOrient)
    /* Orient the principal axes along the coordinate axes */
    orient_mol(&atoms,ftp2fn_null(efNDX,NFILE,fnm),x,bHaveV ? v : NULL);
  
  if ( bScale ) {
    /* scale coordinates and box */
    if (bRho) {
      /* Compute scaling constant */
      real vol,mass,dens;
      
      vol = det(box);
      mass = calc_mass(&atoms,!fn2bTPX(infile));
      dens = (mass*AMU)/(vol*NANO*NANO*NANO);
      fprintf(stderr,"Volume  of input %g (nm^3)\n",vol);
      fprintf(stderr,"Mass    of input %g (a.m.u.)\n",mass);
      fprintf(stderr,"Density of input %g (g/l)\n",dens);
      if (vol==0.0)
	fatal_error(0,"Cannot scale density with zero box\n");
      if (mass==0.0)
	fatal_error(0,"Cannot scale density with zero mass\n");
      
      scale[XX] = scale[YY] = scale[ZZ] = pow(dens/rho,1.0/3.0);
      fprintf(stderr,"Scaling all box vectors by %g\n",scale[XX]);
    }
    scale_conf(atoms.nr,x,box,scale);
  }
  
  if (bRotate) {
    /* Rotate */
    printf("Rotating %g, %g, %g degrees around the X, Y and Z axis respectively\n",rotangles[XX],rotangles[YY],rotangles[ZZ]);
    for(i=0; i<DIM; i++)
      rotangles[i] *= DEG2RAD;
    rotate_conf(natom,x,v,rotangles[XX],rotangles[YY],rotangles[ZZ]);
  }
  
  if (bCalcGeom) {
    /* recalc geometrical center and max and min coordinates and size */
    calc_geom(ssize,sindex,x,gc,min,max,FALSE);
    rvec_sub(max, min, size);
    if (bScale || bOrient || bRotate)
      printf("new system size : %6.3f %6.3f %6.3f\n",
	     size[XX],size[YY],size[ZZ]);
  }

  if (bSetSize || bDist) {
    clear_mat(box);
    /* calculate new boxsize */
    switch(btype[0][0]){
    case 'r':
      if (bSetSize)
	for (i=0; (i<DIM); i++)
	  box[i][i]=newbox[i];
      else 
	for (i=0; (i<DIM); i++)
	  box[i][i]=size[i]+2*dist;
      break;
    case 'c':
    case 'd':
    case 't':
      if (bSetSize)
	d = newbox[0];
      else
	d = diam+2*dist;
      if (btype[0][0] == 'c')
	for(i=0; i<DIM; i++)
	  box[i][i] = d;
      else if (btype[0][0] == 'd') {
	box[XX][XX] = d;
	box[YY][YY] = d;
	box[ZZ][XX] = d/2;
	box[ZZ][YY] = d/2;
	box[ZZ][ZZ] = d*sqrt(2)/2;
      } else {
	box[XX][XX] = d;
	box[YY][XX] = d/3;
	box[YY][YY] = d*sqrt(2)*2/3;
	box[ZZ][XX] = -d/3;
	box[ZZ][YY] = d*sqrt(2)/3;
	box[ZZ][ZZ] = d*sqrt(6)/3;
      }
     break;
    } 
  }

  /* calculate new coords for geometrical center */
  if (!bSetCenter)
    calc_box_center(box,center);

  /* center molecule on 'center' */
  if (bCenter)
    center_conf(natom,x,center,gc);
    
  /* print some */
  if (bCalcGeom) {
    calc_geom(ssize,sindex,x, gc, min, max, FALSE);
    printf("new center      :%7.3f%7.3f%7.3f (nm)\n",gc[XX],gc[YY],gc[ZZ]);
  }
  if (bOrient || bScale || bDist || bSetSize) {
    printf("new box vectors :%7.3f%7.3f%7.3f (nm)\n", 
	   norm(box[XX]), norm(box[YY]), norm(box[ZZ]));
    printf("new box angles  :%7.2f%7.2f%7.2f (degrees)\n",
	   RAD2DEG*acos(cos_angle_no_table(box[YY],box[ZZ])),
	   RAD2DEG*acos(cos_angle_no_table(box[XX],box[ZZ])),
	   RAD2DEG*acos(cos_angle_no_table(box[XX],box[YY])));
    printf("new box volume  :%7.2f               (nm^3)\n",det(box));
  }  

  if (bIndex) {
    fprintf(stderr,"\nSelect a group for output:\n");
    get_index(&atoms,opt2fn_null("-n",NFILE,fnm),
	      1,&isize,&index,&grpname);
    if (opt2bSet("-bf",NFILE,fnm))
      fatal_error(0,"combination not implemented: -bf -n  or -bf -ndef");
    else
      write_sto_conf_indexed(outfile,title,&atoms,x,bHaveV?v:NULL,box,
			     isize,index); 
  }
  else {
    if (outftp != efPDB) {
      write_sto_conf(outfile,title,&atoms,x,bHaveV?v:NULL,box); 
    } 
    else {
      out=ffopen(outfile,"w");
      if (opt2bSet("-bf",NFILE,fnm)) {
	read_bfac(opt2fn("-bf",NFILE,fnm),&n_bfac,&bfac,&bfac_nr);
	set_pdb_conf_bfac(atoms.nr,atoms.nres,&atoms,
			  n_bfac,bfac,bfac_nr,peratom);
      }
      if (opt2parg_bSet("-label",NPA,pa)) {
	for(i=0; (i<atoms.nr); i++) 
	  atoms.atom[i].chain=label[0];
      }
      if (bMead) {
	set_pdb_wide_format(TRUE);
        fprintf(out,"REMARK    The b-factors in this file hold atomic radii\n");
	fprintf(out,"REMARK    The occupancy in this file hold atomic charges\n");
      }
      write_pdbfile(out,title,&atoms,x,box,0,!bLegend && visbox[0]<=0);
      if (bLegend)
	pdb_legend(out,atoms.nr,atoms.nres,&atoms,x);
      if (visbox[0] > 0)
	visualize_box(out,bLegend ? atoms.nr+12 : atoms.nr,
		      bLegend? atoms.nres=12 : atoms.nres,box,visbox);
      fclose(out);
    }  
  }
  
  thanx(stdout);
  
  return 0;
}

