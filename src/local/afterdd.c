#include "typedefs.h"
#include "3dview.h"
#include "statutil.h"
#include "smalloc.h"
#include "copyrite.h"
#include "rdgroup.h"
#include "confio.h"
#include "fatal.h"
#include "vec.h"
#include "physics.h"
#include "random.h"

static void rot_conf(t_atoms *atoms,rvec x[],rvec v[],real trans,real angle,
		     rvec head,rvec tail,matrix box,int isize,atom_id index[])
{
  rvec     arrow,center,xcm;
  real     theta,phi,arrow_len;
  mat4     Rx,Ry,Rz,Rinvy,Rinvz,Mtot,Tcm,Tinvcm,Tx;
  mat4     temp1,temp2,temp3,temp4,temp21,temp43;
  vec4     xv;
  int      i,j,ai;
  
  rvec_sub(tail,head,arrow);
  arrow_len = norm(arrow);
  printf("Arrow vector:   %10.4f  %10.4f  %10.4f\n",
	 arrow[XX],arrow[YY],arrow[ZZ]);
  printf("Effective translation %g nm\n",trans);
  if (arrow_len == 0.0)
    fatal_error(0,"Arrow vector not given");

  /* Compute center of mass and move atoms there */
  clear_rvec(xcm);
  for(i=0; (i<isize); i++)
    rvec_inc(xcm,x[index[i]]);
  for(i=0; (i<DIM); i++)
    xcm[i] /= isize;
  printf("Center of mass: %10.4f  %10.4f  %10.4f\n",xcm[XX],xcm[YY],xcm[ZZ]);
  for(i=0; (i<isize); i++)
    rvec_dec(x[index[i]],xcm);
  
  /* Compute theta and phi that describe the arrow */
  theta = acos(arrow[ZZ]/arrow_len);
  phi   = atan2(arrow[YY]/arrow_len,arrow[XX]/arrow_len);
  printf("Phi = %.1f, Theta = %.1f\n",RAD2DEG*phi,RAD2DEG*theta);

  /* Now the total rotation matrix: */
  /* Rotate a couple of times */
  rotate(ZZ,-phi,Rz);
  rotate(YY,M_PI/2-theta,Ry);
  rotate(XX,angle*DEG2RAD,Rx);
  Rx[WW][XX] = trans;
  rotate(YY,theta-M_PI/2,Rinvy);
  rotate(ZZ,phi,Rinvz);
  
  mult_matrix(temp1,Ry,Rz);
  mult_matrix(temp2,Rinvy,Rx);
  mult_matrix(temp3,temp2,temp1);
  mult_matrix(Mtot,Rinvz,temp3);

  print_m4(debug,"Rz",Rz);
  print_m4(debug,"Ry",Ry);
  print_m4(debug,"Rx",Rx);
  print_m4(debug,"Rinvy",Rinvy);
  print_m4(debug,"Rinvz",Rinvz);
  print_m4(debug,"Mtot",Mtot);

  for(i=0; (i<isize); i++) {
    ai = index[i];
    m4_op(Mtot,x[ai],xv);
    rvec_add(xv,xcm,x[ai]);
    m4_op(Mtot,v[ai],xv);
    copy_rvec(xv,v[ai]);
  }
}

int main(int argc,char *argv[])
{
  char *desc[] = {
    "afterdd reads a pdb file output from DynDom",
    "http://md.chem.rug.nl/~steve/DynDom/dyndom.home.html",
    "It reads the coordinates, and the coordinates of the rotation axis",
    "furthermore it reads an index file containing the domains.",
    "Furthermore it takes the first and last atom of the arrow file",
    "as command line arguments (head and tail) and",
    "finally it takes the translation vector (given in DynDom info file)",
    "and the angle of rotation (also as command line arguments). If the angle",
    "determined by DynDom is given, one should be able to recover the",
    "second structure used for generating the DynDom output.",
    "Because of limited numerical accuracy this should be verified by",
    "computing an all-atom RMSD (using [TT]g_confrms[tt]) rather than by file",
    "comparison (using diff).[PAR]",
    "The purpose of this program is to interpolate and extrapolate the",
    "rotation as found by DynDom. As a result unphysical structures with",
    "long or short bonds, or overlapping atoms may be produced. Visual",
    "inspection, and energy minimization may be necessary to",
    "validate the structure."
  };
  static real trans = 0;
  static rvec head  = { 0,0,0 };
  static rvec tail  = { 0,0,0 };
  static real angle = 0,maxangle = 0;
  static int  label = 0;
  t_pargs pa[] = {
    { "-angle",    FALSE, etREAL, {&angle},
      "Angle of rotation about rotation vector" },
    { "-maxangle", FALSE, etREAL, {&maxangle},
      "DymDom dtermined angle of rotation about rotation vector" },
    { "-trans",    FALSE, etREAL, {&trans},
      "Translation along rotation vector (see DynDom info file)" },
    { "-head",     FALSE, etRVEC, {head},
      "First atom of the arrow vector" },
    { "-tail",     FALSE, etRVEC, {tail},
      "Last atom of the arrow vector" },
    { "-label",    FALSE, etINT,  {&label},
      "Label in the outgoing pdb file is made by adding this value to 'A'" }
  };
  int     i,natoms,isize;
  atom_id *index=NULL;
  char    title[256],*grpname;
  t_atoms atoms;
  rvec    *x,*v;
  matrix  box;
  
  t_filenm fnm[] = {
    { efPDB, "-f", "dyndom",  ffREAD },
    { efPDB, "-o", "rotated", ffWRITE },
    { efNDX, "-n", "domains", ffREAD }
  };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  get_stx_coordnum (opt2fn("-f",NFILE,fnm),&natoms);
  init_t_atoms(&atoms,natoms,TRUE);
  snew(x,natoms);
  snew(v,natoms);
  read_stx_conf(opt2fn("-f",NFILE,fnm),title,&atoms,x,v,box);
  
  printf("Select group to rotate:\n"); 
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
  printf("Going to rotate %s containg %d atoms\n",grpname,isize);

  trans = -trans*0.1*angle/maxangle;
  rot_conf(&atoms,x,v,trans,angle,head,tail,box,isize,index);

  if (label < 0)
    label = 'A';
  else {
    if (label+'A' >= 'Z')
      label = label % 26;
    label += 'A';
  }
  for(i=0; (i<atoms.nr); i++)
    atoms.atom[i].chain = label;
    
  write_sto_conf(opt2fn("-o",NFILE,fnm),title,&atoms,x,v,box);
  
  thanx(stdout);
  
  return 0;
}
