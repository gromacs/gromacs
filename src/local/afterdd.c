#include "typedefs.h"
#include "3dview.h"
#include "statutil.h"
#include "smalloc.h"
#include "copyrite.h"
#include "rdgroup.h"
#include "confio.h"
#include "fatal.h"
#include "vec.h"

static void rot_conf(t_atoms *atoms,rvec x[],rvec v[],rvec trans,real angle,
		     rvec head,rvec tail,matrix box,int isize,atom_id index[])
{
  rvec     dx,center;
  vec4     xv;
  t_3dview *view;
  int      i,j,ai;
  
  rvec_sub(head,tail,dx);
  printf("Arrow vector: %10.4f  %10.4f  %10.4f\n",dx[XX],dx[YY],dx[ZZ]);
  if (norm(dx) == 0.0)
    fatal_error(0,"Arrow vector not given");
    
  /* Fill the projection matrix */
  view = init_view(box);
  for(i=0; (i<DIM); i++) {
    view->eye[i]    = dx[i];
    view->origin[i] = (head[i]+tail[i])/2;
  }
  view->eye[WW] = view->origin[WW] = 0;
  view->sc_x    = view->sc_y       = 1;
  
  calculate_view(view);
  
  for(i=0; (i<isize); i++) {
    ai = index[i];
    m4_op(view->Rot,x[ai],xv);
    copy_rvec(xv,x[ai]);
    m4_op(view->Rot,v[ai],xv);
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
  static rvec trans = { 0,0,0 };
  static rvec head  = { 0,0,0 };
  static rvec tail  = { 0,0,0 };
  static real angle = 0;
  t_pargs pa[] = {
    { "-angle", FALSE, etREAL, {&angle},
      "Angle of rotation about rotation vector" },
    { "-trans", FALSE, etRVEC, {trans},
      "Translation vector (see DynDom info file)" },
    { "-head",  FALSE, etRVEC, {head},
      "First atom of the arrow vector" },
    { "-tail",  FALSE, etRVEC, {tail},
      "Last atom of the arrow vector" }
  };
  int     natoms,isize;
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
  
  rot_conf(&atoms,x,v,trans,angle,head,tail,box,isize,index);
  
  write_sto_conf(opt2fn("-o",NFILE,fnm),title,&atoms,x,v,box);
  
  thanx(stdout);
  
  return 0;
}
