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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_genconf_c = "$Id$";

#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "topexcl.h"
#include "x2top.h"

bool is_bond(char *ai,char *aj,real len)
{
  char atp[6] = "HCNOSX";
#define NATP asize(atp)
  real blen[6][6] = { 
    {  0.00,  0.108, 0.105, 0.10, 0.10, 0.10 },
    {  0.108, 0.15,  0.14,  0.14, 0.16, 0.15 },
    {  0.105, 0.14,  0.14,  0.14, 0.16, 0.15 },
    {  0.10,  0.14,  0.14,  0.14, 0.17, 0.16 },
    {  0.10,  0.16,  0.16,  0.17, 0.20, 0.20 },
    {  0.10,  0.15,  0.15,  0.16, 0.20, 0.20 }
  };
  int i,aai=NATP-1,aaj=NATP-1;
  bool bIsBond;
  
  if (len == 0.0)
    bIsBond = FALSE;
  else {
    /* Check atom type: default unknown (X) */  
    for(i=0; (i<NATP-1); i++) {
      if (ai[0] == atp[i])
	aai=i;
      if (aj[0] == atp[i])
	aaj=i;
    }
    /* There is a bond when the deviation from ideal length is less than
     * 10 %
     */
    bIsBond = (len < blen[aai][aaj]*1.1);
  }
  
  if (debug)
    fprintf(debug,"ai: %5s  aj: %5s  len: %8.3f  bond: %s\n",
	    ai,aj,len,BOOL(bIsBond));
	    
  return bIsBond;
}

void mk_bonds(t_atoms *atoms,rvec x[],t_params *bond,int nbond[])
{
  t_param b;
  int     i,j;
  rvec    dx;
  real    len;
  
  for(i=0; (i<MAXATOMLIST); i++)
    b.a[i] = -1;
  for(i=0; (i<MAXFORCEPARAM); i++)
    b.c[i] = 0.0;
  for(i=0; (i<atoms->nr); i++)
    for(j=i+1; (j<atoms->nr); j++) {
      rvec_sub(x[i],x[j],dx);
      len = norm(dx);
		
      if (is_bond(*atoms->atomname[i],*atoms->atomname[j],len)) {
	b.AI = i;
	b.AJ = j;
	b.C0 = len;
	push_bondnow (bond,&b);
	nbond[i]++;
	nbond[j]++;
      }
    }
}

t_atomtype *set_atom_type(t_atoms *atoms,int nbonds[],
			  int nnm,t_nm2type nm2t[])
{
  static t_symtab symtab;
  t_atomtype *atype;
  char elem[2];
  int  i;
  
  open_symtab(&symtab);
  snew(atype,1);
  atype->nr       = atoms->nr;
  atype->atom     = atoms->atom;
  snew(atype->atomname,atoms->nr);
  elem[1]='\0';
  for(i=0; (i<atoms->nr); i++) {
    elem[0] = (*atoms->atomname[i])[0];
    if ((elem[0] >= 0) && (elem[0] <= 9))
      elem[0] = (*atoms->atomname[i])[1];
    atype->atomname[i] = put_symtab(&symtab,
				    nm2type(nnm,nm2t,elem,nbonds[i]));
  } 
  /* MORE CODE */

  close_symtab(&symtab);
    
  return atype;
}

void lo_set_force_const(t_params *plist,real c[],int nrfp,bool bRound)
{
  int    i,j;
  double cc;
  char   buf[32];
  
  for(i=0; (i<plist->nr); i++) {
    if (bRound) {
      sprintf(buf,"%.3e",plist->param[i].c[0]);
      sscanf(buf,"%lf",&cc);
      c[0] = cc;
    }
    else 
      c[0] = plist->param[i].c[0];
    
    for(j=0; (j<nrfp); j++) {
      plist->param[i].c[j]      = c[j];
      plist->param[i].c[nrfp+j] = c[j];
    }
    plist->param[i].s = strdup("");
  }
}

void set_force_const(t_params plist[],real kb,real kt,real kp,bool bRound)
{
  int i;
  real c[MAXFORCEPARAM];
  
  c[0] = 0;
  c[1] = kb;
  lo_set_force_const(&plist[F_BONDS],c,2,bRound);
  c[1] = kt;
  lo_set_force_const(&plist[F_ANGLES],c,2,bRound);
  c[1] = kp;
  c[2] = 3;
  lo_set_force_const(&plist[F_PDIHS],c,3,bRound);
}

void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[])
{
  int    i,ai,aj,ak,al;
  rvec   r_ij,r_kj,r_kl,m,n;
  real   sign,th,costh,ph,cosph;
  matrix box;
  
  clear_mat(box);
  box[XX][XX] = box[YY][YY] = box[ZZ][ZZ] = 1000;
  for(i=0; (i<ang->nr); i++) {
    ai = ang->param[i].AI;
    aj = ang->param[i].AJ;
    ak = ang->param[i].AK;
    th = RAD2DEG*bond_angle(box,x[ai],x[aj],x[ak],r_ij,r_kj,&costh);
    ang->param[i].C0 = th;
  }
  for(i=0; (i<dih->nr); i++) {
    ai = dih->param[i].AI;
    aj = dih->param[i].AJ;
    ak = dih->param[i].AK;
    al = dih->param[i].AL;
    ph = RAD2DEG*dih_angle(box,x[ai],x[aj],x[ak],x[al],
			   r_ij,r_kj,r_kl,m,n,&cosph,&sign);
    dih->param[i].C0 = ph;
  }
}

void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[])
{
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    fprintf(fp,"Atom %5s has %1d bonds\n",*atoms->atomname[i],nbonds[i]);
  }
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "x2top generates a primitive topology from a coordinate file.",
    "The program assumes all hydrogens are present when defining",
    "the hybridization from the atom name and the number of bonds."
  };
  FILE       *fp;
  t_params   plist[F_NRE];
  t_atoms    *atoms;       /* list with all atoms */
  t_atomtype *atype;
  t_block    excl;
  t_nextnb   nnb;
  t_nm2type  *nm2t;
  int        nnm;
  char       title[STRLEN];
  rvec       *x;        /* coordinates? */
  int        *nbonds,*cgnr;
  int        bts[] = { 1,1,1,2 };
  matrix     box;          /* box length matrix */
  int        natoms;       /* number of atoms in one molecule  */
  int        nres;         /* number of molecules? */
  int        i,j,k,l,m;
  
  t_filenm fnm[] = {
    { efSTX, "-f", "conf", ffREAD  },
    { efTOP, "-o", "out",  ffWRITE }
  };
#define NFILE asize(fnm)
  static real kb    = 4e5,kt = 400,kp = 5;
  static int  nexcl = 3;
  static bool bH14  = FALSE,bAllDih = FALSE, bRound = TRUE;
  t_pargs pa[] = {
    { "-kb",    FALSE, etREAL, &kb,
      "Bonded force constant (kJ/mol/nm^2)" },
    { "-kt",    FALSE, etREAL, &kt,
      "Angle force constant (kJ/mol/rad^2)" },
    { "-kp",    FALSE, etREAL, &kp,
      "Dihedral angle force constant (kJ/mol/rad^2)" },
    { "-nexcl", FALSE, etINT,  &nexcl,
      "Number of exclusions" },
    { "-H14",    FALSE, etBOOL, &bH14, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-alldih", FALSE, etBOOL, {&bAllDih}, 
      "Generate all proper dihedrals" },
    { "-round",  FALSE, etBOOL, &bRound,
      "Round off measured values" }
  };
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  /* Init lookup table for invsqrt */	    
  init_lookup_table(stdout);
  
  /* Init parameter lists */
  init_plist(plist);
  
  /* Read coordinates */
  get_stx_coordnum(opt2fn("-f",NFILE,fnm),&natoms); 
  snew(atoms,1);
  
  /* make space for all the atoms */
  init_t_atoms(atoms,natoms,FALSE);
  snew(x,natoms);              

  read_stx_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,NULL,box);
 
  
  snew(nbonds,atoms->nr);
  
  printf("Generating bonds from distances...\n");
  mk_bonds(atoms,x,&(plist[F_BONDS]),nbonds);

  nm2t = rd_nm2type(&nnm);
  printf("There are %d name to type translations\n",nnm);
  if (debug)
    dump_nm2type(debug,nnm,nm2t);
  atype = set_atom_type(atoms,nbonds,nnm,nm2t);
  
  /* Make Angles and Dihedrals */
  printf("Generating angles and dihedrals from bonds...\n");
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,plist,NULL,bAllDih);
  done_nnb(&nnb);

  fprintf(stderr,
	  "There are %4d dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d atoms\n",
	  plist[F_PDIHS].nr, plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	  plist[F_LJ14].nr, plist[F_BONDS].nr,atoms->nr);

  init_block(&excl);
  snew(cgnr,atoms->nr);

  calc_angles_dihs(&plist[F_ANGLES],&plist[F_PDIHS],x);
  
  set_force_const(plist,kb,kt,kp,bRound);
  
  fp = ftp2FILE(efTOP,NFILE,fnm,"w");
  write_top(fp,NULL,"x2top",atoms,bts,plist,&excl,atype,cgnr,nexcl);

  fclose(fp);
  
  if (debug) {
    dump_hybridization(debug,atoms,nbonds);
  }
  
  thanx(stdout);
  
  return 0;
}
