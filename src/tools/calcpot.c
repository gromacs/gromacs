#include "nrnb.h"
#include "vec.h"
#include "calcpot.h"
#include "mdebin.h"
#include "mshift.h"
#include "smalloc.h"
#include "force.h"
#include "main.h"
#include "filenm.h"
#include "mdrun.h"
#include "ns.h"

static void c_tabpot(int inr,real ix,real iy,real iz,real qi,
		     real pos[],int nj,int type[],int jjnr[],
		     real charge[],real pot[],
		     int ntab,real tabscale,real VFtab[])
{
  int       k,jnr,j3;
  real      rijX,rijY,rijZ;
  real      vijcoul,fijC,fijscal;
  real      fjx,fjy,fjz;
  real      tx,ty,tz;
  real      vctot,vnbtot;
  real      qq,rsq;
  real      r1,r1t,h_1,poti;
  real      eps,eps2,Y,F,Fp,Geps,Heps2,two=2.0,VV,FF;
  int       n0,n1,nnn;
  
  poti   = 0;
  vctot  = 0;
  vnbtot = 0;
  h_1    = tabscale;
  
  /* See comment in c_coultab (inloopc.c) */
  for(k=0; (k<nj); k++) {  
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3];
    rijY           = iy - pos[j3+1];
    rijZ           = iz - pos[j3+2];
         
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    r1             = sqrt(rsq);
    r1t            = r1*tabscale;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
        
#define EXTRACT(nn) { nnn = nn; Y=VFtab[nnn]; F=VFtab[nnn+1]; Geps=VFtab[nnn+2]*eps; Heps2=VFtab[nnn+3]*eps2; Fp=F+Geps+Heps2; VV=Y+eps*Fp; FF=Fp+Geps+two*Heps2; }

    /* Coulomb */
    EXTRACT(n1)
    pot[jnr]      += qi*VV;
    poti          += charge[jnr]*VV;
  }
  
  pot[inr] += poti;
  
}

static void low_calc_pot(FILE *log,int ftype,t_forcerec *fr,
			 rvec x[],t_mdatoms *mdatoms,rvec box_size,real pot[])
{
  int      i,itpA,itpB,gid,m,nj,inr,iinr,nri,k;
  rvec     r_i,f_ip,fw[3],xw[3];
  real     qi,Vnb,Vc,eps;
  t_nblist *nlist;
  t_nl_i   *nl_i;
  t_nl_j   *nl_j;
  int      *typeA,*typeB;
  real     *chargeA,*chargeB;
  rvec     *svec,*fshift;
  int      nr_inter;
  bool     bWater,bTab;
  
  typeA   = mdatoms->typeA;
  chargeA = mdatoms->chargeA;

  svec    = fr->shift_vec;
  fshift  = fr->fshift;
  eps     = fr->epsfac;
  
  nr_inter = 0;

  if (ftype == F_SR) 
    nlist=fr->coul;
  else
    nlist=fr->vdw;
  
  for(gid=0; (gid<fr->nn); gid++) {
    nri  = nlist[gid].nri;
    nl_i = nlist[gid].nl_i;
    Vnb  = 0;
    Vc   = 0;
    
    for (i=0; (i<nri); i++) {
      inr      = nl_i[i].i_atom;
      k        = nl_i[i].shift;
      nj       = nl_i[i].nj;
      nl_j     = &(nlist[gid].nl_j[nl_i[i].j_index]);
      itpA     = typeA[inr];
      if (nl_i[i].bWater)
	fatal_error(0,"Water interaction found...\n");
      nr_inter += nj;
      
      qi = chargeA[inr]*eps;
      rvec_add(x[inr],svec[k],r_i);
      clear_rvec(f_ip);
      
      c_tabpot(inr,r_i[XX],r_i[YY],r_i[ZZ],qi,x[0],nj,typeA,nl_j,
	       chargeA,pot,fr->ntab,fr->tabscale,fr->VFtab);
    }
  }
  fprintf(stderr,"There were %d interactions\n",nr_inter);
}

void calc_pot(FILE *logf,t_nsborder *nsb,t_commrec *cr,t_groups *grps,
	      t_parm *parm,t_topology *top,rvec x[],t_forcerec *fr,
	      t_graph *graph,t_mdatoms *mdatoms,real pot[])
{
  static bool        bFirst=TRUE;
  static t_nrnb      nrnb;
  static rvec        *f,*buf;
  tensor      force_vir,shake_vir;
  real        ener[F_NRE],qsmall;
  rvec        vcm,box_size;
  int         i,j,m,atnr2,fp_ene;

  /* Some dummies */
  int         step,natoms;
  real        t;

  /* Calc the force */
  fprintf(stderr,"Doing single force calculation...\n");

  if (bFirst) {
    snew(f,   nsb->natoms);
    snew(buf, nsb->natoms);
    
    bFirst = FALSE;
  }
  /* Reset long range forces if necessary */
  if (fr->bTwinRange) {
    clear_rvecs(nsb->natoms,fr->flr);
    clear_rvecs(SHIFTS,fr->fshift_lr);
  }
  if (parm->ir.epc != epcNO)
      calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  put_charge_groups_in_box(stdlog,0,top->blocks[ebCGS].nr,FALSE,
			   parm->box,box_size,&(top->blocks[ebCGS]),x,
			   fr->shift_vec,fr->cg_cm);
  mk_mshift(stdlog,graph,parm->box,x);
  /* Do the actual neighbour searching and if twin range electrostatics
   * also do the calculation of long range forces and energies.
   */
  
  ns(logf,fr,x,f,parm->box,grps,&(parm->ir.opts),top,mdatoms,cr,
     &nrnb,nsb,0);
  for(m=0; (m<DIM); m++)
    box_size[m] = parm->box[m][m];
  for(i=0; (i<mdatoms->nr); i++)
    pot[i] = 0;
  low_calc_pot(logf,F_SR,fr,x,mdatoms,box_size,pot); 
  low_calc_pot(logf,F_LJ,fr,x,mdatoms,box_size,pot); 
}

void init_calcpot(int nfile,t_filenm fnm[],t_topology *top,
		  rvec **x,t_parm *parm,t_commrec *cr,
		  t_graph **graph,t_mdatoms **mdatoms,
		  t_nsborder *nsb,t_groups *grps,
		  t_forcerec **fr,real **pot,
		  matrix box)
{
  real     t,t0,lam,lam0,SAfac;
  bool     bTYZ;
  char     *traj,*xtc_traj;
  rvec     *v;
  t_nrnb   nrnb;
  t_mdebin *mdebin;
  int      fp_ene,m;
  rvec     vcm,box_size;
  tensor   force_vir,shake_vir;
  
  /* Initiate */
  cr->nprocs = 1; cr->pid    = 0; cr->left   = 0; cr->right  = 1;
  open_log(ftp2fn(efLOG,nfile,fnm),cr);
#ifdef CINVSQRT
  init_lookup_table(stdlog);
#endif

  init_nrnb(&nrnb);
  init_single(stdlog,parm,ftp2fn(efTPX,nfile,fnm),top,x,&v,mdatoms,nsb);
  init_md(cr,&(parm->ir),&t,&t0,&lam,&lam0,&SAfac,
	  &nrnb,&bTYZ,top,nfile,fnm,&traj,&xtc_traj,&fp_ene,
	  &mdebin,grps,vcm,force_vir,shake_vir,*mdatoms);
  init_groups(stdlog,*mdatoms,&(parm->ir.opts),grps);  

  /* Calculate intramolecular shift vectors to make molecules whole again */
  *graph = mk_graph(&(top->idef),top->atoms.nr,FALSE);
  mk_mshift(stdlog,*graph,parm->box,*x);
  
  /* Turn off watertype optimizations, to ease coding above. */
  parm->ir.solvent_opt = -1;
  
  /* Turn off free energy computation */
  parm->ir.bPert = FALSE;
  
  /* Initiate forcerecord */
  *fr = mk_forcerec();
  init_forcerec(stdlog,*fr,&(parm->ir),&(top->blocks[ebMOLS]),cr,
		&(top->blocks[ebCGS]),&(top->idef),*mdatoms,parm->box,FALSE);

  /* Remove periodicity */  
  for(m=0; (m<DIM); m++)
    box_size[m] = parm->box[m][m];
  if (parm->ir.eBox != ebtNONE)
    do_pbc_first(stdlog,parm,box_size,*fr,*graph,*x);

  copy_mat(parm->box,box);
      
  snew(*pot,nsb->natoms);
  
  sfree(v);
}
