#include <stdio.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "macros.h"
#include "txtdump.h"
#include "tpxio.h"
#include "enxio.h"
#include "names.h"
#include "statutil.h"
#include "copyrite.h"

real ener(matrix P,real e,real e0,int nmol)
{
  real kp = 1;
  real ke = 100;
  
  return (kp*(sqr(P[XX][XX]-1)+sqr(P[YY][YY]-1)+sqr(P[ZZ][ZZ]-1)+
	      sqr(P[XX][YY])+sqr(P[XX][ZZ])+sqr(P[YY][ZZ])) +
	  ke*(e-e0)/nmol);
}

void do_sim(char *enx,
	    t_topology *top,rvec *x,rvec *v,t_inputrec *ir,matrix box)
{
  char *tpx = "optwat.tpx";
  char buf[128];
  
  write_tpx(tpx,0,0.0,0.0,ir,box,top->atoms.nr,x,v,NULL,top);
  sprintf(buf,"xmdrun -s %s -e %s",tpx,enx);
  system(buf);
}

void get_results(char *enx,real P[],real *e,int pstart)
{
  int      fp_ene;
  char     **nms=NULL;
  int      nre,step,ndr,i;
  real     t;
  t_energy *ener;
  
  fp_ene = open_enx(enx,"r");
  
  do_enxnms(fp_ene,&nre,&nms);
  snew(ener,nre);
  
  /* Read until the last frame */
  while (do_enx(fp_ene,&t,&step,&nre,ener,&ndr,NULL));
  
  close_enx(fp_ene);
  
  *e = ener[F_EPOT].e;
  for(i=pstart; (i<pstart+9); i++)
    P[i-pstart] = ener[i].e;
    
  sfree(ener);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "optwat optimizes the force field parameter set of a molecular crystal",
    "to reproduce the pressure tensor and experimental energy.[PAR]",
    "Note that for good results the tpx file must contain input for a",
    "simulated annealing run, or a single point energy calculation at 0 K"
  };
  t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efENX, NULL,      NULL,       ffRW },
    { efLOG, "-g",      NULL,       ffWRITE }
  };
#define NFILE asize(fnm)

  static real e0    = 57,tol = 1;
  static int maxnit = 100,pstart=23;
  static t_pargs pa[] = {
    { "-e0",  FALSE, etREAL, {&e0},
      "Potential energy in kJ/mol" },
    { "-tol", FALSE, etREAL, {&tol},
      "Tolerance for converging" },
    { "-nit", FALSE, etINT, {&maxnit},
      "Max number of iterations" },
    { "-pstart", FALSE, etINT, {&pstart},
      "Index of P[X][X] in the energy file (check with g_energy and subtract 1)" }
  };

  FILE        *fp;  
  t_topology  top;
  t_tpxheader sh;
  t_inputrec  ir;
  rvec        *xx,*vv;
  matrix      box;
  int         i,step,natoms,nmol,nit;
  real        t,lambda,epot,eFF;
  bool        bConverged;
  tensor      P;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  /* Read initial topology and coordaintes etc. */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&sh);
  snew(xx,sh.natoms);
  snew(vv,sh.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,box,&natoms,
	   xx,vv,NULL,&top);

	   
  fp   = ftp2FILE(efLOG,NFILE,fnm,"w");
  nmol = top.blocks[ebMOLS].nr;
	   
  /* Loop over iterations */
  nit = 0;
  do {
    do_sim(ftp2fn(efENX,NFILE,fnm),&top,xx,vv,&ir,box);
    get_results(ftp2fn(efENX,NFILE,fnm),P[0],&epot,pstart);

    eFF = ener(P,epot,e0,nmol);
    bConverged = (eFF < tol);
    
    fprintf(stderr,"Iter %3d eFF = %g, Converged = %s\n",nit,eFF,
	    yesno_names[bConverged]);
    fprintf(fp,"Iter %3d eFF = %g, Converged = %s\n",nit,eFF,
	    yesno_names[bConverged]);
    fprintf(fp,"Epot = %g\n",epot);
    pr_rvecs(fp,0,"Pres",P,DIM);
    fprintf(fp,"\n");
	        
    nit++;
  } while (!bConverged && (nit < maxnit));

  fclose(fp);
    
  thanx(stdout);
  
  return 0;
}
