/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <math.h>
#include "smalloc.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "confio.h"
#include "tpxio.h"
#include "index.h"
#include "gstat.h"
#include "matio.h"
#include "do_fit.h"
#include "princ.h"
#include "gmx_ana.h"

typedef struct t_monomer {
    t_atoms *atoms; // atoms of this monomer
    rvec    *x;     // coordinates
    atom_id *index; // index
    int     isize; // index size
    char    **grp;  // groupname
    rvec    cm; // center of mass
    real    mass; // mass of monomer
} t_monomer;

void done_monomer(t_monomer *monomer) {
    sfree(monomer->index);
    sfree(monomer->grp);
    sfree(monomer->x);
    done_atom(monomer->atoms);
    sfree(monomer);
}

t_atoms *copy_t_atoms_indexed(t_atoms *src, atom_id *index, int isize)
{
  t_atoms *dst;
  int i;

  snew(dst,1);
  /*
   * we will only copy only one index group to new atom structure
   * thus dst->nr = isize
   */
  init_t_atoms(dst,isize,(NULL != src->pdbinfo));
  dst->nr = isize;
  if (NULL != src->atomname)
      snew(dst->atomname,isize);
  if (NULL != src->atomtype)
      snew(dst->atomtype,isize);
  if (NULL != src->atomtypeB)
      snew(dst->atomtypeB,isize);
  for(i=0; (i<isize); i++) {
      dst->atom[i] = src->atom[index[i]];
    if (NULL != src->pdbinfo)
        dst->pdbinfo[i] = src->pdbinfo[index[i]];
    if (NULL != src->atomname)
        dst->atomname[i]  = src->atomname[index[i]];
    if (NULL != src->atomtype)
        dst->atomtype[i] = src->atomtype[index[i]];
    if (NULL != src->atomtypeB)
        dst->atomtypeB[i] = src->atomtypeB[index[i]];
  }
  dst->nres = src->nres;
  for(i=0; (i<isize); i++) {
    dst->resinfo[i] = src->resinfo[index[i]];
  }
  return dst;
}

int gmx_hpgc(int argc,char *argv[])
{
    const char *desc[] = {
        "This tool compute protein filament properties as a quantity of time such as:",
        "[PAR]",
        "helical pitch",
        "[PAR]",
        "number of monomers per turn",
        "[PAR]",
	"helical rotation",
	"[PAR]",
        "As input you need at lest one filament period and any two neighbour monomers"
    };
    static gmx_bool bPBC=TRUE;

#define NPA asize(pa)

    t_pargs pa[] = {
        { "-pbc", FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances" },
    };
    FILE      *fp;
    const char *fnTPS,*fnTRX,*fnNDX;
    t_trxstatus *status;
    t_topology *top=NULL;
    gmx_rmpbc_t  gpbc=NULL;
    t_pbc     pbc;
    gmx_bool  bTPX;
    int        ePBC=-1;
    matrix     box;
    matrix     R;
    char       title[STRLEN];
    rvec       *x;
    rvec       e;
    rvec       dx;
    int       natoms;
    real       t;
    real       *time;
    real       *w_rls;
    real       *theta;
    real       dh;
    real       *npturn;
    real       *pitch;
    int         i,j;
    int        nrframes;
    t_monomer  *mono1=NULL, *mono2=NULL;
    output_env_t oenv;
    static int NDIM=3;

#define NFILE asize(fnm)

  t_filenm   fnm[] = {
      { efTPS,  "-s",         NULL,       ffREAD },
      { efTRX,  "-f",         NULL,       ffREAD },
      { efNDX,  NULL,         NULL,       ffOPTRD },
      { efXVG,  "-pitch",     "pitch",    ffWRITE},
      { efXVG,  "-rot",       "rotation", ffWRITE},
      { efXVG,  "-nmon",      "nmon",     ffWRITE},
  };

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

  /* Try to read files */
  fnTPS = ftp2fn(efTPS,NFILE,fnm);
  fnTRX = ftp2fn(efTRX,NFILE,fnm);

  snew(top,1);
  snew(mono1,1);
  snew(mono2,1);
  snew(mono1->index,1);
  snew(mono2->index,1);
  snew(mono1->grp,1);
  snew(mono2->grp,1);

  bTPX=read_tps_conf(fnTPS,title,top,&ePBC,&x,NULL,box,TRUE);

  printf("\nPlease select two monomers of protein filament:\n");
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),1,&(mono1->isize),&(mono1->index),mono1->grp);
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),1,&(mono2->isize),&(mono2->index),mono2->grp);

  /* Prepare reference frame */
  if (bPBC) {
      gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
      set_pbc(&pbc,ePBC,box);
      gmx_rmpbc(gpbc,top->atoms.nr,box,x);
  }

  natoms=read_first_x(oenv,&status,fnTRX,&t,&x,box);

  if (natoms != top->atoms.nr) {
      fprintf(stderr,"\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n",natoms,top->atoms.nr);
  }

  /*
   * copy atoms records
   */
  mono1->atoms = copy_t_atoms_indexed(&(top->atoms),mono1->index,mono1->isize);
  mono2->atoms = copy_t_atoms_indexed(&(top->atoms),mono2->index,mono2->isize);
  /*
   * Init some memory structures
   */
  snew(time,1);
  snew(theta,1);
  snew(pitch,1);
  snew(npturn,1);
  snew(mono1->x,mono1->isize);
  snew(mono2->x,mono2->isize);

  nrframes = 0;
  /*
   * Init weights for CoM
   */
  snew(w_rls,mono1->isize);
  for(i=0;i<mono1->isize;i++) {
      w_rls[i] = 1.0;
  }

  do {
      /* Fix PBC */
      if (bPBC) {
          set_pbc(&pbc,ePBC,box);
          gmx_rmpbc(gpbc,top->atoms.nr,box,x);
      }
      nrframes++;
      srenew(time,nrframes+1);
      srenew(theta,nrframes+1);
      srenew(pitch,nrframes+1);
      srenew(npturn,nrframes+1);
      time[nrframes-1] = t;
      /*
       * Copy monomers to monomers_t structures
       */
      for(i=0;i<mono1->isize;i++) {
          copy_rvec(x[mono1->index[i]],mono1->x[i]);
      }
      for(i=0;i<mono2->isize;i++) {
          copy_rvec(x[mono2->index[i]],mono2->x[i]);
      }
      /*
       * calculate Center of mass for monomers
       */
      mono1->mass = calc_xcm(mono1->x,mono1->isize,NULL,mono1->atoms->atom,mono1->cm,FALSE);
      mono2->mass = calc_xcm(mono2->x,mono2->isize,NULL,mono2->atoms->atom,mono2->cm,FALSE);
      rvec_sub(mono2->cm,mono1->cm,dx);
      /*
       * reset CoM to coordinate origin
       */
      reset_x(mono1->isize,NULL,mono1->isize,NULL,mono1->x,w_rls);
      reset_x(mono2->isize,NULL,mono2->isize,NULL,mono2->x,w_rls);

      /*
       * need to reset cm to origin before calculating rotation matrix
       */
      calc_fit_R(NDIM,mono1->isize,w_rls,mono1->x,mono2->x,R);
      /*
       * now get theta rotation angle
       */
      theta[nrframes-1] = acos(0.5*(R[XX][XX] + R[YY][YY] + R[ZZ][ZZ] - 1.0));
      /* now we can get rotation axis if angle wasnt nPi */
      e[XX] = (R[ZZ][YY] - R[YY][ZZ])/(2.0*sin(theta[nrframes-1]));
      e[YY] = (R[XX][ZZ] - R[ZZ][XX])/(2.0*sin(theta[nrframes-1]));
      e[ZZ] = (R[YY][XX] - R[XX][YY])/(2.0*sin(theta[nrframes-1]));
      /*
       * Now we want to compute h per monomer
       */
      dh = iprod(dx,e);
      npturn[nrframes-1] = M_2PI/theta[nrframes-1];
      pitch[nrframes-1] = fabs(dh)*npturn[nrframes-1];
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  close_trj(status);
  /*
   * Clean up unneded structures
   */
  done_monomer(mono1);
  done_monomer(mono2);
  done_top(top);
  /*
   * Print data
   */
  fp = xvgropen(opt2fn_null("-pitch",NFILE,fnm),"Helical filament pitch","Time (ps)","Pitch (nm)",oenv);
  for(i=0;i<nrframes;i++) {
      fprintf(fp,"%10.6f%10.6f\n",time[i],pitch[i]);
  }
  xvgrclose(fp);

  fp = xvgropen(opt2fn_null("-rot",NFILE,fnm),"Rotation per monomer","Time (ps)","Rotation (degrees)",oenv);
  for(i=0;i<nrframes;i++) {
      fprintf(fp,"%10.6f%10.6f\n",time[i],theta[i]*RAD2DEG);
  }
  xvgrclose(fp);

  fp = xvgropen(opt2fn_null("-nmon",NFILE,fnm),"Number of monomers per helix turn","Time (ps)","Number of monomers",oenv);
  for(i=0;i<nrframes;i++) {
      fprintf(fp,"%10.6f%10.6f\n",time[i],npturn[i]);
  }
  xvgrclose(fp);

  /*
   * Clean up remaining structures
   */
  sfree(time);
  sfree(pitch);
  sfree(theta);
  sfree(npturn);

  thanx(stderr);

  return 0;
}
