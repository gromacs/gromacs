/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "mshift.h"
#include "xvgr.h"
#include "rmpbc.h"
#include "tpxio.h"
#include "do_fit.h"
#include "gmx_ana.h"


static void get_refx(output_env_t oenv,const char *trxfn,int nfitdim,int skip,
                     int gnx,int *index,
                     gmx_bool bMW,t_topology *top,int ePBC,rvec *x_ref)
{
    int    natoms,nfr_all,nfr,i,j,a,r,c,min_fr;
    t_trxstatus *status;
    real   *ti,min_t;
    double tot_mass,msd,*srmsd,min_srmsd,srmsd_tot;
    rvec   *x,**xi;
    real   xf;
    matrix box,R;
    real   *w_rls;
    gmx_rmpbc_t  gpbc=NULL;


    nfr_all = 0;
    nfr     = 0;
    snew(ti,100);
    snew(xi,100);
    natoms = read_first_x(oenv,&status,trxfn,&ti[nfr],&x,box); 

    snew(w_rls,gnx);
    tot_mass = 0;
    for(a=0; a<gnx; a++)
    {
        if (index[a] >= natoms)
        {
            gmx_fatal(FARGS,"Atom index (%d) is larger than the number of atoms in the trajecory (%d)",index[a]+1,natoms);
        }
        w_rls[a] = (bMW ? top->atoms.atom[index[a]].m : 1.0);
        tot_mass += w_rls[a];
    }
    gpbc = gmx_rmpbc_init(&top->idef,ePBC,natoms,box);

    do
    {
        if (nfr_all % skip == 0)
        {
            gmx_rmpbc(gpbc,natoms,box,x);
            snew(xi[nfr],gnx);
            for(i=0; i<gnx; i++)
            {
                copy_rvec(x[index[i]],xi[nfr][i]);
            }
            reset_x(gnx,NULL,gnx,NULL,xi[nfr],w_rls);
            nfr++;
            if (nfr % 100 == 0)
            {
                srenew(ti,nfr+100);
                srenew(xi,nfr+100);
            }
        }
        nfr_all++;
    }
    while(read_next_x(oenv,status,&ti[nfr],natoms,x,box));
    close_trj(status);
    sfree(x);

    gmx_rmpbc_done(gpbc);

    snew(srmsd,nfr);
    for(i=0; i<nfr; i++)
    {
        printf("\rProcessing frame %d of %d",i,nfr);
        for(j=i+1; j<nfr; j++)
        {
            calc_fit_R(nfitdim,gnx,w_rls,xi[i],xi[j],R);

            msd = 0;
            for(a=0; a<gnx; a++)
            {
                for(r=0; r<DIM; r++)
                {
                    xf = 0;
                    for(c=0; c<DIM; c++)
                    {
                        xf += R[r][c]*xi[j][a][c];
                    }
                    msd += w_rls[a]*sqr(xi[i][a][r] - xf);
                }
            }
            msd /= tot_mass;
            srmsd[i] += sqrt(msd);
            srmsd[j] += sqrt(msd);
        }
        sfree(xi[i]);
    }
    printf("\n");
    sfree(w_rls);

    min_srmsd = GMX_REAL_MAX;
    min_fr    = -1;
    min_t     = -1;
    srmsd_tot = 0;
    for(i=0; i<nfr; i++)
    {
        srmsd[i] /= (nfr - 1);
        if (srmsd[i] < min_srmsd)
        {
            min_srmsd = srmsd[i];
            min_fr    = i;
            min_t     = ti[i];
        }
        srmsd_tot += srmsd[i];
    }
    sfree(srmsd);

    printf("Average RMSD between all structures: %.3f\n",srmsd_tot/nfr);
    printf("Structure with lowest RMSD to all others: time %g, av. RMSD %.3f\n",
           min_t,min_srmsd);

    for(a=0; a<gnx; a++)
    {
        copy_rvec(xi[min_fr][a],x_ref[index[a]]);
    }

    sfree(xi);
}

int gmx_rotmat(int argc,char *argv[])
{
    const char *desc[] = {
        "g_rotmat plots the rotation matrix required for least squares fitting",
        "a conformation onto the reference conformation provided with",
        "[TT]-s[tt]. Translation is removed before fitting.",
        "The output are the three vectors that give the new directions",
        "of the x, y and z directions of the reference conformation,",
        "for example: (zx,zy,zz) is the orientation of the reference",
        "z-axis in the trajectory frame.",
        "[PAR]",
        "This tool is useful for, for instance,",
        "determining the orientation of a molecule",
        "at an interface, possibly on a trajectory produced with",
        "[TT]trjconv -fit rotxy+transxy[tt] to remove the rotation",
        "in the xy-plane.",
        "[PAR]",
        "Option [TT]-ref[tt] determines a reference structure for fitting,",
        "instead of using the structure from [TT]-s[tt]. The structure with",
        "the lowest sum of RMSD's to all other structures is used.",
        "Since the computational cost of this procedure grows with",
        "the square of the number of frames, the [TT]-skip[tt] option",
        "can be useful. A full fit or only a fit in the x/y plane can",
        "be performed.",
        "[PAR]",
        "Option [TT]-fitxy[tt] fits in the x/y plane before determining",
        "the rotation matrix."
    };
    const char *reffit[] = 
        { NULL, "none", "xyz", "xy", NULL }; 
    static int  skip=1;
    static gmx_bool bFitXY=FALSE,bMW=TRUE;
    t_pargs pa[] = {
        { "-ref", FALSE, etENUM, {reffit},
          "Determine the optimal reference structure" },
        { "-skip", FALSE, etINT, {&skip},
          "Use every nr-th frame for -ref" },
        { "-fitxy", FALSE, etBOOL, {&bFitXY},
          "Fit the x/y rotation before determining the rotation" },
        { "-mw", FALSE, etBOOL, {&bMW},
          "Use mass weighted fitting" }
    };
    FILE       *out;
    t_trxstatus *status;
    t_topology top;
    int        ePBC;
    rvec       *x_ref,*x;
    matrix     box,R;
    real       t;
    int        natoms,i;
    char       *grpname,title[256];
    int        gnx;
    gmx_rmpbc_t  gpbc=NULL;
    atom_id    *index;
    output_env_t oenv;
    real       *w_rls;
    const char *leg[]  = { "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz" }; 
#define NLEG asize(leg) 
    t_filenm fnm[] = {
        { efTRX, "-f",   NULL,       ffREAD }, 
        { efTPS, NULL,   NULL,       ffREAD },
        { efNDX, NULL,   NULL,       ffOPTRD },
        { efXVG, NULL,   "rotmat",   ffWRITE }
    }; 
#define NFILE asize(fnm) 
    
    CopyRight(stderr,argv[0]);
    
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv); 
    
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x_ref,NULL,box,bMW);

    gpbc = gmx_rmpbc_init(&top.idef,ePBC,top.atoms.nr,box);
    
    gmx_rmpbc(gpbc,top.atoms.nr,box,x_ref);
    
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

    if (reffit[0][0] != 'n')
    {
        get_refx(oenv,ftp2fn(efTRX,NFILE,fnm),reffit[0][2]=='z' ? 3 : 2,skip,
                 gnx,index,bMW,&top,ePBC,x_ref);
    }

    natoms = read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 

    snew(w_rls,natoms);
    for(i=0; i<gnx; i++)
    {
        if (index[i] >= natoms)
        {
            gmx_fatal(FARGS,"Atom index (%d) is larger than the number of atoms in the trajecory (%d)",index[i]+1,natoms);
        }
        w_rls[index[i]] = (bMW ? top.atoms.atom[index[i]].m : 1.0);
    }

    if (reffit[0][0] == 'n')
    {
        reset_x(gnx,index,natoms,NULL,x_ref,w_rls);
    }
    
    out = xvgropen(ftp2fn(efXVG,NFILE,fnm), 
                   "Fit matrix","Time (ps)","",oenv); 
    xvgr_legend(out,NLEG,leg,oenv);
    
    do
    {
        gmx_rmpbc(gpbc,natoms,box,x);

        reset_x(gnx,index,natoms,NULL,x,w_rls);

        if (bFitXY)
        {
            do_fit_ndim(2,natoms,w_rls,x_ref,x);
        }

        calc_fit_R(DIM,natoms,w_rls,x_ref,x,R);
          
        fprintf(out,
                "%7g %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                t,
                R[XX][XX],R[XX][YY],R[XX][ZZ],
                R[YY][XX],R[YY][YY],R[YY][ZZ],
                R[ZZ][XX],R[ZZ][YY],R[ZZ][ZZ]);
    }
    while(read_next_x(oenv,status,&t,natoms,x,box));

    gmx_rmpbc_done(gpbc);

    close_trj(status);
    
    ffclose(out);
    
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");
    
    thanx(stderr);
    
    return 0;
}
