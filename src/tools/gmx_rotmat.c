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
        "in the xy-plane."
    };
    static bool bMW=TRUE;
    t_pargs pa[] = {
        { "-mw", FALSE, etBOOL, {&bMW},
          "Use mass weighted fitting" }
    };
    FILE       *out;
    int        status;
    t_topology top;
    int        ePBC;
    rvec       *x_top,*x;
    matrix     box,R;
    real       t;
    int        natoms,i;
    char       *grpname,title[256];
    int        gnx;
    atom_id    *index;
    output_env_t oenv;
    real       *w_rls;
    char *leg[]  = { "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz" }; 
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
    
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x_top,NULL,box,bMW);

    rm_pbc(&top.idef,ePBC,top.atoms.nr,box,x_top,x_top);
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

    natoms = read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 

    snew(w_rls,natoms);
    for(i=0; i<gnx; i++)
    {
        if (index[i] >= natoms) {
            gmx_fatal(FARGS,"Atom index (%d) is larger than the number of atoms in the trajecory (%d)",index[i]+1,natoms);
        }
        w_rls[index[i]] = (bMW ? top.atoms.atom[index[i]].m : 1.0);
    }

    reset_x(gnx,index,natoms,NULL,x_top,w_rls);
    
    out = xvgropen(ftp2fn(efXVG,NFILE,fnm), 
                   "Fit matrix","Time (ps)","",oenv); 
    xvgr_legend(out,NLEG,leg,oenv);
    
    do
    {
        rm_pbc(&top.idef,ePBC,natoms,box,x,x);

        reset_x(gnx,index,natoms,NULL,x,w_rls);

        calc_fit_R(DIM,natoms,w_rls,x_top,x,R);
          
        fprintf(out,
                "%7g %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                t,
                R[XX][XX],R[XX][YY],R[XX][ZZ],
                R[YY][XX],R[YY][YY],R[YY][ZZ],
                R[ZZ][XX],R[ZZ][YY],R[ZZ][ZZ]);
    }
    while(read_next_x(oenv,status,&t,natoms,x,box));

    close_trj(status);
    
    fclose(out);
    
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");
    
    thanx(stderr);
    
    return 0;
}
