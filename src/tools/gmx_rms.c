/*
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

#include "smalloc.h"
#include "math.h"
#include "macros.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "vec.h"
#include "index.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "princ.h"
#include "rmpbc.h"
#include "do_fit.h"
#include "matio.h"
#include "tpxio.h"
#include "cmat.h"
#include "viewit.h"
#include "gmx_ana.h"


static void norm_princ(t_atoms *atoms, int isize, atom_id *index, int natoms,
                       rvec *x)
{
    int i, m;
    rvec princ, vec;

    /* equalize principal components: */
    /* orient principal axes, get principal components */
    orient_princ(atoms, isize, index, natoms, x, NULL, princ);

    /* calc our own principal components */
    clear_rvec(vec);
    for (m = 0; m < DIM; m++)
    {
        for (i = 0; i < isize; i++)
            vec[m] += sqr(x[index[i]][m]);
        vec[m] = sqrt(vec[m] / isize);
        /* calculate scaling constants */
        vec[m] = 1 / (sqrt(3) * vec[m]);
    }

    /* scale coordinates */
    for (i = 0; i < natoms; i++)
    {
        for (m = 0; m < DIM; m++)
        {
            x[i][m] *= vec[m];
        }
    }
}

int gmx_rms(int argc, char *argv[])
{
    const char
        *desc[] =
            {
                "g_rms compares two structures by computing the root mean square",
                "deviation (RMSD), the size-independent 'rho' similarity parameter",
                "(rho) or the scaled rho (rhosc), ",
                "see Maiorov & Crippen, PROTEINS [BB]22[bb], 273 (1995).",
                "This is selected by [TT]-what[tt].[PAR]"

                    "Each structure from a trajectory ([TT]-f[tt]) is compared to a",
                "reference structure. The reference structure",
                "is taken from the structure file ([TT]-s[tt]).[PAR]",

                "With option [TT]-mir[tt] also a comparison with the mirror image of",
                "the reference structure is calculated.",
                "This is useful as a reference for 'significant' values, see",
                "Maiorov & Crippen, PROTEINS [BB]22[bb], 273 (1995).[PAR]",

                "Option [TT]-prev[tt] produces the comparison with a previous frame",
                "the specified number of frames ago.[PAR]",

                "Option [TT]-m[tt] produces a matrix in [TT].xpm[tt] format of",
                "comparison values of each structure in the trajectory with respect to",
                "each other structure. This file can be visualized with for instance",
                "[TT]xv[tt] and can be converted to postscript with [TT]xpm2ps[tt].[PAR]",

                "Option [TT]-fit[tt] controls the least-squares fitting of",
                "the structures on top of each other: complete fit (rotation and",
                "translation), translation only, or no fitting at all.[PAR]",

                "Option [TT]-mw[tt] controls whether mass weighting is done or not.",
                "If you select the option (default) and ",
                "supply a valid tpr file masses will be taken from there, ",
                "otherwise the masses will be deduced from the atommass.dat file in",
                "the GROMACS library directory. This is fine for proteins but not",
                "necessarily for other molecules. A default mass of 12.011 amu (Carbon)",
                "is assigned to unknown atoms. You can check whether this happend by",
                "turning on the [TT]-debug[tt] flag and inspecting the log file.[PAR]",

                "With [TT]-f2[tt], the 'other structures' are taken from a second",
                "trajectory, this generates a comparison matrix of one trajectory",
                "versus the other.[PAR]",

                "Option [TT]-bin[tt] does a binary dump of the comparison matrix.[PAR]",

                "Option [TT]-bm[tt] produces a matrix of average bond angle deviations",
                "analogously to the [TT]-m[tt] option. Only bonds between atoms in the",
                "comparison group are considered." };
    static gmx_bool bPBC = TRUE, bFitAll = TRUE, bSplit = FALSE;
    static gmx_bool bDeltaLog = FALSE;
    static int prev = 0, freq = 1, freq2 = 1, nlevels = 80, avl = 0;
    static real rmsd_user_max = -1, rmsd_user_min = -1, bond_user_max = -1,
        bond_user_min = -1, delta_maxy = 0.0;
    /* strings and things for selecting difference method */
    enum
    {
        ewSel, ewRMSD, ewRho, ewRhoSc, ewNR
    };
    int ewhat;
    const char *what[ewNR + 1] =
        { NULL, "rmsd", "rho", "rhosc", NULL };
    const char *whatname[ewNR] =
        { NULL, "RMSD", "Rho", "Rho sc" };
    const char *whatlabel[ewNR] =
        { NULL, "RMSD (nm)", "Rho", "Rho sc" };
    const char *whatxvgname[ewNR] =
        { NULL, "RMSD", "\\8r\\4", "\\8r\\4\\ssc\\N" };
    const char *whatxvglabel[ewNR] =
        { NULL, "RMSD (nm)", "\\8r\\4", "\\8r\\4\\ssc\\N" };
    /* strings and things for fitting methods */
    enum
    {
        efSel, efFit, efReset, efNone, efNR
    };
    int efit;
    const char *fit[efNR + 1] =
        { NULL, "rot+trans", "translation", "none", NULL };
    const char *fitgraphlabel[efNR + 1] =
        { NULL, "lsq fit", "translational fit", "no fit" };
    static int nrms = 1;
    static gmx_bool bMassWeighted = TRUE;
    t_pargs pa[] =
        {
            { "-what", FALSE, etENUM,
                { what }, "Structural difference measure" },
            { "-pbc", FALSE, etBOOL,
                { &bPBC }, "PBC check" },
            { "-fit", FALSE, etENUM,
                { fit }, "Fit to reference structure" },
            { "-prev", FALSE, etINT,
                { &prev }, "Compare with previous frame" },
            { "-split", FALSE, etBOOL,
                { &bSplit }, "Split graph where time is zero" },
            { "-fitall", FALSE, etBOOL,
                { &bFitAll }, "HIDDENFit all pairs of structures in matrix" },
            { "-skip", FALSE, etINT,
                { &freq }, "Only write every nr-th frame to matrix" },
            { "-skip2", FALSE, etINT,
                { &freq2 }, "Only write every nr-th frame to matrix" },
            { "-max", FALSE, etREAL,
                { &rmsd_user_max }, "Maximum level in comparison matrix" },
            { "-min", FALSE, etREAL,
                { &rmsd_user_min }, "Minimum level in comparison matrix" },
            { "-bmax", FALSE, etREAL,
                { &bond_user_max }, "Maximum level in bond angle matrix" },
            { "-bmin", FALSE, etREAL,
                { &bond_user_min }, "Minimum level in bond angle matrix" },
            { "-mw", FALSE, etBOOL,
                { &bMassWeighted }, "Use mass weighting for superposition" },
            { "-nlevels", FALSE, etINT,
                { &nlevels }, "Number of levels in the matrices" },
            { "-ng", FALSE, etINT,
                { &nrms }, "Number of groups to compute RMS between" },
            { "-dlog", FALSE, etBOOL,
                { &bDeltaLog },
                "HIDDENUse a log x-axis in the delta t matrix" },
            { "-dmax", FALSE, etREAL,
                { &delta_maxy }, "HIDDENMaximum level in delta matrix" },
            { "-aver", FALSE, etINT,
                { &avl },
                "HIDDENAverage over this distance in the RMSD matrix" } };
    int natoms_trx, natoms_trx2, natoms;
    int i, j, k, m, teller, teller2, tel_mat, tel_mat2;
#define NFRAME 5000
    int maxframe = NFRAME, maxframe2 = NFRAME;
    real t, *w_rls, *w_rms, *w_rls_m = NULL, *w_rms_m = NULL;
    gmx_bool bNorm, bAv, bFreq2, bFile2, bMat, bBond, bDelta, bMirror, bMass;
    gmx_bool bFit, bReset;
    t_topology top;
    int ePBC;
    t_iatom *iatom = NULL;

    matrix box;
    rvec *x, *xp, *xm = NULL, **mat_x = NULL, **mat_x2, *mat_x2_j = NULL, vec1,
        vec2;
    t_trxstatus *status;
    char buf[256], buf2[256];
    int ncons = 0;
    FILE *fp;
    real rlstot = 0, **rls, **rlsm = NULL, *time, *time2, *rlsnorm = NULL,
        **rmsd_mat = NULL, **bond_mat = NULL, *axis, *axis2, *del_xaxis,
        *del_yaxis, rmsd_max, rmsd_min, rmsd_avg, bond_max, bond_min, ang;
    real **rmsdav_mat = NULL, av_tot, weight, weight_tot;
    real **delta = NULL, delta_max, delta_scalex = 0, delta_scaley = 0,
        *delta_tot;
    int delta_xsize = 0, del_lev = 100, mx, my, abs_my;
    gmx_bool bA1, bA2, bPrev, bTop, *bInMat = NULL;
    int ifit, *irms, ibond = 0, *ind_bond1 = NULL, *ind_bond2 = NULL, n_ind_m =
        0;
    atom_id *ind_fit, **ind_rms, *ind_m = NULL, *rev_ind_m = NULL, *ind_rms_m =
        NULL;
    char *gn_fit, **gn_rms;
    t_rgb rlo, rhi;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc=NULL;

    t_filenm fnm[] =
        {
            { efTPS, NULL, NULL, ffREAD },
            { efTRX, "-f", NULL, ffREAD },
            { efTRX, "-f2", NULL, ffOPTRD },
            { efNDX, NULL, NULL, ffOPTRD },
            { efXVG, NULL, "rmsd", ffWRITE },
            { efXVG, "-mir", "rmsdmir", ffOPTWR },
            { efXVG, "-a", "avgrp", ffOPTWR },
            { efXVG, "-dist", "rmsd-dist", ffOPTWR },
            { efXPM, "-m", "rmsd", ffOPTWR },
            { efDAT, "-bin", "rmsd", ffOPTWR },
            { efXPM, "-bm", "bond", ffOPTWR } };
#define NFILE asize(fnm)

    CopyRight(stderr, argv[0]);
    parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW
        | PCA_BE_NICE, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    /* parse enumerated options: */
    ewhat = nenum(what);
    if (ewhat == ewRho || ewhat == ewRhoSc)
        please_cite(stdout, "Maiorov95");
    efit = nenum(fit);
    bFit = efit == efFit;
    bReset = efit == efReset;
    if (bFit)
    {
        bReset = TRUE; /* for fit, reset *must* be set */
    }
    else
    {
        bFitAll = FALSE;
    }

    /* mark active cmdline options */
    bMirror = opt2bSet("-mir", NFILE, fnm); /* calc RMSD vs mirror of ref. */
    bFile2 = opt2bSet("-f2", NFILE, fnm);
    bMat = opt2bSet("-m", NFILE, fnm);
    bBond = opt2bSet("-bm", NFILE, fnm);
    bDelta = (delta_maxy > 0); /* calculate rmsd vs delta t matrix from *
     *	your RMSD matrix (hidden option       */
    bNorm = opt2bSet("-a", NFILE, fnm);
    bFreq2 = opt2parg_bSet("-skip2", asize(pa), pa);
    if (freq <= 0)
    {
        fprintf(stderr, "The number of frames to skip is <= 0. "
            "Writing out all frames.\n\n");
        freq = 1;
    }
    if (!bFreq2)
    {
        freq2 = freq;
    }
    else if (bFile2 && freq2 <= 0)
    {
        fprintf(stderr,
                "The number of frames to skip in second trajectory is <= 0.\n"
                    "  Writing out all frames.\n\n");
        freq2 = 1;
    }

    bPrev = (prev > 0);
    if (bPrev)
    {
        prev = abs(prev);
        if (freq != 1)
            fprintf(stderr, "WARNING: option -skip also applies to -prev\n");
    }

    if (bFile2 && !bMat && !bBond)
    {
        fprintf(
                stderr,
                "WARNING: second trajectory (-f2) useless when not calculating matrix (-m/-bm),\n"
                    "         will not read from %s\n", opt2fn("-f2", NFILE,
                                                               fnm));
        bFile2 = FALSE;
    }

    if (bDelta)
    {
        bMat = TRUE;
        if (bFile2)
        {
            fprintf(stderr,
                    "WARNING: second trajectory (-f2) useless when making delta matrix,\n"
                        "         will not read from %s\n", opt2fn("-f2",
                                                                   NFILE, fnm));
            bFile2 = FALSE;
        }
    }

    bTop = read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &xp,
                         NULL, box, TRUE);
    snew(w_rls,top.atoms.nr);
    snew(w_rms,top.atoms.nr);

    if (!bTop && bBond)
    {
        fprintf(stderr,
                "WARNING: Need a run input file for bond angle matrix,\n"
                    "         will not calculate bond angle matrix.\n");
        bBond = FALSE;
    }

    if (bReset)
    {
        fprintf(stderr, "Select group for %s fit\n", bFit ? "least squares"
            : "translational");
        get_index(&(top.atoms), ftp2fn_null(efNDX, NFILE, fnm), 1, &ifit,
                  &ind_fit, &gn_fit);
    }
    else
        ifit = 0;

    if (bReset)
    {
        if (bFit && ifit < 3)
            gmx_fatal(FARGS,"Need >= 3 points to fit!\n" );

        bMass = FALSE;
        for(i=0; i<ifit; i++)
        {
            if (bMassWeighted)
                w_rls[ind_fit[i]] = top.atoms.atom[ind_fit[i]].m;
            else
                w_rls[ind_fit[i]] = 1;
            bMass = bMass || (top.atoms.atom[ind_fit[i]].m != 0);
        }
        if (!bMass)
        {
            fprintf(stderr,"All masses in the fit group are 0, using masses of 1\n");
            for(i=0; i<ifit; i++)
            {
                w_rls[ind_fit[i]] = 1;
            }
        }
    }

    if (bMat || bBond)
        nrms=1;

    snew(gn_rms,nrms);
    snew(ind_rms,nrms);
    snew(irms,nrms);

    fprintf(stderr,"Select group%s for %s calculation\n",
            (nrms>1) ? "s" : "",whatname[ewhat]);
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
              nrms,irms,ind_rms,gn_rms);

    if (bNorm)
    {
        snew(rlsnorm,irms[0]);
    }
    snew(rls,nrms);
    for(j=0; j<nrms; j++)
        snew(rls[j],maxframe);
    if (bMirror)
    {
        snew(rlsm,nrms);
        for(j=0; j<nrms; j++)
            snew(rlsm[j],maxframe);
    }
    snew(time,maxframe);
    for(j=0; j<nrms; j++)
    {
        bMass = FALSE;
        for(i=0; i<irms[j]; i++)
        {
            if (bMassWeighted)
                w_rms[ind_rms[j][i]] = top.atoms.atom[ind_rms[j][i]].m;
            else
                w_rms[ind_rms[j][i]] = 1.0;
            bMass = bMass || (top.atoms.atom[ind_rms[j][i]].m != 0);
        }
        if (!bMass) {
            fprintf(stderr,"All masses in group %d are 0, using masses of 1\n",j);
            for(i=0; i<irms[j]; i++)
                w_rms[ind_rms[j][i]] = 1;
        }
    }
    /* Prepare reference frame */
    if (bPBC) {
      gpbc = gmx_rmpbc_init(&top.idef,ePBC,top.atoms.nr,box);
      gmx_rmpbc(gpbc,top.atoms.nr,box,xp);
    }
    if (bReset)
        reset_x(ifit,ind_fit,top.atoms.nr,NULL,xp,w_rls);
    if (bMirror) {
        /* generate reference structure mirror image: */
        snew(xm, top.atoms.nr);
        for(i=0; i<top.atoms.nr; i++) {
            copy_rvec(xp[i],xm[i]);
            xm[i][XX] = -xm[i][XX];
        }
    }
    if (ewhat==ewRhoSc)
        norm_princ(&top.atoms, ifit, ind_fit, top.atoms.nr, xp);

    /* read first frame */
    natoms_trx=read_first_x(oenv,&status,opt2fn("-f",NFILE,fnm),&t,&x,box);
    if (natoms_trx != top.atoms.nr) 
        fprintf(stderr,
                "\nWARNING: topology has %d atoms, whereas trajectory has %d\n",
                top.atoms.nr,natoms_trx);
    natoms = min(top.atoms.nr,natoms_trx);
    if (bMat || bBond || bPrev) {
        snew(mat_x,NFRAME);

        if (bPrev)
            /* With -prev we use all atoms for simplicity */
            n_ind_m = natoms;
        else {
            /* Check which atoms we need (fit/rms) */
            snew(bInMat,natoms);
            for(i=0; i<ifit; i++)
                bInMat[ind_fit[i]] = TRUE;
            n_ind_m=ifit;
            for(i=0; i<irms[0]; i++)
                if (!bInMat[ind_rms[0][i]]) {
                    bInMat[ind_rms[0][i]] = TRUE;
                    n_ind_m++;
                }
        }
        /* Make an index of needed atoms */
        snew(ind_m,n_ind_m);
        snew(rev_ind_m,natoms);
        j = 0;
        for(i=0; i<natoms; i++)
            if (bPrev || bInMat[i]) {
                ind_m[j] = i;
                rev_ind_m[i] = j;
                j++;
            }
        snew(w_rls_m,n_ind_m);
        snew(ind_rms_m,irms[0]);
        snew(w_rms_m,n_ind_m);
        for(i=0; i<ifit; i++)
            w_rls_m[rev_ind_m[ind_fit[i]]] = w_rls[ind_fit[i]];
        for(i=0; i<irms[0]; i++) {
            ind_rms_m[i] = rev_ind_m[ind_rms[0][i]];
            w_rms_m[ind_rms_m[i]] = w_rms[ind_rms[0][i]];
        }
        sfree(bInMat);
    }

    if (bBond) {
        ncons = 0;
        for(k=0; k<F_NRE; k++)
            if (IS_CHEMBOND(k)) {
                iatom  = top.idef.il[k].iatoms;
                ncons += top.idef.il[k].nr/3;
            }
        fprintf(stderr,"Found %d bonds in topology\n",ncons);
        snew(ind_bond1,ncons);
        snew(ind_bond2,ncons);
        ibond=0;
        for(k=0; k<F_NRE; k++)
            if (IS_CHEMBOND(k)) {
                iatom = top.idef.il[k].iatoms;
                ncons = top.idef.il[k].nr/3;
                for (i=0; i<ncons; i++) {
                    bA1=FALSE; 
                    bA2=FALSE;
                    for (j=0; j<irms[0]; j++) {
                        if (iatom[3*i+1] == ind_rms[0][j]) bA1=TRUE; 
                        if (iatom[3*i+2] == ind_rms[0][j]) bA2=TRUE;
                    }
                    if (bA1 && bA2) {
                        ind_bond1[ibond] = rev_ind_m[iatom[3*i+1]];
                        ind_bond2[ibond] = rev_ind_m[iatom[3*i+2]];
                        ibond++;
                    }
                }
            }
        fprintf(stderr,"Using %d bonds for bond angle matrix\n",ibond);
        if (ibond==0)
            gmx_fatal(FARGS,"0 bonds found");
    }

    /* start looping over frames: */
    tel_mat = 0;
    teller = 0;
    do {
        if (bPBC) 
	  gmx_rmpbc(gpbc,natoms,box,x);

        if (bReset)
            reset_x(ifit,ind_fit,natoms,NULL,x,w_rls);
        if (ewhat==ewRhoSc)
            norm_princ(&top.atoms, ifit, ind_fit, natoms, x);

        if (bFit)
            /*do the least squares fit to original structure*/
            do_fit(natoms,w_rls,xp,x);

        if (teller % freq == 0) {
            /* keep frame for matrix calculation */
            if (bMat || bBond || bPrev) {
                if (tel_mat >= NFRAME) 
                    srenew(mat_x,tel_mat+1);
                snew(mat_x[tel_mat],n_ind_m);
                for (i=0;i<n_ind_m;i++)
                    copy_rvec(x[ind_m[i]],mat_x[tel_mat][i]);
            }
            tel_mat++;
        }

        /*calculate energy of root_least_squares*/
        if (bPrev) {
            j=tel_mat-prev-1;
            if (j<0)
                j=0;
            for (i=0;i<n_ind_m;i++)
                copy_rvec(mat_x[j][i],xp[ind_m[i]]);
            if (bReset)
                reset_x(ifit,ind_fit,natoms,NULL,xp,w_rls);
            if (bFit)
                do_fit(natoms,w_rls,x,xp);
        }    
        for(j=0; (j<nrms); j++) 
            rls[j][teller] = 
                calc_similar_ind(ewhat!=ewRMSD,irms[j],ind_rms[j],w_rms,x,xp);
        if (bNorm) {
            for(j=0; (j<irms[0]); j++)
                rlsnorm[j] += 
                    calc_similar_ind(ewhat!=ewRMSD,1,&(ind_rms[0][j]),w_rms,x,xp);
        } 

        if (bMirror) {
            if (bFit)
                /*do the least squares fit to mirror of original structure*/
                do_fit(natoms,w_rls,xm,x);

            for(j=0; j<nrms; j++)
                rlsm[j][teller] = 
                    calc_similar_ind(ewhat!=ewRMSD,irms[j],ind_rms[j],w_rms,x,xm);
        }
        time[teller]=output_env_conv_time(oenv,t);

        teller++;
        if (teller >= maxframe) {
            maxframe +=NFRAME;
            srenew(time,maxframe);
            for(j=0; (j<nrms); j++) 
                srenew(rls[j],maxframe);
            if (bMirror)
                for(j=0; (j<nrms); j++) 
                    srenew(rlsm[j],maxframe);
        }
    } while (read_next_x(oenv,status,&t,natoms_trx,x,box));
    close_trj(status);

    if (bFile2) {
        snew(time2,maxframe2);

        fprintf(stderr,"\nWill read second trajectory file\n");
        snew(mat_x2,NFRAME);
        natoms_trx2 =
	  read_first_x(oenv,&status,opt2fn("-f2",NFILE,fnm),&t,&x,box);
        if ( natoms_trx2 != natoms_trx )
            gmx_fatal(FARGS,
                      "Second trajectory (%d atoms) does not match the first one"
                      " (%d atoms)", natoms_trx2, natoms_trx);
        tel_mat2 = 0;
        teller2 = 0;
        do {
            if (bPBC) 
	      gmx_rmpbc(gpbc,natoms,box,x);

            if (bReset)
                reset_x(ifit,ind_fit,natoms,NULL,x,w_rls);
            if (ewhat==ewRhoSc)
                norm_princ(&top.atoms, ifit, ind_fit, natoms, x);

            if (bFit)
                /*do the least squares fit to original structure*/
                do_fit(natoms,w_rls,xp,x);

            if (teller2 % freq2 == 0) {
                /* keep frame for matrix calculation */
                if (bMat) {
                    if (tel_mat2 >= NFRAME) 
                        srenew(mat_x2,tel_mat2+1);
                    snew(mat_x2[tel_mat2],n_ind_m);
                    for (i=0;i<n_ind_m;i++)
                        copy_rvec(x[ind_m[i]],mat_x2[tel_mat2][i]);
                }
                tel_mat2++;
            }

            time2[teller2]=output_env_conv_time(oenv,t);

            teller2++;
            if (teller2 >= maxframe2) {
                maxframe2 +=NFRAME;
                srenew(time2,maxframe2);
            }
        } while (read_next_x(oenv,status,&t,natoms_trx2,x,box));
        close_trj(status);
    } else {
        mat_x2=mat_x;
        time2=time;
        tel_mat2=tel_mat;
        freq2=freq;
    }
    gmx_rmpbc_done(gpbc);

    if (bMat || bBond) {
        /* calculate RMS matrix */
        fprintf(stderr,"\n");
        if (bMat) {
            fprintf(stderr,"Building %s matrix, %dx%d elements\n",
                    whatname[ewhat],tel_mat,tel_mat2);
            snew(rmsd_mat,tel_mat);
        }
        if (bBond) {
            fprintf(stderr,"Building bond angle matrix, %dx%d elements\n",
                    tel_mat,tel_mat2);
            snew(bond_mat,tel_mat);
        }
        snew(axis,tel_mat);
        snew(axis2,tel_mat2);
        rmsd_max=0;
        if (bFile2)
            rmsd_min=1e10;
        else
            rmsd_min=0;
        rmsd_avg=0;
        bond_max=0;
        bond_min=1e10;
        for(j=0; j<tel_mat2; j++)
            axis2[j]=time2[freq2*j];
        if (bDelta) {
            if (bDeltaLog) {
                delta_scalex=8.0/log(2.0);
                delta_xsize=(int)(log(tel_mat/2)*delta_scalex+0.5)+1;
            }
            else {
                delta_xsize=tel_mat/2;
            }
            delta_scaley=1.0/delta_maxy;
            snew(delta,delta_xsize);
            for(j=0; j<delta_xsize; j++)
                snew(delta[j],del_lev+1);
            if (avl > 0) {
                snew(rmsdav_mat,tel_mat);
                for(j=0; j<tel_mat; j++)
                    snew(rmsdav_mat[j],tel_mat);
            }
        }

        if (bFitAll)
            snew(mat_x2_j,natoms);
        for(i=0; i<tel_mat; i++) {
            axis[i]=time[freq*i];
            fprintf(stderr,"\r element %5d; time %5.2f  ",i,axis[i]);
            if (bMat) snew(rmsd_mat[i],tel_mat2);
            if (bBond) snew(bond_mat[i],tel_mat2); 
            for(j=0; j<tel_mat2; j++) {
                if (bFitAll) {
                    for (k=0;k<n_ind_m;k++)
                        copy_rvec(mat_x2[j][k],mat_x2_j[k]);
                    do_fit(n_ind_m,w_rls_m,mat_x[i],mat_x2_j);
                } else
                    mat_x2_j=mat_x2[j];
                if (bMat) {
                    if (bFile2 || (i<j)) {
                        rmsd_mat[i][j] =
                            calc_similar_ind(ewhat!=ewRMSD,irms[0],ind_rms_m,
                                             w_rms_m,mat_x[i],mat_x2_j);
                        if (rmsd_mat[i][j] > rmsd_max) rmsd_max=rmsd_mat[i][j];
                        if (rmsd_mat[i][j] < rmsd_min) rmsd_min=rmsd_mat[i][j];
                        rmsd_avg += rmsd_mat[i][j];
                    } else
                        rmsd_mat[i][j]=rmsd_mat[j][i];
                }
                if (bBond) {
                    if (bFile2 || (i<=j)) {
                        ang=0.0;
                        for(m=0;m<ibond;m++) {
                            rvec_sub(mat_x[i][ind_bond1[m]],mat_x[i][ind_bond2[m]],vec1);
                            rvec_sub(mat_x2_j[ind_bond1[m]],mat_x2_j[ind_bond2[m]],vec2);
                            ang += acos(cos_angle(vec1,vec2));
                        }
                        bond_mat[i][j]=ang*180.0/(M_PI*ibond);
                        if (bond_mat[i][j] > bond_max) bond_max=bond_mat[i][j];
                        if (bond_mat[i][j] < bond_min) bond_min=bond_mat[i][j];
                    } 
                    else
                        bond_mat[i][j]=bond_mat[j][i];
                }
            }
        }
        if (bFile2)
            rmsd_avg /= tel_mat*tel_mat2;
        else
            rmsd_avg /= tel_mat*(tel_mat - 1)/2;
        if (bMat && (avl > 0)) {
            rmsd_max=0.0;
            rmsd_min=0.0;
            rmsd_avg=0.0;
            for(j=0; j<tel_mat-1; j++) {
                for(i=j+1; i<tel_mat; i++) {
                    av_tot=0;
                    weight_tot=0;
                    for (my=-avl; my<=avl; my++) {
                        if ((j+my>=0) && (j+my<tel_mat)) {
                            abs_my = abs(my);
                            for (mx=-avl; mx<=avl; mx++) {
                                if ((i+mx>=0) && (i+mx<tel_mat)) {
                                    weight = (real)(avl+1-max(abs(mx),abs_my));
                                    av_tot += weight*rmsd_mat[i+mx][j+my];
                                    weight_tot+=weight;
                                }
                            }
                        }
                    }
                    rmsdav_mat[i][j] = av_tot/weight_tot;
                    rmsdav_mat[j][i] = rmsdav_mat[i][j];
                    if (rmsdav_mat[i][j] > rmsd_max) rmsd_max=rmsdav_mat[i][j];
                }
            }
            rmsd_mat=rmsdav_mat;
        }

        if (bMat) {
            fprintf(stderr,"\n%s: Min %f, Max %f, Avg %f\n",
                    whatname[ewhat],rmsd_min,rmsd_max,rmsd_avg);
            rlo.r = 1; rlo.g = 1; rlo.b = 1;
            rhi.r = 0; rhi.g = 0; rhi.b = 0;
            if (rmsd_user_max != -1) rmsd_max=rmsd_user_max;
            if (rmsd_user_min != -1) rmsd_min=rmsd_user_min;
            if ((rmsd_user_max !=  -1) || (rmsd_user_min != -1))
                fprintf(stderr,"Min and Max value set to resp. %f and %f\n",
                        rmsd_min,rmsd_max);
            sprintf(buf,"%s %s matrix",gn_rms[0],whatname[ewhat]);
            write_xpm(opt2FILE("-m",NFILE,fnm,"w"),0,buf,whatlabel[ewhat],
                      output_env_get_time_label(oenv),output_env_get_time_label(oenv),tel_mat,tel_mat2,
                      axis,axis2, rmsd_mat,rmsd_min,rmsd_max,rlo,rhi,&nlevels);
            /* Print the distribution of RMSD values */
            if (opt2bSet("-dist",NFILE,fnm)) 
                low_rmsd_dist(opt2fn("-dist",NFILE,fnm),rmsd_max,tel_mat,rmsd_mat,oenv);

            if (bDelta) {
                snew(delta_tot,delta_xsize);
                for(j=0; j<tel_mat-1; j++) {
                    for(i=j+1; i<tel_mat; i++) {
                        mx=i-j ;
                        if (mx < tel_mat/2) {
                            if (bDeltaLog) 
                                mx=(int)(log(mx)*delta_scalex+0.5);
                            my=(int)(rmsd_mat[i][j]*delta_scaley*del_lev+0.5);
                            delta_tot[mx] += 1.0;
                            if ((rmsd_mat[i][j]>=0) && (rmsd_mat[i][j]<=delta_maxy))
                                delta[mx][my] += 1.0;
                        }
                    }
                }
                delta_max=0;
                for(i=0; i<delta_xsize; i++) {
                    if (delta_tot[i] > 0.0) {
                        delta_tot[i]=1.0/delta_tot[i];
                        for(j=0; j<=del_lev; j++) {
                            delta[i][j] *= delta_tot[i];
                            if (delta[i][j] > delta_max)
                                delta_max=delta[i][j];
                        }
                    }
                }
                fprintf(stderr,"Maximum in delta matrix: %f\n",delta_max);
                snew(del_xaxis,delta_xsize);
                snew(del_yaxis,del_lev+1);
                for (i=0; i<delta_xsize; i++)
                    del_xaxis[i]=axis[i]-axis[0];
                for (i=0; i<del_lev+1; i++)
                    del_yaxis[i]=delta_maxy*i/del_lev;
                sprintf(buf,"%s %s vs. delta t",gn_rms[0],whatname[ewhat]);
                fp = ffopen("delta.xpm","w");
                write_xpm(fp,0,buf,"density",output_env_get_time_label(oenv),whatlabel[ewhat],
                          delta_xsize,del_lev+1,del_xaxis,del_yaxis,
                          delta,0.0,delta_max,rlo,rhi,&nlevels);
                ffclose(fp);
            }
            if (opt2bSet("-bin",NFILE,fnm)) {
                /* NB: File must be binary if we use fwrite */
                fp=ftp2FILE(efDAT,NFILE,fnm,"wb");
                for(i=0;i<tel_mat;i++) 
                    if(fwrite(rmsd_mat[i],sizeof(**rmsd_mat),tel_mat2,fp) != tel_mat2)
                    {
                        gmx_fatal(FARGS,"Error writing to output file");
                    }
                ffclose(fp);
            }
        }
        if (bBond) {
            fprintf(stderr,"\nMin. angle: %f, Max. angle: %f\n",bond_min,bond_max);
            if (bond_user_max != -1) bond_max=bond_user_max;
            if (bond_user_min != -1) bond_min=bond_user_min;
            if ((bond_user_max !=  -1) || (bond_user_min != -1))
                fprintf(stderr,"Bond angle Min and Max set to:\n"
                        "Min. angle: %f, Max. angle: %f\n",bond_min,bond_max);
            rlo.r = 1; rlo.g = 1; rlo.b = 1;
            rhi.r = 0; rhi.g = 0; rhi.b = 0;
            sprintf(buf,"%s av. bond angle deviation",gn_rms[0]);
            write_xpm(opt2FILE("-bm",NFILE,fnm,"w"),0,buf,"degrees",
                      output_env_get_time_label(oenv),output_env_get_time_label(oenv),tel_mat,tel_mat2,
                      axis,axis2, bond_mat,bond_min,bond_max,rlo,rhi,&nlevels);
        }
    }

    bAv=opt2bSet("-a",NFILE,fnm);

    /* Write the RMSD's to file */
    if (!bPrev)
        sprintf(buf,"%s",whatxvgname[ewhat]);
    else
        sprintf(buf,"%s with frame %g %s ago",whatxvgname[ewhat],
                time[prev*freq]-time[0], output_env_get_time_label(oenv));
    fp=xvgropen(opt2fn("-o",NFILE,fnm),buf,output_env_get_xvgr_tlabel(oenv),
                whatxvglabel[ewhat], oenv);
    if (output_env_get_print_xvgr_codes(oenv))
        fprintf(fp,"@ subtitle \"%s%s after %s%s%s\"\n",
                (nrms==1)?"":"of "    , gn_rms[0], fitgraphlabel[efit],
                    bFit     ?" to ":""   , bFit?gn_fit:"");
    if (nrms != 1)
        xvgr_legend(fp,nrms,(const char**)gn_rms,oenv);
    for(i=0; (i<teller); i++) {
        if ( bSplit && i>0 && 
            abs(time[bPrev ? freq*i : i]/output_env_get_time_factor(oenv))<1e-5 ) 
        {
            fprintf(fp,"&\n");
        }
        fprintf(fp,"%12.7f",time[bPrev ? freq*i : i]);
        for(j=0; (j<nrms); j++) {
            fprintf(fp," %12.7f",rls[j][i]);
            if (bAv)
                rlstot+=rls[j][i];
        }
        fprintf(fp,"\n");
    }
    ffclose(fp);
    
    if (bMirror) {
        /* Write the mirror RMSD's to file */
        sprintf(buf,"%s with Mirror",whatxvgname[ewhat]);
        sprintf(buf2,"Mirror %s",whatxvglabel[ewhat]);
        fp=xvgropen(opt2fn("-mir",NFILE,fnm), buf, output_env_get_xvgr_tlabel(oenv), 
                    buf2,oenv);
        if (nrms == 1) {
            if (output_env_get_print_xvgr_codes(oenv))
                fprintf(fp,"@ subtitle \"of %s after lsq fit to mirror of %s\"\n",
                        gn_rms[0],gn_fit);
        }
        else {
            if (output_env_get_print_xvgr_codes(oenv))
                fprintf(fp,"@ subtitle \"after lsq fit to mirror %s\"\n",gn_fit);
            xvgr_legend(fp,nrms,(const char**)gn_rms,oenv);
        }
        for(i=0; (i<teller); i++) {
            if ( bSplit && i>0 && abs(time[i])<1e-5 ) 
                fprintf(fp,"&\n");
            fprintf(fp,"%12.7f",time[i]);
            for(j=0; (j<nrms); j++)
                fprintf(fp," %12.7f",rlsm[j][i]);
            fprintf(fp,"\n");
        }
        ffclose(fp);
    }

    if (bAv) {
        sprintf(buf,"Average %s",whatxvgname[ewhat]);
        sprintf(buf2,"Average %s",whatxvglabel[ewhat]);
        fp = xvgropen(opt2fn("-a",NFILE,fnm), buf, "Residue", buf2,oenv);
        for(j=0; (j<nrms); j++)
            fprintf(fp,"%10d  %10g\n",j,rlstot/teller);
        ffclose(fp);
    }

    if (bNorm) {
        fp = xvgropen("aver.xvg",gn_rms[0],"Residue",whatxvglabel[ewhat],oenv);
        for(j=0; (j<irms[0]); j++)
            fprintf(fp,"%10d  %10g\n",j,rlsnorm[j]/teller);
        ffclose(fp);
    }
    do_view(oenv,opt2fn_null("-a",NFILE,fnm),"-graphtype bar");
    do_view(oenv,opt2fn("-o",NFILE,fnm),NULL);
    do_view(oenv,opt2fn_null("-mir",NFILE,fnm),NULL);
    do_view(oenv,opt2fn_null("-m",NFILE,fnm),NULL);
    do_view(oenv,opt2fn_null("-bm",NFILE,fnm),NULL);
    do_view(oenv,opt2fn_null("-dist",NFILE,fnm),NULL);

    thanx(stderr);

    return 0;
}
