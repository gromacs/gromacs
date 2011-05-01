/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
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
#include <stdio.h>
#include <math.h>

#include "confio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"
#include "correl.h"
#include "gmx_ana.h"
#include "gmx_fft.h"

enum { VACF, MVACF, DOS, DOS_SOLID, DOS_DIFF, DOS_CP, DOS_S, DOS_NR };

static real FD(real Delta,real f)
{
    return (2*pow(Delta,-4.5)*pow(f,7.5) - 
            6*pow(Delta,-3)*pow(f,5) -
            pow(Delta,-1.5)*pow(f,3.5) +
            6*pow(Delta,-1.5)*pow(f,2.5) +
            2*f - 2);
}

static real calc_fluidicity(real Delta,real tol)
{
    real fd0,fd,fd1,f,f0=0,f1=1;
    real tolmin = 1e-6;
    
    /* Since the fluidity is monotonous according to Fig. 2 in Lin2003a, 
       J. Chem. Phys. 112 (2003) p. 11792 we can use a bisection method
       to get it. */
    if (tol < tolmin) 
    {
        fprintf(stderr,"Unrealistic tolerance %g for calc_fluidity. Setting it to %g\n",tol,tolmin);
        tol=1e-6;
    }
    
    do {
        fd0 = FD(Delta,f0);
        fd1 = FD(Delta,f1);
        f = (f0+f1)*0.5;
        fd = FD(Delta,f);
        if (fd < 0)
            f0 = f;
        else if (fd > 0)
            f1 = f;
        else
            return f;
    } while ((f1-f0) > tol);
    
    return f;
}

static real wCsolid(real nu,real beta)
{
    real bhn = beta*PLANCK*nu;
    real ebn,koko;
    
    if (bhn == 0)
        return 1.0;
    else 
    {
        ebn = exp(bhn);
        koko = sqr(1-ebn);
        return sqr(bhn)*ebn/koko;
    }
}

static real wSsolid(real nu,real beta)
{
    real bhn = beta*PLANCK*nu;
    real ebn;
    
    if (bhn == 0) 
    {
        return 1;
    }
    else 
    {
        ebn = exp(-bhn);
        
        return bhn/(exp(bhn)-1) - log(1-exp(-bhn));
    }
}

int gmx_dos(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_dos[tt] computes the Density of States from a simulations.",
        "In order for this to be meaningful the velocities must be saved",
        "in the trajecotry with sufficiently high frequency such as to cover",
        "all vibrations. For flexible systems that would be around a few fs",
        "between saving. Properties based on the DoS are printed on the",
        "standard output."
    };
    const char *bugs[] = {
        "This program needs a lot of memory: total usage equals the number of atoms times 3 times number of frames times 4 (or 8 when run in double precision)."
    };
    FILE       *fp;
    t_topology top;
    int        ePBC=-1;
    t_trxframe fr;
    matrix     box;
    int        gnx;
    char       title[256];
    real       t0,t1,m;
    t_trxstatus *status;
    int        nV,teller,n_alloc,i,j,k,l,fftcode,Nmol,Natom;
    double     dt,V2sum,Vsum,V,tmass,dostot,dos2,dosabs;
    real       **c1,**dos,mi,beta,bfac,*nu,*tt,stddev,c1j;
    output_env_t oenv;
    gmx_fft_t  fft;
    double     cP,S,DiffCoeff,Delta,f,DoS0;
    double     wCdiff,wSdiff;
    
    static     gmx_bool bVerbose=TRUE,bAbsolute=FALSE,bNormalize=FALSE;
    static     gmx_bool bRecip=FALSE;
    static     real Temp=298.15,toler=1e-3;
    t_pargs pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy." },
        { "-recip", FALSE, etBOOL, {&bRecip},
          "Use cm^-1 on X-axis instead of 1/ps for DoS plots." },
        { "-abs", FALSE, etBOOL, {&bAbsolute},
          "Use the absolute value of the Fourier transform of the VACF as the Density of States. Default is to use the real component only" },
        { "-normdos", FALSE, etBOOL, {&bNormalize},
          "Normalize the DoS such that it adds up to 3N. This is a hack that should not be necessary." },
        { "-T", FALSE, etREAL, {&Temp},
          "Temperature in the simulation" },
        { "-toler", FALSE, etREAL, {&toler},
          "[HIDDEN]Tolerance when computing the fluidicity using bisection algorithm" }
    };

    t_filenm  fnm[] = {
        { efTRN, "-f",    NULL,    ffREAD  },
        { efTPX, "-s",    NULL,    ffREAD }, 
        { efNDX, NULL,    NULL,    ffOPTRD },
        { efXVG, "-vacf", "vacf",  ffWRITE },
        { efXVG, "-mvacf","mvacf", ffWRITE },
        { efXVG, "-dos",  "dos",   ffWRITE }
    };
#define NFILE asize(fnm)
    int     npargs;
    t_pargs *ppa;
    const char *DoSlegend[] = {
        "DoS(v)", "DoS(v)[Solid]", "DoS(v)[Diff]" 
    };
    
    CopyRight(stderr,argv[0]);
    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs,pa);
    parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                      NFILE,fnm,npargs,ppa,asize(desc),desc,
                      asize(bugs),bugs,&oenv);
                      
    please_cite(stdout,"Pascal2011a");
    please_cite(stdout,"Caleman2011b");
    
    read_tps_conf(ftp2fn(efTPX,NFILE,fnm),title,&top,&ePBC,NULL,NULL,box,
                  TRUE);
    V = det(box);
    tmass = 0;
    for(i=0; (i<top.atoms.nr); i++)
        tmass += top.atoms.atom[i].m;

    Natom = top.atoms.nr;
    Nmol = top.mols.nr;
    gnx = Natom*DIM;
    
    /* Correlation stuff */
    snew(c1,gnx);
    for(i=0; (i<gnx); i++)
        c1[i]=NULL;
  
    read_first_frame(oenv,&status,ftp2fn(efTRN,NFILE,fnm),&fr,TRX_NEED_V);
    t0=fr.time;
      
    n_alloc=0;
    teller=0;
    Vsum = V2sum = 0;
    nV = 0;
    do {
        if (fr.bBox) 
        {
            V = det(fr.box);
            V2sum += V*V;
            Vsum += V;
            nV++;
        }
        if (teller >= n_alloc) 
        {
            n_alloc+=100;
            for(i=0; i<gnx; i++)
                srenew(c1[i],n_alloc);
        }
        for(i=0; i<gnx; i+=DIM) 
        {
            c1[i+XX][teller] = fr.v[i/DIM][XX];
            c1[i+YY][teller] = fr.v[i/DIM][YY];
            c1[i+ZZ][teller] = fr.v[i/DIM][ZZ];
        }

        t1=fr.time;

        teller ++;
    } while (read_next_frame(oenv,status,&fr));
  
    close_trj(status);

    dt = (t1-t0)/(teller-1);
    if (nV > 0)
    {
        V = Vsum/nV;
    }
    if (bVerbose)
        printf("Going to do %d fourier transforms of length %d. Hang on.\n",
               gnx,teller);
    low_do_autocorr(NULL,oenv,NULL,teller,gnx,teller,c1,dt,eacNormal,0,FALSE,
                    FALSE,FALSE,-1,-1,0,0);
    snew(dos,DOS_NR);
    for(j=0; (j<DOS_NR); j++)
        snew(dos[j],teller+4);

    if (bVerbose)
        printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");
    for(i=0; (i<gnx); i+=DIM) 
    {
        mi = top.atoms.atom[i/DIM].m;
        for(j=0; (j<teller/2); j++) 
        {
            c1j = (c1[i+XX][j] + c1[i+YY][j] + c1[i+ZZ][j]);
            dos[VACF][j]  += c1j/Natom;
            dos[MVACF][j] += mi*c1j;
        }
    }
    fp = xvgropen(opt2fn("-vacf",NFILE,fnm),"Velocity ACF",
                  "Time (ps)","C(t)",oenv);
    snew(tt,teller/2);
    for(j=0; (j<teller/2); j++) 
    {
        tt[j] = j*dt;
        fprintf(fp,"%10g  %10g\n",tt[j],dos[VACF][j]);
    }
    fclose(fp);
    fp = xvgropen(opt2fn("-mvacf",NFILE,fnm),"Mass-weighted velocity ACF",
                  "Time (ps)","C(t)",oenv);
    for(j=0; (j<teller/2); j++) 
    {
        fprintf(fp,"%10g  %10g\n",tt[j],dos[MVACF][j]);
    }
    fclose(fp);
    
    if ((fftcode = gmx_fft_init_1d_real(&fft,teller/2,
                                        GMX_FFT_FLAG_NONE)) != 0) 
    {
        gmx_fatal(FARGS,"gmx_fft_init_1d_real returned %d",fftcode);
    }
    if ((fftcode = gmx_fft_1d_real(fft,GMX_FFT_REAL_TO_COMPLEX,
                                   (void *)dos[MVACF],(void *)dos[DOS])) != 0)
    {
        gmx_fatal(FARGS,"gmx_fft_1d_real returned %d",fftcode);
    }

    /* First compute the DoS */
    beta = 1/(Temp*BOLTZ);
    bfac = dt*beta/2;
    dos2 = 0;
    snew(nu,teller/2);
    for(j=0; (j<teller/2); j++) 
    {
        nu[j] = j/(t1-t0);
        dos2 += sqr(dos[DOS][2*j]) + sqr(dos[DOS][2*j+1]);
        if (bAbsolute)
            dos[DOS][j] = bfac*sqrt(sqr(dos[DOS][2*j]) + sqr(dos[DOS][2*j+1]));
        else
            dos[DOS][j] = bfac*dos[DOS][2*j];
    }
    /* Normalize it */
    dostot = evaluate_integral(teller/2,nu,dos[DOS],NULL,teller/2,&stddev);
    if (bNormalize) 
    {
        for(j=0; (j<teller/2); j++) 
            dos[DOS][j] *= 3*Natom/dostot;
    }
    
    /* Now analyze it */
    DoS0 = dos[DOS][0];
    
    /* Note this eqn. is incorrect in Pascal2011a! */
    Delta = ((2*DoS0/(9*Natom))*sqrt(M_PI*BOLTZ*Temp/tmass)*
             pow((Natom/V),1.0/3.0)*pow(6/M_PI,2.0/3.0));
    f = calc_fluidicity(Delta,toler);
    printf("Nmol = %d, Natom = %d, dt = %g ps, tmass = %g amu\n"
           "V = %g nm^3, Delta = %g, f = %g beta = %g (mol/kJ)\n"
           "DoS0 = %g, Dos2 = %g, DoSTot = %g\n",
           Nmol,Natom,dt,tmass,V,Delta,f,beta,DoS0,dos2,dostot);
           
    /* Now compute solid (2) and diffusive (3) components */
    fp = xvgropen(opt2fn("-dos",NFILE,fnm),"Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})",oenv);
    xvgr_legend(fp,asize(DoSlegend),DoSlegend,oenv);
    for(j=0; (j<teller/2); j++) 
    {
        dos[DOS_DIFF][j] = DoS0/(1+sqr(DoS0*M_PI*nu[j]/(6*f*Natom)));
        dos[DOS_SOLID][j] = dos[DOS][j]-dos[DOS_DIFF][j];
        fprintf(fp,"%10g  %10g  %10g  %10g\n",
                bRecip ? nu[j]/(1e-7*SPEED_OF_LIGHT) : nu[j],
                dos[DOS][j],dos[DOS_SOLID][j],dos[DOS_DIFF][j]);
    }
    fclose(fp);

    /* Finally analyze the results! */    
    wCdiff = 0.5;
    wSdiff = DoS0/(3*BOLTZ); /* Is this correct? */
/*    fp = xvgropen("test.xvg","Cp(v)","v","Cp",oenv);*/
    for(j=0; (j<teller/2); j++) 
    {
        dos[DOS_CP][j] = (dos[DOS_DIFF][j]*wCdiff + 
                          dos[DOS_SOLID][j]*wCsolid(nu[j],beta));
        dos[DOS_S][j]  = (dos[DOS_DIFF][j]*wSdiff + 
                          dos[DOS_SOLID][j]*wSsolid(nu[j],beta));
/*        fprintf(fp,"%10g  %10g  %10g\n",nu[j],
                dos[DOS_CP][j],wCsolid(nu[j],beta));
*/  }
    /*  fclose(fp);*/
    DiffCoeff = evaluate_integral(teller/2,tt,dos[VACF],NULL,teller/2,&stddev);
    DiffCoeff = 1000*DiffCoeff/3.0;
    printf("Diffusion coefficient from VACF %g 10^-5 cm^2/s\n",DiffCoeff);
    printf("Diffusion coefficient from DoS %g 10^-5 cm^2/s\n",
           1000*DoS0/(12*tmass*beta));
    
    cP = BOLTZ * evaluate_integral(teller/2,nu,dos[DOS_CP],NULL,teller/2,&stddev);
    printf("Heat capacity %g J/mol K\n",1000*cP/Nmol);
    
    S  = BOLTZ * evaluate_integral(teller/2,nu,dos[DOS_S],NULL,teller/2,&stddev);
    printf("Entropy %g J/mol K\n",1000*S/Nmol);
    
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
    thanx(stderr);
  
    return 0;
}
