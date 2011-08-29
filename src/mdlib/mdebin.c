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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <float.h>
#include "typedefs.h"
#include "string2.h"
#include "mdebin.h"
#include "smalloc.h"
#include "physics.h"
#include "enxio.h"
#include "vec.h"
#include "disre.h"
#include "main.h"
#include "network.h"
#include "names.h"
#include "orires.h"
#include "constr.h"
#include "mtop_util.h"
#include "xvgr.h"
#include "gmxfio.h"

#include "mdebin_bar.h"


static const char *conrmsd_nm[] = { "Constr. rmsd", "Constr.2 rmsd" };

static const char *boxs_nm[] = { "Box-X", "Box-Y", "Box-Z" };

static const char *tricl_boxs_nm[] = { 
    "Box-XX", "Box-YY", "Box-ZZ",
    "Box-YX", "Box-ZX", "Box-ZY" 
};

static const char *vol_nm[] = { "Volume" };

static const char *dens_nm[] = {"Density" };

static const char *pv_nm[] = {"pV" };

static const char *enthalpy_nm[] = {"Enthalpy" };

static const char *boxvel_nm[] = {
    "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
    "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
};

#define NBOXS asize(boxs_nm)
#define NTRICLBOXS asize(tricl_boxs_nm)


t_mdebin *init_mdebin(ener_file_t fp_ene,
                      const gmx_mtop_t *mtop,
                      const t_inputrec *ir,
                      FILE *fp_dhdl)
{
    const char *ener_nm[F_NRE];
    static const char *vir_nm[] = {
        "Vir-XX", "Vir-XY", "Vir-XZ",
        "Vir-YX", "Vir-YY", "Vir-YZ",
        "Vir-ZX", "Vir-ZY", "Vir-ZZ"
    };
    static const char *sv_nm[] = {
        "ShakeVir-XX", "ShakeVir-XY", "ShakeVir-XZ",
        "ShakeVir-YX", "ShakeVir-YY", "ShakeVir-YZ",
        "ShakeVir-ZX", "ShakeVir-ZY", "ShakeVir-ZZ"
    };
    static const char *fv_nm[] = {
        "ForceVir-XX", "ForceVir-XY", "ForceVir-XZ",
        "ForceVir-YX", "ForceVir-YY", "ForceVir-YZ",
        "ForceVir-ZX", "ForceVir-ZY", "ForceVir-ZZ"
    };
    static const char *pres_nm[] = {
        "Pres-XX","Pres-XY","Pres-XZ",
        "Pres-YX","Pres-YY","Pres-YZ",
        "Pres-ZX","Pres-ZY","Pres-ZZ"
    };
    static const char *surft_nm[] = {
        "#Surf*SurfTen"
    };
    static const char *mu_nm[] = {
        "Mu-X", "Mu-Y", "Mu-Z"
    };
    static const char *vcos_nm[] = {
        "2CosZ*Vel-X"
    };
    static const char *visc_nm[] = {
        "1/Viscosity"
    };
    static const char *baro_nm[] = {
        "Barostat"
    };

    char     **grpnms;
    const gmx_groups_t *groups;
    char     **gnm;
    char     buf[256];
    const char     *bufi;
    t_mdebin *md;
    int      i,j,ni,nj,n,nh,k,kk,ncon,nset;
    gmx_bool     bBHAM,bNoseHoover,b14;

    snew(md,1);

    if (EI_DYNAMICS(ir->eI))
    {
        md->delta_t = ir->delta_t;
    }
    else
    {
        md->delta_t = 0;
    }

    groups = &mtop->groups;

    bBHAM = (mtop->ffparams.functype[0] == F_BHAM);
    b14   = (gmx_mtop_ftype_count(mtop,F_LJ14) > 0 ||
             gmx_mtop_ftype_count(mtop,F_LJC14_Q) > 0);

    ncon = gmx_mtop_ftype_count(mtop,F_CONSTR);
    nset = gmx_mtop_ftype_count(mtop,F_SETTLE);
    md->bConstr    = (ncon > 0 || nset > 0);
    md->bConstrVir = FALSE;
    if (md->bConstr) {
        if (ncon > 0 && ir->eConstrAlg == econtLINCS) {
            if (ir->eI == eiSD2)
                md->nCrmsd = 2;
            else
                md->nCrmsd = 1;
        }
        md->bConstrVir = (getenv("GMX_CONSTRAINTVIR") != NULL);
    } else {
        md->nCrmsd = 0;
    }

    /* Energy monitoring */
    for(i=0;i<egNR;i++)
    {
        md->bEInd[i]=FALSE;
    }

#ifndef GMX_OPENMM
    for(i=0; i<F_NRE; i++)
    {
        md->bEner[i] = FALSE;
        if (i == F_LJ)
            md->bEner[i] = !bBHAM;
        else if (i == F_BHAM)
            md->bEner[i] = bBHAM;
        else if (i == F_EQM)
            md->bEner[i] = ir->bQMMM;
        else if (i == F_COUL_LR)
            md->bEner[i] = (ir->rcoulomb > ir->rlist);
        else if (i == F_LJ_LR)
            md->bEner[i] = (!bBHAM && ir->rvdw > ir->rlist);
        else if (i == F_BHAM_LR)
            md->bEner[i] = (bBHAM && ir->rvdw > ir->rlist);
        else if (i == F_RF_EXCL)
            md->bEner[i] = (EEL_RF(ir->coulombtype) && ir->coulombtype != eelRF_NEC);
        else if (i == F_COUL_RECIP)
            md->bEner[i] = EEL_FULL(ir->coulombtype);
        else if (i == F_LJ14)
            md->bEner[i] = b14;
        else if (i == F_COUL14)
            md->bEner[i] = b14;
        else if (i == F_LJC14_Q || i == F_LJC_PAIRS_NB)
            md->bEner[i] = FALSE;
        else if ((i == F_DVDL) || (i == F_DKDL))
            md->bEner[i] = (ir->efep != efepNO);
        else if (i == F_DHDL_CON)
            md->bEner[i] = (ir->efep != efepNO && md->bConstr);
        else if ((interaction_function[i].flags & IF_VSITE) ||
                 (i == F_CONSTR) || (i == F_CONSTRNC) || (i == F_SETTLE))
            md->bEner[i] = FALSE;
        else if ((i == F_COUL_SR) || (i == F_EPOT) || (i == F_PRES)  || (i==F_EQM))
            md->bEner[i] = TRUE;
        else if ((i == F_GBPOL) && ir->implicit_solvent==eisGBSA)
            md->bEner[i] = TRUE;
        else if ((i == F_NPSOLVATION) && ir->implicit_solvent==eisGBSA && (ir->sa_algorithm != esaNO))
            md->bEner[i] = TRUE;
        else if ((i == F_GB12) || (i == F_GB13) || (i == F_GB14))
            md->bEner[i] = FALSE;
        else if ((i == F_ETOT) || (i == F_EKIN) || (i == F_TEMP))
            md->bEner[i] = EI_DYNAMICS(ir->eI);
        else if (i==F_VTEMP) 
            md->bEner[i] =  (EI_DYNAMICS(ir->eI) && getenv("GMX_VIRIAL_TEMPERATURE"));
        else if (i == F_DISPCORR || i == F_PDISPCORR)
            md->bEner[i] = (ir->eDispCorr != edispcNO);
        else if (i == F_DISRESVIOL)
            md->bEner[i] = (gmx_mtop_ftype_count(mtop,F_DISRES) > 0);
        else if (i == F_ORIRESDEV)
            md->bEner[i] = (gmx_mtop_ftype_count(mtop,F_ORIRES) > 0);
        else if (i == F_CONNBONDS)
            md->bEner[i] = FALSE;
        else if (i == F_COM_PULL)
            md->bEner[i] = (ir->ePull == epullUMBRELLA || ir->ePull == epullCONST_F);
        else if (i == F_ECONSERVED)
            md->bEner[i] = ((ir->etc == etcNOSEHOOVER || ir->etc == etcVRESCALE) &&
                            (ir->epc == epcNO || ir->epc==epcMTTK));
        else
            md->bEner[i] = (gmx_mtop_ftype_count(mtop,i) > 0);
    }
#else
    /* OpenMM always produces only the following 4 energy terms */
    md->bEner[F_EPOT] = TRUE;
    md->bEner[F_EKIN] = TRUE;
    md->bEner[F_ETOT] = TRUE;
    md->bEner[F_TEMP] = TRUE;
#endif

    md->f_nre=0;
    for(i=0; i<F_NRE; i++)
    {
        if (md->bEner[i])
        {
            /* FIXME: The constness should not be cast away */
            /*ener_nm[f_nre]=(char *)interaction_function[i].longname;*/
            ener_nm[md->f_nre]=interaction_function[i].longname;
            md->f_nre++;
        }
    }

    md->epc = ir->epc;
    md->bDiagPres = !TRICLINIC(ir->ref_p);
    md->ref_p = (ir->ref_p[XX][XX]+ir->ref_p[YY][YY]+ir->ref_p[ZZ][ZZ])/DIM;
    md->bTricl = TRICLINIC(ir->compress) || TRICLINIC(ir->deform);
    md->bDynBox = DYNAMIC_BOX(*ir);
    md->etc = ir->etc;
    md->bNHC_trotter = IR_NVT_TROTTER(ir);
    md->bMTTK = IR_NPT_TROTTER(ir);

    md->ebin  = mk_ebin();
    /* Pass NULL for unit to let get_ebin_space determine the units
     * for interaction_function[i].longname
     */
    md->ie    = get_ebin_space(md->ebin,md->f_nre,ener_nm,NULL);
    if (md->nCrmsd)
    {
        /* This should be called directly after the call for md->ie,
         * such that md->iconrmsd follows directly in the list.
         */
        md->iconrmsd = get_ebin_space(md->ebin,md->nCrmsd,conrmsd_nm,"");
    }
    if (md->bDynBox)
    {
        md->ib    = get_ebin_space(md->ebin, 
                                   md->bTricl ? NTRICLBOXS : NBOXS, 
                                   md->bTricl ? tricl_boxs_nm : boxs_nm,
                                   unit_length);
        md->ivol  = get_ebin_space(md->ebin, 1, vol_nm,  unit_volume);
        md->idens = get_ebin_space(md->ebin, 1, dens_nm, unit_density_SI);
        if (md->bDiagPres)
        {
            md->ipv   = get_ebin_space(md->ebin, 1, pv_nm,   unit_energy);
            md->ienthalpy = get_ebin_space(md->ebin, 1, enthalpy_nm,   unit_energy);
        }
    }
    if (md->bConstrVir)
    {
        md->isvir = get_ebin_space(md->ebin,asize(sv_nm),sv_nm,unit_energy);
        md->ifvir = get_ebin_space(md->ebin,asize(fv_nm),fv_nm,unit_energy);
    }
    md->ivir   = get_ebin_space(md->ebin,asize(vir_nm),vir_nm,unit_energy);
    md->ipres  = get_ebin_space(md->ebin,asize(pres_nm),pres_nm,unit_pres_bar);
    md->isurft = get_ebin_space(md->ebin,asize(surft_nm),surft_nm,
                                unit_surft_bar);
    if (md->epc == epcPARRINELLORAHMAN || md->epc == epcMTTK)
    {
        md->ipc = get_ebin_space(md->ebin,md->bTricl ? 6 : 3,
                                 boxvel_nm,unit_vel);
    }
    md->imu    = get_ebin_space(md->ebin,asize(mu_nm),mu_nm,unit_dipole_D);
    if (ir->cos_accel != 0)
    {
        md->ivcos = get_ebin_space(md->ebin,asize(vcos_nm),vcos_nm,unit_vel);
        md->ivisc = get_ebin_space(md->ebin,asize(visc_nm),visc_nm,
                                   unit_invvisc_SI);
    }

    /* Energy monitoring */
    for(i=0;i<egNR;i++)
    {
        md->bEInd[i] = FALSE;
    }
    md->bEInd[egCOULSR] = TRUE;
    md->bEInd[egLJSR  ] = TRUE;

    if (ir->rcoulomb > ir->rlist)
    {
        md->bEInd[egCOULLR] = TRUE;
    }
    if (!bBHAM)
    {
        if (ir->rvdw > ir->rlist)
        {
            md->bEInd[egLJLR]   = TRUE;
        }
    }
    else
    {
        md->bEInd[egLJSR]   = FALSE;
        md->bEInd[egBHAMSR] = TRUE;
        if (ir->rvdw > ir->rlist)
        {
            md->bEInd[egBHAMLR]   = TRUE;
        }
    }
    if (b14)
    {
        md->bEInd[egLJ14] = TRUE;
        md->bEInd[egCOUL14] = TRUE;
    }
    md->nEc=0;
    for(i=0; (i<egNR); i++)
    {
        if (md->bEInd[i])
        {
            md->nEc++;
        }
    }

    n=groups->grps[egcENER].nr;
    md->nEg=n;
    md->nE=(n*(n+1))/2;
    snew(md->igrp,md->nE);
    if (md->nE > 1)
    {
        n=0;
        snew(gnm,md->nEc);
        for(k=0; (k<md->nEc); k++)
        {
            snew(gnm[k],STRLEN);
        }
        for(i=0; (i<groups->grps[egcENER].nr); i++)
        {
            ni=groups->grps[egcENER].nm_ind[i];
            for(j=i; (j<groups->grps[egcENER].nr); j++)
            {
                nj=groups->grps[egcENER].nm_ind[j];
                for(k=kk=0; (k<egNR); k++)
                {
                    if (md->bEInd[k])
                    {
                        sprintf(gnm[kk],"%s:%s-%s",egrp_nm[k],
                                *(groups->grpname[ni]),*(groups->grpname[nj]));
                        kk++;
                    }
                }
                md->igrp[n]=get_ebin_space(md->ebin,md->nEc,
                                           (const char **)gnm,unit_energy);
                n++;
            }
        }
        for(k=0; (k<md->nEc); k++)
        {
            sfree(gnm[k]);
        }
        sfree(gnm);

        if (n != md->nE)
        {
            gmx_incons("Number of energy terms wrong");
        }
    }

    md->nTC=groups->grps[egcTC].nr;
    md->nNHC = ir->opts.nhchainlength; /* shorthand for number of NH chains */ 
    if (md->bMTTK)
    {
        md->nTCP = 1;  /* assume only one possible coupling system for barostat 
                          for now */
    } 
    else 
    {
        md->nTCP = 0;
    }

    if (md->etc == etcNOSEHOOVER) {
        if (md->bNHC_trotter) { 
            md->mde_n = 2*md->nNHC*md->nTC;
        }
        else 
        {
            md->mde_n = 2*md->nTC;
        }
        if (md->epc == epcMTTK)
        {
            md->mdeb_n = 2*md->nNHC*md->nTCP;
        }
    } else { 
        md->mde_n = md->nTC;
        md->mdeb_n = 0;
    }

    snew(md->tmp_r,md->mde_n);
    snew(md->tmp_v,md->mde_n);
    snew(md->grpnms,md->mde_n);
    grpnms = md->grpnms;

    for(i=0; (i<md->nTC); i++)
    {
        ni=groups->grps[egcTC].nm_ind[i];
        sprintf(buf,"T-%s",*(groups->grpname[ni]));
        grpnms[i]=strdup(buf);
    }
    md->itemp=get_ebin_space(md->ebin,md->nTC,(const char **)grpnms,
                             unit_temp_K);

    bNoseHoover = (getenv("GMX_NOSEHOOVER_CHAINS") != NULL); /* whether to print Nose-Hoover chains */

    if (md->etc == etcNOSEHOOVER)
    {
        if (bNoseHoover) 
        {
            if (md->bNHC_trotter) 
            {
                for(i=0; (i<md->nTC); i++) 
                {
                    ni=groups->grps[egcTC].nm_ind[i];
                    bufi = *(groups->grpname[ni]);
                    for(j=0; (j<md->nNHC); j++) 
                    {
                        sprintf(buf,"Xi-%d-%s",j,bufi);
                        grpnms[2*(i*md->nNHC+j)]=strdup(buf);
                        sprintf(buf,"vXi-%d-%s",j,bufi);
                        grpnms[2*(i*md->nNHC+j)+1]=strdup(buf);
                    }
                }
                md->itc=get_ebin_space(md->ebin,md->mde_n,
                                       (const char **)grpnms,unit_invtime);
                if (md->bMTTK) 
                {
                    for(i=0; (i<md->nTCP); i++) 
                    {
                        bufi = baro_nm[0];  /* All barostat DOF's together for now. */
                        for(j=0; (j<md->nNHC); j++) 
                        {
                            sprintf(buf,"Xi-%d-%s",j,bufi);
                            grpnms[2*(i*md->nNHC+j)]=strdup(buf);
                            sprintf(buf,"vXi-%d-%s",j,bufi);
                            grpnms[2*(i*md->nNHC+j)+1]=strdup(buf);
                        }
                    }
                    md->itcb=get_ebin_space(md->ebin,md->mdeb_n,
                                            (const char **)grpnms,unit_invtime);
                }
            } 
            else
            {
                for(i=0; (i<md->nTC); i++) 
                {
                    ni=groups->grps[egcTC].nm_ind[i];
                    bufi = *(groups->grpname[ni]);
                    sprintf(buf,"Xi-%s",bufi);
                    grpnms[2*i]=strdup(buf);
                    sprintf(buf,"vXi-%s",bufi);
                    grpnms[2*i+1]=strdup(buf);
                }
                md->itc=get_ebin_space(md->ebin,md->mde_n,
                                       (const char **)grpnms,unit_invtime);
            }
        }
    }
    else if (md->etc == etcBERENDSEN || md->etc == etcYES || 
             md->etc == etcVRESCALE)
    {
        for(i=0; (i<md->nTC); i++)
        {
            ni=groups->grps[egcTC].nm_ind[i];
            sprintf(buf,"Lamb-%s",*(groups->grpname[ni]));
            grpnms[i]=strdup(buf);
        }
        md->itc=get_ebin_space(md->ebin,md->mde_n,(const char **)grpnms,"");
    }

    sfree(grpnms);


    md->nU=groups->grps[egcACC].nr;
    if (md->nU > 1)
    {
        snew(grpnms,3*md->nU);
        for(i=0; (i<md->nU); i++)
        {
            ni=groups->grps[egcACC].nm_ind[i];
            sprintf(buf,"Ux-%s",*(groups->grpname[ni]));
            grpnms[3*i+XX]=strdup(buf);
            sprintf(buf,"Uy-%s",*(groups->grpname[ni]));
            grpnms[3*i+YY]=strdup(buf);
            sprintf(buf,"Uz-%s",*(groups->grpname[ni]));
            grpnms[3*i+ZZ]=strdup(buf);
        }
        md->iu=get_ebin_space(md->ebin,3*md->nU,(const char **)grpnms,unit_vel);
        sfree(grpnms);
    }

    if ( fp_ene )
    {
        do_enxnms(fp_ene,&md->ebin->nener,&md->ebin->enm);
    }

    md->print_grpnms=NULL;

    /* check whether we're going to write dh histograms */
    md->dhc=NULL; 
    if (ir->separate_dhdl_file == sepdhdlfileNO )
    {
        int i;
        snew(md->dhc, 1);

        mde_delta_h_coll_init(md->dhc, ir);
        md->fp_dhdl = NULL;
    }
    else
    {
        md->fp_dhdl = fp_dhdl;
    }
    md->dhdl_derivatives = (ir->dhdl_derivatives==dhdlderivativesYES);
    return md;
}

FILE *open_dhdl(const char *filename,const t_inputrec *ir,
                const output_env_t oenv)
{
    FILE *fp;
    const char *dhdl="dH/d\\lambda",*deltag="\\DeltaH",*lambda="\\lambda";
    char title[STRLEN],label_x[STRLEN],label_y[STRLEN];
    char **setname;
    char buf[STRLEN];

    sprintf(label_x,"%s (%s)","Time",unit_time);
    if (ir->n_flambda == 0)
    {
        sprintf(title,"%s",dhdl);
        sprintf(label_y,"%s (%s %s)",
                dhdl,unit_energy,"[\\lambda]\\S-1\\N");
    }
    else
    {
        sprintf(title,"%s, %s",dhdl,deltag);
        sprintf(label_y,"(%s)",unit_energy);
    }
    fp = gmx_fio_fopen(filename,"w+");
    xvgr_header(fp,title,label_x,label_y,exvggtXNY,oenv);

    if (ir->delta_lambda == 0)
    {
        sprintf(buf,"T = %g (K), %s = %g",
                ir->opts.ref_t[0],lambda,ir->init_lambda);
    }
    else
    {
        sprintf(buf,"T = %g (K)",
                ir->opts.ref_t[0]);
    }
    xvgr_subtitle(fp,buf,oenv);

    if (ir->n_flambda > 0)
    {
        int nsets,s,nsi=0;
        /* g_bar has to determine the lambda values used in this simulation
         * from this xvg legend.  */
        nsets = ( (ir->dhdl_derivatives==dhdlderivativesYES) ? 1 : 0) + 
                  ir->n_flambda;
        snew(setname,nsets);
        if (ir->dhdl_derivatives == dhdlderivativesYES)
        {
            sprintf(buf,"%s %s %g",dhdl,lambda,ir->init_lambda);
            setname[nsi++] = gmx_strdup(buf);
        }
        for(s=0; s<ir->n_flambda; s++)
        {
            sprintf(buf,"%s %s %g",deltag,lambda,ir->flambda[s]);
            setname[nsi++] = gmx_strdup(buf);
        }
        xvgr_legend(fp,nsets,(const char**)setname,oenv);

        for(s=0; s<nsets; s++)
        {
            sfree(setname[s]);
        }
        sfree(setname);
    }

    return fp;
}

static void copy_energy(t_mdebin *md, real e[],real ecpy[])
{
    int i,j;

    for(i=j=0; (i<F_NRE); i++)
        if (md->bEner[i])
            ecpy[j++] = e[i];
    if (j != md->f_nre) 
        gmx_incons("Number of energy terms wrong");
}

void upd_mdebin(t_mdebin *md, gmx_bool write_dhdl,
                gmx_bool bSum,
                double time,
                real tmass,
                gmx_enerdata_t *enerd,
                t_state *state,
                matrix  box,
                tensor svir,
                tensor fvir,
                tensor vir,
                tensor pres,
                gmx_ekindata_t *ekind,
                rvec mu_tot,
                gmx_constr_t constr)
{
    int    i,j,k,kk,m,n,gid;
    real   crmsd[2],tmp6[6];
    real   bs[NTRICLBOXS],vol,dens,pv,enthalpy;
    real   eee[egNR];
    real   ecopy[F_NRE];
    real   tmp;
    gmx_bool   bNoseHoover;

    /* Do NOT use the box in the state variable, but the separate box provided
     * as an argument. This is because we sometimes need to write the box from
     * the last timestep to match the trajectory frames.
     */
    copy_energy(md, enerd->term,ecopy);
    add_ebin(md->ebin,md->ie,md->f_nre,ecopy,bSum);
    if (md->nCrmsd)
    {
        crmsd[0] = constr_rmsd(constr,FALSE);
        if (md->nCrmsd > 1)
        {
            crmsd[1] = constr_rmsd(constr,TRUE);
        }
        add_ebin(md->ebin,md->iconrmsd,md->nCrmsd,crmsd,FALSE);
    }
    if (md->bDynBox)
    {
        int nboxs;
        if(md->bTricl)
        {
            bs[0] = box[XX][XX];
            bs[1] = box[YY][YY];
            bs[2] = box[ZZ][ZZ];
            bs[3] = box[YY][XX];
            bs[4] = box[ZZ][XX];
            bs[5] = box[ZZ][YY];
            nboxs=NTRICLBOXS;
        }
        else
        {
            bs[0] = box[XX][XX];
            bs[1] = box[YY][YY];
            bs[2] = box[ZZ][ZZ];
            nboxs=NBOXS;
        }
        vol  = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
        dens = (tmass*AMU)/(vol*NANO*NANO*NANO);

        add_ebin(md->ebin,md->ib   ,nboxs,bs   ,bSum);
        add_ebin(md->ebin,md->ivol ,1    ,&vol ,bSum);
        add_ebin(md->ebin,md->idens,1    ,&dens,bSum);

        if (md->bDiagPres)
        {
            /* This is pV (in kJ/mol).  The pressure is the reference pressure,
               not the instantaneous pressure */  
            pv = vol*md->ref_p/PRESFAC;

            add_ebin(md->ebin,md->ipv  ,1    ,&pv  ,bSum);
            enthalpy = pv + enerd->term[F_ETOT];
            add_ebin(md->ebin,md->ienthalpy  ,1    ,&enthalpy  ,bSum);
        }
    }
    if (md->bConstrVir)
    {
        add_ebin(md->ebin,md->isvir,9,svir[0],bSum);
        add_ebin(md->ebin,md->ifvir,9,fvir[0],bSum);
    }
    add_ebin(md->ebin,md->ivir,9,vir[0],bSum);
    add_ebin(md->ebin,md->ipres,9,pres[0],bSum);
    tmp = (pres[ZZ][ZZ]-(pres[XX][XX]+pres[YY][YY])*0.5)*box[ZZ][ZZ];
    add_ebin(md->ebin,md->isurft,1,&tmp,bSum);
    if (md->epc == epcPARRINELLORAHMAN || md->epc == epcMTTK)
    {
        tmp6[0] = state->boxv[XX][XX];
        tmp6[1] = state->boxv[YY][YY];
        tmp6[2] = state->boxv[ZZ][ZZ];
        tmp6[3] = state->boxv[YY][XX];
        tmp6[4] = state->boxv[ZZ][XX];
        tmp6[5] = state->boxv[ZZ][YY];
        add_ebin(md->ebin,md->ipc,md->bTricl ? 6 : 3,tmp6,bSum);
    }
    add_ebin(md->ebin,md->imu,3,mu_tot,bSum);
    if (ekind && ekind->cosacc.cos_accel != 0)
    {
        vol  = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
        dens = (tmass*AMU)/(vol*NANO*NANO*NANO);
        add_ebin(md->ebin,md->ivcos,1,&(ekind->cosacc.vcos),bSum);
        /* 1/viscosity, unit 1/(kg m^-1 s^-1) */
        tmp = 1/(ekind->cosacc.cos_accel/(ekind->cosacc.vcos*PICO)
                 *dens*vol*sqr(box[ZZ][ZZ]*NANO/(2*M_PI)));
        add_ebin(md->ebin,md->ivisc,1,&tmp,bSum);    
    }
    if (md->nE > 1)
    {
        n=0;
        for(i=0; (i<md->nEg); i++)
        {
            for(j=i; (j<md->nEg); j++)
            {
                gid=GID(i,j,md->nEg);
                for(k=kk=0; (k<egNR); k++)
                {
                    if (md->bEInd[k])
                    {
                        eee[kk++] = enerd->grpp.ener[k][gid];
                    }
                }
                add_ebin(md->ebin,md->igrp[n],md->nEc,eee,bSum);
                n++;
            }
        }
    }

    if (ekind)
    {
        for(i=0; (i<md->nTC); i++)
        {
            md->tmp_r[i] = ekind->tcstat[i].T;
        }
        add_ebin(md->ebin,md->itemp,md->nTC,md->tmp_r,bSum);

        /* whether to print Nose-Hoover chains: */
        bNoseHoover = (getenv("GMX_NOSEHOOVER_CHAINS") != NULL); 

        if (md->etc == etcNOSEHOOVER)
        {
            if (bNoseHoover) 
            {
                if (md->bNHC_trotter)
                {
                    for(i=0; (i<md->nTC); i++) 
                    {
                        for (j=0;j<md->nNHC;j++) 
                        {
                            k = i*md->nNHC+j;
                            md->tmp_r[2*k] = state->nosehoover_xi[k];
                            md->tmp_r[2*k+1] = state->nosehoover_vxi[k];
                        }
                    }
                    add_ebin(md->ebin,md->itc,md->mde_n,md->tmp_r,bSum);      

                    if (md->bMTTK) {
                        for(i=0; (i<md->nTCP); i++) 
                        {
                            for (j=0;j<md->nNHC;j++) 
                            {
                                k = i*md->nNHC+j;
                                md->tmp_r[2*k] = state->nhpres_xi[k];
                                md->tmp_r[2*k+1] = state->nhpres_vxi[k];
                            }
                        }
                        add_ebin(md->ebin,md->itcb,md->mdeb_n,md->tmp_r,bSum);      
                    }

                } 
                else 
                {
                    for(i=0; (i<md->nTC); i++)
                    {
                        md->tmp_r[2*i] = state->nosehoover_xi[i];
                        md->tmp_r[2*i+1] = state->nosehoover_vxi[i];
                    }
                    add_ebin(md->ebin,md->itc,md->mde_n,md->tmp_r,bSum);
                }
            }
        }
        else if (md->etc == etcBERENDSEN || md->etc == etcYES || 
                 md->etc == etcVRESCALE)
        {
            for(i=0; (i<md->nTC); i++)
            {
                md->tmp_r[i] = ekind->tcstat[i].lambda;
            }
            add_ebin(md->ebin,md->itc,md->nTC,md->tmp_r,bSum);
        }
    }

    if (ekind && md->nU > 1)
    {
        for(i=0; (i<md->nU); i++)
        {
            copy_rvec(ekind->grpstat[i].u,md->tmp_v[i]);
        }
        add_ebin(md->ebin,md->iu,3*md->nU,md->tmp_v[0],bSum);
    }

    ebin_increase_count(md->ebin,bSum);

    /* BAR + thermodynamic integration values */
    if (write_dhdl)
    {
        if (md->fp_dhdl)
        {
            fprintf(md->fp_dhdl,"%.4f", time);

            if (md->dhdl_derivatives)
            {
                fprintf(md->fp_dhdl," %g", enerd->term[F_DVDL]+ 
                                           enerd->term[F_DKDL]+
                                           enerd->term[F_DHDL_CON]);
            }
            for(i=1; i<enerd->n_lambda; i++)
            {
                fprintf(md->fp_dhdl," %g",
                        enerd->enerpart_lambda[i]-enerd->enerpart_lambda[0]);
            }
            fprintf(md->fp_dhdl,"\n");
        }
        /* and the binary BAR output */
        if (md->dhc)
        {
            mde_delta_h_coll_add_dh(md->dhc, 
                                    enerd->term[F_DVDL]+ enerd->term[F_DKDL]+
                                    enerd->term[F_DHDL_CON],
                                    enerd->enerpart_lambda, time, 
                                    state->lambda);
        }
    }
}

void upd_mdebin_step(t_mdebin *md)
{
    ebin_increase_count(md->ebin,FALSE); 
}

static void npr(FILE *log,int n,char c)
{
    for(; (n>0); n--) fprintf(log,"%c",c);
}

static void pprint(FILE *log,const char *s,t_mdebin *md)
{
    char CHAR='#';
    int  slen;
    char buf1[22],buf2[22];

    slen = strlen(s);
    fprintf(log,"\t<======  ");
    npr(log,slen,CHAR);
    fprintf(log,"  ==>\n");
    fprintf(log,"\t<====  %s  ====>\n",s);
    fprintf(log,"\t<==  ");
    npr(log,slen,CHAR);
    fprintf(log,"  ======>\n\n");

    fprintf(log,"\tStatistics over %s steps using %s frames\n",
            gmx_step_str(md->ebin->nsteps_sim,buf1),
            gmx_step_str(md->ebin->nsum_sim,buf2));
    fprintf(log,"\n");
}

void print_ebin_header(FILE *log,gmx_large_int_t steps,double time,real lamb)
{
    char buf[22];

    fprintf(log,"   %12s   %12s   %12s\n"
            "   %12s   %12.5f   %12.5f\n\n",
            "Step","Time","Lambda",gmx_step_str(steps,buf),time,lamb);
}

void print_ebin(ener_file_t fp_ene,gmx_bool bEne,gmx_bool bDR,gmx_bool bOR,
                FILE *log,
                gmx_large_int_t step,double time,
                int mode,gmx_bool bCompact,
                t_mdebin *md,t_fcdata *fcd,
                gmx_groups_t *groups,t_grpopts *opts)
{
    /*static char **grpnms=NULL;*/
    char        buf[246];
    int         i,j,n,ni,nj,ndr,nor,b;
    int         ndisre=0;
    real        *disre_rm3tav, *disre_rt;

    /* these are for the old-style blocks (1 subblock, only reals), because
       there can be only one per ID for these */
    int         nr[enxNR];
    int         id[enxNR];
    real        *block[enxNR];

    /* temporary arrays for the lambda values to write out */
    double      enxlambda_data[2]; 

    t_enxframe  fr;

    switch (mode)
    {
        case eprNORMAL:
            init_enxframe(&fr);
            fr.t            = time;
            fr.step         = step;
            fr.nsteps       = md->ebin->nsteps;
            fr.dt           = md->delta_t;
            fr.nsum         = md->ebin->nsum;
            fr.nre          = (bEne) ? md->ebin->nener : 0;
            fr.ener         = md->ebin->e;
            ndisre          = bDR ? fcd->disres.npair : 0;
            disre_rm3tav    = fcd->disres.rm3tav;
            disre_rt        = fcd->disres.rt;
            /* Optional additional old-style (real-only) blocks. */
            for(i=0; i<enxNR; i++)
            {
                nr[i] = 0;
            }
            if (fcd->orires.nr > 0 && bOR)
            {
                diagonalize_orires_tensors(&(fcd->orires));
                nr[enxOR]     = fcd->orires.nr;
                block[enxOR]  = fcd->orires.otav;
                id[enxOR]     = enxOR;
                nr[enxORI]    = (fcd->orires.oinsl != fcd->orires.otav) ? 
                          fcd->orires.nr : 0;
                block[enxORI] = fcd->orires.oinsl;
                id[enxORI]    = enxORI;
                nr[enxORT]    = fcd->orires.nex*12;
                block[enxORT] = fcd->orires.eig;
                id[enxORT]    = enxORT;
            }        

            /* whether we are going to wrte anything out: */
            if (fr.nre || ndisre || nr[enxOR] || nr[enxORI])
            {

                /* the old-style blocks go first */
                fr.nblock = 0;
                for(i=0; i<enxNR; i++)
                {
                    if (nr[i] > 0)
                    {
                        fr.nblock = i + 1;
                    }
                }
                add_blocks_enxframe(&fr, fr.nblock);
                for(b=0;b<fr.nblock;b++)
                {
                    add_subblocks_enxblock(&(fr.block[b]), 1);
                    fr.block[b].id=id[b]; 
                    fr.block[b].sub[0].nr = nr[b];
#ifndef GMX_DOUBLE
                    fr.block[b].sub[0].type = xdr_datatype_float;
                    fr.block[b].sub[0].fval = block[b];
#else
                    fr.block[b].sub[0].type = xdr_datatype_double;
                    fr.block[b].sub[0].dval = block[b];
#endif
                }

                /* check for disre block & fill it. */
                if (ndisre>0)
                {
                    int db = fr.nblock;
                    fr.nblock+=1;
                    add_blocks_enxframe(&fr, fr.nblock);

                    add_subblocks_enxblock(&(fr.block[db]), 2);
                    fr.block[db].id=enxDISRE;
                    fr.block[db].sub[0].nr=ndisre;
                    fr.block[db].sub[1].nr=ndisre;
#ifndef GMX_DOUBLE
                    fr.block[db].sub[0].type=xdr_datatype_float;
                    fr.block[db].sub[1].type=xdr_datatype_float;
                    fr.block[db].sub[0].fval=disre_rt;
                    fr.block[db].sub[1].fval=disre_rm3tav;
#else
                    fr.block[db].sub[0].type=xdr_datatype_double;
                    fr.block[db].sub[1].type=xdr_datatype_double;
                    fr.block[db].sub[0].dval=disre_rt;
                    fr.block[db].sub[1].dval=disre_rm3tav;
#endif
                }
                /* here we can put new-style blocks */

                /* Free energy perturbation blocks */
                if (md->dhc)
                {
                    mde_delta_h_coll_handle_block(md->dhc, &fr, fr.nblock);
                }

                /* do the actual I/O */
                do_enx(fp_ene,&fr);
                gmx_fio_check_file_position(enx_file_pointer(fp_ene));
                if (fr.nre)
                {
                    /* We have stored the sums, so reset the sum history */
                    reset_ebin_sums(md->ebin);
                }

                /* we can now free & reset the data in the blocks */
                if (md->dhc)
                    mde_delta_h_coll_reset(md->dhc);
            }
            free_enxframe(&fr);
            break;
        case eprAVER:
            if (log)
            {
                pprint(log,"A V E R A G E S",md);
            }
            break;
        case eprRMS:
            if (log)
            {
                pprint(log,"R M S - F L U C T U A T I O N S",md);
            }
            break;
        default:
            gmx_fatal(FARGS,"Invalid print mode (%d)",mode);
    }

    if (log)
    {
        for(i=0;i<opts->ngtc;i++)
        {
            if(opts->annealing[i]!=eannNO)
            {
                fprintf(log,"Current ref_t for group %s: %8.1f\n",
                        *(groups->grpname[groups->grps[egcTC].nm_ind[i]]),
                        opts->ref_t[i]);
            }
        }
        if (mode==eprNORMAL && fcd->orires.nr>0)
        {
            print_orires_log(log,&(fcd->orires));
        }
        fprintf(log,"   Energies (%s)\n",unit_energy);
        pr_ebin(log,md->ebin,md->ie,md->f_nre+md->nCrmsd,5,mode,TRUE);  
        fprintf(log,"\n");

        if (!bCompact)
        {
            if (md->bDynBox)
            {
                pr_ebin(log,md->ebin,md->ib, md->bTricl ? NTRICLBOXS : NBOXS,5,
                        mode,TRUE);      
                fprintf(log,"\n");
            }
            if (md->bConstrVir)
            {
                fprintf(log,"   Constraint Virial (%s)\n",unit_energy);
                pr_ebin(log,md->ebin,md->isvir,9,3,mode,FALSE);  
                fprintf(log,"\n");
                fprintf(log,"   Force Virial (%s)\n",unit_energy);
                pr_ebin(log,md->ebin,md->ifvir,9,3,mode,FALSE);  
                fprintf(log,"\n");
            }
            fprintf(log,"   Total Virial (%s)\n",unit_energy);
            pr_ebin(log,md->ebin,md->ivir,9,3,mode,FALSE);   
            fprintf(log,"\n");
            fprintf(log,"   Pressure (%s)\n",unit_pres_bar);
            pr_ebin(log,md->ebin,md->ipres,9,3,mode,FALSE);  
            fprintf(log,"\n");
            fprintf(log,"   Total Dipole (%s)\n",unit_dipole_D);
            pr_ebin(log,md->ebin,md->imu,3,3,mode,FALSE);    
            fprintf(log,"\n");

            if (md->nE > 1)
            {
                if (md->print_grpnms==NULL)
                {
                    snew(md->print_grpnms,md->nE);
                    n=0;
                    for(i=0; (i<md->nEg); i++)
                    {
                        ni=groups->grps[egcENER].nm_ind[i];
                        for(j=i; (j<md->nEg); j++)
                        {
                            nj=groups->grps[egcENER].nm_ind[j];
                            sprintf(buf,"%s-%s",*(groups->grpname[ni]),
                                    *(groups->grpname[nj]));
                            md->print_grpnms[n++]=strdup(buf);
                        }
                    }
                }
                sprintf(buf,"Epot (%s)",unit_energy);
                fprintf(log,"%15s   ",buf);
                for(i=0; (i<egNR); i++)
                {
                    if (md->bEInd[i])
                    {
                        fprintf(log,"%12s   ",egrp_nm[i]);
                    }
                }
                fprintf(log,"\n");
                for(i=0; (i<md->nE); i++)
                {
                    fprintf(log,"%15s",md->print_grpnms[i]);
                    pr_ebin(log,md->ebin,md->igrp[i],md->nEc,md->nEc,mode,
                            FALSE);
                }
                fprintf(log,"\n");
            }
            if (md->nTC > 1)
            {
                pr_ebin(log,md->ebin,md->itemp,md->nTC,4,mode,TRUE);
                fprintf(log,"\n");
            }
            if (md->nU > 1)
            {
                fprintf(log,"%15s   %12s   %12s   %12s\n",
                        "Group","Ux","Uy","Uz");
                for(i=0; (i<md->nU); i++)
                {
                    ni=groups->grps[egcACC].nm_ind[i];
                    fprintf(log,"%15s",*groups->grpname[ni]);
                    pr_ebin(log,md->ebin,md->iu+3*i,3,3,mode,FALSE);
                }
                fprintf(log,"\n");
            }
        }
    }

}

void update_energyhistory(energyhistory_t * enerhist,t_mdebin * mdebin)
{
    int i;

    enerhist->nsteps     = mdebin->ebin->nsteps;
    enerhist->nsum       = mdebin->ebin->nsum;
    enerhist->nsteps_sim = mdebin->ebin->nsteps_sim;
    enerhist->nsum_sim   = mdebin->ebin->nsum_sim;
    enerhist->nener      = mdebin->ebin->nener;

    if (mdebin->ebin->nsum > 0)
    {
        /* Check if we need to allocate first */
        if(enerhist->ener_ave == NULL)
        {
            snew(enerhist->ener_ave,enerhist->nener);
            snew(enerhist->ener_sum,enerhist->nener);
        }

        for(i=0;i<enerhist->nener;i++)
        {
            enerhist->ener_ave[i] = mdebin->ebin->e[i].eav;
            enerhist->ener_sum[i] = mdebin->ebin->e[i].esum;
        }
    }

    if (mdebin->ebin->nsum_sim > 0)
    {
        /* Check if we need to allocate first */
        if(enerhist->ener_sum_sim == NULL)
        {
            snew(enerhist->ener_sum_sim,enerhist->nener);
        }

        for(i=0;i<enerhist->nener;i++)
        {
            enerhist->ener_sum_sim[i] = mdebin->ebin->e_sim[i].esum;
        }
    }
    if (mdebin->dhc)
    {
        mde_delta_h_coll_update_energyhistory(mdebin->dhc, enerhist);
    }
}

void restore_energyhistory_from_state(t_mdebin * mdebin,
                                      energyhistory_t * enerhist)
{
    int i;

    if ((enerhist->nsum > 0 || enerhist->nsum_sim > 0) &&
        mdebin->ebin->nener != enerhist->nener)
    {
        gmx_fatal(FARGS,"Mismatch between number of energies in run input (%d) and checkpoint file (%d).",
                  mdebin->ebin->nener,enerhist->nener);
    }

    mdebin->ebin->nsteps     = enerhist->nsteps;
    mdebin->ebin->nsum       = enerhist->nsum;
    mdebin->ebin->nsteps_sim = enerhist->nsteps_sim;
    mdebin->ebin->nsum_sim   = enerhist->nsum_sim;

    for(i=0; i<mdebin->ebin->nener; i++)
    {
        mdebin->ebin->e[i].eav  =
                  (enerhist->nsum > 0 ? enerhist->ener_ave[i] : 0);
        mdebin->ebin->e[i].esum =
                  (enerhist->nsum > 0 ? enerhist->ener_sum[i] : 0);
        mdebin->ebin->e_sim[i].esum =
                  (enerhist->nsum_sim > 0 ? enerhist->ener_sum_sim[i] : 0);
    }
    if (mdebin->dhc)
    {         
        mde_delta_h_coll_restore_energyhistory(mdebin->dhc, enerhist);
    }
}
