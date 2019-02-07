/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
/*! \internal \file
 * \brief Defines code that writes energy-like quantities.
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "energyoutput.h"

#include <cfloat>
#include <cstdlib>
#include <cstring>

#include <string>

#include "gromacs/awh/awh.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/mdebin_bar.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

//! Labels for energy file quantities
//! \{
static const char                 *conrmsd_nm[] = { "Constr. rmsd", "Constr.2 rmsd" };

static std::array<const char *, 3> boxs_nm = { "Box-X", "Box-Y", "Box-Z" };

static std::array<const char *, 6> tricl_boxs_nm = {
    "Box-XX", "Box-YY", "Box-ZZ",
    "Box-YX", "Box-ZX", "Box-ZY"
};

static const char                 *vol_nm[] = { "Volume" };

static const char                 *dens_nm[] = {"Density" };

static const char                 *pv_nm[] = {"pV" };

static const char                 *enthalpy_nm[] = {"Enthalpy" };

static std::array<const char *, 6> boxvel_nm = {
    "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
    "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY"
};

const char *egrp_nm[egNR+1] = {
    "Coul-SR", "LJ-SR", "Buck-SR",
    "Coul-14", "LJ-14", nullptr
};
//! \}

struct t_mde_delta_h_coll;

namespace gmx
{

/*! \brief Implementation class
 *
 * This is the collection of energy averages collected during mdrun, and to
 * be written out to the .edr file. It also writes other legacy output,
 * e.g. to dhdl.xvg file.
 *
 * \todo Use more std containers.
 * \todo Name member variables with underscore suffixes.
 * \todo Remove GMX_CONSTRAINTVIR
 * \todo Write free-energy output also to energy file (after adding more tests)
 */
class EnergyOutput::Impl
{
    public:
        Impl();
        ~Impl();
        void prepare(const gmx_mtop_t *mtop,
                     const t_inputrec *ir,
                     FILE             *fp_dhdl,
                     bool              isRerun);

        void addDataAtEnergyStep(bool                    bDoDHDL,
                                 bool                    bSum,
                                 double                  time,
                                 real                    tmass,
                                 const gmx_enerdata_t   &enerd,
                                 t_state                *state,
                                 t_lambda               *fep,
                                 t_expanded             *expand,
                                 matrix                  box,
                                 tensor                  svir,
                                 tensor                  fvir,
                                 tensor                  vir,
                                 tensor                  pres,
                                 const gmx_ekindata_t   *ekind,
                                 rvec                    mu_tot,
                                 const gmx::Constraints *constr);
        void recordNonEnergyStep();
        void printStepToEnergyFile(ener_file *fp_ene, bool bEne, bool bDR, bool bOR,
                                   FILE *log,
                                   int64_t step, double time,
                                   int mode,
                                   t_fcdata *fcd,
                                   gmx_groups_t *groups, t_grpopts *opts,
                                   gmx::Awh *awh);
        void fillEnergyHistory(energyhistory_t *enerhist) const;
        void restoreFromEnergyHistory(const energyhistory_t &enerhist);

        double  delta_t        = 0;
        t_ebin *ebin           = nullptr;
        //! Indices to use when storing data to ebin
        //! /{
        int     ie             = 0;
        int     iconrmsd       = 0;
        int     ib             = 0;
        int     ivol           = 0;
        int     idens          = 0;
        int     ipv            = 0;
        int     ienthalpy      = 0;
        int     isvir          = 0;
        int     ifvir          = 0;
        int     ipres          = 0;
        int     ivir           = 0;
        int     isurft         = 0;
        int     ipc            = 0;
        int     itemp          = 0;
        int     itc            = 0;
        int     itcb           = 0;
        int     iu             = 0;
        int     imu            = 0;
        int     ivcos          = 0;
        int     ivisc          = 0;
        int     nE             = 0;
        int     nEg            = 0;
        int     nEc            = 0;
        int     nTC            = 0;
        int     nTCP           = 0;
        int     nU             = 0;
        int     nNHC           = 0;
        int    *igrp           = nullptr;
        //! /}

        int     mde_n          = 0;
        int     mdeb_n         = 0;
        real   *tmp_r          = nullptr;
        rvec   *tmp_v          = nullptr;
        bool    bConstr        = false;
        bool    bConstrVir     = false;
        bool    bTricl         = false;
        bool    bDynBox        = false;
        bool    bNHC_trotter   = false;
        bool    bPrintNHChains = false;
        bool    bMTTK          = false;
        bool    bMu            = false; /* true if dipole is calculated */
        bool    bDiagPres      = false;
        bool    bPres          = false;
        int     f_nre          = 0;
        int     epc            = 0;
        real    ref_p          = 0;
        int     etc            = 0;
        int     nCrmsd         = 0;
        bool    bEner[F_NRE]   = {false};
        bool    bEInd[egNR]    = {false};
        char  **print_grpnms   = nullptr;

        FILE   *fp_dhdl      = nullptr; /* the dhdl.xvg output file */
        double *dE           = nullptr; /* energy components for dhdl.xvg output */


        // Extra space for uncrustify
        //! The delta U components (raw data + histogram)
        t_mde_delta_h_coll *dhc          = nullptr;
        real               *temperatures = nullptr;
};


EnergyOutput::Impl::Impl()
    : fp_dhdl(nullptr)
{
    ebin = mk_ebin();
}

void EnergyOutput::Impl::prepare(const gmx_mtop_t *mtop,
                                 const t_inputrec *ir,
                                 FILE             *fp_dhdl_arg,
                                 bool              isRerun)
{
    const char         *ener_nm[F_NRE];
    static const char  *vir_nm[] = {
        "Vir-XX", "Vir-XY", "Vir-XZ",
        "Vir-YX", "Vir-YY", "Vir-YZ",
        "Vir-ZX", "Vir-ZY", "Vir-ZZ"
    };
    static const char  *sv_nm[] = {
        "ShakeVir-XX", "ShakeVir-XY", "ShakeVir-XZ",
        "ShakeVir-YX", "ShakeVir-YY", "ShakeVir-YZ",
        "ShakeVir-ZX", "ShakeVir-ZY", "ShakeVir-ZZ"
    };
    static const char  *fv_nm[] = {
        "ForceVir-XX", "ForceVir-XY", "ForceVir-XZ",
        "ForceVir-YX", "ForceVir-YY", "ForceVir-YZ",
        "ForceVir-ZX", "ForceVir-ZY", "ForceVir-ZZ"
    };
    static const char  *pres_nm[] = {
        "Pres-XX", "Pres-XY", "Pres-XZ",
        "Pres-YX", "Pres-YY", "Pres-YZ",
        "Pres-ZX", "Pres-ZY", "Pres-ZZ"
    };
    static const char  *surft_nm[] = {
        "#Surf*SurfTen"
    };
    static const char  *mu_nm[] = {
        "Mu-X", "Mu-Y", "Mu-Z"
    };
    static const char  *vcos_nm[] = {
        "2CosZ*Vel-X"
    };
    static const char  *visc_nm[] = {
        "1/Viscosity"
    };
    static const char  *baro_nm[] = {
        "Barostat"
    };

    const gmx_groups_t *groups;
    char              **gnm;
    char                buf[256];
    const char         *bufi;
    int                 i, j, ni, nj, n, k, kk, ncon, nset;
    bool                bBHAM, b14;

    if (EI_DYNAMICS(ir->eI))
    {
        delta_t = ir->delta_t;
    }
    else
    {
        delta_t = 0;
    }

    groups = &mtop->groups;

    bBHAM = (mtop->ffparams.numTypes() > 0) && (mtop->ffparams.functype[0] == F_BHAM);
    b14   = (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 ||
             gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0);

    ncon           = gmx_mtop_ftype_count(mtop, F_CONSTR);
    nset           = gmx_mtop_ftype_count(mtop, F_SETTLE);
    bConstr        = (ncon > 0 || nset > 0) && !isRerun;
    bConstrVir     = false;
    if (bConstr)
    {
        if (ncon > 0 && ir->eConstrAlg == econtLINCS)
        {
            nCrmsd = 1;
        }
        bConstrVir = (getenv("GMX_CONSTRAINTVIR") != nullptr);
    }
    else
    {
        nCrmsd = 0;
    }

    /* Energy monitoring */
    for (i = 0; i < egNR; i++)
    {
        bEInd[i] = false;
    }

    for (i = 0; i < F_NRE; i++)
    {
        bEner[i] = false;
        if (isRerun &&
            (i == F_EKIN || i == F_ETOT || i == F_ECONSERVED ||
             i == F_TEMP || i == F_PDISPCORR || i == F_PRES))
        {
            continue;
        }
        if (i == F_LJ)
        {
            bEner[i] = !bBHAM;
        }
        else if (i == F_BHAM)
        {
            bEner[i] = bBHAM;
        }
        else if (i == F_EQM)
        {
            bEner[i] = ir->bQMMM;
        }
        else if (i == F_RF_EXCL)
        {
            bEner[i] = (EEL_RF(ir->coulombtype) && ir->cutoff_scheme == ecutsGROUP);
        }
        else if (i == F_COUL_RECIP)
        {
            bEner[i] = EEL_FULL(ir->coulombtype);
        }
        else if (i == F_LJ_RECIP)
        {
            bEner[i] = EVDW_PME(ir->vdwtype);
        }
        else if (i == F_LJ14)
        {
            bEner[i] = b14;
        }
        else if (i == F_COUL14)
        {
            bEner[i] = b14;
        }
        else if (i == F_LJC14_Q || i == F_LJC_PAIRS_NB)
        {
            bEner[i] = false;
        }
        else if ((i == F_DVDL_COUL && ir->fepvals->separate_dvdl[efptCOUL]) ||
                 (i == F_DVDL_VDW  && ir->fepvals->separate_dvdl[efptVDW]) ||
                 (i == F_DVDL_BONDED && ir->fepvals->separate_dvdl[efptBONDED]) ||
                 (i == F_DVDL_RESTRAINT && ir->fepvals->separate_dvdl[efptRESTRAINT]) ||
                 (i == F_DKDL && ir->fepvals->separate_dvdl[efptMASS]) ||
                 (i == F_DVDL && ir->fepvals->separate_dvdl[efptFEP]))
        {
            bEner[i] = (ir->efep != efepNO);
        }
        else if ((interaction_function[i].flags & IF_VSITE) ||
                 (i == F_CONSTR) || (i == F_CONSTRNC) || (i == F_SETTLE))
        {
            bEner[i] = false;
        }
        else if ((i == F_COUL_SR) || (i == F_EPOT) || (i == F_PRES)  || (i == F_EQM))
        {
            bEner[i] = true;
        }
        else if ((i == F_ETOT) || (i == F_EKIN) || (i == F_TEMP))
        {
            bEner[i] = EI_DYNAMICS(ir->eI);
        }
        else if (i == F_DISPCORR || i == F_PDISPCORR)
        {
            bEner[i] = (ir->eDispCorr != edispcNO);
        }
        else if (i == F_DISRESVIOL)
        {
            bEner[i] = (gmx_mtop_ftype_count(mtop, F_DISRES) > 0);
        }
        else if (i == F_ORIRESDEV)
        {
            bEner[i] = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
        }
        else if (i == F_CONNBONDS)
        {
            bEner[i] = false;
        }
        else if (i == F_COM_PULL)
        {
            bEner[i] = ((ir->bPull && pull_have_potential(ir->pull_work)) ||
                        ir->bRot);
        }
        else if (i == F_ECONSERVED)
        {
            bEner[i] = (integratorHasConservedEnergyQuantity(ir));
        }
        else
        {
            bEner[i] = (gmx_mtop_ftype_count(mtop, i) > 0);
        }
    }

    f_nre = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (bEner[i])
        {
            ener_nm[f_nre] = interaction_function[i].longname;
            f_nre++;
        }
    }

    epc            = isRerun ? epcNO : ir->epc;
    bDiagPres      = !TRICLINIC(ir->ref_p) && !isRerun;
    ref_p          = (ir->ref_p[XX][XX]+ir->ref_p[YY][YY]+ir->ref_p[ZZ][ZZ])/DIM;
    bTricl         = TRICLINIC(ir->compress) || TRICLINIC(ir->deform);
    bDynBox        = inputrecDynamicBox(ir);
    etc            = isRerun ? etcNO : ir->etc;
    bNHC_trotter   = inputrecNvtTrotter(ir) && !isRerun;
    bPrintNHChains = ir->bPrintNHChains && !isRerun;
    bMTTK          = (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)) && !isRerun;
    bMu            = inputrecNeedMutot(ir);
    bPres          = !isRerun;

    /* Pass NULL for unit to let get_ebin_space determine the units
     * for interaction_function[i].longname
     */
    ie    = get_ebin_space(ebin, f_nre, ener_nm, nullptr);
    if (nCrmsd)
    {
        /* This should be called directly after the call for ie,
         * such that iconrmsd follows directly in the list.
         */
        iconrmsd = get_ebin_space(ebin, nCrmsd, conrmsd_nm, "");
    }
    if (bDynBox)
    {
        ib    = get_ebin_space(ebin,
                               bTricl ? tricl_boxs_nm.size() : boxs_nm.size(),
                               bTricl ? tricl_boxs_nm : boxs_nm,
                               unit_length);
        ivol  = get_ebin_space(ebin, 1, vol_nm,  unit_volume);
        idens = get_ebin_space(ebin, 1, dens_nm, unit_density_SI);
        if (bDiagPres)
        {
            ipv       = get_ebin_space(ebin, 1, pv_nm,   unit_energy);
            ienthalpy = get_ebin_space(ebin, 1, enthalpy_nm,   unit_energy);
        }
    }
    if (bConstrVir)
    {
        isvir = get_ebin_space(ebin, asize(sv_nm), sv_nm, unit_energy);
        ifvir = get_ebin_space(ebin, asize(fv_nm), fv_nm, unit_energy);
    }
    if (bPres)
    {
        ivir   = get_ebin_space(ebin, asize(vir_nm), vir_nm, unit_energy);
        ipres  = get_ebin_space(ebin, asize(pres_nm), pres_nm, unit_pres_bar);
        isurft = get_ebin_space(ebin, asize(surft_nm), surft_nm,
                                unit_surft_bar);
    }
    if (epc == epcPARRINELLORAHMAN || epc == epcMTTK)
    {
        ipc = get_ebin_space(ebin, bTricl ? boxvel_nm.size() : DIM,
                             boxvel_nm, unit_vel);
    }
    if (bMu)
    {
        imu    = get_ebin_space(ebin, asize(mu_nm), mu_nm, unit_dipole_D);
    }
    if (ir->cos_accel != 0)
    {
        ivcos = get_ebin_space(ebin, asize(vcos_nm), vcos_nm, unit_vel);
        ivisc = get_ebin_space(ebin, asize(visc_nm), visc_nm,
                               unit_invvisc_SI);
    }

    /* Energy monitoring */
    for (i = 0; i < egNR; i++)
    {
        bEInd[i] = false;
    }
    bEInd[egCOULSR] = true;
    bEInd[egLJSR  ] = true;

    if (bBHAM)
    {
        bEInd[egLJSR]   = false;
        bEInd[egBHAMSR] = true;
    }
    if (b14)
    {
        bEInd[egLJ14]   = true;
        bEInd[egCOUL14] = true;
    }
    nEc = 0;
    for (i = 0; (i < egNR); i++)
    {
        if (bEInd[i])
        {
            nEc++;
        }
    }

    n       = groups->grps[egcENER].nr;
    nEg     = n;
    nE      = (n*(n+1))/2;

    snew(igrp, nE);
    if (nE > 1)
    {
        n = 0;
        snew(gnm, nEc);
        for (k = 0; (k < nEc); k++)
        {
            snew(gnm[k], STRLEN);
        }
        for (i = 0; (i < groups->grps[egcENER].nr); i++)
        {
            ni = groups->grps[egcENER].nm_ind[i];
            for (j = i; (j < groups->grps[egcENER].nr); j++)
            {
                nj = groups->grps[egcENER].nm_ind[j];
                for (k = kk = 0; (k < egNR); k++)
                {
                    if (bEInd[k])
                    {
                        sprintf(gnm[kk], "%s:%s-%s", egrp_nm[k],
                                *(groups->grpname[ni]), *(groups->grpname[nj]));
                        kk++;
                    }
                }
                igrp[n] = get_ebin_space(ebin, nEc,
                                         gnm, unit_energy);
                n++;
            }
        }
        for (k = 0; (k < nEc); k++)
        {
            sfree(gnm[k]);
        }
        sfree(gnm);

        if (n != nE)
        {
            gmx_incons("Number of energy terms wrong");
        }
    }

    nTC  = isRerun ? 0 : groups->grps[egcTC].nr;
    nNHC = ir->opts.nhchainlength; /* shorthand for number of NH chains */
    if (bMTTK)
    {
        nTCP = 1;  /* assume only one possible coupling system for barostat
                          for now */
    }
    else
    {
        nTCP = 0;
    }
    if (etc == etcNOSEHOOVER)
    {
        if (bNHC_trotter)
        {
            mde_n = 2*nNHC*nTC;
        }
        else
        {
            mde_n = 2*nTC;
        }
        if (epc == epcMTTK)
        {
            mdeb_n = 2*nNHC*nTCP;
        }
    }
    else
    {
        mde_n  = nTC;
        mdeb_n = 0;
    }

    snew(tmp_r, mde_n);
    snew(tmp_v, mde_n);
    char **grpnms;
    snew(grpnms, mde_n);

    for (i = 0; (i < nTC); i++)
    {
        ni = groups->grps[egcTC].nm_ind[i];
        sprintf(buf, "T-%s", *(groups->grpname[ni]));
        grpnms[i] = gmx_strdup(buf);
    }
    itemp = get_ebin_space(ebin, nTC, grpnms,
                           unit_temp_K);

    if (etc == etcNOSEHOOVER)
    {
        if (bPrintNHChains)
        {
            if (bNHC_trotter)
            {
                for (i = 0; (i < nTC); i++)
                {
                    ni   = groups->grps[egcTC].nm_ind[i];
                    bufi = *(groups->grpname[ni]);
                    for (j = 0; (j < nNHC); j++)
                    {
                        sprintf(buf, "Xi-%d-%s", j, bufi);
                        grpnms[2*(i*nNHC+j)] = gmx_strdup(buf);
                        sprintf(buf, "vXi-%d-%s", j, bufi);
                        grpnms[2*(i*nNHC+j)+1] = gmx_strdup(buf);
                    }
                }
                itc = get_ebin_space(ebin, mde_n,
                                     grpnms, unit_invtime);
                if (bMTTK)
                {
                    for (i = 0; (i < nTCP); i++)
                    {
                        bufi = baro_nm[0];  /* All barostat DOF's together for now. */
                        for (j = 0; (j < nNHC); j++)
                        {
                            sprintf(buf, "Xi-%d-%s", j, bufi);
                            grpnms[2*(i*nNHC+j)] = gmx_strdup(buf);
                            sprintf(buf, "vXi-%d-%s", j, bufi);
                            grpnms[2*(i*nNHC+j)+1] = gmx_strdup(buf);
                        }
                    }
                    itcb = get_ebin_space(ebin, mdeb_n,
                                          grpnms, unit_invtime);
                }
            }
            else
            {
                for (i = 0; (i < nTC); i++)
                {
                    ni   = groups->grps[egcTC].nm_ind[i];
                    bufi = *(groups->grpname[ni]);
                    sprintf(buf, "Xi-%s", bufi);
                    grpnms[2*i] = gmx_strdup(buf);
                    sprintf(buf, "vXi-%s", bufi);
                    grpnms[2*i+1] = gmx_strdup(buf);
                }
                itc = get_ebin_space(ebin, mde_n,
                                     grpnms, unit_invtime);
            }
        }
    }
    else if (etc == etcBERENDSEN || etc == etcYES ||
             etc == etcVRESCALE)
    {
        for (i = 0; (i < nTC); i++)
        {
            ni = groups->grps[egcTC].nm_ind[i];
            sprintf(buf, "Lamb-%s", *(groups->grpname[ni]));
            grpnms[i] = gmx_strdup(buf);
        }
        itc = get_ebin_space(ebin, mde_n, grpnms, "");
    }

    for (i = 0; i < mde_n; i++)
    {
        sfree(grpnms[i]);
    }
    sfree(grpnms);

    nU = groups->grps[egcACC].nr;
    if (nU > 1)
    {
        snew(grpnms, 3*nU);
        for (i = 0; (i < nU); i++)
        {
            ni = groups->grps[egcACC].nm_ind[i];
            sprintf(buf, "Ux-%s", *(groups->grpname[ni]));
            grpnms[3*i+XX] = gmx_strdup(buf);
            sprintf(buf, "Uy-%s", *(groups->grpname[ni]));
            grpnms[3*i+YY] = gmx_strdup(buf);
            sprintf(buf, "Uz-%s", *(groups->grpname[ni]));
            grpnms[3*i+ZZ] = gmx_strdup(buf);
        }
        iu = get_ebin_space(ebin, 3*nU, grpnms, unit_vel);
        sfree(grpnms);
    }

    print_grpnms = nullptr;

    /* check whether we're going to write dh histograms */
    dhc = nullptr;
    if (ir->fepvals->separate_dhdl_file == esepdhdlfileNO)
    {
        /* Currently dh histograms are only written with dynamics */
        if (EI_DYNAMICS(ir->eI))
        {
            snew(dhc, 1);

            mde_delta_h_coll_init(dhc, ir);
        }
        fp_dhdl = nullptr;
        snew(dE, ir->fepvals->n_lambda);
    }
    else
    {
        fp_dhdl = fp_dhdl_arg;
        snew(dE, ir->fepvals->n_lambda);
    }
    if (ir->bSimTemp)
    {
        int i;
        snew(temperatures, ir->fepvals->n_lambda);
        for (i = 0; i < ir->fepvals->n_lambda; i++)
        {
            temperatures[i] = ir->simtempvals->temperatures[i];
        }
    }
}

EnergyOutput::Impl::~Impl()
{
    sfree(igrp);
    sfree(tmp_r);
    sfree(tmp_v);
    done_ebin(ebin);
    done_mde_delta_h_coll(dhc);
    sfree(dE);
    sfree(temperatures);
}

void EnergyOutput::writeEnergyFileHeader(ener_file *fp_ene)
{
    if (fp_ene)
    {
        do_enxnms(fp_ene, &impl_->ebin->nener, &impl_->ebin->enm);
    }
}

int EnergyOutput::numEnergyTerms() const
{
    return impl_->ebin->nener;
}

} // namespace gmx

/*! Print a lambda vector to a string
 *
 * \param[in] fep                the inputrec's FEP input data
 * \param[in] i                  the index of the lambda vector
 * \param[in] get_native_lambda  whether to print the native lambda
 * \param[in] get_names          whether to print the names rather than the values
 * \param[in] str                the pre-allocated string buffer to print to. */
static void print_lambda_vector(t_lambda *fep, int i,
                                bool get_native_lambda, bool get_names,
                                char *str)
{
    int    j, k = 0;
    int    Nsep = 0;

    for (j = 0; j < efptNR; j++)
    {
        if (fep->separate_dvdl[j])
        {
            Nsep++;
        }
    }
    str[0] = 0; /* reset the string */
    if (Nsep > 1)
    {
        str += sprintf(str, "("); /* set the opening parenthesis*/
    }
    for (j = 0; j < efptNR; j++)
    {
        if (fep->separate_dvdl[j])
        {
            if (!get_names)
            {
                if (get_native_lambda && fep->init_lambda >= 0)
                {
                    str += sprintf(str, "%.4f", fep->init_lambda);
                }
                else
                {
                    str += sprintf(str, "%.4f", fep->all_lambda[j][i]);
                }
            }
            else
            {
                str += sprintf(str, "%s", efpt_singular_names[j]);
            }
            /* print comma for the next item */
            if (k < Nsep-1)
            {
                str += sprintf(str, ", ");
            }
            k++;
        }
    }
    if (Nsep > 1)
    {
        /* and add the closing parenthesis */
        sprintf(str, ")");
    }
}

FILE *open_dhdl(const char *filename, const t_inputrec *ir,
                const gmx_output_env_t *oenv)
{
    FILE       *fp;
    const char *dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda",
    *lambdastate     = "\\lambda state";
    int         i, nsets, nsets_de, nsetsbegin;
    int         n_lambda_terms = 0;
    t_lambda   *fep            = ir->fepvals; /* for simplicity */
    t_expanded *expand         = ir->expandedvals;
    char        lambda_vec_str[STRLEN], lambda_name_str[STRLEN];

    int         nsets_dhdl = 0;
    int         s          = 0;
    int         nsetsextend;
    bool        write_pV = false;

    /* count the number of different lambda terms */
    for (i = 0; i < efptNR; i++)
    {
        if (fep->separate_dvdl[i])
        {
            n_lambda_terms++;
        }
    }

    std::string title, label_x, label_y;
    if (fep->n_lambda == 0)
    {
        title   = gmx::formatString("%s", dhdl);
        label_x = gmx::formatString("Time (ps)");
        label_y = gmx::formatString("%s (%s %s)",
                                    dhdl, unit_energy, "[\\lambda]\\S-1\\N");
    }
    else
    {
        title   = gmx::formatString("%s and %s", dhdl, deltag);
        label_x = gmx::formatString("Time (ps)");
        label_y = gmx::formatString("%s and %s (%s %s)",
                                    dhdl, deltag, unit_energy, "[\\8l\\4]\\S-1\\N");
    }
    fp = gmx_fio_fopen(filename, "w+");
    xvgr_header(fp, title.c_str(), label_x, label_y, exvggtXNY, oenv);

    std::string buf;
    if (!(ir->bSimTemp))
    {
        buf = gmx::formatString("T = %g (K) ", ir->opts.ref_t[0]);
    }
    if ((ir->efep != efepSLOWGROWTH) && (ir->efep != efepEXPANDED))
    {
        if ( (fep->init_lambda >= 0)  && (n_lambda_terms == 1 ))
        {
            /* compatibility output */
            buf += gmx::formatString("%s = %.4f", lambda, fep->init_lambda);
        }
        else
        {
            print_lambda_vector(fep, fep->init_fep_state, true, false,
                                lambda_vec_str);
            print_lambda_vector(fep, fep->init_fep_state, true, true,
                                lambda_name_str);
            buf += gmx::formatString("%s %d: %s = %s",
                                     lambdastate, fep->init_fep_state,
                                     lambda_name_str, lambda_vec_str);
        }
    }
    xvgr_subtitle(fp, buf.c_str(), oenv);


    nsets_dhdl = 0;
    if (fep->dhdl_derivatives == edhdlderivativesYES)
    {
        nsets_dhdl = n_lambda_terms;
    }
    /* count the number of delta_g states */
    nsets_de = fep->lambda_stop_n - fep->lambda_start_n;

    nsets = nsets_dhdl + nsets_de; /* dhdl + fep differences */

    if (fep->n_lambda > 0 && (expand->elmcmove > elmcmoveNO))
    {
        nsets += 1;   /*add fep state for expanded ensemble */
    }

    if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
    {
        nsets += 1;  /* add energy to the dhdl as well */
    }

    nsetsextend = nsets;
    if ((ir->epc != epcNO) && (fep->n_lambda > 0) && (fep->init_lambda < 0))
    {
        nsetsextend += 1; /* for PV term, other terms possible if required for
                             the reduced potential (only needed with foreign
                             lambda, and only output when init_lambda is not
                             set in order to maintain compatibility of the
                             dhdl.xvg file) */
        write_pV     = true;
    }
    std::vector<std::string> setname(nsetsextend);

    if (expand->elmcmove > elmcmoveNO)
    {
        /* state for the fep_vals, if we have alchemical sampling */
        setname[s++] = "Thermodynamic state";
    }

    if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
    {
        std::string energy;
        switch (fep->edHdLPrintEnergy)
        {
            case edHdLPrintEnergyPOTENTIAL:
                energy = gmx::formatString("%s (%s)", "Potential Energy", unit_energy);
                break;
            case edHdLPrintEnergyTOTAL:
            case edHdLPrintEnergyYES:
            default:
                energy = gmx::formatString("%s (%s)", "Total Energy", unit_energy);
        }
        setname[s++] = energy;
    }

    if (fep->dhdl_derivatives == edhdlderivativesYES)
    {
        for (i = 0; i < efptNR; i++)
        {
            if (fep->separate_dvdl[i])
            {
                std::string derivative;
                if ( (fep->init_lambda >= 0)  && (n_lambda_terms == 1 ))
                {
                    /* compatibility output */
                    derivative = gmx::formatString("%s %s %.4f", dhdl, lambda, fep->init_lambda);
                }
                else
                {
                    double lam = fep->init_lambda;
                    if (fep->init_lambda < 0)
                    {
                        lam = fep->all_lambda[i][fep->init_fep_state];
                    }
                    derivative = gmx::formatString("%s %s = %.4f", dhdl, efpt_singular_names[i],
                                                   lam);
                }
                setname[s++] = derivative;
            }
        }
    }

    if (fep->n_lambda > 0)
    {
        /* g_bar has to determine the lambda values used in this simulation
         * from this xvg legend.
         */

        if (expand->elmcmove > elmcmoveNO)
        {
            nsetsbegin = 1;  /* for including the expanded ensemble */
        }
        else
        {
            nsetsbegin = 0;
        }

        if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
        {
            nsetsbegin += 1;
        }
        nsetsbegin += nsets_dhdl;

        for (i = fep->lambda_start_n; i < fep->lambda_stop_n; i++)
        {
            print_lambda_vector(fep, i, false, false, lambda_vec_str);
            std::string buf;
            if ( (fep->init_lambda >= 0)  && (n_lambda_terms == 1 ))
            {
                /* for compatible dhdl.xvg files */
                buf = gmx::formatString("%s %s %s", deltag, lambda, lambda_vec_str);
            }
            else
            {
                buf = gmx::formatString("%s %s to %s", deltag, lambda, lambda_vec_str);
            }

            if (ir->bSimTemp)
            {
                /* print the temperature for this state if doing simulated annealing */
                buf += gmx::formatString("T = %g (%s)",
                                         ir->simtempvals->temperatures[s-(nsetsbegin)],
                                         unit_temp_K);
            }
            setname[s++] = buf;
        }
        if (write_pV)
        {
            setname[s++] = gmx::formatString("pV (%s)", unit_energy);
        }

        xvgrLegend(fp, setname, oenv);
    }

    return fp;
}

namespace gmx
{

void EnergyOutput::Impl::addDataAtEnergyStep(bool                    bDoDHDL,
                                             bool                    bSum,
                                             double                  time,
                                             real                    tmass,
                                             const gmx_enerdata_t   &enerd,
                                             t_state                *state,
                                             t_lambda               *fep,
                                             t_expanded             *expand,
                                             matrix                  box,
                                             tensor                  svir,
                                             tensor                  fvir,
                                             tensor                  vir,
                                             tensor                  pres,
                                             const gmx_ekindata_t   *ekind,
                                             rvec                    mu_tot,
                                             const gmx::Constraints *constr)
{
    int    i, j, k, kk, n, gid;
    real   crmsd[2], tmp6[6];
    real   bs[NTRICLBOXS], vol, dens, pv, enthalpy;
    real   eee[egNR];
    double store_dhdl[efptNR];
    real   store_energy = 0;
    real   tmp;

    /* Do NOT use the box in the state variable, but the separate box provided
     * as an argument. This is because we sometimes need to write the box from
     * the last timestep to match the trajectory frames.
     */
    add_ebin_indexed(ebin, ie, gmx::ArrayRef<bool>(bEner), enerd.term, bSum);
    if (nCrmsd)
    {
        crmsd[0] = constr->rmsd();
        add_ebin(ebin, iconrmsd, nCrmsd, crmsd, false);
    }
    if (bDynBox)
    {
        int nboxs;
        if (bTricl)
        {
            bs[0] = box[XX][XX];
            bs[1] = box[YY][YY];
            bs[2] = box[ZZ][ZZ];
            bs[3] = box[YY][XX];
            bs[4] = box[ZZ][XX];
            bs[5] = box[ZZ][YY];
            nboxs = tricl_boxs_nm.size();
        }
        else
        {
            bs[0] = box[XX][XX];
            bs[1] = box[YY][YY];
            bs[2] = box[ZZ][ZZ];
            nboxs = boxs_nm.size();
        }
        vol  = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
        dens = (tmass*AMU)/(vol*NANO*NANO*NANO);
        add_ebin(ebin, ib, nboxs, bs, bSum);
        add_ebin(ebin, ivol, 1, &vol, bSum);
        add_ebin(ebin, idens, 1, &dens, bSum);

        if (bDiagPres)
        {
            /* This is pV (in kJ/mol).  The pressure is the reference pressure,
               not the instantaneous pressure */
            pv = vol*ref_p/PRESFAC;

            add_ebin(ebin, ipv, 1, &pv, bSum);
            enthalpy = pv + enerd.term[F_ETOT];
            add_ebin(ebin, ienthalpy, 1, &enthalpy, bSum);
        }
    }
    if (bConstrVir)
    {
        add_ebin(ebin, isvir, 9, svir[0], bSum);
        add_ebin(ebin, ifvir, 9, fvir[0], bSum);
    }
    if (bPres)
    {
        add_ebin(ebin, ivir, 9, vir[0], bSum);
        add_ebin(ebin, ipres, 9, pres[0], bSum);
        tmp = (pres[ZZ][ZZ]-(pres[XX][XX]+pres[YY][YY])*0.5)*box[ZZ][ZZ];
        add_ebin(ebin, isurft, 1, &tmp, bSum);
    }
    if (epc == epcPARRINELLORAHMAN || epc == epcMTTK)
    {
        tmp6[0] = state->boxv[XX][XX];
        tmp6[1] = state->boxv[YY][YY];
        tmp6[2] = state->boxv[ZZ][ZZ];
        tmp6[3] = state->boxv[YY][XX];
        tmp6[4] = state->boxv[ZZ][XX];
        tmp6[5] = state->boxv[ZZ][YY];
        add_ebin(ebin, ipc, bTricl ? 6 : 3, tmp6, bSum);
    }
    if (bMu)
    {
        add_ebin(ebin, imu, 3, mu_tot, bSum);
    }
    if (ekind && ekind->cosacc.cos_accel != 0)
    {
        vol  = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
        dens = (tmass*AMU)/(vol*NANO*NANO*NANO);
        add_ebin(ebin, ivcos, 1, &(ekind->cosacc.vcos), bSum);
        /* 1/viscosity, unit 1/(kg m^-1 s^-1) */
        tmp = 1/(ekind->cosacc.cos_accel/(ekind->cosacc.vcos*PICO)
                 *dens*gmx::square(box[ZZ][ZZ]*NANO/(2*M_PI)));
        add_ebin(ebin, ivisc, 1, &tmp, bSum);
    }
    if (nE > 1)
    {
        n = 0;
        for (i = 0; (i < nEg); i++)
        {
            for (j = i; (j < nEg); j++)
            {
                gid = GID(i, j, nEg);
                for (k = kk = 0; (k < egNR); k++)
                {
                    if (bEInd[k])
                    {
                        eee[kk++] = enerd.grpp.ener[k][gid];
                    }
                }
                add_ebin(ebin, igrp[n], nEc, eee, bSum);
                n++;
            }
        }
    }

    if (ekind)
    {
        for (i = 0; (i < nTC); i++)
        {
            tmp_r[i] = ekind->tcstat[i].T;
        }
        add_ebin(ebin, itemp, nTC, tmp_r, bSum);

        if (etc == etcNOSEHOOVER)
        {
            /* whether to print Nose-Hoover chains: */
            if (bPrintNHChains)
            {
                if (bNHC_trotter)
                {
                    for (i = 0; (i < nTC); i++)
                    {
                        for (j = 0; j < nNHC; j++)
                        {
                            k                = i*nNHC+j;
                            tmp_r[2*k]       = state->nosehoover_xi[k];
                            tmp_r[2*k+1]     = state->nosehoover_vxi[k];
                        }
                    }
                    add_ebin(ebin, itc, mde_n, tmp_r, bSum);

                    if (bMTTK)
                    {
                        for (i = 0; (i < nTCP); i++)
                        {
                            for (j = 0; j < nNHC; j++)
                            {
                                k                = i*nNHC+j;
                                tmp_r[2*k]       = state->nhpres_xi[k];
                                tmp_r[2*k+1]     = state->nhpres_vxi[k];
                            }
                        }
                        add_ebin(ebin, itcb, mdeb_n, tmp_r, bSum);
                    }
                }
                else
                {
                    for (i = 0; (i < nTC); i++)
                    {
                        tmp_r[2*i]   = state->nosehoover_xi[i];
                        tmp_r[2*i+1] = state->nosehoover_vxi[i];
                    }
                    add_ebin(ebin, itc, mde_n, tmp_r, bSum);
                }
            }
        }
        else if (etc == etcBERENDSEN || etc == etcYES ||
                 etc == etcVRESCALE)
        {
            for (i = 0; (i < nTC); i++)
            {
                tmp_r[i] = ekind->tcstat[i].lambda;
            }
            add_ebin(ebin, itc, nTC, tmp_r, bSum);
        }
    }

    if (ekind && nU > 1)
    {
        for (i = 0; (i < nU); i++)
        {
            copy_rvec(ekind->grpstat[i].u, tmp_v[i]);
        }
        add_ebin(ebin, iu, 3*nU, tmp_v[0], bSum);
    }

    ebin_increase_count(ebin, bSum);

    /* BAR + thermodynamic integration values */
    if ((fp_dhdl || dhc) && bDoDHDL)
    {
        for (i = 0; i < enerd.n_lambda-1; i++)
        {
            /* zero for simulated tempering */
            dE[i] = enerd.enerpart_lambda[i+1]-enerd.enerpart_lambda[0];
            if (temperatures != nullptr)
            {
                /* MRS: is this right, given the way we have defined the exchange probabilities? */
                /* is this even useful to have at all? */
                dE[i] += (temperatures[i]/
                          temperatures[state->fep_state]-1.0)*
                    enerd.term[F_EKIN];
            }
        }

        if (fp_dhdl)
        {
            fprintf(fp_dhdl, "%.4f", time);
            /* the current free energy state */

            /* print the current state if we are doing expanded ensemble */
            if (expand->elmcmove > elmcmoveNO)
            {
                fprintf(fp_dhdl, " %4d", state->fep_state);
            }
            /* total energy (for if the temperature changes */

            if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
            {
                switch (fep->edHdLPrintEnergy)
                {
                    case edHdLPrintEnergyPOTENTIAL:
                        store_energy = enerd.term[F_EPOT];
                        break;
                    case edHdLPrintEnergyTOTAL:
                    case edHdLPrintEnergyYES:
                    default:
                        store_energy = enerd.term[F_ETOT];
                }
                fprintf(fp_dhdl, " %#.8g", store_energy);
            }

            if (fep->dhdl_derivatives == edhdlderivativesYES)
            {
                for (i = 0; i < efptNR; i++)
                {
                    if (fep->separate_dvdl[i])
                    {
                        /* assumes F_DVDL is first */
                        fprintf(fp_dhdl, " %#.8g", enerd.term[F_DVDL+i]);
                    }
                }
            }
            for (i = fep->lambda_start_n; i < fep->lambda_stop_n; i++)
            {
                fprintf(fp_dhdl, " %#.8g", dE[i]);
            }
            if (bDynBox &&
                bDiagPres &&
                (epc != epcNO) &&
                (enerd.n_lambda > 0) &&
                (fep->init_lambda < 0))
            {
                fprintf(fp_dhdl, " %#.8g", pv);  /* PV term only needed when
                                                        there are alternate state
                                                        lambda and we're not in
                                                        compatibility mode */
            }
            fprintf(fp_dhdl, "\n");
            /* and the binary free energy output */
        }
        if (dhc && bDoDHDL)
        {
            int idhdl = 0;
            for (i = 0; i < efptNR; i++)
            {
                if (fep->separate_dvdl[i])
                {
                    /* assumes F_DVDL is first */
                    store_dhdl[idhdl] = enerd.term[F_DVDL+i];
                    idhdl            += 1;
                }
            }
            store_energy = enerd.term[F_ETOT];
            /* store_dh is dE */
            mde_delta_h_coll_add_dh(dhc,
                                    static_cast<double>(state->fep_state),
                                    store_energy,
                                    pv,
                                    store_dhdl,
                                    dE + fep->lambda_start_n,
                                    time);
        }
    }
}

void EnergyOutput::Impl::recordNonEnergyStep()
{
    ebin_increase_count(ebin, false);
}

namespace
{

//! Legacy output function
void npr(FILE *log, int n, char c)
{
    for (; (n > 0); n--)
    {
        fprintf(log, "%c", c);
    }
}

//! Legacy output function
void pprint(FILE *log, const char *s, t_ebin *ebin)
{
    char CHAR = '#';
    int  slen;
    char buf1[22], buf2[22];

    slen = strlen(s);
    fprintf(log, "\t<======  ");
    npr(log, slen, CHAR);
    fprintf(log, "  ==>\n");
    fprintf(log, "\t<====  %s  ====>\n", s);
    fprintf(log, "\t<==  ");
    npr(log, slen, CHAR);
    fprintf(log, "  ======>\n\n");

    fprintf(log, "\tStatistics over %s steps using %s frames\n",
            gmx_step_str(ebin->nsteps_sim, buf1),
            gmx_step_str(ebin->nsum_sim, buf2));
    fprintf(log, "\n");
}

}   // namespace

void print_ebin_header(FILE *log, int64_t steps, double time)
{
    char buf[22];

    fprintf(log, "   %12s   %12s\n"
            "   %12s   %12.5f\n\n",
            "Step", "Time", gmx_step_str(steps, buf), time);
}

// TODO It is too many responsibilities for this function to handle
// both .edr and .log output for both per-time and time-average data.
void EnergyOutput::Impl::printStepToEnergyFile(ener_file *fp_ene, bool bEne, bool bDR, bool bOR,
                                               FILE *log,
                                               int64_t step, double time,
                                               int mode,
                                               t_fcdata *fcd,
                                               gmx_groups_t *groups, t_grpopts *opts,
                                               gmx::Awh *awh)
{
    /*static char **grpnms=NULL;*/
    char         buf[246];
    int          i, j, n, ni, nj, b;
    int          ndisre = 0;

    /* these are for the old-style blocks (1 subblock, only reals), because
       there can be only one per ID for these */
    int          nr[enxNR];
    int          id[enxNR];
    real        *block[enxNR];

    t_enxframe   fr;

    if (mode == eprAVER && ebin->nsum_sim <= 0)
    {
        if (log)
        {
            fprintf(log, "Not enough data recorded to report energy averages\n");
        }
        return;
    }

    switch (mode)
    {
        case eprNORMAL:
            init_enxframe(&fr);
            fr.t            = time;
            fr.step         = step;
            fr.nsteps       = ebin->nsteps;
            fr.dt           = delta_t;
            fr.nsum         = ebin->nsum;
            fr.nre          = (bEne) ? ebin->nener : 0;
            fr.ener         = ebin->e;
            ndisre          = bDR ? fcd->disres.npair : 0;
            /* Optional additional old-style (real-only) blocks. */
            for (i = 0; i < enxNR; i++)
            {
                nr[i] = 0;
            }
            if (bOR && fcd->orires.nr > 0)
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
                for (i = 0; i < enxNR; i++)
                {
                    if (nr[i] > 0)
                    {
                        fr.nblock = i + 1;
                    }
                }
                add_blocks_enxframe(&fr, fr.nblock);
                for (b = 0; b < fr.nblock; b++)
                {
                    add_subblocks_enxblock(&(fr.block[b]), 1);
                    fr.block[b].id        = id[b];
                    fr.block[b].sub[0].nr = nr[b];
#if !GMX_DOUBLE
                    fr.block[b].sub[0].type = xdr_datatype_float;
                    fr.block[b].sub[0].fval = block[b];
#else
                    fr.block[b].sub[0].type = xdr_datatype_double;
                    fr.block[b].sub[0].dval = block[b];
#endif
                }

                /* check for disre block & fill it. */
                if (ndisre > 0)
                {
                    int db = fr.nblock;
                    fr.nblock += 1;
                    add_blocks_enxframe(&fr, fr.nblock);

                    add_subblocks_enxblock(&(fr.block[db]), 2);
                    fr.block[db].id        = enxDISRE;
                    fr.block[db].sub[0].nr = ndisre;
                    fr.block[db].sub[1].nr = ndisre;
#if !GMX_DOUBLE
                    fr.block[db].sub[0].type = xdr_datatype_float;
                    fr.block[db].sub[1].type = xdr_datatype_float;
                    fr.block[db].sub[0].fval = fcd->disres.rt;
                    fr.block[db].sub[1].fval = fcd->disres.rm3tav;
#else
                    fr.block[db].sub[0].type = xdr_datatype_double;
                    fr.block[db].sub[1].type = xdr_datatype_double;
                    fr.block[db].sub[0].dval = fcd->disres.rt;
                    fr.block[db].sub[1].dval = fcd->disres.rm3tav;
#endif
                }
                /* here we can put new-style blocks */

                /* Free energy perturbation blocks */
                if (dhc)
                {
                    mde_delta_h_coll_handle_block(dhc, &fr, fr.nblock);
                }

                /* we can now free & reset the data in the blocks */
                if (dhc)
                {
                    mde_delta_h_coll_reset(dhc);
                }

                /* AWH bias blocks. */
                if (awh != nullptr)  // TODO: add boolean flag. See in energyoutput.h
                {
                    awh->writeToEnergyFrame(step, &fr);
                }

                /* do the actual I/O */
                do_enx(fp_ene, &fr);
                if (fr.nre)
                {
                    /* We have stored the sums, so reset the sum history */
                    reset_ebin_sums(ebin);
                }
            }
            free_enxframe(&fr);
            break;
        case eprAVER:
            if (log)
            {
                pprint(log, "A V E R A G E S", ebin);
            }
            break;
        case eprRMS:
            if (log)
            {
                pprint(log, "R M S - F L U C T U A T I O N S", ebin);
            }
            break;
        default:
            gmx_fatal(FARGS, "Invalid print mode (%d)", mode);
    }

    if (log)
    {
        if (opts)
        {
            for (i = 0; i < opts->ngtc; i++)
            {
                if (opts->annealing[i] != eannNO)
                {
                    fprintf(log, "Current ref_t for group %s: %8.1f\n",
                            *(groups->grpname[groups->grps[egcTC].nm_ind[i]]),
                            opts->ref_t[i]);
                }
            }
        }
        if (mode == eprNORMAL && bOR && fcd->orires.nr > 0)
        {
            print_orires_log(log, &(fcd->orires));
        }
        fprintf(log, "   Energies (%s)\n", unit_energy);
        pr_ebin(log, ebin, ie, f_nre+nCrmsd, 5, mode, true);
        fprintf(log, "\n");

        if (mode == eprAVER)
        {
            if (bDynBox)
            {
                pr_ebin(log, ebin, ib, bTricl ? tricl_boxs_nm.size() : boxs_nm.size(), 5,
                        mode, true);
                fprintf(log, "\n");
            }
            if (bConstrVir)
            {
                fprintf(log, "   Constraint Virial (%s)\n", unit_energy);
                pr_ebin(log, ebin, isvir, 9, 3, mode, false);
                fprintf(log, "\n");
                fprintf(log, "   Force Virial (%s)\n", unit_energy);
                pr_ebin(log, ebin, ifvir, 9, 3, mode, false);
                fprintf(log, "\n");
            }
            if (bPres)
            {
                fprintf(log, "   Total Virial (%s)\n", unit_energy);
                pr_ebin(log, ebin, ivir, 9, 3, mode, false);
                fprintf(log, "\n");
                fprintf(log, "   Pressure (%s)\n", unit_pres_bar);
                pr_ebin(log, ebin, ipres, 9, 3, mode, false);
                fprintf(log, "\n");
            }
            if (bMu)
            {
                fprintf(log, "   Total Dipole (%s)\n", unit_dipole_D);
                pr_ebin(log, ebin, imu, 3, 3, mode, false);
                fprintf(log, "\n");
            }

            if (nE > 1)
            {
                if (print_grpnms == nullptr)
                {
                    snew(print_grpnms, nE);
                    n = 0;
                    for (i = 0; (i < nEg); i++)
                    {
                        ni = groups->grps[egcENER].nm_ind[i];
                        for (j = i; (j < nEg); j++)
                        {
                            nj = groups->grps[egcENER].nm_ind[j];
                            sprintf(buf, "%s-%s", *(groups->grpname[ni]),
                                    *(groups->grpname[nj]));
                            print_grpnms[n++] = gmx_strdup(buf);
                        }
                    }
                }
                sprintf(buf, "Epot (%s)", unit_energy);
                fprintf(log, "%15s   ", buf);
                for (i = 0; (i < egNR); i++)
                {
                    if (bEInd[i])
                    {
                        fprintf(log, "%12s   ", egrp_nm[i]);
                    }
                }
                fprintf(log, "\n");
                for (i = 0; (i < nE); i++)
                {
                    fprintf(log, "%15s", print_grpnms[i]);
                    pr_ebin(log, ebin, igrp[i], nEc, nEc, mode,
                            false);
                }
                fprintf(log, "\n");
            }
            if (nTC > 1)
            {
                pr_ebin(log, ebin, itemp, nTC, 4, mode, true);
                fprintf(log, "\n");
            }
            if (nU > 1)
            {
                fprintf(log, "%15s   %12s   %12s   %12s\n",
                        "Group", "Ux", "Uy", "Uz");
                for (i = 0; (i < nU); i++)
                {
                    ni = groups->grps[egcACC].nm_ind[i];
                    fprintf(log, "%15s", *groups->grpname[ni]);
                    pr_ebin(log, ebin, iu+3*i, 3, 3, mode, false);
                }
                fprintf(log, "\n");
            }
        }
    }

}

void EnergyOutput::Impl::fillEnergyHistory(energyhistory_t *enerhist) const
{
    enerhist->nsteps     = ebin->nsteps;
    enerhist->nsum       = ebin->nsum;
    enerhist->nsteps_sim = ebin->nsteps_sim;
    enerhist->nsum_sim   = ebin->nsum_sim;

    if (ebin->nsum > 0)
    {
        /* This will only actually resize the first time */
        enerhist->ener_ave.resize(ebin->nener);
        enerhist->ener_sum.resize(ebin->nener);

        for (int i = 0; i < ebin->nener; i++)
        {
            enerhist->ener_ave[i] = ebin->e[i].eav;
            enerhist->ener_sum[i] = ebin->e[i].esum;
        }
    }

    if (ebin->nsum_sim > 0)
    {
        /* This will only actually resize the first time */
        enerhist->ener_sum_sim.resize(ebin->nener);

        for (int i = 0; i < ebin->nener; i++)
        {
            enerhist->ener_sum_sim[i] = ebin->e_sim[i].esum;
        }
    }
    if (dhc)
    {
        mde_delta_h_coll_update_energyhistory(dhc, enerhist);
    }
}

void EnergyOutput::Impl::restoreFromEnergyHistory(const energyhistory_t &enerhist)
{
    std::size_t nener = ebin->nener;

    if ((enerhist.nsum     > 0 && nener != enerhist.ener_sum.size()) ||
        (enerhist.nsum_sim > 0 && nener != enerhist.ener_sum_sim.size()))
    {
        gmx_fatal(FARGS, "Mismatch between number of energies in run input (%zu) and checkpoint file (%zu or %zu).",
                  nener, enerhist.ener_sum.size(), enerhist.ener_sum_sim.size());
    }

    ebin->nsteps     = enerhist.nsteps;
    ebin->nsum       = enerhist.nsum;
    ebin->nsteps_sim = enerhist.nsteps_sim;
    ebin->nsum_sim   = enerhist.nsum_sim;

    for (int i = 0; i < ebin->nener; i++)
    {
        ebin->e[i].eav  =
            (enerhist.nsum > 0 ? enerhist.ener_ave[i] : 0);
        ebin->e[i].esum =
            (enerhist.nsum > 0 ? enerhist.ener_sum[i] : 0);
        ebin->e_sim[i].esum =
            (enerhist.nsum_sim > 0 ? enerhist.ener_sum_sim[i] : 0);
    }
    if (dhc)
    {
        mde_delta_h_coll_restore_energyhistory(dhc, enerhist.deltaHForeignLambdas.get());
    }
}

EnergyOutput::EnergyOutput()
    : impl_(new Impl())
{
}

void EnergyOutput::prepare(const gmx_mtop_t *mtop,
                           const t_inputrec *ir,
                           FILE             *fp_dhdl,
                           bool              isRerun)
{
    impl_->prepare(mtop, ir, fp_dhdl, isRerun);
}

EnergyOutput::~EnergyOutput() = default;

t_ebin *EnergyOutput::getEbin()
{
    return impl_->ebin;
}

void EnergyOutput::addDataAtEnergyStep(bool                    bDoDHDL,
                                       bool                    bSum,
                                       double                  time,
                                       real                    tmass,
                                       const gmx_enerdata_t   &enerd,
                                       t_state                *state,
                                       t_lambda               *fep,
                                       t_expanded             *expand,
                                       matrix                  box,
                                       tensor                  svir,
                                       tensor                  fvir,
                                       tensor                  vir,
                                       tensor                  pres,
                                       const gmx_ekindata_t   *ekind,
                                       rvec                    mu_tot,
                                       const gmx::Constraints *constr)
{
    impl_->addDataAtEnergyStep(bDoDHDL, bSum, time, tmass, enerd, state, fep,
                               expand, box, svir, fvir, vir, pres, ekind, mu_tot, constr);
}

void EnergyOutput::recordNonEnergyStep()
{
    impl_->recordNonEnergyStep();
}

void EnergyOutput::printStepToEnergyFile(ener_file *fp_ene, bool bEne, bool bDR, bool bOR,
                                         FILE *log,
                                         int64_t step, double time,
                                         int mode,
                                         t_fcdata *fcd,
                                         gmx_groups_t *groups, t_grpopts *opts,
                                         gmx::Awh *awh)
{
    impl_->printStepToEnergyFile(fp_ene, bEne, bDR, bOR, log, step, time, mode,
                                 fcd, groups, opts, awh);
}

void EnergyOutput::fillEnergyHistory(energyhistory_t *enerhist) const
{
    impl_->fillEnergyHistory(enerhist);
}

void EnergyOutput::restoreFromEnergyHistory(const energyhistory_t &enerhist)
{
    impl_->restoreFromEnergyHistory(enerhist);
}

} // namespace gmx
