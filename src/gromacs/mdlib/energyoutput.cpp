/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "energyoutput.h"

#include <cfloat>
#include <cstdlib>
#include <cstring>

#include <array>
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
#include "gromacs/mdrunutility/handlerestart.h"
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
#include "gromacs/utility/mdmodulenotification.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

//! Labels for energy file quantities
//! \{
static const char* conrmsd_nm[] = { "Constr. rmsd", "Constr.2 rmsd" };

static std::array<const char*, 3> boxs_nm = { "Box-X", "Box-Y", "Box-Z" };

static std::array<const char*, 6> tricl_boxs_nm = { "Box-XX", "Box-YY", "Box-ZZ",
                                                    "Box-YX", "Box-ZX", "Box-ZY" };

static const char* vol_nm[] = { "Volume" };

static const char* dens_nm[] = { "Density" };

static const char* pv_nm[] = { "pV" };

static const char* enthalpy_nm[] = { "Enthalpy" };

static std::array<const char*, 6> boxvel_nm = { "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
                                                "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY" };

const char* egrp_nm[egNR + 1] = { "Coul-SR", "LJ-SR", "Buck-SR", "Coul-14", "LJ-14", nullptr };
//! \}

namespace gmx
{

/*! \brief Energy output class
 *
 * This is the collection of energy averages collected during mdrun, and to
 * be written out to the .edr file.
 *
 * \todo Use more std containers.
 * \todo Remove GMX_CONSTRAINTVIR
 * \todo Write free-energy output also to energy file (after adding more tests)
 */
EnergyOutput::EnergyOutput(ener_file*               fp_ene,
                           const gmx_mtop_t*        mtop,
                           const t_inputrec*        ir,
                           const pull_t*            pull_work,
                           FILE*                    fp_dhdl,
                           bool                     isRerun,
                           const StartingBehavior   startingBehavior,
                           const MdModulesNotifier& mdModulesNotifier)
{
    const char*        ener_nm[F_NRE];
    static const char* vir_nm[]   = { "Vir-XX", "Vir-XY", "Vir-XZ", "Vir-YX", "Vir-YY",
                                    "Vir-YZ", "Vir-ZX", "Vir-ZY", "Vir-ZZ" };
    static const char* sv_nm[]    = { "ShakeVir-XX", "ShakeVir-XY", "ShakeVir-XZ",
                                   "ShakeVir-YX", "ShakeVir-YY", "ShakeVir-YZ",
                                   "ShakeVir-ZX", "ShakeVir-ZY", "ShakeVir-ZZ" };
    static const char* fv_nm[]    = { "ForceVir-XX", "ForceVir-XY", "ForceVir-XZ",
                                   "ForceVir-YX", "ForceVir-YY", "ForceVir-YZ",
                                   "ForceVir-ZX", "ForceVir-ZY", "ForceVir-ZZ" };
    static const char* pres_nm[]  = { "Pres-XX", "Pres-XY", "Pres-XZ", "Pres-YX", "Pres-YY",
                                     "Pres-YZ", "Pres-ZX", "Pres-ZY", "Pres-ZZ" };
    static const char* surft_nm[] = { "#Surf*SurfTen" };
    static const char* mu_nm[]    = { "Mu-X", "Mu-Y", "Mu-Z" };
    static const char* vcos_nm[]  = { "2CosZ*Vel-X" };
    static const char* visc_nm[]  = { "1/Viscosity" };
    static const char* baro_nm[]  = { "Barostat" };

    const SimulationGroups* groups;
    char**                  gnm;
    char                    buf[256];
    const char*             bufi;
    int                     i, j, ni, nj, n, k, kk, ncon, nset;
    bool                    bBHAM, b14;

    if (EI_DYNAMICS(ir->eI))
    {
        delta_t_ = ir->delta_t;
    }
    else
    {
        delta_t_ = 0;
    }

    groups = &mtop->groups;

    bBHAM = (mtop->ffparams.numTypes() > 0) && (mtop->ffparams.functype[0] == F_BHAM);
    b14   = (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 || gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0);

    ncon         = gmx_mtop_ftype_count(mtop, F_CONSTR);
    nset         = gmx_mtop_ftype_count(mtop, F_SETTLE);
    bool bConstr = (ncon > 0 || nset > 0) && !isRerun;
    bConstrVir_  = false;
    nCrmsd_      = 0;
    if (bConstr)
    {
        if (ncon > 0 && ir->eConstrAlg == econtLINCS)
        {
            nCrmsd_ = 1;
        }
        bConstrVir_ = (getenv("GMX_CONSTRAINTVIR") != nullptr);
    }
    else
    {
        nCrmsd_ = 0;
    }

    /* Energy monitoring */
    for (i = 0; i < egNR; i++)
    {
        bEInd_[i] = false;
    }

    // Setting true only to those energy terms, that have active interactions and
    // are not vsite terms (not VSITE2, VSITE3, VSITE3FD, VSITE3FAD, VSITE3OUT, VSITE4FD, VSITE4FDN, or VSITEN)
    for (i = 0; i < F_NRE; i++)
    {
        bEner_[i] = (gmx_mtop_ftype_count(mtop, i) > 0)
                    && ((interaction_function[i].flags & IF_VSITE) == 0);
    }

    if (!isRerun)
    {
        bEner_[F_EKIN] = EI_DYNAMICS(ir->eI);
        bEner_[F_ETOT] = EI_DYNAMICS(ir->eI);
        bEner_[F_TEMP] = EI_DYNAMICS(ir->eI);

        bEner_[F_ECONSERVED] = integratorHasConservedEnergyQuantity(ir);
        bEner_[F_PDISPCORR]  = (ir->eDispCorr != edispcNO);
        bEner_[F_PRES]       = true;
    }

    bEner_[F_LJ]           = !bBHAM;
    bEner_[F_BHAM]         = bBHAM;
    bEner_[F_EQM]          = ir->bQMMM;
    bEner_[F_RF_EXCL]      = (EEL_RF(ir->coulombtype) && ir->cutoff_scheme == ecutsGROUP);
    bEner_[F_COUL_RECIP]   = EEL_FULL(ir->coulombtype);
    bEner_[F_LJ_RECIP]     = EVDW_PME(ir->vdwtype);
    bEner_[F_LJ14]         = b14;
    bEner_[F_COUL14]       = b14;
    bEner_[F_LJC14_Q]      = false;
    bEner_[F_LJC_PAIRS_NB] = false;


    bEner_[F_DVDL_COUL]      = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptCOUL];
    bEner_[F_DVDL_VDW]       = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptVDW];
    bEner_[F_DVDL_BONDED]    = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptBONDED];
    bEner_[F_DVDL_RESTRAINT] = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptRESTRAINT];
    bEner_[F_DKDL]           = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptMASS];
    bEner_[F_DVDL]           = (ir->efep != efepNO) && ir->fepvals->separate_dvdl[efptFEP];

    bEner_[F_CONSTR]   = false;
    bEner_[F_CONSTRNC] = false;
    bEner_[F_SETTLE]   = false;

    bEner_[F_COUL_SR] = true;
    bEner_[F_EPOT]    = true;

    bEner_[F_DISPCORR]   = (ir->eDispCorr != edispcNO);
    bEner_[F_DISRESVIOL] = (gmx_mtop_ftype_count(mtop, F_DISRES) > 0);
    bEner_[F_ORIRESDEV]  = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
    bEner_[F_COM_PULL]   = ((ir->bPull && pull_have_potential(pull_work)) || ir->bRot);

    MdModulesEnergyOutputToDensityFittingRequestChecker mdModulesAddOutputToDensityFittingFieldRequest;
    mdModulesNotifier.notifier_.notify(&mdModulesAddOutputToDensityFittingFieldRequest);

    bEner_[F_DENSITYFITTING] = mdModulesAddOutputToDensityFittingFieldRequest.energyOutputToDensityFitting_;


    // Counting the energy terms that will be printed and saving their names
    f_nre_ = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (bEner_[i])
        {
            ener_nm[f_nre_] = interaction_function[i].longname;
            f_nre_++;
        }
    }

    epc_            = isRerun ? epcNO : ir->epc;
    bDiagPres_      = !TRICLINIC(ir->ref_p) && !isRerun;
    ref_p_          = (ir->ref_p[XX][XX] + ir->ref_p[YY][YY] + ir->ref_p[ZZ][ZZ]) / DIM;
    bTricl_         = TRICLINIC(ir->compress) || TRICLINIC(ir->deform);
    bDynBox_        = inputrecDynamicBox(ir);
    etc_            = isRerun ? etcNO : ir->etc;
    bNHC_trotter_   = inputrecNvtTrotter(ir) && !isRerun;
    bPrintNHChains_ = ir->bPrintNHChains && !isRerun;
    bMTTK_          = (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)) && !isRerun;
    bMu_            = inputrecNeedMutot(ir);
    bPres_          = !isRerun;

    ebin_ = mk_ebin();
    /* Pass NULL for unit to let get_ebin_space determine the units
     * for interaction_function[i].longname
     */
    ie_ = get_ebin_space(ebin_, f_nre_, ener_nm, nullptr);
    if (nCrmsd_)
    {
        /* This should be called directly after the call for ie_,
         * such that iconrmsd_ follows directly in the list.
         */
        iconrmsd_ = get_ebin_space(ebin_, nCrmsd_, conrmsd_nm, "");
    }
    if (bDynBox_)
    {
        ib_    = get_ebin_space(ebin_, bTricl_ ? tricl_boxs_nm.size() : boxs_nm.size(),
                             bTricl_ ? tricl_boxs_nm.data() : boxs_nm.data(), unit_length);
        ivol_  = get_ebin_space(ebin_, 1, vol_nm, unit_volume);
        idens_ = get_ebin_space(ebin_, 1, dens_nm, unit_density_SI);
        if (bDiagPres_)
        {
            ipv_       = get_ebin_space(ebin_, 1, pv_nm, unit_energy);
            ienthalpy_ = get_ebin_space(ebin_, 1, enthalpy_nm, unit_energy);
        }
    }
    if (bConstrVir_)
    {
        isvir_ = get_ebin_space(ebin_, asize(sv_nm), sv_nm, unit_energy);
        ifvir_ = get_ebin_space(ebin_, asize(fv_nm), fv_nm, unit_energy);
    }
    if (bPres_)
    {
        ivir_   = get_ebin_space(ebin_, asize(vir_nm), vir_nm, unit_energy);
        ipres_  = get_ebin_space(ebin_, asize(pres_nm), pres_nm, unit_pres_bar);
        isurft_ = get_ebin_space(ebin_, asize(surft_nm), surft_nm, unit_surft_bar);
    }
    if (epc_ == epcPARRINELLORAHMAN || epc_ == epcMTTK)
    {
        ipc_ = get_ebin_space(ebin_, bTricl_ ? boxvel_nm.size() : DIM, boxvel_nm.data(), unit_vel);
    }
    if (bMu_)
    {
        imu_ = get_ebin_space(ebin_, asize(mu_nm), mu_nm, unit_dipole_D);
    }
    if (ir->cos_accel != 0)
    {
        ivcos_ = get_ebin_space(ebin_, asize(vcos_nm), vcos_nm, unit_vel);
        ivisc_ = get_ebin_space(ebin_, asize(visc_nm), visc_nm, unit_invvisc_SI);
    }

    /* Energy monitoring */
    for (i = 0; i < egNR; i++)
    {
        bEInd_[i] = false;
    }
    bEInd_[egCOULSR] = true;
    bEInd_[egLJSR]   = true;

    if (bBHAM)
    {
        bEInd_[egLJSR]   = false;
        bEInd_[egBHAMSR] = true;
    }
    if (b14)
    {
        bEInd_[egLJ14]   = true;
        bEInd_[egCOUL14] = true;
    }
    nEc_ = 0;
    for (i = 0; (i < egNR); i++)
    {
        if (bEInd_[i])
        {
            nEc_++;
        }
    }
    n    = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    nEg_ = n;
    nE_  = (n * (n + 1)) / 2;

    snew(igrp_, nE_);
    if (nE_ > 1)
    {
        n = 0;
        snew(gnm, nEc_);
        for (k = 0; (k < nEc_); k++)
        {
            snew(gnm[k], STRLEN);
        }
        for (i = 0; (i < gmx::ssize(groups->groups[SimulationAtomGroupType::EnergyOutput])); i++)
        {
            ni = groups->groups[SimulationAtomGroupType::EnergyOutput][i];
            for (j = i; (j < gmx::ssize(groups->groups[SimulationAtomGroupType::EnergyOutput])); j++)
            {
                nj = groups->groups[SimulationAtomGroupType::EnergyOutput][j];
                for (k = kk = 0; (k < egNR); k++)
                {
                    if (bEInd_[k])
                    {
                        sprintf(gnm[kk], "%s:%s-%s", egrp_nm[k], *(groups->groupNames[ni]),
                                *(groups->groupNames[nj]));
                        kk++;
                    }
                }
                igrp_[n] = get_ebin_space(ebin_, nEc_, gnm, unit_energy);
                n++;
            }
        }
        for (k = 0; (k < nEc_); k++)
        {
            sfree(gnm[k]);
        }
        sfree(gnm);

        if (n != nE_)
        {
            gmx_incons("Number of energy terms wrong");
        }
    }

    nTC_  = isRerun ? 0 : groups->groups[SimulationAtomGroupType::TemperatureCoupling].size();
    nNHC_ = ir->opts.nhchainlength; /* shorthand for number of NH chains */
    if (bMTTK_)
    {
        nTCP_ = 1; /* assume only one possible coupling system for barostat
                           for now */
    }
    else
    {
        nTCP_ = 0;
    }
    if (etc_ == etcNOSEHOOVER)
    {
        if (bNHC_trotter_)
        {
            mde_n_ = 2 * nNHC_ * nTC_;
        }
        else
        {
            mde_n_ = 2 * nTC_;
        }
        if (epc_ == epcMTTK)
        {
            mdeb_n_ = 2 * nNHC_ * nTCP_;
        }
    }
    else
    {
        mde_n_  = nTC_;
        mdeb_n_ = 0;
    }

    snew(tmp_r_, mde_n_);
    // TODO redo the group name memory management to make it more clear
    char** grpnms;
    snew(grpnms, std::max(mde_n_, mdeb_n_)); // Just in case mdeb_n_ > mde_n_

    for (i = 0; (i < nTC_); i++)
    {
        ni = groups->groups[SimulationAtomGroupType::TemperatureCoupling][i];
        sprintf(buf, "T-%s", *(groups->groupNames[ni]));
        grpnms[i] = gmx_strdup(buf);
    }
    itemp_ = get_ebin_space(ebin_, nTC_, grpnms, unit_temp_K);
    for (i = 0; i < nTC_; i++)
    {
        sfree(grpnms[i]);
    }

    int allocated = 0;
    if (etc_ == etcNOSEHOOVER)
    {
        if (bPrintNHChains_)
        {
            if (bNHC_trotter_)
            {
                for (i = 0; (i < nTC_); i++)
                {
                    ni   = groups->groups[SimulationAtomGroupType::TemperatureCoupling][i];
                    bufi = *(groups->groupNames[ni]);
                    for (j = 0; (j < nNHC_); j++)
                    {
                        sprintf(buf, "Xi-%d-%s", j, bufi);
                        grpnms[2 * (i * nNHC_ + j)] = gmx_strdup(buf);
                        sprintf(buf, "vXi-%d-%s", j, bufi);
                        grpnms[2 * (i * nNHC_ + j) + 1] = gmx_strdup(buf);
                    }
                }
                itc_      = get_ebin_space(ebin_, mde_n_, grpnms, unit_invtime);
                allocated = mde_n_;
                if (bMTTK_)
                {
                    for (i = 0; (i < nTCP_); i++)
                    {
                        bufi = baro_nm[0]; /* All barostat DOF's together for now. */
                        for (j = 0; (j < nNHC_); j++)
                        {
                            sprintf(buf, "Xi-%d-%s", j, bufi);
                            grpnms[2 * (i * nNHC_ + j)] = gmx_strdup(buf);
                            sprintf(buf, "vXi-%d-%s", j, bufi);
                            grpnms[2 * (i * nNHC_ + j) + 1] = gmx_strdup(buf);
                        }
                    }
                    itcb_     = get_ebin_space(ebin_, mdeb_n_, grpnms, unit_invtime);
                    allocated = mdeb_n_;
                }
            }
            else
            {
                for (i = 0; (i < nTC_); i++)
                {
                    ni   = groups->groups[SimulationAtomGroupType::TemperatureCoupling][i];
                    bufi = *(groups->groupNames[ni]);
                    sprintf(buf, "Xi-%s", bufi);
                    grpnms[2 * i] = gmx_strdup(buf);
                    sprintf(buf, "vXi-%s", bufi);
                    grpnms[2 * i + 1] = gmx_strdup(buf);
                }
                itc_      = get_ebin_space(ebin_, mde_n_, grpnms, unit_invtime);
                allocated = mde_n_;
            }
        }
    }
    else if (etc_ == etcBERENDSEN || etc_ == etcYES || etc_ == etcVRESCALE)
    {
        for (i = 0; (i < nTC_); i++)
        {
            ni = groups->groups[SimulationAtomGroupType::TemperatureCoupling][i];
            sprintf(buf, "Lamb-%s", *(groups->groupNames[ni]));
            grpnms[i] = gmx_strdup(buf);
        }
        itc_      = get_ebin_space(ebin_, mde_n_, grpnms, "");
        allocated = mde_n_;
    }

    for (i = 0; i < allocated; i++)
    {
        sfree(grpnms[i]);
    }
    sfree(grpnms);

    nU_ = groups->groups[SimulationAtomGroupType::Acceleration].size();
    snew(tmp_v_, nU_);
    if (nU_ > 1)
    {
        snew(grpnms, 3 * nU_);
        for (i = 0; (i < nU_); i++)
        {
            ni = groups->groups[SimulationAtomGroupType::Acceleration][i];
            sprintf(buf, "Ux-%s", *(groups->groupNames[ni]));
            grpnms[3 * i + XX] = gmx_strdup(buf);
            sprintf(buf, "Uy-%s", *(groups->groupNames[ni]));
            grpnms[3 * i + YY] = gmx_strdup(buf);
            sprintf(buf, "Uz-%s", *(groups->groupNames[ni]));
            grpnms[3 * i + ZZ] = gmx_strdup(buf);
        }
        iu_ = get_ebin_space(ebin_, 3 * nU_, grpnms, unit_vel);
        for (i = 0; i < 3 * nU_; i++)
        {
            sfree(grpnms[i]);
        }
        sfree(grpnms);
    }

    /* Note that fp_ene should be valid on the master rank and null otherwise */
    if (fp_ene != nullptr && startingBehavior != StartingBehavior::RestartWithAppending)
    {
        do_enxnms(fp_ene, &ebin_->nener, &ebin_->enm);
    }

    /* check whether we're going to write dh histograms */
    dhc_ = nullptr;
    if (ir->fepvals->separate_dhdl_file == esepdhdlfileNO)
    {
        /* Currently dh histograms are only written with dynamics */
        if (EI_DYNAMICS(ir->eI))
        {
            snew(dhc_, 1);

            mde_delta_h_coll_init(dhc_, ir);
        }
        fp_dhdl_ = nullptr;
        snew(dE_, ir->fepvals->n_lambda);
    }
    else
    {
        fp_dhdl_ = fp_dhdl;
        snew(dE_, ir->fepvals->n_lambda);
    }
    if (ir->bSimTemp)
    {
        int i;
        snew(temperatures_, ir->fepvals->n_lambda);
        numTemperatures_ = ir->fepvals->n_lambda;
        for (i = 0; i < ir->fepvals->n_lambda; i++)
        {
            temperatures_[i] = ir->simtempvals->temperatures[i];
        }
    }
    else
    {
        numTemperatures_ = 0;
    }
}

EnergyOutput::~EnergyOutput()
{
    sfree(igrp_);
    sfree(tmp_r_);
    sfree(tmp_v_);
    done_ebin(ebin_);
    done_mde_delta_h_coll(dhc_);
    sfree(dE_);
    if (numTemperatures_ > 0)
    {
        sfree(temperatures_);
    }
}

} // namespace gmx

/*! \brief Print a lambda vector to a string
 *
 * \param[in] fep                The inputrec's FEP input data
 * \param[in] i                  The index of the lambda vector
 * \param[in] get_native_lambda  Whether to print the native lambda
 * \param[in] get_names          Whether to print the names rather than the values
 * \param[in,out] str            The pre-allocated string buffer to print to.
 */
static void print_lambda_vector(t_lambda* fep, int i, bool get_native_lambda, bool get_names, char* str)
{
    int j, k = 0;
    int Nsep = 0;

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
            if (k < Nsep - 1)
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

FILE* open_dhdl(const char* filename, const t_inputrec* ir, const gmx_output_env_t* oenv)
{
    FILE*       fp;
    const char *dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda",
               *lambdastate = "\\lambda state";
    int         i, nsets, nsets_de, nsetsbegin;
    int         n_lambda_terms = 0;
    t_lambda*   fep            = ir->fepvals; /* for simplicity */
    t_expanded* expand         = ir->expandedvals;
    char        lambda_vec_str[STRLEN], lambda_name_str[STRLEN];

    int  nsets_dhdl = 0;
    int  s          = 0;
    int  nsetsextend;
    bool write_pV = false;

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
        label_y = gmx::formatString("%s (%s %s)", dhdl, unit_energy, "[\\lambda]\\S-1\\N");
    }
    else
    {
        title   = gmx::formatString("%s and %s", dhdl, deltag);
        label_x = gmx::formatString("Time (ps)");
        label_y = gmx::formatString("%s and %s (%s %s)", dhdl, deltag, unit_energy,
                                    "[\\8l\\4]\\S-1\\N");
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
        if ((fep->init_lambda >= 0) && (n_lambda_terms == 1))
        {
            /* compatibility output */
            buf += gmx::formatString("%s = %.4f", lambda, fep->init_lambda);
        }
        else
        {
            print_lambda_vector(fep, fep->init_fep_state, true, false, lambda_vec_str);
            print_lambda_vector(fep, fep->init_fep_state, true, true, lambda_name_str);
            buf += gmx::formatString("%s %d: %s = %s", lambdastate, fep->init_fep_state,
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
        nsets += 1; /*add fep state for expanded ensemble */
    }

    if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
    {
        nsets += 1; /* add energy to the dhdl as well */
    }

    nsetsextend = nsets;
    if ((ir->epc != epcNO) && (fep->n_lambda > 0) && (fep->init_lambda < 0))
    {
        nsetsextend += 1; /* for PV term, other terms possible if required for
                             the reduced potential (only needed with foreign
                             lambda, and only output when init_lambda is not
                             set in order to maintain compatibility of the
                             dhdl.xvg file) */
        write_pV = true;
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
            default: energy = gmx::formatString("%s (%s)", "Total Energy", unit_energy);
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
                if ((fep->init_lambda >= 0) && (n_lambda_terms == 1))
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
                    derivative = gmx::formatString("%s %s = %.4f", dhdl, efpt_singular_names[i], lam);
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
            nsetsbegin = 1; /* for including the expanded ensemble */
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
            if ((fep->init_lambda >= 0) && (n_lambda_terms == 1))
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
                buf += gmx::formatString(
                        "T = %g (%s)", ir->simtempvals->temperatures[s - (nsetsbegin)], unit_temp_K);
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

void EnergyOutput::addDataAtEnergyStep(bool                    bDoDHDL,
                                       bool                    bSum,
                                       double                  time,
                                       real                    tmass,
                                       const gmx_enerdata_t*   enerd,
                                       const t_state*          state,
                                       const t_lambda*         fep,
                                       const t_expanded*       expand,
                                       const matrix            box,
                                       const tensor            svir,
                                       const tensor            fvir,
                                       const tensor            vir,
                                       const tensor            pres,
                                       const gmx_ekindata_t*   ekind,
                                       const rvec              mu_tot,
                                       const gmx::Constraints* constr)
{
    int    j, k, kk, n, gid;
    real   crmsd[2], tmp6[6];
    real   bs[tricl_boxs_nm.size()], vol, dens, pv, enthalpy;
    real   eee[egNR];
    double store_dhdl[efptNR];
    real   store_energy = 0;
    real   tmp;

    /* Do NOT use the box in the state variable, but the separate box provided
     * as an argument. This is because we sometimes need to write the box from
     * the last timestep to match the trajectory frames.
     */
    add_ebin_indexed(ebin_, ie_, gmx::ArrayRef<bool>(bEner_), enerd->term, bSum);
    if (nCrmsd_)
    {
        crmsd[0] = constr->rmsd();
        add_ebin(ebin_, iconrmsd_, nCrmsd_, crmsd, false);
    }
    if (bDynBox_)
    {
        int nboxs;
        if (bTricl_)
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
        vol  = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
        dens = (tmass * AMU) / (vol * NANO * NANO * NANO);
        add_ebin(ebin_, ib_, nboxs, bs, bSum);
        add_ebin(ebin_, ivol_, 1, &vol, bSum);
        add_ebin(ebin_, idens_, 1, &dens, bSum);

        if (bDiagPres_)
        {
            /* This is pV (in kJ/mol).  The pressure is the reference pressure,
               not the instantaneous pressure */
            pv = vol * ref_p_ / PRESFAC;

            add_ebin(ebin_, ipv_, 1, &pv, bSum);
            enthalpy = pv + enerd->term[F_ETOT];
            add_ebin(ebin_, ienthalpy_, 1, &enthalpy, bSum);
        }
    }
    if (bConstrVir_)
    {
        add_ebin(ebin_, isvir_, 9, svir[0], bSum);
        add_ebin(ebin_, ifvir_, 9, fvir[0], bSum);
    }
    if (bPres_)
    {
        add_ebin(ebin_, ivir_, 9, vir[0], bSum);
        add_ebin(ebin_, ipres_, 9, pres[0], bSum);
        tmp = (pres[ZZ][ZZ] - (pres[XX][XX] + pres[YY][YY]) * 0.5) * box[ZZ][ZZ];
        add_ebin(ebin_, isurft_, 1, &tmp, bSum);
    }
    if (epc_ == epcPARRINELLORAHMAN || epc_ == epcMTTK)
    {
        tmp6[0] = state->boxv[XX][XX];
        tmp6[1] = state->boxv[YY][YY];
        tmp6[2] = state->boxv[ZZ][ZZ];
        tmp6[3] = state->boxv[YY][XX];
        tmp6[4] = state->boxv[ZZ][XX];
        tmp6[5] = state->boxv[ZZ][YY];
        add_ebin(ebin_, ipc_, bTricl_ ? 6 : 3, tmp6, bSum);
    }
    if (bMu_)
    {
        add_ebin(ebin_, imu_, 3, mu_tot, bSum);
    }
    if (ekind && ekind->cosacc.cos_accel != 0)
    {
        vol  = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
        dens = (tmass * AMU) / (vol * NANO * NANO * NANO);
        add_ebin(ebin_, ivcos_, 1, &(ekind->cosacc.vcos), bSum);
        /* 1/viscosity, unit 1/(kg m^-1 s^-1) */
        tmp = 1
              / (ekind->cosacc.cos_accel / (ekind->cosacc.vcos * PICO) * dens
                 * gmx::square(box[ZZ][ZZ] * NANO / (2 * M_PI)));
        add_ebin(ebin_, ivisc_, 1, &tmp, bSum);
    }
    if (nE_ > 1)
    {
        n = 0;
        for (int i = 0; (i < nEg_); i++)
        {
            for (j = i; (j < nEg_); j++)
            {
                gid = GID(i, j, nEg_);
                for (k = kk = 0; (k < egNR); k++)
                {
                    if (bEInd_[k])
                    {
                        eee[kk++] = enerd->grpp.ener[k][gid];
                    }
                }
                add_ebin(ebin_, igrp_[n], nEc_, eee, bSum);
                n++;
            }
        }
    }

    if (ekind)
    {
        for (int i = 0; (i < nTC_); i++)
        {
            tmp_r_[i] = ekind->tcstat[i].T;
        }
        add_ebin(ebin_, itemp_, nTC_, tmp_r_, bSum);

        if (etc_ == etcNOSEHOOVER)
        {
            /* whether to print Nose-Hoover chains: */
            if (bPrintNHChains_)
            {
                if (bNHC_trotter_)
                {
                    for (int i = 0; (i < nTC_); i++)
                    {
                        for (j = 0; j < nNHC_; j++)
                        {
                            k                 = i * nNHC_ + j;
                            tmp_r_[2 * k]     = state->nosehoover_xi[k];
                            tmp_r_[2 * k + 1] = state->nosehoover_vxi[k];
                        }
                    }
                    add_ebin(ebin_, itc_, mde_n_, tmp_r_, bSum);

                    if (bMTTK_)
                    {
                        for (int i = 0; (i < nTCP_); i++)
                        {
                            for (j = 0; j < nNHC_; j++)
                            {
                                k                 = i * nNHC_ + j;
                                tmp_r_[2 * k]     = state->nhpres_xi[k];
                                tmp_r_[2 * k + 1] = state->nhpres_vxi[k];
                            }
                        }
                        add_ebin(ebin_, itcb_, mdeb_n_, tmp_r_, bSum);
                    }
                }
                else
                {
                    for (int i = 0; (i < nTC_); i++)
                    {
                        tmp_r_[2 * i]     = state->nosehoover_xi[i];
                        tmp_r_[2 * i + 1] = state->nosehoover_vxi[i];
                    }
                    add_ebin(ebin_, itc_, mde_n_, tmp_r_, bSum);
                }
            }
        }
        else if (etc_ == etcBERENDSEN || etc_ == etcYES || etc_ == etcVRESCALE)
        {
            for (int i = 0; (i < nTC_); i++)
            {
                tmp_r_[i] = ekind->tcstat[i].lambda;
            }
            add_ebin(ebin_, itc_, nTC_, tmp_r_, bSum);
        }
    }

    if (ekind && nU_ > 1)
    {
        for (int i = 0; (i < nU_); i++)
        {
            copy_rvec(ekind->grpstat[i].u, tmp_v_[i]);
        }
        add_ebin(ebin_, iu_, 3 * nU_, tmp_v_[0], bSum);
    }

    ebin_increase_count(1, ebin_, bSum);

    // BAR + thermodynamic integration values
    if ((fp_dhdl_ || dhc_) && bDoDHDL)
    {
        for (gmx::index i = 0; i < static_cast<gmx::index>(enerd->enerpart_lambda.size()) - 1; i++)
        {
            /* zero for simulated tempering */
            dE_[i] = enerd->enerpart_lambda[i + 1] - enerd->enerpart_lambda[0];
            if (numTemperatures_ > 0)
            {
                GMX_RELEASE_ASSERT(numTemperatures_ > state->fep_state,
                                   "Number of lambdas in state is bigger then in input record");
                GMX_RELEASE_ASSERT(
                        numTemperatures_ >= static_cast<gmx::index>(enerd->enerpart_lambda.size()) - 1,
                        "Number of lambdas in energy data is bigger then in input record");
                /* MRS: is this right, given the way we have defined the exchange probabilities? */
                /* is this even useful to have at all? */
                dE_[i] += (temperatures_[i] / temperatures_[state->fep_state] - 1.0) * enerd->term[F_EKIN];
            }
        }

        if (fp_dhdl_)
        {
            fprintf(fp_dhdl_, "%.4f", time);
            /* the current free energy state */

            /* print the current state if we are doing expanded ensemble */
            if (expand->elmcmove > elmcmoveNO)
            {
                fprintf(fp_dhdl_, " %4d", state->fep_state);
            }
            /* total energy (for if the temperature changes */

            if (fep->edHdLPrintEnergy != edHdLPrintEnergyNO)
            {
                switch (fep->edHdLPrintEnergy)
                {
                    case edHdLPrintEnergyPOTENTIAL: store_energy = enerd->term[F_EPOT]; break;
                    case edHdLPrintEnergyTOTAL:
                    case edHdLPrintEnergyYES:
                    default: store_energy = enerd->term[F_ETOT];
                }
                fprintf(fp_dhdl_, " %#.8g", store_energy);
            }

            if (fep->dhdl_derivatives == edhdlderivativesYES)
            {
                for (int i = 0; i < efptNR; i++)
                {
                    if (fep->separate_dvdl[i])
                    {
                        /* assumes F_DVDL is first */
                        fprintf(fp_dhdl_, " %#.8g", enerd->term[F_DVDL + i]);
                    }
                }
            }
            for (int i = fep->lambda_start_n; i < fep->lambda_stop_n; i++)
            {
                fprintf(fp_dhdl_, " %#.8g", dE_[i]);
            }
            if (bDynBox_ && bDiagPres_ && (epc_ != epcNO) && !enerd->enerpart_lambda.empty()
                && (fep->init_lambda < 0))
            {
                fprintf(fp_dhdl_, " %#.8g", pv); /* PV term only needed when
                                                         there are alternate state
                                                         lambda and we're not in
                                                         compatibility mode */
            }
            fprintf(fp_dhdl_, "\n");
            /* and the binary free energy output */
        }
        if (dhc_ && bDoDHDL)
        {
            int idhdl = 0;
            for (int i = 0; i < efptNR; i++)
            {
                if (fep->separate_dvdl[i])
                {
                    /* assumes F_DVDL is first */
                    store_dhdl[idhdl] = enerd->term[F_DVDL + i];
                    idhdl += 1;
                }
            }
            store_energy = enerd->term[F_ETOT];
            /* store_dh is dE */
            mde_delta_h_coll_add_dh(dhc_, static_cast<double>(state->fep_state), store_energy, pv,
                                    store_dhdl, dE_ + fep->lambda_start_n, time);
        }
    }
}

void EnergyOutput::recordNonEnergyStep()
{
    ebin_increase_count(1, ebin_, false);
}

void EnergyOutput::printHeader(FILE* log, int64_t steps, double time)
{
    char buf[22];

    fprintf(log,
            "   %12s   %12s\n"
            "   %12s   %12.5f\n\n",
            "Step", "Time", gmx_step_str(steps, buf), time);
}

void EnergyOutput::printStepToEnergyFile(ener_file* fp_ene,
                                         bool       bEne,
                                         bool       bDR,
                                         bool       bOR,
                                         FILE*      log,
                                         int64_t    step,
                                         double     time,
                                         t_fcdata*  fcd,
                                         gmx::Awh*  awh)
{
    t_enxframe fr;
    init_enxframe(&fr);
    fr.t       = time;
    fr.step    = step;
    fr.nsteps  = ebin_->nsteps;
    fr.dt      = delta_t_;
    fr.nsum    = ebin_->nsum;
    fr.nre     = (bEne) ? ebin_->nener : 0;
    fr.ener    = ebin_->e;
    int ndisre = bDR ? fcd->disres.npair : 0;
    /* these are for the old-style blocks (1 subblock, only reals), because
       there can be only one per ID for these */
    int   nr[enxNR];
    int   id[enxNR];
    real* block[enxNR];
    /* Optional additional old-style (real-only) blocks. */
    for (int i = 0; i < enxNR; i++)
    {
        nr[i] = 0;
    }

    if (bOR && fcd->orires.nr > 0)
    {
        diagonalize_orires_tensors(&(fcd->orires));
        nr[enxOR]     = fcd->orires.nr;
        block[enxOR]  = fcd->orires.otav;
        id[enxOR]     = enxOR;
        nr[enxORI]    = (fcd->orires.oinsl != fcd->orires.otav) ? fcd->orires.nr : 0;
        block[enxORI] = fcd->orires.oinsl;
        id[enxORI]    = enxORI;
        nr[enxORT]    = fcd->orires.nex * 12;
        block[enxORT] = fcd->orires.eig;
        id[enxORT]    = enxORT;
    }

    /* whether we are going to write anything out: */
    if (fr.nre || ndisre || nr[enxOR] || nr[enxORI])
    {
        /* the old-style blocks go first */
        fr.nblock = 0;
        for (int i = 0; i < enxNR; i++)
        {
            if (nr[i] > 0)
            {
                fr.nblock = i + 1;
            }
        }
        add_blocks_enxframe(&fr, fr.nblock);
        for (int b = 0; b < fr.nblock; b++)
        {
            add_subblocks_enxblock(&(fr.block[b]), 1);
            fr.block[b].id        = id[b];
            fr.block[b].sub[0].nr = nr[b];
#if !GMX_DOUBLE
            fr.block[b].sub[0].type = xdr_datatype_float;
            fr.block[b].sub[0].fval = block[b];
#else
            fr.block[b].sub[0].type  = xdr_datatype_double;
            fr.block[b].sub[0].dval  = block[b];
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
        if (dhc_)
        {
            mde_delta_h_coll_handle_block(dhc_, &fr, fr.nblock);
        }

        /* we can now free & reset the data in the blocks */
        if (dhc_)
        {
            mde_delta_h_coll_reset(dhc_);
        }

        /* AWH bias blocks. */
        if (awh != nullptr) // TODO: add boolean flag.
        {
            awh->writeToEnergyFrame(step, &fr);
        }

        /* do the actual I/O */
        do_enx(fp_ene, &fr);
        if (fr.nre)
        {
            /* We have stored the sums, so reset the sum history */
            reset_ebin_sums(ebin_);
        }
    }
    free_enxframe(&fr);
    if (log)
    {
        if (bOR && fcd->orires.nr > 0)
        {
            print_orires_log(log, &(fcd->orires));
        }

        fprintf(log, "   Energies (%s)\n", unit_energy);
        pr_ebin(log, ebin_, ie_, f_nre_ + nCrmsd_, 5, eprNORMAL, true);
        fprintf(log, "\n");
    }
}

void EnergyOutput::printAnnealingTemperatures(FILE* log, SimulationGroups* groups, t_grpopts* opts)
{
    if (log)
    {
        if (opts)
        {
            for (int i = 0; i < opts->ngtc; i++)
            {
                if (opts->annealing[i] != eannNO)
                {
                    fprintf(log, "Current ref_t for group %s: %8.1f\n",
                            *(groups->groupNames[groups->groups[SimulationAtomGroupType::TemperatureCoupling][i]]),
                            opts->ref_t[i]);
                }
            }
            fprintf(log, "\n");
        }
    }
}

void EnergyOutput::printAverages(FILE* log, const SimulationGroups* groups)
{
    if (ebin_->nsum_sim <= 0)
    {
        if (log)
        {
            fprintf(log, "Not enough data recorded to report energy averages\n");
        }
        return;
    }
    if (log)
    {

        char buf1[22], buf2[22];

        fprintf(log, "\t<======  ###############  ==>\n");
        fprintf(log, "\t<====  A V E R A G E S  ====>\n");
        fprintf(log, "\t<==  ###############  ======>\n\n");

        fprintf(log, "\tStatistics over %s steps using %s frames\n",
                gmx_step_str(ebin_->nsteps_sim, buf1), gmx_step_str(ebin_->nsum_sim, buf2));
        fprintf(log, "\n");

        fprintf(log, "   Energies (%s)\n", unit_energy);
        pr_ebin(log, ebin_, ie_, f_nre_ + nCrmsd_, 5, eprAVER, true);
        fprintf(log, "\n");

        if (bDynBox_)
        {
            pr_ebin(log, ebin_, ib_, bTricl_ ? tricl_boxs_nm.size() : boxs_nm.size(), 5, eprAVER, true);
            fprintf(log, "\n");
        }
        if (bConstrVir_)
        {
            fprintf(log, "   Constraint Virial (%s)\n", unit_energy);
            pr_ebin(log, ebin_, isvir_, 9, 3, eprAVER, false);
            fprintf(log, "\n");
            fprintf(log, "   Force Virial (%s)\n", unit_energy);
            pr_ebin(log, ebin_, ifvir_, 9, 3, eprAVER, false);
            fprintf(log, "\n");
        }
        if (bPres_)
        {
            fprintf(log, "   Total Virial (%s)\n", unit_energy);
            pr_ebin(log, ebin_, ivir_, 9, 3, eprAVER, false);
            fprintf(log, "\n");
            fprintf(log, "   Pressure (%s)\n", unit_pres_bar);
            pr_ebin(log, ebin_, ipres_, 9, 3, eprAVER, false);
            fprintf(log, "\n");
        }
        if (bMu_)
        {
            fprintf(log, "   Total Dipole (%s)\n", unit_dipole_D);
            pr_ebin(log, ebin_, imu_, 3, 3, eprAVER, false);
            fprintf(log, "\n");
        }

        if (nE_ > 1)
        {
            int padding = 8 - strlen(unit_energy);
            fprintf(log, "%*sEpot (%s)   ", padding, "", unit_energy);
            for (int i = 0; (i < egNR); i++)
            {
                if (bEInd_[i])
                {
                    fprintf(log, "%12s   ", egrp_nm[i]);
                }
            }
            fprintf(log, "\n");

            int n = 0;
            for (int i = 0; (i < nEg_); i++)
            {
                int ni = groups->groups[SimulationAtomGroupType::EnergyOutput][i];
                for (int j = i; (j < nEg_); j++)
                {
                    int nj = groups->groups[SimulationAtomGroupType::EnergyOutput][j];
                    int padding =
                            14 - (strlen(*(groups->groupNames[ni])) + strlen(*(groups->groupNames[nj])));
                    fprintf(log, "%*s%s-%s", padding, "", *(groups->groupNames[ni]),
                            *(groups->groupNames[nj]));
                    pr_ebin(log, ebin_, igrp_[n], nEc_, nEc_, eprAVER, false);
                    n++;
                }
            }
            fprintf(log, "\n");
        }
        if (nTC_ > 1)
        {
            pr_ebin(log, ebin_, itemp_, nTC_, 4, eprAVER, true);
            fprintf(log, "\n");
        }
        if (nU_ > 1)
        {
            fprintf(log, "%15s   %12s   %12s   %12s\n", "Group", "Ux", "Uy", "Uz");
            for (int i = 0; (i < nU_); i++)
            {
                int ni = groups->groups[SimulationAtomGroupType::Acceleration][i];
                fprintf(log, "%15s", *groups->groupNames[ni]);
                pr_ebin(log, ebin_, iu_ + 3 * i, 3, 3, eprAVER, false);
            }
            fprintf(log, "\n");
        }
    }
}

void EnergyOutput::fillEnergyHistory(energyhistory_t* enerhist) const
{
    const t_ebin* const ebin = ebin_;

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
    if (dhc_)
    {
        mde_delta_h_coll_update_energyhistory(dhc_, enerhist);
    }
}

void EnergyOutput::restoreFromEnergyHistory(const energyhistory_t& enerhist)
{
    unsigned int nener = static_cast<unsigned int>(ebin_->nener);

    if ((enerhist.nsum > 0 && nener != enerhist.ener_sum.size())
        || (enerhist.nsum_sim > 0 && nener != enerhist.ener_sum_sim.size()))
    {
        gmx_fatal(FARGS,
                  "Mismatch between number of energies in run input (%u) and checkpoint file (%zu "
                  "or %zu).",
                  nener, enerhist.ener_sum.size(), enerhist.ener_sum_sim.size());
    }

    ebin_->nsteps     = enerhist.nsteps;
    ebin_->nsum       = enerhist.nsum;
    ebin_->nsteps_sim = enerhist.nsteps_sim;
    ebin_->nsum_sim   = enerhist.nsum_sim;

    for (int i = 0; i < ebin_->nener; i++)
    {
        ebin_->e[i].eav      = (enerhist.nsum > 0 ? enerhist.ener_ave[i] : 0);
        ebin_->e[i].esum     = (enerhist.nsum > 0 ? enerhist.ener_sum[i] : 0);
        ebin_->e_sim[i].esum = (enerhist.nsum_sim > 0 ? enerhist.ener_sum_sim[i] : 0);
    }
    if (dhc_)
    {
        mde_delta_h_coll_restore_energyhistory(dhc_, enerhist.deltaHForeignLambdas.get());
    }
}

int EnergyOutput::numEnergyTerms() const
{
    return ebin_->nener;
}

} // namespace gmx
