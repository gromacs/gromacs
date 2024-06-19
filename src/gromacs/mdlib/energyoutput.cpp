/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xdr_datatype.h"
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
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "energydrifttracker.h"

//! Labels for energy file quantities
//! \{
static const char* const conrmsd_nm[] = { "Constr. rmsd", "Constr.2 rmsd" };

static constexpr std::array<const char*, 3> boxs_nm = { "Box-X", "Box-Y", "Box-Z" };

static constexpr std::array<const char*, 6> tricl_boxs_nm = { "Box-XX", "Box-YY", "Box-ZZ",
                                                              "Box-YX", "Box-ZX", "Box-ZY" };

static const char* const vol_nm[] = { "Volume" };

static const char* const dens_nm[] = { "Density" };

const char* const pvEnergyFieldName = "pV";

const char* const enthalpyEnergyFieldName = "Enthalpy";

const std::array<const char*, 9> virialEnergyFieldNames = { "Vir-XX", "Vir-XY", "Vir-XZ",
                                                            "Vir-YX", "Vir-YY", "Vir-YZ",
                                                            "Vir-ZX", "Vir-ZY", "Vir-ZZ" };

static constexpr std::array<const char*, 6> boxvel_nm = { "Box-Vel-XX", "Box-Vel-YY", "Box-Vel-ZZ",
                                                          "Box-Vel-YX", "Box-Vel-ZX", "Box-Vel-ZY" };

const char* enumValueToString(NonBondedEnergyTerms enumValue)
{
    static constexpr gmx::EnumerationArray<NonBondedEnergyTerms, const char*> nonBondedEnergyTermTypeNames = {
        "Coul-SR", "LJ-SR", "Buck-SR", "Coul-14", "LJ-14"
    };
    return nonBondedEnergyTermTypeNames[enumValue];
}

//! \}

static bool haveFepLambdaMoves(const t_inputrec& inputrec)
{
    return (inputrec.bExpanded && inputrec.expandedvals->elmcmove > LambdaMoveCalculation::No)
           || (inputrec.efep != FreeEnergyPerturbationType::No && inputrec.bDoAwh
               && awhHasFepLambdaDimension(*inputrec.awhParams));
}

namespace gmx
{

/*! \brief Energy output class
 *
 * This is the collection of energy averages collected during mdrun, and to
 * be written out to the .edr file.
 *
 * \todo Use more std containers.
 * \todo Write free-energy output also to energy file (after adding more tests)
 */
EnergyOutput::EnergyOutput(ener_file*                fp_ene,
                           const gmx_mtop_t&         mtop,
                           const t_inputrec&         inputrec,
                           const pull_t*             pull_work,
                           FILE*                     fp_dhdl,
                           bool                      isRerun,
                           const StartingBehavior    startingBehavior,
                           const bool                simulationsShareState,
                           const MDModulesNotifiers& mdModulesNotifiers) :
    haveFepLambdaMoves_(haveFepLambdaMoves(inputrec))
{
    const char*        ener_nm[F_NRE];
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
    int                     i, j, ni, nj, n, ncon, nset;
    bool                    bBHAM, b14;

    if (EI_DYNAMICS(inputrec.eI))
    {
        delta_t_ = inputrec.delta_t;
    }
    else
    {
        delta_t_ = 0;
    }

    groups = &mtop.groups;

    bBHAM = (mtop.ffparams.numTypes() > 0) && (mtop.ffparams.functype[0] == F_BHAM);
    b14   = (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 || gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0);

    ncon         = gmx_mtop_ftype_count(mtop, F_CONSTR);
    nset         = gmx_mtop_ftype_count(mtop, F_SETTLE);
    bool bConstr = (ncon > 0 || nset > 0) && !isRerun;
    nCrmsd_      = 0;
    if (bConstr)
    {
        if (ncon > 0 && inputrec.eConstrAlg == ConstraintAlgorithm::Lincs)
        {
            nCrmsd_ = 1;
        }
    }
    else
    {
        nCrmsd_ = 0;
    }

    /* Energy monitoring */
    for (auto& term : bEInd_)
    {
        term = false;
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
        bEner_[F_EKIN] = EI_DYNAMICS(inputrec.eI);
        bEner_[F_ETOT] = EI_DYNAMICS(inputrec.eI);
        bEner_[F_TEMP] = EI_DYNAMICS(inputrec.eI);

        bEner_[F_ECONSERVED] = integratorHasConservedEnergyQuantity(&inputrec);
        bEner_[F_PDISPCORR]  = (inputrec.eDispCorr != DispersionCorrectionType::No);
        bEner_[F_PRES]       = true;
    }

    bEner_[F_LJ]   = !bBHAM;
    bEner_[F_BHAM] = bBHAM;
    bEner_[F_RF_EXCL] = (usingRF(inputrec.coulombtype) && inputrec.cutoff_scheme == CutoffScheme::Group);
    bEner_[F_COUL_RECIP]   = usingFullElectrostatics(inputrec.coulombtype);
    bEner_[F_LJ_RECIP]     = usingLJPme(inputrec.vdwtype);
    bEner_[F_LJ14]         = b14;
    bEner_[F_COUL14]       = b14;
    bEner_[F_LJC14_Q]      = false;
    bEner_[F_LJC_PAIRS_NB] = false;


    bEner_[F_DVDL_COUL] = (inputrec.efep != FreeEnergyPerturbationType::No)
                          && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Coul];
    bEner_[F_DVDL_VDW] = (inputrec.efep != FreeEnergyPerturbationType::No)
                         && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Vdw];
    bEner_[F_DVDL_BONDED] = (inputrec.efep != FreeEnergyPerturbationType::No)
                            && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Bonded];
    bEner_[F_DVDL_RESTRAINT] =
            (inputrec.efep != FreeEnergyPerturbationType::No)
            && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Restraint];
    bEner_[F_DKDL] = (inputrec.efep != FreeEnergyPerturbationType::No)
                     && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Mass];
    bEner_[F_DVDL] = (inputrec.efep != FreeEnergyPerturbationType::No)
                     && inputrec.fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep];

    bEner_[F_CONSTR]   = false;
    bEner_[F_CONSTRNC] = false;
    bEner_[F_SETTLE]   = false;

    bEner_[F_COUL_SR] = true;
    bEner_[F_EPOT]    = true;

    bEner_[F_DISPCORR]   = (inputrec.eDispCorr != DispersionCorrectionType::No);
    bEner_[F_DISRESVIOL] = (gmx_mtop_ftype_count(mtop, F_DISRES) > 0);
    bEner_[F_ORIRESDEV]  = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
    bEner_[F_COM_PULL]   = ((inputrec.bPull && pull_have_potential(*pull_work)) || inputrec.bRot);

    // Check MDModules for any energy output
    MDModulesEnergyOutputToDensityFittingRequestChecker mdModulesAddOutputToDensityFittingFieldRequest;
    mdModulesNotifiers.simulationSetupNotifier_.notify(&mdModulesAddOutputToDensityFittingFieldRequest);

    bEner_[F_DENSITYFITTING] = mdModulesAddOutputToDensityFittingFieldRequest.energyOutputToDensityFitting_;

    MDModulesEnergyOutputToQMMMRequestChecker mdModulesAddOutputToQMMMFieldRequest;
    mdModulesNotifiers.simulationSetupNotifier_.notify(&mdModulesAddOutputToQMMMFieldRequest);

    bEner_[F_EQM] = mdModulesAddOutputToQMMMFieldRequest.energyOutputToQMMM_;

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

    epc_       = isRerun ? PressureCoupling::No : inputrec.pressureCouplingOptions.epc;
    bDiagPres_ = !TRICLINIC(inputrec.pressureCouplingOptions.ref_p) && !isRerun;
    ref_p_     = (inputrec.pressureCouplingOptions.ref_p[XX][XX]
              + inputrec.pressureCouplingOptions.ref_p[YY][YY]
              + inputrec.pressureCouplingOptions.ref_p[ZZ][ZZ])
             / DIM;
    bTricl_  = TRICLINIC(inputrec.pressureCouplingOptions.compress) || TRICLINIC(inputrec.deform);
    bDynBox_ = inputrecDynamicBox(&inputrec);
    etc_     = isRerun ? TemperatureCoupling::No : inputrec.etc;
    bNHC_trotter_   = inputrecNvtTrotter(&inputrec) && !isRerun;
    bPrintNHChains_ = inputrec.bPrintNHChains && !isRerun;
    bMTTK_          = (inputrecNptTrotter(&inputrec) || inputrecNphTrotter(&inputrec)) && !isRerun;
    bMu_            = inputrecNeedMutot(&inputrec);
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
        ib_    = get_ebin_space(ebin_,
                             bTricl_ ? tricl_boxs_nm.size() : boxs_nm.size(),
                             bTricl_ ? tricl_boxs_nm.data() : boxs_nm.data(),
                             unit_length);
        ivol_  = get_ebin_space(ebin_, 1, vol_nm, unit_volume);
        idens_ = get_ebin_space(ebin_, 1, dens_nm, unit_density_SI);
        if (bDiagPres_)
        {
            ipv_       = get_ebin_space(ebin_, 1, &pvEnergyFieldName, unit_energy);
            ienthalpy_ = get_ebin_space(ebin_, 1, &enthalpyEnergyFieldName, unit_energy);
        }
    }
    if (bPres_)
    {
        ivir_ = get_ebin_space(
                ebin_, virialEnergyFieldNames.size(), virialEnergyFieldNames.data(), unit_energy);
        ipres_  = get_ebin_space(ebin_, asize(pres_nm), pres_nm, unit_pres_bar);
        isurft_ = get_ebin_space(ebin_, asize(surft_nm), surft_nm, unit_surft_bar);
    }
    if (epc_ == PressureCoupling::ParrinelloRahman || epc_ == PressureCoupling::Mttk)
    {
        ipc_ = get_ebin_space(ebin_, bTricl_ ? boxvel_nm.size() : DIM, boxvel_nm.data(), unit_vel);
    }
    if (bMu_)
    {
        imu_ = get_ebin_space(ebin_, asize(mu_nm), mu_nm, unit_dipole_D);
    }
    if (inputrec.cos_accel != 0)
    {
        ivcos_ = get_ebin_space(ebin_, asize(vcos_nm), vcos_nm, unit_vel);
        ivisc_ = get_ebin_space(ebin_, asize(visc_nm), visc_nm, unit_invvisc_SI);
    }

    /* Energy monitoring */
    for (auto& term : bEInd_)
    {
        term = false;
    }
    bEInd_[NonBondedEnergyTerms::CoulombSR] = true;
    bEInd_[NonBondedEnergyTerms::LJSR]      = true;

    if (bBHAM)
    {
        bEInd_[NonBondedEnergyTerms::LJSR]         = false;
        bEInd_[NonBondedEnergyTerms::BuckinghamSR] = true;
    }
    if (b14)
    {
        bEInd_[NonBondedEnergyTerms::LJ14]      = true;
        bEInd_[NonBondedEnergyTerms::Coulomb14] = true;
    }
    nEc_ = 0;
    for (auto term : bEInd_)
    {
        if (term)
        {
            nEc_++;
        }
    }
    n    = groups->groups[SimulationAtomGroupType::EnergyOutput].size();
    nEg_ = n;
    nE_  = (n * (n + 1)) / 2;

    igrp_.resize(nE_);
    if (nE_ > 1)
    {
        n = 0;
        snew(gnm, nEc_);
        for (int k = 0; (k < nEc_); k++)
        {
            snew(gnm[k], STRLEN);
        }
        for (i = 0; (i < gmx::ssize(groups->groups[SimulationAtomGroupType::EnergyOutput])); i++)
        {
            ni = groups->groups[SimulationAtomGroupType::EnergyOutput][i];
            for (j = i; (j < gmx::ssize(groups->groups[SimulationAtomGroupType::EnergyOutput])); j++)
            {
                nj    = groups->groups[SimulationAtomGroupType::EnergyOutput][j];
                int k = 0;
                for (auto key : keysOf(bEInd_))
                {
                    if (bEInd_[key])
                    {
                        sprintf(gnm[k],
                                "%s:%s-%s",
                                enumValueToString(key),
                                *(groups->groupNames[ni]),
                                *(groups->groupNames[nj]));
                        k++;
                    }
                }
                igrp_[n] = get_ebin_space(ebin_, nEc_, gnm, unit_energy);
                n++;
            }
        }
        for (int k = 0; (k < nEc_); k++)
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
    nNHC_ = inputrec.opts.nhchainlength; /* shorthand for number of NH chains */
    if (bMTTK_)
    {
        nTCP_ = 1; /* assume only one possible coupling system for barostat
                           for now */
    }
    else
    {
        nTCP_ = 0;
    }
    if (etc_ == TemperatureCoupling::NoseHoover)
    {
        if (bNHC_trotter_)
        {
            mde_n_ = 2 * nNHC_ * nTC_;
        }
        else
        {
            mde_n_ = 2 * nTC_;
        }
        if (epc_ == PressureCoupling::Mttk)
        {
            mdeb_n_ = 2 * nNHC_ * nTCP_;
        }
    }
    else
    {
        mde_n_  = nTC_;
        mdeb_n_ = 0;
    }

    tmp_r_.resize(mde_n_);
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
    if (etc_ == TemperatureCoupling::NoseHoover)
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
    else if (etc_ == TemperatureCoupling::Berendsen || etc_ == TemperatureCoupling::Yes
             || etc_ == TemperatureCoupling::VRescale)
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

    /* Note that fp_ene should be valid on the main rank and null otherwise */
    if (fp_ene != nullptr && startingBehavior != StartingBehavior::RestartWithAppending)
    {
        do_enxnms(fp_ene, &ebin_->nener, &ebin_->enm);
    }

    /* check whether we're going to write dh histograms */
    dhc_ = nullptr;
    if (inputrec.fepvals->separate_dhdl_file == SeparateDhdlFile::No)
    {
        /* Currently dh histograms are only written with dynamics */
        if (EI_DYNAMICS(inputrec.eI))
        {
            dhc_ = std::make_unique<t_mde_delta_h_coll>(inputrec);
        }
        fp_dhdl_ = nullptr;
        dE_.resize(inputrec.fepvals->n_lambda);
    }
    else
    {
        fp_dhdl_ = fp_dhdl;
        dE_.resize(inputrec.fepvals->n_lambda);
    }
    if (inputrec.bSimTemp)
    {
        temperatures_ = inputrec.simtempvals->temperatures;
    }

    if (EI_MD(inputrec.eI) && !simulationsShareState)
    {
        conservedEnergyTracker_ = std::make_unique<EnergyDriftTracker>(mtop.natoms);
    }
}

EnergyOutput::~EnergyOutput()
{
    done_ebin(ebin_);
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
    int k    = 0;
    int Nsep = 0;

    for (auto j : keysOf(fep->separate_dvdl))
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
    for (auto j : keysOf(fep->separate_dvdl))
    {
        if (fep->separate_dvdl[j])
        {
            if (!get_names)
            {
                if (get_native_lambda && fep->init_lambda_without_states >= 0)
                {
                    str += sprintf(str, "%.4f", fep->init_lambda_without_states);
                }
                else
                {
                    str += sprintf(str, "%.4f", fep->all_lambda[j][i]);
                }
            }
            else
            {
                str += sprintf(str, "%s", enumValueToStringSingular(j));
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
    int       i, nsets, nsets_de, nsetsbegin;
    int       n_lambda_terms = 0;
    t_lambda* fep            = ir->fepvals.get(); /* for simplicity */
    char      lambda_vec_str[STRLEN], lambda_name_str[STRLEN];

    int  nsets_dhdl = 0;
    int  s          = 0;
    int  nsetsextend;
    bool write_pV = false;

    /* count the number of different lambda terms */
    for (auto i : keysOf(fep->separate_dvdl))
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
        label_y = gmx::formatString(
                "%s and %s (%s %s)", dhdl, deltag, unit_energy, "[\\8l\\4]\\S-1\\N");
    }
    fp = gmx_fio_fopen(filename, "w+");
    xvgr_header(fp, title.c_str(), label_x, label_y, exvggtXNY, oenv);

    std::string buf;
    if (!(ir->bSimTemp) && haveConstantEnsembleTemperature(*ir))
    {
        buf = gmx::formatString("T = %g (K) ", constantEnsembleTemperature(*ir));
    }
    if ((ir->efep != FreeEnergyPerturbationType::SlowGrowth)
        && (ir->efep != FreeEnergyPerturbationType::Expanded)
        && !(ir->bDoAwh && awhHasFepLambdaDimension(*ir->awhParams)))
    {
        if ((fep->init_lambda_without_states >= 0) && (n_lambda_terms == 1))
        {
            /* compatibility output */
            buf += gmx::formatString("%s = %.4f", lambda, fep->init_lambda_without_states);
        }
        else
        {
            print_lambda_vector(fep, fep->init_fep_state, true, false, lambda_vec_str);
            print_lambda_vector(fep, fep->init_fep_state, true, true, lambda_name_str);
            buf += gmx::formatString(
                    "%s %d: %s = %s", lambdastate, fep->init_fep_state, lambda_name_str, lambda_vec_str);
        }
    }
    xvgr_subtitle(fp, buf.c_str(), oenv);


    nsets_dhdl = 0;
    if (fep->dhdl_derivatives == DhDlDerivativeCalculation::Yes)
    {
        nsets_dhdl = n_lambda_terms;
    }
    /* count the number of delta_g states */
    nsets_de = fep->lambda_stop_n - fep->lambda_start_n;

    nsets = nsets_dhdl + nsets_de; /* dhdl + fep differences */

    if (haveFepLambdaMoves(*ir))
    {
        nsets += 1; /*add fep state for expanded ensemble */
    }

    if (fep->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
    {
        nsets += 1; /* add energy to the dhdl as well */
    }

    nsetsextend = nsets;
    if ((ir->pressureCouplingOptions.epc != PressureCoupling::No) && (fep->n_lambda > 0)
        && (fep->init_lambda_without_states < 0))
    {
        nsetsextend += 1; /* for PV term, other terms possible if required for
                             the reduced potential (only needed with foreign
                             lambda, and only output when init_lambda is not
                             set in order to maintain compatibility of the
                             dhdl.xvg file) */
        write_pV = true;
    }
    std::vector<std::string> setname(nsetsextend);

    if (haveFepLambdaMoves(*ir))
    {
        /* state for the fep_vals, if we have alchemical sampling */
        setname[s++] = "Thermodynamic state";
    }

    if (fep->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
    {
        std::string energy;
        switch (fep->edHdLPrintEnergy)
        {
            case FreeEnergyPrintEnergy::Potential:
                energy = gmx::formatString("%s (%s)", "Potential Energy", unit_energy);
                break;
            case FreeEnergyPrintEnergy::Total:
            case FreeEnergyPrintEnergy::Yes:
            default: energy = gmx::formatString("%s (%s)", "Total Energy", unit_energy);
        }
        setname[s++] = energy;
    }

    if (fep->dhdl_derivatives == DhDlDerivativeCalculation::Yes)
    {
        for (auto i : keysOf(fep->separate_dvdl))
        {
            if (fep->separate_dvdl[i])
            {
                std::string derivative;
                if ((fep->init_lambda_without_states >= 0) && (n_lambda_terms == 1))
                {
                    /* compatibility output */
                    derivative = gmx::formatString(
                            "%s %s %.4f", dhdl, lambda, fep->init_lambda_without_states);
                }
                else
                {
                    const double lam = fep->initialLambda(i);

                    derivative = gmx::formatString("%s %s = %.4f", dhdl, enumValueToStringSingular(i), lam);
                }
                setname[s++] = derivative;
            }
        }
    }

    if (fep->n_lambda > 0)
    {
        /* gmx bar has to determine the lambda values used in this simulation
         * from this xvg legend.
         */

        if (haveFepLambdaMoves(*ir))
        {
            nsetsbegin = 1; /* for including the expanded ensemble */
        }
        else
        {
            nsetsbegin = 0;
        }

        if (fep->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
        {
            nsetsbegin += 1;
        }
        nsetsbegin += nsets_dhdl;

        for (i = fep->lambda_start_n; i < fep->lambda_stop_n; i++)
        {
            print_lambda_vector(fep, i, false, false, lambda_vec_str);
            std::string buf;
            if ((fep->init_lambda_without_states >= 0) && (n_lambda_terms == 1))
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
                                       const t_lambda*         fep,
                                       const matrix            box,
                                       PTCouplingArrays        ptCouplingArrays,
                                       int                     fep_state,
                                       const tensor            vir,
                                       const tensor            pres,
                                       const gmx_ekindata_t*   ekind,
                                       const rvec              mu_tot,
                                       const gmx::Constraints* constr)
{
    int  j, k, kk, n, gid;
    real crmsd[2], tmp6[6];
    real bs[tricl_boxs_nm.size()], vol, dens, enthalpy;
    real eee[static_cast<int>(NonBondedEnergyTerms::Count)];
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, double> store_dhdl;
    real                                                              store_energy = 0;
    real                                                              tmp;
    real pv = 0.0; // static analyzer warns about uninitialized variable warnings here.

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
        dens = (tmass * gmx::c_amu) / (vol * gmx::c_nano * gmx::c_nano * gmx::c_nano);
        add_ebin(ebin_, ib_, nboxs, bs, bSum);
        add_ebin(ebin_, ivol_, 1, &vol, bSum);
        add_ebin(ebin_, idens_, 1, &dens, bSum);

        if (bDiagPres_)
        {
            /* This is pV (in kJ/mol).  The pressure is the reference pressure,
               not the instantaneous pressure */
            pv = vol * ref_p_ / gmx::c_presfac;

            add_ebin(ebin_, ipv_, 1, &pv, bSum);
            enthalpy = pv + enerd->term[F_ETOT];
            add_ebin(ebin_, ienthalpy_, 1, &enthalpy, bSum);
        }
    }
    if (bPres_)
    {
        add_ebin(ebin_, ivir_, 9, vir[0], bSum);
        add_ebin(ebin_, ipres_, 9, pres[0], bSum);
        tmp = (pres[ZZ][ZZ] - (pres[XX][XX] + pres[YY][YY]) * 0.5) * box[ZZ][ZZ];
        add_ebin(ebin_, isurft_, 1, &tmp, bSum);
    }
    if (epc_ == PressureCoupling::ParrinelloRahman || epc_ == PressureCoupling::Mttk)
    {
        tmp6[0] = ptCouplingArrays.boxv[XX][XX];
        tmp6[1] = ptCouplingArrays.boxv[YY][YY];
        tmp6[2] = ptCouplingArrays.boxv[ZZ][ZZ];
        tmp6[3] = ptCouplingArrays.boxv[YY][XX];
        tmp6[4] = ptCouplingArrays.boxv[ZZ][XX];
        tmp6[5] = ptCouplingArrays.boxv[ZZ][YY];
        add_ebin(ebin_, ipc_, bTricl_ ? 6 : 3, tmp6, bSum);
    }
    if (bMu_)
    {
        add_ebin(ebin_, imu_, 3, mu_tot, bSum);
    }
    if (ekind && ekind->cosacc.cos_accel != 0)
    {
        vol  = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
        dens = (tmass * gmx::c_amu) / (vol * gmx::c_nano * gmx::c_nano * gmx::c_nano);
        add_ebin(ebin_, ivcos_, 1, &(ekind->cosacc.vcos), bSum);
        /* 1/viscosity, unit 1/(kg m^-1 s^-1) */
        tmp = 1
              / (ekind->cosacc.cos_accel / (ekind->cosacc.vcos * gmx::c_pico) * dens
                 * gmx::square(box[ZZ][ZZ] * gmx::c_nano / (2 * M_PI)));
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
                for (k = kk = 0; (k < static_cast<int>(NonBondedEnergyTerms::Count)); k++)
                {
                    if (bEInd_[k])
                    {
                        eee[kk++] = enerd->grpp.energyGroupPairTerms[k][gid];
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
        add_ebin(ebin_, itemp_, nTC_, tmp_r_.data(), bSum);

        if (etc_ == TemperatureCoupling::NoseHoover)
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
                            tmp_r_[2 * k]     = ptCouplingArrays.nosehoover_xi[k];
                            tmp_r_[2 * k + 1] = ptCouplingArrays.nosehoover_vxi[k];
                        }
                    }
                    add_ebin(ebin_, itc_, mde_n_, tmp_r_.data(), bSum);

                    if (bMTTK_)
                    {
                        for (int i = 0; (i < nTCP_); i++)
                        {
                            for (j = 0; j < nNHC_; j++)
                            {
                                k                 = i * nNHC_ + j;
                                tmp_r_[2 * k]     = ptCouplingArrays.nhpres_xi[k];
                                tmp_r_[2 * k + 1] = ptCouplingArrays.nhpres_vxi[k];
                            }
                        }
                        add_ebin(ebin_, itcb_, mdeb_n_, tmp_r_.data(), bSum);
                    }
                }
                else
                {
                    for (int i = 0; (i < nTC_); i++)
                    {
                        tmp_r_[2 * i]     = ptCouplingArrays.nosehoover_xi[i];
                        tmp_r_[2 * i + 1] = ptCouplingArrays.nosehoover_vxi[i];
                    }
                    add_ebin(ebin_, itc_, mde_n_, tmp_r_.data(), bSum);
                }
            }
        }
        else if (etc_ == TemperatureCoupling::Berendsen || etc_ == TemperatureCoupling::Yes
                 || etc_ == TemperatureCoupling::VRescale)
        {
            for (int i = 0; (i < nTC_); i++)
            {
                tmp_r_[i] = ekind->tcstat[i].lambda;
            }
            add_ebin(ebin_, itc_, nTC_, tmp_r_.data(), bSum);
        }
    }

    ebin_increase_count(1, ebin_, bSum);

    // BAR + thermodynamic integration values
    if ((fp_dhdl_ || dhc_) && bDoDHDL)
    {
        const auto& foreignTerms = enerd->foreignLambdaTerms;
        for (int i = 0; i < foreignTerms.numLambdas(); i++)
        {
            /* zero for simulated tempering */
            dE_[i] = foreignTerms.deltaH(i);
            if (!temperatures_.empty())
            {
                GMX_RELEASE_ASSERT(gmx::ssize(temperatures_) > fep_state,
                                   "Number of lambdas in state is bigger then in input record");
                GMX_RELEASE_ASSERT(
                        gmx::ssize(temperatures_) >= foreignTerms.numLambdas(),
                        "Number of lambdas in energy data is bigger then in input record");
                /* MRS: is this right, given the way we have defined the exchange probabilities? */
                /* is this even useful to have at all? */
                dE_[i] += (temperatures_[i] / temperatures_[fep_state] - 1.0) * enerd->term[F_EKIN];
            }
        }

        if (fp_dhdl_)
        {
            fprintf(fp_dhdl_, "%.4f", time);
            /* the current free energy state */

            /* print the current state if we are doing expanded ensemble */
            if (haveFepLambdaMoves_)
            {
                fprintf(fp_dhdl_, " %4d", fep_state);
            }
            /* total energy (for if the temperature changes */

            if (fep->edHdLPrintEnergy != FreeEnergyPrintEnergy::No)
            {
                switch (fep->edHdLPrintEnergy)
                {
                    case FreeEnergyPrintEnergy::Potential:
                        store_energy = enerd->term[F_EPOT];
                        break;
                    case FreeEnergyPrintEnergy::Total:
                    case FreeEnergyPrintEnergy::Yes:
                    default: store_energy = enerd->term[F_ETOT];
                }
                fprintf(fp_dhdl_, " %#.8g", store_energy);
            }

            if (fep->dhdl_derivatives == DhDlDerivativeCalculation::Yes)
            {
                for (auto i : keysOf(fep->separate_dvdl))
                {
                    if (fep->separate_dvdl[i])
                    {
                        /* assumes F_DVDL is first */
                        fprintf(fp_dhdl_, " %#.8g", enerd->term[F_DVDL + static_cast<int>(i)]);
                    }
                }
            }
            for (int i = fep->lambda_start_n; i < fep->lambda_stop_n; i++)
            {
                fprintf(fp_dhdl_, " %#.8g", dE_[i]);
            }
            if (bDynBox_ && bDiagPres_ && (epc_ != PressureCoupling::No)
                && foreignTerms.numLambdas() > 0 && fep->init_lambda_without_states < 0)
            {
                fprintf(fp_dhdl_, " %#.8g", pv); /* PV term only needed when
                                                         there are alternate state
                                                         lambda and we're not in
                                                         compatibility mode */
            }
            fprintf(fp_dhdl_, "\n");
            /* and the binary free energy output */
        }
        if (dhc_)
        {
            int idhdl = 0;
            for (auto i : keysOf(fep->separate_dvdl))
            {
                if (fep->separate_dvdl[i])
                {
                    /* assumes F_DVDL is first */
                    store_dhdl[idhdl] = enerd->term[F_DVDL + static_cast<int>(i)];
                    idhdl += 1;
                }
            }
            store_energy = enerd->term[F_ETOT];
            /* store_dh is dE */
            mde_delta_h_coll_add_dh(dhc_.get(),
                                    static_cast<double>(fep_state),
                                    store_energy,
                                    pv,
                                    store_dhdl,
                                    dE_.data() + fep->lambda_start_n,
                                    time);
        }
    }

    if (conservedEnergyTracker_)
    {
        conservedEnergyTracker_->addPoint(
                time, bEner_[F_ECONSERVED] ? enerd->term[F_ECONSERVED] : enerd->term[F_ETOT]);
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
            "Step",
            "Time",
            gmx_step_str(steps, buf),
            time);
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
    int ndisre = bDR ? fcd->disres->npair : 0;
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

    if (bOR && fcd->orires)
    {
        t_oriresdata& orires = *fcd->orires;
        diagonalize_orires_tensors(&orires);
        nr[enxOR]     = orires.numRestraints;
        block[enxOR]  = orires.orientationsTimeAndEnsembleAv.data();
        id[enxOR]     = enxOR;
        nr[enxORI]    = (orires.orientations.data() != orires.orientationsTimeAndEnsembleAv.data())
                                ? orires.numRestraints
                                : 0;
        block[enxORI] = orires.orientations.data();
        id[enxORI]    = enxORI;
        nr[enxORT]    = gmx::ssize(orires.eigenOutput);
        block[enxORT] = orires.eigenOutput.data();
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
            fr.block[b].sub[0].type = XdrDataType::Float;
            fr.block[b].sub[0].fval = block[b];
#else
            fr.block[b].sub[0].type  = XdrDataType::Double;
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
            const t_disresdata& disres = *fcd->disres;
            fr.block[db].id            = enxDISRE;
            fr.block[db].sub[0].nr     = ndisre;
            fr.block[db].sub[1].nr     = ndisre;
#if !GMX_DOUBLE
            fr.block[db].sub[0].type = XdrDataType::Float;
            fr.block[db].sub[1].type = XdrDataType::Float;
            fr.block[db].sub[0].fval = disres.rt;
            fr.block[db].sub[1].fval = disres.rm3tav;
#else
            fr.block[db].sub[0].type = XdrDataType::Double;
            fr.block[db].sub[1].type = XdrDataType::Double;
            fr.block[db].sub[0].dval = disres.rt;
            fr.block[db].sub[1].dval = disres.rm3tav;
#endif
        }
        /* here we can put new-style blocks */

        /* Free energy perturbation blocks */
        if (dhc_)
        {
            mde_delta_h_coll_handle_block(dhc_.get(), &fr, fr.nblock);
        }

        /* we can now free & reset the data in the blocks */
        if (dhc_)
        {
            mde_delta_h_coll_reset(dhc_.get());
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
        if (bOR && fcd->orires)
        {
            print_orires_log(log, fcd->orires.get());
        }

        fprintf(log, "   Energies (%s)\n", unit_energy);
        pr_ebin(log, ebin_, ie_, f_nre_ + nCrmsd_, 5, eprNORMAL, true);
        fprintf(log, "\n");
    }
}

void EnergyOutput::printAnnealingTemperatures(FILE*                   log,
                                              const SimulationGroups& groups,
                                              const t_grpopts&        opts,
                                              const gmx_ekindata_t&   ekind)
{
    if (log)
    {
        for (int i = 0; i < opts.ngtc; i++)
        {
            if (opts.annealing[i] != SimulatedAnnealing::No)
            {
                fprintf(log,
                        "Current ref_t for group %s: %8.1f\n",
                        *(groups.groupNames[groups.groups[SimulationAtomGroupType::TemperatureCoupling][i]]),
                        ekind.currentReferenceTemperature(i));
            }
        }
        fprintf(log, "\n");
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

        fprintf(log,
                "\tStatistics over %s steps using %s frames\n",
                gmx_step_str(ebin_->nsteps_sim, buf1),
                gmx_step_str(ebin_->nsum_sim, buf2));
        fprintf(log, "\n");

        fprintf(log, "   Energies (%s)\n", unit_energy);
        pr_ebin(log, ebin_, ie_, f_nre_ + nCrmsd_, 5, eprAVER, true);
        fprintf(log, "\n");

        if (bDynBox_)
        {
            pr_ebin(log, ebin_, ib_, bTricl_ ? tricl_boxs_nm.size() : boxs_nm.size(), 5, eprAVER, true);
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
            for (auto key : keysOf(bEInd_))
            {
                if (bEInd_[key])
                {
                    fprintf(log, "%12s   ", enumValueToString(key));
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
                    fprintf(log, "%*s%s-%s", padding, "", *(groups->groupNames[ni]), *(groups->groupNames[nj]));
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
        mde_delta_h_coll_update_energyhistory(dhc_.get(), enerhist);
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
                  nener,
                  enerhist.ener_sum.size(),
                  enerhist.ener_sum_sim.size());
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
        mde_delta_h_coll_restore_energyhistory(dhc_.get(), enerhist.deltaHForeignLambdas.get());
    }
}

int EnergyOutput::numEnergyTerms() const
{
    return ebin_->nener;
}

void EnergyOutput::printEnergyConservation(FILE* fplog, int simulationPart, bool usingMdIntegrator) const
{
    if (fplog == nullptr)
    {
        return;
    }

    if (conservedEnergyTracker_)
    {
        std::string partName = formatString("simulation part #%d", simulationPart);
        fprintf(fplog, "\n%s\n", conservedEnergyTracker_->energyDriftString(partName).c_str());
    }
    else if (usingMdIntegrator)
    {
        fprintf(fplog,
                "\nCannot report drift of the conserved energy quantity because simulations share "
                "state\n\n");
    }
}

} // namespace gmx
