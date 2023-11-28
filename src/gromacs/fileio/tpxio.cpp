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
#include "gmxpre.h"

/* This file is completely threadsafe - keep it that way! */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/txtdump.h"

#define TPX_TAG_RELEASE "release"

/*! \brief Tag string for the file format written to run input files
 * written by this version of the code.
 *
 * Change this if you want to change the run input format in a feature
 * branch. This ensures that there will not be different run input
 * formats around which cannot be distinguished, while not causing
 * problems rebasing the feature branch onto upstream changes. When
 * merging with mainstream GROMACS, set this tag string back to
 * TPX_TAG_RELEASE, and instead add an element to tpxv.
 */
static const std::string tpx_tag = TPX_TAG_RELEASE;

/*! \brief Enum of values that describe the contents of a tpr file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element just above tpxv_Count
 * in this enumeration, and write code below that does the right thing
 * according to the value of file_version.
 */
enum tpxv
{
    tpxv_ComputationalElectrophysiology =
            96, /**< support for ion/water position swaps (computational electrophysiology) */
    tpxv_Use64BitRandomSeed, /**< change ld_seed from int to int64_t */
    tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, /**< potentials for supporting coarse-grained force fields */
    tpxv_InteractiveMolecularDynamics, /**< interactive molecular dynamics (IMD) */
    tpxv_RemoveObsoleteParameters1,    /**< remove optimize_fft, dihre_fc, nstcheckpoint */
    tpxv_PullCoordTypeGeom,            /**< add pull type and geometry per group and flat-bottom */
    tpxv_PullGeomDirRel,               /**< add pull geometry direction-relative */
    tpxv_IntermolecularBondeds, /**< permit inter-molecular bonded interactions in the topology */
    tpxv_CompElWithSwapLayerOffset, /**< added parameters for improved CompEl setups */
    tpxv_CompElPolyatomicIonsAndMultipleIonTypes, /**< CompEl now can handle polyatomic ions and more than two types of ions */
    tpxv_RemoveAdress,                            /**< removed support for AdResS */
    tpxv_PullCoordNGroup,               /**< add ngroup to pull coord */
    tpxv_RemoveTwinRange,               /**< removed support for twin-range interactions */
    tpxv_ReplacePullPrintCOM12,         /**< Replaced print-com-1, 2 with pull-print-com */
    tpxv_PullExternalPotential,         /**< Added pull type external potential */
    tpxv_GenericParamsForElectricField, /**< Introduced KeyValueTree and moved electric field parameters */
    tpxv_AcceleratedWeightHistogram, /**< sampling with accelerated weight histogram method (AWH) */
    tpxv_RemoveImplicitSolvation,    /**< removed support for implicit solvation */
    tpxv_PullPrevStepCOMAsReference, /**< Enabled using the COM of the pull group of the last frame as reference for PBC */
    tpxv_MimicQMMM,   /**< Introduced support for MiMiC QM/MM interface */
    tpxv_PullAverage, /**< Added possibility to output average pull force and position */
    tpxv_GenericInternalParameters, /**< Added internal parameters for mdrun modules*/
    tpxv_VSite2FD,                  /**< Added 2FD type virtual site */
    tpxv_AddSizeField, /**< Added field with information about the size of the serialized tpr file in bytes, excluding the header */
    tpxv_StoreNonBondedInteractionExclusionGroup, /**< Store the non bonded interaction exclusion group in the topology */
    tpxv_VSite1,                                  /**< Added 1 type virtual site */
    tpxv_MTS,                                     /**< Added multiple time stepping */
    tpxv_RemovedConstantAcceleration, /**< Removed support for constant acceleration NEMD. */
    tpxv_TransformationPullCoord,     /**< Support for transformation pull coordinates */
    tpxv_SoftcoreGapsys,              /**< Added gapsys softcore function */
    tpxv_ReaddedConstantAcceleration, /**< Re-added support for constant acceleration NEMD. */
    tpxv_RemoveTholeRfac,             /**< Remove unused rfac parameter from thole listed force */
    tpxv_RemoveAtomtypes,             /**< Remove unused atomtypes parameter from mtop */
    tpxv_EnsembleTemperature,         /**< Add ensemble temperature settings */
    tpxv_AwhGrowthFactor,             /**< Add AWH growth factor */
    tpxv_MassRepartitioning,          /**< Add mass repartitioning */
    tpxv_AwhTargetMetricScaling,      /**< Add AWH friction optimized target distribution */
    tpxv_VerletBufferPressureTol,     /**< Add Verlet buffer pressure tolerance */
    tpxv_Count                        /**< the total number of tpxv versions */
};

/*! \brief Version number of the file format written to run input
 * files by this version of the code.
 *
 * The tpx_version increases whenever the file format in the main
 * development branch changes, due to an extension of the tpxv enum above.
 * Backward compatibility for reading old run input files is maintained
 * by checking this version number against that of the file and then using
 * the correct code path.
 *
 * When developing a feature branch that needs to change the run input
 * file format, change tpx_tag instead. */
static const int tpx_version = tpxv_Count - 1;


/*! \brief
 * Enum keeping track of incompatible changes for older TPR versions.
 *
 * The enum should be updated with a new field when editing the TOPOLOGY
 * or HEADER of the tpx format. In particular, updating ftupd or
 * changing the fields of TprHeaderVersion often trigger such needs.
 *
 * This way we can maintain forward compatibility too for all analysis tools
 * and/or external programs that only need to know the atom/residue names,
 * charges, and bond connectivity.
 *
 * It first appeared in tpx version 26, when I also moved the inputrecord
 * to the end of the tpx file, so we can just skip it if we only
 * want the topology.
 *
 * In particular, it must be increased when adding new elements to
 * ftupd, so that old code can read new .tpr files.
 */
enum class TpxGeneration : int
{
    Initial = 26, //! First version is 26
    AddSizeField, //! TPR header modified for writing as a block.
    AddVSite1,    //! ftupd changed to include VSite1 type.
    Count         //! Number of entries.
};

//! Value of Current TPR generation.
static const int tpx_generation = static_cast<int>(TpxGeneration::Count) - 1;

/* This number should be the most recent backwards incompatible version
 * I.e., if this number is 9, we cannot read tpx version 9 with this code.
 */
static const int tpx_incompatible_version = 57; // GMX4.0 has version 58


/* Struct used to maintain tpx compatibility when function types are added */
typedef struct
{
    int fvnr;  /* file version number in which the function type first appeared */
    int ftype; /* function type */
} t_ftupd;

/*
 * TODO The following three lines make little sense, please clarify if
 * you've had to work out how ftupd works.
 *
 * The entries should be ordered in:
 * 1. ascending function type number
 * 2. ascending file version number
 *
 * Because we support reading of old .tpr file versions (even when
 * mdrun can no longer run the simulation), we need to be able to read
 * obsolete t_interaction_function types. Any data read from such
 * fields is discarded. Their names have _NOLONGERUSED appended to
 * them to make things clear.
 *
 * When adding to or making breaking changes to reading this struct,
 * update TpxGeneration.
 */
static const t_ftupd ftupd[] = {
    { 70, F_RESTRBONDS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES },
    { 76, F_LINEAR_ANGLES },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS },
    { 65, F_CMAP },
    { 60, F_GB12_NOLONGERUSED },
    { 61, F_GB13_NOLONGERUSED },
    { 61, F_GB14_NOLONGERUSED },
    { 72, F_GBPOL_NOLONGERUSED },
    { 72, F_NPSOLVATION_NOLONGERUSED },
    { 93, F_LJ_RECIP },
    { 76, F_ANHARM_POL },
    { 90, F_FBPOSRES },
    { tpxv_VSite1, F_VSITE1 },
    { tpxv_VSite2FD, F_VSITE2FD },
    { tpxv_GenericInternalParameters, F_DENSITYFITTING },
    { 69, F_VTEMP_NOLONGERUSED },
    { 66, F_PDISPCORR },
    { 79, F_DVDL_COUL },
    {
            79,
            F_DVDL_VDW,
    },
    {
            79,
            F_DVDL_BONDED,
    },
    { 79, F_DVDL_RESTRAINT },
    { 79, F_DVDL_TEMPERATURE },
};
#define NFTUPD asize(ftupd)

/**************************************************************
 *
 * Now the higer level routines that do io of the structures and arrays
 *
 **************************************************************/
static void do_pullgrp_tpx_pre95(gmx::ISerializer* serializer, t_pull_group* pgrp, t_pull_coord* pcrd)
{
    rvec tmp;

    int numAtoms = pgrp->ind.size();
    serializer->doInt(&numAtoms);
    pgrp->ind.resize(numAtoms);
    serializer->doIntArray(pgrp->ind.data(), numAtoms);
    int numWeights = pgrp->weight.size();
    serializer->doInt(&numWeights);
    pgrp->weight.resize(numWeights);
    serializer->doRealArray(pgrp->weight.data(), numWeights);
    serializer->doInt(&pgrp->pbcatom);
    serializer->doRvec(&pcrd->vec.as_vec());
    clear_rvec(pcrd->origin);
    serializer->doRvec(&tmp);
    pcrd->init = tmp[0];
    serializer->doReal(&pcrd->rate);
    serializer->doReal(&pcrd->k);
    serializer->doReal(&pcrd->kB);
}

static void do_pull_group(gmx::ISerializer* serializer, t_pull_group* pgrp)
{
    int numAtoms = pgrp->ind.size();
    serializer->doInt(&numAtoms);
    pgrp->ind.resize(numAtoms);
    serializer->doIntArray(pgrp->ind.data(), numAtoms);
    int numWeights = pgrp->weight.size();
    serializer->doInt(&numWeights);
    pgrp->weight.resize(numWeights);
    serializer->doRealArray(pgrp->weight.data(), numWeights);
    serializer->doInt(&pgrp->pbcatom);
}

static void do_pull_coord(gmx::ISerializer* serializer,
                          t_pull_coord*     pcrd,
                          int               file_version,
                          PullingAlgorithm  ePullOld,
                          PullGroupGeometry eGeomOld,
                          ivec              dimOld)
{
    if (file_version >= tpxv_PullCoordNGroup)
    {
        serializer->doEnumAsInt(&pcrd->eType);
        if (file_version >= tpxv_PullExternalPotential)
        {
            if (pcrd->eType == PullingAlgorithm::External)
            {
                serializer->doString(&pcrd->externalPotentialProvider);
            }
            else
            {
                pcrd->externalPotentialProvider.clear();
            }
        }
        else
        {
            if (serializer->reading())
            {
                pcrd->externalPotentialProvider.clear();
            }
        }
        /* Note that we try to support adding new geometries without
         * changing the tpx version. This requires checks when printing the
         * geometry string and a check and fatal_error in init_pull.
         */
        serializer->doEnumAsInt(&pcrd->eGeom);
        serializer->doInt(&pcrd->ngroup);
        if (pcrd->ngroup <= c_pullCoordNgroupMax)
        {
            serializer->doIntArray(pcrd->group.data(), pcrd->ngroup);
        }
        else
        {
            /* More groups in file than supported, this must be a new geometry
             * that is not supported by our current code. Since we will not
             * use the groups for this coord (checks in the pull and WHAM code
             * ensure this), we can ignore the groups and set ngroup=0.
             */
            int* dum;
            snew(dum, pcrd->ngroup);
            serializer->doIntArray(dum, pcrd->ngroup);
            sfree(dum);

            pcrd->ngroup = 0;
        }
        serializer->doIvec(&pcrd->dim.as_vec());
        if (file_version >= tpxv_TransformationPullCoord)
        {
            serializer->doString(&pcrd->expression);
        }
        else
        {
            if (serializer->reading())
            {
                pcrd->expression.clear();
            }
        }
    }
    else
    {
        pcrd->ngroup = 2;
        serializer->doInt(&pcrd->group[0]);
        serializer->doInt(&pcrd->group[1]);
        if (file_version >= tpxv_PullCoordTypeGeom)
        {
            pcrd->ngroup = (pcrd->eGeom == PullGroupGeometry::DirectionRelative ? 4 : 2);
            serializer->doEnumAsInt(&pcrd->eType);
            serializer->doEnumAsInt(&pcrd->eGeom);
            if (pcrd->ngroup == 4)
            {
                serializer->doInt(&pcrd->group[2]);
                serializer->doInt(&pcrd->group[3]);
            }
            serializer->doIvec(&pcrd->dim.as_vec());
        }
        else
        {
            pcrd->eType = ePullOld;
            pcrd->eGeom = eGeomOld;
            copy_ivec(dimOld, pcrd->dim);
        }
    }
    serializer->doRvec(&pcrd->origin.as_vec());
    serializer->doRvec(&pcrd->vec.as_vec());
    if (file_version >= tpxv_PullCoordTypeGeom)
    {
        serializer->doBool(&pcrd->bStart);
    }
    else
    {
        /* This parameter is only printed, but not actually used by mdrun */
        pcrd->bStart = FALSE;
    }
    serializer->doReal(&pcrd->init);
    serializer->doReal(&pcrd->rate);
    serializer->doReal(&pcrd->k);
    serializer->doReal(&pcrd->kB);
}

static void do_expandedvals(gmx::ISerializer* serializer, t_expanded* expand, t_lambda* fepvals, int file_version)
{
    int n_lambda = fepvals->n_lambda;

    /* reset the lambda calculation window */
    fepvals->lambda_start_n = 0;
    fepvals->lambda_stop_n  = n_lambda;
    if (file_version >= 79)
    {
        if (n_lambda > 0)
        {
            expand->init_lambda_weights.resize(n_lambda);
            serializer->doRealArray(expand->init_lambda_weights.data(), n_lambda);
            serializer->doBool(&expand->bInit_weights);
        }

        serializer->doInt(&expand->nstexpanded);
        serializer->doEnumAsInt(&expand->elmcmove);
        serializer->doEnumAsInt(&expand->elamstats);
        serializer->doInt(&expand->lmc_repeats);
        serializer->doInt(&expand->gibbsdeltalam);
        serializer->doInt(&expand->lmc_forced_nstart);
        serializer->doInt(&expand->lmc_seed);
        serializer->doReal(&expand->mc_temp);
        serializer->doBool(&expand->bSymmetrizedTMatrix);
        serializer->doInt(&expand->nstTij);
        serializer->doInt(&expand->minvarmin);
        serializer->doInt(&expand->c_range);
        serializer->doReal(&expand->wl_scale);
        serializer->doReal(&expand->wl_ratio);
        serializer->doReal(&expand->init_wl_delta);
        serializer->doBool(&expand->bWLoneovert);
        serializer->doEnumAsInt(&expand->elmceq);
        serializer->doInt(&expand->equil_steps);
        serializer->doInt(&expand->equil_samples);
        serializer->doInt(&expand->equil_n_at_lam);
        serializer->doReal(&expand->equil_wl_delta);
        serializer->doReal(&expand->equil_ratio);
    }
}

static void do_simtempvals(gmx::ISerializer* serializer, t_simtemp* simtemp, int n_lambda, int file_version)
{
    if (file_version >= 79)
    {
        serializer->doEnumAsInt(&simtemp->eSimTempScale);
        serializer->doReal(&simtemp->simtemp_high);
        serializer->doReal(&simtemp->simtemp_low);
        if (n_lambda > 0)
        {
            if (serializer->reading())
            {
                simtemp->temperatures.resize(n_lambda);
            }
            serializer->doRealArray(simtemp->temperatures.data(), n_lambda);
        }
    }
}

static void do_imd(gmx::ISerializer* serializer, t_IMD* imd)
{
    serializer->doInt(&imd->nat);
    if (serializer->reading())
    {
        snew(imd->ind, imd->nat);
    }
    serializer->doIntArray(imd->ind, imd->nat);
}

static void do_fepvals(gmx::ISerializer* serializer, t_lambda* fepvals, int file_version)
{
    /* i is defined in the ndo_double macro; use g to iterate. */
    real rdum;

    /* free energy values */

    if (file_version >= 79)
    {
        serializer->doInt(&fepvals->init_fep_state);
        serializer->doDouble(&fepvals->init_lambda);
        serializer->doDouble(&fepvals->delta_lambda);
    }
    else if (file_version >= 59)
    {
        serializer->doDouble(&fepvals->init_lambda);
        serializer->doDouble(&fepvals->delta_lambda);
    }
    else
    {
        serializer->doReal(&rdum);
        fepvals->init_lambda = rdum;
        serializer->doReal(&rdum);
        fepvals->delta_lambda = rdum;
    }
    if (file_version >= 79)
    {
        serializer->doInt(&fepvals->n_lambda);
        for (auto g : keysOf(fepvals->all_lambda))
        {
            if (fepvals->n_lambda > 0)
            {
                fepvals->all_lambda[g].resize(fepvals->n_lambda);
                serializer->doDoubleArray(fepvals->all_lambda[g].data(), fepvals->n_lambda);
                serializer->doBoolArray(fepvals->separate_dvdl.begin(), fepvals->separate_dvdl.size());
            }
            else if (fepvals->init_lambda >= 0)
            {
                fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep] = TRUE;
            }
        }
    }
    else if (file_version >= 64)
    {
        serializer->doInt(&fepvals->n_lambda);
        if (serializer->reading())
        {
            /* still allocate the all_lambda array's contents. */
            for (auto g : keysOf(fepvals->all_lambda))
            {
                if (fepvals->n_lambda > 0)
                {
                    fepvals->all_lambda[g].resize(fepvals->n_lambda);
                }
            }
        }
        serializer->doDoubleArray(fepvals->all_lambda[FreeEnergyPerturbationCouplingType::Fep].data(),
                                  fepvals->n_lambda);
        if (fepvals->init_lambda >= 0)
        {
            fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep] = TRUE;

            if (serializer->reading())
            {
                /* copy the contents of the efptFEP lambda component to all
                   the other components */
                for (auto g : keysOf(fepvals->all_lambda))
                {
                    for (int h = 0; h < fepvals->n_lambda; h++)
                    {
                        if (g != FreeEnergyPerturbationCouplingType::Fep)
                        {
                            fepvals->all_lambda[g][h] =
                                    fepvals->all_lambda[FreeEnergyPerturbationCouplingType::Fep][h];
                        }
                    }
                }
            }
        }
    }
    else
    {
        fepvals->n_lambda = 0;
        if (fepvals->init_lambda >= 0)
        {
            fepvals->separate_dvdl[FreeEnergyPerturbationCouplingType::Fep] = TRUE;
        }
    }
    serializer->doReal(&fepvals->sc_alpha);
    serializer->doInt(&fepvals->sc_power);
    if (file_version >= 79)
    {
        serializer->doReal(&fepvals->sc_r_power);
    }
    else
    {
        fepvals->sc_r_power = 6.0;
    }
    if (fepvals->sc_r_power != 6.0)
    {
        gmx_fatal(FARGS, "Only sc-r-power=6 is supported. Value in file is %f", fepvals->sc_r_power);
    }
    serializer->doReal(&fepvals->sc_sigma);
    if (serializer->reading())
    {
        if (file_version >= 71)
        {
            fepvals->sc_sigma_min = fepvals->sc_sigma;
        }
        else
        {
            fepvals->sc_sigma_min = 0;
        }
    }
    if (file_version >= 79)
    {
        serializer->doBool(&fepvals->bScCoul);
    }
    else
    {
        fepvals->bScCoul = TRUE;
    }
    if (file_version >= 64)
    {
        serializer->doInt(&fepvals->nstdhdl);
    }
    else
    {
        fepvals->nstdhdl = 1;
    }

    if (file_version >= 73)
    {
        serializer->doEnumAsInt(&fepvals->separate_dhdl_file);
        serializer->doEnumAsInt(&fepvals->dhdl_derivatives);
    }
    else
    {
        fepvals->separate_dhdl_file = SeparateDhdlFile::Yes;
        fepvals->dhdl_derivatives   = DhDlDerivativeCalculation::Yes;
    }
    if (file_version >= 71)
    {
        serializer->doInt(&fepvals->dh_hist_size);
        serializer->doDouble(&fepvals->dh_hist_spacing);
    }
    else
    {
        fepvals->dh_hist_size    = 0;
        fepvals->dh_hist_spacing = 0.1;
    }
    if (file_version >= 79)
    {
        serializer->doEnumAsInt(&fepvals->edHdLPrintEnergy);
    }
    else
    {
        fepvals->edHdLPrintEnergy = FreeEnergyPrintEnergy::No;
    }
    if (file_version >= tpxv_SoftcoreGapsys)
    {
        serializer->doInt(reinterpret_cast<int*>(&fepvals->softcoreFunction));
        serializer->doReal(&fepvals->scGapsysScaleLinpointLJ);
        serializer->doReal(&fepvals->scGapsysScaleLinpointQ);
        serializer->doReal(&fepvals->scGapsysSigmaLJ);
    }
    else
    {
        fepvals->softcoreFunction        = SoftcoreType::Beutler;
        fepvals->scGapsysScaleLinpointLJ = 0.85;
        fepvals->scGapsysScaleLinpointQ  = 0.3;
        fepvals->scGapsysSigmaLJ         = 0.3;
    }

    /* handle lambda_neighbors */
    if ((file_version >= 83 && file_version < 90) || file_version >= 92)
    {
        serializer->doInt(&fepvals->lambda_neighbors);
        if ((fepvals->lambda_neighbors >= 0) && (fepvals->init_fep_state >= 0)
            && (fepvals->init_lambda < 0))
        {
            fepvals->lambda_start_n = (fepvals->init_fep_state - fepvals->lambda_neighbors);
            fepvals->lambda_stop_n  = (fepvals->init_fep_state + fepvals->lambda_neighbors + 1);
            if (fepvals->lambda_start_n < 0)
            {
                fepvals->lambda_start_n = 0;
            }
            if (fepvals->lambda_stop_n >= fepvals->n_lambda)
            {
                fepvals->lambda_stop_n = fepvals->n_lambda;
            }
        }
        else
        {
            fepvals->lambda_start_n = 0;
            fepvals->lambda_stop_n  = fepvals->n_lambda;
        }
    }
    else
    {
        fepvals->lambda_start_n = 0;
        fepvals->lambda_stop_n  = fepvals->n_lambda;
    }
}

static void do_pull(gmx::ISerializer* serializer, pull_params_t* pull, int file_version, PullingAlgorithm ePullOld)
{
    PullGroupGeometry eGeomOld = PullGroupGeometry::Count;
    ivec              dimOld;
    int               g;

    if (file_version >= 95)
    {
        serializer->doInt(&pull->ngroup);
    }
    serializer->doInt(&pull->ncoord);
    if (file_version < 95)
    {
        pull->ngroup = pull->ncoord + 1;
    }
    if (file_version < tpxv_PullCoordTypeGeom)
    {
        real dum;

        serializer->doEnumAsInt(&eGeomOld);
        serializer->doIvec(&dimOld);
        /* The inner cylinder radius, now removed */
        serializer->doReal(&dum);
    }
    serializer->doReal(&pull->cylinder_r);
    serializer->doReal(&pull->constr_tol);
    if (file_version >= 95)
    {
        serializer->doBool(&pull->bPrintCOM);
        /* With file_version < 95 this value is set below */
    }
    if (file_version >= tpxv_ReplacePullPrintCOM12)
    {
        serializer->doBool(&pull->bPrintRefValue);
        serializer->doBool(&pull->bPrintComp);
    }
    else if (file_version >= tpxv_PullCoordTypeGeom)
    {
        int idum;
        serializer->doInt(&idum); /* used to be bPrintCOM2 */
        serializer->doBool(&pull->bPrintRefValue);
        serializer->doBool(&pull->bPrintComp);
    }
    else
    {
        pull->bPrintRefValue = FALSE;
        pull->bPrintComp     = TRUE;
    }
    serializer->doInt(&pull->nstxout);
    serializer->doInt(&pull->nstfout);
    if (file_version >= tpxv_PullPrevStepCOMAsReference)
    {
        serializer->doBool(&pull->bSetPbcRefToPrevStepCOM);
    }
    else
    {
        pull->bSetPbcRefToPrevStepCOM = FALSE;
    }
    pull->group.resize(pull->ngroup);
    pull->coord.resize(pull->ncoord);
    if (file_version < 95)
    {
        /* epullgPOS for position pulling, before epullgDIRPBC was removed */
        if (eGeomOld == PullGroupGeometry::DirectionPBC)
        {
            gmx_fatal(FARGS, "pull-geometry=position is no longer supported");
        }
        if (eGeomOld > PullGroupGeometry::DirectionPBC)
        {
            switch (eGeomOld)
            {
                case (PullGroupGeometry::DirectionRelative):
                    eGeomOld = PullGroupGeometry::DirectionPBC;
                    break;
                case (PullGroupGeometry::Angle):
                    eGeomOld = PullGroupGeometry::DirectionRelative;
                    break;
                case (PullGroupGeometry::Dihedral): eGeomOld = PullGroupGeometry::Angle; break;
                case (PullGroupGeometry::AngleAxis): eGeomOld = PullGroupGeometry::Dihedral; break;
                case (PullGroupGeometry::Count): eGeomOld = PullGroupGeometry::AngleAxis; break;
                default: GMX_RELEASE_ASSERT(false, "Unhandled old pull type");
            }
        }

        for (g = 0; g < pull->ngroup; g++)
        {
            /* We read and ignore a pull coordinate for group 0 */
            do_pullgrp_tpx_pre95(serializer, &pull->group[g], &pull->coord[std::max(g - 1, 0)]);
            if (g > 0)
            {
                pull->coord[g - 1].group[0] = 0;
                pull->coord[g - 1].group[1] = g;
            }
        }

        pull->bPrintCOM = (!pull->group[0].ind.empty());
    }
    else
    {
        for (g = 0; g < pull->ngroup; g++)
        {
            do_pull_group(serializer, &pull->group[g]);
        }
        for (g = 0; g < pull->ncoord; g++)
        {
            do_pull_coord(serializer, &pull->coord[g], file_version, ePullOld, eGeomOld, dimOld);
            if (serializer->reading())
            {
                pull->coord[g].coordIndex = g;
            }
        }
    }
    if (file_version >= tpxv_PullAverage)
    {
        gmx_bool v;

        v = pull->bXOutAverage;
        serializer->doBool(&v);
        pull->bXOutAverage = v;
        v                  = pull->bFOutAverage;
        serializer->doBool(&v);
        pull->bFOutAverage = v;
    }
}


static void do_rotgrp(gmx::ISerializer* serializer, t_rotgrp* rotg)
{
    serializer->doEnumAsInt(&rotg->eType);
    if (serializer->reading())
    {
        int temp = 0;
        serializer->doInt(&temp);
        rotg->bMassW = static_cast<bool>(temp);
    }
    else
    {
        int temp = static_cast<int>(rotg->bMassW);
        serializer->doInt(&temp);
    }
    serializer->doInt(&rotg->nat);
    if (serializer->reading())
    {
        snew(rotg->ind, rotg->nat);
    }
    serializer->doIntArray(rotg->ind, rotg->nat);
    if (serializer->reading())
    {
        rotg->x_ref_original.resize(rotg->nat);
    }
    for (gmx::RVec& x : rotg->x_ref_original)
    {
        serializer->doRvec(as_rvec_array(&x));
    }
    serializer->doRvec(&rotg->inputVec);
    serializer->doRvec(&rotg->pivot);
    serializer->doReal(&rotg->rate);
    serializer->doReal(&rotg->k);
    serializer->doReal(&rotg->slab_dist);
    serializer->doReal(&rotg->min_gaussian);
    serializer->doReal(&rotg->eps);
    serializer->doEnumAsInt(&rotg->eFittype);
    serializer->doInt(&rotg->PotAngle_nstep);
    serializer->doReal(&rotg->PotAngle_step);
}

static void do_rot(gmx::ISerializer* serializer, t_rot* rot)
{
    int numGroups = rot->grp.size();

    serializer->doInt(&numGroups);
    serializer->doInt(&rot->nstrout);
    serializer->doInt(&rot->nstsout);
    if (serializer->reading())
    {
        rot->grp.resize(numGroups);
    }
    for (auto& grp : rot->grp)
    {
        do_rotgrp(serializer, &grp);
    }
}


static void do_swapgroup(gmx::ISerializer* serializer, t_swapGroup* g)
{

    /* Name of the group or molecule */
    std::string buf;
    if (serializer->reading())
    {
        serializer->doString(&buf);
        g->molname = gmx_strdup(buf.c_str());
    }
    else
    {
        buf = g->molname;
        serializer->doString(&buf);
    }

    /* Number of atoms in the group */
    serializer->doInt(&g->nat);

    /* The group's atom indices */
    if (serializer->reading())
    {
        snew(g->ind, g->nat);
    }
    serializer->doIntArray(g->ind, g->nat);

    /* Requested counts for compartments A and B */
    serializer->doIntArray(g->nmolReq.data(), static_cast<int>(Compartment::Count));
}

static void do_swapcoords_tpx(gmx::ISerializer* serializer, t_swapcoords* swap, int file_version)
{
    /* Enums for better readability of the code */
    enum
    {
        eCompA = 0,
        eCompB
    };
    enum
    {
        eChannel0 = 0,
        eChannel1
    };


    if (file_version >= tpxv_CompElPolyatomicIonsAndMultipleIonTypes)
    {
        /* The total number of swap groups is the sum of the fixed groups
         * (split0, split1, solvent), and the user-defined groups (2+ types of ions)
         */
        serializer->doInt(&swap->ngrp);
        if (serializer->reading())
        {
            snew(swap->grp, swap->ngrp);
        }
        for (int ig = 0; ig < swap->ngrp; ig++)
        {
            do_swapgroup(serializer, &swap->grp[ig]);
        }
        serializer->doBool(&swap->massw_split[eChannel0]);
        serializer->doBool(&swap->massw_split[eChannel1]);
        serializer->doInt(&swap->nstswap);
        serializer->doInt(&swap->nAverage);
        serializer->doReal(&swap->threshold);
        serializer->doReal(&swap->cyl0r);
        serializer->doReal(&swap->cyl0u);
        serializer->doReal(&swap->cyl0l);
        serializer->doReal(&swap->cyl1r);
        serializer->doReal(&swap->cyl1u);
        serializer->doReal(&swap->cyl1l);
    }
    else
    {
        /*** Support reading older CompEl .tpr files ***/

        /* In the original CompEl .tpr files, we always have 5 groups: */
        swap->ngrp = 5;
        snew(swap->grp, swap->ngrp);

        swap->grp[static_cast<int>(SwapGroupSplittingType::Split0)].molname = gmx_strdup("split0"); // group 0: split0
        swap->grp[static_cast<int>(SwapGroupSplittingType::Split1)].molname = gmx_strdup("split1"); // group 1: split1
        swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].molname =
                gmx_strdup("solvent");                // group 2: solvent
        swap->grp[3].molname = gmx_strdup("anions");  // group 3: anions
        swap->grp[4].molname = gmx_strdup("cations"); // group 4: cations

        serializer->doInt(&swap->grp[3].nat);
        serializer->doInt(&swap->grp[static_cast<int>(SwapGroupSplittingType::Solvent)].nat);
        serializer->doInt(&swap->grp[static_cast<int>(SwapGroupSplittingType::Split0)].nat);
        serializer->doBool(&swap->massw_split[eChannel0]);
        serializer->doInt(&swap->grp[static_cast<int>(SwapGroupSplittingType::Split1)].nat);
        serializer->doBool(&swap->massw_split[eChannel1]);
        serializer->doInt(&swap->nstswap);
        serializer->doInt(&swap->nAverage);
        serializer->doReal(&swap->threshold);
        serializer->doReal(&swap->cyl0r);
        serializer->doReal(&swap->cyl0u);
        serializer->doReal(&swap->cyl0l);
        serializer->doReal(&swap->cyl1r);
        serializer->doReal(&swap->cyl1u);
        serializer->doReal(&swap->cyl1l);

        // The order[] array keeps compatibility with older .tpr files
        // by reading in the groups in the classic order
        {
            const int order[4] = { 3,
                                   static_cast<int>(SwapGroupSplittingType::Solvent),
                                   static_cast<int>(SwapGroupSplittingType::Split0),
                                   static_cast<int>(SwapGroupSplittingType::Split1) };

            for (int ig = 0; ig < 4; ig++)
            {
                int g = order[ig];
                snew(swap->grp[g].ind, swap->grp[g].nat);
                serializer->doIntArray(swap->grp[g].ind, swap->grp[g].nat);
            }
        }

        for (int j = eCompA; j <= eCompB; j++)
        {
            serializer->doInt(&swap->grp[3].nmolReq[j]); // group 3 = anions
            serializer->doInt(&swap->grp[4].nmolReq[j]); // group 4 = cations
        }
    } /* End support reading older CompEl .tpr files */

    if (file_version >= tpxv_CompElWithSwapLayerOffset)
    {
        serializer->doReal(&swap->bulkOffset[eCompA]);
        serializer->doReal(&swap->bulkOffset[eCompB]);
    }
}

static void do_legacy_efield(gmx::ISerializer* serializer, gmx::KeyValueTreeObjectBuilder* root)
{
    const char* const dimName[] = { "x", "y", "z" };

    auto appliedForcesObj = root->addObject("applied-forces");
    auto efieldObj        = appliedForcesObj.addObject("electric-field");
    // The content of the tpr file for this feature has
    // been the same since gromacs 4.0 that was used for
    // developing.
    for (int j = 0; j < DIM; ++j)
    {
        int n, nt;
        serializer->doInt(&n);
        serializer->doInt(&nt);
        std::vector<real> aa(n + 1), phi(nt + 1), at(nt + 1), phit(nt + 1);
        serializer->doRealArray(aa.data(), n);
        serializer->doRealArray(phi.data(), n);
        serializer->doRealArray(at.data(), nt);
        serializer->doRealArray(phit.data(), nt);
        if (n > 0)
        {
            if (n > 1 || nt > 1)
            {
                gmx_fatal(FARGS,
                          "Can not handle tpr files with more than one electric field term per "
                          "direction.");
            }
            auto dimObj = efieldObj.addObject(dimName[j]);
            dimObj.addValue<real>("E0", aa[0]);
            dimObj.addValue<real>("omega", at[0]);
            dimObj.addValue<real>("t0", phi[0]);
            dimObj.addValue<real>("sigma", phit[0]);
        }
    }
}


static void do_inputrec(gmx::ISerializer* serializer, t_inputrec* ir, int file_version)
{
    int      i, j, k, idum = 0;
    real     rdum;
    gmx_bool bdum = false;

    ir->tpxFileVersion = file_version;

    if (file_version != tpx_version)
    {
        /* Give a warning about features that are not accessible */
        fprintf(stderr, "Note: file tpx version %d, software tpx version %d\n", file_version, tpx_version);
    }

    if (file_version == 0)
    {
        return;
    }

    gmx::KeyValueTreeBuilder       paramsBuilder;
    gmx::KeyValueTreeObjectBuilder paramsObj = paramsBuilder.rootObject();

    /* Basic inputrec stuff */
    serializer->doEnumAsInt(&ir->eI);
    if (file_version >= 62)
    {
        serializer->doInt64(&ir->nsteps);
    }
    else
    {
        serializer->doInt(&idum);
        ir->nsteps = idum;
    }

    if (file_version >= 62)
    {
        serializer->doInt64(&ir->init_step);
    }
    else
    {
        serializer->doInt(&idum);
        ir->init_step = idum;
    }

    serializer->doInt(&ir->simulation_part);

    if (file_version >= tpxv_MTS)
    {
        serializer->doBool(&ir->useMts);
        int numLevels = ir->mtsLevels.size();
        if (ir->useMts)
        {
            serializer->doInt(&numLevels);
        }
        ir->mtsLevels.resize(numLevels);
        for (auto& mtsLevel : ir->mtsLevels)
        {
            int forceGroups = mtsLevel.forceGroups.to_ulong();
            serializer->doInt(&forceGroups);
            mtsLevel.forceGroups = std::bitset<static_cast<int>(gmx::MtsForceGroups::Count)>(forceGroups);
            serializer->doInt(&mtsLevel.stepFactor);
        }
    }
    else
    {
        ir->useMts = false;
        ir->mtsLevels.clear();
    }

    if (file_version >= tpxv_MassRepartitioning)
    {
        serializer->doReal(&ir->massRepartitionFactor);
    }
    else
    {
        ir->massRepartitionFactor = 1;
    }

    if (file_version >= tpxv_EnsembleTemperature)
    {
        serializer->doEnumAsInt(&ir->ensembleTemperatureSetting);
        serializer->doReal(&ir->ensembleTemperature);
    }

    if (file_version >= 67)
    {
        serializer->doInt(&ir->nstcalcenergy);
    }
    else
    {
        ir->nstcalcenergy = 1;
    }
    if (file_version >= 81)
    {
        serializer->doEnumAsInt(&ir->cutoff_scheme);
        if (file_version < 94)
        {
            // Need to invert the scheme order
            switch (ir->cutoff_scheme)
            {
                case (CutoffScheme::Group): ir->cutoff_scheme = CutoffScheme::Verlet; break;
                case (CutoffScheme::Verlet): ir->cutoff_scheme = CutoffScheme::Group; break;
                default: GMX_RELEASE_ASSERT(false, "Unhandled cutoff scheme type");
            }
        }
    }
    else
    {
        ir->cutoff_scheme = CutoffScheme::Group;
    }
    serializer->doInt(&idum); /* used to be ns_type; not used anymore */
    serializer->doInt(&ir->nstlist);
    serializer->doInt(&idum); /* used to be ndelta; not used anymore */

    serializer->doReal(&ir->rtpi);

    serializer->doInt(&ir->nstcomm);
    serializer->doEnumAsInt(&ir->comm_mode);

    /* ignore nstcheckpoint */
    if (file_version < tpxv_RemoveObsoleteParameters1)
    {
        serializer->doInt(&idum);
    }

    serializer->doInt(&ir->nstcgsteep);

    serializer->doInt(&ir->nbfgscorr);

    serializer->doInt(&ir->nstlog);
    serializer->doInt(&ir->nstxout);
    serializer->doInt(&ir->nstvout);
    serializer->doInt(&ir->nstfout);
    serializer->doInt(&ir->nstenergy);
    serializer->doInt(&ir->nstxout_compressed);
    if (file_version >= 59)
    {
        serializer->doDouble(&ir->init_t);
        serializer->doDouble(&ir->delta_t);
    }
    else
    {
        serializer->doReal(&rdum);
        ir->init_t = rdum;
        serializer->doReal(&rdum);
        ir->delta_t = rdum;
    }
    serializer->doReal(&ir->x_compression_precision);
    if (file_version >= 81)
    {
        serializer->doReal(&ir->verletbuf_tol);
    }
    else
    {
        ir->verletbuf_tol = 0;
    }
    if (file_version >= tpxv_VerletBufferPressureTol)
    {
        serializer->doReal(&ir->verletBufferPressureTolerance);
    }
    else
    {
        ir->verletBufferPressureTolerance = -1;
    }
    serializer->doReal(&ir->rlist);
    if (file_version >= 67 && file_version < tpxv_RemoveTwinRange)
    {
        if (serializer->reading())
        {
            // Reading such a file version could invoke the twin-range
            // scheme, about which mdrun should give a fatal error.
            real dummy_rlistlong = -1;
            serializer->doReal(&dummy_rlistlong);

            ir->useTwinRange = (ir->rlist > 0 && (dummy_rlistlong == 0 || dummy_rlistlong > ir->rlist));
            // When true, this forces mdrun to issue an error (regardless of
            // ir->cutoff_scheme).
            //
            // Otherwise, grompp used to set rlistlong actively. Users
            // were probably also confused and set rlistlong == rlist.
            // However, in all remaining cases, it is safe to let
            // mdrun proceed normally.
        }
    }
    else
    {
        // No need to read or write anything
        ir->useTwinRange = false;
    }
    if (file_version >= 82 && file_version != 90)
    {
        // Multiple time-stepping is no longer enabled, but the old
        // support required the twin-range scheme, for which mdrun
        // already emits a fatal error.
        int dummy_nstcalclr = -1;
        serializer->doInt(&dummy_nstcalclr);
    }
    serializer->doEnumAsInt(&ir->coulombtype);
    if (file_version >= 81)
    {
        serializer->doEnumAsInt(&ir->coulomb_modifier);
    }
    else
    {
        ir->coulomb_modifier = (ir->cutoff_scheme == CutoffScheme::Verlet ? InteractionModifiers::PotShift
                                                                          : InteractionModifiers::None);
    }
    serializer->doReal(&ir->rcoulomb_switch);
    serializer->doReal(&ir->rcoulomb);
    serializer->doEnumAsInt(&ir->vdwtype);
    if (file_version >= 81)
    {
        serializer->doEnumAsInt(&ir->vdw_modifier);
    }
    else
    {
        ir->vdw_modifier = (ir->cutoff_scheme == CutoffScheme::Verlet ? InteractionModifiers::PotShift
                                                                      : InteractionModifiers::None);
    }
    serializer->doReal(&ir->rvdw_switch);
    serializer->doReal(&ir->rvdw);
    serializer->doEnumAsInt(&ir->eDispCorr);
    serializer->doReal(&ir->epsilon_r);
    serializer->doReal(&ir->epsilon_rf);
    serializer->doReal(&ir->tabext);

    // This permits reading a .tpr file that used implicit solvent,
    // and later permitting mdrun to refuse to run it.
    if (serializer->reading())
    {
        if (file_version < tpxv_RemoveImplicitSolvation)
        {
            serializer->doInt(&idum);
            serializer->doInt(&idum);
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            serializer->doInt(&idum);
            ir->implicit_solvent = (idum > 0);
        }
        else
        {
            ir->implicit_solvent = false;
        }
        if (file_version < tpxv_RemoveImplicitSolvation)
        {
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            if (file_version >= 60)
            {
                serializer->doReal(&rdum);
                serializer->doInt(&idum);
            }
            serializer->doReal(&rdum);
        }
    }

    if (file_version >= 81)
    {
        serializer->doReal(&ir->fourier_spacing);
    }
    else
    {
        ir->fourier_spacing = 0.0;
    }
    serializer->doInt(&ir->nkx);
    serializer->doInt(&ir->nky);
    serializer->doInt(&ir->nkz);
    serializer->doInt(&ir->pme_order);
    serializer->doReal(&ir->ewald_rtol);

    if (file_version >= 93)
    {
        serializer->doReal(&ir->ewald_rtol_lj);
    }
    else
    {
        ir->ewald_rtol_lj = ir->ewald_rtol;
    }
    serializer->doEnumAsInt(&ir->ewald_geometry);
    serializer->doReal(&ir->epsilon_surface);

    /* ignore bOptFFT */
    if (file_version < tpxv_RemoveObsoleteParameters1)
    {
        serializer->doBool(&bdum);
    }

    if (file_version >= 93)
    {
        serializer->doEnumAsInt(&ir->ljpme_combination_rule);
    }
    serializer->doBool(&ir->bContinuation);
    serializer->doEnumAsInt(&ir->etc);
    /* before version 18, ir->etc was a gmx_bool (ir->btc),
     * but the values 0 and 1 still mean no and
     * berendsen temperature coupling, respectively.
     */
    if (file_version >= 79)
    {
        serializer->doBool(&ir->bPrintNHChains);
    }
    if (file_version >= 71)
    {
        serializer->doInt(&ir->nsttcouple);
    }
    else
    {
        ir->nsttcouple = ir->nstcalcenergy;
    }
    serializer->doEnumAsInt(&ir->pressureCouplingOptions.epc);
    serializer->doEnumAsInt(&ir->pressureCouplingOptions.epct);
    if (file_version >= 71)
    {
        serializer->doInt(&ir->pressureCouplingOptions.nstpcouple);
    }
    else
    {
        ir->pressureCouplingOptions.nstpcouple = ir->nstcalcenergy;
    }
    serializer->doReal(&ir->pressureCouplingOptions.tau_p);
    serializer->doRvec(&ir->pressureCouplingOptions.ref_p[XX]);
    serializer->doRvec(&ir->pressureCouplingOptions.ref_p[YY]);
    serializer->doRvec(&ir->pressureCouplingOptions.ref_p[ZZ]);
    serializer->doRvec(&ir->pressureCouplingOptions.compress[XX]);
    serializer->doRvec(&ir->pressureCouplingOptions.compress[YY]);
    serializer->doRvec(&ir->pressureCouplingOptions.compress[ZZ]);
    serializer->doEnumAsInt(&ir->pressureCouplingOptions.refcoord_scaling);
    serializer->doRvec(&ir->posres_com);
    serializer->doRvec(&ir->posres_comB);

    if (file_version < 79)
    {
        serializer->doInt(&ir->andersen_seed);
    }
    else
    {
        ir->andersen_seed = 0;
    }

    serializer->doReal(&ir->shake_tol);

    serializer->doEnumAsInt(&ir->efep);
    if (serializer->reading())
    {
        ir->fepvals = std::make_unique<t_lambda>();
    }
    do_fepvals(serializer, ir->fepvals.get(), file_version);

    if (file_version >= 79)
    {
        serializer->doBool(&ir->bSimTemp);
        if (ir->bSimTemp)
        {
            ir->bSimTemp = TRUE;
        }
    }
    else
    {
        ir->bSimTemp = FALSE;
    }
    if (ir->bSimTemp)
    {
        if (serializer->reading())
        {
            ir->simtempvals = std::make_unique<t_simtemp>();
        }
        do_simtempvals(serializer, ir->simtempvals.get(), ir->fepvals->n_lambda, file_version);
    }

    if (file_version >= 79)
    {
        serializer->doBool(&ir->bExpanded);
    }
    if (ir->bExpanded)
    {
        if (serializer->reading())
        {
            ir->expandedvals = std::make_unique<t_expanded>();
        }
        do_expandedvals(serializer, ir->expandedvals.get(), ir->fepvals.get(), file_version);
    }

    serializer->doEnumAsInt(&ir->eDisre);
    serializer->doEnumAsInt(&ir->eDisreWeighting);
    serializer->doBool(&ir->bDisreMixed);
    serializer->doReal(&ir->dr_fc);
    serializer->doReal(&ir->dr_tau);
    serializer->doInt(&ir->nstdisreout);
    serializer->doReal(&ir->orires_fc);
    serializer->doReal(&ir->orires_tau);
    serializer->doInt(&ir->nstorireout);

    /* ignore dihre_fc */
    if (file_version < 79)
    {
        serializer->doReal(&rdum);
    }

    serializer->doReal(&ir->em_stepsize);
    serializer->doReal(&ir->em_tol);
    serializer->doBool(&ir->bShakeSOR);
    serializer->doInt(&ir->niter);
    serializer->doReal(&ir->fc_stepsize);
    serializer->doEnumAsInt(&ir->eConstrAlg);
    serializer->doInt(&ir->nProjOrder);
    serializer->doReal(&ir->LincsWarnAngle);
    serializer->doInt(&ir->nLincsIter);
    serializer->doReal(&ir->bd_fric);
    if (file_version >= tpxv_Use64BitRandomSeed)
    {
        serializer->doInt64(&ir->ld_seed);
    }
    else
    {
        serializer->doInt(&idum);
        ir->ld_seed = idum;
    }

    for (i = 0; i < DIM; i++)
    {
        serializer->doRvec(&ir->deform[i]);
    }
    serializer->doReal(&ir->cos_accel);

    serializer->doInt(&ir->userint1);
    serializer->doInt(&ir->userint2);
    serializer->doInt(&ir->userint3);
    serializer->doInt(&ir->userint4);
    serializer->doReal(&ir->userreal1);
    serializer->doReal(&ir->userreal2);
    serializer->doReal(&ir->userreal3);
    serializer->doReal(&ir->userreal4);

    /* AdResS is removed, but we need to be able to read old files,
       and let mdrun refuse to run them */
    if (file_version >= 77 && file_version < tpxv_RemoveAdress)
    {
        serializer->doBool(&ir->bAdress);
        if (ir->bAdress)
        {
            int  idum, numThermoForceGroups, numEnergyGroups;
            real rdum;
            rvec rvecdum;
            serializer->doInt(&idum);
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            serializer->doReal(&rdum);
            serializer->doInt(&idum);
            serializer->doInt(&idum);
            serializer->doRvec(&rvecdum);
            serializer->doInt(&numThermoForceGroups);
            serializer->doReal(&rdum);
            serializer->doInt(&numEnergyGroups);
            serializer->doInt(&idum);

            if (numThermoForceGroups > 0)
            {
                std::vector<int> idumn(numThermoForceGroups);
                serializer->doIntArray(idumn.data(), idumn.size());
            }
            if (numEnergyGroups > 0)
            {
                std::vector<int> idumn(numEnergyGroups);
                serializer->doIntArray(idumn.data(), idumn.size());
            }
        }
    }
    else
    {
        ir->bAdress = FALSE;
    }

    /* pull stuff */
    {
        PullingAlgorithm ePullOld = PullingAlgorithm::Umbrella;

        if (file_version >= tpxv_PullCoordTypeGeom)
        {
            serializer->doBool(&ir->bPull);
        }
        else
        {
            serializer->doEnumAsInt(&ePullOld);
            ir->bPull = (ePullOld != PullingAlgorithm::Umbrella);
            /* We removed the first ePull=ePullNo for the enum */
            switch (ePullOld)
            {
                case (PullingAlgorithm::Umbrella): break; // this is equal to not using pulling
                case (PullingAlgorithm::Constraint): ePullOld = PullingAlgorithm::Umbrella; break;
                case (PullingAlgorithm::ConstantForce):
                    ePullOld = PullingAlgorithm::Constraint;
                    break;
                case (PullingAlgorithm::FlatBottom):
                    ePullOld = PullingAlgorithm::ConstantForce;
                    break;
                case (PullingAlgorithm::FlatBottomHigh):
                    ePullOld = PullingAlgorithm::FlatBottom;
                    break;
                case (PullingAlgorithm::External):
                    ePullOld = PullingAlgorithm::FlatBottomHigh;
                    break;
                case (PullingAlgorithm::Count): ePullOld = PullingAlgorithm::External; break;
                default: GMX_RELEASE_ASSERT(false, "Unhandled old pull algorithm");
            }
        }
        if (ir->bPull)
        {
            if (serializer->reading())
            {
                ir->pull = std::make_unique<pull_params_t>();
            }
            do_pull(serializer, ir->pull.get(), file_version, ePullOld);
        }
    }

    if (file_version >= tpxv_AcceleratedWeightHistogram)
    {
        serializer->doBool(&ir->bDoAwh);

        if (ir->bDoAwh)
        {
            if (serializer->reading())
            {
                ir->awhParams =
                        std::make_unique<gmx::AwhParams>(serializer,
                                                         file_version < tpxv_AwhGrowthFactor,
                                                         file_version < tpxv_AwhTargetMetricScaling);
            }
            else
            {
                ir->awhParams->serialize(serializer);
            }
        }
    }
    else
    {
        ir->bDoAwh = FALSE;
    }

    /* Enforced rotation */
    if (file_version >= 74)
    {
        serializer->doBool(&ir->bRot);
        if (ir->bRot)
        {
            if (serializer->reading())
            {
                ir->rot = std::make_unique<t_rot>();
            }
            do_rot(serializer, ir->rot.get());
        }
    }
    else
    {
        ir->bRot = FALSE;
    }

    /* Interactive molecular dynamics */
    if (file_version >= tpxv_InteractiveMolecularDynamics)
    {
        serializer->doBool(&ir->bIMD);
        if (ir->bIMD)
        {
            if (serializer->reading())
            {
                snew(ir->imd, 1);
            }
            do_imd(serializer, ir->imd);
        }
    }
    else
    {
        /* We don't support IMD sessions for old .tpr files */
        ir->bIMD = FALSE;
    }

    /* grpopts stuff */
    serializer->doInt(&ir->opts.ngtc);
    if (file_version >= 69)
    {
        serializer->doInt(&ir->opts.nhchainlength);
    }
    else
    {
        ir->opts.nhchainlength = 1;
    }
    if (serializer->reading() && file_version >= tpxv_RemovedConstantAcceleration
        && file_version < tpxv_ReaddedConstantAcceleration)
    {
        ir->opts.ngacc = 0;
    }
    else
    {
        serializer->doInt(&ir->opts.ngacc);
    }
    serializer->doInt(&ir->opts.ngfrz);
    serializer->doInt(&ir->opts.ngener);

    if (serializer->reading())
    {
        snew(ir->opts.nrdf, ir->opts.ngtc);
        snew(ir->opts.ref_t, ir->opts.ngtc);
        snew(ir->opts.annealing, ir->opts.ngtc);
        snew(ir->opts.anneal_npoints, ir->opts.ngtc);
        snew(ir->opts.anneal_time, ir->opts.ngtc);
        snew(ir->opts.anneal_temp, ir->opts.ngtc);
        snew(ir->opts.tau_t, ir->opts.ngtc);
        snew(ir->opts.nFreeze, ir->opts.ngfrz);
        snew(ir->opts.acceleration, ir->opts.ngacc);
        snew(ir->opts.egp_flags, ir->opts.ngener * ir->opts.ngener);
    }
    if (ir->opts.ngtc > 0)
    {
        serializer->doRealArray(ir->opts.nrdf, ir->opts.ngtc);
        serializer->doRealArray(ir->opts.ref_t, ir->opts.ngtc);
        serializer->doRealArray(ir->opts.tau_t, ir->opts.ngtc);
    }
    if (ir->opts.ngfrz > 0)
    {
        serializer->doIvecArray(ir->opts.nFreeze, ir->opts.ngfrz);
    }
    if (ir->opts.ngacc > 0)
    {
        serializer->doRvecArray(ir->opts.acceleration, ir->opts.ngacc);
    }
    if (serializer->reading())
    {
        ir->useConstantAcceleration = false;
        for (int g = 0; g < ir->opts.ngacc; g++)
        {
            if (norm2(ir->opts.acceleration[g]) != 0)
            {
                ir->useConstantAcceleration = true;
            }
        }
    }
    serializer->doIntArray(ir->opts.egp_flags, ir->opts.ngener * ir->opts.ngener);

    /* First read the lists with annealing and npoints for each group */
    serializer->doEnumArrayAsInt(ir->opts.annealing, ir->opts.ngtc);
    serializer->doIntArray(ir->opts.anneal_npoints, ir->opts.ngtc);
    for (j = 0; j < (ir->opts.ngtc); j++)
    {
        k = ir->opts.anneal_npoints[j];
        if (serializer->reading())
        {
            snew(ir->opts.anneal_time[j], k);
            snew(ir->opts.anneal_temp[j], k);
        }
        serializer->doRealArray(ir->opts.anneal_time[j], k);
        serializer->doRealArray(ir->opts.anneal_temp[j], k);
    }
    /* Walls */
    {
        serializer->doInt(&ir->nwall);
        serializer->doEnumAsInt(&ir->wall_type);
        serializer->doReal(&ir->wall_r_linpot);
        serializer->doInt(&ir->wall_atomtype[0]);
        serializer->doInt(&ir->wall_atomtype[1]);
        serializer->doReal(&ir->wall_density[0]);
        serializer->doReal(&ir->wall_density[1]);
        serializer->doReal(&ir->wall_ewald_zfac);
    }

    /* Cosine stuff for electric fields */
    if (file_version < tpxv_GenericParamsForElectricField)
    {
        do_legacy_efield(serializer, &paramsObj);
    }

    /* Swap ions */
    if (file_version >= tpxv_ComputationalElectrophysiology)
    {
        serializer->doEnumAsInt(&ir->eSwapCoords);
        if (ir->eSwapCoords != SwapType::No)
        {
            if (serializer->reading())
            {
                snew(ir->swap, 1);
            }
            do_swapcoords_tpx(serializer, ir->swap, file_version);
        }
    }

    /* QMMM reading - despite defunct we require reading for backwards
     * compatibility and correct serialization
     */
    {
        serializer->doBool(&ir->bQMMM);
        int qmmmScheme;
        serializer->doInt(&qmmmScheme);
        real unusedScalefactor;
        serializer->doReal(&unusedScalefactor);

        // this is still used in Mimic
        serializer->doInt(&ir->opts.ngQM);
        if (ir->opts.ngQM > 0 && ir->bQMMM)
        {
            /* We leave in dummy i/o for removed parameters to avoid
             * changing the tpr format.
             */
            std::vector<int> dummyIntVec(4 * ir->opts.ngQM, 0);
            serializer->doIntArray(dummyIntVec.data(), dummyIntVec.size());
            dummyIntVec.clear();

            // std::vector<bool> has no data()
            std::vector<char> dummyBoolVec(ir->opts.ngQM * sizeof(bool) / sizeof(char));
            serializer->doBoolArray(reinterpret_cast<bool*>(dummyBoolVec.data()), dummyBoolVec.size());
            dummyBoolVec.clear();

            dummyIntVec.resize(2 * ir->opts.ngQM, 0);
            serializer->doIntArray(dummyIntVec.data(), dummyIntVec.size());
            dummyIntVec.clear();

            std::vector<real> dummyRealVec(2 * ir->opts.ngQM, 0);
            serializer->doRealArray(dummyRealVec.data(), dummyRealVec.size());
            dummyRealVec.clear();

            dummyIntVec.resize(3 * ir->opts.ngQM, 0);
            serializer->doIntArray(dummyIntVec.data(), dummyIntVec.size());
            dummyIntVec.clear();
        }
        /* end of QMMM stuff */
    }

    if (file_version >= tpxv_GenericParamsForElectricField)
    {
        if (serializer->reading())
        {
            paramsObj.mergeObject(gmx::deserializeKeyValueTree(serializer));
        }
        else
        {
            GMX_RELEASE_ASSERT(ir->params != nullptr,
                               "Parameters should be present when writing inputrec");
            gmx::serializeKeyValueTree(*ir->params, serializer);
        }
    }
    if (serializer->reading())
    {
        ir->params = new gmx::KeyValueTreeObject(paramsBuilder.build());
        // Initialize internal parameters to an empty kvt for all tpr versions
        ir->internalParameters = std::make_unique<gmx::KeyValueTreeObject>();
    }

    if (file_version >= tpxv_GenericInternalParameters)
    {
        if (serializer->reading())
        {
            ir->internalParameters =
                    std::make_unique<gmx::KeyValueTreeObject>(gmx::deserializeKeyValueTree(serializer));
        }
        else
        {
            GMX_RELEASE_ASSERT(ir->internalParameters != nullptr,
                               "Parameters should be present when writing inputrec");
            gmx::serializeKeyValueTree(*ir->internalParameters, serializer);
        }
    }

    if (file_version < tpxv_EnsembleTemperature)
    {
        if (doSimulatedAnnealing(*ir) || ir->bSimTemp)
        {
            ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::Variable;
        }
        else if (integratorHasReferenceTemperature(*ir))
        {
            ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::Constant;
        }
        else
        {
            ir->ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;
        }

        if (ir->ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant)
        {
            ir->ensembleTemperature = ir->opts.ref_t[0];
        }
        else
        {
            ir->ensembleTemperature = -1;
        }
    }
}


static void do_harm(gmx::ISerializer* serializer, t_iparams* iparams)
{
    serializer->doReal(&iparams->harmonic.rA);
    serializer->doReal(&iparams->harmonic.krA);
    serializer->doReal(&iparams->harmonic.rB);
    serializer->doReal(&iparams->harmonic.krB);
}

static void do_iparams(gmx::ISerializer* serializer, t_functype ftype, t_iparams* iparams, int file_version)
{
    int  idum;
    real rdum;

    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_IDIHS:
            do_harm(serializer, iparams);
            if ((ftype == F_ANGRES || ftype == F_ANGRESZ) && serializer->reading())
            {
                /* Correct incorrect storage of parameters */
                iparams->pdihs.phiB = iparams->pdihs.phiA;
                iparams->pdihs.cpB  = iparams->pdihs.cpA;
            }
            break;
        case F_RESTRANGLES:
            serializer->doReal(&iparams->harmonic.rA);
            serializer->doReal(&iparams->harmonic.krA);
            break;
        case F_LINEAR_ANGLES:
            serializer->doReal(&iparams->linangle.klinA);
            serializer->doReal(&iparams->linangle.aA);
            serializer->doReal(&iparams->linangle.klinB);
            serializer->doReal(&iparams->linangle.aB);
            break;
        case F_FENEBONDS:
            serializer->doReal(&iparams->fene.bm);
            serializer->doReal(&iparams->fene.kb);
            break;

        case F_RESTRBONDS:
            serializer->doReal(&iparams->restraint.lowA);
            serializer->doReal(&iparams->restraint.up1A);
            serializer->doReal(&iparams->restraint.up2A);
            serializer->doReal(&iparams->restraint.kA);
            serializer->doReal(&iparams->restraint.lowB);
            serializer->doReal(&iparams->restraint.up1B);
            serializer->doReal(&iparams->restraint.up2B);
            serializer->doReal(&iparams->restraint.kB);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            serializer->doReal(&iparams->tab.kA);
            serializer->doInt(&iparams->tab.table);
            serializer->doReal(&iparams->tab.kB);
            break;
        case F_CROSS_BOND_BONDS:
            serializer->doReal(&iparams->cross_bb.r1e);
            serializer->doReal(&iparams->cross_bb.r2e);
            serializer->doReal(&iparams->cross_bb.krr);
            break;
        case F_CROSS_BOND_ANGLES:
            serializer->doReal(&iparams->cross_ba.r1e);
            serializer->doReal(&iparams->cross_ba.r2e);
            serializer->doReal(&iparams->cross_ba.r3e);
            serializer->doReal(&iparams->cross_ba.krt);
            break;
        case F_UREY_BRADLEY:
            serializer->doReal(&iparams->u_b.thetaA);
            serializer->doReal(&iparams->u_b.kthetaA);
            serializer->doReal(&iparams->u_b.r13A);
            serializer->doReal(&iparams->u_b.kUBA);
            if (file_version >= 79)
            {
                serializer->doReal(&iparams->u_b.thetaB);
                serializer->doReal(&iparams->u_b.kthetaB);
                serializer->doReal(&iparams->u_b.r13B);
                serializer->doReal(&iparams->u_b.kUBB);
            }
            else
            {
                iparams->u_b.thetaB  = iparams->u_b.thetaA;
                iparams->u_b.kthetaB = iparams->u_b.kthetaA;
                iparams->u_b.r13B    = iparams->u_b.r13A;
                iparams->u_b.kUBB    = iparams->u_b.kUBA;
            }
            break;
        case F_QUARTIC_ANGLES:
            serializer->doReal(&iparams->qangle.theta);
            serializer->doRealArray(iparams->qangle.c, 5);
            break;
        case F_BHAM:
            serializer->doReal(&iparams->bham.a);
            serializer->doReal(&iparams->bham.b);
            serializer->doReal(&iparams->bham.c);
            break;
        case F_MORSE:
            serializer->doReal(&iparams->morse.b0A);
            serializer->doReal(&iparams->morse.cbA);
            serializer->doReal(&iparams->morse.betaA);
            if (file_version >= 79)
            {
                serializer->doReal(&iparams->morse.b0B);
                serializer->doReal(&iparams->morse.cbB);
                serializer->doReal(&iparams->morse.betaB);
            }
            else
            {
                iparams->morse.b0B   = iparams->morse.b0A;
                iparams->morse.cbB   = iparams->morse.cbA;
                iparams->morse.betaB = iparams->morse.betaA;
            }
            break;
        case F_CUBICBONDS:
            serializer->doReal(&iparams->cubic.b0);
            serializer->doReal(&iparams->cubic.kb);
            serializer->doReal(&iparams->cubic.kcub);
            break;
        case F_CONNBONDS: break;
        case F_POLARIZATION: serializer->doReal(&iparams->polarize.alpha); break;
        case F_ANHARM_POL:
            serializer->doReal(&iparams->anharm_polarize.alpha);
            serializer->doReal(&iparams->anharm_polarize.drcut);
            serializer->doReal(&iparams->anharm_polarize.khyp);
            break;
        case F_WATER_POL:
            serializer->doReal(&iparams->wpol.al_x);
            serializer->doReal(&iparams->wpol.al_y);
            serializer->doReal(&iparams->wpol.al_z);
            serializer->doReal(&iparams->wpol.rOH);
            serializer->doReal(&iparams->wpol.rHH);
            serializer->doReal(&iparams->wpol.rOD);
            break;
        case F_THOLE_POL:
            serializer->doReal(&iparams->thole.a);
            serializer->doReal(&iparams->thole.alpha1);
            serializer->doReal(&iparams->thole.alpha2);
            if (file_version < tpxv_RemoveTholeRfac)
            {
                real noRfac = 0;
                serializer->doReal(&noRfac);
            }

            break;
        case F_LJ:
            serializer->doReal(&iparams->lj.c6);
            serializer->doReal(&iparams->lj.c12);
            break;
        case F_LJ14:
            serializer->doReal(&iparams->lj14.c6A);
            serializer->doReal(&iparams->lj14.c12A);
            serializer->doReal(&iparams->lj14.c6B);
            serializer->doReal(&iparams->lj14.c12B);
            break;
        case F_LJC14_Q:
            serializer->doReal(&iparams->ljc14.fqq);
            serializer->doReal(&iparams->ljc14.qi);
            serializer->doReal(&iparams->ljc14.qj);
            serializer->doReal(&iparams->ljc14.c6);
            serializer->doReal(&iparams->ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            serializer->doReal(&iparams->ljcnb.qi);
            serializer->doReal(&iparams->ljcnb.qj);
            serializer->doReal(&iparams->ljcnb.c6);
            serializer->doReal(&iparams->ljcnb.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            serializer->doReal(&iparams->pdihs.phiA);
            serializer->doReal(&iparams->pdihs.cpA);
            serializer->doReal(&iparams->pdihs.phiB);
            serializer->doReal(&iparams->pdihs.cpB);
            serializer->doInt(&iparams->pdihs.mult);
            break;
        case F_RESTRDIHS:
            serializer->doReal(&iparams->pdihs.phiA);
            serializer->doReal(&iparams->pdihs.cpA);
            break;
        case F_DISRES:
            serializer->doInt(&iparams->disres.label);
            serializer->doInt(&iparams->disres.type);
            serializer->doReal(&iparams->disres.low);
            serializer->doReal(&iparams->disres.up1);
            serializer->doReal(&iparams->disres.up2);
            serializer->doReal(&iparams->disres.kfac);
            break;
        case F_ORIRES:
            serializer->doInt(&iparams->orires.ex);
            serializer->doInt(&iparams->orires.label);
            serializer->doInt(&iparams->orires.power);
            serializer->doReal(&iparams->orires.c);
            serializer->doReal(&iparams->orires.obs);
            serializer->doReal(&iparams->orires.kfac);
            break;
        case F_DIHRES:
            if (file_version < 82)
            {
                serializer->doInt(&idum);
                serializer->doInt(&idum);
            }
            serializer->doReal(&iparams->dihres.phiA);
            serializer->doReal(&iparams->dihres.dphiA);
            serializer->doReal(&iparams->dihres.kfacA);
            if (file_version >= 82)
            {
                serializer->doReal(&iparams->dihres.phiB);
                serializer->doReal(&iparams->dihres.dphiB);
                serializer->doReal(&iparams->dihres.kfacB);
            }
            else
            {
                iparams->dihres.phiB  = iparams->dihres.phiA;
                iparams->dihres.dphiB = iparams->dihres.dphiA;
                iparams->dihres.kfacB = iparams->dihres.kfacA;
            }
            break;
        case F_POSRES:
            serializer->doRvec(&iparams->posres.pos0A);
            serializer->doRvec(&iparams->posres.fcA);
            serializer->doRvec(&iparams->posres.pos0B);
            serializer->doRvec(&iparams->posres.fcB);
            break;
        case F_FBPOSRES:
            serializer->doInt(&iparams->fbposres.geom);
            serializer->doRvec(&iparams->fbposres.pos0);
            serializer->doReal(&iparams->fbposres.r);
            serializer->doReal(&iparams->fbposres.k);
            break;
        case F_CBTDIHS: serializer->doRealArray(iparams->cbtdihs.cbtcA, NR_CBTDIHS); break;
        case F_RBDIHS:
            // Fall-through intended
        case F_FOURDIHS:
            /* Fourier dihedrals are internally represented
             * as Ryckaert-Bellemans since those are faster to compute.
             */
            serializer->doRealArray(iparams->rbdihs.rbcA, NR_RBDIHS);
            serializer->doRealArray(iparams->rbdihs.rbcB, NR_RBDIHS);
            break;
        case F_CONSTR:
        case F_CONSTRNC:
            serializer->doReal(&iparams->constr.dA);
            serializer->doReal(&iparams->constr.dB);
            break;
        case F_SETTLE:
            serializer->doReal(&iparams->settle.doh);
            serializer->doReal(&iparams->settle.dhh);
            break;
        case F_VSITE1: break; // VSite1 has 0 parameters
        case F_VSITE2:
        case F_VSITE2FD: serializer->doReal(&iparams->vsite.a); break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            serializer->doReal(&iparams->vsite.a);
            serializer->doReal(&iparams->vsite.b);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            serializer->doReal(&iparams->vsite.a);
            serializer->doReal(&iparams->vsite.b);
            serializer->doReal(&iparams->vsite.c);
            break;
        case F_VSITEN:
            serializer->doInt(&iparams->vsiten.n);
            serializer->doReal(&iparams->vsiten.a);
            break;
        case F_GB12_NOLONGERUSED:
        case F_GB13_NOLONGERUSED:
        case F_GB14_NOLONGERUSED:
            // Implicit solvent parameters can still be read, but never used
            if (serializer->reading())
            {
                if (file_version < 68)
                {
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                }
                if (file_version < tpxv_RemoveImplicitSolvation)
                {
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                    serializer->doReal(&rdum);
                }
            }
            break;
        case F_CMAP:
            serializer->doInt(&iparams->cmap.cmapA);
            serializer->doInt(&iparams->cmap.cmapB);
            break;
        default:
            gmx_fatal(FARGS,
                      "unknown function type %d (%s) in %s line %d",
                      ftype,
                      interaction_function[ftype].name,
                      __FILE__,
                      __LINE__);
    }
}

static void do_ilist(gmx::ISerializer* serializer, InteractionList* ilist)
{
    int nr = ilist->size();
    serializer->doInt(&nr);
    if (serializer->reading())
    {
        ilist->iatoms.resize(nr);
    }
    serializer->doIntArray(ilist->iatoms.data(), ilist->size());
}

static void do_ffparams(gmx::ISerializer* serializer, gmx_ffparams_t* ffparams, int file_version)
{
    serializer->doInt(&ffparams->atnr);
    int numTypes = ffparams->numTypes();
    serializer->doInt(&numTypes);
    if (serializer->reading())
    {
        ffparams->functype.resize(numTypes);
        ffparams->iparams.resize(numTypes);
    }
    /* Read/write all the function types */
    serializer->doIntArray(ffparams->functype.data(), ffparams->functype.size());

    if (file_version >= 66)
    {
        serializer->doDouble(&ffparams->reppow);
    }
    else
    {
        ffparams->reppow = 12.0;
    }

    serializer->doReal(&ffparams->fudgeQQ);

    /* Check whether all these function types are supported by the code.
     * In practice the code is backwards compatible, which means that the
     * numbering may have to be altered from old numbering to new numbering
     */
    for (int i = 0; i < ffparams->numTypes(); i++)
    {
        if (serializer->reading())
        {
            /* Loop over file versions */
            for (int k = 0; k < NFTUPD; k++)
            {
                /* Compare the read file_version to the update table */
                if ((file_version < ftupd[k].fvnr) && (ffparams->functype[i] >= ftupd[k].ftype))
                {
                    ffparams->functype[i] += 1;
                }
            }
        }

        do_iparams(serializer, ffparams->functype[i], &ffparams->iparams[i], file_version);
    }
}

static void add_settle_atoms(InteractionList* ilist)
{
    int i;

    /* Settle used to only store the first atom: add the other two */
    ilist->iatoms.resize(2 * ilist->size());
    for (i = ilist->size() / 4 - 1; i >= 0; i--)
    {
        ilist->iatoms[4 * i + 0] = ilist->iatoms[2 * i + 0];
        ilist->iatoms[4 * i + 1] = ilist->iatoms[2 * i + 1];
        ilist->iatoms[4 * i + 2] = ilist->iatoms[2 * i + 1] + 1;
        ilist->iatoms[4 * i + 3] = ilist->iatoms[2 * i + 1] + 2;
    }
}

static void do_ilists(gmx::ISerializer* serializer, InteractionLists* ilists, int file_version)
{
    GMX_RELEASE_ASSERT(ilists, "Need a valid ilists object");
    GMX_RELEASE_ASSERT(ilists->size() == F_NRE,
                       "The code needs to be in sync with InteractionLists");

    for (int j = 0; j < F_NRE; j++)
    {
        InteractionList& ilist  = (*ilists)[j];
        gmx_bool         bClear = FALSE;
        if (serializer->reading())
        {
            for (int k = 0; k < NFTUPD; k++)
            {
                if ((file_version < ftupd[k].fvnr) && (j == ftupd[k].ftype))
                {
                    bClear = TRUE;
                }
            }
        }
        if (bClear)
        {
            ilist.iatoms.clear();
        }
        else
        {
            do_ilist(serializer, &ilist);
            if (file_version < 78 && j == F_SETTLE && !ilist.empty())
            {
                add_settle_atoms(&ilist);
            }
        }
    }
}

static void do_block(gmx::ISerializer* serializer, t_block* block)
{
    serializer->doInt(&block->nr);
    if (serializer->reading())
    {
        if ((block->nalloc_index > 0) && (nullptr != block->index))
        {
            sfree(block->index);
        }
        block->nalloc_index = block->nr + 1;
        snew(block->index, block->nalloc_index);
    }
    serializer->doIntArray(block->index, block->nr + 1);
}

static void doListOfLists(gmx::ISerializer* serializer, gmx::ListOfLists<int>* listOfLists)
{
    int numLists = listOfLists->ssize();
    serializer->doInt(&numLists);
    int numElements = listOfLists->elementsView().ssize();
    serializer->doInt(&numElements);
    if (serializer->reading())
    {
        std::vector<int> listRanges(numLists + 1);
        serializer->doIntArray(listRanges.data(), numLists + 1);
        std::vector<int> elements(numElements);
        serializer->doIntArray(elements.data(), numElements);
        *listOfLists = gmx::ListOfLists<int>(std::move(listRanges), std::move(elements));
    }
    else
    {
        serializer->doIntArray(const_cast<int*>(listOfLists->listRangesView().data()), numLists + 1);
        serializer->doIntArray(const_cast<int*>(listOfLists->elementsView().data()), numElements);
    }
}

/* This is a primitive routine to make it possible to translate atomic numbers
 * to element names when reading TPR files, without making the Gromacs library
 * directory a dependency on mdrun (which is the case if we need elements.dat).
 */
static const char* atomicnumber_to_element(int atomicnumber)
{
    const char* p;

    /* This does not have to be complete, so we only include elements likely
     * to occur in PDB files.
     */
    switch (atomicnumber)
    {
        case 1: p = "H"; break;
        case 5: p = "B"; break;
        case 6: p = "C"; break;
        case 7: p = "N"; break;
        case 8: p = "O"; break;
        case 9: p = "F"; break;
        case 11: p = "Na"; break;
        case 12: p = "Mg"; break;
        case 15: p = "P"; break;
        case 16: p = "S"; break;
        case 17: p = "Cl"; break;
        case 18: p = "Ar"; break;
        case 19: p = "K"; break;
        case 20: p = "Ca"; break;
        case 25: p = "Mn"; break;
        case 26: p = "Fe"; break;
        case 28: p = "Ni"; break;
        case 29: p = "Cu"; break;
        case 30: p = "Zn"; break;
        case 35: p = "Br"; break;
        case 47: p = "Ag"; break;
        default: p = ""; break;
    }
    return p;
}


static void do_atom(gmx::ISerializer* serializer, t_atom* atom)
{
    serializer->doReal(&atom->m);
    serializer->doReal(&atom->q);
    serializer->doReal(&atom->mB);
    serializer->doReal(&atom->qB);
    serializer->doUShort(&atom->type);
    serializer->doUShort(&atom->typeB);
    serializer->doEnumAsInt(&atom->ptype);
    serializer->doInt(&atom->resind);
    serializer->doInt(&atom->atomnumber);
    if (serializer->reading())
    {
        /* Set element string from atomic number if present.
         * This routine returns an empty string if the name is not found.
         */
        std::strncpy(atom->elem, atomicnumber_to_element(atom->atomnumber), 4);
        /* avoid warnings about potentially unterminated string */
        atom->elem[3] = '\0';
    }
}

static void do_grps(gmx::ISerializer* serializer, gmx::ArrayRef<AtomGroupIndices> grps)
{
    for (auto& group : grps)
    {
        int size = group.size();
        serializer->doInt(&size);
        if (serializer->reading())
        {
            group.resize(size);
        }
        serializer->doIntArray(group.data(), size);
    }
}

static void do_symstr(gmx::ISerializer* serializer, char*** nm, t_symtab* symtab)
{
    int ls;

    if (serializer->reading())
    {
        serializer->doInt(&ls);
        *nm = get_symtab_handle(symtab, ls);
    }
    else
    {
        ls = lookup_symtab(symtab, *nm);
        serializer->doInt(&ls);
    }
}

static void do_strstr(gmx::ISerializer* serializer, int nstr, char*** nm, t_symtab* symtab)
{
    int j;

    for (j = 0; (j < nstr); j++)
    {
        do_symstr(serializer, &(nm[j]), symtab);
    }
}

static void do_resinfo(gmx::ISerializer* serializer, int n, t_resinfo* ri, t_symtab* symtab, int file_version)
{
    int j;

    for (j = 0; (j < n); j++)
    {
        do_symstr(serializer, &(ri[j].name), symtab);
        if (file_version >= 63)
        {
            serializer->doInt(&ri[j].nr);
            serializer->doUChar(&ri[j].ic);
        }
        else
        {
            ri[j].nr = j + 1;
            ri[j].ic = ' ';
        }
    }
}

static void do_atoms(gmx::ISerializer* serializer, t_atoms* atoms, t_symtab* symtab, int file_version)
{
    int i;

    serializer->doInt(&atoms->nr);
    serializer->doInt(&atoms->nres);
    if (serializer->reading())
    {
        /* Since we have always written all t_atom properties in the tpr file
         * (at least for all backward compatible versions), we don't store
         * but simple set the booleans here.
         */
        atoms->haveMass    = TRUE;
        atoms->haveCharge  = TRUE;
        atoms->haveType    = TRUE;
        atoms->haveBState  = TRUE;
        atoms->havePdbInfo = FALSE;

        snew(atoms->atom, atoms->nr);
        snew(atoms->atomname, atoms->nr);
        snew(atoms->atomtype, atoms->nr);
        snew(atoms->atomtypeB, atoms->nr);
        snew(atoms->resinfo, atoms->nres);
        atoms->pdbinfo = nullptr;
    }
    else
    {
        GMX_RELEASE_ASSERT(atoms->haveMass && atoms->haveCharge && atoms->haveType && atoms->haveBState,
                           "Mass, charge, atomtype and B-state parameters should be present in "
                           "t_atoms when writing a tpr file");
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        do_atom(serializer, &atoms->atom[i]);
    }
    do_strstr(serializer, atoms->nr, atoms->atomname, symtab);
    do_strstr(serializer, atoms->nr, atoms->atomtype, symtab);
    do_strstr(serializer, atoms->nr, atoms->atomtypeB, symtab);

    do_resinfo(serializer, atoms->nres, atoms->resinfo, symtab, file_version);
}

static void do_groups(gmx::ISerializer* serializer, SimulationGroups* groups, t_symtab* symtab)
{
    do_grps(serializer, groups->groups);
    int numberOfGroupNames = groups->groupNames.size();
    serializer->doInt(&numberOfGroupNames);
    if (serializer->reading())
    {
        groups->groupNames.resize(numberOfGroupNames);
    }
    do_strstr(serializer, numberOfGroupNames, groups->groupNames.data(), symtab);
    for (auto group : gmx::keysOf(groups->groupNumbers))
    {
        int numberOfGroupNumbers = groups->numberOfGroupNumbers(group);
        serializer->doInt(&numberOfGroupNumbers);
        if (numberOfGroupNumbers != 0)
        {
            if (serializer->reading())
            {
                groups->groupNumbers[group].resize(numberOfGroupNumbers);
            }
            serializer->doUCharArray(groups->groupNumbers[group].data(), numberOfGroupNumbers);
        }
    }
}

static void do_atomtypes(gmx::ISerializer* serializer, int file_version)
{
    int nr;
    serializer->doInt(&nr);
    if (serializer->reading() && file_version < tpxv_RemoveImplicitSolvation)
    {
        std::vector<real> dummy(nr, 0);
        serializer->doRealArray(dummy.data(), dummy.size());
        serializer->doRealArray(dummy.data(), dummy.size());
        serializer->doRealArray(dummy.data(), dummy.size());
    }
    std::vector<int> atomnumbers(nr);
    serializer->doIntArray(atomnumbers.data(), atomnumbers.size());

    if (serializer->reading() && file_version >= 60 && file_version < tpxv_RemoveImplicitSolvation)
    {
        std::vector<real> dummy(nr, 0);
        serializer->doRealArray(dummy.data(), dummy.size());
        serializer->doRealArray(dummy.data(), dummy.size());
    }
}

static void do_symtab(gmx::ISerializer* serializer, t_symtab* symtab)
{
    int       i, nr;
    t_symbuf* symbuf;

    serializer->doInt(&symtab->nr);
    nr = symtab->nr;
    if (serializer->reading())
    {
        snew(symtab->symbuf, 1);
        symbuf          = symtab->symbuf;
        symbuf->bufsize = nr;
        snew(symbuf->buf, nr);
        for (i = 0; (i < nr); i++)
        {
            std::string buf;
            serializer->doString(&buf);
            symbuf->buf[i] = gmx_strdup(buf.c_str());
        }
    }
    else
    {
        symbuf = symtab->symbuf;
        while (symbuf != nullptr)
        {
            for (i = 0; (i < symbuf->bufsize) && (i < nr); i++)
            {
                std::string buf = symbuf->buf[i];
                serializer->doString(&buf);
            }
            nr -= i;
            symbuf = symbuf->next;
        }
        if (nr != 0)
        {
            gmx_fatal(FARGS, "nr of symtab strings left: %d", nr);
        }
    }
}

static void do_cmap(gmx::ISerializer* serializer, gmx_cmap_t* cmap_grid)
{

    int ngrid = cmap_grid->cmapdata.size();
    serializer->doInt(&ngrid);
    serializer->doInt(&cmap_grid->grid_spacing);

    int gs    = cmap_grid->grid_spacing;
    int nelem = gs * gs;

    if (serializer->reading())
    {
        cmap_grid->cmapdata.resize(ngrid);

        for (int i = 0; i < ngrid; i++)
        {
            cmap_grid->cmapdata[i].cmap.resize(4 * nelem);
        }
    }

    for (int i = 0; i < ngrid; i++)
    {
        for (int j = 0; j < nelem; j++)
        {
            serializer->doReal(&cmap_grid->cmapdata[i].cmap[j * 4]);
            serializer->doReal(&cmap_grid->cmapdata[i].cmap[j * 4 + 1]);
            serializer->doReal(&cmap_grid->cmapdata[i].cmap[j * 4 + 2]);
            serializer->doReal(&cmap_grid->cmapdata[i].cmap[j * 4 + 3]);
        }
    }
}


static void do_moltype(gmx::ISerializer* serializer, gmx_moltype_t* molt, t_symtab* symtab, int file_version)
{
    do_symstr(serializer, &(molt->name), symtab);

    do_atoms(serializer, &molt->atoms, symtab, file_version);

    do_ilists(serializer, &molt->ilist, file_version);

    /* TODO: Remove the obsolete charge group index from the file */
    t_block cgs;
    cgs.nr           = molt->atoms.nr;
    cgs.nalloc_index = cgs.nr + 1;
    snew(cgs.index, cgs.nalloc_index);
    for (int i = 0; i < cgs.nr + 1; i++)
    {
        cgs.index[i] = i;
    }
    do_block(serializer, &cgs);
    sfree(cgs.index);

    /* This used to be in the atoms struct */
    doListOfLists(serializer, &molt->excls);
}

static void do_molblock(gmx::ISerializer* serializer, gmx_molblock_t* molb, int numAtomsPerMolecule)
{
    serializer->doInt(&molb->type);
    serializer->doInt(&molb->nmol);
    /* To maintain forward topology reading compatibility, we store #atoms.
     * TODO: Change this to conditional reading of a dummy int when we
     *       increase tpx_generation.
     */
    serializer->doInt(&numAtomsPerMolecule);
    /* Position restraint coordinates */
    int numPosres_xA = molb->posres_xA.size();
    serializer->doInt(&numPosres_xA);
    if (numPosres_xA > 0)
    {
        if (serializer->reading())
        {
            molb->posres_xA.resize(numPosres_xA);
        }
        serializer->doRvecArray(as_rvec_array(molb->posres_xA.data()), numPosres_xA);
    }
    int numPosres_xB = molb->posres_xB.size();
    serializer->doInt(&numPosres_xB);
    if (numPosres_xB > 0)
    {
        if (serializer->reading())
        {
            molb->posres_xB.resize(numPosres_xB);
        }
        serializer->doRvecArray(as_rvec_array(molb->posres_xB.data()), numPosres_xB);
    }
}

static void set_disres_npair(gmx_mtop_t* mtop)
{
    gmx::ArrayRef<t_iparams> ip = mtop->ffparams.iparams;

    for (const auto ilist : IListRange(*mtop))
    {
        const InteractionList& il = ilist.list()[F_DISRES];

        if (!il.empty())
        {
            gmx::ArrayRef<const int> a     = il.iatoms;
            int                      npair = 0;
            for (int i = 0; i < il.size(); i += 3)
            {
                npair++;
                if (i + 3 == il.size() || ip[a[i]].disres.label != ip[a[i + 3]].disres.label)
                {
                    ip[a[i]].disres.npair = npair;
                    npair                 = 0;
                }
            }
        }
    }
}

static void do_mtop(gmx::ISerializer* serializer, gmx_mtop_t* mtop, int file_version)
{
    do_symtab(serializer, &(mtop->symtab));

    do_symstr(serializer, &(mtop->name), &(mtop->symtab));

    do_ffparams(serializer, &mtop->ffparams, file_version);

    int nmoltype = mtop->moltype.size();
    serializer->doInt(&nmoltype);
    if (serializer->reading())
    {
        mtop->moltype.resize(nmoltype);
    }
    for (gmx_moltype_t& moltype : mtop->moltype)
    {
        do_moltype(serializer, &moltype, &mtop->symtab, file_version);
    }

    int nmolblock = mtop->molblock.size();
    serializer->doInt(&nmolblock);
    if (serializer->reading())
    {
        mtop->molblock.resize(nmolblock);
    }
    for (gmx_molblock_t& molblock : mtop->molblock)
    {
        int numAtomsPerMolecule = (serializer->reading() ? 0 : mtop->moltype[molblock.type].atoms.nr);
        do_molblock(serializer, &molblock, numAtomsPerMolecule);
    }
    serializer->doInt(&mtop->natoms);

    if (file_version >= tpxv_IntermolecularBondeds)
    {
        serializer->doBool(&mtop->bIntermolecularInteractions);
        if (mtop->bIntermolecularInteractions)
        {
            if (serializer->reading())
            {
                mtop->intermolecular_ilist = std::make_unique<InteractionLists>();
            }
            do_ilists(serializer, mtop->intermolecular_ilist.get(), file_version);
        }
    }
    else
    {
        mtop->bIntermolecularInteractions = FALSE;
    }

    if (file_version < tpxv_RemoveAtomtypes)
    {
        do_atomtypes(serializer, file_version);
    }


    if (file_version >= 65)
    {
        do_cmap(serializer, &mtop->ffparams.cmap_grid);
    }
    else
    {
        mtop->ffparams.cmap_grid.grid_spacing = 0;
        mtop->ffparams.cmap_grid.cmapdata.clear();
    }

    do_groups(serializer, &mtop->groups, &(mtop->symtab));

    mtop->haveMoleculeIndices = true;

    if (file_version >= tpxv_StoreNonBondedInteractionExclusionGroup)
    {
        std::int64_t intermolecularExclusionGroupSize = gmx::ssize(mtop->intermolecularExclusionGroup);
        serializer->doInt64(&intermolecularExclusionGroupSize);
        GMX_RELEASE_ASSERT(intermolecularExclusionGroupSize >= 0,
                           "Number of atoms with excluded intermolecular non-bonded interactions "
                           "is negative.");
        mtop->intermolecularExclusionGroup.resize(intermolecularExclusionGroupSize); // no effect when writing
        serializer->doIntArray(mtop->intermolecularExclusionGroup.data(),
                               mtop->intermolecularExclusionGroup.size());
    }

    if (serializer->reading())
    {
        close_symtab(&(mtop->symtab));
    }
}

/*! \brief
 * Read the first part of the TPR file to find general system information.
 *
 * If \p TopOnlyOK is true then we can read even future versions
 * of tpx files, provided the \p fileGeneration hasn't changed.
 * If it is false, we need the \p ir too, and bail out
 * if the file is newer than the program.
 *
 * The version and generation of the topology (see top of this file)
 * are returned in the two last arguments, if those arguments are non-nullptr.
 *
 * If possible, we will read the \p ir even when \p TopOnlyOK is true.
 *
 * \param[in,out] serializer The serializer used to handle header processing.
 * \param[in,out] tpx File header datastructure.
 * \param[in]     filename The name of the file being read/written
 * \param[in,out] fio File handle.
 * \param[in] TopOnlyOK If not reading \p ir is fine or not.
 */
static void do_tpxheader(gmx::FileIOXdrSerializer*    serializer,
                         TpxFileHeader*               tpx,
                         const std::filesystem::path& filename,
                         t_fileio*                    fio,
                         bool                         TopOnlyOK)
{
    int  precision;
    int  idum = 0;
    real rdum = 0;

    /* XDR binary topology file */
    precision = sizeof(real);
    std::string buf;
    std::string fileTag;
    if (serializer->reading())
    {
        serializer->doString(&buf);
        if (std::strncmp(buf.c_str(), "VERSION", 7) != 0)
        {
            gmx_fatal(
                    FARGS,
                    "Can not read file %s,\n"
                    "             this file is from a GROMACS version which is older than 2.0\n"
                    "             Make a new one with grompp or use a gro or pdb file, if possible",
                    filename.u8string().c_str());
        }
        // We need to know the precision used to write the TPR file, to match it
        // to the precision of the currently running binary. If the precisions match
        // there is no problem, but mismatching precision needs to be accounted for
        // by reading into temporary variables of the correct precision instead
        // of the desired target datastructures.
        serializer->doInt(&precision);
        tpx->isDouble = (precision == sizeof(double));
        if ((precision != sizeof(float)) && !tpx->isDouble)
        {
            gmx_fatal(FARGS,
                      "Unknown precision in file %s: real is %d bytes "
                      "instead of %zu or %zu",
                      filename.u8string().c_str(),
                      precision,
                      sizeof(float),
                      sizeof(double));
        }
        gmx_fio_setprecision(fio, tpx->isDouble);
        fprintf(stderr,
                "Reading file %s, %s (%s precision)\n",
                filename.u8string().c_str(),
                buf.c_str(),
                tpx->isDouble ? "double" : "single");
    }
    else
    {
        buf = gmx::formatString("VERSION %s", gmx_version());
        serializer->doString(&buf);
        gmx_fio_setprecision(fio, tpx->isDouble);
        serializer->doInt(&precision);
        fileTag = tpx_tag;
    }

    /* Check versions! */
    serializer->doInt(&tpx->fileVersion);

    /* This is for backward compatibility with development versions 77-79
     * where the tag was, mistakenly, placed before the generation,
     * which would cause a segv instead of a proper error message
     * when reading the topology only from tpx with <77 code.
     */
    if (tpx->fileVersion >= 77 && tpx->fileVersion <= 79)
    {
        serializer->doString(&fileTag);
    }

    serializer->doInt(&tpx->fileGeneration);

    if (tpx->fileVersion >= 81)
    {
        serializer->doString(&fileTag);
    }
    if (serializer->reading())
    {
        if (tpx->fileVersion < 77)
        {
            /* Versions before 77 don't have the tag, set it to release */
            fileTag = gmx::formatString("%s", TPX_TAG_RELEASE);
        }

        if (fileTag != tpx_tag)
        {
            fprintf(stderr, "Note: file tpx tag '%s', software tpx tag '%s'\n", fileTag.c_str(), tpx_tag.c_str());

            /* We only support reading tpx files with the same tag as the code
             * or tpx files with the release tag and with lower version number.
             */
            if (fileTag != TPX_TAG_RELEASE && tpx->fileVersion < tpx_version)
            {
                gmx_fatal(FARGS,
                          "tpx tag/version mismatch: reading tpx file (%s) version %d, tag '%s' "
                          "with program for tpx version %d, tag '%s'",
                          filename.u8string().c_str(),
                          tpx->fileVersion,
                          fileTag.c_str(),
                          tpx_version,
                          tpx_tag.c_str());
            }
        }
    }

    if ((tpx->fileVersion <= tpx_incompatible_version)
        || ((tpx->fileVersion > tpx_version) && !TopOnlyOK) || (tpx->fileGeneration > tpx_generation)
        || tpx_version == 80) /*80 was used by both 5.0-dev and 4.6-dev*/
    {
        gmx_fatal(FARGS,
                  "reading tpx file (%s) version %d with version %d program",
                  filename.u8string().c_str(),
                  tpx->fileVersion,
                  tpx_version);
    }

    serializer->doInt(&tpx->natoms);
    serializer->doInt(&tpx->ngtc);

    if (tpx->fileVersion < 62)
    {
        serializer->doInt(&idum);
        serializer->doReal(&rdum);
    }
    if (tpx->fileVersion >= 79)
    {
        serializer->doInt(&tpx->fep_state);
    }
    serializer->doReal(&tpx->lambda);
    serializer->doBool(&tpx->bIr);
    serializer->doBool(&tpx->bTop);
    serializer->doBool(&tpx->bX);
    serializer->doBool(&tpx->bV);
    serializer->doBool(&tpx->bF);
    serializer->doBool(&tpx->bBox);

    if (tpx->fileVersion >= tpxv_AddSizeField && tpx->fileGeneration >= 27)
    {
        if (!serializer->reading())
        {
            GMX_RELEASE_ASSERT(tpx->sizeOfTprBody != 0,
                               "Not possible to write new file with zero TPR body size");
        }
        serializer->doInt64(&tpx->sizeOfTprBody);
    }

    if ((tpx->fileGeneration > tpx_generation))
    {
        /* This can only happen if TopOnlyOK=TRUE */
        tpx->bIr = FALSE;
    }
}

#define do_test(serializer, b, p)                            \
    if ((serializer)->reading() && ((p) != nullptr) && !(b)) \
    gmx_fatal(FARGS, "No %s in input file", #p)

/*! \brief
 * Process the first part of the TPR into the state datastructure.
 *
 * Due to the structure of the legacy code, it is necessary
 * to split up the state reading into two parts, with the
 * box and legacy temperature coupling processed before the
 * topology datastructures.
 *
 * See the documentation for do_tpx_body for the correct order of
 * the operations for reading a tpr file.
 *
 * \param[in] serializer Abstract serializer used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in, out] state Global state data.
 */
static void do_tpx_state_first(gmx::ISerializer* serializer, TpxFileHeader* tpx, t_state* state)
{
    if (serializer->reading())
    {
        state->setFlags(0);
        init_gtc_state(state, tpx->ngtc, 0, 0);
    }
    do_test(serializer, tpx->bBox, state->box);
    if (tpx->bBox)
    {
        serializer->doRvecArray(state->box, DIM);
        if (tpx->fileVersion >= 51)
        {
            serializer->doRvecArray(state->box_rel, DIM);
        }
        else
        {
            /* We initialize box_rel after reading the inputrec */
            clear_mat(state->box_rel);
        }
        serializer->doRvecArray(state->boxv, DIM);
        if (tpx->fileVersion < 56)
        {
            matrix mdum;
            serializer->doRvecArray(mdum, DIM);
        }
    }

    if (state->ngtc > 0)
    {
        real* dumv;
        snew(dumv, state->ngtc);
        if (tpx->fileVersion < 69)
        {
            serializer->doRealArray(dumv, state->ngtc);
        }
        /* These used to be the Berendsen tcoupl_lambda's */
        serializer->doRealArray(dumv, state->ngtc);
        sfree(dumv);
    }
}

/*! \brief
 * Process global topology data.
 *
 * See the documentation for do_tpx_body for the correct order of
 * the operations for reading a tpr file.
 *
 * \param[in] serializer Abstract serializer  used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in,out] mtop Global topology.
 */
static void do_tpx_mtop(gmx::ISerializer* serializer, TpxFileHeader* tpx, gmx_mtop_t* mtop)
{
    do_test(serializer, tpx->bTop, mtop);
    if (tpx->bTop)
    {
        if (mtop)
        {
            do_mtop(serializer, mtop, tpx->fileVersion);
            set_disres_npair(mtop);
            mtop->finalize();
        }
        else
        {
            gmx_mtop_t dum_top;
            do_mtop(serializer, &dum_top, tpx->fileVersion);
        }
    }
}
/*! \brief
 * Process coordinate vectors for state data.
 *
 * Main part of state gets processed here.
 *
 * See the documentation for do_tpx_body for the correct order of
 * the operations for reading a tpr file.
 *
 * \param[in] serializer Abstract serializer used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in,out] state Global state data.
 * \param[in,out] x Individual coordinates for processing, deprecated.
 * \param[in,out] v Individual velocities for processing, deprecated.
 */
static void do_tpx_state_second(gmx::ISerializer* serializer, TpxFileHeader* tpx, t_state* state, rvec* x, rvec* v)
{
    if (!serializer->reading())
    {
        GMX_RELEASE_ASSERT(
                x == nullptr && v == nullptr,
                "Passing separate x and v pointers to do_tpx() is not supported when writing");
    }
    else
    {
        GMX_RELEASE_ASSERT(!(x == nullptr && v != nullptr),
                           "Passing x==NULL and v!=NULL is not supported");
    }

    if (serializer->reading())
    {
        if (x == nullptr)
        {
            // v is also nullptr by the above assertion, so we may
            // need to make memory in state for storing the contents
            // of the tpx file.
            if (tpx->bX)
            {
                state->addEntry(StateEntry::X);
            }
            if (tpx->bV)
            {
                state->addEntry(StateEntry::V);
            }
            state->changeNumAtoms(tpx->natoms);
        }
    }

    if (x == nullptr)
    {
        x = state->x.rvec_array();
        v = state->v.rvec_array();
    }
    do_test(serializer, tpx->bX, x);
    if (tpx->bX)
    {
        if (serializer->reading())
        {
            state->addEntry(StateEntry::X);
        }
        serializer->doRvecArray(x, tpx->natoms);
    }

    // We cannot call do_test() with v as some "integrators" don't use v
    if (tpx->bV)
    {
        if (serializer->reading())
        {
            state->addEntry(StateEntry::V);
        }
        if (!v)
        {
            std::vector<gmx::RVec> dummyVelocities(tpx->natoms);
            serializer->doRvecArray(as_rvec_array(dummyVelocities.data()), tpx->natoms);
        }
        else
        {
            serializer->doRvecArray(v, tpx->natoms);
        }
    }
    else if (v)
    {
        // Velocities are not present in the tpr file, but v has been passed here.
        // We clear v for backward compatibility.
        for (int i = 0; i < tpx->natoms; i++)
        {
            clear_rvec(v[i]);
        }
    }

    // No need to run do_test when the last argument is NULL
    if (tpx->bF)
    {
        std::vector<gmx::RVec> dummyForces(state->numAtoms());
        serializer->doRvecArray(as_rvec_array(dummyForces.data()), tpx->natoms);
    }
}
/*! \brief
 * Process simulation parameters.
 *
 * See the documentation for do_tpx_body for the correct order of
 * the operations for reading a tpr file.
 *
 * \param[in] serializer Abstract serializer used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in,out] ir Datastructure with simulation parameters.
 */
static PbcType do_tpx_ir(gmx::ISerializer* serializer, TpxFileHeader* tpx, t_inputrec* ir)
{
    PbcType  pbcType;
    gmx_bool bPeriodicMols;

    /* Starting with tpx version 26, we have the inputrec
     * at the end of the file, so we can ignore it
     * if the file is never than the software (but still the
     * same generation - see comments at the top of this file.
     *
     *
     */
    pbcType       = PbcType::Unset;
    bPeriodicMols = FALSE;

    do_test(serializer, tpx->bIr, ir);
    if (tpx->bIr)
    {
        if (tpx->fileVersion >= 53)
        {
            /* Removed the pbc info from do_inputrec, since we always want it */
            if (!serializer->reading())
            {
                pbcType       = ir->pbcType;
                bPeriodicMols = ir->bPeriodicMols;
            }
            serializer->doInt(reinterpret_cast<int*>(&pbcType));
            serializer->doBool(&bPeriodicMols);
        }
        if (tpx->fileGeneration <= tpx_generation && ir)
        {
            do_inputrec(serializer, ir, tpx->fileVersion);
            if (tpx->fileVersion < 53)
            {
                pbcType       = ir->pbcType;
                bPeriodicMols = ir->bPeriodicMols;
            }
        }
        if (serializer->reading() && ir && tpx->fileVersion >= 53)
        {
            /* We need to do this after do_inputrec, since that initializes ir */
            ir->pbcType       = pbcType;
            ir->bPeriodicMols = bPeriodicMols;
        }
    }
    return pbcType;
}

/*! \brief
 * Correct and finalize read information.
 *
 * If \p state is nullptr, skip the parts dependent on it.
 *
 * See the documentation for do_tpx_body for the correct order of
 * the operations for reading a tpr file.
 *
 * \param[in] tpx The file header used to check version numbers.
 * \param[out] ir Input rec that needs correction.
 * \param[out] state State needing correction.
 * \param[out] mtop Topology to finalize.
 */
static void do_tpx_finalize(TpxFileHeader* tpx, t_inputrec* ir, t_state* state, gmx_mtop_t* mtop)
{
    if (tpx->fileVersion < 51 && state)
    {
        set_box_rel(ir, state);
    }
    if (tpx->bIr && ir)
    {
        if (state && state->ngtc == 0)
        {
            /* Reading old version without tcoupl state data: set it */
            init_gtc_state(state, ir->opts.ngtc, 0, ir->opts.nhchainlength);
        }
        if (tpx->bTop && mtop)
        {
            if (tpx->fileVersion < 57)
            {
                ir->eDisre = !mtop->moltype[0].ilist[F_DISRES].empty()
                                     ? DistanceRestraintRefinement::Simple
                                     : DistanceRestraintRefinement::None;
            }
        }
    }
}

/*! \brief
 * Process TPR data for file reading/writing.
 *
 * The TPR file gets processed in in four stages due to the organization
 * of the data within it.
 *
 * First, state data for the box is processed in do_tpx_state_first.
 * This is followed by processing the topology in do_tpx_mtop.
 * Coordinate and velocity vectors are handled next in do_tpx_state_second.
 * The last file information processed is the collection of simulation parameters in do_tpx_ir.
 * When reading, a final processing step is undertaken at the end.
 *
 * \param[in] serializer Abstract serializer used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in,out] ir Datastructures with simulation parameters.
 * \param[in,out] state Global state data.
 * \param[in,out] x Individual coordinates for processing, deprecated.
 * \param[in,out] v Individual velocities for processing, deprecated.
 * \param[in,out] mtop Global topology.
 */
static PbcType do_tpx_body(gmx::ISerializer* serializer,
                           TpxFileHeader*    tpx,
                           t_inputrec*       ir,
                           t_state*          state,
                           rvec*             x,
                           rvec*             v,
                           gmx_mtop_t*       mtop)
{
    if (state)
    {
        do_tpx_state_first(serializer, tpx, state);
    }
    do_tpx_mtop(serializer, tpx, mtop);
    if (state)
    {
        do_tpx_state_second(serializer, tpx, state, x, v);
    }
    PbcType pbcType = do_tpx_ir(serializer, tpx, ir);
    if (serializer->reading())
    {
        do_tpx_finalize(tpx, ir, state, mtop);
    }
    return pbcType;
}

/*! \brief
 * Overload for do_tpx_body that defaults to state vectors being nullptr.
 *
 * \param[in] serializer Abstract serializer used to read/write data.
 * \param[in] tpx The file header data.
 * \param[in,out] ir Datastructures with simulation parameters.
 * \param[in,out] mtop Global topology.
 */
static PbcType do_tpx_body(gmx::ISerializer* serializer, TpxFileHeader* tpx, t_inputrec* ir, gmx_mtop_t* mtop)
{
    return do_tpx_body(serializer, tpx, ir, nullptr, nullptr, nullptr, mtop);
}

static t_fileio* open_tpx(const std::filesystem::path& fn, const char* mode)
{
    return gmx_fio_open(fn, mode);
}

static void close_tpx(t_fileio* fio)
{
    gmx_fio_close(fio);
}

/*! \brief
 * Fill information into the header only from state before writing.
 *
 * Populating the header needs to be independent from writing the information
 * to file to allow things like writing the raw byte stream.
 *
 * \param[in] state The current simulation state. Can't write without it.
 * \param[in] ir Parameter and system information.
 * \param[in] mtop Global topology.
 * \returns Fully populated header.
 */
static TpxFileHeader populateTpxHeader(const t_state& state, const t_inputrec* ir, const gmx_mtop_t* mtop)
{
    TpxFileHeader header;
    header.natoms         = state.numAtoms();
    header.ngtc           = state.ngtc;
    header.fep_state      = state.fep_state;
    header.lambda         = state.lambda[FreeEnergyPerturbationCouplingType::Fep];
    header.bIr            = ir != nullptr;
    header.bTop           = mtop != nullptr;
    header.bX             = state.hasEntry(StateEntry::X);
    header.bV             = state.hasEntry(StateEntry::V);
    header.bF             = false;
    header.bBox           = true;
    header.fileVersion    = tpx_version;
    header.fileGeneration = tpx_generation;
    header.isDouble       = (sizeof(real) == sizeof(double));

    return header;
}

/*! \brief
 * Process the body of a TPR file as an opaque data buffer.
 *
 * Reads/writes the information in \p buffer from/to the \p serializer
 * provided to the function. Does not interact with the actual
 * TPR datastructures but with an in memory representation of the
 * data, so that this data can be efficiently read or written from/to
 * an original source.
 *
 * \param[in] serializer The abstract serializer used for reading or writing
 *                       the information in \p buffer.
 * \param[in,out] buffer Information from TPR file as char buffer.
 */
static void doTpxBodyBuffer(gmx::ISerializer* serializer, gmx::ArrayRef<char> buffer)
{
    serializer->doOpaque(buffer.data(), buffer.size());
}

/*! \brief
 * Populates simulation datastructures.
 *
 * Here the information from the serialization interface \p serializer
 * is used to first populate the datastructures containing the simulation
 * information. Depending on the version found in the header \p tpx,
 * this is done using the new reading of the data as one block from disk,
 * followed by complete deserialization of the information read from there.
 * Otherwise, the datastructures are populated as before one by one from disk.
 * The second version is the default for the legacy tools that read the
 * coordinates and velocities separate from the state.
 *
 * After reading in the data, a separate buffer is populated from them
 * containing only \p ir and \p mtop that can be communicated directly
 * to nodes needing the information to set up a simulation.
 *
 * \param[in] tpx The file header.
 * \param[in] serializer The Serialization interface used to read the TPR.
 * \param[out] ir Input rec to populate.
 * \param[out] state State vectors to populate.
 * \param[out] x Coordinates to populate if needed.
 * \param[out] v Velocities to populate if needed.
 * \param[out] mtop Global topology to populate.
 *
 * \returns Partial de-serialized TPR used for communication to nodes.
 */
static PartialDeserializedTprFile readTpxBody(TpxFileHeader*    tpx,
                                              gmx::ISerializer* serializer,
                                              t_inputrec*       ir,
                                              t_state*          state,
                                              rvec*             x,
                                              rvec*             v,
                                              gmx_mtop_t*       mtop)
{
    PartialDeserializedTprFile partialDeserializedTpr;
    if (tpx->fileVersion >= tpxv_AddSizeField && tpx->fileGeneration >= 27)
    {
        partialDeserializedTpr.body.resize(tpx->sizeOfTprBody);
        partialDeserializedTpr.header = *tpx;
        doTpxBodyBuffer(serializer, partialDeserializedTpr.body);

        partialDeserializedTpr.pbcType =
                completeTprDeserialization(&partialDeserializedTpr, ir, state, x, v, mtop);
    }
    else
    {
        partialDeserializedTpr.pbcType = do_tpx_body(serializer, tpx, ir, state, x, v, mtop);
    }
    // Update header to system info for communication to nodes.
    // As we only need to communicate the inputrec and mtop to other nodes,
    // we prepare a new char buffer with the information we have already read
    // in on main.
    partialDeserializedTpr.header = populateTpxHeader(*state, ir, mtop);
    // Long-term we should move to use little endian in files to avoid extra byte swapping,
    // but since we just used the default XDR format (which is big endian) for the TPR
    // header it would cause third-party libraries reading our raw data to tear their hair
    // if we swap the endian in the middle of the file, so we stick to big endian in the
    // TPR file for now - and thus we ask the serializer to swap if this host is little endian.
    gmx::InMemorySerializer tprBodySerializer(gmx::EndianSwapBehavior::SwapIfHostIsLittleEndian);
    do_tpx_body(&tprBodySerializer, &partialDeserializedTpr.header, ir, mtop);
    partialDeserializedTpr.body = tprBodySerializer.finishAndGetBuffer();

    return partialDeserializedTpr;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

TpxFileHeader readTpxHeader(const std::filesystem::path& fileName, bool canReadTopologyOnly)
{
    t_fileio* fio;

    fio = open_tpx(fileName, "r");
    gmx::FileIOXdrSerializer serializer(fio);

    TpxFileHeader tpx;
    do_tpxheader(&serializer, &tpx, fileName, fio, canReadTopologyOnly);
    close_tpx(fio);
    return tpx;
}

void write_tpx_state(const std::filesystem::path& fn,
                     const t_inputrec*            ir,
                     const t_state*               state,
                     const gmx_mtop_t&            mtop)
{
    /* To write a state, we first need to write the state information to a buffer before
     * we append the raw bytes to the file. For this, the header information needs to be
     * populated before we write the main body because it has some information that is
     * otherwise not available.
     */

    t_fileio* fio;

    TpxFileHeader tpx = populateTpxHeader(*state, ir, &mtop);
    // Long-term we should move to use little endian in files to avoid extra byte swapping,
    // but since we just used the default XDR format (which is big endian) for the TPR
    // header it would cause third-party libraries reading our raw data to tear their hair
    // if we swap the endian in the middle of the file, so we stick to big endian in the
    // TPR file for now - and thus we ask the serializer to swap if this host is little endian.
    gmx::InMemorySerializer tprBodySerializer(gmx::EndianSwapBehavior::SwapIfHostIsLittleEndian);

    do_tpx_body(&tprBodySerializer,
                &tpx,
                const_cast<t_inputrec*>(ir),
                const_cast<t_state*>(state),
                nullptr,
                nullptr,
                const_cast<gmx_mtop_t*>(&mtop));

    std::vector<char> tprBody = tprBodySerializer.finishAndGetBuffer();
    tpx.sizeOfTprBody         = tprBody.size();

    fio = open_tpx(fn, "w");
    gmx::FileIOXdrSerializer serializer(fio);
    do_tpxheader(&serializer, &tpx, fn, fio, ir == nullptr);
    doTpxBodyBuffer(&serializer, tprBody);

    close_tpx(fio);
}

PbcType completeTprDeserialization(PartialDeserializedTprFile* partialDeserializedTpr,
                                   t_inputrec*                 ir,
                                   t_state*                    state,
                                   rvec*                       x,
                                   rvec*                       v,
                                   gmx_mtop_t*                 mtop)
{
    // Long-term we should move to use little endian in files to avoid extra byte swapping,
    // but since we just used the default XDR format (which is big endian) for the TPR
    // header it would cause third-party libraries reading our raw data to tear their hair
    // if we swap the endian in the middle of the file, so we stick to big endian in the
    // TPR file for now - and thus we ask the serializer to swap if this host is little endian.
    gmx::InMemoryDeserializer tprBodyDeserializer(partialDeserializedTpr->body,
                                                  partialDeserializedTpr->header.isDouble,
                                                  gmx::EndianSwapBehavior::SwapIfHostIsLittleEndian);
    return do_tpx_body(&tprBodyDeserializer, &partialDeserializedTpr->header, ir, state, x, v, mtop);
}

PbcType completeTprDeserialization(PartialDeserializedTprFile* partialDeserializedTpr,
                                   t_inputrec*                 ir,
                                   gmx_mtop_t*                 mtop)
{
    return completeTprDeserialization(partialDeserializedTpr, ir, nullptr, nullptr, nullptr, mtop);
}

PartialDeserializedTprFile
read_tpx_state(const std::filesystem::path& fn, t_inputrec* ir, t_state* state, gmx_mtop_t* mtop)
{
    t_fileio* fio;
    fio = open_tpx(fn, "r");
    gmx::FileIOXdrSerializer   serializer(fio);
    PartialDeserializedTprFile partialDeserializedTpr;
    do_tpxheader(&serializer, &partialDeserializedTpr.header, fn, fio, ir == nullptr);
    partialDeserializedTpr =
            readTpxBody(&partialDeserializedTpr.header, &serializer, ir, state, nullptr, nullptr, mtop);
    close_tpx(fio);
    return partialDeserializedTpr;
}

PbcType read_tpx(const std::filesystem::path& fn, t_inputrec* ir, matrix box, int* natoms, rvec* x, rvec* v, gmx_mtop_t* mtop)
{
    t_fileio* fio;
    t_state   state;

    TpxFileHeader tpx;
    fio = open_tpx(fn, "r");
    gmx::FileIOXdrSerializer serializer(fio);
    do_tpxheader(&serializer, &tpx, fn, fio, ir == nullptr);
    PartialDeserializedTprFile partialDeserializedTpr =
            readTpxBody(&tpx, &serializer, ir, &state, x, v, mtop);
    close_tpx(fio);
    if (mtop != nullptr && natoms != nullptr)
    {
        *natoms = mtop->natoms;
    }
    if (box)
    {
        copy_mat(state.box, box);
    }
    return partialDeserializedTpr.pbcType;
}

PbcType read_tpx_top(const std::filesystem::path& fn,
                     t_inputrec*                  ir,
                     matrix                       box,
                     int*                         natoms,
                     rvec*                        x,
                     rvec*                        v,
                     t_topology*                  top)
{
    gmx_mtop_t mtop;
    PbcType    pbcType;

    pbcType = read_tpx(fn, ir, box, natoms, x, v, &mtop);

    *top = gmx_mtop_t_to_t_topology(&mtop, true);

    return pbcType;
}

gmx_bool fn2bTPX(const std::filesystem::path& file)
{
    return (efTPR == fn2ftp(file));
}

void pr_tpxheader(FILE* fp, int indent, const char* title, const TpxFileHeader* sh)
{
    if (available(fp, sh, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "bIr    = %spresent\n", sh->bIr ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bBox   = %spresent\n", sh->bBox ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bTop   = %spresent\n", sh->bTop ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bX     = %spresent\n", sh->bX ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bV     = %spresent\n", sh->bV ? "" : "not ");
        pr_indent(fp, indent);
        fprintf(fp, "bF     = %spresent\n", sh->bF ? "" : "not ");

        pr_indent(fp, indent);
        fprintf(fp, "natoms = %d\n", sh->natoms);
        pr_indent(fp, indent);
        fprintf(fp, "lambda = %e\n", sh->lambda);
        pr_indent(fp, indent);
        fprintf(fp, "buffer size = %" PRId64 "\n", sh->sizeOfTprBody);
    }
}
