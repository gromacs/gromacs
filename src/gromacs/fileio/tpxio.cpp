/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

/* This file is completely threadsafe - keep it that way! */

#include "tpxio.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio-xdr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/txtdump.h"

#define TPX_TAG_RELEASE  "release"

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
static const char *tpx_tag = TPX_TAG_RELEASE;

/*! \brief Enum of values that describe the contents of a tpr file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element just above tpxv_Count
 * in this enumeration, and write code below that does the right thing
 * according to the value of file_version.
 */
enum tpxv {
    tpxv_ComputationalElectrophysiology = 96,                /**< support for ion/water position swaps (computational electrophysiology) */
    tpxv_Use64BitRandomSeed,                                 /**< change ld_seed from int to gmx_int64_t */
    tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, /**< potentials for supporting coarse-grained force fields */
    tpxv_InteractiveMolecularDynamics,                       /**< interactive molecular dynamics (IMD) */
    tpxv_RemoveObsoleteParameters1,                          /**< remove optimize_fft, dihre_fc, nstcheckpoint */
    tpxv_PullCoordTypeGeom,                                  /**< add pull type and geometry per group and flat-bottom */
    tpxv_PullGeomDirRel,                                     /**< add pull geometry direction-relative */
    tpxv_IntermolecularBondeds,                              /**< permit inter-molecular bonded interactions in the topology */
    tpxv_CompElWithSwapLayerOffset,                          /**< added parameters for improved CompEl setups */
    tpxv_CompElPolyatomicIonsAndMultipleIonTypes,            /**< CompEl now can handle polyatomic ions and more than two types of ions */
    tpxv_RemoveAdress,                                       /**< removed support for AdResS */
    tpxv_PullCoordNGroup,                                    /**< add ngroup to pull coord */
    tpxv_RemoveTwinRange,                                    /**< removed support for twin-range interactions */
    tpxv_ReplacePullPrintCOM12,                              /**< Replaced print-com-1, 2 with pull-print-com */
    tpxv_PullExternalPotential,                              /**< Added pull type external potential */
    tpxv_GenericParamsForElectricField,                      /**< Introduced KeyValueTree and moved electric field parameters */
    tpxv_AcceleratedWeightHistogram,                         /**< sampling with accelerated weight histogram method (AWH) */
    tpxv_Count                                               /**< the total number of tpxv versions */
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


/* This number should only be increased when you edit the TOPOLOGY section
 * or the HEADER of the tpx format.
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
static const int tpx_generation = 26;

/* This number should be the most recent backwards incompatible version
 * I.e., if this number is 9, we cannot read tpx version 9 with this code.
 */
static const int tpx_incompatible_version = 30; // GMX3.2 has version 31



/* Struct used to maintain tpx compatibility when function types are added */
typedef struct {
    int fvnr;  /* file version number in which the function type first appeared */
    int ftype; /* function type */
} t_ftupd;

/*
 * The entries should be ordered in:
 * 1. ascending function type number
 * 2. ascending file version number
 */
static const t_ftupd ftupd[] = {
    { 34, F_FENEBONDS         },
    { 43, F_TABBONDS          },
    { 43, F_TABBONDSNC        },
    { 70, F_RESTRBONDS        },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES },
    { 76, F_LINEAR_ANGLES     },
    { 34, F_QUARTIC_ANGLES    },
    { 43, F_TABANGLES         },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS },
    { 43, F_TABDIHS           },
    { 65, F_CMAP              },
    { 60, F_GB12              },
    { 61, F_GB13              },
    { 61, F_GB14              },
    { 72, F_GBPOL             },
    { 72, F_NPSOLVATION       },
    { 41, F_LJC14_Q           },
    { 41, F_LJC_PAIRS_NB      },
    { 32, F_BHAM_LR_NOLONGERUSED },
    { 32, F_RF_EXCL           },
    { 32, F_COUL_RECIP        },
    { 93, F_LJ_RECIP          },
    { 46, F_DPD               },
    { 36, F_THOLE_POL         },
    { 90, F_FBPOSRES          },
    { 49, F_VSITE4FDN         },
    { 50, F_VSITEN            },
    { 46, F_COM_PULL          },
    { 46, F_ECONSERVED        },
    { 69, F_VTEMP_NOLONGERUSED},
    { 66, F_PDISPCORR         },
    { 54, F_DVDL_CONSTR       },
    { 76, F_ANHARM_POL        },
    { 79, F_DVDL_COUL         },
    { 79, F_DVDL_VDW,         },
    { 79, F_DVDL_BONDED,      },
    { 79, F_DVDL_RESTRAINT    },
    { 79, F_DVDL_TEMPERATURE  },
};
#define NFTUPD asize(ftupd)

/* Needed for backward compatibility */
#define MAXNODES 256

/**************************************************************
 *
 * Now the higer level routines that do io of the structures and arrays
 *
 **************************************************************/
static void do_pullgrp_tpx_pre95(t_fileio     *fio,
                                 t_pull_group *pgrp,
                                 t_pull_coord *pcrd,
                                 gmx_bool      bRead,
                                 int           file_version)
{
    rvec tmp;

    gmx_fio_do_int(fio, pgrp->nat);
    if (bRead)
    {
        snew(pgrp->ind, pgrp->nat);
    }
    gmx_fio_ndo_int(fio, pgrp->ind, pgrp->nat);
    gmx_fio_do_int(fio, pgrp->nweight);
    if (bRead)
    {
        snew(pgrp->weight, pgrp->nweight);
    }
    gmx_fio_ndo_real(fio, pgrp->weight, pgrp->nweight);
    gmx_fio_do_int(fio, pgrp->pbcatom);
    gmx_fio_do_rvec(fio, pcrd->vec);
    clear_rvec(pcrd->origin);
    gmx_fio_do_rvec(fio, tmp);
    pcrd->init = tmp[0];
    gmx_fio_do_real(fio, pcrd->rate);
    gmx_fio_do_real(fio, pcrd->k);
    if (file_version >= 56)
    {
        gmx_fio_do_real(fio, pcrd->kB);
    }
    else
    {
        pcrd->kB = pcrd->k;
    }
}

static void do_pull_group(t_fileio *fio, t_pull_group *pgrp, gmx_bool bRead)
{
    gmx_fio_do_int(fio, pgrp->nat);
    if (bRead)
    {
        snew(pgrp->ind, pgrp->nat);
    }
    gmx_fio_ndo_int(fio, pgrp->ind, pgrp->nat);
    gmx_fio_do_int(fio, pgrp->nweight);
    if (bRead)
    {
        snew(pgrp->weight, pgrp->nweight);
    }
    gmx_fio_ndo_real(fio, pgrp->weight, pgrp->nweight);
    gmx_fio_do_int(fio, pgrp->pbcatom);
}

static void do_pull_coord(t_fileio *fio, t_pull_coord *pcrd,
                          gmx_bool bRead, int file_version,
                          int ePullOld, int eGeomOld, ivec dimOld)
{
    if (file_version >= tpxv_PullCoordNGroup)
    {
        gmx_fio_do_int(fio,  pcrd->eType);
        if (file_version >= tpxv_PullExternalPotential)
        {
            if (pcrd->eType == epullEXTERNAL)
            {
                if (bRead)
                {
                    char buf[STRLEN];

                    gmx_fio_do_string(fio, buf);
                    pcrd->externalPotentialProvider = gmx_strdup(buf);
                }
                else
                {
                    gmx_fio_do_string(fio, pcrd->externalPotentialProvider);
                }
            }
            else
            {
                pcrd->externalPotentialProvider = nullptr;
            }
        }
        else
        {
            if (bRead)
            {
                pcrd->externalPotentialProvider = nullptr;
            }
        }
        /* Note that we try to support adding new geometries without
         * changing the tpx version. This requires checks when printing the
         * geometry string and a check and fatal_error in init_pull.
         */
        gmx_fio_do_int(fio,  pcrd->eGeom);
        gmx_fio_do_int(fio,  pcrd->ngroup);
        if (pcrd->ngroup <= c_pullCoordNgroupMax)
        {
            gmx_fio_ndo_int(fio, pcrd->group, pcrd->ngroup);
        }
        else
        {
            /* More groups in file than supported, this must be a new geometry
             * that is not supported by our current code. Since we will not
             * use the groups for this coord (checks in the pull and WHAM code
             * ensure this), we can ignore the groups and set ngroup=0.
             */
            int *dum;
            snew(dum, pcrd->ngroup);
            gmx_fio_ndo_int(fio, dum, pcrd->ngroup);
            sfree(dum);

            pcrd->ngroup = 0;
        }
        gmx_fio_do_ivec(fio, pcrd->dim);
    }
    else
    {
        pcrd->ngroup = 2;
        gmx_fio_do_int(fio, pcrd->group[0]);
        gmx_fio_do_int(fio, pcrd->group[1]);
        if (file_version >= tpxv_PullCoordTypeGeom)
        {
            pcrd->ngroup = (pcrd->eGeom == epullgDIRRELATIVE ? 4 : 2);
            gmx_fio_do_int(fio,  pcrd->eType);
            gmx_fio_do_int(fio,  pcrd->eGeom);
            if (pcrd->ngroup == 4)
            {
                gmx_fio_do_int(fio, pcrd->group[2]);
                gmx_fio_do_int(fio, pcrd->group[3]);
            }
            gmx_fio_do_ivec(fio, pcrd->dim);
        }
        else
        {
            pcrd->eType = ePullOld;
            pcrd->eGeom = eGeomOld;
            copy_ivec(dimOld, pcrd->dim);
        }
    }
    gmx_fio_do_rvec(fio, pcrd->origin);
    gmx_fio_do_rvec(fio, pcrd->vec);
    if (file_version >= tpxv_PullCoordTypeGeom)
    {
        gmx_fio_do_gmx_bool(fio, pcrd->bStart);
    }
    else
    {
        /* This parameter is only printed, but not actually used by mdrun */
        pcrd->bStart = FALSE;
    }
    gmx_fio_do_real(fio, pcrd->init);
    gmx_fio_do_real(fio, pcrd->rate);
    gmx_fio_do_real(fio, pcrd->k);
    gmx_fio_do_real(fio, pcrd->kB);
}

static void do_expandedvals(t_fileio *fio, t_expanded *expand, t_lambda *fepvals, gmx_bool bRead, int file_version)
{
    int      n_lambda = fepvals->n_lambda;

    /* reset the lambda calculation window */
    fepvals->lambda_start_n = 0;
    fepvals->lambda_stop_n  = n_lambda;
    if (file_version >= 79)
    {
        if (n_lambda > 0)
        {
            if (bRead)
            {
                snew(expand->init_lambda_weights, n_lambda);
            }
            gmx_fio_ndo_real(fio, expand->init_lambda_weights, n_lambda);
            gmx_fio_do_gmx_bool(fio, expand->bInit_weights);
        }

        gmx_fio_do_int(fio, expand->nstexpanded);
        gmx_fio_do_int(fio, expand->elmcmove);
        gmx_fio_do_int(fio, expand->elamstats);
        gmx_fio_do_int(fio, expand->lmc_repeats);
        gmx_fio_do_int(fio, expand->gibbsdeltalam);
        gmx_fio_do_int(fio, expand->lmc_forced_nstart);
        gmx_fio_do_int(fio, expand->lmc_seed);
        gmx_fio_do_real(fio, expand->mc_temp);
        gmx_fio_do_int(fio, expand->bSymmetrizedTMatrix);
        gmx_fio_do_int(fio, expand->nstTij);
        gmx_fio_do_int(fio, expand->minvarmin);
        gmx_fio_do_int(fio, expand->c_range);
        gmx_fio_do_real(fio, expand->wl_scale);
        gmx_fio_do_real(fio, expand->wl_ratio);
        gmx_fio_do_real(fio, expand->init_wl_delta);
        gmx_fio_do_gmx_bool(fio, expand->bWLoneovert);
        gmx_fio_do_int(fio, expand->elmceq);
        gmx_fio_do_int(fio, expand->equil_steps);
        gmx_fio_do_int(fio, expand->equil_samples);
        gmx_fio_do_int(fio, expand->equil_n_at_lam);
        gmx_fio_do_real(fio, expand->equil_wl_delta);
        gmx_fio_do_real(fio, expand->equil_ratio);
    }
}

static void do_simtempvals(t_fileio *fio, t_simtemp *simtemp, int n_lambda, gmx_bool bRead,
                           int file_version)
{
    if (file_version >= 79)
    {
        gmx_fio_do_int(fio, simtemp->eSimTempScale);
        gmx_fio_do_real(fio, simtemp->simtemp_high);
        gmx_fio_do_real(fio, simtemp->simtemp_low);
        if (n_lambda > 0)
        {
            if (bRead)
            {
                snew(simtemp->temperatures, n_lambda);
            }
            gmx_fio_ndo_real(fio, simtemp->temperatures, n_lambda);
        }
    }
}

static void do_imd(t_fileio *fio, t_IMD *imd, gmx_bool bRead)
{
    gmx_fio_do_int(fio, imd->nat);
    if (bRead)
    {
        snew(imd->ind, imd->nat);
    }
    gmx_fio_ndo_int(fio, imd->ind, imd->nat);
}

static void do_fepvals(t_fileio *fio, t_lambda *fepvals, gmx_bool bRead, int file_version)
{
    /* i is defined in the ndo_double macro; use g to iterate. */
    int      g;
    real     rdum;

    /* free energy values */

    if (file_version >= 79)
    {
        gmx_fio_do_int(fio, fepvals->init_fep_state);
        gmx_fio_do_double(fio, fepvals->init_lambda);
        gmx_fio_do_double(fio, fepvals->delta_lambda);
    }
    else if (file_version >= 59)
    {
        gmx_fio_do_double(fio, fepvals->init_lambda);
        gmx_fio_do_double(fio, fepvals->delta_lambda);
    }
    else
    {
        gmx_fio_do_real(fio, rdum);
        fepvals->init_lambda = rdum;
        gmx_fio_do_real(fio, rdum);
        fepvals->delta_lambda = rdum;
    }
    if (file_version >= 79)
    {
        gmx_fio_do_int(fio, fepvals->n_lambda);
        if (bRead)
        {
            snew(fepvals->all_lambda, efptNR);
        }
        for (g = 0; g < efptNR; g++)
        {
            if (fepvals->n_lambda > 0)
            {
                if (bRead)
                {
                    snew(fepvals->all_lambda[g], fepvals->n_lambda);
                }
                gmx_fio_ndo_double(fio, fepvals->all_lambda[g], fepvals->n_lambda);
                gmx_fio_ndo_int(fio, fepvals->separate_dvdl, efptNR);
            }
            else if (fepvals->init_lambda >= 0)
            {
                fepvals->separate_dvdl[efptFEP] = TRUE;
            }
        }
    }
    else if (file_version >= 64)
    {
        gmx_fio_do_int(fio, fepvals->n_lambda);
        if (bRead)
        {
            int g;

            snew(fepvals->all_lambda, efptNR);
            /* still allocate the all_lambda array's contents. */
            for (g = 0; g < efptNR; g++)
            {
                if (fepvals->n_lambda > 0)
                {
                    snew(fepvals->all_lambda[g], fepvals->n_lambda);
                }
            }
        }
        gmx_fio_ndo_double(fio, fepvals->all_lambda[efptFEP],
                           fepvals->n_lambda);
        if (fepvals->init_lambda >= 0)
        {
            int g, h;

            fepvals->separate_dvdl[efptFEP] = TRUE;

            if (bRead)
            {
                /* copy the contents of the efptFEP lambda component to all
                   the other components */
                for (g = 0; g < efptNR; g++)
                {
                    for (h = 0; h < fepvals->n_lambda; h++)
                    {
                        if (g != efptFEP)
                        {
                            fepvals->all_lambda[g][h] =
                                fepvals->all_lambda[efptFEP][h];
                        }
                    }
                }
            }
        }
    }
    else
    {
        fepvals->n_lambda     = 0;
        fepvals->all_lambda   = nullptr;
        if (fepvals->init_lambda >= 0)
        {
            fepvals->separate_dvdl[efptFEP] = TRUE;
        }
    }
    gmx_fio_do_real(fio, fepvals->sc_alpha);
    if (file_version >= 38)
    {
        gmx_fio_do_int(fio, fepvals->sc_power);
    }
    else
    {
        fepvals->sc_power = 2;
    }
    if (file_version >= 79)
    {
        gmx_fio_do_real(fio, fepvals->sc_r_power);
    }
    else
    {
        fepvals->sc_r_power = 6.0;
    }
    gmx_fio_do_real(fio, fepvals->sc_sigma);
    if (bRead)
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
        gmx_fio_do_int(fio, fepvals->bScCoul);
    }
    else
    {
        fepvals->bScCoul = TRUE;
    }
    if (file_version >= 64)
    {
        gmx_fio_do_int(fio, fepvals->nstdhdl);
    }
    else
    {
        fepvals->nstdhdl = 1;
    }

    if (file_version >= 73)
    {
        gmx_fio_do_int(fio, fepvals->separate_dhdl_file);
        gmx_fio_do_int(fio, fepvals->dhdl_derivatives);
    }
    else
    {
        fepvals->separate_dhdl_file = esepdhdlfileYES;
        fepvals->dhdl_derivatives   = edhdlderivativesYES;
    }
    if (file_version >= 71)
    {
        gmx_fio_do_int(fio, fepvals->dh_hist_size);
        gmx_fio_do_double(fio, fepvals->dh_hist_spacing);
    }
    else
    {
        fepvals->dh_hist_size    = 0;
        fepvals->dh_hist_spacing = 0.1;
    }
    if (file_version >= 79)
    {
        gmx_fio_do_int(fio, fepvals->edHdLPrintEnergy);
    }
    else
    {
        fepvals->edHdLPrintEnergy = edHdLPrintEnergyNO;
    }

    /* handle lambda_neighbors */
    if ((file_version >= 83 && file_version < 90) || file_version >= 92)
    {
        gmx_fio_do_int(fio, fepvals->lambda_neighbors);
        if ( (fepvals->lambda_neighbors >= 0) && (fepvals->init_fep_state >= 0) &&
             (fepvals->init_lambda < 0) )
        {
            fepvals->lambda_start_n = (fepvals->init_fep_state -
                                       fepvals->lambda_neighbors);
            fepvals->lambda_stop_n = (fepvals->init_fep_state +
                                      fepvals->lambda_neighbors + 1);
            if (fepvals->lambda_start_n < 0)
            {
                fepvals->lambda_start_n = 0;;
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

static void do_awhBias(t_fileio *fio, gmx::AwhBiasParams *awhBiasParams, gmx_bool bRead,
                       int gmx_unused file_version)
{
    gmx_fio_do_int(fio, awhBiasParams->eTarget);
    gmx_fio_do_double(fio, awhBiasParams->targetBetaScaling);
    gmx_fio_do_double(fio, awhBiasParams->targetCutoff);
    gmx_fio_do_int(fio, awhBiasParams->eGrowth);
    gmx_fio_do_int(fio, awhBiasParams->bUserData);
    gmx_fio_do_double(fio, awhBiasParams->errorInitial);
    gmx_fio_do_int(fio, awhBiasParams->ndim);
    gmx_fio_do_int(fio, awhBiasParams->shareGroup);
    gmx_fio_do_gmx_bool(fio, awhBiasParams->equilibrateHistogram);

    if (bRead)
    {
        snew(awhBiasParams->dimParams, awhBiasParams->ndim);
    }

    for (int d = 0; d < awhBiasParams->ndim; d++)
    {
        gmx::AwhDimParams *dimParams = &awhBiasParams->dimParams[d];

        gmx_fio_do_int(fio, dimParams->eCoordProvider);
        gmx_fio_do_int(fio, dimParams->coordIndex);
        gmx_fio_do_double(fio, dimParams->origin);
        gmx_fio_do_double(fio, dimParams->end);
        gmx_fio_do_double(fio, dimParams->period);
        gmx_fio_do_double(fio, dimParams->forceConstant);
        gmx_fio_do_double(fio, dimParams->diffusion);
        gmx_fio_do_double(fio, dimParams->coordValueInit);
        gmx_fio_do_double(fio, dimParams->coverDiameter);
    }
}

static void do_awh(t_fileio *fio, gmx::AwhParams *awhParams, gmx_bool bRead,
                   int gmx_unused file_version)
{
    gmx_fio_do_int(fio, awhParams->numBias);
    gmx_fio_do_int(fio, awhParams->nstOut);
    gmx_fio_do_int64(fio, awhParams->seed);
    gmx_fio_do_int(fio, awhParams->nstSampleCoord);
    gmx_fio_do_int(fio, awhParams->numSamplesUpdateFreeEnergy);
    gmx_fio_do_int(fio, awhParams->ePotential);
    gmx_fio_do_gmx_bool(fio, awhParams->shareBiasMultisim);

    if (awhParams->numBias > 0)
    {
        if (bRead)
        {
            snew(awhParams->awhBiasParams, awhParams->numBias);
        }

        for (int k = 0; k < awhParams->numBias; k++)
        {
            do_awhBias(fio, &awhParams->awhBiasParams[k], bRead, file_version);
        }
    }
}

static void do_pull(t_fileio *fio, pull_params_t *pull, gmx_bool bRead,
                    int file_version, int ePullOld)
{
    int  eGeomOld = -1;
    ivec dimOld;
    int  g;

    if (file_version >= 95)
    {
        gmx_fio_do_int(fio, pull->ngroup);
    }
    gmx_fio_do_int(fio, pull->ncoord);
    if (file_version < 95)
    {
        pull->ngroup = pull->ncoord + 1;
    }
    if (file_version < tpxv_PullCoordTypeGeom)
    {
        real dum;

        gmx_fio_do_int(fio, eGeomOld);
        gmx_fio_do_ivec(fio, dimOld);
        /* The inner cylinder radius, now removed */
        gmx_fio_do_real(fio, dum);
    }
    gmx_fio_do_real(fio, pull->cylinder_r);
    gmx_fio_do_real(fio, pull->constr_tol);
    if (file_version >= 95)
    {
        gmx_fio_do_int(fio, pull->bPrintCOM);
        /* With file_version < 95 this value is set below */
    }
    if (file_version >= tpxv_ReplacePullPrintCOM12)
    {
        gmx_fio_do_int(fio, pull->bPrintRefValue);
        gmx_fio_do_int(fio, pull->bPrintComp);
    }
    else if (file_version >= tpxv_PullCoordTypeGeom)
    {
        int idum;
        gmx_fio_do_int(fio, idum); /* used to be bPrintCOM2 */
        gmx_fio_do_int(fio, pull->bPrintRefValue);
        gmx_fio_do_int(fio, pull->bPrintComp);
    }
    else
    {
        pull->bPrintRefValue = FALSE;
        pull->bPrintComp     = TRUE;
    }
    gmx_fio_do_int(fio, pull->nstxout);
    gmx_fio_do_int(fio, pull->nstfout);
    if (bRead)
    {
        snew(pull->group, pull->ngroup);
        snew(pull->coord, pull->ncoord);
    }
    if (file_version < 95)
    {
        /* epullgPOS for position pulling, before epullgDIRPBC was removed */
        if (eGeomOld == epullgDIRPBC)
        {
            gmx_fatal(FARGS, "pull-geometry=position is no longer supported");
        }
        if (eGeomOld > epullgDIRPBC)
        {
            eGeomOld -= 1;
        }

        for (g = 0; g < pull->ngroup; g++)
        {
            /* We read and ignore a pull coordinate for group 0 */
            do_pullgrp_tpx_pre95(fio, &pull->group[g], &pull->coord[std::max(g-1, 0)],
                                 bRead, file_version);
            if (g > 0)
            {
                pull->coord[g-1].group[0] = 0;
                pull->coord[g-1].group[1] = g;
            }
        }

        pull->bPrintCOM = (pull->group[0].nat > 0);
    }
    else
    {
        for (g = 0; g < pull->ngroup; g++)
        {
            do_pull_group(fio, &pull->group[g], bRead);
        }
        for (g = 0; g < pull->ncoord; g++)
        {
            do_pull_coord(fio, &pull->coord[g],
                          bRead, file_version, ePullOld, eGeomOld, dimOld);
        }
    }
}


static void do_rotgrp(t_fileio *fio, t_rotgrp *rotg, gmx_bool bRead)
{
    gmx_fio_do_int(fio, rotg->eType);
    gmx_fio_do_int(fio, rotg->bMassW);
    gmx_fio_do_int(fio, rotg->nat);
    if (bRead)
    {
        snew(rotg->ind, rotg->nat);
    }
    gmx_fio_ndo_int(fio, rotg->ind, rotg->nat);
    if (bRead)
    {
        snew(rotg->x_ref, rotg->nat);
    }
    gmx_fio_ndo_rvec(fio, rotg->x_ref, rotg->nat);
    gmx_fio_do_rvec(fio, rotg->vec);
    gmx_fio_do_rvec(fio, rotg->pivot);
    gmx_fio_do_real(fio, rotg->rate);
    gmx_fio_do_real(fio, rotg->k);
    gmx_fio_do_real(fio, rotg->slab_dist);
    gmx_fio_do_real(fio, rotg->min_gaussian);
    gmx_fio_do_real(fio, rotg->eps);
    gmx_fio_do_int(fio, rotg->eFittype);
    gmx_fio_do_int(fio, rotg->PotAngle_nstep);
    gmx_fio_do_real(fio, rotg->PotAngle_step);
}

static void do_rot(t_fileio *fio, t_rot *rot, gmx_bool bRead)
{
    int g;

    gmx_fio_do_int(fio, rot->ngrp);
    gmx_fio_do_int(fio, rot->nstrout);
    gmx_fio_do_int(fio, rot->nstsout);
    if (bRead)
    {
        snew(rot->grp, rot->ngrp);
    }
    for (g = 0; g < rot->ngrp; g++)
    {
        do_rotgrp(fio, &rot->grp[g], bRead);
    }
}


static void do_swapgroup(t_fileio *fio, t_swapGroup *g, gmx_bool bRead)
{

    /* Name of the group or molecule */
    if (bRead)
    {
        char buf[STRLEN];

        gmx_fio_do_string(fio, buf);
        g->molname = gmx_strdup(buf);
    }
    else
    {
        gmx_fio_do_string(fio, g->molname);
    }

    /* Number of atoms in the group */
    gmx_fio_do_int(fio, g->nat);

    /* The group's atom indices */
    if (bRead)
    {
        snew(g->ind, g->nat);
    }
    gmx_fio_ndo_int(fio, g->ind, g->nat);

    /* Requested counts for compartments A and B */
    gmx_fio_ndo_int(fio, g->nmolReq, eCompNR);
}

static void do_swapcoords_tpx(t_fileio *fio, t_swapcoords *swap, gmx_bool bRead, int file_version)
{
    /* Enums for better readability of the code */
    enum {
        eCompA = 0, eCompB
    };
    enum {
        eChannel0 = 0, eChannel1
    };


    if (file_version >= tpxv_CompElPolyatomicIonsAndMultipleIonTypes)
    {
        /* The total number of swap groups is the sum of the fixed groups
         * (split0, split1, solvent), and the user-defined groups (2+ types of ions)
         */
        gmx_fio_do_int(fio, swap->ngrp);
        if (bRead)
        {
            snew(swap->grp, swap->ngrp);
        }
        for (int ig = 0; ig < swap->ngrp; ig++)
        {
            do_swapgroup(fio, &swap->grp[ig], bRead);
        }
        gmx_fio_do_int(fio, swap->massw_split[eChannel0]);
        gmx_fio_do_int(fio, swap->massw_split[eChannel1]);
        gmx_fio_do_int(fio, swap->nstswap);
        gmx_fio_do_int(fio, swap->nAverage);
        gmx_fio_do_real(fio, swap->threshold);
        gmx_fio_do_real(fio, swap->cyl0r);
        gmx_fio_do_real(fio, swap->cyl0u);
        gmx_fio_do_real(fio, swap->cyl0l);
        gmx_fio_do_real(fio, swap->cyl1r);
        gmx_fio_do_real(fio, swap->cyl1u);
        gmx_fio_do_real(fio, swap->cyl1l);
    }
    else
    {
        /*** Support reading older CompEl .tpr files ***/

        /* In the original CompEl .tpr files, we always have 5 groups: */
        swap->ngrp = 5;
        snew(swap->grp, swap->ngrp);

        swap->grp[eGrpSplit0 ].molname = gmx_strdup("split0" );  // group 0: split0
        swap->grp[eGrpSplit1 ].molname = gmx_strdup("split1" );  // group 1: split1
        swap->grp[eGrpSolvent].molname = gmx_strdup("solvent");  // group 2: solvent
        swap->grp[3          ].molname = gmx_strdup("anions" );  // group 3: anions
        swap->grp[4          ].molname = gmx_strdup("cations");  // group 4: cations

        gmx_fio_do_int(fio, swap->grp[3].nat);
        gmx_fio_do_int(fio, swap->grp[eGrpSolvent].nat);
        gmx_fio_do_int(fio, swap->grp[eGrpSplit0].nat);
        gmx_fio_do_int(fio, swap->massw_split[eChannel0]);
        gmx_fio_do_int(fio, swap->grp[eGrpSplit1].nat);
        gmx_fio_do_int(fio, swap->massw_split[eChannel1]);
        gmx_fio_do_int(fio, swap->nstswap);
        gmx_fio_do_int(fio, swap->nAverage);
        gmx_fio_do_real(fio, swap->threshold);
        gmx_fio_do_real(fio, swap->cyl0r);
        gmx_fio_do_real(fio, swap->cyl0u);
        gmx_fio_do_real(fio, swap->cyl0l);
        gmx_fio_do_real(fio, swap->cyl1r);
        gmx_fio_do_real(fio, swap->cyl1u);
        gmx_fio_do_real(fio, swap->cyl1l);

        // The order[] array keeps compatibility with older .tpr files
        // by reading in the groups in the classic order
        {
            const int order[4] = {3, eGrpSolvent, eGrpSplit0, eGrpSplit1 };

            for (int ig = 0; ig < 4; ig++)
            {
                int g = order[ig];
                snew(swap->grp[g].ind, swap->grp[g].nat);
                gmx_fio_ndo_int(fio, swap->grp[g].ind, swap->grp[g].nat);
            }
        }

        for (int j = eCompA; j <= eCompB; j++)
        {
            gmx_fio_do_int(fio, swap->grp[3].nmolReq[j]); // group 3 = anions
            gmx_fio_do_int(fio, swap->grp[4].nmolReq[j]); // group 4 = cations
        }
    }                                                     /* End support reading older CompEl .tpr files */

    if (file_version >= tpxv_CompElWithSwapLayerOffset)
    {
        gmx_fio_do_real(fio, swap->bulkOffset[eCompA]);
        gmx_fio_do_real(fio, swap->bulkOffset[eCompB]);
    }

}

static void do_legacy_efield(t_fileio *fio, gmx::KeyValueTreeObjectBuilder *root)
{
    const char *const dimName[] = { "x", "y", "z" };

    auto              appliedForcesObj = root->addObject("applied-forces");
    auto              efieldObj        = appliedForcesObj.addObject("electric-field");
    // The content of the tpr file for this feature has
    // been the same since gromacs 4.0 that was used for
    // developing.
    for (int j = 0; j < DIM; ++j)
    {
        int n, nt;
        gmx_fio_do_int(fio, n);
        gmx_fio_do_int(fio, nt);
        std::vector<real> aa(n+1), phi(nt+1), at(nt+1), phit(nt+1);
        gmx_fio_ndo_real(fio, aa.data(),  n);
        gmx_fio_ndo_real(fio, phi.data(), n);
        gmx_fio_ndo_real(fio, at.data(),  nt);
        gmx_fio_ndo_real(fio, phit.data(), nt);
        if (n > 0)
        {
            if (n > 1 || nt > 1)
            {
                gmx_fatal(FARGS, "Can not handle tpr files with more than one electric field term per direction.");
            }
            auto dimObj = efieldObj.addObject(dimName[j]);
            dimObj.addValue<real>("E0", aa[0]);
            dimObj.addValue<real>("omega", at[0]);
            dimObj.addValue<real>("t0", phi[0]);
            dimObj.addValue<real>("sigma", phit[0]);
        }
    }
}


static void do_inputrec(t_fileio *fio, t_inputrec *ir, gmx_bool bRead,
                        int file_version, real *fudgeQQ)
{
    int      i, j, k, idum = 0;
    real     rdum, bd_temp;
    gmx_bool bdum = 0;

    if (file_version != tpx_version)
    {
        /* Give a warning about features that are not accessible */
        fprintf(stderr, "Note: file tpx version %d, software tpx version %d\n",
                file_version, tpx_version);
    }

    if (file_version == 0)
    {
        return;
    }

    gmx::KeyValueTreeBuilder       paramsBuilder;
    gmx::KeyValueTreeObjectBuilder paramsObj = paramsBuilder.rootObject();

    /* Basic inputrec stuff */
    gmx_fio_do_int(fio, ir->eI);
    if (file_version >= 62)
    {
        gmx_fio_do_int64(fio, ir->nsteps);
    }
    else
    {
        gmx_fio_do_int(fio, idum);
        ir->nsteps = idum;
    }

    if (file_version >= 62)
    {
        gmx_fio_do_int64(fio, ir->init_step);
    }
    else
    {
        gmx_fio_do_int(fio, idum);
        ir->init_step = idum;
    }

    if (file_version >= 58)
    {
        gmx_fio_do_int(fio, ir->simulation_part);
    }
    else
    {
        ir->simulation_part = 1;
    }

    if (file_version >= 67)
    {
        gmx_fio_do_int(fio, ir->nstcalcenergy);
    }
    else
    {
        ir->nstcalcenergy = 1;
    }
    if (file_version < 53)
    {
        /* The pbc info has been moved out of do_inputrec,
         * since we always want it, also without reading the inputrec.
         */
        gmx_fio_do_int(fio, ir->ePBC);
        if (file_version >= 45)
        {
            gmx_fio_do_int(fio, ir->bPeriodicMols);
        }
        else
        {
            if (ir->ePBC == 2)
            {
                ir->ePBC          = epbcXYZ;
                ir->bPeriodicMols = TRUE;
            }
            else
            {
                ir->bPeriodicMols = FALSE;
            }
        }
    }
    if (file_version >= 81)
    {
        gmx_fio_do_int(fio, ir->cutoff_scheme);
        if (file_version < 94)
        {
            ir->cutoff_scheme = 1 - ir->cutoff_scheme;
        }
    }
    else
    {
        ir->cutoff_scheme = ecutsGROUP;
    }
    gmx_fio_do_int(fio, ir->ns_type);
    gmx_fio_do_int(fio, ir->nstlist);
    gmx_fio_do_int(fio, idum); /* used to be ndelta; not used anymore */
    if (file_version < 41)
    {
        gmx_fio_do_int(fio, idum);
        gmx_fio_do_int(fio, idum);
    }
    if (file_version >= 45)
    {
        gmx_fio_do_real(fio, ir->rtpi);
    }
    else
    {
        ir->rtpi = 0.05;
    }
    gmx_fio_do_int(fio, ir->nstcomm);
    if (file_version > 34)
    {
        gmx_fio_do_int(fio, ir->comm_mode);
    }
    else if (ir->nstcomm < 0)
    {
        ir->comm_mode = ecmANGULAR;
    }
    else
    {
        ir->comm_mode = ecmLINEAR;
    }
    ir->nstcomm = abs(ir->nstcomm);

    /* ignore nstcheckpoint */
    if (file_version < tpxv_RemoveObsoleteParameters1)
    {
        gmx_fio_do_int(fio, idum);
    }

    gmx_fio_do_int(fio, ir->nstcgsteep);

    gmx_fio_do_int(fio, ir->nbfgscorr);

    gmx_fio_do_int(fio, ir->nstlog);
    gmx_fio_do_int(fio, ir->nstxout);
    gmx_fio_do_int(fio, ir->nstvout);
    gmx_fio_do_int(fio, ir->nstfout);
    gmx_fio_do_int(fio, ir->nstenergy);
    gmx_fio_do_int(fio, ir->nstxout_compressed);
    if (file_version >= 59)
    {
        gmx_fio_do_double(fio, ir->init_t);
        gmx_fio_do_double(fio, ir->delta_t);
    }
    else
    {
        gmx_fio_do_real(fio, rdum);
        ir->init_t = rdum;
        gmx_fio_do_real(fio, rdum);
        ir->delta_t = rdum;
    }
    gmx_fio_do_real(fio, ir->x_compression_precision);
    if (file_version >= 81)
    {
        gmx_fio_do_real(fio, ir->verletbuf_tol);
    }
    else
    {
        ir->verletbuf_tol = 0;
    }
    gmx_fio_do_real(fio, ir->rlist);
    if (file_version >= 67 && file_version < tpxv_RemoveTwinRange)
    {
        if (bRead)
        {
            // Reading such a file version could invoke the twin-range
            // scheme, about which mdrun should give a fatal error.
            real dummy_rlistlong = -1;
            gmx_fio_do_real(fio, dummy_rlistlong);

            if (ir->rlist > 0 && (dummy_rlistlong == 0 || dummy_rlistlong > ir->rlist))
            {
                // Get mdrun to issue an error (regardless of
                // ir->cutoff_scheme).
                ir->useTwinRange = true;
            }
            else
            {
                // grompp used to set rlistlong actively. Users were
                // probably also confused and set rlistlong == rlist.
                // However, in all remaining cases, it is safe to let
                // mdrun proceed normally.
                ir->useTwinRange = false;
            }
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
        gmx_fio_do_int(fio, dummy_nstcalclr);
    }
    gmx_fio_do_int(fio, ir->coulombtype);
    if (file_version < 32 && ir->coulombtype == eelRF)
    {
        ir->coulombtype = eelRF_NEC_UNSUPPORTED;
    }
    if (file_version >= 81)
    {
        gmx_fio_do_int(fio, ir->coulomb_modifier);
    }
    else
    {
        ir->coulomb_modifier = (ir->cutoff_scheme == ecutsVERLET ? eintmodPOTSHIFT : eintmodNONE);
    }
    gmx_fio_do_real(fio, ir->rcoulomb_switch);
    gmx_fio_do_real(fio, ir->rcoulomb);
    gmx_fio_do_int(fio, ir->vdwtype);
    if (file_version >= 81)
    {
        gmx_fio_do_int(fio, ir->vdw_modifier);
    }
    else
    {
        ir->vdw_modifier = (ir->cutoff_scheme == ecutsVERLET ? eintmodPOTSHIFT : eintmodNONE);
    }
    gmx_fio_do_real(fio, ir->rvdw_switch);
    gmx_fio_do_real(fio, ir->rvdw);
    gmx_fio_do_int(fio, ir->eDispCorr);
    gmx_fio_do_real(fio, ir->epsilon_r);
    if (file_version >= 37)
    {
        gmx_fio_do_real(fio, ir->epsilon_rf);
    }
    else
    {
        if (EEL_RF(ir->coulombtype))
        {
            ir->epsilon_rf = ir->epsilon_r;
            ir->epsilon_r  = 1.0;
        }
        else
        {
            ir->epsilon_rf = 1.0;
        }
    }
    gmx_fio_do_real(fio, ir->tabext);

    gmx_fio_do_int(fio, ir->gb_algorithm);
    gmx_fio_do_int(fio, ir->nstgbradii);
    gmx_fio_do_real(fio, ir->rgbradii);
    gmx_fio_do_real(fio, ir->gb_saltconc);
    gmx_fio_do_int(fio, ir->implicit_solvent);
    if (file_version >= 55)
    {
        gmx_fio_do_real(fio, ir->gb_epsilon_solvent);
        gmx_fio_do_real(fio, ir->gb_obc_alpha);
        gmx_fio_do_real(fio, ir->gb_obc_beta);
        gmx_fio_do_real(fio, ir->gb_obc_gamma);
        if (file_version >= 60)
        {
            gmx_fio_do_real(fio, ir->gb_dielectric_offset);
            gmx_fio_do_int(fio, ir->sa_algorithm);
        }
        else
        {
            ir->gb_dielectric_offset = 0.009;
            ir->sa_algorithm         = esaAPPROX;
        }
        gmx_fio_do_real(fio, ir->sa_surface_tension);

        /* Override sa_surface_tension if it is not changed in the mpd-file */
        if (ir->sa_surface_tension < 0)
        {
            if (ir->gb_algorithm == egbSTILL)
            {
                ir->sa_surface_tension = 0.0049 * 100 * CAL2JOULE;
            }
            else if (ir->gb_algorithm == egbHCT || ir->gb_algorithm == egbOBC)
            {
                ir->sa_surface_tension = 0.0054 * 100 * CAL2JOULE;
            }
        }

    }
    else
    {
        /* Better use sensible values than insane (0.0) ones... */
        ir->gb_epsilon_solvent = 80;
        ir->gb_obc_alpha       = 1.0;
        ir->gb_obc_beta        = 0.8;
        ir->gb_obc_gamma       = 4.85;
        ir->sa_surface_tension = 2.092;
    }


    if (file_version >= 81)
    {
        gmx_fio_do_real(fio, ir->fourier_spacing);
    }
    else
    {
        ir->fourier_spacing = 0.0;
    }
    gmx_fio_do_int(fio, ir->nkx);
    gmx_fio_do_int(fio, ir->nky);
    gmx_fio_do_int(fio, ir->nkz);
    gmx_fio_do_int(fio, ir->pme_order);
    gmx_fio_do_real(fio, ir->ewald_rtol);

    if (file_version >= 93)
    {
        gmx_fio_do_real(fio, ir->ewald_rtol_lj);
    }
    else
    {
        ir->ewald_rtol_lj = ir->ewald_rtol;
    }
    gmx_fio_do_int(fio, ir->ewald_geometry);
    gmx_fio_do_real(fio, ir->epsilon_surface);

    /* ignore bOptFFT */
    if (file_version < tpxv_RemoveObsoleteParameters1)
    {
        gmx_fio_do_gmx_bool(fio, bdum);
    }

    if (file_version >= 93)
    {
        gmx_fio_do_int(fio, ir->ljpme_combination_rule);
    }
    gmx_fio_do_gmx_bool(fio, ir->bContinuation);
    gmx_fio_do_int(fio, ir->etc);
    /* before version 18, ir->etc was a gmx_bool (ir->btc),
     * but the values 0 and 1 still mean no and
     * berendsen temperature coupling, respectively.
     */
    if (file_version >= 79)
    {
        gmx_fio_do_gmx_bool(fio, ir->bPrintNHChains);
    }
    if (file_version >= 71)
    {
        gmx_fio_do_int(fio, ir->nsttcouple);
    }
    else
    {
        ir->nsttcouple = ir->nstcalcenergy;
    }
    gmx_fio_do_int(fio, ir->epc);
    gmx_fio_do_int(fio, ir->epct);
    if (file_version >= 71)
    {
        gmx_fio_do_int(fio, ir->nstpcouple);
    }
    else
    {
        ir->nstpcouple = ir->nstcalcenergy;
    }
    gmx_fio_do_real(fio, ir->tau_p);
    gmx_fio_do_rvec(fio, ir->ref_p[XX]);
    gmx_fio_do_rvec(fio, ir->ref_p[YY]);
    gmx_fio_do_rvec(fio, ir->ref_p[ZZ]);
    gmx_fio_do_rvec(fio, ir->compress[XX]);
    gmx_fio_do_rvec(fio, ir->compress[YY]);
    gmx_fio_do_rvec(fio, ir->compress[ZZ]);
    if (file_version >= 47)
    {
        gmx_fio_do_int(fio, ir->refcoord_scaling);
        gmx_fio_do_rvec(fio, ir->posres_com);
        gmx_fio_do_rvec(fio, ir->posres_comB);
    }
    else
    {
        ir->refcoord_scaling = erscNO;
        clear_rvec(ir->posres_com);
        clear_rvec(ir->posres_comB);
    }
    if (file_version < 79)
    {
        gmx_fio_do_int(fio, ir->andersen_seed);
    }
    else
    {
        ir->andersen_seed = 0;
    }

    if (file_version < 37)
    {
        gmx_fio_do_real(fio, rdum);
    }

    gmx_fio_do_real(fio, ir->shake_tol);
    if (file_version < 54)
    {
        // cppcheck-suppress redundantPointerOp
        gmx_fio_do_real(fio, *fudgeQQ);
    }

    gmx_fio_do_int(fio, ir->efep);
    do_fepvals(fio, ir->fepvals, bRead, file_version);

    if (file_version >= 79)
    {
        gmx_fio_do_gmx_bool(fio, ir->bSimTemp);
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
        do_simtempvals(fio, ir->simtempvals, ir->fepvals->n_lambda, bRead, file_version);
    }

    if (file_version >= 79)
    {
        gmx_fio_do_gmx_bool(fio, ir->bExpanded);
        if (ir->bExpanded)
        {
            ir->bExpanded = TRUE;
        }
        else
        {
            ir->bExpanded = FALSE;
        }
    }
    if (ir->bExpanded)
    {
        do_expandedvals(fio, ir->expandedvals, ir->fepvals, bRead, file_version);
    }
    if (file_version >= 57)
    {
        gmx_fio_do_int(fio, ir->eDisre);
    }
    gmx_fio_do_int(fio, ir->eDisreWeighting);
    gmx_fio_do_gmx_bool(fio, ir->bDisreMixed);
    gmx_fio_do_real(fio, ir->dr_fc);
    gmx_fio_do_real(fio, ir->dr_tau);
    gmx_fio_do_int(fio, ir->nstdisreout);
    gmx_fio_do_real(fio, ir->orires_fc);
    gmx_fio_do_real(fio, ir->orires_tau);
    gmx_fio_do_int(fio, ir->nstorireout);

    /* ignore dihre_fc */
    if (file_version < 79)
    {
        gmx_fio_do_real(fio, rdum);
        if (file_version < 56)
        {
            gmx_fio_do_real(fio, rdum);
            gmx_fio_do_int(fio, idum);
        }
    }

    gmx_fio_do_real(fio, ir->em_stepsize);
    gmx_fio_do_real(fio, ir->em_tol);
    gmx_fio_do_gmx_bool(fio, ir->bShakeSOR);
    gmx_fio_do_int(fio, ir->niter);
    gmx_fio_do_real(fio, ir->fc_stepsize);
    gmx_fio_do_int(fio, ir->eConstrAlg);
    gmx_fio_do_int(fio, ir->nProjOrder);
    gmx_fio_do_real(fio, ir->LincsWarnAngle);
    gmx_fio_do_int(fio, ir->nLincsIter);
    if (file_version < 33)
    {
        gmx_fio_do_real(fio, bd_temp);
    }
    gmx_fio_do_real(fio, ir->bd_fric);
    if (file_version >= tpxv_Use64BitRandomSeed)
    {
        gmx_fio_do_int64(fio, ir->ld_seed);
    }
    else
    {
        gmx_fio_do_int(fio, idum);
        ir->ld_seed = idum;
    }
    if (file_version >= 33)
    {
        for (i = 0; i < DIM; i++)
        {
            gmx_fio_do_rvec(fio, ir->deform[i]);
        }
    }
    else
    {
        for (i = 0; i < DIM; i++)
        {
            clear_rvec(ir->deform[i]);
        }
    }
    gmx_fio_do_real(fio, ir->cos_accel);

    gmx_fio_do_int(fio, ir->userint1);
    gmx_fio_do_int(fio, ir->userint2);
    gmx_fio_do_int(fio, ir->userint3);
    gmx_fio_do_int(fio, ir->userint4);
    gmx_fio_do_real(fio, ir->userreal1);
    gmx_fio_do_real(fio, ir->userreal2);
    gmx_fio_do_real(fio, ir->userreal3);
    gmx_fio_do_real(fio, ir->userreal4);

    /* AdResS is removed, but we need to be able to read old files,
       and let mdrun refuse to run them */
    if (file_version >= 77 && file_version < tpxv_RemoveAdress)
    {
        gmx_fio_do_gmx_bool(fio, ir->bAdress);
        if (ir->bAdress)
        {
            int  idum, numThermoForceGroups, numEnergyGroups;
            real rdum;
            rvec rvecdum;
            gmx_fio_do_int(fio, idum);
            gmx_fio_do_real(fio, rdum);
            gmx_fio_do_real(fio, rdum);
            gmx_fio_do_real(fio, rdum);
            gmx_fio_do_int(fio, idum);
            gmx_fio_do_int(fio, idum);
            gmx_fio_do_rvec(fio, rvecdum);
            gmx_fio_do_int(fio, numThermoForceGroups);
            gmx_fio_do_real(fio, rdum);
            gmx_fio_do_int(fio, numEnergyGroups);
            gmx_fio_do_int(fio, idum);

            if (numThermoForceGroups > 0)
            {
                std::vector<int> idumn(numThermoForceGroups);
                gmx_fio_ndo_int(fio, idumn.data(), idumn.size());
            }
            if (numEnergyGroups > 0)
            {
                std::vector<int> idumn(numEnergyGroups);
                gmx_fio_ndo_int(fio, idumn.data(), idumn.size());
            }
        }
    }
    else
    {
        ir->bAdress = FALSE;
    }

    /* pull stuff */
    if (file_version >= 48)
    {
        int ePullOld = 0;

        if (file_version >= tpxv_PullCoordTypeGeom)
        {
            gmx_fio_do_gmx_bool(fio, ir->bPull);
        }
        else
        {
            gmx_fio_do_int(fio, ePullOld);
            ir->bPull = (ePullOld > 0);
            /* We removed the first ePull=ePullNo for the enum */
            ePullOld -= 1;
        }
        if (ir->bPull)
        {
            if (bRead)
            {
                snew(ir->pull, 1);
            }
            do_pull(fio, ir->pull, bRead, file_version, ePullOld);
        }
    }
    else
    {
        ir->bPull = FALSE;
    }

    if (file_version >= tpxv_AcceleratedWeightHistogram)
    {
        gmx_fio_do_gmx_bool(fio, ir->bDoAwh);

        if (ir->bDoAwh)
        {
            if (bRead)
            {
                snew(ir->awhParams, 1);
            }
            do_awh(fio, ir->awhParams, bRead, file_version);
        }
    }
    else
    {
        ir->bDoAwh = FALSE;
    }

    /* Enforced rotation */
    if (file_version >= 74)
    {
        gmx_fio_do_int(fio, ir->bRot);
        if (ir->bRot == TRUE)
        {
            if (bRead)
            {
                snew(ir->rot, 1);
            }
            do_rot(fio, ir->rot, bRead);
        }
    }
    else
    {
        ir->bRot = FALSE;
    }

    /* Interactive molecular dynamics */
    if (file_version >= tpxv_InteractiveMolecularDynamics)
    {
        gmx_fio_do_int(fio, ir->bIMD);
        if (TRUE == ir->bIMD)
        {
            if (bRead)
            {
                snew(ir->imd, 1);
            }
            do_imd(fio, ir->imd, bRead);
        }
    }
    else
    {
        /* We don't support IMD sessions for old .tpr files */
        ir->bIMD = FALSE;
    }

    /* grpopts stuff */
    gmx_fio_do_int(fio, ir->opts.ngtc);
    if (file_version >= 69)
    {
        gmx_fio_do_int(fio, ir->opts.nhchainlength);
    }
    else
    {
        ir->opts.nhchainlength = 1;
    }
    gmx_fio_do_int(fio, ir->opts.ngacc);
    gmx_fio_do_int(fio, ir->opts.ngfrz);
    gmx_fio_do_int(fio, ir->opts.ngener);

    if (bRead)
    {
        snew(ir->opts.nrdf,   ir->opts.ngtc);
        snew(ir->opts.ref_t,  ir->opts.ngtc);
        snew(ir->opts.annealing, ir->opts.ngtc);
        snew(ir->opts.anneal_npoints, ir->opts.ngtc);
        snew(ir->opts.anneal_time, ir->opts.ngtc);
        snew(ir->opts.anneal_temp, ir->opts.ngtc);
        snew(ir->opts.tau_t,  ir->opts.ngtc);
        snew(ir->opts.nFreeze, ir->opts.ngfrz);
        snew(ir->opts.acc,    ir->opts.ngacc);
        snew(ir->opts.egp_flags, ir->opts.ngener*ir->opts.ngener);
    }
    if (ir->opts.ngtc > 0)
    {
        gmx_fio_ndo_real(fio, ir->opts.nrdf, ir->opts.ngtc);
        gmx_fio_ndo_real(fio, ir->opts.ref_t, ir->opts.ngtc);
        gmx_fio_ndo_real(fio, ir->opts.tau_t, ir->opts.ngtc);
        if (file_version < 33 && ir->eI == eiBD)
        {
            for (i = 0; i < ir->opts.ngtc; i++)
            {
                ir->opts.tau_t[i] = bd_temp;
            }
        }
    }
    if (ir->opts.ngfrz > 0)
    {
        gmx_fio_ndo_ivec(fio, ir->opts.nFreeze, ir->opts.ngfrz);
    }
    if (ir->opts.ngacc > 0)
    {
        gmx_fio_ndo_rvec(fio, ir->opts.acc, ir->opts.ngacc);
    }
    gmx_fio_ndo_int(fio, ir->opts.egp_flags,
                    ir->opts.ngener*ir->opts.ngener);

    /* First read the lists with annealing and npoints for each group */
    gmx_fio_ndo_int(fio, ir->opts.annealing, ir->opts.ngtc);
    gmx_fio_ndo_int(fio, ir->opts.anneal_npoints, ir->opts.ngtc);
    for (j = 0; j < (ir->opts.ngtc); j++)
    {
        k = ir->opts.anneal_npoints[j];
        if (bRead)
        {
            snew(ir->opts.anneal_time[j], k);
            snew(ir->opts.anneal_temp[j], k);
        }
        gmx_fio_ndo_real(fio, ir->opts.anneal_time[j], k);
        gmx_fio_ndo_real(fio, ir->opts.anneal_temp[j], k);
    }
    /* Walls */
    if (file_version >= 45)
    {
        gmx_fio_do_int(fio, ir->nwall);
        gmx_fio_do_int(fio, ir->wall_type);
        if (file_version >= 50)
        {
            gmx_fio_do_real(fio, ir->wall_r_linpot);
        }
        else
        {
            ir->wall_r_linpot = -1;
        }
        gmx_fio_do_int(fio, ir->wall_atomtype[0]);
        gmx_fio_do_int(fio, ir->wall_atomtype[1]);
        gmx_fio_do_real(fio, ir->wall_density[0]);
        gmx_fio_do_real(fio, ir->wall_density[1]);
        gmx_fio_do_real(fio, ir->wall_ewald_zfac);
    }
    else
    {
        ir->nwall            = 0;
        ir->wall_type        = 0;
        ir->wall_atomtype[0] = -1;
        ir->wall_atomtype[1] = -1;
        ir->wall_density[0]  = 0;
        ir->wall_density[1]  = 0;
        ir->wall_ewald_zfac  = 3;
    }
    /* Cosine stuff for electric fields */
    if (file_version < tpxv_GenericParamsForElectricField)
    {
        do_legacy_efield(fio, &paramsObj);
    }

    /* Swap ions */
    if (file_version >= tpxv_ComputationalElectrophysiology)
    {
        gmx_fio_do_int(fio, ir->eSwapCoords);
        if (ir->eSwapCoords != eswapNO)
        {
            if (bRead)
            {
                snew(ir->swap, 1);
            }
            do_swapcoords_tpx(fio, ir->swap, bRead, file_version);
        }
    }

    /* QMMM stuff */
    if (file_version >= 39)
    {
        gmx_fio_do_gmx_bool(fio, ir->bQMMM);
        gmx_fio_do_int(fio, ir->QMMMscheme);
        gmx_fio_do_real(fio, ir->scalefactor);
        gmx_fio_do_int(fio, ir->opts.ngQM);
        if (bRead)
        {
            snew(ir->opts.QMmethod,    ir->opts.ngQM);
            snew(ir->opts.QMbasis,     ir->opts.ngQM);
            snew(ir->opts.QMcharge,    ir->opts.ngQM);
            snew(ir->opts.QMmult,      ir->opts.ngQM);
            snew(ir->opts.bSH,         ir->opts.ngQM);
            snew(ir->opts.CASorbitals, ir->opts.ngQM);
            snew(ir->opts.CASelectrons, ir->opts.ngQM);
            snew(ir->opts.SAon,        ir->opts.ngQM);
            snew(ir->opts.SAoff,       ir->opts.ngQM);
            snew(ir->opts.SAsteps,     ir->opts.ngQM);
        }
        if (ir->opts.ngQM > 0)
        {
            gmx_fio_ndo_int(fio, ir->opts.QMmethod, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.QMbasis, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.QMcharge, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.QMmult, ir->opts.ngQM);
            gmx_fio_ndo_gmx_bool(fio, ir->opts.bSH, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.CASorbitals, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.CASelectrons, ir->opts.ngQM);
            gmx_fio_ndo_real(fio, ir->opts.SAon, ir->opts.ngQM);
            gmx_fio_ndo_real(fio, ir->opts.SAoff, ir->opts.ngQM);
            gmx_fio_ndo_int(fio, ir->opts.SAsteps, ir->opts.ngQM);
            /* We leave in dummy i/o for removed parameters to avoid
             * changing the tpr format for every QMMM change.
             */
            std::vector<int> dummy(ir->opts.ngQM, 0);
            gmx_fio_ndo_gmx_bool(fio, dummy.data(), ir->opts.ngQM);
            gmx_fio_ndo_gmx_bool(fio, dummy.data(), ir->opts.ngQM);
        }
        /* end of QMMM stuff */
    }

    if (file_version >= tpxv_GenericParamsForElectricField)
    {
        gmx::FileIOXdrSerializer serializer(fio);
        if (bRead)
        {
            paramsObj.mergeObject(
                    gmx::deserializeKeyValueTree(&serializer));
        }
        else
        {
            GMX_RELEASE_ASSERT(ir->params != nullptr,
                               "Parameters should be present when writing inputrec");
            gmx::serializeKeyValueTree(*ir->params, &serializer);
        }
    }
    if (bRead)
    {
        ir->params = new gmx::KeyValueTreeObject(paramsBuilder.build());
    }
}


static void do_harm(t_fileio *fio, t_iparams *iparams)
{
    gmx_fio_do_real(fio, iparams->harmonic.rA);
    gmx_fio_do_real(fio, iparams->harmonic.krA);
    gmx_fio_do_real(fio, iparams->harmonic.rB);
    gmx_fio_do_real(fio, iparams->harmonic.krB);
}

static void do_iparams(t_fileio *fio, t_functype ftype, t_iparams *iparams,
                       gmx_bool bRead, int file_version)
{
    int      idum;
    real     rdum;

    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_IDIHS:
            do_harm(fio, iparams);
            if ((ftype == F_ANGRES || ftype == F_ANGRESZ) && bRead)
            {
                /* Correct incorrect storage of parameters */
                iparams->pdihs.phiB = iparams->pdihs.phiA;
                iparams->pdihs.cpB  = iparams->pdihs.cpA;
            }
            break;
        case F_RESTRANGLES:
            gmx_fio_do_real(fio, iparams->harmonic.rA);
            gmx_fio_do_real(fio, iparams->harmonic.krA);
            break;
        case F_LINEAR_ANGLES:
            gmx_fio_do_real(fio, iparams->linangle.klinA);
            gmx_fio_do_real(fio, iparams->linangle.aA);
            gmx_fio_do_real(fio, iparams->linangle.klinB);
            gmx_fio_do_real(fio, iparams->linangle.aB);
            break;
        case F_FENEBONDS:
            gmx_fio_do_real(fio, iparams->fene.bm);
            gmx_fio_do_real(fio, iparams->fene.kb);
            break;

        case F_RESTRBONDS:
            gmx_fio_do_real(fio, iparams->restraint.lowA);
            gmx_fio_do_real(fio, iparams->restraint.up1A);
            gmx_fio_do_real(fio, iparams->restraint.up2A);
            gmx_fio_do_real(fio, iparams->restraint.kA);
            gmx_fio_do_real(fio, iparams->restraint.lowB);
            gmx_fio_do_real(fio, iparams->restraint.up1B);
            gmx_fio_do_real(fio, iparams->restraint.up2B);
            gmx_fio_do_real(fio, iparams->restraint.kB);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            gmx_fio_do_real(fio, iparams->tab.kA);
            gmx_fio_do_int(fio, iparams->tab.table);
            gmx_fio_do_real(fio, iparams->tab.kB);
            break;
        case F_CROSS_BOND_BONDS:
            gmx_fio_do_real(fio, iparams->cross_bb.r1e);
            gmx_fio_do_real(fio, iparams->cross_bb.r2e);
            gmx_fio_do_real(fio, iparams->cross_bb.krr);
            break;
        case F_CROSS_BOND_ANGLES:
            gmx_fio_do_real(fio, iparams->cross_ba.r1e);
            gmx_fio_do_real(fio, iparams->cross_ba.r2e);
            gmx_fio_do_real(fio, iparams->cross_ba.r3e);
            gmx_fio_do_real(fio, iparams->cross_ba.krt);
            break;
        case F_UREY_BRADLEY:
            gmx_fio_do_real(fio, iparams->u_b.thetaA);
            gmx_fio_do_real(fio, iparams->u_b.kthetaA);
            gmx_fio_do_real(fio, iparams->u_b.r13A);
            gmx_fio_do_real(fio, iparams->u_b.kUBA);
            if (file_version >= 79)
            {
                gmx_fio_do_real(fio, iparams->u_b.thetaB);
                gmx_fio_do_real(fio, iparams->u_b.kthetaB);
                gmx_fio_do_real(fio, iparams->u_b.r13B);
                gmx_fio_do_real(fio, iparams->u_b.kUBB);
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
            gmx_fio_do_real(fio, iparams->qangle.theta);
            gmx_fio_ndo_real(fio, iparams->qangle.c, 5);
            break;
        case F_BHAM:
            gmx_fio_do_real(fio, iparams->bham.a);
            gmx_fio_do_real(fio, iparams->bham.b);
            gmx_fio_do_real(fio, iparams->bham.c);
            break;
        case F_MORSE:
            gmx_fio_do_real(fio, iparams->morse.b0A);
            gmx_fio_do_real(fio, iparams->morse.cbA);
            gmx_fio_do_real(fio, iparams->morse.betaA);
            if (file_version >= 79)
            {
                gmx_fio_do_real(fio, iparams->morse.b0B);
                gmx_fio_do_real(fio, iparams->morse.cbB);
                gmx_fio_do_real(fio, iparams->morse.betaB);
            }
            else
            {
                iparams->morse.b0B   = iparams->morse.b0A;
                iparams->morse.cbB   = iparams->morse.cbA;
                iparams->morse.betaB = iparams->morse.betaA;
            }
            break;
        case F_CUBICBONDS:
            gmx_fio_do_real(fio, iparams->cubic.b0);
            gmx_fio_do_real(fio, iparams->cubic.kb);
            gmx_fio_do_real(fio, iparams->cubic.kcub);
            break;
        case F_CONNBONDS:
            break;
        case F_POLARIZATION:
            gmx_fio_do_real(fio, iparams->polarize.alpha);
            break;
        case F_ANHARM_POL:
            gmx_fio_do_real(fio, iparams->anharm_polarize.alpha);
            gmx_fio_do_real(fio, iparams->anharm_polarize.drcut);
            gmx_fio_do_real(fio, iparams->anharm_polarize.khyp);
            break;
        case F_WATER_POL:
            gmx_fio_do_real(fio, iparams->wpol.al_x);
            gmx_fio_do_real(fio, iparams->wpol.al_y);
            gmx_fio_do_real(fio, iparams->wpol.al_z);
            gmx_fio_do_real(fio, iparams->wpol.rOH);
            gmx_fio_do_real(fio, iparams->wpol.rHH);
            gmx_fio_do_real(fio, iparams->wpol.rOD);
            break;
        case F_THOLE_POL:
            gmx_fio_do_real(fio, iparams->thole.a);
            gmx_fio_do_real(fio, iparams->thole.alpha1);
            gmx_fio_do_real(fio, iparams->thole.alpha2);
            gmx_fio_do_real(fio, iparams->thole.rfac);
            break;
        case F_LJ:
            gmx_fio_do_real(fio, iparams->lj.c6);
            gmx_fio_do_real(fio, iparams->lj.c12);
            break;
        case F_LJ14:
            gmx_fio_do_real(fio, iparams->lj14.c6A);
            gmx_fio_do_real(fio, iparams->lj14.c12A);
            gmx_fio_do_real(fio, iparams->lj14.c6B);
            gmx_fio_do_real(fio, iparams->lj14.c12B);
            break;
        case F_LJC14_Q:
            gmx_fio_do_real(fio, iparams->ljc14.fqq);
            gmx_fio_do_real(fio, iparams->ljc14.qi);
            gmx_fio_do_real(fio, iparams->ljc14.qj);
            gmx_fio_do_real(fio, iparams->ljc14.c6);
            gmx_fio_do_real(fio, iparams->ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            gmx_fio_do_real(fio, iparams->ljcnb.qi);
            gmx_fio_do_real(fio, iparams->ljcnb.qj);
            gmx_fio_do_real(fio, iparams->ljcnb.c6);
            gmx_fio_do_real(fio, iparams->ljcnb.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            gmx_fio_do_real(fio, iparams->pdihs.phiA);
            gmx_fio_do_real(fio, iparams->pdihs.cpA);
            if ((ftype == F_ANGRES || ftype == F_ANGRESZ) && file_version < 42)
            {
                /* Read the incorrectly stored multiplicity */
                gmx_fio_do_real(fio, iparams->harmonic.rB);
                gmx_fio_do_real(fio, iparams->harmonic.krB);
                iparams->pdihs.phiB = iparams->pdihs.phiA;
                iparams->pdihs.cpB  = iparams->pdihs.cpA;
            }
            else
            {
                gmx_fio_do_real(fio, iparams->pdihs.phiB);
                gmx_fio_do_real(fio, iparams->pdihs.cpB);
                gmx_fio_do_int(fio, iparams->pdihs.mult);
            }
            break;
        case F_RESTRDIHS:
            gmx_fio_do_real(fio, iparams->pdihs.phiA);
            gmx_fio_do_real(fio, iparams->pdihs.cpA);
            break;
        case F_DISRES:
            gmx_fio_do_int(fio, iparams->disres.label);
            gmx_fio_do_int(fio, iparams->disres.type);
            gmx_fio_do_real(fio, iparams->disres.low);
            gmx_fio_do_real(fio, iparams->disres.up1);
            gmx_fio_do_real(fio, iparams->disres.up2);
            gmx_fio_do_real(fio, iparams->disres.kfac);
            break;
        case F_ORIRES:
            gmx_fio_do_int(fio, iparams->orires.ex);
            gmx_fio_do_int(fio, iparams->orires.label);
            gmx_fio_do_int(fio, iparams->orires.power);
            gmx_fio_do_real(fio, iparams->orires.c);
            gmx_fio_do_real(fio, iparams->orires.obs);
            gmx_fio_do_real(fio, iparams->orires.kfac);
            break;
        case F_DIHRES:
            if (file_version < 82)
            {
                gmx_fio_do_int(fio, idum);
                gmx_fio_do_int(fio, idum);
            }
            gmx_fio_do_real(fio, iparams->dihres.phiA);
            gmx_fio_do_real(fio, iparams->dihres.dphiA);
            gmx_fio_do_real(fio, iparams->dihres.kfacA);
            if (file_version >= 82)
            {
                gmx_fio_do_real(fio, iparams->dihres.phiB);
                gmx_fio_do_real(fio, iparams->dihres.dphiB);
                gmx_fio_do_real(fio, iparams->dihres.kfacB);
            }
            else
            {
                iparams->dihres.phiB  = iparams->dihres.phiA;
                iparams->dihres.dphiB = iparams->dihres.dphiA;
                iparams->dihres.kfacB = iparams->dihres.kfacA;
            }
            break;
        case F_POSRES:
            gmx_fio_do_rvec(fio, iparams->posres.pos0A);
            gmx_fio_do_rvec(fio, iparams->posres.fcA);
            gmx_fio_do_rvec(fio, iparams->posres.pos0B);
            gmx_fio_do_rvec(fio, iparams->posres.fcB);
            break;
        case F_FBPOSRES:
            gmx_fio_do_int(fio, iparams->fbposres.geom);
            gmx_fio_do_rvec(fio, iparams->fbposres.pos0);
            gmx_fio_do_real(fio, iparams->fbposres.r);
            gmx_fio_do_real(fio, iparams->fbposres.k);
            break;
        case F_CBTDIHS:
            gmx_fio_ndo_real(fio, iparams->cbtdihs.cbtcA, NR_CBTDIHS);
            break;
        case F_RBDIHS:
            gmx_fio_ndo_real(fio, iparams->rbdihs.rbcA, NR_RBDIHS);
            gmx_fio_ndo_real(fio, iparams->rbdihs.rbcB, NR_RBDIHS);
            break;
        case F_FOURDIHS:
            /* Fourier dihedrals are internally represented
             * as Ryckaert-Bellemans since those are faster to compute.
             */
            gmx_fio_ndo_real(fio, iparams->rbdihs.rbcA, NR_RBDIHS);
            gmx_fio_ndo_real(fio, iparams->rbdihs.rbcB, NR_RBDIHS);
            break;
        case F_CONSTR:
        case F_CONSTRNC:
            gmx_fio_do_real(fio, iparams->constr.dA);
            gmx_fio_do_real(fio, iparams->constr.dB);
            break;
        case F_SETTLE:
            gmx_fio_do_real(fio, iparams->settle.doh);
            gmx_fio_do_real(fio, iparams->settle.dhh);
            break;
        case F_VSITE2:
            gmx_fio_do_real(fio, iparams->vsite.a);
            break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            gmx_fio_do_real(fio, iparams->vsite.a);
            gmx_fio_do_real(fio, iparams->vsite.b);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            gmx_fio_do_real(fio, iparams->vsite.a);
            gmx_fio_do_real(fio, iparams->vsite.b);
            gmx_fio_do_real(fio, iparams->vsite.c);
            break;
        case F_VSITEN:
            gmx_fio_do_int(fio, iparams->vsiten.n);
            gmx_fio_do_real(fio, iparams->vsiten.a);
            break;
        case F_GB12:
        case F_GB13:
        case F_GB14:
            /* We got rid of some parameters in version 68 */
            if (bRead && file_version < 68)
            {
                gmx_fio_do_real(fio, rdum);
                gmx_fio_do_real(fio, rdum);
                gmx_fio_do_real(fio, rdum);
                gmx_fio_do_real(fio, rdum);
            }
            gmx_fio_do_real(fio, iparams->gb.sar);
            gmx_fio_do_real(fio, iparams->gb.st);
            gmx_fio_do_real(fio, iparams->gb.pi);
            gmx_fio_do_real(fio, iparams->gb.gbr);
            gmx_fio_do_real(fio, iparams->gb.bmlt);
            break;
        case F_CMAP:
            gmx_fio_do_int(fio, iparams->cmap.cmapA);
            gmx_fio_do_int(fio, iparams->cmap.cmapB);
            break;
        default:
            gmx_fatal(FARGS, "unknown function type %d (%s) in %s line %d",
                      ftype, interaction_function[ftype].name, __FILE__, __LINE__);
    }
}

static void do_ilist(t_fileio *fio, t_ilist *ilist, gmx_bool bRead, int file_version)
{
    int      i, idum;

    if (file_version < 44)
    {
        for (i = 0; i < MAXNODES; i++)
        {
            gmx_fio_do_int(fio, idum);
        }
    }
    gmx_fio_do_int(fio, ilist->nr);
    if (bRead)
    {
        snew(ilist->iatoms, ilist->nr);
    }
    gmx_fio_ndo_int(fio, ilist->iatoms, ilist->nr);
}

static void do_ffparams(t_fileio *fio, gmx_ffparams_t *ffparams,
                        gmx_bool bRead, int file_version)
{
    int          idum, i;
    unsigned int k;

    gmx_fio_do_int(fio, ffparams->atnr);
    if (file_version < 57)
    {
        gmx_fio_do_int(fio, idum);
    }
    gmx_fio_do_int(fio, ffparams->ntypes);
    if (bRead)
    {
        snew(ffparams->functype, ffparams->ntypes);
        snew(ffparams->iparams, ffparams->ntypes);
    }
    /* Read/write all the function types */
    gmx_fio_ndo_int(fio, ffparams->functype, ffparams->ntypes);

    if (file_version >= 66)
    {
        gmx_fio_do_double(fio, ffparams->reppow);
    }
    else
    {
        ffparams->reppow = 12.0;
    }

    if (file_version >= 57)
    {
        gmx_fio_do_real(fio, ffparams->fudgeQQ);
    }

    /* Check whether all these function types are supported by the code.
     * In practice the code is backwards compatible, which means that the
     * numbering may have to be altered from old numbering to new numbering
     */
    for (i = 0; (i < ffparams->ntypes); i++)
    {
        if (bRead)
        {
            /* Loop over file versions */
            for (k = 0; (k < NFTUPD); k++)
            {
                /* Compare the read file_version to the update table */
                if ((file_version < ftupd[k].fvnr) &&
                    (ffparams->functype[i] >= ftupd[k].ftype))
                {
                    ffparams->functype[i] += 1;
                }
            }
        }

        do_iparams(fio, ffparams->functype[i], &ffparams->iparams[i], bRead,
                   file_version);
    }
}

static void add_settle_atoms(t_ilist *ilist)
{
    int i;

    /* Settle used to only store the first atom: add the other two */
    srenew(ilist->iatoms, 2*ilist->nr);
    for (i = ilist->nr/2-1; i >= 0; i--)
    {
        ilist->iatoms[4*i+0] = ilist->iatoms[2*i+0];
        ilist->iatoms[4*i+1] = ilist->iatoms[2*i+1];
        ilist->iatoms[4*i+2] = ilist->iatoms[2*i+1] + 1;
        ilist->iatoms[4*i+3] = ilist->iatoms[2*i+1] + 2;
    }
    ilist->nr = 2*ilist->nr;
}

static void do_ilists(t_fileio *fio, t_ilist *ilist, gmx_bool bRead,
                      int file_version)
{
    int          j;
    gmx_bool     bClear;
    unsigned int k;

    for (j = 0; (j < F_NRE); j++)
    {
        bClear = FALSE;
        if (bRead)
        {
            for (k = 0; k < NFTUPD; k++)
            {
                if ((file_version < ftupd[k].fvnr) && (j == ftupd[k].ftype))
                {
                    bClear = TRUE;
                }
            }
        }
        if (bClear)
        {
            ilist[j].nr     = 0;
            ilist[j].iatoms = nullptr;
        }
        else
        {
            do_ilist(fio, &ilist[j], bRead, file_version);
            if (file_version < 78 && j == F_SETTLE && ilist[j].nr > 0)
            {
                add_settle_atoms(&ilist[j]);
            }
        }
    }
}

static void do_idef(t_fileio *fio, gmx_ffparams_t *ffparams, gmx_moltype_t *molt,
                    gmx_bool bRead, int file_version)
{
    do_ffparams(fio, ffparams, bRead, file_version);

    if (file_version >= 54)
    {
        gmx_fio_do_real(fio, ffparams->fudgeQQ);
    }

    do_ilists(fio, molt->ilist, bRead, file_version);
}

static void do_block(t_fileio *fio, t_block *block, gmx_bool bRead, int file_version)
{
    int      i, idum, dum_nra, *dum_a;

    if (file_version < 44)
    {
        for (i = 0; i < MAXNODES; i++)
        {
            gmx_fio_do_int(fio, idum);
        }
    }
    gmx_fio_do_int(fio, block->nr);
    if (file_version < 51)
    {
        gmx_fio_do_int(fio, dum_nra);
    }
    if (bRead)
    {
        if ((block->nalloc_index > 0) && (nullptr != block->index))
        {
            sfree(block->index);
        }
        block->nalloc_index = block->nr+1;
        snew(block->index, block->nalloc_index);
    }
    gmx_fio_ndo_int(fio, block->index, block->nr+1);

    if (file_version < 51 && dum_nra > 0)
    {
        snew(dum_a, dum_nra);
        gmx_fio_ndo_int(fio, dum_a, dum_nra);
        sfree(dum_a);
    }
}

static void do_blocka(t_fileio *fio, t_blocka *block, gmx_bool bRead,
                      int file_version)
{
    int      i, idum;

    if (file_version < 44)
    {
        for (i = 0; i < MAXNODES; i++)
        {
            gmx_fio_do_int(fio, idum);
        }
    }
    gmx_fio_do_int(fio, block->nr);
    gmx_fio_do_int(fio, block->nra);
    if (bRead)
    {
        block->nalloc_index = block->nr+1;
        snew(block->index, block->nalloc_index);
        block->nalloc_a = block->nra;
        snew(block->a, block->nalloc_a);
    }
    gmx_fio_ndo_int(fio, block->index, block->nr+1);
    gmx_fio_ndo_int(fio, block->a, block->nra);
}

/* This is a primitive routine to make it possible to translate atomic numbers
 * to element names when reading TPR files, without making the Gromacs library
 * directory a dependency on mdrun (which is the case if we need elements.dat).
 */
static const char *
atomicnumber_to_element(int atomicnumber)
{
    const char * p;

    /* This does not have to be complete, so we only include elements likely
     * to occur in PDB files.
     */
    switch (atomicnumber)
    {
        case 1:  p = "H";  break;
        case 5:  p = "B";  break;
        case 6:  p = "C";  break;
        case 7:  p = "N";  break;
        case 8:  p = "O";  break;
        case 9:  p = "F";  break;
        case 11: p = "Na"; break;
        case 12: p = "Mg"; break;
        case 15: p = "P";  break;
        case 16: p = "S";  break;
        case 17: p = "Cl"; break;
        case 18: p = "Ar"; break;
        case 19: p = "K";  break;
        case 20: p = "Ca"; break;
        case 25: p = "Mn"; break;
        case 26: p = "Fe"; break;
        case 28: p = "Ni"; break;
        case 29: p = "Cu"; break;
        case 30: p = "Zn"; break;
        case 35: p = "Br"; break;
        case 47: p = "Ag"; break;
        default: p = "";   break;
    }
    return p;
}


static void do_atom(t_fileio *fio, t_atom *atom, int ngrp, gmx_bool bRead,
                    int file_version, gmx_groups_t *groups, int atnr)
{
    int    i, myngrp;

    gmx_fio_do_real(fio, atom->m);
    gmx_fio_do_real(fio, atom->q);
    gmx_fio_do_real(fio, atom->mB);
    gmx_fio_do_real(fio, atom->qB);
    gmx_fio_do_ushort(fio, atom->type);
    gmx_fio_do_ushort(fio, atom->typeB);
    gmx_fio_do_int(fio, atom->ptype);
    gmx_fio_do_int(fio, atom->resind);
    if (file_version >= 52)
    {
        gmx_fio_do_int(fio, atom->atomnumber);
        if (bRead)
        {
            /* Set element string from atomic number if present.
             * This routine returns an empty string if the name is not found.
             */
            std::strncpy(atom->elem, atomicnumber_to_element(atom->atomnumber), 4);
            /* avoid warnings about potentially unterminated string */
            atom->elem[3] = '\0';
        }
    }
    else if (bRead)
    {
        atom->atomnumber = 0;
    }
    if (file_version < 39)
    {
        myngrp = 9;
    }
    else
    {
        myngrp = ngrp;
    }

    if (file_version < 57)
    {
        unsigned char uchar[egcNR];
        gmx_fio_ndo_uchar(fio, uchar, myngrp);
        for (i = myngrp; (i < ngrp); i++)
        {
            uchar[i] = 0;
        }
        /* Copy the old data format to the groups struct */
        for (i = 0; i < ngrp; i++)
        {
            groups->grpnr[i][atnr] = uchar[i];
        }
    }
}

static void do_grps(t_fileio *fio, int ngrp, t_grps grps[], gmx_bool bRead,
                    int file_version)
{
    int      j, myngrp;

    if (file_version < 39)
    {
        myngrp = 9;
    }
    else
    {
        myngrp = ngrp;
    }

    for (j = 0; (j < ngrp); j++)
    {
        if (j < myngrp)
        {
            gmx_fio_do_int(fio, grps[j].nr);
            if (bRead)
            {
                snew(grps[j].nm_ind, grps[j].nr);
            }
            gmx_fio_ndo_int(fio, grps[j].nm_ind, grps[j].nr);
        }
        else
        {
            grps[j].nr = 1;
            snew(grps[j].nm_ind, grps[j].nr);
        }
    }
}

static void do_symstr(t_fileio *fio, char ***nm, gmx_bool bRead, t_symtab *symtab)
{
    int ls;

    if (bRead)
    {
        gmx_fio_do_int(fio, ls);
        *nm = get_symtab_handle(symtab, ls);
    }
    else
    {
        ls = lookup_symtab(symtab, *nm);
        gmx_fio_do_int(fio, ls);
    }
}

static void do_strstr(t_fileio *fio, int nstr, char ***nm, gmx_bool bRead,
                      t_symtab *symtab)
{
    int  j;

    for (j = 0; (j < nstr); j++)
    {
        do_symstr(fio, &(nm[j]), bRead, symtab);
    }
}

static void do_resinfo(t_fileio *fio, int n, t_resinfo *ri, gmx_bool bRead,
                       t_symtab *symtab, int file_version)
{
    int  j;

    for (j = 0; (j < n); j++)
    {
        do_symstr(fio, &(ri[j].name), bRead, symtab);
        if (file_version >= 63)
        {
            gmx_fio_do_int(fio, ri[j].nr);
            gmx_fio_do_uchar(fio, ri[j].ic);
        }
        else
        {
            ri[j].nr = j + 1;
            ri[j].ic = ' ';
        }
    }
}

static void do_atoms(t_fileio *fio, t_atoms *atoms, gmx_bool bRead, t_symtab *symtab,
                     int file_version,
                     gmx_groups_t *groups)
{
    int i;

    gmx_fio_do_int(fio, atoms->nr);
    gmx_fio_do_int(fio, atoms->nres);
    if (file_version < 57)
    {
        gmx_fio_do_int(fio, groups->ngrpname);
        for (i = 0; i < egcNR; i++)
        {
            groups->ngrpnr[i] = atoms->nr;
            snew(groups->grpnr[i], groups->ngrpnr[i]);
        }
    }
    if (bRead)
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
        if (file_version < 57)
        {
            snew(groups->grpname, groups->ngrpname);
        }
        atoms->pdbinfo = nullptr;
    }
    else
    {
        GMX_RELEASE_ASSERT(atoms->haveMass && atoms->haveCharge && atoms->haveType && atoms->haveBState, "Mass, charge, atomtype and B-state parameters should be present in t_atoms when writing a tpr file");
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        do_atom(fio, &atoms->atom[i], egcNR, bRead, file_version, groups, i);
    }
    do_strstr(fio, atoms->nr, atoms->atomname, bRead, symtab);
    do_strstr(fio, atoms->nr, atoms->atomtype, bRead, symtab);
    do_strstr(fio, atoms->nr, atoms->atomtypeB, bRead, symtab);

    do_resinfo(fio, atoms->nres, atoms->resinfo, bRead, symtab, file_version);

    if (file_version < 57)
    {
        do_strstr(fio, groups->ngrpname, groups->grpname, bRead, symtab);

        do_grps(fio, egcNR, groups->grps, bRead, file_version);
    }
}

static void do_groups(t_fileio *fio, gmx_groups_t *groups,
                      gmx_bool bRead, t_symtab *symtab,
                      int file_version)
{
    int      g;

    do_grps(fio, egcNR, groups->grps, bRead, file_version);
    gmx_fio_do_int(fio, groups->ngrpname);
    if (bRead)
    {
        snew(groups->grpname, groups->ngrpname);
    }
    do_strstr(fio, groups->ngrpname, groups->grpname, bRead, symtab);
    for (g = 0; g < egcNR; g++)
    {
        gmx_fio_do_int(fio, groups->ngrpnr[g]);
        if (groups->ngrpnr[g] == 0)
        {
            if (bRead)
            {
                groups->grpnr[g] = nullptr;
            }
        }
        else
        {
            if (bRead)
            {
                snew(groups->grpnr[g], groups->ngrpnr[g]);
            }
            gmx_fio_ndo_uchar(fio, groups->grpnr[g], groups->ngrpnr[g]);
        }
    }
}

static void do_atomtypes(t_fileio *fio, t_atomtypes *atomtypes, gmx_bool bRead,
                         int file_version)
{
    int      j;

    gmx_fio_do_int(fio, atomtypes->nr);
    j = atomtypes->nr;
    if (bRead)
    {
        snew(atomtypes->radius, j);
        snew(atomtypes->vol, j);
        snew(atomtypes->surftens, j);
        snew(atomtypes->atomnumber, j);
        snew(atomtypes->gb_radius, j);
        snew(atomtypes->S_hct, j);
    }
    gmx_fio_ndo_real(fio, atomtypes->radius, j);
    gmx_fio_ndo_real(fio, atomtypes->vol, j);
    gmx_fio_ndo_real(fio, atomtypes->surftens, j);
    if (file_version >= 40)
    {
        gmx_fio_ndo_int(fio, atomtypes->atomnumber, j);
    }
    if (file_version >= 60)
    {
        gmx_fio_ndo_real(fio, atomtypes->gb_radius, j);
        gmx_fio_ndo_real(fio, atomtypes->S_hct, j);
    }
}

static void do_symtab(t_fileio *fio, t_symtab *symtab, gmx_bool bRead)
{
    int       i, nr;
    t_symbuf *symbuf;
    char      buf[STRLEN];

    gmx_fio_do_int(fio, symtab->nr);
    nr     = symtab->nr;
    if (bRead)
    {
        snew(symtab->symbuf, 1);
        symbuf          = symtab->symbuf;
        symbuf->bufsize = nr;
        snew(symbuf->buf, nr);
        for (i = 0; (i < nr); i++)
        {
            gmx_fio_do_string(fio, buf);
            symbuf->buf[i] = gmx_strdup(buf);
        }
    }
    else
    {
        symbuf = symtab->symbuf;
        while (symbuf != nullptr)
        {
            for (i = 0; (i < symbuf->bufsize) && (i < nr); i++)
            {
                gmx_fio_do_string(fio, symbuf->buf[i]);
            }
            nr    -= i;
            symbuf = symbuf->next;
        }
        if (nr != 0)
        {
            gmx_fatal(FARGS, "nr of symtab strings left: %d", nr);
        }
    }
}

static void do_cmap(t_fileio *fio, gmx_cmap_t *cmap_grid, gmx_bool bRead)
{
    int i, j, ngrid, gs, nelem;

    gmx_fio_do_int(fio, cmap_grid->ngrid);
    gmx_fio_do_int(fio, cmap_grid->grid_spacing);

    ngrid = cmap_grid->ngrid;
    gs    = cmap_grid->grid_spacing;
    nelem = gs * gs;

    if (bRead)
    {
        snew(cmap_grid->cmapdata, ngrid);

        for (i = 0; i < cmap_grid->ngrid; i++)
        {
            snew(cmap_grid->cmapdata[i].cmap, 4*nelem);
        }
    }

    for (i = 0; i < cmap_grid->ngrid; i++)
    {
        for (j = 0; j < nelem; j++)
        {
            gmx_fio_do_real(fio, cmap_grid->cmapdata[i].cmap[j*4]);
            gmx_fio_do_real(fio, cmap_grid->cmapdata[i].cmap[j*4+1]);
            gmx_fio_do_real(fio, cmap_grid->cmapdata[i].cmap[j*4+2]);
            gmx_fio_do_real(fio, cmap_grid->cmapdata[i].cmap[j*4+3]);
        }
    }
}


static void do_moltype(t_fileio *fio, gmx_moltype_t *molt, gmx_bool bRead,
                       t_symtab *symtab, int file_version,
                       gmx_groups_t *groups)
{
    if (file_version >= 57)
    {
        do_symstr(fio, &(molt->name), bRead, symtab);
    }

    do_atoms(fio, &molt->atoms, bRead, symtab, file_version, groups);

    if (file_version >= 57)
    {
        do_ilists(fio, molt->ilist, bRead, file_version);

        do_block(fio, &molt->cgs, bRead, file_version);
    }

    /* This used to be in the atoms struct */
    do_blocka(fio, &molt->excls, bRead, file_version);
}

static void do_molblock(t_fileio *fio, gmx_molblock_t *molb, gmx_bool bRead)
{
    gmx_fio_do_int(fio, molb->type);
    gmx_fio_do_int(fio, molb->nmol);
    gmx_fio_do_int(fio, molb->natoms_mol);
    /* Position restraint coordinates */
    gmx_fio_do_int(fio, molb->nposres_xA);
    if (molb->nposres_xA > 0)
    {
        if (bRead)
        {
            snew(molb->posres_xA, molb->nposres_xA);
        }
        gmx_fio_ndo_rvec(fio, molb->posres_xA, molb->nposres_xA);
    }
    gmx_fio_do_int(fio, molb->nposres_xB);
    if (molb->nposres_xB > 0)
    {
        if (bRead)
        {
            snew(molb->posres_xB, molb->nposres_xB);
        }
        gmx_fio_ndo_rvec(fio, molb->posres_xB, molb->nposres_xB);
    }

}

static t_block mtop_mols(gmx_mtop_t *mtop)
{
    int     mb, m, a, mol;
    t_block mols;

    mols.nr = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        mols.nr += mtop->molblock[mb].nmol;
    }
    mols.nalloc_index = mols.nr + 1;
    snew(mols.index, mols.nalloc_index);

    a             = 0;
    m             = 0;
    mols.index[m] = a;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        for (mol = 0; mol < mtop->molblock[mb].nmol; mol++)
        {
            a += mtop->molblock[mb].natoms_mol;
            m++;
            mols.index[m] = a;
        }
    }

    return mols;
}

static void add_posres_molblock(gmx_mtop_t *mtop)
{
    t_ilist        *il, *ilfb;
    int             am, i, mol, a;
    gmx_bool        bFE;
    gmx_molblock_t *molb;
    t_iparams      *ip;

    /* posres reference positions are stored in ip->posres (if present) and
       in ip->fbposres (if present). If normal and flat-bottomed posres are present,
       posres.pos0A are identical to fbposres.pos0. */
    il   = &mtop->moltype[0].ilist[F_POSRES];
    ilfb = &mtop->moltype[0].ilist[F_FBPOSRES];
    if (il->nr == 0 && ilfb->nr == 0)
    {
        return;
    }
    am  = 0;
    bFE = FALSE;
    for (i = 0; i < il->nr; i += 2)
    {
        ip = &mtop->ffparams.iparams[il->iatoms[i]];
        am = std::max(am, il->iatoms[i+1]);
        if (ip->posres.pos0B[XX] != ip->posres.pos0A[XX] ||
            ip->posres.pos0B[YY] != ip->posres.pos0A[YY] ||
            ip->posres.pos0B[ZZ] != ip->posres.pos0A[ZZ])
        {
            bFE = TRUE;
        }
    }
    /* This loop is required if we have only flat-bottomed posres:
       - set am
       - bFE == FALSE (no B-state for flat-bottomed posres) */
    if (il->nr == 0)
    {
        for (i = 0; i < ilfb->nr; i += 2)
        {
            am = std::max(am, ilfb->iatoms[i+1]);
        }
    }
    /* Make the posres coordinate block end at a molecule end */
    mol = 0;
    while (am >= mtop->mols.index[mol+1])
    {
        mol++;
    }
    molb             = &mtop->molblock[0];
    molb->nposres_xA = mtop->mols.index[mol+1];
    snew(molb->posres_xA, molb->nposres_xA);
    if (bFE)
    {
        molb->nposres_xB = molb->nposres_xA;
        snew(molb->posres_xB, molb->nposres_xB);
    }
    else
    {
        molb->nposres_xB = 0;
    }
    for (i = 0; i < il->nr; i += 2)
    {
        ip                     = &mtop->ffparams.iparams[il->iatoms[i]];
        a                      = il->iatoms[i+1];
        molb->posres_xA[a][XX] = ip->posres.pos0A[XX];
        molb->posres_xA[a][YY] = ip->posres.pos0A[YY];
        molb->posres_xA[a][ZZ] = ip->posres.pos0A[ZZ];
        if (bFE)
        {
            molb->posres_xB[a][XX] = ip->posres.pos0B[XX];
            molb->posres_xB[a][YY] = ip->posres.pos0B[YY];
            molb->posres_xB[a][ZZ] = ip->posres.pos0B[ZZ];
        }
    }
    if (il->nr == 0)
    {
        /* If only flat-bottomed posres are present, take reference pos from them.
           Here: bFE == FALSE      */
        for (i = 0; i < ilfb->nr; i += 2)
        {
            ip                     = &mtop->ffparams.iparams[ilfb->iatoms[i]];
            a                      = ilfb->iatoms[i+1];
            molb->posres_xA[a][XX] = ip->fbposres.pos0[XX];
            molb->posres_xA[a][YY] = ip->fbposres.pos0[YY];
            molb->posres_xA[a][ZZ] = ip->fbposres.pos0[ZZ];
        }
    }
}

static void set_disres_npair(gmx_mtop_t *mtop)
{
    t_iparams            *ip;
    gmx_mtop_ilistloop_t  iloop;
    t_ilist              *ilist, *il;
    int                   nmol, i, npair;
    t_iatom              *a;

    ip = mtop->ffparams.iparams;

    iloop     = gmx_mtop_ilistloop_init(mtop);
    while (gmx_mtop_ilistloop_next(iloop, &ilist, &nmol))
    {
        il = &ilist[F_DISRES];

        if (il->nr > 0)
        {
            a     = il->iatoms;
            npair = 0;
            for (i = 0; i < il->nr; i += 3)
            {
                npair++;
                if (i+3 == il->nr || ip[a[i]].disres.label != ip[a[i+3]].disres.label)
                {
                    ip[a[i]].disres.npair = npair;
                    npair                 = 0;
                }
            }
        }
    }
}

static void do_mtop(t_fileio *fio, gmx_mtop_t *mtop, gmx_bool bRead,
                    int file_version)
{
    int            mt, mb;
    t_blocka       dumb;

    if (bRead)
    {
        init_mtop(mtop);
    }
    do_symtab(fio, &(mtop->symtab), bRead);

    do_symstr(fio, &(mtop->name), bRead, &(mtop->symtab));

    if (file_version >= 57)
    {
        do_ffparams(fio, &mtop->ffparams, bRead, file_version);

        gmx_fio_do_int(fio, mtop->nmoltype);
    }
    else
    {
        mtop->nmoltype = 1;
    }
    if (bRead)
    {
        snew(mtop->moltype, mtop->nmoltype);
        if (file_version < 57)
        {
            mtop->moltype[0].name = mtop->name;
        }
    }
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        do_moltype(fio, &mtop->moltype[mt], bRead, &mtop->symtab, file_version,
                   &mtop->groups);
    }

    if (file_version >= 57)
    {
        gmx_fio_do_int(fio, mtop->nmolblock);
    }
    else
    {
        mtop->nmolblock = 1;
    }
    if (bRead)
    {
        snew(mtop->molblock, mtop->nmolblock);
    }
    if (file_version >= 57)
    {
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            do_molblock(fio, &mtop->molblock[mb], bRead);
        }
        gmx_fio_do_int(fio, mtop->natoms);
    }
    else
    {
        mtop->molblock[0].type       = 0;
        mtop->molblock[0].nmol       = 1;
        mtop->molblock[0].natoms_mol = mtop->moltype[0].atoms.nr;
        mtop->molblock[0].nposres_xA = 0;
        mtop->molblock[0].nposres_xB = 0;
    }

    if (file_version >= tpxv_IntermolecularBondeds)
    {
        gmx_fio_do_gmx_bool(fio, mtop->bIntermolecularInteractions);
        if (mtop->bIntermolecularInteractions)
        {
            if (bRead)
            {
                snew(mtop->intermolecular_ilist, F_NRE);
            }
            do_ilists(fio, mtop->intermolecular_ilist, bRead, file_version);
        }
    }
    else
    {
        mtop->bIntermolecularInteractions = FALSE;
    }

    do_atomtypes (fio, &(mtop->atomtypes), bRead, file_version);

    if (file_version < 57)
    {
        do_idef (fio, &mtop->ffparams, &mtop->moltype[0], bRead, file_version);
        mtop->natoms = mtop->moltype[0].atoms.nr;
    }

    if (file_version >= 65)
    {
        do_cmap(fio, &mtop->ffparams.cmap_grid, bRead);
    }
    else
    {
        mtop->ffparams.cmap_grid.ngrid        = 0;
        mtop->ffparams.cmap_grid.grid_spacing = 0;
        mtop->ffparams.cmap_grid.cmapdata     = nullptr;
    }

    if (file_version >= 57)
    {
        do_groups(fio, &mtop->groups, bRead, &(mtop->symtab), file_version);
    }

    if (file_version < 57)
    {
        do_block(fio, &mtop->moltype[0].cgs, bRead, file_version);
        do_block(fio, &mtop->mols, bRead, file_version);
        /* Add the posres coordinates to the molblock */
        add_posres_molblock(mtop);
    }
    if (bRead)
    {
        if (file_version >= 57)
        {
            done_block(&mtop->mols);
            mtop->mols = mtop_mols(mtop);
        }
    }

    if (file_version < 51)
    {
        /* Here used to be the shake blocks */
        do_blocka(fio, &dumb, bRead, file_version);
        if (dumb.nr > 0)
        {
            sfree(dumb.index);
        }
        if (dumb.nra > 0)
        {
            sfree(dumb.a);
        }
    }

    if (bRead)
    {
        close_symtab(&(mtop->symtab));
    }
}

/* If TopOnlyOK is TRUE then we can read even future versions
 * of tpx files, provided the fileGeneration hasn't changed.
 * If it is FALSE, we need the inputrecord too, and bail out
 * if the file is newer than the program.
 *
 * The version and generation of the topology (see top of this file)
 * are returned in the two last arguments, if those arguments are non-NULL.
 *
 * If possible, we will read the inputrec even when TopOnlyOK is TRUE.
 */
static void do_tpxheader(t_fileio *fio, gmx_bool bRead, t_tpxheader *tpx,
                         gmx_bool TopOnlyOK, int *fileVersionPointer, int *fileGenerationPointer)
{
    char      buf[STRLEN];
    char      file_tag[STRLEN];
    gmx_bool  bDouble;
    int       precision;
    int       idum = 0;
    real      rdum = 0;
    int       fileVersion;    /* Version number of the code that wrote the file */
    int       fileGeneration; /* Generation version number of the code that wrote the file */

    /* XDR binary topology file */
    precision = sizeof(real);
    if (bRead)
    {
        gmx_fio_do_string(fio, buf);
        if (std::strncmp(buf, "VERSION", 7))
        {
            gmx_fatal(FARGS, "Can not read file %s,\n"
                      "             this file is from a GROMACS version which is older than 2.0\n"
                      "             Make a new one with grompp or use a gro or pdb file, if possible",
                      gmx_fio_getname(fio));
        }
        gmx_fio_do_int(fio, precision);
        bDouble = (precision == sizeof(double));
        if ((precision != sizeof(float)) && !bDouble)
        {
            gmx_fatal(FARGS, "Unknown precision in file %s: real is %d bytes "
                      "instead of %d or %d",
                      gmx_fio_getname(fio), precision, sizeof(float), sizeof(double));
        }
        gmx_fio_setprecision(fio, bDouble);
        fprintf(stderr, "Reading file %s, %s (%s precision)\n",
                gmx_fio_getname(fio), buf, bDouble ? "double" : "single");
    }
    else
    {
        snprintf(buf, STRLEN, "VERSION %s", gmx_version());
        gmx_fio_write_string(fio, buf);
        bDouble = (precision == sizeof(double));
        gmx_fio_setprecision(fio, bDouble);
        gmx_fio_do_int(fio, precision);
        fileVersion = tpx_version;
        sprintf(file_tag, "%s", tpx_tag);
        fileGeneration = tpx_generation;
    }

    /* Check versions! */
    gmx_fio_do_int(fio, fileVersion);

    /* This is for backward compatibility with development versions 77-79
     * where the tag was, mistakenly, placed before the generation,
     * which would cause a segv instead of a proper error message
     * when reading the topology only from tpx with <77 code.
     */
    if (fileVersion >= 77 && fileVersion <= 79)
    {
        gmx_fio_do_string(fio, file_tag);
    }

    gmx_fio_do_int(fio, fileGeneration);

    if (fileVersion >= 81)
    {
        gmx_fio_do_string(fio, file_tag);
    }
    if (bRead)
    {
        if (fileVersion < 77)
        {
            /* Versions before 77 don't have the tag, set it to release */
            sprintf(file_tag, "%s", TPX_TAG_RELEASE);
        }

        if (std::strcmp(file_tag, tpx_tag) != 0)
        {
            fprintf(stderr, "Note: file tpx tag '%s', software tpx tag '%s'\n",
                    file_tag, tpx_tag);

            /* We only support reading tpx files with the same tag as the code
             * or tpx files with the release tag and with lower version number.
             */
            if (std::strcmp(file_tag, TPX_TAG_RELEASE) != 0 && fileVersion < tpx_version)
            {
                gmx_fatal(FARGS, "tpx tag/version mismatch: reading tpx file (%s) version %d, tag '%s' with program for tpx version %d, tag '%s'",
                          gmx_fio_getname(fio), fileVersion, file_tag,
                          tpx_version, tpx_tag);
            }
        }
    }

    if (fileVersionPointer)
    {
        *fileVersionPointer = fileVersion;
    }
    if (fileGenerationPointer)
    {
        *fileGenerationPointer = fileGeneration;
    }

    if ((fileVersion <= tpx_incompatible_version) ||
        ((fileVersion > tpx_version) && !TopOnlyOK) ||
        (fileGeneration > tpx_generation) ||
        tpx_version == 80) /*80 was used by both 5.0-dev and 4.6-dev*/
    {
        gmx_fatal(FARGS, "reading tpx file (%s) version %d with version %d program",
                  gmx_fio_getname(fio), fileVersion, tpx_version);
    }

    gmx_fio_do_int(fio, tpx->natoms);
    gmx_fio_do_int(fio, tpx->ngtc);

    if (fileVersion < 62)
    {
        gmx_fio_do_int(fio, idum);
        gmx_fio_do_real(fio, rdum);
    }
    if (fileVersion >= 79)
    {
        gmx_fio_do_int(fio, tpx->fep_state);
    }
    gmx_fio_do_real(fio, tpx->lambda);
    gmx_fio_do_int(fio, tpx->bIr);
    gmx_fio_do_int(fio, tpx->bTop);
    gmx_fio_do_int(fio, tpx->bX);
    gmx_fio_do_int(fio, tpx->bV);
    gmx_fio_do_int(fio, tpx->bF);
    gmx_fio_do_int(fio, tpx->bBox);

    if ((fileGeneration > tpx_generation))
    {
        /* This can only happen if TopOnlyOK=TRUE */
        tpx->bIr = FALSE;
    }
}

static int do_tpx(t_fileio *fio, gmx_bool bRead,
                  t_inputrec *ir, t_state *state, rvec *x, rvec *v,
                  gmx_mtop_t *mtop)
{
    t_tpxheader     tpx;
    gmx_mtop_t      dum_top;
    gmx_bool        TopOnlyOK;
    int             ePBC;
    gmx_bool        bPeriodicMols;

    if (!bRead)
    {
        GMX_RELEASE_ASSERT(state != nullptr, "Cannot write a null state");
        GMX_RELEASE_ASSERT(x == nullptr && v == nullptr, "Passing separate x and v pointers to do_tpx() is not supported when writing");

        tpx.natoms    = state->natoms;
        tpx.ngtc      = state->ngtc;
        tpx.fep_state = state->fep_state;
        tpx.lambda    = state->lambda[efptFEP];
        tpx.bIr       = (ir       != nullptr);
        tpx.bTop      = (mtop     != nullptr);
        tpx.bX        = (state->flags & (1 << estX));
        tpx.bV        = (state->flags & (1 << estV));
        tpx.bF        = FALSE;
        tpx.bBox      = TRUE;
    }
    else
    {
        GMX_RELEASE_ASSERT(!(x == nullptr && v != nullptr), "Passing x==NULL and v!=NULL is not supported");
    }

    TopOnlyOK = (ir == nullptr);

    int fileVersion;    /* Version number of the code that wrote the file */
    int fileGeneration; /* Generation version number of the code that wrote the file */
    do_tpxheader(fio, bRead, &tpx, TopOnlyOK, &fileVersion, &fileGeneration);

    if (bRead)
    {
        state->flags = 0;
        init_gtc_state(state, tpx.ngtc, 0, 0);
        if (x == nullptr)
        {
            // v is also nullptr by the above assertion, so we may
            // need to make memory in state for storing the contents
            // of the tpx file.
            if (tpx.bX)
            {
                state->flags |= (1 << estX);
            }
            if (tpx.bV)
            {
                state->flags |= (1 << estV);
            }
            state_change_natoms(state, tpx.natoms);
        }
    }

    if (x == nullptr)
    {
        x = as_rvec_array(state->x.data());
        v = as_rvec_array(state->v.data());
    }

#define do_test(fio, b, p) if (bRead && (p != NULL) && !b) gmx_fatal(FARGS, "No %s in %s",#p, gmx_fio_getname(fio))

    do_test(fio, tpx.bBox, state->box);
    if (tpx.bBox)
    {
        gmx_fio_ndo_rvec(fio, state->box, DIM);
        if (fileVersion >= 51)
        {
            gmx_fio_ndo_rvec(fio, state->box_rel, DIM);
        }
        else
        {
            /* We initialize box_rel after reading the inputrec */
            clear_mat(state->box_rel);
        }
        gmx_fio_ndo_rvec(fio, state->boxv, DIM);
        if (fileVersion < 56)
        {
            matrix mdum;
            gmx_fio_ndo_rvec(fio, mdum, DIM);
        }
    }

    if (state->ngtc > 0)
    {
        real *dumv;
        snew(dumv, state->ngtc);
        if (fileVersion < 69)
        {
            gmx_fio_ndo_real(fio, dumv, state->ngtc);
        }
        /* These used to be the Berendsen tcoupl_lambda's */
        gmx_fio_ndo_real(fio, dumv, state->ngtc);
        sfree(dumv);
    }

    do_test(fio, tpx.bTop, mtop);
    if (tpx.bTop)
    {
        if (mtop)
        {
            do_mtop(fio, mtop, bRead, fileVersion);
        }
        else
        {
            do_mtop(fio, &dum_top, bRead, fileVersion);
            done_mtop(&dum_top);
        }
    }
    do_test(fio, tpx.bX, x);
    if (tpx.bX)
    {
        if (bRead)
        {
            state->flags |= (1<<estX);
        }
        gmx_fio_ndo_rvec(fio, x, tpx.natoms);
    }

    do_test(fio, tpx.bV, v);
    if (tpx.bV)
    {
        if (bRead)
        {
            state->flags |= (1<<estV);
        }
        gmx_fio_ndo_rvec(fio, v, tpx.natoms);
    }

    // No need to run do_test when the last argument is NULL
    if (tpx.bF)
    {
        rvec *dummyForces;
        snew(dummyForces, state->natoms);
        gmx_fio_ndo_rvec(fio, dummyForces, tpx.natoms);
        sfree(dummyForces);
    }

    /* Starting with tpx version 26, we have the inputrec
     * at the end of the file, so we can ignore it
     * if the file is never than the software (but still the
     * same generation - see comments at the top of this file.
     *
     *
     */
    ePBC          = -1;
    bPeriodicMols = FALSE;

    do_test(fio, tpx.bIr, ir);
    if (tpx.bIr)
    {
        if (fileVersion >= 53)
        {
            /* Removed the pbc info from do_inputrec, since we always want it */
            if (!bRead)
            {
                ePBC          = ir->ePBC;
                bPeriodicMols = ir->bPeriodicMols;
            }
            gmx_fio_do_int(fio, ePBC);
            gmx_fio_do_gmx_bool(fio, bPeriodicMols);
        }
        if (fileGeneration <= tpx_generation && ir)
        {
            do_inputrec(fio, ir, bRead, fileVersion, mtop ? &mtop->ffparams.fudgeQQ : nullptr);
            if (fileVersion < 51)
            {
                set_box_rel(ir, state);
            }
            if (fileVersion < 53)
            {
                ePBC          = ir->ePBC;
                bPeriodicMols = ir->bPeriodicMols;
            }
        }
        if (bRead && ir && fileVersion >= 53)
        {
            /* We need to do this after do_inputrec, since that initializes ir */
            ir->ePBC          = ePBC;
            ir->bPeriodicMols = bPeriodicMols;
        }
    }

    if (bRead)
    {
        if (tpx.bIr && ir)
        {
            if (state->ngtc == 0)
            {
                /* Reading old version without tcoupl state data: set it */
                init_gtc_state(state, ir->opts.ngtc, 0, ir->opts.nhchainlength);
            }
            if (tpx.bTop && mtop)
            {
                if (fileVersion < 57)
                {
                    if (mtop->moltype[0].ilist[F_DISRES].nr > 0)
                    {
                        ir->eDisre = edrSimple;
                    }
                    else
                    {
                        ir->eDisre = edrNone;
                    }
                }
                set_disres_npair(mtop);
            }
        }

        if (tpx.bTop && mtop)
        {
            gmx_mtop_finalize(mtop);
        }
    }

    return ePBC;
}

static t_fileio *open_tpx(const char *fn, const char *mode)
{
    return gmx_fio_open(fn, mode);
}

static void close_tpx(t_fileio *fio)
{
    gmx_fio_close(fio);
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

void read_tpxheader(const char *fn, t_tpxheader *tpx, gmx_bool TopOnlyOK)
{
    t_fileio *fio;

    fio = open_tpx(fn, "r");
    do_tpxheader(fio, TRUE, tpx, TopOnlyOK, nullptr, nullptr);
    close_tpx(fio);
}

void write_tpx_state(const char *fn,
                     const t_inputrec *ir, const t_state *state, const gmx_mtop_t *mtop)
{
    t_fileio *fio;

    fio = open_tpx(fn, "w");
    do_tpx(fio, FALSE,
           const_cast<t_inputrec *>(ir),
           const_cast<t_state *>(state), nullptr, nullptr,
           const_cast<gmx_mtop_t *>(mtop));
    close_tpx(fio);
}

void read_tpx_state(const char *fn,
                    t_inputrec *ir, t_state *state, gmx_mtop_t *mtop)
{
    t_fileio *fio;

    fio = open_tpx(fn, "r");
    do_tpx(fio, TRUE, ir, state, nullptr, nullptr, mtop);
    close_tpx(fio);
}

int read_tpx(const char *fn,
             t_inputrec *ir, matrix box, int *natoms,
             rvec *x, rvec *v, gmx_mtop_t *mtop)
{
    t_fileio *fio;
    t_state   state;
    int       ePBC;

    fio     = open_tpx(fn, "r");
    ePBC    = do_tpx(fio, TRUE, ir, &state, x, v, mtop);
    close_tpx(fio);
    *natoms = mtop->natoms;
    if (box)
    {
        copy_mat(state.box, box);
    }

    return ePBC;
}

int read_tpx_top(const char *fn,
                 t_inputrec *ir, matrix box, int *natoms,
                 rvec *x, rvec *v, t_topology *top)
{
    gmx_mtop_t  mtop;
    int         ePBC;

    ePBC = read_tpx(fn, ir, box, natoms, x, v, &mtop);

    *top = gmx_mtop_t_to_t_topology(&mtop, true);

    return ePBC;
}

gmx_bool fn2bTPX(const char *file)
{
    return (efTPR == fn2ftp(file));
}

void pr_tpxheader(FILE *fp, int indent, const char *title, const t_tpxheader *sh)
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
    }
}
