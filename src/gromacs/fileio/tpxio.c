/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define TPX_TAG_RELEASE  "release"

/*! \brief Tag string for the file format written to run input files
 * written by this version of the code.
 *
 * Change this if you want to change the run input format in a feature
 * branch. This ensures that there will not be different run input
 * formats around which cannot be distinguished, while not causing
 * problems rebasing the feature branch onto upstream changes. When
 * merging with mainstream GROMACS, set this tag string back to
 * TPX_TAG_RELEASE, and instead add an element to tpxv and set
 * tpx_version to that.
 */
static const char *tpx_tag = TPX_TAG_RELEASE;

/*! \brief Enum of values that describe the contents of a tpr file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element to the end of this
 * enumeration, change the definition of tpx_version, and write code
 * below that does the right thing according to the value of
 * file_version. */
enum tpxv {
    tpxv_ComputationalElectrophysiology = 96,                /**< support for ion/water position swaps (computational electrophysiology) */
    tpxv_Use64BitRandomSeed,                                 /**< change ld_seed from int to gmx_int64_t */
    tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, /**< potentials for supporting coarse-grained force fields */
    tpxv_InteractiveMolecularDynamics,                       /**< interactive molecular dynamics (IMD) */
    tpxv_RemoveObsoleteParameters1,                          /**< remove optimize_fft, dihre_fc, nstcheckpoint */
    tpxv_PullCoordTypeGeom,                                  /**< add pull type and geometry per group and flat-bottom */
    tpxv_PullGeomDirRel,                                     /**< add pull geometry direction-relative */
    tpxv_IntermolecularBondeds                               /**< permit inter-molecular bonded interactions in the topology */
};

/*! \brief Version number of the file format written to run input
 * files by this version of the code.
 *
 * The tpx_version number should be increased whenever the file format
 * in the main development branch changes, generally to the highest
 * value present in tpxv. Backward compatibility for reading old run
 * input files is maintained by checking this version number against
 * that of the file and then using the correct code path.
 *
 * When developing a feature branch that needs to change the run input
 * file format, change tpx_tag instead. */
static const int tpx_version = tpxv_IntermolecularBondeds;


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
static const int tpx_incompatible_version = 9;



/* Struct used to maintain tpx compatibility when function types are added */
typedef struct {
    int fvnr;  /* file version number in which the function type first appeared */
    int ftype; /* function type */
} t_ftupd;

/*
 * The entries should be ordered in:
 * 1. ascending file version number
 * 2. ascending function type number
 */
/*static const t_ftupd ftupd[] = {
   { 20, F_CUBICBONDS        },
   { 20, F_CONNBONDS         },
   { 20, F_HARMONIC          },
   { 20, F_EQM,              },
   { 22, F_DISRESVIOL        },
   { 22, F_ORIRES            },
   { 22, F_ORIRESDEV         },
   { 26, F_FOURDIHS          },
   { 26, F_PIDIHS            },
   { 26, F_DIHRES            },
   { 26, F_DIHRESVIOL        },
   { 30, F_CROSS_BOND_BONDS  },
   { 30, F_CROSS_BOND_ANGLES },
   { 30, F_UREY_BRADLEY      },
   { 30, F_POLARIZATION      },
   { 54, F_DHDL_CON          },
   };*/
/*
 * The entries should be ordered in:
 * 1. ascending function type number
 * 2. ascending file version number
 */
/* question; what is the purpose of the commented code above? */
static const t_ftupd ftupd[] = {
    { 20, F_CUBICBONDS        },
    { 20, F_CONNBONDS         },
    { 20, F_HARMONIC          },
    { 34, F_FENEBONDS         },
    { 43, F_TABBONDS          },
    { 43, F_TABBONDSNC        },
    { 70, F_RESTRBONDS        },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES },
    { 76, F_LINEAR_ANGLES     },
    { 30, F_CROSS_BOND_BONDS  },
    { 30, F_CROSS_BOND_ANGLES },
    { 30, F_UREY_BRADLEY      },
    { 34, F_QUARTIC_ANGLES    },
    { 43, F_TABANGLES         },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS },
    { 26, F_FOURDIHS          },
    { 26, F_PIDIHS            },
    { 43, F_TABDIHS           },
    { 65, F_CMAP              },
    { 60, F_GB12              },
    { 61, F_GB13              },
    { 61, F_GB14              },
    { 72, F_GBPOL             },
    { 72, F_NPSOLVATION       },
    { 41, F_LJC14_Q           },
    { 41, F_LJC_PAIRS_NB      },
    { 32, F_BHAM_LR           },
    { 32, F_RF_EXCL           },
    { 32, F_COUL_RECIP        },
    { 93, F_LJ_RECIP          },
    { 46, F_DPD               },
    { 30, F_POLARIZATION      },
    { 36, F_THOLE_POL         },
    { 90, F_FBPOSRES          },
    { 22, F_DISRESVIOL        },
    { 22, F_ORIRES            },
    { 22, F_ORIRESDEV         },
    { 26, F_DIHRES            },
    { 26, F_DIHRESVIOL        },
    { 49, F_VSITE4FDN         },
    { 50, F_VSITEN            },
    { 46, F_COM_PULL          },
    { 20, F_EQM               },
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
    int  i;
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
    int      i;

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

static void do_pull_coord(t_fileio *fio, t_pull_coord *pcrd, int file_version,
                          int ePullOld, int eGeomOld, ivec dimOld)
{
    int      i;

    gmx_fio_do_int(fio, pcrd->group[0]);
    gmx_fio_do_int(fio, pcrd->group[1]);
    if (file_version >= tpxv_PullCoordTypeGeom)
    {
        gmx_fio_do_int(fio,  pcrd->eType);
        gmx_fio_do_int(fio,  pcrd->eGeom);
        if (pcrd->eGeom == epullgDIRRELATIVE)
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
    /* i is used in the ndo_double macro*/
    int      i;
    real     fv;
    real     rdum;
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
    int      i, g;
    real     fv;
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
        fepvals->all_lambda   = NULL;
        if (fepvals->init_lambda >= 0)
        {
            fepvals->separate_dvdl[efptFEP] = TRUE;
        }
    }
    if (file_version >= 13)
    {
        gmx_fio_do_real(fio, fepvals->sc_alpha);
    }
    else
    {
        fepvals->sc_alpha = 0;
    }
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
    if (file_version >= 15)
    {
        gmx_fio_do_real(fio, fepvals->sc_sigma);
    }
    else
    {
        fepvals->sc_sigma = 0.3;
    }
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

static void do_pull(t_fileio *fio, pull_params_t *pull, gmx_bool bRead,
                    int file_version, int ePullOld)
{
    int  eGeomOld;
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
        gmx_fio_do_int(fio, pull->bPrintCOM1);
        /* With file_version < 95 this value is set below */
    }
    if (file_version >= tpxv_PullCoordTypeGeom)
    {
        gmx_fio_do_int(fio, pull->bPrintCOM2);
        gmx_fio_do_int(fio, pull->bPrintRefValue);
        gmx_fio_do_int(fio, pull->bPrintComp);
    }
    else
    {
        pull->bPrintCOM2     = FALSE;
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
            do_pullgrp_tpx_pre95(fio, &pull->group[g], &pull->coord[max(g-1, 0)],
                                 bRead, file_version);
            if (g > 0)
            {
                pull->coord[g-1].group[0] = 0;
                pull->coord[g-1].group[1] = g;
            }
        }

        pull->bPrintCOM1 = (pull->group[0].nat > 0);
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
                          file_version, ePullOld, eGeomOld, dimOld);
        }
    }
}


static void do_rotgrp(t_fileio *fio, t_rotgrp *rotg, gmx_bool bRead)
{
    int      i;

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


static void do_swapcoords(t_fileio *fio, t_swapcoords *swap, gmx_bool bRead)
{
    int i, j;


    gmx_fio_do_int(fio, swap->nat);
    gmx_fio_do_int(fio, swap->nat_sol);
    for (j = 0; j < 2; j++)
    {
        gmx_fio_do_int(fio, swap->nat_split[j]);
        gmx_fio_do_int(fio, swap->massw_split[j]);
    }
    gmx_fio_do_int(fio, swap->nstswap);
    gmx_fio_do_int(fio, swap->nAverage);
    gmx_fio_do_real(fio, swap->threshold);
    gmx_fio_do_real(fio, swap->cyl0r);
    gmx_fio_do_real(fio, swap->cyl0u);
    gmx_fio_do_real(fio, swap->cyl0l);
    gmx_fio_do_real(fio, swap->cyl1r);
    gmx_fio_do_real(fio, swap->cyl1u);
    gmx_fio_do_real(fio, swap->cyl1l);

    if (bRead)
    {
        snew(swap->ind, swap->nat);
        snew(swap->ind_sol, swap->nat_sol);
        for (j = 0; j < 2; j++)
        {
            snew(swap->ind_split[j], swap->nat_split[j]);
        }
    }

    gmx_fio_ndo_int(fio, swap->ind, swap->nat);
    gmx_fio_ndo_int(fio, swap->ind_sol, swap->nat_sol);
    for (j = 0; j < 2; j++)
    {
        gmx_fio_ndo_int(fio, swap->ind_split[j], swap->nat_split[j]);
    }

    for (j = 0; j < eCompNR; j++)
    {
        gmx_fio_do_int(fio, swap->nanions[j]);
        gmx_fio_do_int(fio, swap->ncations[j]);
    }

}


static void do_inputrec(t_fileio *fio, t_inputrec *ir, gmx_bool bRead,
                        int file_version, real *fudgeQQ)
{
    int      i, j, k, *tmp, idum = 0;
    real     rdum, bd_temp;
    rvec     vdum;
    gmx_bool bSimAnn, bdum = 0;
    real     zerotemptime, finish_t, init_temp, finish_temp;

    if (file_version != tpx_version)
    {
        /* Give a warning about features that are not accessible */
        fprintf(stderr, "Note: file tpx version %d, software tpx version %d\n",
                file_version, tpx_version);
    }

    if (bRead)
    {
        init_inputrec(ir);
    }

    if (file_version == 0)
    {
        return;
    }

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

    if (file_version > 25)
    {
        if (file_version >= 62)
        {
            gmx_fio_do_int64(fio, ir->init_step);
        }
        else
        {
            gmx_fio_do_int(fio, idum);
            ir->init_step = idum;
        }
    }
    else
    {
        ir->init_step = 0;
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
        if ((file_version <= 15) && (ir->ePBC == 2))
        {
            ir->ePBC = epbcNONE;
        }
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
    if (file_version > 25 && file_version < tpxv_RemoveObsoleteParameters1)
    {
        gmx_fio_do_int(fio, idum);
    }

    gmx_fio_do_int(fio, ir->nstcgsteep);

    if (file_version >= 30)
    {
        gmx_fio_do_int(fio, ir->nbfgscorr);
    }
    else if (bRead)
    {
        ir->nbfgscorr = 10;
    }

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
    if (file_version < 19)
    {
        gmx_fio_do_int(fio, idum);
        gmx_fio_do_int(fio, idum);
    }
    if (file_version < 18)
    {
        gmx_fio_do_int(fio, idum);
    }
    if (file_version >= 81)
    {
        gmx_fio_do_real(fio, ir->verletbuf_tol);
    }
    else
    {
        ir->verletbuf_tol = 0;
    }
    gmx_fio_do_real(fio, ir->rlist);
    if (file_version >= 67)
    {
        gmx_fio_do_real(fio, ir->rlistlong);
    }
    if (file_version >= 82 && file_version != 90)
    {
        gmx_fio_do_int(fio, ir->nstcalclr);
    }
    else
    {
        /* Calculate at NS steps */
        ir->nstcalclr = ir->nstlist;
    }
    gmx_fio_do_int(fio, ir->coulombtype);
    if (file_version < 32 && ir->coulombtype == eelRF)
    {
        ir->coulombtype = eelRF_NEC;
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
    if (file_version < 67)
    {
        ir->rlistlong = max_cutoff(ir->rlist, max_cutoff(ir->rvdw, ir->rcoulomb));
    }
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
    if (file_version >= 29)
    {
        gmx_fio_do_real(fio, ir->tabext);
    }
    else
    {
        ir->tabext = 1.0;
    }

    if (file_version > 25)
    {
        gmx_fio_do_int(fio, ir->gb_algorithm);
        gmx_fio_do_int(fio, ir->nstgbradii);
        gmx_fio_do_real(fio, ir->rgbradii);
        gmx_fio_do_real(fio, ir->gb_saltconc);
        gmx_fio_do_int(fio, ir->implicit_solvent);
    }
    else
    {
        ir->gb_algorithm     = egbSTILL;
        ir->nstgbradii       = 1;
        ir->rgbradii         = 1.0;
        ir->gb_saltconc      = 0;
        ir->implicit_solvent = eisNO;
    }
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

    if (file_version >= 24)
    {
        gmx_fio_do_int(fio, ir->ewald_geometry);
    }
    else
    {
        ir->ewald_geometry = eewg3D;
    }

    if (file_version <= 17)
    {
        ir->epsilon_surface = 0;
        if (file_version == 17)
        {
            gmx_fio_do_int(fio, idum);
        }
    }
    else
    {
        gmx_fio_do_real(fio, ir->epsilon_surface);
    }

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
    if (file_version <= 15)
    {
        gmx_fio_do_int(fio, idum);
    }
    if (file_version <= 17)
    {
        gmx_fio_do_int(fio, ir->epct);
        if (file_version <= 15)
        {
            if (ir->epct == 5)
            {
                ir->epct = epctSURFACETENSION;
            }
            gmx_fio_do_int(fio, idum);
        }
        ir->epct -= 1;
        /* we have removed the NO alternative at the beginning */
        if (ir->epct == -1)
        {
            ir->epc  = epcNO;
            ir->epct = epctISOTROPIC;
        }
        else
        {
            ir->epc = epcBERENDSEN;
        }
    }
    else
    {
        gmx_fio_do_int(fio, ir->epc);
        gmx_fio_do_int(fio, ir->epct);
    }
    if (file_version >= 71)
    {
        gmx_fio_do_int(fio, ir->nstpcouple);
    }
    else
    {
        ir->nstpcouple = ir->nstcalcenergy;
    }
    gmx_fio_do_real(fio, ir->tau_p);
    if (file_version <= 15)
    {
        gmx_fio_do_rvec(fio, vdum);
        clear_mat(ir->ref_p);
        for (i = 0; i < DIM; i++)
        {
            ir->ref_p[i][i] = vdum[i];
        }
    }
    else
    {
        gmx_fio_do_rvec(fio, ir->ref_p[XX]);
        gmx_fio_do_rvec(fio, ir->ref_p[YY]);
        gmx_fio_do_rvec(fio, ir->ref_p[ZZ]);
    }
    if (file_version <= 15)
    {
        gmx_fio_do_rvec(fio, vdum);
        clear_mat(ir->compress);
        for (i = 0; i < DIM; i++)
        {
            ir->compress[i][i] = vdum[i];
        }
    }
    else
    {
        gmx_fio_do_rvec(fio, ir->compress[XX]);
        gmx_fio_do_rvec(fio, ir->compress[YY]);
        gmx_fio_do_rvec(fio, ir->compress[ZZ]);
    }
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
    if ((file_version > 25) && (file_version < 79))
    {
        gmx_fio_do_int(fio, ir->andersen_seed);
    }
    else
    {
        ir->andersen_seed = 0;
    }
    if (file_version < 26)
    {
        gmx_fio_do_gmx_bool(fio, bSimAnn);
        gmx_fio_do_real(fio, zerotemptime);
    }

    if (file_version < 37)
    {
        gmx_fio_do_real(fio, rdum);
    }

    gmx_fio_do_real(fio, ir->shake_tol);
    if (file_version < 54)
    {
        gmx_fio_do_real(fio, *fudgeQQ);
    }

    gmx_fio_do_int(fio, ir->efep);
    if (file_version <= 14 && ir->efep != efepNO)
    {
        ir->efep = efepYES;
    }
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
    if (file_version < 22)
    {
        if (ir->eDisreWeighting == 0)
        {
            ir->eDisreWeighting = edrwEqual;
        }
        else
        {
            ir->eDisreWeighting = edrwConservative;
        }
    }
    gmx_fio_do_gmx_bool(fio, ir->bDisreMixed);
    gmx_fio_do_real(fio, ir->dr_fc);
    gmx_fio_do_real(fio, ir->dr_tau);
    gmx_fio_do_int(fio, ir->nstdisreout);
    if (file_version >= 22)
    {
        gmx_fio_do_real(fio, ir->orires_fc);
        gmx_fio_do_real(fio, ir->orires_tau);
        gmx_fio_do_int(fio, ir->nstorireout);
    }
    else
    {
        ir->orires_fc   = 0;
        ir->orires_tau  = 0;
        ir->nstorireout = 0;
    }

    /* ignore dihre_fc */
    if (file_version >= 26 && file_version < 79)
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
    if (file_version >= 22)
    {
        gmx_fio_do_gmx_bool(fio, ir->bShakeSOR);
    }
    else if (bRead)
    {
        ir->bShakeSOR = TRUE;
    }
    if (file_version >= 11)
    {
        gmx_fio_do_int(fio, ir->niter);
    }
    else if (bRead)
    {
        ir->niter = 25;
        fprintf(stderr, "Note: niter not in run input file, setting it to %d\n",
                ir->niter);
    }
    if (file_version >= 21)
    {
        gmx_fio_do_real(fio, ir->fc_stepsize);
    }
    else
    {
        ir->fc_stepsize = 0;
    }
    gmx_fio_do_int(fio, ir->eConstrAlg);
    gmx_fio_do_int(fio, ir->nProjOrder);
    gmx_fio_do_real(fio, ir->LincsWarnAngle);
    if (file_version <= 14)
    {
        gmx_fio_do_int(fio, idum);
    }
    if (file_version >= 26)
    {
        gmx_fio_do_int(fio, ir->nLincsIter);
    }
    else if (bRead)
    {
        ir->nLincsIter = 1;
        fprintf(stderr, "Note: nLincsIter not in run input file, setting it to %d\n",
                ir->nLincsIter);
    }
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
    if (file_version >= 14)
    {
        gmx_fio_do_real(fio, ir->cos_accel);
    }
    else if (bRead)
    {
        ir->cos_accel = 0;
    }
    gmx_fio_do_int(fio, ir->userint1);
    gmx_fio_do_int(fio, ir->userint2);
    gmx_fio_do_int(fio, ir->userint3);
    gmx_fio_do_int(fio, ir->userint4);
    gmx_fio_do_real(fio, ir->userreal1);
    gmx_fio_do_real(fio, ir->userreal2);
    gmx_fio_do_real(fio, ir->userreal3);
    gmx_fio_do_real(fio, ir->userreal4);

    /* AdResS stuff */
    if (file_version >= 77)
    {
        gmx_fio_do_gmx_bool(fio, ir->bAdress);
        if (ir->bAdress)
        {
            if (bRead)
            {
                snew(ir->adress, 1);
            }
            gmx_fio_do_int(fio, ir->adress->type);
            gmx_fio_do_real(fio, ir->adress->const_wf);
            gmx_fio_do_real(fio, ir->adress->ex_width);
            gmx_fio_do_real(fio, ir->adress->hy_width);
            gmx_fio_do_int(fio, ir->adress->icor);
            gmx_fio_do_int(fio, ir->adress->site);
            gmx_fio_do_rvec(fio, ir->adress->refs);
            gmx_fio_do_int(fio, ir->adress->n_tf_grps);
            gmx_fio_do_real(fio, ir->adress->ex_forcecap);
            gmx_fio_do_int(fio, ir->adress->n_energy_grps);
            gmx_fio_do_int(fio, ir->adress->do_hybridpairs);

            if (bRead)
            {
                snew(ir->adress->tf_table_index, ir->adress->n_tf_grps);
            }
            if (ir->adress->n_tf_grps > 0)
            {
                gmx_fio_ndo_int(fio, ir->adress->tf_table_index, ir->adress->n_tf_grps);
            }
            if (bRead)
            {
                snew(ir->adress->group_explicit, ir->adress->n_energy_grps);
            }
            if (ir->adress->n_energy_grps > 0)
            {
                gmx_fio_ndo_int(fio, ir->adress->group_explicit, ir->adress->n_energy_grps);
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
        if (bRead && file_version < 13)
        {
            snew(tmp, ir->opts.ngtc);
            gmx_fio_ndo_int(fio, tmp, ir->opts.ngtc);
            for (i = 0; i < ir->opts.ngtc; i++)
            {
                ir->opts.nrdf[i] = tmp[i];
            }
            sfree(tmp);
        }
        else
        {
            gmx_fio_ndo_real(fio, ir->opts.nrdf, ir->opts.ngtc);
        }
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
    if (file_version >= 12)
    {
        gmx_fio_ndo_int(fio, ir->opts.egp_flags,
                        ir->opts.ngener*ir->opts.ngener);
    }

    if (bRead && file_version < 26)
    {
        for (i = 0; i < ir->opts.ngtc; i++)
        {
            if (bSimAnn)
            {
                ir->opts.annealing[i]      = eannSINGLE;
                ir->opts.anneal_npoints[i] = 2;
                snew(ir->opts.anneal_time[i], 2);
                snew(ir->opts.anneal_temp[i], 2);
                /* calculate the starting/ending temperatures from reft, zerotemptime, and nsteps */
                finish_t                   = ir->init_t + ir->nsteps * ir->delta_t;
                init_temp                  = ir->opts.ref_t[i]*(1-ir->init_t/zerotemptime);
                finish_temp                = ir->opts.ref_t[i]*(1-finish_t/zerotemptime);
                ir->opts.anneal_time[i][0] = ir->init_t;
                ir->opts.anneal_time[i][1] = finish_t;
                ir->opts.anneal_temp[i][0] = init_temp;
                ir->opts.anneal_temp[i][1] = finish_temp;
            }
            else
            {
                ir->opts.annealing[i]      = eannNO;
                ir->opts.anneal_npoints[i] = 0;
            }
        }
    }
    else
    {
        /* file version 26 or later */
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
    for (j = 0; (j < DIM); j++)
    {
        gmx_fio_do_int(fio, ir->ex[j].n);
        gmx_fio_do_int(fio, ir->et[j].n);
        if (bRead)
        {
            snew(ir->ex[j].a,  ir->ex[j].n);
            snew(ir->ex[j].phi, ir->ex[j].n);
            snew(ir->et[j].a,  ir->et[j].n);
            snew(ir->et[j].phi, ir->et[j].n);
        }
        gmx_fio_ndo_real(fio, ir->ex[j].a,  ir->ex[j].n);
        gmx_fio_ndo_real(fio, ir->ex[j].phi, ir->ex[j].n);
        gmx_fio_ndo_real(fio, ir->et[j].a,  ir->et[j].n);
        gmx_fio_ndo_real(fio, ir->et[j].phi, ir->et[j].n);
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
            do_swapcoords(fio, ir->swap, bRead);
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
            snew(ir->opts.bOPT,        ir->opts.ngQM);
            snew(ir->opts.bTS,         ir->opts.ngQM);
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
            gmx_fio_ndo_gmx_bool(fio, ir->opts.bOPT, ir->opts.ngQM);
            gmx_fio_ndo_gmx_bool(fio, ir->opts.bTS, ir->opts.ngQM);
        }
        /* end of QMMM stuff */
    }
}


static void do_harm(t_fileio *fio, t_iparams *iparams)
{
    gmx_fio_do_real(fio, iparams->harmonic.rA);
    gmx_fio_do_real(fio, iparams->harmonic.krA);
    gmx_fio_do_real(fio, iparams->harmonic.rB);
    gmx_fio_do_real(fio, iparams->harmonic.krB);
}

void do_iparams(t_fileio *fio, t_functype ftype, t_iparams *iparams,
                gmx_bool bRead, int file_version)
{
    int      idum;
    real     rdum;

    if (!bRead)
    {
        gmx_fio_set_comment(fio, interaction_function[ftype].name);
    }
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
            if (file_version < 31)
            {
                gmx_fatal(FARGS, "Old tpr files with water_polarization not supported. Make a new.");
            }
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
            if (bRead && file_version < 27)
            {
                copy_rvec(iparams->posres.pos0A, iparams->posres.pos0B);
                copy_rvec(iparams->posres.fcA, iparams->posres.fcB);
            }
            else
            {
                gmx_fio_do_rvec(fio, iparams->posres.pos0B);
                gmx_fio_do_rvec(fio, iparams->posres.fcB);
            }
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
            if (file_version >= 25)
            {
                gmx_fio_ndo_real(fio, iparams->rbdihs.rbcB, NR_RBDIHS);
            }
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
    if (!bRead)
    {
        gmx_fio_unset_comment(fio);
    }
}

static void do_ilist(t_fileio *fio, t_ilist *ilist, gmx_bool bRead, int file_version,
                     int ftype)
{
    int      i, k, idum;

    if (!bRead)
    {
        gmx_fio_set_comment(fio, interaction_function[ftype].name);
    }
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
    if (!bRead)
    {
        gmx_fio_unset_comment(fio);
    }
}

static void do_ffparams(t_fileio *fio, gmx_ffparams_t *ffparams,
                        gmx_bool bRead, int file_version)
{
    int          idum, i, j;
    unsigned int k;

    gmx_fio_do_int(fio, ffparams->atnr);
    if (file_version < 57)
    {
        gmx_fio_do_int(fio, idum);
    }
    gmx_fio_do_int(fio, ffparams->ntypes);
    if (bRead && debug)
    {
        fprintf(debug, "ffparams->atnr = %d, ntypes = %d\n",
                ffparams->atnr, ffparams->ntypes);
    }
    if (bRead)
    {
        snew(ffparams->functype, ffparams->ntypes);
        snew(ffparams->iparams, ffparams->ntypes);
    }
    /* Read/write all the function types */
    gmx_fio_ndo_int(fio, ffparams->functype, ffparams->ntypes);
    if (bRead && debug)
    {
        pr_ivec(debug, 0, "functype", ffparams->functype, ffparams->ntypes, TRUE);
    }

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
                    if (debug)
                    {
                        fprintf(debug, "Incrementing function type %d to %d (due to %s)\n",
                                i, ffparams->functype[i],
                                interaction_function[ftupd[k].ftype].longname);
                        fflush(debug);
                    }
                }
            }
        }

        do_iparams(fio, ffparams->functype[i], &ffparams->iparams[i], bRead,
                   file_version);
        if (bRead && debug)
        {
            pr_iparams(debug, ffparams->functype[i], &ffparams->iparams[i]);
        }
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
    int          i, j, renum[F_NRE];
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
            ilist[j].iatoms = NULL;
        }
        else
        {
            do_ilist(fio, &ilist[j], bRead, file_version, j);
            if (file_version < 78 && j == F_SETTLE && ilist[j].nr > 0)
            {
                add_settle_atoms(&ilist[j]);
            }
        }
        /*
           if (bRead && gmx_debug_at)
           pr_ilist(debug,0,interaction_function[j].longname,
               functype,&ilist[j],TRUE);
         */
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
        if ((block->nalloc_index > 0) && (NULL != block->index))
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
    char * p_elem;

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
            strncpy(atom->elem, atomicnumber_to_element(atom->atomnumber), 4);
            /* avoid warnings about potentially unterminated string */
            atom->elem[3] = '\0';
        }
    }
    else if (bRead)
    {
        atom->atomnumber = NOTSET;
    }
    if (file_version < 23)
    {
        myngrp = 8;
    }
    else if (file_version < 39)
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
    int      i, j, myngrp;

    if (file_version < 23)
    {
        myngrp = 8;
    }
    else if (file_version < 39)
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
        snew(atoms->atom, atoms->nr);
        snew(atoms->atomname, atoms->nr);
        snew(atoms->atomtype, atoms->nr);
        snew(atoms->atomtypeB, atoms->nr);
        snew(atoms->resinfo, atoms->nres);
        if (file_version < 57)
        {
            snew(groups->grpname, groups->ngrpname);
        }
        atoms->pdbinfo = NULL;
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        do_atom(fio, &atoms->atom[i], egcNR, bRead, file_version, groups, i);
    }
    do_strstr(fio, atoms->nr, atoms->atomname, bRead, symtab);
    if (bRead && (file_version <= 20))
    {
        for (i = 0; i < atoms->nr; i++)
        {
            atoms->atomtype[i]  = put_symtab(symtab, "?");
            atoms->atomtypeB[i] = put_symtab(symtab, "?");
        }
    }
    else
    {
        do_strstr(fio, atoms->nr, atoms->atomtype, bRead, symtab);
        do_strstr(fio, atoms->nr, atoms->atomtypeB, bRead, symtab);
    }
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
    int      g, n, i;

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
                groups->grpnr[g] = NULL;
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
    int      i, j;

    if (file_version > 25)
    {
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
    else
    {
        /* File versions prior to 26 cannot do GBSA,
         * so they dont use this structure
         */
        atomtypes->nr         = 0;
        atomtypes->radius     = NULL;
        atomtypes->vol        = NULL;
        atomtypes->surftens   = NULL;
        atomtypes->atomnumber = NULL;
        atomtypes->gb_radius  = NULL;
        atomtypes->S_hct      = NULL;
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
        while (symbuf != NULL)
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


void tpx_make_chain_identifiers(t_atoms *atoms, t_block *mols)
{
    int  m, a, a0, a1, r;
    char c, chainid;
    int  chainnum;

    /* We always assign a new chain number, but save the chain id characters
     * for larger molecules.
     */
#define CHAIN_MIN_ATOMS 15

    chainnum = 0;
    chainid  = 'A';
    for (m = 0; m < mols->nr; m++)
    {
        a0 = mols->index[m];
        a1 = mols->index[m+1];
        if ((a1-a0 >= CHAIN_MIN_ATOMS) && (chainid <= 'Z'))
        {
            c = chainid;
            chainid++;
        }
        else
        {
            c = ' ';
        }
        for (a = a0; a < a1; a++)
        {
            atoms->resinfo[atoms->atom[a].resind].chainnum = chainnum;
            atoms->resinfo[atoms->atom[a].resind].chainid  = c;
        }
        chainnum++;
    }

    /* Blank out the chain id if there was only one chain */
    if (chainid == 'B')
    {
        for (r = 0; r < atoms->nres; r++)
        {
            atoms->resinfo[r].chainid = ' ';
        }
    }
}

static void do_moltype(t_fileio *fio, gmx_moltype_t *molt, gmx_bool bRead,
                       t_symtab *symtab, int file_version,
                       gmx_groups_t *groups)
{
    int i;

    if (file_version >= 57)
    {
        do_symstr(fio, &(molt->name), bRead, symtab);
    }

    do_atoms(fio, &molt->atoms, bRead, symtab, file_version, groups);

    if (bRead && gmx_debug_at)
    {
        pr_atoms(debug, 0, "atoms", &molt->atoms, TRUE);
    }

    if (file_version >= 57)
    {
        do_ilists(fio, molt->ilist, bRead, file_version);

        do_block(fio, &molt->cgs, bRead, file_version);
        if (bRead && gmx_debug_at)
        {
            pr_block(debug, 0, "cgs", &molt->cgs, TRUE);
        }
    }

    /* This used to be in the atoms struct */
    do_blocka(fio, &molt->excls, bRead, file_version);
}

static void do_molblock(t_fileio *fio, gmx_molblock_t *molb, gmx_bool bRead)
{
    int i;

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
        am = max(am, il->iatoms[i+1]);
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
            ip = &mtop->ffparams.iparams[ilfb->iatoms[i]];
            am = max(am, ilfb->iatoms[i+1]);
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
    int            mt, mb, i;
    t_blocka       dumb;

    if (bRead)
    {
        init_mtop(mtop);
    }
    do_symtab(fio, &(mtop->symtab), bRead);
    if (bRead && debug)
    {
        pr_symtab(debug, 0, "symtab", &mtop->symtab);
    }

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
    if (bRead && debug)
    {
        pr_atomtypes(debug, 0, "atomtypes", &mtop->atomtypes, TRUE);
    }

    if (file_version < 57)
    {
        /* Debug statements are inside do_idef */
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
        mtop->ffparams.cmap_grid.cmapdata     = NULL;
    }

    if (file_version >= 57)
    {
        do_groups(fio, &mtop->groups, bRead, &(mtop->symtab), file_version);
    }

    if (file_version < 57)
    {
        do_block(fio, &mtop->moltype[0].cgs, bRead, file_version);
        if (bRead && gmx_debug_at)
        {
            pr_block(debug, 0, "cgs", &mtop->moltype[0].cgs, TRUE);
        }
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
        if (gmx_debug_at)
        {
            pr_block(debug, 0, "mols", &mtop->mols, TRUE);
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
 * of tpx files, provided the file_generation hasn't changed.
 * If it is FALSE, we need the inputrecord too, and bail out
 * if the file is newer than the program.
 *
 * The version and generation if the topology (see top of this file)
 * are returned in the two last arguments.
 *
 * If possible, we will read the inputrec even when TopOnlyOK is TRUE.
 */
static void do_tpxheader(t_fileio *fio, gmx_bool bRead, t_tpxheader *tpx,
                         gmx_bool TopOnlyOK, int *file_version,
                         int *file_generation)
{
    char      buf[STRLEN];
    char      file_tag[STRLEN];
    gmx_bool  bDouble;
    int       precision;
    int       fver, fgen;
    int       idum = 0;
    real      rdum = 0;

    gmx_fio_checktype(fio);
    gmx_fio_setdebug(fio, bDebugMode());

    /* XDR binary topology file */
    precision = sizeof(real);
    if (bRead)
    {
        gmx_fio_do_string(fio, buf);
        if (strncmp(buf, "VERSION", 7))
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
        gmx_fio_write_string(fio, GromacsVersion());
        bDouble = (precision == sizeof(double));
        gmx_fio_setprecision(fio, bDouble);
        gmx_fio_do_int(fio, precision);
        fver = tpx_version;
        sprintf(file_tag, "%s", tpx_tag);
        fgen = tpx_generation;
    }

    /* Check versions! */
    gmx_fio_do_int(fio, fver);

    /* This is for backward compatibility with development versions 77-79
     * where the tag was, mistakenly, placed before the generation,
     * which would cause a segv instead of a proper error message
     * when reading the topology only from tpx with <77 code.
     */
    if (fver >= 77 && fver <= 79)
    {
        gmx_fio_do_string(fio, file_tag);
    }

    if (fver >= 26)
    {
        gmx_fio_do_int(fio, fgen);
    }
    else
    {
        fgen = 0;
    }

    if (fver >= 81)
    {
        gmx_fio_do_string(fio, file_tag);
    }
    if (bRead)
    {
        if (fver < 77)
        {
            /* Versions before 77 don't have the tag, set it to release */
            sprintf(file_tag, "%s", TPX_TAG_RELEASE);
        }

        if (strcmp(file_tag, tpx_tag) != 0)
        {
            fprintf(stderr, "Note: file tpx tag '%s', software tpx tag '%s'\n",
                    file_tag, tpx_tag);

            /* We only support reading tpx files with the same tag as the code
             * or tpx files with the release tag and with lower version number.
             */
            if (strcmp(file_tag, TPX_TAG_RELEASE) != 0 && fver < tpx_version)
            {
                gmx_fatal(FARGS, "tpx tag/version mismatch: reading tpx file (%s) version %d, tag '%s' with program for tpx version %d, tag '%s'",
                          gmx_fio_getname(fio), fver, file_tag,
                          tpx_version, tpx_tag);
            }
        }
    }

    if (file_version != NULL)
    {
        *file_version = fver;
    }
    if (file_generation != NULL)
    {
        *file_generation = fgen;
    }


    if ((fver <= tpx_incompatible_version) ||
        ((fver > tpx_version) && !TopOnlyOK) ||
        (fgen > tpx_generation) ||
        tpx_version == 80) /*80 was used by both 5.0-dev and 4.6-dev*/
    {
        gmx_fatal(FARGS, "reading tpx file (%s) version %d with version %d program",
                  gmx_fio_getname(fio), fver, tpx_version);
    }

    gmx_fio_do_int(fio, tpx->natoms);
    if (fver >= 28)
    {
        gmx_fio_do_int(fio, tpx->ngtc);
    }
    else
    {
        tpx->ngtc = 0;
    }
    if (fver < 62)
    {
        gmx_fio_do_int(fio, idum);
        gmx_fio_do_real(fio, rdum);
    }
    /*a better decision will eventually (5.0 or later) need to be made
       on how to treat the alchemical state of the system, which can now
       vary through a simulation, and cannot be completely described
       though a single lambda variable, or even a single state
       index. Eventually, should probably be a vector. MRS*/
    if (fver >= 79)
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

    if ((fgen > tpx_generation))
    {
        /* This can only happen if TopOnlyOK=TRUE */
        tpx->bIr = FALSE;
    }
}

static int do_tpx(t_fileio *fio, gmx_bool bRead,
                  t_inputrec *ir, t_state *state, rvec *f, gmx_mtop_t *mtop,
                  gmx_bool bXVallocated)
{
    t_tpxheader     tpx;
    t_inputrec      dum_ir;
    gmx_mtop_t      dum_top;
    gmx_bool        TopOnlyOK;
    int             file_version, file_generation;
    int             i;
    rvec           *xptr, *vptr;
    int             ePBC;
    gmx_bool        bPeriodicMols;

    if (!bRead)
    {
        tpx.natoms    = state->natoms;
        tpx.ngtc      = state->ngtc; /* need to add nnhpres here? */
        tpx.fep_state = state->fep_state;
        tpx.lambda    = state->lambda[efptFEP];
        tpx.bIr       = (ir       != NULL);
        tpx.bTop      = (mtop     != NULL);
        tpx.bX        = (state->x != NULL);
        tpx.bV        = (state->v != NULL);
        tpx.bF        = (f        != NULL);
        tpx.bBox      = TRUE;
    }

    TopOnlyOK = (ir == NULL);

    do_tpxheader(fio, bRead, &tpx, TopOnlyOK, &file_version, &file_generation);

    if (bRead)
    {
        state->flags  = 0;
        /* state->lambda = tpx.lambda;*/ /*remove this eventually? */
        /* The init_state calls initialize the Nose-Hoover xi integrals to zero */
        if (bXVallocated)
        {
            xptr = state->x;
            vptr = state->v;
            init_state(state, 0, tpx.ngtc, 0, 0, 0); /* nose-hoover chains */ /* eventually, need to add nnhpres here? */
            state->natoms = tpx.natoms;
            state->nalloc = tpx.natoms;
            state->x      = xptr;
            state->v      = vptr;
        }
        else
        {
            init_state(state, tpx.natoms, tpx.ngtc, 0, 0, 0); /* nose-hoover chains */
        }
    }

#define do_test(fio, b, p) if (bRead && (p != NULL) && !b) gmx_fatal(FARGS, "No %s in %s",#p, gmx_fio_getname(fio))

    do_test(fio, tpx.bBox, state->box);
    if (tpx.bBox)
    {
        gmx_fio_ndo_rvec(fio, state->box, DIM);
        if (file_version >= 51)
        {
            gmx_fio_ndo_rvec(fio, state->box_rel, DIM);
        }
        else
        {
            /* We initialize box_rel after reading the inputrec */
            clear_mat(state->box_rel);
        }
        if (file_version >= 28)
        {
            gmx_fio_ndo_rvec(fio, state->boxv, DIM);
            if (file_version < 56)
            {
                matrix mdum;
                gmx_fio_ndo_rvec(fio, mdum, DIM);
            }
        }
    }

    if (state->ngtc > 0 && file_version >= 28)
    {
        real *dumv;
        /*ndo_double(state->nosehoover_xi,state->ngtc,bDum);*/
        /*ndo_double(state->nosehoover_vxi,state->ngtc,bDum);*/
        /*ndo_double(state->therm_integral,state->ngtc,bDum);*/
        snew(dumv, state->ngtc);
        if (file_version < 69)
        {
            gmx_fio_ndo_real(fio, dumv, state->ngtc);
        }
        /* These used to be the Berendsen tcoupl_lambda's */
        gmx_fio_ndo_real(fio, dumv, state->ngtc);
        sfree(dumv);
    }

    /* Prior to tpx version 26, the inputrec was here.
     * I moved it to enable partial forward-compatibility
     * for analysis/viewer programs.
     */
    if (file_version < 26)
    {
        do_test(fio, tpx.bIr, ir);
        if (tpx.bIr)
        {
            if (ir)
            {
                do_inputrec(fio, ir, bRead, file_version,
                            mtop ? &mtop->ffparams.fudgeQQ : NULL);
                if (bRead && debug)
                {
                    pr_inputrec(debug, 0, "inputrec", ir, FALSE);
                }
            }
            else
            {
                do_inputrec(fio, &dum_ir, bRead, file_version,
                            mtop ? &mtop->ffparams.fudgeQQ : NULL);
                if (bRead && debug)
                {
                    pr_inputrec(debug, 0, "inputrec", &dum_ir, FALSE);
                }
                done_inputrec(&dum_ir);
            }

        }
    }

    do_test(fio, tpx.bTop, mtop);
    if (tpx.bTop)
    {
        if (mtop)
        {
            do_mtop(fio, mtop, bRead, file_version);
        }
        else
        {
            do_mtop(fio, &dum_top, bRead, file_version);
            done_mtop(&dum_top, TRUE);
        }
    }
    do_test(fio, tpx.bX, state->x);
    if (tpx.bX)
    {
        if (bRead)
        {
            state->flags |= (1<<estX);
        }
        gmx_fio_ndo_rvec(fio, state->x, state->natoms);
    }

    do_test(fio, tpx.bV, state->v);
    if (tpx.bV)
    {
        if (bRead)
        {
            state->flags |= (1<<estV);
        }
        gmx_fio_ndo_rvec(fio, state->v, state->natoms);
    }

    do_test(fio, tpx.bF, f);
    if (tpx.bF)
    {
        gmx_fio_ndo_rvec(fio, f, state->natoms);
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
    if (file_version >= 26)
    {
        do_test(fio, tpx.bIr, ir);
        if (tpx.bIr)
        {
            if (file_version >= 53)
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
            if (file_generation <= tpx_generation && ir)
            {
                do_inputrec(fio, ir, bRead, file_version, mtop ? &mtop->ffparams.fudgeQQ : NULL);
                if (bRead && debug)
                {
                    pr_inputrec(debug, 0, "inputrec", ir, FALSE);
                }
                if (file_version < 51)
                {
                    set_box_rel(ir, state);
                }
                if (file_version < 53)
                {
                    ePBC          = ir->ePBC;
                    bPeriodicMols = ir->bPeriodicMols;
                }
            }
            if (bRead && ir && file_version >= 53)
            {
                /* We need to do this after do_inputrec, since that initializes ir */
                ir->ePBC          = ePBC;
                ir->bPeriodicMols = bPeriodicMols;
            }
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
                if (file_version < 57)
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

        if (file_version >= 57)
        {
            char *env;
            int   ienv;
            env = getenv("GMX_NOCHARGEGROUPS");
            if (env != NULL)
            {
                sscanf(env, "%d", &ienv);
                fprintf(stderr, "\nFound env.var. GMX_NOCHARGEGROUPS = %d\n",
                        ienv);
                if (ienv > 0)
                {
                    fprintf(stderr,
                            "Will make single atomic charge groups in non-solvent%s\n",
                            ienv > 1 ? " and solvent" : "");
                    gmx_mtop_make_atomic_charge_groups(mtop, ienv == 1);
                }
                fprintf(stderr, "\n");
            }
        }
    }

    return ePBC;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

t_fileio *open_tpx(const char *fn, const char *mode)
{
    return gmx_fio_open(fn, mode);
}

void close_tpx(t_fileio *fio)
{
    gmx_fio_close(fio);
}

void read_tpxheader(const char *fn, t_tpxheader *tpx, gmx_bool TopOnlyOK,
                    int *file_version, int *file_generation)
{
    t_fileio *fio;

    fio = open_tpx(fn, "r");
    do_tpxheader(fio, TRUE, tpx, TopOnlyOK, file_version, file_generation);
    close_tpx(fio);
}

void write_tpx_state(const char *fn,
                     t_inputrec *ir, t_state *state, gmx_mtop_t *mtop)
{
    t_fileio *fio;

    fio = open_tpx(fn, "w");
    do_tpx(fio, FALSE, ir, state, NULL, mtop, FALSE);
    close_tpx(fio);
}

void read_tpx_state(const char *fn,
                    t_inputrec *ir, t_state *state, rvec *f, gmx_mtop_t *mtop)
{
    t_fileio *fio;

    fio = open_tpx(fn, "r");
    do_tpx(fio, TRUE, ir, state, f, mtop, FALSE);
    close_tpx(fio);
}

int read_tpx(const char *fn,
             t_inputrec *ir, matrix box, int *natoms,
             rvec *x, rvec *v, rvec *f, gmx_mtop_t *mtop)
{
    t_fileio *fio;
    t_state   state;
    int       ePBC;

    state.x = x;
    state.v = v;
    fio     = open_tpx(fn, "r");
    ePBC    = do_tpx(fio, TRUE, ir, &state, f, mtop, TRUE);
    close_tpx(fio);
    *natoms = state.natoms;
    if (box)
    {
        copy_mat(state.box, box);
    }
    state.x = NULL;
    state.v = NULL;
    done_state(&state);

    return ePBC;
}

int read_tpx_top(const char *fn,
                 t_inputrec *ir, matrix box, int *natoms,
                 rvec *x, rvec *v, rvec *f, t_topology *top)
{
    gmx_mtop_t  mtop;
    t_topology *ltop;
    int         ePBC;

    ePBC = read_tpx(fn, ir, box, natoms, x, v, f, &mtop);

    *top = gmx_mtop_t_to_t_topology(&mtop);

    return ePBC;
}

gmx_bool fn2bTPX(const char *file)
{
    return (efTPR == fn2ftp(file));
}

static void done_gmx_groups_t(gmx_groups_t *g)
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        if (NULL != g->grps[i].nm_ind)
        {
            sfree(g->grps[i].nm_ind);
            g->grps[i].nm_ind = NULL;
        }
        if (NULL != g->grpnr[i])
        {
            sfree(g->grpnr[i]);
            g->grpnr[i] = NULL;
        }
    }
    /* The contents of this array is in symtab, don't free it here */
    sfree(g->grpname);
}

gmx_bool read_tps_conf(const char *infile, char *title, t_topology *top, int *ePBC,
                       rvec **x, rvec **v, matrix box, gmx_bool bMass)
{
    t_tpxheader      header;
    int              natoms, i, version, generation;
    gmx_bool         bTop, bXNULL = FALSE;
    gmx_mtop_t      *mtop;
    t_topology      *topconv;
    gmx_atomprop_t   aps;

    bTop  = fn2bTPX(infile);
    *ePBC = -1;
    if (bTop)
    {
        read_tpxheader(infile, &header, TRUE, &version, &generation);
        if (x)
        {
            snew(*x, header.natoms);
        }
        if (v)
        {
            snew(*v, header.natoms);
        }
        snew(mtop, 1);
        *ePBC = read_tpx(infile, NULL, box, &natoms,
                         (x == NULL) ? NULL : *x, (v == NULL) ? NULL : *v, NULL, mtop);
        *top = gmx_mtop_t_to_t_topology(mtop);
        /* In this case we need to throw away the group data too */
        done_gmx_groups_t(&mtop->groups);
        sfree(mtop);
        strcpy(title, *top->name);
        tpx_make_chain_identifiers(&top->atoms, &top->mols);
    }
    else
    {
        get_stx_coordnum(infile, &natoms);
        init_t_atoms(&top->atoms, natoms, (fn2ftp(infile) == efPDB));
        if (x == NULL)
        {
            snew(x, 1);
            bXNULL = TRUE;
        }
        snew(*x, natoms);
        if (v)
        {
            snew(*v, natoms);
        }
        read_stx_conf(infile, title, &top->atoms, *x, (v == NULL) ? NULL : *v, ePBC, box);
        if (bXNULL)
        {
            sfree(*x);
            sfree(x);
        }
        if (bMass)
        {
            aps = gmx_atomprop_init();
            for (i = 0; (i < natoms); i++)
            {
                if (!gmx_atomprop_query(aps, epropMass,
                                        *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                        *top->atoms.atomname[i],
                                        &(top->atoms.atom[i].m)))
                {
                    if (debug)
                    {
                        fprintf(debug, "Can not find mass for atom %s %d %s, setting to 1\n",
                                *top->atoms.resinfo[top->atoms.atom[i].resind].name,
                                top->atoms.resinfo[top->atoms.atom[i].resind].nr,
                                *top->atoms.atomname[i]);
                    }
                }
            }
            gmx_atomprop_destroy(aps);
        }
        top->idef.ntypes = -1;
    }

    return bTop;
}
