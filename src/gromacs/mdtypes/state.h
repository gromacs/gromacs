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
#ifndef GMX_MDTYPES_STATE_H
#define GMX_MDTYPES_STATE_H

#include <array>
#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/*
 * The t_state struct should contain all the (possibly) non-static
 * information required to define the state of the system.
 * Currently the random seeds for SD and BD are missing.
 */

/* These enums are used in flags as (1<<est...).
 * The order of these enums should not be changed,
 * since that affects the checkpoint (.cpt) file format.
 */
enum {
    estLAMBDA,
    estBOX, estBOX_REL, estBOXV, estPRES_PREV, estNH_XI,  estTC_INT,
    estX,   estV,       estSDX_NOTSUPPORTED,  estCGP,
    estLD_RNG_NOTSUPPORTED, estLD_RNGI_NOTSUPPORTED,
    estDISRE_INITF, estDISRE_RM3TAV,
    estORIRE_INITF, estORIRE_DTAV,
    estSVIR_PREV, estNH_VXI, estVETA, estVOL0, estNHPRES_XI, estNHPRES_VXI, estFVIR_PREV,
    estFEPSTATE, estMC_RNG_NOTSUPPORTED, estMC_RNGI_NOTSUPPORTED,
    estNR
};

/* The names of the state entries, defined in src/gmxlib/checkpoint.c */
extern const char *est_names[estNR];

class history_t
{
    public:
        history_t();

        real  disre_initf;  /* The scaling factor for initializing the time av. */
        int   ndisrepairs;  /* The number of distance restraints                */
        real *disre_rm3tav; /* The r^-3 time averaged pair distances            */
        real  orire_initf;  /* The scaling factor for initializing the time av. */
        int   norire_Dtav;  /* The number of matrix element in dtav (npair*5)   */
        real *orire_Dtav;   /* The time averaged orientation tensors            */
};

/* Struct used for checkpointing only.
 * This struct would not be required with unlimited precision.
 * But because of limited precision, the COM motion removal implementation
 * can cause the kinetic energy in the MD loop to differ by a few bits from
 * the kinetic energy one would determine from state.v.
 */
class ekinstate_t
{
    public:
        ekinstate_t();

        gmx_bool             bUpToDate;
        int                  ekin_n;
        tensor              *ekinh;
        tensor              *ekinf;
        tensor              *ekinh_old;
        tensor               ekin_total;
        std::vector<double>  ekinscalef_nhc;
        std::vector<double>  ekinscaleh_nhc;
        std::vector<double>  vscale_nhc;
        real                 dekindl;
        real                 mvcos;
};

typedef struct df_history_t
{
    int      nlambda;        /* total number of lambda states - for history*/

    gmx_bool bEquil;         /* Have we reached equilibration */
    int     *n_at_lam;       /* number of points observed at each lambda */
    real    *wl_histo;       /* histogram for WL flatness determination */
    real     wl_delta;       /* current wang-landau delta */

    real    *sum_weights;    /* weights of the states */
    real    *sum_dg;         /* free energies of the states -- not actually used for weighting, but informational */
    real    *sum_minvar;     /* corrections to weights for minimum variance */
    real    *sum_variance;   /* variances of the states */

    real   **accum_p;        /* accumulated bennett weights for n+1 */
    real   **accum_m;        /* accumulated bennett weights for n-1 */
    real   **accum_p2;       /* accumulated squared bennett weights for n+1 */
    real   **accum_m2;       /* accumulated squared bennett weights for n-1 */

    real   **Tij;            /* transition matrix */
    real   **Tij_empirical;  /* Empirical transition matrix */

} df_history_t;

typedef struct edsamstate_t
{
    /* If one uses essential dynamics or flooding on a group of atoms from
     * more than one molecule, we cannot make this group whole with
     * do_pbc_first_mtop(). We assume that the ED group has the correct PBC
     * representation at the beginning of the simulation and keep track
     * of the shifts to always get it into that representation.
     * For proper restarts from a checkpoint we store the positions of the
     * reference group at the time of checkpoint writing */
    gmx_bool      bFromCpt;     /* Did we start from a checkpoint file?       */
    int           nED;          /* No. of ED/Flooding data sets, if <1 no ED  */
    int          *nref;         /* No. of atoms in i'th reference structure   */
    int          *nav;          /* Same for average structure                 */
    rvec        **old_sref;     /* Positions of the reference atoms
                                   at the last time step (with correct PBC
                                   representation)                            */
    rvec        **old_sref_p;   /* Pointer to these positions                 */
    rvec        **old_sav;      /* Same for the average positions             */
    rvec        **old_sav_p;
}
edsamstate_t;

typedef struct swapstateIons_t
{
    int  nMolReq[eCompNR];                  /* Requested # of molecules per compartment      */
    int *nMolReq_p[eCompNR];                /* Pointer to this data (for .cpt writing)       */
    int  inflow_net[eCompNR];               /* Flux determined from the # of swaps           */
    int *inflow_net_p[eCompNR];             /* Pointer to this data                          */
    int *nMolPast[eCompNR];                 /* Array with nAverage entries for history       */
    int *nMolPast_p[eCompNR];               /* Pointer points to the first entry only        */
    /*                                                                                       */
    /* Channel flux detection, this is counting only and has no influence on whether swaps   */
    /* are performed or not:                                                                 */
    int            fluxfromAtoB[eCompNR];   /* Flux determined from the split cylinders      */
    int           *fluxfromAtoB_p[eCompNR]; /* Pointer to this data                          */
    int            nMol;                    /* # of molecules, size of the following arrays  */
    unsigned char *comp_from;               /* Ion came from which compartment?              */
    unsigned char *channel_label;           /* Through which channel did this ion pass?      */
} swapstateIons_t;

typedef struct swapstate_t
{
    int  eSwapCoords;                       /* Swapping along x, y, or z-direction?          */
    int  nIonTypes;                         /* Number of ion types, this is the size of      */
                                            /* the following arrays                          */
    int  nAverage;                          /* Use average over this many swap attempt       */
                                            /* steps when determining the ion counts         */
    int  fluxleak;                          /* Ions not going through any channel (bad!)     */
    int *fluxleak_p;                        /* Pointer to this data                          */
    /*                                                                                       */
    /* To also make multimeric channel proteins whole, we save the last whole configuration  */
    /* of the channels in the checkpoint file. If we have no checkpoint file, we assume      */
    /* that the starting configuration has the correct PBC representation after making the   */
    /* individual molecules whole                                                            */
    gmx_bool         bFromCpt;                 /* Did we start from a checkpoint file?       */
    int              nat[eChanNR];             /* Size of xc_old_whole, i.e. the number of   */
                                               /* atoms in each channel                      */
    rvec            *xc_old_whole[eChanNR];    /* Last known whole positions of the two      */
                                               /* channels (important for multimeric ch.!)   */
    rvec           **xc_old_whole_p[eChanNR];  /* Pointer to these positions                 */
    swapstateIons_t *ionType;
}
swapstate_t;


class t_state
{
    public:
        // Constructor
        t_state();

        // All things public
        int                      natoms;
        int                      ngtc;
        int                      nnhpres;
        int                      nhchainlength;  /* number of nose-hoover chains               */
        int                      flags;          /* Flags telling which entries are present      */
        int                      fep_state;      /* indicates which of the alchemical states we are in                 */
        std::array<real, efptNR> lambda;         /* lambda vector                               */
        matrix                   box;            /* box vector coordinates                         */
        matrix                   box_rel;        /* Relitaive box vectors to preserve shape        */
        matrix                   boxv;           /* box velocitites for Parrinello-Rahman pcoupl */
        matrix                   pres_prev;      /* Pressure of the previous step for pcoupl  */
        matrix                   svir_prev;      /* Shake virial for previous step for pcoupl */
        matrix                   fvir_prev;      /* Force virial of the previous step for pcoupl  */
        std::vector<double>      nosehoover_xi;  /* for Nose-Hoover tcoupl (ngtc)       */
        std::vector<double>      nosehoover_vxi; /* for N-H tcoupl (ngtc)               */
        std::vector<double>      nhpres_xi;      /* for Nose-Hoover pcoupl for barostat     */
        std::vector<double>      nhpres_vxi;     /* for Nose-Hoover pcoupl for barostat     */
        std::vector<double>      therm_integral; /* for N-H/V-rescale tcoupl (ngtc)     */
        real                     veta;           /* trotter based isotropic P-coupling             */
        real                     vol0;           /* initial volume,required for computing NPT conserverd quantity */
        PaddedRVecVector         x;              /* the coordinates (natoms)                     */
        PaddedRVecVector         v;              /* the velocities (natoms)                      */
        PaddedRVecVector         cg_p;           /* p vector for conjugate gradient minimization */

        ekinstate_t              ekinstate;      /* The state of the kinetic energy data      */

        /* History for special algorithms, should be moved to a history struct */
        history_t               hist;            /* Time history for restraints                  */
        swapstate_t            *swapstate;       /* Position swapping                       */
        df_history_t           *dfhist;          /*Free energy history for free energy analysis  */
        edsamstate_t           *edsamstate;      /* Essential dynamics / flooding history */

        int                     ddp_count;       /* The DD partitioning count for this state  */
        int                     ddp_count_cg_gl; /* The DD part. count for index_gl     */
        std::vector<int>        cg_gl;           /* The global cg number of the local cgs        */
};

typedef struct t_extmass
{
    double *Qinv;  /* inverse mass of thermostat -- computed from inputs, but a good place to store */
    double *QPinv; /* inverse mass of thermostat for barostat -- computed from inputs, but a good place to store */
    double  Winv;  /* Pressure mass inverse -- computed, not input, but a good place to store. Need to make a matrix later */
    tensor  Winvm; /* inverse pressure mass tensor, computed       */
} t_extmass;


typedef struct
{
    real    veta;
    double  rscale;
    double  vscale;
    double  rvscale;
    double  alpha;
    double *vscale_nhc;
} t_vetavars;

//! Allocates memory for temperature coupling
void init_gtc_state(t_state *state, int ngtc, int nnhpres, int nhchainlength);

//! Change the number of atoms represented by this state, allocating memory as needed.
void state_change_natoms(t_state *state, int natoms);

//! Allocates memory for free-energy history
void init_dfhist_state(t_state *state, int dfhistNumLambda);

void comp_state(const t_state *st1, const t_state *st2, gmx_bool bRMSD, real ftol, real abstol);

/*! \brief Allocate an rvec pointer and copy the contents of v to it */
rvec *getRvecArrayFromPaddedRVecVector(const PaddedRVecVector *v,
                                       unsigned int            n);

#endif
