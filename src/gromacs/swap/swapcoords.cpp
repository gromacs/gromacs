/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
/*! \internal \file
 * \brief
 * Implements functions in swapcoords.h.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_swap
 */
#include "gmxpre.h"

#include "swapcoords.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

static const char *SwS      = {"SWAP:"};                                           /**< For output that comes from the swap module */
static const char *SwSEmpty = {"     "};                                           /**< Placeholder for multi-line output */
static const char* IonString[eIonNR] = {"anion", "cation" };                       /**< Type of ion, used for verbose output */
static const char* IonStr[eIonNR]    = {"-", "+"      };                           /**< Type of ion, used for short output */
static const char* CompStr[eCompNR] = {"A", "B" };                                 /**< Compartment name */
static const char *SwapStr[eSwapTypesNR+1] = { "", "X-", "Y-", "Z-", NULL};        /**< Name for the swap types. */
static const char *DimStr[DIM+1] = { "X", "Y", "Z", NULL};                         /**< Name for the swap dimension. */

/* eGrpSplit0 and eGrpSplit1 _must_ be neighbors in this list because
 * we sometimes loop from eGrpSplit0 to eGrpSplit1 */
enum {
    eGrpIons, eGrpSplit0, eGrpSplit1, eGrpSolvent, eGrpNr
};                                                                               /**< Group identifier */
static const char* GrpString[eGrpNr] = { "ion", "split0", "split1", "solvent" }; /**< Group name */

/** Keep track of through which channel the ions have passed */
enum eChannelHistory {
    eChHistPassedNone, eChHistPassedCh0, eChHistPassedCh1, eChHistNr
};
static const char* ChannelString[eChHistNr] = { "none", "channel0", "channel1" };  /**< Name for the channels */

/*! \brief Domain identifier.
 *
 * Keeps track of from which compartment the ions came before passing the
 * channel.
 */
enum eDomain {
    eDomainNotset, eDomainA, eDomainB, eDomainNr
};
static const char* DomainString[eDomainNr] = { "not_assigned", "Domain_A", "Domain_B" }; /**< Name for the domains */



/*! \internal \brief
 * Structure containing compartment-specific data.
 */
typedef struct swap_compartment
{
    int                nat;                   /**< Number of atoms matching the
                                                 compartment conditions.                         */
    int                nat_old;               /**< Number of atoms before swapping.              */
    int                nat_req;               /**< Requested number of atoms.                    */
    real               nat_av;                /**< Time-averaged number of atoms matching
                                                   the compartment conditions.                   */
    int               *nat_past;              /**< Past ion counts for time-averaging.           */
    int               *ind;                   /**< Indices to coll array of atoms.               */
    real              *dist;                  /**< Distance of atom to compartment center.       */
    int                nalloc;                /**< Allocation size for ind array.                */
    int                inflow_netto;          /**< Net inflow of ions into this compartment.     */
} t_compartment;


/*! \internal \brief
 * This structure contains data needed for each of the groups involved in swapping: ions, water,
 * and channels.
 */
typedef struct swap_group
{
    int               nat;                    /**< Number of atoms in the group                    */
    int               apm;                    /**< Number of atoms in each molecule                */
    atom_id          *ind;                    /**< Global atom indices of the group                */
    atom_id          *ind_loc;                /**< Local atom indices of the group                 */
    int               nat_loc;                /**< Number of local group atoms                     */
    int               nalloc_loc;             /**< Allocation size for ind_loc                     */
    rvec             *xc;                     /**< Collective array of group atom positions        */
    ivec             *xc_shifts;              /**< Current (collective) shifts                     */
    ivec             *xc_eshifts;             /**< Extra shifts since last DD step                 */
    rvec             *xc_old;                 /**< Old (collective) positions                      */
    real             *qc;                     /**< Collective array of charges                     */
    int              *c_ind_loc;              /**< Position of local atoms in the
                                                   collective array, [0..nat_loc]                  */
    real             *m;                      /**< Masses (can be omitted)                         */
    unsigned char    *comp_from;              /**< (Collective) Stores from which compartment this
                                                   atom has come. This way we keep track of through
                                                   which channel an ion permeates (only used for
                                                   the ion group)                                  */
    unsigned char    *comp_now;               /**< In which compartment this ion is now            */
    unsigned char    *channel_label;          /**< Which channel was passed at last by this ion?   */
    rvec              center;                 /**< Center of the group; COM if masses are used     */
} t_group;


/*! \internal \brief
 * Main (private) data structure for the position swapping protocol.
 */
typedef struct t_swap
{
    int               swapdim;                       /**< One of XX, YY, ZZ                               */
    t_pbc            *pbc;                           /**< Needed to make molecules whole.                 */
    FILE             *fpout;                         /**< Output file.                                    */
    t_group           group[eGrpNr];                 /**< Ions, solvent or channels?                      */
    t_compartment     comp[eCompNR][eIonNR];         /**< Data for a specific compartment and ion type.   */
    t_compartment     compsol[eCompNR];              /**< Solvent compartments.                           */
    int               fluxfromAtoB[eChanNR][eIonNR]; /**< Net flux per channels and ion type.             */
    int               ncyl0ions;                     /**< Number of ions residing in channel 0.           */
    int               ncyl1ions;                     /**< Same for channel 1.                             */
    int               cyl0and1;                      /**< Ions assigned to cyl0 and cyl1. Not good.       */
    int              *fluxleak;                      /**< Pointer to a single int value holding the
                                                          flux not going through any of the channels.     */
    real              deltaQ;                        /**< The charge imbalance between the compartments.  */
} t_swap;



/*! \brief Check whether point is in channel.
 *
 * A channel is a cylinder defined by a disc
 * with radius r around its center c. The thickness of the cylinder is
 * d_up - d_down.
 *
 * \code
 *               ^  d_up
 *               |
 *     r         |
 *     <---------c--------->
 *               |
 *               v  d_down
 *
 * \endcode
 */
static gmx_bool is_in_channel(
        rvec   point,  /* Point under consideration */
        rvec   center, /* 'Center' of cylinder */
        real   d_up,   /* Upper extension */
        real   d_down, /* Lower extensions */
        real   r_cyl2, /* Cylinder radius squared */
        t_pbc *pbc,
        int    normal) /* The membrane normal direction is typically 3, i.e. ZZ, but can be X or Y also */
{
    rvec dr;
    int  plane1, plane2; /* Directions tangential to membrane */


    plane1 = (normal + 1) % 3; /* typically 0, i.e. XX */
    plane2 = (normal + 2) % 3; /* typically 1, i.e. YY */

    /* Get the distance vector dr between the point and the center of the cylinder */
    pbc_dx(pbc, point, center, dr); /* This puts center in the origin */

    /* Check vertical direction */
    if ( (dr[normal] > d_up) || (dr[normal] < -d_down) )
    {
        return FALSE;
    }

    /* Check radial direction */
    if ( (dr[plane1]*dr[plane1] + dr[plane2]*dr[plane2]) > r_cyl2)
    {
        return FALSE;
    }

    /* All check passed, this point is in the cylinder */
    return TRUE;
}


/*! \brief Prints to swap output file which ions are in which compartment. */
static void print_ionlist(
        t_swap       *s,
        double        time,
        const char    comment[])
{
    int            itype, icomp, i, j;
    t_compartment *comp;


    fprintf(s->fpout, "%12.5e", time);
    for (icomp = 0; icomp < eCompNR; icomp++)
    {
        for (itype = 0; itype < eIonNR; itype++)
        {
            comp = &(s->comp[icomp][itype]);
            fprintf(s->fpout, "%7d%7.1f%7d", comp->nat, comp->nat_av-comp->nat_req, comp->inflow_netto);
        }
    }
    fprintf(s->fpout, "%12.3e%12.3e",
            s->group[eGrpSplit0].center[s->swapdim],
            s->group[eGrpSplit1].center[s->swapdim]);

    for (i = 0; i < eChanNR; i++)
    {
        for (j = 0; j < eIonNR; j++)
        {
            fprintf(s->fpout, "%12d", s->fluxfromAtoB[i][j]);
        }
    }

    /* Also print the number of ions that leaked from A to B: */
    fprintf(s->fpout, "%12d", *s->fluxleak);

    fprintf(s->fpout, "%s\n", comment);
}


/*! \brief Get the center of a group of nat atoms.
 *
 * Since with PBC an atom group might not be whole, use the first atom as the
 * reference atom and determine the center with respect to this reference.
 */
static void get_molecule_center(
        rvec   x[],
        int    nat,
        real  *weights,
        rvec   center,
        t_pbc *pbc)
{
    int  i;
    rvec weightedPBCimage;
    real wi, wsum;
    rvec reference, correctPBCimage, dx;


    /* Use the first atom as the reference and put other atoms near that one */
    /* This does not work for large molecules that span > half of the box! */
    copy_rvec(x[0], reference);

    /* Calculate either the weighted center or simply the center of geometry */
    wsum = 0.0;
    clear_rvec(center);
    for (i = 0; i < nat; i++)
    {
        /* PBC distance between position and reference */
        pbc_dx(pbc, x[i], reference, dx);

        /* Add PBC distance to reference */
        rvec_add(reference, dx, correctPBCimage);

        /* Take weight into account */
        if (NULL == weights)
        {
            wi = 1.0;
        }
        else
        {
            wi = weights[i];
        }
        wsum += wi;
        svmul(wi, correctPBCimage, weightedPBCimage);

        /* Add to center */
        rvec_inc(center, weightedPBCimage);
    }

    /* Normalize */
    svmul(1.0/wsum, center, center);
}



/*! \brief Return TRUE if ion is found in the compartment.
 *
 * Returns TRUE if x is between (w1+gap) and (w2-gap)
 *
 * \code
 *
 * ||-----------|--+--|----------o----------|--+--|---------------------||
 *                w1   ?????????????????????  w2
 *
 * \endcode
 */
static gmx_bool compartment_contains_atom(
        real  w1, /* position of wall atom 1 */
        real  w2, /* position of wall atom 2 */
        real  gap,
        real  x,
        real  l,  /* length of the box, from || to || in the sketch */
        real *distance_from_center)
{
    real m, l_2;


    /* First set the origin in the middle of w1 and w2 */
    m   = 0.5 * (w1 + w2);
    w1 -= m;
    w2 -= m;
    x  -= m;

    /* Now choose the PBC image of x that is closest to the origin: */
    l_2 = 0.5*l;
    while (x  > l_2)
    {
        x -= l;
    }
    while (x <= -l_2)
    {
        x += l;
    }

    *distance_from_center = (real)fabs(x);

    /* Return TRUE if we now are in area "????" */
    if ( (x >= (w1+gap)) &&  (x < (w2-gap)) )
    {
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}


/*! \brief Updates the time-averaged number of ions in a compartment. */
static void update_time_window(t_compartment *comp, int values, int replace)
{
    real average;
    int  i;


    /* Put in the new value */
    if (replace >= 0)
    {
        comp->nat_past[replace] = comp->nat;
    }

    /* Compute the new time-average */
    average = 0.0;
    for (i = 0; i < values; i++)
    {
        average += comp->nat_past[i];
    }
    average     /= values;
    comp->nat_av = average;
}


/*! \brief Add atom with collective index ci to the list 'comp'. */
static void add_to_list(
        int            ci,       /* index of this ion in the collective array xc, qc */
        t_compartment *comp,     /* Compartment to add this atom to */
        real           distance) /* Shortest distance of this atom to the compartment center */
{
    int nr;


    nr = comp->nat;

    if (nr >= comp->nalloc)
    {
        comp->nalloc = over_alloc_dd(nr+1);
        srenew(comp->ind, comp->nalloc);
        srenew(comp->dist, comp->nalloc);
    }
    comp->ind[nr]  = ci;
    comp->dist[nr] = distance;
    comp->nat++;
}


/*! \brief Determine the compartment boundaries from the channel centers. */
static void get_compartment_boundaries(
        int c,
        t_swap *s,
        matrix box,
        real *left, real *right)
{
    real pos0, pos1;
    real leftpos, rightpos, leftpos_orig;


    if (c >= eCompNR)
    {
        gmx_fatal(FARGS, "No compartment %d.", c);
    }

    pos0 = s->group[eGrpSplit0].center[s->swapdim];
    pos1 = s->group[eGrpSplit1].center[s->swapdim];

    if (pos0 < pos1)
    {
        leftpos  = pos0;
        rightpos = pos1;
    }
    else
    {
        leftpos  = pos1;
        rightpos = pos0;
    }

    /* This gets us the other compartment: */
    if (c == eCompB)
    {
        leftpos_orig = leftpos;
        leftpos      = rightpos;
        rightpos     = leftpos_orig + box[s->swapdim][s->swapdim];
    }

    *left  = leftpos;
    *right = rightpos;
}


/*! \brief Determine the per-channel ion flux.
 *
 * To determine the flux through the individual channels, we
 * remember the compartment and channel history of each ion. An ion can be
 * either in channel0 or channel1, or in the remaining volume of compartment
 * A or B.
 *
 * \code
 *    +-----------------+
 *    |                 | B
 *    |                 | B compartment
 *    ||||||||||0|||||||| bilayer with channel 0
 *    |                 | A
 *    |                 | A
 *    |                 | A compartment
 *    |                 | A
 *    |||||1||||||||||||| bilayer with channel 1
 *    |                 | B
 *    |                 | B compartment
 *    +-----------------+
 *
 * \endcode
 */
static void detect_flux_per_channel(
        int             iion,
        int             comp,
        int             iontype,
        rvec            ion_pos,
        unsigned char  *comp_now,
        unsigned char  *comp_from,
        unsigned char  *channel_label,
        t_swapcoords   *sc,
        real            cyl0_r2,
        real            cyl1_r2,
        gmx_int64_t     step,
        gmx_bool        bRerun,
        FILE           *fpout)
{
    gmx_swapcoords_t s;
    int              sd, chan_nr;
    gmx_bool         in_cyl0, in_cyl1;
    char             buf[STRLEN];


    s    = sc->si_priv;
    sd   = s->swapdim;

    /* Check whether ion is inside any of the channels */
    in_cyl0 = is_in_channel(ion_pos, s->group[eGrpSplit0].center, sc->cyl0u, sc->cyl0l, cyl0_r2, s->pbc, sd);
    in_cyl1 = is_in_channel(ion_pos, s->group[eGrpSplit1].center, sc->cyl1u, sc->cyl1l, cyl1_r2, s->pbc, sd);

    if (in_cyl0 && in_cyl1)
    {
        /* Ion appears to be in both channels. Something is severely wrong! */
        s->cyl0and1++;
        *comp_now      = eDomainNotset;
        *comp_from     = eDomainNotset;
        *channel_label = eChHistPassedNone;
    }
    else if (in_cyl0)
    {
        /* Ion is in channel 0 now */
        *channel_label = eChHistPassedCh0;
        *comp_now      = eDomainNotset;
        s->ncyl0ions++;
    }
    else if (in_cyl1)
    {
        /* Ion is in channel 1 now */
        *channel_label = eChHistPassedCh1;
        *comp_now      = eDomainNotset;
        s->ncyl1ions++;
    }
    else
    {
        /* Ion is not in any of the channels, so it must be in domain A or B */
        if (eCompA == comp)
        {
            *comp_now = eDomainA;
        }
        else
        {
            *comp_now = eDomainB;
        }
    }

    /* Only take action, if ion is now in domain A or B, and was before
     * in the other domain!
     */
    if (eDomainNotset == *comp_from)
    {
        /* Maybe we can set the domain now */
        *comp_from = *comp_now;               /* Could still be eDomainNotset, though */
    }
    else if (  (*comp_now  != eDomainNotset ) /* if in channel */
               && (*comp_from != *comp_now)  )
    {
        /* Obviously the ion changed its domain.
         * Count this for the channel through which it has passed. */
        switch (*channel_label)
        {
            case eChHistPassedNone:
                *s->fluxleak = *s->fluxleak + 1;

                fprintf(stderr, " %s Warning! Step %s, ion %d (%s) moved from %s to %s\n",
                        SwS, gmx_step_str(step, buf), iion, IonStr[iontype], DomainString[*comp_from], DomainString[*comp_now]);
                if (bRerun)
                {
                    fprintf(stderr, ", possibly due to a swap in the original simulation.\n");
                }
                else
                {
                    fprintf(stderr, "but did not pass cyl0 or cyl1 as defined in the .mdp file.\n"
                            "Do you have an ion somewhere within the membrane?\n");
                    /* Write this info to the CompEL output file: */
                    fprintf(s->fpout, " # Warning: step %s, ion %d (%s) moved from %s to %s (probably through the membrane)\n",
                            gmx_step_str(step, buf), iion, IonStr[iontype],
                            DomainString[*comp_from], DomainString[*comp_now]);

                }
                break;
            case eChHistPassedCh0:
            case eChHistPassedCh1:
                if (*channel_label == eChHistPassedCh0)
                {
                    chan_nr = 0;
                }
                else
                {
                    chan_nr = 1;
                }

                if (eDomainA == *comp_from)
                {
                    s->fluxfromAtoB[chan_nr][iontype]++;
                }
                else
                {
                    s->fluxfromAtoB[chan_nr][iontype]--;
                }
                fprintf(fpout, "# Atom nr. %d finished passing %s.\n", iion, ChannelString[*channel_label]);
                break;
            default:
                gmx_fatal(FARGS, "%s Unknown channel history entry!\n", SwS);
                break;
        }

        /* This ion has moved to the _other_ compartment ... */
        *comp_from = *comp_now;
        /* ... and it did not pass any channel yet */
        *channel_label = eChHistPassedNone;
    }
}


/*! \brief Get the lists of ions for the two compartments */
static void compartmentalize_ions(
        t_commrec      *cr,
        t_swapcoords   *sc,
        matrix          box,
        gmx_int64_t     step,
        FILE           *fpout,
        gmx_bool        bRerun)
{
    gmx_swapcoords_t s;
    int              i, sd, replace;
    real             left, right;
    t_group         *iong;
    real             dist;
    real             cyl0_r2, cyl1_r2;
    int              comp, type;
    int              sum, not_in_comp[eCompNR]; /* consistency check */
    int              ion_nr_global;


    s    = sc->si_priv;
    iong = &s->group[eGrpIons];
    sd   = s->swapdim;

    cyl0_r2 = sc->cyl0r * sc->cyl0r;
    cyl1_r2 = sc->cyl1r * sc->cyl1r;


    /* Get us a counter that cycles in the range of [0 ... sc->nAverage[ */
    replace = (step/sc->nstswap) % sc->nAverage;

    for (comp = eCompA; comp <= eCompB; comp++)
    {
        /* Get lists of atoms that match criteria for this compartment */
        get_compartment_boundaries(comp, sc->si_priv, box, &left, &right);

        /* First clear the ion lists */
        s->comp[comp][eIonNEG].nat = 0;
        s->comp[comp][eIonPOS].nat = 0;
        not_in_comp[comp]          = 0; /* consistency check */

        /* Loop over the IONS */
        for (i = 0; i < iong->nat; i++)
        {
            /* Anion or cation? */
            type = iong->qc[i] < 0 ? eIonNEG : eIonPOS;

            /* Is this ion in the compartment that we look at? */
            if (compartment_contains_atom(left, right, 0, iong->xc[i][sd], box[sd][sd], &dist) )
            {
                /* Now put it into the list containing only ions of its type */
                add_to_list(i, &s->comp[comp][type], dist);

                /* Master also checks through which channel each ion has passed */
                if (MASTER(cr) && (iong->comp_now != NULL))
                {
                    ion_nr_global = iong->ind[i] + 1; /* PDB index starts at 1 ... */
                    detect_flux_per_channel(ion_nr_global, comp, type, iong->xc[i],
                                            &iong->comp_now[i], &iong->comp_from[i], &iong->channel_label[i],
                                            sc, cyl0_r2, cyl1_r2, step, bRerun, fpout);
                }
            }
            else
            {
                not_in_comp[comp] += 1;
            }
        }
        /* Correct the time-averaged number of ions in both compartments */
        update_time_window(&s->comp[comp][eIonNEG], sc->nAverage, replace);
        update_time_window(&s->comp[comp][eIonPOS], sc->nAverage, replace);
    }

    /* Flux detection warnings */
    if (MASTER(cr) )
    {
        if (s->cyl0and1 > 0)
        {
            fprintf(stderr, "\n"
                    "%s Warning: %d atoms were detected as being in both channels! Probably your split\n"
                    "%s          cylinder is way too large, or one compartment has collapsed (step %" GMX_PRId64 ")\n",
                    SwS, s->cyl0and1, SwS, step);

            fprintf(s->fpout, "Warning: %d atoms were assigned to both channels!\n", s->cyl0and1);

            s->cyl0and1 = 0;
        }
    }


    /* Consistency checks */
    if (not_in_comp[eCompA] + not_in_comp[eCompB] != iong->nat)
    {
        if (NULL != fpout)
        {
            fprintf(fpout, "# Warning: Inconsistency during ion compartmentalization. !inA: %d, !inB: %d, total ions %d\n",
                    not_in_comp[eCompA], not_in_comp[eCompB], iong->nat);
            fflush(fpout);
        }
        else
        {
            fprintf(stderr, "%s rank %d: Inconsistency during ion compartmentalization. !inA: %d, !inB: %d, total ions %d\n",
                    SwS, cr->nodeid, not_in_comp[eCompA], not_in_comp[eCompB], iong->nat);

        }
    }
    sum =   s->comp[eCompA][eIonNEG].nat + s->comp[eCompA][eIonPOS].nat
        + s->comp[eCompB][eIonNEG].nat + s->comp[eCompB][eIonPOS].nat;
    if (sum != iong->nat)
    {
        if (NULL != fpout)
        {
            fprintf(fpout, "# Warning: %d atoms are in the ion group, but altogether %d have been assigned to the compartments.\n",
                    iong->nat, sum);
            fflush(fpout);
        }
        else
        {
            fprintf(stderr, "%s rank %d: %d atoms are in the ion group, but altogether %d have been assigned to the compartments.\n",
                    SwS, cr->nodeid, iong->nat, sum);
        }
    }


}


/*! \brief Set up the compartments and get lists of solvent atoms in each compartment */
static void compartmentalize_solvent(
        t_commrec    *cr,
        t_swapcoords *sc,
        matrix        box,
        FILE         *fpout)
{
    gmx_swapcoords_t s;
    int              apm, i, j, sd;
    real             left, right;
    t_group         *solg;
    real             dist;
    int              comp;
    int              sum, not_in_comp[eCompNR]; /* consistency check */


    s    = sc->si_priv;
    solg = &s->group[eGrpSolvent];
    sd   = s->swapdim;
    apm  = solg->apm;

    for (comp = eCompA; comp <= eCompB; comp++)
    {
        /* Get lists of atoms that match criteria for this compartment */
        get_compartment_boundaries(comp, sc->si_priv, box, &left, &right);

        /* First clear the solvent molecule lists */
        s->compsol[comp].nat = 0;
        not_in_comp[comp]    = 0; /* consistency check */

        /* Loop over the solvent MOLECULES */
        for (i = 0; i < sc->nat_sol; i += apm)
        {
            if (compartment_contains_atom(left, right, 0, solg->xc[i][sd], box[sd][sd], &dist))
            {
                /* Add the whole molecule to the list */
                for (j = 0; j < apm; j++)
                {
                    add_to_list(i+j, &s->compsol[comp], dist);
                }
            }
            else
            {
                not_in_comp[comp] += apm;
            }
        }
    }

    if (NULL != fpout)
    {
        fprintf(fpout, "# Solv. molecules in comp.%s: %d   comp.%s: %d\n",
                CompStr[eCompA], s->compsol[eCompA].nat/apm,
                CompStr[eCompB], s->compsol[eCompB].nat/apm);
    }

    /* Consistency checks */
    if (not_in_comp[eCompA] + not_in_comp[eCompB] != solg->nat)
    {
        if (NULL != fpout)
        {
            fprintf(fpout, "# Warning: Inconsistency during solvent compartmentalization. !inA: %d, !inB: %d, solvent atoms %d\n",
                    not_in_comp[eCompA], not_in_comp[eCompB], solg->nat);
            fflush(fpout);
        }
        else
        {
            fprintf(stderr, "%s rank %d: Inconsistency during solvent compartmentalization. !inA: %d, !inB: %d, solvent atoms %d\n",
                    SwS, cr->nodeid, not_in_comp[eCompA], not_in_comp[eCompB], solg->nat);
        }
    }
    sum = s->compsol[eCompA].nat + s->compsol[eCompB].nat;
    if (sum != solg->nat)
    {
        if (NULL != fpout)
        {
            fprintf(fpout, "# Warning: %d atoms in solvent group, but %d have been assigned to the compartments.\n",
                    solg->nat, sum);
            fflush(fpout);
        }
        else
        {
            fprintf(stderr, "%s rank %d: %d atoms in solvent group, but %d have been assigned to the compartments.\n",
                    SwS, cr->nodeid, solg->nat, sum);
        }
    }
}


/*! \brief Find out how many group atoms are in the compartments initially */
static void get_initial_ioncounts(
        t_inputrec       *ir,
        rvec              x[],   /* the initial positions */
        matrix            box,
        t_commrec        *cr,
        gmx_bool          bRerun)
{
    t_swapcoords *sc;
    t_swap       *s;
    int           i, ii, ind, ic;
    int           req[2], tot[2];


    sc = ir->swap;
    s  = sc->si_priv;

    /* Copy the initial swap group positions to the collective array so
     * that we can compartmentalize */
    for (i = 0; i < sc->nat; i++)
    {
        ind = sc->ind[i];
        copy_rvec(x[ind], s->group[eGrpIons].xc[i]);
    }

    /* Set up the compartments and get lists of atoms in each compartment */
    compartmentalize_ions(cr, sc, box, 0, s->fpout, bRerun);

    /* Set initial concentrations if requested */
    for (ic = 0; ic < eCompNR; ic++)
    {
        s->comp[ic][eIonPOS].nat_req = sc->ncations[ic];
        s->comp[ic][eIonNEG].nat_req = sc->nanions[ic];
    }
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            if (s->comp[ic][ii].nat_req < 0)
            {
                s->comp[ic][ii].nat_req = s->comp[ic][ii].nat;
            }
        }
    }

    /* Check whether the number of requested ions adds up to the total number of ions */
    for (ii = 0; ii < eIonNR; ii++)
    {
        req[ii] = s->comp[eCompA][ii].nat_req + s->comp[eCompB][ii].nat_req;
        tot[ii] = s->comp[eCompA][ii].nat     + s->comp[eCompB][ii].nat;
    }
    if ( (req[eCompA] != tot[eCompA]) || (req[eCompB] != tot[eCompB ]) )
    {
        gmx_fatal(FARGS, "Mismatch of the number of ions summed over both compartments.\n"
                  "You requested a total of %d anions and %d cations,\n"
                  "but there are a total of %d anions and %d cations in the system.\n",
                  req[eIonNEG], req[eIonPOS],
                  tot[eIonNEG], tot[eIonPOS]);
    }

    /* Initialize time-averaging:
     * Write initial concentrations to all time bins to start with */
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            s->comp[ic][ii].nat_av = s->comp[ic][ii].nat;
            for (i = 0; i < sc->nAverage; i++)
            {
                s->comp[ic][ii].nat_past[i] = s->comp[ic][ii].nat;
            }
        }
    }
}


/*! \brief Copy history of ion counts from checkpoint file.
 *
 * When called, the checkpoint file has already been read in. Here we copy
 * over the values from .cpt file to the swap data structure.
 */
static void get_initial_ioncounts_from_cpt(
        t_inputrec *ir, swapstate_t *swapstate,
        t_commrec *cr, gmx_bool bVerbose)
{
    t_swapcoords *sc;
    t_swap       *s;
    int           ic, ii, j;

    sc = ir->swap;
    s  = sc->si_priv;

    if (MASTER(cr))
    {
        /* Copy the past values from the checkpoint values that have been read in already */
        if (bVerbose)
        {
            fprintf(stderr, "%s Copying values from checkpoint\n", SwS);
        }

        for (ic = 0; ic < eCompNR; ic++)
        {
            for (ii = 0; ii < eIonNR; ii++)
            {
                s->comp[ic][ii].nat_req      = swapstate->nat_req[ic][ii];
                s->comp[ic][ii].inflow_netto = swapstate->inflow_netto[ic][ii];

                if (bVerbose)
                {
                    fprintf(stderr, "%s ... Influx netto: %d   Requested: %d   Past values: ", SwS,
                            s->comp[ic][ii].inflow_netto, s->comp[ic][ii].nat_req);
                }

                for (j = 0; j < sc->nAverage; j++)
                {
                    s->comp[ic][ii].nat_past[j] = swapstate->nat_past[ic][ii][j];
                    if (bVerbose)
                    {
                        fprintf(stderr, "%d ", s->comp[ic][ii].nat_past[j]);
                    }
                }
                if (bVerbose)
                {
                    fprintf(stderr, "\n");
                }
            }
        }
    }
}


/*! \brief The master lets all others know about the initial ion counts. */
static void bc_initial_concentrations(
        t_commrec    *cr,
        t_swapcoords *swap)
{
    int     ic, ii;
    t_swap *s;

    s = swap->si_priv;

    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            gmx_bcast(sizeof(s->comp[ic][ii].nat_req), &(s->comp[ic][ii].nat_req), cr);
            gmx_bcast(sizeof(s->comp[ic][ii].nat    ), &(s->comp[ic][ii].nat    ), cr);
            gmx_bcast( swap->nAverage * sizeof(s->comp[ic][ii].nat_past[0]), s->comp[ic][ii].nat_past, cr);
        }
    }
}


/*! \brief Ensure that each atom belongs to at most one of the swap groups. */
static void check_swap_groups(t_swap *s, int nat, gmx_bool bVerbose)
{
    t_group  *g;
    int       i, j;
    atom_id  *nGroup    = NULL; /* This array counts for each atom in the MD system to
                                   how many swap groups it belongs (should be 0 or 1!) */
    int       ind       = -1;
    int       nMultiple = 0;    /* Number of atoms belonging to multiple groups */


    if (bVerbose)
    {
        fprintf(stderr, "%s Making sure each atom belongs to at most one of the swap groups.\n", SwS);
    }

    /* Add one to the group count of atoms belonging to a swap group: */
    snew(nGroup, nat);
    for (i = 0; i < eGrpNr; i++)
    {
        g = &s->group[i];
        for (j = 0; j < g->nat; j++)
        {
            /* Get the global index of this atom of this group: */
            ind = g->ind[j];
            nGroup[ind]++;
        }
    }
    /* Make sure each atom belongs to at most one swap group: */
    for (j = 0; j < g->nat; j++)
    {
        if (nGroup[j] > 1)
        {
            nMultiple++;
        }
    }
    sfree(nGroup);

    if (nMultiple)
    {
        gmx_fatal(FARGS, "%s Cannot perform swapping since %d atom%s allocated to more than one swap index group.\n"
                  "%s Each atom must be allocated to at most one of the split groups, the swap group, or the solvent.\n"
                  "%s Check the .mdp file settings regarding the swap index groups or the index groups themselves.\n",
                  SwS, nMultiple, (1 == nMultiple) ? " is" : "s are", SwSEmpty, SwSEmpty);
    }
}


/*! \brief Get the number of atoms per molecule for this group.
 *
 * Also ensure that all the molecules in this group have this number of atoms.
 */
static int get_group_apm_check(
        int                         group,
        t_swap                     *s,
        gmx_bool                    bVerbose,
        const gmx_mtop_atomlookup_t alook,
        gmx_mtop_t                 *mtop)
{
    int *ind;
    int  nat, apm, i;
    int  molb, molnr, atnr_mol;


    ind = s->group[group].ind;
    nat = s->group[group].nat;

    /* Determine the number of solvent atoms per solvent molecule from the
     * first solvent atom: */
    i = 0;
    gmx_mtop_atomnr_to_molblock_ind(alook, ind[i], &molb, &molnr, &atnr_mol);
    apm = mtop->molblock[molb].natoms_mol;

    if (bVerbose)
    {
        fprintf(stderr, "%s Checking whether all %s molecules consist of %d atom%s\n",
                SwS, GrpString[group], apm, apm > 1 ? "s" : "");
    }

    /* Check whether this is also true for all other solvent atoms */
    for (i = 1; i < nat; i++)
    {
        gmx_mtop_atomnr_to_molblock_ind(alook, ind[i], &molb, &molnr, &atnr_mol);
        if (apm != mtop->molblock[molb].natoms_mol)
        {
            gmx_fatal(FARGS, "Not all %s group molecules consist of %d atoms.",
                      GrpString[group], apm);
        }
    }

    return apm;
}


/*! \brief Print the legend to the swap output file.
 *
 * Also print the initial ion counts
 */
static void print_ionlist_legend(t_inputrec *ir, const output_env_t oenv)
{
    const char **legend;
    int          ic, count, ii;
    char         buf[256];
    t_swap      *s;


    s = ir->swap->si_priv;

    snew(legend, eCompNR*eIonNR*3 + 2 + eChanNR*eIonNR + 1);
    for (ic = count = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            sprintf(buf, "%s %ss", CompStr[ic], IonString[ii]);
            legend[count++] = gmx_strdup(buf);
            sprintf(buf, "%s av. mismatch to %d%s",
                    CompStr[ic], s->comp[ic][ii].nat_req, IonStr[ii]);
            legend[count++] = gmx_strdup(buf);
            sprintf(buf, "%s netto %s influx", CompStr[ic], IonString[ii]);
            legend[count++] = gmx_strdup(buf);
        }
    }
    sprintf(buf, "%scenter of %s of split group 0", SwapStr[ir->eSwapCoords], (NULL != s->group[eGrpSplit0].m) ? "mass" : "geometry");
    legend[count++] = gmx_strdup(buf);
    sprintf(buf, "%scenter of %s of split group 1", SwapStr[ir->eSwapCoords], (NULL != s->group[eGrpSplit1].m) ? "mass" : "geometry");
    legend[count++] = gmx_strdup(buf);

    for (ic = 0; ic < eChanNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            sprintf(buf, "A->ch%d->B %s permeations", ic, IonString[ii]);
            legend[count++] = gmx_strdup(buf);
        }
    }

    sprintf(buf, "leakage");
    legend[count++] = gmx_strdup(buf);

    xvgr_legend(s->fpout, count, legend, oenv);

    fprintf(s->fpout, "# Instantaneous ion counts and time-averaged differences to requested numbers\n");
    fprintf(s->fpout, "#   time[ps]   A_an   diff   t_in  A_cat   diff   t_in   B_an   diff   t_in  B_cat   diff   t_in ");
    fprintf(s->fpout, "   %s-Split0    %s-Split1", DimStr[s->swapdim], DimStr[s->swapdim]);
    fprintf(s->fpout, "  A-ch0-B_an A-ch0-B_cat  A-ch1-B_an A-ch1-B_cat ion_leakage\n");
    fflush(s->fpout);
}


/*! \brief Initialize channel ion flux detection routine.
 *
 * Initialize arrays that keep track of where the ions come from and where
 * they go.
 */
static void detect_flux_per_channel_init(
        t_commrec   *cr,
        t_swap      *s,
        swapstate_t *swapstate,
        gmx_bool     bStartFromCpt)
{
    int      i, ic, ii;
    t_group *g;


    g = &(s->group[eGrpIons]);

    /* All these flux detection routines run on the master only */
    if (!MASTER(cr))
    {
        g->comp_now      = NULL;
        g->comp_from     = NULL;
        g->channel_label = NULL;

        return;
    }

    /******************************************************/
    /* Channel and domain history for the individual ions */
    /******************************************************/
    if (bStartFromCpt) /* set the pointers right */
    {
        g->comp_from     = swapstate->comp_from;
        g->channel_label = swapstate->channel_label;
    }
    else /* allocate memory */
    {
        snew(g->comp_from, g->nat);
        swapstate->comp_from = g->comp_from;
        snew(g->channel_label, g->nat);
        swapstate->channel_label = g->channel_label;
    }
    snew(g->comp_now, g->nat);

    /* Initialize the channel and domain history counters */
    for (i = 0; i < g->nat; i++)
    {
        g->comp_now[i] = eDomainNotset;
        if (!bStartFromCpt)
        {
            g->comp_from[i]     = eDomainNotset;
            g->channel_label[i] = eChHistPassedNone;
        }
    }

    /************************************/
    /* Channel fluxes for both channels */
    /************************************/
    s->ncyl0ions = 0;
    s->ncyl1ions = 0;
    s->cyl0and1  = 0;

    if (bStartFromCpt)
    {
        fprintf(stderr, "%s Copying channel fluxes from checkpoint file data\n", SwS);
    }

    for (ic = 0; ic < eChanNR; ic++)
    {
        fprintf(stderr, "%s Channel %d flux history: ", SwS, ic);
        for (ii = 0; ii < eIonNR; ii++)
        {
            if (bStartFromCpt)
            {
                s->fluxfromAtoB[ic][ii] = swapstate->fluxfromAtoB[ic][ii];
            }
            else
            {
                s->fluxfromAtoB[ic][ii] = 0;
            }

            fprintf(stderr, "%d %s%s   ", s->fluxfromAtoB[ic][ii], IonString[ii], s->fluxfromAtoB[ic][ii] == 1 ? "" : "s");
        }
        fprintf(stderr, "\n");
    }
    if (bStartFromCpt)
    {
        s->fluxleak = swapstate->fluxleak;
    }
    else
    {
        snew(s->fluxleak, 1);
        *s->fluxleak = 0;
        /* Set pointer for checkpoint writing */
        swapstate->fluxleak = s->fluxleak;
    }

    /* Set pointers for checkpoint writing */
    for (ic = 0; ic < eChanNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            swapstate->fluxfromAtoB_p[ic][ii] = &(s->fluxfromAtoB[ic][ii]);
        }
    }
}


/*! \brief Outputs the initial structure to PDB file for debugging reasons.
 *
 * Output the starting structure so that in case of multimeric channels
 * the user can check whether we have the correct PBC image for all atoms.
 * If this is not correct, the ion counts per channel will be very likely
 * wrong.
 */
static void outputStartStructureIfWanted(gmx_mtop_t *mtop, rvec *x, int ePBC, matrix box)
{
    char *env = getenv("GMX_COMPELDUMP");

    if (env != NULL)
    {
        fprintf(stderr, "\n%s Found env.var. GMX_COMPELDUMP, will output CompEL starting structure made whole.\n"
                "%s In case of multimeric channels, please check whether they have the correct PBC representation.\n",
                SwS, SwSEmpty);

        write_sto_conf_mtop("CompELAssumedWholeConfiguration.pdb", *mtop->name, mtop, x, NULL, ePBC, box);
    }
}


/*! \brief Initialize the swapstate structure, used for checkpoint writing.
 *
 * The swapstate struct stores the information we need to make the channels
 * whole again after restarts from a checkpoint file. Here we do the following:\n
 * a) If we did not start from .cpt, we prepare the struct for proper .cpt writing,\n
 * b) if we did start from .cpt, we copy over the last whole structures from .cpt,\n
 * c) in any case, for subsequent checkpoint writing, we set the pointers in\n
 * swapstate to the x_old arrays, which contain the correct PBC representation of
 * multimeric channels at the last time step.
 */
static void init_swapstate(
        swapstate_t      *swapstate,
        t_swapcoords     *sc,
        gmx_mtop_t       *mtop,
        rvec              x[],      /* the initial positions */
        matrix            box,
        int               ePBC)
{
    int                    i, ig;
    rvec                  *x_pbc  = NULL;   /* positions of the whole MD system with molecules made whole */
    t_group               *g;
    t_swap                *s;


    s = sc->si_priv;

    /* We always need the last whole positions such that
     * in the next time step we can make the channels whole again in PBC */
    if (swapstate->bFromCpt)
    {
        /* Copy the last whole positions of each channel from .cpt */
        g = &(s->group[eGrpSplit0]);
        for (i = 0; i <  g->nat; i++)
        {
            copy_rvec(swapstate->xc_old_whole[eChan0][i], g->xc_old[i]);
        }
        g = &(s->group[eGrpSplit1]);
        for (i = 0; i <  g->nat; i++)
        {
            copy_rvec(swapstate->xc_old_whole[eChan1][i], g->xc_old[i]);
        }
    }
    else
    {
        /* Extract the initial split group positions. */

        /* Remove pbc, make molecule whole. */
        snew(x_pbc, mtop->natoms);
        m_rveccopy(mtop->natoms, x, x_pbc);

        /* This can only make individual molecules whole, not multimers */
        do_pbc_mtop(NULL, ePBC, box, mtop, x_pbc);

        /* Output the starting structure? */
        outputStartStructureIfWanted(mtop, x_pbc, ePBC, box);

        /* If this is the first run (i.e. no checkpoint present) we assume
         * that the starting positions give us the correct PBC representation */
        for (ig = eGrpSplit0; ig <= eGrpSplit1; ig++)
        {
            g = &(s->group[ig]);
            for (i = 0; i < g->nat; i++)
            {
                copy_rvec(x_pbc[g->ind[i]], g->xc_old[i]);
            }
        }
        sfree(x_pbc);

        /* Prepare swapstate arrays for later checkpoint writing */
        swapstate->nat[eChan0] = s->group[eGrpSplit0].nat;
        swapstate->nat[eChan1] = s->group[eGrpSplit1].nat;
    }

    /* For subsequent checkpoint writing, set the swapstate pointers to the xc_old
     * arrays that get updated at every swapping step */
    swapstate->xc_old_whole_p[eChan0] = &s->group[eGrpSplit0].xc_old;
    swapstate->xc_old_whole_p[eChan1] = &s->group[eGrpSplit1].xc_old;
}


extern void init_swapcoords(
        FILE              *fplog,
        gmx_bool           bVerbose,
        t_inputrec        *ir,
        const char        *fn,
        gmx_mtop_t        *mtop,
        rvec               x[],
        matrix             box,
        swapstate_t       *swapstate,
        t_commrec         *cr,
        const output_env_t oenv,
        unsigned long      Flags)
{
    int                    i, ic, ig, ii, j;
    t_swapcoords          *sc;
    t_swap                *s;
    t_atom                *atom;
    t_group               *g;
    gmx_bool               bAppend, bStartFromCpt, bRerun;
    gmx_mtop_atomlookup_t  alook = NULL;


    alook = gmx_mtop_atomlookup_init(mtop);

    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
    {
        gmx_fatal(FARGS, "Position swapping is only implemented for domain decomposition!");
    }

    bAppend       = Flags & MD_APPENDFILES;
    bStartFromCpt = Flags & MD_STARTFROMCPT;
    bRerun        = Flags & MD_RERUN;

    sc = ir->swap;
    snew(sc->si_priv, 1);
    s = sc->si_priv;

    if (bRerun)
    {
        if (PAR(cr))
        {
            gmx_fatal(FARGS, "%s This module does not support reruns in parallel\nPlease request a serial run with -nt 1 / -np 1\n", SwS);
        }

        fprintf(stderr, "%s Rerun - using every available frame\n", SwS);
        sc->nstswap  = 1;
        sc->nAverage = 1;  /* averaging makes no sense for reruns */
    }

    if (MASTER(cr) && !bAppend)
    {
        fprintf(fplog, "\nInitializing ion/water position exchanges\n");
        please_cite(fplog, "Kutzner2011b");
    }

    switch (ir->eSwapCoords)
    {
        case eswapX:
            s->swapdim = XX;
            break;
        case eswapY:
            s->swapdim = YY;
            break;
        case eswapZ:
            s->swapdim = ZZ;
            break;
        default:
            s->swapdim = -1;
            break;
    }

    /* Copy some data to the group structures for convenience */
    /* Number of atoms in the group */
    s->group[eGrpIons   ].nat = sc->nat;
    s->group[eGrpSplit0 ].nat = sc->nat_split[0];
    s->group[eGrpSplit1 ].nat = sc->nat_split[1];
    s->group[eGrpSolvent].nat = sc->nat_sol;
    /* Pointer to the indices */
    s->group[eGrpIons   ].ind = sc->ind;
    s->group[eGrpSplit0 ].ind = sc->ind_split[0];
    s->group[eGrpSplit1 ].ind = sc->ind_split[1];
    s->group[eGrpSolvent].ind = sc->ind_sol;

    check_swap_groups(s, mtop->natoms, bVerbose && MASTER(cr));

    /* Allocate space for the collective arrays for all groups */
    for (ig = 0; ig < eGrpNr; ig++)
    {
        g = &(s->group[ig]);
        snew(g->xc, g->nat);
        snew(g->c_ind_loc, g->nat);
        /* For the split groups (the channels) we need some extra memory to
         * be able to make the molecules whole even if they span more than
         * half of the box size. */
        if (eGrpSplit0 == ig || eGrpSplit1 == ig)
        {
            snew(g->xc_shifts, g->nat);
            snew(g->xc_eshifts, g->nat);
            snew(g->xc_old, g->nat);
        }
    }

    if (MASTER(cr))
    {
        init_swapstate(swapstate, sc, mtop, x, box, ir->ePBC);
    }

    /* After init_swapstate we have a set of (old) whole positions for our
     * channels. Now transfer that to all nodes */
    if (PAR(cr))
    {
        for (ig = eGrpSplit0; ig <= eGrpSplit1; ig++)
        {
            g = &(s->group[ig]);
            gmx_bcast((g->nat)*sizeof((g->xc_old)[0]), g->xc_old, (cr));
        }
    }

    /* Make sure that all molecules in the ion and solvent groups contain the
     * same number of atoms each */
    s->group[eGrpIons   ].apm = get_group_apm_check(eGrpIons, s, MASTER(cr) && bVerbose, alook, mtop);
    s->group[eGrpSolvent].apm = get_group_apm_check(eGrpSolvent, s, MASTER(cr) && bVerbose, alook, mtop);

    /* Save masses where needed */
    s->group[eGrpIons   ].m = NULL;
    /* We only need enough space to determine a single solvent molecule's
     * center at at time */
    g = &(s->group[eGrpSolvent]);
    snew(g->m, g->apm);

    /* Need mass-weighted center of split group? */
    for (j = 0, ig = eGrpSplit0; j < eChanNR; ig++, j++)
    {
        g = &(s->group[ig]);
        if (TRUE == sc->massw_split[j])
        {
            /* Save the split group charges if mass-weighting is requested */
            snew(g->m, g->nat);
            for (i = 0; i < g->nat; i++)
            {
                gmx_mtop_atomnr_to_atom(alook, g->ind[i], &atom);
                g->m[i] = atom->m;
            }
        }
        else
        {
            g->m = NULL;
        }
    }

    /* Save the ionic charges */
    g = &(s->group[eGrpIons]);
    snew(g->qc, g->nat);
    for (i = 0; i < g->nat; i++)
    {
        gmx_mtop_atomnr_to_atom(alook, g->ind[i], &atom);
        g->qc[i] = atom->q;
    }

    snew(s->pbc, 1);
    set_pbc(s->pbc, -1, box);


    if (MASTER(cr))
    {
        if (bVerbose)
        {
            fprintf(stderr, "%s Opening output file %s%s\n", SwS, fn, bAppend ? " for appending" : "");
        }

        s->fpout = gmx_fio_fopen(fn, bAppend ? "a" : "w" );

        if (!bAppend)
        {
            xvgr_header(s->fpout, "Ion counts", "Time (ps)", "counts", exvggtXNY, oenv);

            for (ig = 0; ig < eGrpNr; ig++)
            {
                g = &(s->group[ig]);
                fprintf(s->fpout, "# %s group contains %d atom%s", GrpString[ig], g->nat, (g->nat > 1) ? "s" : "");
                if (eGrpSolvent == ig || eGrpIons == ig)
                {
                    fprintf(s->fpout, " with %d atom%s in each molecule", g->apm, (g->apm > 1) ? "s" : "");
                }
                fprintf(s->fpout, ".\n");
            }

            fprintf(s->fpout, "#\n# Initial positions of split groups:\n");
        }

        for (j = 0, ig = eGrpSplit0; j < eChanNR; j++, ig++)
        {
            g = &(s->group[ig]);
            for (i = 0; i < g->nat; i++)
            {
                copy_rvec(x[sc->ind_split[j][i]], g->xc[i]);
            }
            if (eGrpSplit0 == ig || eGrpSplit1 == ig)
            {
                /* xc has the correct PBC representation for the two channels, so we do
                 * not need to correct for that */
                get_center(g->xc, g->m, g->nat, g->center);
            }
            else
            {
                /* For the water molecules, we need to make the molecules whole */
                get_molecule_center(g->xc, g->nat, g->m, g->center, s->pbc);
            }
            if (!bAppend)
            {
                fprintf(s->fpout, "# %s group %s-center %5f nm\n", GrpString[ig],
                        DimStr[s->swapdim], g->center[s->swapdim]);
            }
        }

        if (!bAppend)
        {
            fprintf(s->fpout, "#\n");
            fprintf(s->fpout, "# split0 cylinder radius %f nm, up %f nm, down %f nm\n",
                    sc->cyl0r, sc->cyl0u, sc->cyl0l);
            fprintf(s->fpout, "# split1 cylinder radius %f nm, up %f nm, down %f nm\n",
                    sc->cyl1r, sc->cyl1u, sc->cyl1l);
        }

        if (!bAppend)
        {
            fprintf(s->fpout, "#\n");
            if (!bRerun)
            {
                fprintf(s->fpout, "# Coupling constant (number of swap attempt steps to average over): %d  (translates to %f ps).\n",
                        sc->nAverage, sc->nAverage*sc->nstswap*ir->delta_t);
                fprintf(s->fpout, "# Threshold is %f\n", sc->threshold);
                fprintf(s->fpout, "#\n");
                fprintf(s->fpout, "# Remarks about which atoms passed which channel use global atoms numbers starting at one.\n");
            }
        }
    }
    else
    {
        s->fpout = NULL;
    }

    /* Prepare for parallel or serial run */
    if (PAR(cr))
    {
        for (ig = 0; ig < eGrpNr; ig++)
        {
            g             = &(s->group[ig]);
            g->nat_loc    = 0;
            g->nalloc_loc = 0;
            g->ind_loc    = NULL;
        }
    }
    else
    {
        for (ig = 0; ig < eGrpNr; ig++)
        {
            g          = &(s->group[ig]);
            g->nat_loc = g->nat;
            g->ind_loc = g->ind;
            /* c_ind_loc needs to be set to identity in the serial case */
            for (i = 0; i < g->nat; i++)
            {
                g->c_ind_loc[i] = i;
            }
        }
    }

    /* Allocate memory for the ion counts time window */
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            snew(s->comp[ic][ii].nat_past, sc->nAverage);
        }
    }

    /* Get the initial ion concentrations and let the other nodes know */
    if (MASTER(cr))
    {
        swapstate->nions = s->group[eGrpIons].nat;

        if (bStartFromCpt)
        {
            get_initial_ioncounts_from_cpt(ir, swapstate, cr, bVerbose);
        }
        else
        {
            fprintf(stderr, "%s Determining initial ion counts.\n", SwS);
            get_initial_ioncounts(ir, x, box, cr, bRerun);
        }

        /* Prepare (further) checkpoint writes ... */
        if (bStartFromCpt)
        {
            /* Consistency check */
            if (swapstate->nAverage != sc->nAverage)
            {
                gmx_fatal(FARGS, "%s Ion count averaging steps mismatch! checkpoint: %d, tpr: %d",
                          SwS, swapstate->nAverage, sc->nAverage);
            }
        }
        else
        {
            swapstate->nAverage = sc->nAverage;
        }
        fprintf(stderr, "%s Setting pointers for checkpoint writing\n", SwS);
        for (ic = 0; ic < eCompNR; ic++)
        {
            for (ii = 0; ii < eIonNR; ii++)
            {
                swapstate->nat_req_p[ic][ii]      = &(s->comp[ic][ii].nat_req);
                swapstate->nat_past_p[ic][ii]     = &(s->comp[ic][ii].nat_past[0]);
                swapstate->inflow_netto_p[ic][ii] = &(s->comp[ic][ii].inflow_netto);
            }
        }

        /* Determine the total charge imbalance */
        s->deltaQ =  ( (-1) * s->comp[eCompA][eIonNEG].nat_req + s->comp[eCompA][eIonPOS].nat_req )
            - ( (-1) * s->comp[eCompB][eIonNEG].nat_req + s->comp[eCompB][eIonPOS].nat_req );

        if (bVerbose)
        {
            fprintf(stderr, "%s Requested charge imbalance is Q(A) - Q(B) = %gz.\n", SwS, s->deltaQ);
        }
        if (!bAppend)
        {
            fprintf(s->fpout, "# Requested charge imbalance is Q(A)-Q(B) = %gz.\n", s->deltaQ);
        }
    }

    if (PAR(cr))
    {
        bc_initial_concentrations(cr, ir->swap);
    }

    /* Put the time-averaged number of ions for all compartments */
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            update_time_window(&(s->comp[ic][ii]), sc->nAverage, -1);
        }
    }

    /* Initialize arrays that keep track of through which channel the ions go */
    detect_flux_per_channel_init(cr, s, swapstate, bStartFromCpt);

    /* We need to print the legend if we open this file for the first time. */
    if (MASTER(cr) && !bAppend)
    {
        print_ionlist_legend(ir, oenv);
    }
}


extern void dd_make_local_swap_groups(gmx_domdec_t *dd, t_swapcoords *sc)
{
    t_group *g;
    int      ig;


    /* Make ion group, split groups and solvent group */
    for (ig = 0; ig < eGrpNr; ig++)
    {
        g = &(sc->si_priv->group[ig]);
        dd_make_local_group_indices(dd->ga2la, g->nat, g->ind,
                                    &(g->nat_loc), &(g->ind_loc), &(g->nalloc_loc), g->c_ind_loc);
    }
}


/*! \brief Do we need to swap ions with water molecules at this step?
 *
 * From the requested and average ion counts we determine whether a swap is needed
 * at this time step.
 */
static gmx_bool need_swap(t_swapcoords *sc)
{
    t_swap *s;
    int     ic, ii;


    s = sc->si_priv;
    for (ic = 0; ic < eCompNR; ic++)
    {
        for (ii = 0; ii < eIonNR; ii++)
        {
            if (s->comp[ic][ii].nat_req - s->comp[ic][ii].nat_av >= sc->threshold)
            {
                return TRUE;
            }
        }
    }
    return FALSE;
}


/*! \brief Return index of atom that we can use for swapping.
 *
 * Returns the index of an atom that is far off the compartment boundaries.
 * Other atoms of the molecule (if any) will directly follow the returned index
 */
static int get_index_of_distant_atom(
        t_compartment *comp,
        int            apm) /* Atoms per molecule - just return the first atom index of a molecule */
{
    int  i, ibest = -1;
    real d = GMX_REAL_MAX;


    /* comp->nat contains the original number of atoms in this compartment
     * prior to doing any swaps. Some of these atoms may already have been
     * swapped out, but then they are marked with a distance of GMX_REAL_MAX
     */
    for (i = 0; i < comp->nat_old; i += apm)
    {
        if (comp->dist[i] < d)
        {
            ibest = i;
            d     = comp->dist[ibest];
        }
    }

    if (ibest < 0)
    {
        gmx_fatal(FARGS, "Could not get index of swap atom. Compartment atoms %d before swaps, atoms per molecule %d.",
                  comp->nat_old, apm);
    }

    /* Set the distance of this index to infinity such that it won't get selected again in
     * this time step
     */
    comp->dist[ibest] = GMX_REAL_MAX;

    return comp->ind[ibest];
}


/*! \brief Swaps centers of mass and makes molecule whole if broken */
static void translate_positions(
        rvec  *x,
        int    apm,
        rvec   old_com,
        rvec   new_com,
        t_pbc *pbc)
{
    int  i;
    rvec reference, dx, correctPBCimage;


    /* Use the first atom as the reference for PBC */
    copy_rvec(x[0], reference);

    for (i = 0; i < apm; i++)
    {
        /* PBC distance between position and reference */
        pbc_dx(pbc, x[i], reference, dx);

        /* Add PBC distance to reference */
        rvec_add(reference, dx, correctPBCimage);

        /* Subtract old_com from correct image and add new_com */
        rvec_dec(correctPBCimage, old_com);
        rvec_inc(correctPBCimage, new_com);

        copy_rvec(correctPBCimage, x[i]);
    }
}


/*! \brief Write back the the modified local positions from the collective array to the official positions. */
static void apply_modified_positions(
        t_group *g,
        rvec     x[])
{
    int l, ii, cind;


    for (l = 0; l < g->nat_loc; l++)
    {
        /* Get the right local index to write to */
        ii = g->ind_loc[l];
        /* Where is the local atom in the collective array? */
        cind = g->c_ind_loc[l];

        /* Copy the possibly modified position */
        copy_rvec(g->xc[cind], x[ii]);
    }
}


extern gmx_bool do_swapcoords(
        t_commrec        *cr,
        gmx_int64_t       step,
        double            t,
        t_inputrec       *ir,
        gmx_wallcycle_t   wcycle,
        rvec              x[],
        matrix            box,
        gmx_mtop_t       *mtop,
        gmx_bool          bVerbose,
        gmx_bool          bRerun)
{
    t_swapcoords         *sc;
    t_swap               *s;
    int                   j, ii, ic, ig, im, gmax, nswaps;
    gmx_bool              bSwap = FALSE;
    t_group              *g;
    real                  vacancy[eCompNR][eIonNR];
    int                   isol, iion;
    rvec                  solvent_center, ion_center;
    t_atom               *atom;
    gmx_mtop_atomlookup_t alook = NULL;


    wallcycle_start(wcycle, ewcSWAP);

    sc  = ir->swap;
    s   = sc->si_priv;

    /* Assemble all the positions of the swap group (ig = 0), the split groups
     * (ig = 1,2), and possibly the solvent group (ig = 3) */
    gmax = eGrpNr;

    for (ig = 0; ig < gmax; ig++)
    {
        g = &(s->group[ig]);
        if (eGrpSplit0 == ig || eGrpSplit1 == ig)
        {
            /* The split groups, i.e. the channels. Here we need  the full
             * communicate_group_positions(), so that we can make the molecules
             * whole even in cases where they span more than half of the box in
             * any dimension */
            communicate_group_positions(cr, g->xc, g->xc_shifts, g->xc_eshifts, TRUE,
                                        x, g->nat, g->nat_loc, g->ind_loc, g->c_ind_loc, g->xc_old, box);

            get_center(g->xc, g->m, g->nat, g->center); /* center of split groups == channels */
        }
        else
        {
            /* Swap group (ions), and solvent group. These molecules are small
             * and we can always make them whole with a simple distance check.
             * Therefore we pass NULL as third argument. */
            communicate_group_positions(cr, g->xc, NULL, NULL, FALSE,
                                        x, g->nat, g->nat_loc, g->ind_loc, g->c_ind_loc, NULL, NULL);
        }
    }

    /* Set up the compartments and get lists of atoms in each compartment,
     * determine how many ions each compartment contains */
    compartmentalize_ions(cr, sc, box, step, s->fpout, bRerun);

    /* Output how many ions are in the compartments */
    if (MASTER(cr))
    {
        print_ionlist(s, t, "");
    }

    /* If we are doing a rerun, we are finished here, since we cannot perform
     * swaps anyway */
    if (bRerun)
    {
        return FALSE;
    }

    /* Do we have to perform a swap? */
    bSwap = need_swap(sc);
    if (bSwap)
    {
        g = &(s->group[eGrpSolvent]);
        communicate_group_positions(cr, g->xc, NULL, NULL, FALSE,
                                    x, g->nat, g->nat_loc, g->ind_loc, g->c_ind_loc, NULL, NULL);

        compartmentalize_solvent(cr, sc, box, s->fpout);

        /* Determine where ions are missing and where ions are too many */
        for (ic = 0; ic < eCompNR; ic++)
        {
            for (ii = 0; ii < eIonNR; ii++)
            {
                vacancy[ic][ii] = s->comp[ic][ii].nat_req - s->comp[ic][ii].nat_av;
            }
        }

        /* Remember the original number of ions per compartment */
        for (ic = 0; ic < eCompNR; ic++)
        {
            s->compsol[ic].nat_old = s->compsol[ic].nat;
            for (ii = 0; ii < eIonNR; ii++)
            {
                s->comp[ic][ii].nat_old = s->comp[ic][ii].nat;
            }
        }

        /* Now actually correct the number of ions */
        nswaps = 0;
        alook  = gmx_mtop_atomlookup_init(mtop);
        for (ic = 0; ic < eCompNR; ic++)
        {
            for (ii = 0; ii < eIonNR; ii++)
            {
                while (vacancy[ic][ii] >= sc->threshold)
                {
                    /* Swap in an ion */

                    /* Get the xc-index of the first atom of a solvent molecule of this compartment */
                    isol = get_index_of_distant_atom(&(s->compsol[ic]), s->group[eGrpSolvent].apm );

                    /* Get the xc-index of an ion from the other compartment */
                    iion = get_index_of_distant_atom(&(s->comp[(ic+1)%eCompNR][ii]), s->group[eGrpIons].apm );

                    /* Get the solvent molecule's center of mass */
                    for (im = 0; im < s->group[eGrpSolvent].apm; im++)
                    {
                        gmx_mtop_atomnr_to_atom(alook, s->group[eGrpSolvent].ind[isol+im], &atom);
                        s->group[eGrpSolvent].m[im] = atom->m;
                    }
                    get_molecule_center(&(s->group[eGrpSolvent].xc[isol]), s->group[eGrpSolvent].apm, s->group[eGrpSolvent].m, solvent_center, s->pbc);
                    get_molecule_center(&(s->group[eGrpIons   ].xc[iion]), s->group[eGrpIons   ].apm, NULL, ion_center, s->pbc);

                    /* subtract com_solvent and add com_ion */
                    translate_positions(&(s->group[eGrpSolvent].xc[isol]), s->group[eGrpSolvent].apm, solvent_center, ion_center, s->pbc);
                    /* For the ion, subtract com_ion and add com_solvent */
                    translate_positions(&(s->group[eGrpIons   ].xc[iion]), s->group[eGrpIons   ].apm, ion_center, solvent_center, s->pbc);

                    vacancy[ic              ][ii]--;
                    vacancy[(ic+1) % eCompNR][ii]++;

                    /* Keep track of the changes */
                    s->comp[ic              ][ii].nat++;
                    s->comp[(ic+1) % eCompNR][ii].nat--;
                    s->comp[ic              ][ii].inflow_netto++;
                    s->comp[(ic+1) % eCompNR][ii].inflow_netto--;
                    /* Correct the past time window to still get the right averages from now on */
                    s->comp[ic              ][ii].nat_av++;
                    s->comp[(ic+1) % eCompNR][ii].nat_av--;
                    for (j = 0; j < sc->nAverage; j++)
                    {
                        s->comp[ic              ][ii].nat_past[j]++;
                        s->comp[(ic+1) % eCompNR][ii].nat_past[j]--;
                    }
                    /* Clear ion history */
                    if (MASTER(cr))
                    {
                        s->group[eGrpIons].channel_label[iion] = eChHistPassedNone;
                        s->group[eGrpIons].comp_from[iion]     = eDomainNotset;
                    }
                    /* That was the swap */
                    nswaps++;
                }
            }
        }
        gmx_mtop_atomlookup_destroy(alook);

        if (bVerbose)
        {
            fprintf(stderr, "%s Performed %d swap%s in step %" GMX_PRId64 ".\n", SwS, nswaps, nswaps > 1 ? "s" : "", step);
        }
        if (s->fpout != NULL)
        {
            print_ionlist(s, t, "  # after swap");
        }

        /* Write back the the modified local positions from the collective array to the official coordinates */
        apply_modified_positions(&(s->group[eGrpIons   ]), x);
        apply_modified_positions(&(s->group[eGrpSolvent]), x);
    } /* end of if(bSwap) */

    wallcycle_stop(wcycle, ewcSWAP);

    return bSwap;
}
