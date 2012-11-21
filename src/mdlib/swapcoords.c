#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "futil.h"
#include "swapcoords.h"
#include "gmx_ga2la.h"
#include "groupcoord.h"
#include "mtop_util.h"
#include "macros.h"
#include "vec.h"
#include "names.h"
#include "network.h"
#include "mdrun.h"
#include "xvgr.h"
#include "gmxfio.h"
#include "copyrite.h"


static char* IonString[eIonNr] = {"anion", "cation" };
static char* IonStr[eIonNr] = {"-", "+" };

/* static char* CompString[eCompNr] = {"Compartment A", "Compartment B" }; */
static char* CompStr[eCompNr] = {"A", "B" };

/* eGrpSplit0 and eGrpSplit1 _must_ be neighbors in this list because
 * we sometimes loop from eGrpSplit0 to eGrpSplit1 */
enum {eGrpIons, eGrpSplit0, eGrpSplit1, eGrpSolvent, eGrpNr};
static char* GrpString[eGrpNr] = { "ion", "split0", "split1", "solvent" };

static char *SwapStr[eSwapTypesNR+1] = {
  "", "X-", "Y-", "Z-", "", NULL
};

static char *DimStr[DIM+1] = {
  "X", "Y", "Z", NULL
};

/* Keep track of through which channel the ions have passed */
enum eChannelHistory { eChHistPassedNone, eChHistPassedCh0, eChHistPassedCh1, eChHistNr };
/* static char* ChannelHistString[eChHistNr] = { "none", "channel0", "channel1" }; */

/* Keep track of from which compartment the ions came before passing the channel */
enum eDomain { eDomainNotset, eDomainA, eDomainB, eDomainNr };
static char* DomainString[eDomainNr] = { "not_assigned", "Domain_A", "Domain_B" };


static char *SwS = {"SWAP:"};


/* Helper structure for sorting positions along rotation vector */
typedef struct {
    real xcproj;            /* Projection of xc on the direction vector       */
    int ind;                /* Index of xc                                    */
} sort_along_vec_t;

/* Swap coordinate data  */
typedef struct compartment
{
    int  nat;            /* Number of atoms matching the
                            compartment conditions                    */
    int  nat_old;        /* Number of atoms before swapping           */
    int  nat_req;        /* Requested number of atoms                 */
    real nat_av;         /* Time-averaged number of atoms matching
                            the compartment conditions                */
    int  *nat_past;      /* Past ion counts for time-averaging        */
    int  *ind;           /* Indices to coll array of atoms            */
    real *dist;          /* Distance of atom to compartment center    */
    int  nalloc;         /* Allocation size for ind array             */
    int inflow_netto;
} t_compartment;


typedef struct group
{
    int      nat;           /* Number of atoms in the group                   */
    int      apm;           /* Number of atoms in each molecule               */
    atom_id  *ind;          /* Global atom indices of the group               */
    atom_id  *ind_loc;      /* Local atom indices of the group                */
    int      nat_loc;       /* Number of local group atoms                    */
    int      nalloc_loc;    /* Allocation size for ind_loc                    */
    rvec     *xc;           /* Collective array of group atom positions       */
    real     *qc;           /* Collective array of charges                    */
    int      *c_ind_loc;    /* Position of local atoms in the
                               collective array, [0..nat_loc]                 */
    real     *m;            /* Masses (can be omitted)                        */
    unsigned char *dom_from;/* (Collective) Stores from which compartment this
                               atom has come. This way we keep track of through
                               which channel an ion permeates (only used for
                               the ion group)                                 */
    unsigned char *dom_now; /* Stores in which compartment this ion is now    */
    unsigned char *chan_pass; /* Stores which channel this ion has passed at
                               last                                           */
    rvec     center;        /* Center of the group; COM if masses are used    */
} t_group;


typedef struct swap
{
    int      swapdim;        /* 0=X, 1=Y, 2=Z, -1=auto                        */
    t_pbc    *pbc;           /* Needed for 'auto'                             */
    FILE     *fpout;         /* Output file                                   */
    t_group  group[eGrpNr];
    t_compartment comp[eCompNr][eIonNr];  /* Data for a specific swap volume  *
                                           * and ion type                     */
    t_compartment compsol[eCompNr];       /* Solvent compartments             */
    sort_along_vec_t *data; /* Helper struct for sorting atoms                */
    int       fluxfromAtoB[eChanNr][eIonNr]; /* Netto flux through each of the
                                                channels for each ion type    */
    int       cyl0ions;      /* Ions residing in channel 0                    */
    int       cyl1ions;      /* dito for channel 1                            */
    int       cyl0and1;      /* Ions assigned to cyl0 and cyl1. Not good.     */
    int       *fluxleak;     /* Pointer to a single int value holding the
                                flux not going through any of the channels    */
    real      deltaQ;        /* The charge imbalance between the compartments */
} t_swap;


//static void dump_xcoll(FILE *fp, int nat, rvec x[], real q[], t_swapcoords *swap)
//{
//    int i;
//
////    ddglatnr(gmx_domdec_t *dd,int i);
//
//    for (i=0; i<nat; i++)
//        fprintf(fp, "%d  %9.5f %9.5f %9.5f  %3.1f\n",
//                swap->ind[i], x[i][XX], x[i][YY], x[i][ZZ], q[i]);
//}


/* Check whether point is in channel. Channel is a cylinder defined by a disc
 * with radius r around its center c. The thickness of the cylinder is
 * d_up - d_down.
 *
 *               ^  d_up
 *               |
 *     r         |
 *     <---------c--------->
 *               |
 *               v  d_down
 *
 */
static gmx_bool is_in_channel(
        rvec point,  /* Point under consideration */
        rvec center, /* 'Center' of cylinder */
        real d_up,   /* Upper extension */
        real d_down, /* Lower extensions */
        real r_cyl2, /* Cylinder radius squared */
        t_pbc *pbc,
        int normal)  /* The membrane normal direction is typically 3, i.e. ZZ, but can be X or Y also */
{
    rvec dr;
    int plane1, plane2; /* Directions tangential to membrane */


    plane1 = (normal + 1) % 3; /* typically 0, i.e. XX */
    plane2 = (normal + 2) % 3; /* typically 1, i.e. YY */

    /* Get the distance vector dr between the point and the center of the cylinder */
    pbc_dx(pbc, point, center, dr); /* This puts center in the origin */

    /* Check vertical direction */
    if ( (dr[normal] > d_up) || (dr[normal] < -d_down) )
        return FALSE;

    /* Check radial direction */
    if ( (dr[plane1]*dr[plane1] + dr[plane2]*dr[plane2]) > r_cyl2 )
        return FALSE;

    /* All check passed, this point is in the cylinder */
    return TRUE;
}


static void print_ionlist(
        t_swap *s,
        real time,
        char comment[])
{
    int itype,icomp,i,j;
    t_compartment *comp;


    fprintf(s->fpout, "%12.5e", time);
    for (icomp=0; icomp<eCompNr; icomp++)
    {
        for (itype = 0; itype < eIonNr; itype++)
        {
            comp = &(s->comp[icomp][itype]);
            fprintf(s->fpout, "%7d%7.1f%7d", comp->nat, comp->nat_av-comp->nat_req, comp->inflow_netto);
        }
    }
    if (s->swapdim < 0)
    {
        fprintf(s->fpout, "  %f %f %f   %f %f %f",
                s->group[eGrpSplit0].center[XX], s->group[eGrpSplit0].center[YY],s->group[eGrpSplit0].center[ZZ],
                s->group[eGrpSplit1].center[XX], s->group[eGrpSplit1].center[YY],s->group[eGrpSplit1].center[ZZ]);
    }
    else
    {
        fprintf(s->fpout, "%12.3e%12.3e",
                s->group[eGrpSplit0].center[s->swapdim],
                s->group[eGrpSplit1].center[s->swapdim]);
    }

    for (i=0; i<eChanNr; i++)
        for (j=0; j<eIonNr; j++)
            fprintf(s->fpout, "%12d", s->fluxfromAtoB[i][j]);

    /* Also print the number of ions that leaked from A to B: */
    fprintf(s->fpout, "%12d", *s->fluxleak);

    fprintf(s->fpout, "%s\n", comment);
}


/* Get the center of a group of nat atoms */
static void get_molecule_center(
        rvec x[],
        int nat,
        real *weights,
        rvec center)
{
    int i;
    rvec wpos;
    real wsum;


    clear_rvec(center);

    if (NULL == weights)
    {
        /* Calculate the center of geometry */
        for (i=0; i<nat; i++)
        {
            rvec_inc(center, x[i]);
        }
        wsum = nat;
    }
    else
    {
        /* Calculate weighted center */
        wsum = 0.0;
        for (i=0; i<nat; i++)
        {
            wsum += weights[i];
            svmul(weights[i], x[i], wpos);
            rvec_inc(center, wpos);
        }
    }

    svmul(1.0/wsum, center, center);
}


/* Returns TRUE if x is inbetween (w1+gap) and (w2-gap)
 *
 * ||-----------|--+--|----------o----------|--+--|---------------------||
 *                w1   ?????????????????????  w2
 *
 */
static gmx_bool compartment_contains_atom(
        real w1,  /* position of wall atom 1 */
        real w2,  /* position of wall atom 2 */
        real gap,
        real x,
        real l,   /* length of the box, from || to || in the sketch */
        real *distance_from_center)
{
    real m, l_2;


    /* First set the origin in the middle of w1 and w2 */
    m = 0.5 * (w1 + w2);
    w1 -= m;
    w2 -= m;
    x  -= m;

    /* Now choose the PBC image of x that is closest to the origin: */
    l_2 = 0.5*l;
    while (x  > l_2)
        x -= l;
    while (x <= -l_2)
        x += l;

    *distance_from_center = (real)fabs(x);

    /* Return TRUE if we now are in area "????" */
    if ( (x >= (w1+gap)) &&  (x < (w2-gap)) )
        return TRUE;
    else
        return FALSE;
}


static void update_time_window(t_compartment *comp, int values, int replace)
{
    real average;
    int i;


    /* Put in the new value */
    if (replace >= 0)
        comp->nat_past[replace] = comp->nat;

    /* Compute the new time-average */
    average = 0.0;
    for (i=0; i<values; i++)
        average += comp->nat_past[i];
    average /= values;
    comp->nat_av = average;
}


/* Add atom with collective index ci to the list 'comp' */
static void add_to_list(
        int ci,                  /* index of this ion in the collective array xc, qc */
        t_compartment *comp,     /* Compartment to add this atom to */
        real distance)           /* Shortest distance of this atom to the compartment center */
{
    int nr;


    nr = comp->nat;

    if (nr >= comp->nalloc)
    {
        comp->nalloc = over_alloc_dd(nr+1);
        srenew(comp->ind , comp->nalloc);
        srenew(comp->dist, comp->nalloc);
    }
    comp->ind[nr] = ci;
    comp->dist[nr] = distance;
    comp->nat++;
}


static void get_compartment_boundaries(
        int c,
        t_swap *s,
        matrix box,
        real *left, real *right)
{
    real pos0, pos1;
    real leftpos, rightpos, leftpos_orig;


    if (c >= eCompNr)
        gmx_fatal(FARGS, "No compartment %d.", c);

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
        leftpos  = rightpos;
        rightpos = leftpos_orig + box[s->swapdim][s->swapdim];
    }

    *left = leftpos;
    *right = rightpos;
}


/* To determine the flux through the individual channels, we
 * remember the compartment and channel history of each ion. An ion can be
 * either in channel0 or channel1, or in the remaining volume of compartment
 * A or B.
 *
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
 */
static void detect_flux_per_channel(
        int i,
        int comp,
        int iontype,
        rvec ion_pos,
        unsigned char *dom_now,
        unsigned char *dom_from,
        unsigned char *chan_pass,
        t_swapcoords  *sc,
        real cyl0_r2,
        real cyl1_r2,
        gmx_bool bRerun)
{
    gmx_swapcoords_t s;
    int sd, chan_nr;
    gmx_bool in_cyl0, in_cyl1;


    s    = sc->si_priv;
    sd   = s->swapdim;

    /* Check whether ion is inside any of the channels */
    in_cyl0 = is_in_channel(ion_pos, s->group[eGrpSplit0].center, sc->cyl0u, sc->cyl0l, cyl0_r2, s->pbc, sd);
    in_cyl1 = is_in_channel(ion_pos, s->group[eGrpSplit1].center, sc->cyl1u, sc->cyl1l, cyl1_r2, s->pbc, sd);

    if (in_cyl0 && in_cyl1)
    {
        /* Ion appears to be in both channels. Something is severely wrong! */
        s->cyl0and1++;
        *dom_now   = eDomainNotset;
        *dom_from  = eDomainNotset;
        *chan_pass = eChHistPassedNone;
    }
    else if (in_cyl0)
    {
        /* Ion is in channel 0 now */
        *chan_pass = eChHistPassedCh0;
        *dom_now   = eDomainNotset;
        s->cyl0ions++;
    }
    else if (in_cyl1)
    {
        /* Ion is in channel 1 now */
        *chan_pass = eChHistPassedCh1;
        *dom_now   = eDomainNotset;
        s->cyl1ions++;
    }
    else
    {
        /* Ion is not in any of the channels, so it must be in domain A or B */
        if (eCompA == comp)
            *dom_now = eDomainA;
        else
            *dom_now = eDomainB;
    }

    /* Only take action, if ion is now in domain A or B, and was before
     * in the other domain!
     */
    if (eDomainNotset == *dom_from)
    {
        /* Maybe we can set the domain now */
        *dom_from = *dom_now; /* Could still be eDomainNotset, though */
    }
    else if (  (*dom_now  != eDomainNotset )   /* if in channel */
            && (*dom_from != *dom_now)  )
    {
        /* Obviously the ion changed its domain.
         * Count this for the channel through which it has passed. */
        switch (*chan_pass)
        {
            case eChHistPassedNone:
                *s->fluxleak = *s->fluxleak + 1;

                fprintf(stderr, " %s Ion %d (%s) moved from %s to %s", SwS, i, IonStr[iontype], DomainString[*dom_from],DomainString[*dom_now]);
                if (bRerun)
                    fprintf(stderr, ", possibly due to a swap in the original simulation.\n");
                else
                    fprintf(stderr, " but did not pass any of the split cylinders!\n");
                break;
            case eChHistPassedCh0:
            case eChHistPassedCh1:
                if (*chan_pass == eChHistPassedCh0)
                    chan_nr = 0;
                else
                    chan_nr = 1;

                if (eDomainA == *dom_from)
                    s->fluxfromAtoB[chan_nr][iontype]++;
                else
                    s->fluxfromAtoB[chan_nr][iontype]--;
                break;
            default:
                gmx_fatal(FARGS, "%s Unknown channel history entry!\n", SwS);
                break;
        }

        /* This ion has moved to the _other_ compartment ... */
        *dom_from = *dom_now;
        /* ... and it did not pass any channel yet */
        *chan_pass = eChHistPassedNone;
    }
}


/* Get the lists of ions for the two compartments */
static void compartmentalize_ions(
        t_commrec *cr,
        t_swapcoords *sc,
        matrix box,
        gmx_large_int_t step,
        gmx_bool bVerbose,
        FILE *fpout,
        gmx_bool bRerun)
{
    gmx_swapcoords_t s;
    int i,sd,replace;
    real left, right;
    t_group *iong;
    real dist;
    real cyl0_r2,cyl1_r2;
    int comp,type;
    int sum,not_in_comp[eCompNr]; /* consistency check */


    s    = sc->si_priv;
    iong = &s->group[eGrpIons];
    sd   = s->swapdim;

    cyl0_r2 = sc->cyl0r * sc->cyl0r;
    cyl1_r2 = sc->cyl1r * sc->cyl1r;


    /* Get us a counter that cycles in the range of [0 ... sc->csteps[ */
    replace = (step/sc->nstswap) % sc->csteps;

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
            type = iong->qc[i]<0? eIonNEG : eIonPOS;

            /* Is this ion in the compartment that we look at? */
            if (compartment_contains_atom(left, right, 0, iong->xc[i][sd], box[sd][sd], &dist) )
            {
                /* Now put it into the list containing only ions of its type */
                add_to_list(i, &s->comp[comp][type], dist);

                /* Correct the time-averaged number of ions for this compartment */
                update_time_window(&s->comp[comp][type],sc->csteps,replace);

                /* Master also checks through which channel each ion has passed */
                if (MASTER(cr) && (iong->dom_now != NULL))
                {
                    detect_flux_per_channel(i, comp, type, iong->xc[i],
                            &iong->dom_now[i], &iong->dom_from[i], &iong->chan_pass[i],
                            sc, cyl0_r2, cyl1_r2, bRerun);
                }
            }
            else
            {
                not_in_comp[comp] += 1;
            }
        }
    }

    /* Flux detection warnings */
    if ( MASTER(cr) )
    {
        if (s->cyl0and1 > 0)
        {
            fprintf(stderr, "\n"
                    "%s Warning: %d atoms were detected as being in both channels! Probably your split\n"
                    "%s          cylinder is way too large, or one compartment has collapsed (step ",
                    SwS, s->cyl0and1, SwS);
            fprintf(stderr,gmx_large_int_pfmt,step);
            fprintf(stderr, ")\n");

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
            fprintf(stderr, "%s node %d: Inconsistency during ion compartmentalization. !inA: %d, !inB: %d, total ions %d\n",
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
            fprintf(stderr, "%s node %d: %d atoms are in the ion group, but altogether %d have been assigned to the compartments.\n",
                    SwS, cr->nodeid, iong->nat, sum);
        }
    }


}


/* Set up the compartments and get lists of solvent atoms in each compartment */
static void compartmentalize_solvent(
        t_commrec *cr,
        t_swapcoords *sc,
        matrix box,
        gmx_large_int_t step,
        gmx_bool bVerbose,
        FILE *fpout)
{
    gmx_swapcoords_t s;
    int apm,i,j,sd;
    real left, right;
    t_group *solg;
    real dist;
    int comp;
    int sum, not_in_comp[eCompNr]; /* consistency check */


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
                    add_to_list(i+j, &s->compsol[comp], dist);
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
            fprintf(stderr, "%s node %d: Inconsistency during solvent compartmentalization. !inA: %d, !inB: %d, solvent atoms %d\n",
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
            fprintf(stderr, "%s node %d: %d atoms in solvent group, but %d have been assigned to the compartments.\n",
                    SwS, cr->nodeid, solg->nat, sum);
        }
    }
}


static void mark_molecule(
        unsigned char *mask, int nat, unsigned char domain)
{
    int i;


    for (i=0; i<nat; i++)
        mask[i] = domain;
}


static void compartmentalize_auto(
        t_commrec *cr,
        t_swapcoords *sc,
        gmx_large_int_t step,
        t_forcerec *fr)
{
    int i, sd;
    t_group *g;
    gmx_swapcoords_t s;
    rvec position;
    real cyl0_r2, cyl1_r2;


    s = sc->si_priv;
    sd = s->swapdim;
    cyl0_r2 = sc->cyl0r * sc->cyl0r;
    cyl1_r2 = sc->cyl1r * sc->cyl1r;

    /* First mask all water molecules in the vicinity of the split groups
     * so that we are able to get two distinct water compartments */
    g = &s->group[eGrpSolvent];

    /* Clear solvent mask from last swap step */
    for (i=0; i<g->nat; i++)
        g->dom_now[i] = eChHistPassedNone;

    for (i=0; i<g->nat; i += g->apm)
    {
        /* Get the center of mass of the solvent molecule. Note that
         * g->m only contains masses for a single molecule since all are
         * assumed to be equal */
        get_molecule_center(&g->xc[i], g->apm, g->m, position);

        /* Check whether COM is inside split sphere radius. If yes,
         * mark it as being in the channel */
        if (is_in_channel(position, s->group[eGrpSplit0].center, sc->cyl0u, sc->cyl0l, cyl0_r2, s->pbc, sd))
            mark_molecule(&g->dom_now[i], g->apm, eChHistPassedCh0);
        if (is_in_channel(position, s->group[eGrpSplit1].center, sc->cyl1u, sc->cyl1l, cyl1_r2, s->pbc, sd))
            mark_molecule(&g->dom_now[i], g->apm, eChHistPassedCh1);
    }
}


/* Find out how many group atoms are in the compartments initially */
static void get_initial_ioncounts(
        t_inputrec       *ir,
        swapstate_t      *swapstate,
        rvec             x[],    /* the initial positions */
        matrix           box,
        t_commrec        *cr,
        gmx_bool         bVerbose,
        gmx_bool         bRerun)
{
    t_swapcoords *sc;
    t_swap *s;
    int i, ii, ind, ic;
    int req[2],tot[2];


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
    compartmentalize_ions(cr, sc, box, 0, bVerbose, s->fpout, bRerun);

    /* Set initial concentrations if requested */
    for (ic=0; ic<eCompNr; ic++)
    {
        s->comp[ic][eIonPOS].nat_req = sc->ncations[ic];
        s->comp[ic][eIonNEG].nat_req = sc->nanions[ic];
    }
    for (ic = 0; ic < eCompNr; ic++)
    {
        for (ii = 0; ii < eIonNr; ii++)
        {
            if (s->comp[ic][ii].nat_req < 0)
            {
                s->comp[ic][ii].nat_req = s->comp[ic][ii].nat;
            }
        }
    }

    /* Check whether the number of requested ions adds up to the total number of ions */
    for (ii = 0; ii < eIonNr; ii++)
    {
        req[ii] = s->comp[eCompA][ii].nat_req + s->comp[eCompB][ii].nat_req;
        tot[ii] = s->comp[eCompA][ii].nat     + s->comp[eCompB][ii].nat    ;
    }
    if ( (req[eCompA] != tot[eCompA]) || (req[eCompB] != tot[eCompB ]) )
        gmx_fatal(FARGS, "Mismatch of the number of ions summed over both compartments.\n"
                         "You requested a total of %d anions and %d cations,\n"
                         "but there are a total of %d anions and %d cations in the system.\n",
                  req[eIonNEG], req[eIonPOS],
                  tot[eIonNEG], tot[eIonPOS]);

    /* Initialize time-averaging:
     * Write initial concentrations to all time bins to start with */
    for (ic = 0; ic < eCompNr; ic++)
    {
        for (ii = 0; ii < eIonNr; ii++)
        {
            s->comp[ic][ii].nat_av = s->comp[ic][ii].nat;
            for (i=0; i < sc->csteps; i++)
                s->comp[ic][ii].nat_past[i] = s->comp[ic][ii].nat;
        }
    }
}


static void get_initial_ioncounts_from_cpt(
        t_inputrec *ir, swapstate_t *swapstate,
        t_commrec *cr, gmx_bool bVerbose)
{
    t_swapcoords *sc;
    t_swap *s;
    int ic,ii,j;

    sc = ir->swap;
    s = sc->si_priv;

    if (MASTER(cr))
    {
        /* Copy the past values from the checkpoint values that have been read in already */
        if (bVerbose)
            fprintf(stderr, "%s Copying values from checkpoint\n", SwS);

        for (ic=0; ic<eCompNr; ic++)
        {
            for (ii=0; ii<eIonNr; ii++)
            {
                s->comp[ic][ii].nat_req = swapstate->nat_req[ic][ii];
                s->comp[ic][ii].inflow_netto = swapstate->inflow_netto[ic][ii];

                if (bVerbose)
                    fprintf(stderr, "%s ... Influx netto: %d   Requested: %d   Past values: ", SwS,
                            s->comp[ic][ii].inflow_netto, s->comp[ic][ii].nat_req);

                for (j=0; j < sc->csteps; j++)
                {
                    s->comp[ic][ii].nat_past[j] = swapstate->nat_past[ic][ii][j];
                    if (bVerbose)
                        fprintf(stderr, "%d ", s->comp[ic][ii].nat_past[j]);
                }
                if (bVerbose)
                    fprintf(stderr, "\n");
            }
        }
    }
}


static void bc_initial_concentrations(
        t_commrec *cr,
        t_swapcoords *swap)
{
    int ic, ii;
    t_swap *s;

    s = swap->si_priv;

    for (ic = 0; ic < eCompNr; ic++)
    {
        for (ii = 0; ii < eIonNr; ii++)
        {
            gmx_bcast(sizeof(s->comp[ic][ii].nat_req), &(s->comp[ic][ii].nat_req), cr);
            gmx_bcast(sizeof(s->comp[ic][ii].nat    ), &(s->comp[ic][ii].nat    ), cr);
            gmx_bcast( swap->csteps * sizeof(s->comp[ic][ii].nat_past[0]), s->comp[ic][ii].nat_past, cr);
        }
    }
}


static void ensure_that_groups_differ(t_swap *s, gmx_bool bVerbose)
{
    t_group *ga, *gb;
    int i,j,k;
    gmx_bool bSame;


    if (bVerbose)
        fprintf(stderr, "%s Making shure the groups have no overlapping atoms.\n", SwS);

    for (i = 0; i<eGrpNr; i++)
    {
        ga = &s->group[i];
        for (j = i+1; j<eGrpNr; j++)
        {
            gb = &s->group[j];
            if (bVerbose)
                fprintf(stderr, "%s ... comparing %s and %s group (%d and %d atoms)\n", SwS, GrpString[i], GrpString[j], ga->nat,gb->nat);
            if (ga->nat == gb->nat)
            {
                bSame = TRUE;
                for (k = 0; k<ga->nat; k++)
                {
                    if (ga->ind[k] != gb->ind[k])
                        bSame = FALSE;
                }
                if (TRUE == bSame)
                    gmx_fatal(FARGS, "%s and %s groups are identical. Cannot perform swapping.", GrpString[i], GrpString[j]);
            }
        }
    }
}


static int get_group_apm_check(
        int group,
        t_swap *s,
        gmx_bool bVerbose,
        gmx_mtop_t *mtop)
{
    int *ind;
    int nat, apm, i;
    int molb, molnr, atnr_mol;



    ind = s->group[group].ind;
    nat = s->group[group].nat;

    /* Determine the number of solvent atoms per solvent molecule from the
     * first solvent atom: */
    i=0;
    gmx_mtop_atomnr_to_molblock_ind(mtop, ind[i],&molb,&molnr,&atnr_mol);
    apm = mtop->molblock[molb].natoms_mol;

    if (bVerbose)
        fprintf(stderr, "%s Checking whether all %s molecules consist of %d atom%s\n",
                SwS, GrpString[group], apm, apm>1? "s":"");

    /* Check whether this is also true for all other solvent atoms */
    for (i=1; i<nat; i++)
    {
        gmx_mtop_atomnr_to_molblock_ind(mtop, ind[i],&molb,&molnr,&atnr_mol);
        if (apm != mtop->molblock[molb].natoms_mol)
            gmx_fatal(FARGS, "Not all %s group molecules consist of %d atoms.",
                    GrpString[group], apm);
    }

    return apm;
}


/* Print the legend to the swapcoords output file as well as the
 * initial ion counts */
static void print_ionlist_legend(t_inputrec *ir, const output_env_t oenv)
{
    const char **legend;
    int ic,count,ii;
    char buf[256];
    t_swap *s;


    s = ir->swap->si_priv;

    snew(legend, eCompNr*eIonNr*3 + 2 + eChanNr*eIonNr + 1);
    for (ic=count=0; ic<eCompNr; ic++)
    {
        for (ii=0; ii<eIonNr; ii++)
        {
            sprintf(buf,"%s %ss", CompStr[ic], IonString[ii]);
            legend[count++] = gmx_strdup(buf);
            sprintf(buf,"%s av. mismatch to %d%s",
                    CompStr[ic], s->comp[ic][ii].nat_req, IonStr[ii]);
            legend[count++] = gmx_strdup(buf);
            sprintf(buf,"%s netto %s influx", CompStr[ic], IonString[ii]);
            legend[count++] = gmx_strdup(buf);
        }
    }
    sprintf(buf, "%scenter of %s of split group 0", SwapStr[ir->eSwapCoords], (NULL != s->group[eGrpSplit0].m)? "mass":"geometry");
    legend[count++] = gmx_strdup(buf);
    sprintf(buf, "%scenter of %s of split group 1", SwapStr[ir->eSwapCoords], (NULL != s->group[eGrpSplit1].m)? "mass":"geometry");
    legend[count++] = gmx_strdup(buf);

    for (ic=0; ic<eChanNr; ic++)
    {
        for (ii=0; ii<eIonNr; ii++)
        {
            sprintf(buf, "A->ch%d->B %s permeations", ic, IonString[ii]);
            legend[count++] = gmx_strdup(buf);
        }
    }

    sprintf(buf, "leakage");
    legend[count++] = gmx_strdup(buf);

    xvgr_legend(s->fpout,count,legend, oenv);

    fprintf(s->fpout, "# Instantaneous ion counts and time-averaged differences to requested numbers\n");
    fprintf(s->fpout, "#   time[ps]   A_an   diff   t_in  A_cat   diff   t_in   B_an   diff   t_in  B_cat   diff   t_in ");
    fprintf(s->fpout, "   %s-Split0    %s-Split1", DimStr[s->swapdim], DimStr[s->swapdim]);
    fprintf(s->fpout, "  A-ch0-B_an A-ch0-B_cat  A-ch1-B_an A-ch1-B_cat ion_leakage\n");
    fflush(s->fpout);
}


/* Initialize arrays that keep track of where the ions come from and where
 * they go */
static void detect_flux_per_channel_init(
        t_commrec *cr,
        t_swap *s,
        swapstate_t *swapstate,
        gmx_bool bStartFromCpt)
{
    int i, ic, ii;
    t_group *g;


    g = &(s->group[eGrpIons]);

    /* All these flux detection routines run on the master only */
    if (!MASTER(cr))
    {
        g->dom_now   = NULL;
        g->dom_from  = NULL;
        g->chan_pass = NULL;

        return;
    }

    /******************************************************/
    /* Channel and domain history for the individual ions */
    /******************************************************/
    if (bStartFromCpt) /* set the pointers right */
    {
        g->dom_from  = swapstate->dom_from;
        g->chan_pass = swapstate->chan_pass;
    }
    else /* allocate memory */
    {
        snew(g->dom_from , g->nat);
        swapstate->dom_from = g->dom_from;
        snew(g->chan_pass, g->nat);
        swapstate->chan_pass = g->chan_pass;
    }
    snew(g->dom_now  , g->nat);

    /* Initialize the channel and domain history counters */
    for (i=0; i<g->nat; i++)
    {
        g->dom_now[i] = eDomainNotset;
        if (!bStartFromCpt)
        {
            g->dom_from[i]  = eDomainNotset;
            g->chan_pass[i] = eChHistPassedNone;
        }
    }

    /************************************/
    /* Channel fluxes for both channels */
    /************************************/
    s->cyl0ions                 = 0;
    s->cyl1ions                 = 0;
    s->cyl0and1                 = 0;

    if (bStartFromCpt)
        fprintf(stderr, "%s Copying channel fluxes from checkpoint file data\n", SwS);

    for (ic=0; ic<eChanNr; ic++)
    {
        fprintf(stderr, "%s Channel %d flux history: ", SwS, ic);
        for (ii=0; ii<eIonNr; ii++)
        {
            if (bStartFromCpt)
                s->fluxfromAtoB[ic][ii] = swapstate->fluxfromAtoB[ic][ii];
            else
                s->fluxfromAtoB[ic][ii] = 0;

            fprintf(stderr, "%d %s%s   ", s->fluxfromAtoB[ic][ii], IonString[ii], s->fluxfromAtoB[ic][ii] == 1? "":"s");
        }
        fprintf(stderr, "\n");
    }
    if (bStartFromCpt)
        s->fluxleak = swapstate->fluxleak;
    else
    {
        snew(s->fluxleak, 1);
        *s->fluxleak = 0;
        /* Set pointer for checkpoint writing */
        swapstate->fluxleak = s->fluxleak;
    }


    /* Set pointers for checkpoint writing */
    for (ic=0; ic<eChanNr; ic++)
        for (ii=0; ii<eIonNr; ii++)
            swapstate->fluxfromAtoB_p[ic][ii] = &(s->fluxfromAtoB[ic][ii]);
}


extern void init_swapcoords(
        FILE             *fplog, /* general output file md.log */
        gmx_bool         bVerbose,
        t_inputrec       *ir,
        const char      *fn,    /* output file name for swap data */
        gmx_mtop_t       *mtop,
        rvec             x[],    /* the initial positions */
        matrix           box,
        swapstate_t      *swapstate,
        t_commrec        *cr,
        const output_env_t oenv,
        unsigned long    Flags)
{
    int i,ic,ig,ii,j;
    t_swapcoords *sc;
    t_swap *s;
    t_atom *atom;
    t_group *g;
    gmx_bool bAppend, bStartFromCpt, bRerun;


    if ( (PAR(cr)) && !DOMAINDECOMP(cr) )
        gmx_fatal(FARGS, "Position swapping is only implemented for domain decomposition!");

    bAppend       = Flags & MD_APPENDFILES;
    bStartFromCpt = Flags & MD_STARTFROMCPT;
    bRerun        = Flags & MD_RERUN;

    sc = ir->swap;
    snew(sc->si_priv, 1);
    s = sc->si_priv;

    if (bRerun)
    {
        if (PAR(cr))
            gmx_fatal(FARGS, "%s This module does not support reruns in parallel\nPlease request a serial run with -nt 1 / -np 1\n", SwS);

        fprintf(stderr, "%s Rerun - using every available frame\n", SwS);
        sc->nstswap = 1;
        sc->csteps  = 1;  /* averaging makes no sense for reruns */
    }

    if (MASTER(cr) && !bAppend)
    {
        fprintf(fplog, "\nInitializing ion/water position exchanges\n");
        please_cite(fplog,"Kutzner2011");
    }

    switch (ir->eSwapCoords)
    {
        case eswapX:
            s->swapdim=XX;
            break;
        case eswapY:
            s->swapdim=YY;
            break;
        case eswapZ:
            s->swapdim=ZZ;
            break;
        default:
            s->swapdim=-1;
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

    ensure_that_groups_differ(s, bVerbose && MASTER(cr));

    /* Allocate space for the collective arrays for all groups */
    for (ig=0; ig<eGrpNr; ig++)
    {
        g = &(s->group[ig]);
        snew(g->xc, g->nat);
        snew(g->c_ind_loc, g->nat);
    }

    /* Make shure that all molecules in the ion and solvent groups contain the
     * same number of atoms each */
    s->group[eGrpIons   ].apm = get_group_apm_check(eGrpIons   , s, MASTER(cr) && bVerbose, mtop);
    s->group[eGrpSolvent].apm = get_group_apm_check(eGrpSolvent, s, MASTER(cr) && bVerbose, mtop);

    /* Save masses where needed */
    s->group[eGrpIons   ].m = NULL;
    /* We only need enough space to determine a single solvent molecule's
     * center at at time */
    g = &(s->group[eGrpSolvent]);
    snew(g->m, g->apm);
    if (eswapAuto == ir->eSwapCoords)
    {
        if (MASTER(cr))
            fprintf(stderr, "%s Solvent atom masses used for auto-compartmentalization: ", SwS);
        for (i=0; i<g->apm; i++)
        {
            gmx_mtop_atomnr_to_atom(mtop,g->ind[i],&atom);
            g->m[i] = atom->m;
            if (MASTER(cr))
                fprintf(stderr, "%f ", g->m[i]);
        }
        if (MASTER(cr))
            fprintf(stderr, "u\n");
    }

    /* Need mass-weighted center of split group? */
    for (j=0, ig=eGrpSplit0; j<eChanNr; ig++, j++)
    {
        g = &(s->group[ig]);
        if (TRUE == sc->massw_split[j])
        {
            /* Save the split group charges if mass-weighting is requested */
            snew(g->m, g->nat);
            for (i=0; i<g->nat; i++)
            {
                gmx_mtop_atomnr_to_atom(mtop,g->ind[i],&atom);
                g->m[i] = atom->m;
            }
        }
        else
            g->m = NULL;
    }

    /* Save the ionic charges */
    g = &(s->group[eGrpIons]);
    snew(g->qc, g->nat);
    for (i=0; i<g->nat; i++)
    {
        gmx_mtop_atomnr_to_atom(mtop,g->ind[i],&atom);
        g->qc[i] = atom->q;
    }

    /* Auto-compartmentalization needs some extra initialization */
    if (eswapAuto == ir->eSwapCoords)
    {
        g = &(s->group[eGrpSolvent]);
        snew(g->dom_now, g->nat);
    }
    snew(s->pbc, 1);
    set_pbc(s->pbc, -1, box);


    if (MASTER(cr))
    {
        if (bVerbose)
            fprintf(stderr, "%s Opening output file %s%s\n", SwS, fn, bAppend? " for appending":"");

        s->fpout = gmx_fio_fopen(fn, bAppend ? "a+" : "w+" );

        if (!bAppend)
        {
            xvgr_header(s->fpout, "Ion counts", "Time (ps)", "counts",exvggtXNY,oenv);

            for (ig=0; ig<eGrpNr; ig++)
            {
                g = &(s->group[ig]);
                fprintf(s->fpout, "# %s group contains %d atom%s", GrpString[ig], g->nat,(g->nat>1)?"s":"");
                if (eGrpSolvent==ig || eGrpIons==ig)
                {
                    fprintf(s->fpout, " with %d atom%s in each molecule", g->apm, (g->apm>1)?"s":"");
                }
                fprintf(s->fpout, ".\n");
            }

            fprintf(s->fpout, "#\n# Initial positions of split groups:\n");
        }

        for (j=0, ig=eGrpSplit0; j<eChanNr; j++, ig++)
        {
            g = &(s->group[ig]);
            for (i = 0; i < g->nat; i++)
                copy_rvec(x[sc->ind_split[j][i]],g->xc[i]);
            get_molecule_center(g->xc, g->nat, g->m, g->center);
            if (!bAppend)
            {
                if (eswapAuto == ir->eSwapCoords)
                {
                    fprintf(s->fpout, "# %s group center %5f %5f %5f nm\n", GrpString[ig],
                            g->center[XX],g->center[YY],g->center[ZZ]);
                }
                else
                {
                    fprintf(s->fpout, "# %s group %s-center %5f nm\n", GrpString[ig],
                            DimStr[s->swapdim], g->center[s->swapdim]);
                }
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
                        sc->csteps, sc->csteps*sc->nstswap*ir->delta_t);
                fprintf(s->fpout, "# Threshold is %f\n", sc->threshold);
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
        for (ig=0; ig<eGrpNr; ig++)
        {
            g = &(s->group[ig]);
            g->nat_loc    = 0;
            g->nalloc_loc = 0;
            g->ind_loc = NULL;
        }
     }
     else
     {
         for (ig=0; ig<eGrpNr; ig++)
         {
             g = &(s->group[ig]);
             g->nat_loc = g->nat;
             g->ind_loc = g->ind;
             /* c_ind_loc needs to be set to identity in the serial case */
             for (i=0; i<g->nat; i++)
                 g->c_ind_loc[i] = i;
         }
     }

    /* Allocate memory for the ion counts time window */
    for (ic = 0; ic < eCompNr; ic++)
        for (ii = 0; ii < eIonNr; ii++)
            snew(s->comp[ic][ii].nat_past, sc->csteps);

    /* Get the initial ion concentrations and let the other nodes know */
    if (MASTER(cr))
    {
        swapstate->nions = s->group[eGrpIons].nat;

        if (bStartFromCpt)
            get_initial_ioncounts_from_cpt(ir,swapstate,cr, bVerbose);
        else
        {
            fprintf(stderr, "%s Determining initial ion counts.\n", SwS);
            get_initial_ioncounts(ir,swapstate,x,box,cr,bVerbose, bRerun);
        }

        /* Prepare (further) checkpoint writes ... */
        if (bStartFromCpt)
        {
            /* Consistency check */
            if (swapstate->csteps != sc->csteps)
                gmx_fatal(FARGS, "%s Ion count averaging steps mismatch! checkpoint: %d, tpr: %d",
                        SwS, swapstate->csteps, sc->csteps);
        }
        else
            swapstate->csteps = sc->csteps;
        fprintf(stderr, "%s Setting pointers for checkpoint writing\n", SwS);
        for (ic=0; ic<eCompNr; ic++)
        {
            for (ii=0; ii<eIonNr; ii++)
            {
                swapstate->nat_req_p[ic][ii] = &(s->comp[ic][ii].nat_req);
                swapstate->nat_past_p[ic][ii] = &(s->comp[ic][ii].nat_past[0]);
                swapstate->inflow_netto_p[ic][ii] = &(s->comp[ic][ii].inflow_netto);
            }
        }

        /* Determine the total charge imbalance */
        s->deltaQ =  ( (-1) * s->comp[eCompA][eIonNEG].nat_req + s->comp[eCompA][eIonPOS].nat_req )
                   - ( (-1) * s->comp[eCompB][eIonNEG].nat_req + s->comp[eCompB][eIonPOS].nat_req );

        if (bVerbose)
            fprintf(stderr, "%s Requested charge imbalance is Q(A) - Q(B) = %gz.\n", SwS, s->deltaQ);
        fprintf(s->fpout, "# Requested charge imbalance is Q(A)-Q(B) = %gz.\n", s->deltaQ);
    }

    if (PAR(cr))
        bc_initial_concentrations(cr, ir->swap);

    /* Put the time-averaged number of ions for all compartments */
    for (ic = 0; ic < eCompNr; ic++)
        for (ii = 0; ii < eIonNr; ii++)
            update_time_window(&(s->comp[ic][ii]),sc->csteps,-1);

    /* Initialize arrays that keep track of through which channel the ions go */
    detect_flux_per_channel_init(cr, s, swapstate, bStartFromCpt);

    /* We need to print the legend if we open this file for the first time. */
    if (MASTER(cr) && !bAppend)
        print_ionlist_legend(ir,oenv);
}


extern void dd_make_local_swap_groups(t_commrec *cr, gmx_domdec_t *dd,t_swapcoords *sc,t_mdatoms *md)
{
    t_group *g;
    int ig;


    /* Make ion group, split groups and solvent group */
    for (ig=0; ig<eGrpNr; ig++)
    {
        g = &(sc->si_priv->group[ig]);
        dd_make_local_group_indices(dd->ga2la, g->nat, g->ind,
                &(g->nat_loc), &(g->ind_loc),&(g->nalloc_loc), g->c_ind_loc);
    }
}


static gmx_bool need_swap(t_commrec *cr, t_swapcoords *sc)
{
    t_swap *s;
    int ic, ii;


    s = sc->si_priv;
    for (ic=0; ic<eCompNr; ic++)
    {
        for (ii=0; ii<eIonNr; ii++)
        {
            if (s->comp[ic][ii].nat_req - s->comp[ic][ii].nat_av >= sc->threshold)
                return TRUE;
        }
    }
    return FALSE;
}


/* Returns the index of an atom that is far off the compartment boundaries.
 * Other atoms of the molecule (if any) will directly follow the returned index
 */
static int get_index_of_distant_atom(
        real boxlength,
        t_compartment *comp,
        int apm) /* Atoms per molecule - just return the first atom index of a molecule */
{
    int i,ibest=-1;
    real d = GMX_REAL_MAX;


    /* comp->nat contains the original number of atoms in this compartment
     * prior to doing any swaps. Some of these atoms may already have been
     * swapped out, but then they are marked with a distance of GMX_REAL_MAX
     */
    for (i=0; i<comp->nat_old; i += apm)
    {
        if (comp->dist[i] < d)
        {
            ibest = i;
            d = comp->dist[ibest];
        }
    }

    if (ibest < 0)
        gmx_fatal(FARGS, "Could not get index of swap atom. Compartment atoms %d before swaps, atoms per molecule %d.",
                comp->nat_old, apm);

    /* Set the distance of this index to infinity such that it won't get selected again in
     * this time step
     */
    comp->dist[ibest] = GMX_REAL_MAX;

    return comp->ind[ibest];
}


/* Swaps centers of mass */
static void translate_positions(
        rvec *x,
        int apm,
        rvec old_com,
        rvec new_com)
{
    int i;


    for (i=0; i<apm; i++)
    {
        rvec_dec(x[i], old_com);
        rvec_inc(x[i], new_com);
    }
}


/* Write back the the modified local positions from the collective array to the official coordinates */
static void apply_modified_positions(
        t_group *g,
        rvec x[])
{
    int l,ii,cind;


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
        FILE             *fplog,
        gmx_large_int_t  step,
        real             t,
        t_inputrec       *ir,
        rvec             x[],            /* positions of home particles       */
        rvec             v[],            /* velocities of home particles      */
        matrix           box,
        gmx_mtop_t       *mtop,
        t_forcerec       *fr,
        swapstate_t      *swapstate,
        gmx_bool         bVerbose,
        gmx_bool         bRerun)         /* Is this a rerun ?                 */
{
    t_swapcoords *sc;
    t_swap *s;
    int j, ii, ic, ig, im, gmax, nswaps;
    gmx_bool bSwap = FALSE;
    t_group *g;
    real vacancy[eCompNr][eIonNr];
    int isol,iion;
    rvec solvent_center, ion_center;
    t_atom *atom;


    sc  = ir->swap;
    s = sc->si_priv;

    /* Assemble all the positions of the swap group (ig=0), the split groups
     * (ig=1,2), and possibly the solvent group (ig=3) */
    gmax = eGrpNr;
    if (s->swapdim == eswapAuto)
        gmax--;

    for (ig=0; ig<gmax; ig++)
    {
        g = &(s->group[ig]);
        communicate_group_positions_noshift(cr, g->xc, x, g->nat,
                g->nat_loc, g->ind_loc, g->c_ind_loc);
    }

    /* Set up the compartments and get lists of atoms in each compartment */
    g = &s->group[eGrpSplit0];
    get_molecule_center(g->xc, g->nat, g->m, g->center); /* center of split group 0 */

    g = &s->group[eGrpSplit1];
    get_molecule_center(g->xc, g->nat, g->m, g->center); /* center of split group 1 */

    /* Determine how many ions each compartment contains */
    if (eswapAuto == ir->eSwapCoords)
        compartmentalize_auto(cr, sc, step, fr);
    else
        compartmentalize_ions(cr, sc, box, step, bVerbose, s->fpout, bRerun);

    /* Output how many ions are in the compartments */
    if (MASTER(cr))
        print_ionlist(s, t, "");

    /* If we are doing a rerun, we are finished here, since we cannot perform
     * swaps anyway */
    if (bRerun)
        return FALSE;

    /* Do we have to perform a swap? */
    bSwap = need_swap(cr, sc);
    if (bSwap)
    {
        if (ir->eSwapCoords != eswapAuto)
        {
            g = &(s->group[eGrpSolvent]);
            communicate_group_positions_noshift(cr, g->xc, x, g->nat,
                    g->nat_loc, g->ind_loc, g->c_ind_loc);

            compartmentalize_solvent(cr, sc, box, step, bVerbose, s->fpout);
        }

        /* Determine where ions are missing and where ions are too many */
        for (ic = 0; ic < eCompNr; ic++)
            for (ii = 0; ii < eIonNr; ii++)
                vacancy[ic][ii] = s->comp[ic][ii].nat_req - s->comp[ic][ii].nat_av;

        /* Remember the original number of ions per compartment */
        for (ic = 0; ic < eCompNr; ic++)
        {
            s->compsol[ic].nat_old = s->compsol[ic].nat;
            for (ii = 0; ii < eIonNr; ii++)
                s->comp[ic][ii].nat_old = s->comp[ic][ii].nat;
        }

        /* Now actually correct the number of ions */
        g = &(s->group[eGrpSolvent]);
        nswaps = 0;
        for (ic = 0; ic < eCompNr; ic++)
        {
             for (ii = 0; ii < eIonNr; ii++)
             {
                 while (vacancy[ic][ii] >= sc->threshold)
                 {
                     /* Swap in an ion */

                     /* Get the xc-index of the first atom of a solvent molecule of this compartment */
                     isol = get_index_of_distant_atom(box[s->swapdim][s->swapdim],
                             &(s->compsol[ic]), s->group[eGrpSolvent].apm );

                     /* Get the xc-index of an ion from the other compartment */
                     iion = get_index_of_distant_atom(box[s->swapdim][s->swapdim],
                             &(s->comp[(ic+1)%eCompNr][ii]), s->group[eGrpIons].apm );

                     /* Get the solvent molecule's center of mass */
                     for (im=0; im<s->group[eGrpSolvent].apm; im++)
                     {
                         gmx_mtop_atomnr_to_atom(mtop,s->group[eGrpSolvent].ind[isol+im],&atom);
                         s->group[eGrpSolvent].m[im] = atom->m;
                     }
                     get_molecule_center(&(s->group[eGrpSolvent].xc[isol]), s->group[eGrpSolvent].apm, s->group[eGrpSolvent].m, solvent_center);
                     get_molecule_center(&(s->group[eGrpIons   ].xc[iion]), s->group[eGrpIons   ].apm, NULL                   , ion_center    );

                     /* subtract com_solvent and add com_ion */
                     translate_positions(&(s->group[eGrpSolvent].xc[isol]), s->group[eGrpSolvent].apm, solvent_center, ion_center    );
                     /* For the ion, subtract com_ion and add com_solvent */
                     translate_positions(&(s->group[eGrpIons   ].xc[iion]), s->group[eGrpIons   ].apm, ion_center    , solvent_center);

                     vacancy[ic              ][ii]--;
                     vacancy[(ic+1) % eCompNr][ii]++;

                     /* Keep track of the changes */
                     s->comp[ic              ][ii].nat++;
                     s->comp[(ic+1) % eCompNr][ii].nat--;
                     s->comp[ic              ][ii].inflow_netto++;
                     s->comp[(ic+1) % eCompNr][ii].inflow_netto--;
                     /* Correct the past time window to still get the right averages from now on */
                     s->comp[ic              ][ii].nat_av++;
                     s->comp[(ic+1) % eCompNr][ii].nat_av--;
                     for (j=0; j<sc->csteps; j++)
                     {
                         s->comp[ic              ][ii].nat_past[j]++;
                         s->comp[(ic+1) % eCompNr][ii].nat_past[j]--;
                     }
                     /* Clear ion history */
                     if (MASTER(cr))
                     {
                         s->group[eGrpIons].chan_pass[iion] = eChHistPassedNone;
                         s->group[eGrpIons].dom_from[iion]  = eDomainNotset;
                     }
                     /* That was the swap */
                     nswaps++;
                 }
             }
        }

        if (bVerbose)
        {
            fprintf(stderr, "%s Performed %d swap%s in step ", SwS, nswaps, nswaps>1?"s":"");
            fprintf(stderr, gmx_large_int_pfmt, step);
            fprintf(stderr, "\n");
        }
        if (s->fpout != NULL)
            print_ionlist(s, t, "  # after swap");

        /* Write back the the modified local positions from the collective array to the official coordinates */
        apply_modified_positions(&(s->group[eGrpIons   ]), x);
        apply_modified_positions(&(s->group[eGrpSolvent]), x);
    } /* end of if(bSwap) */

    return bSwap;
}
