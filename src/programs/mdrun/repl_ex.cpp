/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "repl_ex.h"

#include "config.h"

#include <cmath>

#include <random>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"


#define PROBABILITYCUTOFF 100
/* we don't bother evaluating if events are more rare than exp(-100) = 3.7x10^-44 */

//! Rank in the multisimulaiton
#define MSRANK(ms, nodeid)  (nodeid)

enum {
    ereTEMP, ereLAMBDA, ereENDSINGLE, ereTL, ereNR
};
const char *erename[ereNR] = { "temperature", "lambda", "end_single_marker", "temperature and lambda"};
/* end_single_marker merely notes the end of single variable replica exchange. All types higher than
   it are multiple replica exchange methods */
/* Eventually, should add 'pressure', 'temperature and pressure', 'lambda_and_pressure', 'temperature_lambda_pressure'?;
   Let's wait until we feel better about the pressure control methods giving exact ensembles.  Right now, we assume constant pressure  */

typedef struct gmx_repl_ex
{
    int       repl;        /* replica ID */
    int       nrepl;       /* total number of replica */
    real      temp;        /* temperature */
    int       type;        /* replica exchange type from ere enum */
    real    **q;           /* quantity, e.g. temperature or lambda; first index is ere, second index is replica ID */
    gmx_bool  bNPT;        /* use constant pressure and temperature */
    real     *pres;        /* replica pressures */
    int      *ind;         /* replica indices */
    int      *allswaps;    /* used for keeping track of all the replica swaps */
    int       nst;         /* replica exchange interval (number of steps) */
    int       nex;         /* number of exchanges per interval */
    int       seed;        /* random seed */
    int       nattempt[2]; /* number of even and odd replica change attempts */
    real     *prob_sum;    /* sum of probabilities */
    int     **nmoves;      /* number of moves between replicas i and j */
    int      *nexchange;   /* i-th element of the array is the number of exchanges between replica i-1 and i */

    /* these are helper arrays for replica exchange; allocated here so they
       don't have to be allocated each time */
    int      *destinations;
    int     **cyclic;
    int     **order;
    int      *tmpswap;
    gmx_bool *incycle;
    gmx_bool *bEx;

    /* helper arrays to hold the quantities that are exchanged */
    real  *prob;
    real  *Epot;
    real  *beta;
    real  *Vol;
    real **de;

} t_gmx_repl_ex;

static gmx_bool repl_quantity(const gmx_multisim_t *ms,
                              struct gmx_repl_ex *re, int ere, real q)
{
    real    *qall;
    gmx_bool bDiff;
    int      s;

    snew(qall, ms->nsim);
    qall[re->repl] = q;
    gmx_sum_sim(ms->nsim, qall, ms);

    bDiff = FALSE;
    for (s = 1; s < ms->nsim; s++)
    {
        if (qall[s] != qall[0])
        {
            bDiff = TRUE;
        }
    }

    if (bDiff)
    {
        /* Set the replica exchange type and quantities */
        re->type = ere;

        snew(re->q[ere], re->nrepl);
        for (s = 0; s < ms->nsim; s++)
        {
            re->q[ere][s] = qall[s];
        }
    }
    sfree(qall);
    return bDiff;
}

gmx_repl_ex_t
init_replica_exchange(FILE                            *fplog,
                      const gmx_multisim_t            *ms,
                      int                              numAtomsInSystem,
                      const t_inputrec                *ir,
                      const ReplicaExchangeParameters &replExParams)
{
    real                pres;
    int                 i, j, k;
    struct gmx_repl_ex *re;
    gmx_bool            bTemp;
    gmx_bool            bLambda = FALSE;

    fprintf(fplog, "\nInitializing Replica Exchange\n");

    if (ms == nullptr || ms->nsim == 1)
    {
        gmx_fatal(FARGS, "Nothing to exchange with only one replica, maybe you forgot to set the -multi option of mdrun?");
    }
    if (!EI_DYNAMICS(ir->eI))
    {
        gmx_fatal(FARGS, "Replica exchange is only supported by dynamical simulations");
        /* Note that PAR(cr) is defined by cr->nnodes > 1, which is
         * distinct from MULTISIM(cr). A multi-simulation only runs
         * with real MPI parallelism, but this does not imply PAR(cr)
         * is true!
         *
         * Since we are using a dynamical integrator, the only
         * decomposition is DD, so PAR(cr) and DOMAINDECOMP(cr) are
         * synonymous. The only way for cr->nnodes > 1 to be true is
         * if we are using DD. */
    }

    snew(re, 1);

    re->repl     = ms->sim;
    re->nrepl    = ms->nsim;
    snew(re->q, ereENDSINGLE);

    fprintf(fplog, "Repl  There are %d replicas:\n", re->nrepl);

    /* We only check that the number of atoms in the systms match.
     * This, of course, do not guarantee that the systems are the same,
     * but it does guarantee that we can perform replica exchange.
     */
    check_multi_int(fplog, ms, numAtomsInSystem, "the number of atoms", FALSE);
    check_multi_int(fplog, ms, ir->eI, "the integrator", FALSE);
    check_multi_int64(fplog, ms, ir->init_step+ir->nsteps, "init_step+nsteps", FALSE);
    const int nst = replExParams.exchangeInterval;
    check_multi_int64(fplog, ms, (ir->init_step+nst-1)/nst,
                      "first exchange step: init_step/-replex", FALSE);
    check_multi_int(fplog, ms, ir->etc, "the temperature coupling", FALSE);
    check_multi_int(fplog, ms, ir->opts.ngtc,
                    "the number of temperature coupling groups", FALSE);
    check_multi_int(fplog, ms, ir->epc, "the pressure coupling", FALSE);
    check_multi_int(fplog, ms, ir->efep, "free energy", FALSE);
    check_multi_int(fplog, ms, ir->fepvals->n_lambda, "number of lambda states", FALSE);

    re->temp = ir->opts.ref_t[0];
    for (i = 1; (i < ir->opts.ngtc); i++)
    {
        if (ir->opts.ref_t[i] != re->temp)
        {
            fprintf(fplog, "\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
            fprintf(stderr, "\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
        }
    }

    re->type = -1;
    bTemp    = repl_quantity(ms, re, ereTEMP, re->temp);
    if (ir->efep != efepNO)
    {
        bLambda = repl_quantity(ms, re, ereLAMBDA, (real)ir->fepvals->init_fep_state);
    }
    if (re->type == -1)  /* nothing was assigned */
    {
        gmx_fatal(FARGS, "The properties of the %d systems are all the same, there is nothing to exchange", re->nrepl);
    }
    if (bLambda && bTemp)
    {
        re->type = ereTL;
    }

    if (bTemp)
    {
        please_cite(fplog, "Sugita1999a");
        if (ir->epc != epcNO)
        {
            re->bNPT = TRUE;
            fprintf(fplog, "Repl  Using Constant Pressure REMD.\n");
            please_cite(fplog, "Okabe2001a");
        }
        if (ir->etc == etcBERENDSEN)
        {
            gmx_fatal(FARGS, "REMD with the %s thermostat does not produce correct potential energy distributions, consider using the %s thermostat instead",
                      ETCOUPLTYPE(ir->etc), ETCOUPLTYPE(etcVRESCALE));
        }
    }
    if (bLambda)
    {
        if (ir->fepvals->delta_lambda != 0)   /* check this? */
        {
            gmx_fatal(FARGS, "delta_lambda is not zero");
        }
    }
    if (re->bNPT)
    {
        snew(re->pres, re->nrepl);
        if (ir->epct == epctSURFACETENSION)
        {
            pres = ir->ref_p[ZZ][ZZ];
        }
        else
        {
            pres = 0;
            j    = 0;
            for (i = 0; i < DIM; i++)
            {
                if (ir->compress[i][i] != 0)
                {
                    pres += ir->ref_p[i][i];
                    j++;
                }
            }
            pres /= j;
        }
        re->pres[re->repl] = pres;
        gmx_sum_sim(re->nrepl, re->pres, ms);
    }

    /* Make an index for increasing replica order */
    /* only makes sense if one or the other is varying, not both!
       if both are varying, we trust the order the person gave. */
    snew(re->ind, re->nrepl);
    for (i = 0; i < re->nrepl; i++)
    {
        re->ind[i] = i;
    }

    if (re->type < ereENDSINGLE)
    {

        for (i = 0; i < re->nrepl; i++)
        {
            for (j = i+1; j < re->nrepl; j++)
            {
                if (re->q[re->type][re->ind[j]] < re->q[re->type][re->ind[i]])
                {
                    /* Unordered replicas are supposed to work, but there
                     * is still an issues somewhere.
                     * Note that at this point still re->ind[i]=i.
                     */
                    gmx_fatal(FARGS, "Replicas with indices %d < %d have %ss %g > %g, please order your replicas on increasing %s",
                              i, j,
                              erename[re->type],
                              re->q[re->type][i], re->q[re->type][j],
                              erename[re->type]);

                    k          = re->ind[i];
                    re->ind[i] = re->ind[j];
                    re->ind[j] = k;
                }
                else if (re->q[re->type][re->ind[j]] == re->q[re->type][re->ind[i]])
                {
                    gmx_fatal(FARGS, "Two replicas have identical %ss", erename[re->type]);
                }
            }
        }
    }

    /* keep track of all the swaps, starting with the initial placement. */
    snew(re->allswaps, re->nrepl);
    for (i = 0; i < re->nrepl; i++)
    {
        re->allswaps[i] = re->ind[i];
    }

    switch (re->type)
    {
        case ereTEMP:
            fprintf(fplog, "\nReplica exchange in temperature\n");
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5.1f", re->q[re->type][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        case ereLAMBDA:
            fprintf(fplog, "\nReplica exchange in lambda\n");
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %3d", (int)re->q[re->type][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        case ereTL:
            fprintf(fplog, "\nReplica exchange in temperature and lambda state\n");
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5.1f", re->q[ereTEMP][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            for (i = 0; i < re->nrepl; i++)
            {
                fprintf(fplog, " %5d", (int)re->q[ereLAMBDA][re->ind[i]]);
            }
            fprintf(fplog, "\n");
            break;
        default:
            gmx_incons("Unknown replica exchange quantity");
    }
    if (re->bNPT)
    {
        fprintf(fplog, "\nRepl  p");
        for (i = 0; i < re->nrepl; i++)
        {
            fprintf(fplog, " %5.2f", re->pres[re->ind[i]]);
        }

        for (i = 0; i < re->nrepl; i++)
        {
            if ((i > 0) && (re->pres[re->ind[i]] < re->pres[re->ind[i-1]]))
            {
                fprintf(fplog, "\nWARNING: The reference pressures decrease with increasing temperatures\n\n");
                fprintf(stderr, "\nWARNING: The reference pressures decrease with increasing temperatures\n\n");
            }
        }
    }
    re->nst = nst;
    if (replExParams.randomSeed == -1)
    {
        if (MASTERSIM(ms))
        {
            re->seed = static_cast<int>(gmx::makeRandomSeed());
        }
        else
        {
            re->seed = 0;
        }
        gmx_sumi_sim(1, &(re->seed), ms);
    }
    else
    {
        re->seed = replExParams.randomSeed;
    }
    fprintf(fplog, "\nReplica exchange interval: %d\n", re->nst);
    fprintf(fplog, "\nReplica random seed: %d\n", re->seed);

    re->nattempt[0] = 0;
    re->nattempt[1] = 0;

    snew(re->prob_sum, re->nrepl);
    snew(re->nexchange, re->nrepl);
    snew(re->nmoves, re->nrepl);
    for (i = 0; i < re->nrepl; i++)
    {
        snew(re->nmoves[i], re->nrepl);
    }
    fprintf(fplog, "Replica exchange information below: ex and x = exchange, pr = probability\n");

    /* generate space for the helper functions so we don't have to snew each time */

    snew(re->destinations, re->nrepl);
    snew(re->incycle, re->nrepl);
    snew(re->tmpswap, re->nrepl);
    snew(re->cyclic, re->nrepl);
    snew(re->order, re->nrepl);
    for (i = 0; i < re->nrepl; i++)
    {
        snew(re->cyclic[i], re->nrepl+1);
        snew(re->order[i], re->nrepl);
    }
    /* allocate space for the functions storing the data for the replicas */
    /* not all of these arrays needed in all cases, but they don't take
       up much space, since the max size is nrepl**2 */
    snew(re->prob, re->nrepl);
    snew(re->bEx, re->nrepl);
    snew(re->beta, re->nrepl);
    snew(re->Vol, re->nrepl);
    snew(re->Epot, re->nrepl);
    snew(re->de, re->nrepl);
    for (i = 0; i < re->nrepl; i++)
    {
        snew(re->de[i], re->nrepl);
    }
    re->nex = replExParams.numExchanges;
    return re;
}

static void exchange_reals(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, real *v, int n)
{
    real *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#if GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}


static void exchange_doubles(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, double *v, int n)
{
    double *buf;
    int     i;

    if (v)
    {
        snew(buf, n);
#if GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_rvecs(const gmx_multisim_t gmx_unused *ms, int gmx_unused b, rvec *v, int n)
{
    rvec *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#if GMX_MPI
        /*
           MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            copy_rvec(buf[i], v[i]);
        }
        sfree(buf);
    }
}

static void exchange_state(const gmx_multisim_t *ms, int b, t_state *state)
{
    /* When t_state changes, this code should be updated. */
    int ngtc, nnhpres;
    ngtc    = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    exchange_rvecs(ms, b, state->box, DIM);
    exchange_rvecs(ms, b, state->box_rel, DIM);
    exchange_rvecs(ms, b, state->boxv, DIM);
    exchange_reals(ms, b, &(state->veta), 1);
    exchange_reals(ms, b, &(state->vol0), 1);
    exchange_rvecs(ms, b, state->svir_prev, DIM);
    exchange_rvecs(ms, b, state->fvir_prev, DIM);
    exchange_rvecs(ms, b, state->pres_prev, DIM);
    exchange_doubles(ms, b, state->nosehoover_xi.data(), ngtc);
    exchange_doubles(ms, b, state->nosehoover_vxi.data(), ngtc);
    exchange_doubles(ms, b, state->nhpres_xi.data(), nnhpres);
    exchange_doubles(ms, b, state->nhpres_vxi.data(), nnhpres);
    exchange_doubles(ms, b, state->therm_integral.data(), state->ngtc);
    exchange_doubles(ms, b, &state->baros_integral, 1);
    exchange_rvecs(ms, b, as_rvec_array(state->x.data()), state->natoms);
    exchange_rvecs(ms, b, as_rvec_array(state->v.data()), state->natoms);
}

static void copy_state_serial(const t_state *src, t_state *dest)
{
    if (dest != src)
    {
        /* Currently the local state is always a pointer to the global
         * in serial, so we should never end up here.
         * TODO: Implement a (trivial) t_state copy once converted to C++.
         */
        GMX_RELEASE_ASSERT(false, "State copying is currently not implemented in replica exchange");
    }
}

static void scale_velocities(t_state *state, real fac)
{
    int i;

    if (as_rvec_array(state->v.data()))
    {
        for (i = 0; i < state->natoms; i++)
        {
            svmul(fac, state->v[i], state->v[i]);
        }
    }
}

static void print_transition_matrix(FILE *fplog, int n, int **nmoves, int *nattempt)
{
    int   i, j, ntot;
    float Tprint;

    ntot = nattempt[0] + nattempt[1];
    fprintf(fplog, "\n");
    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "    ");  /* put the title closer to the center */
    }
    fprintf(fplog, "Empirical Transition Matrix\n");

    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%8d", (i+1));
    }
    fprintf(fplog, "\n");

    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "Repl");
        for (j = 0; j < n; j++)
        {
            Tprint = 0.0;
            if (nmoves[i][j] > 0)
            {
                Tprint = nmoves[i][j]/(2.0*ntot);
            }
            fprintf(fplog, "%8.4f", Tprint);
        }
        fprintf(fplog, "%3d\n", i);
    }
}

static void print_ind(FILE *fplog, const char *leg, int n, int *ind, gmx_bool *bEx)
{
    int i;

    fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %c %2d", (bEx != nullptr && bEx[i]) ? 'x' : ' ', ind[i]);
    }
    fprintf(fplog, "\n");
}

static void print_allswitchind(FILE *fplog, int n, int *pind, int *allswaps, int *tmpswap)
{
    int i;

    for (i = 0; i < n; i++)
    {
        tmpswap[i] = allswaps[i];
    }
    for (i = 0; i < n; i++)
    {
        allswaps[i] = tmpswap[pind[i]];
    }

    fprintf(fplog, "\nAccepted Exchanges:   ");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%d ", pind[i]);
    }
    fprintf(fplog, "\n");

    /* the "Order After Exchange" is the state label corresponding to the configuration that
       started in state listed in order, i.e.

       3 0 1 2

       means that the:
       configuration starting in simulation 3 is now in simulation 0,
       configuration starting in simulation 0 is now in simulation 1,
       configuration starting in simulation 1 is now in simulation 2,
       configuration starting in simulation 2 is now in simulation 3
     */
    fprintf(fplog, "Order After Exchange: ");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%d ", allswaps[i]);
    }
    fprintf(fplog, "\n\n");
}

static void print_prob(FILE *fplog, const char *leg, int n, real *prob)
{
    int  i;
    char buf[8];

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf, "%4.2f", prob[i]);
            fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog, "     ");
        }
    }
    fprintf(fplog, "\n");
}

static void print_count(FILE *fplog, const char *leg, int n, int *count)
{
    int i;

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %4d", count[i]);
    }
    fprintf(fplog, "\n");
}

static real calc_delta(FILE *fplog, gmx_bool bPrint, struct gmx_repl_ex *re, int a, int b, int ap, int bp)
{

    real   ediff, dpV, delta = 0;
    real  *Epot = re->Epot;
    real  *Vol  = re->Vol;
    real **de   = re->de;
    real  *beta = re->beta;

    /* Two cases; we are permuted and not.  In all cases, setting ap = a and bp = b will reduce
       to the non permuted case */

    switch (re->type)
    {
        case ereTEMP:
            /*
             * Okabe et. al. Chem. Phys. Lett. 335 (2001) 435-439
             */
            ediff = Epot[b] - Epot[a];
            delta = -(beta[bp] - beta[ap])*ediff;
            break;
        case ereLAMBDA:
            /* two cases:  when we are permuted, and not.  */
            /* non-permuted:
               ediff =  E_new - E_old
                     =  [H_b(x_a) + H_a(x_b)] - [H_b(x_b) + H_a(x_a)]
                     =  [H_b(x_a) - H_a(x_a)] + [H_a(x_b) - H_b(x_b)]
                     =  de[b][a] + de[a][b] */

            /* permuted:
               ediff =  E_new - E_old
                     =  [H_bp(x_a) + H_ap(x_b)] - [H_bp(x_b) + H_ap(x_a)]
                     =  [H_bp(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a) + H_a(x_a) - H_ap(x_a)] + [H_ap(x_b) - H_b(x_b) + H_b(x_b) - H_bp(x_b)]
                     =  [H_bp(x_a) - H_a(x_a)] - [H_ap(x_a) - H_a(x_a)] + [H_ap(x_b) - H_b(x_b)] - H_bp(x_b) - H_b(x_b)]
                     =  (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b])    */
            /* but, in the current code implementation, we flip configurations, not indices . . .
               So let's examine that.
                     =  [H_b(x_ap) - H_a(x_a)] - [H_a(x_ap) - H_a(x_a)] + [H_a(x_bp) - H_b(x_b)] - H_b(x_bp) - H_b(x_b)]
                     =  [H_b(x_ap) - H_a(x_ap)]  + [H_a(x_bp) - H_b(x_pb)]
                     = (de[b][ap] - de[a][ap]) + (de[a][bp] - de[b][bp]
                     So, if we exchange b<=> bp and a<=> ap, we return to the same result.
                     So the simple solution is to flip the
                     position of perturbed and original indices in the tests.
             */

            ediff = (de[bp][a] - de[ap][a]) + (de[ap][b] - de[bp][b]);
            delta = ediff*beta[a]; /* assume all same temperature in this case */
            break;
        case ereTL:
            /* not permuted:  */
            /* delta =  reduced E_new - reduced E_old
                     =  [beta_b H_b(x_a) + beta_a H_a(x_b)] - [beta_b H_b(x_b) + beta_a H_a(x_a)]
                     =  [beta_b H_b(x_a) - beta_a H_a(x_a)] + [beta_a H_a(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + beta_b H_a(x_a) - beta_a H_a(x_a)] +
                        [beta_a dH_a(x_b) + beta_a H_b(x_b) - beta_b H_b(x_b)]
                     =  [beta_b dH_b(x_a) + [beta_a dH_a(x_b) +
                        beta_b (H_a(x_a) - H_b(x_b)]) - beta_a (H_a(x_a) - H_b(x_b))
                     =  beta_b dH_b(x_a) + beta_a dH_a(x_b) - (beta_b - beta_a)(H_b(x_b) - H_a(x_a) */
            /* delta = beta[b]*de[b][a] + beta[a]*de[a][b] - (beta[b] - beta[a])*(Epot[b] - Epot[a]; */
            /* permuted (big breath!) */
            /*   delta =  reduced E_new - reduced E_old
                     =  [beta_bp H_bp(x_a) + beta_ap H_ap(x_b)] - [beta_bp H_bp(x_b) + beta_ap H_ap(x_a)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                     =  [beta_bp H_bp(x_a) - beta_ap H_ap(x_a)] + [beta_ap H_ap(x_b) - beta_bp H_bp(x_b)]
                        - beta_pb H_a(x_a) + beta_ap H_a(x_a) + beta_pb H_a(x_a) - beta_ap H_a(x_a)
                        - beta_ap H_b(x_b) + beta_bp H_b(x_b) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [(beta_bp H_bp(x_a) - beta_bp H_a(x_a)) - (beta_ap H_ap(x_a) - beta_ap H_a(x_a))] +
                        [(beta_ap H_ap(x_b)  - beta_ap H_b(x_b)) - (beta_bp H_bp(x_b) - beta_bp H_b(x_b))]
             + beta_pb H_a(x_a) - beta_ap H_a(x_a) + beta_ap H_b(x_b) - beta_bp H_b(x_b)
                     =  [beta_bp (H_bp(x_a) - H_a(x_a)) - beta_ap (H_ap(x_a) - H_a(x_a))] +
                        [beta_ap (H_ap(x_b) - H_b(x_b)) - beta_bp (H_bp(x_b) - H_b(x_b))]
             + beta_pb (H_a(x_a) - H_b(x_b))  - beta_ap (H_a(x_a) - H_b(x_b))
                     =  ([beta_bp de[bp][a] - beta_ap de[ap][a]) + beta_ap de[ap][b]  - beta_bp de[bp][b])
             + (beta_pb-beta_ap)(H_a(x_a) - H_b(x_b))  */
            delta = beta[bp]*(de[bp][a] - de[bp][b]) + beta[ap]*(de[ap][b] - de[ap][a]) - (beta[bp]-beta[ap])*(Epot[b]-Epot[a]);
            break;
        default:
            gmx_incons("Unknown replica exchange quantity");
    }
    if (bPrint)
    {
        fprintf(fplog, "Repl %d <-> %d  dE_term = %10.3e (kT)\n", a, b, delta);
    }
    if (re->bNPT)
    {
        /* revist the calculation for 5.0.  Might be some improvements. */
        dpV = (beta[ap]*re->pres[ap]-beta[bp]*re->pres[bp])*(Vol[b]-Vol[a])/PRESFAC;
        if (bPrint)
        {
            fprintf(fplog, "  dpV = %10.3e  d = %10.3e\n", dpV, delta + dpV);
        }
        delta += dpV;
    }
    return delta;
}

static void
test_for_replica_exchange(FILE                 *fplog,
                          const gmx_multisim_t *ms,
                          struct gmx_repl_ex   *re,
                          const gmx_enerdata_t *enerd,
                          real                  vol,
                          gmx_int64_t           step,
                          real                  time)
{
    int                                  m, i, j, a, b, ap, bp, i0, i1, tmp;
    real                                 delta = 0;
    gmx_bool                             bPrint, bMultiEx;
    gmx_bool                            *bEx      = re->bEx;
    real                                *prob     = re->prob;
    int                                 *pind     = re->destinations; /* permuted index */
    gmx_bool                             bEpot    = FALSE;
    gmx_bool                             bDLambda = FALSE;
    gmx_bool                             bVol     = FALSE;
    gmx::ThreeFry2x64<64>                rng(re->seed, gmx::RandomDomain::ReplicaExchange);
    gmx::UniformRealDistribution<real>   uniformRealDist;
    gmx::UniformIntDistribution<int>     uniformNreplDist(0, re->nrepl-1);

    bMultiEx = (re->nex > 1);  /* multiple exchanges at each state */
    fprintf(fplog, "Replica exchange at step %" GMX_PRId64 " time %.5f\n", step, time);

    if (re->bNPT)
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->Vol[i] = 0;
        }
        bVol               = TRUE;
        re->Vol[re->repl]  = vol;
    }
    if ((re->type == ereTEMP || re->type == ereTL))
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->Epot[i] = 0;
        }
        bEpot              = TRUE;
        re->Epot[re->repl] = enerd->term[F_EPOT];
        /* temperatures of different states*/
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->q[ereTEMP][i]*BOLTZ);
        }
    }
    else
    {
        for (i = 0; i < re->nrepl; i++)
        {
            re->beta[i] = 1.0/(re->temp*BOLTZ);  /* we have a single temperature */
        }
    }
    if (re->type == ereLAMBDA || re->type == ereTL)
    {
        bDLambda = TRUE;
        /* lambda differences. */
        /* de[i][j] is the energy of the jth simulation in the ith Hamiltonian
           minus the energy of the jth simulation in the jth Hamiltonian */
        for (i = 0; i < re->nrepl; i++)
        {
            for (j = 0; j < re->nrepl; j++)
            {
                re->de[i][j] = 0;
            }
        }
        for (i = 0; i < re->nrepl; i++)
        {
            re->de[i][re->repl] = (enerd->enerpart_lambda[(int)re->q[ereLAMBDA][i]+1]-enerd->enerpart_lambda[0]);
        }
    }

    /* now actually do the communication */
    if (bVol)
    {
        gmx_sum_sim(re->nrepl, re->Vol, ms);
    }
    if (bEpot)
    {
        gmx_sum_sim(re->nrepl, re->Epot, ms);
    }
    if (bDLambda)
    {
        for (i = 0; i < re->nrepl; i++)
        {
            gmx_sum_sim(re->nrepl, re->de[i], ms);
        }
    }

    /* make a duplicate set of indices for shuffling */
    for (i = 0; i < re->nrepl; i++)
    {
        pind[i] = re->ind[i];
    }

    rng.restart( step, 0 );

    if (bMultiEx)
    {
        /* multiple random switch exchange */
        int nself = 0;


        for (i = 0; i < re->nex + nself; i++)
        {
            // For now this is superfluous, but just in case we ever add more
            // calls in different branches it is safer to always reset the distribution.
            uniformNreplDist.reset();

            /* randomly select a pair  */
            /* in theory, could reduce this by identifying only which switches had a nonneglibible
               probability of occurring (log p > -100) and only operate on those switches */
            /* find out which state it is from, and what label that state currently has. Likely
               more work that useful. */
            i0 = uniformNreplDist(rng);
            i1 = uniformNreplDist(rng);
            if (i0 == i1)
            {
                nself++;
                continue;  /* self-exchange, back up and do it again */
            }

            a  = re->ind[i0]; /* what are the indices of these states? */
            b  = re->ind[i1];
            ap = pind[i0];
            bp = pind[i1];

            bPrint = FALSE; /* too noisy */
            /* calculate the energy difference */
            /* if the code changes to flip the STATES, rather than the configurations,
               use the commented version of the code */
            /* delta = calc_delta(fplog,bPrint,re,a,b,ap,bp); */
            delta = calc_delta(fplog, bPrint, re, ap, bp, a, b);

            /* we actually only use the first space in the prob and bEx array,
               since there are actually many switches between pairs. */

            if (delta <= 0)
            {
                /* accepted */
                prob[0] = 1;
                bEx[0]  = TRUE;
            }
            else
            {
                if (delta > PROBABILITYCUTOFF)
                {
                    prob[0] = 0;
                }
                else
                {
                    prob[0] = exp(-delta);
                }
                // roll a number to determine if accepted. For now it is superfluous to
                // reset, but just in case we ever add more calls in different branches
                // it is safer to always reset the distribution.
                uniformRealDist.reset();
                bEx[0] = uniformRealDist(rng) < prob[0];
            }
            re->prob_sum[0] += prob[0];

            if (bEx[0])
            {
                /* swap the states */
                tmp      = pind[i0];
                pind[i0] = pind[i1];
                pind[i1] = tmp;
            }
        }
        re->nattempt[0]++;  /* keep track of total permutation trials here */
        print_allswitchind(fplog, re->nrepl, pind, re->allswaps, re->tmpswap);
    }
    else
    {
        /* standard nearest neighbor replica exchange */

        m = (step / re->nst) % 2;
        for (i = 1; i < re->nrepl; i++)
        {
            a = re->ind[i-1];
            b = re->ind[i];

            bPrint = (re->repl == a || re->repl == b);
            if (i % 2 == m)
            {
                delta = calc_delta(fplog, bPrint, re, a, b, a, b);
                if (delta <= 0)
                {
                    /* accepted */
                    prob[i] = 1;
                    bEx[i]  = TRUE;
                }
                else
                {
                    if (delta > PROBABILITYCUTOFF)
                    {
                        prob[i] = 0;
                    }
                    else
                    {
                        prob[i] = exp(-delta);
                    }
                    // roll a number to determine if accepted. For now it is superfluous to
                    // reset, but just in case we ever add more calls in different branches
                    // it is safer to always reset the distribution.
                    uniformRealDist.reset();
                    bEx[i] = uniformRealDist(rng) < prob[i];
                }
                re->prob_sum[i] += prob[i];

                if (bEx[i])
                {
                    /* swap these two */
                    tmp       = pind[i-1];
                    pind[i-1] = pind[i];
                    pind[i]   = tmp;
                    re->nexchange[i]++;  /* statistics for back compatibility */
                }
            }
            else
            {
                prob[i] = -1;
                bEx[i]  = FALSE;
            }
        }
        /* print some statistics */
        print_ind(fplog, "ex", re->nrepl, re->ind, bEx);
        print_prob(fplog, "pr", re->nrepl, prob);
        fprintf(fplog, "\n");
        re->nattempt[m]++;
    }

    /* record which moves were made and accepted */
    for (i = 0; i < re->nrepl; i++)
    {
        re->nmoves[re->ind[i]][pind[i]] += 1;
        re->nmoves[pind[i]][re->ind[i]] += 1;
    }
    fflush(fplog); /* make sure we can see what the last exchange was */
}

static void
cyclic_decomposition(const int *destinations,
                     int      **cyclic,
                     gmx_bool  *incycle,
                     const int  nrepl,
                     int       *nswap)
{

    int i, j, c, p;
    int maxlen = 1;
    for (i = 0; i < nrepl; i++)
    {
        incycle[i] = FALSE;
    }
    for (i = 0; i < nrepl; i++)  /* one cycle for each replica */
    {
        if (incycle[i])
        {
            cyclic[i][0] = -1;
            continue;
        }
        cyclic[i][0] = i;
        incycle[i]   = TRUE;
        c            = 1;
        p            = i;
        for (j = 0; j < nrepl; j++) /* potentially all cycles are part, but we will break first */
        {
            p = destinations[p];    /* start permuting */
            if (p == i)
            {
                cyclic[i][c] = -1;
                if (c > maxlen)
                {
                    maxlen = c;
                }
                break; /* we've reached the original element, the cycle is complete, and we marked the end. */
            }
            else
            {
                cyclic[i][c] = p;  /* each permutation gives a new member of the cycle */
                incycle[p]   = TRUE;
                c++;
            }
        }
    }
    *nswap = maxlen - 1;

    if (debug)
    {
        for (i = 0; i < nrepl; i++)
        {
            fprintf(debug, "Cycle %d:", i);
            for (j = 0; j < nrepl; j++)
            {
                if (cyclic[i][j] < 0)
                {
                    break;
                }
                fprintf(debug, "%2d", cyclic[i][j]);
            }
            fprintf(debug, "\n");
        }
        fflush(debug);
    }
}

static void
compute_exchange_order(int     **cyclic,
                       int     **order,
                       const int nrepl,
                       const int maxswap)
{
    int i, j;

    for (j = 0; j < maxswap; j++)
    {
        for (i = 0; i < nrepl; i++)
        {
            if (cyclic[i][j+1] >= 0)
            {
                order[cyclic[i][j+1]][j] = cyclic[i][j];
                order[cyclic[i][j]][j]   = cyclic[i][j+1];
            }
        }
        for (i = 0; i < nrepl; i++)
        {
            if (order[i][j] < 0)
            {
                order[i][j] = i; /* if it's not exchanging, it should stay this round*/
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Replica Exchange Order\n");
        for (i = 0; i < nrepl; i++)
        {
            fprintf(debug, "Replica %d:", i);
            for (j = 0; j < maxswap; j++)
            {
                if (order[i][j] < 0)
                {
                    break;
                }
                fprintf(debug, "%2d", order[i][j]);
            }
            fprintf(debug, "\n");
        }
        fflush(debug);
    }
}

static void
prepare_to_do_exchange(struct gmx_repl_ex *re,
                       const int           replica_id,
                       int                *maxswap,
                       gmx_bool           *bThisReplicaExchanged)
{
    int i, j;
    /* Hold the cyclic decomposition of the (multiple) replica
     * exchange. */
    gmx_bool bAnyReplicaExchanged = FALSE;
    *bThisReplicaExchanged = FALSE;

    for (i = 0; i < re->nrepl; i++)
    {
        if (re->destinations[i] != re->ind[i])
        {
            /* only mark as exchanged if the index has been shuffled */
            bAnyReplicaExchanged = TRUE;
            break;
        }
    }
    if (bAnyReplicaExchanged)
    {
        /* reinitialize the placeholder arrays */
        for (i = 0; i < re->nrepl; i++)
        {
            for (j = 0; j < re->nrepl; j++)
            {
                re->cyclic[i][j] = -1;
                re->order[i][j]  = -1;
            }
        }

        /* Identify the cyclic decomposition of the permutation (very
         * fast if neighbor replica exchange). */
        cyclic_decomposition(re->destinations, re->cyclic, re->incycle, re->nrepl, maxswap);

        /* Now translate the decomposition into a replica exchange
         * order at each step. */
        compute_exchange_order(re->cyclic, re->order, re->nrepl, *maxswap);

        /* Did this replica do any exchange at any point? */
        for (j = 0; j < *maxswap; j++)
        {
            if (replica_id != re->order[replica_id][j])
            {
                *bThisReplicaExchanged = TRUE;
                break;
            }
        }
    }
}

gmx_bool replica_exchange(FILE *fplog, const t_commrec *cr, struct gmx_repl_ex *re,
                          t_state *state, const gmx_enerdata_t *enerd,
                          t_state *state_local, gmx_int64_t step, real time)
{
    int j;
    int replica_id = 0;
    int exchange_partner;
    int maxswap = 0;
    /* Number of rounds of exchanges needed to deal with any multiple
     * exchanges. */
    /* Where each replica ends up after the exchange attempt(s). */
    /* The order in which multiple exchanges will occur. */
    gmx_bool bThisReplicaExchanged = FALSE;

    if (MASTER(cr))
    {
        replica_id  = re->repl;
        test_for_replica_exchange(fplog, cr->ms, re, enerd, det(state_local->box), step, time);
        prepare_to_do_exchange(re, replica_id, &maxswap, &bThisReplicaExchanged);
    }
    /* Do intra-simulation broadcast so all processors belonging to
     * each simulation know whether they need to participate in
     * collecting the state. Otherwise, they might as well get on with
     * the next thing to do. */
    if (DOMAINDECOMP(cr))
    {
#if GMX_MPI
        MPI_Bcast(&bThisReplicaExchanged, sizeof(gmx_bool), MPI_BYTE, MASTERRANK(cr),
                  cr->mpi_comm_mygroup);
#endif
    }

    if (bThisReplicaExchanged)
    {
        /* Exchange the states */
        /* Collect the global state on the master node */
        if (DOMAINDECOMP(cr))
        {
            dd_collect_state(cr->dd, state_local, state);
        }
        else
        {
            copy_state_serial(state_local, state);
        }

        if (MASTER(cr))
        {
            /* There will be only one swap cycle with standard replica
             * exchange, but there may be multiple swap cycles if we
             * allow multiple swaps. */

            for (j = 0; j < maxswap; j++)
            {
                exchange_partner = re->order[replica_id][j];

                if (exchange_partner != replica_id)
                {
                    /* Exchange the global states between the master nodes */
                    if (debug)
                    {
                        fprintf(debug, "Exchanging %d with %d\n", replica_id, exchange_partner);
                    }
                    exchange_state(cr->ms, exchange_partner, state);
                }
            }
            /* For temperature-type replica exchange, we need to scale
             * the velocities. */
            if (re->type == ereTEMP || re->type == ereTL)
            {
                scale_velocities(state, sqrt(re->q[ereTEMP][replica_id]/re->q[ereTEMP][re->destinations[replica_id]]));
            }

        }

        /* With domain decomposition the global state is distributed later */
        if (!DOMAINDECOMP(cr))
        {
            /* Copy the global state to the local state data structure */
            copy_state_serial(state, state_local);
        }
    }

    return bThisReplicaExchanged;
}

void print_replica_exchange_statistics(FILE *fplog, struct gmx_repl_ex *re)
{
    int  i;

    fprintf(fplog, "\nReplica exchange statistics\n");

    if (re->nex == 0)
    {
        fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
                re->nattempt[0]+re->nattempt[1], re->nattempt[1], re->nattempt[0]);

        fprintf(fplog, "Repl  average probabilities:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "Repl  number of exchanges:\n");
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_count(fplog, "", re->nrepl, re->nexchange);

        fprintf(fplog, "Repl  average number of exchanges:\n");
        for (i = 1; i < re->nrepl; i++)
        {
            if (re->nattempt[i%2] == 0)
            {
                re->prob[i] = 0;
            }
            else
            {
                re->prob[i] =  ((real)re->nexchange[i])/re->nattempt[i%2];
            }
        }
        print_ind(fplog, "", re->nrepl, re->ind, nullptr);
        print_prob(fplog, "", re->nrepl, re->prob);

        fprintf(fplog, "\n");
    }
    /* print the transition matrix */
    print_transition_matrix(fplog, re->nrepl, re->nmoves, re->nattempt);
}
