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

#include "grompp.h"

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include <sys/types.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_maxwell_velocities.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/sortwater.h"
#include "gromacs/gmxpreprocess/tomorse.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/imd/imd.h"
#include "gromacs/legacyheaders/calcgrid.h"
#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/perf_est.h"
#include "gromacs/legacyheaders/splitter.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/warninp.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

static int rm_interactions(int ifunc, int nrmols, t_molinfo mols[])
{
    int  i, n;

    n = 0;
    /* For all the molecule types */
    for (i = 0; i < nrmols; i++)
    {
        n += mols[i].plist[ifunc].nr;
        mols[i].plist[ifunc].nr = 0;
    }
    return n;
}

static int check_atom_names(const char *fn1, const char *fn2,
                            gmx_mtop_t *mtop, t_atoms *at)
{
    int      mb, m, i, j, nmismatch;
    t_atoms *tat;
#define MAXMISMATCH 20

    if (mtop->natoms != at->nr)
    {
        gmx_incons("comparing atom names");
    }

    nmismatch = 0;
    i         = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        tat = &mtop->moltype[mtop->molblock[mb].type].atoms;
        for (m = 0; m < mtop->molblock[mb].nmol; m++)
        {
            for (j = 0; j < tat->nr; j++)
            {
                if (strcmp( *(tat->atomname[j]), *(at->atomname[i]) ) != 0)
                {
                    if (nmismatch < MAXMISMATCH)
                    {
                        fprintf(stderr,
                                "Warning: atom name %d in %s and %s does not match (%s - %s)\n",
                                i+1, fn1, fn2, *(tat->atomname[j]), *(at->atomname[i]));
                    }
                    else if (nmismatch == MAXMISMATCH)
                    {
                        fprintf(stderr, "(more than %d non-matching atom names)\n", MAXMISMATCH);
                    }
                    nmismatch++;
                }
                i++;
            }
        }
    }

    return nmismatch;
}

static void check_eg_vs_cg(gmx_mtop_t *mtop)
{
    int            astart, mb, m, cg, j, firstj;
    unsigned char  firsteg, eg;
    gmx_moltype_t *molt;

    /* Go through all the charge groups and make sure all their
     * atoms are in the same energy group.
     */

    astart = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        for (m = 0; m < mtop->molblock[mb].nmol; m++)
        {
            for (cg = 0; cg < molt->cgs.nr; cg++)
            {
                /* Get the energy group of the first atom in this charge group */
                firstj  = astart + molt->cgs.index[cg];
                firsteg = ggrpnr(&mtop->groups, egcENER, firstj);
                for (j = molt->cgs.index[cg]+1; j < molt->cgs.index[cg+1]; j++)
                {
                    eg = ggrpnr(&mtop->groups, egcENER, astart+j);
                    if (eg != firsteg)
                    {
                        gmx_fatal(FARGS, "atoms %d and %d in charge group %d of molecule type '%s' are in different energy groups",
                                  firstj+1, astart+j+1, cg+1, *molt->name);
                    }
                }
            }
            astart += molt->atoms.nr;
        }
    }
}

static void check_cg_sizes(const char *topfn, t_block *cgs, warninp_t wi)
{
    int  maxsize, cg;
    char warn_buf[STRLEN];

    maxsize = 0;
    for (cg = 0; cg < cgs->nr; cg++)
    {
        maxsize = max(maxsize, cgs->index[cg+1]-cgs->index[cg]);
    }

    if (maxsize > MAX_CHARGEGROUP_SIZE)
    {
        gmx_fatal(FARGS, "The largest charge group contains %d atoms. The maximum is %d.", maxsize, MAX_CHARGEGROUP_SIZE);
    }
    else if (maxsize > 10)
    {
        set_warning_line(wi, topfn, -1);
        sprintf(warn_buf,
                "The largest charge group contains %d atoms.\n"
                "Since atoms only see each other when the centers of geometry of the charge groups they belong to are within the cut-off distance, too large charge groups can lead to serious cut-off artifacts.\n"
                "For efficiency and accuracy, charge group should consist of a few atoms.\n"
                "For all-atom force fields use: CH3, CH2, CH, NH2, NH, OH, CO2, CO, etc.",
                maxsize);
        warning_note(wi, warn_buf);
    }
}

static void check_bonds_timestep(gmx_mtop_t *mtop, double dt, warninp_t wi)
{
    /* This check is not intended to ensure accurate integration,
     * rather it is to signal mistakes in the mdp settings.
     * A common mistake is to forget to turn on constraints
     * for MD after energy minimization with flexible bonds.
     * This check can also detect too large time steps for flexible water
     * models, but such errors will often be masked by the constraints
     * mdp options, which turns flexible water into water with bond constraints,
     * but without an angle constraint. Unfortunately such incorrect use
     * of water models can not easily be detected without checking
     * for specific model names.
     *
     * The stability limit of leap-frog or velocity verlet is 4.44 steps
     * per oscillational period.
     * But accurate bonds distributions are lost far before that limit.
     * To allow relatively common schemes (although not common with Gromacs)
     * of dt=1 fs without constraints and dt=2 fs with only H-bond constraints
     * we set the note limit to 10.
     */
    int            min_steps_warn = 5;
    int            min_steps_note = 10;
    t_iparams     *ip;
    int            molt;
    gmx_moltype_t *moltype, *w_moltype;
    t_atom        *atom;
    t_ilist       *ilist, *ilb, *ilc, *ils;
    int            ftype;
    int            i, a1, a2, w_a1, w_a2, j;
    real           twopi2, limit2, fc, re, m1, m2, period2, w_period2;
    gmx_bool       bFound, bWater, bWarn;
    char           warn_buf[STRLEN];

    ip = mtop->ffparams.iparams;

    twopi2 = sqr(2*M_PI);

    limit2 = sqr(min_steps_note*dt);

    w_a1      = w_a2 = -1;
    w_period2 = -1.0;

    w_moltype = NULL;
    for (molt = 0; molt < mtop->nmoltype; molt++)
    {
        moltype = &mtop->moltype[molt];
        atom    = moltype->atoms.atom;
        ilist   = moltype->ilist;
        ilc     = &ilist[F_CONSTR];
        ils     = &ilist[F_SETTLE];
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (!(ftype == F_BONDS || ftype == F_G96BONDS || ftype == F_HARMONIC))
            {
                continue;
            }

            ilb = &ilist[ftype];
            for (i = 0; i < ilb->nr; i += 3)
            {
                fc = ip[ilb->iatoms[i]].harmonic.krA;
                re = ip[ilb->iatoms[i]].harmonic.rA;
                if (ftype == F_G96BONDS)
                {
                    /* Convert squared sqaure fc to harmonic fc */
                    fc = 2*fc*re;
                }
                a1 = ilb->iatoms[i+1];
                a2 = ilb->iatoms[i+2];
                m1 = atom[a1].m;
                m2 = atom[a2].m;
                if (fc > 0 && m1 > 0 && m2 > 0)
                {
                    period2 = twopi2*m1*m2/((m1 + m2)*fc);
                }
                else
                {
                    period2 = GMX_FLOAT_MAX;
                }
                if (debug)
                {
                    fprintf(debug, "fc %g m1 %g m2 %g period %g\n",
                            fc, m1, m2, sqrt(period2));
                }
                if (period2 < limit2)
                {
                    bFound = FALSE;
                    for (j = 0; j < ilc->nr; j += 3)
                    {
                        if ((ilc->iatoms[j+1] == a1 && ilc->iatoms[j+2] == a2) ||
                            (ilc->iatoms[j+1] == a2 && ilc->iatoms[j+2] == a1))
                        {
                            bFound = TRUE;
                        }
                    }
                    for (j = 0; j < ils->nr; j += 4)
                    {
                        if ((a1 == ils->iatoms[j+1] || a1 == ils->iatoms[j+2] || a1 == ils->iatoms[j+3]) &&
                            (a2 == ils->iatoms[j+1] || a2 == ils->iatoms[j+2] || a2 == ils->iatoms[j+3]))
                        {
                            bFound = TRUE;
                        }
                    }
                    if (!bFound &&
                        (w_moltype == NULL || period2 < w_period2))
                    {
                        w_moltype = moltype;
                        w_a1      = a1;
                        w_a2      = a2;
                        w_period2 = period2;
                    }
                }
            }
        }
    }

    if (w_moltype != NULL)
    {
        bWarn = (w_period2 < sqr(min_steps_warn*dt));
        /* A check that would recognize most water models */
        bWater = ((*w_moltype->atoms.atomname[0])[0] == 'O' &&
                  w_moltype->atoms.nr <= 5);
        sprintf(warn_buf, "The bond in molecule-type %s between atoms %d %s and %d %s has an estimated oscillational period of %.1e ps, which is less than %d times the time step of %.1e ps.\n"
                "%s",
                *w_moltype->name,
                w_a1+1, *w_moltype->atoms.atomname[w_a1],
                w_a2+1, *w_moltype->atoms.atomname[w_a2],
                sqrt(w_period2), bWarn ? min_steps_warn : min_steps_note, dt,
                bWater ?
                "Maybe you asked for fexible water." :
                "Maybe you forgot to change the constraints mdp option.");
        if (bWarn)
        {
            warning(wi, warn_buf);
        }
        else
        {
            warning_note(wi, warn_buf);
        }
    }
}

static void check_vel(gmx_mtop_t *mtop, rvec v[])
{
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     a;

    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &a, &atom))
    {
        if (atom->ptype == eptShell ||
            atom->ptype == eptBond  ||
            atom->ptype == eptVSite)
        {
            clear_rvec(v[a]);
        }
    }
}

static void check_shells_inputrec(gmx_mtop_t *mtop,
                                  t_inputrec *ir,
                                  warninp_t   wi)
{
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     a, nshells = 0;
    char                    warn_buf[STRLEN];

    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &a, &atom))
    {
        if (atom->ptype == eptShell ||
            atom->ptype == eptBond)
        {
            nshells++;
        }
    }
    if (IR_TWINRANGE(*ir) && (nshells > 0))
    {
        snprintf(warn_buf, STRLEN,
                 "The combination of using shells and a twin-range cut-off is not supported");
        warning_error(wi, warn_buf);
    }
    if ((nshells > 0) && (ir->nstcalcenergy != 1))
    {
        set_warning_line(wi, "unknown", -1);
        snprintf(warn_buf, STRLEN,
                 "There are %d shells, changing nstcalcenergy from %d to 1",
                 nshells, ir->nstcalcenergy);
        ir->nstcalcenergy = 1;
        warning(wi, warn_buf);
    }
}

/* TODO Decide whether this function can be consolidated with
 * gmx_mtop_ftype_count */
static gmx_bool nint_ftype(gmx_mtop_t *mtop, t_molinfo *mi, int ftype)
{
    int nint, mb;

    nint = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        nint += mtop->molblock[mb].nmol*mi[mtop->molblock[mb].type].plist[ftype].nr;
    }

    return nint;
}

/* This routine reorders the molecule type array
 * in the order of use in the molblocks,
 * unused molecule types are deleted.
 */
static void renumber_moltypes(gmx_mtop_t *sys,
                              int *nmolinfo, t_molinfo **molinfo)
{
    int       *order, norder, i;
    int        mb, mi;
    t_molinfo *minew;

    snew(order, *nmolinfo);
    norder = 0;
    for (mb = 0; mb < sys->nmolblock; mb++)
    {
        for (i = 0; i < norder; i++)
        {
            if (order[i] == sys->molblock[mb].type)
            {
                break;
            }
        }
        if (i == norder)
        {
            /* This type did not occur yet, add it */
            order[norder] = sys->molblock[mb].type;
            /* Renumber the moltype in the topology */
            norder++;
        }
        sys->molblock[mb].type = i;
    }

    /* We still need to reorder the molinfo structs */
    snew(minew, norder);
    for (mi = 0; mi < *nmolinfo; mi++)
    {
        for (i = 0; i < norder; i++)
        {
            if (order[i] == mi)
            {
                break;
            }
        }
        if (i == norder)
        {
            done_mi(&(*molinfo)[mi]);
        }
        else
        {
            minew[i] = (*molinfo)[mi];
        }
    }
    sfree(*molinfo);

    *nmolinfo = norder;
    *molinfo  = minew;
}

static void molinfo2mtop(int nmi, t_molinfo *mi, gmx_mtop_t *mtop)
{
    int            m;
    gmx_moltype_t *molt;

    mtop->nmoltype = nmi;
    snew(mtop->moltype, nmi);
    for (m = 0; m < nmi; m++)
    {
        molt        = &mtop->moltype[m];
        molt->name  = mi[m].name;
        molt->atoms = mi[m].atoms;
        /* ilists are copied later */
        molt->cgs   = mi[m].cgs;
        molt->excls = mi[m].excls;
    }
}

static void
new_status(const char *topfile, const char *topppfile, const char *confin,
           t_gromppopts *opts, t_inputrec *ir, gmx_bool bZero,
           gmx_bool bGenVel, gmx_bool bVerbose, t_state *state,
           gpp_atomtype_t atype, gmx_mtop_t *sys,
           int *nmi, t_molinfo **mi, t_molinfo **intermolecular_interactions,
           t_params plist[],
           int *comb, double *reppow, real *fudgeQQ,
           gmx_bool bMorse,
           warninp_t wi)
{
    t_molinfo      *molinfo = NULL;
    int             nmolblock;
    gmx_molblock_t *molblock, *molbs;
    t_atoms        *confat;
    int             mb, i, nrmols, nmismatch;
    char            buf[STRLEN];
    gmx_bool        bGB = FALSE;
    char            warn_buf[STRLEN];

    init_mtop(sys);

    /* Set gmx_boolean for GB */
    if (ir->implicit_solvent)
    {
        bGB = TRUE;
    }

    /* TOPOLOGY processing */
    sys->name = do_top(bVerbose, topfile, topppfile, opts, bZero, &(sys->symtab),
                       plist, comb, reppow, fudgeQQ,
                       atype, &nrmols, &molinfo, intermolecular_interactions,
                       ir,
                       &nmolblock, &molblock, bGB,
                       wi);

    sys->nmolblock = 0;
    snew(sys->molblock, nmolblock);

    sys->natoms = 0;
    for (mb = 0; mb < nmolblock; mb++)
    {
        if (sys->nmolblock > 0 &&
            molblock[mb].type == sys->molblock[sys->nmolblock-1].type)
        {
            /* Merge consecutive blocks with the same molecule type */
            sys->molblock[sys->nmolblock-1].nmol += molblock[mb].nmol;
            sys->natoms += molblock[mb].nmol*sys->molblock[sys->nmolblock-1].natoms_mol;
        }
        else if (molblock[mb].nmol > 0)
        {
            /* Add a new molblock to the topology */
            molbs             = &sys->molblock[sys->nmolblock];
            *molbs            = molblock[mb];
            molbs->natoms_mol = molinfo[molbs->type].atoms.nr;
            molbs->nposres_xA = 0;
            molbs->nposres_xB = 0;
            sys->natoms      += molbs->nmol*molbs->natoms_mol;
            sys->nmolblock++;
        }
    }
    if (sys->nmolblock == 0)
    {
        gmx_fatal(FARGS, "No molecules were defined in the system");
    }

    renumber_moltypes(sys, &nrmols, &molinfo);

    if (bMorse)
    {
        convert_harmonics(nrmols, molinfo, atype);
    }

    if (ir->eDisre == edrNone)
    {
        i = rm_interactions(F_DISRES, nrmols, molinfo);
        if (i > 0)
        {
            set_warning_line(wi, "unknown", -1);
            sprintf(warn_buf, "disre = no, removed %d distance restraints", i);
            warning_note(wi, warn_buf);
        }
    }
    if (opts->bOrire == FALSE)
    {
        i = rm_interactions(F_ORIRES, nrmols, molinfo);
        if (i > 0)
        {
            set_warning_line(wi, "unknown", -1);
            sprintf(warn_buf, "orire = no, removed %d orientation restraints", i);
            warning_note(wi, warn_buf);
        }
    }

    /* Copy structures from msys to sys */
    molinfo2mtop(nrmols, molinfo, sys);

    gmx_mtop_finalize(sys);

    /* COORDINATE file processing */
    if (bVerbose)
    {
        fprintf(stderr, "processing coordinates...\n");
    }

    get_stx_coordnum(confin, &state->natoms);
    if (state->natoms != sys->natoms)
    {
        gmx_fatal(FARGS, "number of coordinates in coordinate file (%s, %d)\n"
                  "             does not match topology (%s, %d)",
                  confin, state->natoms, topfile, sys->natoms);
    }
    else
    {
        /* make space for coordinates and velocities */
        char title[STRLEN];
        snew(confat, 1);
        init_t_atoms(confat, state->natoms, FALSE);
        init_state(state, state->natoms, 0, 0, 0, 0);
        read_stx_conf(confin, title, confat, state->x, state->v, NULL, state->box);
        /* This call fixes the box shape for runs with pressure scaling */
        set_box_rel(ir, state);

        nmismatch = check_atom_names(topfile, confin, sys, confat);
        free_t_atoms(confat, TRUE);
        sfree(confat);

        if (nmismatch)
        {
            sprintf(buf, "%d non-matching atom name%s\n"
                    "atom names from %s will be used\n"
                    "atom names from %s will be ignored\n",
                    nmismatch, (nmismatch == 1) ? "" : "s", topfile, confin);
            warning(wi, buf);
        }

        /* Do more checks, mostly related to constraints */
        if (bVerbose)
        {
            fprintf(stderr, "double-checking input for internal consistency...\n");
        }
        {
            int bHasNormalConstraints = 0 < (nint_ftype(sys, molinfo, F_CONSTR) +
                                             nint_ftype(sys, molinfo, F_CONSTRNC));
            int bHasAnyConstraints = bHasNormalConstraints || 0 < nint_ftype(sys, molinfo, F_SETTLE);
            double_check(ir, state->box,
                         bHasNormalConstraints,
                         bHasAnyConstraints,
                         wi);
        }
    }

    if (bGenVel)
    {
        real                   *mass;
        gmx_mtop_atomloop_all_t aloop;
        t_atom                 *atom;
        unsigned int            useed;

        snew(mass, state->natoms);
        aloop = gmx_mtop_atomloop_all_init(sys);
        while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
        {
            mass[i] = atom->m;
        }

        useed = opts->seed;
        if (opts->seed == -1)
        {
            useed = (int)gmx_rng_make_seed();
            fprintf(stderr, "Setting gen_seed to %u\n", useed);
        }
        maxwell_speed(opts->tempi, useed, sys, state->v);

        stop_cm(stdout, state->natoms, mass, state->x, state->v);
        sfree(mass);
    }

    *nmi = nrmols;
    *mi  = molinfo;
}

static void copy_state(const char *slog, t_trxframe *fr,
                       gmx_bool bReadVel, t_state *state,
                       double *use_time)
{
    int i;

    if (fr->not_ok & FRAME_NOT_OK)
    {
        gmx_fatal(FARGS, "Can not start from an incomplete frame");
    }
    if (!fr->bX)
    {
        gmx_fatal(FARGS, "Did not find a frame with coordinates in file %s",
                  slog);
    }

    for (i = 0; i < state->natoms; i++)
    {
        copy_rvec(fr->x[i], state->x[i]);
    }
    if (bReadVel)
    {
        if (!fr->bV)
        {
            gmx_incons("Trajecory frame unexpectedly does not contain velocities");
        }
        for (i = 0; i < state->natoms; i++)
        {
            copy_rvec(fr->v[i], state->v[i]);
        }
    }
    if (fr->bBox)
    {
        copy_mat(fr->box, state->box);
    }

    *use_time = fr->time;
}

static void cont_status(const char *slog, const char *ener,
                        gmx_bool bNeedVel, gmx_bool bGenVel, real fr_time,
                        t_inputrec *ir, t_state *state,
                        gmx_mtop_t *sys,
                        const output_env_t oenv)
/* If fr_time == -1 read the last frame available which is complete */
{
    gmx_bool     bReadVel;
    t_trxframe   fr;
    t_trxstatus *fp;
    int          i;
    double       use_time;

    bReadVel = (bNeedVel && !bGenVel);

    fprintf(stderr,
            "Reading Coordinates%s and Box size from old trajectory\n",
            bReadVel ? ", Velocities" : "");
    if (fr_time == -1)
    {
        fprintf(stderr, "Will read whole trajectory\n");
    }
    else
    {
        fprintf(stderr, "Will read till time %g\n", fr_time);
    }
    if (!bReadVel)
    {
        if (bGenVel)
        {
            fprintf(stderr, "Velocities generated: "
                    "ignoring velocities in input trajectory\n");
        }
        read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X);
    }
    else
    {
        read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X | TRX_NEED_V);

        if (!fr.bV)
        {
            fprintf(stderr,
                    "\n"
                    "WARNING: Did not find a frame with velocities in file %s,\n"
                    "         all velocities will be set to zero!\n\n", slog);
            for (i = 0; i < sys->natoms; i++)
            {
                clear_rvec(state->v[i]);
            }
            close_trj(fp);
            /* Search for a frame without velocities */
            bReadVel = FALSE;
            read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X);
        }
    }

    state->natoms = fr.natoms;

    if (sys->natoms != state->natoms)
    {
        gmx_fatal(FARGS, "Number of atoms in Topology "
                  "is not the same as in Trajectory");
    }
    copy_state(slog, &fr, bReadVel, state, &use_time);

    /* Find the appropriate frame */
    while ((fr_time == -1 || fr.time < fr_time) &&
           read_next_frame(oenv, fp, &fr))
    {
        copy_state(slog, &fr, bReadVel, state, &use_time);
    }

    close_trj(fp);

    /* Set the relative box lengths for preserving the box shape.
     * Note that this call can lead to differences in the last bit
     * with respect to using gmx convert-tpr to create a [REF].tpx[ref] file.
     */
    set_box_rel(ir, state);

    fprintf(stderr, "Using frame at t = %g ps\n", use_time);
    fprintf(stderr, "Starting time for run is %g ps\n", ir->init_t);

    if ((ir->epc != epcNO  || ir->etc == etcNOSEHOOVER) && ener)
    {
        get_enx_state(ener, use_time, &sys->groups, ir, state);
        preserve_box_shape(ir, state->box_rel, state->boxv);
    }
}

static void read_posres(gmx_mtop_t *mtop, t_molinfo *molinfo, gmx_bool bTopB,
                        char *fn,
                        int rc_scaling, int ePBC,
                        rvec com,
                        warninp_t wi)
{
    gmx_bool        bFirst = TRUE, *hadAtom;
    rvec           *x, *v, *xp;
    dvec            sum;
    double          totmass;
    t_atoms         dumat;
    matrix          box, invbox;
    int             natoms, npbcdim = 0;
    char            warn_buf[STRLEN], title[STRLEN];
    int             a, i, ai, j, k, mb, nat_molb;
    gmx_molblock_t *molb;
    t_params       *pr, *prfb;
    t_atom         *atom;

    get_stx_coordnum(fn, &natoms);
    if (natoms != mtop->natoms)
    {
        sprintf(warn_buf, "The number of atoms in %s (%d) does not match the number of atoms in the topology (%d). Will assume that the first %d atoms in the topology and %s match.", fn, natoms, mtop->natoms, min(mtop->natoms, natoms), fn);
        warning(wi, warn_buf);
    }
    snew(x, natoms);
    snew(v, natoms);
    init_t_atoms(&dumat, natoms, FALSE);
    read_stx_conf(fn, title, &dumat, x, v, NULL, box);

    npbcdim = ePBC2npbcdim(ePBC);
    clear_rvec(com);
    if (rc_scaling != erscNO)
    {
        copy_mat(box, invbox);
        for (j = npbcdim; j < DIM; j++)
        {
            clear_rvec(invbox[j]);
            invbox[j][j] = 1;
        }
        m_inv_ur0(invbox, invbox);
    }

    /* Copy the reference coordinates to mtop */
    clear_dvec(sum);
    totmass = 0;
    a       = 0;
    snew(hadAtom, natoms);
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb     = &mtop->molblock[mb];
        nat_molb = molb->nmol*mtop->moltype[molb->type].atoms.nr;
        pr       = &(molinfo[molb->type].plist[F_POSRES]);
        prfb     = &(molinfo[molb->type].plist[F_FBPOSRES]);
        if (pr->nr > 0 || prfb->nr > 0)
        {
            atom = mtop->moltype[molb->type].atoms.atom;
            for (i = 0; (i < pr->nr); i++)
            {
                ai = pr->param[i].AI;
                if (ai >= natoms)
                {
                    gmx_fatal(FARGS, "Position restraint atom index (%d) in moltype '%s' is larger than number of atoms in %s (%d).\n",
                              ai+1, *molinfo[molb->type].name, fn, natoms);
                }
                hadAtom[ai] = TRUE;
                if (rc_scaling == erscCOM)
                {
                    /* Determine the center of mass of the posres reference coordinates */
                    for (j = 0; j < npbcdim; j++)
                    {
                        sum[j] += atom[ai].m*x[a+ai][j];
                    }
                    totmass  += atom[ai].m;
                }
            }
            /* Same for flat-bottomed posres, but do not count an atom twice for COM */
            for (i = 0; (i < prfb->nr); i++)
            {
                ai = prfb->param[i].AI;
                if (ai >= natoms)
                {
                    gmx_fatal(FARGS, "Position restraint atom index (%d) in moltype '%s' is larger than number of atoms in %s (%d).\n",
                              ai+1, *molinfo[molb->type].name, fn, natoms);
                }
                if (rc_scaling == erscCOM && hadAtom[ai] == FALSE)
                {
                    /* Determine the center of mass of the posres reference coordinates */
                    for (j = 0; j < npbcdim; j++)
                    {
                        sum[j] += atom[ai].m*x[a+ai][j];
                    }
                    totmass  += atom[ai].m;
                }
            }
            if (!bTopB)
            {
                molb->nposres_xA = nat_molb;
                snew(molb->posres_xA, molb->nposres_xA);
                for (i = 0; i < nat_molb; i++)
                {
                    copy_rvec(x[a+i], molb->posres_xA[i]);
                }
            }
            else
            {
                molb->nposres_xB = nat_molb;
                snew(molb->posres_xB, molb->nposres_xB);
                for (i = 0; i < nat_molb; i++)
                {
                    copy_rvec(x[a+i], molb->posres_xB[i]);
                }
            }
        }
        a += nat_molb;
    }
    if (rc_scaling == erscCOM)
    {
        if (totmass == 0)
        {
            gmx_fatal(FARGS, "The total mass of the position restraint atoms is 0");
        }
        for (j = 0; j < npbcdim; j++)
        {
            com[j] = sum[j]/totmass;
        }
        fprintf(stderr, "The center of mass of the position restraint coord's is %6.3f %6.3f %6.3f\n", com[XX], com[YY], com[ZZ]);
    }

    if (rc_scaling != erscNO)
    {
        assert(npbcdim <= DIM);

        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            molb     = &mtop->molblock[mb];
            nat_molb = molb->nmol*mtop->moltype[molb->type].atoms.nr;
            if (molb->nposres_xA > 0 || molb->nposres_xB > 0)
            {
                xp = (!bTopB ? molb->posres_xA : molb->posres_xB);
                for (i = 0; i < nat_molb; i++)
                {
                    for (j = 0; j < npbcdim; j++)
                    {
                        if (rc_scaling == erscALL)
                        {
                            /* Convert from Cartesian to crystal coordinates */
                            xp[i][j] *= invbox[j][j];
                            for (k = j+1; k < npbcdim; k++)
                            {
                                xp[i][j] += invbox[k][j]*xp[i][k];
                            }
                        }
                        else if (rc_scaling == erscCOM)
                        {
                            /* Subtract the center of mass */
                            xp[i][j] -= com[j];
                        }
                    }
                }
            }
        }

        if (rc_scaling == erscCOM)
        {
            /* Convert the COM from Cartesian to crystal coordinates */
            for (j = 0; j < npbcdim; j++)
            {
                com[j] *= invbox[j][j];
                for (k = j+1; k < npbcdim; k++)
                {
                    com[j] += invbox[k][j]*com[k];
                }
            }
        }
    }

    free_t_atoms(&dumat, TRUE);
    sfree(x);
    sfree(v);
    sfree(hadAtom);
}

static void gen_posres(gmx_mtop_t *mtop, t_molinfo *mi,
                       char *fnA, char *fnB,
                       int rc_scaling, int ePBC,
                       rvec com, rvec comB,
                       warninp_t wi)
{
    int i, j;

    read_posres  (mtop, mi, FALSE, fnA, rc_scaling, ePBC, com, wi);
    /* It is safer to simply read the b-state posres rather than trying
     * to be smart and copy the positions.
     */
    read_posres(mtop, mi, TRUE, fnB, rc_scaling, ePBC, comB, wi);
}

static void set_wall_atomtype(gpp_atomtype_t at, t_gromppopts *opts,
                              t_inputrec *ir, warninp_t wi)
{
    int  i;
    char warn_buf[STRLEN];

    if (ir->nwall > 0)
    {
        fprintf(stderr, "Searching the wall atom type(s)\n");
    }
    for (i = 0; i < ir->nwall; i++)
    {
        ir->wall_atomtype[i] = get_atomtype_type(opts->wall_atomtype[i], at);
        if (ir->wall_atomtype[i] == NOTSET)
        {
            sprintf(warn_buf, "Specified wall atom type %s is not defined", opts->wall_atomtype[i]);
            warning_error(wi, warn_buf);
        }
    }
}

static int nrdf_internal(t_atoms *atoms)
{
    int i, nmass, nrdf;

    nmass = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        /* Vsite ptype might not be set here yet, so also check the mass */
        if ((atoms->atom[i].ptype == eptAtom ||
             atoms->atom[i].ptype == eptNucleus)
            && atoms->atom[i].m > 0)
        {
            nmass++;
        }
    }
    switch (nmass)
    {
        case 0:  nrdf = 0; break;
        case 1:  nrdf = 0; break;
        case 2:  nrdf = 1; break;
        default: nrdf = nmass*3 - 6; break;
    }

    return nrdf;
}

void
spline1d( double        dx,
          double *      y,
          int           n,
          double *      u,
          double *      y2 )
{
    int    i;
    double p, q;

    y2[0] = 0.0;
    u[0]  = 0.0;

    for (i = 1; i < n-1; i++)
    {
        p     = 0.5*y2[i-1]+2.0;
        y2[i] = -0.5/p;
        q     = (y[i+1]-2.0*y[i]+y[i-1])/dx;
        u[i]  = (3.0*q/dx-0.5*u[i-1])/p;
    }

    y2[n-1] = 0.0;

    for (i = n-2; i >= 0; i--)
    {
        y2[i] = y2[i]*y2[i+1]+u[i];
    }
}


void
interpolate1d( double     xmin,
               double     dx,
               double *   ya,
               double *   y2a,
               double     x,
               double *   y,
               double *   y1)
{
    int    ix;
    double a, b;

    ix = (x-xmin)/dx;

    a = (xmin+(ix+1)*dx-x)/dx;
    b = (x-xmin-ix*dx)/dx;

    *y  = a*ya[ix]+b*ya[ix+1]+((a*a*a-a)*y2a[ix]+(b*b*b-b)*y2a[ix+1])*(dx*dx)/6.0;
    *y1 = (ya[ix+1]-ya[ix])/dx-(3.0*a*a-1.0)/6.0*dx*y2a[ix]+(3.0*b*b-1.0)/6.0*dx*y2a[ix+1];
}


void
setup_cmap (int              grid_spacing,
            int              nc,
            real *           grid,
            gmx_cmap_t *     cmap_grid)
{
    double *tmp_u, *tmp_u2, *tmp_yy, *tmp_y1, *tmp_t2, *tmp_grid;

    int     i, j, k, ii, jj, kk, idx;
    int     offset;
    double  dx, xmin, v, v1, v2, v12;
    double  phi, psi;

    snew(tmp_u, 2*grid_spacing);
    snew(tmp_u2, 2*grid_spacing);
    snew(tmp_yy, 2*grid_spacing);
    snew(tmp_y1, 2*grid_spacing);
    snew(tmp_t2, 2*grid_spacing*2*grid_spacing);
    snew(tmp_grid, 2*grid_spacing*2*grid_spacing);

    dx   = 360.0/grid_spacing;
    xmin = -180.0-dx*grid_spacing/2;

    for (kk = 0; kk < nc; kk++)
    {
        /* Compute an offset depending on which cmap we are using
         * Offset will be the map number multiplied with the
         * grid_spacing * grid_spacing * 2
         */
        offset = kk * grid_spacing * grid_spacing * 2;

        for (i = 0; i < 2*grid_spacing; i++)
        {
            ii = (i+grid_spacing-grid_spacing/2)%grid_spacing;

            for (j = 0; j < 2*grid_spacing; j++)
            {
                jj = (j+grid_spacing-grid_spacing/2)%grid_spacing;
                tmp_grid[i*grid_spacing*2+j] = grid[offset+ii*grid_spacing+jj];
            }
        }

        for (i = 0; i < 2*grid_spacing; i++)
        {
            spline1d(dx, &(tmp_grid[2*grid_spacing*i]), 2*grid_spacing, tmp_u, &(tmp_t2[2*grid_spacing*i]));
        }

        for (i = grid_spacing/2; i < grid_spacing+grid_spacing/2; i++)
        {
            ii  = i-grid_spacing/2;
            phi = ii*dx-180.0;

            for (j = grid_spacing/2; j < grid_spacing+grid_spacing/2; j++)
            {
                jj  = j-grid_spacing/2;
                psi = jj*dx-180.0;

                for (k = 0; k < 2*grid_spacing; k++)
                {
                    interpolate1d(xmin, dx, &(tmp_grid[2*grid_spacing*k]),
                                  &(tmp_t2[2*grid_spacing*k]), psi, &tmp_yy[k], &tmp_y1[k]);
                }

                spline1d(dx, tmp_yy, 2*grid_spacing, tmp_u, tmp_u2);
                interpolate1d(xmin, dx, tmp_yy, tmp_u2, phi, &v, &v1);
                spline1d(dx, tmp_y1, 2*grid_spacing, tmp_u, tmp_u2);
                interpolate1d(xmin, dx, tmp_y1, tmp_u2, phi, &v2, &v12);

                idx = ii*grid_spacing+jj;
                cmap_grid->cmapdata[kk].cmap[idx*4]   = grid[offset+ii*grid_spacing+jj];
                cmap_grid->cmapdata[kk].cmap[idx*4+1] = v1;
                cmap_grid->cmapdata[kk].cmap[idx*4+2] = v2;
                cmap_grid->cmapdata[kk].cmap[idx*4+3] = v12;
            }
        }
    }
}

void init_cmap_grid(gmx_cmap_t *cmap_grid, int ngrid, int grid_spacing)
{
    int i, k, nelem;

    cmap_grid->ngrid        = ngrid;
    cmap_grid->grid_spacing = grid_spacing;
    nelem                   = cmap_grid->grid_spacing*cmap_grid->grid_spacing;

    snew(cmap_grid->cmapdata, ngrid);

    for (i = 0; i < cmap_grid->ngrid; i++)
    {
        snew(cmap_grid->cmapdata[i].cmap, 4*nelem);
    }
}


static int count_constraints(gmx_mtop_t *mtop, t_molinfo *mi, warninp_t wi)
{
    int             count, count_mol, i, mb;
    gmx_molblock_t *molb;
    t_params       *plist;
    char            buf[STRLEN];

    count = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        count_mol = 0;
        molb      = &mtop->molblock[mb];
        plist     = mi[molb->type].plist;

        for (i = 0; i < F_NRE; i++)
        {
            if (i == F_SETTLE)
            {
                count_mol += 3*plist[i].nr;
            }
            else if (interaction_function[i].flags & IF_CONSTRAINT)
            {
                count_mol += plist[i].nr;
            }
        }

        if (count_mol > nrdf_internal(&mi[molb->type].atoms))
        {
            sprintf(buf,
                    "Molecule type '%s' has %d constraints.\n"
                    "For stability and efficiency there should not be more constraints than internal number of degrees of freedom: %d.\n",
                    *mi[molb->type].name, count_mol,
                    nrdf_internal(&mi[molb->type].atoms));
            warning(wi, buf);
        }
        count += molb->nmol*count_mol;
    }

    return count;
}

static void check_gbsa_params_charged(gmx_mtop_t *sys, gpp_atomtype_t atype)
{
    int            i, nmiss, natoms, mt;
    real           q;
    const t_atoms *atoms;

    nmiss = 0;
    for (mt = 0; mt < sys->nmoltype; mt++)
    {
        atoms  = &sys->moltype[mt].atoms;
        natoms = atoms->nr;

        for (i = 0; i < natoms; i++)
        {
            q = atoms->atom[i].q;
            if ((get_atomtype_radius(atoms->atom[i].type, atype)    == 0  ||
                 get_atomtype_vol(atoms->atom[i].type, atype)       == 0  ||
                 get_atomtype_surftens(atoms->atom[i].type, atype)  == 0  ||
                 get_atomtype_gb_radius(atoms->atom[i].type, atype) == 0  ||
                 get_atomtype_S_hct(atoms->atom[i].type, atype)     == 0) &&
                q != 0)
            {
                fprintf(stderr, "\nGB parameter(s) zero for atom type '%s' while charge is %g\n",
                        get_atomtype_name(atoms->atom[i].type, atype), q);
                nmiss++;
            }
        }
    }

    if (nmiss > 0)
    {
        gmx_fatal(FARGS, "Can't do GB electrostatics; the implicit_genborn_params section of the forcefield has parameters with value zero for %d atomtypes that occur as charged atoms.", nmiss);
    }
}


static void check_gbsa_params(gpp_atomtype_t atype)
{
    int  nmiss, i;

    /* If we are doing GBSA, check that we got the parameters we need
     * This checking is to see if there are GBSA paratmeters for all
     * atoms in the force field. To go around this for testing purposes
     * comment out the nerror++ counter temporarily
     */
    nmiss = 0;
    for (i = 0; i < get_atomtype_ntypes(atype); i++)
    {
        if (get_atomtype_radius(i, atype)    < 0 ||
            get_atomtype_vol(i, atype)       < 0 ||
            get_atomtype_surftens(i, atype)  < 0 ||
            get_atomtype_gb_radius(i, atype) < 0 ||
            get_atomtype_S_hct(i, atype)     < 0)
        {
            fprintf(stderr, "\nGB parameter(s) missing or negative for atom type '%s'\n",
                    get_atomtype_name(i, atype));
            nmiss++;
        }
    }

    if (nmiss > 0)
    {
        gmx_fatal(FARGS, "Can't do GB electrostatics; the implicit_genborn_params section of the forcefield is missing parameters for %d atomtypes or they might be negative.", nmiss);
    }

}

static real calc_temp(const gmx_mtop_t *mtop,
                      const t_inputrec *ir,
                      rvec             *v)
{
    double                  sum_mv2;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     a;
    int                     nrdf, g;

    sum_mv2 = 0;

    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &a, &atom))
    {
        sum_mv2 += atom->m*norm2(v[a]);
    }

    nrdf = 0;
    for (g = 0; g < ir->opts.ngtc; g++)
    {
        nrdf += ir->opts.nrdf[g];
    }

    return sum_mv2/(nrdf*BOLTZ);
}

static real get_max_reference_temp(const t_inputrec *ir,
                                   warninp_t         wi)
{
    real     ref_t;
    int      i;
    gmx_bool bNoCoupl;

    ref_t    = 0;
    bNoCoupl = FALSE;
    for (i = 0; i < ir->opts.ngtc; i++)
    {
        if (ir->opts.tau_t[i] < 0)
        {
            bNoCoupl = TRUE;
        }
        else
        {
            ref_t = max(ref_t, ir->opts.ref_t[i]);
        }
    }

    if (bNoCoupl)
    {
        char buf[STRLEN];

        sprintf(buf, "Some temperature coupling groups do not use temperature coupling. We will assume their temperature is not more than %.3f K. If their temperature is higher, the energy error and the Verlet buffer might be underestimated.",
                ref_t);
        warning(wi, buf);
    }

    return ref_t;
}

static void set_verlet_buffer(const gmx_mtop_t *mtop,
                              t_inputrec       *ir,
                              real              buffer_temp,
                              matrix            box,
                              warninp_t         wi)
{
    int                    i;
    verletbuf_list_setup_t ls;
    real                   rlist_1x1;
    int                    n_nonlin_vsite;
    char                   warn_buf[STRLEN];

    printf("Determining Verlet buffer for a tolerance of %g kJ/mol/ps at %g K\n", ir->verletbuf_tol, buffer_temp);

    /* Calculate the buffer size for simple atom vs atoms list */
    ls.cluster_size_i = 1;
    ls.cluster_size_j = 1;
    calc_verlet_buffer_size(mtop, det(box), ir, buffer_temp,
                            &ls, &n_nonlin_vsite, &rlist_1x1);

    /* Set the pair-list buffer size in ir */
    verletbuf_get_list_setup(FALSE, FALSE, &ls);
    calc_verlet_buffer_size(mtop, det(box), ir, buffer_temp,
                            &ls, &n_nonlin_vsite, &ir->rlist);

    if (n_nonlin_vsite > 0)
    {
        sprintf(warn_buf, "There are %d non-linear virtual site constructions. Their contribution to the energy error is approximated. In most cases this does not affect the error significantly.", n_nonlin_vsite);
        warning_note(wi, warn_buf);
    }

    printf("Calculated rlist for %dx%d atom pair-list as %.3f nm, buffer size %.3f nm\n",
           1, 1, rlist_1x1, rlist_1x1-max(ir->rvdw, ir->rcoulomb));

    ir->rlistlong = ir->rlist;
    printf("Set rlist, assuming %dx%d atom pair-list, to %.3f nm, buffer size %.3f nm\n",
           ls.cluster_size_i, ls.cluster_size_j,
           ir->rlist, ir->rlist-max(ir->rvdw, ir->rcoulomb));

    printf("Note that mdrun will redetermine rlist based on the actual pair-list setup\n");

    if (sqr(ir->rlistlong) >= max_cutoff2(ir->ePBC, box))
    {
        gmx_fatal(FARGS, "The pair-list cut-off (%g nm) is longer than half the shortest box vector or longer than the smallest box diagonal element (%g nm). Increase the box size or decrease nstlist or increase verlet-buffer-tolerance.", ir->rlistlong, sqrt(max_cutoff2(ir->ePBC, box)));
    }
}

int gmx_grompp(int argc, char *argv[])
{
    static const char *desc[] = {
        "[THISMODULE] (the gromacs preprocessor)",
        "reads a molecular topology file, checks the validity of the",
        "file, expands the topology from a molecular description to an atomic",
        "description. The topology file contains information about",
        "molecule types and the number of molecules, the preprocessor",
        "copies each molecule as needed. ",
        "There is no limitation on the number of molecule types. ",
        "Bonds and bond-angles can be converted into constraints, separately",
        "for hydrogens and heavy atoms.",
        "Then a coordinate file is read and velocities can be generated",
        "from a Maxwellian distribution if requested.",
        "[THISMODULE] also reads parameters for [gmx-mdrun] ",
        "(eg. number of MD steps, time step, cut-off), and others such as",
        "NEMD parameters, which are corrected so that the net acceleration",
        "is zero.",
        "Eventually a binary file is produced that can serve as the sole input",
        "file for the MD program.[PAR]",

        "[THISMODULE] uses the atom names from the topology file. The atom names",
        "in the coordinate file (option [TT]-c[tt]) are only read to generate",
        "warnings when they do not match the atom names in the topology.",
        "Note that the atom names are irrelevant for the simulation as",
        "only the atom types are used for generating interaction parameters.[PAR]",

        "[THISMODULE] uses a built-in preprocessor to resolve includes, macros, ",
        "etc. The preprocessor supports the following keywords::",
        "",
        "    #ifdef VARIABLE",
        "    #ifndef VARIABLE",
        "    #else",
        "    #endif",
        "    #define VARIABLE",
        "    #undef VARIABLE",
        "    #include \"filename\"",
        "    #include <filename>",
        "",
        "The functioning of these statements in your topology may be modulated by",
        "using the following two flags in your [REF].mdp[ref] file::",
        "",
        "    define = -DVARIABLE1 -DVARIABLE2",
        "    include = -I/home/john/doe",
        "",
        "For further information a C-programming textbook may help you out.",
        "Specifying the [TT]-pp[tt] flag will get the pre-processed",
        "topology file written out so that you can verify its contents.[PAR]",

        "When using position restraints a file with restraint coordinates",
        "can be supplied with [TT]-r[tt], otherwise restraining will be done",
        "with respect to the conformation from the [TT]-c[tt] option.",
        "For free energy calculation the the coordinates for the B topology",
        "can be supplied with [TT]-rb[tt], otherwise they will be equal to",
        "those of the A topology.[PAR]",

        "Starting coordinates can be read from trajectory with [TT]-t[tt].",
        "The last frame with coordinates and velocities will be read,",
        "unless the [TT]-time[tt] option is used. Only if this information",
        "is absent will the coordinates in the [TT]-c[tt] file be used.",
        "Note that these velocities will not be used when [TT]gen_vel = yes[tt]",
        "in your [REF].mdp[ref] file. An energy file can be supplied with",
        "[TT]-e[tt] to read Nose-Hoover and/or Parrinello-Rahman coupling",
        "variables.[PAR]",

        "[THISMODULE] can be used to restart simulations (preserving",
        "continuity) by supplying just a checkpoint file with [TT]-t[tt].",
        "However, for simply changing the number of run steps to extend",
        "a run, using [gmx-convert-tpr] is more convenient than [THISMODULE].",
        "You then supply the old checkpoint file directly to [gmx-mdrun]",
        "with [TT]-cpi[tt]. If you wish to change the ensemble or things",
        "like output frequency, then supplying the checkpoint file to",
        "[THISMODULE] with [TT]-t[tt] along with a new [REF].mdp[ref] file",
        "with [TT]-f[tt] is the recommended procedure. Actually preserving",
        "the ensemble (if possible) still requires passing the checkpoint",
        "file to [gmx-mdrun] [TT]-cpi[tt].[PAR]",

        "By default, all bonded interactions which have constant energy due to",
        "virtual site constructions will be removed. If this constant energy is",
        "not zero, this will result in a shift in the total energy. All bonded",
        "interactions can be kept by turning off [TT]-rmvsbds[tt]. Additionally,",
        "all constraints for distances which will be constant anyway because",
        "of virtual site constructions will be removed. If any constraints remain",
        "which involve virtual sites, a fatal error will result.[PAR]"

        "To verify your run input file, please take note of all warnings",
        "on the screen, and correct where necessary. Do also look at the contents",
        "of the [TT]mdout.mdp[tt] file; this contains comment lines, as well as",
        "the input that [THISMODULE] has read. If in doubt, you can start [THISMODULE]",
        "with the [TT]-debug[tt] option which will give you more information",
        "in a file called [TT]grompp.log[tt] (along with real debug info). You",
        "can see the contents of the run input file with the [gmx-dump]",
        "program. [gmx-check] can be used to compare the contents of two",
        "run input files.[PAR]"

        "The [TT]-maxwarn[tt] option can be used to override warnings printed",
        "by [THISMODULE] that otherwise halt output. In some cases, warnings are",
        "harmless, but usually they are not. The user is advised to carefully",
        "interpret the output messages before attempting to bypass them with",
        "this option."
    };
    t_gromppopts      *opts;
    gmx_mtop_t        *sys;
    int                nmi;
    t_molinfo         *mi, *intermolecular_interactions;
    gpp_atomtype_t     atype;
    t_inputrec        *ir;
    int                natoms, nvsite, comb, mt;
    t_params          *plist;
    t_state           *state;
    matrix             box;
    real               max_spacing, fudgeQQ;
    double             reppow;
    char               fn[STRLEN], fnB[STRLEN];
    const char        *mdparin;
    int                ntype;
    gmx_bool           bNeedVel, bGenVel;
    gmx_bool           have_atomnumber;
    int                n12, n13, n14;
    t_params          *gb_plist = NULL;
    gmx_genborn_t     *born     = NULL;
    output_env_t       oenv;
    gmx_bool           bVerbose = FALSE;
    warninp_t          wi;
    char               warn_buf[STRLEN];
    unsigned int       useed;
    t_atoms            IMDatoms;   /* Atoms to be operated on interactively (IMD) */

    t_filenm           fnm[] = {
        { efMDP, NULL,  NULL,        ffREAD  },
        { efMDP, "-po", "mdout",     ffWRITE },
        { efSTX, "-c",  NULL,        ffREAD  },
        { efSTX, "-r",  NULL,        ffOPTRD },
        { efSTX, "-rb", NULL,        ffOPTRD },
        { efNDX, NULL,  NULL,        ffOPTRD },
        { efTOP, NULL,  NULL,        ffREAD  },
        { efTOP, "-pp", "processed", ffOPTWR },
        { efTPR, "-o",  NULL,        ffWRITE },
        { efTRN, "-t",  NULL,        ffOPTRD },
        { efEDR, "-e",  NULL,        ffOPTRD },
        /* This group is needed by the VMD viewer as the start configuration for IMD sessions: */
        { efGRO, "-imd", "imdgroup", ffOPTWR },
        { efTRN, "-ref", "rotref",   ffOPTRW }
    };
#define NFILE asize(fnm)

    /* Command line options */
    static gmx_bool bRenum   = TRUE;
    static gmx_bool bRmVSBds = TRUE, bZero = FALSE;
    static int      i, maxwarn = 0;
    static real     fr_time = -1;
    t_pargs         pa[]    = {
        { "-v",       FALSE, etBOOL, {&bVerbose},
          "Be loud and noisy" },
        { "-time",    FALSE, etREAL, {&fr_time},
          "Take frame at or first after this time." },
        { "-rmvsbds", FALSE, etBOOL, {&bRmVSBds},
          "Remove constant bonded interactions with virtual sites" },
        { "-maxwarn", FALSE, etINT,  {&maxwarn},
          "Number of allowed warnings during input processing. Not for normal use and may generate unstable systems" },
        { "-zero",    FALSE, etBOOL, {&bZero},
          "Set parameters for bonded interactions without defaults to zero instead of generating an error" },
        { "-renum",   FALSE, etBOOL, {&bRenum},
          "Renumber atomtypes and minimize number of atomtypes" }
    };

    /* Initiate some variables */
    snew(ir, 1);
    snew(opts, 1);
    init_ir(ir, opts);

    /* Parse the command line */
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    wi = init_warning(TRUE, maxwarn);

    /* PARAMETER file processing */
    mdparin = opt2fn("-f", NFILE, fnm);
    set_warning_line(wi, mdparin, -1);
    get_ir(mdparin, opt2fn("-po", NFILE, fnm), ir, opts, wi);

    if (bVerbose)
    {
        fprintf(stderr, "checking input for internal consistency...\n");
    }
    check_ir(mdparin, ir, opts, wi);

    if (ir->ld_seed == -1)
    {
        ir->ld_seed = (gmx_int64_t)gmx_rng_make_seed();
        fprintf(stderr, "Setting the LD random seed to %"GMX_PRId64 "\n", ir->ld_seed);
    }

    if (ir->expandedvals->lmc_seed == -1)
    {
        ir->expandedvals->lmc_seed = (int)gmx_rng_make_seed();
        fprintf(stderr, "Setting the lambda MC random seed to %d\n", ir->expandedvals->lmc_seed);
    }

    bNeedVel = EI_STATE_VELOCITY(ir->eI);
    bGenVel  = (bNeedVel && opts->bGenVel);
    if (bGenVel && ir->bContinuation)
    {
        sprintf(warn_buf,
                "Generating velocities is inconsistent with attempting "
                "to continue a previous run. Choose only one of "
                "gen-vel = yes and continuation = yes.");
        warning_error(wi, warn_buf);
    }

    snew(plist, F_NRE);
    init_plist(plist);
    snew(sys, 1);
    atype = init_atomtype();
    if (debug)
    {
        pr_symtab(debug, 0, "Just opened", &sys->symtab);
    }

    strcpy(fn, ftp2fn(efTOP, NFILE, fnm));
    if (!gmx_fexist(fn))
    {
        gmx_fatal(FARGS, "%s does not exist", fn);
    }
    snew(state, 1);
    new_status(fn, opt2fn_null("-pp", NFILE, fnm), opt2fn("-c", NFILE, fnm),
               opts, ir, bZero, bGenVel, bVerbose, state,
               atype, sys, &nmi, &mi, &intermolecular_interactions,
               plist, &comb, &reppow, &fudgeQQ,
               opts->bMorse,
               wi);

    if (debug)
    {
        pr_symtab(debug, 0, "After new_status", &sys->symtab);
    }

    nvsite = 0;
    /* set parameters for virtual site construction (not for vsiten) */
    for (mt = 0; mt < sys->nmoltype; mt++)
    {
        nvsite +=
            set_vsites(bVerbose, &sys->moltype[mt].atoms, atype, mi[mt].plist);
    }
    /* now throw away all obsolete bonds, angles and dihedrals: */
    /* note: constraints are ALWAYS removed */
    if (nvsite)
    {
        for (mt = 0; mt < sys->nmoltype; mt++)
        {
            clean_vsite_bondeds(mi[mt].plist, sys->moltype[mt].atoms.nr, bRmVSBds);
        }
    }

    if (nvsite && ir->eI == eiNM)
    {
        gmx_fatal(FARGS, "Normal Mode analysis is not supported with virtual sites.\nIf you'd like to help with adding support, we have an open discussion at http://redmine.gromacs.org/issues/879\n");
    }

    if (ir->cutoff_scheme == ecutsVERLET)
    {
        fprintf(stderr, "Removing all charge groups because cutoff-scheme=%s\n",
                ecutscheme_names[ir->cutoff_scheme]);

        /* Remove all charge groups */
        gmx_mtop_remove_chargegroups(sys);
    }

    if (count_constraints(sys, mi, wi) && (ir->eConstrAlg == econtSHAKE))
    {
        if (ir->eI == eiCG || ir->eI == eiLBFGS)
        {
            sprintf(warn_buf, "Can not do %s with %s, use %s",
                    EI(ir->eI), econstr_names[econtSHAKE], econstr_names[econtLINCS]);
            warning_error(wi, warn_buf);
        }
        if (ir->bPeriodicMols)
        {
            sprintf(warn_buf, "Can not do periodic molecules with %s, use %s",
                    econstr_names[econtSHAKE], econstr_names[econtLINCS]);
            warning_error(wi, warn_buf);
        }
    }

    if (EI_SD (ir->eI) &&  ir->etc != etcNO)
    {
        warning_note(wi, "Temperature coupling is ignored with SD integrators.");
    }

    /* If we are doing QM/MM, check that we got the atom numbers */
    have_atomnumber = TRUE;
    for (i = 0; i < get_atomtype_ntypes(atype); i++)
    {
        have_atomnumber = have_atomnumber && (get_atomtype_atomnumber(i, atype) >= 0);
    }
    if (!have_atomnumber && ir->bQMMM)
    {
        warning_error(wi,
                      "\n"
                      "It appears as if you are trying to run a QM/MM calculation, but the force\n"
                      "field you are using does not contain atom numbers fields. This is an\n"
                      "optional field (introduced in GROMACS 3.3) for general runs, but mandatory\n"
                      "for QM/MM. The good news is that it is easy to add - put the atom number as\n"
                      "an integer just before the mass column in ffXXXnb.itp.\n"
                      "NB: United atoms have the same atom numbers as normal ones.\n\n");
    }

    if (ir->bAdress)
    {
        if ((ir->adress->const_wf > 1) || (ir->adress->const_wf < 0))
        {
            warning_error(wi, "AdResS contant weighting function should be between 0 and 1\n\n");
        }
        /** TODO check size of ex+hy width against box size */
    }

    /* Check for errors in the input now, since they might cause problems
     * during processing further down.
     */
    check_warning_error(wi, FARGS);

    if (opt2bSet("-r", NFILE, fnm))
    {
        sprintf(fn, "%s", opt2fn("-r", NFILE, fnm));
    }
    else
    {
        sprintf(fn, "%s", opt2fn("-c", NFILE, fnm));
    }
    if (opt2bSet("-rb", NFILE, fnm))
    {
        sprintf(fnB, "%s", opt2fn("-rb", NFILE, fnm));
    }
    else
    {
        strcpy(fnB, fn);
    }

    if (nint_ftype(sys, mi, F_POSRES) > 0 || nint_ftype(sys, mi, F_FBPOSRES) > 0)
    {
        if (bVerbose)
        {
            fprintf(stderr, "Reading position restraint coords from %s", fn);
            if (strcmp(fn, fnB) == 0)
            {
                fprintf(stderr, "\n");
            }
            else
            {
                fprintf(stderr, " and %s\n", fnB);
            }
        }
        gen_posres(sys, mi, fn, fnB,
                   ir->refcoord_scaling, ir->ePBC,
                   ir->posres_com, ir->posres_comB,
                   wi);
    }

    /* If we are using CMAP, setup the pre-interpolation grid */
    if (plist[F_CMAP].ncmap > 0)
    {
        init_cmap_grid(&sys->ffparams.cmap_grid, plist[F_CMAP].nc, plist[F_CMAP].grid_spacing);
        setup_cmap(plist[F_CMAP].grid_spacing, plist[F_CMAP].nc, plist[F_CMAP].cmap, &sys->ffparams.cmap_grid);
    }

    set_wall_atomtype(atype, opts, ir, wi);
    if (bRenum)
    {
        renum_atype(plist, sys, ir->wall_atomtype, atype, bVerbose);
        ntype = get_atomtype_ntypes(atype);
    }

    if (ir->implicit_solvent != eisNO)
    {
        /* Now we have renumbered the atom types, we can check the GBSA params */
        check_gbsa_params(atype);

        /* Check that all atoms that have charge and/or LJ-parameters also have
         * sensible GB-parameters
         */
        check_gbsa_params_charged(sys, atype);
    }

    /* PELA: Copy the atomtype data to the topology atomtype list */
    copy_atomtype_atomtypes(atype, &(sys->atomtypes));

    if (debug)
    {
        pr_symtab(debug, 0, "After renum_atype", &sys->symtab);
    }

    if (bVerbose)
    {
        fprintf(stderr, "converting bonded parameters...\n");
    }

    ntype = get_atomtype_ntypes(atype);
    convert_params(ntype, plist, mi, intermolecular_interactions,
                   comb, reppow, fudgeQQ, sys);

    if (debug)
    {
        pr_symtab(debug, 0, "After convert_params", &sys->symtab);
    }

    /* set ptype to VSite for virtual sites */
    for (mt = 0; mt < sys->nmoltype; mt++)
    {
        set_vsites_ptype(FALSE, &sys->moltype[mt]);
    }
    if (debug)
    {
        pr_symtab(debug, 0, "After virtual sites", &sys->symtab);
    }
    /* Check velocity for virtual sites and shells */
    if (bGenVel)
    {
        check_vel(sys, state->v);
    }

    /* check for shells and inpurecs */
    check_shells_inputrec(sys, ir, wi);

    /* check masses */
    check_mol(sys, wi);

    for (i = 0; i < sys->nmoltype; i++)
    {
        check_cg_sizes(ftp2fn(efTOP, NFILE, fnm), &sys->moltype[i].cgs, wi);
    }

    if (EI_DYNAMICS(ir->eI) && ir->eI != eiBD)
    {
        check_bonds_timestep(sys, ir->delta_t, wi);
    }

    if (EI_ENERGY_MINIMIZATION(ir->eI) && 0 == ir->nsteps)
    {
        warning_note(wi, "Zero-step energy minimization will alter the coordinates before calculating the energy. If you just want the energy of a single point, try zero-step MD (with unconstrained_start = yes). To do multiple single-point energy evaluations of different configurations of the same topology, use mdrun -rerun.");
    }

    check_warning_error(wi, FARGS);

    if (bVerbose)
    {
        fprintf(stderr, "initialising group options...\n");
    }
    do_index(mdparin, ftp2fn_null(efNDX, NFILE, fnm),
             sys, bVerbose, ir,
             bGenVel ? state->v : NULL,
             wi);

    if (ir->cutoff_scheme == ecutsVERLET && ir->verletbuf_tol > 0 &&
        ir->nstlist > 1)
    {
        if (EI_DYNAMICS(ir->eI) && inputrec2nboundeddim(ir) == 3)
        {
            real buffer_temp;

            if (EI_MD(ir->eI) && ir->etc == etcNO)
            {
                if (bGenVel)
                {
                    buffer_temp = opts->tempi;
                }
                else
                {
                    buffer_temp = calc_temp(sys, ir, state->v);
                }
                if (buffer_temp > 0)
                {
                    sprintf(warn_buf, "NVE simulation: will use the initial temperature of %.3f K for determining the Verlet buffer size", buffer_temp);
                    warning_note(wi, warn_buf);
                }
                else
                {
                    sprintf(warn_buf, "NVE simulation with an initial temperature of zero: will use a Verlet buffer of %d%%. Check your energy drift!",
                            (int)(verlet_buffer_ratio_NVE_T0*100 + 0.5));
                    warning_note(wi, warn_buf);
                }
            }
            else
            {
                buffer_temp = get_max_reference_temp(ir, wi);
            }

            if (EI_MD(ir->eI) && ir->etc == etcNO && buffer_temp == 0)
            {
                /* NVE with initial T=0: we add a fixed ratio to rlist.
                 * Since we don't actually use verletbuf_tol, we set it to -1
                 * so it can't be misused later.
                 */
                ir->rlist         *= 1.0 + verlet_buffer_ratio_NVE_T0;
                ir->verletbuf_tol  = -1;
            }
            else
            {
                /* We warn for NVE simulations with a drift tolerance that
                 * might result in a 1(.1)% drift over the total run-time.
                 * Note that we can't warn when nsteps=0, since we don't
                 * know how many steps the user intends to run.
                 */
                if (EI_MD(ir->eI) && ir->etc == etcNO && ir->nstlist > 1 &&
                    ir->nsteps > 0)
                {
                    const real driftTolerance = 0.01;
                    /* We use 2 DOF per atom = 2kT pot+kin energy,
                     * to be on the safe side with constraints.
                     */
                    const real totalEnergyDriftPerAtomPerPicosecond = 2*BOLTZ*buffer_temp/(ir->nsteps*ir->delta_t);

                    if (ir->verletbuf_tol > 1.1*driftTolerance*totalEnergyDriftPerAtomPerPicosecond)
                    {
                        sprintf(warn_buf, "You are using a Verlet buffer tolerance of %g kJ/mol/ps for an NVE simulation of length %g ps, which can give a final drift of %d%%. For conserving energy to %d%% when using constraints, you might need to set verlet-buffer-tolerance to %.1e.",
                                ir->verletbuf_tol, ir->nsteps*ir->delta_t,
                                (int)(ir->verletbuf_tol/totalEnergyDriftPerAtomPerPicosecond*100 + 0.5),
                                (int)(100*driftTolerance + 0.5),
                                driftTolerance*totalEnergyDriftPerAtomPerPicosecond);
                        warning_note(wi, warn_buf);
                    }
                }

                set_verlet_buffer(sys, ir, buffer_temp, state->box, wi);
            }
        }
    }

    /* Init the temperature coupling state */
    init_gtc_state(state, ir->opts.ngtc, 0, ir->opts.nhchainlength); /* need to add nnhpres here? */

    if (bVerbose)
    {
        fprintf(stderr, "Checking consistency between energy and charge groups...\n");
    }
    check_eg_vs_cg(sys);

    if (debug)
    {
        pr_symtab(debug, 0, "After index", &sys->symtab);
    }

    triple_check(mdparin, ir, sys, wi);
    close_symtab(&sys->symtab);
    if (debug)
    {
        pr_symtab(debug, 0, "After close", &sys->symtab);
    }

    /* make exclusions between QM atoms */
    if (ir->bQMMM)
    {
        if (ir->QMMMscheme == eQMMMschemenormal && ir->ns_type == ensSIMPLE)
        {
            gmx_fatal(FARGS, "electrostatic embedding only works with grid neighboursearching, use ns-type=grid instead\n");
        }
        else
        {
            generate_qmexcl(sys, ir, wi);
        }
    }

    if (ftp2bSet(efTRN, NFILE, fnm))
    {
        if (bVerbose)
        {
            fprintf(stderr, "getting data from old trajectory ...\n");
        }
        cont_status(ftp2fn(efTRN, NFILE, fnm), ftp2fn_null(efEDR, NFILE, fnm),
                    bNeedVel, bGenVel, fr_time, ir, state, sys, oenv);
    }

    if (ir->ePBC == epbcXY && ir->nwall != 2)
    {
        clear_rvec(state->box[ZZ]);
    }

    if (ir->cutoff_scheme != ecutsVERLET && ir->rlist > 0)
    {
        set_warning_line(wi, mdparin, -1);
        check_chargegroup_radii(sys, ir, state->x, wi);
    }

    if (EEL_FULL(ir->coulombtype) || EVDW_PME(ir->vdwtype))
    {
        /* Calculate the optimal grid dimensions */
        copy_mat(state->box, box);
        if (ir->ePBC == epbcXY && ir->nwall == 2)
        {
            svmul(ir->wall_ewald_zfac, box[ZZ], box[ZZ]);
        }
        if (ir->nkx > 0 && ir->nky > 0 && ir->nkz > 0)
        {
            /* Mark fourier_spacing as not used */
            ir->fourier_spacing = 0;
        }
        else if (ir->nkx != 0 && ir->nky != 0 && ir->nkz != 0)
        {
            set_warning_line(wi, mdparin, -1);
            warning_error(wi, "Some of the Fourier grid sizes are set, but all of them need to be set.");
        }
        max_spacing = calc_grid(stdout, box, ir->fourier_spacing,
                                &(ir->nkx), &(ir->nky), &(ir->nkz));
    }

    /* MRS: eventually figure out better logic for initializing the fep
       values that makes declaring the lambda and declaring the state not
       potentially conflict if not handled correctly. */
    if (ir->efep != efepNO)
    {
        state->fep_state = ir->fepvals->init_fep_state;
        for (i = 0; i < efptNR; i++)
        {
            /* init_lambda trumps state definitions*/
            if (ir->fepvals->init_lambda >= 0)
            {
                state->lambda[i] = ir->fepvals->init_lambda;
            }
            else
            {
                if (ir->fepvals->all_lambda[i] == NULL)
                {
                    gmx_fatal(FARGS, "Values of lambda not set for a free energy calculation!");
                }
                else
                {
                    state->lambda[i] = ir->fepvals->all_lambda[i][state->fep_state];
                }
            }
        }
    }

    if (ir->bPull)
    {
        set_pull_init(ir, sys, state->x, state->box, state->lambda[efptMASS], oenv);
    }

    if (ir->bRot)
    {
        set_reference_positions(ir->rot, state->x, state->box,
                                opt2fn("-ref", NFILE, fnm), opt2bSet("-ref", NFILE, fnm),
                                wi);
    }

    /*  reset_multinr(sys); */

    if (EEL_PME(ir->coulombtype))
    {
        float ratio = pme_load_estimate(sys, ir, state->box);
        fprintf(stderr, "Estimate for the relative computational load of the PME mesh part: %.2f\n", ratio);
        /* With free energy we might need to do PME both for the A and B state
         * charges. This will double the cost, but the optimal performance will
         * then probably be at a slightly larger cut-off and grid spacing.
         */
        if ((ir->efep == efepNO && ratio > 1.0/2.0) ||
            (ir->efep != efepNO && ratio > 2.0/3.0))
        {
            warning_note(wi,
                         "The optimal PME mesh load for parallel simulations is below 0.5\n"
                         "and for highly parallel simulations between 0.25 and 0.33,\n"
                         "for higher performance, increase the cut-off and the PME grid spacing.\n");
            if (ir->efep != efepNO)
            {
                warning_note(wi,
                             "For free energy simulations, the optimal load limit increases from 0.5 to 0.667\n");
            }
        }
    }

    {
        char   warn_buf[STRLEN];
        double cio = compute_io(ir, sys->natoms, &sys->groups, F_NRE, 1);
        sprintf(warn_buf, "This run will generate roughly %.0f Mb of data", cio);
        if (cio > 2000)
        {
            set_warning_line(wi, mdparin, -1);
            warning_note(wi, warn_buf);
        }
        else
        {
            printf("%s\n", warn_buf);
        }
    }

    if (bVerbose)
    {
        fprintf(stderr, "writing run input file...\n");
    }

    done_warning(wi, FARGS);
    write_tpx_state(ftp2fn(efTPR, NFILE, fnm), ir, state, sys);

    /* Output IMD group, if bIMD is TRUE */
    write_IMDgroup_to_file(ir->bIMD, ir, state, sys, NFILE, fnm);

    done_state(state);
    sfree(state);
    done_atomtype(atype);
    done_mtop(sys, TRUE);
    done_inputrec_strings();

    return 0;
}
