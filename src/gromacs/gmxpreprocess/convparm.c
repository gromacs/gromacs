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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "convparm.h"

#include <math.h>
#include <string.h>

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static int round_check(real r, int limit, int ftype, const char *name)
{
    int i;

    if (r >= 0)
    {
        i = (int)(r + 0.5);
    }
    else
    {
        i = (int)(r - 0.5);
    }

    if (r-i > 0.01 || r-i < -0.01)
    {
        gmx_fatal(FARGS, "A non-integer value (%f) was supplied for '%s' in %s",
                  r, name, interaction_function[ftype].longname);
    }

    if (i < limit)
    {
        gmx_fatal(FARGS, "Value of '%s' in %s is %d, which is smaller than the minimum of %d",
                  name, interaction_function[ftype].longname, i, limit);
    }

    return i;
}

static void set_ljparams(int comb, double reppow, real v, real w,
                         real *c6, real *c12)
{
    if (comb == eCOMB_ARITHMETIC || comb == eCOMB_GEOM_SIG_EPS)
    {
        if (v >= 0)
        {
            *c6  = 4*w*pow(v, 6.0);
            *c12 = 4*w*pow(v, reppow);
        }
        else
        {
            /* Interpret negative sigma as c6=0 and c12 with -sigma */
            *c6  = 0;
            *c12 = 4*w*pow(-v, reppow);
        }
    }
    else
    {
        *c6  = v;
        *c12 = w;
    }
}

/* A return value of 0 means parameters were assigned successfully,
 * returning -1 means this is an all-zero interaction that should not be added.
 */
static int
assign_param(t_functype ftype, t_iparams *newparam,
             real old[MAXFORCEPARAM], int comb, double reppow)
{
    int      i, j;
    real     tmp;
    gmx_bool all_param_zero = TRUE;

    /* Set to zero */
    for (j = 0; (j < MAXFORCEPARAM); j++)
    {
        newparam->generic.buf[j] = 0.0;
        /* If all parameters are zero we might not add some interaction types (selected below).
         * We cannot apply this to ALL interactions, since many have valid reasons for having
         * zero parameters (e.g. an index to a Cmap interaction, or LJ parameters), but
         * we use it for angles and torsions that are typically generated automatically.
         */
        all_param_zero = (all_param_zero == TRUE) && fabs(old[j]) < GMX_REAL_MIN;
    }

    if (all_param_zero == TRUE)
    {
        if (IS_ANGLE(ftype) || IS_RESTRAINT_TYPE(ftype) || ftype == F_IDIHS ||
            ftype == F_PDIHS || ftype == F_PIDIHS || ftype == F_RBDIHS || ftype == F_FOURDIHS)
        {
            return -1;
        }
    }

    switch (ftype)
    {
        case F_G96ANGLES:
            /* Post processing of input data: store cosine iso angle itself */
            newparam->harmonic.rA  = cos(old[0]*DEG2RAD);
            newparam->harmonic.krA = old[1];
            newparam->harmonic.rB  = cos(old[2]*DEG2RAD);
            newparam->harmonic.krB = old[3];
            break;
        case F_G96BONDS:
            /* Post processing of input data: store square of length itself */
            newparam->harmonic.rA  = sqr(old[0]);
            newparam->harmonic.krA = old[1];
            newparam->harmonic.rB  = sqr(old[2]);
            newparam->harmonic.krB = old[3];
            break;
        case F_FENEBONDS:
            newparam->fene.bm = old[0];
            newparam->fene.kb = old[1];
            break;
        case F_RESTRBONDS:
            newparam->restraint.lowA = old[0];
            newparam->restraint.up1A = old[1];
            newparam->restraint.up2A = old[2];
            newparam->restraint.kA   = old[3];
            newparam->restraint.lowB = old[4];
            newparam->restraint.up1B = old[5];
            newparam->restraint.up2B = old[6];
            newparam->restraint.kB   = old[7];
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            newparam->tab.table = round_check(old[0], 0, ftype, "table index");
            newparam->tab.kA    = old[1];
            newparam->tab.kB    = old[3];
            break;
        case F_CROSS_BOND_BONDS:
            newparam->cross_bb.r1e = old[0];
            newparam->cross_bb.r2e = old[1];
            newparam->cross_bb.krr = old[2];
            break;
        case F_CROSS_BOND_ANGLES:
            newparam->cross_ba.r1e = old[0];
            newparam->cross_ba.r2e = old[1];
            newparam->cross_ba.r3e = old[2];
            newparam->cross_ba.krt = old[3];
            break;
        case F_UREY_BRADLEY:
            newparam->u_b.thetaA  = old[0];
            newparam->u_b.kthetaA = old[1];
            newparam->u_b.r13A    = old[2];
            newparam->u_b.kUBA    = old[3];
            newparam->u_b.thetaB  = old[4];
            newparam->u_b.kthetaB = old[5];
            newparam->u_b.r13B    = old[6];
            newparam->u_b.kUBB    = old[7];
            break;
        case F_QUARTIC_ANGLES:
            newparam->qangle.theta = old[0];
            for (i = 0; i < 5; i++)
            {
                newparam->qangle.c[i] = old[i+1];
            }
            break;
        case F_LINEAR_ANGLES:
            newparam->linangle.aA    = old[0];
            newparam->linangle.klinA = old[1];
            newparam->linangle.aB    = old[2];
            newparam->linangle.klinB = old[3];
            break;
        case F_BONDS:
        case F_ANGLES:
        case F_HARMONIC:
        case F_IDIHS:
            newparam->harmonic.rA  = old[0];
            newparam->harmonic.krA = old[1];
            newparam->harmonic.rB  = old[2];
            newparam->harmonic.krB = old[3];
            break;
        case F_RESTRANGLES:
            newparam->harmonic.rA  = old[0];
            newparam->harmonic.krA = old[1];
            break;
        case F_MORSE:
            newparam->morse.b0A    = old[0];
            newparam->morse.cbA    = old[1];
            newparam->morse.betaA  = old[2];
            newparam->morse.b0B    = old[3];
            newparam->morse.cbB    = old[4];
            newparam->morse.betaB  = old[5];
            break;
        case F_CUBICBONDS:
            newparam->cubic.b0    = old[0];
            newparam->cubic.kb    = old[1];
            newparam->cubic.kcub  = old[2];
            break;
        case F_CONNBONDS:
            break;
        case F_POLARIZATION:
            newparam->polarize.alpha = old[0];
            break;
        case F_ANHARM_POL:
            newparam->anharm_polarize.alpha = old[0];
            newparam->anharm_polarize.drcut = old[1];
            newparam->anharm_polarize.khyp  = old[2];
            break;
        case F_WATER_POL:
            newparam->wpol.al_x   = old[0];
            newparam->wpol.al_y   = old[1];
            newparam->wpol.al_z   = old[2];
            newparam->wpol.rOH    = old[3];
            newparam->wpol.rHH    = old[4];
            newparam->wpol.rOD    = old[5];
            break;
        case F_THOLE_POL:
            newparam->thole.a      = old[0];
            newparam->thole.alpha1 = old[1];
            newparam->thole.alpha2 = old[2];
            if ((old[1] > 0) && (old[2] > 0))
            {
                newparam->thole.rfac = old[0]*pow(old[1]*old[2], -1.0/6.0);
            }
            else
            {
                newparam->thole.rfac = 1;
            }
            break;
        case F_BHAM:
            newparam->bham.a = old[0];
            newparam->bham.b = old[1];
            newparam->bham.c = old[2];
            break;
        case F_LJ14:
            set_ljparams(comb, reppow, old[0], old[1], &newparam->lj14.c6A, &newparam->lj14.c12A);
            set_ljparams(comb, reppow, old[2], old[3], &newparam->lj14.c6B, &newparam->lj14.c12B);
            break;
        case F_LJC14_Q:
            newparam->ljc14.fqq = old[0];
            newparam->ljc14.qi  = old[1];
            newparam->ljc14.qj  = old[2];
            set_ljparams(comb, reppow, old[3], old[4], &newparam->ljc14.c6, &newparam->ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            newparam->ljcnb.qi = old[0];
            newparam->ljcnb.qj = old[1];
            set_ljparams(comb, reppow, old[2], old[3], &newparam->ljcnb.c6, &newparam->ljcnb.c12);
            break;
        case F_LJ:
            set_ljparams(comb, reppow, old[0], old[1], &newparam->lj.c6, &newparam->lj.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            newparam->pdihs.phiA = old[0];
            newparam->pdihs.cpA  = old[1];

            /* Change 20100720: Amber occasionally uses negative multiplicities (mathematically OK),
             * so I have changed the lower limit to -99 /EL
             */
            newparam->pdihs.phiB = old[3];
            newparam->pdihs.cpB  = old[4];
            /* If both force constants are zero there is no interaction. Return -1 to signal
             * this entry should NOT be added.
             */
            if (fabs(newparam->pdihs.cpA) < GMX_REAL_MIN && fabs(newparam->pdihs.cpB) < GMX_REAL_MIN)
            {
                return -1;
            }

            newparam->pdihs.mult = round_check(old[2], -99, ftype, "multiplicity");

            break;
        case F_RESTRDIHS:
            newparam->pdihs.phiA = old[0];
            newparam->pdihs.cpA  = old[1];
            break;
        case F_POSRES:
            newparam->posres.fcA[XX]   = old[0];
            newparam->posres.fcA[YY]   = old[1];
            newparam->posres.fcA[ZZ]   = old[2];
            newparam->posres.fcB[XX]   = old[3];
            newparam->posres.fcB[YY]   = old[4];
            newparam->posres.fcB[ZZ]   = old[5];
            newparam->posres.pos0A[XX] = old[6];
            newparam->posres.pos0A[YY] = old[7];
            newparam->posres.pos0A[ZZ] = old[8];
            newparam->posres.pos0B[XX] = old[9];
            newparam->posres.pos0B[YY] = old[10];
            newparam->posres.pos0B[ZZ] = old[11];
            break;
        case F_FBPOSRES:
            newparam->fbposres.geom     = round_check(old[0], 0, ftype, "geometry");
            if (!(newparam->fbposres.geom > efbposresZERO && newparam->fbposres.geom < efbposresNR))
            {
                gmx_fatal(FARGS, "Invalid geometry for flat-bottomed position restraint.\n"
                          "Expected number between 1 and %d. Found %d\n", efbposresNR-1,
                          newparam->fbposres.geom);
            }
            newparam->fbposres.r        = old[1];
            newparam->fbposres.k        = old[2];
            newparam->fbposres.pos0[XX] = old[3];
            newparam->fbposres.pos0[YY] = old[4];
            newparam->fbposres.pos0[ZZ] = old[5];
            break;
        case F_DISRES:
            newparam->disres.label = round_check(old[0], 0, ftype, "label");
            newparam->disres.type  = round_check(old[1], 1, ftype, "type'");
            newparam->disres.low   = old[2];
            newparam->disres.up1   = old[3];
            newparam->disres.up2   = old[4];
            newparam->disres.kfac  = old[5];
            break;
        case F_ORIRES:
            newparam->orires.ex    = round_check(old[0], 1, ftype, "experiment") - 1;
            newparam->orires.label = round_check(old[1], 1, ftype, "label");
            newparam->orires.power = round_check(old[2], 0, ftype, "power");
            newparam->orires.c     = old[3];
            newparam->orires.obs   = old[4];
            newparam->orires.kfac  = old[5];
            break;
        case F_DIHRES:
            newparam->dihres.phiA   = old[0];
            newparam->dihres.dphiA  = old[1];
            newparam->dihres.kfacA  = old[2];
            newparam->dihres.phiB   = old[3];
            newparam->dihres.dphiB  = old[4];
            newparam->dihres.kfacB  = old[5];
            break;
        case F_RBDIHS:
            for (i = 0; (i < NR_RBDIHS); i++)
            {
                newparam->rbdihs.rbcA[i] = old[i];
                newparam->rbdihs.rbcB[i] = old[NR_RBDIHS+i];
            }
            break;
        case F_CBTDIHS:
            for (i = 0; (i < NR_CBTDIHS); i++)
            {
                newparam->cbtdihs.cbtcA[i] = old[i];
            }
            break;
        case F_FOURDIHS:
            /* Read the dihedral parameters to temporary arrays,
             * and convert them to the computationally faster
             * Ryckaert-Bellemans form.
             */
            /* Use conversion formula for OPLS to Ryckaert-Bellemans: */
            newparam->rbdihs.rbcA[0] = old[1]+0.5*(old[0]+old[2]);
            newparam->rbdihs.rbcA[1] = 0.5*(3.0*old[2]-old[0]);
            newparam->rbdihs.rbcA[2] = 4.0*old[3]-old[1];
            newparam->rbdihs.rbcA[3] = -2.0*old[2];
            newparam->rbdihs.rbcA[4] = -4.0*old[3];
            newparam->rbdihs.rbcA[5] = 0.0;

            newparam->rbdihs.rbcB[0] = old[NR_FOURDIHS+1]+0.5*(old[NR_FOURDIHS+0]+old[NR_FOURDIHS+2]);
            newparam->rbdihs.rbcB[1] = 0.5*(3.0*old[NR_FOURDIHS+2]-old[NR_FOURDIHS+0]);
            newparam->rbdihs.rbcB[2] = 4.0*old[NR_FOURDIHS+3]-old[NR_FOURDIHS+1];
            newparam->rbdihs.rbcB[3] = -2.0*old[NR_FOURDIHS+2];
            newparam->rbdihs.rbcB[4] = -4.0*old[NR_FOURDIHS+3];
            newparam->rbdihs.rbcB[5] = 0.0;
            break;
        case F_CONSTR:
        case F_CONSTRNC:
            newparam->constr.dA = old[0];
            newparam->constr.dB = old[1];
            break;
        case F_SETTLE:
            newparam->settle.doh = old[0];
            newparam->settle.dhh = old[1];
            break;
        case F_VSITE2:
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            newparam->vsite.a = old[0];
            newparam->vsite.b = old[1];
            newparam->vsite.c = old[2];
            newparam->vsite.d = old[3];
            newparam->vsite.e = old[4];
            newparam->vsite.f = old[5];
            break;
        case F_VSITE3FAD:
            newparam->vsite.a = old[1] * cos(DEG2RAD * old[0]);
            newparam->vsite.b = old[1] * sin(DEG2RAD * old[0]);
            newparam->vsite.c = old[2];
            newparam->vsite.d = old[3];
            newparam->vsite.e = old[4];
            newparam->vsite.f = old[5];
            break;
        case F_VSITEN:
            newparam->vsiten.n = round_check(old[0], 1, ftype, "number of atoms");
            newparam->vsiten.a = old[1];
            break;
        case F_CMAP:
            newparam->cmap.cmapA = old[0];
            newparam->cmap.cmapB = old[1];
            break;
        case F_GB12:
        case F_GB13:
        case F_GB14:
            newparam->gb.sar  = old[0];
            newparam->gb.st   = old[1];
            newparam->gb.pi   = old[2];
            newparam->gb.gbr  = old[3];
            newparam->gb.bmlt = old[4];
            break;
        default:
            gmx_fatal(FARGS, "unknown function type %d in %s line %d",
                      ftype, __FILE__, __LINE__);
    }
    return 0;
}

static int enter_params(gmx_ffparams_t *ffparams, t_functype ftype,
                        real forceparams[MAXFORCEPARAM], int comb, real reppow,
                        int start, gmx_bool bAppend)
{
    t_iparams newparam;
    int       type;
    int       rc;

    if ( (rc = assign_param(ftype, &newparam, forceparams, comb, reppow)) < 0)
    {
        /* -1 means this interaction is all-zero and should not be added */
        return rc;
    }

    if (!bAppend)
    {
        for (type = start; (type < ffparams->ntypes); type++)
        {
            if (ffparams->functype[type] == ftype)
            {
                if (F_GB13 == ftype)
                {
                    /* Occasionally, the way the 1-3 reference distance is
                     * computed can lead to non-binary-identical results, but I
                     * don't know why. */
                    if ((gmx_within_tol(newparam.gb.sar,  ffparams->iparams[type].gb.sar,  1e-6)) &&
                        (gmx_within_tol(newparam.gb.st,   ffparams->iparams[type].gb.st,   1e-6)) &&
                        (gmx_within_tol(newparam.gb.pi,   ffparams->iparams[type].gb.pi,   1e-6)) &&
                        (gmx_within_tol(newparam.gb.gbr,  ffparams->iparams[type].gb.gbr,  1e-6)) &&
                        (gmx_within_tol(newparam.gb.bmlt, ffparams->iparams[type].gb.bmlt, 1e-6)))
                    {
                        return type;
                    }
                }
                else
                {
                    if (memcmp(&newparam, &ffparams->iparams[type], (size_t)sizeof(newparam)) == 0)
                    {
                        return type;
                    }
                }
            }
        }
    }
    else
    {
        type = ffparams->ntypes;
    }
    if (debug)
    {
        fprintf(debug, "copying newparam to ffparams->iparams[%d] (ntypes=%d)\n",
                type, ffparams->ntypes);
    }
    memcpy(&ffparams->iparams[type], &newparam, (size_t)sizeof(newparam));

    ffparams->ntypes++;
    ffparams->functype[type] = ftype;

    return type;
}

static void append_interaction(t_ilist *ilist,
                               int type, int nral, atom_id a[MAXATOMLIST])
{
    int i, where1;

    where1     = ilist->nr;
    ilist->nr += nral+1;

    ilist->iatoms[where1++] = type;
    for (i = 0; (i < nral); i++)
    {
        ilist->iatoms[where1++] = a[i];
    }
}

static void enter_function(t_params *p, t_functype ftype, int comb, real reppow,
                           gmx_ffparams_t *ffparams, t_ilist *il,
                           int *maxtypes,
                           gmx_bool bNB, gmx_bool bAppend)
{
    int     k, type, nr, nral, delta, start;

    start = ffparams->ntypes;
    nr    = p->nr;

    for (k = 0; k < nr; k++)
    {
        if (*maxtypes <= ffparams->ntypes)
        {
            *maxtypes += 1000;
            srenew(ffparams->functype, *maxtypes);
            srenew(ffparams->iparams, *maxtypes);
            if (debug)
            {
                fprintf(debug, "%s, line %d: srenewed idef->functype and idef->iparams to %d\n",
                        __FILE__, __LINE__, *maxtypes);
            }
        }
        type = enter_params(ffparams, ftype, p->param[k].c, comb, reppow, start, bAppend);
        /* Type==-1 is used as a signal that this interaction is all-zero and should not be added. */
        if (!bNB && type >= 0)
        {
            nral  = NRAL(ftype);
            delta = nr*(nral+1);
            srenew(il->iatoms, il->nr+delta);
            append_interaction(il, type, nral, p->param[k].a);
        }
    }
}

void convert_params(int atnr, t_params nbtypes[],
                    t_molinfo *mi, t_molinfo *intermolecular_interactions,
                    int comb, double reppow, real fudgeQQ,
                    gmx_mtop_t *mtop)
{
    int             i, j, maxtypes, mt;
    unsigned long   flags;
    gmx_ffparams_t *ffp;
    gmx_moltype_t  *molt;
    t_params       *plist;

    maxtypes = 0;

    ffp           = &mtop->ffparams;
    ffp->ntypes   = 0;
    ffp->atnr     = atnr;
    ffp->functype = NULL;
    ffp->iparams  = NULL;
    ffp->reppow   = reppow;

    enter_function(&(nbtypes[F_LJ]),  (t_functype)F_LJ,    comb, reppow, ffp, NULL,
                   &maxtypes, TRUE, TRUE);
    enter_function(&(nbtypes[F_BHAM]), (t_functype)F_BHAM,  comb, reppow, ffp, NULL,
                   &maxtypes, TRUE, TRUE);

    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        molt = &mtop->moltype[mt];
        for (i = 0; (i < F_NRE); i++)
        {
            molt->ilist[i].nr     = 0;
            molt->ilist[i].iatoms = NULL;

            plist = mi[mt].plist;

            flags = interaction_function[i].flags;
            if ((i != F_LJ) && (i != F_BHAM) && ((flags & IF_BOND) ||
                                                 (flags & IF_VSITE) ||
                                                 (flags & IF_CONSTRAINT)))
            {
                enter_function(&(plist[i]), (t_functype)i, comb, reppow,
                               ffp, &molt->ilist[i],
                               &maxtypes, FALSE, (i == F_POSRES  || i == F_FBPOSRES));
            }
        }
    }

    mtop->bIntermolecularInteractions = FALSE;
    if (intermolecular_interactions != NULL)
    {
        /* Process the intermolecular interaction list */
        snew(mtop->intermolecular_ilist, F_NRE);

        for (i = 0; (i < F_NRE); i++)
        {
            mtop->intermolecular_ilist[i].nr     = 0;
            mtop->intermolecular_ilist[i].iatoms = NULL;

            plist = intermolecular_interactions->plist;

            if (plist[i].nr > 0)
            {
                flags = interaction_function[i].flags;
                /* For intermolecular interactions we (currently)
                 * only support potentials.
                 * Constraints and virtual sites would be possible,
                 * but require a lot of extra (bug-prone) code.
                 */
                if (!(flags & IF_BOND))
                {
                    gmx_fatal(FARGS, "The intermolecular_interaction section may only contain bonded potentials");
                }
                else if (NRAL(i) == 1) /* e.g. position restraints */
                {
                    gmx_fatal(FARGS, "Single atom interactions don't make sense in the intermolecular_interaction section, you can put them in the moleculetype section");
                }
                else if (flags & IF_CHEMBOND)
                {
                    gmx_fatal(FARGS, "The intermolecular_interaction can not contain chemically bonding interactions");
                }
                else
                {
                    enter_function(&(plist[i]), (t_functype)i, comb, reppow,
                                   ffp, &mtop->intermolecular_ilist[i],
                                   &maxtypes, FALSE, FALSE);

                    mtop->bIntermolecularInteractions = TRUE;
                }
            }
        }

        if (!mtop->bIntermolecularInteractions)
        {
            sfree(mtop->intermolecular_ilist);
        }
    }

    if (debug)
    {
        fprintf(debug, "%s, line %d: There are %d functypes in idef\n",
                __FILE__, __LINE__, ffp->ntypes);
    }

    ffp->fudgeQQ = fudgeQQ;
}
