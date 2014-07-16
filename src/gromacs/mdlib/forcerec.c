/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include <assert.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "types/commrec.h"
#include "vec.h"
#include "gromacs/math/utilities.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "force.h"
#include "tables.h"
#include "nonbonded.h"
#include "invblock.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "mshift.h"
#include "txtdump.h"
#include "coulomb.h"
#include "md_support.h"
#include "md_logging.h"
#include "domdec.h"
#include "qmmm.h"
#include "copyrite.h"
#include "mtop_util.h"
#include "nbnxn_simd.h"
#include "nbnxn_search.h"
#include "nbnxn_atomdata.h"
#include "nbnxn_consts.h"
#include "gmx_omp_nthreads.h"
#include "gmx_detect_hardware.h"
#include "inputrec.h"

#include "types/nbnxn_cuda_types_ext.h"
#include "gpu_utils.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

t_forcerec *mk_forcerec(void)
{
    t_forcerec *fr;

    snew(fr, 1);

    return fr;
}

#ifdef DEBUG
static void pr_nbfp(FILE *fp, real *nbfp, gmx_bool bBHAM, int atnr)
{
    int i, j;

    for (i = 0; (i < atnr); i++)
    {
        for (j = 0; (j < atnr); j++)
        {
            fprintf(fp, "%2d - %2d", i, j);
            if (bBHAM)
            {
                fprintf(fp, "  a=%10g, b=%10g, c=%10g\n", BHAMA(nbfp, atnr, i, j),
                        BHAMB(nbfp, atnr, i, j), BHAMC(nbfp, atnr, i, j)/6.0);
            }
            else
            {
                fprintf(fp, "  c6=%10g, c12=%10g\n", C6(nbfp, atnr, i, j)/6.0,
                        C12(nbfp, atnr, i, j)/12.0);
            }
        }
    }
}
#endif

static real *mk_nbfp(const gmx_ffparams_t *idef, gmx_bool bBHAM)
{
    real *nbfp;
    int   i, j, k, atnr;

    atnr = idef->atnr;
    if (bBHAM)
    {
        snew(nbfp, 3*atnr*atnr);
        for (i = k = 0; (i < atnr); i++)
        {
            for (j = 0; (j < atnr); j++, k++)
            {
                BHAMA(nbfp, atnr, i, j) = idef->iparams[k].bham.a;
                BHAMB(nbfp, atnr, i, j) = idef->iparams[k].bham.b;
                /* nbfp now includes the 6.0 derivative prefactor */
                BHAMC(nbfp, atnr, i, j) = idef->iparams[k].bham.c*6.0;
            }
        }
    }
    else
    {
        snew(nbfp, 2*atnr*atnr);
        for (i = k = 0; (i < atnr); i++)
        {
            for (j = 0; (j < atnr); j++, k++)
            {
                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                C6(nbfp, atnr, i, j)   = idef->iparams[k].lj.c6*6.0;
                C12(nbfp, atnr, i, j)  = idef->iparams[k].lj.c12*12.0;
            }
        }
    }

    return nbfp;
}

static real *make_ljpme_c6grid(const gmx_ffparams_t *idef, t_forcerec *fr)
{
    int   i, j, k, atnr;
    real  c6, c6i, c6j, c12i, c12j, epsi, epsj, sigmai, sigmaj;
    real *grid;

    /* For LJ-PME simulations, we correct the energies with the reciprocal space
     * inside of the cut-off. To do this the non-bonded kernels needs to have
     * access to the C6-values used on the reciprocal grid in pme.c
     */

    atnr = idef->atnr;
    snew(grid, 2*atnr*atnr);
    for (i = k = 0; (i < atnr); i++)
    {
        for (j = 0; (j < atnr); j++, k++)
        {
            c6i  = idef->iparams[i*(atnr+1)].lj.c6;
            c12i = idef->iparams[i*(atnr+1)].lj.c12;
            c6j  = idef->iparams[j*(atnr+1)].lj.c6;
            c12j = idef->iparams[j*(atnr+1)].lj.c12;
            c6   = sqrt(c6i * c6j);
            if (fr->ljpme_combination_rule == eljpmeLB
                && !gmx_numzero(c6) && !gmx_numzero(c12i) && !gmx_numzero(c12j))
            {
                sigmai = pow(c12i / c6i, 1.0/6.0);
                sigmaj = pow(c12j / c6j, 1.0/6.0);
                epsi   = c6i * c6i / c12i;
                epsj   = c6j * c6j / c12j;
                c6     = sqrt(epsi * epsj) * pow(0.5*(sigmai+sigmaj), 6);
            }
            /* Store the elements at the same relative positions as C6 in nbfp in order
             * to simplify access in the kernels
             */
            grid[2*(atnr*i+j)] = c6*6.0;
        }
    }
    return grid;
}

static real *mk_nbfp_combination_rule(const gmx_ffparams_t *idef, int comb_rule)
{
    real *nbfp;
    int   i, j, k, atnr;
    real  c6i, c6j, c12i, c12j, epsi, epsj, sigmai, sigmaj;
    real  c6, c12;

    atnr = idef->atnr;
    snew(nbfp, 2*atnr*atnr);
    for (i = 0; i < atnr; ++i)
    {
        for (j = 0; j < atnr; ++j)
        {
            c6i  = idef->iparams[i*(atnr+1)].lj.c6;
            c12i = idef->iparams[i*(atnr+1)].lj.c12;
            c6j  = idef->iparams[j*(atnr+1)].lj.c6;
            c12j = idef->iparams[j*(atnr+1)].lj.c12;
            c6   = sqrt(c6i  * c6j);
            c12  = sqrt(c12i * c12j);
            if (comb_rule == eCOMB_ARITHMETIC
                && !gmx_numzero(c6) && !gmx_numzero(c12))
            {
                sigmai = pow(c12i / c6i, 1.0/6.0);
                sigmaj = pow(c12j / c6j, 1.0/6.0);
                epsi   = c6i * c6i / c12i;
                epsj   = c6j * c6j / c12j;
                c6     = sqrt(epsi * epsj) * pow(0.5*(sigmai+sigmaj), 6);
                c12    = sqrt(epsi * epsj) * pow(0.5*(sigmai+sigmaj), 12);
            }
            C6(nbfp, atnr, i, j)   = c6*6.0;
            C12(nbfp, atnr, i, j)  = c12*12.0;
        }
    }
    return nbfp;
}

/* This routine sets fr->solvent_opt to the most common solvent in the
 * system, e.g. esolSPC or esolTIP4P. It will also mark each charge group in
 * the fr->solvent_type array with the correct type (or esolNO).
 *
 * Charge groups that fulfill the conditions but are not identical to the
 * most common one will be marked as esolNO in the solvent_type array.
 *
 * TIP3p is identical to SPC for these purposes, so we call it
 * SPC in the arrays (Apologies to Bill Jorgensen ;-)
 *
 * NOTE: QM particle should not
 * become an optimized solvent. Not even if there is only one charge
 * group in the Qm
 */

typedef struct
{
    int    model;
    int    count;
    int    vdwtype[4];
    real   charge[4];
} solvent_parameters_t;

static void
check_solvent_cg(const gmx_moltype_t    *molt,
                 int                     cg0,
                 int                     nmol,
                 const unsigned char    *qm_grpnr,
                 const t_grps           *qm_grps,
                 t_forcerec   *          fr,
                 int                    *n_solvent_parameters,
                 solvent_parameters_t  **solvent_parameters_p,
                 int                     cginfo,
                 int                    *cg_sp)
{
    const t_blocka       *excl;
    t_atom               *atom;
    int                   j, k;
    int                   j0, j1, nj;
    gmx_bool              perturbed;
    gmx_bool              has_vdw[4];
    gmx_bool              match;
    real                  tmp_charge[4]  = { 0.0 }; /* init to zero to make gcc4.8 happy */
    int                   tmp_vdwtype[4] = { 0 };   /* init to zero to make gcc4.8 happy */
    int                   tjA;
    gmx_bool              qm;
    solvent_parameters_t *solvent_parameters;

    /* We use a list with parameters for each solvent type.
     * Every time we discover a new molecule that fulfills the basic
     * conditions for a solvent we compare with the previous entries
     * in these lists. If the parameters are the same we just increment
     * the counter for that type, and otherwise we create a new type
     * based on the current molecule.
     *
     * Once we've finished going through all molecules we check which
     * solvent is most common, and mark all those molecules while we
     * clear the flag on all others.
     */

    solvent_parameters = *solvent_parameters_p;

    /* Mark the cg first as non optimized */
    *cg_sp = -1;

    /* Check if this cg has no exclusions with atoms in other charge groups
     * and all atoms inside the charge group excluded.
     * We only have 3 or 4 atom solvent loops.
     */
    if (GET_CGINFO_EXCL_INTER(cginfo) ||
        !GET_CGINFO_EXCL_INTRA(cginfo))
    {
        return;
    }

    /* Get the indices of the first atom in this charge group */
    j0     = molt->cgs.index[cg0];
    j1     = molt->cgs.index[cg0+1];

    /* Number of atoms in our molecule */
    nj     = j1 - j0;

    if (debug)
    {
        fprintf(debug,
                "Moltype '%s': there are %d atoms in this charge group\n",
                *molt->name, nj);
    }

    /* Check if it could be an SPC (3 atoms) or TIP4p (4) water,
     * otherwise skip it.
     */
    if (nj < 3 || nj > 4)
    {
        return;
    }

    /* Check if we are doing QM on this group */
    qm = FALSE;
    if (qm_grpnr != NULL)
    {
        for (j = j0; j < j1 && !qm; j++)
        {
            qm = (qm_grpnr[j] < qm_grps->nr - 1);
        }
    }
    /* Cannot use solvent optimization with QM */
    if (qm)
    {
        return;
    }

    atom = molt->atoms.atom;

    /* Still looks like a solvent, time to check parameters */

    /* If it is perturbed (free energy) we can't use the solvent loops,
     * so then we just skip to the next molecule.
     */
    perturbed = FALSE;

    for (j = j0; j < j1 && !perturbed; j++)
    {
        perturbed = PERTURBED(atom[j]);
    }

    if (perturbed)
    {
        return;
    }

    /* Now it's only a question if the VdW and charge parameters
     * are OK. Before doing the check we compare and see if they are
     * identical to a possible previous solvent type.
     * First we assign the current types and charges.
     */
    for (j = 0; j < nj; j++)
    {
        tmp_vdwtype[j] = atom[j0+j].type;
        tmp_charge[j]  = atom[j0+j].q;
    }

    /* Does it match any previous solvent type? */
    for (k = 0; k < *n_solvent_parameters; k++)
    {
        match = TRUE;


        /* We can only match SPC with 3 atoms and TIP4p with 4 atoms */
        if ( (solvent_parameters[k].model == esolSPC   && nj != 3)  ||
             (solvent_parameters[k].model == esolTIP4P && nj != 4) )
        {
            match = FALSE;
        }

        /* Check that types & charges match for all atoms in molecule */
        for (j = 0; j < nj && match == TRUE; j++)
        {
            if (tmp_vdwtype[j] != solvent_parameters[k].vdwtype[j])
            {
                match = FALSE;
            }
            if (tmp_charge[j] != solvent_parameters[k].charge[j])
            {
                match = FALSE;
            }
        }
        if (match == TRUE)
        {
            /* Congratulations! We have a matched solvent.
             * Flag it with this type for later processing.
             */
            *cg_sp = k;
            solvent_parameters[k].count += nmol;

            /* We are done with this charge group */
            return;
        }
    }

    /* If we get here, we have a tentative new solvent type.
     * Before we add it we must check that it fulfills the requirements
     * of the solvent optimized loops. First determine which atoms have
     * VdW interactions.
     */
    for (j = 0; j < nj; j++)
    {
        has_vdw[j] = FALSE;
        tjA        = tmp_vdwtype[j];

        /* Go through all other tpes and see if any have non-zero
         * VdW parameters when combined with this one.
         */
        for (k = 0; k < fr->ntype && (has_vdw[j] == FALSE); k++)
        {
            /* We already checked that the atoms weren't perturbed,
             * so we only need to check state A now.
             */
            if (fr->bBHAM)
            {
                has_vdw[j] = (has_vdw[j] ||
                              (BHAMA(fr->nbfp, fr->ntype, tjA, k) != 0.0) ||
                              (BHAMB(fr->nbfp, fr->ntype, tjA, k) != 0.0) ||
                              (BHAMC(fr->nbfp, fr->ntype, tjA, k) != 0.0));
            }
            else
            {
                /* Standard LJ */
                has_vdw[j] = (has_vdw[j] ||
                              (C6(fr->nbfp, fr->ntype, tjA, k)  != 0.0) ||
                              (C12(fr->nbfp, fr->ntype, tjA, k) != 0.0));
            }
        }
    }

    /* Now we know all we need to make the final check and assignment. */
    if (nj == 3)
    {
        /* So, is it an SPC?
         * For this we require thatn all atoms have charge,
         * the charges on atom 2 & 3 should be the same, and only
         * atom 1 might have VdW.
         */
        if (has_vdw[1] == FALSE &&
            has_vdw[2] == FALSE &&
            tmp_charge[0]  != 0 &&
            tmp_charge[1]  != 0 &&
            tmp_charge[2]  == tmp_charge[1])
        {
            srenew(solvent_parameters, *n_solvent_parameters+1);
            solvent_parameters[*n_solvent_parameters].model = esolSPC;
            solvent_parameters[*n_solvent_parameters].count = nmol;
            for (k = 0; k < 3; k++)
            {
                solvent_parameters[*n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                solvent_parameters[*n_solvent_parameters].charge[k]  = tmp_charge[k];
            }

            *cg_sp = *n_solvent_parameters;
            (*n_solvent_parameters)++;
        }
    }
    else if (nj == 4)
    {
        /* Or could it be a TIP4P?
         * For this we require thatn atoms 2,3,4 have charge, but not atom 1.
         * Only atom 1 mght have VdW.
         */
        if (has_vdw[1] == FALSE &&
            has_vdw[2] == FALSE &&
            has_vdw[3] == FALSE &&
            tmp_charge[0]  == 0 &&
            tmp_charge[1]  != 0 &&
            tmp_charge[2]  == tmp_charge[1] &&
            tmp_charge[3]  != 0)
        {
            srenew(solvent_parameters, *n_solvent_parameters+1);
            solvent_parameters[*n_solvent_parameters].model = esolTIP4P;
            solvent_parameters[*n_solvent_parameters].count = nmol;
            for (k = 0; k < 4; k++)
            {
                solvent_parameters[*n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                solvent_parameters[*n_solvent_parameters].charge[k]  = tmp_charge[k];
            }

            *cg_sp = *n_solvent_parameters;
            (*n_solvent_parameters)++;
        }
    }

    *solvent_parameters_p = solvent_parameters;
}

static void
check_solvent(FILE  *                fp,
              const gmx_mtop_t  *    mtop,
              t_forcerec  *          fr,
              cginfo_mb_t           *cginfo_mb)
{
    const t_block     *   cgs;
    const t_block     *   mols;
    const gmx_moltype_t  *molt;
    int                   mb, mol, cg_mol, at_offset, cg_offset, am, cgm, i, nmol_ch, nmol;
    int                   n_solvent_parameters;
    solvent_parameters_t *solvent_parameters;
    int                 **cg_sp;
    int                   bestsp, bestsol;

    if (debug)
    {
        fprintf(debug, "Going to determine what solvent types we have.\n");
    }

    mols = &mtop->mols;

    n_solvent_parameters = 0;
    solvent_parameters   = NULL;
    /* Allocate temporary array for solvent type */
    snew(cg_sp, mtop->nmolblock);

    cg_offset = 0;
    at_offset = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        cgs  = &molt->cgs;
        /* Here we have to loop over all individual molecules
         * because we need to check for QMMM particles.
         */
        snew(cg_sp[mb], cginfo_mb[mb].cg_mod);
        nmol_ch = cginfo_mb[mb].cg_mod/cgs->nr;
        nmol    = mtop->molblock[mb].nmol/nmol_ch;
        for (mol = 0; mol < nmol_ch; mol++)
        {
            cgm = mol*cgs->nr;
            am  = mol*cgs->index[cgs->nr];
            for (cg_mol = 0; cg_mol < cgs->nr; cg_mol++)
            {
                check_solvent_cg(molt, cg_mol, nmol,
                                 mtop->groups.grpnr[egcQMMM] ?
                                 mtop->groups.grpnr[egcQMMM]+at_offset+am : 0,
                                 &mtop->groups.grps[egcQMMM],
                                 fr,
                                 &n_solvent_parameters, &solvent_parameters,
                                 cginfo_mb[mb].cginfo[cgm+cg_mol],
                                 &cg_sp[mb][cgm+cg_mol]);
            }
        }
        cg_offset += cgs->nr;
        at_offset += cgs->index[cgs->nr];
    }

    /* Puh! We finished going through all charge groups.
     * Now find the most common solvent model.
     */

    /* Most common solvent this far */
    bestsp = -2;
    for (i = 0; i < n_solvent_parameters; i++)
    {
        if (bestsp == -2 ||
            solvent_parameters[i].count > solvent_parameters[bestsp].count)
        {
            bestsp = i;
        }
    }

    if (bestsp >= 0)
    {
        bestsol = solvent_parameters[bestsp].model;
    }
    else
    {
        bestsol = esolNO;
    }

#ifdef DISABLE_WATER_NLIST
    bestsol = esolNO;
#endif

    fr->nWatMol = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        cgs  = &mtop->moltype[mtop->molblock[mb].type].cgs;
        nmol = (mtop->molblock[mb].nmol*cgs->nr)/cginfo_mb[mb].cg_mod;
        for (i = 0; i < cginfo_mb[mb].cg_mod; i++)
        {
            if (cg_sp[mb][i] == bestsp)
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[i], bestsol);
                fr->nWatMol += nmol;
            }
            else
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[i], esolNO);
            }
        }
        sfree(cg_sp[mb]);
    }
    sfree(cg_sp);

    if (bestsol != esolNO && fp != NULL)
    {
        fprintf(fp, "\nEnabling %s-like water optimization for %d molecules.\n\n",
                esol_names[bestsol],
                solvent_parameters[bestsp].count);
    }

    sfree(solvent_parameters);
    fr->solvent_opt = bestsol;
}

enum {
    acNONE = 0, acCONSTRAINT, acSETTLE
};

static cginfo_mb_t *init_cginfo_mb(FILE *fplog, const gmx_mtop_t *mtop,
                                   t_forcerec *fr, gmx_bool bNoSolvOpt,
                                   gmx_bool *bFEP_NonBonded,
                                   gmx_bool *bExcl_IntraCGAll_InterCGNone)
{
    const t_block        *cgs;
    const t_blocka       *excl;
    const gmx_moltype_t  *molt;
    const gmx_molblock_t *molb;
    cginfo_mb_t          *cginfo_mb;
    gmx_bool             *type_VDW;
    int                  *cginfo;
    int                   cg_offset, a_offset, cgm, am;
    int                   mb, m, ncg_tot, cg, a0, a1, gid, ai, j, aj, excl_nalloc;
    int                  *a_con;
    int                   ftype;
    int                   ia;
    gmx_bool              bId, *bExcl, bExclIntraAll, bExclInter, bHaveVDW, bHaveQ, bHavePerturbedAtoms;

    ncg_tot = ncg_mtop(mtop);
    snew(cginfo_mb, mtop->nmolblock);

    snew(type_VDW, fr->ntype);
    for (ai = 0; ai < fr->ntype; ai++)
    {
        type_VDW[ai] = FALSE;
        for (j = 0; j < fr->ntype; j++)
        {
            type_VDW[ai] = type_VDW[ai] ||
                fr->bBHAM ||
                C6(fr->nbfp, fr->ntype, ai, j) != 0 ||
                C12(fr->nbfp, fr->ntype, ai, j) != 0;
        }
    }

    *bFEP_NonBonded               = FALSE;
    *bExcl_IntraCGAll_InterCGNone = TRUE;

    excl_nalloc = 10;
    snew(bExcl, excl_nalloc);
    cg_offset = 0;
    a_offset  = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        cgs  = &molt->cgs;
        excl = &molt->excls;

        /* Check if the cginfo is identical for all molecules in this block.
         * If so, we only need an array of the size of one molecule.
         * Otherwise we make an array of #mol times #cgs per molecule.
         */
        bId = TRUE;
        am  = 0;
        for (m = 0; m < molb->nmol; m++)
        {
            am = m*cgs->index[cgs->nr];
            for (cg = 0; cg < cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                if (ggrpnr(&mtop->groups, egcENER, a_offset+am+a0) !=
                    ggrpnr(&mtop->groups, egcENER, a_offset   +a0))
                {
                    bId = FALSE;
                }
                if (mtop->groups.grpnr[egcQMMM] != NULL)
                {
                    for (ai = a0; ai < a1; ai++)
                    {
                        if (mtop->groups.grpnr[egcQMMM][a_offset+am+ai] !=
                            mtop->groups.grpnr[egcQMMM][a_offset   +ai])
                        {
                            bId = FALSE;
                        }
                    }
                }
            }
        }

        cginfo_mb[mb].cg_start = cg_offset;
        cginfo_mb[mb].cg_end   = cg_offset + molb->nmol*cgs->nr;
        cginfo_mb[mb].cg_mod   = (bId ? 1 : molb->nmol)*cgs->nr;
        snew(cginfo_mb[mb].cginfo, cginfo_mb[mb].cg_mod);
        cginfo = cginfo_mb[mb].cginfo;

        /* Set constraints flags for constrained atoms */
        snew(a_con, molt->atoms.nr);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_CONSTRAINT)
            {
                int nral;

                nral = NRAL(ftype);
                for (ia = 0; ia < molt->ilist[ftype].nr; ia += 1+nral)
                {
                    int a;

                    for (a = 0; a < nral; a++)
                    {
                        a_con[molt->ilist[ftype].iatoms[ia+1+a]] =
                            (ftype == F_SETTLE ? acSETTLE : acCONSTRAINT);
                    }
                }
            }
        }

        for (m = 0; m < (bId ? 1 : molb->nmol); m++)
        {
            cgm = m*cgs->nr;
            am  = m*cgs->index[cgs->nr];
            for (cg = 0; cg < cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];

                /* Store the energy group in cginfo */
                gid = ggrpnr(&mtop->groups, egcENER, a_offset+am+a0);
                SET_CGINFO_GID(cginfo[cgm+cg], gid);

                /* Check the intra/inter charge group exclusions */
                if (a1-a0 > excl_nalloc)
                {
                    excl_nalloc = a1 - a0;
                    srenew(bExcl, excl_nalloc);
                }
                /* bExclIntraAll: all intra cg interactions excluded
                 * bExclInter:    any inter cg interactions excluded
                 */
                bExclIntraAll       = TRUE;
                bExclInter          = FALSE;
                bHaveVDW            = FALSE;
                bHaveQ              = FALSE;
                bHavePerturbedAtoms = FALSE;
                for (ai = a0; ai < a1; ai++)
                {
                    /* Check VDW and electrostatic interactions */
                    bHaveVDW = bHaveVDW || (type_VDW[molt->atoms.atom[ai].type] ||
                                            type_VDW[molt->atoms.atom[ai].typeB]);
                    bHaveQ  = bHaveQ    || (molt->atoms.atom[ai].q != 0 ||
                                            molt->atoms.atom[ai].qB != 0);

                    bHavePerturbedAtoms = bHavePerturbedAtoms || (PERTURBED(molt->atoms.atom[ai]) != 0);

                    /* Clear the exclusion list for atom ai */
                    for (aj = a0; aj < a1; aj++)
                    {
                        bExcl[aj-a0] = FALSE;
                    }
                    /* Loop over all the exclusions of atom ai */
                    for (j = excl->index[ai]; j < excl->index[ai+1]; j++)
                    {
                        aj = excl->a[j];
                        if (aj < a0 || aj >= a1)
                        {
                            bExclInter = TRUE;
                        }
                        else
                        {
                            bExcl[aj-a0] = TRUE;
                        }
                    }
                    /* Check if ai excludes a0 to a1 */
                    for (aj = a0; aj < a1; aj++)
                    {
                        if (!bExcl[aj-a0])
                        {
                            bExclIntraAll = FALSE;
                        }
                    }

                    switch (a_con[ai])
                    {
                        case acCONSTRAINT:
                            SET_CGINFO_CONSTR(cginfo[cgm+cg]);
                            break;
                        case acSETTLE:
                            SET_CGINFO_SETTLE(cginfo[cgm+cg]);
                            break;
                        default:
                            break;
                    }
                }
                if (bExclIntraAll)
                {
                    SET_CGINFO_EXCL_INTRA(cginfo[cgm+cg]);
                }
                if (bExclInter)
                {
                    SET_CGINFO_EXCL_INTER(cginfo[cgm+cg]);
                }
                if (a1 - a0 > MAX_CHARGEGROUP_SIZE)
                {
                    /* The size in cginfo is currently only read with DD */
                    gmx_fatal(FARGS, "A charge group has size %d which is larger than the limit of %d atoms", a1-a0, MAX_CHARGEGROUP_SIZE);
                }
                if (bHaveVDW)
                {
                    SET_CGINFO_HAS_VDW(cginfo[cgm+cg]);
                }
                if (bHaveQ)
                {
                    SET_CGINFO_HAS_Q(cginfo[cgm+cg]);
                }
                if (bHavePerturbedAtoms && fr->efep != efepNO)
                {
                    SET_CGINFO_FEP(cginfo[cgm+cg]);
                    *bFEP_NonBonded = TRUE;
                }
                /* Store the charge group size */
                SET_CGINFO_NATOMS(cginfo[cgm+cg], a1-a0);

                if (!bExclIntraAll || bExclInter)
                {
                    *bExcl_IntraCGAll_InterCGNone = FALSE;
                }
            }
        }

        sfree(a_con);

        cg_offset += molb->nmol*cgs->nr;
        a_offset  += molb->nmol*cgs->index[cgs->nr];
    }
    sfree(bExcl);

    /* the solvent optimizer is called after the QM is initialized,
     * because we don't want to have the QM subsystemto become an
     * optimized solvent
     */

    check_solvent(fplog, mtop, fr, cginfo_mb);

    if (getenv("GMX_NO_SOLV_OPT"))
    {
        if (fplog)
        {
            fprintf(fplog, "Found environment variable GMX_NO_SOLV_OPT.\n"
                    "Disabling all solvent optimization\n");
        }
        fr->solvent_opt = esolNO;
    }
    if (bNoSolvOpt)
    {
        fr->solvent_opt = esolNO;
    }
    if (!fr->solvent_opt)
    {
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            for (cg = 0; cg < cginfo_mb[mb].cg_mod; cg++)
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[cg], esolNO);
            }
        }
    }

    return cginfo_mb;
}

static int *cginfo_expand(int nmb, cginfo_mb_t *cgi_mb)
{
    int  ncg, mb, cg;
    int *cginfo;

    ncg = cgi_mb[nmb-1].cg_end;
    snew(cginfo, ncg);
    mb = 0;
    for (cg = 0; cg < ncg; cg++)
    {
        while (cg >= cgi_mb[mb].cg_end)
        {
            mb++;
        }
        cginfo[cg] =
            cgi_mb[mb].cginfo[(cg - cgi_mb[mb].cg_start) % cgi_mb[mb].cg_mod];
    }

    return cginfo;
}

static void set_chargesum(FILE *log, t_forcerec *fr, const gmx_mtop_t *mtop)
{
    /*This now calculates sum for q and c6*/
    double         qsum, q2sum, q, c6sum, c6;
    int            mb, nmol, i;
    const t_atoms *atoms;

    qsum   = 0;
    q2sum  = 0;
    c6sum  = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        nmol  = mtop->molblock[mb].nmol;
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        for (i = 0; i < atoms->nr; i++)
        {
            q       = atoms->atom[i].q;
            qsum   += nmol*q;
            q2sum  += nmol*q*q;
            c6      = mtop->ffparams.iparams[atoms->atom[i].type*(mtop->ffparams.atnr+1)].lj.c6;
            c6sum  += nmol*c6;
        }
    }
    fr->qsum[0]   = qsum;
    fr->q2sum[0]  = q2sum;
    fr->c6sum[0]  = c6sum;

    if (fr->efep != efepNO)
    {
        qsum   = 0;
        q2sum  = 0;
        c6sum  = 0;
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            nmol  = mtop->molblock[mb].nmol;
            atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
            for (i = 0; i < atoms->nr; i++)
            {
                q       = atoms->atom[i].qB;
                qsum   += nmol*q;
                q2sum  += nmol*q*q;
                c6      = mtop->ffparams.iparams[atoms->atom[i].typeB*(mtop->ffparams.atnr+1)].lj.c6;
                c6sum  += nmol*c6;
            }
            fr->qsum[1]   = qsum;
            fr->q2sum[1]  = q2sum;
            fr->c6sum[1]  = c6sum;
        }
    }
    else
    {
        fr->qsum[1]   = fr->qsum[0];
        fr->q2sum[1]  = fr->q2sum[0];
        fr->c6sum[1]  = fr->c6sum[0];
    }
    if (log)
    {
        if (fr->efep == efepNO)
        {
            fprintf(log, "System total charge: %.3f\n", fr->qsum[0]);
        }
        else
        {
            fprintf(log, "System total charge, top. A: %.3f top. B: %.3f\n",
                    fr->qsum[0], fr->qsum[1]);
        }
    }
}

void update_forcerec(t_forcerec *fr, matrix box)
{
    if (fr->eeltype == eelGRF)
    {
        calc_rffac(NULL, fr->eeltype, fr->epsilon_r, fr->epsilon_rf,
                   fr->rcoulomb, fr->temp, fr->zsquare, box,
                   &fr->kappa, &fr->k_rf, &fr->c_rf);
    }
}

void set_avcsixtwelve(FILE *fplog, t_forcerec *fr, const gmx_mtop_t *mtop)
{
    const t_atoms  *atoms, *atoms_tpi;
    const t_blocka *excl;
    int             mb, nmol, nmolc, i, j, tpi, tpj, j1, j2, k, n, nexcl, q;
    gmx_int64_t     npair, npair_ij, tmpi, tmpj;
    double          csix, ctwelve;
    int             ntp, *typecount;
    gmx_bool        bBHAM;
    real           *nbfp;
    real           *nbfp_comb = NULL;

    ntp   = fr->ntype;
    bBHAM = fr->bBHAM;
    nbfp  = fr->nbfp;

    /* For LJ-PME, we want to correct for the difference between the
     * actual C6 values and the C6 values used by the LJ-PME based on
     * combination rules. */

    if (EVDW_PME(fr->vdwtype))
    {
        nbfp_comb = mk_nbfp_combination_rule(&mtop->ffparams,
                                             (fr->ljpme_combination_rule == eljpmeLB) ? eCOMB_ARITHMETIC : eCOMB_GEOMETRIC);
        for (tpi = 0; tpi < ntp; ++tpi)
        {
            for (tpj = 0; tpj < ntp; ++tpj)
            {
                C6(nbfp_comb, ntp, tpi, tpj) =
                    C6(nbfp, ntp, tpi, tpj) - C6(nbfp_comb, ntp, tpi, tpj);
                C12(nbfp_comb, ntp, tpi, tpj) = C12(nbfp, ntp, tpi, tpj);
            }
        }
        nbfp = nbfp_comb;
    }
    for (q = 0; q < (fr->efep == efepNO ? 1 : 2); q++)
    {
        csix    = 0;
        ctwelve = 0;
        npair   = 0;
        nexcl   = 0;
        if (!fr->n_tpi)
        {
            /* Count the types so we avoid natoms^2 operations */
            snew(typecount, ntp);
            gmx_mtop_count_atomtypes(mtop, q, typecount);

            for (tpi = 0; tpi < ntp; tpi++)
            {
                for (tpj = tpi; tpj < ntp; tpj++)
                {
                    tmpi = typecount[tpi];
                    tmpj = typecount[tpj];
                    if (tpi != tpj)
                    {
                        npair_ij = tmpi*tmpj;
                    }
                    else
                    {
                        npair_ij = tmpi*(tmpi - 1)/2;
                    }
                    if (bBHAM)
                    {
                        /* nbfp now includes the 6.0 derivative prefactor */
                        csix    += npair_ij*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                    }
                    else
                    {
                        /* nbfp now includes the 6.0/12.0 derivative prefactors */
                        csix    += npair_ij*   C6(nbfp, ntp, tpi, tpj)/6.0;
                        ctwelve += npair_ij*  C12(nbfp, ntp, tpi, tpj)/12.0;
                    }
                    npair += npair_ij;
                }
            }
            sfree(typecount);
            /* Subtract the excluded pairs.
             * The main reason for substracting exclusions is that in some cases
             * some combinations might never occur and the parameters could have
             * any value. These unused values should not influence the dispersion
             * correction.
             */
            for (mb = 0; mb < mtop->nmolblock; mb++)
            {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                excl  = &mtop->moltype[mtop->molblock[mb].type].excls;
                for (i = 0; (i < atoms->nr); i++)
                {
                    if (q == 0)
                    {
                        tpi = atoms->atom[i].type;
                    }
                    else
                    {
                        tpi = atoms->atom[i].typeB;
                    }
                    j1  = excl->index[i];
                    j2  = excl->index[i+1];
                    for (j = j1; j < j2; j++)
                    {
                        k = excl->a[j];
                        if (k > i)
                        {
                            if (q == 0)
                            {
                                tpj = atoms->atom[k].type;
                            }
                            else
                            {
                                tpj = atoms->atom[k].typeB;
                            }
                            if (bBHAM)
                            {
                                /* nbfp now includes the 6.0 derivative prefactor */
                                csix -= nmol*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                            }
                            else
                            {
                                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                                csix    -= nmol*C6 (nbfp, ntp, tpi, tpj)/6.0;
                                ctwelve -= nmol*C12(nbfp, ntp, tpi, tpj)/12.0;
                            }
                            nexcl += nmol;
                        }
                    }
                }
            }
        }
        else
        {
            /* Only correct for the interaction of the test particle
             * with the rest of the system.
             */
            atoms_tpi =
                &mtop->moltype[mtop->molblock[mtop->nmolblock-1].type].atoms;

            npair = 0;
            for (mb = 0; mb < mtop->nmolblock; mb++)
            {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                for (j = 0; j < atoms->nr; j++)
                {
                    nmolc = nmol;
                    /* Remove the interaction of the test charge group
                     * with itself.
                     */
                    if (mb == mtop->nmolblock-1)
                    {
                        nmolc--;

                        if (mb == 0 && nmol == 1)
                        {
                            gmx_fatal(FARGS, "Old format tpr with TPI, please generate a new tpr file");
                        }
                    }
                    if (q == 0)
                    {
                        tpj = atoms->atom[j].type;
                    }
                    else
                    {
                        tpj = atoms->atom[j].typeB;
                    }
                    for (i = 0; i < fr->n_tpi; i++)
                    {
                        if (q == 0)
                        {
                            tpi = atoms_tpi->atom[i].type;
                        }
                        else
                        {
                            tpi = atoms_tpi->atom[i].typeB;
                        }
                        if (bBHAM)
                        {
                            /* nbfp now includes the 6.0 derivative prefactor */
                            csix    += nmolc*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                        }
                        else
                        {
                            /* nbfp now includes the 6.0/12.0 derivative prefactors */
                            csix    += nmolc*C6 (nbfp, ntp, tpi, tpj)/6.0;
                            ctwelve += nmolc*C12(nbfp, ntp, tpi, tpj)/12.0;
                        }
                        npair += nmolc;
                    }
                }
            }
        }
        if (npair - nexcl <= 0 && fplog)
        {
            fprintf(fplog, "\nWARNING: There are no atom pairs for dispersion correction\n\n");
            csix     = 0;
            ctwelve  = 0;
        }
        else
        {
            csix    /= npair - nexcl;
            ctwelve /= npair - nexcl;
        }
        if (debug)
        {
            fprintf(debug, "Counted %d exclusions\n", nexcl);
            fprintf(debug, "Average C6 parameter is: %10g\n", (double)csix);
            fprintf(debug, "Average C12 parameter is: %10g\n", (double)ctwelve);
        }
        fr->avcsix[q]    = csix;
        fr->avctwelve[q] = ctwelve;
    }

    if (EVDW_PME(fr->vdwtype))
    {
        sfree(nbfp_comb);
    }

    if (fplog != NULL)
    {
        if (fr->eDispCorr == edispcAllEner ||
            fr->eDispCorr == edispcAllEnerPres)
        {
            fprintf(fplog, "Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                    fr->avcsix[0], fr->avctwelve[0]);
        }
        else
        {
            fprintf(fplog, "Long Range LJ corr.: <C6> %10.4e\n", fr->avcsix[0]);
        }
    }
}


static void set_bham_b_max(FILE *fplog, t_forcerec *fr,
                           const gmx_mtop_t *mtop)
{
    const t_atoms *at1, *at2;
    int            mt1, mt2, i, j, tpi, tpj, ntypes;
    real           b, bmin;
    real          *nbfp;

    if (fplog)
    {
        fprintf(fplog, "Determining largest Buckingham b parameter for table\n");
    }
    nbfp   = fr->nbfp;
    ntypes = fr->ntype;

    bmin           = -1;
    fr->bham_b_max = 0;
    for (mt1 = 0; mt1 < mtop->nmoltype; mt1++)
    {
        at1 = &mtop->moltype[mt1].atoms;
        for (i = 0; (i < at1->nr); i++)
        {
            tpi = at1->atom[i].type;
            if (tpi >= ntypes)
            {
                gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", i, tpi, ntypes);
            }

            for (mt2 = mt1; mt2 < mtop->nmoltype; mt2++)
            {
                at2 = &mtop->moltype[mt2].atoms;
                for (j = 0; (j < at2->nr); j++)
                {
                    tpj = at2->atom[j].type;
                    if (tpj >= ntypes)
                    {
                        gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", j, tpj, ntypes);
                    }
                    b = BHAMB(nbfp, ntypes, tpi, tpj);
                    if (b > fr->bham_b_max)
                    {
                        fr->bham_b_max = b;
                    }
                    if ((b < bmin) || (bmin == -1))
                    {
                        bmin = b;
                    }
                }
            }
        }
    }
    if (fplog)
    {
        fprintf(fplog, "Buckingham b parameters, min: %g, max: %g\n",
                bmin, fr->bham_b_max);
    }
}

static void make_nbf_tables(FILE *fp, const output_env_t oenv,
                            t_forcerec *fr, real rtab,
                            const t_commrec *cr,
                            const char *tabfn, char *eg1, char *eg2,
                            t_nblists *nbl)
{
    char buf[STRLEN];
    int  i, j;

    if (tabfn == NULL)
    {
        if (debug)
        {
            fprintf(debug, "No table file name passed, can not read table, can not do non-bonded interactions\n");
        }
        return;
    }

    sprintf(buf, "%s", tabfn);
    if (eg1 && eg2)
    {
        /* Append the two energy group names */
        sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1, "_%s_%s.%s",
                eg1, eg2, ftp2ext(efXVG));
    }
    nbl->table_elec_vdw = make_tables(fp, oenv, fr, MASTER(cr), buf, rtab, 0);
    /* Copy the contents of the table to separate coulomb and LJ tables too,
     * to improve cache performance.
     */
    /* For performance reasons we want
     * the table data to be aligned to 16-byte. The pointers could be freed
     * but currently aren't.
     */
    nbl->table_elec.interaction   = GMX_TABLE_INTERACTION_ELEC;
    nbl->table_elec.format        = nbl->table_elec_vdw.format;
    nbl->table_elec.r             = nbl->table_elec_vdw.r;
    nbl->table_elec.n             = nbl->table_elec_vdw.n;
    nbl->table_elec.scale         = nbl->table_elec_vdw.scale;
    nbl->table_elec.scale_exp     = nbl->table_elec_vdw.scale_exp;
    nbl->table_elec.formatsize    = nbl->table_elec_vdw.formatsize;
    nbl->table_elec.ninteractions = 1;
    nbl->table_elec.stride        = nbl->table_elec.formatsize * nbl->table_elec.ninteractions;
    snew_aligned(nbl->table_elec.data, nbl->table_elec.stride*(nbl->table_elec.n+1), 32);

    nbl->table_vdw.interaction   = GMX_TABLE_INTERACTION_VDWREP_VDWDISP;
    nbl->table_vdw.format        = nbl->table_elec_vdw.format;
    nbl->table_vdw.r             = nbl->table_elec_vdw.r;
    nbl->table_vdw.n             = nbl->table_elec_vdw.n;
    nbl->table_vdw.scale         = nbl->table_elec_vdw.scale;
    nbl->table_vdw.scale_exp     = nbl->table_elec_vdw.scale_exp;
    nbl->table_vdw.formatsize    = nbl->table_elec_vdw.formatsize;
    nbl->table_vdw.ninteractions = 2;
    nbl->table_vdw.stride        = nbl->table_vdw.formatsize * nbl->table_vdw.ninteractions;
    snew_aligned(nbl->table_vdw.data, nbl->table_vdw.stride*(nbl->table_vdw.n+1), 32);

    for (i = 0; i <= nbl->table_elec_vdw.n; i++)
    {
        for (j = 0; j < 4; j++)
        {
            nbl->table_elec.data[4*i+j] = nbl->table_elec_vdw.data[12*i+j];
        }
        for (j = 0; j < 8; j++)
        {
            nbl->table_vdw.data[8*i+j] = nbl->table_elec_vdw.data[12*i+4+j];
        }
    }
}

static void count_tables(int ftype1, int ftype2, const gmx_mtop_t *mtop,
                         int *ncount, int **count)
{
    const gmx_moltype_t *molt;
    const t_ilist       *il;
    int                  mt, ftype, stride, i, j, tabnr;

    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        molt = &mtop->moltype[mt];
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (ftype == ftype1 || ftype == ftype2)
            {
                il     = &molt->ilist[ftype];
                stride = 1 + NRAL(ftype);
                for (i = 0; i < il->nr; i += stride)
                {
                    tabnr = mtop->ffparams.iparams[il->iatoms[i]].tab.table;
                    if (tabnr < 0)
                    {
                        gmx_fatal(FARGS, "A bonded table number is smaller than 0: %d\n", tabnr);
                    }
                    if (tabnr >= *ncount)
                    {
                        srenew(*count, tabnr+1);
                        for (j = *ncount; j < tabnr+1; j++)
                        {
                            (*count)[j] = 0;
                        }
                        *ncount = tabnr+1;
                    }
                    (*count)[tabnr]++;
                }
            }
        }
    }
}

static bondedtable_t *make_bonded_tables(FILE *fplog,
                                         int ftype1, int ftype2,
                                         const gmx_mtop_t *mtop,
                                         const char *basefn, const char *tabext)
{
    int            i, ncount, *count;
    char           tabfn[STRLEN];
    bondedtable_t *tab;

    tab = NULL;

    ncount = 0;
    count  = NULL;
    count_tables(ftype1, ftype2, mtop, &ncount, &count);

    if (ncount > 0)
    {
        snew(tab, ncount);
        for (i = 0; i < ncount; i++)
        {
            if (count[i] > 0)
            {
                sprintf(tabfn, "%s", basefn);
                sprintf(tabfn + strlen(basefn) - strlen(ftp2ext(efXVG)) - 1, "_%s%d.%s",
                        tabext, i, ftp2ext(efXVG));
                tab[i] = make_bonded_table(fplog, tabfn, NRAL(ftype1)-2);
            }
        }
        sfree(count);
    }

    return tab;
}

void forcerec_set_ranges(t_forcerec *fr,
                         int ncg_home, int ncg_force,
                         int natoms_force,
                         int natoms_force_constr, int natoms_f_novirsum)
{
    fr->cg0 = 0;
    fr->hcg = ncg_home;

    /* fr->ncg_force is unused in the standard code,
     * but it can be useful for modified code dealing with charge groups.
     */
    fr->ncg_force           = ncg_force;
    fr->natoms_force        = natoms_force;
    fr->natoms_force_constr = natoms_force_constr;

    if (fr->natoms_force_constr > fr->nalloc_force)
    {
        fr->nalloc_force = over_alloc_dd(fr->natoms_force_constr);

        if (fr->bTwinRange)
        {
            srenew(fr->f_twin, fr->nalloc_force);
        }
    }

    if (fr->bF_NoVirSum)
    {
        fr->f_novirsum_n = natoms_f_novirsum;
        if (fr->f_novirsum_n > fr->f_novirsum_nalloc)
        {
            fr->f_novirsum_nalloc = over_alloc_dd(fr->f_novirsum_n);
            srenew(fr->f_novirsum_alloc, fr->f_novirsum_nalloc);
        }
    }
    else
    {
        fr->f_novirsum_n = 0;
    }
}

static real cutoff_inf(real cutoff)
{
    if (cutoff == 0)
    {
        cutoff = GMX_CUTOFF_INF;
    }

    return cutoff;
}

static void make_adress_tf_tables(FILE *fp, const output_env_t oenv,
                                  t_forcerec *fr, const t_inputrec *ir,
                                  const char *tabfn, const gmx_mtop_t *mtop,
                                  matrix     box)
{
    char buf[STRLEN];
    int  i, j;

    if (tabfn == NULL)
    {
        gmx_fatal(FARGS, "No thermoforce table file given. Use -tabletf to specify a file\n");
        return;
    }

    snew(fr->atf_tabs, ir->adress->n_tf_grps);

    sprintf(buf, "%s", tabfn);
    for (i = 0; i < ir->adress->n_tf_grps; i++)
    {
        j = ir->adress->tf_table_index[i]; /* get energy group index */
        sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1, "tf_%s.%s",
                *(mtop->groups.grpname[mtop->groups.grps[egcENER].nm_ind[j]]), ftp2ext(efXVG));
        if (fp)
        {
            fprintf(fp, "loading tf table for energygrp index %d from %s\n", ir->adress->tf_table_index[i], buf);
        }
        fr->atf_tabs[i] = make_atf_table(fp, oenv, fr, buf, box);
    }

}

gmx_bool can_use_allvsall(const t_inputrec *ir, gmx_bool bPrintNote, t_commrec *cr, FILE *fp)
{
    gmx_bool bAllvsAll;

    bAllvsAll =
        (
            ir->rlist == 0            &&
            ir->rcoulomb == 0         &&
            ir->rvdw == 0             &&
            ir->ePBC == epbcNONE      &&
            ir->vdwtype == evdwCUT    &&
            ir->coulombtype == eelCUT &&
            ir->efep == efepNO        &&
            (ir->implicit_solvent == eisNO ||
             (ir->implicit_solvent == eisGBSA && (ir->gb_algorithm == egbSTILL ||
                                                  ir->gb_algorithm == egbHCT   ||
                                                  ir->gb_algorithm == egbOBC))) &&
            getenv("GMX_NO_ALLVSALL") == NULL
        );

    if (bAllvsAll && ir->opts.ngener > 1)
    {
        const char *note = "NOTE: Can not use all-vs-all force loops, because there are multiple energy monitor groups; you might get significantly higher performance when using only a single energy monitor group.\n";

        if (bPrintNote)
        {
            if (MASTER(cr))
            {
                fprintf(stderr, "\n%s\n", note);
            }
            if (fp != NULL)
            {
                fprintf(fp, "\n%s\n", note);
            }
        }
        bAllvsAll = FALSE;
    }

    if (bAllvsAll && fp && MASTER(cr))
    {
        fprintf(fp, "\nUsing SIMD all-vs-all kernels.\n\n");
    }

    return bAllvsAll;
}


static void init_forcerec_f_threads(t_forcerec *fr, int nenergrp)
{
    int t, i;

    /* These thread local data structures are used for bondeds only */
    fr->nthreads = gmx_omp_nthreads_get(emntBonded);

    if (fr->nthreads > 1)
    {
        snew(fr->f_t, fr->nthreads);
        /* Thread 0 uses the global force and energy arrays */
        for (t = 1; t < fr->nthreads; t++)
        {
            fr->f_t[t].f        = NULL;
            fr->f_t[t].f_nalloc = 0;
            snew(fr->f_t[t].fshift, SHIFTS);
            fr->f_t[t].grpp.nener = nenergrp*nenergrp;
            for (i = 0; i < egNR; i++)
            {
                snew(fr->f_t[t].grpp.ener[i], fr->f_t[t].grpp.nener);
            }
        }
    }
}


gmx_bool nbnxn_acceleration_supported(FILE             *fplog,
                                      const t_commrec  *cr,
                                      const t_inputrec *ir,
                                      gmx_bool          bGPU)
{
    if (!bGPU && (ir->vdwtype == evdwPME && ir->ljpme_combination_rule == eljpmeLB))
    {
        md_print_warn(cr, fplog, "LJ-PME with Lorentz-Berthelot is not supported with %s, falling back to %s\n",
                      bGPU ? "GPUs" : "SIMD kernels",
                      bGPU ? "CPU only" : "plain-C kernels");
        return FALSE;
    }

    return TRUE;
}


static void pick_nbnxn_kernel_cpu(const t_inputrec gmx_unused *ir,
                                  int                         *kernel_type,
                                  int                         *ewald_excl)
{
    *kernel_type = nbnxnk4x4_PlainC;
    *ewald_excl  = ewaldexclTable;

#ifdef GMX_NBNXN_SIMD
    {
#ifdef GMX_NBNXN_SIMD_4XN
        *kernel_type = nbnxnk4xN_SIMD_4xN;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
        *kernel_type = nbnxnk4xN_SIMD_2xNN;
#endif

#if defined GMX_NBNXN_SIMD_2XNN && defined GMX_NBNXN_SIMD_4XN
        /* We need to choose if we want 2x(N+N) or 4xN kernels.
         * Currently this is based on the SIMD acceleration choice,
         * but it might be better to decide this at runtime based on CPU.
         *
         * 4xN calculates more (zero) interactions, but has less pair-search
         * work and much better kernel instruction scheduling.
         *
         * Up till now we have only seen that on Intel Sandy/Ivy Bridge,
         * which doesn't have FMA, both the analytical and tabulated Ewald
         * kernels have similar pair rates for 4x8 and 2x(4+4), so we choose
         * 2x(4+4) because it results in significantly fewer pairs.
         * For RF, the raw pair rate of the 4x8 kernel is higher than 2x(4+4),
         * 10% with HT, 50% without HT. As we currently don't detect the actual
         * use of HT, use 4x8 to avoid a potential performance hit.
         * On Intel Haswell 4x8 is always faster.
         */
        *kernel_type = nbnxnk4xN_SIMD_4xN;

#ifndef GMX_SIMD_HAVE_FMA
        if (EEL_PME_EWALD(ir->coulombtype) ||
            EVDW_PME(ir->vdwtype))
        {
            /* We have Ewald kernels without FMA (Intel Sandy/Ivy Bridge).
             * There are enough instructions to make 2x(4+4) efficient.
             */
            *kernel_type = nbnxnk4xN_SIMD_2xNN;
        }
#endif
#endif  /* GMX_NBNXN_SIMD_2XNN && GMX_NBNXN_SIMD_4XN */


        if (getenv("GMX_NBNXN_SIMD_4XN") != NULL)
        {
#ifdef GMX_NBNXN_SIMD_4XN
            *kernel_type = nbnxnk4xN_SIMD_4xN;
#else
            gmx_fatal(FARGS, "SIMD 4xN kernels requested, but Gromacs has been compiled without support for these kernels");
#endif
        }
        if (getenv("GMX_NBNXN_SIMD_2XNN") != NULL)
        {
#ifdef GMX_NBNXN_SIMD_2XNN
            *kernel_type = nbnxnk4xN_SIMD_2xNN;
#else
            gmx_fatal(FARGS, "SIMD 2x(N+N) kernels requested, but Gromacs has been compiled without support for these kernels");
#endif
        }

        /* Analytical Ewald exclusion correction is only an option in
         * the SIMD kernel.
         * Since table lookup's don't parallelize with SIMD, analytical
         * will probably always be faster for a SIMD width of 8 or more.
         * With FMA analytical is sometimes faster for a width if 4 as well.
         * On BlueGene/Q, this is faster regardless of precision.
         * In single precision, this is faster on Bulldozer.
         */
#if GMX_SIMD_REAL_WIDTH >= 8 || \
        (GMX_SIMD_REAL_WIDTH >= 4 && defined GMX_SIMD_HAVE_FMA && !defined GMX_DOUBLE) || \
        defined GMX_SIMD_IBM_QPX
        *ewald_excl = ewaldexclAnalytical;
#endif
        if (getenv("GMX_NBNXN_EWALD_TABLE") != NULL)
        {
            *ewald_excl = ewaldexclTable;
        }
        if (getenv("GMX_NBNXN_EWALD_ANALYTICAL") != NULL)
        {
            *ewald_excl = ewaldexclAnalytical;
        }

    }
#endif /* GMX_NBNXN_SIMD */
}


const char *lookup_nbnxn_kernel_name(int kernel_type)
{
    const char *returnvalue = NULL;
    switch (kernel_type)
    {
        case nbnxnkNotSet:
            returnvalue = "not set";
            break;
        case nbnxnk4x4_PlainC:
            returnvalue = "plain C";
            break;
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
#ifdef GMX_NBNXN_SIMD
#if defined GMX_SIMD_X86_SSE2
            returnvalue = "SSE2";
#elif defined GMX_SIMD_X86_SSE4_1
            returnvalue = "SSE4.1";
#elif defined GMX_SIMD_X86_AVX_128_FMA
            returnvalue = "AVX_128_FMA";
#elif defined GMX_SIMD_X86_AVX_256
            returnvalue = "AVX_256";
#elif defined GMX_SIMD_X86_AVX2_256
            returnvalue = "AVX2_256";
#else
            returnvalue = "SIMD";
#endif
#else  /* GMX_NBNXN_SIMD */
            returnvalue = "not available";
#endif /* GMX_NBNXN_SIMD */
            break;
        case nbnxnk8x8x8_CUDA: returnvalue   = "CUDA"; break;
        case nbnxnk8x8x8_PlainC: returnvalue = "plain C"; break;

        case nbnxnkNR:
        default:
            gmx_fatal(FARGS, "Illegal kernel type selected");
            returnvalue = NULL;
            break;
    }
    return returnvalue;
};

static void pick_nbnxn_kernel(FILE                *fp,
                              const t_commrec     *cr,
                              gmx_bool             use_simd_kernels,
                              gmx_bool             bUseGPU,
                              gmx_bool             bEmulateGPU,
                              const t_inputrec    *ir,
                              int                 *kernel_type,
                              int                 *ewald_excl,
                              gmx_bool             bDoNonbonded)
{
    assert(kernel_type);

    *kernel_type = nbnxnkNotSet;
    *ewald_excl  = ewaldexclTable;

    if (bEmulateGPU)
    {
        *kernel_type = nbnxnk8x8x8_PlainC;

        if (bDoNonbonded)
        {
            md_print_warn(cr, fp, "Emulating a GPU run on the CPU (slow)");
        }
    }
    else if (bUseGPU)
    {
        *kernel_type = nbnxnk8x8x8_CUDA;
    }

    if (*kernel_type == nbnxnkNotSet)
    {
        /* LJ PME with LB combination rule does 7 mesh operations.
         * This so slow that we don't compile SIMD non-bonded kernels for that.
         */
        if (use_simd_kernels &&
            nbnxn_acceleration_supported(fp, cr, ir, FALSE))
        {
            pick_nbnxn_kernel_cpu(ir, kernel_type, ewald_excl);
        }
        else
        {
            *kernel_type = nbnxnk4x4_PlainC;
        }
    }

    if (bDoNonbonded && fp != NULL)
    {
        fprintf(fp, "\nUsing %s %dx%d non-bonded kernels\n\n",
                lookup_nbnxn_kernel_name(*kernel_type),
                nbnxn_kernel_pairlist_simple(*kernel_type) ? NBNXN_CPU_CLUSTER_I_SIZE : NBNXN_GPU_CLUSTER_SIZE,
                nbnxn_kernel_to_cj_size(*kernel_type));

        if (nbnxnk4x4_PlainC == *kernel_type ||
            nbnxnk8x8x8_PlainC == *kernel_type)
        {
            md_print_warn(cr, fp,
                          "WARNING: Using the slow %s kernels. This should\n"
                          "not happen during routine usage on supported platforms.\n\n",
                          lookup_nbnxn_kernel_name(*kernel_type));
        }
    }
}

static void pick_nbnxn_resources(const t_commrec     *cr,
                                 const gmx_hw_info_t *hwinfo,
                                 gmx_bool             bDoNonbonded,
                                 gmx_bool            *bUseGPU,
                                 gmx_bool            *bEmulateGPU,
                                 const gmx_gpu_opt_t *gpu_opt)
{
    gmx_bool bEmulateGPUEnvVarSet;
    char     gpu_err_str[STRLEN];

    *bUseGPU = FALSE;

    bEmulateGPUEnvVarSet = (getenv("GMX_EMULATE_GPU") != NULL);

    /* Run GPU emulation mode if GMX_EMULATE_GPU is defined. Because
     * GPUs (currently) only handle non-bonded calculations, we will
     * automatically switch to emulation if non-bonded calculations are
     * turned off via GMX_NO_NONBONDED - this is the simple and elegant
     * way to turn off GPU initialization, data movement, and cleanup.
     *
     * GPU emulation can be useful to assess the performance one can expect by
     * adding GPU(s) to the machine. The conditional below allows this even
     * if mdrun is compiled without GPU acceleration support.
     * Note that you should freezing the system as otherwise it will explode.
     */
    *bEmulateGPU = (bEmulateGPUEnvVarSet ||
                    (!bDoNonbonded &&
                     gpu_opt->ncuda_dev_use > 0));

    /* Enable GPU mode when GPUs are available or no GPU emulation is requested.
     */
    if (gpu_opt->ncuda_dev_use > 0 && !(*bEmulateGPU))
    {
        /* Each PP node will use the intra-node id-th device from the
         * list of detected/selected GPUs. */
        if (!init_gpu(cr->rank_pp_intranode, gpu_err_str,
                      &hwinfo->gpu_info, gpu_opt))
        {
            /* At this point the init should never fail as we made sure that
             * we have all the GPUs we need. If it still does, we'll bail. */
            gmx_fatal(FARGS, "On rank %d failed to initialize GPU #%d: %s",
                      cr->nodeid,
                      get_gpu_device_id(&hwinfo->gpu_info, gpu_opt,
                                        cr->rank_pp_intranode),
                      gpu_err_str);
        }

        /* Here we actually turn on hardware GPU acceleration */
        *bUseGPU = TRUE;
    }
}

gmx_bool uses_simple_tables(int                 cutoff_scheme,
                            nonbonded_verlet_t *nbv,
                            int                 group)
{
    gmx_bool bUsesSimpleTables = TRUE;
    int      grp_index;

    switch (cutoff_scheme)
    {
        case ecutsGROUP:
            bUsesSimpleTables = TRUE;
            break;
        case ecutsVERLET:
            assert(NULL != nbv && NULL != nbv->grp);
            grp_index         = (group < 0) ? 0 : (nbv->ngrp - 1);
            bUsesSimpleTables = nbnxn_kernel_pairlist_simple(nbv->grp[grp_index].kernel_type);
            break;
        default:
            gmx_incons("unimplemented");
    }
    return bUsesSimpleTables;
}

static void init_ewald_f_table(interaction_const_t *ic,
                               gmx_bool             bUsesSimpleTables,
                               real                 rtab)
{
    real maxr;

    if (bUsesSimpleTables)
    {
        /* Get the Ewald table spacing based on Coulomb and/or LJ
         * Ewald coefficients and rtol.
         */
        ic->tabq_scale = ewald_spline3_table_scale(ic);

        maxr           = (rtab > ic->rcoulomb) ? rtab : ic->rcoulomb;
        ic->tabq_size  = (int)(maxr*ic->tabq_scale) + 2;
    }
    else
    {
        ic->tabq_size = GPU_EWALD_COULOMB_FORCE_TABLE_SIZE;
        /* Subtract 2 iso 1 to avoid access out of range due to rounding */
        ic->tabq_scale = (ic->tabq_size - 2)/ic->rcoulomb;
    }

    sfree_aligned(ic->tabq_coul_FDV0);
    sfree_aligned(ic->tabq_coul_F);
    sfree_aligned(ic->tabq_coul_V);

    sfree_aligned(ic->tabq_vdw_FDV0);
    sfree_aligned(ic->tabq_vdw_F);
    sfree_aligned(ic->tabq_vdw_V);

    if (ic->eeltype == eelEWALD || EEL_PME(ic->eeltype))
    {
        /* Create the original table data in FDV0 */
        snew_aligned(ic->tabq_coul_FDV0, ic->tabq_size*4, 32);
        snew_aligned(ic->tabq_coul_F, ic->tabq_size, 32);
        snew_aligned(ic->tabq_coul_V, ic->tabq_size, 32);
        table_spline3_fill_ewald_lr(ic->tabq_coul_F, ic->tabq_coul_V, ic->tabq_coul_FDV0,
                                    ic->tabq_size, 1/ic->tabq_scale, ic->ewaldcoeff_q, v_q_ewald_lr);
    }

    if (EVDW_PME(ic->vdwtype))
    {
        snew_aligned(ic->tabq_vdw_FDV0, ic->tabq_size*4, 32);
        snew_aligned(ic->tabq_vdw_F, ic->tabq_size, 32);
        snew_aligned(ic->tabq_vdw_V, ic->tabq_size, 32);
        table_spline3_fill_ewald_lr(ic->tabq_vdw_F, ic->tabq_vdw_V, ic->tabq_vdw_FDV0,
                                    ic->tabq_size, 1/ic->tabq_scale, ic->ewaldcoeff_lj, v_lj_ewald_lr);
    }
}

void init_interaction_const_tables(FILE                *fp,
                                   interaction_const_t *ic,
                                   gmx_bool             bUsesSimpleTables,
                                   real                 rtab)
{
    real spacing;

    if (ic->eeltype == eelEWALD || EEL_PME(ic->eeltype) || EVDW_PME(ic->vdwtype))
    {
        init_ewald_f_table(ic, bUsesSimpleTables, rtab);

        if (fp != NULL)
        {
            fprintf(fp, "Initialized non-bonded Ewald correction tables, spacing: %.2e size: %d\n\n",
                    1/ic->tabq_scale, ic->tabq_size);
        }
    }
}

static void clear_force_switch_constants(shift_consts_t *sc)
{
    sc->c2   = 0;
    sc->c3   = 0;
    sc->cpot = 0;
}

static void force_switch_constants(real p,
                                   real rsw, real rc,
                                   shift_consts_t *sc)
{
    /* Here we determine the coefficient for shifting the force to zero
     * between distance rsw and the cut-off rc.
     * For a potential of r^-p, we have force p*r^-(p+1).
     * But to save flops we absorb p in the coefficient.
     * Thus we get:
     * force/p   = r^-(p+1) + c2*r^2 + c3*r^3
     * potential = r^-p + c2/3*r^3 + c3/4*r^4 + cpot
     */
    sc->c2   =  ((p + 1)*rsw - (p + 4)*rc)/(pow(rc, p + 2)*pow(rc - rsw, 2));
    sc->c3   = -((p + 1)*rsw - (p + 3)*rc)/(pow(rc, p + 2)*pow(rc - rsw, 3));
    sc->cpot = -pow(rc, -p) + p*sc->c2/3*pow(rc - rsw, 3) + p*sc->c3/4*pow(rc - rsw, 4);
}

static void potential_switch_constants(real rsw, real rc,
                                       switch_consts_t *sc)
{
    /* The switch function is 1 at rsw and 0 at rc.
     * The derivative and second derivate are zero at both ends.
     * rsw        = max(r - r_switch, 0)
     * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
     * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
     * force      = force*dsw - potential*sw
     * potential *= sw
     */
    sc->c3 = -10*pow(rc - rsw, -3);
    sc->c4 =  15*pow(rc - rsw, -4);
    sc->c5 =  -6*pow(rc - rsw, -5);
}

static void
init_interaction_const(FILE                       *fp,
                       const t_commrec gmx_unused *cr,
                       interaction_const_t       **interaction_const,
                       const t_forcerec           *fr,
                       real                        rtab)
{
    interaction_const_t *ic;
    gmx_bool             bUsesSimpleTables = TRUE;

    snew(ic, 1);

    /* Just allocate something so we can free it */
    snew_aligned(ic->tabq_coul_FDV0, 16, 32);
    snew_aligned(ic->tabq_coul_F, 16, 32);
    snew_aligned(ic->tabq_coul_V, 16, 32);

    ic->rlist           = fr->rlist;
    ic->rlistlong       = fr->rlistlong;

    /* Lennard-Jones */
    ic->vdwtype         = fr->vdwtype;
    ic->vdw_modifier    = fr->vdw_modifier;
    ic->rvdw            = fr->rvdw;
    ic->rvdw_switch     = fr->rvdw_switch;
    ic->ewaldcoeff_lj   = fr->ewaldcoeff_lj;
    ic->ljpme_comb_rule = fr->ljpme_combination_rule;
    ic->sh_lj_ewald     = 0;
    clear_force_switch_constants(&ic->dispersion_shift);
    clear_force_switch_constants(&ic->repulsion_shift);

    switch (ic->vdw_modifier)
    {
        case eintmodPOTSHIFT:
            /* Only shift the potential, don't touch the force */
            ic->dispersion_shift.cpot = -pow(ic->rvdw, -6.0);
            ic->repulsion_shift.cpot  = -pow(ic->rvdw, -12.0);
            if (EVDW_PME(ic->vdwtype))
            {
                real crc2;

                crc2            = sqr(ic->ewaldcoeff_lj*ic->rvdw);
                ic->sh_lj_ewald = (exp(-crc2)*(1 + crc2 + 0.5*crc2*crc2) - 1)*pow(ic->rvdw, -6.0);
            }
            break;
        case eintmodFORCESWITCH:
            /* Switch the force, switch and shift the potential */
            force_switch_constants(6.0, ic->rvdw_switch, ic->rvdw,
                                   &ic->dispersion_shift);
            force_switch_constants(12.0, ic->rvdw_switch, ic->rvdw,
                                   &ic->repulsion_shift);
            break;
        case eintmodPOTSWITCH:
            /* Switch the potential and force */
            potential_switch_constants(ic->rvdw_switch, ic->rvdw,
                                       &ic->vdw_switch);
            break;
        case eintmodNONE:
        case eintmodEXACTCUTOFF:
            /* Nothing to do here */
            break;
        default:
            gmx_incons("unimplemented potential modifier");
    }

    ic->sh_invrc6 = -ic->dispersion_shift.cpot;

    /* Electrostatics */
    ic->eeltype          = fr->eeltype;
    ic->coulomb_modifier = fr->coulomb_modifier;
    ic->rcoulomb         = fr->rcoulomb;
    ic->epsilon_r        = fr->epsilon_r;
    ic->epsfac           = fr->epsfac;
    ic->ewaldcoeff_q     = fr->ewaldcoeff_q;

    if (fr->coulomb_modifier == eintmodPOTSHIFT)
    {
        ic->sh_ewald = gmx_erfc(ic->ewaldcoeff_q*ic->rcoulomb);
    }
    else
    {
        ic->sh_ewald = 0;
    }

    /* Reaction-field */
    if (EEL_RF(ic->eeltype))
    {
        ic->epsilon_rf = fr->epsilon_rf;
        ic->k_rf       = fr->k_rf;
        ic->c_rf       = fr->c_rf;
    }
    else
    {
        /* For plain cut-off we might use the reaction-field kernels */
        ic->epsilon_rf = ic->epsilon_r;
        ic->k_rf       = 0;
        if (fr->coulomb_modifier == eintmodPOTSHIFT)
        {
            ic->c_rf   = 1/ic->rcoulomb;
        }
        else
        {
            ic->c_rf   = 0;
        }
    }

    if (fp != NULL)
    {
        real dispersion_shift;

        dispersion_shift = ic->dispersion_shift.cpot;
        if (EVDW_PME(ic->vdwtype))
        {
            dispersion_shift -= ic->sh_lj_ewald;
        }
        fprintf(fp, "Potential shift: LJ r^-12: %.3e r^-6: %.3e",
                ic->repulsion_shift.cpot, dispersion_shift);

        if (ic->eeltype == eelCUT)
        {
            fprintf(fp, ", Coulomb %.e", -ic->c_rf);
        }
        else if (EEL_PME(ic->eeltype))
        {
            fprintf(fp, ", Ewald %.3e", -ic->sh_ewald);
        }
        fprintf(fp, "\n");
    }

    *interaction_const = ic;

    if (fr->nbv != NULL && fr->nbv->bUseGPU)
    {
        nbnxn_cuda_init_const(fr->nbv->cu_nbv, ic, fr->nbv->grp);

        /* With tMPI + GPUs some ranks may be sharing GPU(s) and therefore
         * also sharing texture references. To keep the code simple, we don't
         * treat texture references as shared resources, but this means that
         * the coulomb_tab and nbfp texture refs will get updated by multiple threads.
         * Hence, to ensure that the non-bonded kernels don't start before all
         * texture binding operations are finished, we need to wait for all ranks
         * to arrive here before continuing.
         *
         * Note that we could omit this barrier if GPUs are not shared (or
         * texture objects are used), but as this is initialization code, there
         * is not point in complicating things.
         */
#ifdef GMX_THREAD_MPI
        if (PAR(cr))
        {
            gmx_barrier(cr);
        }
#endif  /* GMX_THREAD_MPI */
    }

    bUsesSimpleTables = uses_simple_tables(fr->cutoff_scheme, fr->nbv, -1);
    init_interaction_const_tables(fp, ic, bUsesSimpleTables, rtab);
}

static void init_nb_verlet(FILE                *fp,
                           nonbonded_verlet_t **nb_verlet,
                           gmx_bool             bFEP_NonBonded,
                           const t_inputrec    *ir,
                           const t_forcerec    *fr,
                           const t_commrec     *cr,
                           const char          *nbpu_opt)
{
    nonbonded_verlet_t *nbv;
    int                 i;
    char               *env;
    gmx_bool            bEmulateGPU, bHybridGPURun = FALSE;

    nbnxn_alloc_t      *nb_alloc;
    nbnxn_free_t       *nb_free;

    snew(nbv, 1);

    pick_nbnxn_resources(cr, fr->hwinfo,
                         fr->bNonbonded,
                         &nbv->bUseGPU,
                         &bEmulateGPU,
                         fr->gpu_opt);

    nbv->nbs = NULL;

    nbv->ngrp = (DOMAINDECOMP(cr) ? 2 : 1);
    for (i = 0; i < nbv->ngrp; i++)
    {
        nbv->grp[i].nbl_lists.nnbl = 0;
        nbv->grp[i].nbat           = NULL;
        nbv->grp[i].kernel_type    = nbnxnkNotSet;

        if (i == 0) /* local */
        {
            pick_nbnxn_kernel(fp, cr, fr->use_simd_kernels,
                              nbv->bUseGPU, bEmulateGPU, ir,
                              &nbv->grp[i].kernel_type,
                              &nbv->grp[i].ewald_excl,
                              fr->bNonbonded);
        }
        else /* non-local */
        {
            if (nbpu_opt != NULL && strcmp(nbpu_opt, "gpu_cpu") == 0)
            {
                /* Use GPU for local, select a CPU kernel for non-local */
                pick_nbnxn_kernel(fp, cr, fr->use_simd_kernels,
                                  FALSE, FALSE, ir,
                                  &nbv->grp[i].kernel_type,
                                  &nbv->grp[i].ewald_excl,
                                  fr->bNonbonded);

                bHybridGPURun = TRUE;
            }
            else
            {
                /* Use the same kernel for local and non-local interactions */
                nbv->grp[i].kernel_type = nbv->grp[0].kernel_type;
                nbv->grp[i].ewald_excl  = nbv->grp[0].ewald_excl;
            }
        }
    }

    if (nbv->bUseGPU)
    {
        /* init the NxN GPU data; the last argument tells whether we'll have
         * both local and non-local NB calculation on GPU */
        nbnxn_cuda_init(fp, &nbv->cu_nbv,
                        &fr->hwinfo->gpu_info, fr->gpu_opt,
                        cr->rank_pp_intranode,
                        (nbv->ngrp > 1) && !bHybridGPURun);

        if ((env = getenv("GMX_NB_MIN_CI")) != NULL)
        {
            char *end;

            nbv->min_ci_balanced = strtol(env, &end, 10);
            if (!end || (*end != 0) || nbv->min_ci_balanced <= 0)
            {
                gmx_fatal(FARGS, "Invalid value passed in GMX_NB_MIN_CI=%s, positive integer required", env);
            }

            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (passed as env. var.)\n",
                        nbv->min_ci_balanced);
            }
        }
        else
        {
            nbv->min_ci_balanced = nbnxn_cuda_min_ci_balanced(nbv->cu_nbv);
            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (auto-adjusted to the number of GPU multi-processors)\n",
                        nbv->min_ci_balanced);
            }
        }
    }
    else
    {
        nbv->min_ci_balanced = 0;
    }

    *nb_verlet = nbv;

    nbnxn_init_search(&nbv->nbs,
                      DOMAINDECOMP(cr) ? &cr->dd->nc : NULL,
                      DOMAINDECOMP(cr) ? domdec_zones(cr->dd) : NULL,
                      bFEP_NonBonded,
                      gmx_omp_nthreads_get(emntPairsearch));

    for (i = 0; i < nbv->ngrp; i++)
    {
        if (nbv->grp[0].kernel_type == nbnxnk8x8x8_CUDA)
        {
            nb_alloc = &pmalloc;
            nb_free  = &pfree;
        }
        else
        {
            nb_alloc = NULL;
            nb_free  = NULL;
        }

        nbnxn_init_pairlist_set(&nbv->grp[i].nbl_lists,
                                nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type),
                                /* 8x8x8 "non-simple" lists are ATM always combined */
                                !nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type),
                                nb_alloc, nb_free);

        if (i == 0 ||
            nbv->grp[0].kernel_type != nbv->grp[i].kernel_type)
        {
            gmx_bool bSimpleList;
            int      enbnxninitcombrule;

            bSimpleList = nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type);

            if (bSimpleList && (fr->vdwtype == evdwCUT && (fr->vdw_modifier == eintmodNONE || fr->vdw_modifier == eintmodPOTSHIFT)))
            {
                /* Plain LJ cut-off: we can optimize with combination rules */
                enbnxninitcombrule = enbnxninitcombruleDETECT;
            }
            else if (fr->vdwtype == evdwPME)
            {
                /* LJ-PME: we need to use a combination rule for the grid */
                if (fr->ljpme_combination_rule == eljpmeGEOM)
                {
                    enbnxninitcombrule = enbnxninitcombruleGEOM;
                }
                else
                {
                    enbnxninitcombrule = enbnxninitcombruleLB;
                }
            }
            else
            {
                /* We use a full combination matrix: no rule required */
                enbnxninitcombrule = enbnxninitcombruleNONE;
            }


            snew(nbv->grp[i].nbat, 1);
            nbnxn_atomdata_init(fp,
                                nbv->grp[i].nbat,
                                nbv->grp[i].kernel_type,
                                enbnxninitcombrule,
                                fr->ntype, fr->nbfp,
                                ir->opts.ngener,
                                bSimpleList ? gmx_omp_nthreads_get(emntNonbonded) : 1,
                                nb_alloc, nb_free);
        }
        else
        {
            nbv->grp[i].nbat = nbv->grp[0].nbat;
        }
    }
}

void init_forcerec(FILE              *fp,
                   const output_env_t oenv,
                   t_forcerec        *fr,
                   t_fcdata          *fcd,
                   const t_inputrec  *ir,
                   const gmx_mtop_t  *mtop,
                   const t_commrec   *cr,
                   matrix             box,
                   const char        *tabfn,
                   const char        *tabafn,
                   const char        *tabpfn,
                   const char        *tabbfn,
                   const char        *nbpu_opt,
                   gmx_bool           bNoSolvOpt,
                   real               print_force)
{
    int            i, j, m, natoms, ngrp, negp_pp, negptable, egi, egj;
    real           rtab;
    char          *env;
    double         dbl;
    const t_block *cgs;
    gmx_bool       bGenericKernelOnly;
    gmx_bool       bMakeTables, bMakeSeparate14Table, bSomeNormalNbListsAreInUse;
    gmx_bool       bFEP_NonBonded;
    t_nblists     *nbl;
    int           *nm_ind, egp_flags;

    if (fr->hwinfo == NULL)
    {
        /* Detect hardware, gather information.
         * In mdrun, hwinfo has already been set before calling init_forcerec.
         * Here we ignore GPUs, as tools will not use them anyhow.
         */
        fr->hwinfo = gmx_detect_hardware(fp, cr, FALSE);
    }

    /* By default we turn SIMD kernels on, but it might be turned off further down... */
    fr->use_simd_kernels = TRUE;

    fr->bDomDec = DOMAINDECOMP(cr);

    natoms = mtop->natoms;

    if (check_box(ir->ePBC, box))
    {
        gmx_fatal(FARGS, check_box(ir->ePBC, box));
    }

    /* Test particle insertion ? */
    if (EI_TPI(ir->eI))
    {
        /* Set to the size of the molecule to be inserted (the last one) */
        /* Because of old style topologies, we have to use the last cg
         * instead of the last molecule type.
         */
        cgs       = &mtop->moltype[mtop->molblock[mtop->nmolblock-1].type].cgs;
        fr->n_tpi = cgs->index[cgs->nr] - cgs->index[cgs->nr-1];
        if (fr->n_tpi != mtop->mols.index[mtop->mols.nr] - mtop->mols.index[mtop->mols.nr-1])
        {
            gmx_fatal(FARGS, "The molecule to insert can not consist of multiple charge groups.\nMake it a single charge group.");
        }
    }
    else
    {
        fr->n_tpi = 0;
    }

    /* Copy AdResS parameters */
    if (ir->bAdress)
    {
        fr->adress_type           = ir->adress->type;
        fr->adress_const_wf       = ir->adress->const_wf;
        fr->adress_ex_width       = ir->adress->ex_width;
        fr->adress_hy_width       = ir->adress->hy_width;
        fr->adress_icor           = ir->adress->icor;
        fr->adress_site           = ir->adress->site;
        fr->adress_ex_forcecap    = ir->adress->ex_forcecap;
        fr->adress_do_hybridpairs = ir->adress->do_hybridpairs;


        snew(fr->adress_group_explicit, ir->adress->n_energy_grps);
        for (i = 0; i < ir->adress->n_energy_grps; i++)
        {
            fr->adress_group_explicit[i] = ir->adress->group_explicit[i];
        }

        fr->n_adress_tf_grps = ir->adress->n_tf_grps;
        snew(fr->adress_tf_table_index, fr->n_adress_tf_grps);
        for (i = 0; i < fr->n_adress_tf_grps; i++)
        {
            fr->adress_tf_table_index[i] = ir->adress->tf_table_index[i];
        }
        copy_rvec(ir->adress->refs, fr->adress_refs);
    }
    else
    {
        fr->adress_type           = eAdressOff;
        fr->adress_do_hybridpairs = FALSE;
    }

    /* Copy the user determined parameters */
    fr->userint1  = ir->userint1;
    fr->userint2  = ir->userint2;
    fr->userint3  = ir->userint3;
    fr->userint4  = ir->userint4;
    fr->userreal1 = ir->userreal1;
    fr->userreal2 = ir->userreal2;
    fr->userreal3 = ir->userreal3;
    fr->userreal4 = ir->userreal4;

    /* Shell stuff */
    fr->fc_stepsize = ir->fc_stepsize;

    /* Free energy */
    fr->efep        = ir->efep;
    fr->sc_alphavdw = ir->fepvals->sc_alpha;
    if (ir->fepvals->bScCoul)
    {
        fr->sc_alphacoul  = ir->fepvals->sc_alpha;
        fr->sc_sigma6_min = pow(ir->fepvals->sc_sigma_min, 6);
    }
    else
    {
        fr->sc_alphacoul  = 0;
        fr->sc_sigma6_min = 0; /* only needed when bScCoul is on */
    }
    fr->sc_power      = ir->fepvals->sc_power;
    fr->sc_r_power    = ir->fepvals->sc_r_power;
    fr->sc_sigma6_def = pow(ir->fepvals->sc_sigma, 6);

    env = getenv("GMX_SCSIGMA_MIN");
    if (env != NULL)
    {
        dbl = 0;
        sscanf(env, "%lf", &dbl);
        fr->sc_sigma6_min = pow(dbl, 6);
        if (fp)
        {
            fprintf(fp, "Setting the minimum soft core sigma to %g nm\n", dbl);
        }
    }

    fr->bNonbonded = TRUE;
    if (getenv("GMX_NO_NONBONDED") != NULL)
    {
        /* turn off non-bonded calculations */
        fr->bNonbonded = FALSE;
        md_print_warn(cr, fp,
                      "Found environment variable GMX_NO_NONBONDED.\n"
                      "Disabling nonbonded calculations.\n");
    }

    bGenericKernelOnly = FALSE;

    /* We now check in the NS code whether a particular combination of interactions
     * can be used with water optimization, and disable it if that is not the case.
     */

    if (getenv("GMX_NB_GENERIC") != NULL)
    {
        if (fp != NULL)
        {
            fprintf(fp,
                    "Found environment variable GMX_NB_GENERIC.\n"
                    "Disabling all interaction-specific nonbonded kernels, will only\n"
                    "use the slow generic ones in src/gmxlib/nonbonded/nb_generic.c\n\n");
        }
        bGenericKernelOnly = TRUE;
    }

    if (bGenericKernelOnly == TRUE)
    {
        bNoSolvOpt         = TRUE;
    }

    if ( (getenv("GMX_DISABLE_SIMD_KERNELS") != NULL) || (getenv("GMX_NOOPTIMIZEDKERNELS") != NULL) )
    {
        fr->use_simd_kernels = FALSE;
        if (fp != NULL)
        {
            fprintf(fp,
                    "\nFound environment variable GMX_DISABLE_SIMD_KERNELS.\n"
                    "Disabling the usage of any SIMD-specific kernel routines (e.g. SSE2/SSE4.1/AVX).\n\n");
        }
    }

    fr->bBHAM = (mtop->ffparams.functype[0] == F_BHAM);

    /* Check if we can/should do all-vs-all kernels */
    fr->bAllvsAll       = can_use_allvsall(ir, FALSE, NULL, NULL);
    fr->AllvsAll_work   = NULL;
    fr->AllvsAll_workgb = NULL;

    /* All-vs-all kernels have not been implemented in 4.6, and
     * the SIMD group kernels are also buggy in this case. Non-SIMD
     * group kernels are OK. See Redmine #1249. */
    if (fr->bAllvsAll)
    {
        fr->bAllvsAll            = FALSE;
        fr->use_simd_kernels     = FALSE;
        if (fp != NULL)
        {
            fprintf(fp,
                    "\nYour simulation settings would have triggered the efficient all-vs-all\n"
                    "kernels in GROMACS 4.5, but these have not been implemented in GROMACS\n"
                    "4.6. Also, we can't use the accelerated SIMD kernels here because\n"
                    "of an unfixed bug. The reference C kernels are correct, though, so\n"
                    "we are proceeding by disabling all CPU architecture-specific\n"
                    "(e.g. SSE2/SSE4/AVX) routines. If performance is important, please\n"
                    "use GROMACS 4.5.7 or try cutoff-scheme = Verlet.\n\n");
        }
    }

    /* Neighbour searching stuff */
    fr->cutoff_scheme = ir->cutoff_scheme;
    fr->bGrid         = (ir->ns_type == ensGRID);
    fr->ePBC          = ir->ePBC;

    if (fr->cutoff_scheme == ecutsGROUP)
    {
        const char *note = "NOTE: This file uses the deprecated 'group' cutoff_scheme. This will be\n"
            "removed in a future release when 'verlet' supports all interaction forms.\n";

        if (MASTER(cr))
        {
            fprintf(stderr, "\n%s\n", note);
        }
        if (fp != NULL)
        {
            fprintf(fp, "\n%s\n", note);
        }
    }

    /* Determine if we will do PBC for distances in bonded interactions */
    if (fr->ePBC == epbcNONE)
    {
        fr->bMolPBC = FALSE;
    }
    else
    {
        if (!DOMAINDECOMP(cr))
        {
            /* The group cut-off scheme and SHAKE assume charge groups
             * are whole, but not using molpbc is faster in most cases.
             */
            if (fr->cutoff_scheme == ecutsGROUP ||
                (ir->eConstrAlg == econtSHAKE &&
                 (gmx_mtop_ftype_count(mtop, F_CONSTR) > 0 ||
                  gmx_mtop_ftype_count(mtop, F_CONSTRNC) > 0)))
            {
                fr->bMolPBC = ir->bPeriodicMols;
            }
            else
            {
                fr->bMolPBC = TRUE;
                if (getenv("GMX_USE_GRAPH") != NULL)
                {
                    fr->bMolPBC = FALSE;
                    if (fp)
                    {
                        fprintf(fp, "\nGMX_MOLPBC is set, using the graph for bonded interactions\n\n");
                    }
                }
            }
        }
        else
        {
            fr->bMolPBC = dd_bonded_molpbc(cr->dd, fr->ePBC);
        }
    }
    fr->bGB = (ir->implicit_solvent == eisGBSA);

    fr->rc_scaling = ir->refcoord_scaling;
    copy_rvec(ir->posres_com, fr->posres_com);
    copy_rvec(ir->posres_comB, fr->posres_comB);
    fr->rlist                    = cutoff_inf(ir->rlist);
    fr->rlistlong                = cutoff_inf(ir->rlistlong);
    fr->eeltype                  = ir->coulombtype;
    fr->vdwtype                  = ir->vdwtype;
    fr->ljpme_combination_rule   = ir->ljpme_combination_rule;

    fr->coulomb_modifier = ir->coulomb_modifier;
    fr->vdw_modifier     = ir->vdw_modifier;

    /* Electrostatics: Translate from interaction-setting-in-mdp-file to kernel interaction format */
    switch (fr->eeltype)
    {
        case eelCUT:
            fr->nbkernel_elec_interaction = (fr->bGB) ? GMX_NBKERNEL_ELEC_GENERALIZEDBORN : GMX_NBKERNEL_ELEC_COULOMB;
            break;

        case eelRF:
        case eelGRF:
        case eelRF_NEC:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_REACTIONFIELD;
            break;

        case eelRF_ZERO:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_REACTIONFIELD;
            fr->coulomb_modifier          = eintmodEXACTCUTOFF;
            break;

        case eelSWITCH:
        case eelSHIFT:
        case eelUSER:
        case eelENCADSHIFT:
        case eelPMESWITCH:
        case eelPMEUSER:
        case eelPMEUSERSWITCH:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_CUBICSPLINETABLE;
            break;

        case eelPME:
        case eelEWALD:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_EWALD;
            break;

        default:
            gmx_fatal(FARGS, "Unsupported electrostatic interaction: %s", eel_names[fr->eeltype]);
            break;
    }

    /* Vdw: Translate from mdp settings to kernel format */
    switch (fr->vdwtype)
    {
        case evdwCUT:
            if (fr->bBHAM)
            {
                fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_BUCKINGHAM;
            }
            else
            {
                fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_LENNARDJONES;
            }
            break;
        case evdwPME:
            fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_LJEWALD;
            break;

        case evdwSWITCH:
        case evdwSHIFT:
        case evdwUSER:
        case evdwENCADSHIFT:
            fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_CUBICSPLINETABLE;
            break;

        default:
            gmx_fatal(FARGS, "Unsupported vdw interaction: %s", evdw_names[fr->vdwtype]);
            break;
    }

    /* These start out identical to ir, but might be altered if we e.g. tabulate the interaction in the kernel */
    fr->nbkernel_elec_modifier    = fr->coulomb_modifier;
    fr->nbkernel_vdw_modifier     = fr->vdw_modifier;

    fr->rvdw             = cutoff_inf(ir->rvdw);
    fr->rvdw_switch      = ir->rvdw_switch;
    fr->rcoulomb         = cutoff_inf(ir->rcoulomb);
    fr->rcoulomb_switch  = ir->rcoulomb_switch;

    fr->bTwinRange = fr->rlistlong > fr->rlist;
    fr->bEwald     = (EEL_PME(fr->eeltype) || fr->eeltype == eelEWALD);

    fr->reppow     = mtop->ffparams.reppow;

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        fr->bvdwtab    = ((fr->vdwtype != evdwCUT || !gmx_within_tol(fr->reppow, 12.0, 10*GMX_DOUBLE_EPS))
                          && !EVDW_PME(fr->vdwtype));
        /* We have special kernels for standard Ewald and PME, but the pme-switch ones are tabulated above */
        fr->bcoultab   = !(fr->eeltype == eelCUT ||
                           fr->eeltype == eelEWALD ||
                           fr->eeltype == eelPME ||
                           fr->eeltype == eelRF ||
                           fr->eeltype == eelRF_ZERO);

        /* If the user absolutely wants different switch/shift settings for coul/vdw, it is likely
         * going to be faster to tabulate the interaction than calling the generic kernel.
         * However, if generic kernels have been requested we keep things analytically.
         */
        if (fr->nbkernel_elec_modifier == eintmodPOTSWITCH &&
            fr->nbkernel_vdw_modifier == eintmodPOTSWITCH &&
            bGenericKernelOnly == FALSE)
        {
            if ((fr->rcoulomb_switch != fr->rvdw_switch) || (fr->rcoulomb != fr->rvdw))
            {
                fr->bcoultab = TRUE;
                /* Once we tabulate electrostatics, we can use the switch function for LJ,
                 * which would otherwise need two tables.
                 */
            }
        }
        else if ((fr->nbkernel_elec_modifier == eintmodPOTSHIFT && fr->nbkernel_vdw_modifier == eintmodPOTSHIFT) ||
                 ((fr->nbkernel_elec_interaction == GMX_NBKERNEL_ELEC_REACTIONFIELD &&
                   fr->nbkernel_elec_modifier == eintmodEXACTCUTOFF &&
                   (fr->nbkernel_vdw_modifier == eintmodPOTSWITCH || fr->nbkernel_vdw_modifier == eintmodPOTSHIFT))))
        {
            if ((fr->rcoulomb != fr->rvdw) && (bGenericKernelOnly == FALSE))
            {
                fr->bcoultab = TRUE;
            }
        }

        if (fr->nbkernel_elec_modifier == eintmodFORCESWITCH)
        {
            fr->bcoultab = TRUE;
        }
        if (fr->nbkernel_vdw_modifier == eintmodFORCESWITCH)
        {
            fr->bvdwtab = TRUE;
        }

        if (getenv("GMX_REQUIRE_TABLES"))
        {
            fr->bvdwtab  = TRUE;
            fr->bcoultab = TRUE;
        }

        if (fp)
        {
            fprintf(fp, "Table routines are used for coulomb: %s\n", bool_names[fr->bcoultab]);
            fprintf(fp, "Table routines are used for vdw:     %s\n", bool_names[fr->bvdwtab ]);
        }

        if (fr->bvdwtab == TRUE)
        {
            fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_CUBICSPLINETABLE;
            fr->nbkernel_vdw_modifier    = eintmodNONE;
        }
        if (fr->bcoultab == TRUE)
        {
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_CUBICSPLINETABLE;
            fr->nbkernel_elec_modifier    = eintmodNONE;
        }
    }

    if (ir->cutoff_scheme == ecutsVERLET)
    {
        if (!gmx_within_tol(fr->reppow, 12.0, 10*GMX_DOUBLE_EPS))
        {
            gmx_fatal(FARGS, "Cut-off scheme %S only supports LJ repulsion power 12", ecutscheme_names[ir->cutoff_scheme]);
        }
        fr->bvdwtab  = FALSE;
        fr->bcoultab = FALSE;
    }

    /* Tables are used for direct ewald sum */
    if (fr->bEwald)
    {
        if (EEL_PME(ir->coulombtype))
        {
            if (fp)
            {
                fprintf(fp, "Will do PME sum in reciprocal space for electrostatic interactions.\n");
            }
            if (ir->coulombtype == eelP3M_AD)
            {
                please_cite(fp, "Hockney1988");
                please_cite(fp, "Ballenegger2012");
            }
            else
            {
                please_cite(fp, "Essmann95a");
            }

            if (ir->ewald_geometry == eewg3DC)
            {
                if (fp)
                {
                    fprintf(fp, "Using the Ewald3DC correction for systems with a slab geometry.\n");
                }
                please_cite(fp, "In-Chul99a");
            }
        }
        fr->ewaldcoeff_q = calc_ewaldcoeff_q(ir->rcoulomb, ir->ewald_rtol);
        init_ewald_tab(&(fr->ewald_table), ir, fp);
        if (fp)
        {
            fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for Ewald\n",
                    1/fr->ewaldcoeff_q);
        }
    }

    if (EVDW_PME(ir->vdwtype))
    {
        if (fp)
        {
            fprintf(fp, "Will do PME sum in reciprocal space for LJ dispersion interactions.\n");
        }
        please_cite(fp, "Essmann95a");
        fr->ewaldcoeff_lj = calc_ewaldcoeff_lj(ir->rvdw, ir->ewald_rtol_lj);
        if (fp)
        {
            fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for LJ Ewald\n",
                    1/fr->ewaldcoeff_lj);
        }
    }

    /* Electrostatics */
    fr->epsilon_r       = ir->epsilon_r;
    fr->epsilon_rf      = ir->epsilon_rf;
    fr->fudgeQQ         = mtop->ffparams.fudgeQQ;

    /* Parameters for generalized RF */
    fr->zsquare = 0.0;
    fr->temp    = 0.0;

    if (fr->eeltype == eelGRF)
    {
        init_generalized_rf(fp, mtop, ir, fr);
    }

    fr->bF_NoVirSum = (EEL_FULL(fr->eeltype) || EVDW_PME(fr->vdwtype) ||
                       gmx_mtop_ftype_count(mtop, F_POSRES) > 0 ||
                       gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0 ||
                       IR_ELEC_FIELD(*ir) ||
                       (fr->adress_icor != eAdressICOff)
                       );

    if (fr->cutoff_scheme == ecutsGROUP &&
        ncg_mtop(mtop) > fr->cg_nalloc && !DOMAINDECOMP(cr))
    {
        /* Count the total number of charge groups */
        fr->cg_nalloc = ncg_mtop(mtop);
        srenew(fr->cg_cm, fr->cg_nalloc);
    }
    if (fr->shift_vec == NULL)
    {
        snew(fr->shift_vec, SHIFTS);
    }

    if (fr->fshift == NULL)
    {
        snew(fr->fshift, SHIFTS);
    }

    if (fr->nbfp == NULL)
    {
        fr->ntype = mtop->ffparams.atnr;
        fr->nbfp  = mk_nbfp(&mtop->ffparams, fr->bBHAM);
        if (EVDW_PME(fr->vdwtype))
        {
            fr->ljpme_c6grid  = make_ljpme_c6grid(&mtop->ffparams, fr);
        }
    }

    /* Copy the energy group exclusions */
    fr->egp_flags = ir->opts.egp_flags;

    /* Van der Waals stuff */
    if ((fr->vdwtype != evdwCUT) && (fr->vdwtype != evdwUSER) && !fr->bBHAM)
    {
        if (fr->rvdw_switch >= fr->rvdw)
        {
            gmx_fatal(FARGS, "rvdw_switch (%f) must be < rvdw (%f)",
                      fr->rvdw_switch, fr->rvdw);
        }
        if (fp)
        {
            fprintf(fp, "Using %s Lennard-Jones, switch between %g and %g nm\n",
                    (fr->eeltype == eelSWITCH) ? "switched" : "shifted",
                    fr->rvdw_switch, fr->rvdw);
        }
    }

    if (fr->bBHAM && EVDW_PME(fr->vdwtype))
    {
        gmx_fatal(FARGS, "LJ PME not supported with Buckingham");
    }

    if (fr->bBHAM && (fr->vdwtype == evdwSHIFT || fr->vdwtype == evdwSWITCH))
    {
        gmx_fatal(FARGS, "Switch/shift interaction not supported with Buckingham");
    }

    if (fr->bBHAM && fr->cutoff_scheme == ecutsVERLET)
    {
        gmx_fatal(FARGS, "Verlet cutoff-scheme is not supported with Buckingham");
    }

    if (fp)
    {
        fprintf(fp, "Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
                fr->rlist, fr->rcoulomb, fr->bBHAM ? "BHAM" : "LJ", fr->rvdw);
    }

    fr->eDispCorr = ir->eDispCorr;
    if (ir->eDispCorr != edispcNO)
    {
        set_avcsixtwelve(fp, fr, mtop);
    }

    if (fr->bBHAM)
    {
        set_bham_b_max(fp, fr, mtop);
    }

    fr->gb_epsilon_solvent = ir->gb_epsilon_solvent;

    /* Copy the GBSA data (radius, volume and surftens for each
     * atomtype) from the topology atomtype section to forcerec.
     */
    snew(fr->atype_radius, fr->ntype);
    snew(fr->atype_vol, fr->ntype);
    snew(fr->atype_surftens, fr->ntype);
    snew(fr->atype_gb_radius, fr->ntype);
    snew(fr->atype_S_hct, fr->ntype);

    if (mtop->atomtypes.nr > 0)
    {
        for (i = 0; i < fr->ntype; i++)
        {
            fr->atype_radius[i] = mtop->atomtypes.radius[i];
        }
        for (i = 0; i < fr->ntype; i++)
        {
            fr->atype_vol[i] = mtop->atomtypes.vol[i];
        }
        for (i = 0; i < fr->ntype; i++)
        {
            fr->atype_surftens[i] = mtop->atomtypes.surftens[i];
        }
        for (i = 0; i < fr->ntype; i++)
        {
            fr->atype_gb_radius[i] = mtop->atomtypes.gb_radius[i];
        }
        for (i = 0; i < fr->ntype; i++)
        {
            fr->atype_S_hct[i] = mtop->atomtypes.S_hct[i];
        }
    }

    /* Generate the GB table if needed */
    if (fr->bGB)
    {
#ifdef GMX_DOUBLE
        fr->gbtabscale = 2000;
#else
        fr->gbtabscale = 500;
#endif

        fr->gbtabr = 100;
        fr->gbtab  = make_gb_table(oenv, fr);

        init_gb(&fr->born, fr, ir, mtop, ir->gb_algorithm);

        /* Copy local gb data (for dd, this is done in dd_partition_system) */
        if (!DOMAINDECOMP(cr))
        {
            make_local_gb(cr, fr->born, ir->gb_algorithm);
        }
    }

    /* Set the charge scaling */
    if (fr->epsilon_r != 0)
    {
        fr->epsfac = ONE_4PI_EPS0/fr->epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no coulomb interactions */
        fr->epsfac = 0;
    }

    /* Reaction field constants */
    if (EEL_RF(fr->eeltype))
    {
        calc_rffac(fp, fr->eeltype, fr->epsilon_r, fr->epsilon_rf,
                   fr->rcoulomb, fr->temp, fr->zsquare, box,
                   &fr->kappa, &fr->k_rf, &fr->c_rf);
    }

    /*This now calculates sum for q and c6*/
    set_chargesum(fp, fr, mtop);

    /* if we are using LR electrostatics, and they are tabulated,
     * the tables will contain modified coulomb interactions.
     * Since we want to use the non-shifted ones for 1-4
     * coulombic interactions, we must have an extra set of tables.
     */

    /* Construct tables.
     * A little unnecessary to make both vdw and coul tables sometimes,
     * but what the heck... */

    bMakeTables = fr->bcoultab || fr->bvdwtab || fr->bEwald ||
        (ir->eDispCorr != edispcNO && ir_vdw_switched(ir));

    bMakeSeparate14Table = ((!bMakeTables || fr->eeltype != eelCUT || fr->vdwtype != evdwCUT ||
                             fr->coulomb_modifier != eintmodNONE ||
                             fr->vdw_modifier != eintmodNONE ||
                             fr->bBHAM || fr->bEwald) &&
                            (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 ||
                             gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0 ||
                             gmx_mtop_ftype_count(mtop, F_LJC_PAIRS_NB) > 0));

    negp_pp   = ir->opts.ngener - ir->nwall;
    negptable = 0;
    if (!bMakeTables)
    {
        bSomeNormalNbListsAreInUse = TRUE;
        fr->nnblists               = 1;
    }
    else
    {
        bSomeNormalNbListsAreInUse = (ir->eDispCorr != edispcNO);
        for (egi = 0; egi < negp_pp; egi++)
        {
            for (egj = egi; egj < negp_pp; egj++)
            {
                egp_flags = ir->opts.egp_flags[GID(egi, egj, ir->opts.ngener)];
                if (!(egp_flags & EGP_EXCL))
                {
                    if (egp_flags & EGP_TABLE)
                    {
                        negptable++;
                    }
                    else
                    {
                        bSomeNormalNbListsAreInUse = TRUE;
                    }
                }
            }
        }
        if (bSomeNormalNbListsAreInUse)
        {
            fr->nnblists = negptable + 1;
        }
        else
        {
            fr->nnblists = negptable;
        }
        if (fr->nnblists > 1)
        {
            snew(fr->gid2nblists, ir->opts.ngener*ir->opts.ngener);
        }
    }

    if (ir->adress)
    {
        fr->nnblists *= 2;
    }

    snew(fr->nblists, fr->nnblists);

    /* This code automatically gives table length tabext without cut-off's,
     * in that case grompp should already have checked that we do not need
     * normal tables and we only generate tables for 1-4 interactions.
     */
    rtab = ir->rlistlong + ir->tabext;

    if (bMakeTables)
    {
        /* make tables for ordinary interactions */
        if (bSomeNormalNbListsAreInUse)
        {
            make_nbf_tables(fp, oenv, fr, rtab, cr, tabfn, NULL, NULL, &fr->nblists[0]);
            if (ir->adress)
            {
                make_nbf_tables(fp, oenv, fr, rtab, cr, tabfn, NULL, NULL, &fr->nblists[fr->nnblists/2]);
            }
            if (!bMakeSeparate14Table)
            {
                fr->tab14 = fr->nblists[0].table_elec_vdw;
            }
            m = 1;
        }
        else
        {
            m = 0;
        }
        if (negptable > 0)
        {
            /* Read the special tables for certain energy group pairs */
            nm_ind = mtop->groups.grps[egcENER].nm_ind;
            for (egi = 0; egi < negp_pp; egi++)
            {
                for (egj = egi; egj < negp_pp; egj++)
                {
                    egp_flags = ir->opts.egp_flags[GID(egi, egj, ir->opts.ngener)];
                    if ((egp_flags & EGP_TABLE) && !(egp_flags & EGP_EXCL))
                    {
                        nbl = &(fr->nblists[m]);
                        if (fr->nnblists > 1)
                        {
                            fr->gid2nblists[GID(egi, egj, ir->opts.ngener)] = m;
                        }
                        /* Read the table file with the two energy groups names appended */
                        make_nbf_tables(fp, oenv, fr, rtab, cr, tabfn,
                                        *mtop->groups.grpname[nm_ind[egi]],
                                        *mtop->groups.grpname[nm_ind[egj]],
                                        &fr->nblists[m]);
                        if (ir->adress)
                        {
                            make_nbf_tables(fp, oenv, fr, rtab, cr, tabfn,
                                            *mtop->groups.grpname[nm_ind[egi]],
                                            *mtop->groups.grpname[nm_ind[egj]],
                                            &fr->nblists[fr->nnblists/2+m]);
                        }
                        m++;
                    }
                    else if (fr->nnblists > 1)
                    {
                        fr->gid2nblists[GID(egi, egj, ir->opts.ngener)] = 0;
                    }
                }
            }
        }
    }
    else if ((fr->eDispCorr != edispcNO) &&
             ((fr->vdw_modifier == eintmodPOTSWITCH) ||
              (fr->vdw_modifier == eintmodFORCESWITCH) ||
              (fr->vdw_modifier == eintmodPOTSHIFT)))
    {
        /* Tables might not be used for the potential modifier interactions per se, but
         * we still need them to evaluate switch/shift dispersion corrections in this case.
         */
        make_nbf_tables(fp, oenv, fr, rtab, cr, tabfn, NULL, NULL, &fr->nblists[0]);
    }

    if (bMakeSeparate14Table)
    {
        /* generate extra tables with plain Coulomb for 1-4 interactions only */
        fr->tab14 = make_tables(fp, oenv, fr, MASTER(cr), tabpfn, rtab,
                                GMX_MAKETABLES_14ONLY);
    }

    /* Read AdResS Thermo Force table if needed */
    if (fr->adress_icor == eAdressICThermoForce)
    {
        /* old todo replace */

        if (ir->adress->n_tf_grps > 0)
        {
            make_adress_tf_tables(fp, oenv, fr, ir, tabfn, mtop, box);

        }
        else
        {
            /* load the default table */
            snew(fr->atf_tabs, 1);
            fr->atf_tabs[DEFAULT_TF_TABLE] = make_atf_table(fp, oenv, fr, tabafn, box);
        }
    }

    /* Wall stuff */
    fr->nwall = ir->nwall;
    if (ir->nwall && ir->wall_type == ewtTABLE)
    {
        make_wall_tables(fp, oenv, ir, tabfn, &mtop->groups, fr);
    }

    if (fcd && tabbfn)
    {
        fcd->bondtab  = make_bonded_tables(fp,
                                           F_TABBONDS, F_TABBONDSNC,
                                           mtop, tabbfn, "b");
        fcd->angletab = make_bonded_tables(fp,
                                           F_TABANGLES, -1,
                                           mtop, tabbfn, "a");
        fcd->dihtab   = make_bonded_tables(fp,
                                           F_TABDIHS, -1,
                                           mtop, tabbfn, "d");
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "No fcdata or table file name passed, can not read table, can not do bonded interactions\n");
        }
    }

    /* QM/MM initialization if requested
     */
    if (ir->bQMMM)
    {
        fprintf(stderr, "QM/MM calculation requested.\n");
    }

    fr->bQMMM      = ir->bQMMM;
    fr->qr         = mk_QMMMrec();

    /* Set all the static charge group info */
    fr->cginfo_mb = init_cginfo_mb(fp, mtop, fr, bNoSolvOpt,
                                   &bFEP_NonBonded,
                                   &fr->bExcl_IntraCGAll_InterCGNone);
    if (DOMAINDECOMP(cr))
    {
        fr->cginfo = NULL;
    }
    else
    {
        fr->cginfo = cginfo_expand(mtop->nmolblock, fr->cginfo_mb);
    }

    if (!DOMAINDECOMP(cr))
    {
        forcerec_set_ranges(fr, ncg_mtop(mtop), ncg_mtop(mtop),
                            mtop->natoms, mtop->natoms, mtop->natoms);
    }

    fr->print_force = print_force;


    /* coarse load balancing vars */
    fr->t_fnbf    = 0.;
    fr->t_wait    = 0.;
    fr->timesteps = 0;

    /* Initialize neighbor search */
    init_ns(fp, cr, &fr->ns, fr, mtop);

    if (cr->duty & DUTY_PP)
    {
        gmx_nonbonded_setup(fr, bGenericKernelOnly);
        /*
           if (ir->bAdress)
            {
                gmx_setup_adress_kernels(fp,bGenericKernelOnly);
            }
         */
    }

    /* Initialize the thread working data for bonded interactions */
    init_forcerec_f_threads(fr, mtop->groups.grps[egcENER].nr);

    snew(fr->excl_load, fr->nthreads+1);

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        if (ir->rcoulomb != ir->rvdw)
        {
            gmx_fatal(FARGS, "With Verlet lists rcoulomb and rvdw should be identical");
        }

        init_nb_verlet(fp, &fr->nbv, bFEP_NonBonded, ir, fr, cr, nbpu_opt);
    }

    /* fr->ic is used both by verlet and group kernels (to some extent) now */
    init_interaction_const(fp, cr, &fr->ic, fr, rtab);

    if (ir->eDispCorr != edispcNO)
    {
        calc_enervirdiff(fp, ir->eDispCorr, fr);
    }
}

#define pr_real(fp, r) fprintf(fp, "%s: %e\n",#r, r)
#define pr_int(fp, i)  fprintf((fp), "%s: %d\n",#i, i)
#define pr_bool(fp, b) fprintf((fp), "%s: %s\n",#b, bool_names[b])

void pr_forcerec(FILE *fp, t_forcerec *fr)
{
    int i;

    pr_real(fp, fr->rlist);
    pr_real(fp, fr->rcoulomb);
    pr_real(fp, fr->fudgeQQ);
    pr_bool(fp, fr->bGrid);
    pr_bool(fp, fr->bTwinRange);
    /*pr_int(fp,fr->cg0);
       pr_int(fp,fr->hcg);*/
    for (i = 0; i < fr->nnblists; i++)
    {
        pr_int(fp, fr->nblists[i].table_elec_vdw.n);
    }
    pr_real(fp, fr->rcoulomb_switch);
    pr_real(fp, fr->rcoulomb);

    fflush(fp);
}

void forcerec_set_excl_load(t_forcerec           *fr,
                            const gmx_localtop_t *top)
{
    const int *ind, *a;
    int        t, i, j, ntot, n, ntarget;

    ind = top->excls.index;
    a   = top->excls.a;

    ntot = 0;
    for (i = 0; i < top->excls.nr; i++)
    {
        for (j = ind[i]; j < ind[i+1]; j++)
        {
            if (a[j] > i)
            {
                ntot++;
            }
        }
    }

    fr->excl_load[0] = 0;
    n                = 0;
    i                = 0;
    for (t = 1; t <= fr->nthreads; t++)
    {
        ntarget = (ntot*t)/fr->nthreads;
        while (i < top->excls.nr && n < ntarget)
        {
            for (j = ind[i]; j < ind[i+1]; j++)
            {
                if (a[j] > i)
                {
                    n++;
                }
            }
            i++;
        }
        fr->excl_load[t] = i;
    }
}
