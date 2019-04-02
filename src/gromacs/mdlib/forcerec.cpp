/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "forcerec.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/ewald.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/listed_forces/manage_threading.h"
#include "gromacs/listed_forces/pairs.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec_threading.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdlib/wall.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_geometry.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

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
    int        i, j, k, atnr;
    real       c6, c6i, c6j, c12i, c12j, epsi, epsj, sigmai, sigmaj;
    real      *grid;

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
            c6   = std::sqrt(c6i * c6j);
            if (fr->ljpme_combination_rule == eljpmeLB
                && !gmx_numzero(c6) && !gmx_numzero(c12i) && !gmx_numzero(c12j))
            {
                sigmai = gmx::sixthroot(c12i / c6i);
                sigmaj = gmx::sixthroot(c12j / c6j);
                epsi   = c6i * c6i / c12i;
                epsj   = c6j * c6j / c12j;
                c6     = std::sqrt(epsi * epsj) * gmx::power6(0.5*(sigmai+sigmaj));
            }
            /* Store the elements at the same relative positions as C6 in nbfp in order
             * to simplify access in the kernels
             */
            grid[2*(atnr*i+j)] = c6*6.0;
        }
    }
    return grid;
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
                 const AtomGroupIndices *qm_grps,
                 t_forcerec             *fr,
                 int                    *n_solvent_parameters,
                 solvent_parameters_t  **solvent_parameters_p,
                 int                     cginfo,
                 int                    *cg_sp)
{
    t_atom               *atom;
    int                   j, k;
    int                   j0, j1, nj;
    gmx_bool              perturbed;
    gmx_bool              has_vdw[4];
    gmx_bool              match;
    real                  tmp_charge[4]  = { 0.0 }; /* init to zero to make gcc 7 happy */
    int                   tmp_vdwtype[4] = { 0 };   /* init to zero to make gcc 7 happy */
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
    if (qm_grpnr != nullptr)
    {
        for (j = j0; j < j1 && !qm; j++)
        {
            qm = (qm_grpnr[j] < qm_grps->size() - 1);
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
        for (j = 0; j < nj && match; j++)
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
        if (match)
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
        for (k = 0; k < fr->ntype && (!has_vdw[j]); k++)
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
        if (!has_vdw[1] &&
            !has_vdw[2] &&
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
        if (!has_vdw[1] &&
            !has_vdw[2] &&
            !has_vdw[3] &&
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
    const gmx_moltype_t  *molt;
    int                   mol, cg_mol, at_offset, am, cgm, i, nmol_ch, nmol;
    int                   n_solvent_parameters;
    solvent_parameters_t *solvent_parameters;
    int                 **cg_sp;
    int                   bestsp, bestsol;

    if (debug)
    {
        fprintf(debug, "Going to determine what solvent types we have.\n");
    }

    n_solvent_parameters = 0;
    solvent_parameters   = nullptr;
    /* Allocate temporary array for solvent type */
    snew(cg_sp, mtop->molblock.size());

    at_offset = 0;
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
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
                                 mtop->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].empty() ?
                                 nullptr : mtop->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].data()+at_offset+am,
                                 &mtop->groups.groups[SimulationAtomGroupType::QuantumMechanics],
                                 fr,
                                 &n_solvent_parameters, &solvent_parameters,
                                 cginfo_mb[mb].cginfo[cgm+cg_mol],
                                 &cg_sp[mb][cgm+cg_mol]);
            }
        }
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

    fr->nWatMol = 0;
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
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

    if (bestsol != esolNO && fp != nullptr)
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
    int                   cg_offset, a_offset;
    int                   m, cg, a0, a1, gid, ai, j, aj, excl_nalloc;
    int                  *a_con;
    int                   ftype;
    int                   ia;
    gmx_bool              bId, *bExcl, bExclIntraAll, bExclInter, bHaveVDW, bHaveQ, bHavePerturbedAtoms;

    snew(cginfo_mb, mtop->molblock.size());

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
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
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
        for (m = 0; m < molb->nmol; m++)
        {
            int am = m*cgs->index[cgs->nr];
            for (cg = 0; cg < cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                if (getGroupType(mtop->groups, SimulationAtomGroupType::QuantumMechanics, a_offset+am+a0) !=
                    getGroupType(mtop->groups, SimulationAtomGroupType::QuantumMechanics, a_offset   +a0))
                {
                    bId = FALSE;
                }
                if (!mtop->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].empty())
                {
                    for (ai = a0; ai < a1; ai++)
                    {
                        if (mtop->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics][a_offset+am+ai] !=
                            mtop->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics][a_offset   +ai])
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
                for (ia = 0; ia < molt->ilist[ftype].size(); ia += 1+nral)
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
            int cgm = m*cgs->nr;
            int am  = m*cgs->index[cgs->nr];
            for (cg = 0; cg < cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];

                /* Store the energy group in cginfo */
                gid = getGroupType(mtop->groups, SimulationAtomGroupType::EnergyOutput, a_offset+am+a0);
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
    sfree(type_VDW);
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
        for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
        {
            for (cg = 0; cg < cginfo_mb[mb].cg_mod; cg++)
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[cg], esolNO);
            }
        }
    }

    return cginfo_mb;
}

static std::vector<int> cginfo_expand(const int          nmb,
                                      const cginfo_mb_t *cgi_mb)
{
    const int        ncg = cgi_mb[nmb - 1].cg_end;

    std::vector<int> cginfo(ncg);

    int              mb = 0;
    for (int cg = 0; cg < ncg; cg++)
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

static void done_cginfo_mb(cginfo_mb_t *cginfo_mb, int numMolBlocks)
{
    if (cginfo_mb == nullptr)
    {
        return;
    }
    for (int mb = 0; mb < numMolBlocks; ++mb)
    {
        sfree(cginfo_mb[mb].cginfo);
    }
    sfree(cginfo_mb);
}

/* Sets the sum of charges (squared) and C6 in the system in fr.
 * Returns whether the system has a net charge.
 */
static bool set_chargesum(FILE *log, t_forcerec *fr, const gmx_mtop_t *mtop)
{
    /*This now calculates sum for q and c6*/
    double         qsum, q2sum, q, c6sum, c6;

    qsum   = 0;
    q2sum  = 0;
    c6sum  = 0;
    for (const gmx_molblock_t &molb : mtop->molblock)
    {
        int            nmol  = molb.nmol;
        const t_atoms *atoms = &mtop->moltype[molb.type].atoms;
        for (int i = 0; i < atoms->nr; i++)
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
        for (const gmx_molblock_t &molb : mtop->molblock)
        {
            int            nmol  = molb.nmol;
            const t_atoms *atoms = &mtop->moltype[molb.type].atoms;
            for (int i = 0; i < atoms->nr; i++)
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

    /* A cut-off of 1e-4 is used to catch rounding errors due to ascii input */
    return (std::abs(fr->qsum[0]) > 1e-4 ||
            std::abs(fr->qsum[1]) > 1e-4);
}

void update_forcerec(t_forcerec *fr, matrix box)
{
    if (fr->ic->eeltype == eelGRF)
    {
        calc_rffac(nullptr, fr->ic->eeltype, fr->ic->epsilon_r, fr->ic->epsilon_rf,
                   fr->ic->rcoulomb, fr->temp, fr->zsquare, box,
                   &fr->ic->k_rf, &fr->ic->c_rf);
    }
}

static real calcBuckinghamBMax(FILE *fplog, const gmx_mtop_t *mtop)
{
    const t_atoms *at1, *at2;
    int            i, j, tpi, tpj, ntypes;
    real           b, bmin;

    if (fplog)
    {
        fprintf(fplog, "Determining largest Buckingham b parameter for table\n");
    }
    ntypes = mtop->ffparams.atnr;

    bmin            = -1;
    real bham_b_max = 0;
    for (size_t mt1 = 0; mt1 < mtop->moltype.size(); mt1++)
    {
        at1 = &mtop->moltype[mt1].atoms;
        for (i = 0; (i < at1->nr); i++)
        {
            tpi = at1->atom[i].type;
            if (tpi >= ntypes)
            {
                gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", i, tpi, ntypes);
            }

            for (size_t mt2 = mt1; mt2 < mtop->moltype.size(); mt2++)
            {
                at2 = &mtop->moltype[mt2].atoms;
                for (j = 0; (j < at2->nr); j++)
                {
                    tpj = at2->atom[j].type;
                    if (tpj >= ntypes)
                    {
                        gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", j, tpj, ntypes);
                    }
                    b = mtop->ffparams.iparams[tpi*ntypes + tpj].bham.b;
                    if (b > bham_b_max)
                    {
                        bham_b_max = b;
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
                bmin, bham_b_max);
    }

    return bham_b_max;
}

static void make_nbf_tables(FILE *fp,
                            const interaction_const_t *ic, real rtab,
                            const char *tabfn, char *eg1, char *eg2,
                            t_nblists *nbl)
{
    char buf[STRLEN];
    int  i, j;

    if (tabfn == nullptr)
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
    nbl->table_elec_vdw = make_tables(fp, ic, buf, rtab, 0);
    /* Copy the contents of the table to separate coulomb and LJ tables too,
     * to improve cache performance.
     */
    /* For performance reasons we want
     * the table data to be aligned to 16-byte.
     */
    nbl->table_elec                =
        new t_forcetable(GMX_TABLE_INTERACTION_ELEC,
                         nbl->table_elec_vdw->format);
    nbl->table_elec->r             = nbl->table_elec_vdw->r;
    nbl->table_elec->n             = nbl->table_elec_vdw->n;
    nbl->table_elec->scale         = nbl->table_elec_vdw->scale;
    nbl->table_elec->formatsize    = nbl->table_elec_vdw->formatsize;
    nbl->table_elec->ninteractions = 1;
    nbl->table_elec->stride        = nbl->table_elec->formatsize * nbl->table_elec->ninteractions;
    snew_aligned(nbl->table_elec->data, nbl->table_elec->stride*(nbl->table_elec->n+1), 32);

    nbl->table_vdw                =
        new t_forcetable(GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
                         nbl->table_elec_vdw->format);
    nbl->table_vdw->r             = nbl->table_elec_vdw->r;
    nbl->table_vdw->n             = nbl->table_elec_vdw->n;
    nbl->table_vdw->scale         = nbl->table_elec_vdw->scale;
    nbl->table_vdw->formatsize    = nbl->table_elec_vdw->formatsize;
    nbl->table_vdw->ninteractions = 2;
    nbl->table_vdw->stride        = nbl->table_vdw->formatsize * nbl->table_vdw->ninteractions;
    snew_aligned(nbl->table_vdw->data, nbl->table_vdw->stride*(nbl->table_vdw->n+1), 32);

    /* NOTE: Using a single i-loop here leads to mix-up of data in table_vdw
     *       with (at least) gcc 6.2, 6.3 and 6.4 when compiled with -O3 and AVX
     */
    for (i = 0; i <= nbl->table_elec_vdw->n; i++)
    {
        for (j = 0; j < 4; j++)
        {
            nbl->table_elec->data[4*i+j] = nbl->table_elec_vdw->data[12*i+j];
        }
    }
    for (i = 0; i <= nbl->table_elec_vdw->n; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nbl->table_vdw->data[8*i+j] = nbl->table_elec_vdw->data[12*i+4+j];
        }
    }
}

/*!\brief If there's bonded interactions of type \c ftype1 or \c
 * ftype2 present in the topology, build an array of the number of
 * interactions present for each bonded interaction index found in the
 * topology.
 *
 * \c ftype1 or \c ftype2 may be set to -1 to disable seeking for a
 * valid type with that parameter.
 *
 * \c count will be reallocated as necessary to fit the largest bonded
 * interaction index found, and its current size will be returned in
 * \c ncount. It will contain zero for every bonded interaction index
 * for which no interactions are present in the topology.
 */
static void count_tables(int ftype1, int ftype2, const gmx_mtop_t *mtop,
                         int *ncount, int **count)
{
    int                  ftype, i, j, tabnr;

    // Loop over all moleculetypes
    for (const gmx_moltype_t &molt : mtop->moltype)
    {
        // Loop over all interaction types
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            // If the current interaction type is one of the types whose tables we're trying to count...
            if (ftype == ftype1 || ftype == ftype2)
            {
                const InteractionList &il     = molt.ilist[ftype];
                const int              stride = 1 + NRAL(ftype);
                // ... and there are actually some interactions for this type
                for (i = 0; i < il.size(); i += stride)
                {
                    // Find out which table index the user wanted
                    tabnr = mtop->ffparams.iparams[il.iatoms[i]].tab.table;
                    if (tabnr < 0)
                    {
                        gmx_fatal(FARGS, "A bonded table number is smaller than 0: %d\n", tabnr);
                    }
                    // Make room for this index in the data structure
                    if (tabnr >= *ncount)
                    {
                        srenew(*count, tabnr+1);
                        for (j = *ncount; j < tabnr+1; j++)
                        {
                            (*count)[j] = 0;
                        }
                        *ncount = tabnr+1;
                    }
                    // Record that this table index is used and must have a valid file
                    (*count)[tabnr]++;
                }
            }
        }
    }
}

/*!\brief If there's bonded interactions of flavour \c tabext and type
 * \c ftype1 or \c ftype2 present in the topology, seek them in the
 * list of filenames passed to mdrun, and make bonded tables from
 * those files.
 *
 * \c ftype1 or \c ftype2 may be set to -1 to disable seeking for a
 * valid type with that parameter.
 *
 * A fatal error occurs if no matching filename is found.
 */
static bondedtable_t *make_bonded_tables(FILE *fplog,
                                         int ftype1, int ftype2,
                                         const gmx_mtop_t *mtop,
                                         gmx::ArrayRef<const std::string> tabbfnm,
                                         const char *tabext)
{
    int            ncount, *count;
    bondedtable_t *tab;

    tab = nullptr;

    ncount = 0;
    count  = nullptr;
    count_tables(ftype1, ftype2, mtop, &ncount, &count);

    // Are there any relevant tabulated bond interactions?
    if (ncount > 0)
    {
        snew(tab, ncount);
        for (int i = 0; i < ncount; i++)
        {
            // Do any interactions exist that requires this table?
            if (count[i] > 0)
            {
                // This pattern enforces the current requirement that
                // table filenames end in a characteristic sequence
                // before the file type extension, and avoids table 13
                // being recognized and used for table 1.
                std::string patternToFind = gmx::formatString("_%s%d.%s", tabext, i, ftp2ext(efXVG));
                bool        madeTable     = false;
                for (gmx::index j = 0; j < tabbfnm.ssize() && !madeTable; ++j)
                {
                    if (gmx::endsWith(tabbfnm[j], patternToFind))
                    {
                        // Finally read the table from the file found
                        tab[i]    = make_bonded_table(fplog, tabbfnm[j].c_str(), NRAL(ftype1)-2);
                        madeTable = true;
                    }
                }
                if (!madeTable)
                {
                    bool isPlural = (ftype2 != -1);
                    gmx_fatal(FARGS, "Tabulated interaction of type '%s%s%s' with index %d cannot be used because no table file whose name matched '%s' was passed via the gmx mdrun -tableb command-line option.",
                              interaction_function[ftype1].longname,
                              isPlural ? "' or '" : "",
                              isPlural ? interaction_function[ftype2].longname : "",
                              i,
                              patternToFind.c_str());
                }
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
    }

    if (fr->haveDirectVirialContributions)
    {
        fr->forceBufferForDirectVirialContributions->resize(natoms_f_novirsum);
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

/*! \brief Print Coulomb Ewald citations and set ewald coefficients */
static void initCoulombEwaldParameters(FILE *fp, const t_inputrec *ir,
                                       bool systemHasNetCharge,
                                       interaction_const_t *ic)
{
    if (!EEL_PME_EWALD(ir->coulombtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for electrostatic interactions.\n");

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
                fprintf(fp, "Using the Ewald3DC correction for systems with a slab geometry%s.\n",
                        systemHasNetCharge ? " and net charge" : "");
            }
            please_cite(fp, "In-Chul99a");
            if (systemHasNetCharge)
            {
                please_cite(fp, "Ballenegger2009");
            }
        }
    }

    ic->ewaldcoeff_q = calc_ewaldcoeff_q(ir->rcoulomb, ir->ewald_rtol);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for Ewald\n",
                1/ic->ewaldcoeff_q);
    }

    if (ic->coulomb_modifier == eintmodPOTSHIFT)
    {
        GMX_RELEASE_ASSERT(ic->rcoulomb != 0, "Cutoff radius cannot be zero");
        ic->sh_ewald = std::erfc(ic->ewaldcoeff_q*ic->rcoulomb) / ic->rcoulomb;
    }
    else
    {
        ic->sh_ewald = 0;
    }
}

/*! \brief Print Van der Waals Ewald citations and set ewald coefficients */
static void initVdwEwaldParameters(FILE *fp, const t_inputrec *ir,
                                   interaction_const_t *ic)
{
    if (!EVDW_PME(ir->vdwtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for LJ dispersion interactions.\n");
        please_cite(fp, "Essmann95a");
    }
    ic->ewaldcoeff_lj = calc_ewaldcoeff_lj(ir->rvdw, ir->ewald_rtol_lj);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for LJ Ewald\n",
                1/ic->ewaldcoeff_lj);
    }

    if (ic->vdw_modifier == eintmodPOTSHIFT)
    {
        real crc2       = gmx::square(ic->ewaldcoeff_lj*ic->rvdw);
        ic->sh_lj_ewald = (std::exp(-crc2)*(1 + crc2 + 0.5*crc2*crc2) - 1)/gmx::power6(ic->rvdw);
    }
    else
    {
        ic->sh_lj_ewald = 0;
    }
}

gmx_bool uses_simple_tables(int                       cutoff_scheme,
                            const nonbonded_verlet_t *nbv)
{
    gmx_bool bUsesSimpleTables = TRUE;

    switch (cutoff_scheme)
    {
        case ecutsGROUP:
            bUsesSimpleTables = TRUE;
            break;
        case ecutsVERLET:
            GMX_RELEASE_ASSERT(nullptr != nbv, "A non-bonded verlet object is required with the Verlet cutoff-scheme");
            bUsesSimpleTables = nbv->pairlistIsSimple();
            break;
        default:
            gmx_incons("unimplemented");
    }
    return bUsesSimpleTables;
}

static void init_ewald_f_table(interaction_const_t *ic,
                               real                 rtab)
{
    real maxr;

    /* Get the Ewald table spacing based on Coulomb and/or LJ
     * Ewald coefficients and rtol.
     */
    ic->tabq_scale = ewald_spline3_table_scale(ic);

    if (ic->cutoff_scheme == ecutsVERLET)
    {
        maxr = ic->rcoulomb;
    }
    else
    {
        maxr = std::max(ic->rcoulomb, rtab);
    }
    ic->tabq_size  = static_cast<int>(maxr*ic->tabq_scale) + 2;

    sfree_aligned(ic->tabq_coul_FDV0);
    sfree_aligned(ic->tabq_coul_F);
    sfree_aligned(ic->tabq_coul_V);

    sfree_aligned(ic->tabq_vdw_FDV0);
    sfree_aligned(ic->tabq_vdw_F);
    sfree_aligned(ic->tabq_vdw_V);

    if (EEL_PME_EWALD(ic->eeltype))
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
                                   real                 rtab)
{
    if (EEL_PME_EWALD(ic->eeltype) || EVDW_PME(ic->vdwtype))
    {
        init_ewald_f_table(ic, rtab);

        if (fp != nullptr)
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
    sc->c2   =  ((p + 1)*rsw - (p + 4)*rc)/(pow(rc, p + 2)*gmx::square(rc - rsw));
    sc->c3   = -((p + 1)*rsw - (p + 3)*rc)/(pow(rc, p + 2)*gmx::power3(rc - rsw));
    sc->cpot = -pow(rc, -p) + p*sc->c2/3*gmx::power3(rc - rsw) + p*sc->c3/4*gmx::power4(rc - rsw);
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
    sc->c3 = -10/gmx::power3(rc - rsw);
    sc->c4 =  15/gmx::power4(rc - rsw);
    sc->c5 =  -6/gmx::power5(rc - rsw);
}

/*! \brief Construct interaction constants
 *
 * This data is used (particularly) by search and force code for
 * short-range interactions. Many of these are constant for the whole
 * simulation; some are constant only after PME tuning completes.
 */
static void
init_interaction_const(FILE                       *fp,
                       interaction_const_t       **interaction_const,
                       const t_inputrec           *ir,
                       const gmx_mtop_t           *mtop,
                       bool                        systemHasNetCharge)
{
    interaction_const_t *ic;

    snew(ic, 1);

    ic->cutoff_scheme   = ir->cutoff_scheme;

    /* Just allocate something so we can free it */
    snew_aligned(ic->tabq_coul_FDV0, 16, 32);
    snew_aligned(ic->tabq_coul_F, 16, 32);
    snew_aligned(ic->tabq_coul_V, 16, 32);

    /* Lennard-Jones */
    ic->vdwtype         = ir->vdwtype;
    ic->vdw_modifier    = ir->vdw_modifier;
    ic->reppow          = mtop->ffparams.reppow;
    ic->rvdw            = cutoff_inf(ir->rvdw);
    ic->rvdw_switch     = ir->rvdw_switch;
    ic->ljpme_comb_rule = ir->ljpme_combination_rule;
    ic->useBuckingham   = (mtop->ffparams.functype[0] == F_BHAM);
    if (ic->useBuckingham)
    {
        ic->buckinghamBMax = calcBuckinghamBMax(fp, mtop);
    }

    initVdwEwaldParameters(fp, ir, ic);

    clear_force_switch_constants(&ic->dispersion_shift);
    clear_force_switch_constants(&ic->repulsion_shift);

    switch (ic->vdw_modifier)
    {
        case eintmodPOTSHIFT:
            /* Only shift the potential, don't touch the force */
            ic->dispersion_shift.cpot = -1.0/gmx::power6(ic->rvdw);
            ic->repulsion_shift.cpot  = -1.0/gmx::power12(ic->rvdw);
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
    ic->eeltype          = ir->coulombtype;
    ic->coulomb_modifier = ir->coulomb_modifier;
    ic->rcoulomb         = cutoff_inf(ir->rcoulomb);
    ic->rcoulomb_switch  = ir->rcoulomb_switch;
    ic->epsilon_r        = ir->epsilon_r;

    /* Set the Coulomb energy conversion factor */
    if (ic->epsilon_r != 0)
    {
        ic->epsfac = ONE_4PI_EPS0/ic->epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        ic->epsfac = 0;
    }

    /* Reaction-field */
    if (EEL_RF(ic->eeltype))
    {
        ic->epsilon_rf = ir->epsilon_rf;
        /* Generalized reaction field parameters are updated every step */
        if (ic->eeltype != eelGRF)
        {
            calc_rffac(fp, ic->eeltype, ic->epsilon_r, ic->epsilon_rf,
                       ic->rcoulomb, 0, 0, nullptr,
                       &ic->k_rf, &ic->c_rf);
        }

        if (ir->cutoff_scheme == ecutsGROUP && ic->eeltype == eelRF_ZERO)
        {
            /* grompp should have done this, but this scheme is obsolete */
            ic->coulomb_modifier = eintmodEXACTCUTOFF;
        }
    }
    else
    {
        /* For plain cut-off we might use the reaction-field kernels */
        ic->epsilon_rf = ic->epsilon_r;
        ic->k_rf       = 0;
        if (ir->coulomb_modifier == eintmodPOTSHIFT)
        {
            ic->c_rf   = 1/ic->rcoulomb;
        }
        else
        {
            ic->c_rf   = 0;
        }
    }

    initCoulombEwaldParameters(fp, ir, systemHasNetCharge, ic);

    if (fp != nullptr)
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
}

static void
done_interaction_const(interaction_const_t *interaction_const)
{
    sfree_aligned(interaction_const->tabq_coul_FDV0);
    sfree_aligned(interaction_const->tabq_coul_F);
    sfree_aligned(interaction_const->tabq_coul_V);
    sfree(interaction_const);
}

void init_forcerec(FILE                             *fp,
                   const gmx::MDLogger              &mdlog,
                   t_forcerec                       *fr,
                   t_fcdata                         *fcd,
                   const t_inputrec                 *ir,
                   const gmx_mtop_t                 *mtop,
                   const t_commrec                  *cr,
                   matrix                            box,
                   const char                       *tabfn,
                   const char                       *tabpfn,
                   gmx::ArrayRef<const std::string>  tabbfnm,
                   const gmx_hw_info_t              &hardwareInfo,
                   const gmx_device_info_t          *deviceInfo,
                   const bool                        useGpuForBonded,
                   gmx_bool                          bNoSolvOpt,
                   real                              print_force)
{
    int            m, negp_pp, negptable, egi, egj;
    real           rtab;
    char          *env;
    double         dbl;
    const t_block *cgs;
    gmx_bool       bGenericKernelOnly;
    gmx_bool       needGroupSchemeTables, bSomeNormalNbListsAreInUse;
    gmx_bool       bFEP_NonBonded;
    int            egp_flags;


    /* By default we turn SIMD kernels on, but it might be turned off further down... */
    fr->use_simd_kernels = TRUE;

    if (check_box(ir->ePBC, box))
    {
        gmx_fatal(FARGS, "%s", check_box(ir->ePBC, box));
    }

    /* Test particle insertion ? */
    if (EI_TPI(ir->eI))
    {
        /* Set to the size of the molecule to be inserted (the last one) */
        /* Because of old style topologies, we have to use the last cg
         * instead of the last molecule type.
         */
        cgs       = &mtop->moltype[mtop->molblock.back().type].cgs;
        fr->n_tpi = cgs->index[cgs->nr] - cgs->index[cgs->nr-1];
        gmx::RangePartitioning molecules = gmx_mtop_molecules(*mtop);
        if (fr->n_tpi != molecules.block(molecules.numBlocks() - 1).size())
        {
            gmx_fatal(FARGS, "The molecule to insert can not consist of multiple charge groups.\nMake it a single charge group.");
        }
    }
    else
    {
        fr->n_tpi = 0;
    }

    if (ir->coulombtype == eelRF_NEC_UNSUPPORTED)
    {
        gmx_fatal(FARGS, "%s electrostatics is no longer supported",
                  eel_names[ir->coulombtype]);
    }

    if (ir->bAdress)
    {
        gmx_fatal(FARGS, "AdResS simulations are no longer supported");
    }
    if (ir->useTwinRange)
    {
        gmx_fatal(FARGS, "Twin-range simulations are no longer supported");
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
        fr->sc_sigma6_min = gmx::power6(ir->fepvals->sc_sigma_min);
    }
    else
    {
        fr->sc_alphacoul  = 0;
        fr->sc_sigma6_min = 0; /* only needed when bScCoul is on */
    }
    fr->sc_power      = ir->fepvals->sc_power;
    fr->sc_r_power    = ir->fepvals->sc_r_power;
    fr->sc_sigma6_def = gmx::power6(ir->fepvals->sc_sigma);

    env = getenv("GMX_SCSIGMA_MIN");
    if (env != nullptr)
    {
        dbl = 0;
        sscanf(env, "%20lf", &dbl);
        fr->sc_sigma6_min = gmx::power6(dbl);
        if (fp)
        {
            fprintf(fp, "Setting the minimum soft core sigma to %g nm\n", dbl);
        }
    }

    fr->bNonbonded = TRUE;
    if (getenv("GMX_NO_NONBONDED") != nullptr)
    {
        /* turn off non-bonded calculations */
        fr->bNonbonded = FALSE;
        GMX_LOG(mdlog.warning).asParagraph().appendText(
                "Found environment variable GMX_NO_NONBONDED.\n"
                "Disabling nonbonded calculations.");
    }

    bGenericKernelOnly = FALSE;

    /* We now check in the NS code whether a particular combination of interactions
     * can be used with water optimization, and disable it if that is not the case.
     */

    if (getenv("GMX_NB_GENERIC") != nullptr)
    {
        if (fp != nullptr)
        {
            fprintf(fp,
                    "Found environment variable GMX_NB_GENERIC.\n"
                    "Disabling all interaction-specific nonbonded kernels, will only\n"
                    "use the slow generic ones in src/gmxlib/nonbonded/nb_generic.c\n\n");
        }
        bGenericKernelOnly = TRUE;
    }

    if (bGenericKernelOnly)
    {
        bNoSolvOpt         = TRUE;
    }

    if ( (getenv("GMX_DISABLE_SIMD_KERNELS") != nullptr) || (getenv("GMX_NOOPTIMIZEDKERNELS") != nullptr) )
    {
        fr->use_simd_kernels = FALSE;
        if (fp != nullptr)
        {
            fprintf(fp,
                    "\nFound environment variable GMX_DISABLE_SIMD_KERNELS.\n"
                    "Disabling the usage of any SIMD-specific non-bonded & bonded kernel routines\n"
                    "(e.g. SSE2/SSE4.1/AVX).\n\n");
        }
    }

    fr->bBHAM = (mtop->ffparams.functype[0] == F_BHAM);

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
        if (fp != nullptr)
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
            gmx_bool bSHAKE;

            bSHAKE = (ir->eConstrAlg == econtSHAKE &&
                      (gmx_mtop_ftype_count(mtop, F_CONSTR) > 0 ||
                       gmx_mtop_ftype_count(mtop, F_CONSTRNC) > 0));

            /* The group cut-off scheme and SHAKE assume charge groups
             * are whole, but not using molpbc is faster in most cases.
             * With intermolecular interactions we need PBC for calculating
             * distances between atoms in different molecules.
             */
            if ((fr->cutoff_scheme == ecutsGROUP || bSHAKE) &&
                !mtop->bIntermolecularInteractions)
            {
                fr->bMolPBC = ir->bPeriodicMols;

                if (bSHAKE && fr->bMolPBC)
                {
                    gmx_fatal(FARGS, "SHAKE is not supported with periodic molecules");
                }
            }
            else
            {
                /* Not making molecules whole is faster in most cases,
                 * but With orientation restraints we need whole molecules.
                 */
                fr->bMolPBC = (fcd->orires.nr == 0);

                if (getenv("GMX_USE_GRAPH") != nullptr)
                {
                    fr->bMolPBC = FALSE;
                    if (fp)
                    {
                        GMX_LOG(mdlog.warning).asParagraph().appendText("GMX_USE_GRAPH is set, using the graph for bonded interactions");
                    }

                    if (mtop->bIntermolecularInteractions)
                    {
                        GMX_LOG(mdlog.warning).asParagraph().appendText("WARNING: Molecules linked by intermolecular interactions have to reside in the same periodic image, otherwise artifacts will occur!");
                    }
                }

                GMX_RELEASE_ASSERT(fr->bMolPBC || !mtop->bIntermolecularInteractions, "We need to use PBC within molecules with inter-molecular interactions");

                if (bSHAKE && fr->bMolPBC)
                {
                    gmx_fatal(FARGS, "SHAKE is not properly supported with intermolecular interactions. For short simulations where linked molecules remain in the same periodic image, the environment variable GMX_USE_GRAPH can be used to override this check.\n");
                }
            }
        }
        else
        {
            fr->bMolPBC = dd_bonded_molpbc(cr->dd, fr->ePBC);
        }
    }

    fr->rc_scaling = ir->refcoord_scaling;
    copy_rvec(ir->posres_com, fr->posres_com);
    copy_rvec(ir->posres_comB, fr->posres_comB);
    fr->rlist                    = cutoff_inf(ir->rlist);
    fr->ljpme_combination_rule   = ir->ljpme_combination_rule;

    /* This now calculates sum for q and c6*/
    bool systemHasNetCharge = set_chargesum(fp, fr, mtop);

    /* fr->ic is used both by verlet and group kernels (to some extent) now */
    init_interaction_const(fp, &fr->ic, ir, mtop, systemHasNetCharge);
    init_interaction_const_tables(fp, fr->ic, ir->rlist + ir->tabext);

    const interaction_const_t *ic = fr->ic;

    /* TODO: Replace this Ewald table or move it into interaction_const_t */
    if (ir->coulombtype == eelEWALD)
    {
        init_ewald_tab(&(fr->ewald_table), ir, fp);
    }

    /* Electrostatics: Translate from interaction-setting-in-mdp-file to kernel interaction format */
    switch (ic->eeltype)
    {
        case eelCUT:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_COULOMB;
            break;

        case eelRF:
        case eelGRF:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_REACTIONFIELD;
            break;

        case eelRF_ZERO:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_REACTIONFIELD;
            GMX_RELEASE_ASSERT(ic->coulomb_modifier == eintmodEXACTCUTOFF, "With the group scheme RF-zero needs the exact cut-off modifier");
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
        case eelP3M_AD:
        case eelEWALD:
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_EWALD;
            break;

        default:
            gmx_fatal(FARGS, "Unsupported electrostatic interaction: %s", eel_names[ic->eeltype]);
    }
    fr->nbkernel_elec_modifier = ic->coulomb_modifier;

    /* Vdw: Translate from mdp settings to kernel format */
    switch (ic->vdwtype)
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
            gmx_fatal(FARGS, "Unsupported vdw interaction: %s", evdw_names[ic->vdwtype]);
    }
    fr->nbkernel_vdw_modifier = ic->vdw_modifier;

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        fr->bvdwtab    = ((ic->vdwtype != evdwCUT || !gmx_within_tol(ic->reppow, 12.0, 10*GMX_DOUBLE_EPS))
                          && !EVDW_PME(ic->vdwtype));
        /* We have special kernels for standard Ewald and PME, but the pme-switch ones are tabulated above */
        fr->bcoultab   = !(ic->eeltype == eelCUT ||
                           ic->eeltype == eelEWALD ||
                           ic->eeltype == eelPME ||
                           ic->eeltype == eelP3M_AD ||
                           ic->eeltype == eelRF ||
                           ic->eeltype == eelRF_ZERO);

        /* If the user absolutely wants different switch/shift settings for coul/vdw, it is likely
         * going to be faster to tabulate the interaction than calling the generic kernel.
         * However, if generic kernels have been requested we keep things analytically.
         */
        if (fr->nbkernel_elec_modifier == eintmodPOTSWITCH &&
            fr->nbkernel_vdw_modifier == eintmodPOTSWITCH &&
            !bGenericKernelOnly)
        {
            if ((ic->rcoulomb_switch != ic->rvdw_switch) || (ic->rcoulomb != ic->rvdw))
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
            if ((ic->rcoulomb != ic->rvdw) && (!bGenericKernelOnly))
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
            fprintf(fp, "Table routines are used for coulomb: %s\n",
                    gmx::boolToString(fr->bcoultab));
            fprintf(fp, "Table routines are used for vdw:     %s\n",
                    gmx::boolToString(fr->bvdwtab));
        }

        if (fr->bvdwtab)
        {
            fr->nbkernel_vdw_interaction = GMX_NBKERNEL_VDW_CUBICSPLINETABLE;
            fr->nbkernel_vdw_modifier    = eintmodNONE;
        }
        if (fr->bcoultab)
        {
            fr->nbkernel_elec_interaction = GMX_NBKERNEL_ELEC_CUBICSPLINETABLE;
            fr->nbkernel_elec_modifier    = eintmodNONE;
        }
    }

    if (ir->cutoff_scheme == ecutsVERLET)
    {
        if (!gmx_within_tol(ic->reppow, 12.0, 10*GMX_DOUBLE_EPS))
        {
            gmx_fatal(FARGS, "Cut-off scheme %s only supports LJ repulsion power 12", ecutscheme_names[ir->cutoff_scheme]);
        }
        /* Older tpr files can contain Coulomb user tables with the Verlet cutoff-scheme,
         * while mdrun does not (and never did) support this.
         */
        if (EEL_USER(fr->ic->eeltype))
        {
            gmx_fatal(FARGS, "Combination of %s and cutoff scheme %s is not supported",
                      eel_names[ir->coulombtype], ecutscheme_names[ir->cutoff_scheme]);
        }

        fr->bvdwtab  = FALSE;
        fr->bcoultab = FALSE;
    }

    /* 1-4 interaction electrostatics */
    fr->fudgeQQ = mtop->ffparams.fudgeQQ;

    /* Parameters for generalized RF */
    fr->zsquare = 0.0;
    fr->temp    = 0.0;

    if (ic->eeltype == eelGRF)
    {
        init_generalized_rf(fp, mtop, ir, fr);
    }

    fr->haveDirectVirialContributions =
        (EEL_FULL(ic->eeltype) || EVDW_PME(ic->vdwtype) ||
         fr->forceProviders->hasForceProvider() ||
         gmx_mtop_ftype_count(mtop, F_POSRES) > 0 ||
         gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0 ||
         ir->nwall > 0 ||
         ir->bPull ||
         ir->bRot ||
         ir->bIMD);

    if (fr->haveDirectVirialContributions)
    {
        fr->forceBufferForDirectVirialContributions = new std::vector<gmx::RVec>;
    }

    if (fr->shift_vec == nullptr)
    {
        snew(fr->shift_vec, SHIFTS);
    }

    if (fr->fshift == nullptr)
    {
        snew(fr->fshift, SHIFTS);
    }

    if (fr->nbfp == nullptr)
    {
        fr->ntype = mtop->ffparams.atnr;
        fr->nbfp  = mk_nbfp(&mtop->ffparams, fr->bBHAM);
        if (EVDW_PME(ic->vdwtype))
        {
            fr->ljpme_c6grid  = make_ljpme_c6grid(&mtop->ffparams, fr);
        }
    }

    /* Copy the energy group exclusions */
    fr->egp_flags = ir->opts.egp_flags;

    /* Van der Waals stuff */
    if ((ic->vdwtype != evdwCUT) && (ic->vdwtype != evdwUSER) && !fr->bBHAM)
    {
        if (ic->rvdw_switch >= ic->rvdw)
        {
            gmx_fatal(FARGS, "rvdw_switch (%f) must be < rvdw (%f)",
                      ic->rvdw_switch, ic->rvdw);
        }
        if (fp)
        {
            fprintf(fp, "Using %s Lennard-Jones, switch between %g and %g nm\n",
                    (ic->eeltype == eelSWITCH) ? "switched" : "shifted",
                    ic->rvdw_switch, ic->rvdw);
        }
    }

    if (fr->bBHAM && EVDW_PME(ic->vdwtype))
    {
        gmx_fatal(FARGS, "LJ PME not supported with Buckingham");
    }

    if (fr->bBHAM && (ic->vdwtype == evdwSHIFT || ic->vdwtype == evdwSWITCH))
    {
        gmx_fatal(FARGS, "Switch/shift interaction not supported with Buckingham");
    }

    if (fr->bBHAM && fr->cutoff_scheme == ecutsVERLET)
    {
        gmx_fatal(FARGS, "Verlet cutoff-scheme is not supported with Buckingham");
    }

    if (fp && fr->cutoff_scheme == ecutsGROUP)
    {
        fprintf(fp, "Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
                fr->rlist, ic->rcoulomb, fr->bBHAM ? "BHAM" : "LJ", ic->rvdw);
    }

    if (ir->implicit_solvent)
    {
        gmx_fatal(FARGS, "Implict solvation is no longer supported.");
    }

    /* Construct tables for the group scheme. A little unnecessary to
     * make both vdw and coul tables sometimes, but what the
     * heck. Note that both cutoff schemes construct Ewald tables in
     * init_interaction_const_tables. */
    needGroupSchemeTables = (ir->cutoff_scheme == ecutsGROUP &&
                             (fr->bcoultab || fr->bvdwtab));

    negp_pp   = ir->opts.ngener - ir->nwall;
    negptable = 0;
    if (!needGroupSchemeTables)
    {
        bSomeNormalNbListsAreInUse = TRUE;
    }
    else
    {
        bSomeNormalNbListsAreInUse = FALSE;
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
    }

    /* This code automatically gives table length tabext without cut-off's,
     * in that case grompp should already have checked that we do not need
     * normal tables and we only generate tables for 1-4 interactions.
     */
    rtab = ir->rlist + ir->tabext;

    if (needGroupSchemeTables)
    {
        /* make tables for ordinary interactions */
        if (bSomeNormalNbListsAreInUse)
        {
            make_nbf_tables(fp, ic, rtab, tabfn, nullptr, nullptr, nullptr);
            m = 1;
        }
        else
        {
            m = 0;
        }
        if (negptable > 0)
        {
            /* Read the special tables for certain energy group pairs */
            gmx::ArrayRef<const int> nm_ind = mtop->groups.groups[SimulationAtomGroupType::EnergyOutput];
            for (egi = 0; egi < negp_pp; egi++)
            {
                for (egj = egi; egj < negp_pp; egj++)
                {
                    egp_flags = ir->opts.egp_flags[GID(egi, egj, ir->opts.ngener)];
                    if ((egp_flags & EGP_TABLE) && !(egp_flags & EGP_EXCL))
                    {
                        /* Read the table file with the two energy groups names appended */
                        make_nbf_tables(fp, ic, rtab, tabfn,
                                        *mtop->groups.groupNames[nm_ind[egi]],
                                        *mtop->groups.groupNames[nm_ind[egj]],
                                        nullptr);
                        m++;
                    }
                }
            }
        }
    }

    /* We want to use unmodified tables for 1-4 coulombic
     * interactions, so we must in general have an extra set of
     * tables. */
    if (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 ||
        gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0 ||
        gmx_mtop_ftype_count(mtop, F_LJC_PAIRS_NB) > 0)
    {
        fr->pairsTable = make_tables(fp, ic, tabpfn, rtab,
                                     GMX_MAKETABLES_14ONLY);
    }

    /* Wall stuff */
    fr->nwall = ir->nwall;
    if (ir->nwall && ir->wall_type == ewtTABLE)
    {
        make_wall_tables(fp, ir, tabfn, &mtop->groups, fr);
    }

    if (fcd && !tabbfnm.empty())
    {
        // Need to catch std::bad_alloc
        // TODO Don't need to catch this here, when merging with master branch
        try
        {
            fcd->bondtab  = make_bonded_tables(fp,
                                               F_TABBONDS, F_TABBONDSNC,
                                               mtop, tabbfnm, "b");
            fcd->angletab = make_bonded_tables(fp,
                                               F_TABANGLES, -1,
                                               mtop, tabbfnm, "a");
            fcd->dihtab   = make_bonded_tables(fp,
                                               F_TABDIHS, -1,
                                               mtop, tabbfnm, "d");
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    else
    {
        if (debug)
        {
            fprintf(debug, "No fcdata or table file name passed, can not read table, can not do bonded interactions\n");
        }
    }

    // QM/MM initialization if requested
    fr->bQMMM = ir->bQMMM;
    if (fr->bQMMM)
    {
        // Initialize QM/MM if supported
        if (GMX_QMMM)
        {
            GMX_LOG(mdlog.info).asParagraph().
                appendText("Large parts of the QM/MM support is deprecated, and may be removed in a future "
                           "version. Please get in touch with the developers if you find the support useful, "
                           "as help is needed if the functionality is to continue to be available.");
            fr->qr = mk_QMMMrec();
            init_QMMMrec(cr, mtop, ir, fr);
        }
        else
        {
            gmx_incons("QM/MM was requested, but is only available when GROMACS "
                       "is configured with QM/MM support");
        }
    }

    /* Set all the static charge group info */
    fr->cginfo_mb = init_cginfo_mb(fp, mtop, fr, bNoSolvOpt,
                                   &bFEP_NonBonded,
                                   &fr->bExcl_IntraCGAll_InterCGNone);
    if (!DOMAINDECOMP(cr))
    {
        fr->cginfo = cginfo_expand(mtop->molblock.size(), fr->cginfo_mb);
    }

    if (!DOMAINDECOMP(cr))
    {
        forcerec_set_ranges(fr, ncg_mtop(mtop), ncg_mtop(mtop),
                            mtop->natoms, mtop->natoms, mtop->natoms);
    }

    fr->print_force = print_force;

    /* Initialize the thread working data for bonded interactions */
    fr->bondedThreading =
        init_bonded_threading(fp, mtop->groups.groups[SimulationAtomGroupType::EnergyOutput].size());

    fr->nthread_ewc = gmx_omp_nthreads_get(emntBonded);
    snew(fr->ewc_t, fr->nthread_ewc);

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        // We checked the cut-offs in grompp, but double-check here.
        // We have PME+LJcutoff kernels for rcoulomb>rvdw.
        if (EEL_PME_EWALD(ir->coulombtype) && ir->vdwtype == eelCUT)
        {
            GMX_RELEASE_ASSERT(ir->rcoulomb >= ir->rvdw, "With Verlet lists and PME we should have rcoulomb>=rvdw");
        }
        else
        {
            GMX_RELEASE_ASSERT(ir->rcoulomb == ir->rvdw, "With Verlet lists and no PME rcoulomb and rvdw should be identical");
        }

        fr->nbv = Nbnxm::init_nb_verlet(mdlog, bFEP_NonBonded, ir, fr,
                                        cr, hardwareInfo, deviceInfo,
                                        mtop, box);

        if (useGpuForBonded)
        {
            auto stream = DOMAINDECOMP(cr) ?
                Nbnxm::gpu_get_command_stream(fr->nbv->gpu_nbv, Nbnxm::InteractionLocality::NonLocal) :
                Nbnxm::gpu_get_command_stream(fr->nbv->gpu_nbv, Nbnxm::InteractionLocality::Local);
            // TODO the heap allocation is only needed while
            // t_forcerec lacks a constructor.
            fr->gpuBonded = new gmx::GpuBonded(mtop->ffparams,
                                               stream);
        }
    }

    if (ir->eDispCorr != edispcNO)
    {
        fr->dispersionCorrection =
            std::make_unique<DispersionCorrection>(*mtop, *ir, fr->bBHAM,
                                                   fr->ntype,
                                                   gmx::arrayRefFromArray(fr->nbfp, fr->ntype*fr->ntype*2),
                                                   *fr->ic, tabfn);
        fr->dispersionCorrection->print(mdlog);
    }

    if (fp != nullptr)
    {
        /* Here we switch from using mdlog, which prints the newline before
         * the paragraph, to our old fprintf logging, which prints the newline
         * after the paragraph, so we should add a newline here.
         */
        fprintf(fp, "\n");
    }
}

/* Frees GPU memory and sets a tMPI node barrier.
 *
 * Note that this function needs to be called even if GPUs are not used
 * in this run because the PME ranks have no knowledge of whether GPUs
 * are used or not, but all ranks need to enter the barrier below.
 * \todo Remove physical node barrier from this function after making sure
 * that it's not needed anymore (with a shared GPU run).
 */
void free_gpu_resources(t_forcerec                          *fr,
                        const gmx::PhysicalNodeCommunicator &physicalNodeCommunicator)
{
    bool isPPrankUsingGPU = (fr != nullptr) && (fr->nbv != nullptr) && fr->nbv->useGpu();

    /* stop the GPU profiler (only CUDA) */
    stopGpuProfiler();

    if (isPPrankUsingGPU)
    {
        /* Free data in GPU memory and pinned memory before destroying the GPU context */
        fr->nbv.reset();

        delete fr->gpuBonded;
        fr->gpuBonded = nullptr;
    }

    /* With tMPI we need to wait for all ranks to finish deallocation before
     * destroying the CUDA context in free_gpu() as some tMPI ranks may be sharing
     * GPU and context.
     *
     * This is not a concern in OpenCL where we use one context per rank which
     * is freed in nbnxn_gpu_free().
     *
     * Note: it is safe to not call the barrier on the ranks which do not use GPU,
     * but it is easier and more futureproof to call it on the whole node.
     */
    if (GMX_THREAD_MPI)
    {
        physicalNodeCommunicator.barrier();
    }
}

void done_forcerec(t_forcerec *fr, int numMolBlocks)
{
    if (fr == nullptr)
    {
        // PME-only ranks don't have a forcerec
        return;
    }
    done_cginfo_mb(fr->cginfo_mb, numMolBlocks);
    sfree(fr->nbfp);
    done_interaction_const(fr->ic);
    sfree(fr->shift_vec);
    sfree(fr->fshift);
    sfree(fr->ewc_t);
    tear_down_bonded_threading(fr->bondedThreading);
    GMX_RELEASE_ASSERT(fr->gpuBonded == nullptr, "Should have been deleted earlier, when used");
    fr->bondedThreading = nullptr;
    delete fr;
}
