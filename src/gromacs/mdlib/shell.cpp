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

#include "gmxpre.h"

#include "shellfc.h"

#include <stdio.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"


/* Routines to send/recieve coordinates and force
 * of connected atoms.
 */

static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

static int *atom2cg(t_block *cgs)
{
    int *a2cg, cg, i;

    snew(a2cg, cgs->index[cgs->nr]);
    for (cg = 0; cg < cgs->nr; cg++)
    {
        for (i = cgs->index[cg]; i < cgs->index[cg+1]; i++)
        {
            a2cg[i] = cg;
        }
    }

    return a2cg;
}

static int count_intercg_shells(const gmx_mtop_t *mtop,
                                gmx_bool   *bHaveChargeGroups)
{
    int             mb, ftype, nral, i, cg, a;
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;
    int            *a2cg;
    t_ilist        *il;
    t_iatom        *ia;
    int             n_intercg_shells;

    *bHaveChargeGroups = FALSE;

    n_intercg_shells = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];

        if (molt->cgs.nr < molt->atoms.nr)
        {
            *bHaveChargeGroups = TRUE;
        }

        a2cg = atom2cg(&molt->cgs);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (ftype == F_DRUDEBONDS || ftype == F_POLARIZATION || ftype == F_HYPER_POL)
            {
                nral = NRAL(ftype);
                il   = &molt->ilist[ftype];
                ia   = il->iatoms;
                for (i = 0; i < il->nr; i += 1+nral)
                {
                    cg = a2cg[ia[1+i]];
                    for (a = 1; a < nral; a++)
                    {
                        if (a2cg[ia[1+a]] != cg)
                        {
                            n_intercg_shells += molb->nmol;
                            break;
                        }
                    }
                }
            }
        }
        sfree(a2cg);
    }

    return n_intercg_shells;
}

static int **get_shell_pbc(t_ilist *ilist,
                           t_atom *atom, t_mdatoms *md,
                           t_block *cgs, int *a2cg)
{
    int      ftype, nral, i, shi, cg_v, a;
    int      shell = -1;
    t_ilist *il;
    t_iatom *ia;
    int    **shell_pbc, *shell_pbc_f;
    char    *pbc_set;
    gmx_bool bShellOnlyCG_and_FirstAtom;

    if (md == NULL)
    {
        return 0;
    }

    /* Make an array that tells if the pbc of an atom is set */
    snew(pbc_set, cgs->index[cgs->nr]);
    /* PBC is set for all non shells */
    for (a = 0; a < cgs->index[cgs->nr]; a++)
    {
        if ((atom && atom[a].ptype != eptShell) ||
            (md   && md->ptype[a]  != eptShell))
        {
            pbc_set[a] = 1;
        }
    }

    /* There are several possible interactions for shells/Drudes, including
     * everything from F_DRUDEBONDS..F_HYPER_POL, inclusive */
    snew(shell_pbc, F_HYPER_POL-F_DRUDEBONDS+1);

    for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
    {
        nral = NRAL(ftype);
        il   = &ilist[ftype];
        ia   = il->iatoms;

        snew(shell_pbc[ftype-F_DRUDEBONDS], il->nr/(1+nral));
        shell_pbc_f = shell_pbc[ftype-F_DRUDEBONDS];

        i = 0;
        while (i < il->nr)
        {
            shi   = i/(1+nral);
            if (md && (md->ptype[ia[i+1]] == eptShell))
            {
                shell = ia[i+1];
            }
            else if (md && (md->ptype[ia[i+2]] == eptShell))
            {
                shell = ia[i+2];
            }

            if (shell != -1)
            {
                cg_v  = a2cg[shell];
                /* A value of -2 signals that this shell and its bonded 
                 * atoms are all within the same cg, so no pbc is required.
                 */
                shell_pbc_f[shi] = -2;
                /* Check if connected atoms are outside the shell's cg */
                for (a = 1; a < nral; a++)
                {
                    if (a2cg[shell+a] != cg_v)
                    {
                        shell_pbc_f[shi] = -1;
                    }
                }
                if (shell_pbc_f[shi] == -1)
                {
                    /* Check if this is the first processed atom of a shell-only cg */
                    bShellOnlyCG_and_FirstAtom = TRUE;
                    for (a = cgs->index[cg_v]; a < cgs->index[cg_v+1]; a++)
                    {
                        /* Non-shells already have pbc set, so simply check for pbc_set */
                        if (pbc_set[a])
                        {
                            bShellOnlyCG_and_FirstAtom = FALSE;
                            break;
                        }
                    }
                    if (bShellOnlyCG_and_FirstAtom)
                    {
                        /* First processed atom of a shell-only charge group.
                         * The pbc of the input coordinates should be preserved.
                         */
                        shell_pbc_f[shi] = shell;
                    }
                    else if (cg_v != a2cg[shell+1])
                    {
                        /* This shell has a different charge group index
                         * than its first connected atom
                         * and the charge group has more than one atom,
                         * search for the first normal particle
                         * or shell that already had its pbc defined.
                         * If nothing is found, use full pbc for this shell.
                         */
                        for (a = cgs->index[cg_v]; a < cgs->index[cg_v+1]; a++)
                        {
                            if (a != shell && pbc_set[a])
                            {
                                shell_pbc_f[shi] = a;
                                if (gmx_debug_at)
                                {
                                    fprintf(debug, "shell %d match pbc with atom %d\n",
                                            shell+1, a+1);
                                }
                                break;
                            }
                        }
                        if (gmx_debug_at)
                        {
                            fprintf(debug, "shell atom %d  cg %d - %d pbc atom %d\n",
                                    shell+1, cgs->index[cg_v]+1, cgs->index[cg_v+1],
                                    shell_pbc_f[shi]+1);
                        }
                    }
                }
                i += 1+nral;

                /* This shell now has its pbc defined */
                pbc_set[shell] = 1;
            }
        }
    }

    sfree(pbc_set);

    return shell_pbc;
}

gmx_shellfc_t init_shell(const gmx_mtop_t *mtop, t_commrec *cr,
                         gmx_bool bSerial_NoPBC)
{
    int             nshell, i;
    int            *a2cg;
    gmx_shellfc_t   shfc;
    int             mt;
    gmx_moltype_t  *molt;

    /* check if there are shells */
    nshell = 0;
    for (i = 0; i < F_NRE; i++)
    {
        /* One shell/Drude per interaction, so just count them */
        /* All other shell/Drude functions are a subset of one
         * of these, so we do not need to consider all of them explicitly */
        if ((i == F_DRUDEBONDS) || (i == F_POLARIZATION) || (i == F_HYPER_POL))
        {
            nshell += gmx_mtop_ftype_count(mtop, i);
        }
    }

    if (nshell == 0)
    {
        return NULL;
    }

    if (debug)
    {
        fprintf(debug, "INIT SHELL: There are %d shells\n", nshell);
    }

    snew(shfc, 1);

    shfc->n_intercg_shells = count_intercg_shells(mtop, &shfc->bInterCG);

    if (shfc->n_intercg_shells > 0)
    {
        shfc->bInterCG = TRUE;
    }

    if (debug)
    {
        fprintf(debug, "INIT SHELL: found %d inter-cg shells\n", shfc->n_intercg_shells);
    }

    shfc->bHaveChargeGroups = (ncg_mtop(mtop) < mtop->natoms);

    /* If we don't have charge groups, the shell follows its own pbc */
    if (!bSerial_NoPBC &&
        shfc->bInterCG &&
        shfc->bHaveChargeGroups &&
        shfc->n_intercg_shells > 0 && DOMAINDECOMP(cr))
    {
        shfc->nshell_pbc_molt = mtop->nmoltype;
        snew(shfc->shell_pbc_molt, shfc->nshell_pbc_molt);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            molt = &mtop->moltype[mt];
            /* Make an atom to charge group index */
            a2cg = atom2cg(&molt->cgs);
            /* TODO: check */
            shfc->shell_pbc_molt[mt] = get_shell_pbc(molt->ilist,
                                                     molt->atoms.atom, NULL,
                                                     &molt->cgs, a2cg);
            sfree(a2cg);
        }

        snew(shfc->shell_pbc_loc_nalloc, (F_HYPER_POL-F_DRUDEBONDS+1));
        snew(shfc->shell_pbc_loc, (F_HYPER_POL-F_DRUDEBONDS+1));
    }

    if (bSerial_NoPBC)
    {
        shfc->nthreads = 1;
    }
    else
    {
        shfc->nthreads = gmx_omp_nthreads_get(emntSHELL);
    }
    if (!bSerial_NoPBC)
    {
        /* We need one extra thread data structure for the overlap shells */
        snew(shfc->tdata, shfc->nthreads+1);
    }

    shfc->th_ind        = NULL;
    shfc->th_ind_nalloc = 0;

    return shfc;
}

static void prepare_shell_thread(const t_ilist      *ilist,
                                 gmx_shell_thread_t *shell_th)
{
    int ftype;

    for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
    {
        if (ilist[ftype].nr > shell_th->ilist[ftype].nalloc)
        {
            shell_th->ilist[ftype].nalloc = over_alloc_large(ilist[ftype].nr);
            srenew(shell_th->ilist[ftype].iatoms, shell_th->ilist[ftype].nalloc);
        }

        shell_th->ilist[ftype].nr = 0;
    }
}

void split_shells_over_threads(const t_ilist   *ilist,
                               const t_mdatoms *mdatoms,
                               gmx_bool         bLimitRange,
                               gmx_shellfc_t    shfc)
{
    int      th;
    int      shell_atom_range, natperthread;
    int     *th_ind;
    int      ftype;
    t_iatom *iat;
    t_ilist *il_th;
    int      nral1, inc, i, j;

    if (shfc->nthreads == 1)
    {
        /* Nothing to do */
        return;
    }

#pragma omp parallel for num_threads(shfc->nthreads) schedule(static)
    for (th = 0; th < shfc->nthreads; th++)
    {
        prepare_shell_thread(ilist, &shfc->tdata[th]);
    }
    /* Master threads does the (potential) overlap shells */
    prepare_shell_thread(ilist, &shfc->tdata[shfc->nthreads]);

    /* The current way of distributing the shells over threads in primitive.
     * We divide the atom range 0 - natoms_in_shell uniformly over threads,
     * without taking into account how the shells are distributed.
     * Without domain decomposition we bLimitRange=TRUE and we at least
     * tighten the upper bound of the range (useful for common systems)
     */
    if (bLimitRange)
    {
        shell_atom_range = -1;
        for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
        {
            nral1 = 1 + NRAL(ftype);
            iat   = ilist[ftype].iatoms;
            for (i = 0; i < ilist[ftype].nr; i += nral1)
            {
                for (j = i+1; j < i+nral1; j++)
                {
                    shell_atom_range = std::max(shell_atom_range, iat[j]);
                }
            }
        }
        shell_atom_range++;
    }
    else
    {
        shell_atom_range = mdatoms->homenr;
    }
    natperthread = (shell_atom_range + shfc->nthreads - 1)/shfc->nthreads;

    if (debug)
    {
        fprintf(debug, "shell thread dist: natoms %d, range %d, natperthread %d\n", mdatoms->nr, shell_atom_range, natperthread);
    }

    /* To simplify the shell assignment, we make an index which tells us
     * to which thread particles, both non-shells and shells, are assigned.
     */
    if (mdatoms->nr > shfc->th_ind_nalloc)
    {
        shfc->th_ind_nalloc = over_alloc_large(mdatoms->nr);
        srenew(shfc->th_ind, shfc->th_ind_nalloc);
    }
    th_ind = shfc->th_ind;
    th     = 0;
    for (i = 0; i < mdatoms->nr; i++)
    {
        if (mdatoms->ptype[i] == eptShell)
        {
            /* shells are not assigned to a thread yet */
            th_ind[i] = -1;
        }
        else
        {
            /* assign non-shell particles to thread th */
            th_ind[i] = th;
        }
        if (i == (th + 1)*natperthread && th < shfc->nthreads)
        {
            th++;
        }
    }

    for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
    {
        nral1 = 1 + NRAL(ftype);
        inc   = nral1;
        iat   = ilist[ftype].iatoms;
        for (i = 0; i < ilist[ftype].nr; )
        {
            th = iat[1+i]/natperthread;
            for (j = i+2; j < i+nral1; j++)
            {
                if (th_ind[iat[j]] != th)
                {
                    /* Some constructing atoms are not assigned to
                     * thread th, move this shell to a separate batch.
                     */
                    th = shfc->nthreads;
                }
            }
            /* Copy this shell to the thread data struct of thread th */
            il_th = &shfc->tdata[th].ilist[ftype];
            for (j = i; j < i+inc; j++)
            {
                il_th->iatoms[il_th->nr++] = iat[j];
            }
            /* Update this shell's thread index entry */
            th_ind[iat[1+i]] = th;

            i += inc;
        }
    }

    if (debug)
    {
        for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
        {
            fprintf(debug, "%-20s thread dist:",
                    interaction_function[ftype].longname);
            for (th = 0; th < shfc->nthreads+1; th++)
            {
                fprintf(debug, " %4d", shfc->tdata[th].ilist[ftype].nr);
            }
            fprintf(debug, "\n");
        }
    }
}

void set_shell_top(gmx_shellfc_t shfc, gmx_localtop_t *top, t_mdatoms *md,
                   t_commrec *cr)
{
    int *a2cg;

    if (shfc->n_intercg_shells > 0)
    {
        if (shfc->bInterCG)
        {
            /* Make an atom to charge group index */
            a2cg                 = atom2cg(&top->cgs);
            shfc->shell_pbc_loc = get_shell_pbc(top->idef.il, NULL, md,
                                                &top->cgs, a2cg);
            sfree(a2cg);
        }
    }

    if (shfc->nthreads > 1)
    {
        if (shfc->bInterCG)
        {
            gmx_fatal(FARGS, "The combination of threading, shells and charge groups is not implemented");
        }

        split_shells_over_threads(top->idef.il, md, !DOMAINDECOMP(cr), shfc);
    }
}

