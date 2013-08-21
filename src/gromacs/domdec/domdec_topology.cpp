/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 *
 * \brief This file defines functions used by the domdec module
 * while managing the construction, use and error checking for
 * topologies local to a DD rank.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/legacyheaders/chargegroup.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_ga2la.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_vsite.h"

/*! \brief The number of integer item in the local state, used for broadcasting of the state */
#define NITEM_DD_INIT_LOCAL_STATE 5

typedef struct {
    int  *index;  /* Index for each atom into il                  */
    int  *il;     /* ftype|type|a0|...|an|ftype|...               */
} reverse_ilist_t;

typedef struct {
    int  a_start;
    int  a_end;
    int  natoms_mol;
    int  type;
} molblock_ind_t;

/*! \brief Struct for thread local work data for local topology generation */
typedef struct {
    t_idef     idef;             /**< Partial local topology */
    int      **vsite_pbc;        /**< vsite PBC structure */
    int       *vsite_pbc_nalloc; /**< Allocation sizes for vsite_pbc */
    int        nbonded;          /**< The number of bondeds in this struct */
    t_blocka   excl;             /**< List of exclusions */
    int        excl_count;       /**< The total exclusion count for \p excl */
} thread_work_t;

/*! \brief Struct for the reverse topology: links bonded interactions to atomsx */
typedef struct gmx_reverse_top {
    //! @cond Doxygen_Suppress
    gmx_bool         bExclRequired;               /**< Do we require all exclusions to be assigned? */
    int              n_excl_at_max;               /**< The maximum number of exclusions one atom can have */
    gmx_bool         bConstr;                     /**< Are there constraints in this revserse top?  */
    gmx_bool         bSettle;                     /**< Are there settles in this revserse top?  */
    gmx_bool         bBCheck;                     /**< All bonded interactions have to be assigned? */
    gmx_bool         bMultiCGmols;                /**< Are the multi charge-group molecules?        */
    reverse_ilist_t *ril_mt;                      /**< Reverse ilist for all moltypes      */
    int              ril_mt_tot_size;             /**< The size of ril_mt[?].index summed over all entries */
    int              ilsort;                      /**< The sorting state of bondeds for free energy */
    molblock_ind_t  *mbi;                         /**< molblock to global atom index for quick lookup of molblocks on atom index */
    int              nmolblock;                   /**< The number of molblocks, size of \p mbi */

    gmx_bool         bIntermolecularInteractions; /**< Do we have intermolecular interactions? */
    reverse_ilist_t  ril_intermol;                /**< Intermolecular reverse ilist */

    /* Work data structures for multi-threading */
    int            nthread;           /**< The number of threads to be used */
    thread_work_t *th_work;           /**< Thread work array for local topology generation */

    /* Pointers only used for an error message */
    gmx_mtop_t     *err_top_global; /**< Pointer to the global top, only used for error reporting */
    gmx_localtop_t *err_top_local;  /**< Pointer to the local top, only used for error reporting */
    //! @endcond
} gmx_reverse_top_t;

/*! \brief Returns the number of atom entries for il in gmx_reverse_top_t */
static int nral_rt(int ftype)
{
    int nral;

    nral = NRAL(ftype);
    if (interaction_function[ftype].flags & IF_VSITE)
    {
        /* With vsites the reverse topology contains
         * two extra entries for PBC.
         */
        nral += 2;
    }

    return nral;
}

/*! \brief Return whether interactions of type \p ftype need to be assigned exactly once */
static gmx_bool dd_check_ftype(int ftype, gmx_bool bBCheck,
                               gmx_bool bConstr, gmx_bool bSettle)
{
    return (((interaction_function[ftype].flags & IF_BOND) &&
             !(interaction_function[ftype].flags & IF_VSITE) &&
             (bBCheck || !(interaction_function[ftype].flags & IF_LIMZERO))) ||
            (bConstr && (ftype == F_CONSTR || ftype == F_CONSTRNC)) ||
            (bSettle && ftype == F_SETTLE));
}

/*! \brief Print a header on error messages */
static void print_error_header(FILE *fplog, char *moltypename, int nprint)
{
    fprintf(fplog, "\nMolecule type '%s'\n", moltypename);
    fprintf(stderr, "\nMolecule type '%s'\n", moltypename);
    fprintf(fplog,
            "the first %d missing interactions, except for exclusions:\n",
            nprint);
    fprintf(stderr,
            "the first %d missing interactions, except for exclusions:\n",
            nprint);
}

/*! \brief Help print error output when interactions are missing */
static void print_missing_interactions_mb(FILE *fplog, t_commrec *cr,
                                          gmx_reverse_top_t *rt,
                                          char *moltypename,
                                          reverse_ilist_t *ril,
                                          int a_start, int a_end,
                                          int nat_mol, int nmol,
                                          t_idef *idef)
{
    int      nril_mol, *assigned, *gatindex;
    int      ftype, ftype_j, nral, i, j_mol, j, a0, a0_mol, mol, a;
    int      nprint;
    t_ilist *il;
    t_iatom *ia;
    gmx_bool bFound;

    nril_mol = ril->index[nat_mol];
    snew(assigned, nmol*nril_mol);

    gatindex = cr->dd->gatindex;
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rt->bBCheck, rt->bConstr, rt->bSettle))
        {
            nral = NRAL(ftype);
            il   = &idef->il[ftype];
            ia   = il->iatoms;
            for (i = 0; i < il->nr; i += 1+nral)
            {
                a0     = gatindex[ia[1]];
                /* Check if this interaction is in
                 * the currently checked molblock.
                 */
                if (a0 >= a_start && a0 < a_end)
                {
                    mol    = (a0 - a_start)/nat_mol;
                    a0_mol = (a0 - a_start) - mol*nat_mol;
                    j_mol  = ril->index[a0_mol];
                    bFound = FALSE;
                    while (j_mol < ril->index[a0_mol+1] && !bFound)
                    {
                        j       = mol*nril_mol + j_mol;
                        ftype_j = ril->il[j_mol];
                        /* Here we need to check if this interaction has
                         * not already been assigned, since we could have
                         * multiply defined interactions.
                         */
                        if (ftype == ftype_j && ia[0] == ril->il[j_mol+1] &&
                            assigned[j] == 0)
                        {
                            /* Check the atoms */
                            bFound = TRUE;
                            for (a = 0; a < nral; a++)
                            {
                                if (gatindex[ia[1+a]] !=
                                    a_start + mol*nat_mol + ril->il[j_mol+2+a])
                                {
                                    bFound = FALSE;
                                }
                            }
                            if (bFound)
                            {
                                assigned[j] = 1;
                            }
                        }
                        j_mol += 2 + nral_rt(ftype_j);
                    }
                    if (!bFound)
                    {
                        gmx_incons("Some interactions seem to be assigned multiple times");
                    }
                }
                ia += 1 + nral;
            }
        }
    }

    gmx_sumi(nmol*nril_mol, assigned, cr);

    nprint = 10;
    i      = 0;
    for (mol = 0; mol < nmol; mol++)
    {
        j_mol = 0;
        while (j_mol < nril_mol)
        {
            ftype = ril->il[j_mol];
            nral  = NRAL(ftype);
            j     = mol*nril_mol + j_mol;
            if (assigned[j] == 0 &&
                !(interaction_function[ftype].flags & IF_VSITE))
            {
                if (DDMASTER(cr->dd))
                {
                    if (i == 0)
                    {
                        print_error_header(fplog, moltypename, nprint);
                    }
                    fprintf(fplog, "%20s atoms",
                            interaction_function[ftype].longname);
                    fprintf(stderr, "%20s atoms",
                            interaction_function[ftype].longname);
                    for (a = 0; a < nral; a++)
                    {
                        fprintf(fplog, "%5d", ril->il[j_mol+2+a]+1);
                        fprintf(stderr, "%5d", ril->il[j_mol+2+a]+1);
                    }
                    while (a < 4)
                    {
                        fprintf(fplog, "     ");
                        fprintf(stderr, "     ");
                        a++;
                    }
                    fprintf(fplog, " global");
                    fprintf(stderr, " global");
                    for (a = 0; a < nral; a++)
                    {
                        fprintf(fplog, "%6d",
                                a_start+mol*nat_mol+ril->il[j_mol+2+a]+1);
                        fprintf(stderr, "%6d",
                                a_start+mol*nat_mol+ril->il[j_mol+2+a]+1);
                    }
                    fprintf(fplog, "\n");
                    fprintf(stderr, "\n");
                }
                i++;
                if (i >= nprint)
                {
                    break;
                }
            }
            j_mol += 2 + nral_rt(ftype);
        }
    }

    sfree(assigned);
}

/*! \brief Help print error output when interactions are missing */
static void print_missing_interactions_atoms(FILE *fplog, t_commrec *cr,
                                             gmx_mtop_t *mtop, t_idef *idef)
{
    int                mb, a_start, a_end;
    gmx_molblock_t    *molb;
    gmx_reverse_top_t *rt;

    rt = cr->dd->reverse_top;

    /* Print the atoms in the missing interactions per molblock */
    a_end = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb    = &mtop->molblock[mb];
        a_start = a_end;
        a_end   = a_start + molb->nmol*molb->natoms_mol;

        print_missing_interactions_mb(fplog, cr, rt,
                                      *(mtop->moltype[molb->type].name),
                                      &rt->ril_mt[molb->type],
                                      a_start, a_end, molb->natoms_mol,
                                      molb->nmol,
                                      idef);
    }
}

void dd_print_missing_interactions(FILE *fplog, t_commrec *cr, int local_count,  gmx_mtop_t *top_global, t_state *state_local)
{
    int             ndiff_tot, cl[F_NRE], n, ndiff, rest_global, rest_local;
    int             ftype, nral;
    char            buf[STRLEN];
    gmx_domdec_t   *dd;
    gmx_mtop_t     *err_top_global;
    gmx_localtop_t *err_top_local;

    dd = cr->dd;

    err_top_global = dd->reverse_top->err_top_global;
    err_top_local  = dd->reverse_top->err_top_local;

    if (fplog)
    {
        fprintf(fplog, "\nNot all bonded interactions have been properly assigned to the domain decomposition cells\n");
        fflush(fplog);
    }

    ndiff_tot = local_count - dd->nbonded_global;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        nral      = NRAL(ftype);
        cl[ftype] = err_top_local->idef.il[ftype].nr/(1+nral);
    }

    gmx_sumi(F_NRE, cl, cr);

    if (DDMASTER(dd))
    {
        if (fplog)
        {
            fprintf(fplog, "\nA list of missing interactions:\n");
        }
        fprintf(stderr, "\nA list of missing interactions:\n");
        rest_global = dd->nbonded_global;
        rest_local  = local_count;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            /* In the reverse and local top all constraints are merged
             * into F_CONSTR. So in the if statement we skip F_CONSTRNC
             * and add these constraints when doing F_CONSTR.
             */
            if (((interaction_function[ftype].flags & IF_BOND) &&
                 (dd->reverse_top->bBCheck
                  || !(interaction_function[ftype].flags & IF_LIMZERO)))
                || (dd->reverse_top->bConstr && ftype == F_CONSTR)
                || (dd->reverse_top->bSettle && ftype == F_SETTLE))
            {
                n    = gmx_mtop_ftype_count(err_top_global, ftype);
                if (ftype == F_CONSTR)
                {
                    n += gmx_mtop_ftype_count(err_top_global, F_CONSTRNC);
                }
                ndiff = cl[ftype] - n;
                if (ndiff != 0)
                {
                    sprintf(buf, "%20s of %6d missing %6d",
                            interaction_function[ftype].longname, n, -ndiff);
                    if (fplog)
                    {
                        fprintf(fplog, "%s\n", buf);
                    }
                    fprintf(stderr, "%s\n", buf);
                }
                rest_global -= n;
                rest_local  -= cl[ftype];
            }
        }

        ndiff = rest_local - rest_global;
        if (ndiff != 0)
        {
            sprintf(buf, "%20s of %6d missing %6d", "exclusions",
                    rest_global, -ndiff);
            if (fplog)
            {
                fprintf(fplog, "%s\n", buf);
            }
            fprintf(stderr, "%s\n", buf);
        }
    }

    print_missing_interactions_atoms(fplog, cr, err_top_global,
                                     &err_top_local->idef);
    write_dd_pdb("dd_dump_err", 0, "dump", top_global, cr,
                 -1, state_local->x, state_local->box);
    if (DDMASTER(dd))
    {
        if (ndiff_tot > 0)
        {
            gmx_incons("One or more interactions were multiple assigned in the domain decompostion");
        }
        else
        {
            gmx_fatal(FARGS, "%d of the %d bonded interactions could not be calculated because some atoms involved moved further apart than the multi-body cut-off distance (%g nm) or the two-body cut-off distance (%g nm), see option -rdd, for pairs and tabulated bonds also see option -ddcheck", -ndiff_tot, cr->dd->nbonded_global, dd_cutoff_multibody(cr->dd), dd_cutoff_twobody(cr->dd));
        }
    }
}

/*! \brief Return global topology molecule information for global atom index \p i_gl */
static void global_atomnr_to_moltype_ind(gmx_reverse_top_t *rt, int i_gl,
                                         int *mb, int *mt, int *mol, int *i_mol)
{
    molblock_ind_t *mbi   = rt->mbi;
    int             start = 0;
    int             end   = rt->nmolblock; /* exclusive */
    int             mid;

    /* binary search for molblock_ind */
    while (TRUE)
    {
        mid = (start+end)/2;
        if (i_gl >= mbi[mid].a_end)
        {
            start = mid+1;
        }
        else if (i_gl < mbi[mid].a_start)
        {
            end = mid;
        }
        else
        {
            break;
        }
    }

    *mb  = mid;
    mbi += mid;

    *mt    = mbi->type;
    *mol   = (i_gl - mbi->a_start) / mbi->natoms_mol;
    *i_mol = (i_gl - mbi->a_start) - (*mol)*mbi->natoms_mol;
}

/*! \brief Count the exclusions for all atoms in \p cgs */
static void count_excls(const t_block *cgs, const t_blocka *excls,
                        int *n_excl, int *n_intercg_excl, int *n_excl_at_max)
{
    int cg, at0, at1, at, excl, atj;

    *n_excl         = 0;
    *n_intercg_excl = 0;
    *n_excl_at_max  = 0;
    for (cg = 0; cg < cgs->nr; cg++)
    {
        at0 = cgs->index[cg];
        at1 = cgs->index[cg+1];
        for (at = at0; at < at1; at++)
        {
            for (excl = excls->index[at]; excl < excls->index[at+1]; excl++)
            {
                atj = excls->a[excl];
                if (atj > at)
                {
                    (*n_excl)++;
                    if (atj < at0 || atj >= at1)
                    {
                        (*n_intercg_excl)++;
                    }
                }
            }

            *n_excl_at_max = std::max(*n_excl_at_max,
                                      excls->index[at+1] - excls->index[at]);
        }
    }
}

/*! \brief Run the reverse ilist generation and store it when \p bAssign = TRUE */
static int low_make_reverse_ilist(const t_ilist *il_mt,
                                  const t_atom *atom,
                                  int **vsite_pbc, /* should be const */
                                  int *count,
                                  gmx_bool bConstr, gmx_bool bSettle,
                                  gmx_bool bBCheck,
                                  int *r_index, int *r_il,
                                  gmx_bool bLinkToAllAtoms,
                                  gmx_bool bAssign)
{
    int            ftype, nral, i, j, nlink, link;
    const t_ilist *il;
    t_iatom       *ia;
    atom_id        a;
    int            nint;
    gmx_bool       bVSite;

    nint = 0;
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if ((interaction_function[ftype].flags & (IF_BOND | IF_VSITE)) ||
            (bConstr && (ftype == F_CONSTR || ftype == F_CONSTRNC)) ||
            (bSettle && ftype == F_SETTLE))
        {
            bVSite = (interaction_function[ftype].flags & IF_VSITE);
            nral   = NRAL(ftype);
            il     = &il_mt[ftype];
            for (i = 0; i < il->nr; i += 1+nral)
            {
                ia = il->iatoms + i;
                if (bLinkToAllAtoms)
                {
                    if (bVSite)
                    {
                        /* We don't need the virtual sites for the cg-links */
                        nlink = 0;
                    }
                    else
                    {
                        nlink = nral;
                    }
                }
                else
                {
                    /* Couple to the first atom in the interaction */
                    nlink = 1;
                }
                for (link = 0; link < nlink; link++)
                {
                    a = ia[1+link];
                    if (bAssign)
                    {
                        r_il[r_index[a]+count[a]] =
                            (ftype == F_CONSTRNC ? F_CONSTR : ftype);
                        r_il[r_index[a]+count[a]+1] = ia[0];
                        for (j = 1; j < 1+nral; j++)
                        {
                            /* Store the molecular atom number */
                            r_il[r_index[a]+count[a]+1+j] = ia[j];
                        }
                    }
                    if (interaction_function[ftype].flags & IF_VSITE)
                    {
                        if (bAssign)
                        {
                            /* Add an entry to iatoms for storing
                             * which of the constructing atoms are
                             * vsites again.
                             */
                            r_il[r_index[a]+count[a]+2+nral] = 0;
                            for (j = 2; j < 1+nral; j++)
                            {
                                if (atom[ia[j]].ptype == eptVSite)
                                {
                                    r_il[r_index[a]+count[a]+2+nral] |= (2<<j);
                                }
                            }
                            /* Store vsite pbc atom in a second extra entry */
                            r_il[r_index[a]+count[a]+2+nral+1] =
                                (vsite_pbc ? vsite_pbc[ftype-F_VSITE2][i/(1+nral)] : -2);
                        }
                    }
                    else
                    {
                        /* We do not count vsites since they are always
                         * uniquely assigned and can be assigned
                         * to multiple nodes with recursive vsites.
                         */
                        if (bBCheck ||
                            !(interaction_function[ftype].flags & IF_LIMZERO))
                        {
                            nint++;
                        }
                    }
                    count[a] += 2 + nral_rt(ftype);
                }
            }
        }
    }

    return nint;
}

/*! \brief Make the reverse ilist: a list of bonded interactions linked to atoms */
static int make_reverse_ilist(const t_ilist *ilist,
                              const t_atoms *atoms,
                              int **vsite_pbc, /* should be const (C issue) */
                              gmx_bool bConstr, gmx_bool bSettle,
                              gmx_bool bBCheck,
                              gmx_bool bLinkToAllAtoms,
                              reverse_ilist_t *ril_mt)
{
    int nat_mt, *count, i, nint_mt;

    /* Count the interactions */
    nat_mt = atoms->nr;
    snew(count, nat_mt);
    low_make_reverse_ilist(ilist, atoms->atom, vsite_pbc,
                           count,
                           bConstr, bSettle, bBCheck, NULL, NULL,
                           bLinkToAllAtoms, FALSE);

    snew(ril_mt->index, nat_mt+1);
    ril_mt->index[0] = 0;
    for (i = 0; i < nat_mt; i++)
    {
        ril_mt->index[i+1] = ril_mt->index[i] + count[i];
        count[i]           = 0;
    }
    snew(ril_mt->il, ril_mt->index[nat_mt]);

    /* Store the interactions */
    nint_mt =
        low_make_reverse_ilist(ilist, atoms->atom, vsite_pbc,
                               count,
                               bConstr, bSettle, bBCheck,
                               ril_mt->index, ril_mt->il,
                               bLinkToAllAtoms, TRUE);

    sfree(count);

    return nint_mt;
}

/*! \brief Destroys a reverse_ilist_t struct */
static void destroy_reverse_ilist(reverse_ilist_t *ril)
{
    sfree(ril->index);
    sfree(ril->il);
}

/*! \brief Generate the reverse topology */
static gmx_reverse_top_t *make_reverse_top(gmx_mtop_t *mtop, gmx_bool bFE,
                                           int ***vsite_pbc_molt,
                                           gmx_bool bConstr, gmx_bool bSettle,
                                           gmx_bool bBCheck, int *nint)
{
    int                mt, i, mb;
    gmx_reverse_top_t *rt;
    int               *nint_mt;
    gmx_moltype_t     *molt;
    int                thread;

    snew(rt, 1);

    /* Should we include constraints (for SHAKE) in rt? */
    rt->bConstr = bConstr;
    rt->bSettle = bSettle;
    rt->bBCheck = bBCheck;

    rt->bMultiCGmols = FALSE;
    snew(nint_mt, mtop->nmoltype);
    snew(rt->ril_mt, mtop->nmoltype);
    rt->ril_mt_tot_size = 0;
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        molt = &mtop->moltype[mt];
        if (molt->cgs.nr > 1)
        {
            rt->bMultiCGmols = TRUE;
        }

        /* Make the atom to interaction list for this molecule type */
        nint_mt[mt] =
            make_reverse_ilist(molt->ilist, &molt->atoms,
                               vsite_pbc_molt ? vsite_pbc_molt[mt] : NULL,
                               rt->bConstr, rt->bSettle, rt->bBCheck, FALSE,
                               &rt->ril_mt[mt]);

        rt->ril_mt_tot_size += rt->ril_mt[mt].index[molt->atoms.nr];
    }
    if (debug)
    {
        fprintf(debug, "The total size of the atom to interaction index is %d integers\n", rt->ril_mt_tot_size);
    }

    *nint = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        *nint += mtop->molblock[mb].nmol*nint_mt[mtop->molblock[mb].type];
    }
    sfree(nint_mt);

    /* Make an intermolecular reverse top, if necessary */
    rt->bIntermolecularInteractions = mtop->bIntermolecularInteractions;
    if (rt->bIntermolecularInteractions)
    {
        t_atoms atoms_global;

        rt->ril_intermol.index = NULL;
        rt->ril_intermol.il    = NULL;

        atoms_global.nr   = mtop->natoms;
        atoms_global.atom = NULL; /* Only used with virtual sites */

        *nint +=
            make_reverse_ilist(mtop->intermolecular_ilist, &atoms_global,
                               NULL,
                               rt->bConstr, rt->bSettle, rt->bBCheck, FALSE,
                               &rt->ril_intermol);
    }

    if (bFE && gmx_mtop_bondeds_free_energy(mtop))
    {
        rt->ilsort = ilsortFE_UNSORTED;
    }
    else
    {
        rt->ilsort = ilsortNO_FE;
    }

    /* Make a molblock index for fast searching */
    snew(rt->mbi, mtop->nmolblock);
    rt->nmolblock = mtop->nmolblock;
    i             = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        rt->mbi[mb].a_start    = i;
        i += mtop->molblock[mb].nmol*mtop->molblock[mb].natoms_mol;
        rt->mbi[mb].a_end      = i;
        rt->mbi[mb].natoms_mol = mtop->molblock[mb].natoms_mol;
        rt->mbi[mb].type       = mtop->molblock[mb].type;
    }

    rt->nthread = gmx_omp_nthreads_get(emntDomdec);
    snew(rt->th_work, rt->nthread);
    if (vsite_pbc_molt != NULL)
    {
        for (thread = 0; thread < rt->nthread; thread++)
        {
            snew(rt->th_work[thread].vsite_pbc, F_VSITEN-F_VSITE2+1);
            snew(rt->th_work[thread].vsite_pbc_nalloc, F_VSITEN-F_VSITE2+1);
        }
    }

    return rt;
}

void dd_make_reverse_top(FILE *fplog,
                         gmx_domdec_t *dd, gmx_mtop_t *mtop,
                         gmx_vsite_t *vsite,
                         t_inputrec *ir, gmx_bool bBCheck)
{
    if (fplog)
    {
        fprintf(fplog, "\nLinking all bonded interactions to atoms\n");
    }

    /* If normal and/or settle constraints act only within charge groups,
     * we can store them in the reverse top and simply assign them to domains.
     * Otherwise we need to assign them to multiple domains and set up
     * the parallel version constraint algorithm(s).
     */

    dd->reverse_top = make_reverse_top(mtop, ir->efep != efepNO,
                                       vsite ? vsite->vsite_pbc_molt : NULL,
                                       !dd->bInterCGcons, !dd->bInterCGsettles,
                                       bBCheck, &dd->nbonded_global);

    gmx_reverse_top_t *rt = dd->reverse_top;

    if (rt->ril_mt_tot_size >= 200000 &&
        mtop->mols.nr > 1 &&
        mtop->nmolblock == 1 && mtop->molblock[0].nmol == 1)
    {
        /* mtop comes from a pre Gromacs 4 tpr file */
        const char *note = "NOTE: The tpr file used for this simulation is in an old format, for less memory usage and possibly more performance create a new tpr file with an up to date version of grompp";
        if (fplog)
        {
            fprintf(fplog, "\n%s\n\n", note);
        }
        if (DDMASTER(dd))
        {
            fprintf(stderr, "\n%s\n\n", note);
        }
    }

    /* With the Verlet scheme, exclusions are handled in the non-bonded
     * kernels and only exclusions inside the cut-off lead to exclusion
     * forces. Since each atom pair is treated at most once in the non-bonded
     * kernels, it doesn't matter if the exclusions for the same atom pair
     * appear multiple times in the exclusion list. In contrast, the, old,
     * group cut-off scheme loops over a list of exclusions, so there each
     * excluded pair should appear exactly once.
     */
    rt->bExclRequired = (ir->cutoff_scheme == ecutsGROUP &&
                         IR_EXCL_FORCES(*ir));

    int nexcl, mb;

    nexcl              = 0;
    dd->n_intercg_excl = 0;
    rt->n_excl_at_max  = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        gmx_molblock_t *molb;
        gmx_moltype_t  *molt;
        int             n_excl_mol, n_excl_icg, n_excl_at_max;

        molb                = &mtop->molblock[mb];
        molt                = &mtop->moltype[molb->type];
        count_excls(&molt->cgs, &molt->excls,
                    &n_excl_mol, &n_excl_icg, &n_excl_at_max);
        nexcl              += molb->nmol*n_excl_mol;
        dd->n_intercg_excl += molb->nmol*n_excl_icg;
        rt->n_excl_at_max   = std::max(rt->n_excl_at_max, n_excl_at_max);
    }
    if (rt->bExclRequired)
    {
        dd->nbonded_global += nexcl;
        if (EEL_FULL(ir->coulombtype) && dd->n_intercg_excl > 0 && fplog)
        {
            fprintf(fplog, "There are %d inter charge-group exclusions,\n"
                    "will use an extra communication step for exclusion forces for %s\n",
                    dd->n_intercg_excl, eel_names[ir->coulombtype]);
        }
    }

    if (vsite && vsite->n_intercg_vsite > 0)
    {
        if (fplog)
        {
            fprintf(fplog, "There are %d inter charge-group virtual sites,\n"
                    "will an extra communication step for selected coordinates and forces\n",
                    vsite->n_intercg_vsite);
        }
        init_domdec_vsites(dd, vsite->n_intercg_vsite);
    }

    if (dd->bInterCGcons || dd->bInterCGsettles)
    {
        init_domdec_constraints(dd, mtop);
    }
    if (fplog)
    {
        fprintf(fplog, "\n");
    }
}

/*! \brief Store a vsite interaction at the end of \p il
 *
 * This routine is very similar to add_ifunc, but vsites interactions
 * have more work to do than other kinds of interactions, and the
 * complex way nral (and thus vector contents) depends on ftype
 * confuses static analysis tools unless we fuse the vsite
 * atom-indexing organization code with the ifunc-adding code, so that
 * they can see that nral is the same value. */
static gmx_inline void
add_ifunc_for_vsites(t_iatom *tiatoms, gmx_ga2la_t ga2la,
                     int nral, gmx_bool bHomeA,
                     int a, int a_gl, int a_mol,
                     const t_iatom *iatoms,
                     t_ilist *il)
{
    t_iatom *liatoms;

    if (il->nr+1+nral > il->nalloc)
    {
        il->nalloc = over_alloc_large(il->nr+1+nral);
        srenew(il->iatoms, il->nalloc);
    }
    liatoms = il->iatoms + il->nr;
    il->nr += 1 + nral;

    /* Copy the type */
    tiatoms[0] = iatoms[0];

    if (bHomeA)
    {
        /* We know the local index of the first atom */
        tiatoms[1] = a;
    }
    else
    {
        /* Convert later in make_local_vsites */
        tiatoms[1] = -a_gl - 1;
    }

    for (int k = 2; k < 1+nral; k++)
    {
        int ak_gl = a_gl + iatoms[k] - a_mol;
        if (!ga2la_get_home(ga2la, ak_gl, &tiatoms[k]))
        {
            /* Copy the global index, convert later in make_local_vsites */
            tiatoms[k] = -(ak_gl + 1);
        }
        // Note that ga2la_get_home always sets the third parameter if
        // it returns TRUE
    }
    for (int k = 0; k < 1+nral; k++)
    {
        liatoms[k] = tiatoms[k];
    }
}

/*! \brief Store a bonded interaction at the end of \p il */
static gmx_inline void add_ifunc(int nral, t_iatom *tiatoms, t_ilist *il)
{
    t_iatom *liatoms;
    int      k;

    if (il->nr+1+nral > il->nalloc)
    {
        il->nalloc = over_alloc_large(il->nr+1+nral);
        srenew(il->iatoms, il->nalloc);
    }
    liatoms = il->iatoms + il->nr;
    for (k = 0; k <= nral; k++)
    {
        liatoms[k] = tiatoms[k];
    }
    il->nr += 1 + nral;
}

/*! \brief Store a position restraint in idef and iatoms, complex because the parameters are different for each entry */
static void add_posres(int mol, int a_mol, const gmx_molblock_t *molb,
                       t_iatom *iatoms, const t_iparams *ip_in,
                       t_idef *idef)
{
    int        n, a_molb;
    t_iparams *ip;

    /* This position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    n = idef->il[F_POSRES].nr/2;
    if (n+1 > idef->iparams_posres_nalloc)
    {
        idef->iparams_posres_nalloc = over_alloc_dd(n+1);
        srenew(idef->iparams_posres, idef->iparams_posres_nalloc);
    }
    ip = &idef->iparams_posres[n];
    /* Copy the force constants */
    *ip = ip_in[iatoms[0]];

    /* Get the position restraint coordinates from the molblock */
    a_molb = mol*molb->natoms_mol + a_mol;
    if (a_molb >= molb->nposres_xA)
    {
        gmx_incons("Not enough position restraint coordinates");
    }
    ip->posres.pos0A[XX] = molb->posres_xA[a_molb][XX];
    ip->posres.pos0A[YY] = molb->posres_xA[a_molb][YY];
    ip->posres.pos0A[ZZ] = molb->posres_xA[a_molb][ZZ];
    if (molb->nposres_xB > 0)
    {
        ip->posres.pos0B[XX] = molb->posres_xB[a_molb][XX];
        ip->posres.pos0B[YY] = molb->posres_xB[a_molb][YY];
        ip->posres.pos0B[ZZ] = molb->posres_xB[a_molb][ZZ];
    }
    else
    {
        ip->posres.pos0B[XX] = ip->posres.pos0A[XX];
        ip->posres.pos0B[YY] = ip->posres.pos0A[YY];
        ip->posres.pos0B[ZZ] = ip->posres.pos0A[ZZ];
    }
    /* Set the parameter index for idef->iparams_posre */
    iatoms[0] = n;
}

/*! \brief Store a flat-bottomed position restraint in idef and iatoms, complex because the parameters are different for each entry */
static void add_fbposres(int mol, int a_mol, const gmx_molblock_t *molb,
                         t_iatom *iatoms, const t_iparams *ip_in,
                         t_idef *idef)
{
    int        n, a_molb;
    t_iparams *ip;

    /* This flat-bottom position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    n = idef->il[F_FBPOSRES].nr/2;
    if (n+1 > idef->iparams_fbposres_nalloc)
    {
        idef->iparams_fbposres_nalloc = over_alloc_dd(n+1);
        srenew(idef->iparams_fbposres, idef->iparams_fbposres_nalloc);
    }
    ip = &idef->iparams_fbposres[n];
    /* Copy the force constants */
    *ip = ip_in[iatoms[0]];

    /* Get the position restriant coordinats from the molblock */
    a_molb = mol*molb->natoms_mol + a_mol;
    if (a_molb >= molb->nposres_xA)
    {
        gmx_incons("Not enough position restraint coordinates");
    }
    /* Take reference positions from A position of normal posres */
    ip->fbposres.pos0[XX] = molb->posres_xA[a_molb][XX];
    ip->fbposres.pos0[YY] = molb->posres_xA[a_molb][YY];
    ip->fbposres.pos0[ZZ] = molb->posres_xA[a_molb][ZZ];

    /* Note: no B-type for flat-bottom posres */

    /* Set the parameter index for idef->iparams_posre */
    iatoms[0] = n;
}

/*! \brief Store a virtual site interaction, complex because of PBC and recursion */
static void add_vsite(gmx_ga2la_t ga2la, const int *index, const int *rtil,
                      int ftype, int nral,
                      gmx_bool bHomeA, int a, int a_gl, int a_mol,
                      const t_iatom *iatoms,
                      t_idef *idef, int **vsite_pbc, int *vsite_pbc_nalloc)
{
    int     k, vsi, pbc_a_mol;
    t_iatom tiatoms[1+MAXATOMLIST];
    int     j, ftype_r, nral_r;

    /* Add this interaction to the local topology */
    add_ifunc_for_vsites(tiatoms, ga2la, nral, bHomeA, a, a_gl, a_mol, iatoms, &idef->il[ftype]);

    if (vsite_pbc)
    {
        vsi = idef->il[ftype].nr/(1+nral) - 1;
        if (vsi >= vsite_pbc_nalloc[ftype-F_VSITE2])
        {
            vsite_pbc_nalloc[ftype-F_VSITE2] = over_alloc_large(vsi+1);
            srenew(vsite_pbc[ftype-F_VSITE2], vsite_pbc_nalloc[ftype-F_VSITE2]);
        }
        if (bHomeA)
        {
            pbc_a_mol = iatoms[1+nral+1];
            if (pbc_a_mol < 0)
            {
                /* The pbc flag is one of the following two options:
                 * -2: vsite and all constructing atoms are within the same cg, no pbc
                 * -1: vsite and its first constructing atom are in the same cg, do pbc
                 */
                vsite_pbc[ftype-F_VSITE2][vsi] = pbc_a_mol;
            }
            else
            {
                /* Set the pbc atom for this vsite so we can make its pbc
                 * identical to the rest of the atoms in its charge group.
                 * Since the order of the atoms does not change within a charge
                 * group, we do not need the global to local atom index.
                 */
                vsite_pbc[ftype-F_VSITE2][vsi] = a + pbc_a_mol - iatoms[1];
            }
        }
        else
        {
            /* This vsite is non-home (required for recursion),
             * and therefore there is no charge group to match pbc with.
             * But we always turn on full_pbc to assure that higher order
             * recursion works correctly.
             */
            vsite_pbc[ftype-F_VSITE2][vsi] = -1;
        }
    }

    if (iatoms[1+nral])
    {
        /* Check for recursion */
        for (k = 2; k < 1+nral; k++)
        {
            if ((iatoms[1+nral] & (2<<k)) && (tiatoms[k] < 0))
            {
                /* This construction atoms is a vsite and not a home atom */
                if (gmx_debug_at)
                {
                    fprintf(debug, "Constructing atom %d of vsite atom %d is a vsite and non-home\n", iatoms[k]+1, a_mol+1);
                }
                /* Find the vsite construction */

                /* Check all interactions assigned to this atom */
                j = index[iatoms[k]];
                while (j < index[iatoms[k]+1])
                {
                    ftype_r = rtil[j++];
                    nral_r  = NRAL(ftype_r);
                    if (interaction_function[ftype_r].flags & IF_VSITE)
                    {
                        /* Add this vsite (recursion) */
                        add_vsite(ga2la, index, rtil, ftype_r, nral_r,
                                  FALSE, -1, a_gl+iatoms[k]-iatoms[1], iatoms[k],
                                  rtil+j, idef, vsite_pbc, vsite_pbc_nalloc);
                        j += 1 + nral_r + 2;
                    }
                    else
                    {
                        j += 1 + nral_r;
                    }
                }
            }
        }
    }
}

/*! \brief Update the local atom to local charge group index */
static void make_la2lc(gmx_domdec_t *dd)
{
    int *cgindex, *la2lc, cg, a;

    cgindex = dd->cgindex;

    if (dd->nat_tot > dd->la2lc_nalloc)
    {
        dd->la2lc_nalloc = over_alloc_dd(dd->nat_tot);
        snew(dd->la2lc, dd->la2lc_nalloc);
    }
    la2lc = dd->la2lc;

    /* Make the local atom to local cg index */
    for (cg = 0; cg < dd->ncg_tot; cg++)
    {
        for (a = cgindex[cg]; a < cgindex[cg+1]; a++)
        {
            la2lc[a] = cg;
        }
    }
}

/*! \brief Returns the squared distance between charge groups \p i and \p j */
static real dd_dist2(t_pbc *pbc_null, rvec *cg_cm, const int *la2lc, int i, int j)
{
    rvec dx;

    if (pbc_null)
    {
        pbc_dx_aiuc(pbc_null, cg_cm[la2lc[i]], cg_cm[la2lc[j]], dx);
    }
    else
    {
        rvec_sub(cg_cm[la2lc[i]], cg_cm[la2lc[j]], dx);
    }

    return norm2(dx);
}

/*! \brief Append t_blocka block structures 1 to nsrc in src to *dest */
static void combine_blocka(t_blocka *dest, const thread_work_t *src, int nsrc)
{
    int ni, na, s, i;

    ni = src[nsrc-1].excl.nr;
    na = 0;
    for (s = 0; s < nsrc; s++)
    {
        na += src[s].excl.nra;
    }
    if (ni + 1 > dest->nalloc_index)
    {
        dest->nalloc_index = over_alloc_large(ni+1);
        srenew(dest->index, dest->nalloc_index);
    }
    if (dest->nra + na > dest->nalloc_a)
    {
        dest->nalloc_a = over_alloc_large(dest->nra+na);
        srenew(dest->a, dest->nalloc_a);
    }
    for (s = 1; s < nsrc; s++)
    {
        for (i = dest->nr+1; i < src[s].excl.nr+1; i++)
        {
            dest->index[i] = dest->nra + src[s].excl.index[i];
        }
        for (i = 0; i < src[s].excl.nra; i++)
        {
            dest->a[dest->nra+i] = src[s].excl.a[i];
        }
        dest->nr   = src[s].excl.nr;
        dest->nra += src[s].excl.nra;
    }
}

/*! \brief Append t_idef structures 1 to nsrc in src to *dest,
 * virtual sites need special attention, as pbc info differs per vsite.
 */
static void combine_idef(t_idef *dest, const thread_work_t *src, int nsrc,
                         gmx_vsite_t *vsite)
{
    int ftype;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        int n, s;

        n = 0;
        for (s = 1; s < nsrc; s++)
        {
            n += src[s].idef.il[ftype].nr;
        }
        if (n > 0)
        {
            t_ilist *ild;

            ild = &dest->il[ftype];

            if (ild->nr + n > ild->nalloc)
            {
                ild->nalloc = over_alloc_large(ild->nr+n);
                srenew(ild->iatoms, ild->nalloc);
            }

            gmx_bool vpbc;
            int      nral1 = 0, ftv = 0;

            vpbc = ((interaction_function[ftype].flags & IF_VSITE) &&
                    vsite->vsite_pbc_loc != NULL);
            if (vpbc)
            {
                nral1 = 1 + NRAL(ftype);
                ftv   = ftype - F_VSITE2;
                if ((ild->nr + n)/nral1 > vsite->vsite_pbc_loc_nalloc[ftv])
                {
                    vsite->vsite_pbc_loc_nalloc[ftv] =
                        over_alloc_large((ild->nr + n)/nral1);
                    srenew(vsite->vsite_pbc_loc[ftv],
                           vsite->vsite_pbc_loc_nalloc[ftv]);
                }
            }

            for (s = 1; s < nsrc; s++)
            {
                const t_ilist *ils;
                int            i;

                ils = &src[s].idef.il[ftype];
                for (i = 0; i < ils->nr; i++)
                {
                    ild->iatoms[ild->nr+i] = ils->iatoms[i];
                }
                if (vpbc)
                {
                    for (i = 0; i < ils->nr; i += nral1)
                    {
                        vsite->vsite_pbc_loc[ftv][(ild->nr+i)/nral1] =
                            src[s].vsite_pbc[ftv][i/nral1];
                    }
                }

                ild->nr += ils->nr;
            }
        }
    }

    /* Position restraints need an additional treatment */
    if (dest->il[F_POSRES].nr > 0)
    {
        int n, s, i;

        n = dest->il[F_POSRES].nr/2;
        if (n > dest->iparams_posres_nalloc)
        {
            dest->iparams_posres_nalloc = over_alloc_large(n);
            srenew(dest->iparams_posres, dest->iparams_posres_nalloc);
        }
        /* Set n to the number of original position restraints in dest */
        for (s = 1; s < nsrc; s++)
        {
            n -= src[s].idef.il[F_POSRES].nr/2;
        }
        for (s = 1; s < nsrc; s++)
        {
            for (i = 0; i < src[s].idef.il[F_POSRES].nr/2; i++)
            {
                /* Correct the index into iparams_posres */
                dest->il[F_POSRES].iatoms[n*2] = n;
                /* Copy the position restraint force parameters */
                dest->iparams_posres[n] = src[s].idef.iparams_posres[i];
                n++;
            }
        }
    }
}

/*! \brief Check and when available assign bonded interactions for local atom i
 */
static gmx_inline void
check_assign_interactions_atom(int i, int i_gl,
                               int mol, int i_mol,
                               const int *index, const int *rtil,
                               gmx_bool bInterMolInteractions,
                               int ind_start, int ind_end,
                               const gmx_domdec_t *dd,
                               const gmx_domdec_zones_t *zones,
                               const gmx_molblock_t *molb,
                               gmx_bool bRCheckMB, ivec rcheck, gmx_bool bRCheck2B,
                               real rc2,
                               int *la2lc,
                               t_pbc *pbc_null,
                               rvec *cg_cm,
                               const t_iparams *ip_in,
                               t_idef *idef,
                               int **vsite_pbc, int *vsite_pbc_nalloc,
                               int iz,
                               gmx_bool bBCheck,
                               int *nbonded_local)
{
    int j;

    j = ind_start;
    while (j < ind_end)
    {
        int            ftype;
        const t_iatom *iatoms;
        int            nral;
        t_iatom        tiatoms[1 + MAXATOMLIST];

        ftype  = rtil[j++];
        iatoms = rtil + j;
        nral   = NRAL(ftype);
        if (ftype == F_SETTLE)
        {
            /* Settles are only in the reverse top when they
             * operate within a charge group. So we can assign
             * them without checks. We do this only for performance
             * reasons; it could be handled by the code below.
             */
            if (iz == 0)
            {
                /* Home zone: add this settle to the local topology */
                tiatoms[0] = iatoms[0];
                tiatoms[1] = i;
                tiatoms[2] = i + iatoms[2] - iatoms[1];
                tiatoms[3] = i + iatoms[3] - iatoms[1];
                add_ifunc(nral, tiatoms, &idef->il[ftype]);
                (*nbonded_local)++;
            }
            j += 1 + nral;
        }
        else if (interaction_function[ftype].flags & IF_VSITE)
        {
            assert(!bInterMolInteractions);
            /* The vsite construction goes where the vsite itself is */
            if (iz == 0)
            {
                add_vsite(dd->ga2la, index, rtil, ftype, nral,
                          TRUE, i, i_gl, i_mol,
                          iatoms, idef, vsite_pbc, vsite_pbc_nalloc);
            }
            j += 1 + nral + 2;
        }
        else
        {
            gmx_bool bUse;

            /* Copy the type */
            tiatoms[0] = iatoms[0];

            if (nral == 1)
            {
                assert(!bInterMolInteractions);
                /* Assign single-body interactions to the home zone */
                if (iz == 0)
                {
                    bUse       = TRUE;
                    tiatoms[1] = i;
                    if (ftype == F_POSRES)
                    {
                        add_posres(mol, i_mol, molb, tiatoms, ip_in,
                                   idef);
                    }
                    else if (ftype == F_FBPOSRES)
                    {
                        add_fbposres(mol, i_mol, molb, tiatoms, ip_in,
                                     idef);
                    }
                }
                else
                {
                    bUse = FALSE;
                }
            }
            else if (nral == 2)
            {
                /* This is a two-body interaction, we can assign
                 * analogous to the non-bonded assignments.
                 */
                int k_gl, a_loc, kz;

                if (!bInterMolInteractions)
                {
                    /* Get the global index using the offset in the molecule */
                    k_gl = i_gl + iatoms[2] - i_mol;
                }
                else
                {
                    k_gl = iatoms[2];
                }
                if (!ga2la_get(dd->ga2la, k_gl, &a_loc, &kz))
                {
                    bUse = FALSE;
                }
                else
                {
                    if (kz >= zones->n)
                    {
                        kz -= zones->n;
                    }
                    /* Check zone interaction assignments */
                    bUse = ((iz < zones->nizone &&
                             iz <= kz &&
                             kz >= zones->izone[iz].j0 &&
                             kz <  zones->izone[iz].j1) ||
                            (kz < zones->nizone &&
                                  iz > kz &&
                             iz >= zones->izone[kz].j0 &&
                             iz <  zones->izone[kz].j1));
                    if (bUse)
                    {
                        tiatoms[1] = i;
                        tiatoms[2] = a_loc;
                        /* If necessary check the cgcm distance */
                        if (bRCheck2B &&
                            dd_dist2(pbc_null, cg_cm, la2lc,
                                     tiatoms[1], tiatoms[2]) >= rc2)
                        {
                            bUse = FALSE;
                        }
                    }
                }
            }
            else
            {
                /* Assign this multi-body bonded interaction to
                 * the local node if we have all the atoms involved
                 * (local or communicated) and the minimum zone shift
                 * in each dimension is zero, for dimensions
                 * with 2 DD cells an extra check may be necessary.
                 */
                ivec k_zero, k_plus;
                int  k;

                bUse = TRUE;
                clear_ivec(k_zero);
                clear_ivec(k_plus);
                for (k = 1; k <= nral && bUse; k++)
                {
                    gmx_bool bLocal;
                    int      k_gl, a_loc;
                    int      kz;

                    if (!bInterMolInteractions)
                    {
                        /* Get the global index using the offset in the molecule */
                        k_gl = i_gl + iatoms[k] - i_mol;
                    }
                    else
                    {
                        k_gl = iatoms[k];
                    }
                    bLocal = ga2la_get(dd->ga2la, k_gl, &a_loc, &kz);
                    if (!bLocal || kz >= zones->n)
                    {
                        /* We do not have this atom of this interaction
                         * locally, or it comes from more than one cell
                         * away.
                         */
                        bUse = FALSE;
                    }
                    else
                    {
                        int d;

                        tiatoms[k] = a_loc;
                        for (d = 0; d < DIM; d++)
                        {
                            if (zones->shift[kz][d] == 0)
                            {
                                k_zero[d] = k;
                            }
                            else
                            {
                                k_plus[d] = k;
                            }
                        }
                    }
                }
                bUse = (bUse &&
                        k_zero[XX] && k_zero[YY] && k_zero[ZZ]);
                if (bRCheckMB)
                {
                    int d;

                    for (d = 0; (d < DIM && bUse); d++)
                    {
                        /* Check if the cg_cm distance falls within
                         * the cut-off to avoid possible multiple
                         * assignments of bonded interactions.
                         */
                        if (rcheck[d] &&
                            k_plus[d] &&
                            dd_dist2(pbc_null, cg_cm, la2lc,
                                     tiatoms[k_zero[d]], tiatoms[k_plus[d]]) >= rc2)
                        {
                            bUse = FALSE;
                        }
                    }
                }
            }
            if (bUse)
            {
                /* Add this interaction to the local topology */
                add_ifunc(nral, tiatoms, &idef->il[ftype]);
                /* Sum so we can check in global_stat
                 * if we have everything.
                 */
                if (bBCheck ||
                    !(interaction_function[ftype].flags & IF_LIMZERO))
                {
                    (*nbonded_local)++;
                }
            }
            j += 1 + nral;
        }
    }
}

/*! \brief This function looks up and assigns bonded interactions for zone iz.
 *
 * With thread parallelizing each thread acts on a different atom range:
 * at_start to at_end.
 */
static int make_bondeds_zone(gmx_domdec_t *dd,
                             const gmx_domdec_zones_t *zones,
                             const gmx_molblock_t *molb,
                             gmx_bool bRCheckMB, ivec rcheck, gmx_bool bRCheck2B,
                             real rc2,
                             int *la2lc, t_pbc *pbc_null, rvec *cg_cm,
                             const t_iparams *ip_in,
                             t_idef *idef,
                             int **vsite_pbc,
                             int *vsite_pbc_nalloc,
                             int izone,
                             int at_start, int at_end)
{
    int                i, i_gl, mb, mt, mol, i_mol;
    int               *index, *rtil;
    gmx_bool           bBCheck;
    gmx_reverse_top_t *rt;
    int                nbonded_local;

    rt = dd->reverse_top;

    bBCheck = rt->bBCheck;

    nbonded_local = 0;

    for (i = at_start; i < at_end; i++)
    {
        /* Get the global atom number */
        i_gl = dd->gatindex[i];
        global_atomnr_to_moltype_ind(rt, i_gl, &mb, &mt, &mol, &i_mol);
        /* Check all intramolecular interactions assigned to this atom */
        index = rt->ril_mt[mt].index;
        rtil  = rt->ril_mt[mt].il;

        check_assign_interactions_atom(i, i_gl, mol, i_mol,
                                       index, rtil, FALSE,
                                       index[i_mol], index[i_mol+1],
                                       dd, zones,
                                       &molb[mb],
                                       bRCheckMB, rcheck, bRCheck2B, rc2,
                                       la2lc,
                                       pbc_null,
                                       cg_cm,
                                       ip_in,
                                       idef, vsite_pbc, vsite_pbc_nalloc,
                                       izone,
                                       bBCheck,
                                       &nbonded_local);


        if (rt->bIntermolecularInteractions)
        {
            /* Check all intermolecular interactions assigned to this atom */
            index = rt->ril_intermol.index;
            rtil  = rt->ril_intermol.il;

            check_assign_interactions_atom(i, i_gl, mol, i_mol,
                                           index, rtil, TRUE,
                                           index[i_gl], index[i_gl + 1],
                                           dd, zones,
                                           &molb[mb],
                                           bRCheckMB, rcheck, bRCheck2B, rc2,
                                           la2lc,
                                           pbc_null,
                                           cg_cm,
                                           ip_in,
                                           idef, vsite_pbc, vsite_pbc_nalloc,
                                           izone,
                                           bBCheck,
                                           &nbonded_local);
        }
    }

    return nbonded_local;
}

/*! \brief Set the exclusion data for i-zone \p iz for the case of no exclusions */
static void set_no_exclusions_zone(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                                   int iz, t_blocka *lexcls)
{
    int  a0, a1, a;

    a0 = dd->cgindex[zones->cg_range[iz]];
    a1 = dd->cgindex[zones->cg_range[iz+1]];

    for (a = a0+1; a < a1+1; a++)
    {
        lexcls->index[a] = lexcls->nra;
    }
}

/*! \brief Set the exclusion data for i-zone \p iz
 *
 * This is a legacy version for the group scheme of the same routine below.
 * Here charge groups and distance checks to ensure unique exclusions
 * are supported.
 */
static int make_exclusions_zone_cg(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                                   const gmx_moltype_t *moltype,
                                   gmx_bool bRCheck, real rc2,
                                   int *la2lc, t_pbc *pbc_null, rvec *cg_cm,
                                   const int *cginfo,
                                   t_blocka *lexcls,
                                   int iz,
                                   int cg_start, int cg_end)
{
    int             n_excl_at_max, n, count, jla0, jla1, jla;
    int             cg, la0, la1, la, a_gl, mb, mt, mol, a_mol, j, aj_mol;
    const t_blocka *excls;
    gmx_ga2la_t     ga2la;
    int             cell;

    ga2la = dd->ga2la;

    jla0 = dd->cgindex[zones->izone[iz].jcg0];
    jla1 = dd->cgindex[zones->izone[iz].jcg1];

    n_excl_at_max = dd->reverse_top->n_excl_at_max;

    /* We set the end index, but note that we might not start at zero here */
    lexcls->nr = dd->cgindex[cg_end];

    n     = lexcls->nra;
    count = 0;
    for (cg = cg_start; cg < cg_end; cg++)
    {
        if (n + (cg_end - cg_start)*n_excl_at_max > lexcls->nalloc_a)
        {
            lexcls->nalloc_a = over_alloc_large(n + (cg_end - cg_start)*n_excl_at_max);
            srenew(lexcls->a, lexcls->nalloc_a);
        }
        la0 = dd->cgindex[cg];
        la1 = dd->cgindex[cg+1];
        if (GET_CGINFO_EXCL_INTER(cginfo[cg]) ||
            !GET_CGINFO_EXCL_INTRA(cginfo[cg]))
        {
            /* Copy the exclusions from the global top */
            for (la = la0; la < la1; la++)
            {
                lexcls->index[la] = n;
                a_gl              = dd->gatindex[la];
                global_atomnr_to_moltype_ind(dd->reverse_top, a_gl, &mb, &mt, &mol, &a_mol);
                excls = &moltype[mt].excls;
                for (j = excls->index[a_mol]; j < excls->index[a_mol+1]; j++)
                {
                    aj_mol = excls->a[j];
                    /* This computation of jla is only correct intra-cg */
                    jla = la + aj_mol - a_mol;
                    if (jla >= la0 && jla < la1)
                    {
                        /* This is an intra-cg exclusion. We can skip
                         *  the global indexing and distance checking.
                         */
                        /* Intra-cg exclusions are only required
                         * for the home zone.
                         */
                        if (iz == 0)
                        {
                            lexcls->a[n++] = jla;
                            /* Check to avoid double counts */
                            if (jla > la)
                            {
                                count++;
                            }
                        }
                    }
                    else
                    {
                        /* This is a inter-cg exclusion */
                        /* Since exclusions are pair interactions,
                         * just like non-bonded interactions,
                         * they can be assigned properly up
                         * to the DD cutoff (not cutoff_min as
                         * for the other bonded interactions).
                         */
                        if (ga2la_get(ga2la, a_gl+aj_mol-a_mol, &jla, &cell))
                        {
                            if (iz == 0 && cell == 0)
                            {
                                lexcls->a[n++] = jla;
                                /* Check to avoid double counts */
                                if (jla > la)
                                {
                                    count++;
                                }
                            }
                            else if (jla >= jla0 && jla < jla1 &&
                                     (!bRCheck ||
                                      dd_dist2(pbc_null, cg_cm, la2lc, la, jla) < rc2))
                            {
                                /* jla > la, since jla0 > la */
                                lexcls->a[n++] = jla;
                                count++;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            /* There are no inter-cg excls and this cg is self-excluded.
             * These exclusions are only required for zone 0,
             * since other zones do not see themselves.
             */
            if (iz == 0)
            {
                for (la = la0; la < la1; la++)
                {
                    lexcls->index[la] = n;
                    for (j = la0; j < la1; j++)
                    {
                        lexcls->a[n++] = j;
                    }
                }
                count += ((la1 - la0)*(la1 - la0 - 1))/2;
            }
            else
            {
                /* We don't need exclusions for this cg */
                for (la = la0; la < la1; la++)
                {
                    lexcls->index[la] = n;
                }
            }
        }
    }

    lexcls->index[lexcls->nr] = n;
    lexcls->nra               = n;

    return count;
}

/*! \brief Set the exclusion data for i-zone \p iz */
static void make_exclusions_zone(gmx_domdec_t *dd,
                                 gmx_domdec_zones_t *zones,
                                 const gmx_moltype_t *moltype,
                                 const int *cginfo,
                                 t_blocka *lexcls,
                                 int iz,
                                 int at_start, int at_end)
{
    gmx_ga2la_t ga2la;
    int         jla0, jla1;
    int         n_excl_at_max, n, at;

    ga2la = dd->ga2la;

    jla0 = dd->cgindex[zones->izone[iz].jcg0];
    jla1 = dd->cgindex[zones->izone[iz].jcg1];

    n_excl_at_max = dd->reverse_top->n_excl_at_max;

    /* We set the end index, but note that we might not start at zero here */
    lexcls->nr = at_end;

    n = lexcls->nra;
    for (at = at_start; at < at_end; at++)
    {
        if (n + 1000 > lexcls->nalloc_a)
        {
            lexcls->nalloc_a = over_alloc_large(n + 1000);
            srenew(lexcls->a, lexcls->nalloc_a);
        }
        if (GET_CGINFO_EXCL_INTER(cginfo[at]))
        {
            int             a_gl, mb, mt, mol, a_mol, j;
            const t_blocka *excls;

            if (n + n_excl_at_max > lexcls->nalloc_a)
            {
                lexcls->nalloc_a = over_alloc_large(n + n_excl_at_max);
                srenew(lexcls->a, lexcls->nalloc_a);
            }

            /* Copy the exclusions from the global top */
            lexcls->index[at] = n;
            a_gl              = dd->gatindex[at];
            global_atomnr_to_moltype_ind(dd->reverse_top, a_gl,
                                         &mb, &mt, &mol, &a_mol);
            excls = &moltype[mt].excls;
            for (j = excls->index[a_mol]; j < excls->index[a_mol + 1]; j++)
            {
                int aj_mol, at_j, cell;

                aj_mol = excls->a[j];

                if (ga2la_get(ga2la, a_gl + aj_mol - a_mol, &at_j, &cell))
                {
                    /* This check is not necessary, but it can reduce
                     * the number of exclusions in the list, which in turn
                     * can speed up the pair list construction a bit.
                     */
                    if (at_j >= jla0 && at_j < jla1)
                    {
                        lexcls->a[n++] = at_j;
                    }
                }
            }
        }
        else
        {
            /* We don't need exclusions for this atom */
            lexcls->index[at] = n;
        }
    }

    lexcls->index[lexcls->nr] = n;
    lexcls->nra               = n;
}


/*! \brief Ensure we have enough space in \p ba for \p nindex_max indices */
static void check_alloc_index(t_blocka *ba, int nindex_max)
{
    if (nindex_max+1 > ba->nalloc_index)
    {
        ba->nalloc_index = over_alloc_dd(nindex_max+1);
        srenew(ba->index, ba->nalloc_index);
    }
}

/*! \brief Ensure that we have enough space for exclusion storate in \p lexcls */
static void check_exclusions_alloc(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                                   t_blocka *lexcls)
{
    int nr;
    int thread;

    nr = dd->cgindex[zones->izone[zones->nizone-1].cg1];

    check_alloc_index(lexcls, nr);

    for (thread = 1; thread < dd->reverse_top->nthread; thread++)
    {
        check_alloc_index(&dd->reverse_top->th_work[thread].excl, nr);
    }
}

/*! \brief Set the total count indexes for the local exclusions, needed by several functions */
static void finish_local_exclusions(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                                    t_blocka *lexcls)
{
    int la0, la;

    lexcls->nr = dd->cgindex[zones->izone[zones->nizone-1].cg1];

    if (dd->n_intercg_excl == 0)
    {
        /* There are no exclusions involving non-home charge groups,
         * but we need to set the indices for neighborsearching.
         */
        la0 = dd->cgindex[zones->izone[0].cg1];
        for (la = la0; la < lexcls->nr; la++)
        {
            lexcls->index[la] = lexcls->nra;
        }

        /* nr is only used to loop over the exclusions for Ewald and RF,
         * so we can set it to the number of home atoms for efficiency.
         */
        lexcls->nr = dd->cgindex[zones->izone[0].cg1];
    }
}

/*! \brief Clear a t_idef struct */
static void clear_idef(t_idef *idef)
{
    int  ftype;

    /* Clear the counts */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        idef->il[ftype].nr = 0;
    }
}

/*! \brief Generate and store all required local bonded interactions in \p idef and local exclusions in \p lexcls */
static int make_local_bondeds_excls(gmx_domdec_t *dd,
                                    gmx_domdec_zones_t *zones,
                                    const gmx_mtop_t *mtop,
                                    const int *cginfo,
                                    gmx_bool bRCheckMB, ivec rcheck, gmx_bool bRCheck2B,
                                    real rc,
                                    int *la2lc, t_pbc *pbc_null, rvec *cg_cm,
                                    t_idef *idef, gmx_vsite_t *vsite,
                                    t_blocka *lexcls, int *excl_count)
{
    int                nzone_bondeds, nzone_excl;
    int                izone, cg0, cg1;
    real               rc2;
    int                nbonded_local;
    int                thread;
    gmx_reverse_top_t *rt;

    if (dd->reverse_top->bMultiCGmols)
    {
        nzone_bondeds = zones->n;
    }
    else
    {
        /* Only single charge group molecules, so interactions don't
         * cross zone boundaries and we only need to assign in the home zone.
         */
        nzone_bondeds = 1;
    }

    if (dd->n_intercg_excl > 0)
    {
        /* We only use exclusions from i-zones to i- and j-zones */
        nzone_excl = zones->nizone;
    }
    else
    {
        /* There are no inter-cg exclusions and only zone 0 sees itself */
        nzone_excl = 1;
    }

    check_exclusions_alloc(dd, zones, lexcls);

    rt = dd->reverse_top;

    rc2 = rc*rc;

    /* Clear the counts */
    clear_idef(idef);
    nbonded_local = 0;

    lexcls->nr    = 0;
    lexcls->nra   = 0;
    *excl_count   = 0;

    for (izone = 0; izone < nzone_bondeds; izone++)
    {
        cg0 = zones->cg_range[izone];
        cg1 = zones->cg_range[izone + 1];

#pragma omp parallel for num_threads(rt->nthread) schedule(static)
        for (thread = 0; thread < rt->nthread; thread++)
        {
            int       cg0t, cg1t;
            t_idef   *idef_t;
            int     **vsite_pbc;
            int      *vsite_pbc_nalloc;
            t_blocka *excl_t;

            cg0t = cg0 + ((cg1 - cg0)* thread   )/rt->nthread;
            cg1t = cg0 + ((cg1 - cg0)*(thread+1))/rt->nthread;

            if (thread == 0)
            {
                idef_t = idef;
            }
            else
            {
                idef_t = &rt->th_work[thread].idef;
                clear_idef(idef_t);
            }

            if (vsite && vsite->bHaveChargeGroups && vsite->n_intercg_vsite > 0)
            {
                if (thread == 0)
                {
                    vsite_pbc        = vsite->vsite_pbc_loc;
                    vsite_pbc_nalloc = vsite->vsite_pbc_loc_nalloc;
                }
                else
                {
                    vsite_pbc        = rt->th_work[thread].vsite_pbc;
                    vsite_pbc_nalloc = rt->th_work[thread].vsite_pbc_nalloc;
                }
            }
            else
            {
                vsite_pbc        = NULL;
                vsite_pbc_nalloc = NULL;
            }

            rt->th_work[thread].nbonded =
                make_bondeds_zone(dd, zones,
                                  mtop->molblock,
                                  bRCheckMB, rcheck, bRCheck2B, rc2,
                                  la2lc, pbc_null, cg_cm, idef->iparams,
                                  idef_t,
                                  vsite_pbc, vsite_pbc_nalloc,
                                  izone,
                                  dd->cgindex[cg0t], dd->cgindex[cg1t]);

            if (izone < nzone_excl)
            {
                if (thread == 0)
                {
                    excl_t = lexcls;
                }
                else
                {
                    excl_t      = &rt->th_work[thread].excl;
                    excl_t->nr  = 0;
                    excl_t->nra = 0;
                }

                if (dd->cgindex[dd->ncg_tot] == dd->ncg_tot &&
                    !rt->bExclRequired)
                {
                    /* No charge groups and no distance check required */
                    make_exclusions_zone(dd, zones,
                                         mtop->moltype, cginfo,
                                         excl_t,
                                         izone,
                                         cg0t, cg1t);
                }
                else
                {
                    rt->th_work[thread].excl_count =
                        make_exclusions_zone_cg(dd, zones,
                                                mtop->moltype, bRCheck2B, rc2,
                                                la2lc, pbc_null, cg_cm, cginfo,
                                                excl_t,
                                                izone,
                                                cg0t, cg1t);
                }
            }
        }

        if (rt->nthread > 1)
        {
            combine_idef(idef, rt->th_work, rt->nthread, vsite);
        }

        for (thread = 0; thread < rt->nthread; thread++)
        {
            nbonded_local += rt->th_work[thread].nbonded;
        }

        if (izone < nzone_excl)
        {
            if (rt->nthread > 1)
            {
                combine_blocka(lexcls, rt->th_work, rt->nthread);
            }

            for (thread = 0; thread < rt->nthread; thread++)
            {
                *excl_count += rt->th_work[thread].excl_count;
            }
        }
    }

    /* Some zones might not have exclusions, but some code still needs to
     * loop over the index, so we set the indices here.
     */
    for (izone = nzone_excl; izone < zones->nizone; izone++)
    {
        set_no_exclusions_zone(dd, zones, izone, lexcls);
    }

    finish_local_exclusions(dd, zones, lexcls);
    if (debug)
    {
        fprintf(debug, "We have %d exclusions, check count %d\n",
                lexcls->nra, *excl_count);
    }

    return nbonded_local;
}

void dd_make_local_cgs(gmx_domdec_t *dd, t_block *lcgs)
{
    lcgs->nr    = dd->ncg_tot;
    lcgs->index = dd->cgindex;
}

void dd_make_local_top(gmx_domdec_t *dd, gmx_domdec_zones_t *zones,
                       int npbcdim, matrix box,
                       rvec cellsize_min, ivec npulse,
                       t_forcerec *fr,
                       rvec *cgcm_or_x,
                       gmx_vsite_t *vsite,
                       gmx_mtop_t *mtop, gmx_localtop_t *ltop)
{
    gmx_bool bRCheckMB, bRCheck2B;
    real     rc = -1;
    ivec     rcheck;
    int      d, nexcl;
    t_pbc    pbc, *pbc_null = NULL;

    if (debug)
    {
        fprintf(debug, "Making local topology\n");
    }

    dd_make_local_cgs(dd, &ltop->cgs);

    bRCheckMB   = FALSE;
    bRCheck2B   = FALSE;

    if (dd->reverse_top->bMultiCGmols)
    {
        /* We need to check to which cell bondeds should be assigned */
        rc = dd_cutoff_twobody(dd);
        if (debug)
        {
            fprintf(debug, "Two-body bonded cut-off distance is %g\n", rc);
        }

        /* Should we check cg_cm distances when assigning bonded interactions? */
        for (d = 0; d < DIM; d++)
        {
            rcheck[d] = FALSE;
            /* Only need to check for dimensions where the part of the box
             * that is not communicated is smaller than the cut-off.
             */
            if (d < npbcdim && dd->nc[d] > 1 &&
                (dd->nc[d] - npulse[d])*cellsize_min[d] < 2*rc)
            {
                if (dd->nc[d] == 2)
                {
                    rcheck[d] = TRUE;
                    bRCheckMB = TRUE;
                }
                /* Check for interactions between two atoms,
                 * where we can allow interactions up to the cut-off,
                 * instead of up to the smallest cell dimension.
                 */
                bRCheck2B = TRUE;
            }
            if (debug)
            {
                fprintf(debug,
                        "dim %d cellmin %f bonded rcheck[%d] = %d, bRCheck2B = %d\n",
                        d, cellsize_min[d], d, rcheck[d], bRCheck2B);
            }
        }
        if (bRCheckMB || bRCheck2B)
        {
            make_la2lc(dd);
            if (fr->bMolPBC)
            {
                set_pbc_dd(&pbc, fr->ePBC, dd, TRUE, box);
                pbc_null = &pbc;
            }
            else
            {
                pbc_null = NULL;
            }
        }
    }

    dd->nbonded_local =
        make_local_bondeds_excls(dd, zones, mtop, fr->cginfo,
                                 bRCheckMB, rcheck, bRCheck2B, rc,
                                 dd->la2lc,
                                 pbc_null, cgcm_or_x,
                                 &ltop->idef, vsite,
                                 &ltop->excls, &nexcl);

    /* The ilist is not sorted yet,
     * we can only do this when we have the charge arrays.
     */
    ltop->idef.ilsort = ilsortUNKNOWN;

    if (dd->reverse_top->bExclRequired)
    {
        dd->nbonded_local += nexcl;

        forcerec_set_excl_load(fr, ltop);
    }

    ltop->atomtypes  = mtop->atomtypes;

    /* For an error message only */
    dd->reverse_top->err_top_global = mtop;
    dd->reverse_top->err_top_local  = ltop;
}

void dd_sort_local_top(gmx_domdec_t *dd, t_mdatoms *mdatoms,
                       gmx_localtop_t *ltop)
{
    if (dd->reverse_top->ilsort == ilsortNO_FE)
    {
        ltop->idef.ilsort = ilsortNO_FE;
    }
    else
    {
        gmx_sort_ilist_fe(&ltop->idef, mdatoms->chargeA, mdatoms->chargeB);
    }
}

gmx_localtop_t *dd_init_local_top(gmx_mtop_t *top_global)
{
    gmx_localtop_t *top;
    int             i;

    snew(top, 1);

    top->idef.ntypes    = top_global->ffparams.ntypes;
    top->idef.atnr      = top_global->ffparams.atnr;
    top->idef.functype  = top_global->ffparams.functype;
    top->idef.iparams   = top_global->ffparams.iparams;
    top->idef.fudgeQQ   = top_global->ffparams.fudgeQQ;
    top->idef.cmap_grid = top_global->ffparams.cmap_grid;

    for (i = 0; i < F_NRE; i++)
    {
        top->idef.il[i].iatoms = NULL;
        top->idef.il[i].nalloc = 0;
    }
    top->idef.ilsort   = ilsortUNKNOWN;

    return top;
}

void dd_init_local_state(gmx_domdec_t *dd,
                         t_state *state_global, t_state *state_local)
{
    int buf[NITEM_DD_INIT_LOCAL_STATE];

    if (DDMASTER(dd))
    {
        buf[0] = state_global->flags;
        buf[1] = state_global->ngtc;
        buf[2] = state_global->nnhpres;
        buf[3] = state_global->nhchainlength;
        buf[4] = state_global->dfhist.nlambda;
    }
    dd_bcast(dd, NITEM_DD_INIT_LOCAL_STATE*sizeof(int), buf);

    init_state(state_local, 0, buf[1], buf[2], buf[3], buf[4]);
    state_local->flags = buf[0];
}

/*! \brief Check if a link is stored in \p link between charge groups \p cg_gl and \p cg_gl_j and if not so, store a link */
static void check_link(t_blocka *link, int cg_gl, int cg_gl_j)
{
    int      k;
    gmx_bool bFound;

    bFound = FALSE;
    for (k = link->index[cg_gl]; k < link->index[cg_gl+1]; k++)
    {
        GMX_RELEASE_ASSERT(link->a, "Inconsistent NULL pointer while making charge-group links");
        if (link->a[k] == cg_gl_j)
        {
            bFound = TRUE;
        }
    }
    if (!bFound)
    {
        GMX_RELEASE_ASSERT(link->a || link->index[cg_gl+1]+1 > link->nalloc_a,
                           "Inconsistent allocation of link");
        /* Add this charge group link */
        if (link->index[cg_gl+1]+1 > link->nalloc_a)
        {
            link->nalloc_a = over_alloc_large(link->index[cg_gl+1]+1);
            srenew(link->a, link->nalloc_a);
        }
        link->a[link->index[cg_gl+1]] = cg_gl_j;
        link->index[cg_gl+1]++;
    }
}

/*! \brief Allocate and return an array of the charge group index for all atoms */
static int *make_at2cg(t_block *cgs)
{
    int *at2cg, cg, a;

    snew(at2cg, cgs->index[cgs->nr]);
    for (cg = 0; cg < cgs->nr; cg++)
    {
        for (a = cgs->index[cg]; a < cgs->index[cg+1]; a++)
        {
            at2cg[a] = cg;
        }
    }

    return at2cg;
}

t_blocka *make_charge_group_links(gmx_mtop_t *mtop, gmx_domdec_t *dd,
                                  cginfo_mb_t *cginfo_mb)
{
    gmx_bool            bExclRequired;
    int                 mb, cg_offset, cg, cg_gl, a, aj, i, j, ftype, nral, nlink_mol, mol, ncgi;
    gmx_molblock_t     *molb;
    gmx_moltype_t      *molt;
    t_block            *cgs;
    t_blocka           *excls;
    int                *a2c;
    reverse_ilist_t     ril, ril_intermol;
    t_blocka           *link;
    cginfo_mb_t        *cgi_mb;

    /* For each charge group make a list of other charge groups
     * in the system that a linked to it via bonded interactions
     * which are also stored in reverse_top.
     */

    bExclRequired = dd->reverse_top->bExclRequired;

    if (mtop->bIntermolecularInteractions)
    {
        t_atoms atoms;

        atoms.nr   = mtop->natoms;
        atoms.atom = NULL;

        make_reverse_ilist(mtop->intermolecular_ilist, &atoms,
                           NULL, FALSE, FALSE, FALSE, TRUE, &ril_intermol);
    }

    snew(link, 1);
    snew(link->index, ncg_mtop(mtop)+1);
    link->nalloc_a = 0;
    link->a        = NULL;

    link->index[0] = 0;
    cg_offset      = 0;
    ncgi           = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        if (molb->nmol == 0)
        {
            continue;
        }
        molt  = &mtop->moltype[molb->type];
        cgs   = &molt->cgs;
        excls = &molt->excls;
        a2c   = make_at2cg(cgs);
        /* Make a reverse ilist in which the interactions are linked
         * to all atoms, not only the first atom as in gmx_reverse_top.
         * The constraints are discarded here.
         */
        make_reverse_ilist(molt->ilist, &molt->atoms,
                           NULL, FALSE, FALSE, FALSE, TRUE, &ril);

        cgi_mb = &cginfo_mb[mb];

        for (mol = 0; mol < (mtop->bIntermolecularInteractions ? molb->nmol : 1); mol++)
        {
            for (cg = 0; cg < cgs->nr; cg++)
            {
                cg_gl                = cg_offset + cg;
                link->index[cg_gl+1] = link->index[cg_gl];
                for (a = cgs->index[cg]; a < cgs->index[cg+1]; a++)
                {
                    i = ril.index[a];
                    while (i < ril.index[a+1])
                    {
                        ftype = ril.il[i++];
                        nral  = NRAL(ftype);
                        /* Skip the ifunc index */
                        i++;
                        for (j = 0; j < nral; j++)
                        {
                            aj = ril.il[i+j];
                            if (a2c[aj] != cg)
                            {
                                check_link(link, cg_gl, cg_offset+a2c[aj]);
                            }
                        }
                        i += nral_rt(ftype);
                    }
                    if (bExclRequired)
                    {
                        /* Exclusions always go both ways */
                        for (j = excls->index[a]; j < excls->index[a+1]; j++)
                        {
                            aj = excls->a[j];
                            if (a2c[aj] != cg)
                            {
                                check_link(link, cg_gl, cg_offset+a2c[aj]);
                            }
                        }
                    }

                    if (mtop->bIntermolecularInteractions)
                    {
                        i = ril_intermol.index[a];
                        while (i < ril.index[a+1])
                        {
                            ftype = ril_intermol.il[i++];
                            nral  = NRAL(ftype);
                            /* Skip the ifunc index */
                            i++;
                            for (j = 0; j < nral; j++)
                            {
                                aj = ril_intermol.il[i+j];
                                if (a2c[aj] != cg_offset + cg)
                                {
                                    check_link(link, cg_gl, a2c[aj]);
                                }
                            }
                            i += nral_rt(ftype);
                        }
                    }
                }
                if (link->index[cg_gl+1] - link->index[cg_gl] > 0)
                {
                    SET_CGINFO_BOND_INTER(cgi_mb->cginfo[cg]);
                    ncgi++;
                }
            }

            cg_offset += cgs->nr;
        }
        nlink_mol = link->index[cg_offset] - link->index[cg_offset-cgs->nr];

        destroy_reverse_ilist(&ril);
        sfree(a2c);

        if (debug)
        {
            fprintf(debug, "molecule type '%s' %d cgs has %d cg links through bonded interac.\n", *molt->name, cgs->nr, nlink_mol);
        }

        if (molb->nmol > mol)
        {
            /* Copy the data for the rest of the molecules in this block */
            link->nalloc_a += (molb->nmol - mol)*nlink_mol;
            srenew(link->a, link->nalloc_a);
            for (; mol < molb->nmol; mol++)
            {
                for (cg = 0; cg < cgs->nr; cg++)
                {
                    cg_gl                = cg_offset + cg;
                    link->index[cg_gl+1] =
                        link->index[cg_gl+1-cgs->nr] + nlink_mol;
                    for (j = link->index[cg_gl]; j < link->index[cg_gl+1]; j++)
                    {
                        link->a[j] = link->a[j-nlink_mol] + cgs->nr;
                    }
                    if (link->index[cg_gl+1] - link->index[cg_gl] > 0 &&
                        cg_gl - cgi_mb->cg_start < cgi_mb->cg_mod)
                    {
                        SET_CGINFO_BOND_INTER(cgi_mb->cginfo[cg_gl - cgi_mb->cg_start]);
                        ncgi++;
                    }
                }
                cg_offset += cgs->nr;
            }
        }
    }

    if (mtop->bIntermolecularInteractions)
    {
        destroy_reverse_ilist(&ril_intermol);
    }

    if (debug)
    {
        fprintf(debug, "Of the %d charge groups %d are linked via bonded interactions\n", ncg_mtop(mtop), ncgi);
    }

    return link;
}

/*! \brief Set the distance, function type and atom indices for the longest distance between molecules of molecule type \p molt for two-body and multi-body bonded interactions */
static void bonded_cg_distance_mol(gmx_moltype_t *molt, int *at2cg,
                                   gmx_bool bBCheck, gmx_bool bExcl, rvec *cg_cm,
                                   real *r_2b, int *ft2b, int *a2_1, int *a2_2,
                                   real *r_mb, int *ftmb, int *am_1, int *am_2)
{
    int       ftype, nral, i, j, ai, aj, cgi, cgj;
    t_ilist  *il;
    t_blocka *excls;
    real      r2_2b, r2_mb, rij2;

    r2_2b = 0;
    r2_mb = 0;
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, bBCheck, FALSE, FALSE))
        {
            il   = &molt->ilist[ftype];
            nral = NRAL(ftype);
            if (nral > 1)
            {
                for (i = 0; i < il->nr; i += 1+nral)
                {
                    for (ai = 0; ai < nral; ai++)
                    {
                        cgi = at2cg[il->iatoms[i+1+ai]];
                        for (aj = 0; aj < nral; aj++)
                        {
                            cgj = at2cg[il->iatoms[i+1+aj]];
                            if (cgi != cgj)
                            {
                                rij2 = distance2(cg_cm[cgi], cg_cm[cgj]);
                                if (nral == 2 && rij2 > r2_2b)
                                {
                                    r2_2b = rij2;
                                    *ft2b = ftype;
                                    *a2_1 = il->iatoms[i+1+ai];
                                    *a2_2 = il->iatoms[i+1+aj];
                                }
                                if (nral >  2 && rij2 > r2_mb)
                                {
                                    r2_mb = rij2;
                                    *ftmb = ftype;
                                    *am_1 = il->iatoms[i+1+ai];
                                    *am_2 = il->iatoms[i+1+aj];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (bExcl)
    {
        excls = &molt->excls;
        for (ai = 0; ai < excls->nr; ai++)
        {
            cgi = at2cg[ai];
            for (j = excls->index[ai]; j < excls->index[ai+1]; j++)
            {
                cgj = at2cg[excls->a[j]];
                if (cgi != cgj)
                {
                    rij2 = distance2(cg_cm[cgi], cg_cm[cgj]);
                    if (rij2 > r2_2b)
                    {
                        r2_2b = rij2;
                    }
                }
            }
        }
    }

    *r_2b = sqrt(r2_2b);
    *r_mb = sqrt(r2_mb);
}

//! Compute charge group centers of mass for molecule \p molt
static void get_cgcm_mol(gmx_moltype_t *molt, gmx_ffparams_t *ffparams,
                         int ePBC, t_graph *graph, matrix box,
                         gmx_vsite_t *vsite,
                         rvec *x, rvec *xs, rvec *cg_cm)
{
    int n, i;

    if (ePBC != epbcNONE)
    {
        mk_mshift(NULL, graph, ePBC, box, x);

        shift_x(graph, box, x, xs);
        /* By doing an extra mk_mshift the molecules that are broken
         * because they were e.g. imported from another software
         * will be made whole again. Such are the healing powers
         * of GROMACS.
         */
        mk_mshift(NULL, graph, ePBC, box, xs);
    }
    else
    {
        /* We copy the coordinates so the original coordinates remain
         * unchanged, just to be 100% sure that we do not affect
         * binary reproducibility of simulations.
         */
        n = molt->cgs.index[molt->cgs.nr];
        for (i = 0; i < n; i++)
        {
            copy_rvec(x[i], xs[i]);
        }
    }

    if (vsite)
    {
        construct_vsites(vsite, xs, 0.0, NULL,
                         ffparams->iparams, molt->ilist,
                         epbcNONE, TRUE, NULL, NULL);
    }

    calc_cgcm(NULL, 0, molt->cgs.nr, &molt->cgs, xs, cg_cm);
}

//! Returns whether \p molt has a virtual site
static int have_vsite_molt(gmx_moltype_t *molt)
{
    int      i;
    gmx_bool bVSite;

    bVSite = FALSE;
    for (i = 0; i < F_NRE; i++)
    {
        if ((interaction_function[i].flags & IF_VSITE) &&
            molt->ilist[i].nr > 0)
        {
            bVSite = TRUE;
        }
    }

    return bVSite;
}

void dd_bonded_cg_distance(FILE *fplog,
                           gmx_mtop_t *mtop,
                           t_inputrec *ir, rvec *x, matrix box,
                           gmx_bool bBCheck,
                           real *r_2b, real *r_mb)
{
    gmx_bool        bExclRequired;
    int             mb, at_offset, *at2cg, mol;
    t_graph         graph;
    gmx_vsite_t    *vsite;
    gmx_molblock_t *molb;
    gmx_moltype_t  *molt;
    rvec           *xs, *cg_cm;
    real            rmol_2b, rmol_mb;
    int             ft2b  = -1, a_2b_1 = -1, a_2b_2 = -1, ftmb = -1, a_mb_1 = -1, a_mb_2 = -1;
    int             ftm2b = -1, amol_2b_1 = -1, amol_2b_2 = -1, ftmmb = -1, amol_mb_1 = -1, amol_mb_2 = -1;

    bExclRequired = IR_EXCL_FORCES(*ir);

    vsite = init_vsite(mtop, NULL, TRUE);

    *r_2b     = 0;
    *r_mb     = 0;
    at_offset = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        if (molt->cgs.nr == 1 || molb->nmol == 0)
        {
            at_offset += molb->nmol*molt->atoms.nr;
        }
        else
        {
            if (ir->ePBC != epbcNONE)
            {
                mk_graph_ilist(NULL, molt->ilist, 0, molt->atoms.nr, FALSE, FALSE,
                               &graph);
            }

            at2cg = make_at2cg(&molt->cgs);
            snew(xs, molt->atoms.nr);
            snew(cg_cm, molt->cgs.nr);
            for (mol = 0; mol < molb->nmol; mol++)
            {
                get_cgcm_mol(molt, &mtop->ffparams, ir->ePBC, &graph, box,
                             have_vsite_molt(molt) ? vsite : NULL,
                             x+at_offset, xs, cg_cm);

                bonded_cg_distance_mol(molt, at2cg, bBCheck, bExclRequired, cg_cm,
                                       &rmol_2b, &ftm2b, &amol_2b_1, &amol_2b_2,
                                       &rmol_mb, &ftmmb, &amol_mb_1, &amol_mb_2);
                if (rmol_2b > *r_2b)
                {
                    *r_2b  = rmol_2b;
                    ft2b   = ftm2b;
                    a_2b_1 = at_offset + amol_2b_1;
                    a_2b_2 = at_offset + amol_2b_2;
                }
                if (rmol_mb > *r_mb)
                {
                    *r_mb  = rmol_mb;
                    ftmb   = ftmmb;
                    a_mb_1 = at_offset + amol_mb_1;
                    a_mb_2 = at_offset + amol_mb_2;
                }

                at_offset += molt->atoms.nr;
            }
            sfree(cg_cm);
            sfree(xs);
            sfree(at2cg);
            if (ir->ePBC != epbcNONE)
            {
                done_graph(&graph);
            }
        }
    }

    /* We should have a vsite free routine, but here we can simply free */
    sfree(vsite);

    if (fplog && (ft2b >= 0 || ftmb >= 0))
    {
        fprintf(fplog,
                "Initial maximum inter charge-group distances:\n");
        if (ft2b >= 0)
        {
            fprintf(fplog,
                    "    two-body bonded interactions: %5.3f nm, %s, atoms %d %d\n",
                    *r_2b, interaction_function[ft2b].longname,
                    a_2b_1+1, a_2b_2+1);
        }
        if (ftmb >= 0)
        {
            fprintf(fplog,
                    "  multi-body bonded interactions: %5.3f nm, %s, atoms %d %d\n",
                    *r_mb, interaction_function[ftmb].longname,
                    a_mb_1+1, a_mb_2+1);
        }
    }
}
