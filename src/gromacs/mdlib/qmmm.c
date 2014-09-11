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

#include "gromacs/legacyheaders/qmmm.h"

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* declarations of the interfaces to the QM packages. The _SH indicate
 * the QM interfaces can be used for Surface Hopping simulations
 */
#ifdef GMX_QMMM_GAMESS
/* GAMESS interface */

void
init_gamess(t_commrec *cr, t_QMrec *qm, t_MMrec *mm);

real
call_gamess(t_commrec *cr, t_forcerec *fr,
            t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[]);

#elif defined GMX_QMMM_MOPAC
/* MOPAC interface */

void
init_mopac(t_commrec *cr, t_QMrec *qm, t_MMrec *mm);

real
call_mopac(t_commrec *cr, t_forcerec *fr, t_QMrec *qm,
           t_MMrec *mm, rvec f[], rvec fshift[]);

real
call_mopac_SH(t_commrec *cr, t_forcerec *fr, t_QMrec *qm,
              t_MMrec *mm, rvec f[], rvec fshift[]);

#elif defined GMX_QMMM_GAUSSIAN
/* GAUSSIAN interface */

void
init_gaussian(t_commrec *cr, t_QMrec *qm, t_MMrec *mm);

real
call_gaussian_SH(t_commrec *cr, t_forcerec *fr, t_QMrec *qm,
                 t_MMrec *mm, rvec f[], rvec fshift[]);

real
call_gaussian(t_commrec *cr, t_forcerec *fr, t_QMrec *qm,
              t_MMrec *mm, rvec f[], rvec fshift[]);

#elif defined GMX_QMMM_ORCA
/* ORCA interface */

void
init_orca(t_QMrec *qm);

real
call_orca(t_forcerec *fr, t_QMrec *qm,
          t_MMrec *mm, rvec f[], rvec fshift[]);

#endif




/* this struct and these comparison functions are needed for creating
 * a QMMM input for the QM routines from the QMMM neighbor list.
 */

typedef struct {
    int      j;
    int      shift;
} t_j_particle;

static int struct_comp(const void *a, const void *b)
{

    return (int)(((t_j_particle *)a)->j)-(int)(((t_j_particle *)b)->j);

} /* struct_comp */

static int int_comp(const void *a, const void *b)
{

    return (*(int *)a) - (*(int *)b);

} /* int_comp */

static int QMlayer_comp(const void *a, const void *b)
{

    return (int)(((t_QMrec *)a)->nrQMatoms)-(int)(((t_QMrec *)b)->nrQMatoms);

} /* QMlayer_comp */

real call_QMroutine(t_commrec gmx_unused *cr, t_forcerec gmx_unused *fr, t_QMrec gmx_unused *qm,
                    t_MMrec gmx_unused *mm, rvec gmx_unused f[], rvec gmx_unused fshift[])
{
    /* makes a call to the requested QM routine (qm->QMmethod)
     * Note that f is actually the gradient, i.e. -f
     */
    real
        QMener = 0.0;

    /* do a semi-empiprical calculation */

    if (qm->QMmethod < eQMmethodRHF && !(mm->nrMMatoms))
    {
#ifdef GMX_QMMM_MOPAC
        if (qm->bSH)
        {
            QMener = call_mopac_SH(cr, fr, qm, mm, f, fshift);
        }
        else
        {
            QMener = call_mopac(cr, fr, qm, mm, f, fshift);
        }
#else
        gmx_fatal(FARGS, "Semi-empirical QM only supported with Mopac.");
#endif
    }
    else
    {
        /* do an ab-initio calculation */
        if (qm->bSH && qm->QMmethod == eQMmethodCASSCF)
        {
#ifdef GMX_QMMM_GAUSSIAN
            QMener = call_gaussian_SH(cr, fr, qm, mm, f, fshift);
#else
            gmx_fatal(FARGS, "Ab-initio Surface-hopping only supported with Gaussian.");
#endif
        }
        else
        {
#ifdef GMX_QMMM_GAMESS
            QMener = call_gamess(cr, fr, qm, mm, f, fshift);
#elif defined GMX_QMMM_GAUSSIAN
            QMener = call_gaussian(cr, fr, qm, mm, f, fshift);
#elif defined GMX_QMMM_ORCA
            QMener = call_orca(fr, qm, mm, f, fshift);
#else
            gmx_fatal(FARGS, "Ab-initio calculation only supported with Gamess, Gaussian or ORCA.");
#endif
        }
    }
    return (QMener);
}

void init_QMroutine(t_commrec gmx_unused *cr, t_QMrec gmx_unused *qm, t_MMrec gmx_unused *mm)
{
    /* makes a call to the requested QM routine (qm->QMmethod)
     */
    if (qm->QMmethod < eQMmethodRHF)
    {
#ifdef GMX_QMMM_MOPAC
        /* do a semi-empiprical calculation */
        init_mopac(cr, qm, mm);
#else
        gmx_fatal(FARGS, "Semi-empirical QM only supported with Mopac.");
#endif
    }
    else
    {
        /* do an ab-initio calculation */
#ifdef GMX_QMMM_GAMESS
        init_gamess(cr, qm, mm);
#elif defined GMX_QMMM_GAUSSIAN
        init_gaussian(cr, qm, mm);
#elif defined GMX_QMMM_ORCA
        init_orca(qm);
#else
        gmx_fatal(FARGS, "Ab-initio calculation only supported with Gamess, Gaussian or ORCA.");
#endif
    }
} /* init_QMroutine */

void update_QMMM_coord(rvec x[], t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
    /* shifts the QM and MM particles into the central box and stores
     * these shifted coordinates in the coordinate arrays of the
     * QMMMrec. These coordinates are passed on the QM subroutines.
     */
    int
        i;

    /* shift the QM atoms into the central box
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        rvec_sub(x[qm->indexQM[i]], fr->shift_vec[qm->shiftQM[i]], qm->xQM[i]);
    }
    /* also shift the MM atoms into the central box, if any
     */
    for (i = 0; i < mm->nrMMatoms; i++)
    {
        rvec_sub(x[mm->indexMM[i]], fr->shift_vec[mm->shiftMM[i]], mm->xMM[i]);
    }
} /* update_QMMM_coord */

static void punch_QMMM_excl(t_QMrec *qm, t_MMrec *mm, t_blocka *excls)
{
    /* punch a file containing the bonded interactions of each QM
     * atom with MM atoms. These need to be excluded in the QM routines
     * Only needed in case of QM/MM optimizations
     */
    FILE
       *out = NULL;
    int
        i, j, k, nrexcl = 0, *excluded = NULL, max = 0;


    out = fopen("QMMMexcl.dat", "w");

    /* this can be done more efficiently I think
     */
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        nrexcl = 0;
        for (j = excls->index[qm->indexQM[i]];
             j < excls->index[qm->indexQM[i]+1];
             j++)
        {
            for (k = 0; k < mm->nrMMatoms; k++)
            {
                if (mm->indexMM[k] == excls->a[j]) /* the excluded MM atom */
                {
                    if (nrexcl >= max)
                    {
                        max += 1000;
                        srenew(excluded, max);
                    }
                    excluded[nrexcl++] = k;
                    continue;
                }
            }
        }
        /* write to file: */
        fprintf(out, "%5d %5d\n", i+1, nrexcl);
        for (j = 0; j < nrexcl; j++)
        {
            fprintf(out, "%5d ", excluded[j]);
        }
        fprintf(out, "\n");
    }
    free(excluded);
    fclose(out);
} /* punch_QMMM_excl */


/* end of QMMM subroutines */

/* QMMM core routines */

t_QMrec *mk_QMrec(void)
{
    t_QMrec *qm;
    snew(qm, 1);
    return qm;
} /* mk_QMrec */

t_MMrec *mk_MMrec(void)
{
    t_MMrec *mm;
    snew(mm, 1);
    return mm;
} /* mk_MMrec */

static void init_QMrec(int grpnr, t_QMrec *qm, int nr, int *atomarray,
                       gmx_mtop_t *mtop, t_inputrec *ir)
{
    /* fills the t_QMrec struct of QM group grpnr
     */
    int                   i;
    gmx_mtop_atomlookup_t alook;
    t_atom               *atom;


    qm->nrQMatoms = nr;
    snew(qm->xQM, nr);
    snew(qm->indexQM, nr);
    snew(qm->shiftQM, nr); /* the shifts */
    for (i = 0; i < nr; i++)
    {
        qm->indexQM[i] = atomarray[i];
    }

    alook = gmx_mtop_atomlookup_init(mtop);

    snew(qm->atomicnumberQM, nr);
    for (i = 0; i < qm->nrQMatoms; i++)
    {
        gmx_mtop_atomnr_to_atom(alook, qm->indexQM[i], &atom);
        qm->nelectrons       += mtop->atomtypes.atomnumber[atom->type];
        qm->atomicnumberQM[i] = mtop->atomtypes.atomnumber[atom->type];
    }

    gmx_mtop_atomlookup_destroy(alook);

    qm->QMcharge       = ir->opts.QMcharge[grpnr];
    qm->multiplicity   = ir->opts.QMmult[grpnr];
    qm->nelectrons    -= ir->opts.QMcharge[grpnr];

    qm->QMmethod       = ir->opts.QMmethod[grpnr];
    qm->QMbasis        = ir->opts.QMbasis[grpnr];
    /* trajectory surface hopping setup (Gaussian only) */
    qm->bSH            = ir->opts.bSH[grpnr];
    qm->CASorbitals    = ir->opts.CASorbitals[grpnr];
    qm->CASelectrons   = ir->opts.CASelectrons[grpnr];
    qm->SAsteps        = ir->opts.SAsteps[grpnr];
    qm->SAon           = ir->opts.SAon[grpnr];
    qm->SAoff          = ir->opts.SAoff[grpnr];
    /* hack to prevent gaussian from reinitializing all the time */
    qm->nQMcpus        = 0; /* number of CPU's to be used by g01, is set
                             * upon initializing gaussian
                             * (init_gaussian()
                             */
    /* print the current layer to allow users to check their input */
    fprintf(stderr, "Layer %d\nnr of QM atoms %d\n", grpnr, nr);
    fprintf(stderr, "QMlevel: %s/%s\n\n",
            eQMmethod_names[qm->QMmethod], eQMbasis_names[qm->QMbasis]);

    /* frontier atoms */
    snew(qm->frontatoms, nr);
    /* Lennard-Jones coefficients */
    snew(qm->c6, nr);
    snew(qm->c12, nr);
    /* do we optimize the QM separately using the algorithms of the QM program??
     */
    qm->bTS      = ir->opts.bTS[grpnr];
    qm->bOPT     = ir->opts.bOPT[grpnr];

} /* init_QMrec */

t_QMrec *copy_QMrec(t_QMrec *qm)
{
    /* copies the contents of qm into a new t_QMrec struct */
    t_QMrec
       *qmcopy;
    int
        i;

    qmcopy            = mk_QMrec();
    qmcopy->nrQMatoms = qm->nrQMatoms;
    snew(qmcopy->xQM, qmcopy->nrQMatoms);
    snew(qmcopy->indexQM, qmcopy->nrQMatoms);
    snew(qmcopy->atomicnumberQM, qm->nrQMatoms);
    snew(qmcopy->shiftQM, qmcopy->nrQMatoms); /* the shifts */
    for (i = 0; i < qmcopy->nrQMatoms; i++)
    {
        qmcopy->shiftQM[i]        = qm->shiftQM[i];
        qmcopy->indexQM[i]        = qm->indexQM[i];
        qmcopy->atomicnumberQM[i] = qm->atomicnumberQM[i];
    }
    qmcopy->nelectrons   = qm->nelectrons;
    qmcopy->multiplicity = qm->multiplicity;
    qmcopy->QMcharge     = qm->QMcharge;
    qmcopy->nelectrons   = qm->nelectrons;
    qmcopy->QMmethod     = qm->QMmethod;
    qmcopy->QMbasis      = qm->QMbasis;
    /* trajectory surface hopping setup (Gaussian only) */
    qmcopy->bSH          = qm->bSH;
    qmcopy->CASorbitals  = qm->CASorbitals;
    qmcopy->CASelectrons = qm->CASelectrons;
    qmcopy->SAsteps      = qm->SAsteps;
    qmcopy->SAon         = qm->SAon;
    qmcopy->SAoff        = qm->SAoff;
    qmcopy->bOPT         = qm->bOPT;

    /* Gaussian init. variables */
    qmcopy->nQMcpus      = qm->nQMcpus;
    for (i = 0; i < DIM; i++)
    {
        qmcopy->SHbasis[i] = qm->SHbasis[i];
    }
    qmcopy->QMmem        = qm->QMmem;
    qmcopy->accuracy     = qm->accuracy;
    qmcopy->cpmcscf      = qm->cpmcscf;
    qmcopy->SAstep       = qm->SAstep;
    snew(qmcopy->frontatoms, qm->nrQMatoms);
    snew(qmcopy->c12, qmcopy->nrQMatoms);
    snew(qmcopy->c6, qmcopy->nrQMatoms);
    if (qmcopy->bTS || qmcopy->bOPT)
    {
        for (i = 1; i < qmcopy->nrQMatoms; i++)
        {
            qmcopy->frontatoms[i] = qm->frontatoms[i];
            qmcopy->c12[i]        = qm->c12[i];
            qmcopy->c6[i]         = qm->c6[i];
        }
    }

    return(qmcopy);

} /*copy_QMrec */

t_QMMMrec *mk_QMMMrec(void)
{

    t_QMMMrec *qr;

    snew(qr, 1);

    return qr;

} /* mk_QMMMrec */

void init_QMMMrec(t_commrec  *cr,
                  gmx_mtop_t *mtop,
                  t_inputrec *ir,
                  t_forcerec *fr)
{
    /* we put the atomsnumbers of atoms that belong to the QMMM group in
     * an array that will be copied later to QMMMrec->indexQM[..]. Also
     * it will be used to create an QMMMrec->bQMMM index array that
     * simply contains true/false for QM and MM (the other) atoms.
     */

    gmx_groups_t            *groups;
    atom_id                 *qm_arr = NULL, vsite, ai, aj;
    int                      qm_max = 0, qm_nr = 0, i, j, jmax, k, l, nrvsite2 = 0;
    t_QMMMrec               *qr;
    t_MMrec                 *mm;
    t_iatom                 *iatoms;
    real                     c12au, c6au;
    gmx_mtop_atomloop_all_t  aloop;
    t_atom                  *atom;
    gmx_mtop_ilistloop_all_t iloop;
    int                      a_offset;
    t_ilist                 *ilist_mol;
    gmx_mtop_atomlookup_t    alook;

    c6au  = (HARTREE2KJ*AVOGADRO*pow(BOHR2NM, 6));
    c12au = (HARTREE2KJ*AVOGADRO*pow(BOHR2NM, 12));
    /* issue a fatal if the user wants to run with more than one node */
    if (PAR(cr))
    {
        gmx_fatal(FARGS, "QM/MM does not work in parallel, use a single rank instead\n");
    }

    /* Make a local copy of the QMMMrec */
    qr = fr->qr;

    /* bQMMM[..] is an array containing TRUE/FALSE for atoms that are
     * QM/not QM. We first set all elemenst at false. Afterwards we use
     * the qm_arr (=MMrec->indexQM) to changes the elements
     * corresponding to the QM atoms at TRUE.  */

    qr->QMMMscheme     = ir->QMMMscheme;

    /* we take the possibility into account that a user has
     * defined more than one QM group:
     */
    /* an ugly work-around in case there is only one group In this case
     * the whole system is treated as QM. Otherwise the second group is
     * always the rest of the total system and is treated as MM.
     */

    /* small problem if there is only QM.... so no MM */

    jmax = ir->opts.ngQM;

    if (qr->QMMMscheme == eQMMMschemeoniom)
    {
        qr->nrQMlayers = jmax;
    }
    else
    {
        qr->nrQMlayers = 1;
    }

    groups = &mtop->groups;

    /* there are jmax groups of QM atoms. In case of multiple QM groups
     * I assume that the users wants to do ONIOM. However, maybe it
     * should also be possible to define more than one QM subsystem with
     * independent neighbourlists. I have to think about
     * that.. 11-11-2003
     */
    snew(qr->qm, jmax);
    for (j = 0; j < jmax; j++)
    {
        /* new layer */
        aloop = gmx_mtop_atomloop_all_init(mtop);
        while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
        {
            if (qm_nr >= qm_max)
            {
                qm_max += 1000;
                srenew(qm_arr, qm_max);
            }
            if (ggrpnr(groups, egcQMMM, i) == j)
            {
                /* hack for tip4p */
                qm_arr[qm_nr++] = i;
            }
        }
        if (qr->QMMMscheme == eQMMMschemeoniom)
        {
            /* add the atoms to the bQMMM array
             */

            /* I assume that users specify the QM groups from small to
             * big(ger) in the mdp file
             */
            qr->qm[j] = mk_QMrec();
            /* we need to throw out link atoms that in the previous layer
             * existed to separate this QMlayer from the previous
             * QMlayer. We use the iatoms array in the idef for that
             * purpose. If all atoms defining the current Link Atom (Dummy2)
             * are part of the current QM layer it needs to be removed from
             * qm_arr[].  */

            iloop = gmx_mtop_ilistloop_all_init(mtop);
            while (gmx_mtop_ilistloop_all_next(iloop, &ilist_mol, &a_offset))
            {
                nrvsite2 = ilist_mol[F_VSITE2].nr;
                iatoms   = ilist_mol[F_VSITE2].iatoms;

                for (k = 0; k < nrvsite2; k += 4)
                {
                    vsite = a_offset + iatoms[k+1]; /* the vsite         */
                    ai    = a_offset + iatoms[k+2]; /* constructing atom */
                    aj    = a_offset + iatoms[k+3]; /* constructing atom */
                    if (ggrpnr(groups, egcQMMM, vsite) == ggrpnr(groups, egcQMMM, ai)
                        &&
                        ggrpnr(groups, egcQMMM, vsite) == ggrpnr(groups, egcQMMM, aj))
                    {
                        /* this dummy link atom needs to be removed from the qm_arr
                         * before making the QMrec of this layer!
                         */
                        for (i = 0; i < qm_nr; i++)
                        {
                            if (qm_arr[i] == vsite)
                            {
                                /* drop the element */
                                for (l = i; l < qm_nr; l++)
                                {
                                    qm_arr[l] = qm_arr[l+1];
                                }
                                qm_nr--;
                            }
                        }
                    }
                }
            }

            /* store QM atoms in this layer in the QMrec and initialise layer
             */
            init_QMrec(j, qr->qm[j], qm_nr, qm_arr, mtop, ir);

            /* we now store the LJ C6 and C12 parameters in QM rec in case
             * we need to do an optimization
             */
            if (qr->qm[j]->bOPT || qr->qm[j]->bTS)
            {
                for (i = 0; i < qm_nr; i++)
                {
                    /* nbfp now includes the 6.0/12.0 derivative prefactors */
                    qr->qm[j]->c6[i]  =  C6(fr->nbfp, mtop->ffparams.atnr, atom->type, atom->type)/c6au/6.0;
                    qr->qm[j]->c12[i] = C12(fr->nbfp, mtop->ffparams.atnr, atom->type, atom->type)/c12au/12.0;
                }
            }
            /* now we check for frontier QM atoms. These occur in pairs that
             * construct the vsite
             */
            iloop = gmx_mtop_ilistloop_all_init(mtop);
            while (gmx_mtop_ilistloop_all_next(iloop, &ilist_mol, &a_offset))
            {
                nrvsite2 = ilist_mol[F_VSITE2].nr;
                iatoms   = ilist_mol[F_VSITE2].iatoms;

                for (k = 0; k < nrvsite2; k += 4)
                {
                    vsite = a_offset + iatoms[k+1]; /* the vsite         */
                    ai    = a_offset + iatoms[k+2]; /* constructing atom */
                    aj    = a_offset + iatoms[k+3]; /* constructing atom */
                    if (ggrpnr(groups, egcQMMM, ai) < (groups->grps[egcQMMM].nr-1) &&
                        (ggrpnr(groups, egcQMMM, aj) >= (groups->grps[egcQMMM].nr-1)))
                    {
                        /* mark ai as frontier atom */
                        for (i = 0; i < qm_nr; i++)
                        {
                            if ( (qm_arr[i] == ai) || (qm_arr[i] == vsite) )
                            {
                                qr->qm[j]->frontatoms[i] = TRUE;
                            }
                        }
                    }
                    else if (ggrpnr(groups, egcQMMM, aj) < (groups->grps[egcQMMM].nr-1) &&
                             (ggrpnr(groups, egcQMMM, ai) >= (groups->grps[egcQMMM].nr-1)))
                    {
                        /* mark aj as frontier atom */
                        for (i = 0; i < qm_nr; i++)
                        {
                            if ( (qm_arr[i] == aj) || (qm_arr[i] == vsite))
                            {
                                qr->qm[j]->frontatoms[i] = TRUE;
                            }
                        }
                    }
                }
            }
        }
    }
    if (qr->QMMMscheme != eQMMMschemeoniom)
    {

        /* standard QMMM, all layers are merged together so there is one QM
         * subsystem and one MM subsystem.
         * Also we set the charges to zero in the md->charge arrays to prevent
         * the innerloops from doubly counting the electostatic QM MM interaction
         */

        alook = gmx_mtop_atomlookup_init(mtop);

        for (k = 0; k < qm_nr; k++)
        {
            gmx_mtop_atomnr_to_atom(alook, qm_arr[k], &atom);
            atom->q  = 0.0;
            atom->qB = 0.0;
        }
        qr->qm[0] = mk_QMrec();
        /* store QM atoms in the QMrec and initialise
         */
        init_QMrec(0, qr->qm[0], qm_nr, qm_arr, mtop, ir);
        if (qr->qm[0]->bOPT || qr->qm[0]->bTS)
        {
            for (i = 0; i < qm_nr; i++)
            {
                gmx_mtop_atomnr_to_atom(alook, qm_arr[i], &atom);
                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                qr->qm[0]->c6[i]  =  C6(fr->nbfp, mtop->ffparams.atnr, atom->type, atom->type)/c6au/6.0;
                qr->qm[0]->c12[i] = C12(fr->nbfp, mtop->ffparams.atnr, atom->type, atom->type)/c12au/12.0;
            }
        }

        /* find frontier atoms and mark them true in the frontieratoms array.
         */
        for (i = 0; i < qm_nr; i++)
        {
            gmx_mtop_atomnr_to_ilist(alook, qm_arr[i], &ilist_mol, &a_offset);
            nrvsite2 = ilist_mol[F_VSITE2].nr;
            iatoms   = ilist_mol[F_VSITE2].iatoms;

            for (k = 0; k < nrvsite2; k += 4)
            {
                vsite = a_offset + iatoms[k+1]; /* the vsite         */
                ai    = a_offset + iatoms[k+2]; /* constructing atom */
                aj    = a_offset + iatoms[k+3]; /* constructing atom */
                if (ggrpnr(groups, egcQMMM, ai) < (groups->grps[egcQMMM].nr-1) &&
                    (ggrpnr(groups, egcQMMM, aj) >= (groups->grps[egcQMMM].nr-1)))
                {
                    /* mark ai as frontier atom */
                    if ( (qm_arr[i] == ai) || (qm_arr[i] == vsite) )
                    {
                        qr->qm[0]->frontatoms[i] = TRUE;
                    }
                }
                else if (ggrpnr(groups, egcQMMM, aj) < (groups->grps[egcQMMM].nr-1) &&
                         (ggrpnr(groups, egcQMMM, ai) >= (groups->grps[egcQMMM].nr-1)))
                {
                    /* mark aj as frontier atom */
                    if ( (qm_arr[i] == aj) || (qm_arr[i] == vsite) )
                    {
                        qr->qm[0]->frontatoms[i] = TRUE;
                    }
                }
            }
        }

        gmx_mtop_atomlookup_destroy(alook);

        /* MM rec creation */
        mm               = mk_MMrec();
        mm->scalefactor  = ir->scalefactor;
        mm->nrMMatoms    = (mtop->natoms)-(qr->qm[0]->nrQMatoms); /* rest of the atoms */
        qr->mm           = mm;
    }
    else /* ONIOM */
    {    /* MM rec creation */
        mm               = mk_MMrec();
        mm->scalefactor  = ir->scalefactor;
        mm->nrMMatoms    = 0;
        qr->mm           = mm;
    }

    /* these variables get updated in the update QMMMrec */

    if (qr->nrQMlayers == 1)
    {
        /* with only one layer there is only one initialisation
         * needed. Multilayer is a bit more complicated as it requires
         * re-initialisation at every step of the simulation. This is due
         * to the use of COMMON blocks in the fortran QM subroutines.
         */
        if (qr->qm[0]->QMmethod < eQMmethodRHF)
        {
#ifdef GMX_QMMM_MOPAC
            /* semi-empiprical 1-layer ONIOM calculation requested (mopac93) */
            init_mopac(cr, qr->qm[0], qr->mm);
#else
            gmx_fatal(FARGS, "Semi-empirical QM only supported with Mopac.");
#endif
        }
        else
        {
            /* ab initio calculation requested (gamess/gaussian/ORCA) */
#ifdef GMX_QMMM_GAMESS
            init_gamess(cr, qr->qm[0], qr->mm);
#elif defined GMX_QMMM_GAUSSIAN
            init_gaussian(cr, qr->qm[0], qr->mm);
#elif defined GMX_QMMM_ORCA
            init_orca(qr->qm[0]);
#else
            gmx_fatal(FARGS, "Ab-initio calculation only supported with Gamess, Gaussian or ORCA.");
#endif
        }
    }
} /* init_QMMMrec */

void update_QMMMrec(t_commrec      *cr,
                    t_forcerec     *fr,
                    rvec            x[],
                    t_mdatoms      *md,
                    matrix          box,
                    gmx_localtop_t *top)
{
    /* updates the coordinates of both QM atoms and MM atoms and stores
     * them in the QMMMrec.
     *
     * NOTE: is NOT yet working if there are no PBC. Also in ns.c, simple
     * ns needs to be fixed!
     */
    int
        mm_max = 0, mm_nr = 0, mm_nr_new, i, j, is, k, shift;
    t_j_particle
       *mm_j_particles = NULL, *qm_i_particles = NULL;
    t_QMMMrec
       *qr;
    t_nblist
        QMMMlist;
    rvec
        dx, crd;
    int
       *MMatoms;
    t_QMrec
       *qm;
    t_MMrec
       *mm;
    t_pbc
        pbc;
    int
       *parallelMMarray = NULL;
    real
        c12au, c6au;

    c6au  = (HARTREE2KJ*AVOGADRO*pow(BOHR2NM, 6));
    c12au = (HARTREE2KJ*AVOGADRO*pow(BOHR2NM, 12));

    /* every cpu has this array. On every processor we fill this array
     * with 1's and 0's. 1's indicate the atoms is a QM atom on the
     * current cpu in a later stage these arrays are all summed. indexes
     * > 0 indicate the atom is a QM atom. Every node therefore knows
     * whcih atoms are part of the QM subsystem.
     */
    /* copy some pointers */
    qr          = fr->qr;
    mm          = qr->mm;
    QMMMlist    = fr->QMMMlist;



    /*  init_pbc(box);  needs to be called first, see pbc.h */
    set_pbc_dd(&pbc, fr->ePBC, DOMAINDECOMP(cr) ? cr->dd : NULL, FALSE, box);
    /* only in standard (normal) QMMM we need the neighbouring MM
     * particles to provide a electric field of point charges for the QM
     * atoms.
     */
    if (qr->QMMMscheme == eQMMMschemenormal) /* also implies 1 QM-layer */
    {
        /* we NOW create/update a number of QMMMrec entries:
         *
         * 1) the shiftQM, containing the shifts of the QM atoms
         *
         * 2) the indexMM array, containing the index of the MM atoms
         *
         * 3) the shiftMM, containing the shifts of the MM atoms
         *
         * 4) the shifted coordinates of the MM atoms
         *
         * the shifts are used for computing virial of the QM/MM particles.
         */
        qm = qr->qm[0]; /* in case of normal QMMM, there is only one group */
        snew(qm_i_particles, QMMMlist.nri);
        if (QMMMlist.nri)
        {
            qm_i_particles[0].shift = XYZ2IS(0, 0, 0);
            for (i = 0; i < QMMMlist.nri; i++)
            {
                qm_i_particles[i].j     = QMMMlist.iinr[i];

                if (i)
                {
                    qm_i_particles[i].shift = pbc_dx_aiuc(&pbc, x[QMMMlist.iinr[0]],
                                                          x[QMMMlist.iinr[i]], dx);

                }
                /* However, since nri >= nrQMatoms, we do a quicksort, and throw
                 * out double, triple, etc. entries later, as we do for the MM
                 * list too.
                 */

                /* compute the shift for the MM j-particles with respect to
                 * the QM i-particle and store them.
                 */

                crd[0] = IS2X(QMMMlist.shift[i]) + IS2X(qm_i_particles[i].shift);
                crd[1] = IS2Y(QMMMlist.shift[i]) + IS2Y(qm_i_particles[i].shift);
                crd[2] = IS2Z(QMMMlist.shift[i]) + IS2Z(qm_i_particles[i].shift);
                is     = XYZ2IS(crd[0], crd[1], crd[2]);
                for (j = QMMMlist.jindex[i];
                     j < QMMMlist.jindex[i+1];
                     j++)
                {
                    if (mm_nr >= mm_max)
                    {
                        mm_max += 1000;
                        srenew(mm_j_particles, mm_max);
                    }

                    mm_j_particles[mm_nr].j     = QMMMlist.jjnr[j];
                    mm_j_particles[mm_nr].shift = is;
                    mm_nr++;
                }
            }

            /* quicksort QM and MM shift arrays and throw away multiple entries */



            qsort(qm_i_particles, QMMMlist.nri,
                  (size_t)sizeof(qm_i_particles[0]),
                  struct_comp);
            qsort(mm_j_particles, mm_nr,
                  (size_t)sizeof(mm_j_particles[0]),
                  struct_comp);
            /* remove multiples in the QM shift array, since in init_QMMM() we
             * went through the atom numbers from 0 to md.nr, the order sorted
             * here matches the one of QMindex already.
             */
            j = 0;
            for (i = 0; i < QMMMlist.nri; i++)
            {
                if (i == 0 || qm_i_particles[i].j != qm_i_particles[i-1].j)
                {
                    qm_i_particles[j++] = qm_i_particles[i];
                }
            }
            mm_nr_new = 0;
            if (qm->bTS || qm->bOPT)
            {
                /* only remove double entries for the MM array */
                for (i = 0; i < mm_nr; i++)
                {
                    if ((i == 0 || mm_j_particles[i].j != mm_j_particles[i-1].j)
                        && !md->bQM[mm_j_particles[i].j])
                    {
                        mm_j_particles[mm_nr_new++] = mm_j_particles[i];
                    }
                }
            }
            /* we also remove mm atoms that have no charges!
             * actually this is already done in the ns.c
             */
            else
            {
                for (i = 0; i < mm_nr; i++)
                {
                    if ((i == 0 || mm_j_particles[i].j != mm_j_particles[i-1].j)
                        && !md->bQM[mm_j_particles[i].j]
                        && (md->chargeA[mm_j_particles[i].j]
                            || (md->chargeB && md->chargeB[mm_j_particles[i].j])))
                    {
                        mm_j_particles[mm_nr_new++] = mm_j_particles[i];
                    }
                }
            }
            mm_nr = mm_nr_new;
            /* store the data retrieved above into the QMMMrec
             */
            k = 0;
            /* Keep the compiler happy,
             * shift will always be set in the loop for i=0
             */
            shift = 0;
            for (i = 0; i < qm->nrQMatoms; i++)
            {
                /* not all qm particles might have appeared as i
                 * particles. They might have been part of the same charge
                 * group for instance.
                 */
                if (qm->indexQM[i] == qm_i_particles[k].j)
                {
                    shift = qm_i_particles[k++].shift;
                }
                /* use previous shift, assuming they belong the same charge
                 * group anyway,
                 */

                qm->shiftQM[i] = shift;
            }
        }
        /* parallel excecution */
        if (PAR(cr))
        {
            snew(parallelMMarray, 2*(md->nr));
            /* only MM particles have a 1 at their atomnumber. The second part
             * of the array contains the shifts. Thus:
             * p[i]=1/0 depending on wether atomnumber i is a MM particle in the QM
             * step or not. p[i+md->nr] is the shift of atomnumber i.
             */
            for (i = 0; i < 2*(md->nr); i++)
            {
                parallelMMarray[i] = 0;
            }

            for (i = 0; i < mm_nr; i++)
            {
                parallelMMarray[mm_j_particles[i].j]          = 1;
                parallelMMarray[mm_j_particles[i].j+(md->nr)] = mm_j_particles[i].shift;
            }
            gmx_sumi(md->nr, parallelMMarray, cr);
            mm_nr = 0;

            mm_max = 0;
            for (i = 0; i < md->nr; i++)
            {
                if (parallelMMarray[i])
                {
                    if (mm_nr >= mm_max)
                    {
                        mm_max += 1000;
                        srenew(mm->indexMM, mm_max);
                        srenew(mm->shiftMM, mm_max);
                    }
                    mm->indexMM[mm_nr]   = i;
                    mm->shiftMM[mm_nr++] = parallelMMarray[i+md->nr]/parallelMMarray[i];
                }
            }
            mm->nrMMatoms = mm_nr;
            free(parallelMMarray);
        }
        /* serial execution */
        else
        {
            mm->nrMMatoms = mm_nr;
            srenew(mm->shiftMM, mm_nr);
            srenew(mm->indexMM, mm_nr);
            for (i = 0; i < mm_nr; i++)
            {
                mm->indexMM[i] = mm_j_particles[i].j;
                mm->shiftMM[i] = mm_j_particles[i].shift;
            }

        }
        /* (re) allocate memory for the MM coordiate array. The QM
         * coordinate array was already allocated in init_QMMM, and is
         * only (re)filled in the update_QMMM_coordinates routine
         */
        srenew(mm->xMM, mm->nrMMatoms);
        /* now we (re) fill the array that contains the MM charges with
         * the forcefield charges. If requested, these charges will be
         * scaled by a factor
         */
        srenew(mm->MMcharges, mm->nrMMatoms);
        for (i = 0; i < mm->nrMMatoms; i++) /* no free energy yet */
        {
            mm->MMcharges[i] = md->chargeA[mm->indexMM[i]]*mm->scalefactor;
        }
        if (qm->bTS || qm->bOPT)
        {
            /* store (copy) the c6 and c12 parameters into the MMrec struct
             */
            srenew(mm->c6, mm->nrMMatoms);
            srenew(mm->c12, mm->nrMMatoms);
            for (i = 0; i < mm->nrMMatoms; i++)
            {
                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                mm->c6[i]  = C6(fr->nbfp, top->idef.atnr, md->typeA[mm->indexMM[i]], md->typeA[mm->indexMM[i]])/c6au/6.0;
                mm->c12[i] = C12(fr->nbfp, top->idef.atnr, md->typeA[mm->indexMM[i]], md->typeA[mm->indexMM[i]])/c12au/12.0;
            }
            punch_QMMM_excl(qr->qm[0], mm, &(top->excls));
        }
        /* the next routine fills the coordinate fields in the QMMM rec of
         * both the qunatum atoms and the MM atoms, using the shifts
         * calculated above.
         */

        update_QMMM_coord(x, fr, qr->qm[0], qr->mm);
        free(qm_i_particles);
        free(mm_j_particles);
    }
    else /* ONIOM */ /* ????? */
    {
        mm->nrMMatoms = 0;
        /* do for each layer */
        for (j = 0; j < qr->nrQMlayers; j++)
        {
            qm             = qr->qm[j];
            qm->shiftQM[0] = XYZ2IS(0, 0, 0);
            for (i = 1; i < qm->nrQMatoms; i++)
            {
                qm->shiftQM[i] = pbc_dx_aiuc(&pbc, x[qm->indexQM[0]], x[qm->indexQM[i]],
                                             dx);
            }
            update_QMMM_coord(x, fr, qm, mm);
        }
    }
} /* update_QMMM_rec */


real calculate_QMMM(t_commrec *cr,
                    rvec x[], rvec f[],
                    t_forcerec *fr)
{
    real
        QMener = 0.0;
    /* a selection for the QM package depending on which is requested
     * (Gaussian, GAMESS-UK, MOPAC or ORCA) needs to be implemented here. Now
     * it works through defines.... Not so nice yet
     */
    t_QMMMrec
    *qr;
    t_QMrec
    *qm, *qm2;
    t_MMrec
    *mm = NULL;
    rvec
    *forces  = NULL, *fshift = NULL,
    *forces2 = NULL, *fshift2 = NULL; /* needed for multilayer ONIOM */
    int
        i, j, k;
    /* make a local copy the QMMMrec pointer
     */
    qr = fr->qr;
    mm = qr->mm;

    /* now different procedures are carried out for one layer ONION and
     * normal QMMM on one hand and multilayer oniom on the other
     */
    if (qr->QMMMscheme == eQMMMschemenormal || qr->nrQMlayers == 1)
    {
        qm = qr->qm[0];
        snew(forces, (qm->nrQMatoms+mm->nrMMatoms));
        snew(fshift, (qm->nrQMatoms+mm->nrMMatoms));
        QMener = call_QMroutine(cr, fr, qm, mm, forces, fshift);
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                f[qm->indexQM[i]][j]          -= forces[i][j];
                fr->fshift[qm->shiftQM[i]][j] += fshift[i][j];
            }
        }
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                f[mm->indexMM[i]][j]          -= forces[qm->nrQMatoms+i][j];
                fr->fshift[mm->shiftMM[i]][j] += fshift[qm->nrQMatoms+i][j];
            }

        }
        free(forces);
        free(fshift);
    }
    else                                       /* Multi-layer ONIOM */
    {
        for (i = 0; i < qr->nrQMlayers-1; i++) /* last layer is special */
        {
            qm  = qr->qm[i];
            qm2 = copy_QMrec(qr->qm[i+1]);

            qm2->nrQMatoms = qm->nrQMatoms;

            for (j = 0; j < qm2->nrQMatoms; j++)
            {
                for (k = 0; k < DIM; k++)
                {
                    qm2->xQM[j][k]       = qm->xQM[j][k];
                }
                qm2->indexQM[j]        = qm->indexQM[j];
                qm2->atomicnumberQM[j] = qm->atomicnumberQM[j];
                qm2->shiftQM[j]        = qm->shiftQM[j];
            }

            qm2->QMcharge = qm->QMcharge;
            /* this layer at the higher level of theory */
            srenew(forces, qm->nrQMatoms);
            srenew(fshift, qm->nrQMatoms);
            /* we need to re-initialize the QMroutine every step... */
            init_QMroutine(cr, qm, mm);
            QMener += call_QMroutine(cr, fr, qm, mm, forces, fshift);

            /* this layer at the lower level of theory */
            srenew(forces2, qm->nrQMatoms);
            srenew(fshift2, qm->nrQMatoms);
            init_QMroutine(cr, qm2, mm);
            QMener -= call_QMroutine(cr, fr, qm2, mm, forces2, fshift2);
            /* E = E1high-E1low The next layer includes the current layer at
             * the lower level of theory, which provides + E2low
             * this is similar for gradients
             */
            for (i = 0; i < qm->nrQMatoms; i++)
            {
                for (j = 0; j < DIM; j++)
                {
                    f[qm->indexQM[i]][j]          -= (forces[i][j]-forces2[i][j]);
                    fr->fshift[qm->shiftQM[i]][j] += (fshift[i][j]-fshift2[i][j]);
                }
            }
            free(qm2);
        }
        /* now the last layer still needs to be done: */
        qm      = qr->qm[qr->nrQMlayers-1]; /* C counts from 0 */
        init_QMroutine(cr, qm, mm);
        srenew(forces, qm->nrQMatoms);
        srenew(fshift, qm->nrQMatoms);
        QMener += call_QMroutine(cr, fr, qm, mm, forces, fshift);
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                f[qm->indexQM[i]][j]          -= forces[i][j];
                fr->fshift[qm->shiftQM[i]][j] += fshift[i][j];
            }
        }
        free(forces);
        free(fshift);
        free(forces2);
        free(fshift2);
    }
    if (qm->bTS || qm->bOPT)
    {
        /* qm[0] still contains the largest ONIOM QM subsystem
         * we take the optimized coordiates and put the in x[]
         */
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                x[qm->indexQM[i]][j] = qm->xQM[i][j];
            }
        }
    }
    return(QMener);
} /* calculate_QMMM */

/* end of QMMM core routines */
