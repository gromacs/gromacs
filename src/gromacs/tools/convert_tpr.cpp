/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "convert_tpr.h"

#include <cmath>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/random/seed.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#define RANGECHK(i, n) if ((i) >= (n)) gmx_fatal(FARGS, "Your index file contains atomnumbers (e.g. %d)\nthat are larger than the number of atoms in the tpr file (%d)", (i), (n))

static gmx_bool *bKeepIt(int gnx, int natoms, int index[])
{
    gmx_bool *b;
    int       i;

    snew(b, natoms);
    for (i = 0; (i < gnx); i++)
    {
        RANGECHK(index[i], natoms);
        b[index[i]] = TRUE;
    }

    return b;
}

static int *invind(int gnx, int natoms, int index[])
{
    int     *inv;
    int      i;

    snew(inv, natoms);
    for (i = 0; (i < gnx); i++)
    {
        RANGECHK(index[i], natoms);
        inv[index[i]] = i;
    }

    return inv;
}

static void reduce_block(gmx_bool bKeep[], t_block *block,
                         const char *name)
{
    int     *index;
    int      i, j, newi, newj;

    snew(index, block->nr);

    newi = newj = 0;
    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            if (bKeep[j])
            {
                newj++;
            }
        }
        if (newj > index[newi])
        {
            newi++;
            index[newi] = newj;
        }
    }

    fprintf(stderr, "Reduced block %8s from %6d to %6d index-, %6d to %6d a-entries\n",
            name, block->nr, newi, block->index[block->nr], newj);
    block->index = index;
    block->nr    = newi;
}

static void reduce_blocka(int invindex[], gmx_bool bKeep[], t_blocka *block,
                          const char *name)
{
    int     *index, *a;
    int      i, j, k, newi, newj;

    snew(index, block->nr);
    snew(a, block->nra);

    newi = newj = 0;
    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            k = block->a[j];
            if (bKeep[k])
            {
                a[newj] = invindex[k];
                newj++;
            }
        }
        if (newj > index[newi])
        {
            newi++;
            index[newi] = newj;
        }
    }

    fprintf(stderr, "Reduced block %8s from %6d to %6d index-, %6d to %6d a-entries\n",
            name, block->nr, newi, block->nra, newj);
    block->index = index;
    block->a     = a;
    block->nr    = newi;
    block->nra   = newj;
}

static void reduce_rvec(int gnx, int index[], rvec vv[])
{
    rvec *ptr;
    int   i;

    snew(ptr, gnx);
    for (i = 0; (i < gnx); i++)
    {
        copy_rvec(vv[index[i]], ptr[i]);
    }
    for (i = 0; (i < gnx); i++)
    {
        copy_rvec(ptr[i], vv[i]);
    }
    sfree(ptr);
}

static void reduce_atom(int gnx, int index[], t_atom atom[], char ***atomname,
                        int *nres, t_resinfo *resinfo)
{
    t_atom    *ptr;
    char    ***aname;
    t_resinfo *rinfo;
    int        i, nr;

    snew(ptr, gnx);
    snew(aname, gnx);
    snew(rinfo, atom[index[gnx-1]].resind+1);
    for (i = 0; (i < gnx); i++)
    {
        ptr[i]   = atom[index[i]];
        aname[i] = atomname[index[i]];
    }
    nr = -1;
    for (i = 0; (i < gnx); i++)
    {
        atom[i]     = ptr[i];
        atomname[i] = aname[i];
        if ((i == 0) || (atom[i].resind != atom[i-1].resind))
        {
            nr++;
            rinfo[nr] = resinfo[atom[i].resind];
        }
        atom[i].resind = nr;
    }
    nr++;
    for (i = 0; (i < nr); i++)
    {
        resinfo[i] = rinfo[i];
    }
    *nres = nr;

    sfree(aname);
    sfree(ptr);
    sfree(rinfo);
}

static void reduce_ilist(int invindex[], gmx_bool bKeep[],
                         t_ilist *il, int nratoms, const char *name)
{
    t_iatom *ia;
    int      i, j, newnr;
    gmx_bool bB;

    if (il->nr)
    {
        snew(ia, il->nr);
        newnr = 0;
        for (i = 0; (i < il->nr); i += nratoms+1)
        {
            bB = TRUE;
            for (j = 1; (j <= nratoms); j++)
            {
                bB = bB && bKeep[il->iatoms[i+j]];
            }
            if (bB)
            {
                ia[newnr++] = il->iatoms[i];
                for (j = 1; (j <= nratoms); j++)
                {
                    ia[newnr++] = invindex[il->iatoms[i+j]];
                }
            }
        }
        fprintf(stderr, "Reduced ilist %8s from %6d to %6d entries\n",
                name, il->nr/(nratoms+1),
                newnr/(nratoms+1));

        il->nr = newnr;
        for (i = 0; (i < newnr); i++)
        {
            il->iatoms[i] = ia[i];
        }

        sfree(ia);
    }
}

static void reduce_topology_x(int gnx, int index[],
                              gmx_mtop_t *mtop, rvec x[], rvec v[])
{
    t_topology   top;
    gmx_bool    *bKeep;
    int         *invindex;
    int          i;

    top      = gmx_mtop_t_to_t_topology(mtop, false);
    bKeep    = bKeepIt(gnx, top.atoms.nr, index);
    invindex = invind(gnx, top.atoms.nr, index);

    reduce_block(bKeep, &(top.cgs), "cgs");
    reduce_block(bKeep, &(top.mols), "mols");
    reduce_blocka(invindex, bKeep, &(top.excls), "excls");
    reduce_rvec(gnx, index, x);
    reduce_rvec(gnx, index, v);
    reduce_atom(gnx, index, top.atoms.atom, top.atoms.atomname,
                &(top.atoms.nres), top.atoms.resinfo);

    for (i = 0; (i < F_NRE); i++)
    {
        reduce_ilist(invindex, bKeep, &(top.idef.il[i]),
                     interaction_function[i].nratoms,
                     interaction_function[i].name);
    }

    top.atoms.nr = gnx;

    mtop->nmoltype = 1;
    snew(mtop->moltype, mtop->nmoltype);
    mtop->moltype[0].name  = mtop->name;
    mtop->moltype[0].atoms = top.atoms;
    for (i = 0; i < F_NRE; i++)
    {
        mtop->moltype[0].ilist[i] = top.idef.il[i];
    }
    mtop->moltype[0].atoms = top.atoms;
    mtop->moltype[0].cgs   = top.cgs;
    mtop->moltype[0].excls = top.excls;

    mtop->nmolblock = 1;
    snew(mtop->molblock, mtop->nmolblock);
    mtop->molblock[0].type       = 0;
    mtop->molblock[0].nmol       = 1;
    mtop->molblock[0].natoms_mol = top.atoms.nr;
    mtop->molblock[0].nposres_xA = 0;
    mtop->molblock[0].nposres_xB = 0;

    mtop->natoms                 = top.atoms.nr;
}

static void zeroq(int index[], gmx_mtop_t *mtop)
{
    int mt, i;

    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        for (i = 0; (i < mtop->moltype[mt].atoms.nr); i++)
        {
            mtop->moltype[mt].atoms.atom[index[i]].q  = 0;
            mtop->moltype[mt].atoms.atom[index[i]].qB = 0;
        }
    }
}

int gmx_convert_tpr(int argc, char *argv[])
{
    const char       *desc[] = {
        "[THISMODULE] can edit run input files in three ways.[PAR]",
        "[BB]1.[bb] by modifying the number of steps in a run input file",
        "with options [TT]-extend[tt], [TT]-until[tt] or [TT]-nsteps[tt]",
        "(nsteps=-1 means unlimited number of steps)[PAR]",
        "[BB]2.[bb] by creating a [REF].tpx[ref] file for a subset of your original",
        "tpx file, which is useful when you want to remove the solvent from",
        "your [REF].tpx[ref] file, or when you want to make e.g. a pure C[GRK]alpha[grk] [REF].tpx[ref] file.",
        "Note that you may need to use [TT]-nsteps -1[tt] (or similar) to get",
        "this to work.",
        "[BB]WARNING: this [REF].tpx[ref] file is not fully functional[bb].[PAR]",
        "[BB]3.[bb] by setting the charges of a specified group",
        "to zero. This is useful when doing free energy estimates",
        "using the LIE (Linear Interaction Energy) method."
    };

    const char       *top_fn;
    int               i;
    gmx_int64_t       nsteps_req, run_step;
    double            run_t, state_t;
    gmx_bool          bSel;
    gmx_bool          bNsteps, bExtend, bUntil;
    gmx_mtop_t        mtop;
    t_atoms           atoms;
    t_state           state;
    int               gnx;
    char             *grpname;
    int              *index = nullptr;
    char              buf[200], buf2[200];
    gmx_output_env_t *oenv;
    t_filenm          fnm[] = {
        { efTPR, nullptr,  nullptr,    ffREAD  },
        { efNDX, nullptr,  nullptr,    ffOPTRD },
        { efTPR, "-o",  "tprout", ffWRITE }
    };
#define NFILE asize(fnm)

    /* Command line options */
    static int      nsteps_req_int = 0;
    static real     extend_t       = 0.0, until_t = 0.0;
    static gmx_bool bZeroQ         = FALSE;
    static t_pargs  pa[]           = {
        { "-extend",        FALSE, etREAL, {&extend_t},
          "Extend runtime by this amount (ps)" },
        { "-until",         FALSE, etREAL, {&until_t},
          "Extend runtime until this ending time (ps)" },
        { "-nsteps",        FALSE, etINT,  {&nsteps_req_int},
          "Change the number of steps" },
        { "-zeroq",         FALSE, etBOOL, {&bZeroQ},
          "Set the charges of a group (from the index) to zero" }
    };

    /* Parse the command line */
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    /* Convert int to gmx_int64_t */
    nsteps_req = nsteps_req_int;
    bNsteps    = opt2parg_bSet("-nsteps", asize(pa), pa);
    bExtend    = opt2parg_bSet("-extend", asize(pa), pa);
    bUntil     = opt2parg_bSet("-until", asize(pa), pa);

    top_fn = ftp2fn(efTPR, NFILE, fnm);
    fprintf(stderr, "Reading toplogy and stuff from %s\n", top_fn);

    t_inputrec  irInstance;
    t_inputrec *ir = &irInstance;
    read_tpx_state(top_fn, ir, &state, &mtop);
    run_step = ir->init_step;
    run_t    = ir->init_step*ir->delta_t + ir->init_t;

    if (bNsteps)
    {
        fprintf(stderr, "Setting nsteps to %s\n", gmx_step_str(nsteps_req, buf));
        ir->nsteps = nsteps_req;
    }
    else
    {
        /* Determine total number of steps remaining */
        if (bExtend)
        {
            ir->nsteps = ir->nsteps - (run_step - ir->init_step) + (gmx_int64_t)(extend_t/ir->delta_t + 0.5);
            printf("Extending remaining runtime of by %g ps (now %s steps)\n",
                   extend_t, gmx_step_str(ir->nsteps, buf));
        }
        else if (bUntil)
        {
            printf("nsteps = %s, run_step = %s, current_t = %g, until = %g\n",
                   gmx_step_str(ir->nsteps, buf),
                   gmx_step_str(run_step, buf2),
                   run_t, until_t);
            ir->nsteps = (gmx_int64_t)((until_t - run_t)/ir->delta_t + 0.5);
            printf("Extending remaining runtime until %g ps (now %s steps)\n",
                   until_t, gmx_step_str(ir->nsteps, buf));
        }
        else
        {
            ir->nsteps -= run_step - ir->init_step;
            /* Print message */
            printf("%s steps (%g ps) remaining from first run.\n",
                   gmx_step_str(ir->nsteps, buf), ir->nsteps*ir->delta_t);
        }
    }

    if (bNsteps || bZeroQ || (ir->nsteps > 0))
    {
        ir->init_step = run_step;

        if (ftp2bSet(efNDX, NFILE, fnm) ||
            !(bNsteps || bExtend || bUntil))
        {
            atoms = gmx_mtop_global_atoms(&mtop);
            get_index(&atoms, ftp2fn_null(efNDX, NFILE, fnm), 1,
                      &gnx, &index, &grpname);
            if (!bZeroQ)
            {
                bSel = (gnx != state.natoms);
                for (i = 0; ((i < gnx) && (!bSel)); i++)
                {
                    bSel = (i != index[i]);
                }
            }
            else
            {
                bSel = FALSE;
            }
            if (bSel)
            {
                fprintf(stderr, "Will write subset %s of original tpx containing %d "
                        "atoms\n", grpname, gnx);
                reduce_topology_x(gnx, index, &mtop, as_rvec_array(state.x.data()), as_rvec_array(state.v.data()));
                state.natoms = gnx;
            }
            else if (bZeroQ)
            {
                zeroq(index, &mtop);
                fprintf(stderr, "Zero-ing charges for group %s\n", grpname);
            }
            else
            {
                fprintf(stderr, "Will write full tpx file (no selection)\n");
            }
        }

        state_t = ir->init_t + ir->init_step*ir->delta_t;
        sprintf(buf,   "Writing statusfile with starting step %s%s and length %s%s steps...\n", "%10", GMX_PRId64, "%10", GMX_PRId64);
        fprintf(stderr, buf, ir->init_step, ir->nsteps);
        fprintf(stderr, "                                 time %10.3f and length %10.3f ps\n",
                state_t, ir->nsteps*ir->delta_t);
        write_tpx_state(opt2fn("-o", NFILE, fnm), ir, &state, &mtop);
    }
    else
    {
        printf("You've simulated long enough. Not writing tpr file\n");
    }

    return 0;
}
