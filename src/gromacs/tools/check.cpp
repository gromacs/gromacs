/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "check.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

typedef struct
{
    int bStep;
    int bTime;
    int bLambda;
    int bX;
    int bV;
    int bF;
    int bBox;
} t_count;

typedef struct
{
    float bStep;
    float bTime;
    float bLambda;
    float bX;
    float bV;
    float bF;
    float bBox;
} t_fr_time;

static void comp_tpx(const char* fn1, const char* fn2, gmx_bool bRMSD, real ftol, real abstol)
{
    const char*                 ff[2];
    std::unique_ptr<t_inputrec> ir[2];
    t_state                     state[2];
    gmx_mtop_t                  mtop[2];
    int                         i;

    ff[0] = fn1;
    ff[1] = fn2;
    for (i = 0; i < (fn2 ? 2 : 1); i++)
    {
        ir[i] = std::make_unique<t_inputrec>();
        read_tpx_state(ff[i], ir[i].get(), &state[i], &(mtop[i]));
        gmx::MDModules().adjustInputrecBasedOnModules(ir[i].get());
    }
    if (fn2)
    {
        cmp_inputrec(stdout, ir[0].get(), ir[1].get(), ftol, abstol);
        compareMtop(stdout, mtop[0], mtop[1], ftol, abstol);
        comp_state(&state[0], &state[1], bRMSD, ftol, abstol);
    }
    else
    {
        if (ir[0]->efep == FreeEnergyPerturbationType::No)
        {
            fprintf(stdout, "inputrec->efep = %s\n", enumValueToString(ir[0]->efep));
        }
        else
        {
            if (ir[0]->bPull)
            {
                comp_pull_AB(stdout, *ir[0]->pull, ftol, abstol);
            }
            compareMtopAB(stdout, mtop[0], ftol, abstol);
        }
    }
}

static void comp_trx(const gmx_output_env_t* oenv, const char* fn1, const char* fn2, gmx_bool bRMSD, real ftol, real abstol)
{
    int          i;
    const char*  fn[2];
    t_trxframe   fr[2];
    t_trxstatus* status[2];
    gmx_bool     b[2];

    fn[0] = fn1;
    fn[1] = fn2;
    fprintf(stderr, "Comparing trajectory files %s and %s\n", fn1, fn2);
    for (i = 0; i < 2; i++)
    {
        b[i] = read_first_frame(oenv, &status[i], fn[i], &fr[i], TRX_READ_X | TRX_READ_V | TRX_READ_F);
    }

    if (b[0] && b[1])
    {
        do
        {
            comp_frame(stdout, &(fr[0]), &(fr[1]), bRMSD, ftol, abstol);

            for (i = 0; i < 2; i++)
            {
                b[i] = read_next_frame(oenv, status[i], &fr[i]);
            }
        } while (b[0] && b[1]);

        for (i = 0; i < 2; i++)
        {
            if (b[i] && !b[1 - i])
            {
                fprintf(stdout, "\nEnd of file on %s but not on %s\n", fn[1 - i], fn[i]);
            }
            close_trx(status[i]);
        }
    }
    if (!b[0] && !b[1])
    {
        fprintf(stdout, "\nBoth files read correctly\n");
    }
}

static void chk_coords(int frame, int natoms, rvec* x, matrix box, real fac, real tol)
{
    int  i, j;
    int  nNul = 0;
    real vol  = det(box);

    for (i = 0; (i < natoms); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            if ((vol > 0) && (fabs(x[i][j]) > fac * box[j][j]))
            {
                printf("Warning at frame %d: coordinates for atom %d are large (%g)\n", frame, i, x[i][j]);
            }
        }
        if ((fabs(x[i][XX]) < tol) && (fabs(x[i][YY]) < tol) && (fabs(x[i][ZZ]) < tol))
        {
            nNul++;
        }
    }
    if (nNul > 0)
    {
        printf("Warning at frame %d: there are %d particles with all coordinates zero\n", frame, nNul);
    }
}

static void chk_vels(int frame, int natoms, rvec* v)
{
    int i, j;

    for (i = 0; (i < natoms); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            if (fabs(v[i][j]) > 500)
            {
                printf("Warning at frame %d. Velocities for atom %d are large (%g)\n", frame, i, v[i][j]);
            }
        }
    }
}

static void chk_forces(int frame, int natoms, rvec* f)
{
    int i, j;

    for (i = 0; (i < natoms); i++)
    {
        for (j = 0; (j < DIM); j++)
        {
            if (fabs(f[i][j]) > 10000)
            {
                printf("Warning at frame %d. Forces for atom %d are large (%g)\n", frame, i, f[i][j]);
            }
        }
    }
}

static void chk_bonds(const InteractionDefinitions* idef, PbcType pbcType, rvec* x, matrix box, real tol)
{
    int   ftype, k, ai, aj, type;
    real  b0, blen, deviation;
    t_pbc pbc;
    rvec  dx;

    gmx::ArrayRef<const t_iparams> iparams = idef->iparams;

    set_pbc(&pbc, pbcType, box);
    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if ((interaction_function[ftype].flags & IF_CHEMBOND) == IF_CHEMBOND)
        {
            for (k = 0; (k < idef->il[ftype].size());)
            {
                type = idef->il[ftype].iatoms[k++];
                ai   = idef->il[ftype].iatoms[k++];
                aj   = idef->il[ftype].iatoms[k++];
                b0   = 0;
                switch (ftype)
                {
                    case F_BONDS: b0 = iparams[type].harmonic.rA; break;
                    case F_G96BONDS: b0 = std::sqrt(iparams[type].harmonic.rA); break;
                    case F_MORSE: b0 = iparams[type].morse.b0A; break;
                    case F_CUBICBONDS: b0 = iparams[type].cubic.b0; break;
                    case F_CONSTR: b0 = iparams[type].constr.dA; break;
                    default: break;
                }
                if (b0 != 0)
                {
                    pbc_dx(&pbc, x[ai], x[aj], dx);
                    blen      = norm(dx);
                    deviation = gmx::square(blen - b0);
                    if (std::sqrt(deviation / gmx::square(b0)) > tol)
                    {
                        fprintf(stderr,
                                "Distance between atoms %d and %d is %.3f, should be %.3f\n",
                                ai + 1,
                                aj + 1,
                                blen,
                                b0);
                    }
                }
            }
        }
    }
}

static void chk_trj(const gmx_output_env_t* oenv, const char* fn, const char* tpr, real tol)
{
    t_trxframe   fr;
    t_count      count;
    t_fr_time    first, last;
    int          j = -1, new_natoms, natoms;
    real         old_t1, old_t2;
    gmx_bool     bShowTimestep = TRUE, newline = FALSE;
    t_trxstatus* status;
    gmx_mtop_t   mtop;
    t_state      state;
    t_inputrec   ir;

    std::unique_ptr<gmx_localtop_t> top;
    if (tpr)
    {
        read_tpx_state(tpr, &ir, &state, &mtop);
        top = std::make_unique<gmx_localtop_t>(mtop.ffparams);
        gmx_mtop_generate_local_top(mtop, top.get(), ir.efep != FreeEnergyPerturbationType::No);
    }
    new_natoms = -1;
    natoms     = -1;

    printf("Checking file %s\n", fn);

    j      = 0;
    old_t2 = -2.0;
    old_t1 = -1.0;

    count.bStep   = 0;
    count.bTime   = 0;
    count.bLambda = 0;
    count.bX      = 0;
    count.bV      = 0;
    count.bF      = 0;
    count.bBox    = 0;

    first.bStep   = 0;
    first.bTime   = 0;
    first.bLambda = 0;
    first.bX      = 0;
    first.bV      = 0;
    first.bF      = 0;
    first.bBox    = 0;

    last.bStep   = 0;
    last.bTime   = 0;
    last.bLambda = 0;
    last.bX      = 0;
    last.bV      = 0;
    last.bF      = 0;
    last.bBox    = 0;

    read_first_frame(oenv, &status, fn, &fr, TRX_READ_X | TRX_READ_V | TRX_READ_F);

    do
    {
        if (j == 0)
        {
            fprintf(stderr, "\n# Atoms  %d\n", fr.natoms);
            if (fr.bPrec)
            {
                fprintf(stderr, "Precision %g (nm)\n", 1 / fr.prec);
            }
        }
        newline = TRUE;
        if ((natoms > 0) && (new_natoms != natoms))
        {
            fprintf(stderr, "\nNumber of atoms at t=%g don't match (%d, %d)\n", old_t1, natoms, new_natoms);
            newline = FALSE;
        }
        if (j >= 2)
        {
            if (std::fabs((fr.time - old_t1) - (old_t1 - old_t2))
                > 0.1 * (std::fabs(fr.time - old_t1) + std::fabs(old_t1 - old_t2)))
            {
                bShowTimestep = FALSE;
                fprintf(stderr,
                        "%sTimesteps at t=%g don't match (%g, %g)\n",
                        newline ? "\n" : "",
                        old_t1,
                        old_t1 - old_t2,
                        fr.time - old_t1);
            }
        }
        natoms = new_natoms;
        if (tpr)
        {
            chk_bonds(&top->idef, ir.pbcType, fr.x, fr.box, tol);
        }
        if (fr.bX)
        {
            chk_coords(j, natoms, fr.x, fr.box, 1e5, tol);
        }
        if (fr.bV)
        {
            chk_vels(j, natoms, fr.v);
        }
        if (fr.bF)
        {
            chk_forces(j, natoms, fr.f);
        }

        old_t2 = old_t1;
        old_t1 = fr.time;
        j++;
        new_natoms = fr.natoms;
#define INC(s, n, f, l, item)     \
    if ((s).item != 0)            \
    {                             \
        if ((n).item == 0)        \
        {                         \
            first.item = fr.time; \
        }                         \
        last.item = fr.time;      \
        (n).item++;               \
    }
        INC(fr, count, first, last, bStep)
        INC(fr, count, first, last, bTime)
        INC(fr, count, first, last, bLambda)
        INC(fr, count, first, last, bX)
        INC(fr, count, first, last, bV)
        INC(fr, count, first, last, bF)
        INC(fr, count, first, last, bBox)
#undef INC
    } while (read_next_frame(oenv, status, &fr));

    fprintf(stderr, "\n");

    close_trx(status);

    fprintf(stderr, "\nItem        #frames");
    if (bShowTimestep)
    {
        fprintf(stderr, " Timestep (ps)");
    }
    fprintf(stderr, "\n");
#define PRINTITEM(label, item)                                                    \
    fprintf(stderr, "%-10s  %6d", label, count.item);                             \
    if ((bShowTimestep) && (count.item > 1))                                      \
    {                                                                             \
        fprintf(stderr, "    %g\n", (last.item - first.item) / (count.item - 1)); \
    }                                                                             \
    else                                                                          \
        fprintf(stderr, "\n")
    PRINTITEM("Step", bStep);
    PRINTITEM("Time", bTime);
    PRINTITEM("Lambda", bLambda);
    PRINTITEM("Coords", bX);
    PRINTITEM("Velocities", bV);
    PRINTITEM("Forces", bF);
    PRINTITEM("Box", bBox);
}

static void chk_tps(const char* fn, real vdw_fac, real bon_lo, real bon_hi)
{
    int        natom, i, j, k;
    t_topology top;
    PbcType    pbcType;
    t_atoms*   atoms;
    rvec *     x, *v;
    rvec       dx;
    matrix     box;
    t_pbc      pbc;
    gmx_bool   bV, bX, bB, bFirst, bOut;
    real       r2, ekin, temp1, temp2, dist2, vdwfac2, bonlo2, bonhi2;
    real*      atom_vdw;

    fprintf(stderr, "Checking coordinate file %s\n", fn);
    read_tps_conf(fn, &top, &pbcType, &x, &v, box, TRUE);
    atoms = &top.atoms;
    natom = atoms->nr;
    fprintf(stderr, "%d atoms in file\n", atoms->nr);

    /* check coordinates and box */
    bV = FALSE;
    bX = FALSE;
    for (i = 0; (i < natom) && !(bV && bX); i++)
    {
        for (j = 0; (j < DIM) && !(bV && bX); j++)
        {
            bV = bV || (v[i][j] != 0);
            bX = bX || (x[i][j] != 0);
        }
    }
    bB = FALSE;
    for (i = 0; (i < DIM) && !bB; i++)
    {
        for (j = 0; (j < DIM) && !bB; j++)
        {
            bB = bB || (box[i][j] != 0);
        }
    }

    fprintf(stderr, "coordinates %s\n", bX ? "found" : "absent");
    fprintf(stderr, "box         %s\n", bB ? "found" : "absent");
    fprintf(stderr, "velocities  %s\n", bV ? "found" : "absent");
    fprintf(stderr, "\n");

    /* check velocities */
    if (bV)
    {
        ekin = 0.0;
        for (i = 0; (i < natom); i++)
        {
            for (j = 0; (j < DIM); j++)
            {
                ekin += 0.5 * atoms->atom[i].m * v[i][j] * v[i][j];
            }
        }
        temp1 = (2.0 * ekin) / (natom * DIM * gmx::c_boltz);
        temp2 = (2.0 * ekin) / (natom * (DIM - 1) * gmx::c_boltz);
        fprintf(stderr, "Kinetic energy: %g (kJ/mol)\n", ekin);
        fprintf(stderr,
                "Assuming the number of degrees of freedom to be "
                "Natoms * %d or Natoms * %d,\n"
                "the velocities correspond to a temperature of the system\n"
                "of %g K or %g K respectively.\n\n",
                DIM,
                DIM - 1,
                temp1,
                temp2);
    }

    /* check coordinates */
    if (bX)
    {
        vdwfac2 = gmx::square(vdw_fac);
        bonlo2  = gmx::square(bon_lo);
        bonhi2  = gmx::square(bon_hi);

        fprintf(stderr,
                "Checking for atoms closer than %g and not between %g and %g,\n"
                "relative to sum of Van der Waals distance:\n",
                vdw_fac,
                bon_lo,
                bon_hi);
        snew(atom_vdw, natom);
        AtomProperties aps;
        for (i = 0; (i < natom); i++)
        {
            aps.setAtomProperty(epropVDW,
                                *(atoms->resinfo[atoms->atom[i].resind].name),
                                *(atoms->atomname[i]),
                                &(atom_vdw[i]));
            if (debug)
            {
                fprintf(debug,
                        "%5d %4s %4s %7g\n",
                        i + 1,
                        *(atoms->resinfo[atoms->atom[i].resind].name),
                        *(atoms->atomname[i]),
                        atom_vdw[i]);
            }
        }
        if (bB)
        {
            set_pbc(&pbc, pbcType, box);
        }

        bFirst = TRUE;
        for (i = 0; (i < natom); i++)
        {
            if (((i + 1) % 10) == 0)
            {
                fprintf(stderr, "\r%5d", i + 1);
                fflush(stderr);
            }
            for (j = i + 1; (j < natom); j++)
            {
                if (bB)
                {
                    pbc_dx(&pbc, x[i], x[j], dx);
                }
                else
                {
                    rvec_sub(x[i], x[j], dx);
                }
                r2    = iprod(dx, dx);
                dist2 = gmx::square(atom_vdw[i] + atom_vdw[j]);
                if ((r2 <= dist2 * bonlo2) || ((r2 >= dist2 * bonhi2) && (r2 <= dist2 * vdwfac2)))
                {
                    if (bFirst)
                    {
                        fprintf(stderr,
                                "\r%5s %4s %8s %5s  %5s %4s %8s %5s  %6s\n",
                                "atom#",
                                "name",
                                "residue",
                                "r_vdw",
                                "atom#",
                                "name",
                                "residue",
                                "r_vdw",
                                "distance");
                        bFirst = FALSE;
                    }
                    fprintf(stderr,
                            "\r%5d %4s %4s%4d %-5.3g  %5d %4s %4s%4d %-5.3g  %-6.4g\n",
                            i + 1,
                            *(atoms->atomname[i]),
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            atoms->resinfo[atoms->atom[i].resind].nr,
                            atom_vdw[i],
                            j + 1,
                            *(atoms->atomname[j]),
                            *(atoms->resinfo[atoms->atom[j].resind].name),
                            atoms->resinfo[atoms->atom[j].resind].nr,
                            atom_vdw[j],
                            std::sqrt(r2));
                }
            }
        }
        if (bFirst)
        {
            fprintf(stderr, "\rno close atoms found\n");
        }
        fprintf(stderr, "\r      \n");

        if (bB)
        {
            /* check box */
            bFirst = TRUE;
            k      = 0;
            for (i = 0; (i < natom) && (k < 10); i++)
            {
                bOut = FALSE;
                for (j = 0; (j < DIM) && !bOut; j++)
                {
                    bOut = bOut || (x[i][j] < 0) || (x[i][j] > box[j][j]);
                }
                if (bOut)
                {
                    k++;
                    if (bFirst)
                    {
                        fprintf(stderr, "Atoms outside box ( ");
                        for (j = 0; (j < DIM); j++)
                        {
                            fprintf(stderr, "%g ", box[j][j]);
                        }
                        fprintf(stderr,
                                "):\n"
                                "(These may occur often and are normally not a problem)\n"
                                "%5s %4s %8s %5s  %s\n",
                                "atom#",
                                "name",
                                "residue",
                                "r_vdw",
                                "coordinate");
                        bFirst = FALSE;
                    }
                    fprintf(stderr,
                            "%5d %4s %4s%4d %-5.3g",
                            i,
                            *(atoms->atomname[i]),
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            atoms->resinfo[atoms->atom[i].resind].nr,
                            atom_vdw[i]);
                    for (j = 0; (j < DIM); j++)
                    {
                        fprintf(stderr, " %6.3g", x[i][j]);
                    }
                    fprintf(stderr, "\n");
                }
            }
            if (k == 10)
            {
                fprintf(stderr, "(maybe more)\n");
            }
            if (bFirst)
            {
                fprintf(stderr, "no atoms found outside box\n");
            }
            fprintf(stderr, "\n");
        }
    }
}

static void chk_ndx(const char* fn)
{
    const auto grps = init_index(fn);
    if (debug)
    {
        pr_blocka(debug, 0, fn, grps, false);
    }
    else
    {
        printf("Contents of index file %s\n", fn);
        printf("--------------------------------------------------\n");
        printf("Nr.   Group               #Entries   First    Last\n");
        for (gmx::Index i = 0; i < gmx::ssize(grps); i++)
        {
            printf("%4td  %-20s%8td%8d%8d\n",
                   i,
                   grps[i].name.c_str(),
                   gmx::ssize(grps[i].particleIndices),
                   grps[i].particleIndices[0] + 1,
                   grps[i].particleIndices.back() + 1);
        }
    }
}

static void chk_enx(const char* fn)
{
    int          nre, fnr;
    ener_file_t  in;
    gmx_enxnm_t* enm = nullptr;
    t_enxframe*  fr;
    gmx_bool     bShowTStep;
    gmx_bool     timeSet;
    real         t0, old_t1, old_t2;
    char         buf[22];

    fprintf(stderr, "Checking energy file %s\n\n", fn);

    in = open_enx(fn, "r");
    do_enxnms(in, &nre, &enm);
    fprintf(stderr, "%d groups in energy file", nre);
    snew(fr, 1);
    old_t2     = -2.0;
    old_t1     = -1.0;
    fnr        = 0;
    t0         = 0;
    timeSet    = FALSE;
    bShowTStep = TRUE;

    while (do_enx(in, fr))
    {
        if (fnr >= 2)
        {
            if (fabs((fr->t - old_t1) - (old_t1 - old_t2))
                > 0.1 * (fabs(fr->t - old_t1) + std::fabs(old_t1 - old_t2)))
            {
                bShowTStep = FALSE;
                fprintf(stderr, "\nTimesteps at t=%g don't match (%g, %g)\n", old_t1, old_t1 - old_t2, fr->t - old_t1);
            }
        }
        old_t2 = old_t1;
        old_t1 = fr->t;
        if (!timeSet)
        {
            t0      = fr->t;
            timeSet = TRUE;
        }
        if (fnr == 0)
        {
            fprintf(stderr, "\rframe: %6s (index %6d), t: %10.3f\n", gmx_step_str(fr->step, buf), fnr, fr->t);
        }
        fnr++;
    }
    fprintf(stderr, "\n\nFound %d frames", fnr);
    if (bShowTStep && fnr > 1)
    {
        fprintf(stderr, " with a timestep of %g ps", (old_t1 - t0) / (fnr - 1));
    }
    fprintf(stderr, ".\n");

    free_enxframe(fr);
    free_enxnms(nre, enm);
    sfree(fr);
}

int gmx_check(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] reads a trajectory ([REF].tng[ref], [REF].trr[ref] or ",
        "[REF].xtc[ref]), an energy file ([REF].edr[ref])",
        "or an index file ([REF].ndx[ref])",
        "and prints out useful information about them.[PAR]",
        "Option [TT]-c[tt] checks for presence of coordinates,",
        "velocities and box in the file, for close contacts (smaller than",
        "[TT]-vdwfac[tt] and not bonded, i.e. not between [TT]-bonlo[tt]",
        "and [TT]-bonhi[tt], all relative to the sum of both Van der Waals",
        "radii) and atoms outside the box (these may occur often and are",
        "no problem). If velocities are present, an estimated temperature",
        "will be calculated from them.[PAR]",
        "If an index file, is given its contents will be summarized.[PAR]",
        "If both a trajectory and a [REF].tpr[ref] file are given (with [TT]-s1[tt])",
        "the program will check whether the bond lengths defined in the tpr",
        "file are indeed correct in the trajectory. If not you may have",
        "non-matching files due to e.g. deshuffling or due to problems with",
        "virtual sites. With these flags, [TT]gmx check[tt] provides a quick check for ",
        "such problems.[PAR]",
        "The program can compare two run input ([REF].tpr[ref])",
        "files",
        "when both [TT]-s1[tt] and [TT]-s2[tt] are supplied. When comparing",
        "run input files this way, the default relative tolerance is reduced",
        "to 0.000001 and the absolute tolerance set to zero to find any differences",
        "not due to minor compiler optimization differences, although you can",
        "of course still set any other tolerances through the options.",
        "Similarly a pair of trajectory files can be compared (using the [TT]-f2[tt]",
        "option), or a pair of energy files (using the [TT]-e2[tt] option).[PAR]",
        "For free energy simulations the A and B state topology from one",
        "run input file can be compared with options [TT]-s1[tt] and [TT]-ab[tt].[PAR]"
    };
    t_filenm fnm[] = { { efTRX, "-f", nullptr, ffOPTRD },  { efTRX, "-f2", nullptr, ffOPTRD },
                       { efTPR, "-s1", "top1", ffOPTRD },  { efTPR, "-s2", "top2", ffOPTRD },
                       { efTPS, "-c", nullptr, ffOPTRD },  { efEDR, "-e", nullptr, ffOPTRD },
                       { efEDR, "-e2", "ener2", ffOPTRD }, { efNDX, "-n", nullptr, ffOPTRD },
                       { efTEX, "-m", nullptr, ffOPTWR } };
#define NFILE asize(fnm)
    const char *fn1 = nullptr, *fn2 = nullptr, *tex = nullptr;

    gmx_output_env_t* oenv;
    real              vdw_fac  = 0.8;
    real              bon_lo   = 0.4;
    real              bon_hi   = 0.7;
    gmx_bool          bRMSD    = FALSE;
    real              ftol     = 0.001;
    real              abstol   = 0.001;
    gmx_bool          bCompAB  = FALSE;
    char*             lastener = nullptr;
    t_pargs           pa[]     = {
        { "-vdwfac",
          FALSE,
          etREAL,
          { &vdw_fac },
          "Fraction of sum of VdW radii used as warning cutoff" },
        { "-bonlo", FALSE, etREAL, { &bon_lo }, "Min. fract. of sum of VdW radii for bonded atoms" },
        { "-bonhi", FALSE, etREAL, { &bon_hi }, "Max. fract. of sum of VdW radii for bonded atoms" },
        { "-rmsd", FALSE, etBOOL, { &bRMSD }, "Print RMSD for x, v and f" },
        { "-tol",
          FALSE,
          etREAL,
          { &ftol },
          "Relative tolerance for comparing real values defined as "
          "[MATH]2*(a-b)/([MAG]a[mag]+[MAG]b[mag])[math]" },
        { "-abstol",
          FALSE,
          etREAL,
          { &abstol },
          "Absolute tolerance, useful when sums are close to zero." },
        { "-ab", FALSE, etBOOL, { &bCompAB }, "Compare the A and B topology from one file" },
        { "-lastener",
          FALSE,
          etSTR,
          { &lastener },
          "Last energy term to compare (if not given all are tested). It makes sense to go up "
          "until the Pressure." }
    };

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    fn1 = opt2fn_null("-f", NFILE, fnm);
    fn2 = opt2fn_null("-f2", NFILE, fnm);
    tex = opt2fn_null("-m", NFILE, fnm);

    if (tex)
    {
        fprintf(stderr,
                "LaTeX file writing has been removed from gmx check. "
                "Please use gmx report-methods instead for it.\n");
    }
    if (fn1 && fn2)
    {
        comp_trx(oenv, fn1, fn2, bRMSD, ftol, abstol);
    }
    else if (fn1)
    {
        chk_trj(oenv, fn1, opt2fn_null("-s1", NFILE, fnm), ftol);
    }
    else if (fn2)
    {
        fprintf(stderr, "Please give me TWO trajectory (.xtc/.trr/.tng) files!\n");
    }
    output_env_done(oenv);

    fn1 = opt2fn_null("-s1", NFILE, fnm);
    fn2 = opt2fn_null("-s2", NFILE, fnm);
    if ((fn1 && fn2) || bCompAB)
    {
        if (bCompAB)
        {
            if (fn1 == nullptr)
            {
                gmx_fatal(FARGS, "With -ab you need to set the -s1 option");
            }
            fn2 = nullptr;
        }

        fprintf(stderr, "Note: When comparing run input files, default tolerances are reduced.\n");
        if (!opt2parg_bSet("-tol", asize(pa), pa))
        {
            ftol = 0.000001;
        }
        if (!opt2parg_bSet("-abstol", asize(pa), pa))
        {
            abstol = 0;
        }
        comp_tpx(fn1, fn2, bRMSD, ftol, abstol);
    }
    else if ((fn1 && !opt2fn_null("-f", NFILE, fnm)) || (!fn1 && fn2))
    {
        fprintf(stderr, "Please give me TWO run input (.tpr) files\n");
    }

    fn1 = opt2fn_null("-e", NFILE, fnm);
    fn2 = opt2fn_null("-e2", NFILE, fnm);
    if (fn1 && fn2)
    {
        comp_enx(fn1, fn2, ftol, abstol, lastener);
    }
    else if (fn1)
    {
        chk_enx(ftp2fn(efEDR, NFILE, fnm));
    }
    else if (fn2)
    {
        fprintf(stderr, "Please give me TWO energy (.edr/.ene) files!\n");
    }

    if (ftp2bSet(efTPS, NFILE, fnm))
    {
        chk_tps(ftp2fn(efTPS, NFILE, fnm), vdw_fac, bon_lo, bon_hi);
    }

    if (ftp2bSet(efNDX, NFILE, fnm))
    {
        chk_ndx(ftp2fn(efNDX, NFILE, fnm));
    }

    return 0;
}
