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

#include "gromacs/topology/idef.h"

#include <cstdio>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/txtdump.h"

static void printHarmonicInteraction(gmx::TextWriter* writer,
                                     const t_iparams& iparams,
                                     const char*      r,
                                     const char*      kr)
{
    writer->writeLineFormatted("%sA=%12.5e, %sA=%12.5e, %sB=%12.5e, %sB=%12.5e",
                               r,
                               iparams.harmonic.rA,
                               kr,
                               iparams.harmonic.krA,
                               r,
                               iparams.harmonic.rB,
                               kr,
                               iparams.harmonic.krB);
}

void pr_iparams(FILE* fp, t_functype ftype, const t_iparams& iparams)
{
    gmx::StringOutputStream stream;
    {
        gmx::TextWriter writer(&stream);
        printInteractionParameters(&writer, ftype, iparams);
    }
    fputs(stream.toString().c_str(), fp);
}

void printInteractionParameters(gmx::TextWriter* writer, t_functype ftype, const t_iparams& iparams)
{
    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES: printHarmonicInteraction(writer, iparams, "th", "ct"); break;
        case F_CROSS_BOND_BONDS:
            writer->writeLineFormatted("r1e=%15.8e, r2e=%15.8e, krr=%15.8e",
                                       iparams.cross_bb.r1e,
                                       iparams.cross_bb.r2e,
                                       iparams.cross_bb.krr);
            break;
        case F_CROSS_BOND_ANGLES:
            writer->writeLineFormatted("r1e=%15.8e, r1e=%15.8e, r3e=%15.8e, krt=%15.8e",
                                       iparams.cross_ba.r1e,
                                       iparams.cross_ba.r2e,
                                       iparams.cross_ba.r3e,
                                       iparams.cross_ba.krt);
            break;
        case F_LINEAR_ANGLES:
            writer->writeLineFormatted("klinA=%15.8e, aA=%15.8e, klinB=%15.8e, aB=%15.8e",
                                       iparams.linangle.klinA,
                                       iparams.linangle.aA,
                                       iparams.linangle.klinB,
                                       iparams.linangle.aB);
            break;
        case F_UREY_BRADLEY:
            writer->writeLineFormatted(
                    "thetaA=%15.8e, kthetaA=%15.8e, r13A=%15.8e, kUBA=%15.8e, thetaB=%15.8e, "
                    "kthetaB=%15.8e, r13B=%15.8e, kUBB=%15.8e",
                    iparams.u_b.thetaA,
                    iparams.u_b.kthetaA,
                    iparams.u_b.r13A,
                    iparams.u_b.kUBA,
                    iparams.u_b.thetaB,
                    iparams.u_b.kthetaB,
                    iparams.u_b.r13B,
                    iparams.u_b.kUBB);
            break;
        case F_QUARTIC_ANGLES:
            writer->writeStringFormatted("theta=%15.8e", iparams.qangle.theta);
            for (int i = 0; i < 5; i++)
            {
                writer->writeStringFormatted(", c%c=%15.8e", '0' + i, iparams.qangle.c[i]);
            }
            writer->ensureLineBreak();
            break;
        case F_BHAM:
            writer->writeLineFormatted(
                    "a=%15.8e, b=%15.8e, c=%15.8e", iparams.bham.a, iparams.bham.b, iparams.bham.c);
            break;
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC: printHarmonicInteraction(writer, iparams, "b0", "cb"); break;
        case F_IDIHS: printHarmonicInteraction(writer, iparams, "xi", "cx"); break;
        case F_MORSE:
            writer->writeLineFormatted(
                    "b0A=%15.8e, cbA=%15.8e, betaA=%15.8e, b0B=%15.8e, cbB=%15.8e, betaB=%15.8e",
                    iparams.morse.b0A,
                    iparams.morse.cbA,
                    iparams.morse.betaA,
                    iparams.morse.b0B,
                    iparams.morse.cbB,
                    iparams.morse.betaB);
            break;
        case F_CUBICBONDS:
            writer->writeLineFormatted("b0=%15.8e, kb=%15.8e, kcub=%15.8e",
                                       iparams.cubic.b0,
                                       iparams.cubic.kb,
                                       iparams.cubic.kcub);
            break;
        case F_CONNBONDS: writer->ensureEmptyLine(); break;
        case F_FENEBONDS:
            writer->writeLineFormatted("bm=%15.8e, kb=%15.8e", iparams.fene.bm, iparams.fene.kb);
            break;
        case F_RESTRBONDS:
            writer->writeLineFormatted(
                    "lowA=%15.8e, up1A=%15.8e, up2A=%15.8e, kA=%15.8e, lowB=%15.8e, up1B=%15.8e, "
                    "up2B=%15.8e, kB=%15.8e,",
                    iparams.restraint.lowA,
                    iparams.restraint.up1A,
                    iparams.restraint.up2A,
                    iparams.restraint.kA,
                    iparams.restraint.lowB,
                    iparams.restraint.up1B,
                    iparams.restraint.up2B,
                    iparams.restraint.kB);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            writer->writeLineFormatted(
                    "tab=%d, kA=%15.8e, kB=%15.8e", iparams.tab.table, iparams.tab.kA, iparams.tab.kB);
            break;
        case F_POLARIZATION:
            writer->writeLineFormatted("alpha=%15.8e", iparams.polarize.alpha);
            break;
        case F_ANHARM_POL:
            writer->writeLineFormatted("alpha=%15.8e drcut=%15.8e khyp=%15.8e",
                                       iparams.anharm_polarize.alpha,
                                       iparams.anharm_polarize.drcut,
                                       iparams.anharm_polarize.khyp);
            break;
        case F_THOLE_POL:
            writer->writeLineFormatted("a=%15.8e, alpha1=%15.8e, alpha2=%15.8e",
                                       iparams.thole.a,
                                       iparams.thole.alpha1,
                                       iparams.thole.alpha2);
            break;
        case F_WATER_POL:
            writer->writeLineFormatted(
                    "al_x=%15.8e, al_y=%15.8e, al_z=%15.8e, rOH=%9.6f, rHH=%9.6f, rOD=%9.6f",
                    iparams.wpol.al_x,
                    iparams.wpol.al_y,
                    iparams.wpol.al_z,
                    iparams.wpol.rOH,
                    iparams.wpol.rHH,
                    iparams.wpol.rOD);
            break;
        case F_LJ:
            writer->writeLineFormatted("c6=%15.8e, c12=%15.8e", iparams.lj.c6, iparams.lj.c12);
            break;
        case F_LJ14:
            writer->writeLineFormatted("c6A=%15.8e, c12A=%15.8e, c6B=%15.8e, c12B=%15.8e",
                                       iparams.lj14.c6A,
                                       iparams.lj14.c12A,
                                       iparams.lj14.c6B,
                                       iparams.lj14.c12B);
            break;
        case F_LJC14_Q:
            writer->writeLineFormatted("fqq=%15.8e, qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e",
                                       iparams.ljc14.fqq,
                                       iparams.ljc14.qi,
                                       iparams.ljc14.qj,
                                       iparams.ljc14.c6,
                                       iparams.ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            writer->writeLineFormatted("qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e",
                                       iparams.ljcnb.qi,
                                       iparams.ljcnb.qj,
                                       iparams.ljcnb.c6,
                                       iparams.ljcnb.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            writer->writeLineFormatted("phiA=%15.8e, cpA=%15.8e, phiB=%15.8e, cpB=%15.8e, mult=%d",
                                       iparams.pdihs.phiA,
                                       iparams.pdihs.cpA,
                                       iparams.pdihs.phiB,
                                       iparams.pdihs.cpB,
                                       iparams.pdihs.mult);
            break;
        case F_DISRES:
            writer->writeLineFormatted(
                    "label=%4d, type=%1d, low=%15.8e, up1=%15.8e, up2=%15.8e, fac=%15.8e)",
                    iparams.disres.label,
                    iparams.disres.type,
                    iparams.disres.low,
                    iparams.disres.up1,
                    iparams.disres.up2,
                    iparams.disres.kfac);
            break;
        case F_ORIRES:
            writer->writeLineFormatted(
                    "ex=%4d, label=%d, power=%4d, c=%15.8e, obs=%15.8e, kfac=%15.8e)",
                    iparams.orires.ex,
                    iparams.orires.label,
                    iparams.orires.power,
                    iparams.orires.c,
                    iparams.orires.obs,
                    iparams.orires.kfac);
            break;
        case F_DIHRES:
            writer->writeLineFormatted(
                    "phiA=%15.8e, dphiA=%15.8e, kfacA=%15.8e, phiB=%15.8e, dphiB=%15.8e, "
                    "kfacB=%15.8e",
                    iparams.dihres.phiA,
                    iparams.dihres.dphiA,
                    iparams.dihres.kfacA,
                    iparams.dihres.phiB,
                    iparams.dihres.dphiB,
                    iparams.dihres.kfacB);
            break;
        case F_POSRES:
            writer->writeLineFormatted(
                    "pos0A=(%15.8e,%15.8e,%15.8e), fcA=(%15.8e,%15.8e,%15.8e), "
                    "pos0B=(%15.8e,%15.8e,%15.8e), fcB=(%15.8e,%15.8e,%15.8e)",
                    iparams.posres.pos0A[XX],
                    iparams.posres.pos0A[YY],
                    iparams.posres.pos0A[ZZ],
                    iparams.posres.fcA[XX],
                    iparams.posres.fcA[YY],
                    iparams.posres.fcA[ZZ],
                    iparams.posres.pos0B[XX],
                    iparams.posres.pos0B[YY],
                    iparams.posres.pos0B[ZZ],
                    iparams.posres.fcB[XX],
                    iparams.posres.fcB[YY],
                    iparams.posres.fcB[ZZ]);
            break;
        case F_FBPOSRES:
            writer->writeLineFormatted(
                    "pos0=(%15.8e,%15.8e,%15.8e), geometry=%d, r=%15.8e, k=%15.8e",
                    iparams.fbposres.pos0[XX],
                    iparams.fbposres.pos0[YY],
                    iparams.fbposres.pos0[ZZ],
                    iparams.fbposres.geom,
                    iparams.fbposres.r,
                    iparams.fbposres.k);
            break;
        case F_RBDIHS:
            for (int i = 0; i < NR_RBDIHS; i++)
            {
                writer->writeStringFormatted(
                        "%srbcA[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams.rbdihs.rbcA[i]);
            }
            writer->ensureLineBreak();
            for (int i = 0; i < NR_RBDIHS; i++)
            {
                writer->writeStringFormatted(
                        "%srbcB[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams.rbdihs.rbcB[i]);
            }
            writer->ensureLineBreak();
            break;
        case F_FOURDIHS:
        {
            /* Use the OPLS -> Ryckaert-Bellemans formula backwards to get
             * the OPLS potential constants back.
             */
            const real* rbcA = iparams.rbdihs.rbcA;
            const real* rbcB = iparams.rbdihs.rbcB;
            real        VA[4], VB[4];

            VA[3] = -0.25 * rbcA[4];
            VA[2] = -0.5 * rbcA[3];
            VA[1] = 4.0 * VA[3] - rbcA[2];
            VA[0] = 3.0 * VA[2] - 2.0 * rbcA[1];

            VB[3] = -0.25 * rbcB[4];
            VB[2] = -0.5 * rbcB[3];
            VB[1] = 4.0 * VB[3] - rbcB[2];
            VB[0] = 3.0 * VB[2] - 2.0 * rbcB[1];

            for (int i = 0; i < NR_FOURDIHS; i++)
            {
                writer->writeStringFormatted("%sFourA[%d]=%15.8e", i == 0 ? "" : ", ", i, VA[i]);
            }
            writer->ensureLineBreak();
            for (int i = 0; i < NR_FOURDIHS; i++)
            {
                writer->writeStringFormatted("%sFourB[%d]=%15.8e", i == 0 ? "" : ", ", i, VB[i]);
            }
            writer->ensureLineBreak();
            break;
        }

        case F_CONSTR:
        case F_CONSTRNC:
            writer->writeLineFormatted("dA=%15.8e, dB=%15.8e", iparams.constr.dA, iparams.constr.dB);
            break;
        case F_SETTLE:
            writer->writeLineFormatted("doh=%15.8e, dhh=%15.8e", iparams.settle.doh, iparams.settle.dhh);
            break;
        case F_VSITE1: writer->ensureEmptyLine(); break;
        case F_VSITE2:
        case F_VSITE2FD: writer->writeLineFormatted("a=%15.8e", iparams.vsite.a); break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            writer->writeLineFormatted("a=%15.8e, b=%15.8e", iparams.vsite.a, iparams.vsite.b);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            writer->writeLineFormatted(
                    "a=%15.8e, b=%15.8e, c=%15.8e", iparams.vsite.a, iparams.vsite.b, iparams.vsite.c);
            break;
        case F_VSITEN:
            writer->writeLineFormatted("n=%2d, a=%15.8e", iparams.vsiten.n, iparams.vsiten.a);
            break;
        case F_GB12_NOLONGERUSED:
        case F_GB13_NOLONGERUSED:
        case F_GB14_NOLONGERUSED:
            // These could only be generated by grompp, not written in
            // a .top file. Now that implicit solvent is not
            // supported, they can't be generated, and the values are
            // ignored if read from an old .tpr file. So there is
            // nothing to print.
            break;
        case F_CMAP:
            writer->writeLineFormatted("cmapA=%1d, cmapB=%1d", iparams.cmap.cmapA, iparams.cmap.cmapB);
            break;
        case F_RESTRANGLES: printHarmonicInteraction(writer, iparams, "costheta0", "ktheta"); break;
        case F_RESTRDIHS:
            writer->writeLineFormatted("phiA=%15.8e, cpA=%15.8e", iparams.pdihs.phiA, iparams.pdihs.cpA);
            break;
        case F_CBTDIHS:
            writer->writeLineFormatted("kphi=%15.8e", iparams.cbtdihs.cbtcA[0]);
            for (int i = 1; i < NR_CBTDIHS; i++)
            {
                writer->writeStringFormatted(", cbtcA[%d]=%15.8e", i - 1, iparams.cbtdihs.cbtcA[i]);
            }
            writer->ensureLineBreak();
            break;
        default:
            gmx_fatal(FARGS,
                      "unknown function type %d (%s) in %s line %d",
                      ftype,
                      interaction_function[ftype].name,
                      __FILE__,
                      __LINE__);
    }
}

template<typename T>
static void printIlist(FILE*             fp,
                       int               indent,
                       const char*       title,
                       const t_functype* functype,
                       const T&          ilist,
                       gmx_bool          bShowNumbers,
                       gmx_bool          bShowParameters,
                       const t_iparams*  iparams)
{
    indent = pr_title(fp, indent, title);
    pr_indent(fp, indent);
    fprintf(fp, "nr: %d\n", ilist.size());
    if (!ilist.empty())
    {
        pr_indent(fp, indent);
        fprintf(fp, "iatoms:\n");
        int j = 0;
        for (int i = 0; i < ilist.size();)
        {
            pr_indent(fp, indent + INDENT);
            const int type  = ilist.iatoms[i];
            const int ftype = functype[type];
            if (bShowNumbers)
            {
                fprintf(fp, "%d type=%d ", j, type);
            }
            j++;
            printf("(%s)", interaction_function[ftype].name);
            for (int k = 0; k < interaction_function[ftype].nratoms; k++)
            {
                fprintf(fp, " %3d", ilist.iatoms[i + 1 + k]);
            }
            if (bShowParameters)
            {
                fprintf(fp, "  ");
                pr_iparams(fp, ftype, iparams[type]);
            }
            fprintf(fp, "\n");
            i += 1 + interaction_function[ftype].nratoms;
        }
    }
}

void pr_ilist(FILE*                  fp,
              int                    indent,
              const char*            title,
              const t_functype*      functype,
              const InteractionList& ilist,
              gmx_bool               bShowNumbers,
              gmx_bool               bShowParameters,
              const t_iparams*       iparams)
{
    printIlist(fp, indent, title, functype, ilist, bShowNumbers, bShowParameters, iparams);
}

void pr_idef(FILE* fp, int indent, const char* title, const t_idef* idef, gmx_bool bShowNumbers, gmx_bool bShowParameters)
{
    if (available(fp, idef, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "atnr=%d\n", idef->atnr);
        pr_indent(fp, indent);
        fprintf(fp, "ntypes=%d\n", idef->ntypes);
        for (int i = 0; i < idef->ntypes; i++)
        {
            pr_indent(fp, indent + INDENT);
            fprintf(fp,
                    "functype[%d]=%s, ",
                    bShowNumbers ? i : -1,
                    interaction_function[idef->functype[i]].name);
            pr_iparams(fp, idef->functype[i], idef->iparams[i]);
        }
        pr_real(fp, indent, "fudgeQQ", idef->fudgeQQ);

        for (int j = 0; (j < F_NRE); j++)
        {
            printIlist(fp,
                       indent,
                       interaction_function[j].longname,
                       idef->functype,
                       idef->il[j],
                       bShowNumbers,
                       bShowParameters,
                       idef->iparams);
        }
    }
}

void init_idef(t_idef* idef)
{
    idef->ntypes           = 0;
    idef->atnr             = 0;
    idef->functype         = nullptr;
    idef->iparams          = nullptr;
    idef->fudgeQQ          = 0.0;
    idef->iparams_posres   = nullptr;
    idef->iparams_fbposres = nullptr;
    for (int f = 0; f < F_NRE; ++f)
    {
        idef->il[f].iatoms = nullptr;
        idef->il[f].nalloc = 0;
        idef->il[f].nr     = 0;
    }
}

InteractionDefinitions::InteractionDefinitions(const gmx_ffparams_t& ffparams) :
    iparams(ffparams.iparams), functype(ffparams.functype), cmap_grid(ffparams.cmap_grid)
{
}

void InteractionDefinitions::clear()
{
    /* Clear the counts */
    for (auto& ilist : il)
    {
        ilist.clear();
    }
    iparams_posres.clear();
    iparams_fbposres.clear();
}

void done_idef(t_idef* idef)
{
    sfree(idef->functype);
    sfree(idef->iparams);
    sfree(idef->iparams_posres);
    sfree(idef->iparams_fbposres);
    for (int f = 0; f < F_NRE; ++f)
    {
        sfree(idef->il[f].iatoms);
    }

    init_idef(idef);
}
