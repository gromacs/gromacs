/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2019,2020,2021, by the GROMACS development team, led by
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
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/listed_forces/listed_forces.h"
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
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdlib/wall.h"
#include "gromacs/mdlib/wholemoleculetransform.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/nbnxm/nbnxm.h"
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

#include "gpuforcereduction.h"

ForceHelperBuffers::ForceHelperBuffers(bool haveDirectVirialContributions) :
    haveDirectVirialContributions_(haveDirectVirialContributions)
{
    shiftForces_.resize(SHIFTS);
}

void ForceHelperBuffers::resize(int numAtoms)
{
    if (haveDirectVirialContributions_)
    {
        forceBufferForDirectVirialContributions_.resize(numAtoms);
    }
}

std::vector<real> makeNonBondedParameterLists(const gmx_ffparams_t& forceFieldParams, bool useBuckinghamPotential)
{
    std::vector<real> nbfp;
    int               atnr;

    atnr = forceFieldParams.atnr;
    if (useBuckinghamPotential)
    {
        nbfp.resize(3 * atnr * atnr);
        int k = 0;
        for (int i = 0; (i < atnr); i++)
        {
            for (int j = 0; (j < atnr); j++, k++)
            {
                BHAMA(nbfp, atnr, i, j) = forceFieldParams.iparams[k].bham.a;
                BHAMB(nbfp, atnr, i, j) = forceFieldParams.iparams[k].bham.b;
                /* nbfp now includes the 6.0 derivative prefactor */
                BHAMC(nbfp, atnr, i, j) = forceFieldParams.iparams[k].bham.c * 6.0;
            }
        }
    }
    else
    {
        nbfp.resize(2 * atnr * atnr);
        int k = 0;
        for (int i = 0; (i < atnr); i++)
        {
            for (int j = 0; (j < atnr); j++, k++)
            {
                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                C6(nbfp, atnr, i, j)  = forceFieldParams.iparams[k].lj.c6 * 6.0;
                C12(nbfp, atnr, i, j) = forceFieldParams.iparams[k].lj.c12 * 12.0;
            }
        }
    }

    return nbfp;
}

std::vector<real> makeLJPmeC6GridCorrectionParameters(const gmx_ffparams_t& forceFieldParams,
                                                      const t_forcerec&     forceRec)
{
    int  i, j, k, atnr;
    real c6, c6i, c6j, c12i, c12j, epsi, epsj, sigmai, sigmaj;

    /* For LJ-PME simulations, we correct the energies with the reciprocal space
     * inside of the cut-off. To do this the non-bonded kernels needs to have
     * access to the C6-values used on the reciprocal grid in pme.c
     */

    atnr = forceFieldParams.atnr;
    std::vector<real> grid(2 * atnr * atnr, 0.0);
    for (i = k = 0; (i < atnr); i++)
    {
        for (j = 0; (j < atnr); j++, k++)
        {
            c6i  = forceFieldParams.iparams[i * (atnr + 1)].lj.c6;
            c12i = forceFieldParams.iparams[i * (atnr + 1)].lj.c12;
            c6j  = forceFieldParams.iparams[j * (atnr + 1)].lj.c6;
            c12j = forceFieldParams.iparams[j * (atnr + 1)].lj.c12;
            c6   = std::sqrt(c6i * c6j);
            if (forceRec.ljpme_combination_rule == LongRangeVdW::LB && !gmx_numzero(c6)
                && !gmx_numzero(c12i) && !gmx_numzero(c12j))
            {
                sigmai = gmx::sixthroot(c12i / c6i);
                sigmaj = gmx::sixthroot(c12j / c6j);
                epsi   = c6i * c6i / c12i;
                epsj   = c6j * c6j / c12j;
                c6     = std::sqrt(epsi * epsj) * gmx::power6(0.5 * (sigmai + sigmaj));
            }
            /* Store the elements at the same relative positions as C6 in nbfp in order
             * to simplify access in the kernels
             */
            grid[2 * (atnr * i + j)] = c6 * 6.0;
        }
    }
    return grid;
}

enum
{
    acNONE = 0,
    acCONSTRAINT,
    acSETTLE
};

static std::vector<cginfo_mb_t> init_cginfo_mb(const gmx_mtop_t& mtop, const t_forcerec* fr)
{
    gmx_bool* type_VDW;
    int*      a_con;

    snew(type_VDW, fr->ntype);
    for (int ai = 0; ai < fr->ntype; ai++)
    {
        type_VDW[ai] = FALSE;
        for (int j = 0; j < fr->ntype; j++)
        {
            type_VDW[ai] = type_VDW[ai] || fr->bBHAM || C6(fr->nbfp, fr->ntype, ai, j) != 0
                           || C12(fr->nbfp, fr->ntype, ai, j) != 0;
        }
    }

    std::vector<cginfo_mb_t> cginfoPerMolblock;
    int                      a_offset = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        const gmx_moltype_t&  molt = mtop.moltype[molb.type];
        const auto&           excl = molt.excls;

        /* Check if the cginfo is identical for all molecules in this block.
         * If so, we only need an array of the size of one molecule.
         * Otherwise we make an array of #mol times #cgs per molecule.
         */
        gmx_bool bId = TRUE;
        for (int m = 0; m < molb.nmol; m++)
        {
            const int am = m * molt.atoms.nr;
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                if (getGroupType(mtop.groups, SimulationAtomGroupType::QuantumMechanics, a_offset + am + a)
                    != getGroupType(mtop.groups, SimulationAtomGroupType::QuantumMechanics, a_offset + a))
                {
                    bId = FALSE;
                }
                if (!mtop.groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].empty())
                {
                    if (mtop.groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics][a_offset + am + a]
                        != mtop.groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics][a_offset + a])
                    {
                        bId = FALSE;
                    }
                }
            }
        }

        cginfo_mb_t cginfo_mb;
        cginfo_mb.cg_start = a_offset;
        cginfo_mb.cg_end   = a_offset + molb.nmol * molt.atoms.nr;
        cginfo_mb.cg_mod   = (bId ? 1 : molb.nmol) * molt.atoms.nr;
        cginfo_mb.cginfo.resize(cginfo_mb.cg_mod);
        gmx::ArrayRef<int> cginfo = cginfo_mb.cginfo;

        /* Set constraints flags for constrained atoms */
        snew(a_con, molt.atoms.nr);
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_CONSTRAINT)
            {
                const int nral = NRAL(ftype);
                for (int ia = 0; ia < molt.ilist[ftype].size(); ia += 1 + nral)
                {
                    int a;

                    for (a = 0; a < nral; a++)
                    {
                        a_con[molt.ilist[ftype].iatoms[ia + 1 + a]] =
                                (ftype == F_SETTLE ? acSETTLE : acCONSTRAINT);
                    }
                }
            }
        }

        for (int m = 0; m < (bId ? 1 : molb.nmol); m++)
        {
            const int molculeOffsetInBlock = m * molt.atoms.nr;
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                const t_atom& atom     = molt.atoms.atom[a];
                int&          atomInfo = cginfo[molculeOffsetInBlock + a];

                /* Store the energy group in cginfo */
                int gid = getGroupType(mtop.groups,
                                       SimulationAtomGroupType::EnergyOutput,
                                       a_offset + molculeOffsetInBlock + a);
                SET_CGINFO_GID(atomInfo, gid);

                bool bHaveVDW = (type_VDW[atom.type] || type_VDW[atom.typeB]);
                bool bHaveQ   = (atom.q != 0 || atom.qB != 0);

                bool haveExclusions = false;
                /* Loop over all the exclusions of atom ai */
                for (const int j : excl[a])
                {
                    if (j != a)
                    {
                        haveExclusions = true;
                        break;
                    }
                }

                switch (a_con[a])
                {
                    case acCONSTRAINT: SET_CGINFO_CONSTR(atomInfo); break;
                    case acSETTLE: SET_CGINFO_SETTLE(atomInfo); break;
                    default: break;
                }
                if (haveExclusions)
                {
                    SET_CGINFO_EXCL_INTER(atomInfo);
                }
                if (bHaveVDW)
                {
                    SET_CGINFO_HAS_VDW(atomInfo);
                }
                if (bHaveQ)
                {
                    SET_CGINFO_HAS_Q(atomInfo);
                }
                if (fr->efep != FreeEnergyPerturbationType::No && PERTURBED(atom))
                {
                    SET_CGINFO_FEP(atomInfo);
                }
            }
        }

        sfree(a_con);

        cginfoPerMolblock.push_back(cginfo_mb);

        a_offset += molb.nmol * molt.atoms.nr;
    }
    sfree(type_VDW);

    return cginfoPerMolblock;
}

static std::vector<int> cginfo_expand(const int nmb, gmx::ArrayRef<const cginfo_mb_t> cgi_mb)
{
    const int ncg = cgi_mb[nmb - 1].cg_end;

    std::vector<int> cginfo(ncg);

    int mb = 0;
    for (int cg = 0; cg < ncg; cg++)
    {
        while (cg >= cgi_mb[mb].cg_end)
        {
            mb++;
        }
        cginfo[cg] = cgi_mb[mb].cginfo[(cg - cgi_mb[mb].cg_start) % cgi_mb[mb].cg_mod];
    }

    return cginfo;
}

/* Sets the sum of charges (squared) and C6 in the system in fr.
 * Returns whether the system has a net charge.
 */
static bool set_chargesum(FILE* log, t_forcerec* fr, const gmx_mtop_t& mtop)
{
    /*This now calculates sum for q and c6*/
    double qsum, q2sum, q, c6sum, c6;

    qsum  = 0;
    q2sum = 0;
    c6sum = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int            nmol  = molb.nmol;
        const t_atoms* atoms = &mtop.moltype[molb.type].atoms;
        for (int i = 0; i < atoms->nr; i++)
        {
            q = atoms->atom[i].q;
            qsum += nmol * q;
            q2sum += nmol * q * q;
            c6 = mtop.ffparams.iparams[atoms->atom[i].type * (mtop.ffparams.atnr + 1)].lj.c6;
            c6sum += nmol * c6;
        }
    }
    fr->qsum[0]  = qsum;
    fr->q2sum[0] = q2sum;
    fr->c6sum[0] = c6sum;

    if (fr->efep != FreeEnergyPerturbationType::No)
    {
        qsum  = 0;
        q2sum = 0;
        c6sum = 0;
        for (const gmx_molblock_t& molb : mtop.molblock)
        {
            int            nmol  = molb.nmol;
            const t_atoms* atoms = &mtop.moltype[molb.type].atoms;
            for (int i = 0; i < atoms->nr; i++)
            {
                q = atoms->atom[i].qB;
                qsum += nmol * q;
                q2sum += nmol * q * q;
                c6 = mtop.ffparams.iparams[atoms->atom[i].typeB * (mtop.ffparams.atnr + 1)].lj.c6;
                c6sum += nmol * c6;
            }
            fr->qsum[1]  = qsum;
            fr->q2sum[1] = q2sum;
            fr->c6sum[1] = c6sum;
        }
    }
    else
    {
        fr->qsum[1]  = fr->qsum[0];
        fr->q2sum[1] = fr->q2sum[0];
        fr->c6sum[1] = fr->c6sum[0];
    }
    if (log)
    {
        if (fr->efep == FreeEnergyPerturbationType::No)
        {
            fprintf(log, "System total charge: %.3f\n", fr->qsum[0]);
        }
        else
        {
            fprintf(log, "System total charge, top. A: %.3f top. B: %.3f\n", fr->qsum[0], fr->qsum[1]);
        }
    }

    /* A cut-off of 1e-4 is used to catch rounding errors due to ascii input */
    return (std::abs(fr->qsum[0]) > 1e-4 || std::abs(fr->qsum[1]) > 1e-4);
}

static real calcBuckinghamBMax(FILE* fplog, const gmx_mtop_t& mtop)
{
    const t_atoms *at1, *at2;
    int            i, j, tpi, tpj, ntypes;
    real           b, bmin;

    if (fplog)
    {
        fprintf(fplog, "Determining largest Buckingham b parameter for table\n");
    }
    ntypes = mtop.ffparams.atnr;

    bmin            = -1;
    real bham_b_max = 0;
    for (size_t mt1 = 0; mt1 < mtop.moltype.size(); mt1++)
    {
        at1 = &mtop.moltype[mt1].atoms;
        for (i = 0; (i < at1->nr); i++)
        {
            tpi = at1->atom[i].type;
            if (tpi >= ntypes)
            {
                gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", i, tpi, ntypes);
            }

            for (size_t mt2 = mt1; mt2 < mtop.moltype.size(); mt2++)
            {
                at2 = &mtop.moltype[mt2].atoms;
                for (j = 0; (j < at2->nr); j++)
                {
                    tpj = at2->atom[j].type;
                    if (tpj >= ntypes)
                    {
                        gmx_fatal(FARGS, "Atomtype[%d] = %d, maximum = %d", j, tpj, ntypes);
                    }
                    b = mtop.ffparams.iparams[tpi * ntypes + tpj].bham.b;
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
        fprintf(fplog, "Buckingham b parameters, min: %g, max: %g\n", bmin, bham_b_max);
    }

    return bham_b_max;
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
static void count_tables(int ftype1, int ftype2, const gmx_mtop_t& mtop, int* ncount, int** count)
{
    int ftype, i, j, tabnr;

    // Loop over all moleculetypes
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        // Loop over all interaction types
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            // If the current interaction type is one of the types whose tables we're trying to count...
            if (ftype == ftype1 || ftype == ftype2)
            {
                const InteractionList& il     = molt.ilist[ftype];
                const int              stride = 1 + NRAL(ftype);
                // ... and there are actually some interactions for this type
                for (i = 0; i < il.size(); i += stride)
                {
                    // Find out which table index the user wanted
                    tabnr = mtop.ffparams.iparams[il.iatoms[i]].tab.table;
                    if (tabnr < 0)
                    {
                        gmx_fatal(FARGS, "A bonded table number is smaller than 0: %d\n", tabnr);
                    }
                    // Make room for this index in the data structure
                    if (tabnr >= *ncount)
                    {
                        srenew(*count, tabnr + 1);
                        for (j = *ncount; j < tabnr + 1; j++)
                        {
                            (*count)[j] = 0;
                        }
                        *ncount = tabnr + 1;
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
static std::vector<bondedtable_t> make_bonded_tables(FILE*                            fplog,
                                                     int                              ftype1,
                                                     int                              ftype2,
                                                     const gmx_mtop_t&                mtop,
                                                     gmx::ArrayRef<const std::string> tabbfnm,
                                                     const char*                      tabext)
{
    std::vector<bondedtable_t> tab;

    int  ncount = 0;
    int* count  = nullptr;
    count_tables(ftype1, ftype2, mtop, &ncount, &count);

    // Are there any relevant tabulated bond interactions?
    if (ncount > 0)
    {
        tab.resize(ncount);
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
                        tab[i]    = make_bonded_table(fplog, tabbfnm[j].c_str(), NRAL(ftype1) - 2);
                        madeTable = true;
                    }
                }
                if (!madeTable)
                {
                    bool isPlural = (ftype2 != -1);
                    gmx_fatal(FARGS,
                              "Tabulated interaction of type '%s%s%s' with index %d cannot be used "
                              "because no table file whose name matched '%s' was passed via the "
                              "gmx mdrun -tableb command-line option.",
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

void forcerec_set_ranges(t_forcerec* fr, int natoms_force, int natoms_force_constr, int natoms_f_novirsum)
{
    fr->natoms_force        = natoms_force;
    fr->natoms_force_constr = natoms_force_constr;

    for (auto& forceHelperBuffers : fr->forceHelperBuffers)
    {
        forceHelperBuffers.resize(natoms_f_novirsum);
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
static void initCoulombEwaldParameters(FILE*                fp,
                                       const t_inputrec&    ir,
                                       bool                 systemHasNetCharge,
                                       interaction_const_t* ic)
{
    if (!EEL_PME_EWALD(ir.coulombtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for electrostatic interactions.\n");

        if (ir.coulombtype == CoulombInteractionType::P3mAD)
        {
            please_cite(fp, "Hockney1988");
            please_cite(fp, "Ballenegger2012");
        }
        else
        {
            please_cite(fp, "Essmann95a");
        }

        if (ir.ewald_geometry == EwaldGeometry::ThreeDC)
        {
            if (fp)
            {
                fprintf(fp,
                        "Using the Ewald3DC correction for systems with a slab geometry%s.\n",
                        systemHasNetCharge ? " and net charge" : "");
            }
            please_cite(fp, "In-Chul99a");
            if (systemHasNetCharge)
            {
                please_cite(fp, "Ballenegger2009");
            }
        }
    }

    ic->ewaldcoeff_q = calc_ewaldcoeff_q(ir.rcoulomb, ir.ewald_rtol);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for Ewald\n", 1 / ic->ewaldcoeff_q);
    }

    if (ic->coulomb_modifier == InteractionModifiers::PotShift)
    {
        GMX_RELEASE_ASSERT(ic->rcoulomb != 0, "Cutoff radius cannot be zero");
        ic->sh_ewald = std::erfc(ic->ewaldcoeff_q * ic->rcoulomb) / ic->rcoulomb;
    }
    else
    {
        ic->sh_ewald = 0;
    }
}

/*! \brief Print Van der Waals Ewald citations and set ewald coefficients */
static void initVdwEwaldParameters(FILE* fp, const t_inputrec& ir, interaction_const_t* ic)
{
    if (!EVDW_PME(ir.vdwtype))
    {
        return;
    }

    if (fp)
    {
        fprintf(fp, "Will do PME sum in reciprocal space for LJ dispersion interactions.\n");
        please_cite(fp, "Essmann95a");
    }
    ic->ewaldcoeff_lj = calc_ewaldcoeff_lj(ir.rvdw, ir.ewald_rtol_lj);
    if (fp)
    {
        fprintf(fp, "Using a Gaussian width (1/beta) of %g nm for LJ Ewald\n", 1 / ic->ewaldcoeff_lj);
    }

    if (ic->vdw_modifier == InteractionModifiers::PotShift)
    {
        real crc2       = gmx::square(ic->ewaldcoeff_lj * ic->rvdw);
        ic->sh_lj_ewald = (std::exp(-crc2) * (1 + crc2 + 0.5 * crc2 * crc2) - 1) / gmx::power6(ic->rvdw);
    }
    else
    {
        ic->sh_lj_ewald = 0;
    }
}

/* Generate Coulomb and/or Van der Waals Ewald long-range correction tables
 *
 * Tables are generated for one or both, depending on if the pointers
 * are non-null. The spacing for both table sets is the same and obeys
 * both accuracy requirements, when relevant.
 */
static void init_ewald_f_table(const interaction_const_t& ic,
                               const real                 rlist,
                               const real                 tabext,
                               EwaldCorrectionTables*     coulombTables,
                               EwaldCorrectionTables*     vdwTables)
{
    const bool useCoulombTable = (EEL_PME_EWALD(ic.eeltype) && coulombTables != nullptr);
    const bool useVdwTable     = (EVDW_PME(ic.vdwtype) && vdwTables != nullptr);

    /* Get the Ewald table spacing based on Coulomb and/or LJ
     * Ewald coefficients and rtol.
     */
    const real tableScale = ewald_spline3_table_scale(ic, useCoulombTable, useVdwTable);

    const bool havePerturbedNonbondeds = (ic.softCoreParameters != nullptr);

    real tableLen = ic.rcoulomb;
    if ((useCoulombTable || useVdwTable) && havePerturbedNonbondeds && rlist + tabext > 0.0)
    {
        /* TODO: Ideally this should also check if couple-intramol == no, but that isn't
         * stored in ir. Grompp puts that info into an opts structure that doesn't make it into the tpr.
         * The alternative is to look through all the exclusions and check if they come from
         * couple-intramol == no. Meanwhile, always having larger tables should only affect
         * memory consumption, not speed (barring cache issues).
         */
        tableLen = rlist + tabext;
    }
    const int tableSize = static_cast<int>(tableLen * tableScale) + 2;

    if (useCoulombTable)
    {
        *coulombTables =
                generateEwaldCorrectionTables(tableSize, tableScale, ic.ewaldcoeff_q, v_q_ewald_lr);
    }

    if (useVdwTable)
    {
        *vdwTables = generateEwaldCorrectionTables(tableSize, tableScale, ic.ewaldcoeff_lj, v_lj_ewald_lr);
    }
}

void init_interaction_const_tables(FILE* fp, interaction_const_t* ic, const real rlist, const real tableExtensionLength)
{
    if (EEL_PME_EWALD(ic->eeltype) || EVDW_PME(ic->vdwtype))
    {
        init_ewald_f_table(
                *ic, rlist, tableExtensionLength, ic->coulombEwaldTables.get(), ic->vdwEwaldTables.get());
        if (fp != nullptr)
        {
            fprintf(fp,
                    "Initialized non-bonded Ewald tables, spacing: %.2e size: %zu\n\n",
                    1 / ic->coulombEwaldTables->scale,
                    ic->coulombEwaldTables->tableF.size());
        }
    }
}

static void clear_force_switch_constants(shift_consts_t* sc)
{
    sc->c2   = 0;
    sc->c3   = 0;
    sc->cpot = 0;
}

static void force_switch_constants(real p, real rsw, real rc, shift_consts_t* sc)
{
    /* Here we determine the coefficient for shifting the force to zero
     * between distance rsw and the cut-off rc.
     * For a potential of r^-p, we have force p*r^-(p+1).
     * But to save flops we absorb p in the coefficient.
     * Thus we get:
     * force/p   = r^-(p+1) + c2*r^2 + c3*r^3
     * potential = r^-p + c2/3*r^3 + c3/4*r^4 + cpot
     */
    sc->c2   = ((p + 1) * rsw - (p + 4) * rc) / (pow(rc, p + 2) * gmx::square(rc - rsw));
    sc->c3   = -((p + 1) * rsw - (p + 3) * rc) / (pow(rc, p + 2) * gmx::power3(rc - rsw));
    sc->cpot = -pow(rc, -p) + p * sc->c2 / 3 * gmx::power3(rc - rsw)
               + p * sc->c3 / 4 * gmx::power4(rc - rsw);
}

static void potential_switch_constants(real rsw, real rc, switch_consts_t* sc)
{
    /* The switch function is 1 at rsw and 0 at rc.
     * The derivative and second derivate are zero at both ends.
     * rsw        = max(r - r_switch, 0)
     * sw         = 1 + c3*rsw^3 + c4*rsw^4 + c5*rsw^5
     * dsw        = 3*c3*rsw^2 + 4*c4*rsw^3 + 5*c5*rsw^4
     * force      = force*dsw - potential*sw
     * potential *= sw
     */
    sc->c3 = -10 / gmx::power3(rc - rsw);
    sc->c4 = 15 / gmx::power4(rc - rsw);
    sc->c5 = -6 / gmx::power5(rc - rsw);
}

/*! \brief Construct interaction constants
 *
 * This data is used (particularly) by search and force code for
 * short-range interactions. Many of these are constant for the whole
 * simulation; some are constant only after PME tuning completes.
 */
static interaction_const_t init_interaction_const(FILE*             fp,
                                                  const t_inputrec& ir,
                                                  const gmx_mtop_t& mtop,
                                                  bool              systemHasNetCharge)
{
    interaction_const_t interactionConst;

    interactionConst.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
    interactionConst.vdwEwaldTables     = std::make_unique<EwaldCorrectionTables>();

    /* Lennard-Jones */
    interactionConst.vdwtype         = ir.vdwtype;
    interactionConst.vdw_modifier    = ir.vdw_modifier;
    interactionConst.reppow          = mtop.ffparams.reppow;
    interactionConst.rvdw            = cutoff_inf(ir.rvdw);
    interactionConst.rvdw_switch     = ir.rvdw_switch;
    interactionConst.ljpme_comb_rule = ir.ljpme_combination_rule;
    interactionConst.useBuckingham   = (mtop.ffparams.functype[0] == F_BHAM);
    if (interactionConst.useBuckingham)
    {
        interactionConst.buckinghamBMax = calcBuckinghamBMax(fp, mtop);
    }

    initVdwEwaldParameters(fp, ir, &interactionConst);

    clear_force_switch_constants(&interactionConst.dispersion_shift);
    clear_force_switch_constants(&interactionConst.repulsion_shift);

    switch (interactionConst.vdw_modifier)
    {
        case InteractionModifiers::PotShift:
            /* Only shift the potential, don't touch the force */
            interactionConst.dispersion_shift.cpot = -1.0 / gmx::power6(interactionConst.rvdw);
            interactionConst.repulsion_shift.cpot  = -1.0 / gmx::power12(interactionConst.rvdw);
            break;
        case InteractionModifiers::ForceSwitch:
            /* Switch the force, switch and shift the potential */
            force_switch_constants(
                    6.0, interactionConst.rvdw_switch, interactionConst.rvdw, &interactionConst.dispersion_shift);
            force_switch_constants(
                    12.0, interactionConst.rvdw_switch, interactionConst.rvdw, &interactionConst.repulsion_shift);
            break;
        case InteractionModifiers::PotSwitch:
            /* Switch the potential and force */
            potential_switch_constants(
                    interactionConst.rvdw_switch, interactionConst.rvdw, &interactionConst.vdw_switch);
            break;
        case InteractionModifiers::None:
        case InteractionModifiers::ExactCutoff:
            /* Nothing to do here */
            break;
        default: gmx_incons("unimplemented potential modifier");
    }

    /* Electrostatics */
    interactionConst.eeltype          = ir.coulombtype;
    interactionConst.coulomb_modifier = ir.coulomb_modifier;
    interactionConst.rcoulomb         = cutoff_inf(ir.rcoulomb);
    interactionConst.rcoulomb_switch  = ir.rcoulomb_switch;
    interactionConst.epsilon_r        = ir.epsilon_r;

    /* Set the Coulomb energy conversion factor */
    if (interactionConst.epsilon_r != 0)
    {
        interactionConst.epsfac = ONE_4PI_EPS0 / interactionConst.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        interactionConst.epsfac = 0;
    }

    /* Reaction-field */
    if (EEL_RF(interactionConst.eeltype))
    {
        GMX_RELEASE_ASSERT(interactionConst.eeltype != CoulombInteractionType::GRFNotused,
                           "GRF is no longer supported");
        interactionConst.reactionFieldPermitivity = ir.epsilon_rf;
        calc_rffac(fp,
                   interactionConst.epsilon_r,
                   interactionConst.reactionFieldPermitivity,
                   interactionConst.rcoulomb,
                   &interactionConst.reactionFieldCoefficient,
                   &interactionConst.reactionFieldShift);
    }
    else
    {
        /* For plain cut-off we might use the reaction-field kernels */
        interactionConst.reactionFieldPermitivity = interactionConst.epsilon_r;
        interactionConst.reactionFieldCoefficient = 0;
        if (ir.coulomb_modifier == InteractionModifiers::PotShift)
        {
            interactionConst.reactionFieldShift = 1 / interactionConst.rcoulomb;
        }
        else
        {
            interactionConst.reactionFieldShift = 0;
        }
    }

    initCoulombEwaldParameters(fp, ir, systemHasNetCharge, &interactionConst);

    if (fp != nullptr)
    {
        real dispersion_shift;

        dispersion_shift = interactionConst.dispersion_shift.cpot;
        if (EVDW_PME(interactionConst.vdwtype))
        {
            dispersion_shift -= interactionConst.sh_lj_ewald;
        }
        fprintf(fp,
                "Potential shift: LJ r^-12: %.3e r^-6: %.3e",
                interactionConst.repulsion_shift.cpot,
                dispersion_shift);

        if (interactionConst.eeltype == CoulombInteractionType::Cut)
        {
            fprintf(fp, ", Coulomb %.e", -interactionConst.reactionFieldShift);
        }
        else if (EEL_PME(interactionConst.eeltype))
        {
            fprintf(fp, ", Ewald %.3e", -interactionConst.sh_ewald);
        }
        fprintf(fp, "\n");
    }

    if (ir.efep != FreeEnergyPerturbationType::No)
    {
        GMX_RELEASE_ASSERT(ir.fepvals, "ir.fepvals should be set wth free-energy");
        interactionConst.softCoreParameters =
                std::make_unique<interaction_const_t::SoftCoreParameters>(*ir.fepvals);
    }

    return interactionConst;
}

void init_forcerec(FILE*                            fplog,
                   const gmx::MDLogger&             mdlog,
                   t_forcerec*                      forcerec,
                   const t_inputrec&                inputrec,
                   const gmx_mtop_t&                mtop,
                   const t_commrec*                 commrec,
                   matrix                           box,
                   const char*                      tabfn,
                   const char*                      tabpfn,
                   gmx::ArrayRef<const std::string> tabbfnm,
                   real                             print_force)
{
    /* The CMake default turns SIMD kernels on, but it might be turned off further down... */
    forcerec->use_simd_kernels = GMX_USE_SIMD_KERNELS;

    if (check_box(inputrec.pbcType, box))
    {
        gmx_fatal(FARGS, "%s", check_box(inputrec.pbcType, box));
    }

    /* Test particle insertion ? */
    if (EI_TPI(inputrec.eI))
    {
        /* Set to the size of the molecule to be inserted (the last one) */
        gmx::RangePartitioning molecules = gmx_mtop_molecules(mtop);
        forcerec->n_tpi                  = molecules.block(molecules.numBlocks() - 1).size();
    }
    else
    {
        forcerec->n_tpi = 0;
    }

    if (inputrec.coulombtype == CoulombInteractionType::RFNecUnsupported
        || inputrec.coulombtype == CoulombInteractionType::GRFNotused)
    {
        gmx_fatal(FARGS, "%s electrostatics is no longer supported", enumValueToString(inputrec.coulombtype));
    }

    if (inputrec.bAdress)
    {
        gmx_fatal(FARGS, "AdResS simulations are no longer supported");
    }
    if (inputrec.useTwinRange)
    {
        gmx_fatal(FARGS, "Twin-range simulations are no longer supported");
    }
    /* Copy the user determined parameters */
    forcerec->userint1  = inputrec.userint1;
    forcerec->userint2  = inputrec.userint2;
    forcerec->userint3  = inputrec.userint3;
    forcerec->userint4  = inputrec.userint4;
    forcerec->userreal1 = inputrec.userreal1;
    forcerec->userreal2 = inputrec.userreal2;
    forcerec->userreal3 = inputrec.userreal3;
    forcerec->userreal4 = inputrec.userreal4;

    /* Shell stuff */
    forcerec->fc_stepsize = inputrec.fc_stepsize;

    /* Free energy */
    forcerec->efep = inputrec.efep;

    if ((getenv("GMX_DISABLE_SIMD_KERNELS") != nullptr) || (getenv("GMX_NOOPTIMIZEDKERNELS") != nullptr))
    {
        forcerec->use_simd_kernels = FALSE;
        if (fplog != nullptr)
        {
            fprintf(fplog,
                    "\nFound environment variable GMX_DISABLE_SIMD_KERNELS.\n"
                    "Disabling the usage of any SIMD-specific non-bonded & bonded kernel routines\n"
                    "(e.g. SSE2/SSE4.1/AVX).\n\n");
        }
    }

    forcerec->bBHAM = (mtop.ffparams.functype[0] == F_BHAM);

    /* Neighbour searching stuff */
    forcerec->pbcType = inputrec.pbcType;

    /* Determine if we will do PBC for distances in bonded interactions */
    if (forcerec->pbcType == PbcType::No)
    {
        forcerec->bMolPBC = FALSE;
    }
    else
    {
        const bool useEwaldSurfaceCorrection =
                (EEL_PME_EWALD(inputrec.coulombtype) && inputrec.epsilon_surface != 0);
        const bool haveOrientationRestraints = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
        if (!DOMAINDECOMP(commrec))
        {
            forcerec->bMolPBC = true;

            if (useEwaldSurfaceCorrection || haveOrientationRestraints)
            {
                forcerec->wholeMoleculeTransform =
                        std::make_unique<gmx::WholeMoleculeTransform>(mtop, inputrec.pbcType);
            }
        }
        else
        {
            forcerec->bMolPBC = dd_bonded_molpbc(*commrec->dd, forcerec->pbcType);

            /* With Ewald surface correction it is useful to support e.g. running water
             * in parallel with update groups.
             * With orientation restraints there is no sensible use case supported with DD.
             */
            if ((useEwaldSurfaceCorrection && !dd_moleculesAreAlwaysWhole(*commrec->dd))
                || haveOrientationRestraints)
            {
                gmx_fatal(FARGS,
                          "You requested Ewald surface correction or orientation restraints, "
                          "but molecules are broken "
                          "over periodic boundary conditions by the domain decomposition. "
                          "Run without domain decomposition instead.");
            }
        }

        if (useEwaldSurfaceCorrection)
        {
            GMX_RELEASE_ASSERT(!DOMAINDECOMP(commrec) || dd_moleculesAreAlwaysWhole(*commrec->dd),
                               "Molecules can not be broken by PBC with epsilon_surface > 0");
        }
    }

    forcerec->rc_scaling = inputrec.refcoord_scaling;
    copy_rvec(inputrec.posres_com, forcerec->posres_com);
    copy_rvec(inputrec.posres_comB, forcerec->posres_comB);
    forcerec->rlist                  = cutoff_inf(inputrec.rlist);
    forcerec->ljpme_combination_rule = inputrec.ljpme_combination_rule;

    /* This now calculates sum for q and c6*/
    bool systemHasNetCharge = set_chargesum(fplog, forcerec, mtop);

    /* Make data structure used by kernels */
    forcerec->ic = std::make_unique<interaction_const_t>(
            init_interaction_const(fplog, inputrec, mtop, systemHasNetCharge));
    init_interaction_const_tables(fplog, forcerec->ic.get(), forcerec->rlist, inputrec.tabext);

    const interaction_const_t* interactionConst = forcerec->ic.get();

    /* TODO: Replace this Ewald table or move it into interaction_const_t */
    if (inputrec.coulombtype == CoulombInteractionType::Ewald)
    {
        forcerec->ewald_table = std::make_unique<gmx_ewald_tab_t>(inputrec, fplog);
    }

    /* Electrostatics: Translate from interaction-setting-in-mdp-file to kernel interaction format */
    switch (interactionConst->eeltype)
    {
        case CoulombInteractionType::Cut:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::Coulomb;
            break;

        case CoulombInteractionType::RF:
        case CoulombInteractionType::RFZero:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::ReactionField;
            break;

        case CoulombInteractionType::Switch:
        case CoulombInteractionType::Shift:
        case CoulombInteractionType::User:
        case CoulombInteractionType::PmeSwitch:
        case CoulombInteractionType::PmeUser:
        case CoulombInteractionType::PmeUserSwitch:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::CubicSplineTable;
            break;

        case CoulombInteractionType::Pme:
        case CoulombInteractionType::P3mAD:
        case CoulombInteractionType::Ewald:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::Ewald;
            break;

        default:
            gmx_fatal(FARGS,
                      "Unsupported electrostatic interaction: %s",
                      enumValueToString(interactionConst->eeltype));
    }
    forcerec->nbkernel_elec_modifier = interactionConst->coulomb_modifier;

    /* Vdw: Translate from mdp settings to kernel format */
    switch (interactionConst->vdwtype)
    {
        case VanDerWaalsType::Cut:
            if (forcerec->bBHAM)
            {
                forcerec->nbkernel_vdw_interaction = NbkernelVdwType::Buckingham;
            }
            else
            {
                forcerec->nbkernel_vdw_interaction = NbkernelVdwType::LennardJones;
            }
            break;
        case VanDerWaalsType::Pme:
            forcerec->nbkernel_vdw_interaction = NbkernelVdwType::LJEwald;
            break;

        case VanDerWaalsType::Switch:
        case VanDerWaalsType::Shift:
        case VanDerWaalsType::User:
            forcerec->nbkernel_vdw_interaction = NbkernelVdwType::CubicSplineTable;
            break;

        default:
            gmx_fatal(FARGS, "Unsupported vdw interaction: %s", enumValueToString(interactionConst->vdwtype));
    }
    forcerec->nbkernel_vdw_modifier = interactionConst->vdw_modifier;

    if (!gmx_within_tol(interactionConst->reppow, 12.0, 10 * GMX_DOUBLE_EPS))
    {
        gmx_fatal(FARGS, "Only LJ repulsion power 12 is supported");
    }
    /* Older tpr files can contain Coulomb user tables with the Verlet cutoff-scheme,
     * while mdrun does not (and never did) support this.
     */
    if (EEL_USER(forcerec->ic->eeltype))
    {
        gmx_fatal(FARGS,
                  "Electrostatics type %s is currently not supported",
                  enumValueToString(inputrec.coulombtype));
    }

    forcerec->bvdwtab  = FALSE;
    forcerec->bcoultab = FALSE;

    /* 1-4 interaction electrostatics */
    forcerec->fudgeQQ = mtop.ffparams.fudgeQQ;

    // Multiple time stepping
    forcerec->useMts = inputrec.useMts;

    if (forcerec->useMts)
    {
        GMX_ASSERT(gmx::checkMtsRequirements(inputrec).empty(),
                   "All MTS requirements should be met here");
    }

    const bool haveDirectVirialContributionsFast =
            forcerec->forceProviders->hasForceProvider() || gmx_mtop_ftype_count(mtop, F_POSRES) > 0
            || gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0 || inputrec.nwall > 0 || inputrec.bPull
            || inputrec.bRot || inputrec.bIMD;
    const bool haveDirectVirialContributionsSlow =
            EEL_FULL(interactionConst->eeltype) || EVDW_PME(interactionConst->vdwtype);
    for (int i = 0; i < (forcerec->useMts ? 2 : 1); i++)
    {
        bool haveDirectVirialContributions =
                (((!forcerec->useMts || i == 0) && haveDirectVirialContributionsFast)
                 || ((!forcerec->useMts || i == 1) && haveDirectVirialContributionsSlow));
        forcerec->forceHelperBuffers.emplace_back(haveDirectVirialContributions);
    }

    if (forcerec->shift_vec == nullptr)
    {
        snew(forcerec->shift_vec, SHIFTS);
    }

    if (forcerec->nbfp.empty())
    {
        forcerec->ntype = mtop.ffparams.atnr;
        forcerec->nbfp  = makeNonBondedParameterLists(mtop.ffparams, forcerec->bBHAM);
        if (EVDW_PME(interactionConst->vdwtype))
        {
            forcerec->ljpme_c6grid = makeLJPmeC6GridCorrectionParameters(mtop.ffparams, *forcerec);
        }
    }

    /* Copy the energy group exclusions */
    forcerec->egp_flags = inputrec.opts.egp_flags;

    /* Van der Waals stuff */
    if ((interactionConst->vdwtype != VanDerWaalsType::Cut)
        && (interactionConst->vdwtype != VanDerWaalsType::User) && !forcerec->bBHAM)
    {
        if (interactionConst->rvdw_switch >= interactionConst->rvdw)
        {
            gmx_fatal(FARGS,
                      "rvdw_switch (%f) must be < rvdw (%f)",
                      interactionConst->rvdw_switch,
                      interactionConst->rvdw);
        }
        if (fplog)
        {
            fprintf(fplog,
                    "Using %s Lennard-Jones, switch between %g and %g nm\n",
                    (interactionConst->eeltype == CoulombInteractionType::Switch) ? "switched" : "shifted",
                    interactionConst->rvdw_switch,
                    interactionConst->rvdw);
        }
    }

    if (forcerec->bBHAM && EVDW_PME(interactionConst->vdwtype))
    {
        gmx_fatal(FARGS, "LJ PME not supported with Buckingham");
    }

    if (forcerec->bBHAM
        && (interactionConst->vdwtype == VanDerWaalsType::Shift
            || interactionConst->vdwtype == VanDerWaalsType::Switch))
    {
        gmx_fatal(FARGS, "Switch/shift interaction not supported with Buckingham");
    }

    if (forcerec->bBHAM)
    {
        gmx_fatal(FARGS, "The Verlet cutoff-scheme does not (yet) support Buckingham");
    }

    if (inputrec.implicit_solvent)
    {
        gmx_fatal(FARGS, "Implict solvation is no longer supported.");
    }


    /* This code automatically gives table length tabext without cut-off's,
     * in that case grompp should already have checked that we do not need
     * normal tables and we only generate tables for 1-4 interactions.
     */
    real rtab = inputrec.rlist + inputrec.tabext;

    /* We want to use unmodified tables for 1-4 coulombic
     * interactions, so we must in general have an extra set of
     * tables. */
    if (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 || gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0
        || gmx_mtop_ftype_count(mtop, F_LJC_PAIRS_NB) > 0)
    {
        forcerec->pairsTable = make_tables(fplog, interactionConst, tabpfn, rtab, GMX_MAKETABLES_14ONLY);
    }

    /* Wall stuff */
    forcerec->nwall = inputrec.nwall;
    if (inputrec.nwall && inputrec.wall_type == WallType::Table)
    {
        make_wall_tables(fplog, inputrec, tabfn, &mtop.groups, forcerec);
    }

    forcerec->fcdata = std::make_unique<t_fcdata>();

    if (!tabbfnm.empty())
    {
        t_fcdata& fcdata = *forcerec->fcdata;
        // Need to catch std::bad_alloc
        // TODO Don't need to catch this here, when merging with master branch
        try
        {
            // TODO move these tables into a separate struct and store reference in ListedForces
            fcdata.bondtab = make_bonded_tables(fplog, F_TABBONDS, F_TABBONDSNC, mtop, tabbfnm, "b");
            fcdata.angletab = make_bonded_tables(fplog, F_TABANGLES, -1, mtop, tabbfnm, "a");
            fcdata.dihtab   = make_bonded_tables(fplog, F_TABDIHS, -1, mtop, tabbfnm, "d");
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    else
    {
        if (debug)
        {
            fprintf(debug,
                    "No fcdata or table file name passed, can not read table, can not do bonded "
                    "interactions\n");
        }
    }

    /* Initialize the thread working data for bonded interactions */
    if (forcerec->useMts)
    {
        // Add one ListedForces object for each MTS level
        bool isFirstLevel = true;
        for (const auto& mtsLevel : inputrec.mtsLevels)
        {
            ListedForces::InteractionSelection interactionSelection;
            const auto&                        forceGroups = mtsLevel.forceGroups;
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Pair)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Pairs));
            }
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Dihedral)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Dihedrals));
            }
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Angle)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Angles));
            }
            if (isFirstLevel)
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Rest));
                isFirstLevel = false;
            }
            forcerec->listedForces.emplace_back(
                    mtop.ffparams,
                    mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                    gmx_omp_nthreads_get(emntBonded),
                    interactionSelection,
                    fplog);
        }
    }
    else
    {
        // Add one ListedForces object with all listed interactions
        forcerec->listedForces.emplace_back(
                mtop.ffparams,
                mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                gmx_omp_nthreads_get(emntBonded),
                ListedForces::interactionSelectionAll(),
                fplog);
    }

    // QM/MM initialization if requested
    if (inputrec.bQMMM)
    {
        gmx_incons("QM/MM was requested, but is no longer available in GROMACS");
    }

    /* Set all the static charge group info */
    forcerec->cginfo_mb = init_cginfo_mb(mtop, forcerec);
    if (!DOMAINDECOMP(commrec))
    {
        forcerec->cginfo = cginfo_expand(mtop.molblock.size(), forcerec->cginfo_mb);
    }

    if (!DOMAINDECOMP(commrec))
    {
        forcerec_set_ranges(forcerec, mtop.natoms, mtop.natoms, mtop.natoms);
    }

    forcerec->print_force = print_force;

    forcerec->nthread_ewc = gmx_omp_nthreads_get(emntBonded);
    snew(forcerec->ewc_t, forcerec->nthread_ewc);

    if (inputrec.eDispCorr != DispersionCorrectionType::No)
    {
        forcerec->dispersionCorrection = std::make_unique<DispersionCorrection>(
                mtop, inputrec, forcerec->bBHAM, forcerec->ntype, forcerec->nbfp, *forcerec->ic, tabfn);
        forcerec->dispersionCorrection->print(mdlog);
    }

    if (fplog != nullptr)
    {
        /* Here we switch from using mdlog, which prints the newline before
         * the paragraph, to our old fprintf logging, which prints the newline
         * after the paragraph, so we should add a newline here.
         */
        fprintf(fplog, "\n");
    }
}

t_forcerec::t_forcerec() = default;

t_forcerec::~t_forcerec()
{
    /* Note: This code will disappear when types are converted to C++ */
    sfree(shift_vec);
    sfree(ewc_t);
}
