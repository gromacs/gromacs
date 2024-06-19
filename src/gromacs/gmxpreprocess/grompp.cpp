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

#include "grompp.h"

#include <cerrno>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <sys/types.h>

#include "gromacs/applied_forces/awh/read_params.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fft/calcgrid.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/add_par.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gen_maxwell_velocities.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/massrepartitioning.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/tomorse.h"
#include "gromacs/gmxpreprocess/topio.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/imd/imd.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/perf_est.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/seed.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;
struct pull_t;

/* TODO The implementation details should move to their own source file. */
InteractionOfType::InteractionOfType(gmx::ArrayRef<const int>  atoms,
                                     gmx::ArrayRef<const real> params,
                                     const std::string&        name) :
    atoms_(atoms.begin(), atoms.end()), interactionTypeName_(name)
{
    GMX_RELEASE_ASSERT(
            params.size() <= forceParam_.size(),
            gmx::formatString("Cannot have more parameters than the maximum number possible (%d)", MAXFORCEPARAM)
                    .c_str());
    std::array<real, MAXFORCEPARAM>::iterator forceParamIt = forceParam_.begin();
    for (const auto param : params)
    {
        *forceParamIt++ = param;
    }
    while (forceParamIt != forceParam_.end())
    {
        *forceParamIt++ = NOTSET;
    }
}

const int& InteractionOfType::ai() const
{
    GMX_RELEASE_ASSERT(!atoms_.empty(), "Need to have at least one atom set");
    return atoms_[0];
}

const int& InteractionOfType::aj() const
{
    GMX_RELEASE_ASSERT(atoms_.size() > 1, "Need to have at least two atoms set");
    return atoms_[1];
}

const int& InteractionOfType::ak() const
{
    GMX_RELEASE_ASSERT(atoms_.size() > 2, "Need to have at least three atoms set");
    return atoms_[2];
}

const int& InteractionOfType::al() const
{
    GMX_RELEASE_ASSERT(atoms_.size() > 3, "Need to have at least four atoms set");
    return atoms_[3];
}

const int& InteractionOfType::am() const
{
    GMX_RELEASE_ASSERT(atoms_.size() > 4, "Need to have at least five atoms set");
    return atoms_[4];
}

const real& InteractionOfType::c0() const
{
    return forceParam_[0];
}

const real& InteractionOfType::c1() const
{
    return forceParam_[1];
}

const real& InteractionOfType::c2() const
{
    return forceParam_[2];
}

const std::string& InteractionOfType::interactionTypeName() const
{
    return interactionTypeName_;
}

void InteractionOfType::sortBondAtomIds()
{
    if (aj() < ai())
    {
        std::swap(atoms_[0], atoms_[1]);
    }
}

void InteractionOfType::sortAngleAtomIds()
{
    if (ak() < ai())
    {
        std::swap(atoms_[0], atoms_[2]);
    }
}

void InteractionOfType::sortDihedralAtomIds()
{
    if (al() < ai())
    {
        std::swap(atoms_[0], atoms_[3]);
        std::swap(atoms_[1], atoms_[2]);
    }
}

void InteractionOfType::sortAtomIds()
{
    if (isBond())
    {
        sortBondAtomIds();
    }
    else if (isAngle())
    {
        sortAngleAtomIds();
    }
    else if (isDihedral())
    {
        sortDihedralAtomIds();
    }
    else
    {
        GMX_THROW(gmx::InternalError(
                "Can not sort parameters other than bonds, angles or dihedrals"));
    }
};

void InteractionOfType::setForceParameter(int pos, real value)
{
    GMX_RELEASE_ASSERT(pos < MAXFORCEPARAM,
                       "Can't set parameter beyond the maximum number of parameters");
    forceParam_[pos] = value;
}

void MoleculeInformation::initMolInfo()
{
    init_block(&mols);
    excls.clear();
    init_t_atoms(&atoms, 0, FALSE);
}

void MoleculeInformation::partialCleanUp()
{
    done_block(&mols);
}

void MoleculeInformation::fullCleanUp()
{
    done_atom(&atoms);
    done_block(&mols);
}

static int rm_interactions(int ifunc, gmx::ArrayRef<MoleculeInformation> mols)
{
    int n = 0;
    /* For all the molecule types */
    for (auto& mol : mols)
    {
        n += mol.interactions[ifunc].size();
        mol.interactions[ifunc].interactionTypes.clear();
    }
    return n;
}

static int check_atom_names(const char*          fn1,
                            const char*          fn2,
                            gmx_mtop_t*          mtop,
                            const t_atoms*       at,
                            const gmx::MDLogger& logger)
{
    int      m, i, j, nmismatch;
    t_atoms* tat;

    constexpr int c_maxNumberOfMismatches = 20;

    if (mtop->natoms != at->nr)
    {
        gmx_incons("comparing atom names");
    }

    nmismatch = 0;
    i         = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        tat = &mtop->moltype[molb.type].atoms;
        for (m = 0; m < molb.nmol; m++)
        {
            for (j = 0; j < tat->nr; j++)
            {
                if (strcmp(*(tat->atomname[j]), *(at->atomname[i])) != 0)
                {
                    if (nmismatch < c_maxNumberOfMismatches)
                    {
                        GMX_LOG(logger.warning)
                                .asParagraph()
                                .appendTextFormatted(
                                        "atom name %d in %s and %s does not match (%s - %s)",
                                        i + 1,
                                        fn1,
                                        fn2,
                                        *(tat->atomname[j]),
                                        *(at->atomname[i]));
                    }
                    else if (nmismatch == c_maxNumberOfMismatches)
                    {
                        GMX_LOG(logger.warning)
                                .asParagraph()
                                .appendTextFormatted("(more than %d non-matching atom names)",
                                                     c_maxNumberOfMismatches);
                    }
                    nmismatch++;
                }
                i++;
            }
        }
    }

    return nmismatch;
}

static void check_bonds_timestep(const gmx_mtop_t* mtop, double dt, WarningHandler* wi)
{
    /* This check is not intended to ensure accurate integration,
     * rather it is to signal mistakes in the mdp settings.
     * A common mistake is to forget to turn on constraints
     * for MD after energy minimization with flexible bonds.
     * This check can also detect too large time steps for flexible water
     * models, but such errors will often be masked by the constraints
     * mdp options, which turns flexible water into water with bond constraints,
     * but without an angle constraint. Unfortunately such incorrect use
     * of water models can not easily be detected without checking
     * for specific model names.
     *
     * The stability limit of leap-frog or velocity verlet is 4.44 steps
     * per oscillational period.
     * But accurate bonds distributions are lost far before that limit.
     * To allow relatively common schemes (although not common with Gromacs)
     * of dt=1 fs without constraints and dt=2 fs with only H-bond constraints
     * we set the note limit to 10.
     */
    int  min_steps_warn = 5;
    int  min_steps_note = 10;
    int  ftype;
    int  i, a1, a2, w_a1, w_a2, j;
    real twopi2, limit2, fc, re, m1, m2, period2, w_period2;
    bool bFound, bWater, bWarn;

    /* Get the interaction parameters */
    gmx::ArrayRef<const t_iparams> ip = mtop->ffparams.iparams;

    twopi2 = gmx::square(2 * M_PI);

    limit2 = gmx::square(min_steps_note * dt);

    w_a1 = w_a2 = -1;
    w_period2   = -1.0;

    const gmx_moltype_t* w_moltype = nullptr;
    for (const gmx_moltype_t& moltype : mtop->moltype)
    {
        const t_atom*           atom  = moltype.atoms.atom;
        const InteractionLists& ilist = moltype.ilist;
        const InteractionList&  ilc   = ilist[F_CONSTR];
        const InteractionList&  ils   = ilist[F_SETTLE];
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (!(ftype == F_BONDS || ftype == F_G96BONDS || ftype == F_HARMONIC))
            {
                continue;
            }

            const InteractionList& ilb = ilist[ftype];
            for (i = 0; i < ilb.size(); i += 3)
            {
                fc = ip[ilb.iatoms[i]].harmonic.krA;
                re = ip[ilb.iatoms[i]].harmonic.rA;
                if (ftype == F_G96BONDS)
                {
                    /* Convert squared sqaure fc to harmonic fc */
                    fc = 2 * fc * re;
                }
                a1 = ilb.iatoms[i + 1];
                a2 = ilb.iatoms[i + 2];
                m1 = atom[a1].m;
                m2 = atom[a2].m;
                if (fc > 0 && m1 > 0 && m2 > 0)
                {
                    period2 = twopi2 * m1 * m2 / ((m1 + m2) * fc);
                }
                else
                {
                    period2 = GMX_FLOAT_MAX;
                }
                if (debug)
                {
                    fprintf(debug, "fc %g m1 %g m2 %g period %g\n", fc, m1, m2, std::sqrt(period2));
                }
                if (period2 < limit2)
                {
                    bFound = false;
                    for (j = 0; j < ilc.size(); j += 3)
                    {
                        if ((ilc.iatoms[j + 1] == a1 && ilc.iatoms[j + 2] == a2)
                            || (ilc.iatoms[j + 1] == a2 && ilc.iatoms[j + 2] == a1))
                        {
                            bFound = true;
                        }
                    }
                    for (j = 0; j < ils.size(); j += 4)
                    {
                        if ((a1 == ils.iatoms[j + 1] || a1 == ils.iatoms[j + 2] || a1 == ils.iatoms[j + 3])
                            && (a2 == ils.iatoms[j + 1] || a2 == ils.iatoms[j + 2]
                                || a2 == ils.iatoms[j + 3]))
                        {
                            bFound = true;
                        }
                    }
                    if (!bFound && (w_moltype == nullptr || period2 < w_period2))
                    {
                        w_moltype = &moltype;
                        w_a1      = a1;
                        w_a2      = a2;
                        w_period2 = period2;
                    }
                }
            }
        }
    }

    if (w_moltype != nullptr)
    {
        bWarn = (w_period2 < gmx::square(min_steps_warn * dt));
        /* A check that would recognize most water models */
        bWater = ((*w_moltype->atoms.atomname[0])[0] == 'O' && w_moltype->atoms.nr <= 5);
        std::string warningMessage = gmx::formatString(
                "The bond in molecule-type %s between atoms %d %s and %d %s has an estimated "
                "oscillational period of %.1e ps, which is less than %d times the time step of "
                "%.1e ps.\n"
                "%s",
                *w_moltype->name,
                w_a1 + 1,
                *w_moltype->atoms.atomname[w_a1],
                w_a2 + 1,
                *w_moltype->atoms.atomname[w_a2],
                std::sqrt(w_period2),
                bWarn ? min_steps_warn : min_steps_note,
                dt,
                bWater ? "Maybe you asked for fexible water."
                       : "Maybe you forgot to change the constraints mdp option.");
        if (bWarn)
        {
            wi->addWarning(warningMessage);
        }
        else
        {
            wi->addNote(warningMessage);
        }
    }
}

static void check_vel(gmx_mtop_t* mtop, rvec v[])
{
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        if (local.ptype == ParticleType::Shell || local.ptype == ParticleType::Bond
            || local.ptype == ParticleType::VSite)
        {
            clear_rvec(v[i]);
        }
    }
}

static void check_shells_inputrec(gmx_mtop_t* mtop, t_inputrec* ir, WarningHandler* wi)
{
    int nshells = 0;

    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        if (local.ptype == ParticleType::Shell || local.ptype == ParticleType::Bond)
        {
            nshells++;
        }
    }
    if ((nshells > 0) && (ir->nstcalcenergy != 1))
    {
        wi->setFileAndLineNumber("unknown", -1);
        std::string warningMessage = gmx::formatString(
                "There are %d shells, changing nstcalcenergy from %d to 1", nshells, ir->nstcalcenergy);
        ir->nstcalcenergy = 1;
        wi->addWarning(warningMessage);
    }
}

/* TODO Decide whether this function can be consolidated with
 * gmx_mtop_ftype_count */
static int nint_ftype(gmx_mtop_t* mtop, gmx::ArrayRef<const MoleculeInformation> mi, int ftype)
{
    int nint = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        nint += molb.nmol * mi[molb.type].interactions[ftype].size();
    }

    return nint;
}

/* This routine reorders the molecule type array
 * in the order of use in the molblocks,
 * unused molecule types are deleted.
 */
static void renumber_moltypes(gmx_mtop_t* sys, std::vector<MoleculeInformation>* molinfo)
{

    std::vector<int> order;
    for (gmx_molblock_t& molblock : sys->molblock)
    {
        const auto found = std::find_if(order.begin(), order.end(), [&molblock](const auto& entry) {
            return molblock.type == entry;
        });
        if (found == order.end())
        {
            /* This type did not occur yet, add it */
            order.push_back(molblock.type);
            molblock.type = order.size() - 1;
        }
        else
        {
            molblock.type = std::distance(order.begin(), found);
        }
    }

    /* We still need to reorder the molinfo structs */
    std::vector<MoleculeInformation> minew(order.size());
    int                              index = 0;
    for (auto& mi : *molinfo)
    {
        const auto found = std::find(order.begin(), order.end(), index);
        if (found != order.end())
        {
            int pos    = std::distance(order.begin(), found);
            minew[pos] = mi;
        }
        else
        {
            // Need to manually clean up memory ....
            mi.fullCleanUp();
        }
        index++;
    }

    *molinfo = minew;
}

static void molinfo2mtop(gmx::ArrayRef<const MoleculeInformation> mi, gmx_mtop_t* mtop)
{
    mtop->moltype.resize(mi.size());
    int pos = 0;
    for (const auto& mol : mi)
    {
        gmx_moltype_t& molt = mtop->moltype[pos];
        molt.name           = mol.name;
        molt.atoms          = mol.atoms;
        /* ilists are copied later */
        molt.excls = mol.excls;
        pos++;
    }
}

static void new_status(const char*                                 topfile,
                       const std::optional<std::filesystem::path>& topppfile,
                       const char*                                 confin,
                       t_gromppopts*                               opts,
                       t_inputrec*                                 ir,
                       gmx_bool                                    bZero,
                       bool                                        bGenVel,
                       bool                                        bVerbose,
                       t_state*                                    state,
                       PreprocessingAtomTypes*                     atypes,
                       gmx_mtop_t*                                 sys,
                       std::vector<MoleculeInformation>*           mi,
                       std::unique_ptr<MoleculeInformation>*       intermolecular_interactions,
                       gmx::ArrayRef<InteractionsOfType>           interactions,
                       CombinationRule*                            comb,
                       double*                                     reppow,
                       real*                                       fudgeQQ,
                       gmx_bool                                    bMorse,
                       WarningHandler*                             wi,
                       const gmx::MDLogger&                        logger)
{
    std::vector<gmx_molblock_t> molblock;
    int                         i, nmismatch;
    bool                        ffParametrizedWithHBondConstraints;

    /* TOPOLOGY processing */
    sys->name = do_top(bVerbose,
                       topfile,
                       topppfile,
                       opts,
                       bZero,
                       &(sys->symtab),
                       interactions,
                       comb,
                       reppow,
                       fudgeQQ,
                       atypes,
                       mi,
                       intermolecular_interactions,
                       ir,
                       &molblock,
                       &ffParametrizedWithHBondConstraints,
                       wi,
                       logger);

    sys->molblock.clear();

    sys->natoms = 0;
    for (const gmx_molblock_t& molb : molblock)
    {
        if (!sys->molblock.empty() && molb.type == sys->molblock.back().type)
        {
            /* Merge consecutive blocks with the same molecule type */
            sys->molblock.back().nmol += molb.nmol;
        }
        else if (molb.nmol > 0)
        {
            /* Add a new molblock to the topology */
            sys->molblock.push_back(molb);
        }
        sys->natoms += molb.nmol * (*mi)[sys->molblock.back().type].atoms.nr;
    }
    if (sys->molblock.empty())
    {
        gmx_fatal(FARGS, "No molecules were defined in the system");
    }

    renumber_moltypes(sys, mi);

    if (bMorse)
    {
        convert_harmonics(*mi, atypes);
    }

    if (ir->eDisre == DistanceRestraintRefinement::None)
    {
        i = rm_interactions(F_DISRES, *mi);
        if (i > 0)
        {
            wi->setFileAndLineNumber("unknown", -1);
            std::string warningMessage =
                    gmx::formatString("disre = no, removed %d distance restraints", i);
            wi->addNote(warningMessage);
        }
    }
    if (!opts->bOrire)
    {
        i = rm_interactions(F_ORIRES, *mi);
        if (i > 0)
        {
            wi->setFileAndLineNumber("unknown", -1);
            std::string warningMessage =
                    gmx::formatString("orire = no, removed %d orientation restraints", i);
            wi->addWarning(warningMessage);
        }
    }

    /* Copy structures from msys to sys */
    molinfo2mtop(*mi, sys);

    sys->finalize();

    /* Discourage using the, previously common, setup of applying constraints
     * to all bonds with force fields that have been parametrized with H-bond
     * constraints only. Do not print note with large timesteps or vsites.
     */
    if (opts->nshake == eshALLBONDS && ffParametrizedWithHBondConstraints && ir->delta_t < 0.0026
        && gmx_mtop_ftype_count(*sys, F_VSITE3FD) == 0)
    {
        wi->setFileAndLineNumber("unknown", -1);
        wi->addNote(
                "You are using constraints on all bonds, whereas the forcefield "
                "has been parametrized only with constraints involving hydrogen atoms. "
                "We suggest using constraints = h-bonds instead, this will also improve "
                "performance.");
    }

    /* COORDINATE file processing */
    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("processing coordinates...");
    }

    t_topology* conftop;
    rvec*       x = nullptr;
    rvec*       v = nullptr;
    snew(conftop, 1);
    // Note that all components in v are set to zero when no v is present in confin
    read_tps_conf(confin, conftop, nullptr, &x, EI_DYNAMICS(ir->eI) ? &v : nullptr, state->box, FALSE);
    if (conftop->atoms.nr != sys->natoms)
    {
        gmx_fatal(FARGS,
                  "number of coordinates in coordinate file (%s, %d)\n"
                  "             does not match topology (%s, %d)",
                  confin,
                  conftop->atoms.nr,
                  topfile,
                  sys->natoms);
    }
    /* It would be nice to get rid of the copies below, but we don't know
     * a priori if the number of atoms in confin matches what we expect.
     */
    state->addEntry(StateEntry::X);
    if (EI_DYNAMICS(ir->eI))
    {
        state->addEntry(StateEntry::V);
    }
    state->changeNumAtoms(sys->natoms);
    std::copy(x, x + state->numAtoms(), state->x.data());
    sfree(x);
    if (EI_DYNAMICS(ir->eI))
    {
        GMX_RELEASE_ASSERT(v, "With dynamics we expect a velocity vector");
        std::copy(v, v + state->numAtoms(), state->v.data());
        sfree(v);
    }
    /* This call fixes the box shape for runs with pressure scaling */
    set_box_rel(ir, state);

    nmismatch = check_atom_names(topfile, confin, sys, &conftop->atoms, logger);
    done_top(conftop);
    sfree(conftop);

    if (nmismatch)
    {
        std::string warningMessage = gmx::formatString(
                "%d non-matching atom name%s\n"
                "atom names from %s will be used\n"
                "atom names from %s will be ignored\n",
                nmismatch,
                (nmismatch == 1) ? "" : "s",
                topfile,
                confin);
        wi->addWarning(warningMessage);
    }

    /* Do more checks, mostly related to constraints */
    if (bVerbose)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("double-checking input for internal consistency...");
    }
    {
        bool bHasNormalConstraints =
                0 < (nint_ftype(sys, *mi, F_CONSTR) + nint_ftype(sys, *mi, F_CONSTRNC));
        bool bHasAnyConstraints = bHasNormalConstraints || 0 < nint_ftype(sys, *mi, F_SETTLE);
        double_check(ir, state->box, bHasNormalConstraints, bHasAnyConstraints, wi);
    }

    if (bGenVel)
    {
        std::vector<real> mass(state->numAtoms());

        for (const AtomProxy atomP : AtomRange(*sys))
        {
            const t_atom& local = atomP.atom();
            int           i     = atomP.globalAtomNumber();
            mass[i]             = local.m;
        }

        if (opts->bMadeSeed)
        {
            GMX_LOG(logger.info).asParagraph().appendTextFormatted("Setting gen_seed to %d", opts->seed);
        }
        GMX_RELEASE_ASSERT(state->hasEntry(StateEntry::V),
                           "Generate velocities only makes sense when they are used");
        maxwell_speed(opts->tempi, opts->seed, sys, state->v.rvec_array(), logger);

        stop_cm(logger, state->numAtoms(), mass.data(), state->x.rvec_array(), state->v.rvec_array());
    }
}

static void copy_state(const char* slog, t_trxframe* fr, bool bReadVel, t_state* state, double* use_time)
{
    if (fr->not_ok & FRAME_NOT_OK)
    {
        gmx_fatal(FARGS, "Can not start from an incomplete frame");
    }
    if (!fr->bX)
    {
        gmx_fatal(FARGS, "Did not find a frame with coordinates in file %s", slog);
    }

    std::copy(fr->x, fr->x + state->numAtoms(), state->x.data());
    if (bReadVel)
    {
        if (!fr->bV)
        {
            gmx_incons("Trajecory frame unexpectedly does not contain velocities");
        }
        std::copy(fr->v, fr->v + state->numAtoms(), state->v.data());
    }
    if (fr->bBox)
    {
        copy_mat(fr->box, state->box);
    }

    *use_time = fr->time;
}

static void cont_status(const char*                                 slog,
                        const std::optional<std::filesystem::path>& ener,
                        bool                                        bNeedVel,
                        bool                                        bGenVel,
                        real                                        fr_time,
                        t_inputrec*                                 ir,
                        t_state*                                    state,
                        gmx_mtop_t*                                 sys,
                        const gmx_output_env_t*                     oenv,
                        const gmx::MDLogger&                        logger)
/* If fr_time == -1 read the last frame available which is complete */
{
    bool         bReadVel;
    t_trxframe   fr;
    t_trxstatus* fp;
    double       use_time;

    bReadVel = (bNeedVel && !bGenVel);

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted("Reading Coordinates%s and Box size from old trajectory",
                                 bReadVel ? ", Velocities" : "");
    if (fr_time == -1)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("Will read whole trajectory");
    }
    else
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("Will read till time %g", fr_time);
    }
    if (!bReadVel)
    {
        if (bGenVel)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted(
                            "Velocities generated: "
                            "ignoring velocities in input trajectory");
        }
        read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X);
    }
    else
    {
        read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X | TRX_NEED_V);

        if (!fr.bV)
        {
            GMX_LOG(logger.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "WARNING: Did not find a frame with velocities in file %s,\n"
                            "         all velocities will be set to zero!",
                            slog);
            for (auto& vi : makeArrayRef(state->v))
            {
                vi = { 0, 0, 0 };
            }
            close_trx(fp);
            /* Search for a frame without velocities */
            bReadVel = false;
            read_first_frame(oenv, &fp, slog, &fr, TRX_NEED_X);
        }
    }

    if (sys->natoms != fr.natoms)
    {
        gmx_fatal(FARGS,
                  "Number of atoms in Topology "
                  "is not the same as in Trajectory");
    }
    state->changeNumAtoms(sys->natoms);
    copy_state(slog, &fr, bReadVel, state, &use_time);

    /* Find the appropriate frame */
    while ((fr_time == -1 || fr.time < fr_time) && read_next_frame(oenv, fp, &fr))
    {
        copy_state(slog, &fr, bReadVel, state, &use_time);
    }

    close_trx(fp);

    /* Set the relative box lengths for preserving the box shape.
     * Note that this call can lead to differences in the last bit
     * with respect to using gmx convert-tpr to create a [REF].tpx[ref] file.
     */
    set_box_rel(ir, state);

    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Using frame at t = %g ps", use_time);
    GMX_LOG(logger.info).asParagraph().appendTextFormatted("Starting time for run is %g ps", ir->init_t);

    if ((ir->pressureCouplingOptions.epc != PressureCoupling::No || ir->etc == TemperatureCoupling::NoseHoover)
        && ener)
    {
        get_enx_state(ener.value(), use_time, sys->groups, ir, state);
        preserveBoxShape(ir->pressureCouplingOptions, ir->deform, state->box_rel, state->boxv);
    }
}

static void read_posres(gmx_mtop_t*                              mtop,
                        gmx::ArrayRef<const MoleculeInformation> molinfo,
                        gmx_bool                                 bTopB,
                        const char*                              fn,
                        RefCoordScaling                          rc_scaling,
                        PbcType                                  pbcType,
                        rvec                                     com,
                        WarningHandler*                          wi,
                        const gmx::MDLogger&                     logger)
{
    gmx_bool*   hadAtom;
    rvec *      x, *v;
    dvec        sum;
    double      totmass;
    t_topology* top;
    matrix      box, invbox;
    int         natoms, npbcdim = 0;
    int         a, nat_molb;
    t_atom*     atom;

    snew(top, 1);
    read_tps_conf(fn, top, nullptr, &x, &v, box, FALSE);
    natoms = top->atoms.nr;
    done_top(top);
    sfree(top);
    if (natoms != mtop->natoms)
    {
        std::string warningMessage = gmx::formatString(
                "The number of atoms in %s (%d) does not match the number of atoms in the topology "
                "(%d). Will assume that the first %d atoms in the topology and %s match.",
                fn,
                natoms,
                mtop->natoms,
                std::min(mtop->natoms, natoms),
                fn);
        wi->addWarning(warningMessage);
    }

    npbcdim = numPbcDimensions(pbcType);
    GMX_RELEASE_ASSERT(npbcdim <= DIM, "Invalid npbcdim");
    clear_rvec(com);
    if (rc_scaling != RefCoordScaling::No)
    {
        copy_mat(box, invbox);
        for (int j = npbcdim; j < DIM; j++)
        {
            clear_rvec(invbox[j]);
            invbox[j][j] = 1;
        }
        gmx::invertBoxMatrix(invbox, invbox);
    }

    /* Copy the reference coordinates to mtop */
    clear_dvec(sum);
    totmass = 0;
    a       = 0;
    snew(hadAtom, natoms);
    for (gmx_molblock_t& molb : mtop->molblock)
    {
        nat_molb                       = molb.nmol * mtop->moltype[molb.type].atoms.nr;
        const InteractionsOfType* pr   = &(molinfo[molb.type].interactions[F_POSRES]);
        const InteractionsOfType* prfb = &(molinfo[molb.type].interactions[F_FBPOSRES]);
        if (pr->size() > 0 || prfb->size() > 0)
        {
            atom = mtop->moltype[molb.type].atoms.atom;
            for (const auto& restraint : pr->interactionTypes)
            {
                int ai = restraint.ai();
                if (ai >= natoms)
                {
                    gmx_fatal(FARGS,
                              "Position restraint atom index (%d) in moltype '%s' is larger than "
                              "number of atoms in %s (%d).\n",
                              ai + 1,
                              *molinfo[molb.type].name,
                              fn,
                              natoms);
                }
                hadAtom[ai] = TRUE;
                if (rc_scaling == RefCoordScaling::Com)
                {
                    /* Determine the center of mass of the posres reference coordinates */
                    for (int j = 0; j < npbcdim; j++)
                    {
                        sum[j] += atom[ai].m * x[a + ai][j];
                    }
                    totmass += atom[ai].m;
                }
            }
            /* Same for flat-bottomed posres, but do not count an atom twice for COM */
            for (const auto& restraint : prfb->interactionTypes)
            {
                int ai = restraint.ai();
                if (ai >= natoms)
                {
                    gmx_fatal(FARGS,
                              "Position restraint atom index (%d) in moltype '%s' is larger than "
                              "number of atoms in %s (%d).\n",
                              ai + 1,
                              *molinfo[molb.type].name,
                              fn,
                              natoms);
                }
                if (rc_scaling == RefCoordScaling::Com && !hadAtom[ai])
                {
                    /* Determine the center of mass of the posres reference coordinates */
                    for (int j = 0; j < npbcdim; j++)
                    {
                        sum[j] += atom[ai].m * x[a + ai][j];
                    }
                    totmass += atom[ai].m;
                }
            }
            if (!bTopB)
            {
                molb.posres_xA.resize(nat_molb);
                for (int i = 0; i < nat_molb; i++)
                {
                    copy_rvec(x[a + i], molb.posres_xA[i]);
                }
            }
            else
            {
                molb.posres_xB.resize(nat_molb);
                for (int i = 0; i < nat_molb; i++)
                {
                    copy_rvec(x[a + i], molb.posres_xB[i]);
                }
            }
        }
        a += nat_molb;
    }
    if (rc_scaling == RefCoordScaling::Com)
    {
        if (totmass == 0)
        {
            gmx_fatal(FARGS, "The total mass of the position restraint atoms is 0");
        }
        for (int j = 0; j < npbcdim; j++)
        {
            com[j] = sum[j] / totmass;
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "The center of mass of the position restraint coord's is %6.3f %6.3f %6.3f",
                        com[XX],
                        com[YY],
                        com[ZZ]);
    }

    if (rc_scaling != RefCoordScaling::No)
    {
        GMX_ASSERT(npbcdim <= DIM, "Only DIM dimensions can have PBC");

        for (gmx_molblock_t& molb : mtop->molblock)
        {
            nat_molb = molb.nmol * mtop->moltype[molb.type].atoms.nr;
            if (!molb.posres_xA.empty() || !molb.posres_xB.empty())
            {
                std::vector<gmx::RVec>& xp = (!bTopB ? molb.posres_xA : molb.posres_xB);
                for (int i = 0; i < nat_molb; i++)
                {
                    for (int j = 0; j < npbcdim; j++)
                    {
                        if (rc_scaling == RefCoordScaling::All)
                        {
                            /* Convert from Cartesian to crystal coordinates */
                            xp[i][j] *= invbox[j][j];
                            for (int k = j + 1; k < npbcdim; k++)
                            {
                                xp[i][j] += invbox[k][j] * xp[i][k];
                            }
                        }
                        else if (rc_scaling == RefCoordScaling::Com)
                        {
                            /* Subtract the center of mass */
                            xp[i][j] -= com[j];
                        }
                    }
                }
            }
        }

        if (rc_scaling == RefCoordScaling::Com)
        {
            /* Convert the COM from Cartesian to crystal coordinates */
            for (int j = 0; j < npbcdim; j++)
            {
                com[j] *= invbox[j][j];
                for (int k = j + 1; k < npbcdim; k++)
                {
                    com[j] += invbox[k][j] * com[k];
                }
            }
        }
    }

    sfree(x);
    sfree(v);
    sfree(hadAtom);
}

static void gen_posres(gmx_mtop_t*                              mtop,
                       gmx::ArrayRef<const MoleculeInformation> mi,
                       const char*                              fnA,
                       const char*                              fnB,
                       RefCoordScaling                          rc_scaling,
                       PbcType                                  pbcType,
                       rvec                                     com,
                       rvec                                     comB,
                       WarningHandler*                          wi,
                       const gmx::MDLogger&                     logger)
{
    read_posres(mtop, mi, FALSE, fnA, rc_scaling, pbcType, com, wi, logger);
    /* It is safer to simply read the b-state posres rather than trying
     * to be smart and copy the positions.
     */
    read_posres(mtop, mi, TRUE, fnB, rc_scaling, pbcType, comB, wi, logger);
}

static void set_wall_atomtype(PreprocessingAtomTypes* at,
                              t_gromppopts*           opts,
                              t_inputrec*             ir,
                              WarningHandler*         wi,
                              const gmx::MDLogger&    logger)
{
    int i;

    if (ir->nwall > 0)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("Searching the wall atom type(s)");
    }
    for (i = 0; i < ir->nwall; i++)
    {
        auto atomType = at->atomTypeFromName(opts->wall_atomtype[i]);
        if (!atomType.has_value())
        {
            std::string warningMessage = gmx::formatString(
                    "Specified wall atom type %s is not defined", opts->wall_atomtype[i]);
            wi->addError(warningMessage);
        }
        else
        {
            ir->wall_atomtype[i] = *atomType;
        }
    }
}

static int nrdf_internal(const t_atoms* atoms)
{
    int i, nmass, nrdf;

    nmass = 0;
    for (i = 0; i < atoms->nr; i++)
    {
        /* Vsite ptype might not be set here yet, so also check the mass */
        if ((atoms->atom[i].ptype == ParticleType::Atom || atoms->atom[i].ptype == ParticleType::Nucleus)
            && atoms->atom[i].m > 0)
        {
            nmass++;
        }
    }
    switch (nmass)
    {
        case 0: // Fall through intended
        case 1: nrdf = 0; break;
        case 2: nrdf = 1; break;
        default: nrdf = nmass * 3 - 6; break;
    }

    return nrdf;
}

static void spline1d(double dx, const double* y, int n, double* u, double* y2)
{
    int    i;
    double p, q;

    y2[0] = 0.0;
    u[0]  = 0.0;

    for (i = 1; i < n - 1; i++)
    {
        p     = 0.5 * y2[i - 1] + 2.0;
        y2[i] = -0.5 / p;
        q     = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx;
        u[i]  = (3.0 * q / dx - 0.5 * u[i - 1]) / p;
    }

    y2[n - 1] = 0.0;

    for (i = n - 2; i >= 0; i--)
    {
        y2[i] = y2[i] * y2[i + 1] + u[i];
    }
}


static void
interpolate1d(double xmin, double dx, const double* ya, const double* y2a, double x, double* y, double* y1)
{
    int    ix;
    double a, b;

    ix = static_cast<int>((x - xmin) / dx);

    a = (xmin + (ix + 1) * dx - x) / dx;
    b = (x - xmin - ix * dx) / dx;

    *y = a * ya[ix] + b * ya[ix + 1]
         + ((a * a * a - a) * y2a[ix] + (b * b * b - b) * y2a[ix + 1]) * (dx * dx) / 6.0;
    *y1 = (ya[ix + 1] - ya[ix]) / dx - (3.0 * a * a - 1.0) / 6.0 * dx * y2a[ix]
          + (3.0 * b * b - 1.0) / 6.0 * dx * y2a[ix + 1];
}


static void setup_cmap(int grid_spacing, int nc, gmx::ArrayRef<const real> grid, gmx_cmap_t* cmap_grid)
{
    int    i, j, k, ii, jj, kk, idx;
    int    offset;
    double dx, xmin, v, v1, v2, v12;
    double phi, psi;

    std::vector<double> tmp_u(2 * grid_spacing, 0.0);
    std::vector<double> tmp_u2(2 * grid_spacing, 0.0);
    std::vector<double> tmp_yy(2 * grid_spacing, 0.0);
    std::vector<double> tmp_y1(2 * grid_spacing, 0.0);
    std::vector<double> tmp_t2(2 * grid_spacing * 2 * grid_spacing, 0.0);
    std::vector<double> tmp_grid(2 * grid_spacing * 2 * grid_spacing, 0.0);

    dx   = 360.0 / grid_spacing;
    xmin = -180.0 - dx * grid_spacing / 2;

    for (kk = 0; kk < nc; kk++)
    {
        /* Compute an offset depending on which cmap we are using
         * Offset will be the map number multiplied with the
         * grid_spacing * grid_spacing * 2
         */
        offset = kk * grid_spacing * grid_spacing * 2;

        for (i = 0; i < 2 * grid_spacing; i++)
        {
            ii = (i + grid_spacing - grid_spacing / 2) % grid_spacing;

            for (j = 0; j < 2 * grid_spacing; j++)
            {
                jj = (j + grid_spacing - grid_spacing / 2) % grid_spacing;
                tmp_grid[i * grid_spacing * 2 + j] = grid[offset + ii * grid_spacing + jj];
            }
        }

        for (i = 0; i < 2 * grid_spacing; i++)
        {
            spline1d(dx,
                     &(tmp_grid[2 * grid_spacing * i]),
                     2 * grid_spacing,
                     tmp_u.data(),
                     &(tmp_t2[2 * grid_spacing * i]));
        }

        for (i = grid_spacing / 2; i < grid_spacing + grid_spacing / 2; i++)
        {
            ii  = i - grid_spacing / 2;
            phi = ii * dx - 180.0;

            for (j = grid_spacing / 2; j < grid_spacing + grid_spacing / 2; j++)
            {
                jj  = j - grid_spacing / 2;
                psi = jj * dx - 180.0;

                for (k = 0; k < 2 * grid_spacing; k++)
                {
                    interpolate1d(xmin,
                                  dx,
                                  &(tmp_grid[2 * grid_spacing * k]),
                                  &(tmp_t2[2 * grid_spacing * k]),
                                  psi,
                                  &tmp_yy[k],
                                  &tmp_y1[k]);
                }

                spline1d(dx, tmp_yy.data(), 2 * grid_spacing, tmp_u.data(), tmp_u2.data());
                interpolate1d(xmin, dx, tmp_yy.data(), tmp_u2.data(), phi, &v, &v1);
                spline1d(dx, tmp_y1.data(), 2 * grid_spacing, tmp_u.data(), tmp_u2.data());
                interpolate1d(xmin, dx, tmp_y1.data(), tmp_u2.data(), phi, &v2, &v12);

                idx                                       = ii * grid_spacing + jj;
                cmap_grid->cmapdata[kk].cmap[idx * 4]     = grid[offset + ii * grid_spacing + jj];
                cmap_grid->cmapdata[kk].cmap[idx * 4 + 1] = v1;
                cmap_grid->cmapdata[kk].cmap[idx * 4 + 2] = v2;
                cmap_grid->cmapdata[kk].cmap[idx * 4 + 3] = v12;
            }
        }
    }
}

static void init_cmap_grid(gmx_cmap_t* cmap_grid, int ngrid, int grid_spacing)
{
    int i, nelem;

    cmap_grid->grid_spacing = grid_spacing;
    nelem                   = cmap_grid->grid_spacing * cmap_grid->grid_spacing;

    cmap_grid->cmapdata.resize(ngrid);

    for (i = 0; i < ngrid; i++)
    {
        cmap_grid->cmapdata[i].cmap.resize(4 * nelem);
    }
}


static int count_constraints(const gmx_mtop_t*                        mtop,
                             gmx::ArrayRef<const MoleculeInformation> mi,
                             WarningHandler*                          wi)
{
    int count, count_mol;

    count = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        count_mol                                            = 0;
        gmx::ArrayRef<const InteractionsOfType> interactions = mi[molb.type].interactions;

        for (int i = 0; i < F_NRE; i++)
        {
            if (i == F_SETTLE)
            {
                count_mol += 3 * interactions[i].size();
            }
            else if (interaction_function[i].flags & IF_CONSTRAINT)
            {
                count_mol += interactions[i].size();
            }
        }

        if (count_mol > nrdf_internal(&mi[molb.type].atoms))
        {
            std::string warningMessage = gmx::formatString(
                    "Molecule type '%s' has %d constraints.\n"
                    "For stability and efficiency there should not be more constraints than "
                    "internal number of degrees of freedom: %d.\n",
                    *mi[molb.type].name,
                    count_mol,
                    nrdf_internal(&mi[molb.type].atoms));
            wi->addWarning(warningMessage);
        }
        count += molb.nmol * count_mol;
    }

    return count;
}

static real calc_temp(const gmx_mtop_t* mtop, const t_inputrec* ir, rvec* v)
{
    double sum_mv2 = 0;
    for (const AtomProxy atomP : AtomRange(*mtop))
    {
        const t_atom& local = atomP.atom();
        int           i     = atomP.globalAtomNumber();
        sum_mv2 += local.m * norm2(v[i]);
    }

    double nrdf = 0;
    for (int g = 0; g < ir->opts.ngtc; g++)
    {
        nrdf += ir->opts.nrdf[g];
    }

    return sum_mv2 / (nrdf * gmx::c_boltz);
}

static real get_max_reference_temp(const t_inputrec* ir, WarningHandler* wi)
{
    real ref_t;
    int  i;
    bool bNoCoupl;

    ref_t    = 0;
    bNoCoupl = false;
    for (i = 0; i < ir->opts.ngtc; i++)
    {
        if (ir->opts.tau_t[i] < 0)
        {
            bNoCoupl = true;
        }
        else
        {
            ref_t = std::max(ref_t, ir->opts.ref_t[i]);
        }
    }

    if (bNoCoupl)
    {
        std::string warningMessage = gmx::formatString(
                "Some temperature coupling groups do not use temperature coupling. We will assume "
                "their temperature is not more than %.3f K. If their temperature is higher, the "
                "energy error and the Verlet buffer might be underestimated.",
                ref_t);
        wi->addWarning(warningMessage);
    }

    return ref_t;
}

/* Checks if there are unbound atoms in moleculetype molt.
 * Prints a note for each unbound atoms and a warning if any is present.
 */
static void checkForUnboundAtoms(const gmx_moltype_t* molt,
                                 gmx_bool             bVerbose,
                                 WarningHandler*      wi,
                                 const gmx::MDLogger& logger)
{
    const t_atoms* atoms = &molt->atoms;

    if (atoms->nr == 1)
    {
        /* Only one atom, there can't be unbound atoms */
        return;
    }

    std::vector<int> count(atoms->nr, 0);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (((interaction_function[ftype].flags & IF_BOND) && NRAL(ftype) == 2 && ftype != F_CONNBONDS)
            || (interaction_function[ftype].flags & IF_CONSTRAINT) || ftype == F_SETTLE)
        {
            const InteractionList& il   = molt->ilist[ftype];
            const int              nral = NRAL(ftype);

            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                for (int j = 0; j < nral; j++)
                {
                    count[il.iatoms[i + 1 + j]]++;
                }
            }
        }
    }

    int numDanglingAtoms = 0;
    for (int a = 0; a < atoms->nr; a++)
    {
        if (atoms->atom[a].ptype != ParticleType::VSite && count[a] == 0)
        {
            if (bVerbose)
            {
                GMX_LOG(logger.warning)
                        .asParagraph()
                        .appendTextFormatted(
                                "Atom %d '%s' in moleculetype '%s' is not bound by a potential or "
                                "constraint to any other atom in the same moleculetype.",
                                a + 1,
                                *atoms->atomname[a],
                                *molt->name);
            }
            numDanglingAtoms++;
        }
    }

    if (numDanglingAtoms > 0)
    {
        std::string warningMessage = gmx::formatString(
                "In moleculetype '%s' %d atoms are not bound by a potential or constraint to any "
                "other atom in the same moleculetype. Although technically this might not cause "
                "issues in a simulation, this often means that the user forgot to add a "
                "bond/potential/constraint or put multiple molecules in the same moleculetype "
                "definition by mistake. Run with -v to get information for each atom.",
                *molt->name,
                numDanglingAtoms);
        wi->addNote(warningMessage);
    }
}

/* Checks all moleculetypes for unbound atoms */
static void checkForUnboundAtoms(const gmx_mtop_t*    mtop,
                                 gmx_bool             bVerbose,
                                 WarningHandler*      wi,
                                 const gmx::MDLogger& logger)
{
    for (const gmx_moltype_t& molt : mtop->moltype)
    {
        checkForUnboundAtoms(&molt, bVerbose, wi, logger);
    }
}

/*! \brief Checks if there are decoupled modes in moleculetype \p molt.
 *
 * The specific decoupled modes this routine check for are angle modes
 * where the two bonds are constrained and the atoms a both ends are only
 * involved in a single constraint; the mass of the two atoms needs to
 * differ by more than \p massFactorThreshold.
 */
static bool haveDecoupledModeInMol(const gmx_moltype_t&           molt,
                                   gmx::ArrayRef<const t_iparams> iparams,
                                   real                           massFactorThreshold)
{
    if (molt.ilist[F_CONSTR].empty() && molt.ilist[F_CONSTRNC].empty())
    {
        return false;
    }

    const t_atom* atom = molt.atoms.atom;

    const auto atomToConstraints =
            gmx::make_at2con(molt, iparams, gmx::FlexibleConstraintTreatment::Exclude);

    bool haveDecoupledMode = false;
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_ATYPE)
        {
            const int              nral = NRAL(ftype);
            const InteractionList& il   = molt.ilist[ftype];
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                /* Here we check for the mass difference between the atoms
                 * at both ends of the angle, that the atoms at the ends have
                 * 1 contraint and the atom in the middle at least 3; we check
                 * that the 3 atoms are linked by constraints below.
                 * We check for at least three constraints for the middle atom,
                 * since with only the two bonds in the angle, we have 3-atom
                 * molecule, which has much more thermal exhange in this single
                 * angle mode than molecules with more atoms.
                 * Note that this check also catches molecules where more atoms
                 * are connected to one or more atoms in the angle, but by
                 * bond potentials instead of angles. But such cases will not
                 * occur in "normal" molecules and it doesn't hurt running
                 * those with higher accuracy settings as well.
                 */
                int a0 = il.iatoms[1 + i];
                int a1 = il.iatoms[1 + i + 1];
                int a2 = il.iatoms[1 + i + 2];
                if ((atom[a0].m > atom[a2].m * massFactorThreshold || atom[a2].m > atom[a0].m * massFactorThreshold)
                    && atomToConstraints[a0].ssize() == 1 && atomToConstraints[a2].ssize() == 1
                    && atomToConstraints[a1].ssize() >= 3)
                {
                    int constraint0 = atomToConstraints[a0][0];
                    int constraint2 = atomToConstraints[a2][0];

                    bool foundAtom0 = false;
                    bool foundAtom2 = false;
                    for (const int constraint : atomToConstraints[a1])
                    {
                        if (constraint == constraint0)
                        {
                            foundAtom0 = true;
                        }
                        if (constraint == constraint2)
                        {
                            foundAtom2 = true;
                        }
                    }
                    if (foundAtom0 && foundAtom2)
                    {
                        haveDecoupledMode = true;
                    }
                }
            }
        }
    }

    return haveDecoupledMode;
}

/*! \brief Checks if the Verlet buffer and constraint accuracy is sufficient for decoupled dynamic modes.
 *
 * When decoupled modes are present and the accuracy in insufficient,
 * this routine issues a warning if the accuracy is insufficient.
 */
static void checkDecoupledModeAccuracy(const gmx_mtop_t* mtop, const t_inputrec* ir, WarningHandler* wi)
{
    /* We only have issues with decoupled modes with normal MD.
     * With stochastic dynamics equipartitioning is enforced strongly.
     */
    if (!EI_MD(ir->eI))
    {
        return;
    }

    /* When atoms of very different mass are involved in an angle potential
     * and both bonds in the angle are constrained, the dynamic modes in such
     * angles have very different periods and significant energy exchange
     * takes several nanoseconds. Thus even a small amount of error in
     * different algorithms can lead to violation of equipartitioning.
     * The parameters below are mainly based on an all-atom chloroform model
     * with all bonds constrained. Equipartitioning is off by more than 1%
     * (need to run 5-10 ns) when the difference in mass between H and Cl
     * is more than a factor 13 and the accuracy is less than the thresholds
     * given below. This has been verified on some other molecules.
     *
     * Note that the buffer and shake parameters have unit length and
     * energy/time, respectively, so they will "only" work correctly
     * for atomistic force fields using MD units.
     */
    const real massFactorThreshold      = 13.0;
    const real bufferToleranceThreshold = 1e-4;
    const int  lincsIterationThreshold  = 2;
    const int  lincsOrderThreshold      = 4;
    const real shakeToleranceThreshold  = 0.005 * ir->delta_t;

    bool lincsWithSufficientTolerance =
            (ir->eConstrAlg == ConstraintAlgorithm::Lincs
             && ir->nLincsIter >= lincsIterationThreshold && ir->nProjOrder >= lincsOrderThreshold);
    bool shakeWithSufficientTolerance = (ir->eConstrAlg == ConstraintAlgorithm::Shake
                                         && ir->shake_tol <= 1.1 * shakeToleranceThreshold);
    if (ir->cutoff_scheme == CutoffScheme::Verlet && ir->verletbuf_tol <= 1.1 * bufferToleranceThreshold
        && (lincsWithSufficientTolerance || shakeWithSufficientTolerance))
    {
        return;
    }

    bool haveDecoupledMode = false;
    for (const gmx_moltype_t& molt : mtop->moltype)
    {
        if (haveDecoupledModeInMol(molt, mtop->ffparams.iparams, massFactorThreshold))
        {
            haveDecoupledMode = true;
        }
    }

    if (haveDecoupledMode)
    {
        std::string message = gmx::formatString(
                "There are atoms at both ends of an angle, connected by constraints "
                "and with masses that differ by more than a factor of %g. This means "
                "that there are likely dynamic modes that are only very weakly coupled.",
                massFactorThreshold);
        if (ir->cutoff_scheme == CutoffScheme::Verlet)
        {
            message += gmx::formatString(
                    " To ensure good equipartitioning, you need to either not use "
                    "constraints on all bonds (but, if possible, only on bonds involving "
                    "hydrogens) or use integrator = %s or decrease one or more tolerances: "
                    "verlet-buffer-tolerance <= %g, LINCS iterations >= %d, LINCS order "
                    ">= %d or SHAKE tolerance <= %g",
                    enumValueToString(IntegrationAlgorithm::SD1),
                    bufferToleranceThreshold,
                    lincsIterationThreshold,
                    lincsOrderThreshold,
                    shakeToleranceThreshold);
        }
        else
        {
            message += gmx::formatString(
                    " To ensure good equipartitioning, we suggest to switch to the %s "
                    "cutoff-scheme, since that allows for better control over the Verlet "
                    "buffer size and thus over the energy drift.",
                    enumValueToString(CutoffScheme::Verlet));
        }
        wi->addWarning(message);
    }
}

static void set_verlet_buffer(const gmx_mtop_t*              mtop,
                              t_inputrec*                    ir,
                              real                           buffer_temp,
                              gmx::ArrayRef<const gmx::RVec> coordinates,
                              matrix                         box,
                              WarningHandler*                wi,
                              const gmx::MDLogger&           logger)
{
    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "Determining Verlet buffer for a tolerance of %g kJ/mol/ps at %g K",
                    ir->verletbuf_tol,
                    buffer_temp);

    const real effectiveAtomDensity = computeEffectiveAtomDensity(
            coordinates, box, std::max(ir->rcoulomb, ir->rvdw), MPI_COMM_NULL);

    /* Calculate the buffer size for simple atom vs atoms list */
    VerletbufListSetup listSetup1x1;
    listSetup1x1.cluster_size_i = 1;
    listSetup1x1.cluster_size_j = 1;
    const real rlist_1x1        = calcVerletBufferSize(*mtop,
                                                effectiveAtomDensity,
                                                *ir,
                                                ir->verletBufferPressureTolerance,
                                                ir->nstlist,
                                                ir->nstlist - 1,
                                                buffer_temp,
                                                listSetup1x1);

    /* Set the pair-list buffer size in ir */
    VerletbufListSetup listSetup4x4 = verletbufGetSafeListSetup(ListSetupType::CpuNoSimd);
    ir->rlist                       = calcVerletBufferSize(*mtop,
                                     effectiveAtomDensity,
                                     *ir,
                                     ir->verletBufferPressureTolerance,
                                     ir->nstlist,
                                     ir->nstlist - 1,
                                     buffer_temp,
                                     listSetup4x4);

    const int n_nonlin_vsite = gmx::countNonlinearVsites(*mtop);
    if (n_nonlin_vsite > 0)
    {
        std::string warningMessage = gmx::formatString(
                "There are %d non-linear virtual site constructions. Their contribution to the "
                "energy error is approximated. In most cases this does not affect the error "
                "significantly.",
                n_nonlin_vsite);
        wi->addNote(warningMessage);
    }

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "Calculated rlist for %dx%d atom pair-list as %.3f nm, buffer size %.3f nm",
                    1,
                    1,
                    rlist_1x1,
                    rlist_1x1 - std::max(ir->rvdw, ir->rcoulomb));

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "Set rlist, assuming %dx%d atom pair-list, to %.3f nm, buffer size %.3f nm",
                    listSetup4x4.cluster_size_i,
                    listSetup4x4.cluster_size_j,
                    ir->rlist,
                    ir->rlist - std::max(ir->rvdw, ir->rcoulomb));

    GMX_LOG(logger.info)
            .asParagraph()
            .appendTextFormatted(
                    "Note that mdrun will redetermine rlist based on the actual pair-list setup");

    if (gmx::square(ir->rlist) >= max_cutoff2(ir->pbcType, box))
    {
        gmx_fatal(FARGS,
                  "The pair-list cut-off (%g nm) is longer than half the shortest box vector or "
                  "longer than the smallest box diagonal element (%g nm). Increase the box size or "
                  "decrease nstlist or increase verlet-buffer-tolerance.",
                  ir->rlist,
                  std::sqrt(max_cutoff2(ir->pbcType, box)));
    }
}

// Computes and returns that largest distance between non-perturbed excluded atom pairs
static std::tuple<real, int, int> maxNonPerturbedExclusionDistance(const gmx_mtop_t& mtop,
                                                                   const bool        useFep,
                                                                   const PbcType     pbcType,
                                                                   gmx::ArrayRef<const gmx::RVec> x,
                                                                   const matrix box)
{
    t_pbc pbc;

    set_pbc(&pbc, pbcType, box);

    real dx2Max = 0;
    int  atom0  = -1;
    int  atom1  = -1;

    int moleculeOffset = 0;
    for (const auto& mb : mtop.molblock)
    {
        const gmx::ListOfLists<int>& excls = mtop.moltype[mb.type].excls;
        const t_atoms&               atoms = mtop.moltype[mb.type].atoms;
        GMX_RELEASE_ASSERT(gmx::ssize(excls) == atoms.nr,
                           "There should be one exclusion list per atom");

        for (int mol = 0; mol < mb.nmol; mol++)
        {
            for (int iAtom = 0; iAtom < atoms.nr; iAtom++)
            {
                if (useFep && PERTURBED(atoms.atom[iAtom]))
                {
                    // This atom is perturbed, so all its exclusions are perturbed
                    continue;
                }

                const auto& jAtoms = excls[iAtom];
                for (int jAtom : jAtoms)
                {
                    if (jAtom != iAtom && !(useFep && PERTURBED(atoms.atom[jAtom])))
                    {
                        rvec dx;

                        pbc_dx(&pbc, x[moleculeOffset + iAtom], x[moleculeOffset + jAtom], dx);
                        const real distanceSquared = norm2(dx);
                        if (distanceSquared > dx2Max)
                        {
                            dx2Max = distanceSquared;
                            atom0  = moleculeOffset + iAtom;
                            atom1  = moleculeOffset + jAtom;
                        }
                    }
                }
            }

            moleculeOffset += atoms.nr;
        }
    }

    return std::make_tuple(std::sqrt(dx2Max), atom0, atom1);
}

// Computes and logs the maximum exclusion distance. Checks whether (non-perturbed) excluded
// atom pairs are close to the cut-off distance and if so, generates a warning/error
static void checkExclusionDistances(const gmx_mtop_t&              mtop,
                                    const t_inputrec&              ir,
                                    gmx::ArrayRef<const gmx::RVec> x,
                                    const matrix                   box,
                                    const gmx::MDLogger&           logger,
                                    WarningHandler*                wi)
{
    // Check the maximum distance for (non-perturbed) excluded pairs here,
    // as it should not be longer than the cut-off distance, but we can't
    // easily ensure that during the run.
    const bool useFep = (ir.efep != FreeEnergyPerturbationType::No);
    const auto [maxExclusionDistance, atom0, atom1] =
            maxNonPerturbedExclusionDistance(mtop, useFep, ir.pbcType, x, box);

    const std::string distanceString = gmx::formatString(
            "The largest distance between%s excluded atoms is %.3f nm between atom %d and %d",
            useFep ? " non-perturbed" : "",
            maxExclusionDistance,
            atom0 + 1,
            atom1 + 1);

    GMX_LOG(logger.info).asParagraph().appendText(distanceString);

    const real cutoffDistance = std::max(ir.rvdw, ir.rcoulomb);
    if (maxExclusionDistance >= cutoffDistance)
    {
        std::string text = gmx::formatString(
                "%s, which is larger than the cut-off distance. This will "
                "lead to missing long-range corrections in the forces and energies.",
                distanceString.c_str());
        if (EI_ENERGY_MINIMIZATION(ir.eI))
        {
            text += " If you expect that minimization will bring such distances within the "
                    "cut-off, you can ignore this warning.";
            wi->addWarning(text);
        }
        else
        {
            wi->addError(text);
        }
    }
    else if (maxExclusionDistance > 0.9_real * cutoffDistance)
    {
        std::string text = gmx::formatString(
                "%s, which is larger than 90%% of the cut-off distance. "
                "When excluded pairs go beyond the cut-off distance, this leads to missing "
                "long-range corrections in the forces and energies.",
                distanceString.c_str());
        if (EI_DYNAMICS(ir.eI))
        {
            wi->addWarning(text);
        }
        else
        {
            wi->addNote(text);
        }
    }
}

//! Add the velocity profile of \p deform to the velocities in \p state
static void deformInitFlow(t_state* state, const matrix deform)
{
    // Deform gives the speed of box vector elements, we need to scale relative to the box size
    matrix coordToVelocity;
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            coordToVelocity[d1][d2] = deform[d1][d2] / state->box[d1][d1];
        }
    }

    for (int i = 0; i < state->numAtoms(); i++)
    {
        for (int d1 = 0; d1 < DIM; d1++)
        {
            for (int d2 = 0; d2 < DIM; d2++)
            {
                state->v[i][d2] += coordToVelocity[d1][d2] * state->x[i][d1];
            }
        }
    }
}

int gmx_grompp(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] (the gromacs preprocessor)",
        "reads a molecular topology file, checks the validity of the",
        "file, expands the topology from a molecular description to an atomic",
        "description. The topology file contains information about",
        "molecule types and the number of molecules, the preprocessor",
        "copies each molecule as needed. ",
        "There is no limitation on the number of molecule types. ",
        "Bonds and bond-angles can be converted into constraints, separately",
        "for hydrogens and heavy atoms.",
        "Then a coordinate file is read and velocities can be generated",
        "from a Maxwellian distribution if requested.",
        "[THISMODULE] also reads parameters for [gmx-mdrun] ",
        "(eg. number of MD steps, time step, cut-off).",
        "Eventually a binary file is produced that can serve as the sole input",
        "file for the MD program.[PAR]",

        "[THISMODULE] uses the atom names from the topology file. The atom names",
        "in the coordinate file (option [TT]-c[tt]) are only read to generate",
        "warnings when they do not match the atom names in the topology.",
        "Note that the atom names are irrelevant for the simulation as",
        "only the atom types are used for generating interaction parameters.[PAR]",

        "[THISMODULE] uses a built-in preprocessor to resolve includes, macros, ",
        "etc. The preprocessor supports the following keywords::",
        "",
        "    #ifdef VARIABLE",
        "    #ifndef VARIABLE",
        "    #else",
        "    #endif",
        "    #define VARIABLE",
        "    #undef VARIABLE",
        "    #include \"filename\"",
        "    #include <filename>",
        "",
        "The functioning of these statements in your topology may be modulated by",
        "using the following two flags in your [REF].mdp[ref] file::",
        "",
        "    define = -DVARIABLE1 -DVARIABLE2",
        "    include = -I/home/john/doe",
        "",
        "For further information a C-programming textbook may help you out.",
        "Specifying the [TT]-pp[tt] flag will get the pre-processed",
        "topology file written out so that you can verify its contents.[PAR]",

        "When using position restraints, a file with restraint coordinates",
        "must be supplied with [TT]-r[tt] (can be the same file as supplied",
        "for [TT]-c[tt]). For free energy calculations, separate reference",
        "coordinates for the B topology can be supplied with [TT]-rb[tt],",
        "otherwise they will be equal to those of the A topology.[PAR]",

        "Starting coordinates can be read from trajectory with [TT]-t[tt].",
        "The last frame with coordinates and velocities will be read,",
        "unless the [TT]-time[tt] option is used. Only if this information",
        "is absent will the coordinates in the [TT]-c[tt] file be used.",
        "Note that these velocities will not be used when [TT]gen_vel = yes[tt]",
        "in your [REF].mdp[ref] file. An energy file can be supplied with",
        "[TT]-e[tt] to read Nose-Hoover and/or Parrinello-Rahman coupling",
        "variables.[PAR]",

        "[THISMODULE] can be used to restart simulations (preserving",
        "continuity) by supplying just a checkpoint file with [TT]-t[tt].",
        "However, for simply changing the number of run steps to extend",
        "a run, using [gmx-convert-tpr] is more convenient than [THISMODULE].",
        "You then supply the old checkpoint file directly to [gmx-mdrun]",
        "with [TT]-cpi[tt]. If you wish to change the ensemble or things",
        "like output frequency, then supplying the checkpoint file to",
        "[THISMODULE] with [TT]-t[tt] along with a new [REF].mdp[ref] file",
        "with [TT]-f[tt] is the recommended procedure. Actually preserving",
        "the ensemble (if possible) still requires passing the checkpoint",
        "file to [gmx-mdrun] [TT]-cpi[tt].[PAR]",

        "By default, all bonded interactions which have constant energy due to",
        "virtual site constructions will be removed. If this constant energy is",
        "not zero, this will result in a shift in the total energy. All bonded",
        "interactions can be kept by turning off [TT]-rmvsbds[tt]. Additionally,",
        "all constraints for distances which will be constant anyway because",
        "of virtual site constructions will be removed. If any constraints remain",
        "which involve virtual sites, a fatal error will result.[PAR]",

        "To verify your run input file, please take note of all warnings",
        "on the screen, and correct where necessary. Do also look at the contents",
        "of the [TT]mdout.mdp[tt] file; this contains comment lines, as well as",
        "the input that [THISMODULE] has read. If in doubt, you can start [THISMODULE]",
        "with the [TT]-debug[tt] option which will give you more information",
        "in a file called [TT]grompp.log[tt] (along with real debug info). You",
        "can see the contents of the run input file with the [gmx-dump]",
        "program. [gmx-check] can be used to compare the contents of two",
        "run input files.[PAR]",

        "The [TT]-maxwarn[tt] option can be used to override warnings printed",
        "by [THISMODULE] that otherwise halt output. In some cases, warnings are",
        "harmless, but usually they are not. The user is advised to carefully",
        "interpret the output messages before attempting to bypass them with",
        "this option."
    };
    std::vector<MoleculeInformation>     mi;
    std::unique_ptr<MoleculeInformation> intermolecular_interactions;
    int                                  nvsite;
    CombinationRule                      comb;
    real                                 fudgeQQ;
    double                               reppow;
    const char*                          mdparin;
    bool                                 bNeedVel, bGenVel;
    gmx_output_env_t*                    oenv;
    gmx_bool                             bVerbose = FALSE;

    t_filenm fnm[] = { { efMDP, nullptr, nullptr, ffREAD },
                       { efMDP, "-po", "mdout", ffWRITE },
                       { efSTX, "-c", nullptr, ffREAD },
                       { efSTX, "-r", "restraint", ffOPTRD },
                       { efSTX, "-rb", "restraint", ffOPTRD },
                       { efNDX, nullptr, nullptr, ffOPTRD },
                       { efTOP, nullptr, nullptr, ffREAD },
                       { efTOP, "-pp", "processed", ffOPTWR },
                       { efTPR, "-o", nullptr, ffWRITE },
                       { efTRN, "-t", nullptr, ffOPTRD },
                       { efEDR, "-e", nullptr, ffOPTRD },
                       /* This group is needed by the VMD viewer as the start configuration for IMD sessions: */
                       { efGRO, "-imd", "imdgroup", ffOPTWR },
                       { efTRN, "-ref", "rotref", ffOPTRW | ffALLOW_MISSING },
                       /* This group is needed by the QMMM MDModule: */
                       { efQMI, "-qmi", nullptr, ffOPTRD } };
#define NFILE asize(fnm)

    /* Command line options */
    gmx_bool bRenum   = TRUE;
    gmx_bool bRmVSBds = TRUE, bZero = FALSE;
    int      maxwarn = 0;
    real     fr_time = -1;
    t_pargs  pa[]    = {
        { "-v", FALSE, etBOOL, { &bVerbose }, "Be loud and noisy" },
        { "-time", FALSE, etREAL, { &fr_time }, "Take frame at or first after this time." },
        { "-rmvsbds",
          FALSE,
          etBOOL,
          { &bRmVSBds },
          "Remove constant bonded interactions with virtual sites" },
        { "-maxwarn",
          FALSE,
          etINT,
          { &maxwarn },
          "Number of allowed warnings during input processing. Not for normal use and may "
          "generate unstable systems" },
        { "-zero",
          FALSE,
          etBOOL,
          { &bZero },
          "Set parameters for bonded interactions without defaults to zero instead of "
          "generating an error" },
        { "-renum",
          FALSE,
          etBOOL,
          { &bRenum },
          "Renumber atomtypes and minimize number of atomtypes" }
    };

    /* Parse the command line */
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    /* Initiate some variables */
    gmx::MDModules mdModules;
    t_inputrec     irInstance;
    t_inputrec*    ir = &irInstance;
    t_gromppopts   optsInstance;
    t_gromppopts*  opts = &optsInstance;
    snew(opts->include, STRLEN);
    snew(opts->define, STRLEN);

    gmx::LoggerBuilder builder;
    builder.addTargetStream(gmx::MDLogger::LogLevel::Info, &gmx::TextOutputFile::standardOutput());
    builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &gmx::TextOutputFile::standardError());
    gmx::LoggerOwner    logOwner(builder.build());
    const gmx::MDLogger logger(logOwner.logger());


    WarningHandler wi{ true, maxwarn };

    /* PARAMETER file processing */
    mdparin = opt2fn("-f", NFILE, fnm);
    wi.setFileAndLineNumber(mdparin, -1);
    try
    {
        get_ir(mdparin, opt2fn("-po", NFILE, fnm), &mdModules, ir, opts, WriteMdpHeader::yes, &wi);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    // Now that the MDModules have their options assigned from get_ir, subscribe
    // to eventual notifications during pre-processing their data
    mdModules.subscribeToPreProcessingNotifications();

    // Notify MDModules of existing logger
    mdModules.notifiers().preProcessingNotifier_.notify(logger);

    // Notify MDModules of existing warninp
    mdModules.notifiers().preProcessingNotifier_.notify(&wi);

    // Notify QMMM MDModule of external QM input file command-line option
    {
        gmx::QMInputFileName qmInputFileName = { ftp2bSet(efQMI, NFILE, fnm), ftp2fn(efQMI, NFILE, fnm) };
        mdModules.notifiers().preProcessingNotifier_.notify(qmInputFileName);
    }

    if (bVerbose)
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("checking input for internal consistency...");
    }
    check_ir(mdparin, mdModules.notifiers(), ir, opts, &wi);

    if (ir->ld_seed == -1)
    {
        ir->ld_seed = static_cast<int>(gmx::makeRandomSeed());
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Setting the LD random seed to %" PRId64 "", ir->ld_seed);
    }

    if (ir->expandedvals->lmc_seed == -1)
    {
        ir->expandedvals->lmc_seed = static_cast<int>(gmx::makeRandomSeed());
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("Setting the lambda MC random seed to %d", ir->expandedvals->lmc_seed);
    }

    bNeedVel = EI_STATE_VELOCITY(ir->eI);
    bGenVel  = (bNeedVel && opts->bGenVel);
    if (bGenVel && ir->bContinuation)
    {
        std::string warningMessage = gmx::formatString(
                "Generating velocities is inconsistent with attempting "
                "to continue a previous run. Choose only one of "
                "gen-vel = yes and continuation = yes.");
        wi.addError(warningMessage);
    }

    std::array<InteractionsOfType, F_NRE> interactions;
    gmx_mtop_t                            sys;
    PreprocessingAtomTypes                atypes;
    if (debug)
    {
        pr_symtab(debug, 0, "Just opened", &sys.symtab);
    }

    const char* fn = ftp2fn(efTOP, NFILE, fnm);
    if (!gmx_fexist(fn))
    {
        gmx_fatal(FARGS, "%s does not exist", fn);
    }

    t_state state;
    new_status(fn,
               opt2path_optional("-pp", NFILE, fnm),
               opt2fn("-c", NFILE, fnm),
               opts,
               ir,
               bZero,
               bGenVel,
               bVerbose,
               &state,
               &atypes,
               &sys,
               &mi,
               &intermolecular_interactions,
               interactions,
               &comb,
               &reppow,
               &fudgeQQ,
               opts->bMorse,
               &wi,
               logger);

    if (debug)
    {
        pr_symtab(debug, 0, "After new_status", &sys.symtab);
    }

    nvsite = 0;
    /* set parameters for virtual site construction (not for vsiten) */
    for (size_t mt = 0; mt < sys.moltype.size(); mt++)
    {
        nvsite += set_vsites(bVerbose, &sys.moltype[mt].atoms, &atypes, mi[mt].interactions, logger);
    }
    /* now throw away all obsolete bonds, angles and dihedrals: */
    /* note: constraints are ALWAYS removed */
    if (nvsite)
    {
        for (size_t mt = 0; mt < sys.moltype.size(); mt++)
        {
            clean_vsite_bondeds(mi[mt].interactions, sys.moltype[mt].atoms.nr, bRmVSBds, logger);
        }
    }

    if ((count_constraints(&sys, mi, &wi) != 0) && (ir->eConstrAlg == ConstraintAlgorithm::Shake))
    {
        if (ir->eI == IntegrationAlgorithm::CG || ir->eI == IntegrationAlgorithm::LBFGS)
        {
            std::string warningMessage =
                    gmx::formatString("Can not do %s with %s, use %s",
                                      enumValueToString(ir->eI),
                                      enumValueToString(ConstraintAlgorithm::Shake),
                                      enumValueToString(ConstraintAlgorithm::Lincs));
            wi.addError(warningMessage);
        }
        if (ir->bPeriodicMols)
        {
            std::string warningMessage =
                    gmx::formatString("Can not do periodic molecules with %s, use %s",
                                      enumValueToString(ConstraintAlgorithm::Shake),
                                      enumValueToString(ConstraintAlgorithm::Lincs));
            wi.addError(warningMessage);
        }
    }

    if (EI_SD(ir->eI) && ir->etc != TemperatureCoupling::No)
    {
        wi.addNote("Temperature coupling is ignored with SD integrators.");
    }

    /* Check for errors in the input now, since they might cause problems
     * during processing further down.
     */
    check_warning_error(wi, FARGS);

    if (nint_ftype(&sys, mi, F_POSRES) > 0 || nint_ftype(&sys, mi, F_FBPOSRES) > 0)
    {
        if (ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
            || ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            std::string warningMessage = gmx::formatString(
                    "You are combining position restraints with %s pressure coupling, which can "
                    "lead to instabilities. If you really want to combine position restraints with "
                    "pressure coupling, we suggest to use %s pressure coupling instead.",
                    enumValueToString(ir->pressureCouplingOptions.epc),
                    enumValueToString(PressureCoupling::CRescale));
            wi.addNote(warningMessage);
        }

        const char* fn = opt2fn("-r", NFILE, fnm);
        const char* fnB;

        if (!gmx_fexist(fn))
        {
            gmx_fatal(FARGS,
                      "Cannot find position restraint file %s (option -r).\n"
                      "From GROMACS-2018, you need to specify the position restraint "
                      "coordinate files explicitly to avoid mistakes, although you can "
                      "still use the same file as you specify for the -c option.",
                      fn);
        }

        if (opt2bSet("-rb", NFILE, fnm))
        {
            fnB = opt2fn("-rb", NFILE, fnm);
            if (!gmx_fexist(fnB))
            {
                gmx_fatal(FARGS,
                          "Cannot find B-state position restraint file %s (option -rb).\n"
                          "From GROMACS-2018, you need to specify the position restraint "
                          "coordinate files explicitly to avoid mistakes, although you can "
                          "still use the same file as you specify for the -c option.",
                          fnB);
            }
        }
        else
        {
            fnB = fn;
        }

        if (bVerbose)
        {
            std::string message = gmx::formatString("Reading position restraint coords from %s", fn);
            if (strcmp(fn, fnB) != 0)
            {
                message += gmx::formatString(" and %s", fnB);
            }
            GMX_LOG(logger.info).asParagraph().appendText(message);
        }
        gen_posres(&sys,
                   mi,
                   fn,
                   fnB,
                   ir->pressureCouplingOptions.refcoord_scaling,
                   ir->pbcType,
                   ir->posres_com,
                   ir->posres_comB,
                   &wi,
                   logger);
    }

    /* If we are using CMAP, setup the pre-interpolation grid */
    if (interactions[F_CMAP].ncmap() > 0)
    {
        init_cmap_grid(&sys.ffparams.cmap_grid,
                       interactions[F_CMAP].numCmaps_,
                       interactions[F_CMAP].cmapGridSpacing_);
        setup_cmap(interactions[F_CMAP].cmapGridSpacing_,
                   interactions[F_CMAP].numCmaps_,
                   interactions[F_CMAP].cmap,
                   &sys.ffparams.cmap_grid);
    }

    set_wall_atomtype(&atypes, opts, ir, &wi, logger);
    if (bRenum)
    {
        atypes.renumberTypes(interactions, &sys, ir->wall_atomtype, bVerbose);
    }

    if (ir->implicit_solvent)
    {
        gmx_fatal(FARGS, "Implicit solvation is no longer supported");
    }

    if (debug)
    {
        pr_symtab(debug, 0, "After atype.renumberTypes", &sys.symtab);
    }

    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("converting bonded parameters...");
    }

    const int ntype = atypes.size();
    convertInteractionsOfType(
            ntype, interactions, mi, intermolecular_interactions.get(), comb, reppow, fudgeQQ, &sys);

    if (debug)
    {
        pr_symtab(debug, 0, "After converInteractionsOfType", &sys.symtab);
    }

    /* set ptype to VSite for virtual sites */
    for (gmx_moltype_t& moltype : sys.moltype)
    {
        set_vsites_ptype(FALSE, &moltype, logger);
    }
    if (debug)
    {
        pr_symtab(debug, 0, "After virtual sites", &sys.symtab);
    }
    /* Check velocity for virtual sites and shells */
    if (bGenVel)
    {
        check_vel(&sys, state.v.rvec_array());
    }

    /* check for shells and inpurecs */
    check_shells_inputrec(&sys, ir, &wi);

    /* check masses */
    check_mol(&sys, &wi);

    checkRBDihedralSum(sys, *ir, &wi);

    if (haveFepPerturbedMassesInSettles(sys))
    {
        wi.addError(
                "SETTLE is not implemented for atoms whose mass is perturbed. "
                "You might instead use normal constraints.");
    }

    checkForUnboundAtoms(&sys, bVerbose, &wi, logger);

    // Now that we have the topology finalized and checked, we can repartition masses
    if (ir->massRepartitionFactor > 1)
    {
        const bool useFep = (ir->efep != FreeEnergyPerturbationType::No);

        gmx::repartitionAtomMasses(&sys, useFep, ir->massRepartitionFactor, &wi);
    }
    else if (ir->massRepartitionFactor < 1)
    {
        wi.addError("The mass repartitioning factor should be >= 1");
    }

    if (EI_DYNAMICS(ir->eI) && ir->eI != IntegrationAlgorithm::BD)
    {
        check_bonds_timestep(&sys, ir->delta_t, &wi);
    }

    checkDecoupledModeAccuracy(&sys, ir, &wi);

    if (EI_ENERGY_MINIMIZATION(ir->eI) && 0 == ir->nsteps)
    {
        wi.addNote(
                "Zero-step energy minimization will alter the coordinates before calculating the "
                "energy. If you just want the energy of a single point, try zero-step MD (with "
                "unconstrained_start = yes). To do multiple single-point energy evaluations of "
                "different configurations of the same topology, use mdrun -rerun.");
    }

    check_warning_error(wi, FARGS);

    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("initialising group options...");
    }
    do_index(mdparin, ftp2path_optional(efNDX, NFILE, fnm), &sys, bVerbose, mdModules.notifiers(), ir, &wi);

    // Notify topology to MdModules for pre-processing after all indexes were built
    mdModules.notifiers().preProcessingNotifier_.notify(&sys);

    if (usingFullElectrostatics(ir->coulombtype) || usingLJPme(ir->vdwtype))
    {
        // We may have exclusion forces beyond the cut-off distance.
        checkExclusionDistances(sys, *ir, state.x, state.box, logger, &wi);
    }

    if (ir->cutoff_scheme == CutoffScheme::Verlet && ir->verletbuf_tol > 0)
    {
        if (EI_DYNAMICS(ir->eI) && inputrec2nboundeddim(ir) == 3)
        {
            real buffer_temp;

            if (EI_MD(ir->eI) && ir->etc == TemperatureCoupling::No)
            {
                if (bGenVel)
                {
                    buffer_temp = opts->tempi;
                }
                else
                {
                    buffer_temp = calc_temp(&sys, ir, state.v.rvec_array());
                }
                if (buffer_temp > 0)
                {
                    std::string warningMessage = gmx::formatString(
                            "NVE simulation: will use the initial temperature of %.3f K for "
                            "determining the Verlet buffer size",
                            buffer_temp);
                    wi.addNote(warningMessage);
                }
                else
                {
                    std::string warningMessage = gmx::formatString(
                            "NVE simulation with an initial temperature of zero: will use a Verlet "
                            "buffer of %d%%. Check your energy drift!",
                            gmx::roundToInt(verlet_buffer_ratio_NVE_T0 * 100));
                    wi.addNote(warningMessage);
                }
            }
            else
            {
                buffer_temp = get_max_reference_temp(ir, &wi);
            }

            if (EI_MD(ir->eI) && ir->etc == TemperatureCoupling::No && buffer_temp == 0)
            {
                /* NVE with initial T=0: we add a fixed ratio to rlist.
                 * Since we don't actually use verletbuf_tol, we set it to -1
                 * so it can't be misused later.
                 */
                ir->rlist *= 1.0 + verlet_buffer_ratio_NVE_T0;
                ir->verletbuf_tol = -1;
            }
            else
            {
                /* We warn for NVE simulations with a drift tolerance that
                 * might result in a 1(.1)% drift over the total run-time.
                 * Note that we can't warn when nsteps=0, since we don't
                 * know how many steps the user intends to run.
                 */
                if (EI_MD(ir->eI) && ir->etc == TemperatureCoupling::No && ir->nstlist > 1 && ir->nsteps > 0)
                {
                    const real driftTolerance = 0.01;
                    /* We use 2 DOF per atom = 2kT pot+kin energy,
                     * to be on the safe side with constraints.
                     */
                    const real totalEnergyDriftPerAtomPerPicosecond =
                            2 * gmx::c_boltz * buffer_temp / (ir->nsteps * ir->delta_t);

                    if (ir->verletbuf_tol > 1.1 * driftTolerance * totalEnergyDriftPerAtomPerPicosecond)
                    {
                        std::string warningMessage = gmx::formatString(
                                "You are using a Verlet buffer tolerance of %g kJ/mol/ps for an "
                                "NVE simulation of length %g ps, which can give a final drift of "
                                "%d%%. For conserving energy to %d%% when using constraints, you "
                                "might need to set verlet-buffer-tolerance to %.1e.",
                                ir->verletbuf_tol,
                                ir->nsteps * ir->delta_t,
                                gmx::roundToInt(ir->verletbuf_tol / totalEnergyDriftPerAtomPerPicosecond * 100),
                                gmx::roundToInt(100 * driftTolerance),
                                driftTolerance * totalEnergyDriftPerAtomPerPicosecond);
                        wi.addNote(warningMessage);
                    }
                }

                set_verlet_buffer(&sys, ir, buffer_temp, state.x, state.box, &wi, logger);
            }
        }
    }

    /* Init the temperature coupling state */
    init_gtc_state(&state, ir->opts.ngtc, 0, ir->opts.nhchainlength); /* need to add nnhpres here? */

    /* After we are done with all checks on the state, we can add the flow profile */
    if (opts->deformInitFlow)
    {
        deformInitFlow(&state, ir->deform);
    }

    if (debug)
    {
        pr_symtab(debug, 0, "After index", &sys.symtab);
    }

    triple_check(mdparin, ir, &sys, &wi);
    close_symtab(&sys.symtab);
    if (debug)
    {
        pr_symtab(debug, 0, "After close", &sys.symtab);
    }

    if (ir->eI == IntegrationAlgorithm::Mimic)
    {
        generate_qmexcl(&sys, ir, logger);
    }

    if (ftp2bSet(efTRN, NFILE, fnm))
    {
        if (bVerbose)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted("getting data from old trajectory ...");
        }
        cont_status(ftp2fn(efTRN, NFILE, fnm),
                    ftp2path_optional(efEDR, NFILE, fnm),
                    bNeedVel,
                    bGenVel,
                    fr_time,
                    ir,
                    &state,
                    &sys,
                    oenv,
                    logger);
    }

    if (ir->pbcType == PbcType::XY && ir->nwall != 2)
    {
        clear_rvec(state.box[ZZ]);
    }

    if (usingFullElectrostatics(ir->coulombtype) || usingLJPme(ir->vdwtype))
    {
        /* Calculate the optimal grid dimensions */
        matrix          scaledBox;
        EwaldBoxZScaler boxScaler(inputrecPbcXY2Walls(ir), ir->wall_ewald_zfac);
        boxScaler.scaleBox(state.box, scaledBox);

        if (ir->nkx > 0 && ir->nky > 0 && ir->nkz > 0)
        {
            /* Mark fourier_spacing as not used */
            ir->fourier_spacing = 0;
        }
        else if (ir->nkx != 0 && ir->nky != 0 && ir->nkz != 0)
        {
            wi.setFileAndLineNumber(mdparin, -1);
            wi.addError("Some of the Fourier grid sizes are set, but all of them need to be set.");
        }
        const int minGridSize = minimalPmeGridSize(ir->pme_order);
        calcFftGrid(stdout, scaledBox, ir->fourier_spacing, minGridSize, &(ir->nkx), &(ir->nky), &(ir->nkz));
        if (ir->nkx < minGridSize || ir->nky < minGridSize || ir->nkz < minGridSize)
        {
            wi.addError(
                    "The PME grid size should be >= 2*(pme-order - 1); either manually "
                    "increase the grid size or decrease pme-order");
        }
    }

    /* MRS: eventually figure out better logic for initializing the fep
       values that makes declaring the lambda and declaring the state not
       potentially conflict if not handled correctly. */
    if (ir->efep != FreeEnergyPerturbationType::No)
    {
        state.fep_state = ir->fepvals->init_fep_state;
        for (const auto couplingType : gmx::EnumerationWrapper<FreeEnergyPerturbationCouplingType>{})
        {
            state.lambda[static_cast<int>(couplingType)] = ir->fepvals->initialLambda(couplingType);
        }
    }

    pull_t* pull = nullptr;

    if (ir->bPull)
    {
        pull = set_pull_init(
                ir, sys, state.x, state.box, state.lambda[FreeEnergyPerturbationCouplingType::Mass], &wi);
    }

    /* Modules that supply external potential for pull coordinates
     * should register those potentials here. finish_pull() will check
     * that providers have been registerd for all external potentials.
     */

    if (ir->bDoAwh)
    {
        tensor compressibility = { { 0 } };
        if (ir->pressureCouplingOptions.epc != PressureCoupling::No)
        {
            copy_mat(ir->pressureCouplingOptions.compress, compressibility);
        }
        real initialLambda = 0;
        if (ir->efep != FreeEnergyPerturbationType::No)
        {
            initialLambda = ir->fepvals->initialLambda(FreeEnergyPerturbationCouplingType::Fep);
        }
        setStateDependentAwhParams(
                ir->awhParams.get(), *ir->pull, pull, state.box, ir->pbcType, compressibility, *ir, initialLambda, sys, &wi);
    }

    if (ir->bPull)
    {
        finish_pull(pull);
    }

    if (ir->bRot)
    {
        set_reference_positions(ir->rot.get(),
                                state.x.rvec_array(),
                                state.box,
                                opt2fn("-ref", NFILE, fnm),
                                opt2bSet("-ref", NFILE, fnm),
                                &wi);
    }

    /*  reset_multinr(sys); */

    if (usingPme(ir->coulombtype))
    {
        float ratio = pme_load_estimate(sys, *ir, state.box);
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "Estimate for the relative computational load of the PME mesh part: %.2f", ratio);
        /* With free energy we might need to do PME both for the A and B state
         * charges. This will double the cost, but the optimal performance will
         * then probably be at a slightly larger cut-off and grid spacing.
         */
        if ((ir->efep == FreeEnergyPerturbationType::No && ratio > 1.0 / 2.0)
            || (ir->efep != FreeEnergyPerturbationType::No && ratio > 2.0 / 3.0))
        {
            wi.addNote(
                    "The optimal PME mesh load for parallel simulations is below 0.5\n"
                    "and for highly parallel simulations between 0.25 and 0.33,\n"
                    "for higher performance, increase the cut-off and the PME grid spacing.\n");
            if (ir->efep != FreeEnergyPerturbationType::No)
            {
                wi.addNote(
                        "For free energy simulations, the optimal load limit increases from "
                        "0.5 to 0.667\n");
            }
        }
    }

    {
        double      cio = compute_io(ir, sys.natoms, sys.groups, F_NRE, 1);
        std::string warningMessage =
                gmx::formatString("This run will generate roughly %.0f Mb of data", cio);
        if (cio > 2000)
        {
            wi.setFileAndLineNumber(mdparin, -1);
            wi.addNote(warningMessage);
        }
        else
        {
            GMX_LOG(logger.info).asParagraph().appendText(warningMessage);
        }
    }

    // Hand over box and coordiantes to MdModules before they evaluate their final parameters
    {
        gmx::CoordinatesAndBoxPreprocessed coordinatesAndBoxPreprocessed;
        coordinatesAndBoxPreprocessed.coordinates_ = state.x.arrayRefWithPadding();
        copy_mat(state.box, coordinatesAndBoxPreprocessed.box_);
        coordinatesAndBoxPreprocessed.pbc_ = ir->pbcType;
        mdModules.notifiers().preProcessingNotifier_.notify(coordinatesAndBoxPreprocessed);

        // Send also the constant ensemble temperature if available.
        gmx::EnsembleTemperature ensembleTemperature(*ir);
        mdModules.notifiers().preProcessingNotifier_.notify(ensembleTemperature);
    }

    // Add the md modules internal parameters that are not mdp options
    // e.g., atom indices

    {
        gmx::KeyValueTreeBuilder internalParameterBuilder;
        mdModules.notifiers().preProcessingNotifier_.notify(internalParameterBuilder.rootObject());
        ir->internalParameters =
                std::make_unique<gmx::KeyValueTreeObject>(internalParameterBuilder.build());
    }

    if (ir->comm_mode != ComRemovalAlgorithm::No)
    {
        const int nstglobalcomm = computeGlobalCommunicationPeriod(ir);
        if (ir->nstcomm % nstglobalcomm != 0)
        {
            wi.addNote(

                    gmx::formatString(
                            "COM removal frequency is set to (%d).\n"
                            "Other settings require a global communication frequency of %d.\n"
                            "Note that this will require additional global communication steps,\n"
                            "which will reduce performance when using multiple ranks.\n"
                            "Consider setting nstcomm to a multiple of %d.",
                            ir->nstcomm,
                            nstglobalcomm,
                            nstglobalcomm));
        }
    }

    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("writing run input file...");
    }

    done_warning(wi, FARGS);
    write_tpx_state(ftp2fn(efTPR, NFILE, fnm), ir, &state, sys);

    /* Output IMD group, if bIMD is TRUE */
    gmx::write_IMDgroup_to_file(ir->bIMD, ir, &state, sys, NFILE, fnm);

    sfree(opts->define);
    sfree(opts->wall_atomtype[0]);
    sfree(opts->wall_atomtype[1]);
    sfree(opts->include);
    sfree(opts->couple_moltype);

    for (auto& mol : mi)
    {
        // Some of the contents of molinfo have been stolen, so
        // fullCleanUp can't be called.
        mol.partialCleanUp();
    }
    done_inputrec_strings();
    output_env_done(oenv);

    return 0;
}
