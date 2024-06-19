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

#include "insert_molecules.h"

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/conformation_utilities.h"
#include "gromacs/gmxpreprocess/makeexclusiondistances.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/atomsbuilder.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
class CommandLineModuleSettings;
} // namespace gmx

using gmx::RVec;

/* enum for random rotations of inserted solutes */
enum class RotationType : int
{
    XYZ,
    Z,
    None,
    Count
};
static const gmx::EnumerationArray<RotationType, const char*> c_rotationTypeNames = {
    { "xyz", "z", "none" }
};

static void center_molecule(gmx::ArrayRef<RVec> x)
{
    rvec center = { 0 };
    for (auto& xi : x)
    {
        rvec_inc(center, xi);
    }
    svmul(1.0 / x.size(), center, center);
    for (auto& xi : x)
    {
        rvec_dec(xi, center);
    }
}

static void generate_trial_conf(gmx::ArrayRef<RVec>       xin,
                                const rvec                offset,
                                RotationType              enum_rot,
                                gmx::DefaultRandomEngine* rng,
                                std::vector<RVec>*        xout)
{
    gmx::UniformRealDistribution<real> dist(0, 2.0 * M_PI);
    xout->assign(xin.begin(), xin.end());

    real alfa = 0.0, beta = 0.0, gamma = 0.0;
    switch (enum_rot)
    {
        case RotationType::XYZ:
            alfa  = dist(*rng);
            beta  = dist(*rng);
            gamma = dist(*rng);
            break;
        case RotationType::Z:
            alfa = beta = 0.;
            gamma       = dist(*rng);
            break;
        case RotationType::None: alfa = beta = gamma = 0.; break;
        default: GMX_THROW(gmx::InternalError("Invalid RotationType"));
    }
    if (enum_rot == RotationType::XYZ || enum_rot == RotationType::Z)
    {
        rotate_conf(xout->size(), as_rvec_array(xout->data()), nullptr, alfa, beta, gamma);
    }
    for (size_t i = 0; i < xout->size(); ++i)
    {
        rvec_inc((*xout)[i], offset);
    }
}

static bool isInsertionAllowed(gmx::AnalysisNeighborhoodSearch* search,
                               const std::vector<real>&         exclusionDistances,
                               const std::vector<RVec>&         x,
                               const std::vector<real>&         exclusionDistances_insrt,
                               const t_atoms&                   atoms,
                               const std::set<int>&             removableAtoms,
                               gmx::AtomsRemover*               remover)
{
    gmx::AnalysisNeighborhoodPositions  pos(x);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = search->startPairSearch(pos);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        const real r1 = exclusionDistances[pair.refIndex()];
        const real r2 = exclusionDistances_insrt[pair.testIndex()];
        if (pair.distance2() < gmx::square(r1 + r2))
        {
            if (removableAtoms.count(pair.refIndex()) == 0)
            {
                return false;
            }
            // TODO: If molecule information is available, this should ideally
            // use it to remove whole molecules.
            remover->markResidue(atoms, pair.refIndex(), true);
        }
    }
    return true;
}

static void insert_mols(int                  nmol_insrt,
                        int                  ntry,
                        int                  seed,
                        real                 defaultDistance,
                        real                 scaleFactor,
                        t_atoms*             atoms,
                        t_symtab*            symtab,
                        std::vector<RVec>*   x,
                        const std::set<int>& removableAtoms,
                        const t_atoms&       atoms_insrt,
                        gmx::ArrayRef<RVec>  x_insrt,
                        PbcType              pbcType,
                        matrix               box,
                        const std::string&   posfn,
                        const rvec           deltaR,
                        RotationType         enum_rot)
{
    fprintf(stderr, "Initialising inter-atomic distances...\n");
    AtomProperties aps;
    std::vector<real> exclusionDistances(makeExclusionDistances(atoms, &aps, defaultDistance, scaleFactor));
    const std::vector<real> exclusionDistances_insrt(
            makeExclusionDistances(&atoms_insrt, &aps, defaultDistance, scaleFactor));

    const real maxInsertRadius =
            *std::max_element(exclusionDistances_insrt.begin(), exclusionDistances_insrt.end());
    real maxRadius = maxInsertRadius;
    if (!exclusionDistances.empty())
    {
        const real maxExistingRadius =
                *std::max_element(exclusionDistances.begin(), exclusionDistances.end());
        maxRadius = std::max(maxInsertRadius, maxExistingRadius);
    }

    // TODO: Make all of this exception-safe.
    gmx::AnalysisNeighborhood nb;
    nb.setCutoff(maxInsertRadius + maxRadius);


    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
    }
    fprintf(stderr, "Using random seed %d\n", seed);

    gmx::DefaultRandomEngine rng(seed);

    t_pbc pbc;
    set_pbc(&pbc, pbcType, box);

    /* With -ip, take nmol_insrt from file posfn */
    double**   rpos              = nullptr;
    const bool insertAtPositions = !posfn.empty();
    if (insertAtPositions)
    {
        int ncol;
        nmol_insrt = read_xvg(posfn.c_str(), &rpos, &ncol);
        if (ncol != 3)
        {
            gmx_fatal(FARGS, "Expected 3 columns (x/y/z coordinates) in file %s\n", posfn.c_str());
        }
        fprintf(stderr, "Read %d positions from file %s\n\n", nmol_insrt, posfn.c_str());
    }

    gmx::AtomsBuilder builder(atoms, symtab);
    gmx::AtomsRemover remover(*atoms);
    {
        const int finalAtomCount    = atoms->nr + nmol_insrt * atoms_insrt.nr;
        const int finalResidueCount = atoms->nres + nmol_insrt * atoms_insrt.nres;
        builder.reserve(finalAtomCount, finalResidueCount);
        x->reserve(finalAtomCount);
        exclusionDistances.reserve(finalAtomCount);
    }

    std::vector<RVec> x_n(x_insrt.size());

    int                                mol        = 0;
    int                                trial      = 0;
    int                                firstTrial = 0;
    int                                failed     = 0;
    gmx::UniformRealDistribution<real> dist;

    while (mol < nmol_insrt && trial < ntry * nmol_insrt)
    {
        rvec offset_x;
        if (!insertAtPositions)
        {
            // Insert at random positions.
            offset_x[XX] = box[XX][XX] * dist(rng);
            offset_x[YY] = box[YY][YY] * dist(rng);
            offset_x[ZZ] = box[ZZ][ZZ] * dist(rng);
        }
        else
        {
            // Skip a position if ntry trials were not successful.
            if (trial >= firstTrial + ntry)
            {
                fprintf(stderr,
                        " skipped position (%.3f, %.3f, %.3f)\n",
                        rpos[XX][mol],
                        rpos[YY][mol],
                        rpos[ZZ][mol]);
                ++mol;
                ++failed;
                firstTrial = trial;
                continue;
            }
            // Insert at positions taken from option -ip file.
            offset_x[XX] = rpos[XX][mol] + deltaR[XX] * (2 * dist(rng) - 1);
            offset_x[YY] = rpos[YY][mol] + deltaR[YY] * (2 * dist(rng) - 1);
            offset_x[ZZ] = rpos[ZZ][mol] + deltaR[ZZ] * (2 * dist(rng) - 1);
        }
        fprintf(stderr, "\rTry %d", ++trial);
        fflush(stderr);

        generate_trial_conf(x_insrt, offset_x, enum_rot, &rng, &x_n);
        gmx::AnalysisNeighborhoodPositions pos(*x);
        gmx::AnalysisNeighborhoodSearch    search = nb.initSearch(&pbc, pos);
        if (isInsertionAllowed(
                    &search, exclusionDistances, x_n, exclusionDistances_insrt, *atoms, removableAtoms, &remover))
        {
            x->insert(x->end(), x_n.begin(), x_n.end());
            exclusionDistances.insert(exclusionDistances.end(),
                                      exclusionDistances_insrt.begin(),
                                      exclusionDistances_insrt.end());
            builder.mergeAtoms(atoms_insrt);
            ++mol;
            firstTrial = trial;
            fprintf(stderr, " success (now %d atoms)!\n", builder.currentAtomCount());
        }
    }

    fprintf(stderr, "\n");
    /* print number of molecules added */
    fprintf(stderr, "Added %d molecules (out of %d requested)\n", mol - failed, nmol_insrt);

    const int originalAtomCount    = atoms->nr;
    const int originalResidueCount = atoms->nres;
    remover.refreshAtomCount(*atoms);
    remover.removeMarkedElements(x);
    remover.removeMarkedAtoms(atoms);
    if (atoms->nr < originalAtomCount)
    {
        fprintf(stderr,
                "Replaced %d residues (%d atoms)\n",
                originalResidueCount - atoms->nres,
                originalAtomCount - atoms->nr);
    }

    if (rpos != nullptr)
    {
        for (int i = 0; i < DIM; ++i)
        {
            sfree(rpos[i]);
        }
        sfree(rpos);
    }
}

namespace gmx
{

namespace
{

class InsertMolecules : public ICommandLineOptionsModule, public ITopologyProvider
{
public:
    InsertMolecules() :
        bBox_(false),
        nmolIns_(0),
        concIns_(0.0),
        concInsIsSet_(false),
        nmolTry_(10),
        seed_(0),
        defaultDistance_(0.105),
        scaleFactor_(0.57),
        enumRot_(RotationType::XYZ)
    {
        clear_rvec(newBox_);
        clear_rvec(deltaR_);
        clear_mat(box_);
    }

    // From ITopologyProvider
    gmx_mtop_t* getTopology(bool /*required*/) override { return &top_; }
    int         getAtomCount() override { return 0; }

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

private:
    SelectionCollection selections_;

    std::string  inputConfFile_;
    std::string  insertConfFile_;
    std::string  positionFile_;
    std::string  outputConfFile_;
    rvec         newBox_;
    bool         bBox_;
    int          nmolIns_;
    real         concIns_;
    bool         concInsIsSet_;
    int          nmolTry_;
    int          seed_;
    real         defaultDistance_;
    real         scaleFactor_;
    rvec         deltaR_;
    RotationType enumRot_;
    Selection    replaceSel_;

    gmx_mtop_t        top_;
    std::vector<RVec> x_;
    matrix            box_;
    PbcType           pbcType_;
};

void InsertMolecules::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    const char* const desc[] = {
        "[THISMODULE] inserts [TT]-nmol[tt] copies of the system specified in",
        "the [TT]-ci[tt] input file. The number of copies can also be determined",
        "by the concentration [TT]-conc[tt] in mol/liter and box volume.",
        "The insertions take place either into vacant space in the solute conformation",
        "given with [TT]-f[tt], or into an empty box given by [TT]-box[tt].",
        "Specifying both [TT]-f[tt] and [TT]-box[tt] behaves like [TT]-f[tt],",
        "but places a new box around the solute before insertions.",
        "Any velocities present are discarded.",
        "",
        "It is possible to also insert into a solvated configuration and",
        "replace solvent atoms with the inserted atoms. To do this, use",
        "[TT]-replace[tt] to specify a selection that identifies the atoms",
        "that can be replaced. The tool assumes that all molecules in this",
        "selection consist of single residues: each residue from this",
        "selection that overlaps with the inserted molecules will be removed",
        "instead of preventing insertion.",
        "",
        "By default, the insertion positions are random (with initial seed",
        "specified by [TT]-seed[tt]). The program iterates until [TT]-nmol[tt]",
        "molecules have been inserted in the box. Molecules are not inserted",
        "where the distance between any existing atom and any atom of the",
        "inserted molecule is less than the sum based on the van der Waals",
        "radii of both atoms. A database ([TT]vdwradii.dat[tt]) of van der",
        "Waals radii is read by the program, and the resulting radii scaled",
        "by [TT]-scale[tt]. If radii are not found in the database, those",
        "atoms are assigned the (pre-scaled) distance [TT]-radius[tt].",
        "Note that the usefulness of those radii depends on the atom names,",
        "and thus varies widely with force field.",
        "",
        "A total of [TT]-nmol[tt] * [TT]-try[tt] insertion attempts are made",
        "before giving up. Increase [TT]-try[tt] if you have several small",
        "holes to fill. Option [TT]-rot[tt] specifies whether the insertion",
        "molecules are randomly oriented before insertion attempts.",
        "",
        "Alternatively, the molecules can be inserted only at positions defined in",
        "positions.dat ([TT]-ip[tt]). That file should have 3 columns (x,y,z),",
        "that give the displacements compared to the input molecule position",
        "([TT]-ci[tt]). Hence, if that file should contain the absolute",
        "positions, the molecule must be centered on (0,0,0) before using",
        "[THISMODULE] (e.g. from [gmx-editconf] [TT]-center[tt]).",
        "Comments in that file starting with # are ignored. Option [TT]-dr[tt]",
        "defines the maximally allowed displacements during insertial trials.",
        "[TT]-try[tt] and [TT]-rot[tt] work as in the default mode (see above)."
    };

    settings->setHelpText(desc);

    std::shared_ptr<SelectionOptionBehavior> selectionOptionBehavior(
            new SelectionOptionBehavior(&selections_, this));
    settings->addOptionsBehavior(selectionOptionBehavior);

    // TODO: Replace use of legacyType.
    options->addOption(FileNameOption("f")
                               .legacyType(efSTX)
                               .inputFile()
                               .store(&inputConfFile_)
                               .defaultBasename("protein")
                               .description("Existing configuration to insert into"));
    options->addOption(FileNameOption("ci")
                               .legacyType(efSTX)
                               .inputFile()
                               .required()
                               .store(&insertConfFile_)
                               .defaultBasename("insert")
                               .description("Configuration to insert"));
    options->addOption(FileNameOption("ip")
                               .filetype(OptionFileType::GenericData)
                               .inputFile()
                               .store(&positionFile_)
                               .defaultBasename("positions")
                               .description("Predefined insertion trial positions"));
    options->addOption(FileNameOption("o")
                               .legacyType(efSTO)
                               .outputFile()
                               .required()
                               .store(&outputConfFile_)
                               .defaultBasename("out")
                               .description("Output configuration after insertion"));

    options->addOption(
            SelectionOption("replace").onlyAtoms().store(&replaceSel_).description("Atoms that can be removed if overlapping"));
    selectionOptionBehavior->initOptions(options);

    options->addOption(RealOption("box").vector().store(newBox_).storeIsSet(&bBox_).description(
            "Box size (in nm)"));
    options->addOption(IntegerOption("nmol").store(&nmolIns_).description(
            "Number of extra molecules to insert"));
    options->addOption(
            RealOption("conc")
                    .store(&concIns_)
                    .storeIsSet(&concInsIsSet_)
                    .description("Concentration (in mol/liter) of extra molecules to insert. "
                                 "This overrides [TT]-nmol[tt]"));
    options->addOption(IntegerOption("try").store(&nmolTry_).description(
            "Try inserting [TT]-nmol[tt] times [TT]-try[tt] times"));
    options->addOption(IntegerOption("seed").store(&seed_).description(
            "Random generator seed (0 means generate)"));
    options->addOption(
            RealOption("radius").store(&defaultDistance_).description("Default van der Waals distance"));
    options->addOption(
            RealOption("scale").store(&scaleFactor_).description("Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The default value of 0.57 yields density close to 1000 g/l for proteins in water."));
    options->addOption(RealOption("dr").vector().store(deltaR_).description(
            "Allowed displacement in x/y/z from positions in [TT]-ip[tt] file"));
    options->addOption(EnumOption<RotationType>("rot")
                               .enumValue(c_rotationTypeNames)
                               .store(&enumRot_)
                               .description("Rotate inserted molecules randomly"));
}

void InsertMolecules::optionsFinished()
{
    if ((nmolIns_ <= 0 && concIns_ <= 0.0) && positionFile_.empty())
    {
        GMX_THROW(
                InconsistentInputError("Either -nmol or -conc must be larger than 0, "
                                       "or positions must be given with -ip."));
    }
    if (inputConfFile_.empty() && !bBox_)
    {
        GMX_THROW(
                InconsistentInputError("When no solute (-f) is specified, "
                                       "a box size (-box) must be specified."));
    }
    if (replaceSel_.isValid() && inputConfFile_.empty())
    {
        GMX_THROW(
                InconsistentInputError("Replacement (-replace) only makes sense "
                                       "together with an existing configuration (-f)."));
    }

    if (!inputConfFile_.empty())
    {
        bool  bTprFileWasRead;
        rvec* temporaryX = nullptr;
        fprintf(stderr, "Reading solute configuration\n");
        readConfAndTopology(
                inputConfFile_.c_str(), &bTprFileWasRead, &top_, &pbcType_, &temporaryX, nullptr, box_);
        x_.assign(temporaryX, temporaryX + top_.natoms);
        sfree(temporaryX);
        if (top_.natoms == 0)
        {
            fprintf(stderr, "Note: no atoms in %s\n", inputConfFile_.c_str());
        }
    }
}

int InsertMolecules::run()
{
    std::set<int> removableAtoms;
    if (replaceSel_.isValid())
    {
        t_pbc pbc;
        set_pbc(&pbc, pbcType_, box_);
        t_trxframe* fr;
        snew(fr, 1);
        fr->natoms = x_.size();
        fr->bX     = TRUE;
        fr->x      = as_rvec_array(x_.data());
        selections_.evaluate(fr, &pbc);
        sfree(fr);
        removableAtoms.insert(replaceSel_.atomIndices().begin(), replaceSel_.atomIndices().end());
        // TODO: It could be nice to check that removableAtoms contains full
        // residues, since we anyways remove whole residues instead of
        // individual atoms.
    }

    PbcType pbcTypeForOutput = pbcType_;
    if (bBox_)
    {
        pbcTypeForOutput = PbcType::Xyz;
        clear_mat(box_);
        box_[XX][XX] = newBox_[XX];
        box_[YY][YY] = newBox_[YY];
        box_[ZZ][ZZ] = newBox_[ZZ];
    }
    const real vol = det(box_);
    if (!(vol > 0))
    {
        gmx_fatal(FARGS,
                  "Undefined solute box.\nCreate one with gmx editconf "
                  "or give explicit -box command line option");
    }
    if (concInsIsSet_)
    {
        if (concIns_ > 0)
        {
            nmolIns_ = gmx::roundToInt(concIns_ * vol * gmx::c_avogadro / 1e24);
            if (!nmolIns_)
            {
                gmx_fatal(FARGS, "Very low concentration is set, nothing to insert");
            }
        }
        else
        {
            gmx_fatal(FARGS, "Concentration is defined, but it is equal to zero or negative");
        }
    }

    gmx_mtop_t        topInserted;
    t_atoms           atomsInserted;
    std::vector<RVec> xInserted;
    {
        bool    bTprFileWasRead;
        PbcType pbcType_dummy;
        matrix  box_dummy;
        rvec*   temporaryX;
        readConfAndTopology(
                insertConfFile_.c_str(), &bTprFileWasRead, &topInserted, &pbcType_dummy, &temporaryX, nullptr, box_dummy);
        xInserted.assign(temporaryX, temporaryX + topInserted.natoms);
        sfree(temporaryX);
        atomsInserted = gmx_mtop_global_atoms(topInserted);
        if (atomsInserted.nr == 0)
        {
            gmx_fatal(FARGS, "No molecule in %s, please check your input", insertConfFile_.c_str());
        }
        if (top_.name == nullptr)
        {
            top_.name = topInserted.name;
        }
        if (positionFile_.empty())
        {
            center_molecule(xInserted);
        }
    }

    t_atoms atoms = gmx_mtop_global_atoms(top_);

    /* add nmol_ins molecules of atoms_ins
       in random orientation at random place */
    insert_mols(nmolIns_,
                nmolTry_,
                seed_,
                defaultDistance_,
                scaleFactor_,
                &atoms,
                &top_.symtab,
                &x_,
                removableAtoms,
                atomsInserted,
                xInserted,
                pbcTypeForOutput,
                box_,
                positionFile_,
                deltaR_,
                enumRot_);

    /* write new configuration to file confout */
    fprintf(stderr, "Writing generated configuration to %s\n", outputConfFile_.c_str());
    write_sto_conf(
            outputConfFile_.c_str(), *top_.name, &atoms, as_rvec_array(x_.data()), nullptr, pbcTypeForOutput, box_);

    /* print size of generated configuration */
    fprintf(stderr, "\nOutput configuration contains %d atoms in %d residues\n", atoms.nr, atoms.nres);

    done_atom(&atoms);
    done_atom(&atomsInserted);

    return 0;
}

} // namespace


const char* InsertMoleculesInfo::name()
{
    static const char* name = "insert-molecules";
    return name;
}

const char* InsertMoleculesInfo::shortDescription()
{
    static const char* shortDescription = "Insert molecules into existing vacancies";
    return shortDescription;
}

ICommandLineOptionsModulePointer InsertMoleculesInfo::create()
{
    return ICommandLineOptionsModulePointer(new InsertMolecules());
}

} // namespace gmx
