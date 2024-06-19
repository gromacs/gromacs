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
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Sasa.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com> (C++ conversion)
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "sasa.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "surfacearea.h"

enum class PbcType : int;

namespace gmx
{
class AnalysisDataParallelOptions;
class SelectionCollection;

namespace analysismodules
{

namespace
{

//! \addtogroup module_trajectoryanalysis
//! \{

//! Tracks information on two nearest neighbors of a single surface dot.
struct t_conect
{
    //! Index of the second nearest neighbor dot.
    int aa;
    //! Index of the nearest neighbor dot.
    int ab;
    //! Squared distance to `aa`.
    real d2a;
    //! Squared distance to `ab`.
    real d2b;
};

/*! \brief
 * Updates nearest neighbor information for a surface dot.
 *
 * \param[in,out] c  Nearest neighbor information array to update.
 * \param[in]     i  Index in `c` to update.
 * \param[in]     j  Index of the other surface dot to add to the array.
 * \param[in]     d2 Squared distance between `i` and `j`.
 */
void add_rec(t_conect c[], int i, int j, real d2)
{
    if (c[i].aa == -1)
    { // NOLINT bugprone-branch-clone
        c[i].aa  = j;
        c[i].d2a = d2;
    }
    else if (c[i].ab == -1)
    { // NOLINT bugprone-branch-clone
        c[i].ab  = j;
        c[i].d2b = d2;
    }
    else if (d2 < c[i].d2a)
    {
        c[i].aa  = j;
        c[i].d2a = d2;
    }
    else if (d2 < c[i].d2b)
    {
        c[i].ab  = j;
        c[i].d2b = d2;
    }
    /* Swap them if necessary: a must be larger than b */
    if (c[i].d2a < c[i].d2b)
    {
        j        = c[i].ab;
        c[i].ab  = c[i].aa;
        c[i].aa  = j;
        d2       = c[i].d2b;
        c[i].d2b = c[i].d2a;
        c[i].d2a = d2;
    }
}

/*! \brief
 * Adds CONECT records for surface dots.
 *
 * \param[in] fn  PDB file to append the CONECT records to.
 * \param[in] n   Number of dots in `x`.
 * \param[in] x   Array of surface dot positions.
 *
 * Adds a CONECT record that connects each surface dot to its two nearest
 * neighbors.  The function is copied verbatim from the old gmx_sas.c
 * implementation.
 */
void do_conect(const char* fn, int n, rvec x[])
{
    t_conect* c = nullptr;

    fprintf(stderr, "Building CONECT records\n");
    snew(c, n);
    for (int i = 0; (i < n); i++)
    {
        c[i].aa = c[i].ab = -1;
    }

    for (int i = 0; (i < n); i++)
    {
        for (int j = i + 1; (j < n); j++)
        {
            rvec dx;
            rvec_sub(x[i], x[j], dx);
            const real d2 = iprod(dx, dx);
            add_rec(c, i, j, d2);
            add_rec(c, j, i, d2);
        }
    }
    FILE* fp = gmx_ffopen(fn, "a");
    for (int i = 0; (i < n); i++)
    {
        if ((c[i].aa == -1) || (c[i].ab == -1))
        {
            fprintf(stderr, "Warning dot %d has no connections\n", i + 1);
        }
        fprintf(fp, "CONECT%5d%5d%5d\n", i + 1, c[i].aa + 1, c[i].ab + 1);
    }
    gmx_ffclose(fp);
    sfree(c);
}

/*! \brief
 * Plots the surface into a PDB file, optionally including the original atoms.
 */
void connolly_plot(const char*  fn,
                   int          ndots,
                   const real   dots[],
                   rvec         x[],
                   t_atoms*     atoms,
                   t_symtab*    symtab,
                   PbcType      pbcType,
                   const matrix box,
                   gmx_bool     bIncludeSolute)
{
    const char* const atomnm = "DOT";
    const char* const resnm  = "DOT";
    const char* const title  = "Connolly Dot Surface Generated by gmx sasa";

    rvec* xnew = nullptr;

    if (bIncludeSolute)
    {
        int i0 = atoms->nr;
        int r0 = atoms->nres;
        srenew(atoms->atom, atoms->nr + ndots);
        memset(&atoms->atom[i0], 0, sizeof(*atoms->atom) * ndots);
        srenew(atoms->atomname, atoms->nr + ndots);
        srenew(atoms->resinfo, r0 + 1);
        atoms->atom[i0].resind = r0;
        t_atoms_set_resinfo(atoms, i0, symtab, resnm, r0 + 1, ' ', 0, ' ');
        if (atoms->pdbinfo != nullptr)
        {
            srenew(atoms->pdbinfo, atoms->nr + ndots);
        }
        snew(xnew, atoms->nr + ndots);
        for (int i = 0; (i < atoms->nr); i++)
        {
            copy_rvec(x[i], xnew[i]);
        }
        int k = 0;
        for (int i = 0; (i < ndots); i++)
        {
            int ii0                 = i0 + i;
            atoms->atomname[ii0]    = put_symtab(symtab, atomnm);
            atoms->atom[ii0].resind = r0;
            xnew[ii0][XX]           = dots[k++];
            xnew[ii0][YY]           = dots[k++];
            xnew[ii0][ZZ]           = dots[k++];
            if (atoms->pdbinfo != nullptr)
            {
                atoms->pdbinfo[ii0].type   = PdbRecordType::Atom;
                atoms->pdbinfo[ii0].atomnr = ii0;
                atoms->pdbinfo[ii0].bfac   = 0.0;
                atoms->pdbinfo[ii0].occup  = 0.0;
            }
        }
        atoms->nr   = i0 + ndots;
        atoms->nres = r0 + 1;
        write_sto_conf(fn, title, atoms, xnew, nullptr, pbcType, const_cast<rvec*>(box));
        atoms->nres = r0;
        atoms->nr   = i0;
    }
    else
    {
        t_atoms aaa;
        init_t_atoms(&aaa, ndots, TRUE);
        aaa.atom[0].resind = 0;
        t_atoms_set_resinfo(&aaa, 0, symtab, resnm, 1, ' ', 0, ' ');
        snew(xnew, ndots);
        int k = 0;
        for (int i = 0; (i < ndots); i++)
        {
            int ii0                 = i;
            aaa.atomname[ii0]       = put_symtab(symtab, atomnm);
            aaa.pdbinfo[ii0].type   = PdbRecordType::Atom;
            aaa.pdbinfo[ii0].atomnr = ii0;
            aaa.atom[ii0].resind    = 0;
            xnew[ii0][XX]           = dots[k++];
            xnew[ii0][YY]           = dots[k++];
            xnew[ii0][ZZ]           = dots[k++];
            aaa.pdbinfo[ii0].bfac   = 0.0;
            aaa.pdbinfo[ii0].occup  = 0.0;
        }
        aaa.nr = ndots;
        write_sto_conf(fn, title, &aaa, xnew, nullptr, pbcType, const_cast<rvec*>(box));
        do_conect(fn, ndots, xnew);
        done_atom(&aaa);
    }
    sfree(xnew);
}

/********************************************************************
 * Actual analysis module
 */

/*! \brief
 * Implements `gmx sas` trajectory analysis module.
 */
class Sasa : public TrajectoryAnalysisModule
{
public:
    Sasa();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;

    TrajectoryAnalysisModuleDataPointer startFrames(const AnalysisDataParallelOptions& opt,
                                                    const SelectionCollection& selections) override;
    void                                analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    /*! \brief
     * Surface areas as a function of time.
     *
     * First column is for the calculation group, and the rest for the
     * output groups.  This data is always produced.
     */
    AnalysisData area_;
    /*! \brief
     * Per-atom surface areas as a function of time.
     *
     * Contains one data set for each column in `area_`.
     * Each column corresponds to a selection position in `surfaceSel_`.
     * This data is only produced if atom or residue areas have been
     * requested.
     */
    AnalysisData atomArea_;
    /*! \brief
     * Per-residue surface areas as a function of time.
     *
     * Contains one data set for each column in `area_`.
     * Each column corresponds to a distinct residue `surfaceSel_`.
     * For example, if `surfaceSel_` selects residues 2, 5, and 7, there
     * will be three columns here.
     * This data is only produced if atom or residue areas have been
     * requested.
     */
    AnalysisData residueArea_;
    /*! \brief
     * Free energy estimates as a function of time.
     *
     * Column layout is the same as for `area_`.
     * This data is only produced if the output is requested.
     */
    AnalysisData dgSolv_;
    /*! \brief
     * Total volume and density of the calculation group as a function of
     * time.
     *
     * The first column is the volume and the second column is the density.
     * This data is only produced if the output is requested.
     */
    AnalysisData volume_;

    /*! \brief
     * The selection to calculate the surface for.
     *
     * Selection::originalId() and Selection::mappedId() store the mapping
     * from the positions to the columns of `residueArea_`.
     * The selection is computed with SelectionOption::dynamicMask(), i.e.,
     * even in the presence of a dynamic selection, the number of returned
     * positions is fixed, and SelectionPosition::selected() is used.
     */
    Selection surfaceSel_;
    /*! \brief
     * List of optional additional output groups.
     *
     * Each of these must be a subset of the `surfaceSel_`.
     * Selection::originalId() and Selection::mappedId() store the mapping
     * from the positions to the corresponsing positions in `surfaceSel_`.
     */
    SelectionList outputSel_;

    std::string fnArea_;
    std::string fnAtomArea_;
    std::string fnResidueArea_;
    std::string fnDGSolv_;
    std::string fnVolume_;
    std::string fnConnolly_;

    double solsize_;
    int    ndots_;
    // double                  minarea_;
    double dgsDefault_;
    bool   bIncludeSolute_;

    //! Global topology corresponding to the input.
    gmx_mtop_t* mtop_;
    //! Per-atom data corresponding to the input.
    AtomsDataPtr atoms_;
    //! Combined VdW and probe radii for each atom in the calculation group.
    std::vector<real> radii_;
    /*! \brief
     * Solvation free energy coefficients for each atom in the calculation
     * group.
     *
     * Empty if the free energy output has not been requested.
     */
    std::vector<real> dgsFactor_;
    //! Calculation algorithm.
    SurfaceAreaCalculator calculator_;

    // Copy and assign disallowed by base.
};

Sasa::Sasa() :
    solsize_(0.14), ndots_(24), dgsDefault_(0), bIncludeSolute_(true), mtop_(nullptr), atoms_(nullptr)
{
    // minarea_ = 0.5;
    registerAnalysisDataset(&area_, "area");
    registerAnalysisDataset(&atomArea_, "atomarea");
    registerAnalysisDataset(&residueArea_, "resarea");
    registerAnalysisDataset(&dgSolv_, "dgsolv");
    registerAnalysisDataset(&volume_, "volume");
}

void Sasa::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] computes solvent accessible surface areas.",
        "See Eisenhaber F, Lijnzaad P, Argos P, Sander C, & Scharf M",
        "(1995) J. Comput. Chem. 16, 273-284 for the algorithm used.",
        "With [TT]-q[tt], the Connolly surface can be generated as well",
        "in a [REF].pdb[ref] file where the nodes are represented as atoms",
        "and the edges connecting the nearest nodes as CONECT records.",
        "[TT]-odg[tt] allows for estimation of solvation free energies",
        "from per-atom solvation energies per exposed surface area.[PAR]",

        "The program requires a selection for the surface calculation to be",
        "specified with [TT]-surface[tt]. This should always consist of all",
        "non-solvent atoms in the system. The area of this group is always",
        "calculated. Optionally, [TT]-output[tt] can specify additional",
        "selections, which should be subsets of the calculation group.",
        "The solvent-accessible areas for these groups are also extracted",
        "from the full surface.[PAR]",

        "The average and standard deviation of the area over the trajectory",
        "can be calculated per residue and atom (options [TT]-or[tt] and",
        "[TT]-oa[tt]).[PAR]",
        //"In combination with the latter option an [REF].itp[ref] file can be",
        //"generated (option [TT]-i[tt])",
        //"which can be used to restrain surface atoms.[PAR]",

        "With the [TT]-tv[tt] option the total volume and density of the",
        "molecule can be computed. With [TT]-pbc[tt] (the default), you",
        "must ensure that your molecule/surface group is not split across PBC.",
        "Otherwise, you will get non-sensical results.",
        "Please also consider whether the normal probe radius is appropriate",
        "in this case or whether you would rather use, e.g., 0. It is good",
        "to keep in mind that the results for volume and density are very",
        "approximate. For example, in ice Ih, one can easily fit water molecules in the",
        "pores which would yield a volume that is too low, and surface area and density",
        "that are both too high."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .required()
                               .store(&fnArea_)
                               .defaultBasename("area")
                               .description("Total area as a function of time"));
    options->addOption(
            FileNameOption("odg")
                    .filetype(OptionFileType::Plot)
                    .outputFile()
                    .store(&fnDGSolv_)
                    .defaultBasename("dgsolv")
                    .description("Estimated solvation free energy as a function of time"));
    options->addOption(FileNameOption("or")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnResidueArea_)
                               .defaultBasename("resarea")
                               .description("Average area per residue"));
    options->addOption(FileNameOption("oa")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAtomArea_)
                               .defaultBasename("atomarea")
                               .description("Average area per atom"));
    options->addOption(FileNameOption("tv")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnVolume_)
                               .defaultBasename("volume")
                               .description("Total volume and density as a function of time"));
    options->addOption(FileNameOption("q")
                               .filetype(OptionFileType::PDB)
                               .outputFile()
                               .store(&fnConnolly_)
                               .defaultBasename("connolly")
                               .description("PDB file for Connolly surface"));

    options->addOption(
            DoubleOption("probe").store(&solsize_).description("Radius of the solvent probe (nm)"));
    options->addOption(IntegerOption("ndots").store(&ndots_).description(
            "Number of dots per sphere, more dots means more accuracy"));
    options->addOption(
            BooleanOption("prot").store(&bIncludeSolute_).description("Output the protein to the Connolly [REF].pdb[ref] file too"));
    options->addOption(
            DoubleOption("dgs").store(&dgsDefault_).description("Default value for solvation free energy per area (kJ/mol/nm^2)"));

    // Selections must select atoms for the VdW radii lookup to work.
    // The calculation group uses dynamicMask() so that the coordinates
    // match a static array of VdW radii.
    options->addOption(SelectionOption("surface")
                               .store(&surfaceSel_)
                               .required()
                               .onlySortedAtoms()
                               .dynamicMask()
                               .description("Surface calculation selection"));
    options->addOption(
            SelectionOption("output").storeVector(&outputSel_).onlySortedAtoms().multiValue().description("Output selection(s)"));

    // Atom names etc. are required for the VdW radii lookup.
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void Sasa::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    mtop_  = top.mtop();
    atoms_ = top.copyAtoms();

    // bITP   = opt2bSet("-i", nfile, fnm);
    const bool bResAt = !fnResidueArea_.empty() || !fnAtomArea_.empty(); // || bITP;
    const bool bDGsol = !fnDGSolv_.empty();

    if (solsize_ < 0)
    {
        solsize_ = 1e-3;
        fprintf(stderr, "Probe size too small, setting it to %g\n", solsize_);
    }
    if (ndots_ < 20)
    {
        ndots_ = 20;
        fprintf(stderr, "Ndots too small, setting it to %d\n", ndots_);
    }

    please_cite(stderr, "Eisenhaber95");

    if (bDGsol)
    {
        if (!top.hasFullTopology())
        {
            GMX_THROW(InconsistentInputError(
                    "Cannot compute Delta G of solvation without a tpr file"));
        }
        else
        {
            if (strcmp(*(atoms_->atomtype[0]), "?") == 0)
            {
                GMX_THROW(InconsistentInputError(
                        "Your input tpr file is too old (does not contain atom types). Cannot not "
                        "compute Delta G of solvation"));
            }
            else
            {
                printf("Free energy of solvation predictions:\n");
                please_cite(stdout, "Eisenberg86a");
            }
        }
    }

    // Now compute atomic radii including solvent probe size.
    // Also, fetch solvation free energy coefficients and
    // compute the residue indices that map the calculation atoms
    // to the columns of residueArea_.
    radii_.reserve(surfaceSel_.posCount());
    if (bDGsol)
    {
        dgsFactor_.reserve(surfaceSel_.posCount());
    }

    const int resCount = surfaceSel_.initOriginalIdsToGroup(top.mtop(), INDEX_RES);

    // TODO: Not exception-safe, but nice solution would be to have a C++
    // atom properties class...
    AtomProperties aps;

    ArrayRef<const int> atomIndices = surfaceSel_.atomIndices();
    int                 ndefault    = 0;
    for (int i = 0; i < surfaceSel_.posCount(); i++)
    {
        const int ii     = atomIndices[i];
        const int resind = atoms_->atom[ii].resind;
        real      radius = 0;
        if (!aps.setAtomProperty(epropVDW, *(atoms_->resinfo[resind].name), *(atoms_->atomname[ii]), &radius))
        {
            ndefault++;
        }
        radii_.push_back(radius + solsize_);
        if (bDGsol)
        {
            real dgsFactor = 0;
            if (!aps.setAtomProperty(
                        epropDGsol, *(atoms_->resinfo[resind].name), *(atoms_->atomtype[ii]), &dgsFactor))
            {
                dgsFactor = dgsDefault_;
            }
            dgsFactor_.push_back(dgsFactor);
        }
    }
    if (ndefault > 0)
    {
        fprintf(stderr, "WARNING: could not find a Van der Waals radius for %d atoms\n", ndefault);
    }

    // Pre-compute mapping from the output groups to the calculation group,
    // and store it in the selection ID map for easy lookup.
    for (size_t g = 0; g < outputSel_.size(); ++g)
    {
        ArrayRef<const int> outputIndices = outputSel_[g].atomIndices();
        for (int i = 0, j = 0; i < outputSel_[g].posCount(); ++i)
        {
            while (j < surfaceSel_.posCount() && outputIndices[i] > atomIndices[j])
            {
                ++j;
            }
            if (j == surfaceSel_.posCount() || outputIndices[i] != atomIndices[j])
            {
                const std::string message = formatString(
                        "Output selection '%s' is not a subset of "
                        "the surface selection (atom %d is the first "
                        "atom not in the surface selection)",
                        outputSel_[g].name(),
                        outputIndices[i] + 1);
                GMX_THROW(InconsistentInputError(message));
            }
            outputSel_[g].setOriginalId(i, j);
        }
    }

    calculator_.setDotCount(ndots_);
    calculator_.setRadii(radii_);

    // Initialize all the output data objects and initialize the output plotters.

    area_.setColumnCount(0, 1 + outputSel_.size());
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnArea_);
        plotm->setTitle("Solvent Accessible Surface");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Area (nm\\S2\\N)");
        plotm->appendLegend("Total");
        for (size_t i = 0; i < outputSel_.size(); ++i)
        {
            plotm->appendLegend(outputSel_[i].name());
        }
        area_.addModule(plotm);
    }

    if (bResAt)
    {
        atomArea_.setDataSetCount(1 + outputSel_.size());
        residueArea_.setDataSetCount(1 + outputSel_.size());
        for (size_t i = 0; i <= outputSel_.size(); ++i)
        {
            atomArea_.setColumnCount(i, surfaceSel_.posCount());
            residueArea_.setColumnCount(i, resCount);
        }
        {
            AnalysisDataAverageModulePointer avem(new AnalysisDataAverageModule);
            for (int i = 0; i < surfaceSel_.posCount(); ++i)
            {
                avem->setXAxisValue(i, surfaceSel_.position(i).atomIndices()[0] + 1);
            }
            atomArea_.addModule(avem);
            if (!fnAtomArea_.empty())
            {
                AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
                plotm->setFileName(fnAtomArea_);
                plotm->setTitle("Area per atom over the trajectory");
                plotm->setXLabel("Atom");
                plotm->setXFormat(8, 0);
                plotm->setYLabel("Area (nm\\S2\\N)");
                plotm->setErrorsAsSeparateColumn(true);
                plotm->appendLegend("Average (nm\\S2\\N)");
                plotm->appendLegend("Standard deviation (nm\\S2\\N)");
                avem->addModule(plotm);
            }
        }
        {
            AnalysisDataAverageModulePointer avem(new AnalysisDataAverageModule);
            int                              nextRow = 0;
            for (int i = 0; i < surfaceSel_.posCount(); ++i)
            {
                const int residueGroup = surfaceSel_.position(i).mappedId();
                if (residueGroup >= nextRow)
                {
                    GMX_ASSERT(residueGroup == nextRow,
                               "Inconsistent (non-uniformly increasing) residue grouping");
                    const int atomIndex    = surfaceSel_.position(i).atomIndices()[0];
                    const int residueIndex = atoms_->atom[atomIndex].resind;
                    avem->setXAxisValue(nextRow, atoms_->resinfo[residueIndex].nr);
                    ++nextRow;
                }
            }
            residueArea_.addModule(avem);
            if (!fnResidueArea_.empty())
            {
                AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
                plotm->setFileName(fnResidueArea_);
                plotm->setTitle("Area per residue over the trajectory");
                plotm->setXLabel("Residue");
                plotm->setXFormat(8, 0);
                plotm->setYLabel("Area (nm\\S2\\N)");
                plotm->setErrorsAsSeparateColumn(true);
                plotm->appendLegend("Average (nm\\S2\\N)");
                plotm->appendLegend("Standard deviation (nm\\S2\\N)");
                avem->addModule(plotm);
            }
        }
    }

    if (!fnDGSolv_.empty())
    {
        dgSolv_.setColumnCount(0, 1 + outputSel_.size());
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDGSolv_);
        plotm->setTitle("Free Energy of Solvation");
        plotm->setXAxisIsTime();
        plotm->setYLabel("D Gsolv");
        plotm->appendLegend("Total");
        for (size_t i = 0; i < outputSel_.size(); ++i)
        {
            plotm->appendLegend(outputSel_[i].name());
        }
        dgSolv_.addModule(plotm);
    }

    if (!fnVolume_.empty())
    {
        volume_.setColumnCount(0, 2);
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnVolume_);
        plotm->setTitle("Volume and Density");
        plotm->setXAxisIsTime();
        plotm->appendLegend("Volume (nm\\S3\\N)");
        plotm->appendLegend("Density (g/l)");
        volume_.addModule(plotm);
    }
}

/*! \brief
 * Temporary memory for use within a single-frame calculation.
 */
class SasaModuleData : public TrajectoryAnalysisModuleData
{
public:
    /*! \brief
     * Reserves memory for the frame-local data.
     *
     * `residueCount` will be zero if per-residue data is not being
     * calculated.
     */
    SasaModuleData(TrajectoryAnalysisModule*          module,
                   const AnalysisDataParallelOptions& opt,
                   const SelectionCollection&         selections,
                   int                                atomCount,
                   int                                residueCount) :
        TrajectoryAnalysisModuleData(module, opt, selections)
    {
        index_.reserve(atomCount);
        // If the calculation group is not dynamic, pre-calculate
        // the index, since it is not going to change.
        for (int i = 0; i < atomCount; ++i)
        {
            index_.push_back(i);
        }
        atomAreas_.resize(atomCount);
        res_a_.resize(residueCount);
    }

    void finish() override { finishDataHandles(); }

    //! Indices of the calculation selection positions selected for the frame.
    std::vector<int> index_;
    /*! \brief
     * Atom areas for each calculation selection position for the frame.
     *
     * One entry for each position in the calculation group.
     * Values for atoms not selected are set to zero.
     */
    std::vector<real> atomAreas_;
    /*! \brief
     * Working array to accumulate areas for each residue.
     *
     * One entry for each distinct residue in the calculation group;
     * indices are not directly residue numbers or residue indices.
     *
     * This vector is empty if residue area calculations are not being
     * performed.
     */
    std::vector<real> res_a_;
};

TrajectoryAnalysisModuleDataPointer Sasa::startFrames(const AnalysisDataParallelOptions& opt,
                                                      const SelectionCollection&         selections)
{
    return TrajectoryAnalysisModuleDataPointer(new SasaModuleData(
            this, opt, selections, surfaceSel_.posCount(), residueArea_.columnCount(0)));
}

/*! \brief
 * Helper method to compute the areas for a single selection.
 *
 * \param[in]  surfaceSel     The calculation selection.
 * \param[in]  sel            The selection to compute the areas for (can be
 *     `surfaceSel` or one of the output selections).
 * \param[in]  atomAreas      Atom areas for each position in `surfaceSel`.
 * \param[in]  dgsFactor      Free energy coefficients for each position in
 *     `surfaceSel`. If empty, free energies are not calculated.
 * \param[out] totalAreaOut Total area of `sel` (sum of atom areas it selects).
 * \param[out] dgsolvOut      Solvation free energy.
 *     Will be zero of `dgsFactor` is empty.
 * \param      atomAreaHandle Data handle to use for storing atom areas for `sel`.
 * \param      resAreaHandle  Data handle to use for storing residue areas for `sel`.
 * \param      resAreaWork    Work array for accumulating the residue areas.
 *     If empty, atom and residue areas are not calculated.
 *
 * `atomAreaHandle` and `resAreaHandle` are not used if `resAreaWork` is empty.
 */
void computeAreas(const Selection&         surfaceSel,
                  const Selection&         sel,
                  const std::vector<real>& atomAreas,
                  const std::vector<real>& dgsFactor,
                  real*                    totalAreaOut,
                  real*                    dgsolvOut,
                  AnalysisDataHandle       atomAreaHandle,
                  AnalysisDataHandle       resAreaHandle,
                  std::vector<real>*       resAreaWork)
{
    const bool bResAt    = !resAreaWork->empty();
    const bool bDGsolv   = !dgsFactor.empty();
    real       totalArea = 0;
    real       dgsolv    = 0;

    if (bResAt)
    {
        std::fill(resAreaWork->begin(), resAreaWork->end(), 0.0_real);
    }
    for (int i = 0; i < sel.posCount(); ++i)
    {
        // Get the index of the atom in the calculation group.
        // For the output groups, the mapping has been precalculated in
        // initAnalysis().
        const int ii = (sel != surfaceSel ? sel.position(i).mappedId() : i);
        if (!surfaceSel.position(ii).selected())
        {
            // For the calculation group, skip unselected atoms.
            if (sel == surfaceSel)
            {
                continue;
            }
            GMX_THROW(InconsistentInputError(
                    "Output selection is not a subset of the surface selection"));
        }
        // Get the internal index of the matching residue.
        // These have been precalculated in initAnalysis().
        const int  ri       = surfaceSel.position(ii).mappedId();
        const real atomArea = atomAreas[ii];
        totalArea += atomArea;
        if (bResAt)
        {
            atomAreaHandle.setPoint(ii, atomArea);
            (*resAreaWork)[ri] += atomArea;
        }
        if (bDGsolv)
        {
            dgsolv += atomArea * dgsFactor[ii];
        }
    }
    if (bResAt)
    {
        for (size_t i = 0; i < (*resAreaWork).size(); ++i)
        {
            resAreaHandle.setPoint(i, (*resAreaWork)[i]);
        }
    }
    *totalAreaOut = totalArea;
    *dgsolvOut    = dgsolv;
}

void Sasa::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   ah         = pdata->dataHandle(area_);
    AnalysisDataHandle   dgh        = pdata->dataHandle(dgSolv_);
    AnalysisDataHandle   aah        = pdata->dataHandle(atomArea_);
    AnalysisDataHandle   rah        = pdata->dataHandle(residueArea_);
    AnalysisDataHandle   vh         = pdata->dataHandle(volume_);
    const Selection&     surfaceSel = TrajectoryAnalysisModuleData::parallelSelection(surfaceSel_);
    const SelectionList& outputSel  = TrajectoryAnalysisModuleData::parallelSelections(outputSel_);
    SasaModuleData&      frameData  = *static_cast<SasaModuleData*>(pdata);

    const bool bResAt    = !frameData.res_a_.empty();
    const bool bDGsol    = !dgsFactor_.empty();
    const bool bConnolly = (frnr == 0 && !fnConnolly_.empty());

    // Update indices of selected atoms in the work array.
    if (surfaceSel.isDynamic())
    {
        frameData.index_.clear();
        for (int i = 0; i < surfaceSel.posCount(); ++i)
        {
            if (surfaceSel.position(i).selected())
            {
                frameData.index_.push_back(i);
            }
        }
    }

    // Determine what needs to be calculated.
    int flag = 0;
    if (bResAt || bDGsol || !outputSel.empty())
    {
        flag |= FLAG_ATOM_AREA;
    }
    if (bConnolly)
    {
        flag |= FLAG_DOTS;
    }
    if (volume_.columnCount() > 0)
    {
        flag |= FLAG_VOLUME;
    }

    // Do the low-level calculation.
    // totarea and totvolume receive the values for the calculation group.
    // area array contains the per-atom areas for the selected positions.
    // surfacedots contains nsurfacedots entries, and contains the actual
    // points.
    real  totarea   = 0;
    real  totvolume = 0;
    real *area = nullptr, *surfacedots = nullptr;
    int   nsurfacedots = 0;
    calculator_.calculate(surfaceSel.coordinates().data(),
                          pbc,
                          frameData.index_.size(),
                          frameData.index_.data(),
                          flag,
                          &totarea,
                          &totvolume,
                          &area,
                          &surfacedots,
                          &nsurfacedots);
    // Unpack the atomwise areas into the frameData.atomAreas_ array for easier
    // indexing in the case of dynamic surfaceSel.
    if (area != nullptr)
    {
        if (surfaceSel.isDynamic())
        {
            std::fill(frameData.atomAreas_.begin(), frameData.atomAreas_.end(), 0.0_real);
            for (size_t i = 0; i < frameData.index_.size(); ++i)
            {
                frameData.atomAreas_[frameData.index_[i]] = area[i];
            }
        }
        else
        {
            std::copy(area, area + surfaceSel.posCount(), frameData.atomAreas_.begin());
        }
        sfree(area);
    }
    const sfree_guard dotsGuard(surfacedots);

    if (bConnolly)
    {
        if (fr.natoms != mtop_->natoms)
        {
            GMX_THROW(
                    InconsistentInputError("Connolly plot (-q) is only supported for trajectories "
                                           "that contain all the atoms"));
        }
        // This is somewhat nasty, as it modifies the atoms and symtab
        // structures.  But since it is only used in the first frame, and no
        // one else uses the topology after initialization, it may just work
        // even with future parallelization.
        connolly_plot(fnConnolly_.c_str(),
                      nsurfacedots,
                      surfacedots,
                      fr.x,
                      atoms_.get(),
                      &mtop_->symtab,
                      fr.pbcType,
                      fr.box,
                      bIncludeSolute_);
    }

    ah.startFrame(frnr, fr.time);
    if (bResAt)
    {
        aah.startFrame(frnr, fr.time);
        rah.startFrame(frnr, fr.time);
    }
    if (bDGsol)
    {
        dgh.startFrame(frnr, fr.time);
    }

    ah.setPoint(0, totarea);

    real totalArea = 0;
    real dgsolv    = 0;
    if (bResAt || bDGsol)
    {
        computeAreas(surfaceSel,
                     surfaceSel,
                     frameData.atomAreas_,
                     dgsFactor_,
                     &totalArea,
                     &dgsolv,
                     aah,
                     rah,
                     &frameData.res_a_);
        if (bDGsol)
        {
            dgh.setPoint(0, dgsolv);
        }
    }
    for (size_t g = 0; g < outputSel.size(); ++g)
    {
        if (bResAt)
        {
            aah.selectDataSet(g + 1);
            rah.selectDataSet(g + 1);
        }
        computeAreas(surfaceSel,
                     outputSel[g],
                     frameData.atomAreas_,
                     dgsFactor_,
                     &totalArea,
                     &dgsolv,
                     aah,
                     rah,
                     &frameData.res_a_);
        ah.setPoint(g + 1, totalArea);
        if (bDGsol)
        {
            dgh.setPoint(g + 1, dgsolv);
        }
    }

    ah.finishFrame();
    if (bResAt)
    {
        aah.finishFrame();
        rah.finishFrame();
    }
    if (bDGsol)
    {
        dgh.finishFrame();
    }

    if (vh.isValid())
    {
        real totmass = 0;
        for (int i = 0; i < surfaceSel.posCount(); ++i)
        {
            totmass += surfaceSel.position(i).mass();
        }
        const real density = totmass * gmx::c_amu / (totvolume * gmx::c_nano * gmx::c_nano * gmx::c_nano);
        vh.startFrame(frnr, fr.time);
        vh.setPoint(0, totvolume);
        vh.setPoint(1, density);
        vh.finishFrame();
    }
}

void Sasa::finishAnalysis(int /*nframes*/) {}

void Sasa::writeOutput() {}

//! \}

} // namespace

const char SasaInfo::name[]             = "sasa";
const char SasaInfo::shortDescription[] = "Compute solvent accessible surface area";

TrajectoryAnalysisModulePointer SasaInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Sasa);
}

} // namespace analysismodules

} // namespace gmx
