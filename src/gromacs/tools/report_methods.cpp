/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include "report_methods.h"

#include <filesystem>
#include <memory>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileredirector.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class CommandLineModuleSettings;

void writeHeader(TextWriter* writer, const std::string& text, const std::string& section, bool writeFormattedText)
{
    std::string formattedText;
    if (writeFormattedText)
    {
        formattedText = "\\" + section + "{" + text + "}\n";
    }
    else
    {
        formattedText = section + ": " + text + "\n";
    }
    writer->writeString(formattedText);
}

void writeSystemInformation(TextWriter* writer, const gmx_mtop_t& top, bool writeFormattedText)
{
    int                       nmol, nvsite = 0;
    gmx_mtop_atomloop_block_t aloop;
    const t_atom*             atom;

    writeHeader(writer, "Simulation system", "subsection", writeFormattedText);
    aloop = gmx_mtop_atomloop_block_init(top);
    while (gmx_mtop_atomloop_block_next(aloop, &atom, &nmol))
    {
        if (atom->ptype == ParticleType::VSite)
        {
            nvsite += nmol;
        }
    }
    {
        writer->writeLine(formatString("A system of %d molecules (%d atoms) was simulated.",
                                       gmx_mtop_num_molecules(top),
                                       top.natoms - nvsite));
    }
    if (nvsite)
    {
        writer->writeLine(formatString("Virtual sites were used in some of the molecules."));
    }
    writer->ensureEmptyLine();
}

void writeParameterInformation(TextWriter* writer, const t_inputrec& ir, bool writeFormattedText)
{
    writeHeader(writer, "Simulation settings", "subsection", writeFormattedText);
    writer->writeLine(formatString("A total of %g ns were simulated with a time step of %g fs.",
                                   ir.nsteps * ir.delta_t * 0.001,
                                   1000 * ir.delta_t));
    writer->writeLine(formatString("Neighbor searching was performed every %d steps.", ir.nstlist));
    writer->writeLine(formatString("The %s algorithm was used for electrostatic interactions.",
                                   enumValueToString(ir.coulombtype)));
    writer->writeLine(formatString("with a cut-off of %g nm.", ir.rcoulomb));
    if (ir.coulombtype == CoulombInteractionType::Pme)
    {
        writer->writeLine(
                formatString("A reciprocal grid of %d x %d x %d cells was used with %dth order "
                             "B-spline interpolation.",
                             ir.nkx,
                             ir.nky,
                             ir.nkz,
                             ir.pme_order));
    }
    writer->writeLine(formatString(
            "A single cut-off of %g nm was used for Van der Waals interactions.", ir.rlist));
    if (ir.etc != TemperatureCoupling::No)
    {
        writer->writeLine(formatString("Temperature coupling was done with the %s algorithm.",
                                       enumValueToString(ir.etc)));
    }
    if (ir.pressureCouplingOptions.epc != PressureCoupling::No)
    {
        writer->writeLine(formatString("Pressure coupling was done with the %s algorithm.",
                                       enumValueToString(ir.pressureCouplingOptions.epc)));
    }
    writer->ensureEmptyLine();
}

void writeInformation(TextOutputFile*   outputStream,
                      const t_inputrec& ir,
                      const gmx_mtop_t& top,
                      bool              writeFormattedText,
                      bool              notStdout)
{
    TextWriter writer(outputStream);
    writer.ensureEmptyLine();
    writeHeader(&writer, "Methods", "section", writeFormattedText);
    writeSystemInformation(&writer, top, writeFormattedText);
    writeParameterInformation(&writer, ir, writeFormattedText);
    writer.ensureEmptyLine();

    if (notStdout)
    {
        writer.close();
    }
}

namespace
{

class ReportMethods : public ICommandLineOptionsModule
{
public:
    ReportMethods() : writeLatex_(false), writePlainText_(false) {}

    // From ICommandLineOptionsModule
    void init(CommandLineModuleSettings* /*settings*/) override {}
    void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings) override;
    void optionsFinished() override;
    int  run() override;

private:
    //! File name for the output LaTeX file or empty.
    std::string outputFileLatex_;
    //! File name for the unformatted output file or empty.
    std::string outputFileUnformatted_;
    //! File name of the run input file with full topology.
    std::string inputTopology_;
    //! Boolean reporting if writing to the LaTeX output file is requested.
    bool writeLatex_;
    //! Boolean reporting if writing to unformatted output is requested.
    bool writePlainText_;
};

void ReportMethods::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings* settings)
{
    const char* const desc[] = { "[THISMODULE] reports basic system information for the run input",
                                 "file specified with [TT]-s[tt] either to the",
                                 "terminal, to a LaTeX formatted output file if run with",
                                 "the [TT]-m[tt] option or to an unformatted file with",
                                 "the [TT]-o[tt] option.",
                                 "The functionality has been moved here from its previous",
                                 "place in [gmx-check]." };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("s")
                               .filetype(OptionFileType::Topology)
                               .inputFile()
                               .required()
                               .store(&inputTopology_)
                               .defaultBasename("topol")
                               .description("Run input file for report"));

    // TODO: Replace use of legacyType.
    options->addOption(FileNameOption("m")
                               .legacyType(efTEX)
                               .outputFile()
                               .store(&outputFileLatex_)
                               .storeIsSet(&writeLatex_)
                               .defaultBasename("report")
                               .description("LaTeX formatted report output"));
    options->addOption(FileNameOption("o")
                               .legacyType(efOUT)
                               .outputFile()
                               .store(&outputFileUnformatted_)
                               .storeIsSet(&writePlainText_)
                               .defaultBasename("report")
                               .description("Unformatted report output to file"));
}

void ReportMethods::optionsFinished() {}

int ReportMethods::run()
{
    t_state    state;
    t_inputrec ir;
    gmx_mtop_t top;
    read_tpx_state(inputTopology_, &ir, &state, &top);
    if (writeLatex_)
    {
        TextOutputFile file(outputFileLatex_);
        writeInformation(&file, ir, top, true, true);
    }
    if (writePlainText_)
    {
        TextOutputFile file(outputFileUnformatted_);
        writeInformation(&file, ir, top, false, true);
    }
    TextOutputFile& stdoutFile = TextOutputFile::standardOutput();
    writeInformation(&stdoutFile, ir, top, false, false);

    return 0;
}

} // namespace

LIBGROMACS_EXPORT const char ReportMethodsInfo::name[] = "report-methods";
LIBGROMACS_EXPORT const char ReportMethodsInfo::shortDescription[] =
        "Write short summary about the simulation setup to a text file "
        "and/or to the standard output.";
ICommandLineOptionsModulePointer ReportMethodsInfo::create()
{
    return ICommandLineOptionsModulePointer(std::make_unique<ReportMethods>());
}

} // namespace gmx
