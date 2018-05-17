/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "report.h"

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selectionoptionbehavior.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/fileredirector.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

namespace
{

class Report : public ICommandLineOptionsModule, public ITopologyProvider
{
    public:
        Report()
            : bLatex_(false), bWrite_(false), bFullTopology_(false), top_()
        {
        }

        // From ITopologyProvider
        virtual gmx_mtop_t *getTopology(bool /*required*/) { return &top_; }
        virtual int getAtomCount() {return 0; }

        // From ICommandLineOptionsModule
        virtual void init(CommandLineModuleSettings * /*settings*/)
        {
        }
        virtual void initOptions(IOptionsContainer                 *options,
                                 ICommandLineOptionsModuleSettings *settings);
        virtual void optionsFinished();
        virtual int run();

    private:

        /*! \brief
         * Write appropiate Header to output stream.
         *
         * \param[in] writer TextWriter object for writing information.
         * \param[in] text String with the header before writing.
         * \param[in] section String with section text for header.
         */
        void writeHeader(TextWriter *writer, std::string text, std::string section);

        /*! \brief
         * Wrapper for writing out information.
         *
         * This function is actual called from within the run method
         * to write the information to the terminal or to file.
         * New write out methods should be added to it instead of adding them in run.
         */
        void writeInformation();

        /*! \brief
         * Write information about system parameters.
         *
         * This method writes the basic information for the system parameters
         * and simulation settings as reported in the \p ir.
         *
         * \param[in] writer TextWriter object for writing information.
         * \param[in] ir Pointer to inputrec of the run input.
         */
        void writeParameterInformation(TextWriter *writer, const t_inputrec *ir);

        /*! \brief
         * Write information about the molecules in the system.
         *
         * This method should write all possible information about
         * the molecular composition of the system.
         *
         * \param[in] writer TextWriter object for writing information.
         *
         */
        void writeSystemInformation(TextWriter *writer);

        //! File name for the output LaTeX file or empty.
        std::string         outputFileLatex_;
        //! File name for the unformatted output file or empty.
        std::string         outputFIleUnformatted_;
        //! File name of the run input file with full topology.
        std::string         inputTopology_;
        //! Boolean reporting if writing to the LaTeX output file is requested.
        bool                bLatex_;
        //! Boolean reporting if writing to unformatted output is requested.
        bool                bWrite_;
        //! Boolean to check if the run input has a full topology.
        bool                bFullTopology_;
        //! Topology loaded from run input.
        gmx_mtop_t          top_;

};

void Report::initOptions(IOptionsContainer                 *options,
                         ICommandLineOptionsModuleSettings *settings)
{
    const char *const desc[] = {
        "[THISMODULE] reports basic system information for the run input",
        "file specfied with [TT]-s[tt] either to the",
        "terminal, to a LaTeX formatted output file if run with",
        "the [TT]-m[tt] option or to an unformatted file with",
        "the [TT]-o[tt] option.",
        "The functionality has been moved here from its previous",
        "place in [gmx-check]."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("s")
                           .filetype(eftTopology).inputFile().required()
                           .store(&inputTopology_)
                           .defaultBasename("topol")
                           .description("Run input file for report"));

    // TODO: Replace use of legacyType.
    options->addOption(FileNameOption("m")
                           .legacyType(efTEX).outputFile()
                           .store(&outputFileLatex_).storeIsSet(&bLatex_)
                           .defaultBasename("report")
                           .description("LaTeX formatted report output"));
    options->addOption(FileNameOption("o")
                           .legacyType(efOUT).outputFile()
                           .store(&outputFIleUnformatted_).storeIsSet(&bWrite_)
                           .defaultBasename("report")
                           .description("Unformatted report output to file"));
}

void Report::optionsFinished()
{
    if (inputTopology_.empty())
    {
        GMX_THROW(InconsistentInputError("Need a run input file to write report."));
    }
    {
        matrix box;
        readConfAndTopology(inputTopology_.c_str(), &bFullTopology_, &top_, nullptr,
                            nullptr, nullptr, box);
        if (!bFullTopology_)
        {
            GMX_THROW(InconsistentInputError("Need a complete run input file to write report"));
        }
    }
    if (bLatex_ && bWrite_)
    {
        GMX_THROW(InconsistentInputError("Can't write both formatted and unformatted output"));
    }
}

void Report::writeHeader(TextWriter *writer, std::string text, std::string section)
{
    std::string formattedText;
    if (bLatex_)
    {
        formattedText = "\\" + section + "{" + text + "}\n";
    }
    else
    {
        formattedText = section + ": " + text + "\n";
    }
    writer->writeString(formattedText);
}

void Report::writeSystemInformation(TextWriter *writer)
{
    int                       nmol, nvsite = 0;
    gmx_mtop_atomloop_block_t aloop;
    const t_atom             *atom;

    writeHeader(writer, "Simulation system", "subsection");
    aloop = gmx_mtop_atomloop_block_init(&top_);
    while (gmx_mtop_atomloop_block_next(aloop, &atom, &nmol))
    {
        if (atom->ptype == eptVSite)
        {
            nvsite += nmol;
        }
    }
    {
        std::string text = formatString("A system of %d molecules (%d atoms) was simulated.\n",
                                        gmx_mtop_num_molecules(top_), top_.natoms-nvsite);
        writer->writeLine(text);
    }
    if (nvsite)
    {
        std::string text = formatString("Virtual sites were used in some of the molecules.\n");
        writer->writeLine(text);
    }
    writer->ensureEmptyLine();
}

void Report::writeParameterInformation(TextWriter *writer, const t_inputrec *ir)
{
    writeHeader(writer, "Simulation settings", "subsection");
    std::string text;
    text = formatString("A total of %g ns were simulated with a time step of %g fs.\n",
                        ir->nsteps*ir->delta_t*0.001, 1000*ir->delta_t);
    writer->writeLine(text);
    text = formatString("Neighbor searching was performed every %d steps.\n", ir->nstlist);
    writer->writeLine(text);
    text = formatString("The %s algorithm was used for electrostatic interactions.\n",
                        EELTYPE(ir->coulombtype));
    writer->writeLine(text);
    text = formatString("with a cut-off of %g nm.\n", ir->rcoulomb);
    writer->writeLine(text);
    if (ir->coulombtype == eelPME)
    {
        text = formatString("A reciprocal grid of %d x %d x %d cells was used with %dth order B-spline interpolation.\n", ir->nkx, ir->nky, ir->nkz, ir->pme_order);
        writer->writeLine(text);
    }
    if (ir->rvdw > ir->rlist)
    {
        text = formatString("A twin-range Van der Waals cut-off (%g/%g nm) was used, where the long range forces were updated during neighborsearching.\n", ir->rlist, ir->rvdw);
        writer->writeLine(text);
    }
    else
    {
        text = formatString("A single cut-off of %g was used for Van der Waals interactions.\n", ir->rlist);
        writer->writeLine(text);
    }
    if (ir->etc != 0)
    {
        text = formatString("Temperature coupling was done with the %s algorithm.\n",
                            etcoupl_names[ir->etc]);
        writer->writeLine(text);
    }
    if (ir->epc != 0)
    {
        text = formatString("Pressure coupling was done with the %s algorithm.\n",
                            epcoupl_names[ir->epc]);
        writer->writeLine(text);
    }
    writer->ensureEmptyLine();
}

void Report::writeInformation()
{
    t_state        state;
    gmx_mtop_t     mtop;

    t_inputrec     ir;
    read_tpx_state(inputTopology_.c_str(), &ir, &state, &mtop);

    TextOutputFile         &stdoutFile = TextOutputFile::standardOutput();
    TextOutputStreamPointer test;
    TextOutputFile         *file              = nullptr;
    bool                    writeOutputToFile = (bLatex_ || bWrite_);
    std::string             outputFile;
    if (bLatex_)
    {
        outputFile = outputFileLatex_;
    }
    else if (bWrite_)
    {
        outputFile = outputFIleUnformatted_;
    }
    TextWriter              writer(writeOutputToFile ? file = (new TextOutputFile(outputFile)) : &stdoutFile);
    writer.ensureEmptyLine();
    writeHeader(&writer, "Methods", "section");
    writeSystemInformation(&writer);
    writeParameterInformation(&writer, &ir);
    writer.ensureEmptyLine();
    if (bLatex_)
    {
        writer.close();
        delete(file);
    }
}

int Report::run()
{
    writeInformation();

    return 0;
}

}   // namespace

const char ReportInfo::name[]             = "report";
const char ReportInfo::shortDescription[] =
    "Write short summary about the simulation setup to a text file";
ICommandLineOptionsModulePointer ReportInfo::create()
{
    return ICommandLineOptionsModulePointer(new Report());
}

} // namespace gmx
