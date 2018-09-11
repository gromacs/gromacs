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
/*! \libinternal \file
 * \brief Define the filename options to MD simulations.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrun
 */

#ifndef GROMACS_MDFILENAMES_H
#define GROMACS_MDFILENAMES_H

#include <array>

#include "gromacs/commandline/filenm.h"

namespace gmx
{

/*! \libinternal
 * \brief MD filename options.
 *
 * Default construction provides the initial container of filename options used
 * by gmx::Mdrunner.
 *
 * \ingroup module_mdrun
 * \inlibraryapi
 */
class MdFilenames
{
    public:

        /*!
         * \brief Access stored data.
         *
         * \return reference to internal data structure.
         * \{
         */
        std::array<t_filenm, 34> &operator()() noexcept { return filenames_; }
        const std::array<t_filenm, 34> &operator()() const noexcept { return filenames_; }
        //! \}

    private:
        //! Filenames and properties from command-line argument values.
        std::array<t_filenm, 34> filenames_ =
        {{{ efTPR, nullptr,     nullptr,     ffREAD },
          { efTRN, "-o",        nullptr,     ffWRITE },
          { efCOMPRESSED, "-x", nullptr,     ffOPTWR },
          { efCPT, "-cpi",      nullptr,     ffOPTRD | ffALLOW_MISSING },
          { efCPT, "-cpo",      nullptr,     ffOPTWR },
          { efSTO, "-c",        "confout",   ffWRITE },
          { efEDR, "-e",        "ener",      ffWRITE },
          { efLOG, "-g",        "md",        ffWRITE },
          { efXVG, "-dhdl",     "dhdl",      ffOPTWR },
          { efXVG, "-field",    "field",     ffOPTWR },
          { efXVG, "-table",    "table",     ffOPTRD },
          { efXVG, "-tablep",   "tablep",    ffOPTRD },
          { efXVG, "-tableb",   "table",     ffOPTRDMULT },
          { efTRX, "-rerun",    "rerun",     ffOPTRD },
          { efXVG, "-tpi",      "tpi",       ffOPTWR },
          { efXVG, "-tpid",     "tpidist",   ffOPTWR },
          { efEDI, "-ei",       "sam",       ffOPTRD },
          { efXVG, "-eo",       "edsam",     ffOPTWR },
          { efXVG, "-devout",   "deviatie",  ffOPTWR },
          { efXVG, "-runav",    "runaver",   ffOPTWR },
          { efXVG, "-px",       "pullx",     ffOPTWR },
          { efXVG, "-pf",       "pullf",     ffOPTWR },
          { efXVG, "-ro",       "rotation",  ffOPTWR },
          { efLOG, "-ra",       "rotangles", ffOPTWR },
          { efLOG, "-rs",       "rotslabs",  ffOPTWR },
          { efLOG, "-rt",       "rottorque", ffOPTWR },
          { efMTX, "-mtx",      "nm",        ffOPTWR },
          { efRND, "-multidir", nullptr,     ffOPTRDMULT},
          { efXVG, "-awh",      "awhinit",   ffOPTRD },
          { efDAT, "-membed",   "membed",    ffOPTRD },
          { efTOP, "-mp",       "membed",    ffOPTRD },
          { efNDX, "-mn",       "membed",    ffOPTRD },
          { efXVG, "-if",       "imdforces", ffOPTWR },
          { efXVG, "-swap",     "swapions",  ffOPTWR }}};

};

}      // end namespace gmx

#endif //GROMACS_MDFILENAMES_H
