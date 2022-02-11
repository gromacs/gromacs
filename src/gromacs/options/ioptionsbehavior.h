/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \file
 * \brief
 * Declares gmx::IOptionsBehavior.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_IOPTIONSBEHAVIOR_H
#define GMX_OPTIONS_IOPTIONSBEHAVIOR_H

namespace gmx
{

class Options;

/*! \brief
 * Interface to provide extension points for options parsing.
 *
 * Currently, this is only used in the context of ICommandLineOptionsModule and
 * some other command-line handling, but it is declared in the options module
 * for the lack of a better place: most implementations of the interface are in
 * modules that do not otherwise depend on the commandline module.
 *
 * \if libapi
 * Any code that wants to support these extension points needs to use
 * OptionsBehaviorCollection and call the methods there at appropriate points.
 * This is not (at least, not currently) integrated in any automatic way to the
 * actual Options object.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_options
 */
class IOptionsBehavior
{
public:
    virtual ~IOptionsBehavior();

    /*! \brief
     * Called when the behavior is associated with an options object.
     *
     * This method can, e.g., use Options::addManager() to associate
     * managers with the options object.
     */
    virtual void initBehavior(Options* options) = 0;
    /*! \brief
     * Called when all option values have been assigned.
     *
     * This is called just before Options::finish(), and can, e.g., do
     * operations that still influence the option values.
     */
    virtual void optionsFinishing(Options* options) = 0;
    /*! \brief
     * Called when all option values have been processed.
     *
     * This is called after Options::finish() (and typically after
     * higher-level optionsFinished() methods, such as that in
     * ICommandLineOptionsModule).  This can add behavior that performs
     * tasks based on the option values provided.
     */
    virtual void optionsFinished() = 0;
};

} // namespace gmx

#endif
