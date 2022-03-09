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

/*! \file
 * \brief Tools for converting simulation input data to and from TPR files.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi_compat
 */

#ifndef GMXAPICOMPAT_TPR_H
#define GMXAPICOMPAT_TPR_H

#include <memory>
#include <vector>

#include "gmxapi/compat/mdparams.h"
#include "gmxapi/gmxapicompat.h"

namespace gmxapicompat
{

/*!
 * \brief Manager for TPR file resources.
 *
 * To avoid copies, this resource-owning object is shared by consumers of its
 * resources, even when different resources are consumed.
 *
 * Multiple read-only handles may be issued if there are no write-handles.
 * One write handle may be issued if there are no other open handles.
 *
 * A const TprFile may only issue read file-handles, allowing handles to be
 * issued more quickly by avoiding atomic resource locking.
 *
 * \note Shared ownership of file manager could be avoided if owned by a Context.
 * It is appropriate for a Context to own and mediate access to the manager because
 * the Context should provide the filesystem abstraction to more intelligently
 * map named file paths to resources. For now, handles and other consumers share ownership
 * of the TprContents manager object via shared_ptr.
 */
class TprContents;

class TprReadHandle
{
public:
    explicit TprReadHandle(std::shared_ptr<TprContents> tprFile);
    explicit TprReadHandle(TprContents&& tprFile);
    TprReadHandle(const TprReadHandle&) = default;
    TprReadHandle& operator=(const TprReadHandle&) = default;
    TprReadHandle(TprReadHandle&&) noexcept        = default;
    TprReadHandle& operator=(TprReadHandle&&) noexcept = default;
    ~TprReadHandle();

    /*!
     * \brief Allow API functions to access data resources.
     *
     * Used internally. The entire TPR contents are never extracted to the
     * client, but API implementation details need to be
     * able to access some or all entire contents in later operations.
     *
     * \return Reference-counted handle to data container.
     */
    std::shared_ptr<TprContents> get() const;

private:
    std::shared_ptr<TprContents> tprContents_;
};

/*!
 * \brief Helper function for early implementation.
 *
 * Allows extraction of TPR file information from special params objects.
 *
 * \todo This is a very temporary shim! Find a better way to construct simulation input.
 */
TprReadHandle getSourceFileHandle(const GmxMdParams& params);

class StructureSource
{
public:
    std::shared_ptr<TprContents> tprFile_;
};

class TopologySource
{
public:
    std::shared_ptr<TprContents> tprFile_;
};

class SimulationState
{
public:
    std::shared_ptr<TprContents> tprFile_;
};

/*!
 * \brief Copy TPR file.
 *
 * \param input TPR source to copy from
 * \param outFile output TPR file name
 * \return true if successful. else false.
 */
bool copy_tprfile(const gmxapicompat::TprReadHandle& input, const std::string& outFile);

/*!
 * \brief Copy and possibly update TPR file by name.
 *
 * \param inFile Input file name
 * \param outFile Output file name
 * \param endTime Replace `nsteps` in infile with `endTime/dt`
 * \return true if successful, else false
 */
bool rewrite_tprfile(const std::string& inFile, const std::string& outFile, double endTime);

} // end namespace gmxapicompat

#endif // GMXAPICOMPAT_TPR_H
