/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \internal
 * \brief Encapsulates membed methods
 *
 * \author Joe Jordan <ejjordan@kth.se>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_MEMBEDHOLDER_H
#define GMX_MDRUN_MEMBEDHOLDER_H

#include <cstdio>

#include "gromacs/utility/real.h"

struct gmx_membed_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
class t_state;

namespace gmx
{

/*! \brief Membed SimulatorBuilder parameter type.
 *
 * Does not (yet) encapsulate ownership semantics of resources. Simulator is
 * not (necessarily) granted ownership of resources. Client is responsible for
 * maintaining the validity of resources for the life time of the Simulator,
 * then for cleaning up those resources.
 */
class MembedHolder
{
public:
    //! Build holder from input information.
    explicit MembedHolder(int nfile, const t_filenm fnm[]);
    //! Move is possible.
    MembedHolder(MembedHolder&& holder) noexcept;
    //! Move assignment is possible.
    MembedHolder& operator=(MembedHolder&& holder) noexcept;
    //! Copy is not allowed.
    MembedHolder(const MembedHolder&) = delete;
    //! Copy assignment is not allowed.
    MembedHolder& operator=(const MembedHolder&) = delete;

    ~MembedHolder();

    //! Get information about membed being used.
    [[nodiscard]] bool doMembed() const { return doMembed_; }

    /*! \brief
     * Fully initialize underlying datastructure.
     *
     * \param[in] fplog Handle to log file.
     * \param[in] nfile How many input files are there.
     * \param[in] fnm   Input file collection datastructure.
     * \param[in,out] mtop  Handle to mtop, can be modified.
     * \param[in,out] inputrec Handle to inputrec, can be modified.
     * \param[in,out] state    Simulation state information, can be modified.
     * \param[in,out] cr       Communication information.
     * \param[out]    cpt      Some kind of checkpoint information.
     */
    void initializeMembed(FILE*          fplog,
                          int            nfile,
                          const t_filenm fnm[],
                          gmx_mtop_t*    mtop,
                          t_inputrec*    inputrec,
                          t_state*       state,
                          t_commrec*     cr,
                          real*          cpt);

    //! Get handle to membed object.
    gmx_membed_t* membed();

private:
    //! Pointer to membed object. TODO replace with unique_ptr or get rid of this.
    gmx_membed_t* membed_ = nullptr;
    //! Whether membed is being used.
    bool doMembed_ = false;
};

} // namespace gmx

#endif // GMX_MDRUN_MEMBEDHOLDER_H
