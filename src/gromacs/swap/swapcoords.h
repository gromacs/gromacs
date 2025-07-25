/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
/*! \libinternal
 * \defgroup module_swap "Computational Electrophysiology" position swapping (swap)
 * \ingroup group_mdrun
 * \brief
 * Implements the "Computational Electrophysiology" protocol.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 */
/*! \libinternal \file
 * \brief
 * The "Computational Electrophysiology" protocol for ion/water position swapping.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \inlibraryapi
 * \ingroup module_swap
 */
#ifndef GMX_SWAP_SWAPCOORDS_H
#define GMX_SWAP_SWAPCOORDS_H

#include <cstdint>
#include <cstdio>

#include <memory>
#include <string_view>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/vectypes.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct swaphistory_t;
struct t_inputrec;
class t_state;
struct t_swapcoords;
struct ObservablesHistory;

namespace gmx
{
template<typename T>
class ArrayRef;
class MpiComm;
enum class StartingBehavior;
class IMDModule;
class LocalAtomSetManager;
struct MdrunOptions;

/*! \internal
    \brief Information about the computational-electrophysiology module.
 *
 * Provides name and method to create a computational-electrophysiology module.
 */
struct SwapCoordinatesModuleInfo
{
    /*! \brief
     * Creates a module for Computational Electrophysiology swapping.
     */
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static constexpr std::string_view sc_name = "swap-coordinates";
};

} // namespace gmx

class SwapCoords
{
public:
    SwapCoords();
    ~SwapCoords();
    //! Impl class, currently public while module evolves
    class Impl;
    //! Impl object, currently public while module evolves
    std::unique_ptr<Impl> impl_;
};

/*! \brief Initialize ion / water position swapping ("Computational Electrophysiology").
 *
 * This routine does the memory allocation for various helper arrays, opens
 * the output file, sets up swap data checkpoint writing, etc. and returns it.
 *
 * \param[in] fplog         General output file, normally md.log.
 * \param[in] ir            Structure containing MD input parameters, among those
 *                          also the structure needed for position swapping.
 * \param[in] fn            Output file name for swap data.
 * \param[in] mtop          Molecular topology.
 * \param[in] globalState   The global state, only used on the main rank.
 * \param[in] oh            Contains struct with swap data that is read from or written to checkpoint.
 * \param[in] mpiComm       Communicator object for my group.
 * \param[in] dd            Domain decomposition object, is nullptr when DD is not active.
 * \param[in] atomSets      Manager tending to swap atom indices.
 * \param[in] oenv          Needed to open the swap output XVGR file.
 * \param[in] mdrunOptions  Options for mdrun.
 * \param[in] startingBehavior  Describes whether this is a restart appending to output files
 */
std::unique_ptr<SwapCoords> init_swapcoords(FILE*                     fplog,
                                            const t_inputrec*         ir,
                                            const char*               fn,
                                            const gmx_mtop_t&         mtop,
                                            const t_state*            globalState,
                                            ObservablesHistory*       oh,
                                            const gmx::MpiComm&       mpiComm,
                                            const gmx_domdec_t*       dd,
                                            gmx::LocalAtomSetManager* atomSets,
                                            const gmx_output_env_t*   oenv,
                                            const gmx::MdrunOptions&  mdrunOptions,
                                            gmx::StartingBehavior     startingBehavior);


/*! \brief "Computational Electrophysiology" main routine within MD loop.
 *
 * \param[in] mpiComm  Communicator object for my group.
 * \param[in] step     The number of the MD time step.
 * \param[in] t        The time.
 * \param[in] ir       Structure containing MD input parameters
 * \param[in,out] s    The structure needed for position swapping.
 * \param[in] wcycle   Count wallcycles of swap routines for diagnostic output.
 * \param[in] x        Positions of home particles this node owns.
 * \param[in] box      The simulation box.
 * \param[in] bVerbose Should we be quiet or verbose?
 * \param[in] bRerun   Are we doing a rerun?
 *
 * \returns Whether at least one pair of molecules was swapped.
 */
gmx_bool do_swapcoords(const gmx::MpiComm&      mpiComm,
                       int64_t                  step,
                       double                   t,
                       const t_inputrec*        ir,
                       SwapCoords*              s,
                       gmx_wallcycle*           wcycle,
                       gmx::ArrayRef<gmx::RVec> x,
                       matrix                   box,
                       gmx_bool                 bVerbose,
                       gmx_bool                 bRerun);

#endif
