/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_CONFIO_H
#define GMX_FILEIO_CONFIO_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

/* For reading coordinate files it is assumed that enough memory
 * has been allocated beforehand.
 */
struct gmx_mtop_t;
struct t_atoms;
struct t_symtab;
struct t_topology;
enum class PbcType : int;

void write_sto_conf_indexed(const char*    outfile,
                            const char*    title,
                            const t_atoms* atoms,
                            const rvec     x[],
                            const rvec*    v,
                            PbcType        pbcType,
                            const matrix   box,
                            int            nindex,
                            int            index[]);
/* like write_sto_conf, but indexed */

void write_sto_conf(const char*    outfile,
                    const char*    title,
                    const t_atoms* atoms,
                    const rvec     x[],
                    const rvec*    v,
                    PbcType        pbcType,
                    const matrix   box);
/* write atoms, x, v (if .gro and not NULL) and box (if not NULL)
 * to an STO (.gro or .pdb) file */

void write_sto_conf_mtop(const char*       outfile,
                         const char*       title,
                         const gmx_mtop_t* mtop,
                         const rvec        x[],
                         const rvec*       v,
                         PbcType           pbcType,
                         const matrix      box);
/* As write_sto_conf, but uses a gmx_mtop_t struct */

/*! \brief Read a configuration and, when available, a topology from a tpr or structure file.
 *
 * When reading from a tpr file, the complete topology is returned in \p mtop.
 * When reading from a structure file, only the atoms struct in \p mtop contains data.
 *
 * \param[in]     infile        Input file name
 * \param[out]    haveTopology  true when a topology was read and stored in mtop
 * \param[out]    mtop          The topology, either complete or only atom data
 * \param[out]    pbcType       Enum reporting the type of PBC
 * \param[in,out] x             Coordinates will be stored when *x!=NULL
 * \param[in,out] v             Velocities will be stored when *v!=NULL
 * \param[out]    box           Box dimensions
 */
void readConfAndTopology(const char* infile,
                         bool*       haveTopology,
                         gmx_mtop_t* mtop,
                         PbcType*    pbcType,
                         rvec**      x,
                         rvec**      v,
                         matrix      box);

/*! \brief Read a configuration from a structure file.
 *
 * This should eventually be superseded by TopologyInformation
 *
 * \param[in]     infile        Input file name
 * \param[out]    symtab        The symbol table
 * \param[out]    name          The title of the molecule, e.g. from pdb TITLE record
 * \param[out]    atoms         The global t_atoms struct
 * \param[out]    pbcType       Enum reporting the type of PBC
 * \param[in,out] x             Coordinates will be stored when *x!=NULL
 * \param[in,out] v             Velocities will be stored when *v!=NULL
 * \param[out]    box           Box dimensions
 */
void readConfAndAtoms(const char* infile,
                      t_symtab*   symtab,
                      char**      name,
                      t_atoms*    atoms,
                      PbcType*    pbcType,
                      rvec**      x,
                      rvec**      v,
                      matrix      box);

/*! \brief Read a configuration and, when available, a topology from a tpr or structure file.
 *
 * Deprecated, superseded by readConfAndTopology().
 * When \p requireMasses = TRUE, this routine must return a topology with
 * mass data. Masses are either read from a tpr input file, or otherwise
 * looked up from the mass database, and when such lookup fails a fatal error
 * results.
 * When \p requireMasses = FALSE, masses will still be read from tpr input and
 * their presence is signaled with the \p haveMass flag in t_atoms of \p top.
 *
 * \param[in]     infile        Input file name
 * \param[out]    top           The topology, either complete or only atom data. Caller is
 *                              responsible for calling done_top().
 * \param[out]    pbcType       Enum reporting the type of PBC
 * \param[in,out] x             Coordinates will be stored when *x!=NULL
 * \param[in,out] v             Velocities will be stored when *v!=NULL
 * \param[out]    box           Box dimensions
 * \param[in]     requireMasses Require masses to be present, either from tpr or from the mass
 * database \returns if a topology is available
 */
gmx_bool read_tps_conf(const char*        infile,
                       struct t_topology* top,
                       PbcType*           pbcType,
                       rvec**             x,
                       rvec**             v,
                       matrix             box,
                       gmx_bool           requireMasses);

#endif
