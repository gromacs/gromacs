/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Reading and writing of macromolecular structures in the mmCIF/PDBx format.
 *
 * This module wraps the RCSB mmCIF parser (bundled under
 * src/external/mmcif-dict-suite-*) to import/export the \c atom_site records
 * of an mmCIF/PDBx file into/from a GROMACS ::t_atoms structure, and to extract
 * the (possibly extended, > 4 character) PDB entry ID.
 *
 * The interface deliberately mirrors the legacy pdbio.h free functions and is
 * kept self-contained (no mmCIF-library types appear here) so that it can serve
 * as a drop-in structure-format plugin for a future unified I/O framework.
 *
 * \author Amr Alhossary
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_CIFIO_H
#define GMX_FILEIO_CIFIO_H

#include <cstdio>

#include <filesystem>

#include "gromacs/utility/vectypes.h"

enum class PbcType : int;
struct t_atoms;
struct t_symtab;
struct t_trxframe;
//! Opaque parsed-mmCIF trajectory handle; defined in cifio.cpp.
struct CifTrajectory;
/*! \brief Bond list, as used for PDB CONECT records; defined in pdbio.cpp.
 *
 * Declared here rather than including pdbio.h so that this header stays
 * self-contained. The pdbio.h typedef \c gmx_conect is a pointer to this type,
 * so callers can pass their gmx_conect unchanged. */
struct gmx_conect_t;

/*! \brief Reader/writer for mmCIF/PDBx structure files.
 *
 * All members are static; the class only groups the CIF I/O entry points and
 * gives them a natural home to grow into a polymorphic format handler later.
 */
class CifIO
{
public:
    /*! \brief Read \c atom_site (and cell) records from an mmCIF/PDBx file.
     *
     * Fills \p atoms, \p x, \p box and \p pbcType exactly like
     * gmx_pdb_read_conf() does for PDB. Only the first model is read. If
     * \p name is not nullptr, the extended PDB entry ID (\c _entry.id, falling
     * back to \c _struct.title / \c _database_2.database_code) is duplicated
     * into it with gmx_strdup(). The caller must have sized \p atoms and \p x
     * for the atom count returned by countAtoms().
     *
     * If \p conect is not nullptr, the bonds listed in \c struct_conn are added
     * to it, mirroring how read_pdbfile() collects CONECT records.
     */
    static void read(const std::filesystem::path& infile,
                     t_symtab*                    symtab,
                     char**                       name,
                     t_atoms*                     atoms,
                     rvec                         x[],
                     PbcType*                     pbcType,
                     matrix                       box,
                     gmx_conect_t*                conect = nullptr);

    //! Count the atoms (first model only) in an mmCIF/PDBx file.
    static void countAtoms(const std::filesystem::path& infile, int* natoms);

    /*! \brief Write all atoms of \p atoms to an mmCIF/PDBx file.
     *
     * \p conect may be nullptr; see writeIndexed() for what it adds.
     */
    static void write(const std::filesystem::path& outfile,
                      const char*                  title,
                      const t_atoms*               atoms,
                      const rvec                   x[],
                      PbcType                      pbcType,
                      const matrix                 box,
                      gmx_conect_t*                conect = nullptr);

    /*! \brief Write the \p nindex atoms selected by \p index to an mmCIF/PDBx file.
     *
     * When \p conect is not nullptr, its bonds are written as a \c struct_conn
     * loop (the mmCIF counterpart of the PDB LINK and SSBOND records). Bonds to
     * atoms outside \p index are skipped. Cis peptide bonds found in the written
     * coordinates are annotated in \c struct_mon_prot_cis independently of
     * \p conect.
     */
    static void writeIndexed(const std::filesystem::path& outfile,
                             const char*                  title,
                             const t_atoms*               atoms,
                             const rvec                   x[],
                             PbcType                      pbcType,
                             const matrix                 box,
                             int                          nindex,
                             const int                    index[],
                             gmx_conect_t*                conect = nullptr);

    /*! \name Trajectory (multi-model) I/O
     *
     * mmCIF files can hold several models (e.g. an NMR ensemble or a saved
     * trajectory) in one atom_site loop, distinguished by pdbx_PDB_model_num.
     * These entry points let the trajectory framework (trxio) treat each model
     * as a frame. Reading parses the whole file once (the bundled parser is not
     * line-streamable) and then iterates models; writing streams one model per
     * frame to an already-open FILE*.
     * \{
     */

    //! Parse \p infile and return an opaque multi-model handle (nullptr on empty).
    static CifTrajectory* openTrajectory(const std::filesystem::path& infile);
    //! Number of atoms per model in an opened trajectory handle.
    static int trajectoryNumAtoms(CifTrajectory* handle);
    /*! \brief Fill \p fr with the next model's coordinates/box/time.
     * \p fr->x must already be allocated for trajectoryNumAtoms() atoms.
     * Returns false when all models have been consumed. */
    static bool readTrajectoryFrame(CifTrajectory* handle, t_trxframe* fr);
    //! Release a handle from openTrajectory().
    static void closeTrajectory(CifTrajectory* handle);

    /*! \brief Stream one model (frame) of \p nindex atoms to \p out.
     * When \p writeHeader is true the data block, optional cell, the optional
     * \c struct_conn connectivity of \p conect and the atom_site loop header are
     * written first; \p modelNr is the 1-based pdbx_PDB_model_num for this frame.
     * Connectivity does not change between frames, so \p conect is only used
     * when \p writeHeader is true. */
    static void writeTrajectoryFrame(FILE*          out,
                                     const char*    title,
                                     const t_atoms* atoms,
                                     const rvec     x[],
                                     const matrix   box,
                                     int            nindex,
                                     const int      index[],
                                     int            modelNr,
                                     bool           writeHeader,
                                     gmx_conect_t*  conect = nullptr);
    //! \}
};

#endif // GMX_FILEIO_CIFIO_H
