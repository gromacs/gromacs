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
#include "gmxpre.h"

#include "cifio.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

// RCSB mmCIF library (bundled under src/external/mmcif-dict-suite-*). Kept
// strictly inside this translation unit so that cifio.h stays free of any
// third-party types.
#include "CifFile.h"
#include "CifFileUtil.h"
#include "ISTable.h"
#include "TableFile.h"

namespace
{

//! Values used by mmCIF to denote an unknown ('?') or inapplicable ('.') item.
bool isCifNull(const std::string& s)
{
    return s.empty() || s == "?" || s == ".";
}

//! Fetch a whole atom_site column by name; returns an empty vector if absent.
std::vector<std::string> getColumn(ISTable* table, const std::string& name)
{
    std::vector<std::string> col;
    if (table != nullptr && table->IsColumnPresent(name))
    {
        table->GetColumn(col, name);
    }
    return col;
}

//! Value of column \p col at row \p row, or empty string if the column is absent.
const std::string& cell(const std::vector<std::string>& col, size_t row)
{
    static const std::string empty;
    return (row < col.size()) ? col[row] : empty;
}

/*! \brief Row indices of the atom_site table that belong to the first model.
 *
 * mmCIF files may hold several NMR/ensemble models in one atom_site loop; like
 * the PDB reader we keep only the first one. If there is no model column, all
 * rows are returned.
 */
std::vector<unsigned int> firstModelRows(ISTable* atomSite)
{
    const unsigned int        nRows = atomSite->GetNumRows();
    std::vector<unsigned int> rows;
    std::vector<std::string>  model = getColumn(atomSite, "pdbx_PDB_model_num");
    const std::string         first = model.empty() ? std::string() : model[0];
    for (unsigned int r = 0; r < nRows; r++)
    {
        if (model.empty() || model[r] == first)
        {
            rows.push_back(r);
        }
    }
    return rows;
}

//! Parse \p infile and return the atom_site table of its first data block (or nullptr).
ISTable* parseFirstBlock(const std::filesystem::path& infile,
                         std::unique_ptr<CifFile>&    cifOwner,
                         std::string*                 blockName)
{
    CifFile* raw = ParseCifSimple(infile.string());
    if (raw == nullptr)
    {
        gmx_fatal(FARGS, "Could not parse mmCIF file %s", infile.string().c_str());
    }
    cifOwner.reset(raw);
    if (raw->GetNumBlocks() == 0)
    {
        return nullptr;
    }
    const std::string block = raw->GetFirstBlockName();
    if (blockName != nullptr)
    {
        *blockName = block;
    }
    return raw->GetBlock(block).GetTablePtr("atom_site");
}

//! Convert unit-cell parameters (in Angstrom / degrees) to a GROMACS box (nm).
void cellToBox(double a, double b, double c, double alpha, double beta, double gamma, matrix box)
{
    // Same construction as read_cryst1() in pdbio.cpp, Angstrom -> nm.
    const double fa = a * 0.1;
    const double fb = b * 0.1;
    const double fc = c * 0.1;
    clear_mat(box);
    box[XX][XX] = fa;
    if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0))
    {
        const double cosa = (alpha != 90.0) ? std::cos(alpha * gmx::c_deg2Rad) : 0.0;
        const double cosb = (beta != 90.0) ? std::cos(beta * gmx::c_deg2Rad) : 0.0;
        const double cosg = (gamma != 90.0) ? std::cos(gamma * gmx::c_deg2Rad) : 0.0;
        const double sing = (gamma != 90.0) ? std::sin(gamma * gmx::c_deg2Rad) : 1.0;
        box[YY][XX]       = fb * cosg;
        box[YY][YY]       = fb * sing;
        box[ZZ][XX]       = fc * cosb;
        box[ZZ][YY]       = fc * (cosa - cosb * cosg) / sing;
        box[ZZ][ZZ] = std::sqrt(fc * fc - box[ZZ][XX] * box[ZZ][XX] - box[ZZ][YY] * box[ZZ][YY]);
    }
    else
    {
        box[YY][YY] = fb;
        box[ZZ][ZZ] = fc;
    }
}

//! Read the unit cell of \p block into \p box / \p pbcType (no-op if absent).
void readCell(CifFile* cif, const std::string& block, PbcType* pbcType, matrix box)
{
    ISTable* cellTable = cif->GetBlock(block).GetTablePtr("cell");
    if (cellTable == nullptr || cellTable->GetNumRows() == 0)
    {
        return;
    }
    auto get = [cellTable](const char* name) -> std::string
    { return cellTable->IsColumnPresent(name) ? (*cellTable)(0, name) : std::string(); };
    const std::string a = get("length_a");
    const std::string b = get("length_b");
    const std::string c = get("length_c");
    if (isCifNull(a) || isCifNull(b) || isCifNull(c))
    {
        return;
    }
    auto angle = [&get](const char* name) -> double
    {
        const std::string v = get(name);
        return isCifNull(v) ? 90.0 : std::strtod(v.c_str(), nullptr);
    };
    if (box != nullptr)
    {
        cellToBox(std::strtod(a.c_str(), nullptr),
                  std::strtod(b.c_str(), nullptr),
                  std::strtod(c.c_str(), nullptr),
                  angle("angle_alpha"),
                  angle("angle_beta"),
                  angle("angle_gamma"),
                  box);
    }
    if (pbcType != nullptr)
    {
        *pbcType = PbcType::Xyz;
    }
}

//! Extract the (possibly extended) PDB entry ID / title of \p block.
std::string readEntryId(CifFile* cif, const std::string& block)
{
    Block& b = cif->GetBlock(block);
    // _entry.id is the canonical (extended, up to 12 char) identifier.
    if (ISTable* entry = b.GetTablePtr("entry");
        entry != nullptr && entry->GetNumRows() > 0 && entry->IsColumnPresent("id"))
    {
        const std::string id = (*entry)(0, "id");
        if (!isCifNull(id))
        {
            return id;
        }
    }
    // Fall back to a human-readable title, then to the deposition database code.
    if (ISTable* st = b.GetTablePtr("struct");
        st != nullptr && st->GetNumRows() > 0 && st->IsColumnPresent("title"))
    {
        const std::string title = (*st)(0, "title");
        if (!isCifNull(title))
        {
            return title;
        }
    }
    if (ISTable* db = b.GetTablePtr("database_2");
        db != nullptr && db->GetNumRows() > 0 && db->IsColumnPresent("database_code"))
    {
        const std::string code = (*db)(0, "database_code");
        if (!isCifNull(code))
        {
            return code;
        }
    }
    return std::string();
}

//! Format a coordinate (nm) as an Angstrom mmCIF value.
std::string coordToStr(real xNm)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3f", 10.0 * xNm);
    return buf;
}

//! Format a real with the given precision.
std::string realToStr(real value, const char* fmt)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), fmt, value);
    return buf;
}

//! Build the 16 atom_site field values for atom \p i (column order as written below).
std::vector<std::string> buildAtomSiteRow(const t_atoms* atoms, int i, const rvec x[], int modelNr)
{
    const int         resind = atoms->atom[i].resind;
    const std::string resnm  = *atoms->resinfo[resind].name;
    const std::string atomnm = *atoms->atomname[i];
    const std::string seqStr = std::to_string(atoms->resinfo[resind].nr);

    char chainCh = atoms->resinfo[resind].chainid;
    if (chainCh == 0 || chainCh == ' ')
    {
        chainCh = 'A';
    }
    const std::string chainStr(1, chainCh);

    const char*       elem   = atoms->atom[i].elem;
    const std::string symbol = (elem[0] != '\0') ? std::string(elem) : std::string("?");

    std::string groupPdb = "ATOM";
    real        occup    = 1.0;
    real        bfac     = 0.0;
    if (atoms->pdbinfo != nullptr)
    {
        groupPdb = (atoms->pdbinfo[i].type == PdbRecordType::Hetatm) ? "HETATM" : "ATOM";
        occup    = atoms->pdbinfo[i].occup;
        bfac     = atoms->pdbinfo[i].bfac;
    }

    return std::vector<std::string>{ groupPdb,
                                     std::to_string(i + 1),
                                     symbol,
                                     atomnm,
                                     resnm,
                                     chainStr,
                                     seqStr,
                                     coordToStr(x[i][XX]),
                                     coordToStr(x[i][YY]),
                                     coordToStr(x[i][ZZ]),
                                     realToStr(occup, "%.2f"),
                                     realToStr(bfac, "%.2f"),
                                     seqStr,
                                     resnm,
                                     chainStr,
                                     std::to_string(modelNr) };
}

//! The atom_site column names, in the order produced by buildAtomSiteRow().
const std::vector<const char*>& atomSiteColumns()
{
    static const std::vector<const char*> columns = { "group_PDB",     "id",
                                                      "type_symbol",   "label_atom_id",
                                                      "label_comp_id", "label_asym_id",
                                                      "label_seq_id",  "Cartn_x",
                                                      "Cartn_y",       "Cartn_z",
                                                      "occupancy",     "B_iso_or_equiv",
                                                      "auth_seq_id",   "auth_comp_id",
                                                      "auth_asym_id",  "pdbx_PDB_model_num" };
    return columns;
}

} // namespace

void CifIO::read(const std::filesystem::path& infile,
                 t_symtab*                    symtab,
                 char**                       name,
                 t_atoms*                     atoms,
                 rvec                         x[],
                 PbcType*                     pbcType,
                 matrix                       box)
{
    if (pbcType != nullptr)
    {
        *pbcType = PbcType::No;
    }
    if (box != nullptr)
    {
        clear_mat(box);
    }

    atoms->haveMass    = FALSE;
    atoms->haveCharge  = FALSE;
    atoms->haveType    = FALSE;
    atoms->haveBState  = FALSE;
    atoms->havePdbInfo = (atoms->pdbinfo != nullptr);

    std::unique_ptr<CifFile> cif;
    std::string              block;
    ISTable*                 atomSite = parseFirstBlock(infile, cif, &block);
    if (atomSite == nullptr)
    {
        if (name != nullptr)
        {
            *name = gmx_strdup("");
        }
        return;
    }

    // Read the columns we care about once each; author (auth_*) items match the
    // PDB semantics GROMACS users expect, with the label_* items as fallback.
    const std::vector<std::string> group    = getColumn(atomSite, "group_PDB");
    const std::vector<std::string> serial   = getColumn(atomSite, "id");
    const std::vector<std::string> element  = getColumn(atomSite, "type_symbol");
    std::vector<std::string>       atomName = getColumn(atomSite, "auth_atom_id");
    if (atomName.empty())
    {
        atomName = getColumn(atomSite, "label_atom_id");
    }
    std::vector<std::string> resName = getColumn(atomSite, "auth_comp_id");
    if (resName.empty())
    {
        resName = getColumn(atomSite, "label_comp_id");
    }
    std::vector<std::string> chain = getColumn(atomSite, "auth_asym_id");
    if (chain.empty())
    {
        chain = getColumn(atomSite, "label_asym_id");
    }
    std::vector<std::string> seq = getColumn(atomSite, "auth_seq_id");
    if (seq.empty())
    {
        seq = getColumn(atomSite, "label_seq_id");
    }
    const std::vector<std::string> insCode = getColumn(atomSite, "pdbx_PDB_ins_code");
    const std::vector<std::string> altLoc  = getColumn(atomSite, "label_alt_id");
    const std::vector<std::string> cartnX  = getColumn(atomSite, "Cartn_x");
    const std::vector<std::string> cartnY  = getColumn(atomSite, "Cartn_y");
    const std::vector<std::string> cartnZ  = getColumn(atomSite, "Cartn_z");
    const std::vector<std::string> occCol  = getColumn(atomSite, "occupancy");
    const std::vector<std::string> bfacCol = getColumn(atomSite, "B_iso_or_equiv");

    const std::vector<unsigned int> rows = firstModelRows(atomSite);

    int         natom     = 0;
    int         chainnum  = 0;
    std::string prevChain = cell(chain, rows.empty() ? 0 : rows[0]);
    for (unsigned int r : rows)
    {
        if (natom >= atoms->nr)
        {
            gmx_fatal(FARGS,
                      "\nFound more atoms in mmCIF file %s than expected (%d)",
                      infile.string().c_str(),
                      atoms->nr);
        }

        const std::string& chainStr = cell(chain, r);
        if (natom > 0 && chainStr != prevChain)
        {
            chainnum++;
        }
        prevChain = chainStr;

        const std::string&  resnm   = cell(resName, r);
        const int           resnr   = std::atoi(cell(seq, r).c_str());
        const std::string&  ins     = cell(insCode, r);
        const unsigned char resic   = isCifNull(ins) ? ' ' : static_cast<unsigned char>(ins[0]);
        const char          chainid = chainStr.empty() ? ' ' : chainStr[0];

        // Start a new residue on any change of number, insertion code or name.
        if (natom == 0 || atoms->resinfo[atoms->atom[natom - 1].resind].nr != resnr
            || atoms->resinfo[atoms->atom[natom - 1].resind].ic != resic
            || std::strcmp(*atoms->resinfo[atoms->atom[natom - 1].resind].name, resnm.c_str()) != 0)
        {
            atoms->atom[natom].resind = (natom == 0) ? 0 : atoms->atom[natom - 1].resind + 1;
            atoms->nres               = atoms->atom[natom].resind + 1;
            t_atoms_set_resinfo(atoms, natom, symtab, resnm.c_str(), resnr, resic, chainnum, chainid);
        }
        else
        {
            atoms->atom[natom].resind = atoms->atom[natom - 1].resind;
        }

        const std::string& elem       = cell(element, r);
        atoms->atomname[natom]        = put_symtab(symtab, cell(atomName, r).c_str());
        atoms->atom[natom].m          = 0.0;
        atoms->atom[natom].q          = 0.0;
        atoms->atom[natom].atomnumber = 0;
        std::strncpy(atoms->atom[natom].elem, elem.c_str(), sizeof(atoms->atom[natom].elem) - 1);
        atoms->atom[natom].elem[sizeof(atoms->atom[natom].elem) - 1] = '\0';

        x[natom][XX] = std::strtod(cell(cartnX, r).c_str(), nullptr) * 0.1;
        x[natom][YY] = std::strtod(cell(cartnY, r).c_str(), nullptr) * 0.1;
        x[natom][ZZ] = std::strtod(cell(cartnZ, r).c_str(), nullptr) * 0.1;

        if (atoms->pdbinfo != nullptr)
        {
            const std::string& group_r = cell(group, r);
            const std::string& alt     = cell(altLoc, r);
            atoms->pdbinfo[natom].type =
                    (group_r == "HETATM") ? PdbRecordType::Hetatm : PdbRecordType::Atom;
            atoms->pdbinfo[natom].atomnr = std::atoi(cell(serial, r).c_str());
            atoms->pdbinfo[natom].altloc = isCifNull(alt) ? ' ' : alt[0];
            std::strncpy(atoms->pdbinfo[natom].atomnm,
                         cell(atomName, r).c_str(),
                         sizeof(atoms->pdbinfo[natom].atomnm) - 1);
            atoms->pdbinfo[natom].atomnm[sizeof(atoms->pdbinfo[natom].atomnm) - 1] = '\0';
            atoms->pdbinfo[natom].occup        = std::strtod(cell(occCol, r).c_str(), nullptr);
            atoms->pdbinfo[natom].bfac         = std::strtod(cell(bfacCol, r).c_str(), nullptr);
            atoms->pdbinfo[natom].bAnisotropic = FALSE;
        }
        natom++;
    }

    readCell(cif.get(), block, pbcType, box);

    if (name != nullptr)
    {
        *name = gmx_strdup(readEntryId(cif.get(), block).c_str());
    }
}

void CifIO::countAtoms(const std::filesystem::path& infile, int* natoms)
{
    std::unique_ptr<CifFile> cif;
    ISTable*                 atomSite = parseFirstBlock(infile, cif, nullptr);
    *natoms                           = (atomSite == nullptr) ? 0 : firstModelRows(atomSite).size();
}

void CifIO::writeIndexed(const std::filesystem::path& outfile,
                         const char*                  title,
                         const t_atoms*               atoms,
                         const rvec                   x[],
                         PbcType                      pbcType,
                         const matrix                 box,
                         int                          nindex,
                         const int                    index[])
{
    const std::string entryId = (title != nullptr && title[0] != '\0') ? title : "gromacs";

    CifFile cif(/*verbose*/ false);
    // mmCIF data-block names cannot contain whitespace.
    std::string blockName;
    for (char ch : entryId)
    {
        blockName.push_back(std::isspace(static_cast<unsigned char>(ch)) ? '_' : ch);
    }
    const std::string actualBlock = cif.AddBlock(blockName);
    Block&            block       = cif.GetBlock(actualBlock);

    // Tables are created with AddTable() so the block owns them; the returned
    // references stay valid until cif.Write() below (a local ISTable would be
    // destroyed too early and leave the CifFile with dangling table data).

    // _entry.id carries the (possibly extended) identifier; _struct.title the text.
    ISTable& entry = block.AddTable("entry");
    entry.AddColumn("id");
    entry.AddRow(std::vector<std::string>{ entryId });

    ISTable& st = block.AddTable("struct");
    st.AddColumn("entry_id");
    st.AddColumn("title");
    st.AddRow(std::vector<std::string>{ entryId, entryId });

    // Unit cell + symmetry (CRYST1 equivalent), only when a box is present.
    const bool haveBox = box != nullptr && pbcType != PbcType::No
                         && (norm2(box[XX]) != 0.0F || norm2(box[YY]) != 0.0F || norm2(box[ZZ]) != 0.0F);
    if (haveBox)
    {
        const real a     = norm(box[XX]);
        const real b     = norm(box[YY]);
        const real c     = norm(box[ZZ]);
        const real alpha = (b * c != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[YY], box[ZZ]) : 90.0F;
        const real beta  = (a * c != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[XX], box[ZZ]) : 90.0F;
        const real gamma = (a * b != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[XX], box[YY]) : 90.0F;

        ISTable& cell = block.AddTable("cell");
        for (const char* col : { "entry_id",
                                 "length_a",
                                 "length_b",
                                 "length_c",
                                 "angle_alpha",
                                 "angle_beta",
                                 "angle_gamma" })
        {
            cell.AddColumn(col);
        }
        cell.AddRow(std::vector<std::string>{ entryId,
                                              realToStr(10.0 * a, "%.3f"),
                                              realToStr(10.0 * b, "%.3f"),
                                              realToStr(10.0 * c, "%.3f"),
                                              realToStr(alpha, "%.2f"),
                                              realToStr(beta, "%.2f"),
                                              realToStr(gamma, "%.2f") });

        ISTable& sym = block.AddTable("symmetry");
        sym.AddColumn("entry_id");
        sym.AddColumn("space_group_name_H-M");
        sym.AddColumn("Int_Tables_number");
        sym.AddRow(std::vector<std::string>{ entryId, "P 1", "1" });
    }

    // atom_site loop.
    ISTable& atomSite = block.AddTable("atom_site");
    for (const char* col : atomSiteColumns())
    {
        atomSite.AddColumn(col);
    }
    for (int ii = 0; ii < nindex; ii++)
    {
        atomSite.AddRow(buildAtomSiteRow(atoms, index[ii], x, /*modelNr*/ 1));
    }

    cif.Write(outfile.string());
}

void CifIO::write(const std::filesystem::path& outfile,
                  const char*                  title,
                  const t_atoms*               atoms,
                  const rvec                   x[],
                  PbcType                      pbcType,
                  const matrix                 box)
{
    std::vector<int> index(atoms->nr);
    for (int i = 0; i < atoms->nr; i++)
    {
        index[i] = i;
    }
    writeIndexed(outfile, title, atoms, x, pbcType, box, atoms->nr, index.data());
}

//! Internal state behind CifIO's opaque CifTrajectory handle (see cifio.h).
struct CifTrajectory
{
    std::unique_ptr<CifFile>               cif;
    std::vector<std::string>               cartnX;
    std::vector<std::string>               cartnY;
    std::vector<std::string>               cartnZ;
    std::vector<std::vector<unsigned int>> modelRows; //!< row indices per model, in file order
    int                                    natoms  = 0;
    bool                                   haveBox = false;
    matrix                                 box     = { { 0 } };
    size_t                                 cursor  = 0; //!< next model to hand out
};

CifTrajectory* CifIO::openTrajectory(const std::filesystem::path& infile)
{
    auto        traj = std::make_unique<CifTrajectory>();
    std::string block;
    ISTable*    atomSite = parseFirstBlock(infile, traj->cif, &block);
    if (atomSite == nullptr)
    {
        return nullptr;
    }
    traj->cartnX = getColumn(atomSite, "Cartn_x");
    traj->cartnY = getColumn(atomSite, "Cartn_y");
    traj->cartnZ = getColumn(atomSite, "Cartn_z");

    const std::vector<std::string> model = getColumn(atomSite, "pdbx_PDB_model_num");
    const unsigned int             nRows = atomSite->GetNumRows();
    std::vector<std::string>       order; // model ids in first-seen order
    for (unsigned int r = 0; r < nRows; r++)
    {
        const std::string id = model.empty() ? std::string("1") : model[r];
        auto              it = std::find(order.begin(), order.end(), id);
        size_t            mi;
        if (it == order.end())
        {
            mi = order.size();
            order.push_back(id);
            traj->modelRows.emplace_back();
        }
        else
        {
            mi = static_cast<size_t>(std::distance(order.begin(), it));
        }
        traj->modelRows[mi].push_back(r);
    }
    traj->natoms = traj->modelRows.empty() ? 0 : static_cast<int>(traj->modelRows[0].size());

    PbcType pt = PbcType::No;
    readCell(traj->cif.get(), block, &pt, traj->box);
    traj->haveBox = (pt != PbcType::No);

    return traj.release();
}

int CifIO::trajectoryNumAtoms(CifTrajectory* traj)
{
    return (traj != nullptr) ? traj->natoms : 0;
}

bool CifIO::readTrajectoryFrame(CifTrajectory* traj, t_trxframe* fr)
{
    if (traj == nullptr || traj->cursor >= traj->modelRows.size())
    {
        return false;
    }
    const std::vector<unsigned int>& rows = traj->modelRows[traj->cursor];
    fr->natoms                            = traj->natoms;
    const int n                           = std::min(traj->natoms, static_cast<int>(rows.size()));
    for (int i = 0; i < n; i++)
    {
        const unsigned int r = rows[i];
        fr->x[i][XX]         = std::strtod(cell(traj->cartnX, r).c_str(), nullptr) * 0.1;
        fr->x[i][YY]         = std::strtod(cell(traj->cartnY, r).c_str(), nullptr) * 0.1;
        fr->x[i][ZZ]         = std::strtod(cell(traj->cartnZ, r).c_str(), nullptr) * 0.1;
    }
    fr->bX = TRUE;
    if (traj->haveBox)
    {
        copy_mat(traj->box, fr->box);
    }
    fr->bBox  = traj->haveBox;
    fr->step  = static_cast<int64_t>(traj->cursor);
    fr->bStep = TRUE;
    fr->time  = static_cast<real>(traj->cursor);
    fr->bTime = TRUE;
    fr->bPrec = TRUE;
    fr->prec  = 10000;
    traj->cursor++;
    return true;
}

void CifIO::closeTrajectory(CifTrajectory* traj)
{
    delete traj;
}

void CifIO::writeTrajectoryFrame(FILE*          out,
                                 const char*    title,
                                 const t_atoms* atoms,
                                 const rvec     x[],
                                 const matrix   box,
                                 int            nindex,
                                 const int      index[],
                                 int            modelNr,
                                 bool           writeHeader)
{
    if (writeHeader)
    {
        const std::string entryId = (title != nullptr && title[0] != '\0') ? title : "gromacs";
        std::string       blockName;
        for (char ch : entryId)
        {
            blockName.push_back(std::isspace(static_cast<unsigned char>(ch)) ? '_' : ch);
        }
        std::fprintf(out, "data_%s\n#\n", blockName.c_str());
        std::fprintf(out, "_entry.id   %s\n#\n", entryId.c_str());

        const bool haveBox =
                box != nullptr
                && (norm2(box[XX]) != 0.0F || norm2(box[YY]) != 0.0F || norm2(box[ZZ]) != 0.0F);
        if (haveBox)
        {
            const real a = norm(box[XX]);
            const real b = norm(box[YY]);
            const real c = norm(box[ZZ]);
            const real alpha = (b * c != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[YY], box[ZZ]) : 90.0F;
            const real beta = (a * c != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[XX], box[ZZ]) : 90.0F;
            const real gamma = (a * b != 0.0F) ? gmx::c_rad2Deg * gmx_angle(box[XX], box[YY]) : 90.0F;
            std::fprintf(out,
                         "_cell.length_a   %.3f\n_cell.length_b   %.3f\n_cell.length_c   %.3f\n",
                         10.0 * a,
                         10.0 * b,
                         10.0 * c);
            std::fprintf(out,
                         "_cell.angle_alpha   %.2f\n_cell.angle_beta   %.2f\n_cell.angle_gamma   "
                         "%.2f\n#\n",
                         alpha,
                         beta,
                         gamma);
        }
        std::fprintf(out, "loop_\n");
        for (const char* col : atomSiteColumns())
        {
            std::fprintf(out, "_atom_site.%s\n", col);
        }
    }
    for (int ii = 0; ii < nindex; ii++)
    {
        const std::vector<std::string> row = buildAtomSiteRow(atoms, index[ii], x, modelNr);
        for (size_t k = 0; k < row.size(); k++)
        {
            std::string tok = row[k].empty() ? std::string(".") : row[k];
            if (tok.find(' ') != std::string::npos)
            {
                tok = "'" + tok + "'";
            }
            std::fprintf(out, "%s%s", tok.c_str(), (k + 1 < row.size()) ? " " : "\n");
        }
    }
}
