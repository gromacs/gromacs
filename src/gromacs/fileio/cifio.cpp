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
#include <initializer_list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/units.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"

// RCSB mmCIF library (bundled under src/external/mmcif-dict-suite). Its headers
// are included only here; cifio.h merely declares the few library types that
// CifTrajectory names.
//
// Note that most of its accessors are not const-qualified: CifFile::GetBlock(),
// Block::GetTablePtr() and ISTable::IsColumnPresent() all need a non-const
// object even though none of them modifies the file.
// The reading code below therefore has to hold non-const pointers throughout,
// however read-only its intent.
#include "CifFile.h"
#include "CifFileUtil.h"
#include "ISTable.h"
#include "TableFile.h"

namespace
{

//! Values used by mmCIF to denote an unknown ('?') or inapplicable ('.') item.
bool isCifNull(std::string_view s)
{
    return s.empty() || s == "?" || s == ".";
}

/*! \brief The first of \p candidates that is a column of \p table, or an empty string if none is.
 *
 * Many items are optional, so the reading loops resolve each column name once
 * with this and then address cells with ISTable::operator() through valueAt().
 */
std::string findColName(ISTable* table, std::initializer_list<std::string> candidates)
{
    for (const std::string& name : candidates)
    {
        if (table->IsColumnPresent(name))
        {
            return name;
        }
    }
    return std::string();
}

//! Value at \p rowNum of \p colName as resolved by findColName(); "" when there is none.
const std::string& valueAt(const ISTable& table, unsigned int rowNum, const std::string& colName)
{
    static const std::string c_absent;
    return colName.empty() ? c_absent : table(rowNum, colName);
}

/*! \brief Row numbers of the atom_site table that belong to the first model.
 *
 * mmCIF files may hold several NMR/ensemble models in one atom_site loop; like
 * the PDB reader we keep only the first one.
 */
std::vector<unsigned int> firstModelRows(ISTable* atomSite)
{
    const unsigned int        numRows = atomSite->GetNumRows();
    std::vector<unsigned int> rowNums;
    rowNums.reserve(numRows);
    const std::string modelColName = "pdbx_PDB_model_num";
    if (numRows == 0 || !atomSite->IsColumnPresent(modelColName))
    {
        // Without a model column, every row belongs to the one model there is.
        for (unsigned int rowNum = 0; rowNum < numRows; rowNum++)
        {
            rowNums.push_back(rowNum);
        }
        return rowNums;
    }
    const std::string& firstModel = (*atomSite)(0, modelColName);
    for (unsigned int rowNum = 0; rowNum < numRows; rowNum++)
    {
        if ((*atomSite)(rowNum, modelColName) == firstModel)
        {
            rowNums.push_back(rowNum);
        }
    }
    return rowNums;
}

//! Parse \p infile, throwing when the parser reports an error.
std::unique_ptr<CifFile> parseCifFile(const std::filesystem::path& infile)
{
    std::unique_ptr<CifFile> cif(ParseCifSimple(infile.string()));
    // ParseCifSimple() always hands back a file object, even for input it could
    // not make sense of: what went wrong is reported in the parsing diagnostics.
    GMX_RELEASE_ASSERT(cif != nullptr, "ParseCifSimple() should always return a file");
    const std::string& diagnostics = cif->GetParsingDiags();
    const std::string  fileName    = infile.string();
    if (diagnostics.find("ERROR") != std::string::npos)
    {
        GMX_THROW(gmx::InvalidInputError("Could not parse mmCIF file " + fileName + "\n" + diagnostics));
    }
    if (!diagnostics.empty())
    {
        // Not fatal, but the file is not clean and that should not pass unnoticed.
        std::fprintf(stderr, "WARNING: while parsing %s:\n%s", fileName.c_str(), diagnostics.c_str());
    }
    return cif;
}

//! The first data block of \p cif, or nullptr when it has none.
Block* firstBlock(CifFile* cif)
{
    return (cif->GetNumBlocks() == 0) ? nullptr : &cif->GetBlock(cif->GetFirstBlockName());
}

//! The atom_site table of \p block, or nullptr when there is no block or no such table.
ISTable* atomSiteTable(Block* block)
{
    return (block == nullptr) ? nullptr : block->GetTablePtr("atom_site");
}

//! Convert unit-cell parameters (in Angstrom / degrees) to a GROMACS box (nm).
void unitCellToBox(double a, double b, double c, double alpha, double beta, double gamma, matrix box)
{
    // Same construction as read_cryst1() in pdbio.cpp, Angstrom -> nm.
    const double fa = a / gmx::c_nm2A;
    const double fb = b / gmx::c_nm2A;
    const double fc = c / gmx::c_nm2A;
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
        box[ZZ][ZZ] = std::sqrt((fc * fc) - (box[ZZ][XX] * box[ZZ][XX]) - (box[ZZ][YY] * box[ZZ][YY]));
    }
    else
    {
        box[YY][YY] = fb;
        box[ZZ][ZZ] = fc;
    }
}

/*! \brief Read the unit cell of \p block into \p box / \p pbcType.
 *
 * Returns whether a usable unit cell was found, and touches \p box and \p pbcType
 * only when it was. Callers must therefore have set them to their no-box
 * defaults beforehand: an mmCIF file legitimately carries no unit cell at all (an
 * NMR ensemble, for one), so that is reported by returning false rather than
 * treated as an error. A unit cell that is present but unusable does warn, since
 * unlike an absent one it means the file is malformed.
 */
bool readUnitCell(Block& block, PbcType* pbcType, matrix box)
{
    ISTable* unitCellTable = block.GetTablePtr("cell");
    if (unitCellTable == nullptr || unitCellTable->GetNumRows() == 0)
    {
        return false;
    }
    auto get = [unitCellTable](const char* colName) -> std::string {
        return unitCellTable->IsColumnPresent(colName) ? (*unitCellTable)(0, colName) : std::string();
    };

    const std::string lengths[DIM] = { get("length_a"), get("length_b"), get("length_c") };
    double            edge[DIM];
    for (int d = 0; d < DIM; d++)
    {
        if (isCifNull(lengths[d]))
        {
            // '?' and '.' say the edge is unknown or does not apply, which is
            // the file declining to give a cell rather than getting it wrong.
            return false;
        }
        // Otherwise a value that is not a positive number means the file is
        // broken, so reject it instead of letting strtod() turn it into a
        // zero-volume box that later code would take at face value.
        char* endOfNumber = nullptr;
        edge[d]           = std::strtod(lengths[d].c_str(), &endOfNumber);
        if (endOfNumber == lengths[d].c_str() || edge[d] <= 0.0)
        {
            std::fprintf(stderr,
                         "WARNING: mmCIF cell gives no usable length_%c (\"%s\"), so no box is "
                         "read from this file\n",
                         static_cast<char>('a' + d),
                         lengths[d].c_str());
            return false;
        }
    }

    auto angle = [&get](const char* name) -> double
    {
        const std::string v = get(name);
        return isCifNull(v) ? 90.0 : std::strtod(v.c_str(), nullptr);
    };
    if (box != nullptr)
    {
        unitCellToBox(
                edge[XX], edge[YY], edge[ZZ], angle("angle_alpha"), angle("angle_beta"), angle("angle_gamma"), box);
    }
    if (pbcType != nullptr)
    {
        *pbcType = PbcType::Xyz;
    }
    return true;
}

//! Extract the (possibly extended) PDB entry ID / title of \p block.
std::string readEntryId(Block& block)
{
    // _entry.id is the canonical (extended, up to 12 char) identifier.
    if (ISTable* entry = block.GetTablePtr("entry");
        entry != nullptr && entry->GetNumRows() > 0 && entry->IsColumnPresent("id"))
    {
        const std::string id = (*entry)(0, "id");
        if (!isCifNull(id))
        {
            return id;
        }
    }
    // Fall back to a human-readable title, then to the deposition database code.
    if (ISTable* st = block.GetTablePtr("struct");
        st != nullptr && st->GetNumRows() > 0 && st->IsColumnPresent("title"))
    {
        const std::string title = (*st)(0, "title");
        if (!isCifNull(title))
        {
            return title;
        }
    }
    if (ISTable* db = block.GetTablePtr("database_2");
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
    return gmx::formatString("%.3f", gmx::c_nm2A * xNm);
}

//! The chain identifier written for residue \p resind ('A' when it is unset).
char chainId(const t_atoms* atoms, int resind)
{
    const char chainCh = atoms->resinfo[resind].chainid;
    return (chainCh == 0 || chainCh == ' ') ? 'A' : chainCh;
}

/*! \brief Fill row \p rowNum of \p atomSite, which has the atomSiteColNames() columns, with atom \p i.
 *
 * Each value is put in its cell by column name, so this is where a GROMACS
 * field meets the mmCIF item that holds it.
 */
void fillAtomSiteRow(ISTable& atomSite, unsigned int rowNum, const t_atoms* atoms, int i, const rvec x[], int modelNr)
{
    const int         resind = atoms->atom[i].resind;
    const std::string resnm  = *atoms->resinfo[resind].name;
    const std::string atomnm = *atoms->atomname[i];
    const std::string seqStr = std::to_string(atoms->resinfo[resind].nr);
    const std::string chainStr(1, chainId(atoms, resind));

    const char*       elem   = atoms->atom[i].elem;
    const std::string symbol = (elem[0] != '\0') ? std::string(elem) : std::string("?");

    std::string groupPdb = "ATOM";
    std::string altLoc   = ".";
    real        occup    = 1.0;
    real        bfac     = 0.0;
    if (atoms->pdbinfo != nullptr)
    {
        groupPdb = (atoms->pdbinfo[i].type == PdbRecordType::Hetatm) ? "HETATM" : "ATOM";
        occup    = atoms->pdbinfo[i].occup;
        bfac     = atoms->pdbinfo[i].bfac;
        if (atoms->pdbinfo[i].altloc != 0 && atoms->pdbinfo[i].altloc != ' ')
        {
            altLoc = std::string(1, atoms->pdbinfo[i].altloc);
        }
    }

    atomSite.UpdateCell(rowNum, "group_PDB", groupPdb);
    atomSite.UpdateCell(rowNum, "id", std::to_string(i + 1));
    atomSite.UpdateCell(rowNum, "type_symbol", symbol);
    atomSite.UpdateCell(rowNum, "label_atom_id", atomnm);
    atomSite.UpdateCell(rowNum, "label_alt_id", altLoc);
    atomSite.UpdateCell(rowNum, "label_comp_id", resnm);
    atomSite.UpdateCell(rowNum, "label_asym_id", chainStr);
    // GROMACS has no notion of an mmCIF entity and writes no entity category,
    // so there is nothing to point at.
    atomSite.UpdateCell(rowNum, "label_entity_id", ".");
    atomSite.UpdateCell(rowNum, "label_seq_id", seqStr);
    atomSite.UpdateCell(rowNum, "Cartn_x", coordToStr(x[i][XX]));
    atomSite.UpdateCell(rowNum, "Cartn_y", coordToStr(x[i][YY]));
    atomSite.UpdateCell(rowNum, "Cartn_z", coordToStr(x[i][ZZ]));
    atomSite.UpdateCell(rowNum, "occupancy", gmx::formatString("%.2f", occup));
    atomSite.UpdateCell(rowNum, "B_iso_or_equiv", gmx::formatString("%.2f", bfac));
    atomSite.UpdateCell(rowNum, "auth_seq_id", seqStr);
    atomSite.UpdateCell(rowNum, "auth_comp_id", resnm);
    atomSite.UpdateCell(rowNum, "auth_asym_id", chainStr);
    atomSite.UpdateCell(rowNum, "pdbx_PDB_model_num", std::to_string(modelNr));
}

/*! \brief The atom_site column names, in the order they are written.
 *
 * label_alt_id and label_entity_id carry nothing GROMACS knows about, but readers
 * are entitled to expect them: gemmi, for one, reads no atoms at all from a file
 * that leaves them out.
 */
const std::vector<const char*>& atomSiteColNames()
{
    static const std::vector<const char*> colNames = { "group_PDB",     "id",
                                                       "type_symbol",   "label_atom_id",
                                                       "label_alt_id",  "label_comp_id",
                                                       "label_asym_id", "label_entity_id",
                                                       "label_seq_id",  "Cartn_x",
                                                       "Cartn_y",       "Cartn_z",
                                                       "occupancy",     "B_iso_or_equiv",
                                                       "auth_seq_id",   "auth_comp_id",
                                                       "auth_asym_id",  "pdbx_PDB_model_num" };
    return colNames;
}

//! Format \p value as a CIF token, quoting it when it holds a space.
std::string cifToken(const std::string& value)
{
    std::string tok = value.empty() ? std::string(".") : value;
    if (tok.find(' ') != std::string::npos)
    {
        tok.insert(tok.begin(), '\'');
        tok.push_back('\'');
    }
    return tok;
}

//! Add \p rows to \p block as a loop named \p category; does nothing when empty.
void addLoopTable(Block&                                        block,
                  const char*                                   category,
                  gmx::ArrayRef<const char* const>              colNames,
                  gmx::ArrayRef<const std::vector<std::string>> rows)
{
    if (rows.empty())
    {
        return;
    }
    ISTable& table = block.AddTable(category);
    for (const char* colName : colNames)
    {
        table.AddColumn(colName);
    }
    for (const std::vector<std::string>& row : rows)
    {
        table.AddRow(row);
    }
}

//! Stream \p rows to \p out as a loop named \p category; does nothing when empty.
void writeLoop(FILE*                                         out,
               const char*                                   category,
               gmx::ArrayRef<const char* const>              colNames,
               gmx::ArrayRef<const std::vector<std::string>> rows)
{
    if (rows.empty())
    {
        return;
    }
    std::fprintf(out, "loop_\n");
    for (const char* colName : colNames)
    {
        std::fprintf(out, "_%s.%s\n", category, colName);
    }
    for (const std::vector<std::string>& row : rows)
    {
        for (size_t colNum = 0; colNum < row.size(); colNum++)
        {
            std::fprintf(out, "%s%s", cifToken(row[colNum]).c_str(), (colNum + 1 < row.size()) ? " " : "\n");
        }
    }
    std::fprintf(out, "#\n");
}

/*! \name Connectivity and cis-peptide annotation
 *
 * The PDB CONECT record has no mmCIF counterpart at all (every field maps to
 * "n.a." in the wwPDB PDB/PDBx correspondence). Connectivity is instead expressed
 * by struct_conn, onto which the PDB LINK and SSBOND records map, so that is what
 * the bonds of a gmx_conect are written as. Writing every bond, rather than only
 * those not implied by the chemical component dictionary, goes beyond how the PDB
 * uses struct_conn; it is done deliberately so that -conect means the same thing
 * for mmCIF as for PDB output, where GROMACS likewise writes every bond.
 * \{
 */

//! The struct_conn column names, in the order produced by buildStructConnRows().
const std::vector<const char*>& structConnColNames()
{
    static const std::vector<const char*> colNames = { "id",
                                                       "conn_type_id",
                                                       "ptnr1_label_asym_id",
                                                       "ptnr1_label_comp_id",
                                                       "ptnr1_label_seq_id",
                                                       "ptnr1_label_atom_id",
                                                       "pdbx_ptnr1_label_alt_id",
                                                       "pdbx_ptnr1_PDB_ins_code",
                                                       "ptnr1_auth_asym_id",
                                                       "ptnr1_auth_comp_id",
                                                       "ptnr1_auth_seq_id",
                                                       "ptnr1_symmetry",
                                                       "ptnr2_label_asym_id",
                                                       "ptnr2_label_comp_id",
                                                       "ptnr2_label_seq_id",
                                                       "ptnr2_label_atom_id",
                                                       "pdbx_ptnr2_label_alt_id",
                                                       "pdbx_ptnr2_PDB_ins_code",
                                                       "ptnr2_auth_asym_id",
                                                       "ptnr2_auth_comp_id",
                                                       "ptnr2_auth_seq_id",
                                                       "ptnr2_symmetry",
                                                       "pdbx_dist_value" };
    return colNames;
}

//! The struct_conn_type column names.
const std::vector<const char*>& structConnTypeColNames()
{
    static const std::vector<const char*> colNames = { "id", "criteria", "reference" };
    return colNames;
}

/*! \brief Append the eleven fields describing bond partner \p i to \p row.
 *
 * The LINK/SSBOND correspondence points at the auth_* chain, residue name and
 * residue number and at the label_ atom name; the label_ chain, residue name and
 * residue number are written as well because deposited entries carry both and
 * readers differ in which they use. Both sets repeat what fillAtomSiteRow()
 * writes for the same atom, so the pointers resolve by construction.
 */
void appendConnPartner(std::vector<std::string>* row, const t_atoms* atoms, int i)
{
    const int         resind = atoms->atom[i].resind;
    const std::string comp   = *atoms->resinfo[resind].name;
    const std::string asym(1, chainId(atoms, resind));
    const std::string seq = std::to_string(atoms->resinfo[resind].nr);

    row->push_back(asym);
    row->push_back(comp);
    row->push_back(seq);
    row->push_back(*atoms->atomname[i]);
    // atom_site holds neither alternate locations nor insertion codes, so there
    // is nothing for these two to point at.
    row->emplace_back(".");
    row->emplace_back("?");
    row->push_back(asym);
    row->push_back(comp);
    row->push_back(seq);
    row->emplace_back("1_555");
}

//! Whether atom \p i is a cysteine-like SG sulfur, i.e. a disulfide candidate.
bool isThiolSulfur(const t_atoms* atoms, int i)
{
    if (std::strcmp(*atoms->atomname[i], "SG") != 0)
    {
        return false;
    }
    // The element is only known when the caller filled it in; when it is, it has
    // to agree, so that an atom merely named SG is not mistaken for a sulfur.
    const char* elem = atoms->atom[i].elem;
    return elem[0] == '\0' || std::strcmp(elem, "S") == 0;
}

/*! \brief Whether some residue holds two written atoms of the same name.
 *
 * struct_conn addresses atoms by name rather than by atom_site.id, so such a
 * name makes every connection to it ambiguous. Atoms of one residue are
 * contiguous, so one pass over the written atoms finds these.
 */
bool hasAmbiguousAtomNames(const t_atoms* atoms, const std::vector<bool>& written)
{
    std::vector<std::string> names;
    int                      currentResind = -1;
    for (int i = 0; i < atoms->nr; i++)
    {
        if (!written[i])
        {
            continue;
        }
        const int resind = atoms->atom[i].resind;
        if (resind != currentResind)
        {
            names.clear();
            currentResind = resind;
        }
        const std::string name = *atoms->atomname[i];
        if (std::find(names.begin(), names.end(), name) != names.end())
        {
            return true;
        }
        names.push_back(name);
    }
    return false;
}

/*! \brief The struct_conn rows for the bonds of \p conect.
 *
 * Bonds with an end outside the \p nindex atoms selected by \p index are left
 * out, so that the records only ever point at atoms that were written. The
 * connection types actually used are returned through \p connTypes, for the
 * accompanying struct_conn_type table.
 */
std::vector<std::vector<std::string>> buildStructConnRows(const t_atoms*            atoms,
                                                          const rvec                x[],
                                                          int                       nindex,
                                                          const int                 index[],
                                                          gmx_conect_t*             conect,
                                                          std::vector<std::string>* connTypes)
{
    std::vector<std::vector<std::string>> rows;
    const int                             numBonds = gmx_conect_nconect(conect);
    if (numBonds == 0)
    {
        return rows;
    }

    std::vector<bool> written(atoms->nr, false);
    for (int ii = 0; ii < nindex; ii++)
    {
        written[index[ii]] = true;
    }

    int numCovale = 0;
    int numDisulf = 0;
    for (int b = 0; b < numBonds; b++)
    {
        int ai = 0;
        int aj = 0;
        gmx_conect_get(conect, b, &ai, &aj);
        if (ai < 0 || ai >= atoms->nr || aj < 0 || aj >= atoms->nr || !written[ai] || !written[aj])
        {
            continue;
        }

        const bool disulfide = isThiolSulfur(atoms, ai) && isThiolSulfur(atoms, aj);

        std::vector<std::string> row;
        row.reserve(structConnColNames().size());
        row.push_back(disulfide ? "disulf" + std::to_string(++numDisulf)
                                : "covale" + std::to_string(++numCovale));
        row.emplace_back(disulfide ? "disulf" : "covale");
        appendConnPartner(&row, atoms, ai);
        appendConnPartner(&row, atoms, aj);
        row.push_back(gmx::formatString("%.3f", gmx::c_nm2A * std::sqrt(distance2(x[ai], x[aj]))));
        rows.push_back(std::move(row));
    }

    if (numCovale > 0)
    {
        connTypes->emplace_back("covale");
    }
    if (numDisulf > 0)
    {
        connTypes->emplace_back("disulf");
    }
    if (!rows.empty() && hasAmbiguousAtomNames(atoms, written))
    {
        std::fprintf(stderr,
                     "WARNING: some residues hold several atoms with the same name. mmCIF "
                     "struct_conn records identify atoms by name, so the connectivity written "
                     "for those residues is ambiguous.\n");
    }
    return rows;
}

//! How struct_conn addresses one atom: chain, residue number and atom name.
std::string connectionKey(std::string_view asym, std::string_view seq, std::string_view atom)
{
    std::string key;
    key.reserve(asym.size() + seq.size() + atom.size() + 2);
    key.append(asym).append(1, '|').append(seq).append(1, '|').append(atom);
    return key;
}

//! The names of the struct_conn columns that address one bond partner, each with a fallback.
struct PartnerColNames
{
    std::string asym, asymFallback;
    std::string seq, seqFallback;
    std::string atom, atomFallback;
};

//! Value of \p colName at \p rowNum of \p table, or of \p fallbackColName when that is null.
const std::string& preferredValue(const ISTable&     table,
                                  unsigned int       rowNum,
                                  const std::string& colName,
                                  const std::string& fallbackColName)
{
    const std::string& value = valueAt(table, rowNum, colName);
    return isCifNull(value) ? valueAt(table, rowNum, fallbackColName) : value;
}

/*! \brief Add the bonds listed in the struct_conn table of \p block to \p conect.
 *
 * The LINK and SSBOND correspondence points at the auth_ chain and residue number
 * and at the label_ atom name, but files in the wild do not all fill the same
 * subset, so the other set of items is used as a fallback. Records naming atoms
 * that were not read, such as those of a further model, are counted and reported.
 */
void readStructConn(Block& block, const t_atoms* atoms, int numAtoms, gmx_conect_t* conect)
{
    ISTable* structConn = block.GetTablePtr("struct_conn");
    if (structConn == nullptr)
    {
        return;
    }

    // Where the atoms that were read sit, keyed as struct_conn addresses them.
    std::map<std::string, int> atomIndex;
    for (int i = 0; i < numAtoms; i++)
    {
        const t_resinfo& resinfo = atoms->resinfo[atoms->atom[i].resind];
        atomIndex.emplace(connectionKey(std::string(1, resinfo.chainid),
                                        std::to_string(resinfo.nr),
                                        *atoms->atomname[i]),
                          i);
    }

    // Prepare the names of the columns so they can be
    // efficiently reused while reading every row.
    const PartnerColNames partner1ColNames = { findColName(structConn, { "ptnr1_auth_asym_id" }),
                                               findColName(structConn, { "ptnr1_label_asym_id" }),
                                               findColName(structConn, { "ptnr1_auth_seq_id" }),
                                               findColName(structConn, { "ptnr1_label_seq_id" }),
                                               findColName(structConn, { "ptnr1_label_atom_id" }),
                                               findColName(structConn, { "ptnr1_auth_atom_id" }) };
    const PartnerColNames partner2ColNames = { findColName(structConn, { "ptnr2_auth_asym_id" }),
                                               findColName(structConn, { "ptnr2_label_asym_id" }),
                                               findColName(structConn, { "ptnr2_auth_seq_id" }),
                                               findColName(structConn, { "ptnr2_label_seq_id" }),
                                               findColName(structConn, { "ptnr2_label_atom_id" }),
                                               findColName(structConn, { "ptnr2_auth_atom_id" }) };
    // The atom a struct_conn row names as one partner, or atomIndex.end() if it was not read.
    auto findPartner = [&atomIndex, structConn](unsigned int rowNum, const PartnerColNames& colNames)
    {
        return atomIndex.find(connectionKey(
                preferredValue(*structConn, rowNum, colNames.asym, colNames.asymFallback),
                preferredValue(*structConn, rowNum, colNames.seq, colNames.seqFallback),
                preferredValue(*structConn, rowNum, colNames.atom, colNames.atomFallback)));
    };

    int numUnresolved = 0;
    for (unsigned int rowNum = 0; rowNum < structConn->GetNumRows(); rowNum++)
    {
        const auto partner1 = findPartner(rowNum, partner1ColNames);
        if (partner1 == atomIndex.end())
        {
            numUnresolved++;
            continue;
        }
        const auto partner2 = findPartner(rowNum, partner2ColNames);
        if (partner2 == atomIndex.end())
        {
            numUnresolved++;
            continue;
        }
        gmx_conect_add(conect, partner1->second, partner2->second);
    }
    if (numUnresolved > 0)
    {
        std::fprintf(stderr,
                     "WARNING: %d struct_conn record(s) name atoms that were not read; they are "
                     "ignored\n",
                     numUnresolved);
    }
}
//! \}

} // namespace

void CifIO::read(const std::filesystem::path& infile,
                 t_symtab*                    symtab,
                 char**                       name,
                 t_atoms*                     atoms,
                 rvec                         x[],
                 PbcType*                     pbcType,
                 matrix                       box,
                 gmx_conect_t*                conect)
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

    const std::unique_ptr<CifFile> cif      = parseCifFile(infile);
    Block* const                   block    = firstBlock(cif.get());
    ISTable* const                 atomSite = atomSiteTable(block);
    if (atomSite == nullptr)
    {
        if (name != nullptr)
        {
            *name = gmx_strdup("");
        }
        return;
    }

    // Resolve the column of each item we read once; author (auth_*) items match
    // the PDB semantics GROMACS users expect, with the label_* items as fallback.
    // This resolution is done once before entering the loop over rows
    // so that the loop runs faster.
    const std::string groupColName    = findColName(atomSite, { "group_PDB" });
    const std::string serialColName   = findColName(atomSite, { "id" });
    const std::string elementColName  = findColName(atomSite, { "type_symbol" });
    const std::string atomNameColName = findColName(atomSite, { "auth_atom_id", "label_atom_id" });
    const std::string resNameColName  = findColName(atomSite, { "auth_comp_id", "label_comp_id" });
    const std::string chainColName    = findColName(atomSite, { "auth_asym_id", "label_asym_id" });
    const std::string seqColName      = findColName(atomSite, { "auth_seq_id", "label_seq_id" });
    const std::string insCodeColName  = findColName(atomSite, { "pdbx_PDB_ins_code" });
    const std::string altLocColName   = findColName(atomSite, { "label_alt_id" });
    const std::string xColName        = findColName(atomSite, { "Cartn_x" });
    const std::string yColName        = findColName(atomSite, { "Cartn_y" });
    const std::string zColName        = findColName(atomSite, { "Cartn_z" });
    const std::string occupColName    = findColName(atomSite, { "occupancy" });
    const std::string bfacColName     = findColName(atomSite, { "B_iso_or_equiv" });

    const ISTable&                  site    = *atomSite;
    const std::vector<unsigned int> rowNums = firstModelRows(atomSite);

    int              natom    = 0;
    int              chainnum = 0;
    std::string_view prevChain =
            rowNums.empty() ? std::string_view() : valueAt(site, rowNums[0], chainColName);
    // mmCIF permits chain identifiers of more than one character, while t_resinfo
    // stores a single char, so collect whatever has to be truncated and report it
    // once rather than per atom.
    std::set<std::string> truncatedChains;
    for (unsigned int rowNum : rowNums)
    {
        if (natom >= atoms->nr)
        {
            GMX_THROW(gmx::InconsistentInputError(
                    gmx::formatString("Found more atoms in mmCIF file %s than the %d expected",
                                      infile.string().c_str(),
                                      atoms->nr)));
        }

        const std::string_view chainStr     = valueAt(site, rowNum, chainColName);
        const bool             chainChanged = (natom > 0 && chainStr != prevChain);
        if (chainChanged)
        {
            chainnum++;
        }
        prevChain = chainStr;

        const char*            resnm   = valueAt(site, rowNum, resNameColName).c_str();
        const int              resnr   = std::atoi(valueAt(site, rowNum, seqColName).c_str());
        const std::string_view ins     = valueAt(site, rowNum, insCodeColName);
        const unsigned char    resic   = isCifNull(ins) ? ' ' : static_cast<unsigned char>(ins[0]);
        const char             chainid = chainStr.empty() ? ' ' : chainStr[0];
        if (chainStr.size() > 1)
        {
            truncatedChains.emplace(chainStr);
        }

        // Start a new residue on any change of chain, number, insertion code or name.
        // Chains numbered independently can put two residues with the same name and
        // number next to each other, e.g. a ligand 401 on chain A, then on chain B.
        if (natom == 0 || chainChanged || atoms->resinfo[atoms->atom[natom - 1].resind].nr != resnr
            || atoms->resinfo[atoms->atom[natom - 1].resind].ic != resic
            || std::strcmp(*atoms->resinfo[atoms->atom[natom - 1].resind].name, resnm) != 0)
        {
            atoms->atom[natom].resind = (natom == 0) ? 0 : atoms->atom[natom - 1].resind + 1;
            atoms->nres               = atoms->atom[natom].resind + 1;
            t_atoms_set_resinfo(atoms, natom, symtab, resnm, resnr, resic, chainnum, chainid);
        }
        else
        {
            atoms->atom[natom].resind = atoms->atom[natom - 1].resind;
        }

        const std::string& atomName   = valueAt(site, rowNum, atomNameColName);
        const char*        elem       = valueAt(site, rowNum, elementColName).c_str();
        atoms->atomname[natom]        = put_symtab(symtab, atomName.c_str());
        atoms->atom[natom].m          = 0.0;
        atoms->atom[natom].q          = 0.0;
        atoms->atom[natom].atomnumber = 0;
        std::strncpy(atoms->atom[natom].elem, elem, sizeof(atoms->atom[natom].elem) - 1);
        atoms->atom[natom].elem[sizeof(atoms->atom[natom].elem) - 1] = '\0';

        x[natom][XX] = std::strtod(valueAt(site, rowNum, xColName).c_str(), nullptr) / gmx::c_nm2A;
        x[natom][YY] = std::strtod(valueAt(site, rowNum, yColName).c_str(), nullptr) / gmx::c_nm2A;
        x[natom][ZZ] = std::strtod(valueAt(site, rowNum, zColName).c_str(), nullptr) / gmx::c_nm2A;

        if (atoms->pdbinfo != nullptr)
        {
            const std::string_view group = valueAt(site, rowNum, groupColName);
            const std::string_view alt   = valueAt(site, rowNum, altLocColName);
            atoms->pdbinfo[natom].type = (group == "HETATM") ? PdbRecordType::Hetatm : PdbRecordType::Atom;
            atoms->pdbinfo[natom].atomnr = std::atoi(valueAt(site, rowNum, serialColName).c_str());
            atoms->pdbinfo[natom].altloc = isCifNull(alt) ? ' ' : alt[0];
            std::strncpy(atoms->pdbinfo[natom].atomnm,
                         atomName.c_str(),
                         sizeof(atoms->pdbinfo[natom].atomnm) - 1);
            atoms->pdbinfo[natom].atomnm[sizeof(atoms->pdbinfo[natom].atomnm) - 1] = '\0';
            atoms->pdbinfo[natom].occup =
                    std::strtod(valueAt(site, rowNum, occupColName).c_str(), nullptr);
            atoms->pdbinfo[natom].bfac = std::strtod(valueAt(site, rowNum, bfacColName).c_str(), nullptr);
            atoms->pdbinfo[natom].bAnisotropic = FALSE;
        }
        natom++;
    }

    if (!truncatedChains.empty())
    {
        std::string list;
        for (const std::string& id : truncatedChains)
        {
            list += (list.empty() ? "" : ", ") + id;
        }
        std::fprintf(stderr,
                     "WARNING: mmCIF chain identifiers longer than one character (%s) were "
                     "truncated to their first character, because that is all GROMACS stores. "
                     "Chains that differ only beyond it can no longer be told apart by "
                     "identifier.\n",
                     list.c_str());
    }

    readUnitCell(*block, pbcType, box);

    if (conect != nullptr)
    {
        readStructConn(*block, atoms, natom, conect);
    }

    if (name != nullptr)
    {
        *name = gmx_strdup(readEntryId(*block).c_str());
    }
}

void CifIO::countAtoms(const std::filesystem::path& infile, int* natoms)
{
    const std::unique_ptr<CifFile> cif      = parseCifFile(infile);
    ISTable* const                 atomSite = atomSiteTable(firstBlock(cif.get()));
    *natoms = (atomSite == nullptr) ? 0 : firstModelRows(atomSite).size();
}

void CifIO::writeIndexed(const std::filesystem::path& outfile,
                         const char*                  title,
                         const t_atoms*               atoms,
                         const rvec                   x[],
                         PbcType                      pbcType,
                         const matrix                 box,
                         int                          nindex,
                         const int                    index[],
                         gmx_conect_t*                conect)
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

        ISTable& unitCellTable = block.AddTable("cell");
        for (const char* colName : { "entry_id",
                                     "length_a",
                                     "length_b",
                                     "length_c",
                                     "angle_alpha",
                                     "angle_beta",
                                     "angle_gamma" })
        {
            unitCellTable.AddColumn(colName);
        }
        unitCellTable.AddRow(std::vector<std::string>{ entryId,
                                                       gmx::formatString("%.3f", gmx::c_nm2A * a),
                                                       gmx::formatString("%.3f", gmx::c_nm2A * b),
                                                       gmx::formatString("%.3f", gmx::c_nm2A * c),
                                                       gmx::formatString("%.2f", alpha),
                                                       gmx::formatString("%.2f", beta),
                                                       gmx::formatString("%.2f", gamma) });

        ISTable& sym = block.AddTable("symmetry");
        sym.AddColumn("entry_id");
        sym.AddColumn("space_group_name_H-M");
        sym.AddColumn("Int_Tables_number");
        sym.AddRow(std::vector<std::string>{ entryId, "P 1", "1" });
    }

    // Connectivity and cis peptide annotation precede the coordinates, as they
    // do in deposited entries.
    std::vector<std::string> connTypes;
    addLoopTable(block,
                 "struct_conn",
                 structConnColNames(),
                 buildStructConnRows(atoms, x, nindex, index, conect, &connTypes));
    std::vector<std::vector<std::string>> connTypeRows;
    connTypeRows.reserve(connTypes.size());
    for (const std::string& type : connTypes)
    {
        connTypeRows.push_back(std::vector<std::string>{ type, ".", "." });
    }
    addLoopTable(block, "struct_conn_type", structConnTypeColNames(), connTypeRows);

    // atom_site loop.
    ISTable& atomSite = block.AddTable("atom_site");
    for (const char* colName : atomSiteColNames())
    {
        atomSite.AddColumn(colName);
    }
    for (int ii = 0; ii < nindex; ii++)
    {
        atomSite.AddRow();
        fillAtomSiteRow(atomSite, atomSite.GetNumRows() - 1, atoms, index[ii], x, /*modelNr*/ 1);
    }

    cif.Write(outfile.string());
}

void CifIO::write(const std::filesystem::path& outfile,
                  const char*                  title,
                  const t_atoms*               atoms,
                  const rvec                   x[],
                  PbcType                      pbcType,
                  const matrix                 box,
                  gmx_conect_t*                conect)
{
    std::vector<int> index(atoms->nr);
    for (int i = 0; i < atoms->nr; i++)
    {
        index[i] = i;
    }
    writeIndexed(outfile, title, atoms, x, pbcType, box, atoms->nr, index.data(), conect);
}

CifTrajectory::CifTrajectory(std::unique_ptr<CifFile> cif, Block& block, ISTable* atomSite) :
    cif_(std::move(cif)),
    atomSite_(atomSite),
    xColName_(findColName(atomSite, { "Cartn_x" })),
    yColName_(findColName(atomSite, { "Cartn_y" })),
    zColName_(findColName(atomSite, { "Cartn_z" }))
{
    // Without a model column every row reads as model "", so all go to one model.
    const std::string        modelColName = findColName(atomSite_, { "pdbx_PDB_model_num" });
    const unsigned int       numRows      = atomSite_->GetNumRows();
    std::vector<std::string> order; // model ids in first-seen order
    for (unsigned int rowNum = 0; rowNum < numRows; rowNum++)
    {
        const std::string& id = valueAt(*atomSite_, rowNum, modelColName);
        auto               it = std::find(order.begin(), order.end(), id);
        size_t             mi;
        if (it == order.end())
        {
            mi = order.size();
            order.push_back(id);
            modelRowNums_.emplace_back();
        }
        else
        {
            mi = static_cast<size_t>(std::distance(order.begin(), it));
        }
        modelRowNums_[mi].push_back(rowNum);
    }
    numAtoms_ = modelRowNums_.empty() ? 0 : static_cast<int>(modelRowNums_[0].size());
    haveBox_  = readUnitCell(block, nullptr, box_);
}

CifTrajectory::~CifTrajectory() = default;

int CifTrajectory::numAtoms() const
{
    return numAtoms_;
}

bool CifTrajectory::readNextFrame(t_trxframe* fr)
{
    if (cursor_ >= modelRowNums_.size())
    {
        return false;
    }
    const std::vector<unsigned int>& rowNums = modelRowNums_[cursor_];
    const ISTable&                   site    = *atomSite_;
    fr->natoms                               = numAtoms_;
    const int n = std::min(numAtoms_, static_cast<int>(rowNums.size()));
    for (int i = 0; i < n; i++)
    {
        const unsigned int rowNum = rowNums[i];
        fr->x[i][XX] = std::strtod(valueAt(site, rowNum, xColName_).c_str(), nullptr) / gmx::c_nm2A;
        fr->x[i][YY] = std::strtod(valueAt(site, rowNum, yColName_).c_str(), nullptr) / gmx::c_nm2A;
        fr->x[i][ZZ] = std::strtod(valueAt(site, rowNum, zColName_).c_str(), nullptr) / gmx::c_nm2A;
    }
    fr->bX = TRUE;
    if (haveBox_)
    {
        copy_mat(box_, fr->box);
    }
    fr->bBox  = haveBox_;
    fr->step  = static_cast<int64_t>(cursor_);
    fr->bStep = TRUE;
    fr->time  = static_cast<real>(cursor_);
    fr->bTime = TRUE;
    fr->bPrec = TRUE;
    fr->prec  = 10000;
    cursor_++;
    return true;
}

std::unique_ptr<CifTrajectory> CifIO::openTrajectory(const std::filesystem::path& infile)
{
    std::unique_ptr<CifFile> cif      = parseCifFile(infile);
    Block* const             block    = firstBlock(cif.get());
    ISTable* const           atomSite = atomSiteTable(block);
    if (atomSite == nullptr)
    {
        return nullptr;
    }
    return std::make_unique<CifTrajectory>(std::move(cif), *block, atomSite);
}

void CifIO::writeTrajectoryFrame(FILE*          out,
                                 const char*    title,
                                 const t_atoms* atoms,
                                 const rvec     x[],
                                 const matrix   box,
                                 int            nindex,
                                 const int      index[],
                                 int            modelNr,
                                 bool           writeHeader,
                                 gmx_conect_t*  conect)
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
        std::fprintf(out, "_entry.id   %s\n#\n", cifToken(entryId).c_str());

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
                         gmx::c_nm2A * a,
                         gmx::c_nm2A * b,
                         gmx::c_nm2A * c);
            std::fprintf(out,
                         "_cell.angle_alpha   %.2f\n_cell.angle_beta   %.2f\n_cell.angle_gamma   "
                         "%.2f\n#\n",
                         alpha,
                         beta,
                         gamma);
        }
        // Connectivity does not change from frame to frame, so it is written once,
        // before the atom_site loop that the remaining frames are appended to.
        std::vector<std::string> connTypes;
        writeLoop(out,
                  "struct_conn",
                  structConnColNames(),
                  buildStructConnRows(atoms, x, nindex, index, conect, &connTypes));
        std::vector<std::vector<std::string>> connTypeRows;
        connTypeRows.reserve(connTypes.size());
        for (const std::string& type : connTypes)
        {
            connTypeRows.push_back(std::vector<std::string>{ type, ".", "." });
        }
        writeLoop(out, "struct_conn_type", structConnTypeColNames(), connTypeRows);

        std::fprintf(out, "loop_\n");
        for (const char* colName : atomSiteColNames())
        {
            std::fprintf(out, "_atom_site.%s\n", colName);
        }
    }
    // Each atom is filled into one reusable row by the same code as the table
    // writer uses, and then streamed out in column order.
    const std::vector<const char*>& colNames = atomSiteColNames();
    ISTable                         atomSiteRow;
    for (const char* colName : colNames)
    {
        atomSiteRow.AddColumn(colName);
    }
    atomSiteRow.AddRow();
    for (int ii = 0; ii < nindex; ii++)
    {
        fillAtomSiteRow(atomSiteRow, 0, atoms, index[ii], x, modelNr);
        for (size_t colNum = 0; colNum < colNames.size(); colNum++)
        {
            std::fprintf(out,
                         "%s%s",
                         cifToken(atomSiteRow(0, colNames[colNum])).c_str(),
                         (colNum + 1 < colNames.size()) ? " " : "\n");
        }
    }
}
