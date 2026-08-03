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
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/fileio/pdbio.h"
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

//! The chain identifier written for residue \p resind ('A' when it is unset).
char chainId(const t_atoms* atoms, int resind)
{
    const char chainCh = atoms->resinfo[resind].chainid;
    return (chainCh == 0 || chainCh == ' ') ? 'A' : chainCh;
}

//! Build the 18 atom_site field values for atom \p i (column order as written below).
std::vector<std::string> buildAtomSiteRow(const t_atoms* atoms, int i, const rvec x[], int modelNr)
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

    return std::vector<std::string>{ groupPdb,
                                     std::to_string(i + 1),
                                     symbol,
                                     atomnm,
                                     altLoc,
                                     resnm,
                                     chainStr,
                                     // GROMACS has no notion of an mmCIF entity and writes
                                     // no entity category, so there is nothing to point at.
                                     ".",
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

/*! \brief The atom_site column names, in the order produced by buildAtomSiteRow().
 *
 * label_alt_id and label_entity_id carry nothing GROMACS knows about, but readers
 * are entitled to expect them: gemmi, for one, reads no atoms at all from a file
 * that leaves them out.
 */
const std::vector<const char*>& atomSiteColumns()
{
    static const std::vector<const char*> columns = { "group_PDB",     "id",
                                                      "type_symbol",   "label_atom_id",
                                                      "label_alt_id",  "label_comp_id",
                                                      "label_asym_id", "label_entity_id",
                                                      "label_seq_id",  "Cartn_x",
                                                      "Cartn_y",       "Cartn_z",
                                                      "occupancy",     "B_iso_or_equiv",
                                                      "auth_seq_id",   "auth_comp_id",
                                                      "auth_asym_id",  "pdbx_PDB_model_num" };
    return columns;
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
void addLoopTable(Block&                                       block,
                  const char*                                  category,
                  const std::vector<const char*>&              columns,
                  const std::vector<std::vector<std::string>>& rows)
{
    if (rows.empty())
    {
        return;
    }
    ISTable& table = block.AddTable(category);
    for (const char* col : columns)
    {
        table.AddColumn(col);
    }
    for (const std::vector<std::string>& row : rows)
    {
        table.AddRow(row);
    }
}

//! Stream \p rows to \p out as a loop named \p category; does nothing when empty.
void writeLoop(FILE*                                        out,
               const char*                                  category,
               const std::vector<const char*>&              columns,
               const std::vector<std::vector<std::string>>& rows)
{
    if (rows.empty())
    {
        return;
    }
    std::fprintf(out, "loop_\n");
    for (const char* col : columns)
    {
        std::fprintf(out, "_%s.%s\n", category, col);
    }
    for (const std::vector<std::string>& row : rows)
    {
        for (size_t k = 0; k < row.size(); k++)
        {
            std::fprintf(out, "%s%s", cifToken(row[k]).c_str(), (k + 1 < row.size()) ? " " : "\n");
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
 *
 * The related CISPEP record maps onto struct_mon_prot_cis, which is written from
 * the coordinates alone and therefore needs no connectivity.
 * \{
 */

//! The struct_conn column names, in the order produced by buildStructConnRows().
const std::vector<const char*>& structConnColumns()
{
    static const std::vector<const char*> columns = { "id",
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
    return columns;
}

//! The struct_conn_type column names.
const std::vector<const char*>& structConnTypeColumns()
{
    static const std::vector<const char*> columns = { "id", "criteria", "reference" };
    return columns;
}

/*! \brief Append the eleven fields describing bond partner \p i to \p row.
 *
 * The LINK/SSBOND correspondence points at the auth_* chain, residue name and
 * residue number and at the label_ atom name; the label_ chain, residue name and
 * residue number are written as well because deposited entries carry both and
 * readers differ in which they use. Both sets repeat what buildAtomSiteRow()
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
        row.reserve(structConnColumns().size());
        row.push_back(disulfide ? "disulf" + std::to_string(++numDisulf)
                                : "covale" + std::to_string(++numCovale));
        row.emplace_back(disulfide ? "disulf" : "covale");
        appendConnPartner(&row, atoms, ai);
        appendConnPartner(&row, atoms, aj);
        row.push_back(realToStr(10.0 * std::sqrt(distance2(x[ai], x[aj])), "%.3f"));
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

//! The struct_mon_prot_cis column names, in the order produced by buildCisPeptideRows().
const std::vector<const char*>& structMonProtCisColumns()
{
    static const std::vector<const char*> columns = { "pdbx_id",
                                                      "label_comp_id",
                                                      "label_asym_id",
                                                      "label_seq_id",
                                                      "label_alt_id",
                                                      "pdbx_PDB_ins_code",
                                                      "auth_comp_id",
                                                      "auth_asym_id",
                                                      "auth_seq_id",
                                                      "pdbx_label_comp_id_2",
                                                      "pdbx_label_asym_id_2",
                                                      "pdbx_label_seq_id_2",
                                                      "pdbx_PDB_ins_code_2",
                                                      "pdbx_auth_comp_id_2",
                                                      "pdbx_auth_asym_id_2",
                                                      "pdbx_auth_seq_id_2",
                                                      "pdbx_PDB_model_num",
                                                      "pdbx_omega_angle" };
    return columns;
}

/*! \brief The dihedral angle p1-p2-p3-p4 in degrees, in (-180, 180].
 *
 * Computed here rather than with dih_angle(), which lives in the listed_forces
 * module that fileio does not depend on.
 */
real dihedralDegrees(const rvec p1, const rvec p2, const rvec p3, const rvec p4)
{
    rvec b1, b2, b3, n1, n2, m1, unitB2;
    rvec_sub(p2, p1, b1);
    rvec_sub(p3, p2, b2);
    rvec_sub(p4, p3, b3);
    cprod(b1, b2, n1);
    cprod(b2, b3, n2);
    svmul(1.0 / norm(b2), b2, unitB2);
    cprod(n1, unitB2, m1);
    return gmx::c_rad2Deg * std::atan2(iprod(m1, n2), iprod(n1, n2));
}

/*! \brief The struct_mon_prot_cis rows for the peptide bonds in cis conformation.
 *
 * A peptide bond is recognised by residues r and r+1 of one chain holding the
 * backbone atoms CA and C, resp. N and CA, at a plausible C-N bond distance;
 * requiring those names is what distinguishes a peptide from anything else, so
 * no residue type database is needed. Only atoms among the \p nindex atoms
 * selected by \p index are considered, so a bond cannot be annotated across a
 * gap in the selection.
 */
std::vector<std::vector<std::string>>
buildCisPeptideRows(const t_atoms* atoms, const rvec x[], int nindex, const int index[])
{
    //! Largest |omega| (degrees) still called cis, as used by the wwPDB.
    constexpr real c_maxCisOmega = 30.0;
    //! Largest plausible peptide C-N bond length (nm); longer means a chain break.
    constexpr real c_maxPeptideBondLength = 0.2;

    std::vector<std::vector<std::string>> rows;
    if (atoms->nres < 2)
    {
        return rows;
    }

    // Backbone atoms of each residue, restricted to the written selection.
    std::vector<int> backboneN(atoms->nres, -1);
    std::vector<int> backboneCa(atoms->nres, -1);
    std::vector<int> backboneC(atoms->nres, -1);
    for (int ii = 0; ii < nindex; ii++)
    {
        const int   i      = index[ii];
        const int   resind = atoms->atom[i].resind;
        const char* name   = *atoms->atomname[i];
        if (std::strcmp(name, "N") == 0 && backboneN[resind] < 0)
        {
            backboneN[resind] = i;
        }
        else if (std::strcmp(name, "CA") == 0 && backboneCa[resind] < 0)
        {
            backboneCa[resind] = i;
        }
        else if (std::strcmp(name, "C") == 0 && backboneC[resind] < 0)
        {
            backboneC[resind] = i;
        }
    }

    for (int r = 0; r + 1 < atoms->nres; r++)
    {
        const int caPrev = backboneCa[r];
        const int cPrev  = backboneC[r];
        const int nNext  = backboneN[r + 1];
        const int caNext = backboneCa[r + 1];
        if (caPrev < 0 || cPrev < 0 || nNext < 0 || caNext < 0
            || chainId(atoms, r) != chainId(atoms, r + 1)
            || distance2(x[cPrev], x[nNext]) > c_maxPeptideBondLength * c_maxPeptideBondLength)
        {
            continue;
        }
        const real omega = dihedralDegrees(x[caPrev], x[cPrev], x[nNext], x[caNext]);
        if (std::abs(omega) >= c_maxCisOmega)
        {
            continue;
        }

        const std::string comp1 = *atoms->resinfo[r].name;
        const std::string asym1(1, chainId(atoms, r));
        const std::string seq1  = std::to_string(atoms->resinfo[r].nr);
        const std::string comp2 = *atoms->resinfo[r + 1].name;
        const std::string asym2(1, chainId(atoms, r + 1));
        const std::string seq2 = std::to_string(atoms->resinfo[r + 1].nr);

        rows.push_back(std::vector<std::string>{ std::to_string(rows.size() + 1),
                                                 comp1,
                                                 asym1,
                                                 seq1,
                                                 ".",
                                                 "?",
                                                 comp1,
                                                 asym1,
                                                 seq1,
                                                 comp2,
                                                 asym2,
                                                 seq2,
                                                 "?",
                                                 comp2,
                                                 asym2,
                                                 seq2,
                                                 "1",
                                                 realToStr(omega, "%.2f") });
    }
    return rows;
}

//! How struct_conn addresses one atom: chain, residue number and atom name.
std::string connectionKey(const std::string& asym, const std::string& seq, const std::string& atom)
{
    return asym + '|' + seq + '|' + atom;
}

//! Value of \p primary at row \p row, falling back to \p fallback when it is null.
const std::string& preferredCell(const std::vector<std::string>& primary,
                                 const std::vector<std::string>& fallback,
                                 size_t                          row)
{
    const std::string& value = cell(primary, row);
    return isCifNull(value) ? cell(fallback, row) : value;
}

/*! \brief Add the bonds listed in the struct_conn table of \p block to \p conect.
 *
 * The LINK and SSBOND correspondence points at the auth_ chain and residue number
 * and at the label_ atom name, but files in the wild do not all fill the same
 * subset, so the other set of items is used as a fallback. Records naming atoms
 * that were not read, such as those of a further model, are counted and reported.
 */
void readStructConn(CifFile* cif, const std::string& block, const t_atoms* atoms, int numAtoms, gmx_conect_t* conect)
{
    ISTable* structConn = cif->GetBlock(block).GetTablePtr("struct_conn");
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

    const std::vector<std::string> authAsym1  = getColumn(structConn, "ptnr1_auth_asym_id");
    const std::vector<std::string> labelAsym1 = getColumn(structConn, "ptnr1_label_asym_id");
    const std::vector<std::string> authSeq1   = getColumn(structConn, "ptnr1_auth_seq_id");
    const std::vector<std::string> labelSeq1  = getColumn(structConn, "ptnr1_label_seq_id");
    const std::vector<std::string> labelAtom1 = getColumn(structConn, "ptnr1_label_atom_id");
    const std::vector<std::string> authAtom1  = getColumn(structConn, "ptnr1_auth_atom_id");
    const std::vector<std::string> authAsym2  = getColumn(structConn, "ptnr2_auth_asym_id");
    const std::vector<std::string> labelAsym2 = getColumn(structConn, "ptnr2_label_asym_id");
    const std::vector<std::string> authSeq2   = getColumn(structConn, "ptnr2_auth_seq_id");
    const std::vector<std::string> labelSeq2  = getColumn(structConn, "ptnr2_label_seq_id");
    const std::vector<std::string> labelAtom2 = getColumn(structConn, "ptnr2_label_atom_id");
    const std::vector<std::string> authAtom2  = getColumn(structConn, "ptnr2_auth_atom_id");

    int numUnresolved = 0;
    for (unsigned int r = 0; r < structConn->GetNumRows(); r++)
    {
        const auto first  = atomIndex.find(connectionKey(preferredCell(authAsym1, labelAsym1, r),
                                                        preferredCell(authSeq1, labelSeq1, r),
                                                        preferredCell(labelAtom1, authAtom1, r)));
        const auto second = atomIndex.find(connectionKey(preferredCell(authAsym2, labelAsym2, r),
                                                         preferredCell(authSeq2, labelSeq2, r),
                                                         preferredCell(labelAtom2, authAtom2, r)));
        if (first == atomIndex.end() || second == atomIndex.end())
        {
            numUnresolved++;
            continue;
        }
        gmx_conect_add(conect, first->second, second->second);
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

    if (conect != nullptr)
    {
        readStructConn(cif.get(), block, atoms, natom, conect);
    }

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

    // Connectivity and cis peptide annotation precede the coordinates, as they
    // do in deposited entries.
    std::vector<std::string> connTypes;
    addLoopTable(block,
                 "struct_conn",
                 structConnColumns(),
                 buildStructConnRows(atoms, x, nindex, index, conect, &connTypes));
    std::vector<std::vector<std::string>> connTypeRows;
    connTypeRows.reserve(connTypes.size());
    for (const std::string& type : connTypes)
    {
        connTypeRows.push_back(std::vector<std::string>{ type, ".", "." });
    }
    addLoopTable(block, "struct_conn_type", structConnTypeColumns(), connTypeRows);
    addLoopTable(block,
                 "struct_mon_prot_cis",
                 structMonProtCisColumns(),
                 buildCisPeptideRows(atoms, x, nindex, index));

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
        // Connectivity does not change from frame to frame, so it is written once,
        // before the atom_site loop that the remaining frames are appended to.
        // The cis peptide annotation of struct_mon_prot_cis is per model and is
        // deliberately not written here: the atom_site loop stays open until the
        // last frame, so no further loop can follow it.
        std::vector<std::string> connTypes;
        writeLoop(out,
                  "struct_conn",
                  structConnColumns(),
                  buildStructConnRows(atoms, x, nindex, index, conect, &connTypes));
        std::vector<std::vector<std::string>> connTypeRows;
        connTypeRows.reserve(connTypes.size());
        for (const std::string& type : connTypes)
        {
            connTypeRows.push_back(std::vector<std::string>{ type, ".", "." });
        }
        writeLoop(out, "struct_conn_type", structConnTypeColumns(), connTypeRows);

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
            std::fprintf(out, "%s%s", cifToken(row[k]).c_str(), (k + 1 < row.size()) ? " " : "\n");
        }
    }
}
