/*
 * FDASettings.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <set>
#include <sstream>
#include "FDASettings.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"

using namespace fda;

const std::string FDASettings::compat_fm_version = "1.5";

const int FDASettings::compat_new_entry = -280480;

FDASettings::FDASettings(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop, bool parallel_execution)
 : atom_based_result_type(ResultType::NO),
   residue_based_result_type(ResultType::NO),
   one_pair(OnePair::DETAILED),
   v2s(Vector2Scalar::NORM),
   residues_renumber(ResiduesRenumber::AUTO),
   no_end_zeros(false),
   syslen_atoms(mtop->natoms),
   syslen_residues(0),
   time_averaging_period(1),
   sys_in_group1(syslen_atoms, 0),
   sys_in_group2(syslen_atoms, 0),
   type(InteractionType_NONE),
   nonbonded_exclusion_on(true),
   bonded_exclusion_on(true),
   index_group1(-1),
   index_group2(-1),
   groups(nullptr),
   groupnames(nullptr)
{
    /// Parallel execution not implemented yet
    if (parallel_execution)
        gmx_fatal(FARGS, "FDA parallel execution not implemented yet! Please start with '-nt 1' on the mdrun command line\n");

    // check for the pf configuration file (specified with -pfi option);
    // if it doesn't exist, return NULL to specify that no pf handling is done;
    // otherwise, check also for specification of the index file (-pfn)
    if (opt2bSet("-pfi", nfile, fnm)) {
        if (!opt2bSet("-pfn", nfile, fnm))
            gmx_fatal(FARGS, "-pfi option (pairwise forces configuration) specified, an index file (-pfn) is also needed.\n");
    } else {
        gmx_fatal(FARGS, "No pairwise forces input file, no sense to compute pairwise forces.\n");
    }

    // Use GROMACS function to read lines of form key = value from input file
    const char *pf_file_in = opt2fn("-pfi", nfile, fnm);
    warninp_t wi = init_warning(FALSE, 0);
    int ninp;
    gmx::TextInputFile stream(pf_file_in);
    t_inpfile *inp = read_inpfile(&stream, pf_file_in, &ninp, wi);

    // Check for deprecated keywords
    std::vector<std::pair<std::string, std::string>> deprecated_keywords{
        {"pf_onepair", "onepair"},
        {"pf_group1", "group1"},
        {"pf_group2", "group2"},
        {"pf_type", "type"},
        {"pf_atombased", "atombased"},
        {"pf_residuebased", "residuebased"},
        {"pf_vector2scalar", "vector2scalar"},
        {"pf_residuesrenumber", "residuesrenumber"},
        {"pf_time_averages_period", "time_averages_period"},
        {"pf_no_end_zeros", "no_end_zeros"}
    };
    for (auto const& pair : deprecated_keywords) {
        if (search_einp(ninp, inp, pair.first.c_str()) != -1) {
                std::string message = "Deprecated keyword '" + pair.first + "' was used, please use '" + pair.second + "' instead.\n";
            gmx_fatal(FARGS, message.c_str());
        }
    }

    // Read result types
    std::stringstream(get_estr(&ninp, &inp, "atombased", "no")) >> atom_based_result_type;
    std::stringstream(get_estr(&ninp, &inp, "residuebased", "no")) >> residue_based_result_type;

    // OnePair has to be initialized before the atoms/residues are initialized
    // because the data structures used for storing atoms/residues depend on it
    std::stringstream(get_estr(&ninp, &inp, "onepair", "detailed")) >> one_pair;

    std::stringstream(get_estr(&ninp, &inp, "residuesrenumber", "auto")) >> residues_renumber;
    std::cout << "ResidueRenumber: " << residues_renumber << std::endl;

    std::stringstream(get_estr(&ninp, &inp, "vector2scalar", "norm")) >> v2s;
    std::cout << "Vector2Scalar: " << v2s << std::endl;

    no_end_zeros = strcasecmp(get_estr(&ninp, &inp, "no_end_zeros", "no"), "no");

    if ((compatibility_mode(atom_based_result_type) or compatibility_mode(residue_based_result_type)) and v2s != Vector2Scalar::NORM)
        gmx_fatal(FARGS, "When using compat mode, pf_vector2scalar should be set to norm.\n");

    std::string type_string;
    std::stringstream(get_estr(&ninp, &inp, "type", "all")) >> type_string;
    type = from_string(type_string);
    std::cout << "Pairwise interactions selected: " << type_string << std::endl;

    if (type == InteractionType_NONE)
        gmx_fatal(FARGS, "No interactions selected, no sense to compute pairwise forces.\n");

    // Read group names
    std::string name_group1, name_group2;
    std::stringstream(get_estr(&ninp, &inp, "group1", "Protein")) >> name_group1;
    std::stringstream(get_estr(&ninp, &inp, "group2", "Protein")) >> name_group2;
    std::cout << "Pairwise forces for groups: " << name_group1 << " and " << name_group2 << std::endl;

    // Get atom indices of groups defined in pfn-file
    groups = init_index(opt2fn("-pfn", nfile, fnm), &groupnames);

    // Get group number
    for (int i = 0; i < groups->nr; ++i) {
        if (gmx_strcasecmp(groupnames[i], name_group1.c_str()) == 0) index_group1 = i;
        if (gmx_strcasecmp(groupnames[i], name_group2.c_str()) == 0) index_group2 = i;
    }

    // Check if using compatibility mode, there should be only one group
    if (compatibility_mode(atom_based_result_type) or compatibility_mode(residue_based_result_type)) {
        if (name_group1 != name_group2)
            gmx_fatal(FARGS, "When using compat mode, the two group names should the identical.\n");
        else
            groupname = name_group1;
    }

    // Read group information
    char** group_names;
    t_blocka *groups = init_index(opt2fn("-pfn", nfile, fnm), &group_names);
    if(groups->nr == 0) gmx_fatal(FARGS, "No groups found in the indexfile.\n");

    // Map atoms to residues
    fill_atom2residue(mtop);

    // Set sys_in_group arrays
    for (int i = 0; i != groups->nr; ++i) {
        std::vector<int> group_atoms;
        if (group_names[i] == name_group1) {
            group_atoms = std::vector<int>(groups->a + groups->index[i], groups->a + groups->index[i + 1]);
            for (auto g : group_atoms) sys_in_group1[g] = 1;
        }
        if (group_names[i] == name_group2) {
            group_atoms = std::vector<int>(groups->a + groups->index[i], groups->a + groups->index[i + 1]);
            for (auto g : group_atoms) sys_in_group2[g] = 1;
        }
        if (PF_or_PS_mode(atom_based_result_type)) {
            for (auto g : group_atoms) sys2pf_atoms[g] = sys2pf_atoms.size();
        }
        if (PF_or_PS_mode(residue_based_result_type)) {
            std::vector<int> group_residues = groupatoms2residues(group_atoms);
            for (auto g : group_residues) sys2pf_residues[g] = sys2pf_residues.size();
        }
    }

    // Read time averaging period
    time_averaging_period = get_eint(&ninp, &inp, "time_averages_period", 1, wi);
    if (time_averaging_period < 0)
        gmx_fatal(FARGS, "Invalid value for time_averages_period: %d\n", time_averaging_period);

    // Check for valid input options using time averaging
    if (time_averaging_period != 1) {
        if (one_pair != OnePair::SUMMED)
            gmx_fatal(FARGS, "Can only save scalar time averages from summed interactions.\n");
        if (PF_or_PS_mode(atom_based_result_type)) {
            if (!(compatibility_mode(atom_based_result_type) or atom_based_result_type == ResultType::PAIRWISE_FORCES_SCALAR))
                gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
        }
        if (PF_or_PS_mode(residue_based_result_type)) {
            if (!(compatibility_mode(residue_based_result_type) or residue_based_result_type == ResultType::PAIRWISE_FORCES_SCALAR))
                gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
        }
    }

    // Check if groups are defined for PF/PF mode
    if (PF_or_PS_mode(atom_based_result_type) or PF_or_PS_mode(residue_based_result_type)) {
        if (sys_in_group1.empty() or sys_in_group2.empty())
            gmx_fatal(FARGS, "No atoms in one or both groups.\n");
    }

    // Check that there is an index file
    if (!opt2bSet("-pfn", nfile, fnm))
        gmx_fatal(FARGS, "No index file (-pfn) for pairwise forces.\n");

    // Get output file names
    if (PF_or_PS_mode(atom_based_result_type)) {
        if (atom_based_result_type == ResultType::PUNCTUAL_STRESS) {
            atom_based_result_filename = gmx_strdup(opt2fn("-psa", nfile, fnm));
        } else {
            atom_based_result_filename = gmx_strdup(opt2fn("-pfa", nfile, fnm));
        }
        if (atom_based_result_filename.empty())
            gmx_fatal(FARGS, "No file for writing out atom-based values.\n");
    }

    if (PF_or_PS_mode(residue_based_result_type)) {
        if (residue_based_result_type == ResultType::PUNCTUAL_STRESS) {
            residue_based_result_filename = gmx_strdup(opt2fn("-psr", nfile, fnm));
        } else {
            residue_based_result_filename = gmx_strdup(opt2fn("-pfr", nfile, fnm));
        }
        if (residue_based_result_filename.empty())
            gmx_fatal(FARGS, "No file for writing out residue-based values.\n");
    }

    if (VS_mode(atom_based_result_type)) {
        if (atom_based_result_type == ResultType::VIRIAL_STRESS)
            atom_based_result_filename = gmx_strdup(opt2fn("-vsa", nfile, fnm));
        if (atom_based_result_type == ResultType::VIRIAL_STRESS_VON_MISES)
            atom_based_result_filename = gmx_strdup(opt2fn("-vma", nfile, fnm));
        if (atom_based_result_filename.empty())
            gmx_fatal(FARGS, "No file for writing out virial stress.\n");
    }

    if ((stress_mode(atom_based_result_type) or stress_mode(residue_based_result_type)) and (one_pair != OnePair::SUMMED))
        gmx_fatal(FARGS, "Per atom data can only be computed from summed interactions.\n");

    // Energy groups exclusions
    nonbonded_exclusion_on = strcasecmp(get_estr(&ninp, &inp, "energy_grp_exclusion", "yes"), "no");
    bonded_exclusion_on = strcasecmp(get_estr(&ninp, &inp, "bonded_exclusion", "yes"), "no");
}

std::vector<int> FDASettings::groupatoms2residues(std::vector<int> const& group_atoms) const
{
    std::set<int> group_residues;
    for (auto atom : group_atoms) group_residues.insert(atom_2_residue[atom]);
    return std::vector<int>(group_residues.begin(), group_residues.end());
}

void FDASettings::fill_atom2residue(gmx_mtop_t *mtop)
{
    t_atoms *atoms;
    t_atom *atom_info;
    gmx_molblock_t *mb;

    // atom 2 residue correspondence tables, both are filled, one will be used in the end
    std::vector<int> a2r_resnr(syslen_atoms);
    std::vector<int> a2r_renum(syslen_atoms);

    // get the maximum number of residues in any molecule block;
    // any residue nr. obtained via atoms->resinfo[].nr is guaranteed to be smaller than this;
    // this will be used for residue nr. collision detection - initialized to -1,
    // will be set to the renumbered value on the first encounter, then a new renumbered value means a collision
    int resnrmax = gmx_mtop_maxresnr(mtop, 0);
    std::vector<int> resnr2renum(resnrmax + 1, -1);

    int atom_global_index = 0;
    int renum = 0; //< renumbered residue nr.; increased monotonically, so could theoretically be as large as the nr. of atoms => type int
    bool bResnrCollision = false;
    for (int moltype_index = 0; moltype_index < mtop->nmolblock; ++moltype_index) {
        mb = &mtop->molblock[moltype_index];
        for (int mol_index = 0; mol_index < mb->nmol; ++mol_index) {
            atoms = &mtop->moltype[mb->type].atoms;
            for (int atom_index = 0; atom_index < atoms->nr; ++atom_index) {
                atom_info = &atoms->atom[atom_index];
                int resnr = atoms->resinfo[atom_info->resind].nr;
                renum = get_global_residue_number(mtop, atom_global_index);
                if ((resnr2renum[resnr] != renum) && (resnr2renum[resnr] != -1)) {
                    bResnrCollision = true;
                }
                resnr2renum[resnr] = renum;
                a2r_resnr[atom_global_index] = resnr;
                a2r_renum[atom_global_index] = renum;
                atom_global_index++;
            }
        }
    }

    // renum is set to the residue number of the last atom, so should be largest value so far
    switch (residues_renumber) {
        case ResiduesRenumber::AUTO:
            if (bResnrCollision) {
                std::cerr << "Residue number collision detected, residues will be renumbered." << std::endl;
                atom_2_residue = a2r_renum;
                syslen_residues = renum + 1;
            }
            else {
                std::cerr << "Residues will not be renumbered." << std::endl;
                atom_2_residue = a2r_resnr;
                syslen_residues = resnrmax + 1;
            }
            break;
        case ResiduesRenumber::DO:
            std::cerr << "Residues will be renumbered." << std::endl;
            atom_2_residue = a2r_renum;
            syslen_residues = renum + 1;
            break;
        case ResiduesRenumber::DONT:
            if (bResnrCollision) {
                std::cerr << "Residue number collision detected, residues will NOT be renumbered." << std::endl;
            } else {
                std::cerr << "Residues will not be renumbered." << std::endl;
                atom_2_residue = a2r_resnr;
                syslen_residues = resnrmax + 1;
            }
            break;
    }
}

int FDASettings::get_global_residue_number(gmx_mtop_t *mtop, int atnr_global) const
{
    int mb;
    int a_start, a_end, maxresnr, at_loc;
    t_atoms *atoms=NULL;

    mb = -1;
    a_end = 0;
    maxresnr = 0;

    // find the molecule block to which this atom belongs
    // maxresnr counts the number of residues encountered along the way
    do
    {
        if (mb >= 0) maxresnr += mtop->molblock[mb].nmol*atoms->nres;
        mb++;
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        a_start = a_end;
        a_end = a_start + mtop->molblock[mb].nmol*atoms->nr;
    }
    while (atnr_global >= a_end);

    // get the location of the atom inside the current block
    // as the block could appear repeated multiple times (repetition count stored in mtop->molblock[mb].nmol),
    // a simple substraction is not enough, modulo the total nr. of atoms in this block needs to be performed afterwards
    at_loc = (atnr_global - a_start) % atoms->nr;

    // (atnr_global - a_start)/atoms = in which repetition of the molblock this atom is found
    // above *atoms->nres = the residues corresponding to these repetitions
    // atoms->atom[at_loc].resind = residue index from the start of the molblock
    // original function had maxresnr + 1 + (...), to make residue numbering start from 1
    return maxresnr + (atnr_global - a_start)/atoms->nr*atoms->nres + atoms->atom[at_loc].resind;
}
