/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Declares gmx::analysismodules::BasePair.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_BASEPAIRDB_H
#define GMX_TRAJECTORYANALYSIS_BASEPAIRDB_H

#include "gmxpre.h"

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"

namespace gmx
{

namespace analysismodules
{

//! Possible nucleotides including unknown
enum Nucleotide
{
    Cytosine, Guanine, Adenine, Uracil, Thymine, Nucnull
};

//! Designation of isomers
enum Isomerism
{
    Cis, Trans, Isonull
};

//! Return string corresponding to isomerism
const char *isomerismName(int iso);

//! The different dimer types
enum BondType
{
    Ww, Ss, Hh, Bh, Bs, bondSize, Bondnull
};

/*! \brief
 * Comparing basepairs residue names to a template.
 */
enum BasePairComparison
{
    BasePairComparison_MisMatch,
    BasePairComparison_Match,
    BasePairComparison_Swap,
    BasePairComparison_Both
};

//! Information about residues
class ResInfo
{
    public:
        /*! \brief
         * Constructor
         */
        ResInfo(char type, int resNr, int start) { type_ = type; resNr_ = resNr; start_ = start; }
        //! Destructor
        ~ResInfo() {}

        char residueType() const { return type_; }

        int residueNumber() const { return resNr_; }

        int atomStart() const { return start_; }

        unsigned int searchAtom(const char *name) const;

        unsigned int nAtoms() const { return atomName_.size(); }

        void addResInfoAtom(std::string name, bool invalid)
        {
            atomName_.push_back(name);
            invalid_.push_back(invalid);
        }

        std::string atomName(unsigned int i) const
        {
            GMX_ASSERT(i < atomName_.size(), "Atom number out of range");
            return atomName_[i];
        }

        bool invalid(unsigned int i) const
        {
            GMX_ASSERT(i < atomName_.size(), "Atom number out of range");
            return invalid_[i];
        }
    private:
        //! Residue type
        char                     type_;
        //! Residue number
        int                      resNr_;
        //! Start in the global atoms array
        int                      start_;
        //! Whether or not this particle is invalid
        std::vector<bool>        invalid_;
        //! Names of the particles
        std::vector<std::string> atomName_;
};


/*! \libinternal \brief
 * Utility class to store pairs of atom names
 */
class BondAtom
{
    public:
        BondAtom(std::string ai, std::string aj)
        {
            ai_ = ai;
            aj_ = aj;
        }

        ~BondAtom() {}

        std::string atomI() { return ai_; }

        std::string atomJ() { return aj_; }
    private:
        std::string ai_, aj_;
};

/*! \libinternal \brief
 * Utility class to store and operate on a nucleotide base pair
 */
class BasePair
{
    public:
        // Constructor
        BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd);

        // Destructor
        virtual ~BasePair();

        // Set the types of both nucleotide encoded in a two character string
        void setNucleotides(const std::string &nucleotides);

        // Return number of atoms
        unsigned int nrAtoms() const { return atomName_.size(); }

        // Add an atom
        void addBasePairAtom(rvec               x,
                             int                atomnum,
                             const std::string &atomname,
                             const std::string &resname,
                             real               m);

        void resetToOrigin()
        {
            reset_x_ndim(3, nrAtoms(),  NULL, nrAtoms(), NULL, x_, atomMass());
        }

        // Return the atom positions
        ConstArrayRef<rvec> atomPosition() const
        {
            return constArrayRefFromArray(x_, nrAtoms());
        }

        // Add information about a pair of bonded atoms
        void addBondAtoms(int ai, int aj);

        // Return enum describing the match with the template
        BasePairComparison checkTemplate(char nuc0, char nuc1);

        // Calculate the RMSD
        real computeRootMeanSquareDeviation(rvec * vec, int baseSize);

        // Return the maximun distance in nm TODO: between what?
        real maximumDistance() { return dist; }

        // Set the maximum distance to another atom
        void setAtomDist(int index);

        // Return phosphate
        bool phosphateRmsd() { return phosphateRmsd_; }

        // Check if an atom is valid
        bool checkAtom (const std::string &name);

        // Return nucleotide type. Index should be 0 or 1 or the routine will throw.
        char nucleotideType (unsigned int index);

        // Set bond type
        void setBondType (const std::string &bondType);

        // Get bond type
        int bondType();

        // Get bond name
        std::string bondName(unsigned int type);

        // Check whether a bond exist by calculating the H bond distance
        bool isHBond(const ResInfo *tempID1,
                     const ResInfo *tempID2,
                     const t_pbc   *pbc,
                     real           bondDist,
                     Selection      sel);

        // Set the isomerism type
        void setIsomerism (std::string name);

        // Get the isomerism type
        int isomerismType() { return iso; };

        // Get the mass
        real *atomMass();

        // Returns the bond index nr bond from pair
        const char * bondIndex(unsigned int bond, int pair) const;

        // Returns the number of bonded atoms
        unsigned int nrBondAtoms() const { return bondAtoms.size(); }

        /*! \brief
         * Set the periodic boundary condition variables
         */
        void setPBC(int ePBC, matrix box) { ePBC_ = ePBC; copy_mat(box_, box); }

        /*! \brief
         * \return a pointer to the box
         */
        matrix *box() { return &box_; }

        // Get ePBC
        int  getEpbc()
        {
            return ePBC_;
        }

        // Get the RMSD of the template
        real getTemplateRMSD() { return templateRMSD; }

        // set the RMSD of the template
        void setTemplateRMSD(real rmsd)
        {
            templateRMSD = rmsd;
            return;
        }

    private:
        real                      dist;
        Nucleotide                type[2];
        std::vector <std::string> atomName_;
        std::vector <int>         atomNumber_;
        std::vector<real>         atomMass_;
        rvec                     *x_;
        Isomerism                 iso;
        BondType                  bondtype[2];
        bool                      hydrogenRmsd_;
        bool                      sugarRmsd_;
        bool                      phosphateRmsd_;
        std::vector <BondAtom *>  bondAtoms;
        int                       resId[2];
        real                      templateRMSD;
        //! PBC type
        int                       ePBC_;
        //! Box matrix
        matrix                    box_;

        // Return type of atom
        char atomType(const char name[]);

        // Check whether the C1' distance exceed the template threshold
        bool isC1Contact(const ResInfo *tempID1,
                         const ResInfo *tempID2,
                         Selection      sel);

        // Get the tyoe of the bond name
        BondType searchBondType(std::string name);

};

} // namespace analysismodules

} // namespace gmx

#endif
