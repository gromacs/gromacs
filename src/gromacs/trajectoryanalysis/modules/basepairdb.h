/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::analysismodules::BasePair and help classes
 *
 * This class is to manage the nucleic acid base pair database
 * entries.
 *
 * \author Nina Fischer <nina.fischer@icm.uu.se>
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_BASEPAIRDB_H
#define GMX_TRAJECTORYANALYSIS_BASEPAIRDB_H

#include "gmxpre.h"

#include <gromacs/math/do_fit.h>
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
    Cytosine, Guanine, Adenine, Uracil, Thymine, NucleotideSize = Thymine
};

//! Convert one letter nucleotide to enum
Nucleotide searchNucleotide(char c);

//! Convert one enum nucleotide to letter
char nucleotideToChar(Nucleotide nuc);

const unsigned int numberOfNucleotideAtoms[] = {31, 34, 33, 30};

enum EndDirectionality
{
    NotEnd, FivePrimeEnd, ThreePrimeEnd
};

//! Designation of isomers
enum Isomerism
{
    Cis, Trans, IsomerismSize = Trans
};

//! Return string corresponding to isomerism
const char *isomerismName(Isomerism iso);

//! The different dimer types
enum BondType
{
    Ww, Ss, Hh, Bh, Bs, BondTypeSize = Bs
};

//! Return string corresponding to bondtype
const char *bondTypeName(BondType bt);
//! Get the type of the bond name
BondType searchBondType(std::string name);

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
         * Constructor convert character to enum
         */
        ResInfo(char type, int resNr, int start)
        {
            type_  = searchNucleotide(type);
            resNr_ = resNr;
            start_ = start;
            endD_  = NotEnd;
        }

        //! Return the residue type
        Nucleotide residueType() const { return type_; }

        //! Return the charactor of above.
        char residueTypeChar() const { return nucleotideToChar(type_); }

        //! Return the residue number
        int residueNumber() const { return resNr_; }

        //! Return the index of the first atom in the coordinates
        int atomStart() const { return start_; }

        //! Return index of a particular atom
        unsigned int searchAtom(const char *name) const;

        //! Return number of atoms
        unsigned int nAtoms() const { return atomName_.size(); }

        //! Add an atom and set the validity flag
        void addResInfoAtom(std::string name, bool invalid)
        {
            atomName_.push_back(name);
            if (name.compare("H5T") == 0)
            {
                endD_   = FivePrimeEnd;
                invalid = true;
            }
            if (name.compare("H3T") == 0)
            {
                endD_   = ThreePrimeEnd;
                invalid = true;
            }
            invalid_.push_back(invalid);
        }

        //! Return the atom name corresponding to atom i
        std::string atomName(unsigned int i) const
        {
            GMX_ASSERT(i < atomName_.size(), "Atom number out of range");
            return atomName_[i];
        }

        //! Return whether an atom number is within range
        bool invalid(unsigned int i) const
        {
            GMX_ASSERT(i < atomName_.size(), "Atom number out of range");
            return invalid_[i];
        }

        EndDirectionality endDirectionality() const { return endD_; }

    private:
        //! Residue type
        Nucleotide               type_;
        //! Residue number
        int                      resNr_;
        //! Start in the global atoms array
        int                      start_;
        //! Whether or not this particle is invalid
        std::vector<bool>        invalid_;
        //! Names of the particles
        std::vector<std::string> atomName_;
        //! Whether or not the residue is at an end.
        EndDirectionality        endD_;
};


/*! \libinternal \brief
 * Utility class to store pairs of atom names
 */
class BondAtom
{
    public:
        BondAtom(std::string ai, std::string aj) : ai_(ai), aj_(aj) {}

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
        //! Constructor
        BasePair(bool hydrogenRmsd, bool sugarRmsd, bool phosphateRmsd);

        //! Destructor
        virtual ~BasePair();

        //! Set the types of both nucleotide encoded in a two character string
        void setNucleotides(const std::string &nucleotides);

        //! Return number of atoms
        unsigned int nrAtoms() const { return atomName_.size(); }

        /*! \brief Add an atom
         *
         * Adds an atom to the base pair when reading the database.
         * Only atoms taking part in the weighted RMSD comparison
         * are stored (as determined with command line flags in the
         * rnaAnalysis tool).
         *
         * Throws if the residue name of the first atom does not
         * match the base pair type set earlier.
         *
         * \param[in] x        The coordinates
         * \param[in] atomnum  The atomnumber as stored in the database
         * \param[in] atomname The atomname as stored in the database
         * \param[in] resname  The residue name as stored in the database
         * \param[in] resnum   The residue number
         * \param[in] m        The mass of the atom
         */
        void addBasePairAtom(rvec               x,
                             unsigned int       atomnum,
                             const std::string &atomname,
                             const std::string &resname,
                             unsigned int       resnum,
                             real               m);

        //! Reset base pair to be centered on the origin
        void resetToOrigin()
        {
            reset_x_ndim(3, nrAtoms(),  NULL, nrAtoms(), NULL, x_, atomMass());
        }

        //! Return a pointer to the atom positions
        ConstArrayRef<rvec> atomPosition() const
        {
            return constArrayRefFromArray(x_, nrAtoms());
        }

        /*! \brief Add information about a connection
         *
         * The CONECT records in the database are used to
         * store the hydrogen bonds between bases. These
         * are stored here.
         *
         * Will throw if the atom indices do not correspond to
         * known atoms in the nucleotides, or if the atoms belong
         * to the same residue.
         *
         * \param[in] ai The first atom number (in first nucleotide)
         * \param[in] aj The second atom number (in second nucleotide)
         */
        void addBondAtoms(unsigned int ai, unsigned int aj);

        /*! \brief Return enum describing the match with the template
         *
         * Checks whether the nucleotides passed match those in
         * the template stored here.
         * \param[in] nuc0 The first nucleotide
         * \param[in] nuc1 The second nucleotide
         * \return Enum describing the match
         */
        BasePairComparison checkTemplate(Nucleotide nuc0, Nucleotide nuc1);

        /*! \brief Calculate the RMSD between template and pass coordinates
         *
         * \param[in] vec The coordinate array
         * \param[in] baseSize Number of elements in the coordinates
         * \return The root mean square deviation
         */
        real computeRootMeanSquareDeviation(const ResInfo *tempID1,
                                            const ResInfo *tempID2,
                                            rvec *vec, unsigned int baseSize);

        // Return the maximun distance in nm TODO: between what?
        real maximumDistance() { return dist; }

        // Set the maximum distance to another atom
        void setAtomDist(int index);

        // Return phosphate
        bool phosphateRmsd() { return phosphateRmsd_; }

        // Check if an atom is valid
        bool checkAtom (const std::string &name);

        // Return nucleotide type. Index should be 0 or 1 or the routine will throw.
        Nucleotide nucleotideType (unsigned int index);

        // Return the charactor of above.
        char nucleotideTypeChar (unsigned int index);

        //! Set bond types
        void setBondTypes(const std::string &bondTypes);

        //! Get bond pair index
        unsigned int bondPairIndex()
        {
            return (static_cast<unsigned int>(bondtype_[0])*(1+BondTypeSize) +
                    static_cast<unsigned int>(bondtype_[1]));
        }
        //! Return the bond name
        std::string bondName();

        //! Return a string formatted in the same way as the template database
        std::string templateName();

        // Check whether a bond exist by calculating the H bond distance
        bool isHBond(const ResInfo   *tempID1,
                     const ResInfo   *tempID2,
                     const t_pbc     *pbc,
                     real             bondDist,
                     const Selection &sel);

        // Set the isomerism type
        void setIsomerism (std::string name);
        // Get the isomerism type
        Isomerism isomerismType() { return iso; };
        // Get the split point of two residues.
        unsigned int splitPoint() { return splitPoint_; };

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
        real                        dist;
        Nucleotide                  type[2];
        std::vector <std::string>   atomName_;
        std::vector <unsigned int>  atomNumber_;
        std::vector <unsigned int>  resNumber_;
        std::vector <real>          atomMass_;
        rvec                       *x_;
        Isomerism                   iso;
        BondType                    bondtype_[2];
        bool                        hydrogenRmsd_;
        bool                        sugarRmsd_;
        bool                        phosphateRmsd_;
        std::vector <BondAtom *>    bondAtoms;
        int                         resId[2];
        real                        templateRMSD;
        //! PBC type
        int                         ePBC_;
        //! Box matrix
        matrix                      box_;
        unsigned int                splitPoint_;

        // Return type of atom
        char atomType(const char name[]);

        // Check whether the C1' distance exceed the template threshold
        bool isC1Contact(const ResInfo *tempID1,
                         const ResInfo *tempID2,
                         Selection      sel);


};

} // namespace analysismodules

} // namespace gmx

#endif
