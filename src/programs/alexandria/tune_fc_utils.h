/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef TUNE_FC_UTILS_H
#define TUNE_FC_UTILS_H

#include <string>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "communication.h"
#include "mymol.h"
#include "plistwrapper.h"
#include "poldata.h"

struct t_commrec;

namespace alexandria
{

class AtomTypes
{
    public:

        AtomTypes () {}

        AtomTypes(int                ncopies,
                  const std::string &name,
                  const std::string &vdwParams,
                  int                index)
            :
              ncopies_(ncopies),
              name_(name),
              vdwParams_(vdwParams),
              poldataIndex_(index)
        {
            extractParams();
        }

        void inc() { ncopies_++; }

        int nCopies() const { return ncopies_; }

        void setParamString(const std::string &params);

        int poldataIndex() const { return poldataIndex_; }

        const std::string &paramString() const { return vdwParams_; }

        const std::vector<double> &paramValues() const { return p_; }

        const std::string &name() const { return name_; }

        size_t nParams() const { return p_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        int                 ncopies_;
        std::string         name_;
        std::string         vdwParams_;
        std::vector<double> p_;
        int                 poldataIndex_;
        void extractParams();
};

using AtomTypesIterator = typename std::vector<AtomTypes>::iterator;

class NonBondParams
{
    public:

        NonBondParams () {}

        NonBondParams(bool bOpt, InteractionType  itype)

            :
              bOpt_(bOpt),
              itype_(itype)
        {}


        void addNonBonded(AtomTypes at) { at_.push_back(std::move(at)); }

        void analyzeIdef(const std::vector<MyMol> &mm,
                         const Poldata            *pd);

        void makeReverseIndex();

        int reverseIndex(int poldataIndex)
        {
            GMX_RELEASE_ASSERT(poldataIndex >= 0 && poldataIndex < static_cast<int>(reverseIndex_.size()), "Incorrect poldataIndex");
            GMX_RELEASE_ASSERT(reverseIndex_[poldataIndex] != -1, "The reverseIndex is incorrect");

            return reverseIndex_[poldataIndex];
        }

        AtomTypesIterator beginAT() { return at_.begin(); }

        AtomTypesIterator endAT() { return at_.end(); }

        InteractionType interactionType() const { return itype_; }

        size_t nAT() const { return at_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        bool                   bOpt_;
        InteractionType        itype_;
        std::vector<AtomTypes> at_;
        std::vector<int>       reverseIndex_;
        std::vector<double>    params_;

};

/*! \brief Helper class storing bond/angle/dihedral names
 *
 * For one bond/angle/dihedral here the name of the bondtypes
 * are stored as in e.g. c c h for an angle, along with the number
 * of occurrences in the force field.
 */
class BondNames
{
    public:

        BondNames () {}

        BondNames(int                ncopies,
                  int                ftype,
                  const std::string &name,
                  double             geometry,
                  const std::string &params,
                  int                index,
                  double             bondorder = 0)

            :
              ncopies_(ncopies),
              ftype_(ftype),
              name_(name),
              geometry_(geometry),
              params_(params),
              bondorder_(bondorder),
              poldataIndex_(index)
        {
            extractParams();
        }

        void inc() { ncopies_++; }

        int nCopies() const { return ncopies_; }

        void setParamString(const std::string &params);

        const std::string &name() const { return name_; }

        double geometry() const { return geometry_; }
        
        void setGeometry(double geometry) { geometry_ = geometry; }
        
        double bondorder() const { return bondorder_; }

        int poldataIndex() const { return poldataIndex_; }

        const std::string &paramString() const { return params_; }

        const std::vector<double> &paramValues() const { return p_; }

        size_t nParams() const { return p_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:

        //! Number of copies in the molecule data set
        int                 ncopies_;
        //! Function type for this particular bond
        int                 ftype_;
        //! Name of this bond/angle/dihedral
        std::string         name_;
        //! Reference bond length or (dihedral) angle
        double              geometry_;
        //! String holding all the parameters
        std::string         params_;
        //! Vector containing all the parameters
        std::vector<double> p_;
        //! The bond order in case this is a bond
        double              bondorder_;
        //! Index in Poldata structure
        int                 poldataIndex_;
        //! Internal routine to extract the parameters
        void extractParams();
};

using BondNamesIterator = typename std::vector<BondNames>::iterator;

/*! \brief Class holding for one type of interactions all names
 *
 * Class holding the OptNames for each interaction type.
 */
class ForceConstants
{

    public:

        ForceConstants () {}

    ForceConstants(int ftype, 
                   InteractionType itype, 
                   bool bOpt)
            :
              ftype_(ftype),
              itype_(itype),
              bOpt_(bOpt)
        {
        }

        void addForceConstant(BondNames bn) { bn_.push_back(std::move(bn)); }

        void analyzeIdef(const std::vector<MyMol> &mm,
                         const Poldata            *pd);

        /*! \brief Make reverse index from Poldata to BondNames
         *
         * The BondNames structure stores the Poldata index for
         * all interactions. This routine makes an index to convert
         * the Poldata index to the index in BondNames.
         */
        void makeReverseIndex();

        int reverseIndex(int poldataIndex)
        {
            GMX_RELEASE_ASSERT(poldataIndex >= 0 && poldataIndex < static_cast<int>(reverseIndex_.size()), "Incorrect poldataIndex");
            GMX_RELEASE_ASSERT(reverseIndex_[poldataIndex] != -1, "The reverseIndex is incorrect");

            return reverseIndex_[poldataIndex];
        }

        int ftype() const { return ftype_; }

        InteractionType interactionType() const { return itype_; }

        void dump(FILE *fp) const;

        BondNamesIterator beginBN() { return bn_.begin(); }

        BondNamesIterator endBN() { return bn_.end(); }

        size_t nbad() const { return bn_.size(); }

        CommunicationStatus Send(t_commrec *cr, int dest);

        CommunicationStatus Receive(t_commrec *cr, int src);

    private:
        int                    ftype_;
        InteractionType        itype_;
        bool                   bOpt_;
        std::vector<BondNames> bn_;
        std::vector<int>       reverseIndex_;
        std::vector<double>    params_;
};


class PoldataUpdate
{
public:
    PoldataUpdate() {}
    PoldataUpdate(InteractionType     iType,
                  int                 index,
                  double              geometry,
                  std::string         paramString) : iType_(iType), index_(index), geometry_(geometry), paramString_(paramString)
    {}
    
    /*! \brief
     * Implement the changes in Poldata.
     *
     * \param[inout] pd The Poldata structure
     */
    void execute(Poldata *pd);
    /*! \brief
     * Dump the contents of the structure to a file.
     * \param[in] fp File pointer to dumpt to if not nullptr
     */
    void dump(FILE *fp) const;
    
    CommunicationStatus Send(t_commrec *cr, int dest);
    
    CommunicationStatus Receive(t_commrec *cr, int src);

private:
    InteractionType     iType_;
    int                 index_;
    double              geometry_;
    std::string         paramString_;
};


} // namespace alexandria
#endif
