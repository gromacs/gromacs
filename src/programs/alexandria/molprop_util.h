/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include "molprop.h"
#include "molselect.h"

struct t_topology;

/*! \brief
 * Enumerated type for MolPropSort function
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropSortAlgorithm {
    MPSA_MOLNAME     = 0,
    MPSA_FORMULA     = 1,
    MPSA_COMPOSITION = 2,
    MPSA_SELECTION   = 3,
    MPSA_NR          = 4
};

namespace alexandria
{
class QmCalc
{
    private:
        std::string method_, basis_, type_, lot_;
        int         count_;
    public:
        QmCalc(const std::string &method,
               const std::string &basis,
               const std::string &type) :
            method_(method), basis_(basis), type_(type), lot_(method), count_(1)
        {
            lot_.append("/");
            lot_.append(basis);
        };

        const std::string &method() const { return method_; }
        const std::string &basis() const { return basis_; }
        const std::string &type() const { return type_; }
        const std::string &lot() const { return lot_; }
        int count() const { return count_; }
        void increment() { count_++; }
        void decrement()
        {
            GMX_RELEASE_ASSERT(count_ > 0, "Trying to reduce count below zero");
            count_--;
        }
};

class QmCount
{
    private:
        std::vector<QmCalc>      qmc_;
        std::vector<std::string> conf_;

        std::vector<QmCalc>::const_iterator findCalc(const std::string &method,
                                                     const std::string &basis,
                                                     const std::string &type) const
        {
            return std::find_if(qmc_.begin(), qmc_.end(),
                                [method, basis, type](const QmCalc &qmc)
                                { return (qmc.method().compare(method) == 0 &&
                                          qmc.basis().compare(basis) == 0 &&
                                          qmc.type().compare(type) == 0); });
        }
        std::vector<QmCalc>::iterator findCalc(const std::string &method,
                                               const std::string &basis,
                                               const std::string &type)
        {
            return std::find_if(qmc_.begin(), qmc_.end(),
                                [method, basis, type](const QmCalc &qmc)
                                { return (qmc.method().compare(method) == 0 &&
                                          qmc.basis().compare(basis) == 0 &&
                                          qmc.type().compare(type) == 0); });
        }

    public:
        QmCount() {};

        std::vector<QmCalc>::iterator beginCalc() { return qmc_.begin(); }
        std::vector<QmCalc>::iterator endCalc() { return qmc_.end(); }
        std::vector<QmCalc>::const_iterator beginCalc() const { return qmc_.begin(); }
        std::vector<QmCalc>::const_iterator endCalc() const { return qmc_.end(); }
        size_t nCalc() const { return qmc_.size(); }
        int qmCalcCount(const std::string &method,
                        const std::string &basis,
                        const std::string &type) const;

        void addConf(const std::string &conformation);

        void addCalc(const std::string &method,
                     const std::string &basis,
                     const std::string &type);

};



void generate_composition(std::vector<MolProp> &mp,
                          const Poldata        &pd);

void generate_formula(std::vector<MolProp> &mp,
                      gmx_atomprop_t        ap);

int merge_doubles(std::vector<alexandria::MolProp> &mp,
                  char *doubles, bool bForceMerge);

int merge_xml(int nfile, char **infiles,
              std::vector<alexandria::MolProp> &mp,
              char *outf, char *sorted, char *doubles,
              gmx_atomprop_t ap,
              const Poldata &pd,
              bool bForceMerge);

/* Check the available molprops to see what kind of calculations are stored in there */
void find_calculations(std::vector<alexandria::MolProp> &mp,
                       MolPropObservable                 mpo,
                       const char                       *fc_str,
                       QmCount                          *qmc);

/*! \brief
 * Sorts a vector of molprops
 *
 * Function that uses the std::sort routine and can apply different sorting
 * keys.
 *
 * \param[inout]  mp        The vector of MolProp
 * \param[in]     mpsa      The algorithm used for sorting
 * \param[in]     apt       Database of atom properties
 * \param[in]     mgs       Optional structure containing selection criteria
 * \ingroup module_alexandria
 */
void MolPropSort(std::vector<MolProp> &mp,
                 MolPropSortAlgorithm  mpsa,
                 gmx_atomprop_t        apt,
                 const MolSelect      &gms);

} // namespace alexandria

#endif
