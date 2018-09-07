/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
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
 
 
#ifndef MOLPROP_H
#define MOLPROP_H

#include <string.h>

#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/real.h"

#include "communication.h"
#include "phase.h"
#include "poldata.h"

/*! \brief
 * Enumerated type holding the types of observables stored in MolProp
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropObservable {
    MPO_POTENTIAL,
    MPO_DIPOLE,
    MPO_QUADRUPOLE,
    MPO_POLARIZABILITY,
    MPO_ENERGY,
    MPO_ENTROPY,
    MPO_CHARGE,
    MPO_NR
};

/*! \brief
 * Enum to select either QM or Experimental data or either
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum iqmType {
    iqmExp, iqmBoth, iqmQM, iqmNR
};

/*! \breif
 * Strings describing the MolPropObservable enum elements
 */
extern const char *mpo_name[MPO_NR];

/*! \brief
 * Strings describing the MolPropObservable enum units
 */
extern const char *mpo_unit[MPO_NR];

/*! \brief
 * Contains all classes related to alexandria force field tools
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
namespace alexandria
{

/*! \brief
 * Enum describing the type of the QM job computed by the Gaussian software
 */
enum jobType {
    JOB_OPT       = 0,
    JOB_POP       = 1,
    JOB_POLAR     = 2,
    JOB_G2        = 3,
    JOB_G3        = 4,
    JOB_G4        = 5,
    JOB_CBSQB3    = 6,
    JOB_W1U       = 7,
    JOB_W1BD      = 8,
    JOB_SP        = 9,
    JOB_UNKNOWN   = 10,
    JOB_NR        = 11,
};

/*! \brief
 * Return string corresponding to the job type
 */
const char *jobType2string(alexandria::jobType jType);

/*! \brief
 * Strings describing the job type
 */
jobType string2jobType(const std::string &str);

/*! \brief
 * Enum describing the source of the data
 */
enum DataSource {
    dsExperiment, dsTheory
};

/*! \brief
 * Return string corresponding to the data source
 */
const char *dataSourceName(DataSource ds);

/*! \brief
 * Return DataSource corresponding to the data source string
 */
DataSource dataSourceFromName(const std::string &name);

/*! \brief
 * Specifies the name of an atom type and the number in a molecular composition
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomNum
{
    private:
        /*! \brief
         * Atom name
         */
        std::string catom_;
        /*! \brief
         * Atom number
         */
        int         cnumber_;
    public:
        //! Default constructor
        AtomNum() {};

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(const char *catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        /*! \brief
         * Creates a new AtomNum object.
         *
         * \param[in] catom   Atom name
         * \param[in] cnumber Number of copies of this atom
         */
        AtomNum(std::string catom, int cnumber) { SetAtom(catom); SetNumber(cnumber); }

        /*! \brief
         * Return the name of the atom for this AtomNum
         */
        const std::string &getAtom() const { return catom_; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(std::string catom) { catom_ = catom; }

        /*! \brief
         * Set the name of the atom for this AtomNum
         */
        void SetAtom(const char *catom) { catom_.assign(catom); }

        /*! \brief
         * Return the number of atoms for this AtomNum
         */
        int getNumber() const { return cnumber_; }

        /*! \brief
         * Set the number of atoms for this AtomNum
         */
        void SetNumber(int cnumber) { cnumber_ = cnumber; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over a vector of AtomNum
using  AtomNumIterator      = typename std::vector<AtomNum>::iterator;
using  AtomNumConstIterator = typename std::vector<AtomNum>::const_iterator;

/*! \brief
 * Contains the molecular composition in terms of atoms
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularComposition
{
    private:
        /*! \brief
         * Composition name
         */
        std::string          compname_;
        /*! \brief
         * A vector of AtomNum object
         */
        std::vector<AtomNum> atomnum_;
    public:
        //! Defult constructor
        MolecularComposition() {}

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(const char *compname) { compname_.assign(compname); }

        /*! \brief
         * Creates a new MolecularComposition object.
         *
         * \param[in] compname  Name of the composition type
         */
        MolecularComposition(std::string compname) { compname_ = compname; }

        /*! \brief 
         * Return the composition name
         */
        const std::string &getCompName() const { return compname_; }

        /*! \brief
         * Set the composition name
         */
        void SetCompName(std::string compname) { compname_ = compname; }

        /*! \brief
         * Set the composition name
         */
        void SetCompName(char *compname) { compname_.assign(compname); }

        /*! \brief
         * Add an AtomNum object to the composition
         *
         * \param[in] an  Atom number
         */
        void AddAtom(AtomNum an);

        /*! \brief
         * Remove the atom with name catom from the composition
         *
         *\param[in] catom   Atom name
         */
        void DeleteAtom(const char *catom)
        {
            std::string _str(catom); DeleteAtom(_str);
        }

        /*! \brief 
         * Remove the atom with name catom from the composition
         *
         * \param[in] catom Atom name
         */
        void DeleteAtom(std::string catom);

        /*! \brief 
         * Replace the oldatom by newatom
         *
         * \param[in] oldatom   Name of the old atom
         * \param[in] newatom   Name of the new atom
         */
        void ReplaceAtom(const char *oldatom, const char *newatom)
        {
            std::string so(oldatom), sn(newatom); ReplaceAtom(so, sn);
        }
        
        /*! \brief 
         * Replace the oldatom by newatom
         *
         * \param[in] oldatom   Name of the old atom
         * \param[in] newatom   Name of the new atom
         */
        void ReplaceAtom(std::string oldatom, std::string newatom);

        /*! \brief
         * Return iterator to begin looping over AtomNum
         */
        AtomNumIterator BeginAtomNum() { return atomnum_.begin(); }
        
        /*! \brief
         * Return iterator to end looping over AtomNum
         */
        AtomNumIterator EndAtomNum() { return atomnum_.end(); }

        /*! \brief
         * Return iterator pointing to a specific atom or EndAtomNum if not found
         *
         * \param[in] an Atom number
         */
        AtomNumIterator SearchAtom(std::string an);

        /*! \brief
         * Return the number of atoms of a certain type
         *
         *\param[in] atom Atom name
         */
        int CountAtoms(const char *atom);

        /*! \brief
         * Return the number of atoms of a certain type
         *
         *\param[in] atom Atom name
         */
        int CountAtoms(std::string atom);

        /*! \brief
         * Return the total number of atoms
         */
        int CountAtoms();

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularComposition items
using MolecularCompositionIterator      = typename std::vector<MolecularComposition>::iterator;
using MolecularCompositionConstIterator = typename std::vector<MolecularComposition>::const_iterator;

/*! \brief
 * Generic molecular property base clase
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class GenericProperty
{
    private:
        /*! \brief
         * Type of the property
         */
        std::string type_;
        /*! \brief
         * Unit of the property
         */
        std::string unit_;
        /*! \brief
         * Phase in which the property is measured or computed: e.x. Gas, Liquid, and Solid
         */
        ePhase      eP_;
        /*! \brief
         * Temperature at which the property is measured or computed.
         */
        double      T_;
    public:
        //! Default constructor
        GenericProperty() { T_ = 0; eP_ = epNR; };

        /*! \brief
         * Creates a new GenericProperty object.
         *
         * \param[in] type  Type of the property
         * \param[in] unit  Unit of the property
         * \param[in] T     Temperature
         */
        GenericProperty(std::string type, std::string unit, double T, ePhase ep)
        { SetType(type); SetUnit(unit); setTemperature(T); setPhase(ep); }

        /*! \brief
         * Return the property type
         */
        const std::string &getType() const { return type_; }
        
        /*! \brief
         * Return the unit of the property
         */
        const std::string &getUnit() const { return unit_; }
        
        /*! \brief
         * Return the temperature
         */
        double getTemperature() const { return T_; }

        /*! \brief 
         * Return the phase
         */
        ePhase getPhase() const { return eP_; }

        /*! \brief
         * Set the type of the property
         *
         *\param[in] type  Type of property
         */
        void SetType(std::string type);

        /*! \brief
         * Set the unit of the property
         *
         *\param[in] unit Unit of the property
         */
        void SetUnit(std::string unit);

        /*! \brief
         * Set the temperature of the property
         *
         *\param[in] T Temperature
         */
        void setTemperature(double T) { T_ = T; }

        /*! \brief
         * Set the phase of the property
         *
         *\param[in] ep Phase of the property
         */
        void setPhase(ePhase ep) { eP_ = ep; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};

/*! \brief
 * Contains the elements of the molecular quadrupole
 *
 * The six elements of the upper diagonal of a quadrupole tensor are stored.
 * The values are dependent on the orientation of the molecule and are relative
 * to the center of charge in case that the total charge or total dipole of
 * the molecules are non-zero.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularQuadrupole : public GenericProperty
{
    private:
        double xx_, yy_, zz_, xy_, xz_, yz_;
    public:
        //! Default constructor
        MolecularQuadrupole() {}

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularQuadrupole(std::string type, std::string unit, double T,
                            double xx, double yy, double zz,
                            double xy, double xz, double yz) : 
                            GenericProperty(type, unit, T, epGAS) 
                            { Set(xx, yy, zz, xy, xz, yz); };

        //! Set all the elements of the qudrupole tensor
        void Set(double xx, double yy, double zz, double xy, double xz, double yz) 
        { xx_ = xx; yy_ = yy; zz_ = zz; xy_ = xy; xz_ = xz; yz_ = yz; };

        //! get all the elements of the qudrupole tensor
        void get(double *xx, double *yy, double *zz, double *xy, double *xz, double *yz) 
        { *xx = xx_; *yy = yy_; *zz = zz_; *xy = xy_; *xz = xz_; *yz = yz_; };

        //! Return the XX component of the quadrupole tensor
        double getXX() const { return xx_; }

        //! Return the YY component of the quadrupole tensor
        double getYY() const { return yy_; }

        //! Return the ZZ component of the quadrupole tensor
        double getZZ() const { return zz_; }

        //! Return the XY component of the quadrupole tensor
        double getXY() const { return xy_; }

        //! Return the XZ component of the quadrupole tensor
        double getXZ() const { return xz_; }

        //! Return the YZ component of the quadrupole tensor
        double getYZ() const { return yz_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularQuadrupole items
using MolecularQuadrupoleIterator      = typename std::vector<MolecularQuadrupole>::iterator;
 using MolecularQuadrupoleConstIterator = typename std::vector<MolecularQuadrupole>::const_iterator;

/*! \brief
 * Contains the elements of the molecular polarizability tensor
 *
 * The six elements of the upper diagonal of a polarizability tensor are stored
 * along with the average molecular polarizability and the error if known.
 * The values are dependent on the orientation of the molecule.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularPolarizability : public GenericProperty
{
    private:
        double xx_, yy_, zz_, xy_, xz_, yz_, average_, error_;
    public:
        //! Default constructor
        MolecularPolarizability() {}

        //! Constructor initiating all elements of the quadrupole tensor
        MolecularPolarizability(std::string type, std::string unit, double T,
                                double xx, double yy, double zz,
                                double xy, double xz, double yz,
                                double average, double error) 
                                : GenericProperty(type, unit, T, epGAS) 
                                { Set(xx, yy, zz, xy, xz, yz, average, error); };

        //! Set all the elements of the polarizability tensor
        void Set(double xx, double yy, double zz,
                 double xy, double xz, double yz,
                 double average, double error);

        //! get all the elements of the polarizability tensor
        void get(double *xx, double *yy, double *zz,
                 double *xy, double *xz, double *yz,
                 double *average, double *error) const
        {
            *xx      = xx_; *yy = yy_; *zz = zz_;
            *xy      = xy_; *xz = xz_; *yz = yz_;
            *average = average_; *error = error_;
        }

        //! Return the XX component of the polarizability tensor
        double getXX() const { return xx_; }

        //! Return the YY component of the polarizability tensor
        double getYY() const { return yy_; }

        //! Return the ZZ component of the polarizability tensor
        double getZZ() const { return zz_; }

        //! Return the XY component of the polarizability tensor
        double getXY() const { return xy_; }

        //! Return the XZ component of the polarizability tensor
        double getXZ() const { return xz_; }

        //! Return the YZ component of the polarizability tensor
        double getYZ() const { return yz_; }

        //! Return the average of the polarizability tensor
        double getAverage() const { return average_; }

        //! Return the error in the polarizability tensor
        double getError() const { return error_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularPolarizability items
using  MolecularPolarizabilityIterator      = typename std::vector<MolecularPolarizability>::iterator;
using  MolecularPolarizabilityConstIterator = typename std::vector<MolecularPolarizability>::const_iterator;

/*! \brief
 * Contains a molecular energy
 *
 * Different energy terms associated with a molecule can be stored based
 * on quantum chemistry calculations or experimental data.
 * For example the Enthalpy of formation at different temperatures.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularEnergy : public GenericProperty
{
    private:
        double _value, error_;
    public:
        //! Default constructor
        MolecularEnergy() {};

        //! Constructor storing all properties related to this energy term
        MolecularEnergy(std::string type, std::string unit, 
                        double T, ePhase ep, double value, double error) 
                        : GenericProperty(type, unit, T, ep) 
                        { Set(value, error); };

        //! Set the value and error for the energy
        void Set(double value, double error) { _value = value; error_ = error; };

        //! get the value and error for this energy
        void get(double *value, double *error) const { *value = _value; *error = error_; };

        //! Set the value for the energy
        void SetValue(double value) { _value = value; };

        //! Return the energy value
        double getValue() const { return _value; };

        //! Set the error in the energy value
        void SetError(double error) { _value = error; };

        //! Return the error in the energy
        double getError() const { return error_; };

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolecularEnergy items
using  MolecularEnergyIterator      = typename std::vector<MolecularEnergy>::iterator;
using  MolecularEnergyConstIterator = typename std::vector<MolecularEnergy>::const_iterator;

/*! \brief
 * Contains the dipole vector
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolecularDipole : public GenericProperty
{
    private:
        double _x, _y, _z;
        double _aver, error_;
    public:
        //! Default constructor
        MolecularDipole() {}

        //! Constructor storing all properties related to this dipole
        MolecularDipole(std::string type, std::string unit, double T,
                        double x, double y, double z, double aver, double error) 
                        : GenericProperty(type, unit, T, epGAS) 
                        { Set(x, y, z, aver, error); }

        //! Set all properties related to this dipole
        void Set(double x, double y, double z, double aver, double error) 
        { _x = x; _y = y; _z = z; _aver = aver; error_ = error; };

        //! Return all properties of this dipole
        void get(double *x, double *y, double *z, double *aver, double *error) const
        { *x = _x; *y = _y; *z = _z; *aver = _aver; *error = error_; };

        //! Return the average dipole value
        double getAver() const { return _aver; }

        //! Return the error in the average dipole
        double getError() const { return error_; }

        //! Return the X component of the dipole
        double getX() const { return _x; }

        //! Return the Y component of the dipole
        double getY() const { return _y; }

        //! Return the Z component of the dipole
        double getZ() const { return _z; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over a vector of MolecularDipole
using MolecularDipoleIterator      = typename std::vector<MolecularDipole>::iterator;
using MolecularDipoleConstIterator = typename std::vector<MolecularDipole>::const_iterator;

/*! \brief
 * Contains the electrostatic potential in a coordinate close to a molecule.
 *
 * The electrostatic potential (ESP) can be computed using quantum chemistry and
 * stored in this class.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class ElectrostaticPotential
{
    private:
        std::string xyzUnit_, vUnit_;
        int         espID_;
        double      x_, y_, z_, V_;
    public:
        //! Default constructor
        ElectrostaticPotential() {}

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(std::string xyz_unit, std::string V_unit, 
                               int espid, double x, double y, double z, double V) 
                               { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Constructor that set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        ElectrostaticPotential(const char *xyz_unit, const char *V_unit, 
                               int espid, double x, double y, double z, double V) 
                               { Set(xyz_unit, V_unit, espid, x, y, z, V); };

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(std::string xyz_unit, std::string V_unit, 
                 int espid, double x, double y, double z, double V) 
                 { xyzUnit_ = xyz_unit; vUnit_ = V_unit; espID_ = espid; x_ = x; y_ = y; z_ = z; V_ = V; };

        //! Set the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void Set(const char *xyz_unit, const char *V_unit, int espid, double x, double y, double z, double V)
        {
            std::string x_(xyz_unit), V_(V_unit); Set(x_, V_, espid, x, y, z, V);
        }

        //! Return the units of coordinates and potential, the ESP id, the coordinates and the potential itself
        void get(char **xyz_unit, char **V_unit, int *espid, double *x, double *y, double *z, double *V) const
        {
            *xyz_unit = strdup(xyzUnit_.c_str()); *V_unit = strdup(vUnit_.c_str()); *espid = espID_; *x = x_; *y = y_; *z = z_; *V = V_;
        };

        //! Return the unit of the coordinates
        const std::string &getXYZunit() const { return xyzUnit_; }

        //! Return the unit of the potential
        const std::string &getVunit() const { return vUnit_; }

        //! Return the ESP id from the original calculation
        int getEspid() const { return espID_; }

        //! Return the X coordinate of the ESP point
        double getX() const { return x_; }

        //! Return the Y coordinate of the ESP point
        double getY() const { return y_; }

        //! Return the Z coordinate of the ESP point
        double getZ() const { return z_; }

        //! Return the electrostatic potential at this point in space
        double getV() const { return V_; }

        //! Set the potential for this instance of the class
        void SetV(double V) { V_ = V; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over ElectrostaticPotential items
using ElectrostaticPotentialIterator      = typename std::vector<ElectrostaticPotential>::iterator;
using ElectrostaticPotentialConstIterator = typename std::vector<ElectrostaticPotential>::const_iterator;

/*! \brief
 * Chemical bond in a molecule with associated bond order.
 *
 * The chemical bonds in a molecule are stored here along with the bond order
 * which together can be used for determining the atom types in a force field
 * calculation. In the present implementation this data is generated by OpenBabel.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Bond
{
    private:
        int ai_, aj_, bondorder_;
    public:
        //! Default constructor
        Bond() {}

        //! Constructor setting the ids of the atoms and the bondorder
        Bond(int ai, int aj, int bondorder) { Set(ai, aj, bondorder); }

        //! Sets the ids of the atoms and the bondorder
        void Set(int ai, int aj, int bondorder) {ai_ = ai; aj_ = aj; bondorder_ = bondorder; };

        //! Returns the ids of the atoms and the bondorder
        void get(int *ai, int *aj, int *bondorder) const
        { *ai = ai_; *aj = aj_; *bondorder = bondorder_; };

        //! Returns the first atom id
        int getAi() const { return ai_; }

        //! Returns the second atom id
        int getAj() const { return aj_; }

        //! Returns the bondorder
        int getBondOrder() const { return bondorder_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over Bond items
using BondIterator      = typename std::vector<Bond>::iterator;
using BondConstIterator = typename std::vector<Bond>::const_iterator;

/*! \brief
 * Contains the charge of an atom
 *
 * The charge of an atom can be derived from quantum chemistry codes in different
 * ways, despite it not being a physical observable.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomicCharge : public GenericProperty
{
    private:
        double q_;
    public:
        //! Default constructor
        AtomicCharge() {}

        //! Constructor setting type, unit and charge itself
        AtomicCharge(std::string type, std::string unit, double T, double q) 
            : GenericProperty(type, unit, T, epGAS) { SetQ(q); };

        //! Set the charge to q
        void SetQ(double q) { q_ = q; };

        //! Return the charge
        double getQ() const { return q_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over AtomicCharge items
using  AtomicChargeIterator      = typename std::vector<AtomicCharge>::iterator;
using  AtomicChargeConstIterator = typename std::vector<AtomicCharge>::const_iterator;

/*! \brief
 * Contains data on an atom based on a calculation.
 *
 * This class coordinates, name, atom type and an array of
 * AtomicCharge values based on different methods.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class CalcAtom
{
    private:
        std::string               name_, obType_, unit_;
        double                    x_, y_, z_;
        int                       atomID_;
        std::vector<AtomicCharge> q_;
    public:
        //! Default constructor
        CalcAtom() {}

        //! Constructor initiating the name, type and atomid
        CalcAtom(const char *name, const char *obtype, int atomid)
        {
            name_.assign(name); obType_.assign(obtype); atomID_ = atomid;
        };

        //! Constructor initiating the name, type and atomid
        CalcAtom(std::string name, std::string obtype, int atomid)
        {
            name_ = name; obType_ = obtype; atomID_ = atomid;
        };


        //! Function returning true if the two atoms are equal
        bool Equal(CalcAtom ca);

        //! Add an AtomicCharge element to the atom
        void AddCharge(AtomicCharge aq);

        //! Begin Iterator over AtomicCharge items
        AtomicChargeIterator BeginQ() { return q_.begin(); }

        //! End Iterator over AtomicCharge items
        AtomicChargeIterator EndQ() { return q_.end(); }

        //! Return the atom id of the atom
        int getAtomid() const { return atomID_; }

        //! Return the name of the atom
        const std::string &getName() const { return name_; }

        //! Return the OpenBabel type of the atom
        const std::string &getObtype() const { return obType_; }

        //! Return the unit of the coordinates of the atom
        const std::string &getUnit() const { return unit_; }

        //! Set the unit of the coordinates of the atom
        void SetUnit(std::string unit);

        //! Set the unit of the coordinates of the atom
        void SetUnit(const char *unit) { std::string s(unit); SetUnit(s); }

        //! Set the coordinates of the atom
        void SetCoords(double x, double y, double z) { x_ = x; y_ = y; z_ = z; }

        //! Return all the coordinates of the atom
        void getCoords(double *x, double *y, double *z) const
        { *x = x_; *y = y_; *z = z_; }

        //! Return the X coordinate of the atom
        double getX() const { return x_; }

        //! Return the Y coordinate of the atom
        double getY() const { return y_; }

        //! Return the Z coordinate of the atom
        double getZ() const { return z_; }

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over CalcAtom items
using  CalcAtomIterator      = typename std::vector<CalcAtom>::iterator;
using  CalcAtomConstIterator = typename std::vector<CalcAtom>::const_iterator;

/*! \brief
 * Contains molecular data based on experiments
 *
 * This is a composite class holding the results from experiment, either
 * the dipole, polarizability, energy or the quadrupole. A reference to the
 * publication (or handbook) containing the data is stored, as well as the
 * conformation of the molecule (if known).
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Experiment
{
    public:
        //! Default constructor
        Experiment() { }

        //! Constructor initiating an Experiment with reference and conformation
        Experiment(std::string reference, std::string conformation) :
            dataSource_(dsExperiment), reference_(reference),
            conformation_(conformation), jobtype_(JOB_UNKNOWN)
        {}

        //! Constructor initiating a Calculation
        Experiment(std::string program, std::string method,
                   std::string basisset, std::string reference,
                   std::string conformation, std::string datafile,
                   jobType jtype);

        //! Return the type of data
        DataSource dataSource() const { return dataSource_; }

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Add a MolecularDipole element containing polarizability data
        void AddPolar(MolecularPolarizability mdp) { polar_.push_back(mdp); }

        //! Return Begin Iterator over polarizability elements
        MolecularPolarizabilityIterator BeginPolar() { return polar_.begin(); }

        //! Return End Iterator over polarizability elements
        MolecularPolarizabilityIterator EndPolar()   { return polar_.end(); }

        //! Number of polarizability values
        int NPolar() { return polar_.size(); }

        //! Number of dipole values
        int NDipole() { return dipole_.size(); }

        //! Add a MolecularDipole element containing dipole data
        void AddDipole(MolecularDipole mdp) { dipole_.push_back(mdp); }

        //! Return Begin Iterator over dipole elements
        MolecularDipoleIterator BeginDipole() { return dipole_.begin(); }

        //! Return End Iterator over dipole elements
        MolecularDipoleIterator EndDipole()   { return dipole_.end(); }

        //! Number of quadrupole values
        int NQuadrupole() { return quadrupole_.size(); }

        //! Add a MolecularQuadrupole element
        void AddQuadrupole(MolecularQuadrupole mq) { quadrupole_.push_back(mq); }

        //! Return Begin Iterator over quadrupole elements
        MolecularQuadrupoleIterator BeginQuadrupole() { return quadrupole_.begin(); }

        //! Return End Iterator over quadrupole elements
        MolecularQuadrupoleIterator EndQuadrupole() { return quadrupole_.end(); }

        //! Add a MolecularEnergy element
        void AddEnergy(MolecularEnergy me) { energy_.push_back(me); }

        //! Return Begin Iterator over energy elements
        MolecularEnergyIterator BeginEnergy() { return energy_.begin(); }

        //! Return End Iterator over energy elements
        MolecularEnergyIterator EndEnergy()   { return energy_.end(); }

        //! Return the molecular conformation
        const std::string &getConformation() const { return conformation_; }

        //! Return the literature reference
        const std::string &getReference() const { return reference_; }

        //! Return the type of calculation
        const jobType &getJobtype() const { return jobtype_; }

        //! Add a CalcAtom object to the list of atoms
        void AddAtom(CalcAtom ca);

        //! Return the number of atoms
        int NAtom() { return catom_.size(); }

        //! Iterator Begin over CalcAtom objects
        CalcAtomIterator BeginAtom() { return catom_.begin(); }

        //! Iterator End over CalcAtom objects
        CalcAtomIterator EndAtom() { return catom_.end(); }

        //! Return iterator pointint to this particular atom or EndAtom() if not found
        CalcAtomIterator SearchAtom(CalcAtom ca);

        //! Add ElectrostaticPotential element to the array
        void AddPotential(ElectrostaticPotential ep) { potential_.push_back(ep); }

        //! Return the number of potential points
        int NPotential() { return potential_.size(); };

        //! Iterator Begin over ElectrostaticPotential objects
        ElectrostaticPotentialIterator BeginPotential() { return potential_.begin(); }

        //! Iterator End over ElectrostaticPotential objects
        ElectrostaticPotentialIterator EndPotential() { return potential_.end(); }

        //! Return the program used to perform the calculation
        const std::string &getProgram() const { return program_; }

        //! Return the basis set used to perform the calculation
        const std::string &getBasisset() const { return basisset_; }

        //! Return the method used to perform the calculation
        const std::string &getMethod() const { return method_; }

        //! Return the datafile from which the calculation output was extracted
        const std::string &getDatafile() const { return datafile_; }

        /*! \brief
         * Convenience function that fetches a value from this experiment
         *
         * \param[in]  type   The type of the data (dependent on whether it is dipole, energy etc.)
         *                    If size is 0, then this is ignored.
         * \param[in]  mpo    Enum selecting the type of data to fetch
         * \param[out] value  The value of e.g. the energy
         * \param[out] error  The error in the value of e.g. the energy
         * \param[out] T      Temperature
         * \param[out] vec    Vector data to be output, dipole or diagonal element of polarizability
         * \param[out] quadrupole The quadrupole tensor
         * \return true on success
         */
        bool getVal(const std::string  type, 
                    MolPropObservable  mpo,
                    double            *value,
                    double            *error, 
                    double            *T,
                    double             vec[3],
                    tensor             quadrupole);

        /*! \brief
         * Return the HF energy of the molecule
         */
        bool getHF(double *value);

        //! Merge in another object. Return number of warnings.
        int Merge(std::vector<Experiment>::iterator src);

        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);

    private:
        DataSource                           dataSource_;
        std::string                          reference_;
        std::string                          conformation_;
        std::string                          program_;
        std::string                          method_;
        std::string                          basisset_;
        std::string                          datafile_;
        jobType                              jobtype_;
        std::vector<CalcAtom>                catom_;
        std::vector<ElectrostaticPotential>  potential_;
        std::vector<MolecularDipole>         dipole_;
        std::vector<MolecularEnergy>         energy_;
        std::vector<MolecularQuadrupole>     quadrupole_;
        std::vector<MolecularPolarizability> polar_;
};
//! Iterates over Experiment items
using  ExperimentIterator = typename std::vector<Experiment>::iterator;
using  ExperimentConstIterator = typename std::vector<Experiment>::const_iterator;

/*! \brief
 * Contains molecular properties from a range of sources.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class MolProp
{
    private:
        int                               index_;
        double                            mass_;
        int                               charge_, multiplicity_;
        std::string                       formula_, texform_, molname_, iupac_, cas_, cid_, inchi_;
        std::vector<std::string>          category_;
        std::vector<MolecularComposition> mol_comp_;
        std::vector<Experiment>           exper_;
        std::vector<Bond>                 bond_;
    public:
        //! Construct a number MolProp object
        MolProp() { index_ = -1; mass_ = 0; charge_ = 0; multiplicity_ = 0; }

        /*! \brief
         * Check the internal consistency of this object
         *
         * \todo Implement this
         */
        void CheckConsistency();

        //! Set the index number for sorting
        void SetIndex(int index) { index_ = index; }

        //! Return the index number for sorting
        int getIndex() const { return index_; }

        //! Set the molecular mass
        void SetMass(double mass) { mass_ = mass; }

        //! Return the molecular mass
        double getMass() const { return mass_; }

        //! Set the total charge of the molecule
        void SetCharge(double charge) { charge_ = charge; }

        //! Return the total charge of the molecule
        int getCharge() const { return charge_; }

        //! Set the multiplicity of the molecule
        void SetMultiplicity(int multiplicity) { multiplicity_ = multiplicity; }

        //! Return the multiplicity  of the molecule
        int getMultiplicity() const { return multiplicity_; }

        /*! \brief
         * Merge the content of another MolProp into this one
         *
         * \param[in] mpi The object to be merged into the present one
         * \return Number of warnings
         * \todo Check and double check
         */
        int Merge(std::vector<MolProp>::iterator mpi);

        //! Dump the contents of this object to a file
        void Dump(FILE *fp);

        //! Set the LaTeX formula
        void SetTexFormula(const std::string &formula) { texform_.assign(formula); }

        //! Set the formula
        void SetFormula(const std::string &formula) { formula_.assign(formula); }

        //! Return the formula
        const std::string &formula() const { return formula_; }

        //! Return the LaTeX formula
        const std::string &getTexFormula() const;

        /*! \brief
         * Generate the chemical formula for this molecule based on atoms
         * present in a calculation
         *
         * \param[in] ap Data structure containing information about atoms
         * \todo Check and double check. If there is no calculation data no
         * formula can be generated
         */
        bool GenerateFormula(gmx_atomprop_t ap);

        //! Set the molname
        void SetMolname(const char *molname) { molname_.assign(molname); }

        //! Set the molname
        void SetMolname(std::string molname) { molname_ = molname; }

        //! Return the molname
        const std::string &getMolname() const { return molname_; }

        //! Set the IUPAC name
        void SetIupac(const char *iupac) { iupac_.assign(iupac); }

        //! Set the IUPAC name
        void SetIupac(std::string iupac) { iupac_ = iupac; }

        //! Return IUPAC name or, if not found, the molname
        const std::string &getIupac() const
        {
            if (iupac_.size() > 0)
            {
                return iupac_;
            }
            else
            {
                return molname_;
            }
        }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(const char *cas) { cas_.assign(cas); }

        //! Set the CAS (Chemical Abstract Service) identifier, see http://www.cas.org/
        void SetCas(std::string cas) { cas_.assign(cas); }

        //! Return the CAS (Chemical Abstract Service) identifier, see http:://www.cas.org
        const std::string &getCas() const { return cas_; }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(const char *cid) { cid_.assign(cid); }

        //! Set the CID (Chemspider identifier) see http:://www.chemspider.com
        void SetCid(std::string cid) { cid_ = cid; }

        //! Return the CID (Chemspider identifier) see http:://www.chemspider.com
        const std::string &getCid() const { return cid_; }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(const char *inchi) { inchi_.assign(inchi); }

        //! Set the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        void SetInchi(std::string inchi) { inchi_ = inchi; }

        //! Return the IUPAC International Chemical Identifier (InChI) see http://www.iupac.org/home/publications/e-resources/inchi.html
        const std::string &getInchi() const { return inchi_; }

        //! Convenience function
        bool getPropRef(MolPropObservable mpo, iqmType iQM,
                        const std::string &lot,
                        const std::string &conf,
                        const std::string &type,
                        double *value, double *error, double *T,
                        std::string &ref, std::string &mylot,
                        double vec[3], tensor quadrupole);

        //! And another one
        bool getProp(MolPropObservable mpo, iqmType iQM,
                     const std::string &lot,
                     const std::string &conf,
                     const std::string &type,
                     double *value, double *error, double *T);

        //! Returns true if the HF energy of the optimized geometry exists and returns the HF
        bool getOptHF(double *value);
        
        //! Returns the number of Opt and SP experiments for a molecule in allmols.dat
        int NOptSP();

        //! Add a classification category for this molecule
        void AddCategory(const std::string &category)
        {
            if (!SearchCategory(category)) { category_.push_back(category); }
        }

        //! Return the number of categories
        int NCategory() { return category_.size(); }

        //! Return iterator Begin over categories
        std::vector<std::string>::iterator BeginCategory() { return category_.begin(); }

        //! Return iterator End over categories
        std::vector<std::string>::iterator EndCategory() { return category_.end(); }

        //! Return true if catname is an existing category
        bool SearchCategory(const std::string &catname) const;

        //! Delete a composition type if present
        void DeleteComposition(const std::string &compname);

        //! Add a composition entry
        void AddComposition(MolecularComposition mc);

        std::vector<MolecularComposition> &MolComp()
        { return mol_comp_; }

        const std::vector<MolecularComposition> &MolComp() const
        { return mol_comp_; }

        //! Begin Iterator over MolecularCompostion items
        MolecularCompositionConstIterator BeginMolecularComposition() const { return mol_comp_.begin(); }

        //! End Iterator over MolecularCompostion items
        MolecularCompositionConstIterator EndMolecularComposition() const { return mol_comp_.end(); }

        //! Begin Iterator over MolecularCompostion items
        MolecularCompositionIterator BeginMolecularComposition() { return mol_comp_.begin(); }

        //! End Iterator over MolecularCompostion items
        MolecularCompositionIterator EndMolecularComposition()   { return mol_comp_.end(); }

        //! Last Iterator over MolecularCompostion items
        MolecularComposition *LastMolecularComposition()   { return &(mol_comp_.back()); }

        //! Search for particular MolecularCompostion item or return EndMolecularComposition if not found
        MolecularCompositionIterator SearchMolecularComposition(std::string str);

        //! Return number of atoms in the first composition if present, or 0 otherwise
        int NAtom();

        //! Routine to generate compositions based on calculation data
        bool GenerateComposition(const Poldata &pd);

        //! Returns boolean stating whether a particular composition is present
        bool HasComposition(const std::string &composition) const;

        //! Add a Bond element
        void AddBond(Bond b);

        //! Check whether a Bond element is present already
        bool BondExists(Bond b);

        //! Return the number of Bond elements
        int NBond() const { return bond_.size(); }

        //! Begin Iterator over Bond elements
        BondIterator BeginBond() { return bond_.begin(); }

        //! End Iterator over Bond elements
        BondIterator EndBond() { return bond_.end(); }

        //! Add an experiment
        void AddExperiment(Experiment myexp) { exper_.push_back(myexp); }

        void Stats()
        {
            printf("%s - %s - %d experiments\n",
                   molname_.c_str(), formula_.c_str(),
                   (int)exper_.size());
        }
        //! Return the number of experiments
        int NExperiment() const { return exper_.size(); }

        //! Iterator Begin over experiments
        ExperimentIterator BeginExperiment() { return exper_.begin(); }

        //! Iterator End over experiments
        ExperimentIterator EndExperiment() { return exper_.end(); }

        //! Return pointer to the last inserted experiment or nullptr if the number of experiments is zero
        Experiment *LastExperiment()
        {
            if (NExperiment() > 0)
            {
                return &(exper_.back());
            }
            else
            {
                return nullptr;
            }
        }

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndExperiment in case it is not found
        ExperimentIterator getLot(const char *lot);

        //! Return a calculation iterator corresponding to the level of theory (lot) parameter, or EndCalculation in case it is not found
        //! The operator should hold the requested observable of the type (can be nullptr)
        ExperimentIterator getLotPropType(const char       *lot,
                                          MolPropObservable mpo,
                                          const char       *type);
        /*! \brief
         * Sends this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] dest      Destination processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receives this object over an MPI connection
         *
         * \param[in] commrec   GROMACS data structure for MPI communication
         * \param[in] src       Source processor
         * \return the CommunicationStatus of the operation
         */
        CommunicationStatus Receive(t_commrec *cr, int src);
};
//! Iterates over MolProp items
using  MolPropIterator      = typename std::vector<MolProp>::iterator;
using  MolPropConstIterator = typename std::vector<MolProp>::const_iterator;

/*! \brief Utility to compare temperatures
 *
 * Compares two temperatures
 * \param[in] Tref The reference temperature
 * \param[in] T    The temperature to be tested
 * \return true if the reference T < 0, or the difference between the two
 *              T is negligable.
 */
bool bCheckTemperature(double Tref, double T);

}

#endif
