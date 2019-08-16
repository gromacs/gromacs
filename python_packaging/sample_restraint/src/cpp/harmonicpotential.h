//
// Created by Eric Irrgang on 10/13/17.
//

#ifndef GROMACS_HARMONICPOTENTIAL_H
#define GROMACS_HARMONICPOTENTIAL_H

#include <iostream>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

/*! \file
 * \brief Implement a harmonic pair force.
 *
 * Calculations and additional behavior is defined in harmonicpotential.cpp
 *
 * \todo This code has not been updated in a while...
 * This needs to be updated and tested more rigorously.
 *
 * Ref. https://github.com/kassonlab/gmxapi/issues/55
 *      https://github.com/kassonlab/gmxapi/issues/77
 *      https://github.com/kassonlab/gmxapi/issues/78
 */

namespace plugin
{

class Harmonic
{
    public:
        Harmonic(real equilibrium,
                 real springconstant) :
            R0_{equilibrium},
            k_{springconstant}
        {};

        Harmonic() :
            Harmonic{0.0, 0.0}
        {};

        // Allow easier automatic generation of bindings.
        struct input_param_type
        {
//             not yet used
        };

        /*!
         * \brief Calculate harmonic force on particle at position v in reference to position v0.
         *
         * \param v position at which to evaluate force
         * \param v0 position of harmonic bond reference
         * \return F = -k ((v - v0)/|v - v0| - R0);
         *
         * R0 == 1.0 is the equilibrium distance in the harmonic potential.
         * k == 1.0 is the spring constant.
         *
         * In the case of a pair of harmonically bonded particles, the force on particle i is evaluated with particle j as
         * the reference point with
         * \code
         * auto force = calculateForce(r_i, r_j);
         * \endcode
         *
         * The force on particle j is the opposite as the force vector for particle i. E.g.
         * \code
         * assert(-1 * force, calculateForce(r_j, r_i));
         * \endcode
         */
        gmx::PotentialPointData calculate(gmx::Vector v,
                                          gmx::Vector v0,
                                          gmx_unused double t);

        // The class will either be inherited as a mix-in or inherit a CRTP base class. Either way, it probably needs
        // proper virtual destructor management.
        virtual ~Harmonic()
        {
        }

    private:
        // set equilibrium separation distance in GROMACS units.
        // TODO: be clearer about units
        real R0_;
        // set spring constant in native GROMACS units.
        // TODO: be clearer about units
        real k_;
};

// implement IRestraintPotential in terms of Harmonic
// To be templated and moved.
class HarmonicRestraint : public ::gmx::IRestraintPotential, private Harmonic
{
    public:
        /*!
         * \brief Create an instance of the restraint (used in libgromacs)
         *
         * Each pair restraint instance operates on one pair of atomic sites.
         *
         * \param site1 first atomic site in the pair.
         * \param site2 second atomic site in the pair.
         * \param R0 targeted equilibrium pair separation.
         * \param k spring constant.
         */
        HarmonicRestraint(int site1,
                          int site2,
                          real R0,
                          real k) :
            Harmonic{R0, k},
            site1_{site1},
            site2_{site2}
        {};

        ~HarmonicRestraint() override = default;

        /*!
         * \brief Implement required interface of gmx::IRestraintPotential
         *
         * \return list of configured site indices.
         *
         * \todo remove to template header
         * \todo abstraction of site references
         */
        std::vector<int> sites() const override;

        /*!
         * \brief Implement the interface gmx::IRestraintPotential
         *
         * Dispatch to calculate() method.
         *
         * \param r1 coordinate of first site
         * \param r2 reference coordinate (second site)
         * \param t simulation time
         * \return calculated force and energy
         *
         * \todo remove to template header.
         */
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override;

    private:
        int site1_{0};
        int site2_{0};
};

/*!
 * \brief Wraps HarmonicPotential with a gmxapi compatible "module".
 *
 * Objects of this type allow the potential class to be instantiated as the simulation is launched.
 */
class HarmonicModule : public gmxapi::MDModule
{
    public:
        using param_t = Harmonic::input_param_type;

        HarmonicModule(int site1,
                       int site2,
                       real R0,
                       real k)
        {
            site1_ = site1;
            site2_ = site2;
            R0_ = R0;
            k_ = k;
        }


        const char* name() const override
        {
            return "HarmonicModule";
        }

        /*!
         * \brief implement gmxapi::MDModule::getRestraint()
         *
         * \return Handle to configured library object.
         */
        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<HarmonicRestraint>(site1_,
                                                                 site2_,
                                                                 R0_,
                                                                 k_);
            return restraint;
        }

        /*!
         * \brief Set restraint parameters.
         *
         * \todo generalize this
         * \param site1
         * \param site2
         * \param k
         * \param R0
         */
        void setParams(int site1,
                       int site2,
                       real R0,
                       real k)
        {
            site1_ = site1;
            site2_ = site2;
            R0_ = R0;
            k_ = k;
        }

    private:
        int site1_;
        int site2_;
        real R0_;
        real k_;
};

} // end namespace plugin

#endif //GROMACS_HARMONICPOTENTIAL_H
