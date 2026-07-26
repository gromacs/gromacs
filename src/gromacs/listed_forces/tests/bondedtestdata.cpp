/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*! \internal \file
 * \brief Implementations for shared bonded test data and utilities.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "bondedtestdata.h"

#include <algorithm>
#include <iterator>
#include <unordered_map>

#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"

namespace gmx
{
namespace test
{

// iListInput setter implementations
iListInput iListInput::setHarmonic(InteractionFunction ft, real rA, real krA, real rB, real krB)
{
    iparams.harmonic.rA  = rA;
    iparams.harmonic.rB  = rB;
    iparams.harmonic.krA = krA;
    iparams.harmonic.krB = krB;
    ftype                = ft;
    fep                  = (rA != rB || krA != krB);
    return *this;
}

iListInput iListInput::setHarmonic(InteractionFunction ft, real rA, real krA)
{
    return setHarmonic(ft, rA, krA, rA, krA);
}

iListInput iListInput::setCubic(real b0, real kb, real kcub)
{
    ftype              = InteractionFunction::CubicBonds;
    iparams.cubic.b0   = b0;
    iparams.cubic.kb   = kb;
    iparams.cubic.kcub = kcub;
    return *this;
}

iListInput iListInput::setMorse(real b0A, real cbA, real betaA, real b0B, real cbB, real betaB)
{
    ftype               = InteractionFunction::MorsePotential;
    iparams.morse.b0A   = b0A;
    iparams.morse.cbA   = cbA;
    iparams.morse.betaA = betaA;
    iparams.morse.b0B   = b0B;
    iparams.morse.cbB   = cbB;
    iparams.morse.betaB = betaB;
    fep                 = (b0A != b0B || cbA != cbB || betaA != betaB);
    return *this;
}

iListInput iListInput::setMorse(real b0A, real cbA, real betaA)
{
    return setMorse(b0A, cbA, betaA, b0A, cbA, betaA);
}

iListInput iListInput::setFene(real bm, real kb)
{
    ftype           = InteractionFunction::FENEBonds;
    iparams.fene.bm = bm;
    iparams.fene.kb = kb;
    return *this;
}

iListInput iListInput::setLinearAngle(real klinA, real aA, real klinB, real aB)
{
    ftype                  = InteractionFunction::LinearAngles;
    iparams.linangle.klinA = klinA;
    iparams.linangle.aA    = aA;
    iparams.linangle.klinB = klinB;
    iparams.linangle.aB    = aB;
    fep                    = (klinA != klinB || aA != aB);
    return *this;
}

iListInput iListInput::setLinearAngle(real klinA, real aA)
{
    return setLinearAngle(klinA, aA, klinA, aA);
}

iListInput iListInput::setUreyBradley(real thetaA,
                                      real kthetaA,
                                      real r13A,
                                      real kUBA,
                                      real thetaB,
                                      real kthetaB,
                                      real r13B,
                                      real kUBB)
{
    ftype               = InteractionFunction::UreyBradleyPotential;
    iparams.u_b.thetaA  = thetaA;
    iparams.u_b.kthetaA = kthetaA;
    iparams.u_b.r13A    = r13A;
    iparams.u_b.kUBA    = kUBA;
    iparams.u_b.thetaB  = thetaB;
    iparams.u_b.kthetaB = kthetaB;
    iparams.u_b.r13B    = r13B;
    iparams.u_b.kUBB    = kUBB;
    fep                 = (thetaA != thetaB || kthetaA != kthetaB || r13A != r13B || kUBA != kUBB);
    return *this;
}

iListInput iListInput::setUreyBradley(real thetaA, real kthetaA, real r13A, real kUBA)
{
    return setUreyBradley(thetaA, kthetaA, r13A, kUBA, thetaA, kthetaA, r13A, kUBA);
}

iListInput iListInput::setCrossBondBonds(real r1e, real r2e, real krr)
{
    ftype                = InteractionFunction::CrossBondBonds;
    iparams.cross_bb.r1e = r1e;
    iparams.cross_bb.r2e = r2e;
    iparams.cross_bb.krr = krr;
    return *this;
}

iListInput iListInput::setCrossBondAngles(real r1e, real r2e, real r3e, real krt)
{
    ftype                = InteractionFunction::CrossBondAngles;
    iparams.cross_ba.r1e = r1e;
    iparams.cross_ba.r2e = r2e;
    iparams.cross_ba.r3e = r3e;
    iparams.cross_ba.krt = krt;
    return *this;
}

iListInput iListInput::setQuarticAngles(real theta, const real c[5])
{
    ftype                = InteractionFunction::QuarticAngles;
    iparams.qangle.theta = theta;
    iparams.qangle.c[0]  = c[0];
    iparams.qangle.c[1]  = c[1];
    iparams.qangle.c[2]  = c[2];
    iparams.qangle.c[3]  = c[3];
    iparams.qangle.c[4]  = c[4];
    return *this;
}

iListInput iListInput::setPDihedrals(InteractionFunction ft, real phiA, real cpA, int mult, real phiB, real cpB)
{
    ftype              = ft;
    iparams.pdihs.phiA = phiA;
    iparams.pdihs.cpA  = cpA;
    iparams.pdihs.phiB = phiB;
    iparams.pdihs.cpB  = cpB;
    iparams.pdihs.mult = mult;
    fep                = (phiA != phiB || cpA != cpB);
    return *this;
}

iListInput iListInput::setPDihedrals(InteractionFunction ft, real phiA, real cpA, int mult)
{
    return setPDihedrals(ft, phiA, cpA, mult, phiA, cpA);
}

iListInput iListInput::setRbDihedrals(const real rbcA[NR_RBDIHS], const real rbcB[NR_RBDIHS])
{
    ftype = InteractionFunction::RyckaertBellemansDihedrals;
    fep   = false;
    for (int i = 0; i < NR_RBDIHS; i++)
    {
        iparams.rbdihs.rbcA[i] = rbcA[i];
        iparams.rbdihs.rbcB[i] = rbcB[i];
        fep                    = fep || (rbcA[i] != rbcB[i]);
    }
    return *this;
}

iListInput iListInput::setRbDihedrals(const real rbc[NR_RBDIHS])
{
    return setRbDihedrals(rbc, rbc);
}

iListInput iListInput::setPolarization(real alpha)
{
    ftype                  = InteractionFunction::Polarization;
    fep                    = false;
    iparams.polarize.alpha = alpha;
    return *this;
}

iListInput iListInput::setAnharmPolarization(real alpha, real drcut, real khyp)
{
    ftype                         = InteractionFunction::AnharmonicPolarization;
    fep                           = false;
    iparams.anharm_polarize.alpha = alpha;
    iparams.anharm_polarize.drcut = drcut;
    iparams.anharm_polarize.khyp  = khyp;
    return *this;
}

iListInput iListInput::setTholePolarization(real a, real alpha1, real alpha2)
{
    ftype                = InteractionFunction::TholePolarization;
    fep                  = false;
    iparams.thole.a      = a;
    iparams.thole.alpha1 = alpha1;
    iparams.thole.alpha2 = alpha2;
    return *this;
}

iListInput iListInput::setWaterPolarization(real alpha_x, real alpha_y, real alpha_z, real rOH, real rHH, real rOD)
{
    ftype             = InteractionFunction::WaterPolarization;
    fep               = false;
    iparams.wpol.al_x = alpha_x;
    iparams.wpol.al_y = alpha_y;
    iparams.wpol.al_z = alpha_z;
    iparams.wpol.rOH  = rOH;
    iparams.wpol.rHH  = rHH;
    iparams.wpol.rOD  = rOD;
    return *this;
}

std::ostream& operator<<(std::ostream& out, const iListInput& input)
{
    using std::endl;
    const InteractionFunction ftype = input.ftype.value();
    out << "Function type " << static_cast<int>(ftype) << " called " << interaction_function[ftype].name
        << " ie. labelled '" << interaction_function[ftype].longname << "' in an energy file" << endl;

    StringOutputStream stream;
    {
        TextWriter writer(&stream);
        printInteractionParameters(&writer, input.ftype.value(), input.iparams);
    }
    out << "Function parameters " << stream.toString();
    out << "Parameters trigger FEP? " << (input.fep ? "true" : "false") << endl;
    return out;
}

void fillIatoms(std::optional<InteractionFunction> ftype, std::vector<t_iatom>* iatoms)
{
    std::unordered_map<int, std::vector<int>> ia   = { { 2, { 0, 0, 1, 0, 1, 2, 0, 2, 3 } },
                                                       { 3, { 0, 0, 1, 2, 0, 1, 2, 3 } },
                                                       { 4, { 0, 0, 1, 2, 3 } },
                                                       { 5, { 0, 0, 1, 2, 3, 0 } } };
    int                                       nral = interaction_function[ftype.value()].nratoms;
    for (auto& i : ia[nral])
    {
        iatoms->push_back(i);
    }
}

void checkOutput(TestReferenceChecker* checker, const OutputQuantities& output, BondedKernelFlavor bondedKernelFlavor)
{
    if (computeEnergy(bondedKernelFlavor))
    {
        checker->checkReal(output.energy, "Epot ");
        checker->checkReal(output.dvdlambda, "dVdlambda ");
    }
    checker->checkSequence(std::begin(output.f), std::end(output.f), "Forces");
}

} // namespace test
} // namespace gmx
