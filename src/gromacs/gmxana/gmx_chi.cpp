/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

static gmx_bool bAllowed(real phi, real psi)
{
    static const char* map[] = { "1100000000000000001111111000000000001111111111111111111111111",
                                 "1100000000000000001111110000000000011111111111111111111111111",
                                 "1100000000000000001111110000000000011111111111111111111111111",
                                 "1100000000000000001111100000000000111111111111111111111111111",
                                 "1100000000000000001111100000000000111111111111111111111111111",
                                 "1100000000000000001111100000000001111111111111111111111111111",
                                 "1100000000000000001111100000000001111111111111111111111111111",
                                 "1100000000000000001111100000000011111111111111111111111111111",
                                 "1110000000000000001111110000000111111111111111111111111111111",
                                 "1110000000000000001111110000001111111111111111111111111111111",
                                 "1110000000000000001111111000011111111111111111111111111111111",
                                 "1110000000000000001111111100111111111111111111111111111111111",
                                 "1110000000000000001111111111111111111111111111111111111111111",
                                 "1110000000000000001111111111111111111111111111111111111111111",
                                 "1110000000000000001111111111111111111111111111111111111111111",
                                 "1110000000000000001111111111111111111111111111111111111111111",
                                 "1110000000000000001111111111111110011111111111111111111111111",
                                 "1110000000000000001111111111111100000111111111111111111111111",
                                 "1110000000000000001111111111111000000000001111111111111111111",
                                 "1100000000000000001111111111110000000000000011111111111111111",
                                 "1100000000000000001111111111100000000000000011111111111111111",
                                 "1000000000000000001111111111000000000000000001111111111111110",
                                 "0000000000000000001111111110000000000000000000111111111111100",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000111111111111000000000000000",
                                 "1100000000000000000000000000000001111111111111100000000000111",
                                 "1100000000000000000000000000000001111111111111110000000000111",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000",
                                 "0000000000000000000000000000000000000000000000000000000000000" };
    int                x, y;

#define INDEX(ppp) (((static_cast<int>(360 + (ppp)*gmx::c_rad2Deg)) % 360) / 6)
    x = INDEX(phi);
    y = INDEX(psi);
#undef INDEX

    return map[x][y] == '1';
}

static std::vector<int> make_chi_ind(gmx::ArrayRef<t_dlist> dlist)
{
    /* There are dlist.size() residues with max edMax dihedrals with 4
     * atoms each.  This may be an over-allocation, which is reduced
     * later. */
    std::vector<int> id(dlist.size() * edMax * 4);

    int n = 0;
    for (auto& dihedral : dlist)
    {
        /* Phi, fake for N-terminal residue */
        dihedral.j0[edPhi] = n / 4;
        if (dihedral.atm.minC >= 0)
        {
            id[n++] = dihedral.atm.minC;
        }
        else
        {
            id[n++] = dihedral.atm.H;
        }
        id[n++] = dihedral.atm.N;
        id[n++] = dihedral.atm.Cn[1];
        id[n++] = dihedral.atm.C;
    }
    for (auto& dihedral : dlist)
    {
        /* Psi, fake for C-terminal residue */
        dihedral.j0[edPsi] = n / 4;
        id[n++]            = dihedral.atm.N;
        id[n++]            = dihedral.atm.Cn[1];
        id[n++]            = dihedral.atm.C;
        if (dihedral.atm.maxN >= 0)
        {
            id[n++] = dihedral.atm.maxN;
        }
        else
        {
            id[n++] = dihedral.atm.O;
        }
    }
    for (auto& dihedral : dlist)
    {
        /* Omega */
        if (has_dihedral(edOmega, dihedral))
        {
            dihedral.j0[edOmega] = n / 4;
            id[n++]              = dihedral.atm.minCalpha;
            id[n++]              = dihedral.atm.minC;
            id[n++]              = dihedral.atm.N;
            id[n++]              = dihedral.atm.Cn[1];
        }
    }
    for (int Xi = 0; (Xi < MAXCHI); Xi++)
    {
        /* Chi# */
        for (auto& dihedral : dlist)
        {
            if (dihedral.atm.Cn[Xi + 3] != -1)
            {
                dihedral.j0[edChi1 + Xi] = n / 4;
                id[n++]                  = dihedral.atm.Cn[Xi];
                id[n++]                  = dihedral.atm.Cn[Xi + 1];
                id[n++]                  = dihedral.atm.Cn[Xi + 2];
                id[n++]                  = dihedral.atm.Cn[Xi + 3];
            }
        }
    }
    id.resize(n);

    return id;
}

static void do_dihcorr(const char*                  fn,
                       int                          nf,
                       int                          ndih,
                       real**                       dih,
                       real                         dt,
                       gmx::ArrayRef<const t_dlist> dlist,
                       real                         time[],
                       int                          maxchi,
                       gmx_bool                     bPhi,
                       gmx_bool                     bPsi,
                       gmx_bool                     bChi,
                       gmx_bool                     bOmega,
                       const gmx_output_env_t*      oenv)
{
    char name1[256], name2[256];
    int  j, Xi;

    do_autocorr(fn, oenv, "Dihedral Autocorrelation Function", nf, ndih, dih, dt, eacCos, FALSE);
    /* Dump em all */
    j = 0;
    for (const auto& dihedral : dlist)
    {
        if (bPhi)
        {
            print_one(oenv, "corrphi", dihedral.name, "Phi ACF for", "C(t)", nf / 2, time, dih[j]);
        }
        j++;
    }
    for (const auto& dihedral : dlist)
    {
        if (bPsi)
        {
            print_one(oenv, "corrpsi", dihedral.name, "Psi ACF for", "C(t)", nf / 2, time, dih[j]);
        }
        j++;
    }
    for (const auto& dihedral : dlist)
    {
        if (has_dihedral(edOmega, dihedral))
        {
            if (bOmega)
            {
                print_one(oenv, "corromega", dihedral.name, "Omega ACF for", "C(t)", nf / 2, time, dih[j]);
            }
            j++;
        }
    }
    for (Xi = 0; (Xi < maxchi); Xi++)
    {
        sprintf(name1, "corrchi%d", Xi + 1);
        sprintf(name2, "Chi%d ACF for", Xi + 1);
        for (const auto& dihedral : dlist)
        {
            if (dihedral.atm.Cn[Xi + 3] != -1)
            {
                if (bChi)
                {
                    print_one(oenv, name1, dihedral.name, name2, "C(t)", nf / 2, time, dih[j]);
                }
                j++;
            }
        }
    }
    fprintf(stderr, "\n");
}

static void copy_dih_data(const real in[], real out[], int nf, gmx_bool bLEAVE)
{
    /* if bLEAVE, do nothing to data in copying to out
     * otherwise multiply by 180/pi to convert rad to deg */
    int  i;
    real mult;
    if (bLEAVE)
    {
        mult = 1;
    }
    else
    {
        mult = (180.0 / M_PI);
    }
    for (i = 0; (i < nf); i++)
    {
        out[i] = in[i] * mult;
    }
}

static void dump_em_all(gmx::ArrayRef<const t_dlist> dlist,
                        int                          nf,
                        real                         time[],
                        real**                       dih,
                        int                          maxchi,
                        gmx_bool                     bPhi,
                        gmx_bool                     bPsi,
                        gmx_bool                     bChi,
                        gmx_bool                     bOmega,
                        gmx_bool                     bRAD,
                        const gmx_output_env_t*      oenv)
{
    char  name[256], titlestr[256], ystr[256];
    real* data;

    snew(data, nf);
    gmx::sfree_guard dataGuard(data);
    if (bRAD)
    {
        std::strcpy(ystr, "Angle (rad)");
    }
    else
    {
        std::strcpy(ystr, "Angle (degrees)");
    }

    /* Dump em all */
    int j = 0;
    for (const auto& dihedral : dlist)
    {
        if (bPhi)
        {
            copy_dih_data(dih[j], data, nf, bRAD);
            print_one(oenv, "phi", dihedral.name, "\\xf\\f{}", ystr, nf, time, data);
        }
        j++;
    }
    for (const auto& dihedral : dlist)
    {
        if (bPsi)
        {
            copy_dih_data(dih[j], data, nf, bRAD);
            print_one(oenv, "psi", dihedral.name, "\\xy\\f{}", ystr, nf, time, data);
        }
        j++;
    }
    for (const auto& dihedral : dlist)
    {
        if (has_dihedral(edOmega, dihedral))
        {
            if (bOmega)
            {
                copy_dih_data(dih[j], data, nf, bRAD);
                print_one(oenv, "omega", dihedral.name, "\\xw\\f{}", ystr, nf, time, data);
            }
            j++;
        }
    }

    for (int Xi = 0; (Xi < maxchi); Xi++)
    {
        for (const auto& dihedral : dlist)
        {
            if (dihedral.atm.Cn[Xi + 3] != -1)
            {
                if (bChi)
                {
                    sprintf(name, "chi%d", Xi + 1);
                    sprintf(titlestr, "\\xc\\f{}\\s%d\\N", Xi + 1);
                    copy_dih_data(dih[j], data, nf, bRAD);
                    print_one(oenv, name, dihedral.name, titlestr, ystr, nf, time, data);
                }
                j++;
            }
        }
    }
    fprintf(stderr, "\n");
}

static void reset_one(real dih[], int nf, real phase)
{
    int j;

    for (j = 0; (j < nf); j++)
    {
        dih[j] += phase;
        while (dih[j] < -M_PI)
        {
            dih[j] += 2 * M_PI;
        }
        while (dih[j] >= M_PI)
        {
            dih[j] -= 2 * M_PI;
        }
    }
}

static int reset_em_all(gmx::ArrayRef<const t_dlist> dlist, int nf, real** dih, int maxchi)
{
    /* Reset em all */
    int j = 0;
    /* Phi */
    for (const auto& dihedral : dlist)
    {
        if (dihedral.atm.minC == -1)
        {
            reset_one(dih[j++], nf, M_PI);
        }
        else
        {
            reset_one(dih[j++], nf, 0);
        }
    }
    /* Psi */
    for (const auto& dihedral : dlist)
    {
        if (dihedral.atm.maxN == -1)
        {
            reset_one(dih[j++], nf, M_PI);
        }
        else
        {
            reset_one(dih[j++], nf, 0);
        }
    }
    /* Omega */
    for (const auto& dihedral : dlist)
    {
        if (has_dihedral(edOmega, dihedral))
        {
            reset_one(dih[j++], nf, 0);
        }
    }
    /* Chi 1 thru maxchi */
    for (int Xi = 0; (Xi < maxchi); Xi++)
    {
        for (const auto& dihedral : dlist)
        {
            if (dihedral.atm.Cn[Xi + 3] != -1)
            {
                reset_one(dih[j], nf, 0);
                j++;
            }
        }
    }
    fprintf(stderr, "j after resetting (nr. active dihedrals) = %d\n", j);
    return j;
}

static void histogramming(FILE*                    log,
                          int                      nbin,
                          int                      nf,
                          int                      maxchi,
                          real**                   dih,
                          gmx::ArrayRef<t_dlist>   dlist,
                          gmx::ArrayRef<const int> index,
                          gmx_bool                 bPhi,
                          gmx_bool                 bPsi,
                          gmx_bool                 bOmega,
                          gmx_bool                 bChi,
                          gmx_bool                 bNormalize,
                          gmx_bool                 bSSHisto,
                          const char*              ssdump,
                          real                     bfac_max,
                          const t_atoms*           atoms,
                          gmx_bool                 bDo_jc,
                          const char*              fn,
                          const gmx_output_env_t*  oenv)
{
    /* also gets 3J couplings and order parameters S2 */
    // Avoid warnings about narrowing conversions from double to real
#ifdef _MSC_VER
#    pragma warning(disable : 4838)
#endif
    t_karplus kkkphi[]  = { { "J_NHa1", 6.51, -1.76, 1.6, -M_PI / 3, 0.0, 0.0 },
                           { "J_NHa2", 6.51, -1.76, 1.6, M_PI / 3, 0.0, 0.0 },
                           { "J_HaC'", 4.0, 1.1, 0.1, 0.0, 0.0, 0.0 },
                           { "J_NHCb", 4.7, -1.5, -0.2, M_PI / 3, 0.0, 0.0 },
                           { "J_Ci-1Hai", 4.5, -1.3, -1.2, 2 * M_PI / 3, 0.0, 0.0 } };
    t_karplus kkkpsi[]  = { { "J_HaN", -0.88, -0.61, -0.27, M_PI / 3, 0.0, 0.0 } };
    t_karplus kkkchi1[] = { { "JHaHb2", 9.5, -1.6, 1.8, -M_PI / 3, 0, 0.0 },
                            { "JHaHb3", 9.5, -1.6, 1.8, 0, 0, 0.0 } };
#ifdef _MSC_VER
#    pragma warning(default : 4838)
#endif
#define NKKKPHI asize(kkkphi)
#define NKKKPSI asize(kkkpsi)
#define NKKKCHI asize(kkkchi1)
#define NJC (NKKKPHI + NKKKPSI + NKKKCHI)

    FILE *                   fp, *ssfp[3] = { nullptr, nullptr, nullptr };
    const char*              sss[3] = { "sheet", "helix", "coil" };
    real                     S2;
    real **                  Jc, **Jcsig;
    int*                     histmp;
    int                      m, n, nn, nres, hindex, angle;
    gmx_bool                 bBfac, bOccup;
    char                     hisfile[256], hhisfile[256], title[256], *ss_str = nullptr;
    std::vector<std::string> leg;

    if (bSSHisto)
    {
        fp = gmx_ffopen(ssdump, "r");
        if (1 != fscanf(fp, "%d", &nres))
        {
            gmx_fatal(FARGS, "Error reading from file %s", ssdump);
        }

        snew(ss_str, nres + 1);
        if (1 != fscanf(fp, "%s", ss_str))
        {
            gmx_fatal(FARGS, "Error reading from file %s", ssdump);
        }

        gmx_ffclose(fp);
    }

    // Build a list of unique residue names found in the dihedral
    // list, so we can loop over those unique names conveniently.
    std::unordered_set<std::string> uniqueResidueNames;
    for (const auto& dihedral : dlist)
    {
        uniqueResidueNames.emplace(dihedral.residueName);
    }
    // Build the lookup tables for data relating to the all dihedrals
    // from each unique residue name represented in the dihedral list.
    std::array<std::map<std::string, std::vector<std::vector<int>>>, 3> his_aa_ss;
    std::vector<std::map<std::string, std::vector<int>>>                his_aa(edMax);
    for (const auto& residueName : uniqueResidueNames)
    {
        if (bSSHisto)
        {
            for (auto& secondaryStructure : his_aa_ss)
            {
                secondaryStructure[residueName] =
                        std::vector<std::vector<int>>(edMax, std::vector<int>(nbin, 0));
            }
        }
        for (auto& dihedraltype : his_aa)
        {
            dihedraltype[residueName] = std::vector<int>(nbin, 0);
        }
    }
    snew(histmp, nbin);

    snew(Jc, dlist.size());
    snew(Jcsig, dlist.size());
    for (size_t i = 0; i < dlist.size(); i++)
    {
        snew(Jc[i], NJC);
        snew(Jcsig[i], NJC);
    }

    int j = 0;
    n     = 0;
    for (int Dih = 0; (Dih < NONCHI + maxchi); Dih++)
    {
        int i = 0;
        for (auto& dihedral : dlist)
        {
            if (((Dih < edOmega)) || ((Dih == edOmega) && (has_dihedral(edOmega, dihedral)))
                || ((Dih > edOmega) && (dihedral.atm.Cn[Dih - NONCHI + 3] != -1)))
            {
                make_histo(log, nf, dih[j], nbin, histmp, -M_PI, M_PI);

                if (bSSHisto)
                {
                    /* Assume there is only one structure, the first.
                     * Compute index in histogram.
                     */
                    /* Check the atoms to see whether their B-factors are low enough
                     * Check atoms to see their occupancy is 1.
                     */
                    bBfac = bOccup = TRUE;
                    for (nn = 0; (nn < 4); nn++, n++)
                    {
                        bBfac  = bBfac && (atoms->pdbinfo[index[n]].bfac <= bfac_max);
                        bOccup = bOccup && (atoms->pdbinfo[index[n]].occup == 1);
                    }
                    if (bOccup && ((bfac_max <= 0) || bBfac))
                    {
                        hindex = static_cast<int>(((dih[j][0] + M_PI) * nbin) / (2 * M_PI));
                        range_check(hindex, 0, nbin);

                        /* Assign dihedral to either of the structure determined
                         * histograms
                         */
                        switch (ss_str[dihedral.resnr])
                        {
                            case 'E': his_aa_ss[0][dihedral.residueName][Dih][hindex]++; break;
                            case 'H': his_aa_ss[1][dihedral.residueName][Dih][hindex]++; break;
                            default: his_aa_ss[2][dihedral.residueName][Dih][hindex]++; break;
                        }
                    }
                    else if (debug)
                    {
                        fprintf(debug, "Res. %d has incomplete occupancy or bfacs > %g\n", dihedral.resnr, bfac_max);
                    }
                }
                else
                {
                    n += 4;
                }

                switch (Dih)
                {
                    case edPhi:
                        calc_distribution_props(nbin, histmp, -M_PI, NKKKPHI, kkkphi, &S2);

                        for (m = 0; (m < NKKKPHI); m++)
                        {
                            Jc[i][m]    = kkkphi[m].Jc;
                            Jcsig[i][m] = kkkphi[m].Jcsig;
                        }
                        break;
                    case edPsi:
                        calc_distribution_props(nbin, histmp, -M_PI, NKKKPSI, kkkpsi, &S2);

                        for (m = 0; (m < NKKKPSI); m++)
                        {
                            Jc[i][NKKKPHI + m]    = kkkpsi[m].Jc;
                            Jcsig[i][NKKKPHI + m] = kkkpsi[m].Jcsig;
                        }
                        break;
                    case edChi1:
                        calc_distribution_props(nbin, histmp, -M_PI, NKKKCHI, kkkchi1, &S2);
                        for (m = 0; (m < NKKKCHI); m++)
                        {
                            Jc[i][NKKKPHI + NKKKPSI + m]    = kkkchi1[m].Jc;
                            Jcsig[i][NKKKPHI + NKKKPSI + m] = kkkchi1[m].Jcsig;
                        }
                        break;
                    default: /* covers edOmega and higher Chis than Chi1 */
                        calc_distribution_props(nbin, histmp, -M_PI, 0, nullptr, &S2);
                        break;
                }
                dihedral.S2[Dih] = S2;

                /* Sum distribution per amino acid type as well */
                for (int k = 0; (k < nbin); k++)
                {
                    his_aa[Dih][dihedral.residueName][k] += histmp[k];
                    histmp[k] = 0;
                }
                j++;
            }
            else /* dihed not defined */
            {
                dihedral.S2[Dih] = 0.0;
            }
            ++i;
        }
    }
    sfree(histmp);

    /* Print out Jcouplings */
    fprintf(log, "\n *** J-Couplings from simulation (plus std. dev.) ***\n\n");
    fprintf(log, "Residue   ");
    for (int i = 0; (i < NKKKPHI); i++)
    {
        fprintf(log, "%7s   SD", kkkphi[i].name);
    }
    for (int i = 0; (i < NKKKPSI); i++)
    {
        fprintf(log, "%7s   SD", kkkpsi[i].name);
    }
    for (int i = 0; (i < NKKKCHI); i++)
    {
        fprintf(log, "%7s   SD", kkkchi1[i].name);
    }
    fprintf(log, "\n");
    for (int i = 0; (i < NJC + 1); i++)
    {
        fprintf(log, "------------");
    }
    fprintf(log, "\n");
    {
        int i = 0;
        for (const auto& dihedral : dlist)
        {
            fprintf(log, "%-10s", dihedral.name);
            for (int j = 0; (j < NJC); j++)
            {
                fprintf(log, "  %5.2f %4.2f", Jc[i][j], Jcsig[i][j]);
            }
            fprintf(log, "\n");
            ++i;
        }
        fprintf(log, "\n");
    }

    /* and to -jc file... */
    if (bDo_jc)
    {
        fp = xvgropen(fn, "\\S3\\NJ-Couplings from Karplus Equation", "Residue", "Coupling", oenv);
        for (int i = 0; (i < NKKKPHI); i++)
        {
            leg.emplace_back(kkkphi[i].name);
        }
        for (int i = 0; (i < NKKKPSI); i++)
        {
            leg.emplace_back(kkkpsi[i].name);
        }
        for (int i = 0; (i < NKKKCHI); i++)
        {
            leg.emplace_back(kkkchi1[i].name);
        }
        xvgrLegend(fp, leg, oenv);
        fprintf(fp, "%5s ", "#Res.");
        for (int i = 0; (i < NJC); i++)
        {
            fprintf(fp, "%10s ", leg[i].c_str());
        }
        fprintf(fp, "\n");
        {
            int i = 0;
            for (const auto& dihedral : dlist)
            {
                fprintf(fp, "%5d ", dihedral.resnr);
                for (int j = 0; (j < NJC); j++)
                {
                    fprintf(fp, "  %8.3f", Jc[i][j]);
                }
                fprintf(fp, "\n");
                ++i;
            }
        }
        xvgrclose(fp);
    }
    /* finished -jc stuff */

    std::vector<real> normhisto(nbin);
    for (const auto& residueName : uniqueResidueNames)
    {
        for (int Dih = 0; (Dih < edMax); Dih++)
        {
            /* First check whether something is in there */
            int j;
            for (j = 0; (j < nbin); j++)
            {
                if (his_aa[Dih][residueName][j] != 0)
                {
                    break;
                }
            }
            if ((j < nbin)
                && ((bPhi && (Dih == edPhi)) || (bPsi && (Dih == edPsi))
                    || (bOmega && (Dih == edOmega)) || (bChi && (Dih >= edChi1))))
            {
                if (bNormalize)
                {
                    normalize_histo(his_aa[Dih][residueName], (360.0 / nbin), normhisto);
                }

                switch (Dih)
                {
                    case edPhi:
                        sprintf(hisfile, "histo-phi%s", residueName.c_str());
                        sprintf(title, "\\xf\\f{} Distribution for %s", residueName.c_str());
                        break;
                    case edPsi:
                        sprintf(hisfile, "histo-psi%s", residueName.c_str());
                        sprintf(title, "\\xy\\f{} Distribution for %s", residueName.c_str());
                        break;
                    case edOmega:
                        sprintf(hisfile, "histo-omega%s", residueName.c_str());
                        sprintf(title, "\\xw\\f{} Distribution for %s", residueName.c_str());
                        break;
                    default:
                        sprintf(hisfile, "histo-chi%d%s", Dih - NONCHI + 1, residueName.c_str());
                        sprintf(title,
                                "\\xc\\f{}\\s%d\\N Distribution for %s",
                                Dih - NONCHI + 1,
                                residueName.c_str());
                }
                std::strcpy(hhisfile, hisfile);
                std::strcat(hhisfile, ".xvg");
                fp = xvgropen(hhisfile, title, "Degrees", "", oenv);
                if (output_env_get_print_xvgr_codes(oenv))
                {
                    fprintf(fp, "@ with g0\n");
                }
                xvgr_world(fp, -180, 0, 180, 0.1, oenv);
                if (output_env_get_print_xvgr_codes(oenv))
                {
                    fprintf(fp,
                            "# this effort to set graph size fails unless you run with -autoscale "
                            "none or -autoscale y flags\n");
                    fprintf(fp, "@ xaxis tick on\n");
                    fprintf(fp, "@ xaxis tick major 90\n");
                    fprintf(fp, "@ xaxis tick minor 30\n");
                    fprintf(fp, "@ xaxis ticklabel prec 0\n");
                    fprintf(fp, "@ yaxis tick off\n");
                    fprintf(fp, "@ yaxis ticklabel off\n");
                    fprintf(fp, "@ type xy\n");
                }
                if (bSSHisto)
                {
                    for (int k = 0; (k < 3); k++)
                    {
                        std::string sshisfile = gmx::formatString("%s-%s.xvg", hisfile, sss[k]);
                        ssfp[k]               = gmx_ffopen(sshisfile, "w");
                    }
                }
                for (int j = 0; (j < nbin); j++)
                {
                    angle = -180 + (360 / nbin) * j;
                    if (bNormalize)
                    {
                        fprintf(fp, "%5d  %10g\n", angle, normhisto[j]);
                    }
                    else
                    {
                        fprintf(fp, "%5d  %10d\n", angle, his_aa[Dih][residueName][j]);
                    }
                    if (bSSHisto)
                    {
                        for (int k = 0; (k < 3); k++)
                        {
                            fprintf(ssfp[k], "%5d  %10d\n", angle, his_aa_ss[k][residueName][Dih][j]);
                        }
                    }
                }
                fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
                xvgrclose(fp);
                if (bSSHisto)
                {
                    for (int k = 0; (k < 3); k++)
                    {
                        fprintf(ssfp[k], "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
                        gmx_ffclose(ssfp[k]);
                    }
                }
            }
        }
    }

    if (bSSHisto)
    {
        sfree(ss_str);
    }
    for (size_t i = 0; i < dlist.size(); i++)
    {
        sfree(Jc[i]);
        sfree(Jcsig[i]);
    }
    sfree(Jc);
    sfree(Jcsig);
}

static FILE* rama_file(const char* fn, const char* title, const char* xaxis, const char* yaxis, const gmx_output_env_t* oenv)
{
    FILE* fp;

    fp = xvgropen(fn, title, xaxis, yaxis, oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ with g0\n");
    }
    xvgr_world(fp, -180, -180, 180, 180, oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ xaxis tick on\n");
        fprintf(fp, "@ xaxis tick major 90\n");
        fprintf(fp, "@ xaxis tick minor 30\n");
        fprintf(fp, "@ xaxis ticklabel prec 0\n");
        fprintf(fp, "@ yaxis tick on\n");
        fprintf(fp, "@ yaxis tick major 90\n");
        fprintf(fp, "@ yaxis tick minor 30\n");
        fprintf(fp, "@ yaxis ticklabel prec 0\n");
        fprintf(fp, "@    s0 type xy\n");
        fprintf(fp, "@    s0 symbol 2\n");
        fprintf(fp, "@    s0 symbol size 0.410000\n");
        fprintf(fp, "@    s0 symbol fill 1\n");
        fprintf(fp, "@    s0 symbol color 1\n");
        fprintf(fp, "@    s0 symbol linewidth 1\n");
        fprintf(fp, "@    s0 symbol linestyle 1\n");
        fprintf(fp, "@    s0 symbol center false\n");
        fprintf(fp, "@    s0 symbol char 0\n");
        fprintf(fp, "@    s0 skip 0\n");
        fprintf(fp, "@    s0 linestyle 0\n");
        fprintf(fp, "@    s0 linewidth 1\n");
        fprintf(fp, "@ type xy\n");
    }
    return fp;
}

static void do_rama(int                          nf,
                    gmx::ArrayRef<const t_dlist> dlist,
                    real**                       dih,
                    gmx_bool                     bViol,
                    gmx_bool                     bRamOmega,
                    const gmx_output_env_t*      oenv)
{
    FILE *        fp, *gp = nullptr;
    gmx_bool      bOm;
    char          fn[256];
    int           Xi1, Xi2, Phi, Psi, Om = 0, nlevels;
    constexpr int NMAT = 120;
    real **       mat  = nullptr, phi, psi, omega, axis[NMAT], lo, hi;
    t_rgb         rlo  = { 1.0, 0.0, 0.0 };
    t_rgb         rmid = { 1.0, 1.0, 1.0 };
    t_rgb         rhi  = { 0.0, 0.0, 1.0 };

    for (const auto& dihedral : dlist)
    {
        if ((has_dihedral(edPhi, dihedral)) && has_dihedral(edPsi, dihedral))
        {
            sprintf(fn, "ramaPhiPsi%s.xvg", dihedral.name);
            fp  = rama_file(fn, "Ramachandran Plot", "\\8f\\4 (deg)", "\\8y\\4 (deg)", oenv);
            bOm = bRamOmega && has_dihedral(edOmega, dihedral);
            if (bOm)
            {
                Om = dihedral.j0[edOmega];
                snew(mat, NMAT);
                for (int j = 0; (j < NMAT); j++)
                {
                    snew(mat[j], NMAT);
                    axis[j] = -180 + gmx::exactDiv(360 * j, NMAT);
                }
            }
            if (bViol)
            {
                sprintf(fn, "violPhiPsi%s.xvg", dihedral.name);
                gp = gmx_ffopen(fn, "w");
            }
            Phi = dihedral.j0[edPhi];
            Psi = dihedral.j0[edPsi];
            for (int j = 0; (j < nf); j++)
            {
                phi = gmx::c_rad2Deg * dih[Phi][j];
                psi = gmx::c_rad2Deg * dih[Psi][j];
                fprintf(fp, "%10g  %10g\n", phi, psi);
                if (bViol)
                {
                    fprintf(gp,
                            "%d\n",
                            static_cast<int>(!bAllowed(dih[Phi][j], gmx::c_rad2Deg * dih[Psi][j])));
                }
                if (bOm)
                {
                    omega = gmx::c_rad2Deg * dih[Om][j];
                    mat[static_cast<int>(((phi * NMAT) / 360) + gmx::exactDiv(NMAT, 2))]
                       [static_cast<int>(((psi * NMAT) / 360) + gmx::exactDiv(NMAT, 2))] += omega;
                }
            }
            if (bViol)
            {
                gmx_ffclose(gp);
            }
            xvgrclose(fp);
            if (bOm)
            {
                sprintf(fn, "ramomega%s.xpm", dihedral.name);
                fp = gmx_ffopen(fn, "w");
                lo = hi = 0;
                for (int j = 0; (j < NMAT); j++)
                {
                    for (int k = 0; (k < NMAT); k++)
                    {
                        mat[j][k] /= nf;
                        lo = std::min(mat[j][k], lo);
                        hi = std::max(mat[j][k], hi);
                    }
                }
                /* Symmetrise */
                if (std::abs(lo) > std::abs(hi))
                {
                    hi = -lo;
                }
                else
                {
                    lo = -hi;
                }
                /* Add 180 */
                for (int j = 0; (j < NMAT); j++)
                {
                    for (int k = 0; (k < NMAT); k++)
                    {
                        mat[j][k] += 180;
                    }
                }
                lo += 180;
                hi += 180;
                nlevels = 20;
                write_xpm3(fp,
                           0,
                           "Omega/Ramachandran Plot",
                           "Deg",
                           "Phi",
                           "Psi",
                           NMAT,
                           NMAT,
                           axis,
                           axis,
                           mat,
                           lo,
                           180.0,
                           hi,
                           rlo,
                           rmid,
                           rhi,
                           &nlevels);
                gmx_ffclose(fp);
                for (int j = 0; (j < NMAT); j++)
                {
                    sfree(mat[j]);
                }
                sfree(mat);
            }
        }
        if (has_dihedral(edChi1, dihedral) && has_dihedral(edChi2, dihedral))
        {
            sprintf(fn, "ramaX1X2%s.xvg", dihedral.name);
            fp  = rama_file(fn,
                           "\\8c\\4\\s1\\N-\\8c\\4\\s2\\N Ramachandran Plot",
                           "\\8c\\4\\s1\\N (deg)",
                           "\\8c\\4\\s2\\N (deg)",
                           oenv);
            Xi1 = dihedral.j0[edChi1];
            Xi2 = dihedral.j0[edChi2];
            for (int j = 0; (j < nf); j++)
            {
                fprintf(fp, "%10g  %10g\n", gmx::c_rad2Deg * dih[Xi1][j], gmx::c_rad2Deg * dih[Xi2][j]);
            }
            xvgrclose(fp);
        }
        else
        {
            fprintf(stderr, "No chi1 & chi2 angle for %s\n", dihedral.name);
        }
    }
}


static void print_transitions(const char*                  fn,
                              int                          maxchi,
                              gmx::ArrayRef<const t_dlist> dlist,
                              real                         dt,
                              const gmx_output_env_t*      oenv)
{
    /* based on order_params below */
    FILE* fp;

    /*  must correspond with enum in gstat.h */
    std::array<std::string, 9> leg = {
        "Phi", "Psi", "Omega", "Chi1", "Chi2", "Chi3", "Chi4", "Chi5", "Chi6",
    };
    /* Print order parameters */
    fp = xvgropen(fn, "Dihedral Rotamer Transitions", "Residue", "Transitions/ns", oenv);
    xvgrLegend(fp, gmx::makeArrayRef(leg).subArray(0, NONCHI + maxchi), oenv);

    fprintf(fp, "%5s ", "#Res.");
    fprintf(fp, "%10s %10s %10s ", leg[edPhi].c_str(), leg[edPsi].c_str(), leg[edOmega].c_str());
    for (int Xi = 0; Xi < maxchi; Xi++)
    {
        fprintf(fp, "%10s ", leg[NONCHI + Xi].c_str());
    }
    fprintf(fp, "\n");

    for (const auto& dihedral : dlist)
    {
        fprintf(fp, "%5d ", dihedral.resnr);
        for (int Dih = 0; (Dih < NONCHI + maxchi); Dih++)
        {
            fprintf(fp, "%10.3f ", dihedral.ntr[Dih] / dt);
        }
        /* fprintf(fp,"%12s\n",dihedral.name);  this confuses xmgrace */
        fprintf(fp, "\n");
    }
    xvgrclose(fp);
}

static void order_params(FILE*                        log,
                         const char*                  fn,
                         int                          maxchi,
                         gmx::ArrayRef<const t_dlist> dlist,
                         const char*                  pdbfn,
                         real                         bfac_init,
                         t_atoms*                     atoms,
                         const rvec                   x[],
                         PbcType                      pbcType,
                         matrix                       box,
                         gmx_bool                     bPhi,
                         gmx_bool                     bPsi,
                         gmx_bool                     bChi,
                         const gmx_output_env_t*      oenv)
{
    FILE* fp;
    int   nh[edMax];
    real  S2Max, S2Min;

    /* except for S2Min/Max, must correspond with enum in pp2shift.h:38 */
    std::array<std::string, 11> const_leg = { "S2Min", "S2Max", "Phi",  "Psi",  "Omega", "Chi1",
                                              "Chi2",  "Chi3",  "Chi4", "Chi5", "Chi6" };
#define NLEG asize(leg)

    /* Print order parameters */
    fp = xvgropen(fn, "Dihedral Order Parameters", "Residue", "S2", oenv);
    xvgrLegend(fp, gmx::makeArrayRef(const_leg).subArray(0, 2 + NONCHI + maxchi), oenv);

    for (int Dih = 0; (Dih < edMax); Dih++)
    {
        nh[Dih] = 0;
    }

    fprintf(fp, "%5s ", "#Res.");
    fprintf(fp, "%10s %10s ", const_leg[0].c_str(), const_leg[1].c_str());
    fprintf(fp,
            "%10s %10s %10s ",
            const_leg[2 + edPhi].c_str(),
            const_leg[2 + edPsi].c_str(),
            const_leg[2 + edOmega].c_str());
    for (int Xi = 0; Xi < maxchi; Xi++)
    {
        fprintf(fp, "%10s ", const_leg[2 + NONCHI + Xi].c_str());
    }
    fprintf(fp, "\n");

    for (const auto& dihedral : dlist)
    {
        S2Max = -10;
        S2Min = 10;
        for (int Dih = 0; (Dih < NONCHI + maxchi); Dih++)
        {
            if (dihedral.S2[Dih] != 0)
            {
                if (dihedral.S2[Dih] > S2Max)
                {
                    S2Max = dihedral.S2[Dih];
                }
                if (dihedral.S2[Dih] < S2Min)
                {
                    S2Min = dihedral.S2[Dih];
                }
            }
            if (dihedral.S2[Dih] > 0.8)
            {
                nh[Dih]++;
            }
        }
        fprintf(fp, "%5d ", dihedral.resnr);
        fprintf(fp, "%10.3f %10.3f ", S2Min, S2Max);
        for (int Dih = 0; (Dih < NONCHI + maxchi); Dih++)
        {
            fprintf(fp, "%10.3f ", dihedral.S2[Dih]);
        }
        fprintf(fp, "\n");
        /* fprintf(fp,"%12s\n",dihedral.name);  this confuses xmgrace */
    }
    xvgrclose(fp);

    if (nullptr != pdbfn)
    {
        real x0, y0, z0;

        atoms->havePdbInfo = TRUE;

        if (nullptr == atoms->pdbinfo)
        {
            snew(atoms->pdbinfo, atoms->nr);
        }
        for (int i = 0; (i < atoms->nr); i++)
        {
            atoms->pdbinfo[i].bfac = bfac_init;
        }

        for (const auto& dihedral : dlist)
        {
            atoms->pdbinfo[dihedral.atm.N].bfac = -dihedral.S2[0]; /* Phi */
            atoms->pdbinfo[dihedral.atm.H].bfac = -dihedral.S2[0]; /* Phi */
            atoms->pdbinfo[dihedral.atm.C].bfac = -dihedral.S2[1]; /* Psi */
            atoms->pdbinfo[dihedral.atm.O].bfac = -dihedral.S2[1]; /* Psi */
            for (int Xi = 0; (Xi < maxchi); Xi++)                  /* Chi's */
            {
                if (dihedral.atm.Cn[Xi + 3] != -1)
                {
                    atoms->pdbinfo[dihedral.atm.Cn[Xi + 1]].bfac = -dihedral.S2[NONCHI + Xi];
                }
            }
        }

        fp = gmx_ffopen(pdbfn, "w");
        fprintf(fp, "REMARK generated by gmx chi\n");
        fprintf(fp,
                "REMARK "
                "B-factor field contains negative of dihedral order parameters\n");
        write_pdbfile(fp, nullptr, atoms, x, pbcType, box, ' ', 0, nullptr);
        x0 = y0 = z0 = 1000.0;
        for (int i = 0; (i < atoms->nr); i++)
        {
            x0 = std::min(x0, x[i][XX]);
            y0 = std::min(y0, x[i][YY]);
            z0 = std::min(z0, x[i][ZZ]);
        }
        x0 *= 10.0; /* nm -> angstrom */
        y0 *= 10.0; /* nm -> angstrom */
        z0 *= 10.0; /* nm -> angstrom */
        for (int i = 0; (i < 10); i++)
        {
            gmx_fprintf_pdb_atomline(fp,
                                     PdbRecordType::Atom,
                                     atoms->nr + 1 + i,
                                     "CA",
                                     ' ',
                                     "LEG",
                                     ' ',
                                     atoms->nres + 1,
                                     ' ',
                                     x0,
                                     y0,
                                     z0 + (1.2 * i),
                                     0.0,
                                     -0.1 * i,
                                     "");
        }
        gmx_ffclose(fp);
    }

    fprintf(log, "Dihedrals with S2 > 0.8\n");
    fprintf(log, "Dihedral: ");
    if (bPhi)
    {
        fprintf(log, " Phi  ");
    }
    if (bPsi)
    {
        fprintf(log, " Psi ");
    }
    if (bChi)
    {
        for (int Xi = 0; (Xi < maxchi); Xi++)
        {
            fprintf(log, " %s ", const_leg[2 + NONCHI + Xi].c_str());
        }
    }
    fprintf(log, "\nNumber:   ");
    if (bPhi)
    {
        fprintf(log, "%4d  ", nh[0]);
    }
    if (bPsi)
    {
        fprintf(log, "%4d  ", nh[1]);
    }
    if (bChi)
    {
        for (int Xi = 0; (Xi < maxchi); Xi++)
        {
            fprintf(log, "%4d  ", nh[NONCHI + Xi]);
        }
    }
    fprintf(log, "\n");
}

int gmx_chi(int argc, char* argv[])
{
    const char* desc[] = {
        "[THISMODULE] computes [GRK]phi[grk], [GRK]psi[grk], [GRK]omega[grk],",
        "and [GRK]chi[grk] dihedrals for all your",
        "amino acid backbone and sidechains.",
        "It can compute dihedral angle as a function of time, and as",
        "histogram distributions.",
        "The distributions [TT](histo-(dihedral)(RESIDUE).xvg[tt]) are cumulative over all ",
        "residues of each type.[PAR]",
        "If option [TT]-corr[tt] is given, the program will",
        "calculate dihedral autocorrelation functions. The function used",
        "is C(t) = [CHEVRON][COS][GRK]chi[grk]([GRK]tau[grk])[cos] ",
        "[COS][GRK]chi[grk]([GRK]tau[grk]+t)[cos][chevron]. The use of cosines",
        "rather than angles themselves, resolves the problem of periodicity.",
        "(Van der Spoel & Berendsen (1997), Biophys. J. 72, 2032-2041).",
        "Separate files for each dihedral of each residue",
        "[TT](corr(dihedral)(RESIDUE)(nresnr).xvg[tt]) are output, as well as a",
        "file containing the information for all residues (argument of [TT]-corr[tt]).[PAR]",
        "With option [TT]-all[tt], the angles themselves as a function of time for",
        "each residue are printed to separate files [TT](dihedral)(RESIDUE)(nresnr).xvg[tt].",
        "These can be in radians or degrees.[PAR]",
        "A log file (argument [TT]-g[tt]) is also written. This contains",
        "",
        " * information about the number of residues of each type.",
        " * The NMR ^3J coupling constants from the Karplus equation.",
        " * a table for each residue of the number of transitions between ",
        "   rotamers per nanosecond,  and the order parameter S^2 of each dihedral.",
        " * a table for each residue of the rotamer occupancy.",
        "",
        "All rotamers are taken as 3-fold, except for [GRK]omega[grk] and [GRK]chi[grk] dihedrals",
        "to planar groups (i.e. [GRK]chi[grk][SUB]2[sub] of aromatics, Asp and Asn; ",
        "[GRK]chi[grk][SUB]3[sub] of Glu",
        "and Gln; and [GRK]chi[grk][SUB]4[sub] of Arg), which are 2-fold. \"rotamer 0\" means ",
        "that the dihedral was not in the core region of each rotamer. ",
        "The width of the core region can be set with [TT]-core_rotamer[tt][PAR]",

        "The S^2 order parameters are also output to an [REF].xvg[ref] file",
        "(argument [TT]-o[tt] ) and optionally as a [REF].pdb[ref] file with",
        "the S^2 values as B-factor (argument [TT]-p[tt]). ",
        "The total number of rotamer transitions per timestep",
        "(argument [TT]-ot[tt]), the number of transitions per rotamer",
        "(argument [TT]-rt[tt]), and the ^3J couplings (argument [TT]-jc[tt]), ",
        "can also be written to [REF].xvg[ref] files. Note that the analysis",
        "of rotamer transitions assumes that the supplied trajectory frames",
        "are equally spaced in time.[PAR]",

        "If [TT]-chi_prod[tt] is set (and [TT]-maxchi[tt] > 0), cumulative rotamers, e.g.",
        "1+9([GRK]chi[grk][SUB]1[sub]-1)+3([GRK]chi[grk][SUB]2[sub]-1)+",
        "([GRK]chi[grk][SUB]3[sub]-1) (if the residue has three 3-fold ",
        "dihedrals and [TT]-maxchi[tt] >= 3)",
        "are calculated. As before, if any dihedral is not in the core region,",
        "the rotamer is taken to be 0. The occupancies of these cumulative ",
        "rotamers (starting with rotamer 0) are written to the file",
        "that is the argument of [TT]-cp[tt], and if the [TT]-all[tt] flag",
        "is given, the rotamers as functions of time",
        "are written to [TT]chiproduct(RESIDUE)(nresnr).xvg[tt] ",
        "and their occupancies to [TT]histo-chiproduct(RESIDUE)(nresnr).xvg[tt].[PAR]",

        "The option [TT]-r[tt] generates a contour plot of the average [GRK]omega[grk] angle",
        "as a function of the [GRK]phi[grk] and [GRK]psi[grk] angles, that is, in a Ramachandran ",
        "plot the average [GRK]omega[grk] angle is plotted using color coding.",

    };

    const char* bugs[] = {
        "N-terminal [GRK]phi[grk] and C-terminal [GRK]psi[grk] dihedrals are calculated in a "
        "non-standard way, using H-N-CA-C for [GRK]phi[grk] instead of "
        "C(-)-N-CA-C, and N-CA-C-O for [GRK]psi[grk] instead of N-CA-C-N(+). "
        "This causes (usually small) discrepancies with the output of other "
        "tools like [gmx-rama].",
        "Rotamers with multiplicity 2 are printed in [TT]chi.log[tt] as if they had ",
        "multiplicity 3, with the 3rd (g(+)) always having probability 0"
    };

    /* defaults */
    static int         r0 = 1, rN = -1, ndeg = 1, maxchi = 2;
    static gmx_bool    bAll = FALSE;
    static gmx_bool    bPhi = FALSE, bPsi = FALSE, bOmega = FALSE;
    static real        bfac_init = -1.0, bfac_max = 0;
    static const char* maxchistr[] = { nullptr, "0", "1", "2", "3", "4", "5", "6", nullptr };
    static gmx_bool    bRama = FALSE, bShift = FALSE, bViol = FALSE, bRamOmega = FALSE;
    static gmx_bool bNormHisto = TRUE, bChiProduct = FALSE, bHChi = FALSE, bRAD = FALSE, bPBC = TRUE;
    static real     core_frac = 0.5;
    t_pargs         pa[]      = {
        { "-r0", FALSE, etINT, { &r0 }, "starting residue" },
        { "-rN", FALSE, etINT, { &rN }, "last residue" },
        { "-phi", FALSE, etBOOL, { &bPhi }, "Output for [GRK]phi[grk] dihedral angles" },
        { "-psi", FALSE, etBOOL, { &bPsi }, "Output for [GRK]psi[grk] dihedral angles" },
        { "-omega",
          FALSE,
          etBOOL,
          { &bOmega },
          "Output for [GRK]omega[grk] dihedrals (peptide bonds)" },
        { "-rama",
          FALSE,
          etBOOL,
          { &bRama },
          "Generate [GRK]phi[grk]/[GRK]psi[grk] and "
          "[GRK]chi[grk][SUB]1[sub]/[GRK]chi[grk][SUB]2[sub] Ramachandran plots" },
        { "-viol",
          FALSE,
          etBOOL,
          { &bViol },
          "Write a file that gives 0 or 1 for violated Ramachandran angles" },
        { "-periodic", FALSE, etBOOL, { &bPBC }, "Print dihedral angles modulo 360 degrees" },
        { "-all", FALSE, etBOOL, { &bAll }, "Output separate files for every dihedral." },
        { "-rad",
          FALSE,
          etBOOL,
          { &bRAD },
          "in angle vs time files, use radians rather than degrees." },
        { "-shift",
          FALSE,
          etBOOL,
          { &bShift },
          "Compute chemical shifts from [GRK]phi[grk]/[GRK]psi[grk] angles" },
        { "-binwidth", FALSE, etINT, { &ndeg }, "bin width for histograms (degrees)" },
        { "-core_rotamer",
          FALSE,
          etREAL,
          { &core_frac },
          "only the central [TT]-core_rotamer[tt]\\*(360/multiplicity) belongs to each rotamer "
          "(the rest is assigned to rotamer 0)" },
        { "-maxchi", FALSE, etENUM, { maxchistr }, "calculate first ndih [GRK]chi[grk] dihedrals" },
        { "-normhisto", FALSE, etBOOL, { &bNormHisto }, "Normalize histograms" },
        { "-ramomega",
          FALSE,
          etBOOL,
          { &bRamOmega },
          "compute average omega as a function of [GRK]phi[grk]/[GRK]psi[grk] and plot it in an "
          "[REF].xpm[ref] plot" },
        { "-bfact",
          FALSE,
          etREAL,
          { &bfac_init },
          "B-factor value for [REF].pdb[ref] file for atoms with no calculated dihedral order "
          "parameter" },
        { "-chi_prod",
          FALSE,
          etBOOL,
          { &bChiProduct },
          "compute a single cumulative rotamer for each residue" },
        { "-HChi", FALSE, etBOOL, { &bHChi }, "Include dihedrals to sidechain hydrogens" },
        { "-bmax",
          FALSE,
          etREAL,
          { &bfac_max },
          "Maximum B-factor on any of the atoms that make up a dihedral, for the dihedral angle to "
          "be considered in the statistics. Applies to database work where a number of X-Ray "
          "structures is analyzed. [TT]-bmax[tt] <= 0 means no limit." }
    };

    FILE*             log;
    int               idum, nbin;
    rvec*             x;
    PbcType           pbcType;
    matrix            box;
    char              grpname[256];
    gmx_bool          bChi, bCorr, bSSHisto;
    gmx_bool          bDo_rt, bDo_oh, bDo_ot, bDo_jc;
    real              dt = 0, traj_t_ns;
    gmx_output_env_t* oenv;

    int    nactdih, nf;
    real **dih, *trans_frac, *aver_angle, *time;
    int ** chi_lookup, *multiplicity;

    t_filenm fnm[] = { { efSTX, "-s", nullptr, ffREAD },
                       { efTRX, "-f", nullptr, ffREAD },
                       { efXVG, "-o", "order", ffWRITE },
                       { efPDB, "-p", "order", ffOPTWR },
                       { efDAT, "-ss", "ssdump", ffOPTRD },
                       { efXVG, "-jc", "Jcoupling", ffWRITE },
                       { efXVG, "-corr", "dihcorr", ffOPTWR },
                       { efLOG, "-g", "chi", ffWRITE },
                       /* add two more arguments copying from gmx angle */
                       { efXVG, "-ot", "dihtrans", ffOPTWR },
                       { efXVG, "-oh", "trhisto", ffOPTWR },
                       { efXVG, "-rt", "restrans", ffOPTWR },
                       { efXVG, "-cp", "chiprodhisto", ffOPTWR } };
#define NFILE asize(fnm)
    int      npargs;
    t_pargs* ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);
    gmx::sfree_guard ppaGuard(ppa);
    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE, fnm, npargs, ppa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        return 0;
    }
    gmx::unique_cptr<gmx_output_env_t, output_env_done> oenvGuard(oenv);

    /* Handle result from enumerated type */
    sscanf(maxchistr[0], "%d", &maxchi);
    bChi = (maxchi > 0);

    log = gmx_ffopen(ftp2fn(efLOG, NFILE, fnm), "w");

    if (bRamOmega)
    {
        bOmega = TRUE;
        bPhi   = TRUE;
        bPsi   = TRUE;
    }

    /* set some options */
    bDo_rt = (opt2bSet("-rt", NFILE, fnm));
    bDo_oh = (opt2bSet("-oh", NFILE, fnm));
    bDo_ot = (opt2bSet("-ot", NFILE, fnm));
    bDo_jc = (opt2bSet("-jc", NFILE, fnm));
    bCorr  = (opt2bSet("-corr", NFILE, fnm));
    if (bCorr)
    {
        fprintf(stderr, "Will calculate autocorrelation\n");
    }

    if (core_frac > 1.0)
    {
        fprintf(stderr, "core_rotamer fraction > 1.0 ; will use 1.0\n");
        core_frac = 1.0;
    }
    if (core_frac < 0.0)
    {
        fprintf(stderr, "core_rotamer fraction < 0.0 ; will use 0.0\n");
        core_frac = 0.0;
    }

    if (maxchi > MAXCHI)
    {
        fprintf(stderr, "Will only calculate first %d Chi dihedrals instead of %d.\n", MAXCHI, maxchi);
        maxchi = MAXCHI;
    }
    bSSHisto = ftp2bSet(efDAT, NFILE, fnm);
    nbin     = 360 / ndeg;

    /* Find the chi angles using atoms struct and a list of amino acids */
    t_symtab symtab;
    char*    name;
    t_atoms  atoms;
    open_symtab(&symtab);
    gmx::unique_cptr<t_symtab, done_symtab> symtabGuard(&symtab);
    readConfAndAtoms(ftp2fn(efSTX, NFILE, fnm), &symtab, &name, &atoms, &pbcType, &x, nullptr, box);
    gmx::sfree_guard                     nameGuard(name);
    gmx::sfree_guard                     xGuard(x);
    gmx::unique_cptr<t_atoms, done_atom> atomsGuard(&atoms);
    if (atoms.pdbinfo == nullptr)
    {
        snew(atoms.pdbinfo, atoms.nr);
    }
    fprintf(log, "Title: %s\n", name);

    std::vector<t_dlist> dlist = mk_dlist(log, &atoms, bPhi, bPsi, bChi, bHChi, maxchi, r0, rN);
    fprintf(stderr, "%zu residues with dihedrals found\n", dlist.size());

    if (dlist.empty())
    {
        gmx_fatal(FARGS, "No dihedrals in your structure!\n");
    }

    /* Make a linear index for reading all dihedral atoms (4 per dihedral). */
    std::vector<int> index = make_chi_ind(dlist);
    int              ndih  = index.size() / 4; // 4 atoms per dihedral
    fprintf(stderr, "%d dihedrals found\n", ndih);

    snew(dih, ndih);

    /* COMPUTE ALL DIHEDRALS! */
    read_ang_dih(ftp2fn(efTRX, NFILE, fnm),
                 FALSE,
                 TRUE,
                 FALSE,
                 bPBC,
                 1,
                 &idum,
                 &nf,
                 &time,
                 index.size(),
                 index.data(),
                 &trans_frac,
                 &aver_angle,
                 dih,
                 oenv);
    gmx::sfree_guard timeGuard(time);
    gmx::sfree_guard transFracGuard(trans_frac);
    gmx::sfree_guard averAngleGuard(aver_angle);

    dt = (time[nf - 1] - time[0]) / (nf - 1); /* might want this for corr or n. transit*/
    if (bCorr)
    {
        if (nf < 2)
        {
            gmx_fatal(FARGS, "Need at least 2 frames for correlation");
        }
    }

    /* put angles in -M_PI to M_PI ! and correct phase factor for phi and psi
     * pass nactdih instead of ndih to low_ana_dih_trans
     * to prevent accessing off end of arrays when maxchi < 5 or 6. */
    nactdih = reset_em_all(dlist, nf, dih, maxchi);

    if (bAll)
    {
        dump_em_all(dlist, nf, time, dih, maxchi, bPhi, bPsi, bChi, bOmega, bRAD, oenv);
    }

    /* Histogramming & J coupling constants & calc of S2 order params */
    histogramming(log,
                  nbin,
                  nf,
                  maxchi,
                  dih,
                  dlist,
                  index,
                  bPhi,
                  bPsi,
                  bOmega,
                  bChi,
                  bNormHisto,
                  bSSHisto,
                  ftp2fn(efDAT, NFILE, fnm),
                  bfac_max,
                  &atoms,
                  bDo_jc,
                  opt2fn("-jc", NFILE, fnm),
                  oenv);

    /* transitions
     *
     * added multiplicity */

    snew(multiplicity, ndih);
    gmx::sfree_guard multiplicityGuard(multiplicity);
    mk_multiplicity_lookup(multiplicity, maxchi, dlist, ndih);

    std::strcpy(grpname, "All residues, ");
    if (bPhi)
    {
        std::strcat(grpname, "Phi ");
    }
    if (bPsi)
    {
        std::strcat(grpname, "Psi ");
    }
    if (bOmega)
    {
        std::strcat(grpname, "Omega ");
    }
    if (bChi)
    {
        std::strcat(grpname, "Chi 1-");
        sprintf(grpname + std::strlen(grpname), "%i", maxchi);
    }


    low_ana_dih_trans(bDo_ot,
                      opt2fn("-ot", NFILE, fnm),
                      bDo_oh,
                      opt2fn("-oh", NFILE, fnm),
                      maxchi,
                      dih,
                      dlist,
                      nf,
                      nactdih,
                      grpname,
                      multiplicity,
                      time,
                      FALSE,
                      core_frac,
                      oenv);

    /* Order parameters */
    order_params(log,
                 opt2fn("-o", NFILE, fnm),
                 maxchi,
                 dlist,
                 ftp2fn_null(efPDB, NFILE, fnm),
                 bfac_init,
                 &atoms,
                 x,
                 pbcType,
                 box,
                 bPhi,
                 bPsi,
                 bChi,
                 oenv);

    /* Print ramachandran maps! */
    if (bRama)
    {
        do_rama(nf, dlist, dih, bViol, bRamOmega, oenv);
    }

    if (bShift)
    {
        do_pp2shifts(log, nf, dlist, dih);
    }

    /* rprint S^2, transitions, and rotamer occupancies to log */
    traj_t_ns = 0.001 * (time[nf - 1] - time[0]);
    pr_dlist(log, dlist, traj_t_ns, edPrintST, bPhi, bPsi, bChi, bOmega, maxchi);
    pr_dlist(log, dlist, traj_t_ns, edPrintRO, bPhi, bPsi, bChi, bOmega, maxchi);
    gmx_ffclose(log);
    /* transitions to xvg */
    if (bDo_rt)
    {
        print_transitions(opt2fn("-rt", NFILE, fnm), maxchi, dlist, traj_t_ns, oenv);
    }

    /* chi_product trajectories (ie one "rotamer number" for each residue) */
    if (bChiProduct && bChi)
    {
        snew(chi_lookup, dlist.size());
        for (size_t i = 0; i < dlist.size(); i++)
        {
            snew(chi_lookup[i], maxchi);
        }
        mk_chi_lookup(chi_lookup, maxchi, dlist);

        get_chi_product_traj(dih,
                             nf,
                             maxchi,
                             dlist,
                             time,
                             chi_lookup,
                             multiplicity,
                             FALSE,
                             bNormHisto,
                             core_frac,
                             bAll,
                             opt2fn("-cp", NFILE, fnm),
                             oenv);

        for (size_t i = 0; i < dlist.size(); i++)
        {
            sfree(chi_lookup[i]);
        }
        sfree(chi_lookup);
    }

    /* Correlation comes last because it messes up the angles */
    if (bCorr)
    {
        do_dihcorr(opt2fn("-corr", NFILE, fnm), nf, ndih, dih, dt, dlist, time, maxchi, bPhi, bPsi, bChi, bOmega, oenv);
    }


    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");
    do_view(oenv, opt2fn("-jc", NFILE, fnm), "-nxy");
    if (bCorr)
    {
        do_view(oenv, opt2fn("-corr", NFILE, fnm), "-nxy");
    }

    for (int i = 0; (i < ndih); i++)
    {
        sfree(dih[i]);
    }
    sfree(dih);

    return 0;
}
