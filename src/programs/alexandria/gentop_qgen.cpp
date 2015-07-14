/*
 * This source file is part of the Aleandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"
#include <ctype.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/units.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/linearalgebra/matrix.h"
#include "coulombintegrals/coulombintegrals.h"
#include "molprop.h"
#include "poldata.h"
#include "gentop_qgen.h"
#include "gmx_resp.h"

namespace alexandria
{


  GentopQgen::GentopQgen(Poldata * pd, t_atoms *atoms, gmx_atomprop_t aps,
			 rvec *x,
			 ChargeDistributionModel   iChargeDistributionModel,
			 ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
			 real hfac, int qtotal, real epsr)
  {
    bWarned = false;  
    bAllocSave = false;  
    natom = 0;
    eQGEN = 0;
    qtotal = 0; 
    chieq = 0;
    hfac = 0;
    epsr = 0;
    char        *atp;
    gmx_bool     bSup = TRUE;
    int          i, j, k, atm, nz;

   
    iChargeDistributionModel   = iChargeDistributionModel;
    iChargeGenerationAlgorithm = iChargeGenerationAlgorithm;
    hfac                       = hfac;
    qtotal                     = qtotal;
    if (epsr <= 1)
      {
        epsr = 1;
      }
    epsr   = epsr;
    for (i = j = 0; (i < atoms->nr); i++)
      {
        if (atoms->atom[i].ptype == eptAtom)
	  {
            natom++;
	  }
      }
    
    chi0.resize(natom);
    
    rhs.resize(natom+1);
    
    elem.resize(natom);
   
    atomnr.resize(natom);
    row.resize(natom);
       
    Jab.resize(natom);
    zeta.resize(natom);
    j00.resize(natom);
    q.resize(natom);
    //this->x.resize(natom);

    bAllocSave = FALSE;
   
    nZeta.resize(natom);

    /* Special case for chi_eq */
    nZeta[natom] = 1;
   
    q[natom].resize(natom);

    for (i = j = 0; (i < atoms->nr) && bSup; i++)
      {
        if (atoms->atom[i].ptype == eptAtom)
	  {
            
	    Jab[natom].resize(natom+1);
            atm = atoms->atom[i].atomnumber;
            if (atm == NOTSET)
	      {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
	      }
            atp = *atoms->atomtype[j];
            if (pd->haveEemSupport(iChargeDistributionModel, atp, TRUE) == 0)
	      {
                atp = gmx_atomprop_element(aps, atm);
                if (pd->haveEemSupport(iChargeDistributionModel, atp, TRUE) == 0)
		  {
                    fprintf(stderr, "No charge distribution support for atom %s (element %s), model %s\n",
                            *atoms->atomtype[j], atp, Poldata::getEemtypeName(iChargeDistributionModel));
                    bSup = FALSE;
		  }
	      }
            if (bSup)
	      {
                elem[j]   = strdup(atp);
                atomnr[j] = atm;
                nz              = pd->getNzeta(iChargeDistributionModel, atp);
                nZeta[j]  = nz;
                
		q[j].resize(nz);
                
		zeta[j].resize(nz);
                
		row[j].resize(nz);
                for (k = 0; (k < nz); k++)
		  {
                    q[j][k]    = pd->getQ(iChargeDistributionModel, *atoms->atomtype[j], k);
                    zeta[j][k] = pd->getZeta(iChargeDistributionModel, *atoms->atomtype[j], k);
                    row[j][k]  = pd->getRow(iChargeDistributionModel, *atoms->atomtype[j], k);
                    if (row[j][k] > SLATER_MAX)
		      {
                        if (debug)
			  {
                            fprintf(debug, "Can not handle higher slaters than %d for atom %s %s\n",
                                    SLATER_MAX,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
			  }
                        row[j][k] = SLATER_MAX;
		      }
		  }
                chi0[j]  = 0;
                j00[j]   = 0;
                copy_rvec(x[i], x[j]);
                j++;
	      }
	  }
      }
    if (bSup)
      {
      }
    else
      {
	//    done();
      }
  }




  GentopQgen::~GentopQgen()
  {
    //    int  i;

    /*   sfree(chi0);
	 sfree(rhs);
	 sfree(atomnr);
	 sfree(j00);
	 sfree(x);
	 for (i = 0; (i < natom); i++)
	 {
         sfree(row[i]);
	 sfree(q[i]);
	 sfree(zeta[i]);
	 sfree(Jab[i]);
	 sfree(elem[i]);
	 if (bAllocSave)
	 {
	 sfree(qsave[i]);
	 sfree(zetasave[i]);
	 }
	 }
	 sfree(row);
	 sfree(zeta);
	 sfree(elem);
	 sfree(q);
	 sfree(Jab);
	 sfree(nZeta);
	 if (bAllocSave)
	 {
	 sfree(qsave);
	 sfree(zetasave);
	 bAllocSave = FALSE;
	 }*/
  }


  void GentopQgen::saveParams( Resp * gr)
  {
    int i, j;

    if (!bAllocSave)
      {
	qsave.resize(natom);
	zetasave.resize(natom);
      }
    for (i = 0; (i < natom); i++)
      {
        if (!bAllocSave)
	  {
	    qsave[i].resize(nZeta[i]);
	    zetasave[i].resize(nZeta[i]);
	  }
        for (j = 0; (j < nZeta[i]); j++)
	  {
            if (NULL != gr)
	      {
		q[i][j]    = (real)gr->get_q( i, j);
                zeta[i][j] = gr->get_zeta( i, j);
	      }
            qsave[i][j]    = q[i][j];
            zetasave[i][j] = zeta[i][j];
	  }
      }
    bAllocSave = TRUE;
  }

  void GentopQgen::getParams( Resp * gr)
  {
    int i, j;

    if (bAllocSave)
      {
        for (i = 0; (i < natom); i++)
	  {
            for (j = 0; (j < nZeta[i]); j++)
	      {
                q[i][j]    = qsave[i][j];
                zeta[i][j] = zetasave[i][j];
                if (NULL != gr)
		  {
                    gr->set_q( i, j, q[i][j]);
                    gr->set_zeta( i, j, zeta[i][j]);
		  }
	      }
	  }
      }
    else
      {
        fprintf(stderr, "WARNING: no ESP charges generated.\n");
      }
  }

  int GentopQgen::getNzeta( int atom)
  {
    if ((0 <= atom) && (atom < natom))
      {
        return nZeta[atom];
      }
    return NOTSET;
  
  }

  int GentopQgen::getRow( int atom, int z)
  {
    if ((0 <= atom) && (atom < natom) &&
        (0 <= z) && (z <= nZeta[atom]))
      {
        return row[atom][z];
      }
    return NOTSET;
  
  }
  double GentopQgen::getQ(int atom, int z)
  {
    if ((0 <= atom) && (atom < natom) &&
        (0 <= z) && (z <= nZeta[atom]))
      {
        return q[atom][z];
      }
    return NOTSET;

  }


  double GentopQgen::getZeta(int atom, int z)
  {
    if ((0 <= atom) && (atom < natom) &&
        (0 <= z) && (z <= nZeta[atom]))
      {
        return zeta[atom][z];
      }
    return NOTSET;
  }

  real CoulombNN(real r)
  {
    return 1/r;
  }

  real GentopQgen::calcJab(ChargeDistributionModel iChargeDistributionModel,
			   rvec xi, rvec xj,
			   int nZi, int nZj,
			   std::vector<real> zeta_i, std::vector<real> zeta_j,
			   std::vector<int> rowi, std::vector<int> rowj)
  {
    int  i, j;
    rvec dx;
    real r;
    real eTot = 0;

    rvec_sub(xi, xj, dx);
    r = norm(dx);
    if (r == 0)
      {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
      }
    //    if ((*zeta_i <= 0) || (*zeta_j <= 0)) correct TODO
    if ((zeta_i[0] <= 0) || (zeta_j[0] <= 0))
      {
        iChargeDistributionModel = eqdAXp;
      }
    switch (iChargeDistributionModel)
      {
      case eqdAXp:
	eTot = CoulombNN(r);
	break;
      case eqdAXs:
      case eqdRappe:
      case eqdYang:
	eTot = 0;
	for (i = nZi-1; (i < nZi); i++)
	  {
	    for (j = nZj-1; (j < nZj); j++)
	      {
		eTot += Coulomb_SS(r, rowi[i], rowj[j], zeta_i[i], zeta_j[j]);
	      }
	  }
	break;
      case eqdAXg:
	eTot = 0;
	for (i = nZi-1; (i < nZi); i++)
	  {
	    for (j = nZj-1; (j < nZj); j++)
	      {
		eTot += Coulomb_GG(r, zeta_i[i], zeta_j[j]);
	      }
	  }
	break;
      default:
	gmx_fatal(FARGS, "Unsupported model %d in calc_jab", iChargeDistributionModel);
      }

    return ONE_4PI_EPS0*(eTot)/ELECTRONVOLT;
  }

  void GentopQgen::solveQEem(FILE *fp,  real hardness_factor)
  {
    double **a, qtot, q;
    int      i, j, n;

    n = natom+1;
    a = alloc_matrix(n, n);
    for (i = 0; (i < n-1); i++)
      {
        for (j = 0; (j < n-1); j++)
	  {
            a[i][j] = Jab[i][j];
	  }
        a[i][i] = hardness_factor*Jab[i][i];
      }
    for (j = 0; (j < n-1); j++)
      {
        a[n-1][j] = 1;
      }
    for (i = 0; (i < n-1); i++)
      {
        a[i][n-1] = -1;
      }
    a[n-1][n-1] = 0;

    if (matrix_invert(fp, n, a) == 0)
      {
        for (i = 0; (i < n); i++)
	  {
            q = 0;
            for (j = 0; (j < n); j++)
	      {
                q += a[i][j]*rhs[j];
	      }
            this->q[i][nZeta[i]-1] = q;
            if (fp)
	      {
                fprintf(fp, "%2d RHS = %10g Charge= %10g\n", i, rhs[i], q);
	      }
	  }
      }
    else
      {
        for (i = 0; (i < n); i++)
	  {
            this->q[i][nZeta[i]] = 1;
	  }
      }
    chieq = this->q[n-1][nZeta[n-1]-1];
    qtot        = 0;
    for (i = 0; (i < n-1); i++)
      {
        for (j = 0; (j < nZeta[i]); j++)
	  {
            qtot += this->q[i][j];
	  }
      }

    if (fp && (fabs(qtot - qtotal) > 1e-2))
      {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, qtotal);
      }
    free_matrix(a);
  }

  void GentopQgen::updateJ00()
  {
    int    i;
    double j0, qq;
    double zetaH = 1.0698;

    for (i = 0; (i < natom); i++)
      {
        j0 = j00[i]/epsr;
        if (((iChargeDistributionModel == eqdYang) ||
             (iChargeDistributionModel == eqdRappe)) &&
            (atomnr[i] == 1))
	  {
            qq = q[i][nZeta[i]-1];
            j0 = (1+qq/zetaH)*j0;

            if (debug && (j0 < 0) && !bWarned)
	      {
                fprintf(debug, "WARNING: J00 = %g for atom %d. The equations will be instable.\n", j0, i+1);
                bWarned = TRUE;
	      }
	  }
        Jab[i][i] = (j0 > 0) ? j0 : 0;
      }
  }

  void GentopQgen::debugFun(FILE *fp)
  {
    int i, j;

    for (i = 0; (i < natom); i++)
      {
        fprintf(fp, "THIS: i: %2d chi0: %8g J0: %8g q:",
                i+1, chi0[i], Jab[i][i]);
        for (j = 0; (j < nZeta[i]); j++)
	  {
            fprintf(fp, " %8g", q[i][j]);
	  }
        fprintf(fp, "\n");
      }
    fprintf(fp, "qgen Jab matrix:\n");
    for (i = 0; (i < natom); i++)
      {
        for (j = 0; (j <= i); j++)
	  {
            fprintf(fp, "  %6.2f", Jab[i][j]);
	  }
        fprintf(fp, "\n");
      }
    fprintf(fp, "\n");
  }

  real GentopQgen::calcSij(int i, int j)
  {
    real dist, dism, Sij = 1.0;
    rvec dx;
    int  l, m, tag;

    rvec_sub(x[i], x[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (atomnr[i] != 1) && (atomnr[j] != 1))
      {
        Sij = Sij*1.64;
      }
    else if ((dist < 0.122) && (atomnr[i] != 1) && (atomnr[j] != 1))
      {
        if ((atomnr[i] != 8) && (atomnr[j] != 8))
	  {
            Sij = Sij*2.23;
	  }
        else
	  {
            Sij = Sij*2.23;
	  }
      }
    else if (dist < 0.125)
      {
        tag = 0;
        if ((atomnr[i] == 6) && (atomnr[j] == 8))
	  {
            tag = i;
	  }
        else if ((atomnr[i] == 8) && (atomnr[j] == 6))
	  {
            tag = j;
	  }
        if (tag != 0)
	  {
            printf("found CO\n");
            for (l = 0; (l < natom); l++)
	      {
                if (atomnr[l] == 1)
		  {
                    printf("found H\n");
                    dism = 0.0;
                    for (m = 0; (m < DIM); m++)
		      {
                        dism = dism+sqr(x[tag][m]-x[l][m]);
		      }

                    printf("dist: %8.3f\n", sqrt(dism));
                    if (sqrt(dism) < 0.105)
		      {
                        printf("dist %5d %5d %5s  %5s %8.3f\n",
                               i, l, elem[tag], elem[l], sqrt(dism));
                        Sij = Sij*1.605;
		      }
		  }
	      }
	  }
      }
    else if ((atomnr[i] == 6) && (atomnr[j] == 8))
      {
        Sij = Sij*1.03;
      }
    else if (((atomnr[j] == 6) && (atomnr[i] == 7) && (dist < 0.15)) ||
             ((atomnr[i] == 6) && (atomnr[j] == 7) && (dist < 0.15)))
      {
        if (atomnr[i] == 6)
	  {
            tag = i;
	  }
        else
	  {
            tag = j;
	  }
        for (l = 0; (l < natom); l++)
	  {
            if (atomnr[l] == 8)
	      {
                printf("found Oxy\n");
                dism = 0.0;
                for (m = 0; (m < DIM); m++)
		  {
                    dism = dism+sqr(x[tag][m]-x[l][m]);
		  }
                if (sqrt(dism) < 0.130)
		  {
                    printf("found peptide bond\n");
                    Sij = Sij*0.66;
		  }
                else
		  {
                    Sij = Sij*1.1;
		  }
	      }
	  }
      }
    return Sij;
  }

  void GentopQgen::calcJab()
  {
    int    i, j;
    double Jab;

    for (i = 0; (i < natom); i++)
      {
        for (j = i+1; (j < natom); j++)
	  {
            Jab = calcJab(iChargeDistributionModel,
			  x[i], x[j],
			  nZeta[i], nZeta[j],
			  zeta[i], zeta[j],
			  row[i], row[j]);
            if (iChargeDistributionModel == eqdYang)
	      {
                Jab = Jab*calcSij(i, j);
	      }
            this->Jab[j][i] = this->Jab[i][j] = Jab/epsr;
	  }
      }
  }

  void GentopQgen::calcRhs()
  {
    int    i, j, k, l;
    rvec   dx;
    real   r, j1, j1q, qcore;

    /* This right hand side is for all models */
    for (i = 0; (i < natom); i++)
      {
        rhs[i] = -chi0[i];
      }
    rhs[natom] = qtotal;

    /* In case the charge is split in nuclear charge and electronic charge
     * we need to add some more stuff. See paper for details.
     */
    for (i = 0; (i < natom); i++)
      {
        j1q   = 0;
        qcore = 0;
        for (k = 0; (k < nZeta[i]-1); k++)
	  {
            j1q   += j00[i]*q[i][k];
            qcore += q[i][k];
	  }
        j1 = 0;
        /* This assignment for k is superfluous because of the previous loop,
         * but if I take it out it will at some stage break the loop below where
         * exactly this value of k is needed.
         */
        k  = nZeta[i]-1;
        for (j = 0; (j < natom); j++)
	  {
            if (i != j)
	      {
                rvec_sub(x[i], x[j], dx);
                r = norm(dx);
                switch (iChargeDistributionModel)
		  {
		  case eqdAXs:
		  case eqdRappe:
		  case eqdYang:
		    for (l = 0; (l < nZeta[j]-1); l++)
		      {
			j1 += q[j][l]*Coulomb_SS(r, k, l, zeta[i][k], zeta[j][l]);
		      }
		    break;
		  case eqdAXg:
		    for (l = 0; (l < nZeta[j]-1); l++)
		      {
			j1 += q[j][l]*Coulomb_GG(r, zeta[i][k], zeta[j][l]);
		      }
		    break;
		  default:
		    break;
		  }
	      }
	  }
        rhs[i]           -= j1q + ONE_4PI_EPS0*j1/ELECTRONVOLT;
        rhs[natom] -= qcore;
      }
  }

  int atomicnumber2rowXX(int elem)
  {
    int row;

    /* Compute which row in the periodic table is this element */
    if (elem <= 2)
      {
        row = 1;
      }
    else if (elem <= 10)
      {
        row = 2;
      }
    else if (elem <= 18)
      {
        row = 3;
      }
    else if (elem <= 36)
      {
        row = 4;
      }
    else if (elem <= 54)
      {
        row = 5;
      }
    else if (elem <= 86)
      {
        row = 6;
      }
    else
      {
        row = 7;
      }

    return row;
  }


  void GentopQgen::print(FILE *fp, t_atoms *atoms)
  {
    int  i, j, k, m;
    rvec mu = { 0, 0, 0 };
    real qq;

    if (eQGEN == eQGEN_OK)
      {
        if (fp)
	  {
            fprintf(fp, "Res  Atom   Nr       J0     chi0 row        q zeta (1/nm)\n");
	  }
        for (i = j = 0; (i < atoms->nr); i++)
	  {
            if (atoms->atom[i].ptype == eptAtom)
	      {
                qq = 0;
                for (k = 0; (k < nZeta[j]); k++)
		  {
                    qq += q[j][k];
		  }

                atoms->atom[i].q = qq;
                for (m = 0; (m < DIM); m++)
		  {
                    mu[m] += qq* x[i][m] * ENM2DEBYE;
		  }
                if (fp)
		  {
                    fprintf(fp, "%4s %4s%5d %8g %8g",
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            *(atoms->atomname[i]), i+1, j00[j], chi0[j]);
                    for (k = 0; (k < nZeta[j]); k++)
		      {
                        fprintf(fp, " %3d %8.5f %8.4f", row[j][k], q[j][k],
                                zeta[j][k]);
		      }
                    fprintf(fp, "\n");
		  }
                j++;
	      }
	  }
        if (fp)
	  {
            fprintf(fp, "<chieq> = %10g\n|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n",
                    chieq, norm(mu), mu[XX], mu[YY], mu[ZZ]);
	  }
      }
  }

  void GentopQgen::message( int len, char buf[], Resp * gr)
  {
    switch (eQGEN)
      {
      case eQGEN_OK:
	if (NULL != gr)
	  {
	    gr->calc_pot();
	    gr->calc_rms();
	    gr->statistics( len, buf);
	  }
	else
	  {
	    sprintf(buf, "Charge generation finished correctly.\n");
	  }
	break;
      case eQGEN_NOTCONVERGED:
	sprintf(buf, "Charge generation did not converge.\n");
	break;
      case eQGEN_NOSUPPORT:
	sprintf(buf, "No charge generation support for (some of) the atomtypes.\n");
	break;
      case eQGEN_ERROR:
      default:
	sprintf(buf, "Unknown status %d in charge generation.\n", eQGEN);
      }
  }

  void GentopQgen::checkSupport(Poldata * pd, gmx_atomprop_t aps)
  {
    int      i;
    gmx_bool bSup = TRUE;

    for (i = 0; (i < natom); i++)
      {
        if (pd->haveEemSupport(iChargeDistributionModel, elem[i], TRUE) == 0)
	  {
            /*sfree(elem[i]);*/
            elem[i] = strdup(gmx_atomprop_element(aps, atomnr[i]));
            if (pd->haveEemSupport(iChargeDistributionModel, elem[i], TRUE) == 0)
	      {
                fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                        elem[i], Poldata::getEemtypeName(iChargeDistributionModel));
                bSup = FALSE;
	      }
	  }
      }
    if (bSup)
      {
        eQGEN = eQGEN_OK;
      }
    else
      {
        eQGEN = eQGEN_NOSUPPORT;
      }
  }

  void GentopQgen::updatePd(t_atoms *atoms, Poldata * pd)
  {
    int i, j, n, nz;

    for (i = j = 0; (i < atoms->nr); i++)
      {
        if (atoms->atom[i].ptype == eptAtom)
	  {
            chi0[j]  = pd->getChi0(iChargeDistributionModel, elem[j]);
            j00[j]   = pd->getJ00(iChargeDistributionModel, elem[j]);
            nz             = pd->getNzeta(iChargeDistributionModel, elem[j]);
            for (n = 0; (n < nz); n++)
	      {
                zeta[j][n] = pd->getZeta(iChargeDistributionModel,elem[j], n);
                q[j][n]    = pd->getQ(iChargeDistributionModel,elem[j], n);
                row[j][n]  = pd->getRow(iChargeDistributionModel,elem[j], n);
	      }
            j++;
	  }
      }
  }

  int GentopQgen::generateChargesSm(FILE *fp,
				    Poldata * pd,
				    t_atoms *atoms,
				    real tol, int maxiter, gmx_atomprop_t aps,
				    real *chieq)
  {
    real       *qq = NULL;
    int         i, j, iter;
    real        rms;

    checkSupport(pd, aps);
    if (eQGEN_OK == eQGEN)
      {

        updatePd(atoms, pd);

        snew(qq, atoms->nr+1);
        for (i = j = 0; (i < atoms->nr); i++)
	  {
            if (atoms->atom[i].ptype != eptShell)
	      {
                qq[j] = q[j][nZeta[j]-1];
                j++;
	      }
	  }
        iter = 0;
        calcJab();
        calcRhs();
        do
	  {
            updateJ00();
            if (debug)
	      {
                debugFun(debug);
	      }
            solveQEem(debug, 1.0);
            rms = 0;
            for (i = j = 0; (i < atoms->nr); i++)
	      {
                if (atoms->atom[i].ptype != eptShell)
		  {
                    rms  += sqr(qq[j] - q[j][nZeta[j]-1]);
                    qq[j] = q[j][nZeta[j]-1];
                    j++;
		  }
	      }
            rms = sqrt(rms/atoms->nr);
            iter++;
	  }
        while ((rms > tol) && (iter < maxiter));

        if (iter < maxiter)
	  {
            eQGEN = eQGEN_OK;
	  }
        else
	  {
            eQGEN = eQGEN_NOTCONVERGED;
	  }

        if (fp)
	  {
            if (eQGEN == eQGEN_OK)
	      {
                fprintf(fp, "Converged to tolerance %g after %d iterations\n",
                        tol, iter);
	      }
            else
	      {
                fprintf(fp, "Did not converge within %d iterations. RMS = %g\n",
                        maxiter, rms);
	      }
	  }
        *chieq = this->chieq;
        sfree(qq);
      }

    if (eQGEN_OK == eQGEN)
      {
        print(fp, atoms);
      }

    return eQGEN;
  }

  int GentopQgen::generateChargesBultinck(FILE *fp,
					  Poldata * pd, t_atoms *atoms,
					  gmx_atomprop_t aps)
  {
    checkSupport(pd, aps);
    if (eQGEN_OK == eQGEN)
      {
        updatePd(atoms, pd);

        calcJab();
        calcRhs();
        updateJ00();
        solveQEem(debug, 2.0);

        print(fp, atoms);
      }

    return eQGEN;
  }

  int GentopQgen::generateCharges(FILE *fp,
				  Resp * gr,
				  const char *molname, Poldata * pd,
				  t_atoms *atoms,
				  real tol, int maxiter, int maxcycle,
				  gmx_atomprop_t aps)
  {
    int  cc, eQGEN_min = eQGEN_NOTCONVERGED;
    real chieq, chi2, chi2min = GMX_REAL_MAX;

    /* Generate charges */
    switch (iChargeGenerationAlgorithm)
      {
      case eqgRESP:
	if (NULL == gr)
	  {
	    gmx_incons("No RESP data structure");
	  }
	if (fp)
	  {
	    fprintf(fp, "Generating %s charges for %s using RESP algorithm\n",
		    Poldata::getEemtypeName(iChargeDistributionModel), molname);
	  }
	for (cc = 0; (cc < maxcycle); cc++)
	  {
	    if (fp)
	      {
		fprintf(fp, "Cycle %d/%d\n", cc+1, maxcycle);
	      }
	    /* Fit charges to electrostatic potential */
	    eQGEN = gr->optimize_charges(fp, maxiter, tol, &chi2);
	    if (eQGEN == eQGEN_OK)
	      {
		eQGEN_min = eQGEN;
		if (chi2 <= chi2min)
		  {
		    saveParams(gr);
		    chi2min = chi2;
		  }

		if (NULL != fp)
		  {
		    fprintf(fp, "chi2 = %g kJ/mol e\n", chi2);
		  }
		print(fp, atoms);
	      }
	  }
	if (maxcycle > 1)
	  {
	    if (fp)
	      {
		fprintf(fp, "---------------------------------\nchi2 at minimum is %g\n", chi2min);
	      }
	    getParams(gr);
	    print(fp, atoms);
	  }
	eQGEN = eQGEN_min;
	break;
      default:
	/* Use empirical algorithms */
	if (fp)
	  {
	    fprintf(fp, "Generating charges for %s using %s algorithm\n",
		    molname, Poldata::getEemtypeName(iChargeDistributionModel));
	  }
	if (iChargeDistributionModel == eqdBultinck)
	  {
	    (void) generateChargesBultinck(fp, pd, atoms, aps);
	  }
	else
	  {
	    (void) generateChargesSm(fp, pd, atoms, tol, maxiter, aps, &chieq);
	  }
	saveParams(gr);
      }
    return eQGEN;
  }

}
