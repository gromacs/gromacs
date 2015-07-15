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
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/legacyheaders/nmsimplex.h"
#include "coulombintegrals/coulombintegrals.h"
#include "poldata.h"
#include "gmx_resp.h"
#include "gentop_qgen.h"
#include "stringutil.h"




namespace alexandria
{
  Resp::Resp(ChargeDistributionModel iDistributionModel,
	     bool bAXpRESP, real qfac, real b_hyper, real qtot,
	     real zmin, real zmax, real deltaZ, bool bZatype,
	     real watoms, real rDecrZeta, bool bRandZeta,
	     bool bRandQ, real penalty_fac, bool bFitZeta, bool bEntropy,
	     const char *dzatoms, unsigned int seed)
  {



    _qtot                  = qtot;
    _qsum                  = qtot;
    _bAXpRESP              = bAXpRESP;
    _qfac                  = qfac;
    _bHyper               = b_hyper;
    _wtot                  = 0;
    _iDistributionModel    = iDistributionModel;
    _zmin                  = zmin;
    _zmax                  = zmax;
    _deltaZ               = deltaZ;
    std::vector<std::string> ptr = split(dzatoms, ' ');
    for (unsigned int i = 0; (i < ptr.size()); ++i)
      {
        _dzatoms[i].assign(ptr[i].c_str());
      }
    _pfac      = penalty_fac;
    _qmin      = -2;
    _qmax      = 2; /* e */
    _nesp      = 0;
    _natom     = 0;
    _natype    = 0;
    _seed      = seed;
    _bEntropy  = bEntropy;
    _bZatype   = bZatype;
    _rDecrZeta = rDecrZeta;
    _bRandZeta = bRandZeta;
    _bRandQ    = bRandQ;
    _bFitZeta  = bFitZeta;
    _watoms    = watoms;
    _nparam    = 0;

  
  }

  Resp::~Resp()
  {
    //int i;

    sfree(_x);
    sfree(_esp);
    //  sfree(_pot);
    //    i = 0;
    /* for (i = 0; (i < _natom); i++)
      {
        delete _ra[i];
	}*/
    //    sfree(_ra);
  }

  void Resp::getAtomInfo( t_atoms *atoms,
                            t_symtab *symtab, rvec **x)
  {
    int          i;
    const char  *rnm;

    init_t_atoms(atoms, _natom, true);
    if (NULL == (*x))
      {
        snew((*x), atoms->nr);
      }
    if (NULL == atoms->atomtype)
      {
        snew(atoms->atomtype, atoms->nr);
      }
    for (i = 0; (i < _natom); i++)
      {
        atoms->atom[i].atomnumber = _ra[i]->getAtomnumber();
        atoms->atom[i].q          = _ra[i]->getQ();
        atoms->atomname[i]        = put_symtab(symtab, _ra[i]->getAtomtype().c_str());
        atoms->atomtype[i]        = put_symtab(symtab, _ra[i]->getAtomtype().c_str());
        atoms->atom[i].resind     = 0;

        strncpy(atoms->atom[i].elem, _ra[i]->getAtomtype().c_str(),
                sizeof(atoms->atom[i].elem)-1);
        copy_rvec(_x[i], (*x)[i]);
      }
    rnm = ("" != _stoichiometry) ? _stoichiometry.c_str() : (const char *)"BOE";
    t_atoms_set_resinfo(atoms, 0, symtab, rnm, 1, ' ', 1, ' ');

    atoms->nres = 1;
  }

  void Resp::updateAtomtypes( t_atoms *atoms)
  {
    int i, j;

    for (i = 0; (i < _natom); i++)
      {
        _ra[i]->getAtomtype() = strdup(*atoms->atomtype[i]);
        for (j = 0; (j < i); j++)
	  {
            if (0 == strcmp(*atoms->atomtype[i], *atoms->atomtype[j]))
	      {
                break;
	      }
	  }
        if (j == i)
	  {
            _natype++;
	  }
        _ra[i]->setAtype(j);
      }
  }

  void Resp::addAtomCoords( rvec *x)
  {
    int        i;

    srenew(_x, _natom);
    for (i = 0; (i < _natom); i++)
      {
        copy_rvec(x[i], _x[i]);
      }
  }

  void Resp::fillZeta( alexandria::Poldata * pd)
  {
    int i, zz;

    for (i = 0; (i < _natom); i++)
      {
        for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	  {
            setZeta( i, zz, pd->getZeta( _iDistributionModel,
					 _ra[i]->getAtomtype().c_str(), zz));
	  }
      }
  }

  void Resp::fillQ( t_atoms *atoms)
  {
    int    i, zz;
    double q;

    for (i = 0; (i < _natom); i++)
      {
        q = 0;
        for (zz = 0; (zz < _ra[i]->getNZeta()-1); zz++)
	  {
            q -= _ra[i]->getQ(zz);
	  }
        q += atoms->atom[i].q;
        setQ( i, _ra[i]->getNZeta()-1, q);
      }
  }

  bool Resp::addAtomInfo( t_atoms *atoms, alexandria::Poldata * pd)
  {
    int  i;

    _natom    = atoms->nr;
    _ra.resize(_natom);

    for (i = 0; (i < _natom); i++)
      {
	_ra[i] =  new Ra(atoms->atom[i].atomnumber, atoms->atom[i].type,
			 *(atoms->atomtype[i]), pd, _iDistributionModel, _dzatoms);
        if (_ra[i]->setUpcorrectly() == false)
	  {
            return false;
	  }
      }
    return true;
  }

  void Resp::summary(FILE *fp, 
		     std::vector<int> &symmetric_atoms)
  {
    int i;

    if (NULL != fp)
      {
        fprintf(fp, "There are %d atoms, %d atomtypes %d parameters for (R)ESP fitting.\n",
                _natom, _natype, _nparam);
        for (i = 0; (i < _natom); i++)
	  {
            fprintf(fp, " %d", symmetric_atoms[i]);
	  }
        fprintf(fp, "\n");
      }
  }



  void Resp::addParam( int atom, eParm eparm, int zz)
  {
    range_check(atom, 0, _natom);
    if ((zz >= 0) && (zz < _ra[atom]->getNZeta()))
      {
        if (eparm == eparmQ)
	  {
            _ra[atom]->setIq(zz, _nparam++);
            if (debug)
	      {
                fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
	      }
	  }
        else if (_bFitZeta)
	  {
            if (_ra[atom]->getZeta(zz) != 0)
	      {
                _ra[atom]->setIz(zz, _nparam++);
                if (debug)
		  {
                    fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
		  }
	      }
            else
	      {
                _ra[atom]->setIz(zz, -1);
	      }
	  }
      }
  }

  void Resp::addAtomSymmetry(std::vector<int> &symmetric_atoms)
  {
    int        i, k, zz;

    if (_ra.empty())
      {
        gmx_fatal(FARGS, "resp_atom struct not initialized");
      }

    /* Map the symmetric atoms */
    for (i = 0; (i < _natom); i++)
      {
        if (0 == i)
	  {
            /* The first charge is not a free variable, it follows from the total charge.
             * Only add the zeta values here.
             */
            for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	      {
                addParam( i, eparmZ, zz);
	      }
	  }
        else if (symmetric_atoms[i] == i)
	  {
            addParam( i, eparmQ, _ra[i]->getNZeta()-1);

            if (_bZatype)
	      {
                for (k = 0; (k < i); k++)
		  {
                    if (_ra[i]->getAtype() == _ra[k]->getAtype())
		      {
                        break;
		      }
		  }
                if (k == i)
		  {
                    for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		      {
                        addParam( i, eparmZ, zz);
		      }
		  }
                else
		  {
                    for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		      {
                        _ra[i]->setIz(zz, _ra[k]->getIz(zz));
		      }
		  }
	      }
            else
	      {
                for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		  {
                    addParam( i, eparmZ, zz);
		  }
	      }
	  }
        else if (symmetric_atoms[i] > i)
	  {
            gmx_fatal(FARGS, "The symmetric_atoms array can not point to larger atom numbers");
	  }
        else if (_ra[i]->getNZeta() > 0)
	  {
            _ra[i]->setIq(_ra[i]->getNZeta()-1,
			  _ra[symmetric_atoms[i]]->getIq(_ra[i]->getNZeta()-1));
            for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	      {
                _ra[i]->setIz(zz,_ra[symmetric_atoms[i]]->getIz(zz));
	      }
            if (debug)
	      {
                fprintf(debug, "Atom %d is a copy of atom %d\n",
                        i+1, symmetric_atoms[i]+1);
	      }
	  }
        else
	  {
            if (0 == i)
	      {
                _ra[i]->setIq(_ra[i]->getNZeta()-1, -1);
	      }
            else
	      {
                addParam( i, eparmQ, _ra[i]->getNZeta()-1);
	      }

            for (k = 0; (k < _ra[i]->getNZeta()); k++)
	      {
                addParam( i, eparmZ, k);
	      }
	  }
      }
    if (debug)
      {
        int maxz = 0;
        for (i = 0; (i < _natom); i++)
	  {
            if (_ra[i]->getNZeta() > maxz)
	      {
                maxz = _ra[i]->getNZeta();
	      }
	  }

        fprintf(debug, "GRQ: %3s %5s", "nr", "type");
        for (i = 0; (i < maxz); i++)
	  {
            fprintf(debug, " %8s %4s %8s %4s\n", "q", "iq", "zeta", "iz");
	  }
        fprintf(debug, "\n");
        for (i = 0; (i < _natom); i++)
	  {
            fprintf(debug, "GRQ: %3d %5s", i+1, _ra[i]->getAtomtype().c_str());
            for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	      {
                fprintf(debug, " %8.4f %4d %8.4f %4d\n",
                        _ra[i]->getQ(zz), _ra[i]->getIq(zz),
                        _ra[i]->getZeta(zz), _ra[i]->getIz(zz));
	      }
            fprintf(debug, "\n");
	  }
        fprintf(debug, "_qsum = %g\n", _qsum);
      }
  }

  void Resp::writeHisto( const char *fn, char *title, output_env_t oenv)
  {
    FILE       *fp;
    gmx_stats_t gs;
    real       *x, *y;
    int         i, nbin = 100;

    if (NULL == fn)
      {
        return;
      }
    gs = gmx_stats_init();
    for (i = 0; (i < _nesp); i++)
      {
        gmx_stats_add_point(gs, i, gmx2convert(_potCalc[i], eg2cHartree_e), 0, 0);
      }

    gmx_stats_make_histogram(gs, 0, &nbin, ehistoY, 1, &x, &y);

    fp = xvgropen(fn, title, "Pot (1/a.u.)", "()", oenv);
    for (i = 0; (i < nbin); i++)
      {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
      }
    sfree(x);
    sfree(y);
    fclose(fp);
    gmx_stats_done(gs);
  }


  //TODO  is src right or shood this be src
  void Resp::writeDiffCube(Resp * src,const char *cube_fn,
			     const char *hist_fn, char *title, output_env_t oenv,
			     int rho)
  {
    FILE       *fp;
    int         i, m, ix, iy, iz, zz;
    real        pp, q, r, rmin;
    rvec        dx;
    gmx_stats_t gst = NULL, ppcorr = NULL;

    if (NULL != hist_fn)
      {
        gst    = gmx_stats_init();
        ppcorr = gmx_stats_init();
      }
    if (NULL != cube_fn)
      {
        fp = gmx_ffopen(cube_fn, "w");
        fprintf(fp, "%s\n", title);
        fprintf(fp, "POTENTIAL\n");
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n",
                _natom,
                gmx2convert(_origin[XX], eg2cBohr),
                gmx2convert(_origin[YY], eg2cBohr),
                gmx2convert(_origin[ZZ], eg2cBohr));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[XX],
                gmx2convert(_space[XX], eg2cBohr), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[YY],
                0.0, gmx2convert(_space[YY], eg2cBohr), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[ZZ],
                0.0, 0.0, gmx2convert(_space[ZZ], eg2cBohr));

        for (m = 0; (m < _natom); m++)
	  {
            q = 0;
            for (zz = 0; (zz < _ra[m]->getNZeta()); zz++)
	      {
                q += _ra[m]->getQ(zz);
	      }
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    _ra[m]->getAtomnumber(), q,
                    gmx2convert(_x[m][XX], eg2cBohr),
                    gmx2convert(_x[m][YY], eg2cBohr),
                    gmx2convert(_x[m][ZZ], eg2cBohr));
	  }

        for (ix = m = 0; ix < _nxyz[XX]; ix++)
	  {
            for (iy = 0; iy < _nxyz[YY]; iy++)
	      {
                for (iz = 0; iz < _nxyz[ZZ]; iz++, m++)
		  {
                    if (NULL != src)
		      {
                        pp = _potCalc[m] - src->_pot[m];
                        if (NULL != ppcorr)
			  {
                            gmx_stats_add_point(ppcorr,
                                                gmx2convert(src->_pot[m], eg2cHartree_e),
                                                gmx2convert(_potCalc[m], eg2cHartree_e), 0, 0);
			  }
		      }
                    else
		      {
                        if (rho == 0)
			  {
                            pp = gmx2convert(_potCalc[m], eg2cHartree_e);
			  }
                        else
			  {
                            pp = _rho[m]*pow(BOHR2NM, 3);
			  }
		      }
                    fprintf(fp, "%13.5e", pp);
                    if (iz % 6 == 5)
		      {
                        fprintf(fp, "\n");
		      }
                    if (NULL != gst)
		      {
                        rmin = 1000;
                        /* Add point to histogram! */
                        for (i = 0; (i < _natom); i++)
			  {
                            rvec_sub(_x[i], _esp[m], dx);
                            r = norm(dx);
                            if (r < rmin)
			      {
                                rmin = r;
			      }
			  }
                        gmx_stats_add_point(gst, rmin, pp, 0, 0);
		      }
		  }
                if ((iz % 6) != 0)
		  {
                    fprintf(fp, "\n");
		  }
	      }
	  }
        fclose(fp);
      }
    if (NULL != gst)
      {
        int   nb = 0;
        real *x  = NULL, *y = NULL;

        fp = xvgropen(hist_fn, "Absolute deviation from QM", "Distance (nm)",
                      "Potential", oenv);
        gmx_stats_dump_xy(gst, fp);
        if (0)
	  {
            gmx_stats_make_histogram(gst, 0.01, &nb, ehistoX, 0, &x, &y);
            gmx_stats_done(gst);
            for (i = 0; (i < nb); i++)
	      {
                fprintf(fp, "%10g  %10g\n", x[i], y[i]);
	      }
            sfree(x);
            sfree(y);
	  }
        fclose(fp);
        fp = xvgropen("diff-pot.xvg", "Correlation between QM and Calc", "Pot (QM)",
                      "Pot (Calc)", oenv);
        gmx_stats_dump_xy(ppcorr, fp);
        fclose(fp);
      }
  }

  void Resp::writeCube(const char *fn, char *title)
  {
    writeDiffCube(NULL,  fn, NULL, title, NULL, 0);
  }

  void Resp::writeRho( const char *fn, char *title)
  {
    writeDiffCube(NULL,  fn, NULL, title, NULL, 1);
  }

  void Resp::readCube( const char *fn, bool bESPonly)
  {
    char         **strings;
    bool           bOK;
    double         lx, ly, lz, pp, qq;
    int            nlines, line = 0, m, ix, iy, iz, n, anr, nxyz[DIM];
    double         origin[DIM], space[DIM];
    const  char   *forms[] = {
      "%lf", "%*s%lf", "%*s%*s%lf", "%*s%*s%*s%lf",
      "%*s%*s%*s%*s%lf", "%*s%*s%*s%*s%*s%lf"
    };
    if (NULL == fn)
      {
        return;
      }

    nlines = get_file(fn, &strings);
    bOK    = (nlines > 100);
    if (bOK)
      {
        printf("%s\n", strings[line++]);
      }
    bOK = (line < nlines) && (strcmp(strings[line++], "POTENTIAL") != 0);
    if (bOK)
      {
        bOK = (line < nlines) && (4 == sscanf(strings[line++], "%d%lf%lf%lf",
                                              &n, &origin[XX], &origin[YY], &origin[ZZ]));
      }
    if (bOK && !bESPonly)
      {
        _natom      = n;
        _origin[XX] = origin[XX];
        _origin[YY] = origin[YY];
        _origin[ZZ] = origin[ZZ];
      }
    if (bOK)
      {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%lf",
                                              &nxyz[XX], &space[XX]));
      }
    if (bOK)
      {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%*s%lf",
                                              &nxyz[YY], &space[YY]));
      }
    if (bOK)
      {
        bOK = (line < nlines) && (2 == sscanf(strings[line++], "%d%*s%*s%lf",
                                              &nxyz[ZZ], &space[ZZ]));
      }
    if (bOK)
      {
        for (m = 0; (m < DIM); m++)
	  {
            _nxyz[m]  = nxyz[m];
            _space[m] = space[m];
	  }
        for (m = 0; (m < DIM); m++)
	  {
            _origin[m] = convert2gmx(_origin[m], eg2cBohr);
            _space[m]  = convert2gmx(_space[m], eg2cBohr);
	  }
      }
    if (bOK && ((line+_natom) < nlines))
      {
        snew(_x, _natom);
        for (m = 0; (m < _natom); m++)
	  {
            bOK = (5 == sscanf(strings[line++], "%d%lf%lf%lf%lf",
                               &anr, &qq, &lx, &ly, &lz));
            if (bOK)
	      {
                if (!bESPonly)
		  {
                    _ra[m]->setAtomnumber(anr);
                    if (_ra[m]->getNZeta() > 0)
		      {
                        _ra[m]->setQ(_ra[m]->getNZeta()-1, qq);
		      }
		  }
                _x[m][XX] = convert2gmx(lx, eg2cBohr);
                _x[m][YY] = convert2gmx(ly, eg2cBohr);
                _x[m][ZZ] = convert2gmx(lz, eg2cBohr);
	      }
	  }
      }
    if (bOK)
      {
        _nesp = _nxyz[XX]*_nxyz[YY]*_nxyz[ZZ];
        _pot.resize(_nesp);
	snew(_esp,_nesp);
        for (ix = m = 0; ix < _nxyz[XX]; ix++)
	  {
            for (iy = 0; iy < _nxyz[YY]; iy++)
	      {
                for (iz = 0; iz < _nxyz[ZZ]; iz++, m++)
		  {
                    _esp[m][XX] = _origin[XX] + ix*_space[XX];
                    _esp[m][YY] = _origin[YY] + iy*_space[YY];
                    _esp[m][ZZ] = _origin[ZZ] + iz*_space[ZZ];
                    bOK            = (1 == sscanf(strings[line], forms[iz % 6], &pp));
                    if (bOK)
		      {
                        _pot[m] = convert2gmx(pp, eg2cHartree_e);
		      }
                    if (iz % 6 == 5)
		      {
                        line++;
		      }
		  }
                if ((iz % 6) != 0)
		  {
                    line++;
		  }
	      }
	  }
      }
    bOK = (line == nlines);
    if (!bOK)
      {
        gmx_fatal(FARGS, "Error reading %s, line %d out of %d", fn, line, nlines);
      }

    for (m = 0; (m < nlines); m++)
      {
        sfree(strings[m]);
      }
    sfree(strings);
  }

  void Resp::copyGrid(Resp * src)
  {
    int m;

    for (m = 0; (m < DIM); m++)
      {
        _origin[m] = src->_origin[m];
        _space[m]  = src->_space[m];
        _nxyz[m]   = src->_nxyz[m];
      }
    _nesp = src->_nesp;
    snew(_esp, _nesp);
    _pot.resize(_nesp);
    _potCalc.resize(_nesp);
    for (m = 0; (m < _nesp); m++)
      {
        copy_rvec(src->_esp[m], _esp[m]);
      }
  }

  Resp * Resp::copy()
  {
    Resp * dest = new Resp();

    memcpy(dest, this, sizeof(*this));

    return dest;
  }

  void Resp::makeGrid( real spacing, matrix box, rvec x[])
  {
    int  i, j, k, m, n;
    rvec xyz;

    if (0 != _nesp)
      {
        fprintf(stderr, "Overwriting existing ESP grid\n");
      }
    if (0 <= spacing)
      {
        spacing = 0.1;
        fprintf(stderr, "spacing too small, setting it to %g\n", spacing);
      }
    snew(_x, _natom);
    for (i = 0; (i < _natom); i++)
      {
        copy_rvec(x[i], _x[i]);
      }
    _nesp = 1;
    for (m = 0; (m < DIM); m++)
      {
        _nxyz[m]  = 1+(int) (box[m][m]/spacing);
        _space[m] = box[m][m]/_nxyz[m];
        _nesp    *= _nxyz[m];
      }
    n = 0;
    snew(_esp, _nesp);
    _potCalc.resize(_nesp);
    for (i = 0; (i < _nxyz[XX]); i++)
      {
        xyz[XX] = (i-0.5*_nxyz[XX])*_space[XX];
        for (j = 0; (j < _nxyz[YY]); j++)
	  {
            xyz[YY] = (j-0.5*_nxyz[YY])*_space[YY];
            for (k = 0; (k < _nxyz[ZZ]); k++)
	      {
                xyz[ZZ] = (k-0.5*_nxyz[ZZ])*_space[ZZ];
                copy_rvec(xyz, _esp[n]);
                n++;
	      }
	  }
      }
  }

  void Resp::calcRho()
  {
    unsigned int  i, j, k;
    real r, z, V, vv, pi32;
    rvec dx;

    pi32 = pow(M_PI, -1.5);
    if ((int)_rho.size() < _nesp)
      {
        _rho.resize(_nesp);
      }
    for (i = 0; (i < _rho.size()); i++)
      {
        V = 0;
        for (j = 0; ((int)j < _natom); j++)
	  {
            vv = 0;
            rvec_sub(_esp[i], _x[j], dx);
            r = norm(dx);
            switch (_iDistributionModel)
	      {
	      case eqdBultinck:
	      case eqdAXp:
		return;
	      case eqdAXs:
		vv = 0;
		break;
	      case eqdYang:
	      case eqdRappe:
		vv = _ra[j]->getQ(0)*Nuclear_SS(r, _ra[j]->getRow(0),
					     _ra[j]->getZeta(0));
		break;
	      case eqdAXg:
		vv = 0;
		for (k = 0; ((int)k < _ra[j]->getNZeta()); k++)
		  {
		    z = _ra[j]->getZeta(k);
		    if (z > 0)
		      {
			vv -= (_ra[j]->getQ(k)*pi32*exp(-sqr(r*z))*
			       pow(z, 3));
		      }
		  }
		break;
	      default:
		gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %d!",
			  _iDistributionModel);
	      }
            V  += vv;
	  }
        _rho[i] = V;
      }
  }

  void Resp::calcPot()
  {
    int    i, j, k, m;
    double r, r2, dx, V, vv;

    for (i = 0; (i < _nesp); i++)
      {
        V = 0;
        for (j = 0; (j < _natom); j++)
	  {
            vv = 0;
            r2 = 0;
            for (m = 0; (m < DIM); m++)
	      {
                dx  = _esp[i][m]-_x[j][m];
                r2 += dx*dx;
	      }
            r = sqrt(r2);
            switch (_iDistributionModel)
	      {
	      case eqdBultinck:
	      case eqdAXp:
		if (r > 0.01)
		  {
		    vv = _ra[j]->getQ(0)/r;
		  }
		break;
	      case eqdAXs:
		vv = 0;
		for (k = 0; (k < _ra[j]->getNZeta()); k++)
		  {
		    vv += _ra[j]->getQ(k)*Nuclear_SS(r, _ra[j]->getRow(k),
						     _ra[j]->getZeta(k));
		  }
		break;
	      case eqdYang:
	      case eqdRappe:
		vv = _ra[j]->getQ(0)*Nuclear_SS(r, _ra[j]->getRow(0),
					     _ra[j]->getZeta(0));
		break;
	      case eqdAXg:
		vv = 0;
		for (k = 0; (k < _ra[j]->getNZeta()); k++)
		  {
		    vv += _ra[j]->getQ(k)*Nuclear_GG(r, _ra[j]->getZeta(k));
		  }
		break;
	      default:
		gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %s!",
			  alexandria::Poldata::getEemtypeName(_iDistributionModel));
	      }
            V  += vv;
	  }
        _potCalc[i] = V*ONE_4PI_EPS0;
      }
  }

  void Resp::warning(const char *fn, int line)
  {
    fprintf(stderr, "WARNING: It seems like you have two sets of ESP data in your file\n         %s\n", fn);
    fprintf(stderr, "         using the second set, starting at line %d\n", line);
  }

  const char *Resp::getStoichiometry()
  {
    return _stoichiometry.c_str();
  }

  void Resp::getSetVector(bool         bSet,
			    bool         bRandQ,
			    bool         bRandZeta,
			    unsigned int seed,
			    double      *nmx)
  {
    int       i, n, zz, zzz, nrest;
    double    qtot, dq, qi, zeta;
    gmx_rng_t rnd = NULL;

    if (bSet && (bRandQ || bRandZeta))
      {
        rnd = gmx_rng_init(seed);
      }
    _penalty = 0;
    n           = 0;
    qtot        = 0;
    nrest       = 0;
    for (i = 0; (i < _natom); i++)
      {
        if (bSet)
	  {
            /* First do charges */
            qi = 0;
            for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	      {
                if (_ra[i]->getIq(zz) == n)
		  {
                    if (_ra[i]->getQ(zz) == 0)
		      {
                        nmx[n] = -qi;
                        if (bRandQ)
			  {
                            nmx[n] += 0.2*(gmx_rng_uniform_real(rnd)-0.5);
			  }
		      }
                    else
		      {
                        nmx[n] = _ra[i]->getQ(zz);
		      }
                    n++;
		  }
                qi += _ra[i]->getQ(zz);
	      }
            /* Then do zeta */
            if (_bFitZeta)
	      {
                for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		  {
                    if (_ra[i]->getIz(zz) == n)
		      {
                        real zmin = _zmin;
                        real zmax = _zmax;

                        if ((_deltaZ > 0) && (_ra[i]->getBRestrained()))
			  {
                            zmax = _ra[i]->getZetaRef(zz)+_deltaZ;
                            zmin = _ra[i]->getZetaRef(zz)-_deltaZ;
			  }
                        if ((zz > 1) && (_rDecrZeta >= 0))
			  {
                            zmax = _ra[i]->getZeta(zz-1)-_rDecrZeta;
                            if (zmax < zmin)
			      {
                                zmax = (zmin+_ra[i]->getZeta(zz-1))/2;
			      }
			  }
                        if (bRandZeta)
			  {
                            nmx[n] = zmin + (zmax-zmin)*gmx_rng_uniform_real(rnd);
			  }
                        else
			  {
                            nmx[n] = _ra[i]->getZeta(zz);
			  }
                        _ra[i]->setZeta(zz, nmx[n]);
                        n++;
		      }
		  }
	      }
	  }
        else
	  {
            /* Initialize to something strange */
            if (_bFitZeta)
	      {
                for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		  {
                    if (_ra[i]->getZeta(zz) != 0)
		      {
                        _ra[i]->setZeta(zz, NOTSET);
		      }
		  }
	      }
            _ra[i]->setQ(_ra[i]->getNZeta()-1, NOTSET);

            /* First do charges */
            for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	      {
                if (_ra[i]->getIq(zz) == n)
		  {
                    _ra[i]->setQ(zz, nmx[n]);
                    qtot           += _ra[i]->getQ(zz);
                    n++;
		  }
                else if ((_ra[i]->getIq(zz) < n) && (_ra[i]->getIq(zz) >= 0))
		  {
                    for (zzz = 0; (zzz < i); zzz++)
		      {
                        if (_ra[zzz]->getIq(zz) == _ra[i]->getIq(zz))
			  {
                            _ra[i]->setQ(zz, _ra[zzz]->getQ(zz));
                            break;
			  }
		      }
                    if (zzz == i)
		      {
                        gmx_fatal(FARGS, "Can not find a previous atom with iq[%d] = %d", zz, n);
		      }

                    /* Only sum those atoms to qtot, that are not part of
                       the "rest" charge */
                    if (_ra[i]->getIq(zz) != -1)
		      {
                        qtot += _ra[i]->getQ(zz);
		      }
		  }
                else if (zz == _ra[i]->getNZeta()-1)
		  {
                    nrest++;
		  }
                else
		  {
                    qtot += _ra[i]->getQ(zz);
		  }
	      }

            if (_bFitZeta)
	      {
                /* Then do zeta */
                for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
		  {
                    if (_ra[i]->getIz(zz) == n)
		      {
                        zeta               = nmx[n];
                        _ra[i]->setZeta(zz, zeta);
                        if (_deltaZ >= 0)
			  {
                            real zmin = _ra[i]->getZetaRef(zz)-_deltaZ;
                            real zmax = _ra[i]->getZetaRef(zz)+_deltaZ;
                            if (zeta <= zmin)
			      {
                                _penalty += sqr(zeta-zmin);
			      }
                            else if (zeta >= zmax)
			      {
                                _penalty += sqr(zmax-zeta);
			      }
			  }
                        else
			  {
                            if (zeta <= _zmin)
			      {
                                _penalty += sqr(_zmin-zeta);
			      }
                            else if (zeta >= _zmax)
			      {
                                _penalty += sqr(_zmax-zeta);
			      }
                            if ((_rDecrZeta >= 0) && (zz > 0) &&
                                (_ra[i]->getZeta(zz-1) != 0) &&
                                ((_ra[i]->getZeta(zz-1) - zeta) < _rDecrZeta))
			      {
                                _penalty += sqr(_ra[i]->getZeta(zz-1) - zeta - _rDecrZeta);
			      }
			  }
                        n++;
		      }
                    else if ((_ra[i]->getIz(zz) < n) && (_ra[i]->getIz(zz) >= 0))
		      {
                        for (zzz = 0; (zzz < i); zzz++)
			  {
                            if (_ra[zzz]->getIz(zz) == _ra[i]->getIz(zz))
			      {
                                _ra[i]->setZeta(zz, _ra[zzz]->getZeta(zz));
                                break;
			      }
			  }
                        if (zzz == i)
			  {
                            gmx_fatal(FARGS, "Can not find a previous atom with iz[%d] = %d", zz, n);
			  }
		      }
                    else if ((_ra[i]->getIz(zz) == -1) && (_ra[i]->getZeta(zz) != 0))
		      {
                        gmx_fatal(FARGS, "ra[%d]->iz[%d] = %d whereas ra[%d]->zeta[%d] = %g", i, zz, _ra[i]->getIz(zz), i, zz, _ra[i]->getZeta(zz));
		      }
		  }
	      }
	  }
      }
    if (NULL != rnd)
      {
        gmx_rng_destroy(rnd);
      }
    if (n != _nparam)
      {
        gmx_fatal(FARGS, "Whoopsydaisies! n = %d, should be %d. bSet = %d", n, _nparam, bSet);
      }

    if (nrest > 0)
      {
        dq = (_qtot-qtot)/nrest;
        if (debug)
	  {
            fprintf(debug, "_qtot = %g, qtot = %g, nrest = %d, dq = %g\n",
                    _qtot, qtot, nrest, dq);
	  }
        for (i = 0; (i < _natom); i++)
	  {
            if (_ra[i]->getIq(_ra[i]->getNZeta()-1) == -1)
	      {
                _ra[i]->setQ(_ra[i]->getNZeta()-1,dq);
	      }
	  }
      }
    /* Check for excessive charges */
    for (i = 0; (i < _natom); i++)
      {
        qi = 0;
        for (zz = 0; (zz < _ra[i]->getNZeta()); zz++)
	  {
            qi += _ra[i]->getQ(zz);
	  }
        if (qi < _qmin)
	  {
            _penalty += sqr(_qmin-qi);
	  }
        else if (qi > _qmax)
	  {
            _penalty += sqr(_qmax-qi);
	  }
        else if ((qi < -0.02) && (_ra[i]->getAtomnumber() == 1))
	  {
            _penalty += qi*qi;
	  }
      }
    _penalty *= _pfac;
  }

  void Resp::addPoint( double x, double y,
                        double z, double V)
  {
    int i;

    i = _nesp++;
    snew(_esp,_nesp);
    _pot.resize(_nesp);
    _potCalc.resize(_nesp);
    _esp[i][XX]  = x;
    _esp[i][YY]  = y;
    _esp[i][ZZ]  = z;
    _pot[i]      = V;
    _potCalc[i] = 0;
  }

  real Resp::myWeight( int iatom)
  {
    if (iatom < _natom)
      {
        return _watoms;
      }
    else
      {
        return 1.0;
      }
  }

  void Resp::potLsq( gmx_stats_t lsq)
  {
    int    i;
    double w;

    for (i = 0; (i < _nesp); i++)
      {
        w = myWeight( i);
        if (w > 0)
	  {
            gmx_stats_add_point(lsq,
                                gmx2convert(_pot[i], eg2cHartree_e),
                                gmx2convert(_potCalc[i], eg2cHartree_e), 0, 0);
	  }
      }
  }

  void Resp::calcRms()
  {
    int    i;
    double pot2, s2, sum2, w, wtot, entropy;
    char   buf[STRLEN];

    pot2 = sum2 = wtot = entropy = 0;
    sprintf(buf, " - weight %g in fit", _watoms);
    for (i = 0; (i < _nesp); i++)
      {
        w = myWeight( i);
        if ((NULL != debug) && (i < 2*_natom))
	  {
            fprintf(debug, "ESP %d QM: %g EEM: %g DIFF: %g%s\n",
                    i, _pot[i], _potCalc[i],
                    _pot[i]-_potCalc[i],
                    (i < _natom)  ? buf : "");
	  }
        s2    = w*sqr(_pot[i]-_potCalc[i]);
        if ((s2 > 0) && (_bEntropy))
	  {
            entropy += s2*log(s2);
	  }
        sum2 += s2;
        pot2 += w*sqr(_pot[i]);
        wtot += w;
      }
    _wtot = wtot;
    if (wtot > 0)
      {
        _rms     = gmx2convert(sqrt(sum2/wtot), eg2cHartree_e);
        _entropy = gmx2convert(entropy/wtot, eg2cHartree_e);
      }
    else
      {
        _rms     = 0;
        _entropy = 0;
      }
    _rrms = sqrt(sum2/pot2);
  }

  double Resp::getRms( real *wtot)
  {
    calcRms();
    *wtot = _wtot;
    if (_bEntropy)
      {
        return _entropy;
      }
    else
      {
        return _rms;
      }
  }

  void Resp::calcPenalty()
  {
    int    i;
    double p, b2;

    p = 0;
    if (_bAXpRESP && (_iDistributionModel == eqdAXp))
      {
        b2 = sqr(_bHyper);
        for (i = 0; (i < _natom); i++)
	  {
            p += sqrt(sqr(_ra[i]->getQ(0)) + b2) - _bHyper;
	  }
        p = (_qfac * p);
      }
    _penalty += p;
  }

  //Writen in c stile, needed as an function argument
  double Resp::chargeFunction(void * gr,double v[])
  {
    Resp * resp = (Resp *)gr;
    double     rms = 0;
    real       wtot;

    resp->getSetVector( false, false, false, resp->_seed, v);
    resp->calcPot();
    resp->calcPenalty();
    rms = resp->getRms( &wtot);

    return rms; // + _penalty;
  }

  void Resp::statistics( int len, char buf[])
  {
    if (len >= 100)
      {
        sprintf(buf, "RMS: %10e [Hartree/e] RRMS: %10e Entropy: %10e Penalty: %10e",
                _rms, _rrms, _entropy, _penalty);
      }
    else
      {
        fprintf(stderr, "buflen too small (%d) in gmx_resp_statistics\n", len);
      }
  }

  int Resp::optimizeCharges(FILE *fp,  int maxiter,
			     real toler, real *rms)
  {
    double *param;
    double  ccc;
    int     bConv;
    char    buf[STRLEN];

    snew(param, _nparam);

    getSetVector( true, _bRandQ, _bRandZeta, _seed, param);

    bConv = nmsimplex(fp, (void *)this, chargeFunction, param, _nparam,
                      toler, 1, maxiter, &ccc);
    if (bConv)
      {
        statistics( STRLEN-1, buf);
      }
    else
      {
        printf("NM Simplex did not converge\n\n");
      }

    if (_bEntropy)
      {
        *rms = _entropy;
      }
    else
      {
        *rms = _rms;
      }

    getSetVector( false, false, false, _seed, param);

    sfree(param);

    if (bConv)
      {
        return eQGEN_OK;
      }
    else
      {
        return eQGEN_NOTCONVERGED;
      }
  }


  void Resp::potcomp( const char *potcomp,
                      const char *pdbdiff, output_env_t oenv)
  {
    int     i;
    double  pp, exp, eem;
    FILE   *fp;
    int     unit = eg2cHartree_e;

    if (NULL != potcomp)
      {
        const char *pcleg[2] = { "Atoms", "ESP points" };
        fp = xvgropen(potcomp, "Electrostatic potential", unit2string(unit), unit2string(unit), oenv);
        xvgr_legend(fp, 2, pcleg, oenv);
        fprintf(fp, "@type xy\n");
        for (i = 0; (i < _nesp); i++)
	  {
            /* Conversion may or may not be in vain depending on unit */
            exp = gmx2convert(_pot[i], unit);
            eem = gmx2convert(_potCalc[i], unit);
            if (i == _natom)
	      {
                fprintf(fp, "&\n");
                fprintf(fp, "@type xy\n");
	      }
            fprintf(fp, "%10.5e  %10.5e\n", exp, eem);
	  }
        fprintf(fp, "&\n");
        fclose(fp);
      }
    if (NULL != pdbdiff)
      {
        fp = fopen(pdbdiff, "w");
        fprintf(fp, "REMARK All distances are scaled by a factor of two.\n");
        for (i = 0; (i < _nesp); i++)
	  {
            exp = gmx2convert(_pot[i], eg2cHartree_e);
            eem = gmx2convert(_potCalc[i], eg2cHartree_e);
            pp  = _pot[i]-_potCalc[i];
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, "HE", "HE", ' ', i+1, ' ', 20*_esp[i][XX],
                    20*_esp[i][YY], 20*_esp[i][ZZ], 0.0, pp);
	  }
        fclose(fp);
      }
  }

  double Resp::getQtot( int atom)
  {
    int    i;
    double q = 0;

    range_check(atom, 0, _natom);
    for (i = 0; (i < _ra[atom]->getNZeta()); i++)
      {
        q += _ra[atom]->getQ(i);
      }
    return q;
  }

  double Resp::getQ( int atom, int zz)
  {
    range_check(atom, 0, _natom);
    range_check(zz, 0, _ra[atom]->getNZeta());

    return _ra[atom]->getQ(zz);
  }

  double Resp::getZeta( int atom, int zz)
  {
    range_check(atom, 0, _natom);
    range_check(zz, 0, _ra[atom]->getNZeta());

    return _ra[atom]->getZeta(zz);
  }

  void Resp::setQ( int atom, int zz, double q)
  {
    range_check(atom, 0, _natom);
    range_check(zz, 0, _ra[atom]->getNZeta());

    _ra[atom]->setQ(zz,q);
  }

  void Resp::setZeta( int atom, int zz, double zeta)
  {
    range_check(atom, 0, _natom);
    range_check(zz, 0, _ra[atom]->getNZeta());

    _ra[atom]->setZeta(zz,zeta);
  }
}
