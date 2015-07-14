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
	     real zmin, real zmax, real delta_z, bool bZatype,
	     real watoms, real rDecrZeta, bool bRandZeta,
	     bool bRandQ, real penalty_fac, bool bFitZeta, bool bEntropy,
	     const char *dzatoms, unsigned int seed)
  {



    this->qtot                  = qtot;
    this->qsum                  = qtot;
    this->bAXpRESP              = bAXpRESP;
    this->qfac                  = qfac;
    this->b_hyper               = b_hyper;
    this->wtot                  = 0;
    this->iDistributionModel    = iDistributionModel;
    this->zmin                  = zmin;
    this->zmax                  = zmax;
    this->delta_z               = delta_z;
    std::vector<std::string> ptr = split(dzatoms, ' ');
    snew(this->dzatoms, ptr.size());
    for (unsigned int i = 0; (i < ptr.size()); ++i)
      {
        this->dzatoms[i] = strdup(ptr[i].c_str());
      }
    this->pfac      = penalty_fac;
    this->qmin      = -2;
    this->qmax      = 2; /* e */
    this->nesp      = 0;
    this->nrho      = 0;
    this->natom     = 0;
    this->natype    = 0;
    this->seed      = seed;
    this->bEntropy  = bEntropy;
    this->bZatype   = bZatype;
    this->rDecrZeta = rDecrZeta;
    this->bRandZeta = bRandZeta;
    this->bRandQ    = bRandQ;
    this->bFitZeta  = bFitZeta;
    this->watoms    = watoms;
    this->nparam    = 0;

  
  }

  void Resp::get_atom_info( t_atoms *atoms,
                            t_symtab *symtab, rvec **x)
  {
    int          i;
    const char  *rnm;

    init_t_atoms(atoms, this->natom, true);
    if (NULL == (*x))
      {
        snew((*x), atoms->nr);
      }
    if (NULL == atoms->atomtype)
      {
        snew(atoms->atomtype, atoms->nr);
      }
    for (i = 0; (i < this->natom); i++)
      {
        atoms->atom[i].atomnumber = this->ra[i]->atomnumber;
        atoms->atom[i].q          = this->ra[i]->get_q();
        atoms->atomname[i]        = put_symtab(symtab, this->ra[i]->atomtype);
        atoms->atomtype[i]        = put_symtab(symtab, this->ra[i]->atomtype);
        atoms->atom[i].resind     = 0;

        strncpy(atoms->atom[i].elem, this->ra[i]->atomtype,
                sizeof(atoms->atom[i].elem)-1);
        copy_rvec(this->x[i], (*x)[i]);
      }
    rnm = (NULL != this->stoichiometry) ? this->stoichiometry : (const char *)"BOE";
    t_atoms_set_resinfo(atoms, 0, symtab, rnm, 1, ' ', 1, ' ');

    atoms->nres = 1;
  }

  void Resp::update_atomtypes( t_atoms *atoms)
  {
    int i, j;

    for (i = 0; (i < this->natom); i++)
      {
        this->ra[i]->atomtype = strdup(*atoms->atomtype[i]);
        for (j = 0; (j < i); j++)
	  {
            if (0 == strcmp(*atoms->atomtype[i], *atoms->atomtype[j]))
	      {
                break;
	      }
	  }
        if (j == i)
	  {
            this->natype++;
	  }
        this->ra[i]->atype = j;
      }
  }

  void Resp::add_atom_coords( rvec *x)
  {
    int        i;

    srenew(this->x, this->natom);
    for (i = 0; (i < this->natom); i++)
      {
        copy_rvec(x[i], this->x[i]);
      }
  }

  void Resp::fill_zeta( alexandria::Poldata * pd)
  {
    int i, zz;

    for (i = 0; (i < this->natom); i++)
      {
        for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	  {
            Resp::set_zeta( i, zz, pd->getZeta( this->iDistributionModel,
						 this->ra[i]->atomtype, zz));
	  }
      }
  }

  void Resp::fill_q( t_atoms *atoms)
  {
    int    i, zz;
    double q;

    for (i = 0; (i < this->natom); i++)
      {
        q = 0;
        for (zz = 0; (zz < this->ra[i]->nZeta-1); zz++)
	  {
            q -= this->ra[i]->q[zz];
	  }
        q += atoms->atom[i].q;
        Resp::set_q( i, this->ra[i]->nZeta-1, q);
      }
  }

  bool Resp::add_atom_info( t_atoms *atoms, alexandria::Poldata * pd)
  {
    int  i;

    this->natom    = atoms->nr;
    snew(this->ra, this->natom);

    for (i = 0; (i < this->natom); i++)
      {
	this->ra[i] =  new Ra(atoms->atom[i].atomnumber, atoms->atom[i].type,
			     *(atoms->atomtype[i]), pd, this->iDistributionModel, this->dzatoms);
        if (this->ra[i]->setUpcorrectly() == false)
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
                this->natom, this->natype, this->nparam);
        for (i = 0; (i < this->natom); i++)
	  {
            fprintf(fp, " %d", symmetric_atoms[i]);
	  }
        fprintf(fp, "\n");
      }
  }



  void Resp::add_param( int atom, eParm eparm, int zz)
  {
    range_check(atom, 0, this->natom);
    if ((zz >= 0) && (zz < this->ra[atom]->nZeta))
      {
        if (eparm == eparmQ)
	  {
            this->ra[atom]->iq[zz] = this->nparam++;
            if (debug)
	      {
                fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
	      }
	  }
        else if (this->bFitZeta)
	  {
            if (this->ra[atom]->zeta[zz] != 0)
	      {
                this->ra[atom]->iz[zz] = this->nparam++;
                if (debug)
		  {
                    fprintf(debug, "GRESP: Adding parameter %d for atom %d zz %d\n", eparm, atom, zz);
		  }
	      }
            else
	      {
                this->ra[atom]->iz[zz] = -1;
	      }
	  }
      }
  }

  void Resp::add_atom_symmetry(std::vector<int> &symmetric_atoms)
  {
    int        i, k, zz;

    if (NULL == this->ra)
      {
        gmx_fatal(FARGS, "resp_atom struct not initialized");
      }

    /* Map the symmetric atoms */
    for (i = 0; (i < this->natom); i++)
      {
        if (0 == i)
	  {
            /* The first charge is not a free variable, it follows from the total charge.
             * Only add the zeta values here.
             */
            for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	      {
                add_param( i, eparmZ, zz);
	      }
	  }
        else if (symmetric_atoms[i] == i)
	  {
            add_param( i, eparmQ, this->ra[i]->nZeta-1);

            if (this->bZatype)
	      {
                for (k = 0; (k < i); k++)
		  {
                    if (this->ra[i]->atype == this->ra[k]->atype)
		      {
                        break;
		      }
		  }
                if (k == i)
		  {
                    for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		      {
                        add_param( i, eparmZ, zz);
		      }
		  }
                else
		  {
                    for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		      {
                        this->ra[i]->iz[zz] = this->ra[k]->iz[zz];
		      }
		  }
	      }
            else
	      {
                for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		  {
                    add_param( i, eparmZ, zz);
		  }
	      }
	  }
        else if (symmetric_atoms[i] > i)
	  {
            gmx_fatal(FARGS, "The symmetric_atoms array can not point to larger atom numbers");
	  }
        else if (this->ra[i]->nZeta > 0)
	  {
            this->ra[i]->iq[this->ra[i]->nZeta-1] =
	      this->ra[symmetric_atoms[i]]->iq[this->ra[i]->nZeta-1];
            for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	      {
                this->ra[i]->iz[zz] = this->ra[symmetric_atoms[i]]->iz[zz];
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
                this->ra[i]->iq[this->ra[i]->nZeta-1] = -1;
	      }
            else
	      {
                add_param( i, eparmQ, this->ra[i]->nZeta-1);
	      }

            for (k = 0; (k < this->ra[i]->nZeta); k++)
	      {
                add_param( i, eparmZ, k);
	      }
	  }
      }
    if (debug)
      {
        int maxz = 0;
        for (i = 0; (i < this->natom); i++)
	  {
            if (this->ra[i]->nZeta > maxz)
	      {
                maxz = this->ra[i]->nZeta;
	      }
	  }

        fprintf(debug, "GRQ: %3s %5s", "nr", "type");
        for (i = 0; (i < maxz); i++)
	  {
            fprintf(debug, " %8s %4s %8s %4s\n", "q", "iq", "zeta", "iz");
	  }
        fprintf(debug, "\n");
        for (i = 0; (i < this->natom); i++)
	  {
            fprintf(debug, "GRQ: %3d %5s", i+1, this->ra[i]->atomtype);
            for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	      {
                fprintf(debug, " %8.4f %4d %8.4f %4d\n",
                        this->ra[i]->q[zz], this->ra[i]->iq[zz],
                        this->ra[i]->zeta[zz], this->ra[i]->iz[zz]);
	      }
            fprintf(debug, "\n");
	  }
        fprintf(debug, "this->qsum = %g\n", this->qsum);
      }
  }

  void Resp::write_histo( const char *fn, char *title, output_env_t oenv)
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
    for (i = 0; (i < this->nesp); i++)
      {
        gmx_stats_add_point(gs, i, gmx2convert(this->pot_calc[i], eg2cHartree_e), 0, 0);
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
  void Resp::write_diff_cube(Resp * src,const char *cube_fn,
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
                this->natom,
                gmx2convert(this->origin[XX], eg2cBohr),
                gmx2convert(this->origin[YY], eg2cBohr),
                gmx2convert(this->origin[ZZ], eg2cBohr));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", this->nxyz[XX],
                gmx2convert(this->space[XX], eg2cBohr), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", this->nxyz[YY],
                0.0, gmx2convert(this->space[YY], eg2cBohr), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", this->nxyz[ZZ],
                0.0, 0.0, gmx2convert(this->space[ZZ], eg2cBohr));

        for (m = 0; (m < this->natom); m++)
	  {
            q = 0;
            for (zz = 0; (zz < this->ra[m]->nZeta); zz++)
	      {
                q += this->ra[m]->q[zz];
	      }
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    this->ra[m]->atomnumber, q,
                    gmx2convert(this->x[m][XX], eg2cBohr),
                    gmx2convert(this->x[m][YY], eg2cBohr),
                    gmx2convert(this->x[m][ZZ], eg2cBohr));
	  }

        for (ix = m = 0; ix < this->nxyz[XX]; ix++)
	  {
            for (iy = 0; iy < this->nxyz[YY]; iy++)
	      {
                for (iz = 0; iz < this->nxyz[ZZ]; iz++, m++)
		  {
                    if (NULL != src)
		      {
                        pp = this->pot_calc[m] - src->pot[m];
                        if (NULL != ppcorr)
			  {
                            gmx_stats_add_point(ppcorr,
                                                gmx2convert(src->pot[m], eg2cHartree_e),
                                                gmx2convert(this->pot_calc[m], eg2cHartree_e), 0, 0);
			  }
		      }
                    else
		      {
                        if (rho == 0)
			  {
                            pp = gmx2convert(this->pot_calc[m], eg2cHartree_e);
			  }
                        else
			  {
                            pp = this->rho[m]*pow(BOHR2NM, 3);
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
                        for (i = 0; (i < this->natom); i++)
			  {
                            rvec_sub(this->x[i], this->esp[m], dx);
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

  void Resp::write_cube(const char *fn, char *title)
  {
    write_diff_cube(NULL,  fn, NULL, title, NULL, 0);
  }

  void Resp::write_rho( const char *fn, char *title)
  {
    write_diff_cube(NULL,  fn, NULL, title, NULL, 1);
  }

  void Resp::read_cube( const char *fn, bool bESPonly)
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
        this->natom      = n;
        this->origin[XX] = origin[XX];
        this->origin[YY] = origin[YY];
        this->origin[ZZ] = origin[ZZ];
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
            this->nxyz[m]  = nxyz[m];
            this->space[m] = space[m];
	  }
        for (m = 0; (m < DIM); m++)
	  {
            this->origin[m] = convert2gmx(this->origin[m], eg2cBohr);
            this->space[m]  = convert2gmx(this->space[m], eg2cBohr);
	  }
      }
    if (bOK && ((line+this->natom) < nlines))
      {
        snew(this->x, this->natom);
        for (m = 0; (m < this->natom); m++)
	  {
            bOK = (5 == sscanf(strings[line++], "%d%lf%lf%lf%lf",
                               &anr, &qq, &lx, &ly, &lz));
            if (bOK)
	      {
                if (!bESPonly)
		  {
                    this->ra[m]->atomnumber = anr;
                    if (this->ra[m]->nZeta > 0)
		      {
                        this->ra[m]->q[this->ra[m]->nZeta-1] = qq;
		      }
		  }
                this->x[m][XX] = convert2gmx(lx, eg2cBohr);
                this->x[m][YY] = convert2gmx(ly, eg2cBohr);
                this->x[m][ZZ] = convert2gmx(lz, eg2cBohr);
	      }
	  }
      }
    if (bOK)
      {
        this->nesp = this->nxyz[XX]*this->nxyz[YY]*this->nxyz[ZZ];
        snew(this->pot, this->nesp);
        snew(this->esp, this->nesp);
        for (ix = m = 0; ix < this->nxyz[XX]; ix++)
	  {
            for (iy = 0; iy < this->nxyz[YY]; iy++)
	      {
                for (iz = 0; iz < this->nxyz[ZZ]; iz++, m++)
		  {
                    this->esp[m][XX] = this->origin[XX] + ix*this->space[XX];
                    this->esp[m][YY] = this->origin[YY] + iy*this->space[YY];
                    this->esp[m][ZZ] = this->origin[ZZ] + iz*this->space[ZZ];
                    bOK            = (1 == sscanf(strings[line], forms[iz % 6], &pp));
                    if (bOK)
		      {
                        this->pot[m] = convert2gmx(pp, eg2cHartree_e);
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

  void Resp::copy_grid(Resp * src)
  {
    int m;

    for (m = 0; (m < DIM); m++)
      {
        this->origin[m] = src->origin[m];
        this->space[m]  = src->space[m];
        this->nxyz[m]   = src->nxyz[m];
      }
    this->nesp = src->nesp;
    snew(this->esp, this->nesp);
    snew(this->pot, this->nesp);
    snew(this->pot_calc, this->nesp);
    for (m = 0; (m < this->nesp); m++)
      {
        copy_rvec(src->esp[m], this->esp[m]);
      }
  }

  Resp * Resp::copy()
  {
    Resp * dest = new Resp();

    memcpy(dest, this, sizeof(*this));

    return dest;
  }

  void Resp::make_grid( real spacing, matrix box, rvec x[])
  {
    int  i, j, k, m, n;
    rvec xyz;

    if (0 != this->nesp)
      {
        fprintf(stderr, "Overwriting existing ESP grid\n");
      }
    if (0 <= spacing)
      {
        spacing = 0.1;
        fprintf(stderr, "spacing too small, setting it to %g\n", spacing);
      }
    snew(this->x, this->natom);
    for (i = 0; (i < this->natom); i++)
      {
        copy_rvec(x[i], this->x[i]);
      }
    this->nesp = 1;
    for (m = 0; (m < DIM); m++)
      {
        this->nxyz[m]  = 1+(int) (box[m][m]/spacing);
        this->space[m] = box[m][m]/this->nxyz[m];
        this->nesp    *= this->nxyz[m];
      }
    n = 0;
    snew(this->esp, this->nesp);
    snew(this->pot_calc, this->nesp);
    for (i = 0; (i < this->nxyz[XX]); i++)
      {
        xyz[XX] = (i-0.5*this->nxyz[XX])*this->space[XX];
        for (j = 0; (j < this->nxyz[YY]); j++)
	  {
            xyz[YY] = (j-0.5*this->nxyz[YY])*this->space[YY];
            for (k = 0; (k < this->nxyz[ZZ]); k++)
	      {
                xyz[ZZ] = (k-0.5*this->nxyz[ZZ])*this->space[ZZ];
                copy_rvec(xyz, this->esp[n]);
                n++;
	      }
	  }
      }
  }

  void Resp::calc_rho()
  {
    int  i, j, k;
    real r, z, V, vv, pi32;
    rvec dx;

    pi32 = pow(M_PI, -1.5);
    if (this->nrho < this->nesp)
      {
        srenew(this->rho, this->nesp);
        this->nrho = this->nesp;
      }
    for (i = 0; (i < this->nrho); i++)
      {
        V = 0;
        for (j = 0; (j < this->natom); j++)
	  {
            vv = 0;
            rvec_sub(this->esp[i], this->x[j], dx);
            r = norm(dx);
            switch (this->iDistributionModel)
	      {
	      case eqdBultinck:
	      case eqdAXp:
		return;
	      case eqdAXs:
		vv = 0;
		break;
	      case eqdYang:
	      case eqdRappe:
		vv = this->ra[j]->q[0]*Nuclear_SS(r, this->ra[j]->row[0],
						 this->ra[j]->zeta[0]);
		break;
	      case eqdAXg:
		vv = 0;
		for (k = 0; (k < this->ra[j]->nZeta); k++)
		  {
		    z = this->ra[j]->zeta[k];
		    if (z > 0)
		      {
			vv -= (this->ra[j]->q[k]*pi32*exp(-sqr(r*z))*
			       pow(z, 3));
		      }
		  }
		break;
	      default:
		gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %d!",
			  this->iDistributionModel);
	      }
            V  += vv;
	  }
        this->rho[i] = V;
      }
  }

  void Resp::calc_pot()
  {
    int    i, j, k, m;
    double r, r2, dx, V, vv;

    for (i = 0; (i < this->nesp); i++)
      {
        V = 0;
        for (j = 0; (j < this->natom); j++)
	  {
            vv = 0;
            r2 = 0;
            for (m = 0; (m < DIM); m++)
	      {
                dx  = this->esp[i][m]-this->x[j][m];
                r2 += dx*dx;
	      }
            r = sqrt(r2);
            switch (this->iDistributionModel)
	      {
	      case eqdBultinck:
	      case eqdAXp:
		if (r > 0.01)
		  {
		    vv = this->ra[j]->q[0]/r;
		  }
		break;
	      case eqdAXs:
		vv = 0;
		for (k = 0; (k < this->ra[j]->nZeta); k++)
		  {
		    vv += this->ra[j]->q[k]*Nuclear_SS(r, this->ra[j]->row[k],
						      this->ra[j]->zeta[k]);
		  }
		break;
	      case eqdYang:
	      case eqdRappe:
		vv = this->ra[j]->q[0]*Nuclear_SS(r, this->ra[j]->row[0],
						 this->ra[j]->zeta[0]);
		break;
	      case eqdAXg:
		vv = 0;
		for (k = 0; (k < this->ra[j]->nZeta); k++)
		  {
		    vv += this->ra[j]->q[k]*Nuclear_GG(r, this->ra[j]->zeta[k]);
		  }
		break;
	      default:
		gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %s!",
			  alexandria::Poldata::getEemtypeName(this->iDistributionModel));
	      }
            V  += vv;
	  }
        this->pot_calc[i] = V*ONE_4PI_EPS0;
      }
  }

  void Resp::warning(const char *fn, int line)
  {
    fprintf(stderr, "WARNING: It seems like you have two sets of ESP data in your file\n         %s\n", fn);
    fprintf(stderr, "         using the second set, starting at line %d\n", line);
  }

  const char *Resp::get_stoichiometry()
  {
    return this->stoichiometry;
  }

  void Resp::get_set_vector(bool         bSet,
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
    this->penalty = 0;
    n           = 0;
    qtot        = 0;
    nrest       = 0;
    for (i = 0; (i < this->natom); i++)
      {
        if (bSet)
	  {
            /* First do charges */
            qi = 0;
            for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	      {
                if (this->ra[i]->iq[zz] == n)
		  {
                    if (this->ra[i]->q[zz] == 0)
		      {
                        nmx[n] = -qi;
                        if (bRandQ)
			  {
                            nmx[n] += 0.2*(gmx_rng_uniform_real(rnd)-0.5);
			  }
		      }
                    else
		      {
                        nmx[n] = this->ra[i]->q[zz];
		      }
                    n++;
		  }
                qi += this->ra[i]->q[zz];
	      }
            /* Then do zeta */
            if (this->bFitZeta)
	      {
                for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		  {
                    if (this->ra[i]->iz[zz] == n)
		      {
                        real zmin = this->zmin;
                        real zmax = this->zmax;

                        if ((this->delta_z > 0) && (this->ra[i]->bRestrained))
			  {
                            zmax = this->ra[i]->zeta_ref[zz]+this->delta_z;
                            zmin = this->ra[i]->zeta_ref[zz]-this->delta_z;
			  }
                        if ((zz > 1) && (this->rDecrZeta >= 0))
			  {
                            zmax = this->ra[i]->zeta[zz-1]-this->rDecrZeta;
                            if (zmax < zmin)
			      {
                                zmax = (zmin+this->ra[i]->zeta[zz-1])/2;
			      }
			  }
                        if (bRandZeta)
			  {
                            nmx[n] = zmin + (zmax-zmin)*gmx_rng_uniform_real(rnd);
			  }
                        else
			  {
                            nmx[n] = this->ra[i]->zeta[zz];
			  }
                        this->ra[i]->zeta[zz] = nmx[n];
                        n++;
		      }
		  }
	      }
	  }
        else
	  {
            /* Initialize to something strange */
            if (this->bFitZeta)
	      {
                for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		  {
                    if (this->ra[i]->zeta[zz] != 0)
		      {
                        this->ra[i]->zeta[zz] = NOTSET;
		      }
		  }
	      }
            this->ra[i]->q[this->ra[i]->nZeta-1] = NOTSET;

            /* First do charges */
            for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	      {
                if (this->ra[i]->iq[zz] == n)
		  {
                    this->ra[i]->q[zz] = nmx[n];
                    qtot           += this->ra[i]->q[zz];
                    n++;
		  }
                else if ((this->ra[i]->iq[zz] < n) && (this->ra[i]->iq[zz] >= 0))
		  {
                    for (zzz = 0; (zzz < i); zzz++)
		      {
                        if (this->ra[zzz]->iq[zz] == this->ra[i]->iq[zz])
			  {
                            this->ra[i]->q[zz] = this->ra[zzz]->q[zz];
                            break;
			  }
		      }
                    if (zzz == i)
		      {
                        gmx_fatal(FARGS, "Can not find a previous atom with iq[%d] = %d", zz, n);
		      }

                    /* Only sum those atoms to qtot, that are not part of
                       the "rest" charge */
                    if (this->ra[i]->iq[zz] != -1)
		      {
                        qtot += this->ra[i]->q[zz];
		      }
		  }
                else if (zz == this->ra[i]->nZeta-1)
		  {
                    nrest++;
		  }
                else
		  {
                    qtot += this->ra[i]->q[zz];
		  }
	      }

            if (this->bFitZeta)
	      {
                /* Then do zeta */
                for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
		  {
                    if (this->ra[i]->iz[zz] == n)
		      {
                        zeta               = nmx[n];
                        this->ra[i]->zeta[zz] = zeta;
                        if (this->delta_z >= 0)
			  {
                            real zmin = this->ra[i]->zeta_ref[zz]-this->delta_z;
                            real zmax = this->ra[i]->zeta_ref[zz]+this->delta_z;
                            if (zeta <= zmin)
			      {
                                this->penalty += sqr(zeta-zmin);
			      }
                            else if (zeta >= zmax)
			      {
                                this->penalty += sqr(zmax-zeta);
			      }
			  }
                        else
			  {
                            if (zeta <= this->zmin)
			      {
                                this->penalty += sqr(this->zmin-zeta);
			      }
                            else if (zeta >= this->zmax)
			      {
                                this->penalty += sqr(this->zmax-zeta);
			      }
                            if ((this->rDecrZeta >= 0) && (zz > 0) &&
                                (this->ra[i]->zeta[zz-1] != 0) &&
                                ((this->ra[i]->zeta[zz-1] - zeta) < this->rDecrZeta))
			      {
                                this->penalty += sqr(this->ra[i]->zeta[zz-1] - zeta - this->rDecrZeta);
			      }
			  }
                        n++;
		      }
                    else if ((this->ra[i]->iz[zz] < n) && (this->ra[i]->iz[zz] >= 0))
		      {
                        for (zzz = 0; (zzz < i); zzz++)
			  {
                            if (this->ra[zzz]->iz[zz] == this->ra[i]->iz[zz])
			      {
                                this->ra[i]->zeta[zz] = this->ra[zzz]->zeta[zz];
                                break;
			      }
			  }
                        if (zzz == i)
			  {
                            gmx_fatal(FARGS, "Can not find a previous atom with iz[%d] = %d", zz, n);
			  }
		      }
                    else if ((this->ra[i]->iz[zz] == -1) && (this->ra[i]->zeta[zz] != 0))
		      {
                        gmx_fatal(FARGS, "ra[%d]->iz[%d] = %d whereas ra[%d]->zeta[%d] = %g", i, zz, this->ra[i]->iz[zz], i, zz, this->ra[i]->zeta[zz]);
		      }
		  }
	      }
	  }
      }
    if (NULL != rnd)
      {
        gmx_rng_destroy(rnd);
      }
    if (n != this->nparam)
      {
        gmx_fatal(FARGS, "Whoopsydaisies! n = %d, should be %d. bSet = %d", n, this->nparam, bSet);
      }

    if (nrest > 0)
      {
        dq = (this->qtot-qtot)/nrest;
        if (debug)
	  {
            fprintf(debug, "this->qtot = %g, qtot = %g, nrest = %d, dq = %g\n",
                    this->qtot, qtot, nrest, dq);
	  }
        for (i = 0; (i < this->natom); i++)
	  {
            if (this->ra[i]->iq[this->ra[i]->nZeta-1] == -1)
	      {
                this->ra[i]->q[this->ra[i]->nZeta-1] = dq;
	      }
	  }
      }
    /* Check for excessive charges */
    for (i = 0; (i < this->natom); i++)
      {
        qi = 0;
        for (zz = 0; (zz < this->ra[i]->nZeta); zz++)
	  {
            qi += this->ra[i]->q[zz];
	  }
        if (qi < this->qmin)
	  {
            this->penalty += sqr(this->qmin-qi);
	  }
        else if (qi > this->qmax)
	  {
            this->penalty += sqr(this->qmax-qi);
	  }
        else if ((qi < -0.02) && (this->ra[i]->atomnumber == 1))
	  {
            this->penalty += qi*qi;
	  }
      }
    this->penalty *= this->pfac;
  }

  void Resp::add_point( double x, double y,
                        double z, double V)
  {
    int i;

    i = this->nesp++;
    srenew(this->esp, this->nesp);
    srenew(this->pot, this->nesp);
    srenew(this->pot_calc, this->nesp);
    this->esp[i][XX]  = x;
    this->esp[i][YY]  = y;
    this->esp[i][ZZ]  = z;
    this->pot[i]      = V;
    this->pot_calc[i] = 0;
  }

  real Resp::my_weight( int iatom)
  {
    if (iatom < this->natom)
      {
        return this->watoms;
      }
    else
      {
        return 1.0;
      }
  }

  void Resp::pot_lsq( gmx_stats_t lsq)
  {
    int    i;
    double w;

    for (i = 0; (i < this->nesp); i++)
      {
        w = my_weight( i);
        if (w > 0)
	  {
            gmx_stats_add_point(lsq,
                                gmx2convert(this->pot[i], eg2cHartree_e),
                                gmx2convert(this->pot_calc[i], eg2cHartree_e), 0, 0);
	  }
      }
  }

  void Resp::calc_rms()
  {
    int    i;
    double pot2, s2, sum2, w, wtot, entropy;
    char   buf[STRLEN];

    pot2 = sum2 = wtot = entropy = 0;
    sprintf(buf, " - weight %g in fit", this->watoms);
    for (i = 0; (i < this->nesp); i++)
      {
        w = my_weight( i);
        if ((NULL != debug) && (i < 2*this->natom))
	  {
            fprintf(debug, "ESP %d QM: %g EEM: %g DIFF: %g%s\n",
                    i, this->pot[i], this->pot_calc[i],
                    this->pot[i]-this->pot_calc[i],
                    (i < this->natom)  ? buf : "");
	  }
        s2    = w*sqr(this->pot[i]-this->pot_calc[i]);
        if ((s2 > 0) && (this->bEntropy))
	  {
            entropy += s2*log(s2);
	  }
        sum2 += s2;
        pot2 += w*sqr(this->pot[i]);
        wtot += w;
      }
    this->wtot = wtot;
    if (wtot > 0)
      {
        this->rms     = gmx2convert(sqrt(sum2/wtot), eg2cHartree_e);
        this->entropy = gmx2convert(entropy/wtot, eg2cHartree_e);
      }
    else
      {
        this->rms     = 0;
        this->entropy = 0;
      }
    this->rrms = sqrt(sum2/pot2);
  }

  double Resp::get_rms( real *wtot)
  {
    calc_rms();
    *wtot = this->wtot;
    if (this->bEntropy)
      {
        return this->entropy;
      }
    else
      {
        return this->rms;
      }
  }

  void Resp::calc_penalty()
  {
    int    i;
    double p, b2;

    p = 0;
    if (this->bAXpRESP && (this->iDistributionModel == eqdAXp))
      {
        b2 = sqr(this->b_hyper);
        for (i = 0; (i < this->natom); i++)
	  {
            p += sqrt(sqr(this->ra[i]->q[0]) + b2) - this->b_hyper;
	  }
        p = (this->qfac * p);
      }
    this->penalty += p;
  }

  //Writen in c stile, needed as an function argument
  double Resp::charge_function(void * gr,double v[])
  {
    Resp * resp = (Resp *)gr;
    double     rms = 0;
    real       wtot;

    resp->get_set_vector( false, false, false, resp->seed, v);
    resp->calc_pot();
    resp->calc_penalty();
    rms = resp->get_rms( &wtot);

    return rms; // + this->penalty;
  }

  void Resp::statistics( int len, char buf[])
  {
    if (len >= 100)
      {
        sprintf(buf, "RMS: %10e [Hartree/e] RRMS: %10e Entropy: %10e Penalty: %10e",
                this->rms, this->rrms, this->entropy, this->penalty);
      }
    else
      {
        fprintf(stderr, "buflen too small (%d) in gmx_resp_statistics\n", len);
      }
  }

  int Resp::optimize_charges(FILE *fp,  int maxiter,
			     real toler, real *rms)
  {
    double *param;
    double  ccc;
    int     bConv;
    char    buf[STRLEN];

    snew(param, this->nparam);

    get_set_vector( true, this->bRandQ, this->bRandZeta, this->seed, param);

    bConv = nmsimplex(fp, (void *)this, charge_function, param, this->nparam,
                      toler, 1, maxiter, &ccc);
    if (bConv)
      {
        statistics( STRLEN-1, buf);
      }
    else
      {
        printf("NM Simplex did not converge\n\n");
      }

    if (this->bEntropy)
      {
        *rms = this->entropy;
      }
    else
      {
        *rms = this->rms;
      }

    get_set_vector( false, false, false, this->seed, param);

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

  Resp::~Resp()
  {
    int i;

    sfree(this->x);
    sfree(this->esp);
    sfree(this->pot);
    i = 0;
    while (NULL != this->dzatoms[i])
      {
        sfree(this->dzatoms[i]);
        i++;
      }
    if (NULL != this->dzatoms)
      {
        sfree(this->dzatoms);
      }
    for (i = 0; (i < this->natom); i++)
      {
        delete this->ra[i];
      }
    sfree(this->ra);
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
        for (i = 0; (i < this->nesp); i++)
	  {
            /* Conversion may or may not be in vain depending on unit */
            exp = gmx2convert(this->pot[i], unit);
            eem = gmx2convert(this->pot_calc[i], unit);
            if (i == this->natom)
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
        for (i = 0; (i < this->nesp); i++)
	  {
            exp = gmx2convert(this->pot[i], eg2cHartree_e);
            eem = gmx2convert(this->pot_calc[i], eg2cHartree_e);
            pp  = this->pot[i]-this->pot_calc[i];
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, "HE", "HE", ' ', i+1, ' ', 20*this->esp[i][XX],
                    20*this->esp[i][YY], 20*this->esp[i][ZZ], 0.0, pp);
	  }
        fclose(fp);
      }
  }

  double Resp::get_qtot( int atom)
  {
    int    i;
    double q = 0;

    range_check(atom, 0, this->natom);
    for (i = 0; (i < this->ra[atom]->nZeta); i++)
      {
        q += this->ra[atom]->q[i];
      }
    return q;
  }

  double Resp::get_q( int atom, int zz)
  {
    range_check(atom, 0, this->natom);
    range_check(zz, 0, this->ra[atom]->nZeta);

    return this->ra[atom]->q[zz];
  }

  double Resp::get_zeta( int atom, int zz)
  {
    range_check(atom, 0, this->natom);
    range_check(zz, 0, this->ra[atom]->nZeta);

    return this->ra[atom]->zeta[zz];
  }

  void Resp::set_q( int atom, int zz, double q)
  {
    range_check(atom, 0, this->natom);
    range_check(zz, 0, this->ra[atom]->nZeta);

    this->ra[atom]->q[zz] = q;
  }

  void Resp::set_zeta( int atom, int zz, double zeta)
  {
    range_check(atom, 0, this->natom);
    range_check(zz, 0, this->ra[atom]->nZeta);

    this->ra[atom]->zeta[zz] = zeta;
  }
}
