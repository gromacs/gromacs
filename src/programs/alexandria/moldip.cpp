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
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/shellfc.h"
#include "stringutil.h"

// Alexandria stuff
#include "poldata_xml.h"
#include "molprop_xml.h"
#include "gmx_simple_comm.h"
#include "moldip.h"

#define STRLEN 256

static void add_index_count(t_index_count *ic, const char *name, gmx_bool bConst)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            if (ic->bConst[i] != bConst)
            {
                gmx_fatal(FARGS, "Trying to add atom %s as both constant and optimized",
                          name);
            }
            else
            {
                fprintf(stderr, "Trying to add %s twice\n", name);
            }
        }
    }
    if (i == ic->n)
    {
        ic->n++;
        srenew(ic->name, ic->n);
        srenew(ic->tot_count, ic->n);
        srenew(ic->count, ic->n);
        srenew(ic->bConst, ic->n);
        ic->name[i]      = strdup(name);
        ic->tot_count[i] = 0;
        ic->count[i]     = 0;
        ic->bConst[i]    = bConst;
        if (bConst)
        {
            ic->nconst++;
        }
        else
        {
            ic->nopt++;
        }
    }
}

static void sum_index_count(t_index_count *ic, t_commrec *cr)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        ic->tot_count[i] = ic->count[i];
    }
    if (cr->nnodes > 1)
    {
        gmx_sumi(ic->n, ic->tot_count, cr);
    }
}

static void inc_index_count(t_index_count *ic, char *name)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            ic->count[i]++;
            break;
        }
    }
    if (i == ic->n)
    {
        gmx_fatal(FARGS, "No such atom %s", name);
    }
}

static void dec_index_count(t_index_count *ic, char *name)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            if (ic->count[i] > 0)
            {
                ic->count[i]--;
                break;
            }
            else
            {
                fprintf(stderr, "Trying to decrease number of atoms %s below zero\n",
                        name);
            }
        }
    }
    if (i == ic->n)
    {
        fprintf(stderr, "No such atom %s ic->n = %d\n", name, ic->n);
    }
}

static int n_index_count(t_index_count *ic, char *name)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            return i;
        }
    }
    return -1;
}

static int c_index_count(t_index_count *ic, char *name)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            return ic->tot_count[i];
        }
    }
    return 0;
}

char *opt_index_count(t_index_count *ic)
{
    for (; (ic->nopt_c < ic->n); ic->nopt_c++)
    {
        if (!ic->bConst[ic->nopt_c])
        {
            return ic->name[ic->nopt_c++];
        }
    }
    ic->nopt_c = 0;

    return NULL;
}

static gmx_bool const_index_count(t_index_count *ic, char *name)
{
    int i;

    for (i = 0; (i < ic->n); i++)
    {
        if (strcasecmp(ic->name[i], name) == 0)
        {
            return ic->bConst[i];
        }
    }
    return FALSE;
}

static void dump_index_count(t_index_count *ic, FILE *fp,
                             ChargeGenerationModel iModel, gmx_poldata_t pd,
                             gmx_bool bFitZeta)
{
    int    i, j, nZeta, nZopt;
    double zz;
    if (fp)
    {
        fprintf(fp, "Atom index for this optimization.\n");
        fprintf(fp, "Name  Number  Action   #Zeta\n");
        for (i = 0; (i < ic->n); i++)
        {
            nZeta = gmx_poldata_get_nzeta(pd, iModel, ic->name[i]);
            nZopt = 0;
            for (j = 0; (j < nZeta); j++)
            {
                zz = gmx_poldata_get_zeta(pd, iModel, ic->name[i], j);
                if (zz > 0)
                {
                    nZopt++;
                }
            }
            if (ic->bConst[i])
            {
                fprintf(fp, "%-4s  %6d  Constant\n", ic->name[i], ic->count[i]);
            }
            else
            {
                fprintf(fp, "%-4s  %6d  Optimized %4d%s\n",
                        ic->name[i], ic->count[i], nZopt,
                        bFitZeta ? " optimized" : " constant");
            }
        }
        fprintf(fp, "\n");
        fflush(fp);
    }
}

static int clean_index_count(t_index_count *ic, int minimum_data, FILE *fp)
{
    int i, j, nremove = 0;

    for (i = 0; (i < ic->n); )
    {
        if (!ic->bConst[i] && (ic->tot_count[i] < minimum_data))
        {
            if (fp)
            {
                fprintf(fp, "Not enough support in data set for optimizing %s\n",
                        ic->name[i]);
            }
            sfree(ic->name[i]);
            for (j = i; (j < ic->n-1); j++)
            {
                ic->name[j]       = ic->name[j+1];
                ic->count[j]      = ic->count[j+1];
                ic->tot_count[j]  = ic->tot_count[j+1];
                ic->bConst[j]     = ic->bConst[j+1];
            }
            nremove++;
            ic->n--;
            ic->nopt--;
        }
        else
        {
            i++;
        }
    }
    return nremove;
}

static void update_index_count_bool(t_index_count *ic, gmx_poldata_t pd,
                                    const char *string, gmx_bool bSet,
                                    gmx_bool bAllowZero, ChargeGenerationModel iModel)
{
    std::vector<std::string> ptr = split(string, ' ');
    for (std::vector<std::string>::iterator k = ptr.begin();
         (k < ptr.end()); ++k)
    {
        if (gmx_poldata_have_eem_support(pd, iModel, k->c_str(), bAllowZero))
        {
            add_index_count(ic, k->c_str(), bSet);
        }
    }
}

static int check_data_sufficiency(FILE *fp,
                                  std::vector<alexandria::MyMol> mol,
                                  int minimum_data, gmx_poldata_t pd,
                                  t_index_count *ic, ChargeGenerationModel iModel, char *opt_elem, char *const_elem,
                                  t_commrec *cr, gmx_bool bPol,
                                  gmx_bool bFitZeta)
{
    std::vector<alexandria::MyMol>::iterator mmi;

    int                                      j, nremove, nsupported;
    gmx_mtop_atomloop_all_t                  aloop;
    t_atom                                  *atom;
    char                                    *myname;
    int                                      k, at_global;
    ChargeGenerationModel                    mymodel;

    /* Parse opt_elem list to test which elements to optimize */
    if (NULL != const_elem)
    {
        update_index_count_bool(ic, pd, const_elem, TRUE, FALSE, iModel);
    }
    if (NULL != opt_elem)
    {
        update_index_count_bool(ic, pd, opt_elem, FALSE, TRUE, iModel);
    }
    else
    {
        while (gmx_poldata_get_eemprops(pd, &mymodel, &myname, NULL, NULL, NULL, NULL, NULL) != 0)
        {
            if ((mymodel == iModel) &&
                !const_index_count(ic, myname) &&
                gmx_poldata_have_eem_support(pd, iModel, myname, FALSE))
            {
                add_index_count(ic, myname, FALSE);
            }
        }
    }
    sum_index_count(ic, cr);
    dump_index_count(ic, debug, iModel, pd, bFitZeta);
    for (mmi = mol.begin(); (mmi < mol.end()); mmi++)
    {
        if (mmi->eSupp != eSupportNo)
        {
            aloop = gmx_mtop_atomloop_all_init(mmi->mtop_);
            k     = 0;
            while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom) &&
                   (mmi->eSupp != eSupportNo))
            {
                if ((atom->atomnumber > 0) || !bPol)
                {
                    if (n_index_count(ic, *(mmi->topology_->atoms.atomtype[k])) == -1)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Removing %s because of lacking support for atom %s\n",
                                    mmi->GetMolname().c_str(),
                                    *(mmi->topology_->atoms.atomtype[k]));
                        }
                        mmi->eSupp = eSupportNo;
                    }
                }
                k++;
            }
            if (mmi->eSupp != eSupportNo)
            {
                gmx_assert(k, mmi->topology_->atoms.nr);
                aloop = gmx_mtop_atomloop_all_init(mmi->mtop_);
                k     = 0;
                while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
                {
                    if ((atom->atomnumber > 0) || !bPol)
                    {
                        inc_index_count(ic, *(mmi->topology_->atoms.atomtype[k]));
                    }
                    k++;
                }
                gmx_assert(k, mmi->topology_->atoms.nr);
            }
        }
    }
    do
    {
        sum_index_count(ic, cr);
        dump_index_count(ic, debug, iModel, pd, bFitZeta);
        nremove = 0;
        for (mmi = mol.begin(); (mmi < mol.end()); mmi++)
        {
            if (mmi->eSupp != eSupportNo)
            {
                j     = 0;
                aloop = gmx_mtop_atomloop_all_init(mmi->mtop_);
                while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
                {
                    if (c_index_count(ic, *(mmi->topology_->atoms.atomtype[j])) < minimum_data)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Removing %s because of no support for name %s\n",
                                    mmi->GetMolname().c_str(),
                                    *(mmi->topology_->atoms.atomtype[j]));
                        }
                        break;
                    }
                    j++;
                }
                if (j < mmi->mtop_->natoms)
                {
                    aloop = gmx_mtop_atomloop_all_init(mmi->mtop_);
                    k     = 0;
                    while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
                    {
                        dec_index_count(ic, *(mmi->topology_->atoms.atomtype[k++]));
                    }
                    mmi->eSupp = eSupportNo;
                    nremove++;
                }
            }
        }
        if (cr->nnodes > 1)
        {
            /* Sum nremove */
            gmx_sumi(1, &nremove, cr);
        }
    }
    while (nremove > 0);
    nremove = clean_index_count(ic, minimum_data, debug);
    sum_index_count(ic, cr);
    dump_index_count(ic, fp, iModel, pd, bFitZeta);

    nsupported = 0;
    for (mmi = mol.begin(); (mmi < mol.end()); mmi++)
    {
        if (mmi->eSupp == eSupportLocal)
        {
            if (NULL != debug)
            {
                fprintf(debug, "Supported molecule %s on CPU %d\n",
                        mmi->GetMolname().c_str(), cr->nodeid);
            }
            nsupported++;
        }
    }
    if (cr->nnodes > 1)
    {
        gmx_sumi(1, &nsupported, cr);
    }
    if (fp)
    {
        fprintf(fp, "Removed %d atomtypes\n", nremove);
        fprintf(fp, "There are %d supported molecules left.\n\n", nsupported);
    }
    return nsupported;
}

namespace alexandria
{

MolDip::MolDip()
{
    _ic     = NULL;
    _cr     = NULL;
    _fixchi = NULL;
    for (int i = 0; (i < ermsNR); i++)
    {
        _fc[i]   = 0;
        _ener[i] = 0;
    }
}

void MolDip::Init(t_commrec *cr, gmx_bool bQM, gmx_bool bGaussianBug,
                  ChargeGenerationModel iModel, real rDecrZeta, real epsr,
                  real J0_0, real Chi0_0, real w_0,
                  real J0_1, real Chi0_1, real w_1,
                  real fc_bound, real fc_mu, real fc_quad, real fc_charge,
                  real fc_esp, real fc_epot, real fc_force, char *fixchi,
                  gmx_bool bOptHfac, real hfac,
                  gmx_bool bPol, gmx_bool bFitZeta)
{
    _cr             = cr;
    _bQM            = bQM;
    _bDone          = FALSE;
    _bFinal         = FALSE;
    _bGaussianBug   = bGaussianBug;
    _bFitZeta       = bFitZeta;
    _iModel         = iModel;
    _decrzeta       = rDecrZeta;
    _epsr           = epsr;
    _J0_0           = J0_0;
    _Chi0_0         = Chi0_0;
    _w_0            = w_0;
    _J0_1           = J0_1;
    _Chi0_1         = Chi0_1;
    _w_1            = w_1;
    _fc[ermsMU]     = fc_mu;
    _fc[ermsBOUNDS] = fc_bound;
    _fc[ermsQUAD]   = fc_quad;
    _fc[ermsCHARGE] = fc_charge;
    _fc[ermsESP]    = fc_esp;
    _fc[ermsForce2] = fc_force;
    _fc[ermsEPOT]   = fc_epot;
    _fixchi         = strdup(fixchi);
    _hfac           = hfac;
    _hfac0          = hfac;
    _bOptHfac       = bOptHfac;
    _bPol           = bPol;
}

void MolDip::Read(FILE *fp, const char *fn, const char *pd_fn,
                  int minimum_data,
                  gmx_bool bZero,
                  char *opt_elem, char *const_elem,
                  char *lot,
                  output_env_t oenv, gmx_molselect_t gms,
                  real watoms, gmx_bool bCheckSupport,
                  unsigned int seed)
{
    int                              i, n, nwarn = 0, nmol_cpu;
    int                              nexcl, imm_count[immNR];
    immStatus                        imm;
    std::vector<alexandria::MolProp> mp;
    alexandria::GaussAtomProp        gap;

    for (i = 0; (i < immNR); i++)
    {
        imm_count[i] = 0;
    }

    /* Read the EEM parameters */
    _atomprop   = gmx_atomprop_init();

    /* Force field data */
    if ((_pd = gmx_poldata_read(pd_fn, _atomprop)) == NULL)
    {
        gmx_fatal(FARGS, "Can not read the force field information. File %s missing or incorrect.", pd_fn);
    }

    if ((n = gmx_poldata_get_numprops(_pd, _iModel)) == 0)
    {
        gmx_fatal(FARGS, "File %s does not contain the requested parameters for model %d", pd_fn, _iModel);
    }

    nexcl = gmx_poldata_get_nexcl(_pd);

    if (NULL != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s:\n---\n",
                n, pd_fn);
        fprintf(fp, "---\n\n");
    }

    /* Read other stuff */
    if (MASTER(_cr))
    {
        /* Now read the molecules */
        MolPropRead(fn, mp);
        for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            if (false == mpi->GenerateComposition(_pd))
            {
                mpi = mp.erase(mpi);
            }
            mpi->CheckConsistency();
        }
        nmol_cpu = mp.size()/_cr->nnodes + 1;
    }
    else
    {
        nmol_cpu = 0;
    }
    if (PAR(_cr))
    {
        gmx_sumi(1, &nmol_cpu, _cr);
    }

    if (MASTER(_cr))
    {
        i = -1;
        for (alexandria::MolPropIterator mpi = mp.begin(); (mpi < mp.end()); mpi++)
        {
            if (imsTrain == gmx_molselect_status(gms, mpi->GetIupac().c_str()))
            {
                int               dest = (n % _cr->nnodes);
                alexandria::MyMol mpnew;

                mpnew.Merge(*mpi);

                imm = mpnew.GenerateTopology(_atomprop, _pd, lot, _iModel, _bPol, nexcl);

                if (_iModel != eqgNone)
                {
                    if (immOK == imm)
                    {
                        mpnew.gr_ = gmx_resp_init(_iModel, TRUE, 0.001, 0.1, mpnew.GetCharge(),
                                                  1, 100, 5,
                                                  TRUE, watoms, 5, TRUE, TRUE,
                                                  1, TRUE,
                                                  TRUE, NULL, seed);
                        if (NULL == mpnew.gr_)
                        {
                            imm = immRespInit;
                        }
                    }

                    if (immOK == imm)
                    {
                        imm = mpnew.GenerateCharges(_pd, _atomprop, _iModel, _hfac, _epsr,
                                                    lot, TRUE, NULL);
                    }
                }
                if (immOK == imm)
                {
                    imm = mpnew.GenerateChargeGroups(ecgAtom, FALSE, NULL, 1);
                }
                if (immOK == imm)
                {
                    imm = mpnew.GenerateGromacs(oenv, _cr);
                }
                if (immOK == imm)
                {
                    imm = mpnew.GetExpProps(_bQM, bZero, lot, gap);
                }

                if (immOK == imm)
                {
                    if (dest > 0)
                    {
                        mpnew.eSupp = eSupportRemote;
                        /* Send another molecule */
                        if (NULL != debug)
                        {
                            fprintf(debug, "Going to send %s to cpu %d\n",
                                    mpi->GetMolname().c_str(), dest);
                        }
                        gmx_send_int(_cr, dest, 1);
                        CommunicationStatus cs = mpi->Send(_cr, dest);
                        if (CS_OK != cs)
                        {
                            imm = immCommProblem;
                        }
                        else
                        {
                            imm = (immStatus)gmx_recv_int(_cr, dest);
                        }
                        if (imm != immOK)
                        {
                            fprintf(stderr, "Molecule %s was not accepted on node %d - error %s\n",
                                    mpnew.GetMolname().c_str(), dest, alexandria::immsg(imm));
                        }
                        else if (NULL != debug)
                        {
                            fprintf(debug, "Succesfully beamed over %s\n", mpi->GetMolname().c_str());
                        }

                    }
                    else
                    {
                        mpnew.eSupp = eSupportLocal;
                    }
                    if (immOK == imm)
                    {
                        _mymol.push_back(mpnew);
                        n++;
                        if (NULL != debug)
                        {
                            fprintf(debug, "Added %s, n = %d\n",
                                    mpnew.GetMolname().c_str(), n);
                        }
                    }
                }
                if ((immOK != imm) && (NULL != debug))
                {
                    fprintf(debug, "IMM: Dest: %d %s - %s\n",
                            dest, mpi->GetMolname().c_str(), immsg(imm));
                }
            }
            else
            {
                imm = immTest;
            }
            imm_count[imm]++;
        }
        /* Send signal done with transferring molecules */
        for (i = 1; (i < _cr->nnodes); i++)
        {
            gmx_send_int(_cr, i, 0);
        }
    }
    else
    {
        /***********************************************
         *
         *           S L A V E   N O D E S
         *
         ***********************************************/
        n = 0;
        while (gmx_recv_int(_cr, 0) == 1)
        {
            /* Receive another molecule */
            alexandria::MyMol mpnew;

            if (NULL != debug)
            {
                fprintf(debug, "Going to retrieve new molecule\n");
            }
            CommunicationStatus cs = mpnew.Receive(_cr, 0);
            if (CS_OK != cs)
            {
                imm = immCommProblem;
            }
            else if (NULL != debug)
            {
                fprintf(debug, "Succesfully retrieved %s\n", mpnew.GetMolname().c_str());
                fflush(debug);
            }

            imm = mpnew.GenerateTopology(_atomprop, _pd, lot, _iModel, _bPol, nexcl);

            if (immOK == imm)
            {
                mpnew.gr_ = gmx_resp_init(_iModel, TRUE, 0.001, 0.1, mpnew.GetCharge(),
                                          1, 100, 5,
                                          TRUE, watoms, 5, TRUE, TRUE,
                                          1, TRUE,
                                          TRUE, NULL, seed);
                if (NULL == mpnew.gr_)
                {
                    imm = immRespInit;
                }
            }

            if (immOK == imm)
            {
                imm = mpnew.GenerateCharges(_pd, _atomprop, _iModel, _hfac, _epsr,
                                            lot, TRUE, NULL);
            }
            if (immOK == imm)
            {
                imm = mpnew.GenerateChargeGroups(ecgAtom, FALSE, NULL, 1);
            }
            if (immOK == imm)
            {
                imm = mpnew.GenerateGromacs(oenv, _cr);
            }
            if (immOK == imm)
            {
                imm = mpnew.GetExpProps(_bQM, bZero, lot, gap);
            }

            mpnew.eSupp = eSupportLocal;
            imm_count[imm]++;
            if (immOK == imm)
            {
                _mymol.push_back(mpnew);
                if (NULL != debug)
                {
                    fprintf(debug, "Added molecule %s\n", mpnew.GetMolname().c_str());
                }
            }
            gmx_send_int(_cr, 0, imm);
        }
    }
    int  nnn = nmol_cpu;
    int *nmolpar;

    if (PAR(_cr))
    {
        snew(nmolpar, _cr->nnodes);
        nmolpar[_cr->nodeid] = nnn;
        gmx_sumi(_cr->nnodes, nmolpar, _cr);
    }

    if (fp)
    {
        fprintf(fp, "There were %d warnings because of zero error bars.\n", nwarn);
        int nmoltot = 0;
        if (PAR(_cr))
        {
            for (i = 0; (i < _cr->nnodes); i++)
            {
                fprintf(fp, "node %d has %d molecules\n", i, nmolpar[i]);
                nmoltot += nmolpar[i];
            }
        }
        else
        {
            nmoltot = mp.size();
        }
        fprintf(fp, "Made topologies for %d out of %d molecules.\n", n,
                (MASTER(_cr)) ? nmoltot : nmol_cpu);

        for (i = 0; (i < immNR); i++)
        {
            if (imm_count[i] > 0)
            {
                fprintf(fp, "%d molecules - %s.\n", imm_count[i], alexandria::immsg((immStatus)i));
            }
        }
        if (imm_count[immOK] != (int)mp.size())
        {
            fprintf(fp, "Check %s.debug for more information.\nYou may have to use the -debug 1 flag.\n\n", ShortProgram());
        }
    }
    sfree(nmolpar);
    snew(_ic, 1);
    if (bCheckSupport)
    {
        _nmol_support =
            check_data_sufficiency(MASTER(_cr) ? fp : NULL, _mymol,
                                   minimum_data, _pd, _ic,
                                   _iModel, opt_elem, const_elem, _cr,
                                   _bPol, _bFitZeta);
        if (_nmol_support == 0)
        {
            gmx_fatal(FARGS, "No support for any molecule!");
        }
    }
    else
    {
        _nmol_support = _mymol.size();
    }
}

static void split_shell_charges(gmx_mtop_t *mtop, t_idef *idef)
{
    int                     k, ai, aj;
    real                    q, Z;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom, *atom_i, *atom_j;
    int                     at_global;
    gmx_mtop_atomlookup_t   alook;

    alook = gmx_mtop_atomlookup_init(mtop);

    for (k = 0; (k < idef->il[F_POLARIZATION].nr); )
    {
        k++; // Skip over the type.
        ai = idef->il[F_POLARIZATION].iatoms[k++];
        aj = idef->il[F_POLARIZATION].iatoms[k++];

        gmx_mtop_atomnr_to_atom(alook, ai, &atom_i);
        gmx_mtop_atomnr_to_atom(alook, aj, &atom_j);

        if ((atom_i->ptype == eptAtom) &&
            (atom_j->ptype == eptShell))
        {
            q         = atom_i->q;
            Z         = atom_i->atomnumber;
            atom_i->q = Z;
            atom_j->q = q-Z;
        }
        else if ((atom_i->ptype == eptAtom) &&
                 (atom_j->ptype == eptShell))
        {
            q         = atom_j->q;
            Z         = atom_j->atomnumber;
            atom_j->q = Z;
            atom_i->q = q-Z;
        }
        else
        {
            gmx_incons("Polarization entry does not have one atom and one shell");
        }
    }
    q     = 0;
    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
    {
        q += atom->q;
    }
    Z = gmx_nint(q);
    if (fabs(q-Z) > 1e-3)
    {
        gmx_fatal(FARGS, "Total charge in molecule is not zero, but %f", q-Z);
    }
    gmx_mtop_atomlookup_destroy(alook);
}

void MolDip::CalcDeviation()
{
    int                     j, atomnr;
    double                  qq, qtot;
    real                    etot[ermsNR];
    real                    t      = 0;
    rvec                    mu_tot = {0, 0, 0};
    gmx_enerdata_t         *epot;
    tensor                  force_vir = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    t_nrnb                  my_nrnb;
    gmx_wallcycle_t         wcycle;
    gmx_bool                bConverged;
    int                     eQ;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                 *atom;
    int                     at_global;

    if (PAR(_cr))
    {
        gmx_bcast(sizeof(_bDone), &_bDone, _cr);
        gmx_bcast(sizeof(_bFinal), &_bFinal, _cr);
    }
    if (_bDone)
    {
        return;
    }
    if (PAR(_cr))
    {
        gmx_poldata_comm_eemprops(_pd, _cr);
    }
    init_nrnb(&my_nrnb);
    snew(epot, 1);

    wcycle  = wallcycle_init(stdout, 0, _cr, 1, 0);
    for (j = 0; (j < ermsNR); j++)
    {
        etot[j]  = 0;
        _ener[j] = 0;
    }
    for (std::vector<MyMol>::iterator mymol = _mymol.begin(); (mymol < _mymol.end()); mymol++)
    {
        if ((mymol->eSupp == eSupportLocal) ||
            (_bFinal && (mymol->eSupp == eSupportRemote)))
        {
            /* Reset energies */
            for (j = 0; (j < ermsNR); j++)
            {
                _ener[j] = 0;
            }

            if (NULL == mymol->qgen_)
            {
                mymol->qgen_ =
                    gentop_qgen_init(_pd, &(mymol->topology_->atoms), _atomprop,
                                     mymol->x_, _iModel, _hfac,
                                     mymol->GetCharge(), _epsr);
            }
            /*if (strcmp(mymol->molname,"1-butene") == 0)
               fprintf(stderr,"Ready for %s\n",mymol->molname);*/
            eQ = generate_charges_sm(debug, mymol->qgen_,
                                     _pd, &(mymol->topology_->atoms),
                                     1e-4, 100, _atomprop,
                                     &(mymol->chieq));
            if (eQ != eQGEN_OK)
            {
                char buf[STRLEN];
                qgen_message(mymol->qgen_, STRLEN, buf, NULL);
                fprintf(stderr, "%s\n", buf);
            }
            else
            {
                aloop = gmx_mtop_atomloop_all_init(mymol->mtop_);
                j     = 0;
                while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
                {
                    atom->q = mymol->topology_->atoms.atom[j].q;
                    j++;
                }
                gmx_assert(j, mymol->topology_->atoms.nr);
            }

            /* Now optimize the shell positions */
            if (mymol->shell_)
            {
                split_shell_charges(mymol->mtop_, &mymol->ltop_->idef);
                atoms2md(mymol->mtop_, mymol->inputrec_, 0, NULL, 0,
                         mymol->md_);
                (void)
                relax_shell_flexcon(debug, _cr, FALSE, 0,
                                    mymol->inputrec_, TRUE,
                                    GMX_FORCE_ALLFORCES,
                                    mymol->ltop_, NULL, NULL, NULL,
                                    &(mymol->state_),
                                    mymol->f_, force_vir, mymol->md_,
                                    &my_nrnb, wcycle, NULL,
                                    &(mymol->mtop_->groups),
                                    mymol->shell_, mymol->fr_, FALSE, t, mu_tot,
                                    &bConverged, NULL, NULL);
            }
            /* Compute the molecular dipole */
            mymol->CalcMultipoles();

            /* Compute the ESP on the points */
            if ((NULL != mymol->gr_) && _bQM)
            {
                /*gmx_resp_add_atom_info(mymol->gr,&(mymol->atoms),_pd);*/
                gmx_resp_fill_zeta(mymol->gr_, _pd);
                gmx_resp_fill_q(mymol->gr_, &(mymol->topology_->atoms));
                gmx_resp_calc_pot(mymol->gr_);
            }
            qtot = 0;
            for (j = 0; (j < mymol->topology_->atoms.nr); j++)
            {
                atomnr = mymol->topology_->atoms.atom[j].atomnumber;
                qq     = mymol->topology_->atoms.atom[j].q;
                qtot  += qq;
                if (((qq < 0) && (atomnr == 1)) ||
                    ((qq > 0) && ((atomnr == 8)  || (atomnr == 9) ||
                                  (atomnr == 16) || (atomnr == 17) ||
                                  (atomnr == 35) || (atomnr == 53))))
                {
                    _ener[ermsBOUNDS] += fabs(qq);
                }
                if (_bQM)
                {
                    _ener[ermsCHARGE] += sqr(qq-mymol->qESP[j]);
                }
            }
            if (0 && (fabs(qtot-mymol->GetCharge()) > 1e-2))
            {
                fprintf(stderr, "Warning qtot for %s is %g, should be %d\n",
                        mymol->GetMolname().c_str(),
                        qtot, mymol->GetCharge());
            }
            if (_bQM)
            {
                int  mm, nn;
                rvec dmu;
                real wtot;

                rvec_sub(mymol->mu_calc, mymol->mu_exp, dmu);
                _ener[ermsMU]   = iprod(dmu, dmu);
                for (mm = 0; (mm < DIM); mm++)
                {
                    if (0)
                    {
                        for (nn = 0; (nn < DIM); nn++)
                        {
                            _ener[ermsQUAD] += sqr(mymol->Q_exp[mm][nn] - mymol->Q_calc[mm][nn]);
                        }
                    }
                    else
                    {
                        _ener[ermsQUAD] += sqr(mymol->Q_exp[mm][mm] - mymol->Q_calc[mm][mm]);
                    }
                }
                if (NULL != mymol->gr_)
                {
                    _ener[ermsESP] += gmx_resp_get_rms(mymol->gr_, &wtot);
                    if (NULL != debug)
                    {
                        fprintf(debug, "RMS %s = %g\n",
                                mymol->GetMolname().c_str(), _ener[ermsESP]);
                    }
                }
            }
            else
            {
                _ener[ermsMU]     = sqr(mymol->dip_calc - mymol->dip_exp);
            }
            for (j = 0; (j < ermsNR); j++)
            {
                etot[j] += _ener[j];
            }
        }
    }
    for (j = 0; (j < ermsTOT); j++)
    {
        _ener[j]       += _fc[j]*etot[j]/_nmol_support;
        _ener[ermsTOT] += _ener[j];
    }
    sfree(epot);
    if (debug)
    {
        fprintf(debug, "ENER:");
        for (j = 0; (j < ermsNR); j++)
        {
            fprintf(debug, "  %8.3f", etot[j]);
        }
        fprintf(debug, "\n");
    }
    /* Global sum energies */
    if (PAR(_cr))
    {
#ifdef GMX_DOUBLE
        gmx_sumd(ermsNR, _ener, _cr);
#else
        gmx_sumf(ermsNR, _ener, _cr);
#endif
    }
}


}
