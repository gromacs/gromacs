/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROtesk MACabre and Sinister
 */

#ifndef	_recio_h
#define	_recio_h

#ifdef HAVE_IDENT
#ident	"@(#) recio.h 1.30 9/25/97"
#endif /* HAVE_IDENT */

/*
 * This module defines macro's for the binary write and read operation on 
 * some structures:
 *
 * header_blockio(fp,io,sh)
 *   FILE *       fp
 *   io           identifier: write / read
 *   t_statheader sh
 *   Writes or reads an status header to or from the file specified by fp.
 * user_blockio(fp,io,user)
 *   FILE *       fp
 *   io           identifier: write / read
 *   t_user       user
 *   Writes or reads an user record to or from the file specified by fp.
 * emin_blockio(fp,io,emin)
 *   FILE *       fp
 *   io           identifier: write / read
 *   t_emin       emin
 *   Writes or reads an energy minimization record to or from the file 
 *   specified by fp.
 * inputrec_blockio(fp,io,ir)
 *   FILE *       fp
 *   io           identifier: write / read
 *   t_inputrec   ir
 *   Writes or reads an input record to or from the file specified by fp.
 *
 * Example:
 * 
 *    FILE *fp;
 *    t_inputrec ir;
 *    fp=fopen("ir.dat","r");
 *    inputrec_blockio(fp,read,ir);
 *    fclose(fp);
 *
 * Notice that no address operator is needed for the structure.
 */
 
#include "binio.h"
#include "smalloc.h"

#define header_blockio(fp,io,sh) \
  do \
    { \
      block##io(fp,(sh).ir_size); \
      block##io(fp,(sh).e_size); \
      block##io(fp,(sh).box_size); \
      block##io(fp,(sh).vir_size); \
      block##io(fp,(sh).pres_size); \
      block##io(fp,(sh).top_size); \
      block##io(fp,(sh).sym_size); \
      block##io(fp,(sh).x_size); \
      block##io(fp,(sh).v_size); \
      block##io(fp,(sh).f_size); \
      \
      block##io(fp,(sh).natoms); \
      block##io(fp,(sh).step); \
      block##io(fp,(sh).nre); \
      block##io(fp,(sh).t); \
      block##io(fp,(sh).lambda); \
    } \
  while (0)

#define grpopts_blockio(fp,io,opts) \
do { \
  block##io(fp,(opts).ngtc); \
  block##io(fp,(opts).ngacc); \
  block##io(fp,(opts).ngfrz); \
  block##io(fp,(opts).ngener); \
  if (strcmp("read",#io) == 0) { \
    snew((opts).nrdf,(opts).ngtc); \
    snew((opts).ref_t,(opts).ngtc); \
    snew((opts).tau_t,(opts).ngtc); \
    snew((opts).nFreeze,(opts).ngfrz); \
    snew((opts).acc,(opts).ngacc); \
  } \
  if ((opts).ngtc > 0) {\
    nblock##io(fp,(opts).ngtc,(opts).nrdf); \
    nblock##io(fp,(opts).ngtc,(opts).ref_t); \
    nblock##io(fp,(opts).ngtc,(opts).tau_t); \
		      }\
  if ((opts).ngfrz > 0) nblock##io(fp,DIM*(opts).ngfrz,(opts).nFreeze[0]); \
  if ((opts).ngacc > 0) nblock##io(fp,DIM*(opts).ngacc,(opts).acc[0]); \
} while(0)

#define cosine_blockio(fp,io,ef) \
do { \
  block##io(fp,(ef).n); \
  if (strcmp("read",#io) == 0) { \
    snew((ef).a,(ef).n); \
    snew((ef).phi,(ef).n); \
  } \
  nblock##io(fp,(ef).n,(ef).a); \
  nblock##io(fp,(ef).n,(ef).phi); \
} while (0); 
  
#define inputrec_blockio(fp,io,ir) \
  do \
    { \
      int i; \
      block##io(fp,(ir).eI); \
      block##io(fp,(ir).nsteps); \
      block##io(fp,(ir).ns_type); \
      block##io(fp,(ir).nstlist); \
      block##io(fp,(ir).ndelta); \
      block##io(fp,(ir).nstcomm); \
      block##io(fp,(ir).nstprint); \
      block##io(fp,(ir).nstxout); \
      block##io(fp,(ir).nstvout); \
      block##io(fp,(ir).nstfout); \
      block##io(fp,(ir).nstgrp); \
      block##io(fp,(ir).nstxtcout); \
      block##io(fp,(ir).init_t); \
      block##io(fp,(ir).delta_t); \
      block##io(fp,(ir).xtcprec); \
      block##io(fp,(ir).niter); \
      block##io(fp,(ir).watertype); \
      block##io(fp,(ir).nwatoms); \
      block##io(fp,(ir).gausswidth); \
      block##io(fp,(ir).nkx); \
      block##io(fp,(ir).nky); \
      block##io(fp,(ir).nkz); \
      block##io(fp,(ir).eBox); \
      block##io(fp,(ir).bShakeFirst); \
      block##io(fp,(ir).btc); \
      block##io(fp,(ir).ntcmemory); \
      block##io(fp,(ir).bpc); \
      block##io(fp,(ir).npcmemory); \
      block##io(fp,(ir).tau_p); \
      block##io(fp,(ir).ref_p); \
      block##io(fp,(ir).compress); \
      block##io(fp,(ir).bSimAnn); \
      block##io(fp,(ir).zero_temp_time); \
      block##io(fp,(ir).eeltype); \
      block##io(fp,(ir).rshort); \
      block##io(fp,(ir).rlong); \
      block##io(fp,(ir).bLJcorr); \
      block##io(fp,(ir).epsilon_r); \
      block##io(fp,(ir).tol); \
      block##io(fp,(ir).fudgeLJ); \
      block##io(fp,(ir).fudgeQQ); \
      block##io(fp,(ir).bPert); \
      block##io(fp,(ir).init_lambda); \
      block##io(fp,(ir).delta_lambda); \
      block##io(fp,(ir).dr_fc); \
      block##io(fp,(ir).dr_tau); \
      block##io(fp,(ir).dihr_fc); \
      block##io(fp,(ir).em_stepsize); \
      block##io(fp,(ir).em_tol); \
      block##io(fp,(ir).eShakeType); \
      block##io(fp,(ir).nProjOrder); \
      block##io(fp,(ir).LincsWarnAngle); \
      block##io(fp,(ir).nstLincsout); \
      block##io(fp,(ir).ld_temp); \
      block##io(fp,(ir).ld_fric); \
      block##io(fp,(ir).ld_seed); \
      block##io(fp,(ir).userint1); \
      block##io(fp,(ir).userint2); \
      block##io(fp,(ir).userint3); \
      block##io(fp,(ir).userint4); \
      block##io(fp,(ir).userreal1); \
      block##io(fp,(ir).userreal2); \
      block##io(fp,(ir).userreal3); \
      block##io(fp,(ir).userreal4); \
      grpopts_blockio(fp,io,(ir).opts); \
      for(i=0; (i<DIM); i++) \
        cosine_blockio(fp,io,(ir).ex[i]); \
      for(i=0; (i<DIM); i++) \
        cosine_blockio(fp,io,(ir).et[i]); \
    } \
  while (0)

#define iparams_blockio(fp,io,ftype,iparams) \
  do \
    { \
      switch (ftype) \
        { \
        case F_BHAM: \
          block##io(fp,(iparams).bham.a); \
          block##io(fp,(iparams).bham.b); \
          block##io(fp,(iparams).bham.c); \
          break; \
        case F_ANGLES: \
        case F_BONDS: \
        case F_IDIHS: \
          block##io(fp,(iparams).harmonic.rA); \
          block##io(fp,(iparams).harmonic.krA); \
          block##io(fp,(iparams).harmonic.rB); \
          block##io(fp,(iparams).harmonic.krB); \
          break; \
        case F_LJ: \
          block##io(fp,(iparams).lj.c6); \
          block##io(fp,(iparams).lj.c12); \
          break; \
        case F_LJ14: \
          block##io(fp,(iparams).lj14.c6A); \
          block##io(fp,(iparams).lj14.c12A); \
          block##io(fp,(iparams).lj14.c6B); \
          block##io(fp,(iparams).lj14.c12B); \
          break; \
        case F_PDIHS: \
          block##io(fp,(iparams).pdihs.phiA); \
          block##io(fp,(iparams).pdihs.cpA); \
          block##io(fp,(iparams).pdihs.phiB); \
          block##io(fp,(iparams).pdihs.cpB); \
          block##io(fp,(iparams).pdihs.mult); \
          break; \
        case F_DISRES: \
          block##io(fp,(iparams).disres.type); \
          block##io(fp,(iparams).disres.index); \
          block##io(fp,(iparams).disres.rx0); \
          block##io(fp,(iparams).disres.rx1); \
          block##io(fp,(iparams).disres.rx2); \
          block##io(fp,(iparams).disres.rx3); \
          break; \
        case F_POSRES: \
          nblock##io(fp,DIM,(iparams).posres.pos0); \
          nblock##io(fp,DIM,(iparams).posres.fc); \
          break; \
        case F_MORSE: \
	  block##io(fp,(iparams).morse.b0); \
	  block##io(fp,(iparams).morse.cb); \
	  block##io(fp,(iparams).morse.beta); \
          break; \
        case F_RBDIHS: \
          nblock##io(fp,NR_RBDIHS,(iparams).rbdihs.rbc); \
          break; \
        case F_SHAKE: \
          block##io(fp,(iparams).shake.dA); \
          block##io(fp,(iparams).shake.dB); \
          break; \
        case F_SETTLE: \
          block##io(fp,(iparams).settle.doh); \
          block##io(fp,(iparams).settle.dhh); \
          break; \
        case F_DUMMY1: \
        case F_DUMMY2: \
        case F_DUMMY3: \
          block##io(fp,(iparams).dummy.a); \
          block##io(fp,(iparams).dummy.b); \
          block##io(fp,(iparams).dummy.c); \
          block##io(fp,(iparams).dummy.d); \
          block##io(fp,(iparams).dummy.e); \
          block##io(fp,(iparams).dummy.f); \
          break; \
        default: \
          fatal_error(0,"unknown function type %d in %s line %d", \
                      ftype,__FILE__,__LINE__); \
        } \
    } \
  while (0)

#endif	/* _recio_h */
