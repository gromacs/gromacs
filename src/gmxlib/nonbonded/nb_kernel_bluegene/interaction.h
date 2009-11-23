/*
 * Copyright (c) 2006, International Business Machines (IBM) Inc.
 *
 * Author: Mathias Puetz
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 4. Neither the name of IBM nor the names of its contributors may be used
 * to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY IBM AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL IBM OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This include file defines some macros which help to write SIMD code on BlueGene
 */

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

/* Some defines to hold constants describing table dimensions. */

#if COULOMB == COULOMB_TAB

 #define LJTAB_OFS 4

 #if VDW == VDW_TAB
  #define TAB_MULT 12
 #else
  #define TAB_MULT 4
 #endif

#else

 #define LJTAB_OFS 0
 #define TAB_MULT 8

#endif


/*****************************************************************************
 *                        scalar interaction calculations
 *
 * That is, those not using SIMD dual forms. Such code can be
 * identified by macro names that end in an underscore, and variable
 * names that begin with a single underscore.
 *
 *****************************************************************************/

/* This code takes care of calculating eps in non-dual variables.
*/

#if (! defined MAKE_EPS_FAST_SCALAR)

#define make_eps_(r,r0,eps)  \
    r0  = (double)((int) r); \
    eps = r - r0

#else

/* The following fragment doesn't work at -O3 on BlueGene/L, because
 * the compiler optimizes away the addition and subtraction of _round,
 * even though it does not optimize away the same operations in the
 * dual case! Depending how these compilers evolve, one or other form
 * may be useful in the future.
 */
#define make_eps_(r,r0,eps) \
    r0  = r - _half;        \
    r0  = r0 + _round;      \
    r0  = r0 - _round;      \
    eps = r - r0

#endif


/* Macros that will evaluate Coulomb energies and forces for non-dual
 * variables. The code has to change slightly in cases where there are
 * no VDW interactions to compute.
 */

#if COULOMB == COULOMB_TAB

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr)          \
    _rt     = (_rsq * _rinv) *_tabscale;                 \
    make_eps_(_rt,_n0,_eps);                             \
    nnn1    = TAB_MULT*((int)_rt);                       \
    _Y      = VFtab[nnn1+0];                             \
    _F      = VFtab[nnn1+1];                             \
    _G      = VFtab[nnn1+2];                             \
    _H      = VFtab[nnn1+3];                             \
    _GHeps  = _G + _H * _eps;                            \
    _VV     = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF     = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    _fijC   = _qq * _FF;                                 \
    vctot  += _qq * _VV

  #define calc_coulomb_only_force_(_qq,_rinv,_rsq) \
    _fscal  = _fijC * (_rinv * _tabscale)

#if VDW == VDW_NONE
  #define calc_coulomb_force_(_qq,_rinv,_rsq) calc_coulomb_only_force_(_qq,_rinv,_rsq)
#else
  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal += _fijC * (_rinv * _tabscale)
#endif

#elif COULOMB == GENERALIZED_BORN

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr)              \
    _r          = _rinv * _rsq;                              \
    _isaprod    = _isai * invsqrta[jnr];                     \
     _iqq       = _qq * _isaprod;                            \
    _gbscale    = _isaprod * _gbtabscale;                    \
    _rt         = _r * _gbscale;                             \
    make_eps_(_rt,_n0,_eps);                                 \
    nnn1        = 4*((int)_rt);                              \
    _Y          = GBtab[nnn1+0];                             \
    _F          = GBtab[nnn1+1];                             \
    _G          = GBtab[nnn1+2];                             \
    _H          = GBtab[nnn1+3];                             \
    _GHeps      = _G + _H * _eps;                            \
    _VV         = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF         = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    _fijC       = _iqq * _FF;                                \
    _vcoul      = _iqq * _VV;                                \
    vctot      += _vcoul;                                    \
    _dvdatmp    = _vcoul + _fijC * _rt;                      \
    dvdasum    -= _dvdatmp;                                  \
    dvda[jnr]  -= _dvdatmp;

  #define calc_coulomb_only_force_(_qq,_rinv,_rsq) \
    _fscal      = _fijC * (_rinv * _gbscale)

#if VDW == VDW_NONE
  #define calc_coulomb_force_(_qq,_rinv,_rsq) calc_coulomb_only_force_(_qq,_rinv,_rsq)
#else
  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal     += _fijC * (_rinv * _gbscale)
#endif

#elif COULOMB == REACTION_FIELD

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) \
    _krsq   = _krf * _rsq;                      \
    vctot  += _qq * (_rinv + _krsq - _crf)

  #define calc_coulomb_only_force_(_qq,_rinv,_rsq) \
    _fscal  = -_qq * (_rinv - 2.0 * _krsq) * _rinvsq

#if VDW == VDW_NONE
  #define calc_coulomb_force_(_qq,_rinv,_rsq) calc_coulomb_only_force_(_qq,_rinv,_rsq)
#else
  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal -= _qq * (_rinv - 2.0 * _krsq) * _rinvsq
#endif

#elif COULOMB == COULOMB_CUTOFF

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) \
    vctot += _qq * _rinv 

  #define calc_coulomb_only_force_(_qq,_rinv,_rsq) \
    _fscal  = -_qq * _rinv * _rinvsq

#if VDW == VDW_NONE
  #define calc_coulomb_force_(_qq,_rinv,_rsq) calc_coulomb_only_force_(_qq,_rinv,_rsq)
#else
  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal -= _qq * _rinv * _rinvsq
#endif

#else

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) vctot = vctot

  #define calc_coulomb_only_force_(_qq,_rinv,_rsq) \
    _fscal = 0.0

#if VDW == VDW_NONE
  #define calc_coulomb_force_(_qq,_rinv,_rsq) calc_coulomb_only_force_(_qq,_rinv,_rsq)
#else
  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    {}
#endif

#endif


/* Macros that will evaluate VDW energies and forces for non-dual
 * variables. In cases where the VDW parameters are constant
 * (e.g. water loops), different forms of these codes are faster.
 */

#if VDW == LENNARD_JONES

 #if CONST_LJ

  #define calc_vdw_pot_(_rinv,_rsq,jnr)       \
    _rinvsix = _rinvsq * _rinvsq * _rinvsq;   \
    _Vvdw6   = _rinvsix * _c6;                \
    _Vvdw12  = _rinvsix * _rinvsix *  _c12;   \
    Vvdwtot += _Vvdw12 - _Vvdw6

 #else

  #define calc_vdw_pot_(_rinv,_rsq,jnr)       \
    tj       = nti+2*type[jnr];               \
    _c6      = vdwparam[tj];                  \
    _c12     = vdwparam[tj+1];                \
    _rinvsix = _rinvsq * _rinvsq * _rinvsq;   \
    _Vvdw6   = _rinvsix * _c6;                \
    _Vvdw12  = _rinvsix * _rinvsix *  _c12;   \
    Vvdwtot += _Vvdw12 - _Vvdw6

 #endif

  #define calc_vdw_only_force_(_rinv,_rsq)          \
    _fscal = (6.0 * _Vvdw6 - 12.0 * _Vvdw12) * _rinvsq

  #define calc_vdw_force_(_rinv,_rsq) calc_vdw_only_force_(_rinv,_rsq)

#elif VDW == BUCKINGHAM

 #ifdef CONST_LJ

  #define calc_vdw_pot_(_rinv,_rsq,jnr)             \
    _Vvdw6   = _rinvsq * _rinvsq * _rinvsq  * _c6;  \
    _r       = _cexp2 * _rsq * _rinv;               \
    _Vvdwexp = _cexp1 * exp(- _r);                  \
    Vvdwtot += _Vvdwexp - _Vvdw6

 #else

  #define calc_vdw_pot_(_rinv,_rsq,jnr)             \
    tj       = nti+3*type[jnr];                     \
    _c6      = vdwparam[tj];                        \
    _cexp1   = vdwparam[tj+1];                      \
    _cexp2   = vdwparam[tj+2];                      \
    _Vvdw6   = _rinvsq * _rinvsq * _rinvsq  * _c6;  \
    _r       = _cexp2 * _rsq * _rinv;               \
    _Vvdwexp = _cexp1 * exp(- _r);                  \
    Vvdwtot += _Vvdwexp - _Vvdw6

 #endif

  #define calc_vdw_only_force_(_rinv,_rsq)          \
    _fscal = (6.0 * _Vvdw6 - _r * _Vvdwexp) * _rinvsq

  #define calc_vdw_force_(_rinv,_rsq) calc_vdw_only_force_(_rinv,_rsq)

#elif VDW == VDW_TAB

 #ifdef CONST_LJ

  #define calc_vdw_pot_(_rinv,_rsq,jnr)                   \
    _rt      = _rsq * (_rinv *_tabscale);                 \
    make_eps_(_rt,_n0,_eps);                             \
    nnn1     = TAB_MULT*((int)_rt) + LJTAB_OFS;           \
    _Y       = VFtab[nnn1+0];                             \
    _F       = VFtab[nnn1+1];                             \
    _G       = VFtab[nnn1+2];                             \
    _H       = VFtab[nnn1+3];                             \
    _GHeps   = _G + _H * _eps;                            \
    _VV      = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF      = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    Vvdwtot += _c6 * _VV;                                 \
    _fijD    = _c6 * _FF;                                 \
    _Y       = VFtab[nnn1+4];                             \
    _F       = VFtab[nnn1+5];                             \
    _G       = VFtab[nnn1+6];                             \
    _H       = VFtab[nnn1+7];                             \
    _GHeps   = _G + _H * _eps;                            \
    _VV      = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF      = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    Vvdwtot += _c12 * _VV;                                \
    _fijD   += _c12 * _FF

 #else

  #define calc_vdw_pot_(_rinv,_rsq,jnr)                   \
    tj       = nti+2*type[jnr];                           \
    _c6      = vdwparam[tj];                              \
    _c12     = vdwparam[tj+1];                            \
    _rt      = _rsq * (_rinv *_tabscale);                 \
    make_eps_(_rt,_n0,_eps);                             \
    nnn1     = TAB_MULT*((int)_rt) + LJTAB_OFS;           \
    _Y       = VFtab[nnn1+0];                             \
    _F       = VFtab[nnn1+1];                             \
    _G       = VFtab[nnn1+2];                             \
    _H       = VFtab[nnn1+3];                             \
    _GHeps   = _G + _H * _eps;                            \
    _VV      = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF      = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    Vvdwtot += _c6 * _VV;                                 \
    _fijD    = _c6 * _FF;                                 \
    _Y       = VFtab[nnn1+4];                             \
    _F       = VFtab[nnn1+5];                             \
    _G       = VFtab[nnn1+6];                             \
    _H       = VFtab[nnn1+7];                             \
    _GHeps   = _G + _H * _eps;                            \
    _VV      = _Y + _eps * (_F + _eps * _GHeps);          \
    _FF      = _F + _eps * (_GHeps + _GHeps + _eps * _H); \
    Vvdwtot += _c12 * _VV;                                \
    _fijD   += _c12 * _FF

 #endif

  #define calc_vdw_only_force_(_rinv,_rsq)     \
    _fscal  = _fijD * (_rinv * _tabscale)

  #define calc_vdw_force_(_rinv,_rsq) calc_vdw_only_force_(_rinv,_rsq)

#else

  #define calc_vdw_pot_(_rinv,_rsq,jnr) \
    Vvdwtot  = Vvdwtot                  \

  #define calc_vdw_only_force_(_rinv,_rsq) \
    _fscal   = 0.0

#if COULOMB == COULOMB_NONE
  #define calc_vdw_force_(_rinv,_rsq) calc_vdw_only_force_(_rinv,_rsq)
#else
  #define calc_vdw_force_(_rinv,_rsq) \
    {}
#endif

#endif


/* Based on the above conditionally-defined non-dual-form macros, we
 * can now paste together other macros that will actually be useful
 * and succinct to call from kernel functions to compute energies
 * and/or forces for pairs of interacting atoms.
 */

#define nfcalc_interaction_(_qq,_rinv,_rsq,jnr) \
  _rinvsq = _rinv * _rinv;                      \
  calc_coulomb_pot_(_qq,_rinv,_rsq,jnr);        \
  calc_vdw_pot_(_rinv,_rsq,jnr)


#define calc_interaction_(_qq,_rinv,_rsq,jnr) \
                                              \
  nfcalc_interaction_(_qq,_rinv,_rsq,jnr);    \
  calc_vdw_force_(_rinv,_rsq);                \
  calc_coulomb_force_(_qq,_rinv,_rsq)


#define nfcalc_coul_interaction_(_qq,_rinv,_rsq,jnr) \
  _rinvsq = _rinv * _rinv;                           \
  calc_coulomb_pot_(_qq,_rinv,_rsq,jnr)


#define calc_coul_interaction_(_qq,_rinv,_rsq,jnr) \
                                                   \
  nfcalc_coul_interaction_(_qq,_rinv,_rsq,jnr);    \
  calc_coulomb_only_force_(_qq,_rinv,_rsq)


#define nfcalc_vdw_interaction_(_rinv,_rsq,jnr) \
  _rinvsq = _rinv * _rinv;                      \
  calc_vdw_pot_(_rinv,_rsq,jnr)


#define calc_vdw_interaction_(_rinv,_rsq,jnr) \
                                              \
  nfcalc_vdw_interaction_(_rinv,_rsq,jnr);    \
  calc_vdw_only_force_(_rinv,_rsq)


/*****************************************************************************/
/*             BLUE GENE double hummer interaction calculations              */
/*****************************************************************************/

/* Macros for optimized table lookups */
#ifdef GMX_DOUBLE

#define dual_load  __lfpd
#define cross_load __lfxd

#else

#define dual_load  __lfps
#define cross_load __lfxs

#endif

/*
   Old comments that were left here by Mathias Puetz:

   trying to use quad loads for table lookups, but it's slower than
   double load version (probably because of new dependencies
   
    YF1      = __lfps(&VFtab[nnn1]);                         \
    YF2      = __lfxs(&VFtab[nnn2]);                         \
    GH1      = __lfps(&VFtab[nnn1+2]);                       \
    GH2      = __lfxs(&VFtab[nnn2+2]);                       \
    Y        = __fpsel(lr,YF1,YF2);                          \
    F        = __fxmr(__fpsel(rl,YF1,YF2));                  \
    G        = __fpsel(lr,GH1,GH2);                          \
    H        = __fpsel(rl,GH1,GH2);                          \
    GHeps    = __fxmadd(G,H,eps);                            \
    VV       = __fpmadd(Y,__fpmadd(F,GHeps,eps),eps);        \
    vctot    = __fpmadd(vctot,qq,VV)

    FF      = __fpmadd(F,__fxmadd(__fpadd(GHeps,GHeps),H,eps),eps); \

*/

/* This code takes care of calculating dual forms of eps. It requires
 * truncating floats to ints and converting the result back to
 * float. The int and float rounded forms have to occupy different
 * registers in different CPU units. There are three ways to do this.
 *
 * By default, do the integer form in the usual way, and use a
 * rounding trick that is known to be faster on BlueGene/L and
 * BlueGene. Otherwise, if the user configures the slow version,
 * because of the bug report in IssueID #429, use a different
 * slow version according to whether we are on BlueGene/L or
 * BlueGene/P.
 */

#if defined MAKE_EPS_SLOW_DUAL

#if defined __blrts__
/* We're on BlueGene/L, use moderately fast code */

#define make_eps(r,ri,conv,i1,i2,r0,eps)                   \
    ri      = __fpctiwz(r);                                \
    __stfpd(conv,ri);                                      \
    i1      = ((int *)conv)[1];                            \
    i2      = ((int *)conv)[3];                            \
    {};                                                    \
    {};                                                    \
    r0      = __cmplx((double)i1,(double)i2);              \
    eps     = __fpsub(r,r0)

#else /* We're not on BlueGene/L, use safe code */

#define make_eps(r,ri,conv,i1,i2,r0,eps)                   \
    i1      = (int)__creal(r);                             \
    i2      = (int)__cimag(r);                             \
    r0      = __cmplx((double)i1,(double)i2);              \
    eps     = __fpsub(r,r0)

#endif

#else /* use the rounding trick for speed */

#define make_eps(r,ri,conv,i1,i2,r0,eps)                   \
    ri      = __fpctiwz(r);                                \
    __stfpd(conv,ri);                                      \
    i1      = ((int *)conv)[1];                            \
    i2      = ((int *)conv)[3];                            \
    r0      = __fpsub(r,half);                             \
    r0      = __fpadd(r0,round);                           \
    r0      = __fpsub(r0,round);                           \
    eps     = __fpsub(r,r0)

#endif


/* Macros that will evaluate Coulomb energies and forces for dual
 * variables. The code has to change slightly in cases where there are
 * no VDW interactions to compute.
 */

#if COULOMB == COULOMB_TAB

/* Various alternative forms for this code have been developed
 * and tested. However, the BlueGene compilers are works-in-progress,
 * so retaining the old versions may be useful in the future.
 */

#if defined COULOMB_TAB_PUETZ

#define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)                    \
                                                                        \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                             \
    nnn1    *= TAB_MULT;                                                \
    nnn2    *= TAB_MULT;                                                \
    {};                                                                 \
    {};                                                                 \
    {};                                                                 \
    {};                                                                 \
    Y        = __cmplx(VFtab[nnn1+0],VFtab[nnn2+0]);                    \
    F        = __cmplx(VFtab[nnn1+1],VFtab[nnn2+1]);                    \
    G        = __cmplx(VFtab[nnn1+2],VFtab[nnn2+2]);                    \
    H        = __cmplx(VFtab[nnn1+3],VFtab[nnn2+3]);                    \
    {};                                                                 \
    GHeps    = __fpmadd(G,eps,H);                                       \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));                   \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H));    \
    fijC     = __fpmul(qq,FF);                                          \
    vctot    = __fpmadd(vctot,VV,qq)

#elif defined COULOMB_TAB_FPSEL

/* The dual variable H is nonlocal to the others in the code below. */

#define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)                    \
                                                                        \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                             \
    nnn1    *= TAB_MULT;                                                \
    nnn2    *= TAB_MULT;                                                \
    buf1     = dual_load(VFtab+nnn1);                                   \
    buf2     = cross_load(VFtab+nnn2);                                  \
    buf3     = dual_load(VFtab+nnn1+2);                                 \
    buf4     = cross_load(VFtab+nnn2+2);                                \
    Y        = __fpsel(lr,buf1,buf2);                                   \
    F        = __fxmr(__fpsel(rl,buf1,buf2));                           \
    G        = __fpsel(lr,buf3,buf4);                                   \
    H        = __fpsel(rl,buf3,buf4);                                   \
    {};                                                                 \
    GHeps    = __fxmadd(G,H,eps);                                       \
    VV       = __fpmadd(Y,__fpmadd(F,GHeps,eps),eps);                   \
    FF       = __fpmadd(F,__fxmadd(__fpadd(GHeps,GHeps),H,eps),eps);    \
    fijC     = __fpmul(qq,FF);                                          \
    vctot    = __fpmadd(vctot,qq,VV)

#elif defined COULOMB_TAB_EPSX

/* The dual variables F, H and FF are nonlocal to the others in the
 * code below, and buf1 has an copy of eps of opposite locality.
 */

#define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)                    \
                                                                        \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                             \
    nnn1    *= TAB_MULT;                                                \
    nnn2    *= TAB_MULT;                                                \
    buf1     = dual_load(VFtab+nnn1+0);                                 \
    buf2     = cross_load(VFtab+nnn2+0);                                \
    buf3     = dual_load(VFtab+nnn1+2);                                 \
    buf4     = cross_load(VFtab+nnn2+2);                                \
    Y        = __cmplx(__creal(buf1), __cimag(buf2));                   \
    F        = __cmplx(__creal(buf2), __cimag(buf1));                   \
    G        = __cmplx(__creal(buf3), __cimag(buf4));                   \
    H        = __cmplx(__creal(buf4), __cimag(buf3));                   \
    buf1     = __fxmr(eps);                                             \
    GHeps    = __fxmadd(G,H,eps);                                       \
    VV       = __fxmadd(Y,__fxmadd(F,GHeps,buf1),eps);                  \
    FF       = __fxmadd(F,__fxmadd(__fpadd(GHeps,GHeps),H,eps),buf1);   \
    fijC     = __fxmul(FF,qq);                                          \
    vctot    = __fpmadd(vctot,VV,qq)

#elif defined COULOMB_TAB_FX

/* The dual variables F and H are nonlocal to the others in the code
 * below, and buf1 has an copy of F of opposite locality. */

/* This code breaks on BlueGene/L at -O3. */

#define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)                    \
                                                                        \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                             \
    nnn1    *= TAB_MULT;                                                \
    nnn2    *= TAB_MULT;                                                \
    buf1      = dual_load(VFtab+nnn1+0);                                \
    buf2      = cross_load(VFtab+nnn2+0);                               \
    buf3      = dual_load(VFtab+nnn1+2);                                \
    buf4      = cross_load(VFtab+nnn2+2);                               \
    Y        = __cmplx(__creal(buf1), __cimag(buf2));                   \
    F        = __cmplx(__creal(buf2), __cimag(buf1));                   \
    G        = __cmplx(__creal(buf3), __cimag(buf4));                   \
    H        = __cmplx(__creal(buf4), __cimag(buf3));                   \
    buf1     = __fxmr(F);                                               \
    GHeps    = __fxmadd(G,H,eps);                                       \
    VV       = __fpmadd(Y,__fpmadd(buf1,GHeps,eps),eps);                \
    FF       = __fpmadd(buf1,__fxmadd(__fpadd(GHeps,GHeps),H,eps),eps); \
    fijC     = __fpmulqq(qq,FF);                                        \
    vctot    = __fpmaddqq(vctot,VV,qq)

#else /* default to the fastest known code */

/* The dual variable H is nonlocal to the others in the code below. */

#define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)                    \
                                                                        \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                             \
    nnn1    *= TAB_MULT;                                                \
    nnn2    *= TAB_MULT;                                                \
    {};                                                                 \
    {};                                                                 \
    buf3     = dual_load(VFtab+nnn1+2);                                 \
    buf4     = cross_load(VFtab+nnn2+2);                                \
    Y        = __cmplx(VFtab[nnn1+0],VFtab[nnn2+0]);                    \
    F        = __cmplx(VFtab[nnn1+1],VFtab[nnn2+1]);                    \
    G        = __cmplx(__creal(buf3), __cimag(buf4));                   \
    H        = __cmplx(__creal(buf4), __cimag(buf3));                   \
    {};                                                                 \
    GHeps    = __fxmadd(G,H,eps);                                       \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));                   \
    FF       = __fpmadd(F,eps,__fxmadd(__fpadd(GHeps,GHeps),H,eps));    \
    fijC     = __fpmul(qq,FF);                                          \
    vctot    = __fpmadd(vctot,VV,qq)

#endif

  #define calc_coulomb_only_force(qq,rinv,rsq)                      \
    fscal   = __fpmul(fijC,__fxpmul(rinv,_tabscale))

#if VDW == VDW_NONE
  #define calc_coulomb_force(qq,rinv,rsq) calc_coulomb_only_force(qq,rinv,rsq)
#else
  #define calc_coulomb_force(qq,rinv,rsq)                           \
    fscal   = __fpmadd(fscal,fijC,__fxpmul(rinv,_tabscale))
#endif

#elif COULOMB == GENERALIZED_BORN

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)               \
                                                                     \
    r        = __fpmul(rinv,rsq);                                    \
    isaprod  = __cmplx(invsqrta[jnr1],invsqrta[jnr2]);               \
    isaprod  = __fxpmul(isaprod,_isai);                              \
    iqq      = __fpmul(isaprod,qq);                                  \
    gbscale  = __fxpmul(isaprod,_gbtabscale);                        \
    rt       = __fpmul(r,gbscale);                                   \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                          \
    nnn1    *= 4;                                                    \
    nnn2    *= 4;                                                    \
    Y        = __cmplx(GBtab[nnn1+0],GBtab[nnn2+0]);                 \
    F        = __cmplx(GBtab[nnn1+1],GBtab[nnn2+1]);                 \
    G        = __cmplx(GBtab[nnn1+2],GBtab[nnn2+2]);                 \
    H        = __cmplx(GBtab[nnn1+3],GBtab[nnn2+3]);                 \
    r        = __fpmul(rinv,rsq);                                    \
    GHeps    = __fpmadd(G,eps,H);                                    \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));                \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    vcoul    = __fpmul(iqq,VV);                                      \
    vctot    = __fpmadd(vctot,iqq,VV);                               \
    fijC     = __fpmul(iqq,FF);                                      \
    dvdatmp  = __fpmadd(vcoul,fijC,rt);                              \
    dvdasum  = __fpsub(dvdasum,dvdatmp);                             \
    dvdaj    = __fpsub(dvdaj,dvdatmp)

  #define calc_coulomb_only_force(qq,rinv,rsq)                      \
    fscal   = __fpmul(fijC,__fpmul(rinv,gbscale))

#if VDW == VDW_NONE
  #define calc_coulomb_force(qq,rinv,rsq) calc_coulomb_only_force(qq,rinv,rsq)
#else
  #define calc_coulomb_force(qq,rinv,rsq)                           \
    fscal   = __fpmadd(fscal,fijC,__fpmul(rinv,gbscale))
#endif

#elif COULOMB == REACTION_FIELD

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)            \
                                                                  \
    krsq    = __fxpmul(rsq,_krf);                                 \
    vctot   = __fpmadd(vctot,qq,__fpadd(rinv,krsq));              \
    vctot   = __fxcpnmsub(vctot,qq,_crf)

  #define calc_coulomb_only_force(qq,rinv,rsq)                      \
    fscal   = __fpmul(__fpmul(qq,__fxcpmsub(rinv,krsq,2.0)),rinvsq)

#if VDW == VDW_NONE
  #define calc_coulomb_force(qq,rinv,rsq) calc_coulomb_only_force(qq,rinv,rsq)
#else
  #define calc_coulomb_force(qq,rinv,rsq)                           \
    fscal   = __fpmadd(fscal,__fpmul(qq,__fxcpmsub(rinv,krsq,2.0)),rinvsq)
#endif

#elif COULOMB == COULOMB_CUTOFF

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                       \
    vctot   = __fpmadd(vctot,qq,rinv)

  #define calc_coulomb_only_force(qq,rinv,rsq)                      \
    fscal   = __fpneg(__fpmul(__fpmul(qq,rinv),rinvsq))

#if VDW == VDW_NONE
  #define calc_coulomb_force(qq,rinv,rsq) calc_coulomb_only_force(qq,rinv,rsq)
#else
  #define calc_coulomb_force(qq,rinv,rsq)                           \
    fscal   = __fpnmsub(fscal,__fpmul(qq,rinv),rinvsq)
#endif

#else

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                       \
    vctot = vctot

  #define calc_coulomb_only_force(qq,rinv,rsq)                      \
    fscal   = zero

#if VDW == VDW_NONE
  #define calc_coulomb_force(qq,rinv,rsq) calc_coulomb_only_force(qq,rinv,rsq)
#else
  #define calc_coulomb_force(qq,rinv,rsq)                           \
    {}
#endif

#endif


/* Macros that will evaluate VDW energies and forces for non-dual
 * variables. In cases where the VDW parameters are constant
 * (e.g. water loops), different forms of these codes are faster.
 */

#if VDW == LENNARD_JONES

 #ifdef CONST_LJ

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)             \
                                                            \
    rinvsix  = __fpmul(rinvsq,__fpmul(rinvsq,rinvsq));      \
    Vvdw6    = __fxpmul(rinvsix,_c6);                       \
    Vvdw12   = __fxpmul(__fpmul(rinvsix,rinvsix),_c12);     \
    Vvdwtot  = __fpadd(Vvdwtot,__fpsub(Vvdw12,Vvdw6))

 #else

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)             \
                                                            \
    tj1      = nti+2*type[jnr1];                            \
    tj2      = nti+2*type[jnr2];                            \
    c6       = __cmplx(vdwparam[tj1],vdwparam[tj2]);        \
    c12      = __cmplx(vdwparam[tj1+1],vdwparam[tj2+1]);    \
    rinvsix  = __fpmul(rinvsq,__fpmul(rinvsq,rinvsq));      \
    Vvdw6    = __fpmul(rinvsix,c6);                         \
    Vvdw12   = __fpmul(__fpmul(rinvsix,rinvsix),c12);       \
    Vvdwtot  = __fpadd(Vvdwtot,__fpsub(Vvdw12,Vvdw6))

 #endif

  #define calc_vdw_only_force(rinv,rsq)                                   \
    fscal   = __fpmul(__fxcpmsub(__fxpmul(Vvdw12,12.0),Vvdw6,6.0),rinvsq)

  #define calc_vdw_force(rinv,rsq) calc_vdw_only_force(rinv,rsq)

#elif VDW == BUCKINGHAM

 #ifdef CONST_LJ

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)                             \
                                                                            \
    Vvdw6    = __fpmul(__fxpmul(rinvsq,_c6),__fpmul(rinvsq,rinvsq));        \
    r        = __fxpmul(__fpmul(rinv,rsq),_cexp2);                          \
    Vvdwexp  = __fxpmul(__cmplx(exp(__creal(-r)),exp(__cimag(-r))),_cexp1); \
    Vvdwtot  = __fpadd(Vvdwtot,__fpsub(Vvdwexp,Vvdw6))

 #else

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)                           \
                                                                          \
    tj1      = nti+3*type[jnr1];                                          \
    tj2      = nti+3*type[jnr2];                                          \
    c6       = __cmplx(vdwparam[tj1],vdwparam[tj2]);                      \
    cexp1    = __cmplx(vdwparam[tj1+1],vdwparam[tj2+1]);                  \
    cexp2    = __cmplx(vdwparam[tj1+2],vdwparam[tj2+2]);                  \
    Vvdw6    = __fpmul(__fpmul(rinvsq,c6),__fpmul(rinvsq,rinvsq));        \
    r        = __fpmul(__fpmul(rinv,rsq),cexp2);                          \
    Vvdwexp  = __fpmul(cexp1,__cmplx(exp(__creal(-r)),exp(__cimag(-r)))); \
    Vvdwtot  = __fpadd(Vvdwtot,__fpsub(Vvdwexp,Vvdw6))

 #endif

  #define calc_vdw_only_force(rinv,rsq)                                   \
    fscal   = __fpmul(__fxcpmsub(__fpmul(Vvdwexp,r),Vvdw6,6.0),rinvsq)

  #define calc_vdw_force(rinv,rsq) calc_vdw_only_force(rinv,rsq)

#elif VDW == VDW_TAB

 #ifdef CONST_LJ

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)              \
                                                             \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));        \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                  \
    nnn1     = nnn1*TAB_MULT + LJTAB_OFS;                    \
    nnn2     = nnn2*TAB_MULT + LJTAB_OFS;                    \
    Y        = __cmplx(VFtab[nnn1+0],VFtab[nnn2+0]);         \
    F        = __cmplx(VFtab[nnn1+1],VFtab[nnn2+1]);         \
    G        = __cmplx(VFtab[nnn1+2],VFtab[nnn2+2]);         \
    H        = __cmplx(VFtab[nnn1+3],VFtab[nnn2+3]);         \
    GHeps    = __fpmadd(G,eps,H);                            \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));        \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    Vvdwtot  = __fxcpmadd(Vvdwtot,VV,_c6);                   \
    fijD     = __fxpmul(FF,_c6);                             \
    Y        = __cmplx(VFtab[nnn1+4],VFtab[nnn2+4]);         \
    F        = __cmplx(VFtab[nnn1+5],VFtab[nnn2+5]);         \
    G        = __cmplx(VFtab[nnn1+6],VFtab[nnn2+6]);         \
    H        = __cmplx(VFtab[nnn1+7],VFtab[nnn2+7]);         \
    GHeps    = __fpmadd(G,eps,H);                            \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));        \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    Vvdwtot  = __fxcpmadd(Vvdwtot,VV,_c12);                  \
    fijD     = __fxcpmadd(fijD,FF,_c12);

 #else

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)              \
                                                             \
    tj1      = nti+2*type[jnr1];                             \
    tj2      = nti+2*type[jnr2];                             \
    c6       = __cmplx(vdwparam[tj1],vdwparam[tj2]);         \
    c12      = __cmplx(vdwparam[tj1+1],vdwparam[tj2+1]);     \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));        \
    make_eps(rt,rti,conv,nnn1,nnn2,n0,eps);                  \
    nnn1     = nnn1*TAB_MULT + LJTAB_OFS;                    \
    nnn2     = nnn2*TAB_MULT + LJTAB_OFS;                    \
    Y        = __cmplx(VFtab[nnn1+0],VFtab[nnn2+0]);         \
    F        = __cmplx(VFtab[nnn1+1],VFtab[nnn2+1]);         \
    G        = __cmplx(VFtab[nnn1+2],VFtab[nnn2+2]);         \
    H        = __cmplx(VFtab[nnn1+3],VFtab[nnn2+3]);         \
    GHeps    = __fpmadd(G,eps,H);                            \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));        \
    Vvdwtot  = __fpmadd(Vvdwtot,VV,c6);                      \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    fijD     = __fpmul(FF,c6);                               \
    Y        = __cmplx(VFtab[nnn1+4],VFtab[nnn2+4]);         \
    F        = __cmplx(VFtab[nnn1+5],VFtab[nnn2+5]);         \
    G        = __cmplx(VFtab[nnn1+6],VFtab[nnn2+6]);         \
    H        = __cmplx(VFtab[nnn1+7],VFtab[nnn2+7]);         \
    GHeps    = __fpmadd(G,eps,H);                            \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));        \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    Vvdwtot  = __fpmadd(Vvdwtot,c12,VV);                       \
    fijD     = __fpmadd(fijD,FF,c12)

 #endif

  #define calc_vdw_only_force(rinv,rsq)              \
    fscal   = __fpmul(fijD,__fxpmul(rinv,_tabscale))

  #define calc_vdw_force(rinv,rsq) calc_vdw_only_force(rinv,rsq)

#else /* VDW == VDW_NONE */

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2) \
                                                \
    Vvdwtot = Vvdwtot

  #define calc_vdw_only_force(rinv,rsq) \
    fscal   = zero

#if COULOMB == COULOMB_NONE
  #define calc_vdw_force(rinv,rsq) calc_vdw_only_force(rinv,rsq)
#else
  #define calc_vdw_force(rinv,rsq) \
    {}
#endif

#endif


/* Based on the above conditionally-defined dual-form macros, we can
 * now paste together other macros that will actually be useful and
 * succinct to call from kernel functions to compute energies and/or
 * forces for pairs of interacting atoms.
 */

#define nfcalc_interaction(qq,rinv,rsq,conv,jnr1,jnr2) \
  rinvsq   = __fpmul(rinv,rinv);                       \
  calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2);        \
  calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)


#define calc_interaction(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                     \
  nfcalc_interaction(qq,rinv,rsq,conv,jnr1,jnr2);    \
  calc_vdw_force(rinv,rsq);                          \
  calc_coulomb_force(qq,rinv,rsq)


#define nfcalc_coul_interaction(qq,rinv,rsq,conv,jnr1,jnr2) \
  rinvsq   = __fpmul(rinv,rinv);                            \
  calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)


#define calc_coul_interaction(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                          \
  nfcalc_coul_interaction(qq,rinv,rsq,conv,jnr1,jnr2);    \
  calc_coulomb_only_force(qq,rinv,rsq)


#define nfcalc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2) \
  rinvsq   = __fpmul(rinv,rinv);                        \
  calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)


#define calc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2) \
                                                      \
  nfcalc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2);    \
  calc_vdw_only_force(rinv,rsq)


/* Now some other helper macros. */

/* calulates the sqare of the euclidian norms of a two 3-dimensional vectors in parallel
 *
 * INPUT
 *
 *  x,y,z : components of the 3-vector   (vector)
 *
 * RETURNS
 *
 *  square of the euclidian norm         (vector)
 *
 */
#define vdist2(x,y,z) __fpmadd(__fpmadd(__fpmul(x,x),y,y),z,z)


/* refines two reciprocal square roots rinv of rsq with one Newton-Raphson iteration in parallel
 */
#define sqrt_newton(rinv,rsq) __fpmul(__fpmul(rinv,half),__fpnmsub(three,rsq,__fpmul(rinv,rinv)))


/* refines the reciprocal square root rinv of rsq with one Newton-Raphson iteration (scalar version)
 */
#define sqrt_newton_scalar(rinv,rsq) ((_half * rinv) * (3.0 - rsq * (rinv * rinv)));


/* refines two reciprocal estimates rinv of r with one Newton-Raphson iteration in parallel
 */
#define reci_newton(rinv,r) __fpmadd(rinv,rinv,__fpnmsub(one,r,rinv))


/* refines the reciprocal estimate rinv of r with one Newton-Raphson iteration (scalar version)
 */
#define reci_newton_scalar(rinv,r) (rinv + rinv*(1.0 - r*rinv))


#define load6_and_merge(v,i,j,x,y,z) \
  x = __cmplx(v[i]  ,v[j  ]);      \
  y = __cmplx(v[i+1],v[j+1]);      \
  z = __cmplx(v[i+2],v[j+2])


#define load3_and_clone(v,i,x,y,z) \
  x = __cmplx(v[i  ],v[i  ]);      \
  y = __cmplx(v[i+1],v[i+1]);      \
  z = __cmplx(v[i+2],v[i+2])


#define load3_real(v,i,x,y,z)   \
  x = __cmplx(v[i  ],0.0);      \
  y = __cmplx(v[i+1],0.0);      \
  z = __cmplx(v[i+2],0.0)


#define load3_imag(v,i,x,y,z)   \
  x = __cmplx(0.0,v[i  ]);      \
  y = __cmplx(0.0,v[i+1]);      \
  z = __cmplx(0.0,v[i+2])


#define real_store3(v,i,x,y,z) \
      v[i  ] = __creal(x);     \
      v[i+1] = __creal(y);     \
      v[i+2] = __creal(z)


#define imag_store3(v,i,x,y,z) \
      v[i]   = __cimag(x);     \
      v[i+1] = __cimag(y);     \
      v[i+2] = __cimag(z)


#define sum_and_store3(v,i,x,y,z)   \
  v[i  ] = __creal(x) + __cimag(x); \
  v[i+1] = __creal(y) + __cimag(y); \
  v[i+2] = __creal(z) + __cimag(z)


#define sum_and_add3(v,i,x,y,z)   \
  v[i  ] += __creal(x) + __cimag(x); \
  v[i+1] += __creal(y) + __cimag(y); \
  v[i+2] += __creal(z) + __cimag(z)


#define split_and_store6(v,i,j,x,y,z)  \
  real_store3(v,i,x,y,z);              \
  imag_store3(v,j,x,y,z)

#endif /* _INTERACTION_H_ */

