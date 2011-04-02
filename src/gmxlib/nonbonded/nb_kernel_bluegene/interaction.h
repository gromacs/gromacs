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


#ifndef _INTERACTION_H_
#define _INTERACTION_H_

/*****************************************************************************/
/*                        scalar interaction calculations                    */
/*****************************************************************************/

#if COULOMB == COULOMB_TAB

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr)          \
    _rt     = (_rsq * _rinv) *_tabscale;                 \
    _n0     = (real)((int)_rt);                          \
    _eps    = _rt-_n0;                                   \
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

  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal += _fijC * (_rinv * _tabscale)

#elif COULOMB == GENERALIZED_BORN

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr)              \
    _r          = _rinv * _rsq;                              \
    _isaprod    = _isai * invsqrta[jnr];                     \
     _iqq       = _qq * _isaprod;                            \
    _gbscale    = _isaprod * _gbtabscale;                    \
    _rt         = _r * _gbscale;                             \
    _n0         = (real)((int)_rt);                          \
    _eps        = _rt-_n0;                                   \
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

  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal += _fijC * (_rinv * _gbscale)

#elif COULOMB == REACTION_FIELD

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) \
    _krsq   = _krf * _rsq;                      \
    vctot  += _qq * (_rinv + _krsq - _crf)

  #define calc_coulomb_force_(_qq,_rinv,_rsq)       \
    _fscal -= _qq * (_rinv - 2.0 * _krsq) * _rinvsq

#elif COULOMB == COULOMB_CUTOFF

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) \
    vctot += _qq * _rinv 

  #define calc_coulomb_force_(_qq,_rinv,_rsq) \
    _fscal -= _qq * _rinv * _rinvsq

#else

  #define calc_coulomb_pot_(_qq,_rinv,_rsq,jnr) vctot = vctot
  #define calc_coulomb_force_(_qq,_rinv,_rsq) _fscal = _fscal

#endif

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

  #define calc_vdw_force_(_rinv,_rsq)               \
    _fscal = (6.0 * _Vvdw6 - 12.0 * _Vvdw12) * _rinvsq

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

  #define calc_vdw_force_(_rinv,_rsq)                \
    _fscal = (6.0 * _Vvdw6 - _r * _Vvdwexp) * _rinvsq

#elif VDW == VDW_TAB

 #ifdef CONST_LJ

  #define calc_vdw_pot_(_rinv,_rsq,jnr)                   \
    _rt      = _rsq * (_rinv *_tabscale);                 \
    _n0      = (real)((int)_rt);                          \
    _eps     = _rt-_n0;                                   \
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
    _n0      = (real)((int)_rt);                          \
    _eps     = _rt-_n0;                                   \
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

  #define calc_vdw_force_(_rinv,_rsq)     \
    _fscal  = _fijD * (_rinv * _tabscale)

#else

  #define calc_vdw_pot_(_rinv,_rsq,jnr) \
    Vvdwtot  = Vvdwtot                  \

  #define calc_vdw_force_(_rinv,_rsq) \
    _fscal = 0.0

#endif


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
  _fscal = 0.0;                                    \
  calc_coulomb_force_(_qq,_rinv,_rsq)


#define nfcalc_vdw_interaction_(_rinv,_rsq,jnr) \
  _rinvsq = _rinv * _rinv;                      \
  calc_vdw_pot_(_rinv,_rsq,jnr)


#define calc_vdw_interaction_(_rinv,_rsq,jnr) \
                                              \
  nfcalc_vdw_interaction_(_rinv,_rsq,jnr);    \
  calc_vdw_force_(_rinv,_rsq)


/*****************************************************************************/
/*             BLUE GENE double hummer interaction calculations              */
/*****************************************************************************/

/*
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

/* The optimized version of converts2ints is disabled on BG/P
 * because of issues on BG/P reported in bugzilla 429
 */
#if defined __blrts__

#define convert2ints(x,xi,conv,i1,i2)                      \
    xi      = __fpctiwz(x);                                \
    __stfpd(conv,xi);                                      \
    i1     = ((int *)conv)[1];                             \
    i2     = ((int *)conv)[3]

#else

#define convert2ints(x,xi,conv,i1,i2)                     \
    i1     = (int)__creal(x);                             \
    i2     = (int)__cimag(x)

#endif

#if COULOMB == COULOMB_TAB

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)       \
                                                             \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));        \
    convert2ints(rt,rti,conv,nnn1,nnn2);                     \
    n0       = __cmplx((double)nnn1,(double)nnn2);           \
    eps      = __fpsub(rt,n0);                               \
    nnn1    *= TAB_MULT;                                     \
    nnn2    *= TAB_MULT;                                     \
    Y        = __cmplx(VFtab[nnn1+0],VFtab[nnn2+0]);         \
    F        = __cmplx(VFtab[nnn1+1],VFtab[nnn2+1]);         \
    G        = __cmplx(VFtab[nnn1+2],VFtab[nnn2+2]);         \
    H        = __cmplx(VFtab[nnn1+3],VFtab[nnn2+3]);         \
    GHeps    = __fpmadd(G,eps,H);                            \
    VV       = __fpmadd(Y,eps,__fpmadd(F,eps,GHeps));        \
    FF       = __fpmadd(F,eps,__fpmadd(__fpadd(GHeps,GHeps),eps,H)); \
    fijC    = __fpmul(qq,FF);                                       \
    vctot    = __fpmadd(vctot,qq,VV)


  #define calc_coulomb_force(qq,rinv,rsq)                           \
                                                                    \
    fscal   = __fpmadd(fscal,fijC,__fxpmul(rinv,_tabscale))

#elif COULOMB == GENERALIZED_BORN

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)               \
                                                                     \
    r        = __fpmul(rinv,rsq);                                    \
    isaprod  = __cmplx(invsqrta[jnr1],invsqrta[jnr2]);               \
    isaprod  = __fxpmul(isaprod,_isai);                              \
    iqq      = __fpmul(isaprod,qq);                                  \
    gbscale  = __fxpmul(isaprod,_gbtabscale);                        \
    rt       = __fpmul(r,gbscale);                                   \
    convert2ints(rt,rti,conv,nnn1,nnn2);                             \
    n0       = __cmplx((double)nnn1,(double)nnn2);                   \
    eps      = __fpsub(rt,n0);                                       \
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

  #define calc_coulomb_force(qq,rinv,rsq) \
                                          \
    fscal   = __fpmadd(fscal,fijC,__fpmul(rinv,gbscale))

#elif COULOMB == REACTION_FIELD

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2)            \
                                                                  \
    krsq    = __fxpmul(rsq,_krf);                                 \
    vctot   = __fpmadd(vctot,qq,__fpadd(rinv,krsq));              \
    vctot   = __fxcpnmsub(vctot,qq,_crf)

  #define calc_coulomb_force(qq,rinv,rsq)                                  \
                                                                           \
    fscal   = __fpmadd(fscal,__fpmul(qq,__fxcpmsub(rinv,krsq,2.0)),rinvsq)

#elif COULOMB == COULOMB_CUTOFF

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                       \
    vctot   = __fpmadd(vctot,qq,rinv)

  #define calc_coulomb_force(qq,rinv,rsq)              \
                                                       \
    fscal   = __fpnmsub(fscal,__fpmul(qq,rinv),rinvsq)

#else

  #define calc_coulomb_pot(qq,rinv,rsq,conv,jnr1,jnr2) \
                                                       \
    vctot = vctot

  #define calc_coulomb_force(qq,rinv,rsq) \
                                          \
    fscal = fscal


#endif


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

  #define calc_vdw_force(rinv,rsq)                                        \
                                                                          \
    fscal   = __fpmul(__fxcpmsub(__fxpmul(Vvdw12,12.0),Vvdw6,6.0),rinvsq)

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

  #define calc_vdw_force(rinv,rsq)                                        \
                                                                          \
    fscal   = __fpmul(__fxcpmsub(__fpmul(Vvdwexp,r),Vvdw6,6.0),rinvsq)

#elif VDW == VDW_TAB

 #ifdef CONST_LJ

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)              \
                                                             \
    rt       = __fpmul(rsq,__fxpmul(rinv,_tabscale));        \
    convert2ints(rt,rti,conv,nnn1,nnn2);                     \
    n0       = __cmplx((double)nnn1,(double)nnn2);           \
    eps      = __fpsub(rt,n0);                               \
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
    convert2ints(rt,rti,conv,nnn1,nnn2);                     \
    n0       = __cmplx((double)nnn1,(double)nnn2);           \
    eps      = __fpsub(rt,n0);                               \
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

  #define calc_vdw_force(rinv,rsq)                   \
                                                     \
    fscal   = __fpmul(fijD,__fxpmul(rinv,_tabscale))

#else

  #define calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2) \
                                                \
    Vvdwtot = Vvdwtot

  #define calc_vdw_force(rinv,rsq) \
                                   \
    fscal   = zero

#endif



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
  fscal = zero;                                           \
  calc_coulomb_force(qq,rinv,rsq)


#define nfcalc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2) \
  rinvsq   = __fpmul(rinv,rinv);                        \
  calc_vdw_pot(rinv,rsq,conv,jnr1,jnr2)


#define calc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2) \
                                                      \
  nfcalc_vdw_interaction(rinv,rsq,conv,jnr1,jnr2);    \
  fscal = zero;                                       \
  calc_vdw_force(rinv,rsq)



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
#define sqrt_newton(rinv,rsq) __fpmul(__fxpmul(rinv,0.5),__fpnmsub(three,rsq,__fpmul(rinv,rinv)))


/* refines the reciprocal square root rinv of rsq with one Newton-Raphson iteration (scalar version)
 */
#define sqrt_newton_scalar(rinv,rsq) ((0.5 * rinv) * (3.0 - rsq * (rinv * rinv)));


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

