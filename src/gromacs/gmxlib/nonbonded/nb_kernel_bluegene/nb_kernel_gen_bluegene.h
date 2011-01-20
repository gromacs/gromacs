#include "types/simple.h"

#undef FULL_INTERACTION
#undef FULL_INTERACTION_
#undef COUL_INTERACTION
#undef COUL_INTERACTION_

#ifdef NO_FORCE

 #define FULL_INTERACTION  nfcalc_interaction
 #define FULL_INTERACTION_ nfcalc_interaction_
 #define COUL_INTERACTION  nfcalc_coul_interaction
 #define COUL_INTERACTION_ nfcalc_coul_interaction_

#else

 #define FULL_INTERACTION  calc_interaction
 #define FULL_INTERACTION_ calc_interaction_
 #define COUL_INTERACTION  calc_coul_interaction
 #define COUL_INTERACTION_ calc_coul_interaction_

#endif

void NB_KERNEL (
                    int *    p_nri,
                    int *    iinr,
                    int *    jindex,
                    int *    jjnr,
                    int *    shift,
                    real *   shiftvec,
                    real *   fshift,
                    int *    gid,
                    real *   pos,
                    real *   faction,
                    real *   charge,
                    real *   p_facel,
                    real *   p_krf,
                    real *   p_crf,
                    real *   Vc,
                    int *    type,
                    int *    p_ntype,
                    real *   vdwparam,
                    real *   Vvdw,
                    real *   p_tabscale,
                    real *   VFtab,
                    real *   invsqrta,
                    real *   dvda,
                    real *   p_gbtabscale,
                    real *   GBtab,
                    int *    p_nthreads,
                    int *    count,
                    void *   mtx,
                    int *    outeriter,
                    int *    inneriter,
                    real *   work)
{
    double _Complex zero  = __cmplx(0.0,0.0);
    double _Complex rone  = __cmplx(1.0,0.0);
    double _Complex three = __cmplx(3.0,3.0);
    double conv1[2],conv2[2],conv3[2];

    real   _facel,_tabscale,_gbtabscale,_krf,_crf;
    int    nri,ntype,nthreads,n,ii,nj0,nj1;

#pragma disjoint(*shiftvec,*fshift,*pos,*faction,*charge,*p_facel,*p_krf,*p_crf,*Vc, \
                 *vdwparam,*Vvdw,*p_tabscale,*VFtab,*invsqrta,*dvda,*p_gbtabscale,*GBtab,*work)

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    _facel           = *p_facel;       
    _tabscale        = *p_tabscale;    
    _gbtabscale      = *p_gbtabscale;    
    _krf             = *p_krf;    
    _crf             = *p_crf;    
    nj1              = 0;              

    for(n=0; (n<nri); n++)
    {
	double _Complex ix,iy,iz;
	
	// initialize potential energies and forces for this paricle

        double _Complex vctot    = zero;              
        double _Complex Vvdwtot  = zero;
	double _Complex dvdasum  = zero;
	double _Complex fix      = zero;
	double _Complex fiy      = zero;
	double _Complex fiz      = zero;

        // shift is the position of the n-th water group

        int is3  = 3*shift[n];
	
	// shiftvec is the center of a water group

        real _shX  = shiftvec[is3];  
        real _shY  = shiftvec[is3+1];
        real _shZ  = shiftvec[is3+2];

        int ii  = iinr[n];        
        int ii3 = 3*ii;
#if VDW == BUCKINGHAM
        int nti = 3*ntype*type[ii];
#else
	int nti = 2*ntype*type[ii];
#endif
	int k,ggid;

	real _iq   = _facel * charge[ii];
	real _isai = invsqrta[ii];

	// add the shift vector to all water atoms

        real _ix = _shX + pos[ii3+0];
        real _iy = _shY + pos[ii3+1];
        real _iz = _shZ + pos[ii3+2];
	
	ix = __cmplx(_ix,_ix);
	iy = __cmplx(_iy,_iy);
	iz = __cmplx(_iz,_iz);

        nj0 = jindex[n];      
        nj1 = jindex[n+1];    

	/* unrolled 6 times : unrolled twice for SIMDization and three times to hide dependencies */

        for(k=nj0; (k<nj1-5); k+=6)
        {
	    double _Complex dx1,dy1,dz1,fjx1,fjy1,fjz1;
	    double _Complex rsq1,rinv1;
	    double _Complex dx2,dy2,dz2,fjx2,fjy2,fjz2;
	    double _Complex rsq2,rinv2;
	    double _Complex dx3,dy3,dz3,fjx3,fjy3,fjz3;
	    double _Complex rsq3,rinv3;

	    double _Complex rinvsq,krsq,fscal;
	    double _Complex rt,rti,r,eps,Y,F,G,H,GHeps,VV,FF,fijC,fijD,n0,qq;
	    double _Complex isaprod,isaj,iqq,gbscale,dvdaj,vcoul,dvdatmp;
	    double _Complex c6,c12,cexp1,cexp2,rinvsix,Vvdw6,Vvdw12,Vvdwexp;

	    int nnn1,nnn2;
	    int jnr11,jnr21,jnr12,jnr22,jnr13,jnr23,j11,j21,j12,j22,j13,j23,tj1,tj2;

	    jnr11  = jjnr[k+0];
	    jnr21  = jjnr[k+1];
	    jnr12  = jjnr[k+2];
	    jnr22  = jjnr[k+3];
	    jnr13  = jjnr[k+4];
	    jnr23  = jjnr[k+5];

            j11    = 3*jnr11;        
            j21    = 3*jnr21;
            j12    = 3*jnr12;        
            j22    = 3*jnr22;
            j13    = 3*jnr13;        
            j23    = 3*jnr23;

	    load6_and_merge(pos,j11,j21,dx1,dy1,dz1);
	    load6_and_merge(pos,j12,j22,dx2,dy2,dz2);
	    load6_and_merge(pos,j13,j23,dx3,dy3,dz3);

	    dx1   = __fpsub(ix,dx1);
	    dy1   = __fpsub(iy,dy1);
	    dz1   = __fpsub(iz,dz1);

	    dx2   = __fpsub(ix,dx2);
	    dy2   = __fpsub(iy,dy2);
	    dz2   = __fpsub(iz,dz2);

	    dx3   = __fpsub(ix,dx3);
	    dy3   = __fpsub(iy,dy3);
	    dz3   = __fpsub(iz,dz3);

#if COULOMB != COULOMB_NONE
	    qq   = __fxpmul(__cmplx(charge[jnr11],charge[jnr21]),_iq);

#if COULOMB == GENERALIZED_BORN
	    dvdaj = __cmplx(dvda[jnr11],dvda[jnr21]);
#endif
#endif

	    rsq1  = vdist2(dx1,dy1,dz1);
	    rsq2  = vdist2(dx2,dy2,dz2);
	    rsq3  = vdist2(dx3,dy3,dz3);

	    rinv1 = __fprsqrte(rsq1);
	    rinv2 = __fprsqrte(rsq2);
	    rinv3 = __fprsqrte(rsq3);

	    rinv1 = sqrt_newton(rinv1,rsq1);
	    rinv2 = sqrt_newton(rinv2,rsq2);
	    rinv3 = sqrt_newton(rinv3,rsq3);

#ifdef GMX_DOUBLE

	    rinv1 = sqrt_newton(rinv1,rsq1);
	    rinv2 = sqrt_newton(rinv2,rsq2);
	    rinv3 = sqrt_newton(rinv3,rsq3);

#endif

#ifndef NO_FORCE
	    load6_and_merge(faction,j11,j21,fjx1,fjy1,fjz1);
#endif
	    
	    FULL_INTERACTION(qq,rinv1,rsq1,conv1,jnr11,jnr21);

#if COULOMB != COULOMB_NONE
	    qq    = __fxpmul(__cmplx(charge[jnr12],charge[jnr22]),_iq);

#if COULOMB == GENERALIZED_BORN
	    dvda[jnr11] -= __creal(dvdaj);
	    dvda[jnr21] -= __creal(dvdaj);

	    dvdaj = __cmplx(dvda[jnr12],dvda[jnr22]);
#endif
#endif

#ifndef NO_FORCE
            fix   = __fpnmsub(fix,dx1,fscal);
            fiy   = __fpnmsub(fiy,dy1,fscal);      
            fiz   = __fpnmsub(fiz,dz1,fscal);

            fjx1  = __fpmadd(fjx1,dx1,fscal);
            fjy1  = __fpmadd(fjy1,dy1,fscal);
            fjz1  = __fpmadd(fjz1,dz1,fscal);

	    load6_and_merge(faction,j12,j22,fjx2,fjy2,fjz2);
#endif

	    FULL_INTERACTION(qq,rinv2,rsq2,conv2,jnr12,jnr22);

#if COULOMB != COULOMB_NONE
		qq    = __fxpmul(__cmplx(charge[jnr13],charge[jnr23]),_iq);
			 
#if COULOMB == GENERALIZED_BORN
	    dvda[jnr12] -= __creal(dvdaj);
	    dvda[jnr22] -= __creal(dvdaj);

	    dvdaj = __cmplx(dvda[jnr13],dvda[jnr23]);
#endif
#endif
			
#ifndef NO_FORCE
	    split_and_store6(faction,j11,j21,fjx1,fjy1,fjz1);

            fix   = __fpnmsub(fix,dx2,fscal);
            fiy   = __fpnmsub(fiy,dy2,fscal);      
            fiz   = __fpnmsub(fiz,dz2,fscal);

            fjx2  = __fpmadd(fjx2,dx2,fscal);
            fjy2  = __fpmadd(fjy2,dy2,fscal);
            fjz2  = __fpmadd(fjz2,dz2,fscal);

	    load6_and_merge(faction,j13,j23,fjx3,fjy3,fjz3);
#endif

	    FULL_INTERACTION(qq,rinv3,rsq3,conv3,jnr13,jnr23);

#if COULOMB == GENERALIZED_BORN
	    dvda[jnr13] -= __creal(dvdaj);
	    dvda[jnr23] -= __creal(dvdaj);
#endif

#ifndef NO_FORCE
	    split_and_store6(faction,j12,j22,fjx2,fjy2,fjz2);

            fix   = __fpnmsub(fix,dx3,fscal);
            fiy   = __fpnmsub(fiy,dy3,fscal);      
            fiz   = __fpnmsub(fiz,dz3,fscal);

            fjx3  = __fpmadd(fjx3,dx3,fscal);
            fjy3  = __fpmadd(fjy3,dy3,fscal);
            fjz3  = __fpmadd(fjz3,dz3,fscal);

	    split_and_store6(faction,j13,j23,fjx3,fjy3,fjz3);
#endif
        }

        for(;(k<nj1); k++)
        {
	    real _dx,_dy,_dz,_rsq,_rinv;
	    real _rinvsq,_krsq,_fscal;
	    real _r,_rt,_eps,_Y,_F,_G,_H,_GHeps,_VV,_FF,_fijC,_fijD,_n0,_qq;
	    real _isaprod,_isaj,_iqq,_dvdatmp,_vcoul,_gbscale;
	    real _c6,_c12,_cexp1,_cexp2,_rinvsix,_Vvdw6,_Vvdw12,_Vvdwexp;

	    int nnn1;
	    int jnr,j3,tj;

            jnr    = jjnr[k];        
            j3     = 3*jnr;          

            _dx    = __creal(ix) - pos[j3+0];
            _dy    = __creal(iy) - pos[j3+1];
            _dz    = __creal(iz) - pos[j3+2];

            _qq    = _iq * charge[jnr];

            _rsq   = _dx*_dx + _dy*_dy + _dz*_dz;

            _rinv  = __frsqrte(_rsq);

            _rinv  = sqrt_newton_scalar(_rinv,_rsq);

#ifdef GMX_DOUBLE
	    
            _rinv  = sqrt_newton_scalar(_rinv,_rsq);

#endif
	    FULL_INTERACTION_(_qq,_rinv,_rsq,jnr);

#ifndef NO_FORCE
            fix             -= _fscal*_dx; 
            fiy             -= _fscal*_dy; 
            fiz             -= _fscal*_dz; 
            faction[j3+0]   += _fscal*_dx;
            faction[j3+1]   += _fscal*_dy;
            faction[j3+2]   += _fscal*_dz;
#endif
        }
        
	ggid = gid[n];

#ifndef NO_FORCE
        sum_and_add3(faction,ii3,fix,fiy,fiz);
        sum_and_add3(fshift ,is3,fix,fiy,fiz);
#endif

#if COULOMB != COULOMB_NONE
        Vc[ggid]    += __creal(vctot)   + __cimag(vctot);
#if COULOMB == GENERALIZED_BORN
        dvda[ii]    += __creal(dvdasum) + __cimag(dvdasum);
#endif
#endif

#if VDW != VDW_NONE
        Vvdw[ggid]  += __creal(Vvdwtot) + __cimag(Vvdwtot);
#endif
    }

    *outeriter  = nri;            
    *inneriter  = nj1;         
}

