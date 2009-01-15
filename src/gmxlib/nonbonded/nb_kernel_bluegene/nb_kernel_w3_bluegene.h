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
    double _Complex qO,qH;
    double _Complex zero  = __cmplx(0.0,0.0);
    double _Complex rone  = __cmplx(1.0,0.0);
    double _Complex three = __cmplx(3.0,3.0);
    double conv1[2],conv2[2],conv3[2];

    real            _qO,_qH,_facel,_tabscale,_gbtabscale,_krf,_crf;

    int             nri,ntype,nthreads,n,ii,nj0,nj1,nti;

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
    ii               = iinr[0];        

    _qO              = _facel * charge[ii];     
    _qH              = _facel * charge[ii+1];   
    qO               = __cmplx(_qO,_qO);
    qH               = __cmplx(_qH,_qH);
#if VDW == BUCKINGHAM
    nti              = 3*ntype*type[ii];
#else
    nti              = 2*ntype*type[ii];
#endif
    
    nj1 = 0;

    for(n=0; (n<nri); n++)
    {
	double _Complex ix1,ix2,ix3,iy1,iy2,iy3,iz1,iz2,iz3,ix23,iy23,iz23,fix,fiy,fiz;
	
	// initialize potential energies and forces for this paricle

        double _Complex vctot     = zero;              
        double _Complex Vvdwtot   = zero;              
	double _Complex fix1      = zero;
	double _Complex fix2      = zero;
	double _Complex fix3      = zero;
	double _Complex fiy1      = zero;
	double _Complex fiy2      = zero;
	double _Complex fiy3      = zero;
	double _Complex fiz1      = zero;
	double _Complex fiz2      = zero;
	double _Complex fiz3      = zero;

        // shift is the position of the n-th water group

        int is3  = 3*shift[n];
	
	// shiftvec is the center of a water group

        real  _shX  = shiftvec[is3+0];  
        real  _shY  = shiftvec[is3+1];
        real  _shZ  = shiftvec[is3+2];

        int ii  = iinr[n];        
        int ii3 = 3*ii;
	int k,ggid;

	// add the shift vector to all water atoms

        real  _ix1 = _shX + pos[ii3+0];
        real  _iy1 = _shY + pos[ii3+1];
        real  _iz1 = _shZ + pos[ii3+2];
        real  _ix2 = _shX + pos[ii3+3];
        real  _iy2 = _shY + pos[ii3+4];
        real  _iz2 = _shZ + pos[ii3+5];
        real  _ix3 = _shX + pos[ii3+6];
        real  _iy3 = _shY + pos[ii3+7];
        real  _iz3 = _shZ + pos[ii3+8];

	// clone all positions in complex variables
	
	ix1 = __cmplx(_ix1,_ix1);
	iy1 = __cmplx(_iy1,_iy1);
	iz1 = __cmplx(_iz1,_iz1);
              
	ix2 = __cmplx(_ix2,_ix2);
	iy2 = __cmplx(_iy2,_iy2);
	iz2 = __cmplx(_iz2,_iz2);

	ix3 = __cmplx(_ix3,_ix3);
	iy3 = __cmplx(_iy3,_iy3);
	iz3 = __cmplx(_iz3,_iz3);

        nj0 = jindex[n];
        nj1 = jindex[n+1];

	// unrolled twice for SIMDization

        for(k=nj0; (k<nj1-1); k+=2)
        {
	    double _Complex dx11,dx21,dx31,jx,fjx;
	    double _Complex dy11,dy21,dy31,jy,fjy;
	    double _Complex dz11,dz21,dz31,jz,fjz;
	    double _Complex rsq11,rsq21,rsq31;
	    double _Complex rinv11,rinv21,rinv31;
	    double _Complex rinvsq,krsq,fscal;
	    double _Complex rt,rti,r,n0,eps,Y,F,G,H,GHeps,VV,FF,fijC,fijD,qj,qq;
	    double _Complex c6,c12,cexp1,cexp2,Vvdwexp,Vvdw6,Vvdw12,rinvsix;

	    int nnn1,nnn2;
	    int jnr1,jnr2,j1,j2,tj1,tj2;

	    jnr1  = jjnr[k];
	    jnr2  = jjnr[k+1];

            j1    = 3*jnr1;        
            j2    = 3*jnr2;

	    load6_and_merge(pos,j1,j2,jx,jy,jz);

	    dx11   = __fpsub(ix1,jx);
	    dy11   = __fpsub(iy1,jy);
	    dz11   = __fpsub(iz1,jz);
	    dx21   = __fpsub(ix2,jx);
	    dy21   = __fpsub(iy2,jy);
	    dz21   = __fpsub(iz2,jz);
	    dx31   = __fpsub(ix3,jx);
	    dy31   = __fpsub(iy3,jy);
	    dz31   = __fpsub(iz3,jz);

	    qj = __cmplx(charge[jnr1],charge[jnr2]);

	    rsq11  = vdist2(dx11,dy11,dz11);
	    rsq21  = vdist2(dx21,dy21,dz21);
	    rsq31  = vdist2(dx31,dy31,dz31);

	    rinv11 = __fprsqrte(rsq11);
	    rinv21 = __fprsqrte(rsq21);
	    rinv31 = __fprsqrte(rsq31);

	    rinv11 = sqrt_newton(rinv11,rsq11);
	    rinv21 = sqrt_newton(rinv21,rsq21);
            rinv31 = sqrt_newton(rinv31,rsq31);

#ifdef GMX_DOUBLE

	    rinv11 = sqrt_newton(rinv11,rsq11);
	    rinv21 = sqrt_newton(rinv21,rsq21);
            rinv31 = sqrt_newton(rinv31,rsq31);

#endif
	    qq = __fpmul(qO,qj);

#ifndef NO_FORCE
	    load6_and_merge(faction,j1,j2,fjx,fjy,fjz);
#endif

	    FULL_INTERACTION(qq,rinv11,rsq11,conv1,jnr1,jnr2);

	    qq = __fpmul(qH,qj);

#ifndef NO_FORCE
            fix1  = __fpnmsub(fix1,dx11,fscal);
            fiy1  = __fpnmsub(fiy1,dy11,fscal);      
            fiz1  = __fpnmsub(fiz1,dz11,fscal);

            fjx   = __fpmadd(fjx,dx11,fscal);
            fjy   = __fpmadd(fjy,dy11,fscal);
            fjz   = __fpmadd(fjz,dz11,fscal);
#endif

	    COUL_INTERACTION(qq,rinv21,rsq21,conv2,jnr1,jnr2);

#ifndef NO_FORCE
            fix2  = __fpnmsub(fix2,dx21,fscal);
            fiy2  = __fpnmsub(fiy2,dy21,fscal);      
            fiz2  = __fpnmsub(fiz2,dz21,fscal);

            fjx   = __fpmadd(fjx,dx21,fscal);
            fjy   = __fpmadd(fjy,dy21,fscal);
            fjz   = __fpmadd(fjz,dz21,fscal);
#endif

	    COUL_INTERACTION(qq,rinv31,rsq31,conv3,jnr1,jnr2);

#ifndef NO_FORCE
            fix3  = __fpnmsub(fix3,dx31,fscal);
            fiy3  = __fpnmsub(fiy3,dy31,fscal);      
            fiz3  = __fpnmsub(fiz3,dz31,fscal);

            fjx   = __fpmadd(fjx,dx31,fscal);
            fjy   = __fpmadd(fjy,dy31,fscal);
            fjz   = __fpmadd(fjz,dz31,fscal);

	    split_and_store6(faction,j1,j2,fjx,fjy,fjz);
#endif
        }

	ix23 = __cmplx(_ix2,_ix3);
	iy23 = __cmplx(_iy2,_iy3);
	iz23 = __cmplx(_iz2,_iz3);


	// actually we should not simdize the remainder loop, because it's slower

        for(;(k<nj1); k++)
        {
	    double _Complex dx23,dy23,dz23,tx,ty,tz,jx,jy,jz;
	    double _Complex fjx,fjy,fjz;
	    double _Complex rsq23,rinv23;
	    double _Complex rinvsq,krsq,fscal;
	    double _Complex rt,rti,r,eps,Y,F,G,H,GHeps,VV,FF,fijC,fijD,n0,qq;
	    
	    real _dx11,_dy11,_dz11,_rsq11,_rinv11;
	    real _rinvsq,_krsq,_fscal;
	    real _rt,_r,_eps,_Y,_F,_H,_G,_GHeps,_VV,_FF,_fijC,_fijD,_n0,_qq;
	    real _qj,_c6,_c12,_cexp1,_cexp2,_Vvdwexp,_Vvdw6,_Vvdw12,_rinvsix;

	    int nnn1,nnn2;
	    int jnr,j3,tj;

	    jnr    = jjnr[k];
            j3     = 3*jnr;

	    load3_and_clone(pos,j3,jx,jy,jz);

            _dx11    = __creal(ix1) - __creal(jx);
            _dy11    = __creal(iy1) - __creal(jy);
	    _dz11    = __creal(iz1) - __creal(jz);
	    dx23     = __fpsub(ix23,jx);
	    dy23     = __fpsub(iy23,jy);
 	    dz23     = __fpsub(iz23,jz);

            _rsq11   = _dx11 * _dx11 + _dy11 * _dy11 + _dz11 * _dz11;
	    rsq23    = vdist2(dx23,dy23,dz23);

	    _rinv11  = __frsqrte(_rsq11);
	    rinv23   = __fprsqrte(rsq23 );

	    _rinv11  = sqrt_newton_scalar(_rinv11,_rsq11);
            rinv23   = sqrt_newton(rinv23,rsq23);

#ifdef GMX_DOUBLE

	    _rinv11  = sqrt_newton_scalar(_rinv11,_rsq11);
            rinv23   = sqrt_newton(rinv23,rsq23);

#endif

	    _qq       = _qO * charge[jnr];

#ifndef NO_FORCE
	    load3_real(faction,j3,fjx,fjy,fjz);
#endif

	    FULL_INTERACTION_(_qq,_rinv11,_rsq11,jnr);

	    qq        = __fxpmul(qH,charge[jnr]);

#ifndef NO_FORCE
            fix1     -= _fscal * _dx11;
            fiy1     -= _fscal * _dy11;
            fiz1     -= _fscal * _dz11;

            fjx      += _fscal * _dx11;
            fjy      += _fscal * _dy11;
            fjz      += _fscal * _dz11;
#endif

	    COUL_INTERACTION(qq,rinv23,rsq23,conv1,jnr,jnr);

#ifndef NO_FORCE
	    tx    = __fpmul(fscal,dx23);
	    ty    = __fpmul(fscal,dy23);
	    tz    = __fpmul(fscal,dz23);

	    fix2  = __fpnmsub(fix2,rone,tx);
	    fiy2  = __fpnmsub(fiy2,rone,ty);
	    fiz2  = __fpnmsub(fiz2,rone,tz);
	    fix3  = __fxnmsub(fix3,rone,tx);
	    fiy3  = __fxnmsub(fiy3,rone,ty);
	    fiz3  = __fxnmsub(fiz3,rone,tz);

            fjx   = __fpadd(fjx,tx);
            fjy   = __fpadd(fjy,ty);
            fjz   = __fpadd(fjz,tz);

	    sum_and_store3(faction,j3,fjx,fjy,fjz);
#endif
	}

	ggid = gid[n];

#ifndef NO_FORCE
        sum_and_add3(faction,ii3  ,fix1,fiy1,fiz1);
        sum_and_add3(faction,ii3+3,fix2,fiy2,fiz2);
        sum_and_add3(faction,ii3+6,fix3,fiy3,fiz3);

	fix = __fpadd(fix1,__fpadd(fix2,fix3));
	fiy = __fpadd(fiy1,__fpadd(fiy2,fiy3));
	fiz = __fpadd(fiz1,__fpadd(fiz2,fiz3));

	sum_and_add3(fshift,is3,fix,fiy,fiz);
#endif

#if COULOMB != COULOMB_NONE
        Vc[ggid]   += __creal(vctot)   + __cimag(vctot);
#endif

#if VDW != VDW_NONE
        Vvdw[ggid] += __creal(Vvdwtot)   + __cimag(Vvdwtot);
#endif
    }
         
    *outeriter       = nri;            
    *inneriter       = nj1;            
}
