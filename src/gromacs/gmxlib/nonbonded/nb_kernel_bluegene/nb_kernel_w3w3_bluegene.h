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

                    int *   p_nri,
                    int *   iinr,
                    int *   jindex,
                    int *   jjnr,
                    int *   shift,
                    real *  shiftvec,
                    real *  fshift,
                    int *   gid,
                    real *  pos,
                    real *  faction,
                    real *  charge,
                    real *  p_facel,
                    real *  p_krf,
                    real *  p_crf,
                    real *  Vc,
                    int *   type,
                    int *   p_ntype,
                    real *  vdwparam,
                    real *  Vvdw,
                    real *  p_tabscale,
                    real *  VFtab,
                    real *  invsqrta,
                    real *  dvda,
                    real *  p_gbtabscale,
                    real *  GBtab,
                    int *   p_nthreads,
                    int *   count,
                    void *  mtx,
                    int *   outeriter,
                    int *   inneriter,
                    real *  work)
{
    double _Complex qqOO_OO,qqOO_OH,qqOH_OH,qqOH_HH,qqHH_HH;
    double _Complex zero  = __cmplx(0.0,0.0);
    double _Complex rone  = __cmplx(1.0,0.0);
    double _Complex two   = __cmplx(2.0,2.0);
    double _Complex lr    = __cmplx(-1.0,1.0);
    double _Complex rl    = __cmplx(1.0,-1.0);
    double _Complex three = __cmplx(3.0,3.0);

    double conv1[2],conv2[2],conv3[2],conv4[2];

    real  _qO,_qH,_qqOO,_qqOH,_qqHH;
    real  _facel,_tabscale,_gbtabscale,_krf,_crf,_c6,_c12,_cexp1,_cexp2;

    int  nri,nti,ntype,nthreads,n,ii,tj,nj0,nj1;

#pragma disjoint(*shiftvec,*fshift,*pos,*faction,*charge,*p_facel,*p_krf,*p_crf,*Vc, \
                 *vdwparam,*Vvdw,*p_tabscale,*VFtab,*invsqrta,*dvda,*p_gbtabscale,*GBtab,*work)

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    _facel           = *p_facel;
    _tabscale        = *p_tabscale;
    _krf             = *p_krf;
    _crf             = *p_crf;
    _gbtabscale      = *p_gbtabscale;    

    ii               = iinr[0];        

    _qO              = charge[ii];     
    _qH              = charge[ii+1];   
    _qqOO            = _facel*_qO*_qO;    
    _qqOH            = _facel*_qO*_qH;
    _qqHH            = _facel*_qH*_qH;    
    qqOO_OO          = __cmplx(_qqOO,_qqOO);
    qqOH_OH          = __cmplx(_qqOH,_qqOH);
    qqOO_OH          = __cmplx(_qqOO,_qqOH);
    qqOH_HH          = __cmplx(_qqOH,_qqHH);
    qqHH_HH          = __cmplx(_qqHH,_qqHH);

#if VDW == LENNARD_JONES || VDW == VDW_TAB
    tj               = 2*(ntype+1)*type[ii];
    _c6              = vdwparam[tj];   
    _c12             = vdwparam[tj+1]; 
#elif VDW == BUCKINGHAM
    tj               = 3*(ntype+1)*type[ii];
    _c6              = vdwparam[tj];
    _cexp1           = vdwparam[tj+1];
    _cexp2           = vdwparam[tj+2];
#endif

    nj1              = 0;              
    
    for(n=0; (n<nri); n++)
    {
	double _Complex ix1,ix2,ix3,iy1,iy2,iy3,iz1,iz2,iz3,ix23,iy23,iz23;
	
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

        real _shX  = shiftvec[is3];  
        real _shY  = shiftvec[is3+1];
        real _shZ  = shiftvec[is3+2];

        int ii  = iinr[n];        
        int ii3 = 3*ii;
	int k,ggid;

	// add the shift vector to all water atoms

        real _ix1 = _shX + pos[ii3+0];
        real _iy1 = _shY + pos[ii3+1];
        real _iz1 = _shZ + pos[ii3+2];
        real _ix2 = _shX + pos[ii3+3];
        real _iy2 = _shY + pos[ii3+4];
        real _iz2 = _shZ + pos[ii3+5];
        real _ix3 = _shX + pos[ii3+6];
        real _iy3 = _shY + pos[ii3+7];
        real _iz3 = _shZ + pos[ii3+8];

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

        for(k=nj0; (k<nj1-1); k+=2)
        {
	    double _Complex jx1,jx2,jx3,fjx1,fjx2,fjx3,dx11,dx12,dx13,dx21,dx22,dx23,dx31,dx32,dx33,tx;
	    double _Complex jy1,jy2,jy3,fjy1,fjy2,fjy3,dy11,dy12,dy13,dy21,dy22,dy23,dy31,dy32,dy33,ty;
	    double _Complex jz1,jz2,jz3,fjz1,fjz2,fjz3,dz11,dz12,dz13,dz21,dz22,dz23,dz31,dz32,dz33,tz;
	    double _Complex rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
	    double _Complex rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
	    double _Complex rinvsq,krsq,gbscale,fscal;
	    double _Complex rinvsix,Vvdw6,Vvdw12,Vvdwexp;
	    double _Complex rt,r,rti,n0,eps,Y,F,G,H,GHeps,VV,FF,fijC,fijD;
	    
	    int nnn1,nnn2;
	    int jnr1,jnr2,j31,j32;

	    jnr1   = jjnr[k];
	    jnr2   = jjnr[k+1];

            j31    = 3*jnr1;        
            j32    = 3*jnr2;

	    // Coulomb and possibly VdW between i-water and Oxygens of j1- and j2-water

	    load6_and_merge(pos,j31,j32,jx1,jy1,jz1);

	    dx11   = __fpsub(ix1,jx1);
	    dy11   = __fpsub(iy1,jy1);
	    dz11   = __fpsub(iz1,jz1);

	    dx21   = __fpsub(ix2,jx1);
	    dy21   = __fpsub(iy2,jy1);
	    dz21   = __fpsub(iz2,jz1);

	    dx31   = __fpsub(ix3,jx1);
	    dy31   = __fpsub(iy3,jy1);
	    dz31   = __fpsub(iz3,jz1);

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

	    // only the oxygens will have a LJ-interaction

#ifndef NO_FORCE
	    load6_and_merge(faction,j31,j32,fjx1,fjy1,fjz1);
#endif

	    FULL_INTERACTION(qqOO_OO,rinv11,rsq11,conv1,jnr1,jnr2);

#ifndef NO_FORCE
            fix1   = __fpnmsub(fix1,dx11,fscal);
            fiy1   = __fpnmsub(fiy1,dy11,fscal);      
            fiz1   = __fpnmsub(fiz1,dz11,fscal);

            fjx1   = __fpmadd(fjx1,dx11,fscal);
            fjy1   = __fpmadd(fjy1,dy11,fscal);
            fjz1   = __fpmadd(fjz1,dz11,fscal);
#endif

	    COUL_INTERACTION(qqOH_OH,rinv21,rsq21,conv2,jnr1,jnr2);

#ifndef NO_FORCE
            fix2   = __fpnmsub(fix2,dx21,fscal);
            fiy2   = __fpnmsub(fiy2,dy21,fscal);      
            fiz2   = __fpnmsub(fiz2,dz21,fscal);

            fjx1   = __fpmadd(fjx1,dx21,fscal);
            fjy1   = __fpmadd(fjy1,dy21,fscal);
            fjz1   = __fpmadd(fjz1,dz21,fscal);
#endif

	    COUL_INTERACTION(qqOH_OH,rinv31,rsq31,conv3,jnr1,jnr2);

#ifndef NO_FORCE
            fix3   = __fpnmsub(fix3,dx31,fscal);
            fiy3   = __fpnmsub(fiy3,dy31,fscal);      
            fiz3   = __fpnmsub(fiz3,dz31,fscal);

            fjx1   = __fpmadd(fjx1,dx31,fscal);
            fjy1   = __fpmadd(fjy1,dy31,fscal);
            fjz1   = __fpmadd(fjz1,dz31,fscal);

	    split_and_store6(faction,j31,j32,fjx1,fjy1,fjz1);
#endif

	    // Coulomb between i-water and first hydrogen in j1- and j2-water

	    load6_and_merge(pos,j31+3,j32+3,jx2,jy2,jz2);

	    dx12   = __fpsub(ix1,jx2);
	    dy12   = __fpsub(iy1,jy2);
	    dz12   = __fpsub(iz1,jz2);
	    dx22   = __fpsub(ix2,jx2);
	    dy22   = __fpsub(iy2,jy2);
	    dz22   = __fpsub(iz2,jz2);
	    dx32   = __fpsub(ix3,jx2);
	    dy32   = __fpsub(iy3,jy2);
	    dz32   = __fpsub(iz3,jz2);

	    rsq12  = vdist2(dx12,dy12,dz12);
	    rsq22  = vdist2(dx22,dy22,dz22);
	    rsq32  = vdist2(dx32,dy32,dz32);

	    rinv12 = __fprsqrte(rsq12);
	    rinv22 = __fprsqrte(rsq22);
	    rinv32 = __fprsqrte(rsq32);

            rinv12 = sqrt_newton(rinv12,rsq12);
            rinv22 = sqrt_newton(rinv22,rsq22);
            rinv32 = sqrt_newton(rinv32,rsq32);

#ifdef GMX_DOUBLE

            rinv12 = sqrt_newton(rinv12,rsq12);
            rinv22 = sqrt_newton(rinv22,rsq22);
            rinv32 = sqrt_newton(rinv32,rsq32);

#endif

#ifndef NO_FORCE
	    load6_and_merge(faction,j31+3,j32+3,fjx2,fjy2,fjz2);
#endif

	    COUL_INTERACTION(qqOH_OH,rinv12,rsq12,conv1,jnr1,jnr2);

#ifndef NO_FORCE
            fix1   = __fpnmsub(fix1,dx12,fscal);
            fiy1   = __fpnmsub(fiy1,dy12,fscal);      
            fiz1   = __fpnmsub(fiz1,dz12,fscal);

            fjx2   = __fpmadd(fjx2,dx12,fscal);
            fjy2   = __fpmadd(fjy2,dy12,fscal);
            fjz2   = __fpmadd(fjz2,dz12,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv22,rsq22,conv2,jnr1,jnr2);

#ifndef NO_FORCE
            fix2   = __fpnmsub(fix2,dx22,fscal);
            fiy2   = __fpnmsub(fiy2,dy22,fscal);      
            fiz2   = __fpnmsub(fiz2,dz22,fscal);

            fjx2   = __fpmadd(fjx2,dx22,fscal);
            fjy2   = __fpmadd(fjy2,dy22,fscal);
            fjz2   = __fpmadd(fjz2,dz22,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv32,rsq32,conv3,jnr1,jnr2);

#ifndef NO_FORCE
            fix3   = __fpnmsub(fix3,dx32,fscal);
            fiy3   = __fpnmsub(fiy3,dy32,fscal);      
            fiz3   = __fpnmsub(fiz3,dz32,fscal);

            fjx2   = __fpmadd(fjx2,dx32,fscal);
            fjy2   = __fpmadd(fjy2,dy32,fscal);
            fjz2   = __fpmadd(fjz2,dz32,fscal);

	    split_and_store6(faction,j31+3,j32+3,fjx2,fjy2,fjz2);
#endif


	    // Coulomb between i-water and second hydrogen in j1- and j2-water

	    load6_and_merge(pos,j31+6,j32+6,jx3,jy3,jz3);

	    dx13   = __fpsub(ix1,jx3);
	    dy13   = __fpsub(iy1,jy3);
	    dz13   = __fpsub(iz1,jz3);
	    dx23   = __fpsub(ix2,jx3);
	    dy23   = __fpsub(iy2,jy3);
	    dz23   = __fpsub(iz2,jz3);
	    dx33   = __fpsub(ix3,jx3);
	    dy33   = __fpsub(iy3,jy3);
	    dz33   = __fpsub(iz3,jz3);

	    rsq13  = vdist2(dx13,dy13,dz13);
	    rsq23  = vdist2(dx23,dy23,dz23);
	    rsq33  = vdist2(dx33,dy33,dz33);

	    rinv13 = __fprsqrte(rsq13);
	    rinv23 = __fprsqrte(rsq23);
	    rinv33 = __fprsqrte(rsq33);

            rinv13 = sqrt_newton(rinv13,rsq13);
            rinv23 = sqrt_newton(rinv23,rsq23);
            rinv33 = sqrt_newton(rinv33,rsq33);

#ifdef GMX_DOUBLE

            rinv13 = sqrt_newton(rinv13,rsq13);
            rinv23 = sqrt_newton(rinv23,rsq23);
            rinv33 = sqrt_newton(rinv33,rsq33);

#endif

#ifndef NO_FORCE
	    load6_and_merge(faction,j31+6,j32+6,fjx3,fjy3,fjz3);
#endif

	    COUL_INTERACTION(qqOH_OH,rinv13,rsq13,conv1,jnr1,jnr2);

#ifndef NO_FORCE
            fix1   = __fpnmsub(fix1,dx13,fscal);
            fiy1   = __fpnmsub(fiy1,dy13,fscal);      
            fiz1   = __fpnmsub(fiz1,dz13,fscal);

            fjx3   = __fpmadd(fjx3,dx13,fscal);
            fjy3   = __fpmadd(fjy3,dy13,fscal);
            fjz3   = __fpmadd(fjz3,dz13,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv23,rsq23,conv2,jnr1,jnr2);

#ifndef NO_FORCE
            fix2   = __fpnmsub(fix2,dx23,fscal);
            fiy2   = __fpnmsub(fiy2,dy23,fscal);      
            fiz2   = __fpnmsub(fiz2,dz23,fscal);

            fjx3   = __fpmadd(fjx3,dx23,fscal);
            fjy3   = __fpmadd(fjy3,dy23,fscal);
            fjz3   = __fpmadd(fjz3,dz23,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv33,rsq33,conv3,jnr1,jnr2);

#ifndef NO_FORCE
            fix3   = __fpnmsub(fix3,dx33,fscal);
            fiy3   = __fpnmsub(fiy3,dy33,fscal);      
            fiz3   = __fpnmsub(fiz3,dz33,fscal);

            fjx3   = __fpmadd(fjx3,dx33,fscal);
            fjy3   = __fpmadd(fjy3,dy33,fscal);
            fjz3   = __fpmadd(fjz3,dz33,fscal);

	    split_and_store6(faction,j31+6,j32+6,fjx3,fjy3,fjz3);
#endif
        }

	ix23 = __cmplx(_ix2,_ix3);
	iy23 = __cmplx(_iy2,_iy3);
	iz23 = __cmplx(_iz2,_iz3);
        
        for(; (k<nj1); k++)
        {
	    double _Complex dx123,dx223,dx323,dx231,tx;
	    double _Complex dy123,dy223,dy323,dy231,ty;
	    double _Complex dz123,dz223,dz323,dz231,tz;
	    double _Complex jx23,jy23,jz23,jx1,jy1,jz1;
	    double _Complex fjx23,fjy23,fjz23,fjx1,fjy1,fjz1;
	    double _Complex rsq123,rsq223,rsq323,rsq231;
	    double _Complex rinv123,rinv223,rinv323,rinv231;
	    double _Complex rinvsq,krsq,fscal;
	    double _Complex rt,r,rti,eps,YF1,YF2,GH1,GH2,Y,F,G,H,GHeps,VV,FF,fijC,fijD,n0;
	    real _dx11,_dy11,_dz11,_rsq11,_rinv11;
	    real _rinvsq,_krsq,_fscal;
	    real _rt,_r,_eps,_Y,_F,_G,_H,_GHeps,_VV,_FF,_fijC,_fijD,_n0;
	    real _rinvsix,_Vvdw6,_Vvdw12,_Vvdwexp;

	    int nnn1,nnn2;
	    int jnr,j3;

	    jnr    = jjnr[k];
            j3     = 3*jnr;

	    load3_and_clone(pos,j3,jx1,jy1,jz1);

            _dx11   = ix1 - __creal(jx1);
            _dy11   = iy1 - __creal(jy1);
	    _dz11   = iz1 - __creal(jz1);
	    dx231   = __fpsub(ix23,jx1);
	    dy231   = __fpsub(iy23,jy1);
	    dz231   = __fpsub(iz23,jz1);

            _rsq11  = _dx11 * _dx11 + _dy11 * _dy11 + _dz11 * _dz11;
	    rsq231  = vdist2(dx231,dy231,dz231);
	    _rinv11 = __frsqrte(_rsq11);
	    rinv231 = __fprsqrte(rsq231);

	    _rinv11 = sqrt_newton_scalar(_rinv11,_rsq11);
            rinv231 = sqrt_newton(rinv231,rsq231);

#ifdef GMX_DOUBLE

	    _rinv11 = sqrt_newton_scalar(_rinv11,_rsq11);
            rinv231 = sqrt_newton(rinv231,rsq231);

#endif

	    FULL_INTERACTION_(_qqOO,_rinv11,_rsq11,jnr);

#ifndef NO_FORCE

	    load3_real(faction,j3,fjx1,fjy1,fjz1);

            fix1            -= _dx11 * _fscal;
            fiy1            -= _dy11 * _fscal;
            fiz1            -= _dz11 * _fscal;

            fjx1            += _dx11 * _fscal;
            fjy1            += _dy11 * _fscal;
            fjz1            += _dz11 * _fscal;
#endif

	    COUL_INTERACTION(qqOH_OH,rinv231,rsq231,conv4,jnr,jnr);

#ifndef NO_FORCE
	    tx     = __fpmul(dx231,fscal);
	    ty     = __fpmul(dy231,fscal);
	    tz     = __fpmul(dz231,fscal);

            fix2    = __fpnmsub(fix2,rone,tx);      
            fiy2    = __fpnmsub(fiy2,rone,ty); 
            fiz2    = __fpnmsub(fiz2,rone,tz);
            fix3    = __fxnmsub(fix3,rone,tx);      
            fiy3    = __fxnmsub(fiy3,rone,ty); 
            fiz3    = __fxnmsub(fiz3,rone,tz);

            fjx1   = __fpadd(fjx1,tx);
            fjy1   = __fpadd(fjy1,ty);
            fjz1   = __fpadd(fjz1,tz);

	    sum_and_store3(faction,j3,fjx1,fjy1,fjz1);
#endif

	    load6_and_merge(pos,j3+3,j3+6,jx23,jy23,jz23);

	    dx123   = __fpsub(ix1,jx23);
	    dy123   = __fpsub(iy1,jy23);
	    dz123   = __fpsub(iz1,jz23);
	    dx223   = __fpsub(ix2,jx23);
	    dy223   = __fpsub(iy2,jy23);
	    dz223   = __fpsub(iz2,jz23);
	    dx323   = __fpsub(ix3,jx23);
	    dy323   = __fpsub(iy3,jy23);
	    dz323   = __fpsub(iz3,jz23);

	    rsq123  = vdist2(dx123,dy123,dz123);
	    rsq223  = vdist2(dx223,dy223,dz223);
	    rsq323  = vdist2(dx323,dy323,dz323);

	    rinv123 = __fprsqrte(rsq123);
	    rinv223 = __fprsqrte(rsq223);
	    rinv323 = __fprsqrte(rsq323);

	    rinv123 = sqrt_newton(rinv123,rsq123);
	    rinv223 = sqrt_newton(rinv223,rsq223);
	    rinv323 = sqrt_newton(rinv323,rsq323);

#ifdef GMX_DOUBLE

	    rinv123 = sqrt_newton(rinv123,rsq123);
	    rinv223 = sqrt_newton(rinv223,rsq223);
	    rinv323 = sqrt_newton(rinv323,rsq323);

#endif
	    COUL_INTERACTION(qqOH_OH,rinv123,rsq123,conv1,jnr,jnr);

#ifndef NO_FORCE
	    load6_and_merge(faction,j3+3,j3+6,fjx23,fjy23,fjz23);

            fix1   = __fpnmsub(fix1,dx123,fscal);
            fiy1   = __fpnmsub(fiy1,dy123,fscal);      
            fiz1   = __fpnmsub(fiz1,dz123,fscal);

            fjx23  = __fpmadd(fjx23,dx123,fscal);
            fjy23  = __fpmadd(fjy23,dy123,fscal);
            fjz23  = __fpmadd(fjz23,dz123,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv223,rsq223,conv2,jnr,jnr);

#ifndef NO_FORCE
            fix2   = __fpnmsub(fix2,dx223,fscal);
            fiy2   = __fpnmsub(fiy2,dy223,fscal);      
            fiz2   = __fpnmsub(fiz2,dz223,fscal);

            fjx23  = __fpmadd(fjx23,dx223,fscal);
            fjy23  = __fpmadd(fjy23,dy223,fscal);
            fjz23  = __fpmadd(fjz23,dz223,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv323,rsq323,conv3,jnr,jnr);

#ifndef NO_FORCE
            fix3   = __fpnmsub(fix3,dx323,fscal);
            fiy3   = __fpnmsub(fiy3,dy323,fscal);      
            fiz3   = __fpnmsub(fiz3,dz323,fscal);

            fjx23  = __fpmadd(fjx23,dx323,fscal);
            fjy23  = __fpmadd(fjy23,dy323,fscal);
            fjz23  = __fpmadd(fjz23,dz323,fscal);

	    split_and_store6(faction,j3+3,j3+6,fjx23,fjy23,fjz23);
#endif
	}

	ggid = gid[n];

#ifndef NO_FORCE
        sum_and_add3(faction,ii3  ,fix1,fiy1,fiz1);
        sum_and_add3(faction,ii3+3,fix2,fiy2,fiz2);
        sum_and_add3(faction,ii3+6,fix3,fiy3,fiz3);

	fshift[is3]      = fshift[is3]
	                 + (__creal(fix1) + __cimag(fix1)) 
	                 + (__creal(fix2) + __cimag(fix2))
	                 + (__creal(fix3) + __cimag(fix3));
        fshift[is3+1]    = fshift[is3+1]
	                 + (__creal(fiy1) + __cimag(fiy1)) 
	                 + (__creal(fiy2) + __cimag(fiy2))
	                 + (__creal(fiy3) + __cimag(fiy3));
        fshift[is3+2]    = fshift[is3+2]
	                 + (__creal(fiz1) + __cimag(fiz1)) 
	                 + (__creal(fiz2) + __cimag(fiz2))
	                 + (__creal(fiz3) + __cimag(fiz3));
#endif

#if COULOMB != COULOMB_NONE
        Vc[ggid]         += __creal(vctot)   + __cimag(vctot);
#endif

#if VDW != VDW_NONE
        Vvdw[ggid]       += __creal(Vvdwtot) + __cimag(Vvdwtot);
#endif
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}
