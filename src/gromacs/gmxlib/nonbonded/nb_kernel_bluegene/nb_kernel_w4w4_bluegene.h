#include "types/simple.h"

#undef FULL_INTERACTION
#undef FULL_INTERACTION_
#undef COUL_INTERACTION
#undef COUL_INTERACTION_
#undef VDW_INTERACTION
#undef VDW_INTERACTION_

#ifdef NO_FORCE

 #define FULL_INTERACTION  nfcalc_interaction
 #define FULL_INTERACTION_ nfcalc_interaction_
 #define COUL_INTERACTION  nfcalc_coul_interaction
 #define COUL_INTERACTION_ nfcalc_coul_interaction_
 #define VDW_INTERACTION   nfcalc_vdw_interaction
 #define VDW_INTERACTION_  nfcalc_vdw_interaction_

#else

 #define FULL_INTERACTION  calc_interaction
 #define FULL_INTERACTION_ calc_interaction_
 #define COUL_INTERACTION  calc_coul_interaction
 #define COUL_INTERACTION_ calc_coul_interaction_
 #define VDW_INTERACTION   calc_vdw_interaction
 #define VDW_INTERACTION_  calc_vdw_interaction_

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
    double _Complex qqMM_MM,qqMM_MH,qqMH_MH,qqMH_HH,qqHH_HH;
    double _Complex zero  = __cmplx(0.0,0.0);
    double _Complex one   = __cmplx(1.0,1.0);
    double _Complex rone  = __cmplx(1.0,0.0);
    double _Complex lr    = __cmplx(1.0,-1.0);
    double _Complex rl    = __cmplx(-1.0,1.0);
    double _Complex two   = __cmplx(2.0,2.0);
    double _Complex three = __cmplx(3.0,3.0);

    double conv1[2],conv2[2],conv3[2],conv4[2];

    real  _qM,_qH,_qqMM,_qqMH,_qqHH;
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

    _qH              = charge[ii+1];   
    _qM              = charge[ii+3];

    _qqMM            = _facel*_qM*_qM;    
    _qqMH            = _facel*_qM*_qH;
    _qqHH            = _facel*_qH*_qH;    
    qqMM_MM          = __cmplx(_qqMM,_qqMM);
    qqMH_MH          = __cmplx(_qqMH,_qqMH);
    qqMM_MH          = __cmplx(_qqMM,_qqMH);
    qqMH_HH          = __cmplx(_qqMH,_qqHH);
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
	double _Complex ix1,ix2,ix3,ix4,iy1,iy2,iy3,iy4,iz1,iz2,iz3,iz4,ix14,iy14,iz14,ix23,iy23,iz23;
	
	// initialize potential energies and forces for this paricle

        double _Complex vctot     = zero;
        double _Complex Vvdwtot   = zero;
	double _Complex fix1      = zero;
	double _Complex fix2      = zero;
	double _Complex fix3      = zero;
	double _Complex fix4      = zero;
	double _Complex fiy1      = zero;
	double _Complex fiy2      = zero;
	double _Complex fiy3      = zero;
	double _Complex fiy4      = zero;
	double _Complex fiz1      = zero;
	double _Complex fiz2      = zero;
	double _Complex fiz3      = zero;
	double _Complex fiz4      = zero;

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
        real _ix4 = _shX + pos[ii3+9];
        real _iy4 = _shY + pos[ii3+10];
        real _iz4 = _shZ + pos[ii3+11];

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

	ix4 = __cmplx(_ix4,_ix4);
	iy4 = __cmplx(_iy4,_iy4);
	iz4 = __cmplx(_iz4,_iz4);

        nj0 = jindex[n];      
        nj1 = jindex[n+1];    

        for(k=nj0; (k<nj1-1); k+=2)
        {
	    double _Complex jx1,jx2,jx3,jx4,fjx1,fjx2,fjx3,fjx4,tx;
	    double _Complex dx11,dx12,dx13,dx22,dx23,dx24,dx32,dx33,dx34,dx42,dx43,dx44;
	    double _Complex jy1,jy2,jy3,jy4,fjy1,fjy2,fjy3,fjy4,ty;
	    double _Complex dy11,dy12,dy13,dy22,dy23,dy24,dy32,dy33,dy34,dy42,dy43,dy44;
	    double _Complex jz1,jz2,jz3,jz4,fjz1,fjz2,fjz3,fjz4,tz;
	    double _Complex dz11,dz12,dz13,dz22,dz23,dz24,dz32,dz33,dz34,dz42,dz43,dz44;
	    double _Complex rsq11,rsq12,rsq13,rsq22,rsq23,rsq24,rsq32,rsq33,rsq34,rsq42,rsq43,rsq44;
	    double _Complex rinv11,rinv12,rinv13,rinv22,rinv23,rinv24,rinv32,rinv33,rinv34,rinv42,rinv43,rinv44;
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

	    load6_and_merge(pos,j31  ,j32  ,jx1,jy1,jz1);
	    load6_and_merge(pos,j31+9,j32+9,jx4,jy4,jz4);
	    
	    /*
	    a1 = __lfpd(&pos[j31]);
	    a2 = __lfxd(&pos[j32]);
	    b1 = __lfpd(&pos[j31+2]);
	    b2 = __lfxd(&pos[j32+2]);
	    c1 = __lfpd(&pos[j31+4]);
	    c2 = __lfxd(&pos[j32+4]);
	    
	    jx1 = __fpsel(a1,a2,rone);
	    jx2 = __fpsel(a1,a2,cone);
	    jy1 = __fpsel(b1,b2,rone);
	    jy2 = __fpsel(b1,b2,cone);
	    jz1 = __fpsel(c1,c2,rone);
	    jz2 = __fpsel(c1,c2,cone);
	    */

	    dx11   = __fpsub(ix1,jx1);
	    dy11   = __fpsub(iy1,jy1);
	    dz11   = __fpsub(iz1,jz1);

	    dx24   = __fpsub(ix2,jx4);
	    dy24   = __fpsub(iy2,jy4);
	    dz24   = __fpsub(iz2,jz4);

	    dx34   = __fpsub(ix3,jx4);
	    dy34   = __fpsub(iy3,jy4);
	    dz34   = __fpsub(iz3,jz4);

	    dx44   = __fpsub(ix4,jx4);
	    dy44   = __fpsub(iy4,jy4);
	    dz44   = __fpsub(iz4,jz4);

	    rsq11  = vdist2(dx11,dy11,dz11);
	    rsq24  = vdist2(dx24,dy24,dz24);
	    rsq34  = vdist2(dx34,dy34,dz34);
	    rsq44  = vdist2(dx44,dy44,dz44);

	    rinv11 = __fprsqrte(rsq11);
	    rinv24 = __fprsqrte(rsq24);
	    rinv34 = __fprsqrte(rsq34);
	    rinv44 = __fprsqrte(rsq44);

	    rinv11 = sqrt_newton(rinv11,rsq11);
	    rinv24 = sqrt_newton(rinv24,rsq24);
            rinv34 = sqrt_newton(rinv34,rsq34);
            rinv44 = sqrt_newton(rinv44,rsq44);

#ifdef GMX_DOUBLE

	    rinv11 = sqrt_newton(rinv11,rsq11);
	    rinv24 = sqrt_newton(rinv24,rsq24);
            rinv34 = sqrt_newton(rinv34,rsq34);
            rinv44 = sqrt_newton(rinv44,rsq44);

#endif

	    // only the oxygens will have a VDW-interaction

#ifndef NO_FORCE
	    load6_and_merge(faction,j31,j32,fjx1,fjy1,fjz1);
#endif

	    VDW_INTERACTION(rinv11,rsq11,conv1,jnr1,jnr2);

#ifndef NO_FORCE
            fix1   = __fpnmsub(fix1,dx11,fscal);
            fiy1   = __fpnmsub(fiy1,dy11,fscal);      
            fiz1   = __fpnmsub(fiz1,dz11,fscal);

            fjx1   = __fpmadd(fjx1,dx11,fscal);
            fjy1   = __fpmadd(fjy1,dy11,fscal);
            fjz1   = __fpmadd(fjz1,dz11,fscal);
#endif

	    COUL_INTERACTION(qqMM_MM,rinv44,rsq44,conv4,jnr1,jnr2);

#ifndef NO_FORCE
	    load6_and_merge (faction,j31+9,j32+9,fjx4,fjy4,fjz4);
	    split_and_store6(faction,j31  ,j32  ,fjx1,fjy1,fjz1);

            fix4   = __fpnmsub(fix4,dx44,fscal);
            fiy4   = __fpnmsub(fiy4,dy44,fscal);      
            fiz4   = __fpnmsub(fiz4,dz44,fscal);

            fjx4   = __fpmadd(fjx4,dx44,fscal);
            fjy4   = __fpmadd(fjy4,dy44,fscal);
            fjz4   = __fpmadd(fjz4,dz44,fscal);
#endif

	    COUL_INTERACTION(qqMH_MH,rinv24,rsq24,conv2,jnr1,jnr2);

#ifndef NO_FORCE
            fix2   = __fpnmsub(fix2,dx24,fscal);
            fiy2   = __fpnmsub(fiy2,dy24,fscal);      
            fiz2   = __fpnmsub(fiz2,dz24,fscal);

            fjx4   = __fpmadd(fjx4,dx24,fscal);
            fjy4   = __fpmadd(fjy4,dy24,fscal);
            fjz4   = __fpmadd(fjz4,dz24,fscal);
#endif

	    COUL_INTERACTION(qqMH_MH,rinv34,rsq34,conv3,jnr1,jnr2);

#ifndef NO_FORCE
            fix3   = __fpnmsub(fix3,dx34,fscal);
            fiy3   = __fpnmsub(fiy3,dy34,fscal);      
            fiz3   = __fpnmsub(fiz3,dz34,fscal);

            fjx4   = __fpmadd(fjx4,dx34,fscal);
            fjy4   = __fpmadd(fjy4,dy34,fscal);
            fjz4   = __fpmadd(fjz4,dz34,fscal);

	    split_and_store6(faction,j31+9,j32+9,fjx4,fjy4,fjz4);
#endif
	    
	    // Coulomb between i-water and first hydrogen in j1- and j2-water

	    load6_and_merge(pos,j31+3,j32+3,jx2,jy2,jz2);

	    dx22   = __fpsub(ix2,jx2);
	    dy22   = __fpsub(iy2,jy2);
	    dz22   = __fpsub(iz2,jz2);
	    dx32   = __fpsub(ix3,jx2);
	    dy32   = __fpsub(iy3,jy2);
	    dz32   = __fpsub(iz3,jz2);
	    dx42   = __fpsub(ix4,jx2);
	    dy42   = __fpsub(iy4,jy2);
	    dz42   = __fpsub(iz4,jz2);

	    rsq22  = vdist2(dx22,dy22,dz22);
	    rsq32  = vdist2(dx32,dy32,dz32);
	    rsq42  = vdist2(dx42,dy42,dz42);

	    rinv22 = __fprsqrte(rsq22);
	    rinv32 = __fprsqrte(rsq32);
	    rinv42 = __fprsqrte(rsq42);

            rinv22 = sqrt_newton(rinv22,rsq22);
            rinv32 = sqrt_newton(rinv32,rsq32);
            rinv42 = sqrt_newton(rinv42,rsq42);

#ifdef GMX_DOUBLE

            rinv22 = sqrt_newton(rinv22,rsq22);
            rinv32 = sqrt_newton(rinv32,rsq32);
            rinv42 = sqrt_newton(rinv42,rsq42);

#endif

#ifndef NO_FORCE
	    load6_and_merge(faction,j31+3,j32+3,fjx2,fjy2,fjz2);
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
#endif

	    COUL_INTERACTION(qqMH_MH,rinv42,rsq42,conv4,jnr1,jnr2);

#ifndef NO_FORCE
            fix4   = __fpnmsub(fix4,dx42,fscal);
            fiy4   = __fpnmsub(fiy4,dy42,fscal);      
            fiz4   = __fpnmsub(fiz4,dz42,fscal);

            fjx2   = __fpmadd(fjx2,dx42,fscal);
            fjy2   = __fpmadd(fjy2,dy42,fscal);
            fjz2   = __fpmadd(fjz2,dz42,fscal);

	    split_and_store6(faction,j31+3,j32+3,fjx2,fjy2,fjz2);
#endif

	    // Coulomb between i-water and second hydrogen in j1- and j2-water

	    load6_and_merge(pos,j31+6,j32+6,jx3,jy3,jz3);

	    dx43   = __fpsub(ix4,jx3);
	    dy43   = __fpsub(iy4,jy3);
	    dz43   = __fpsub(iz4,jz3);

	    dx23   = __fpsub(ix2,jx3);
	    dy23   = __fpsub(iy2,jy3);
	    dz23   = __fpsub(iz2,jz3);

	    dx33   = __fpsub(ix3,jx3);
	    dy33   = __fpsub(iy3,jy3);
	    dz33   = __fpsub(iz3,jz3);

	    rsq43  = vdist2(dx43,dy43,dz43);
	    rsq23  = vdist2(dx23,dy23,dz23);
	    rsq33  = vdist2(dx33,dy33,dz33);

	    rinv43 = __fprsqrte(rsq43);
	    rinv23 = __fprsqrte(rsq23);
	    rinv33 = __fprsqrte(rsq33);

            rinv43 = sqrt_newton(rinv43,rsq43);
            rinv23 = sqrt_newton(rinv23,rsq23);
            rinv33 = sqrt_newton(rinv33,rsq33);

#ifdef GMX_DOUBLE

            rinv43 = sqrt_newton(rinv43,rsq43);
            rinv23 = sqrt_newton(rinv23,rsq23);
            rinv33 = sqrt_newton(rinv33,rsq33);

#endif

#ifndef NO_FORCE
	    load6_and_merge(faction,j31+6,j32+6,fjx3,fjy3,fjz3);
#endif

	    COUL_INTERACTION(qqMH_MH,rinv43,rsq43,conv1,jnr1,jnr2);

#ifndef NO_FORCE
            fix4   = __fpnmsub(fix4,dx43,fscal);
            fiy4   = __fpnmsub(fiy4,dy43,fscal);      
            fiz4   = __fpnmsub(fiz4,dz43,fscal);

            fjx3   = __fpmadd(fjx3,dx43,fscal);
            fjy3   = __fpmadd(fjy3,dy43,fscal);
            fjz3   = __fpmadd(fjz3,dz43,fscal);
#endif

	    COUL_INTERACTION(qqHH_HH,rinv23,rsq23,conv2,jnr1,jnr2);

#ifndef NO_FORCE
	    split_and_store6(faction,j31+9,j32+9,fjx4,fjy4,fjz4);

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

	ix14 = __cmplx(_ix1,_ix4);
	iy14 = __cmplx(_iy1,_iy4);
	iz14 = __cmplx(_iz1,_iz4);
	ix23 = __cmplx(_ix2,_ix3);
	iy23 = __cmplx(_iy2,_iy3);
	iz23 = __cmplx(_iz2,_iz3);
        
        for(; (k<nj1); k++)
        {
	    double _Complex dx14,dx423,dx223,dx323,dx234,tx;
	    double _Complex dy14,dy423,dy223,dy323,dy234,ty;
	    double _Complex dz14,dz423,dz223,dz323,dz234,tz;
	    double _Complex jx23,jy23,jz23,jx14,jy14,jz14;
	    double _Complex fjx23,fjy23,fjz23,fjx14,fjy14,fjz14;
	    double _Complex rsq14,rsq423,rsq223,rsq323,rsq234;
	    double _Complex rinv14,rinv423,rinv223,rinv323,rinv234;
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

	    load6_and_merge(pos,j3,j3+9,jx14,jy14,jz14);

	    dx14    = __fpsub(ix14,jx14);
	    dy14    = __fpsub(iy14,jy14);
	    dz14    = __fpsub(iz14,jz14);

	    dx234   = __fxcpnmsub(ix23,one,__cimag(jx14));
	    dy234   = __fxcpnmsub(iy23,one,__cimag(jy14));
	    dz234   = __fxcpnmsub(iz23,one,__cimag(jz14));

	    rsq14   = vdist2(dx14,dy14,dz14);
	    rsq234  = vdist2(dx234,dy234,dz234);

	    rinv14  = __fprsqrte(rsq14);
	    rinv234 = __fprsqrte(rsq234);

	    rinv14  = sqrt_newton(rinv14,rsq14);
            rinv234 = sqrt_newton(rinv234,rsq234);

#ifdef GMX_DOUBLE

	    rinv14  = sqrt_newton(rinv14,rsq14);
            rinv234 = sqrt_newton(rinv234,rsq234);

#endif

	    VDW_INTERACTION_(__creal(rinv14),__creal(rsq14),jnr);
	    fscal = _fscal;
	    COUL_INTERACTION_(_qqMM,__cimag(rinv14),__cimag(rsq14),jnr);
	    fscal = __cmplx(__creal(fscal),_fscal);

#ifndef NO_FORCE

	    load6_and_merge(faction,j3,j3+9,fjx14,fjy14,fjz14);

	    tx = __fpmul(dx14,fscal);
	    ty = __fpmul(dy14,fscal);
	    tz = __fpmul(dz14,fscal);

            fix1   -= __creal(tx);
            fiy1   -= __creal(ty);
            fiz1   -= __creal(tz);
            fix4   -= __cimag(tx);
            fiy4   -= __cimag(ty);
            fiz4   -= __cimag(tz);

            fjx14   = __fpadd(fjx14,tx);
            fjy14   = __fpadd(fjy14,ty);
            fjz14   = __fpadd(fjz14,tz);

#endif

	    COUL_INTERACTION(qqMH_MH,rinv234,rsq234,conv4,jnr,jnr);

	    load6_and_merge(pos,j3+3,j3+6,jx23,jy23,jz23);

#ifndef NO_FORCE
	    tx     = __fpmul(dx234,fscal);
	    ty     = __fpmul(dy234,fscal);
	    tz     = __fpmul(dz234,fscal);

            fix2    = __fpnmsub(fix2,rone,tx);      
            fiy2    = __fpnmsub(fiy2,rone,ty); 
            fiz2    = __fpnmsub(fiz2,rone,tz);
            fix3    = __fxnmsub(fix3,rone,tx);      
            fiy3    = __fxnmsub(fiy3,rone,ty); 
            fiz3    = __fxnmsub(fiz3,rone,tz);

	    tx      = __fxcxma(tx,tx,rone);
	    ty      = __fxcxma(ty,ty,rone);
	    tz      = __fxcxma(tz,tz,rone);

	    fjx14  = __fxmadd(fjx14,rone,tx);
	    fjy14  = __fxmadd(fjy14,rone,ty);
	    fjz14  = __fxmadd(fjz14,rone,tz);

	    split_and_store6(faction,j3,j3+9,fjx14,fjy14,fjz14);
#endif

	    dx423   = __fpsub(ix4,jx23);
	    dy423   = __fpsub(iy4,jy23);
	    dz423   = __fpsub(iz4,jz23);

	    dx223   = __fpsub(ix2,jx23);
	    dy223   = __fpsub(iy2,jy23);
	    dz223   = __fpsub(iz2,jz23);

	    dx323   = __fpsub(ix3,jx23);
	    dy323   = __fpsub(iy3,jy23);
	    dz323   = __fpsub(iz3,jz23);

	    rsq423  = vdist2(dx423,dy423,dz423);
	    rsq223  = vdist2(dx223,dy223,dz223);
	    rsq323  = vdist2(dx323,dy323,dz323);

	    rinv423 = __fprsqrte(rsq423);
	    rinv223 = __fprsqrte(rsq223);
	    rinv323 = __fprsqrte(rsq323);

	    rinv423 = sqrt_newton(rinv423,rsq423);
	    rinv223 = sqrt_newton(rinv223,rsq223);
	    rinv323 = sqrt_newton(rinv323,rsq323);

#ifdef GMX_DOUBLE

	    rinv423 = sqrt_newton(rinv423,rsq423);
	    rinv223 = sqrt_newton(rinv223,rsq223);
	    rinv323 = sqrt_newton(rinv323,rsq323);

#endif
	    COUL_INTERACTION(qqMH_MH,rinv423,rsq423,conv1,jnr,jnr);

#ifndef NO_FORCE
	    load6_and_merge(faction,j3+3,j3+6,fjx23,fjy23,fjz23);

            fix4   = __fpnmsub(fix4,dx423,fscal);
            fiy4   = __fpnmsub(fiy4,dy423,fscal);      
            fiz4   = __fpnmsub(fiz4,dz423,fscal);

            fjx23  = __fpmadd(fjx23,dx423,fscal);
            fjy23  = __fpmadd(fjy23,dy423,fscal);
            fjz23  = __fpmadd(fjz23,dz423,fscal);
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
        sum_and_add3(faction,ii3+9,fix4,fiy4,fiz4);

	fshift[is3]      = fshift[is3]
	                 + (__creal(fix1) + __cimag(fix1)) 
	                 + (__creal(fix2) + __cimag(fix2))
	                 + (__creal(fix3) + __cimag(fix3))
	                 + (__creal(fix4) + __cimag(fix4));
        fshift[is3+1]    = fshift[is3+1]
	                 + (__creal(fiy1) + __cimag(fiy1)) 
	                 + (__creal(fiy2) + __cimag(fiy2))
	                 + (__creal(fiy3) + __cimag(fiy3))
	                 + (__creal(fiy4) + __cimag(fiy4));
        fshift[is3+2]    = fshift[is3+2]
	                 + (__creal(fiz1) + __cimag(fiz1)) 
	                 + (__creal(fiz2) + __cimag(fiz2))
	                 + (__creal(fiz3) + __cimag(fiz3))
	                 + (__creal(fiz4) + __cimag(fiz4));
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
