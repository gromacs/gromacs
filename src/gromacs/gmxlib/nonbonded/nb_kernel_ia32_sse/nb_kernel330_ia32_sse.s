##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##


.globl nb_kernel330_ia32_sse
.globl _nb_kernel330_ia32_sse
nb_kernel330_ia32_sse:  
_nb_kernel330_ia32_sse: 
.set nb330_p_nri, 8
.set nb330_iinr, 12
.set nb330_jindex, 16
.set nb330_jjnr, 20
.set nb330_shift, 24
.set nb330_shiftvec, 28
.set nb330_fshift, 32
.set nb330_gid, 36
.set nb330_pos, 40
.set nb330_faction, 44
.set nb330_charge, 48
.set nb330_p_facel, 52
.set nb330_argkrf, 56
.set nb330_argcrf, 60
.set nb330_Vc, 64
.set nb330_type, 68
.set nb330_p_ntype, 72
.set nb330_vdwparam, 76
.set nb330_Vvdw, 80
.set nb330_p_tabscale, 84
.set nb330_VFtab, 88
.set nb330_invsqrta, 92
.set nb330_dvda, 96
.set nb330_p_gbtabscale, 100
.set nb330_GBtab, 104
.set nb330_p_nthreads, 108
.set nb330_count, 112
.set nb330_mtx, 116
.set nb330_outeriter, 120
.set nb330_inneriter, 124
.set nb330_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330_ix, 0
.set nb330_iy, 16
.set nb330_iz, 32
.set nb330_iq, 48
.set nb330_dx, 64
.set nb330_dy, 80
.set nb330_dz, 96
.set nb330_two, 112
.set nb330_tsc, 128
.set nb330_qq, 144
.set nb330_c6, 160
.set nb330_c12, 176
.set nb330_fscal, 192
.set nb330_vctot, 208
.set nb330_Vvdwtot, 224
.set nb330_fix, 240
.set nb330_fiy, 256
.set nb330_fiz, 272
.set nb330_half, 288
.set nb330_three, 304
.set nb330_is3, 320
.set nb330_ii3, 324
.set nb330_ntia, 328
.set nb330_innerjjnr, 332
.set nb330_innerk, 336
.set nb330_n, 340
.set nb330_nn1, 344
.set nb330_nri, 348
.set nb330_facel, 352
.set nb330_ntype, 356
.set nb330_nouter, 360
.set nb330_ninner, 364
.set nb330_salign, 368
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $372,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb330_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb330_p_nri(%ebp),%ecx
        movl nb330_p_facel(%ebp),%esi
        movl nb330_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb330_nri(%esp)
        movl %esi,nb330_facel(%esp)
        movl %edi,nb330_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330_nouter(%esp)
        movl %eax,nb330_ninner(%esp)


        movl nb330_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb330_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb330_half(%esp)
        movss nb330_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb330_half(%esp)
        movaps %xmm2,nb330_two(%esp)
        movaps %xmm3,nb330_three(%esp)

_nb_kernel330_ia32_sse.nb330_threadloop: 
        movl  nb330_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330_ia32_sse.nb330_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330_ia32_sse.nb330_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330_n(%esp)
        movl %ebx,nb330_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330_ia32_sse.nb330_outerstart
        jmp _nb_kernel330_ia32_sse.nb330_end

_nb_kernel330_ia32_sse.nb330_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330_nouter(%esp),%ebx
        movl %ebx,nb330_nouter(%esp)

_nb_kernel330_ia32_sse.nb330_outer: 
        movl  nb330_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb330_is3(%esp)      ## store is3 

        movl  nb330_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb330_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb330_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb330_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb330_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb330_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb330_ix(%esp)
        movaps %xmm1,nb330_iy(%esp)
        movaps %xmm2,nb330_iz(%esp)

        movl  %ebx,nb330_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb330_vctot(%esp)
        movaps %xmm4,nb330_Vvdwtot(%esp)
        movaps %xmm4,nb330_fix(%esp)
        movaps %xmm4,nb330_fiy(%esp)
        movaps %xmm4,nb330_fiz(%esp)

        movl  nb330_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330_pos(%ebp),%esi
        movl  nb330_faction(%ebp),%edi
        movl  nb330_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb330_ninner(%esp),%ecx
        movl  %ecx,nb330_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb330_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel330_ia32_sse.nb330_unroll_loop
        jmp   _nb_kernel330_ia32_sse.nb330_finish_inner
_nb_kernel330_ia32_sse.nb330_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb330_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb330_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb330_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb330_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb330_qq(%esp)

        movl nb330_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb330_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb330_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb330_c6(%esp)
        movaps %xmm6,nb330_c12(%esp)

        movl nb330_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb330_ix(%esp),%xmm4
        movaps nb330_iy(%esp),%xmm5
        movaps nb330_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb330_dx(%esp)
        movaps %xmm5,nb330_dy(%esp)
        movaps %xmm6,nb330_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330_tsc(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb330_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        movaps nb330_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330_vctot(%esp),%xmm5
        movaps %xmm5,nb330_vctot(%esp)

        ## put scalar force on stack temporarily 
        movaps %xmm3,nb330_fscal(%esp)

        ## dispersion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        addps  nb330_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb330_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb330_fscal(%esp)
        movaps %xmm5,nb330_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330_fscal(%esp),%xmm7

        addps  nb330_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb330_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb330_dx(%esp),%xmm0
        movaps nb330_dy(%esp),%xmm1
        movaps nb330_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb330_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb330_fix(%esp),%xmm3
        movaps nb330_fiy(%esp),%xmm4
        movaps nb330_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb330_fix(%esp)
        movaps %xmm4,nb330_fiy(%esp)
        movaps %xmm5,nb330_fiz(%esp)
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%edi,%eax,4),%xmm4
        movlps (%edi,%ecx,4),%xmm6
        movhps (%edi,%ebx,4),%xmm4
        movhps (%edi,%edx,4),%xmm6

        movaps %xmm4,%xmm3
        shufps $136,%xmm6,%xmm3 ## constant 10001000
        shufps $221,%xmm6,%xmm4 ## constant 11011101                          

        ## now xmm3-xmm5 contains fjx, fjy, fjz 
        subps  %xmm0,%xmm3
        subps  %xmm1,%xmm4

        ## unpack them back so we can store them - first x & y in xmm3/xmm4 

        movaps %xmm3,%xmm6
        unpcklps %xmm4,%xmm6
        unpckhps %xmm4,%xmm3
        ## xmm6(l)=x & y for j1, (h) for j2 
        ## xmm3(l)=x & y for j3, (h) for j4 
        movlps %xmm6,(%edi,%eax,4)
        movlps %xmm3,(%edi,%ecx,4)

        movhps %xmm6,(%edi,%ebx,4)
        movhps %xmm3,(%edi,%edx,4)

        ## and the z forces 
        movss  8(%edi,%eax,4),%xmm4
        movss  8(%edi,%ebx,4),%xmm5
        movss  8(%edi,%ecx,4),%xmm6
        movss  8(%edi,%edx,4),%xmm7
        subss  %xmm2,%xmm4
        shufps $229,%xmm2,%xmm2 ## constant 11100101
        subss  %xmm2,%xmm5
        shufps $234,%xmm2,%xmm2 ## constant 11101010
        subss  %xmm2,%xmm6
        shufps $255,%xmm2,%xmm2 ## constant 11111111
        subss  %xmm2,%xmm7
        movss  %xmm4,8(%edi,%eax,4)
        movss  %xmm5,8(%edi,%ebx,4)
        movss  %xmm6,8(%edi,%ecx,4)
        movss  %xmm7,8(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb330_innerk(%esp)
        jl    _nb_kernel330_ia32_sse.nb330_finish_inner
        jmp   _nb_kernel330_ia32_sse.nb330_unroll_loop
_nb_kernel330_ia32_sse.nb330_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb330_innerk(%esp)
        movl  nb330_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel330_ia32_sse.nb330_dopair
        jmp   _nb_kernel330_ia32_sse.nb330_checksingle
_nb_kernel330_ia32_sse.nb330_dopair: 
        movl nb330_charge(%ebp),%esi

    movl  nb330_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb330_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb330_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb330_qq(%esp)

        movl nb330_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb330_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb330_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb330_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb330_c6(%esp)
        movaps %xmm6,nb330_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb330_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb330_ix(%esp),%xmm4
        movaps nb330_iy(%esp),%xmm5
        movaps nb330_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb330_dx(%esp)
        movaps %xmm5,nb330_dy(%esp)
        movaps %xmm6,nb330_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb330_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%ecx,4),%xmm5
        movhps (%esi,%edx,4),%xmm5 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ecx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## constant 10001000
        shufps $221,%xmm7,%xmm7 ## constant 11011101
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        movaps nb330_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330_vctot(%esp),%xmm5
        movaps %xmm5,nb330_vctot(%esp)

        ## put scalar force on stack temporarily 
        movaps %xmm3,nb330_fscal(%esp)

        ## dispersion 
        movlps 16(%esi,%ecx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm5,%xmm5 ## constant 11011101

        movlps 24(%esi,%ecx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm7    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## constant 10001000
        shufps $221,%xmm7,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        addps  nb330_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb330_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb330_fscal(%esp)
        movaps %xmm5,nb330_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%ecx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%ecx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm7    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330_fscal(%esp),%xmm7

        addps  nb330_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb330_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb330_dx(%esp),%xmm0
        movaps nb330_dy(%esp),%xmm1
        movaps nb330_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb330_fix(%esp),%xmm3
        movaps nb330_fiy(%esp),%xmm4
        movaps nb330_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb330_fix(%esp)
        movaps %xmm4,nb330_fiy(%esp)
        movaps %xmm5,nb330_fiz(%esp)
        ## update the fj's 
        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)

        shufps $225,%xmm0,%xmm0 ## constant 11100001
        shufps $225,%xmm1,%xmm1 ## constant 11100001
        shufps $225,%xmm2,%xmm2 ## constant 11100001

        movss   (%edi,%ebx,4),%xmm3
        movss   4(%edi,%ebx,4),%xmm4
        movss   8(%edi,%ebx,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%ebx,4)
        movss   %xmm4,4(%edi,%ebx,4)
        movss   %xmm5,8(%edi,%ebx,4)

_nb_kernel330_ia32_sse.nb330_checksingle:       
        movl  nb330_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330_ia32_sse.nb330_dosingle
        jmp    _nb_kernel330_ia32_sse.nb330_updateouterdata
_nb_kernel330_ia32_sse.nb330_dosingle: 
        movl nb330_charge(%ebp),%esi
        movl nb330_pos(%ebp),%edi
        movl  nb330_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb330_iq(%esp),%xmm6
        movaps %xmm6,nb330_qq(%esp)

        movl nb330_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb330_vdwparam(%ebp),%esi
        shll %ecx
        addl nb330_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb330_c6(%esp)
        movaps %xmm6,nb330_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb330_ix(%esp),%xmm4
        movaps nb330_iy(%esp),%xmm5
        movaps nb330_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb330_dx(%esp)
        movaps %xmm5,nb330_dy(%esp)
        movaps %xmm6,nb330_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb330_VFtab(%ebp),%esi
        movd %mm6,%ebx

        leal (%ebx,%ebx,2),%ebx

        movlps (%esi,%ebx,4),%xmm4
        movlps 8(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        movaps nb330_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb330_vctot(%esp),%xmm5
        movss %xmm5,nb330_vctot(%esp)

        ## put scalar force on stack temporarily 
        movaps %xmm3,nb330_fscal(%esp)

        ## dispersion 
        movlps 16(%esi,%ebx,4),%xmm4
        movlps 24(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        addps  nb330_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb330_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb330_fscal(%esp)
        movss %xmm5,nb330_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%ebx,4),%xmm4
        movlps 40(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb330_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330_fscal(%esp),%xmm7

        addss  nb330_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb330_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb330_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4
        movl   nb330_faction(%ebp),%edi

        movaps nb330_dx(%esp),%xmm0
        movaps nb330_dy(%esp),%xmm1
        movaps nb330_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb330_fix(%esp),%xmm3
        movaps nb330_fiy(%esp),%xmm4
        movaps nb330_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb330_fix(%esp)
        movaps %xmm4,nb330_fiy(%esp)
        movaps %xmm5,nb330_fiz(%esp)
        ## update fj 

        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)
_nb_kernel330_ia32_sse.nb330_updateouterdata: 
        movl  nb330_ii3(%esp),%ecx
        movl  nb330_faction(%ebp),%edi
        movl  nb330_fshift(%ebp),%esi
        movl  nb330_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb330_fix(%esp),%xmm0
        movaps nb330_fiy(%esp),%xmm1
        movaps nb330_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## increment fshift force  
        movss  (%esi,%edx,4),%xmm3
        movss  4(%esi,%edx,4),%xmm4
        movss  8(%esi,%edx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%esi,%edx,4)
        movss  %xmm4,4(%esi,%edx,4)
        movss  %xmm5,8(%esi,%edx,4)

        ## get n from stack
        movl nb330_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb330_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb330_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb330_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb330_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb330_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330_ia32_sse.nb330_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330_n(%esp)
        jmp _nb_kernel330_ia32_sse.nb330_outer
_nb_kernel330_ia32_sse.nb330_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330_ia32_sse.nb330_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330_ia32_sse.nb330_threadloop
_nb_kernel330_ia32_sse.nb330_end: 
        emms

        movl nb330_nouter(%esp),%eax
        movl nb330_ninner(%esp),%ebx
        movl nb330_outeriter(%ebp),%ecx
        movl nb330_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb330_salign(%esp),%eax
        addl %eax,%esp
        addl $372,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret








.globl nb_kernel330nf_ia32_sse
.globl _nb_kernel330nf_ia32_sse
nb_kernel330nf_ia32_sse:        
_nb_kernel330nf_ia32_sse:       
.set nb330nf_p_nri, 8
.set nb330nf_iinr, 12
.set nb330nf_jindex, 16
.set nb330nf_jjnr, 20
.set nb330nf_shift, 24
.set nb330nf_shiftvec, 28
.set nb330nf_fshift, 32
.set nb330nf_gid, 36
.set nb330nf_pos, 40
.set nb330nf_faction, 44
.set nb330nf_charge, 48
.set nb330nf_p_facel, 52
.set nb330nf_argkrf, 56
.set nb330nf_argcrf, 60
.set nb330nf_Vc, 64
.set nb330nf_type, 68
.set nb330nf_p_ntype, 72
.set nb330nf_vdwparam, 76
.set nb330nf_Vvdw, 80
.set nb330nf_p_tabscale, 84
.set nb330nf_VFtab, 88
.set nb330nf_invsqrta, 92
.set nb330nf_dvda, 96
.set nb330nf_p_gbtabscale, 100
.set nb330nf_GBtab, 104
.set nb330nf_p_nthreads, 108
.set nb330nf_count, 112
.set nb330nf_mtx, 116
.set nb330nf_outeriter, 120
.set nb330nf_inneriter, 124
.set nb330nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb330nf_ix, 0
.set nb330nf_iy, 16
.set nb330nf_iz, 32
.set nb330nf_iq, 48
.set nb330nf_tsc, 64
.set nb330nf_qq, 80
.set nb330nf_c6, 96
.set nb330nf_c12, 112
.set nb330nf_vctot, 128
.set nb330nf_Vvdwtot, 144
.set nb330nf_half, 160
.set nb330nf_three, 176
.set nb330nf_is3, 192
.set nb330nf_ii3, 196
.set nb330nf_ntia, 200
.set nb330nf_innerjjnr, 204
.set nb330nf_innerk, 208
.set nb330nf_n, 212
.set nb330nf_nn1, 216
.set nb330nf_nri, 220
.set nb330nf_facel, 224
.set nb330nf_ntype, 228
.set nb330nf_nouter, 232
.set nb330nf_ninner, 236
.set nb330nf_salign, 240
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $244,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb330nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb330nf_p_nri(%ebp),%ecx
        movl nb330nf_p_facel(%ebp),%esi
        movl nb330nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb330nf_nri(%esp)
        movl %esi,nb330nf_facel(%esp)
        movl %edi,nb330nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb330nf_nouter(%esp)
        movl %eax,nb330nf_ninner(%esp)


        movl nb330nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb330nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb330nf_half(%esp)
        movss nb330nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb330nf_half(%esp)
        movaps %xmm3,nb330nf_three(%esp)

_nb_kernel330nf_ia32_sse.nb330nf_threadloop: 
        movl  nb330nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel330nf_ia32_sse.nb330nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel330nf_ia32_sse.nb330nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb330nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb330nf_n(%esp)
        movl %ebx,nb330nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel330nf_ia32_sse.nb330nf_outerstart
        jmp _nb_kernel330nf_ia32_sse.nb330nf_end

_nb_kernel330nf_ia32_sse.nb330nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb330nf_nouter(%esp),%ebx
        movl %ebx,nb330nf_nouter(%esp)

_nb_kernel330nf_ia32_sse.nb330nf_outer: 
        movl  nb330nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb330nf_is3(%esp)            ## store is3 

        movl  nb330nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb330nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb330nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb330nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb330nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb330nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb330nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb330nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb330nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb330nf_ix(%esp)
        movaps %xmm1,nb330nf_iy(%esp)
        movaps %xmm2,nb330nf_iz(%esp)

        movl  %ebx,nb330nf_ii3(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb330nf_vctot(%esp)
        movaps %xmm4,nb330nf_Vvdwtot(%esp)

        movl  nb330nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb330nf_pos(%ebp),%esi
        movl  nb330nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb330nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb330nf_ninner(%esp),%ecx
        movl  %ecx,nb330nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb330nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel330nf_ia32_sse.nb330nf_unroll_loop
        jmp   _nb_kernel330nf_ia32_sse.nb330nf_finish_inner
_nb_kernel330nf_ia32_sse.nb330nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb330nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb330nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb330nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb330nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb330nf_qq(%esp)

        movl nb330nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb330nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb330nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%esi,%eax,4),%xmm6
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm6
        movhps (%esi,%edx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb330nf_c6(%esp)
        movaps %xmm6,nb330nf_c12(%esp)

        movl nb330nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move four coordinates to xmm0-xmm2   

        movlps (%esi,%eax,4),%xmm4
        movlps (%esi,%ecx,4),%xmm5
        movss 8(%esi,%eax,4),%xmm2
        movss 8(%esi,%ecx,4),%xmm6

        movhps (%esi,%ebx,4),%xmm4
        movhps (%esi,%edx,4),%xmm5

        movss 8(%esi,%ebx,4),%xmm0
        movss 8(%esi,%edx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ix-iz to xmm4-xmm6 
        movaps nb330nf_ix(%esp),%xmm4
        movaps nb330nf_iy(%esp),%xmm5
        movaps nb330nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%esp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb330nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330nf_vctot(%esp),%xmm5
        movaps %xmm5,nb330nf_vctot(%esp)

        ## dispersion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 
        ## put scalar force on stack 
        addps  nb330nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%eax,4),%xmm5
        movlps 32(%esi,%ecx,4),%xmm7
        movhps 32(%esi,%ebx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%eax,4),%xmm7
        movlps 40(%esi,%ecx,4),%xmm3
        movhps 40(%esi,%ebx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $4,nb330nf_innerk(%esp)
        jl    _nb_kernel330nf_ia32_sse.nb330nf_finish_inner
        jmp   _nb_kernel330nf_ia32_sse.nb330nf_unroll_loop
_nb_kernel330nf_ia32_sse.nb330nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb330nf_innerk(%esp)
        movl  nb330nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel330nf_ia32_sse.nb330nf_dopair
        jmp   _nb_kernel330nf_ia32_sse.nb330nf_checksingle
_nb_kernel330nf_ia32_sse.nb330nf_dopair: 
        movl nb330nf_charge(%ebp),%esi

        movl  nb330nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb330nf_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb330nf_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb330nf_qq(%esp)

        movl nb330nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb330nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb330nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb330nf_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb330nf_c6(%esp)
        movaps %xmm6,nb330nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101


        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb330nf_ix(%esp),%xmm4
        movaps nb330nf_iy(%esp),%xmm5
        movaps nb330nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb330nf_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx
        leal  (%ecx,%ecx,2),%ecx
        leal  (%edx,%edx,2),%edx

        movlps (%esi,%ecx,4),%xmm5
        movhps (%esi,%edx,4),%xmm5 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%ecx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## constant 10001000
        shufps $221,%xmm7,%xmm7 ## constant 11011101
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb330nf_vctot(%esp),%xmm5
        movaps %xmm5,nb330nf_vctot(%esp)

        ## dispersion 
        movlps 16(%esi,%ecx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm5,%xmm5 ## constant 11011101

        movlps 24(%esi,%ecx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm7    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## constant 10001000
        shufps $221,%xmm7,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 
        ## put scalar force on stack 
        addps  nb330nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%ecx,4),%xmm5
        movhps 32(%esi,%edx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 40(%esi,%ecx,4),%xmm7
        movhps 40(%esi,%edx,4),%xmm7    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb330nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb330nf_Vvdwtot(%esp)

_nb_kernel330nf_ia32_sse.nb330nf_checksingle:   
        movl  nb330nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel330nf_ia32_sse.nb330nf_dosingle
        jmp    _nb_kernel330nf_ia32_sse.nb330nf_updateouterdata
_nb_kernel330nf_ia32_sse.nb330nf_dosingle: 
        movl nb330nf_charge(%ebp),%esi
        movl nb330nf_pos(%ebp),%edi
        movl  nb330nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb330nf_iq(%esp),%xmm6
        movaps %xmm6,nb330nf_qq(%esp)

        movl nb330nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb330nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb330nf_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb330nf_c6(%esp)
        movaps %xmm6,nb330nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb330nf_ix(%esp),%xmm4
        movaps nb330nf_iy(%esp),%xmm5
        movaps nb330nf_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb330nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb330nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb330nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb330nf_VFtab(%ebp),%esi
        movd %mm6,%ebx

        leal (%ebx,%ebx,2),%ebx

        movlps (%esi,%ebx,4),%xmm4
        movlps 8(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb330nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb330nf_vctot(%esp),%xmm5
        movss %xmm5,nb330nf_vctot(%esp)

        ## dispersion 
        movlps 16(%esi,%ebx,4),%xmm4
        movlps 24(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack  
        addss  nb330nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb330nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 32(%esi,%ebx,4),%xmm4
        movlps 40(%esi,%ebx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb330nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addss  nb330nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb330nf_Vvdwtot(%esp)

_nb_kernel330nf_ia32_sse.nb330nf_updateouterdata: 
        ## get n from stack
        movl nb330nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb330nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb330nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb330nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb330nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb330nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb330nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel330nf_ia32_sse.nb330nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb330nf_n(%esp)
        jmp _nb_kernel330nf_ia32_sse.nb330nf_outer
_nb_kernel330nf_ia32_sse.nb330nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb330nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel330nf_ia32_sse.nb330nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel330nf_ia32_sse.nb330nf_threadloop
_nb_kernel330nf_ia32_sse.nb330nf_end: 
        emms

        movl nb330nf_nouter(%esp),%eax
        movl nb330nf_ninner(%esp),%ebx
        movl nb330nf_outeriter(%ebp),%ecx
        movl nb330nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb330nf_salign(%esp),%eax
        addl %eax,%esp
        addl $244,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




