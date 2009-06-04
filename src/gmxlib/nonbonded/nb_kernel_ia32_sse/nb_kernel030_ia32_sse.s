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



.globl nb_kernel030_ia32_sse
.globl _nb_kernel030_ia32_sse
nb_kernel030_ia32_sse:  
_nb_kernel030_ia32_sse: 
.set nb030_p_nri, 8
.set nb030_iinr, 12
.set nb030_jindex, 16
.set nb030_jjnr, 20
.set nb030_shift, 24
.set nb030_shiftvec, 28
.set nb030_fshift, 32
.set nb030_gid, 36
.set nb030_pos, 40
.set nb030_faction, 44
.set nb030_charge, 48
.set nb030_p_facel, 52
.set nb030_p_krf, 56
.set nb030_p_crf, 60
.set nb030_Vc, 64
.set nb030_type, 68
.set nb030_p_ntype, 72
.set nb030_vdwparam, 76
.set nb030_Vvdw, 80
.set nb030_p_tabscale, 84
.set nb030_VFtab, 88
.set nb030_invsqrta, 92
.set nb030_dvda, 96
.set nb030_p_gbtabscale, 100
.set nb030_GBtab, 104
.set nb030_p_nthreads, 108
.set nb030_count, 112
.set nb030_mtx, 116
.set nb030_outeriter, 120
.set nb030_inneriter, 124
.set nb030_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030_ix, 0
.set nb030_iy, 16
.set nb030_iz, 32
.set nb030_dx, 48
.set nb030_dy, 64
.set nb030_dz, 80
.set nb030_two, 96
.set nb030_tsc, 112
.set nb030_c6, 128
.set nb030_c12, 144
.set nb030_fscal, 160
.set nb030_Vvdwtot, 176
.set nb030_fix, 192
.set nb030_fiy, 208
.set nb030_fiz, 224
.set nb030_half, 240
.set nb030_three, 256
.set nb030_is3, 272
.set nb030_ii3, 276
.set nb030_ntia, 280
.set nb030_innerjjnr, 284
.set nb030_innerk, 288
.set nb030_n, 292
.set nb030_nn1, 296
.set nb030_nri, 300
.set nb030_ntype, 304
.set nb030_nouter, 308
.set nb030_ninner, 312
.set nb030_salign, 316
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $320,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb030_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb030_p_nri(%ebp),%ecx
        movl nb030_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb030_nri(%esp)
        movl %edi,nb030_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030_nouter(%esp)
        movl %eax,nb030_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb030_half(%esp)
        movss nb030_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movl nb030_p_tabscale(%ebp),%eax
        movss (%eax),%xmm0
        movaps %xmm1,nb030_half(%esp)
        movaps %xmm2,nb030_two(%esp)
        movaps %xmm3,nb030_three(%esp)
        shufps $0,%xmm0,%xmm0
        movaps %xmm0,nb030_tsc(%esp)

_nb_kernel030_ia32_sse.nb030_threadloop: 
        movl  nb030_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030_ia32_sse.nb030_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030_ia32_sse.nb030_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030_n(%esp)
        movl %ebx,nb030_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030_ia32_sse.nb030_outerstart
        jmp _nb_kernel030_ia32_sse.nb030_end

_nb_kernel030_ia32_sse.nb030_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030_nouter(%esp),%ebx
        movl %ebx,nb030_nouter(%esp)

_nb_kernel030_ia32_sse.nb030_outer: 
        movl  nb030_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb030_is3(%esp)      ## store is3 

        movl  nb030_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb030_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb030_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb030_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb030_ix(%esp)
        movaps %xmm1,nb030_iy(%esp)
        movaps %xmm2,nb030_iz(%esp)

        movl  %ebx,nb030_ii3(%esp)

        ## clear tot potential and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb030_Vvdwtot(%esp)
        movaps %xmm4,nb030_fix(%esp)
        movaps %xmm4,nb030_fiy(%esp)
        movaps %xmm4,nb030_fiz(%esp)

        movl  nb030_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb030_pos(%ebp),%esi
        movl  nb030_faction(%ebp),%edi
        movl  nb030_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb030_ninner(%esp),%ecx
        movl  %ecx,nb030_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb030_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel030_ia32_sse.nb030_unroll_loop
        jmp   _nb_kernel030_ia32_sse.nb030_finish_inner
_nb_kernel030_ia32_sse.nb030_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb030_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb030_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb030_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb030_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb030_ntia(%esp),%edi
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

        movaps %xmm4,nb030_c6(%esp)
        movaps %xmm6,nb030_c12(%esp)

        movl nb030_pos(%ebp),%esi        ## base of pos[] 

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

        ## move nb030_ix-iz to xmm4-xmm6 
        movaps nb030_ix(%esp),%xmm4
        movaps nb030_iy(%esp),%xmm5
        movaps nb030_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb030_dx(%esp)
        movaps %xmm5,nb030_dy(%esp)
        movaps %xmm6,nb030_dz(%esp)
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
        movaps nb030_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030_tsc(%esp),%xmm4

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
        pslld $3,%mm6
        pslld $3,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb030_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps 0(%esi,%eax,4),%xmm5
        movlps 0(%esi,%ecx,4),%xmm7
        movhps 0(%esi,%ebx,4),%xmm5
        movhps 0(%esi,%edx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb030_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb030_fscal(%esp)
        movaps %xmm5,nb030_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7
        mulps  %xmm4,%xmm5
        addps  nb030_fscal(%esp),%xmm7

        addps  nb030_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb030_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb030_dx(%esp),%xmm0
        movaps nb030_dy(%esp),%xmm1
        movaps nb030_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb030_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb030_fix(%esp),%xmm3
        movaps nb030_fiy(%esp),%xmm4
        movaps nb030_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb030_fix(%esp)
        movaps %xmm4,nb030_fiy(%esp)
        movaps %xmm5,nb030_fiz(%esp)
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
        subl $4,nb030_innerk(%esp)
        jl    _nb_kernel030_ia32_sse.nb030_finish_inner
        jmp   _nb_kernel030_ia32_sse.nb030_unroll_loop
_nb_kernel030_ia32_sse.nb030_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb030_innerk(%esp)
        movl  nb030_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel030_ia32_sse.nb030_dopair
        jmp   _nb_kernel030_ia32_sse.nb030_checksingle
_nb_kernel030_ia32_sse.nb030_dopair: 
        movl  nb030_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb030_innerjjnr(%esp)
        xorps %xmm7,%xmm7

        movl nb030_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb030_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb030_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb030_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## shuffle constant  00001000      
        shufps $13,%xmm6,%xmm6 ## suffle constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb030_c6(%esp)
        movaps %xmm6,nb030_c12(%esp)

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

        movl   nb030_faction(%ebp),%edi
        ## move nb030_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb030_ix(%esp),%xmm4
        movaps nb030_iy(%esp),%xmm5
        movaps nb030_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb030_dx(%esp)
        movaps %xmm5,nb030_dy(%esp)
        movaps %xmm6,nb030_dz(%esp)
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
        movaps nb030_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movl nb030_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movlps 0(%esi,%ecx,4),%xmm5
        movhps 0(%esi,%edx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm5,%xmm5 ## constant 11011101

        movlps 8(%esi,%ecx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm7    ## other half of dispersion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## constant 10001000
        shufps $221,%xmm7,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb030_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb030_fscal(%esp)
        movaps %xmm5,nb030_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%ecx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%ecx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm7    ## other half of repulsion table 
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb030_fscal(%esp),%xmm7

        addps  nb030_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb030_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4

        movaps nb030_dx(%esp),%xmm0
        movaps nb030_dy(%esp),%xmm1
        movaps nb030_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb030_fix(%esp),%xmm3
        movaps nb030_fiy(%esp),%xmm4
        movaps nb030_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb030_fix(%esp)
        movaps %xmm4,nb030_fiy(%esp)
        movaps %xmm5,nb030_fiz(%esp)
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

_nb_kernel030_ia32_sse.nb030_checksingle:       
        movl  nb030_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030_ia32_sse.nb030_dosingle
        jmp    _nb_kernel030_ia32_sse.nb030_updateouterdata
_nb_kernel030_ia32_sse.nb030_dosingle: 
        movl nb030_pos(%ebp),%edi
        movl  nb030_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6

        movl nb030_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb030_vdwparam(%ebp),%esi
        shll %ecx
        addl nb030_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb030_c6(%esp)
        movaps %xmm6,nb030_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb030_ix(%esp),%xmm4
        movaps nb030_iy(%esp),%xmm5
        movaps nb030_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb030_dx(%esp)
        movaps %xmm5,nb030_dy(%esp)
        movaps %xmm6,nb030_dz(%esp)
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
        movaps nb030_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movl nb030_VFtab(%ebp),%esi
        movd %mm6,%ebx

        ## dispersion 
        movlps 0(%esi,%ebx,4),%xmm4
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
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb030_Vvdwtot(%esp),%xmm5
        movaps %xmm7,nb030_fscal(%esp)
        movss %xmm5,nb030_Vvdwtot(%esp)

        ## repulsion 
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
        mulps  nb030_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        addps  nb030_fscal(%esp),%xmm7

        addss  nb030_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb030_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps nb030_tsc(%esp),%xmm7
        mulps %xmm0,%xmm7
        subps  %xmm7,%xmm4
        movl   nb030_faction(%ebp),%edi

        movaps nb030_dx(%esp),%xmm0
        movaps nb030_dy(%esp),%xmm1
        movaps nb030_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb030_fix(%esp),%xmm3
        movaps nb030_fiy(%esp),%xmm4
        movaps nb030_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb030_fix(%esp)
        movaps %xmm4,nb030_fiy(%esp)
        movaps %xmm5,nb030_fiz(%esp)
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
_nb_kernel030_ia32_sse.nb030_updateouterdata: 
        movl  nb030_ii3(%esp),%ecx
        movl  nb030_faction(%ebp),%edi
        movl  nb030_fshift(%ebp),%esi
        movl  nb030_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb030_fix(%esp),%xmm0
        movaps nb030_fiy(%esp),%xmm1
        movaps nb030_fiz(%esp),%xmm2

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
        movl nb030_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb030_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb030_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb030_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030_ia32_sse.nb030_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030_n(%esp)
        jmp _nb_kernel030_ia32_sse.nb030_outer
_nb_kernel030_ia32_sse.nb030_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030_ia32_sse.nb030_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030_ia32_sse.nb030_threadloop
_nb_kernel030_ia32_sse.nb030_end: 
        emms

        movl nb030_nouter(%esp),%eax
        movl nb030_ninner(%esp),%ebx
        movl nb030_outeriter(%ebp),%ecx
        movl nb030_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb030_salign(%esp),%eax
        addl %eax,%esp
        addl $320,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel030nf_ia32_sse
.globl _nb_kernel030nf_ia32_sse
nb_kernel030nf_ia32_sse:        
_nb_kernel030nf_ia32_sse:       
.set nb030nf_p_nri, 8
.set nb030nf_iinr, 12
.set nb030nf_jindex, 16
.set nb030nf_jjnr, 20
.set nb030nf_shift, 24
.set nb030nf_shiftvec, 28
.set nb030nf_fshift, 32
.set nb030nf_gid, 36
.set nb030nf_pos, 40
.set nb030nf_faction, 44
.set nb030nf_charge, 48
.set nb030nf_p_facel, 52
.set nb030nf_p_krf, 56
.set nb030nf_p_crf, 60
.set nb030nf_Vc, 64
.set nb030nf_type, 68
.set nb030nf_p_ntype, 72
.set nb030nf_vdwparam, 76
.set nb030nf_Vvdw, 80
.set nb030nf_p_tabscale, 84
.set nb030nf_VFtab, 88
.set nb030nf_invsqrta, 92
.set nb030nf_dvda, 96
.set nb030nf_p_gbtabscale, 100
.set nb030nf_GBtab, 104
.set nb030nf_p_nthreads, 108
.set nb030nf_count, 112
.set nb030nf_mtx, 116
.set nb030nf_outeriter, 120
.set nb030nf_inneriter, 124
.set nb030nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb030nf_ix, 0
.set nb030nf_iy, 16
.set nb030nf_iz, 32
.set nb030nf_tsc, 48
.set nb030nf_c6, 64
.set nb030nf_c12, 80
.set nb030nf_Vvdwtot, 96
.set nb030nf_half, 112
.set nb030nf_three, 128
.set nb030nf_is3, 144
.set nb030nf_ii3, 148
.set nb030nf_ntia, 152
.set nb030nf_innerjjnr, 156
.set nb030nf_innerk, 160
.set nb030nf_n, 164
.set nb030nf_nn1, 168
.set nb030nf_nri, 172
.set nb030nf_ntype, 176
.set nb030nf_nouter, 180
.set nb030nf_ninner, 184
.set nb030nf_salign, 188
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $192,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb030nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb030nf_p_nri(%ebp),%ecx
        movl nb030nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb030nf_nri(%esp)
        movl %edi,nb030nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb030nf_nouter(%esp)
        movl %eax,nb030nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb030nf_half(%esp)
        movss nb030nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movl nb030nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm0
        movaps %xmm1,nb030nf_half(%esp)
        movaps %xmm3,nb030nf_three(%esp)
        shufps $0,%xmm0,%xmm0
        movaps %xmm0,nb030nf_tsc(%esp)

_nb_kernel030nf_ia32_sse.nb030nf_threadloop: 
        movl  nb030nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel030nf_ia32_sse.nb030nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock  
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel030nf_ia32_sse.nb030nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb030nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb030nf_n(%esp)
        movl %ebx,nb030nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel030nf_ia32_sse.nb030nf_outerstart
        jmp _nb_kernel030nf_ia32_sse.nb030nf_end

_nb_kernel030nf_ia32_sse.nb030nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb030nf_nouter(%esp),%ebx
        movl %ebx,nb030nf_nouter(%esp)

_nb_kernel030nf_ia32_sse.nb030nf_outer: 
        movl  nb030nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb030nf_is3(%esp)            ## store is3 

        movl  nb030nf_shiftvec(%ebp),%eax       ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb030nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        movl  nb030nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb030nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb030nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx                ## ebx = 3*ii=ii3 
        movl  nb030nf_pos(%ebp),%eax            ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb030nf_ix(%esp)
        movaps %xmm1,nb030nf_iy(%esp)
        movaps %xmm2,nb030nf_iz(%esp)

        movl  %ebx,nb030nf_ii3(%esp)

        ## clear tot potential and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb030nf_Vvdwtot(%esp)

        movl  nb030nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movl  nb030nf_pos(%ebp),%esi
        movl  nb030nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb030nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb030nf_ninner(%esp),%ecx
        movl  %ecx,nb030nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb030nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel030nf_ia32_sse.nb030nf_unroll_loop
        jmp   _nb_kernel030nf_ia32_sse.nb030nf_finish_inner
_nb_kernel030nf_ia32_sse.nb030nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb030nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb030nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb030nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb030nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb030nf_ntia(%esp),%edi
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

        movaps %xmm4,nb030nf_c6(%esp)
        movaps %xmm6,nb030nf_c12(%esp)

        movl nb030nf_pos(%ebp),%esi        ## base of pos[] 

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

        ## move nb030nf_ix-iz to xmm4-xmm6 
        movaps nb030nf_ix(%esp),%xmm4
        movaps nb030nf_iy(%esp),%xmm5
        movaps nb030nf_iz(%esp),%xmm6

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
        movaps nb030nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%esp),%xmm4

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
        pslld $3,%mm6
        pslld $3,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movl nb030nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps 0(%esi,%eax,4),%xmm5
        movlps 0(%esi,%ecx,4),%xmm7
        movhps 0(%esi,%ebx,4),%xmm5
        movhps 0(%esi,%edx,4),%xmm7    ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movlps 8(%esi,%ecx,4),%xmm3
        movhps 8(%esi,%ebx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm3    ## other half of dispersion table 
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

        movaps nb030nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot 
        addps  nb030nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movlps 16(%esi,%ecx,4),%xmm7
        movhps 16(%esi,%ebx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movlps 24(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ebx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm3    ## other half of repulsion table 
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

        movaps nb030nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5
        addps  nb030nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        ## should we do one more iteration? 
        subl $4,nb030nf_innerk(%esp)
        jl    _nb_kernel030nf_ia32_sse.nb030nf_finish_inner
        jmp   _nb_kernel030nf_ia32_sse.nb030nf_unroll_loop
_nb_kernel030nf_ia32_sse.nb030nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb030nf_innerk(%esp)
        movl  nb030nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel030nf_ia32_sse.nb030nf_dopair
        jmp   _nb_kernel030nf_ia32_sse.nb030nf_checksingle
_nb_kernel030nf_ia32_sse.nb030nf_dopair: 
        movl  nb030nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb030nf_innerjjnr(%esp)
        xorps %xmm7,%xmm7

        movl nb030nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb030nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb030nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb030nf_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## shuffle constant 00001000       
        shufps $13,%xmm6,%xmm6 ## shuffle constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb030nf_c6(%esp)
        movaps %xmm6,nb030nf_c12(%esp)

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

        ## move nb030nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb030nf_ix(%esp),%xmm4
        movaps nb030nf_iy(%esp),%xmm5
        movaps nb030nf_iz(%esp),%xmm6

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
        movaps nb030nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movl nb030nf_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## dispersion 
        movlps 0(%esi,%ecx,4),%xmm5
        movhps 0(%esi,%edx,4),%xmm5   ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## constant 10001000
        shufps $221,%xmm5,%xmm5 ## constant 11011101

        movlps 8(%esi,%ecx,4),%xmm7
        movhps 8(%esi,%edx,4),%xmm7    ## other half of dispersion table 
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

        movaps nb030nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ##  Update Vvdwtot  
        addps  nb030nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%ecx,4),%xmm5
        movhps 16(%esi,%edx,4),%xmm5    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%ecx,4),%xmm7
        movhps 24(%esi,%edx,4),%xmm7    ## other half of repulsion table 
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

        movaps nb030nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb030nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb030nf_Vvdwtot(%esp)

_nb_kernel030nf_ia32_sse.nb030nf_checksingle:   
        movl  nb030nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel030nf_ia32_sse.nb030nf_dosingle
        jmp    _nb_kernel030nf_ia32_sse.nb030nf_updateouterdata
_nb_kernel030nf_ia32_sse.nb030nf_dosingle: 
        movl nb030nf_pos(%ebp),%edi
        movl  nb030nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6

        movl nb030nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb030nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb030nf_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb030nf_c6(%esp)
        movaps %xmm6,nb030nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb030nf_ix(%esp),%xmm4
        movaps nb030nf_iy(%esp),%xmm5
        movaps nb030nf_iz(%esp),%xmm6

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
        movaps nb030nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb030nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb030nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6

        movl nb030nf_VFtab(%ebp),%esi
        movd %mm6,%ebx

        ## dispersion 
        movlps 0(%esi,%ebx,4),%xmm4
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
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb030nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb030nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb030nf_Vvdwtot(%esp)

        ## repulsion 
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

        movaps nb030nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addss  nb030nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb030nf_Vvdwtot(%esp)

_nb_kernel030nf_ia32_sse.nb030nf_updateouterdata: 
        ## get n from stack
        movl nb030nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb030nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb030nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb030nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb030nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel030nf_ia32_sse.nb030nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb030nf_n(%esp)
        jmp _nb_kernel030nf_ia32_sse.nb030nf_outer
_nb_kernel030nf_ia32_sse.nb030nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb030nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel030nf_ia32_sse.nb030nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel030nf_ia32_sse.nb030nf_threadloop
_nb_kernel030nf_ia32_sse.nb030nf_end: 
        emms

        movl nb030nf_nouter(%esp),%eax
        movl nb030nf_ninner(%esp),%ebx
        movl nb030nf_outeriter(%ebp),%ecx
        movl nb030nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb030nf_salign(%esp),%eax
        addl %eax,%esp
        addl $192,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

