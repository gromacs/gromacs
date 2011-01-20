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


.globl nb_kernel300_ia32_sse
.globl _nb_kernel300_ia32_sse
nb_kernel300_ia32_sse:  
_nb_kernel300_ia32_sse: 
.set nb300_p_nri, 8
.set nb300_iinr, 12
.set nb300_jindex, 16
.set nb300_jjnr, 20
.set nb300_shift, 24
.set nb300_shiftvec, 28
.set nb300_fshift, 32
.set nb300_gid, 36
.set nb300_pos, 40
.set nb300_faction, 44
.set nb300_charge, 48
.set nb300_p_facel, 52
.set nb300_argkrf, 56
.set nb300_argcrf, 60
.set nb300_Vc, 64
.set nb300_type, 68
.set nb300_p_ntype, 72
.set nb300_vdwparam, 76
.set nb300_Vvdw, 80
.set nb300_p_tabscale, 84
.set nb300_VFtab, 88
.set nb300_invsqrta, 92
.set nb300_dvda, 96
.set nb300_p_gbtabscale, 100
.set nb300_GBtab, 104
.set nb300_p_nthreads, 108
.set nb300_count, 112
.set nb300_mtx, 116
.set nb300_outeriter, 120
.set nb300_inneriter, 124
.set nb300_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb300_ix, 0
.set nb300_iy, 16
.set nb300_iz, 32
.set nb300_iq, 48
.set nb300_dx, 64
.set nb300_dy, 80
.set nb300_dz, 96
.set nb300_two, 112
.set nb300_tsc, 128
.set nb300_qq, 144
.set nb300_fs, 160
.set nb300_vctot, 176
.set nb300_fix, 192
.set nb300_fiy, 208
.set nb300_fiz, 224
.set nb300_half, 240
.set nb300_three, 256
.set nb300_is3, 272
.set nb300_ii3, 276
.set nb300_innerjjnr, 280
.set nb300_innerk, 284
.set nb300_n, 288
.set nb300_nn1, 292
.set nb300_nri, 296
.set nb300_facel, 300
.set nb300_nouter, 304
.set nb300_ninner, 308
.set nb300_salign, 312
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $316,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb300_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb300_p_nri(%ebp),%ecx
        movl nb300_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb300_nri(%esp)
        movl %esi,nb300_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300_nouter(%esp)
        movl %eax,nb300_ninner(%esp)


        movl nb300_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb300_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb300_half(%esp)
        movss nb300_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb300_half(%esp)
        movaps %xmm2,nb300_two(%esp)
        movaps %xmm3,nb300_three(%esp)

_nb_kernel300_ia32_sse.nb300_threadloop: 
        movl  nb300_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300_ia32_sse.nb300_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300_ia32_sse.nb300_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300_n(%esp)
        movl %ebx,nb300_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300_ia32_sse.nb300_outerstart
        jmp _nb_kernel300_ia32_sse.nb300_end

_nb_kernel300_ia32_sse.nb300_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300_nouter(%esp),%ebx
        movl %ebx,nb300_nouter(%esp)

_nb_kernel300_ia32_sse.nb300_outer: 
        movl  nb300_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb300_is3(%esp)      ## store is3 

        movl  nb300_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb300_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb300_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb300_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb300_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb300_ix(%esp)
        movaps %xmm1,nb300_iy(%esp)
        movaps %xmm2,nb300_iz(%esp)

        movl  %ebx,nb300_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb300_vctot(%esp)
        movaps %xmm4,nb300_fix(%esp)
        movaps %xmm4,nb300_fiy(%esp)
        movaps %xmm4,nb300_fiz(%esp)

        movl  nb300_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300_pos(%ebp),%esi
        movl  nb300_faction(%ebp),%edi
        movl  nb300_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb300_ninner(%esp),%ecx
        movl  %ecx,nb300_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb300_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel300_ia32_sse.nb300_unroll_loop
        jmp   _nb_kernel300_ia32_sse.nb300_finish_inner
_nb_kernel300_ia32_sse.nb300_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb300_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb300_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb300_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb300_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3

        movaps %xmm3,nb300_qq(%esp)

        movl nb300_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb300_ix(%esp),%xmm4
        movaps nb300_iy(%esp),%xmm5
        movaps nb300_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb300_dx(%esp)
        movaps %xmm5,nb300_dy(%esp)
        movaps %xmm6,nb300_dz(%esp)
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
        movaps nb300_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300_tsc(%esp),%xmm4

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

        movl nb300_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

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
        mulps  nb300_two(%esp),%xmm7    ## two*Heps2 
        movaps nb300_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300_vctot(%esp),%xmm5
        movaps %xmm5,nb300_vctot(%esp)

        xorps  %xmm4,%xmm4

        mulps nb300_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

        movaps nb300_dx(%esp),%xmm0
        movaps nb300_dy(%esp),%xmm1
        movaps nb300_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb300_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb300_fix(%esp),%xmm3
        movaps nb300_fiy(%esp),%xmm4
        movaps nb300_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb300_fix(%esp)
        movaps %xmm4,nb300_fiy(%esp)
        movaps %xmm5,nb300_fiz(%esp)
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
        subl $4,nb300_innerk(%esp)
        jl    _nb_kernel300_ia32_sse.nb300_finish_inner
        jmp   _nb_kernel300_ia32_sse.nb300_unroll_loop
_nb_kernel300_ia32_sse.nb300_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb300_innerk(%esp)
        movl  nb300_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel300_ia32_sse.nb300_dopair
        jmp   _nb_kernel300_ia32_sse.nb300_checksingle
_nb_kernel300_ia32_sse.nb300_dopair: 
        movl nb300_charge(%ebp),%esi

    movl  nb300_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb300_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb300_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb300_qq(%esp)

        movl nb300_pos(%ebp),%edi

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

        movl   nb300_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb300_ix(%esp),%xmm4
        movaps nb300_iy(%esp),%xmm5
        movaps nb300_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb300_dx(%esp)
        movaps %xmm5,nb300_dy(%esp)
        movaps %xmm6,nb300_dz(%esp)
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
        movaps nb300_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb300_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

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
        mulps  nb300_two(%esp),%xmm7    ## two*Heps2 
        movaps nb300_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300_vctot(%esp),%xmm5
        movaps %xmm5,nb300_vctot(%esp)

        xorps  %xmm4,%xmm4

        mulps nb300_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

        movaps nb300_dx(%esp),%xmm0
        movaps nb300_dy(%esp),%xmm1
        movaps nb300_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb300_fix(%esp),%xmm3
        movaps nb300_fiy(%esp),%xmm4
        movaps nb300_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb300_fix(%esp)
        movaps %xmm4,nb300_fiy(%esp)
        movaps %xmm5,nb300_fiz(%esp)
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

_nb_kernel300_ia32_sse.nb300_checksingle:       
        movl  nb300_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300_ia32_sse.nb300_dosingle
        jmp    _nb_kernel300_ia32_sse.nb300_updateouterdata
_nb_kernel300_ia32_sse.nb300_dosingle: 
        movl nb300_charge(%ebp),%esi
        movl nb300_pos(%ebp),%edi
        movl  nb300_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb300_iq(%esp),%xmm6
        movaps %xmm6,nb300_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb300_ix(%esp),%xmm4
        movaps nb300_iy(%esp),%xmm5
        movaps nb300_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb300_dx(%esp)
        movaps %xmm5,nb300_dy(%esp)
        movaps %xmm6,nb300_dz(%esp)
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
        movaps nb300_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb300_VFtab(%ebp),%esi
        movd %mm6,%ebx

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
        mulps  nb300_two(%esp),%xmm7    ## two*Heps2 
        movaps nb300_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb300_vctot(%esp),%xmm5
        movss %xmm5,nb300_vctot(%esp)

        xorps %xmm4,%xmm4

        mulps nb300_tsc(%esp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4
        movl   nb300_faction(%ebp),%edi

        movaps nb300_dx(%esp),%xmm0
        movaps nb300_dy(%esp),%xmm1
        movaps nb300_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb300_fix(%esp),%xmm3
        movaps nb300_fiy(%esp),%xmm4
        movaps nb300_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb300_fix(%esp)
        movaps %xmm4,nb300_fiy(%esp)
        movaps %xmm5,nb300_fiz(%esp)
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
_nb_kernel300_ia32_sse.nb300_updateouterdata: 
        movl  nb300_ii3(%esp),%ecx
        movl  nb300_faction(%ebp),%edi
        movl  nb300_fshift(%ebp),%esi
        movl  nb300_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb300_fix(%esp),%xmm0
        movaps nb300_fiy(%esp),%xmm1
        movaps nb300_fiz(%esp),%xmm2

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
        movl nb300_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb300_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb300_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb300_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300_ia32_sse.nb300_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300_n(%esp)
        jmp _nb_kernel300_ia32_sse.nb300_outer
_nb_kernel300_ia32_sse.nb300_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300_ia32_sse.nb300_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300_ia32_sse.nb300_threadloop
_nb_kernel300_ia32_sse.nb300_end: 
        emms

        movl nb300_nouter(%esp),%eax
        movl nb300_ninner(%esp),%ebx
        movl nb300_outeriter(%ebp),%ecx
        movl nb300_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb300_salign(%esp),%eax
        addl %eax,%esp
        addl $316,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




.globl nb_kernel300nf_ia32_sse
.globl _nb_kernel300nf_ia32_sse
nb_kernel300nf_ia32_sse:        
_nb_kernel300nf_ia32_sse:       
.set nb300nf_p_nri, 8
.set nb300nf_iinr, 12
.set nb300nf_jindex, 16
.set nb300nf_jjnr, 20
.set nb300nf_shift, 24
.set nb300nf_shiftvec, 28
.set nb300nf_fshift, 32
.set nb300nf_gid, 36
.set nb300nf_pos, 40
.set nb300nf_faction, 44
.set nb300nf_charge, 48
.set nb300nf_p_facel, 52
.set nb300nf_argkrf, 56
.set nb300nf_argcrf, 60
.set nb300nf_Vc, 64
.set nb300nf_type, 68
.set nb300nf_p_ntype, 72
.set nb300nf_vdwparam, 76
.set nb300nf_Vvdw, 80
.set nb300nf_p_tabscale, 84
.set nb300nf_VFtab, 88
.set nb300nf_invsqrta, 92
.set nb300nf_dvda, 96
.set nb300nf_p_gbtabscale, 100
.set nb300nf_GBtab, 104
.set nb300nf_p_nthreads, 108
.set nb300nf_count, 112
.set nb300nf_mtx, 116
.set nb300nf_outeriter, 120
.set nb300nf_inneriter, 124
.set nb300nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb300nf_ix, 0
.set nb300nf_iy, 16
.set nb300nf_iz, 32
.set nb300nf_iq, 48
.set nb300nf_tsc, 64
.set nb300nf_qq, 80
.set nb300nf_vctot, 96
.set nb300nf_half, 112
.set nb300nf_three, 128
.set nb300nf_is3, 144
.set nb300nf_ii3, 148
.set nb300nf_innerjjnr, 152
.set nb300nf_innerk, 156
.set nb300nf_n, 160
.set nb300nf_nn1, 164
.set nb300nf_nri, 168
.set nb300nf_facel, 172
.set nb300nf_nouter, 176
.set nb300nf_ninner, 180
.set nb300nf_salign, 184
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $188,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb300nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb300nf_p_nri(%ebp),%ecx
        movl nb300nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl %ecx,nb300nf_nri(%esp)
        movl %esi,nb300nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300nf_nouter(%esp)
        movl %eax,nb300nf_ninner(%esp)


        movl nb300nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb300nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb300nf_half(%esp)
        movss nb300nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb300nf_half(%esp)
        movaps %xmm3,nb300nf_three(%esp)

_nb_kernel300nf_ia32_sse.nb300nf_threadloop: 
        movl  nb300nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300nf_ia32_sse.nb300nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300nf_ia32_sse.nb300nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300nf_n(%esp)
        movl %ebx,nb300nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300nf_ia32_sse.nb300nf_outerstart
        jmp _nb_kernel300nf_ia32_sse.nb300nf_end
_nb_kernel300nf_ia32_sse.nb300nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300nf_nouter(%esp),%ebx
        movl %ebx,nb300nf_nouter(%esp)

_nb_kernel300nf_ia32_sse.nb300nf_outer: 
        movl  nb300nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb300nf_is3(%esp)            ## store is3 

        movl  nb300nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb300nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb300nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb300nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb300nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb300nf_ix(%esp)
        movaps %xmm1,nb300nf_iy(%esp)
        movaps %xmm2,nb300nf_iz(%esp)

        movl  %ebx,nb300nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb300nf_vctot(%esp)

        movl  nb300nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300nf_pos(%ebp),%esi
        movl  nb300nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb300nf_ninner(%esp),%ecx
        movl  %ecx,nb300nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb300nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel300nf_ia32_sse.nb300nf_unroll_loop
        jmp   _nb_kernel300nf_ia32_sse.nb300nf_finish_inner
_nb_kernel300nf_ia32_sse.nb300nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb300nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb300nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb300nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb300nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3

        movaps %xmm3,nb300nf_qq(%esp)

        movl nb300nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb300nf_ix(%esp),%xmm4
        movaps nb300nf_iy(%esp),%xmm5
        movaps nb300nf_iz(%esp),%xmm6

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
        movaps nb300nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%esp),%xmm4

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

        movl nb300nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

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
        movaps nb300nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  

        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300nf_vctot(%esp),%xmm5
        movaps %xmm5,nb300nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb300nf_innerk(%esp)
        jl    _nb_kernel300nf_ia32_sse.nb300nf_finish_inner
        jmp   _nb_kernel300nf_ia32_sse.nb300nf_unroll_loop
_nb_kernel300nf_ia32_sse.nb300nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb300nf_innerk(%esp)
        movl  nb300nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel300nf_ia32_sse.nb300nf_dopair
        jmp   _nb_kernel300nf_ia32_sse.nb300nf_checksingle
_nb_kernel300nf_ia32_sse.nb300nf_dopair: 
        movl nb300nf_charge(%ebp),%esi

    movl  nb300nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb300nf_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb300nf_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb300nf_qq(%esp)

        movl nb300nf_pos(%ebp),%edi

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

        movaps nb300nf_ix(%esp),%xmm4
        movaps nb300nf_iy(%esp),%xmm5
        movaps nb300nf_iz(%esp),%xmm6

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
        movaps nb300nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb300nf_VFtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

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
        movaps nb300nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300nf_vctot(%esp),%xmm5
        movaps %xmm5,nb300nf_vctot(%esp)

_nb_kernel300nf_ia32_sse.nb300nf_checksingle:   
        movl  nb300nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300nf_ia32_sse.nb300nf_dosingle
        jmp    _nb_kernel300nf_ia32_sse.nb300nf_updateouterdata
_nb_kernel300nf_ia32_sse.nb300nf_dosingle: 
        movl nb300nf_charge(%ebp),%esi
        movl nb300nf_pos(%ebp),%edi
        movl  nb300nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb300nf_iq(%esp),%xmm6
        movaps %xmm6,nb300nf_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb300nf_ix(%esp),%xmm4
        movaps nb300nf_iy(%esp),%xmm5
        movaps nb300nf_iz(%esp),%xmm6

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
        movaps nb300nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb300nf_VFtab(%ebp),%esi
        movd %mm6,%ebx

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
        movaps nb300nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb300nf_vctot(%esp),%xmm5
        movss %xmm5,nb300nf_vctot(%esp)

_nb_kernel300nf_ia32_sse.nb300nf_updateouterdata: 
        ## get n from stack
        movl nb300nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb300nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb300nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb300nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300nf_ia32_sse.nb300nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300nf_n(%esp)
        jmp _nb_kernel300nf_ia32_sse.nb300nf_outer
_nb_kernel300nf_ia32_sse.nb300nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300nf_ia32_sse.nb300nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300nf_ia32_sse.nb300nf_threadloop
_nb_kernel300nf_ia32_sse.nb300nf_end: 
        emms

        movl nb300nf_nouter(%esp),%eax
        movl nb300nf_ninner(%esp),%ebx
        movl nb300nf_outeriter(%ebp),%ecx
        movl nb300nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb300nf_salign(%esp),%eax
        addl %eax,%esp
        addl $188,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


