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



.globl nb_kernel310_ia32_sse
.globl _nb_kernel310_ia32_sse
nb_kernel310_ia32_sse:  
_nb_kernel310_ia32_sse: 
.set nb310_p_nri, 8
.set nb310_iinr, 12
.set nb310_jindex, 16
.set nb310_jjnr, 20
.set nb310_shift, 24
.set nb310_shiftvec, 28
.set nb310_fshift, 32
.set nb310_gid, 36
.set nb310_pos, 40
.set nb310_faction, 44
.set nb310_charge, 48
.set nb310_p_facel, 52
.set nb310_argkrf, 56
.set nb310_argcrf, 60
.set nb310_Vc, 64
.set nb310_type, 68
.set nb310_p_ntype, 72
.set nb310_vdwparam, 76
.set nb310_Vvdw, 80
.set nb310_p_tabscale, 84
.set nb310_VFtab, 88
.set nb310_invsqrta, 92
.set nb310_dvda, 96
.set nb310_p_gbtabscale, 100
.set nb310_GBtab, 104
.set nb310_p_nthreads, 108
.set nb310_count, 112
.set nb310_mtx, 116
.set nb310_outeriter, 120
.set nb310_inneriter, 124
.set nb310_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb310_ix, 0
.set nb310_iy, 16
.set nb310_iz, 32
.set nb310_iq, 48
.set nb310_dx, 64
.set nb310_dy, 80
.set nb310_dz, 96
.set nb310_two, 112
.set nb310_six, 128
.set nb310_twelve, 144
.set nb310_tsc, 160
.set nb310_qq, 176
.set nb310_c6, 192
.set nb310_c12, 208
.set nb310_fscal, 224
.set nb310_vctot, 240
.set nb310_Vvdwtot, 256
.set nb310_fix, 272
.set nb310_fiy, 288
.set nb310_fiz, 304
.set nb310_half, 320
.set nb310_three, 336
.set nb310_is3, 352
.set nb310_ii3, 356
.set nb310_ntia, 360
.set nb310_innerjjnr, 364
.set nb310_innerk, 368
.set nb310_n, 372
.set nb310_nn1, 376
.set nb310_nri, 380
.set nb310_facel, 384
.set nb310_ntype, 388
.set nb310_nouter, 392
.set nb310_ninner, 396
.set nb310_salign, 400
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $404,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb310_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb310_p_nri(%ebp),%ecx
        movl nb310_p_facel(%ebp),%esi
        movl nb310_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb310_nri(%esp)
        movl %esi,nb310_facel(%esp)
        movl %edi,nb310_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb310_nouter(%esp)
        movl %eax,nb310_ninner(%esp)


        movl nb310_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb310_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb310_half(%esp)
        movss nb310_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## 6.0
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## constant 12.0
        movaps %xmm1,nb310_half(%esp)
        movaps %xmm2,nb310_two(%esp)
        movaps %xmm3,nb310_three(%esp)
        movaps %xmm4,nb310_six(%esp)
        movaps %xmm5,nb310_twelve(%esp)

_nb_kernel310_ia32_sse.nb310_threadloop: 
        movl  nb310_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel310_ia32_sse.nb310_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310_ia32_sse.nb310_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310_n(%esp)
        movl %ebx,nb310_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310_ia32_sse.nb310_outerstart
        jmp _nb_kernel310_ia32_sse.nb310_end

_nb_kernel310_ia32_sse.nb310_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310_nouter(%esp),%ebx
        movl %ebx,nb310_nouter(%esp)

_nb_kernel310_ia32_sse.nb310_outer: 
        movl  nb310_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb310_is3(%esp)      ## store is3 

        movl  nb310_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb310_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb310_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb310_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb310_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb310_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb310_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb310_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb310_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb310_ix(%esp)
        movaps %xmm1,nb310_iy(%esp)
        movaps %xmm2,nb310_iz(%esp)

        movl  %ebx,nb310_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb310_vctot(%esp)
        movaps %xmm4,nb310_Vvdwtot(%esp)
        movaps %xmm4,nb310_fix(%esp)
        movaps %xmm4,nb310_fiy(%esp)
        movaps %xmm4,nb310_fiz(%esp)

        movl  nb310_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb310_pos(%ebp),%esi
        movl  nb310_faction(%ebp),%edi
        movl  nb310_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb310_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb310_ninner(%esp),%ecx
        movl  %ecx,nb310_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb310_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel310_ia32_sse.nb310_unroll_loop
        jmp   _nb_kernel310_ia32_sse.nb310_finish_inner
_nb_kernel310_ia32_sse.nb310_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb310_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb310_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb310_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb310_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb310_qq(%esp)

        movl nb310_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb310_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb310_ntia(%esp),%edi
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

        movaps %xmm4,nb310_c6(%esp)
        movaps %xmm6,nb310_c12(%esp)

        movl nb310_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb310_ix(%esp),%xmm4
        movaps nb310_iy(%esp),%xmm5
        movaps nb310_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb310_dx(%esp)
        movaps %xmm5,nb310_dy(%esp)
        movaps %xmm6,nb310_dz(%esp)
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
        movaps nb310_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310_tsc(%esp),%xmm4

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

        movl nb310_VFtab(%ebp),%esi
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
        mulps  nb310_two(%esp),%xmm7    ## two*Heps2 
        movaps nb310_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310_vctot(%esp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movaps %xmm5,nb310_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310_c6(%esp),%xmm6
        mulps  nb310_c12(%esp),%xmm4
        movaps nb310_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        mulps  nb310_twelve(%esp),%xmm4
        subps  %xmm6,%xmm7
        mulps  nb310_tsc(%esp),%xmm3
        mulps  nb310_six(%esp),%xmm6
        movaps %xmm7,nb310_Vvdwtot(%esp)
        subps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        subps  %xmm3,%xmm4
        mulps  %xmm0,%xmm4

        movaps nb310_dx(%esp),%xmm0
        movaps nb310_dy(%esp),%xmm1
        movaps nb310_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx

        movl   nb310_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb310_fix(%esp),%xmm3
        movaps nb310_fiy(%esp),%xmm4
        movaps nb310_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb310_fix(%esp)
        movaps %xmm4,nb310_fiy(%esp)
        movaps %xmm5,nb310_fiz(%esp)
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
        subl $4,nb310_innerk(%esp)
        jl    _nb_kernel310_ia32_sse.nb310_finish_inner
        jmp   _nb_kernel310_ia32_sse.nb310_unroll_loop
_nb_kernel310_ia32_sse.nb310_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb310_innerk(%esp)
        movl  nb310_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel310_ia32_sse.nb310_dopair
        jmp   _nb_kernel310_ia32_sse.nb310_checksingle
_nb_kernel310_ia32_sse.nb310_dopair: 
        movl nb310_charge(%ebp),%esi
    movl  nb310_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb310_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb310_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb310_qq(%esp)

        movl nb310_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb310_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb310_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb310_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb310_c6(%esp)
        movaps %xmm6,nb310_c12(%esp)

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

        movl   nb310_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb310_ix(%esp),%xmm4
        movaps nb310_iy(%esp),%xmm5
        movaps nb310_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb310_dx(%esp)
        movaps %xmm5,nb310_dy(%esp)
        movaps %xmm6,nb310_dz(%esp)
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
        movaps nb310_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb310_VFtab(%ebp),%esi
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
        mulps  nb310_two(%esp),%xmm7    ## two*Heps2 
        movaps nb310_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310_vctot(%esp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movaps %xmm5,nb310_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310_c6(%esp),%xmm6
        mulps  nb310_c12(%esp),%xmm4
        movaps nb310_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        mulps  nb310_twelve(%esp),%xmm4
        subps  %xmm6,%xmm7
        mulps  nb310_tsc(%esp),%xmm3
        mulps  nb310_six(%esp),%xmm6
        movaps %xmm7,nb310_Vvdwtot(%esp)
        subps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        subps  %xmm3,%xmm4
        mulps  %xmm0,%xmm4

        movaps nb310_dx(%esp),%xmm0
        movaps nb310_dy(%esp),%xmm1
        movaps nb310_dz(%esp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb310_fix(%esp),%xmm3
        movaps nb310_fiy(%esp),%xmm4
        movaps nb310_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb310_fix(%esp)
        movaps %xmm4,nb310_fiy(%esp)
        movaps %xmm5,nb310_fiz(%esp)
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

_nb_kernel310_ia32_sse.nb310_checksingle:       
        movl  nb310_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel310_ia32_sse.nb310_dosingle
        jmp    _nb_kernel310_ia32_sse.nb310_updateouterdata
_nb_kernel310_ia32_sse.nb310_dosingle: 
        movl nb310_charge(%ebp),%esi
        movl nb310_pos(%ebp),%edi
        movl  nb310_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb310_iq(%esp),%xmm6
        movaps %xmm6,nb310_qq(%esp)

        movl nb310_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb310_vdwparam(%ebp),%esi
        shll %ecx
        addl nb310_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb310_c6(%esp)
        movaps %xmm6,nb310_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb310_ix(%esp),%xmm4
        movaps nb310_iy(%esp),%xmm5
        movaps nb310_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb310_dx(%esp)
        movaps %xmm5,nb310_dy(%esp)
        movaps %xmm6,nb310_dz(%esp)
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
        movaps nb310_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb310_VFtab(%ebp),%esi
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
        mulps  nb310_two(%esp),%xmm7    ## two*Heps2 
        movaps nb310_qq(%esp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb310_vctot(%esp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movss %xmm5,nb310_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310_c6(%esp),%xmm6
        mulps  nb310_c12(%esp),%xmm4
        movss nb310_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        mulps  nb310_twelve(%esp),%xmm4
        subps  %xmm6,%xmm7
        mulps  nb310_tsc(%esp),%xmm3
        mulps  nb310_six(%esp),%xmm6
        movss %xmm7,nb310_Vvdwtot(%esp)
        subps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        subps  %xmm3,%xmm4
        mulps  %xmm0,%xmm4

        movaps nb310_dx(%esp),%xmm0
        movaps nb310_dy(%esp),%xmm1
        movaps nb310_dz(%esp),%xmm2

        movl   nb310_faction(%ebp),%edi
        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb310_fix(%esp),%xmm3
        movaps nb310_fiy(%esp),%xmm4
        movaps nb310_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb310_fix(%esp)
        movaps %xmm4,nb310_fiy(%esp)
        movaps %xmm5,nb310_fiz(%esp)
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
_nb_kernel310_ia32_sse.nb310_updateouterdata: 
        movl  nb310_ii3(%esp),%ecx
        movl  nb310_faction(%ebp),%edi
        movl  nb310_fshift(%ebp),%esi
        movl  nb310_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb310_fix(%esp),%xmm0
        movaps nb310_fiy(%esp),%xmm1
        movaps nb310_fiz(%esp),%xmm2

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
        movl nb310_n(%esp),%esi
        ## get group index for i particle 
        movl  nb310_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb310_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb310_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb310_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb310_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb310_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310_ia32_sse.nb310_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310_n(%esp)
        jmp _nb_kernel310_ia32_sse.nb310_outer
_nb_kernel310_ia32_sse.nb310_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310_ia32_sse.nb310_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310_ia32_sse.nb310_threadloop
_nb_kernel310_ia32_sse.nb310_end: 
        emms

        movl nb310_nouter(%esp),%eax
        movl nb310_ninner(%esp),%ebx
        movl nb310_outeriter(%ebp),%ecx
        movl nb310_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb310_salign(%esp),%eax
        addl %eax,%esp
        addl $404,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel310nf_ia32_sse
.globl _nb_kernel310nf_ia32_sse
nb_kernel310nf_ia32_sse:        
_nb_kernel310nf_ia32_sse:       
.set nb310nf_p_nri, 8
.set nb310nf_iinr, 12
.set nb310nf_jindex, 16
.set nb310nf_jjnr, 20
.set nb310nf_shift, 24
.set nb310nf_shiftvec, 28
.set nb310nf_fshift, 32
.set nb310nf_gid, 36
.set nb310nf_pos, 40
.set nb310nf_faction, 44
.set nb310nf_charge, 48
.set nb310nf_p_facel, 52
.set nb310nf_argkrf, 56
.set nb310nf_argcrf, 60
.set nb310nf_Vc, 64
.set nb310nf_type, 68
.set nb310nf_p_ntype, 72
.set nb310nf_vdwparam, 76
.set nb310nf_Vvdw, 80
.set nb310nf_p_tabscale, 84
.set nb310nf_VFtab, 88
.set nb310nf_invsqrta, 92
.set nb310nf_dvda, 96
.set nb310nf_p_gbtabscale, 100
.set nb310nf_GBtab, 104
.set nb310nf_p_nthreads, 108
.set nb310nf_count, 112
.set nb310nf_mtx, 116
.set nb310nf_outeriter, 120
.set nb310nf_inneriter, 124
.set nb310nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb310nf_ix, 0
.set nb310nf_iy, 16
.set nb310nf_iz, 32
.set nb310nf_iq, 48
.set nb310nf_tsc, 64
.set nb310nf_qq, 80
.set nb310nf_c6, 96
.set nb310nf_c12, 112
.set nb310nf_vctot, 128
.set nb310nf_Vvdwtot, 144
.set nb310nf_half, 160
.set nb310nf_three, 176
.set nb310nf_is3, 192
.set nb310nf_ii3, 196
.set nb310nf_ntia, 200
.set nb310nf_innerjjnr, 204
.set nb310nf_innerk, 208
.set nb310nf_n, 212
.set nb310nf_nn1, 216
.set nb310nf_nri, 220
.set nb310nf_facel, 224
.set nb310nf_ntype, 228
.set nb310nf_nouter, 232
.set nb310nf_ninner, 236
.set nb310nf_salign, 240
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
        movl %eax,nb310nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb310nf_p_nri(%ebp),%ecx
        movl nb310nf_p_facel(%ebp),%esi
        movl nb310nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb310nf_nri(%esp)
        movl %esi,nb310nf_facel(%esp)
        movl %edi,nb310nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb310nf_nouter(%esp)
        movl %eax,nb310nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb310nf_half(%esp)
        movss nb310nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb310nf_half(%esp)
        movaps %xmm3,nb310nf_three(%esp)
        movl nb310nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm5
        shufps $0,%xmm5,%xmm5
        movaps %xmm5,nb310nf_tsc(%esp)

_nb_kernel310nf_ia32_sse.nb310nf_threadloop: 
        movl  nb310nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel310nf_ia32_sse.nb310nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310nf_ia32_sse.nb310nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310nf_n(%esp)
        movl %ebx,nb310nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310nf_ia32_sse.nb310nf_outerstart
        jmp _nb_kernel310nf_ia32_sse.nb310nf_end

_nb_kernel310nf_ia32_sse.nb310nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310nf_nouter(%esp),%ebx
        movl %ebx,nb310nf_nouter(%esp)

_nb_kernel310nf_ia32_sse.nb310nf_outer: 
        movl  nb310nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb310nf_is3(%esp)            ## store is3 

        movl  nb310nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb310nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb310nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb310nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb310nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb310nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb310nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb310nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb310nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb310nf_ix(%esp)
        movaps %xmm1,nb310nf_iy(%esp)
        movaps %xmm2,nb310nf_iz(%esp)

        movl  %ebx,nb310nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb310nf_vctot(%esp)
        movaps %xmm4,nb310nf_Vvdwtot(%esp)

        movl  nb310nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb310nf_pos(%ebp),%esi
        movl  nb310nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb310nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb310nf_ninner(%esp),%ecx
        movl  %ecx,nb310nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb310nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel310nf_ia32_sse.nb310nf_unroll_loop
        jmp   _nb_kernel310nf_ia32_sse.nb310nf_finish_inner
_nb_kernel310nf_ia32_sse.nb310nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb310nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb310nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb310nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb310nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb310nf_qq(%esp)

        movl nb310nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb310nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb310nf_ntia(%esp),%edi
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

        movaps %xmm4,nb310nf_c6(%esp)
        movaps %xmm6,nb310nf_c12(%esp)

        movl nb310nf_pos(%ebp),%esi        ## base of pos[] 

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
        movaps nb310nf_ix(%esp),%xmm4
        movaps nb310nf_iy(%esp),%xmm5
        movaps nb310nf_iz(%esp),%xmm6

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
        movaps nb310nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%esp),%xmm4

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

        movl nb310nf_VFtab(%ebp),%esi
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
        movaps nb310nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310nf_vctot(%esp),%xmm5
        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6
        movaps %xmm5,nb310nf_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%esp),%xmm6
        mulps  nb310nf_c12(%esp),%xmm4
        movaps nb310nf_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb310nf_Vvdwtot(%esp)


        ## should we do one more iteration? 
        subl $4,nb310nf_innerk(%esp)
        jl    _nb_kernel310nf_ia32_sse.nb310nf_finish_inner
        jmp   _nb_kernel310nf_ia32_sse.nb310nf_unroll_loop
_nb_kernel310nf_ia32_sse.nb310nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb310nf_innerk(%esp)
        movl  nb310nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel310nf_ia32_sse.nb310nf_dopair
        jmp   _nb_kernel310nf_ia32_sse.nb310nf_checksingle
_nb_kernel310nf_ia32_sse.nb310nf_dopair: 
        movl nb310nf_charge(%ebp),%esi
    movl  nb310nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb310nf_innerjjnr(%esp)
        xorps %xmm7,%xmm7
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## constant 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb310nf_iq(%esp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb310nf_qq(%esp)

        movl nb310nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb310nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb310nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb310nf_pos(%ebp),%edi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb310nf_c6(%esp)
        movaps %xmm6,nb310nf_c12(%esp)

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

        movaps nb310nf_ix(%esp),%xmm4
        movaps nb310nf_iy(%esp),%xmm5
        movaps nb310nf_iz(%esp),%xmm6

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
        movaps nb310nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb310nf_VFtab(%ebp),%esi
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
        movaps nb310nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310nf_vctot(%esp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movaps %xmm5,nb310nf_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%esp),%xmm6
        mulps  nb310nf_c12(%esp),%xmm4
        movaps nb310nf_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb310nf_Vvdwtot(%esp)

_nb_kernel310nf_ia32_sse.nb310nf_checksingle:   
        movl  nb310nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel310nf_ia32_sse.nb310nf_dosingle
        jmp    _nb_kernel310nf_ia32_sse.nb310nf_updateouterdata
_nb_kernel310nf_ia32_sse.nb310nf_dosingle: 
        movl nb310nf_charge(%ebp),%esi
        movl nb310nf_pos(%ebp),%edi
        movl  nb310nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorps  %xmm6,%xmm6
        movss (%esi,%eax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb310nf_iq(%esp),%xmm6
        movaps %xmm6,nb310nf_qq(%esp)

        movl nb310nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb310nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb310nf_ntia(%esp),%ecx
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb310nf_c6(%esp)
        movaps %xmm6,nb310nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        movaps nb310nf_ix(%esp),%xmm4
        movaps nb310nf_iy(%esp),%xmm5
        movaps nb310nf_iz(%esp),%xmm6

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
        movaps nb310nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%esp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%esp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movl nb310nf_VFtab(%ebp),%esi
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
        movaps nb310nf_qq(%esp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb310nf_vctot(%esp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movss %xmm5,nb310nf_vctot(%esp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%esp),%xmm6
        mulps  nb310nf_c12(%esp),%xmm4
        movss nb310nf_Vvdwtot(%esp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movss %xmm7,nb310nf_Vvdwtot(%esp)

_nb_kernel310nf_ia32_sse.nb310nf_updateouterdata: 
        ## get n from stack
        movl nb310nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb310nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb310nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb310nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb310nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb310nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb310nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310nf_ia32_sse.nb310nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310nf_n(%esp)
        jmp _nb_kernel310nf_ia32_sse.nb310nf_outer
_nb_kernel310nf_ia32_sse.nb310nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310nf_ia32_sse.nb310nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310nf_ia32_sse.nb310nf_threadloop
_nb_kernel310nf_ia32_sse.nb310nf_end: 
        emms

        movl nb310nf_nouter(%esp),%eax
        movl nb310nf_ninner(%esp),%ebx
        movl nb310nf_outeriter(%ebp),%ecx
        movl nb310nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb310nf_salign(%esp),%eax
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


