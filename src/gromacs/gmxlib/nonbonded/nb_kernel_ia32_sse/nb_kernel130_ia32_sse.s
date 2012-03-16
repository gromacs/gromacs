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



.globl nb_kernel130_ia32_sse
.globl _nb_kernel130_ia32_sse
nb_kernel130_ia32_sse:  
_nb_kernel130_ia32_sse: 
.set nb130_p_nri, 8
.set nb130_iinr, 12
.set nb130_jindex, 16
.set nb130_jjnr, 20
.set nb130_shift, 24
.set nb130_shiftvec, 28
.set nb130_fshift, 32
.set nb130_gid, 36
.set nb130_pos, 40
.set nb130_faction, 44
.set nb130_charge, 48
.set nb130_p_facel, 52
.set nb130_argkrf, 56
.set nb130_argcrf, 60
.set nb130_Vc, 64
.set nb130_type, 68
.set nb130_p_ntype, 72
.set nb130_vdwparam, 76
.set nb130_Vvdw, 80
.set nb130_p_tabscale, 84
.set nb130_VFtab, 88
.set nb130_invsqrta, 92
.set nb130_dvda, 96
.set nb130_p_gbtabscale, 100
.set nb130_GBtab, 104
.set nb130_p_nthreads, 108
.set nb130_count, 112
.set nb130_mtx, 116
.set nb130_outeriter, 120
.set nb130_inneriter, 124
.set nb130_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb130_ix, 0
.set nb130_iy, 16
.set nb130_iz, 32
.set nb130_iq, 48
.set nb130_dx, 64
.set nb130_dy, 80
.set nb130_dz, 96
.set nb130_c6, 112
.set nb130_c12, 128
.set nb130_tsc, 144
.set nb130_fstmp, 160
.set nb130_vctot, 176
.set nb130_Vvdwtot, 192
.set nb130_fix, 208
.set nb130_fiy, 224
.set nb130_fiz, 240
.set nb130_half, 256
.set nb130_three, 272
.set nb130_two, 288
.set nb130_krf, 304
.set nb130_crf, 320
.set nb130_is3, 336
.set nb130_ii3, 340
.set nb130_ntia, 344
.set nb130_innerjjnr, 348
.set nb130_innerk, 352
.set nb130_n, 356
.set nb130_nn1, 360
.set nb130_nri, 364
.set nb130_facel, 368
.set nb130_ntype, 372
.set nb130_nouter, 376
.set nb130_ninner, 380
.set nb130_salign, 384
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $400,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb130_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb130_p_nri(%ebp),%ecx
        movl nb130_p_facel(%ebp),%esi
        movl nb130_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb130_nri(%esp)
        movl %esi,nb130_facel(%esp)
        movl %edi,nb130_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb130_nouter(%esp)
        movl %eax,nb130_ninner(%esp)

        movl nb130_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb130_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb130_half(%esp)
        movss nb130_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb130_half(%esp)
        movaps %xmm2,nb130_two(%esp)
        movaps %xmm3,nb130_three(%esp)

_nb_kernel130_ia32_sse.nb130_threadloop: 
        movl  nb130_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel130_ia32_sse.nb130_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130_ia32_sse.nb130_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130_n(%esp)
        movl %ebx,nb130_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130_ia32_sse.nb130_outerstart
        jmp _nb_kernel130_ia32_sse.nb130_end

_nb_kernel130_ia32_sse.nb130_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130_nouter(%esp),%ebx
        movl %ebx,nb130_nouter(%esp)

_nb_kernel130_ia32_sse.nb130_outer: 
        movl  nb130_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb130_is3(%esp)      ## store is3 

        movl  nb130_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb130_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb130_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb130_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb130_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb130_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb130_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb130_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb130_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb130_ix(%esp)
        movaps %xmm1,nb130_iy(%esp)
        movaps %xmm2,nb130_iz(%esp)

        movl  %ebx,nb130_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb130_vctot(%esp)
        movaps %xmm4,nb130_Vvdwtot(%esp)
        movaps %xmm4,nb130_fix(%esp)
        movaps %xmm4,nb130_fiy(%esp)
        movaps %xmm4,nb130_fiz(%esp)

        movl  nb130_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb130_pos(%ebp),%esi
        movl  nb130_faction(%ebp),%edi
        movl  nb130_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb130_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb130_ninner(%esp),%ecx
        movl  %ecx,nb130_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb130_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel130_ia32_sse.nb130_unroll_loop
        jmp   _nb_kernel130_ia32_sse.nb130_finish_inner
_nb_kernel130_ia32_sse.nb130_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb130_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb130_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb130_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb130_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb130_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb130_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb130_ntia(%esp),%edi
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

        movaps %xmm4,nb130_c6(%esp)
        movaps %xmm6,nb130_c12(%esp)

        movl nb130_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
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
        movaps nb130_ix(%esp),%xmm4
        movaps nb130_iy(%esp),%xmm5
        movaps nb130_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb130_dx(%esp)
        movaps %xmm5,nb130_dy(%esp)
        movaps %xmm6,nb130_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb130_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb130_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb130_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        mulps %xmm0,%xmm3   ## vcoul
        mulps %xmm3,%xmm2       ## vcoul*rinv

        addps  nb130_vctot(%esp),%xmm3
        movaps %xmm3,nb130_vctot(%esp)

        movaps %xmm2,nb130_fstmp(%esp)

        ## LJ table
        mulps  %xmm1,%xmm4 ## r
        mulps  nb130_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
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

        movl nb130_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half dispersion table 
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
        mulps  nb130_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb130_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        movaps nb130_fstmp(%esp),%xmm3
        mulps  nb130_tsc(%esp),%xmm7
        subps  %xmm7,%xmm3

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb130_Vvdwtot(%esp),%xmm5
        movaps %xmm3,nb130_fstmp(%esp)
        movaps %xmm5,nb130_Vvdwtot(%esp)

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
        mulps  nb130_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb130_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        movaps nb130_fstmp(%esp),%xmm3
        mulps  nb130_tsc(%esp),%xmm7
        subps  %xmm7,%xmm3

        addps  nb130_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps %xmm0,%xmm3

        movd %mm0,%eax
        movd %mm1,%ebx
        movd %mm2,%ecx
        movd %mm3,%edx


        movaps nb130_dx(%esp),%xmm0
        movaps nb130_dy(%esp),%xmm1
        movaps nb130_dz(%esp),%xmm2

        movl   nb130_faction(%ebp),%edi
        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb130_fix(%esp),%xmm3
        movaps nb130_fiy(%esp),%xmm4
        movaps nb130_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb130_fix(%esp)
        movaps %xmm4,nb130_fiy(%esp)
        movaps %xmm5,nb130_fiz(%esp)
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
        subl $4,nb130_innerk(%esp)
        jl    _nb_kernel130_ia32_sse.nb130_finish_inner
        jmp   _nb_kernel130_ia32_sse.nb130_unroll_loop
_nb_kernel130_ia32_sse.nb130_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb130_innerk(%esp)
        movl  nb130_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel130_ia32_sse.nb130_dopair
        jmp   _nb_kernel130_ia32_sse.nb130_checksingle
_nb_kernel130_ia32_sse.nb130_dopair: 
        movl nb130_charge(%ebp),%esi

    movl  nb130_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb130_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb130_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb130_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb130_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb130_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb130_c6(%esp)
        movaps %xmm6,nb130_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb130_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        movl   nb130_faction(%ebp),%edi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb130_ix(%esp),%xmm4
        movaps nb130_iy(%esp),%xmm5
        movaps nb130_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb130_dx(%esp)
        movaps %xmm5,nb130_dy(%esp)
        movaps %xmm6,nb130_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movaps nb130_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb130_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb130_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2
        mulps %xmm0,%xmm3   ## vcoul
        mulps %xmm3,%xmm2       ## vcoul*rinv

        addps  nb130_vctot(%esp),%xmm3
        movaps %xmm3,nb130_vctot(%esp)

        movaps %xmm2,nb130_fstmp(%esp)

        ## LJ table
        mulps  %xmm1,%xmm4 ## r
        mulps  nb130_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb130_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movhps (%esi,%ebx,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb130_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb130_c6(%esp),%xmm4
        mulps  %xmm4,%xmm7       ## fijD 
        mulps  %xmm4,%xmm5       ## Vvdw6 
        movaps nb130_fstmp(%esp),%xmm3
        mulps  nb130_tsc(%esp),%xmm7
        subps  %xmm7,%xmm3

        ## put scalar force on stack Update Vvdwtot directly 
        addps  nb130_Vvdwtot(%esp),%xmm5
        movaps %xmm3,nb130_fstmp(%esp)
        movaps %xmm5,nb130_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movhps 16(%esi,%ebx,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movhps 24(%esi,%ebx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        mulps  nb130_two(%esp),%xmm7    ## two*Heps2 
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb130_c12(%esp),%xmm4
        mulps  %xmm4,%xmm7 ## fijR 
        mulps  %xmm4,%xmm5 ## Vvdw12 
        movaps nb130_fstmp(%esp),%xmm3
        mulps  nb130_tsc(%esp),%xmm7
        subps  %xmm7,%xmm3

        addps  nb130_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130_Vvdwtot(%esp)
        xorps  %xmm4,%xmm4

        mulps %xmm0,%xmm3

        movd %mm0,%eax
        movd %mm1,%ebx

        movaps nb130_dx(%esp),%xmm0
        movaps nb130_dy(%esp),%xmm1
        movaps nb130_dz(%esp),%xmm2

        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb130_fix(%esp),%xmm3
        movaps nb130_fiy(%esp),%xmm4
        movaps nb130_fiz(%esp),%xmm5
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        addps  %xmm2,%xmm5
        movaps %xmm3,nb130_fix(%esp)
        movaps %xmm4,nb130_fiy(%esp)
        movaps %xmm5,nb130_fiz(%esp)
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

_nb_kernel130_ia32_sse.nb130_checksingle:       
        movl  nb130_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130_ia32_sse.nb130_dosingle
        jmp    _nb_kernel130_ia32_sse.nb130_updateouterdata
_nb_kernel130_ia32_sse.nb130_dosingle:  
        movl nb130_charge(%ebp),%esi
        movl nb130_pos(%ebp),%edi
        movl  nb130_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb130_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb130_vdwparam(%ebp),%esi
        shll %ecx
        addl nb130_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb130_c6(%esp)
        movaps %xmm6,nb130_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb130_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb130_ix(%esp),%xmm4
        movaps nb130_iy(%esp),%xmm5
        movaps nb130_iz(%esp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## store dr 
        movaps %xmm4,nb130_dx(%esp)
        movaps %xmm5,nb130_dy(%esp)
        movaps %xmm6,nb130_dz(%esp)
        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        movss nb130_krf(%esp),%xmm7
        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movss %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb130_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb130_half(%esp),%xmm0
        mulss  %xmm4,%xmm7      ## xmm7=krsq 
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 
        movss %xmm0,%xmm1
        movss %xmm0,%xmm2
        mulss %xmm0,%xmm3   ## vcoul
        mulss %xmm3,%xmm2       ## vcoul*rinv

        addss  nb130_vctot(%esp),%xmm3
        movss %xmm3,nb130_vctot(%esp)

        movss %xmm2,nb130_fstmp(%esp)

        ## LJ table
        mulss  %xmm1,%xmm4 ## r
        mulss  nb130_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movl nb130_VFtab(%ebp),%esi
        movd %mm6,%eax

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb130_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss  nb130_c6(%esp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        movss  nb130_fstmp(%esp),%xmm3
        mulss  nb130_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb130_Vvdwtot(%esp),%xmm5
        movss %xmm3,nb130_fstmp(%esp)
        movss %xmm5,nb130_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb130_two(%esp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss  nb130_c12(%esp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        movss  nb130_fstmp(%esp),%xmm3
        mulss  nb130_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3

        addss  nb130_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb130_Vvdwtot(%esp)

        mulss %xmm0,%xmm3

        movd %mm0,%eax

        movaps nb130_dx(%esp),%xmm0
        movaps nb130_dy(%esp),%xmm1
        movaps nb130_dz(%esp),%xmm2

        mulss  %xmm3,%xmm0
        mulss  %xmm3,%xmm1
        mulss  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movaps nb130_fix(%esp),%xmm3
        movaps nb130_fiy(%esp),%xmm4
        movaps nb130_fiz(%esp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movaps %xmm3,nb130_fix(%esp)
        movaps %xmm4,nb130_fiy(%esp)
        movaps %xmm5,nb130_fiz(%esp)
        ## update fj 
        movl nb130_faction(%ebp),%edi

        movss   (%edi,%eax,4),%xmm3
        movss   4(%edi,%eax,4),%xmm4
        movss   8(%edi,%eax,4),%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movss   %xmm3,(%edi,%eax,4)
        movss   %xmm4,4(%edi,%eax,4)
        movss   %xmm5,8(%edi,%eax,4)
_nb_kernel130_ia32_sse.nb130_updateouterdata: 
        movl  nb130_ii3(%esp),%ecx
        movl  nb130_faction(%ebp),%edi
        movl  nb130_fshift(%ebp),%esi
        movl  nb130_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movaps nb130_fix(%esp),%xmm0
        movaps nb130_fiy(%esp),%xmm1
        movaps nb130_fiz(%esp),%xmm2

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
        movl nb130_n(%esp),%esi
        ## get group index for i particle 
        movl  nb130_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb130_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb130_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb130_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb130_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb130_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130_ia32_sse.nb130_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130_n(%esp)
        jmp _nb_kernel130_ia32_sse.nb130_outer
_nb_kernel130_ia32_sse.nb130_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130_ia32_sse.nb130_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130_ia32_sse.nb130_threadloop
_nb_kernel130_ia32_sse.nb130_end: 
        emms

        movl nb130_nouter(%esp),%eax
        movl nb130_ninner(%esp),%ebx
        movl nb130_outeriter(%ebp),%ecx
        movl nb130_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb130_salign(%esp),%eax
        addl %eax,%esp
        addl $400,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret







.globl nb_kernel130nf_ia32_sse
.globl _nb_kernel130nf_ia32_sse
nb_kernel130nf_ia32_sse:        
_nb_kernel130nf_ia32_sse:       
.set nb130nf_p_nri, 8
.set nb130nf_iinr, 12
.set nb130nf_jindex, 16
.set nb130nf_jjnr, 20
.set nb130nf_shift, 24
.set nb130nf_shiftvec, 28
.set nb130nf_fshift, 32
.set nb130nf_gid, 36
.set nb130nf_pos, 40
.set nb130nf_faction, 44
.set nb130nf_charge, 48
.set nb130nf_p_facel, 52
.set nb130nf_argkrf, 56
.set nb130nf_argcrf, 60
.set nb130nf_Vc, 64
.set nb130nf_type, 68
.set nb130nf_p_ntype, 72
.set nb130nf_vdwparam, 76
.set nb130nf_Vvdw, 80
.set nb130nf_p_tabscale, 84
.set nb130nf_VFtab, 88
.set nb130nf_invsqrta, 92
.set nb130nf_dvda, 96
.set nb130nf_p_gbtabscale, 100
.set nb130nf_GBtab, 104
.set nb130nf_p_nthreads, 108
.set nb130nf_count, 112
.set nb130nf_mtx, 116
.set nb130nf_outeriter, 120
.set nb130nf_inneriter, 124
.set nb130nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb130nf_ix, 0
.set nb130nf_iy, 16
.set nb130nf_iz, 32
.set nb130nf_iq, 48
.set nb130nf_dx, 64
.set nb130nf_dy, 80
.set nb130nf_dz, 96
.set nb130nf_c6, 112
.set nb130nf_c12, 128
.set nb130nf_tsc, 144
.set nb130nf_vctot, 160
.set nb130nf_Vvdwtot, 176
.set nb130nf_half, 192
.set nb130nf_three, 208
.set nb130nf_krf, 224
.set nb130nf_crf, 240
.set nb130nf_is3, 256
.set nb130nf_ii3, 260
.set nb130nf_ntia, 264
.set nb130nf_innerjjnr, 268
.set nb130nf_innerk, 272
.set nb130nf_n, 276
.set nb130nf_nn1, 280
.set nb130nf_nri, 284
.set nb130nf_facel, 288
.set nb130nf_ntype, 292
.set nb130nf_nouter, 296
.set nb130nf_ninner, 300
.set nb130nf_salign, 304
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $304,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb130nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb130nf_p_nri(%ebp),%ecx
        movl nb130nf_p_facel(%ebp),%esi
        movl nb130nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%esi),%esi
        movl (%edi),%edi
        movl %ecx,nb130nf_nri(%esp)
        movl %esi,nb130nf_facel(%esp)
        movl %edi,nb130nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb130nf_nouter(%esp)
        movl %eax,nb130nf_ninner(%esp)

        movl nb130nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb130nf_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb130nf_half(%esp)
        movss nb130nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb130nf_half(%esp)
        movaps %xmm3,nb130nf_three(%esp)

_nb_kernel130nf_ia32_sse.nb130nf_threadloop: 
        movl  nb130nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel130nf_ia32_sse.nb130nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel130nf_ia32_sse.nb130nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb130nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb130nf_n(%esp)
        movl %ebx,nb130nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel130nf_ia32_sse.nb130nf_outerstart
        jmp _nb_kernel130nf_ia32_sse.nb130nf_end

_nb_kernel130nf_ia32_sse.nb130nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb130nf_nouter(%esp),%ebx
        movl %ebx,nb130nf_nouter(%esp)

_nb_kernel130nf_ia32_sse.nb130nf_outer: 
        movl  nb130nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb130nf_is3(%esp)            ## store is3 

        movl  nb130nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb130nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        movl  nb130nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        mulss nb130nf_facel(%esp),%xmm3
        shufps $0,%xmm3,%xmm3

        movl  nb130nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb130nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb130nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb130nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addss (%eax,%ebx,4),%xmm0
        addss 4(%eax,%ebx,4),%xmm1
        addss 8(%eax,%ebx,4),%xmm2

        movaps %xmm3,nb130nf_iq(%esp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb130nf_ix(%esp)
        movaps %xmm1,nb130nf_iy(%esp)
        movaps %xmm2,nb130nf_iz(%esp)

        movl  %ebx,nb130nf_ii3(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb130nf_vctot(%esp)
        movaps %xmm4,nb130nf_Vvdwtot(%esp)

        movl  nb130nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb130nf_pos(%ebp),%esi
        movl  nb130nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb130nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb130nf_ninner(%esp),%ecx
        movl  %ecx,nb130nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb130nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel130nf_ia32_sse.nb130nf_unroll_loop
        jmp   _nb_kernel130nf_ia32_sse.nb130nf_finish_inner
_nb_kernel130nf_ia32_sse.nb130nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb130nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 
        addl $16,nb130nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb130nf_charge(%ebp),%esi     ## base of charge[] 

        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ecx,4),%xmm4
        movss (%esi,%ebx,4),%xmm6
        movss (%esi,%edx,4),%xmm7

        movaps nb130nf_iq(%esp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movl nb130nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb130nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb130nf_ntia(%esp),%edi
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

        movaps %xmm4,nb130nf_c6(%esp)
        movaps %xmm6,nb130nf_c12(%esp)

        movl nb130nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        mulps %xmm2,%xmm3
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
        movaps nb130nf_ix(%esp),%xmm4
        movaps nb130nf_iy(%esp),%xmm5
        movaps nb130nf_iz(%esp),%xmm6

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

        movaps nb130nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb130nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb130nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv    
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm3
        addps  nb130nf_vctot(%esp),%xmm3
        movaps %xmm3,nb130nf_vctot(%esp)

        ## LJ table
        mulps  %xmm1,%xmm4 ## r
        mulps  nb130nf_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
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

        movl nb130nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movlps (%esi,%ecx,4),%xmm7
        movhps (%esi,%ebx,4),%xmm5
        movhps (%esi,%edx,4),%xmm7 ## got half dispersion table 
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

        movaps nb130nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addps  nb130nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130nf_Vvdwtot(%esp)

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

        movaps nb130nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb130nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $4,nb130nf_innerk(%esp)
        jl    _nb_kernel130nf_ia32_sse.nb130nf_finish_inner
        jmp   _nb_kernel130nf_ia32_sse.nb130nf_unroll_loop
_nb_kernel130nf_ia32_sse.nb130nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb130nf_innerk(%esp)
        movl  nb130nf_innerk(%esp),%edx
        andl  $2,%edx
        jnz   _nb_kernel130nf_ia32_sse.nb130nf_dopair
        jmp   _nb_kernel130nf_ia32_sse.nb130nf_checksingle
_nb_kernel130nf_ia32_sse.nb130nf_dopair: 
        movl nb130nf_charge(%ebp),%esi

    movl  nb130nf_innerjjnr(%esp),%ecx

        movl  (%ecx),%eax
        movl  4(%ecx),%ebx
        addl $8,nb130nf_innerjjnr(%esp)

        xorps %xmm3,%xmm3
        movss (%esi,%eax,4),%xmm3
        movss (%esi,%ebx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## constant 00001100 
        shufps $88,%xmm3,%xmm3 ## constant 01011000 ;# xmm3(0,1) has the charges 

        movl nb130nf_type(%ebp),%esi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%esi,%ecx,4),%ecx
        movl (%esi,%edx,4),%edx
        movl nb130nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb130nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%esi,%ecx,4),%xmm6
        movhps (%esi,%edx,4),%xmm6
        movl nb130nf_pos(%ebp),%edi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## constant 00001000       
        shufps $13,%xmm6,%xmm6 ## constant 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb130nf_c6(%esp)
        movaps %xmm6,nb130nf_c12(%esp)

        leal  (%eax,%eax,2),%eax
        leal  (%ebx,%ebx,2),%ebx
        ## move coordinates to xmm0-xmm2 
        movlps (%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2
        movhps (%edi,%ebx,4),%xmm1
        movss 8(%edi,%ebx,4),%xmm0

        mulps  nb130nf_iq(%esp),%xmm3

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## constant 10001000

        shufps $136,%xmm0,%xmm0 ## constant 10001000
        shufps $221,%xmm1,%xmm1 ## constant 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb130nf_ix(%esp),%xmm4
        movaps nb130nf_iy(%esp),%xmm5
        movaps nb130nf_iz(%esp),%xmm6

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

        movaps nb130nf_krf(%esp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb130nf_three(%esp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb130nf_half(%esp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv    
        movaps %xmm0,%xmm1
        mulps %xmm0,%xmm3
        addps  nb130nf_vctot(%esp),%xmm3
        movaps %xmm3,nb130nf_vctot(%esp)

        ## LJ table
        mulps  %xmm1,%xmm4 ## r
        mulps  nb130nf_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movl nb130nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movhps (%esi,%ebx,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movhps 8(%esi,%ebx,4),%xmm7
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

        movaps nb130nf_c6(%esp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addps  nb130nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movhps 16(%esi,%ebx,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movhps 24(%esi,%ebx,4),%xmm7
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

        movaps nb130nf_c12(%esp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb130nf_Vvdwtot(%esp),%xmm5
        movaps %xmm5,nb130nf_Vvdwtot(%esp)

_nb_kernel130nf_ia32_sse.nb130nf_checksingle:   
        movl  nb130nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel130nf_ia32_sse.nb130nf_dosingle
        jmp    _nb_kernel130nf_ia32_sse.nb130nf_updateouterdata
_nb_kernel130nf_ia32_sse.nb130nf_dosingle: 
        movl nb130nf_charge(%ebp),%esi
        movl nb130nf_pos(%ebp),%edi
        movl  nb130nf_innerjjnr(%esp),%ecx
        xorps %xmm3,%xmm3
        movl  (%ecx),%eax
        movss (%esi,%eax,4),%xmm3       ## xmm3(0) has the charge       

        movl nb130nf_type(%ebp),%esi
        movl %eax,%ecx
        movl (%esi,%ecx,4),%ecx
        movl nb130nf_vdwparam(%ebp),%esi
        shll %ecx
        addl nb130nf_ntia(%esp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%esi,%ecx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## constant 11111100    
        shufps $253,%xmm6,%xmm6 ## constant 11111101    

        movaps %xmm4,nb130nf_c6(%esp)
        movaps %xmm6,nb130nf_c12(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movss (%edi,%eax,4),%xmm0
        movss 4(%edi,%eax,4),%xmm1
        movss 8(%edi,%eax,4),%xmm2

        mulps  nb130nf_iq(%esp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb130nf_ix(%esp),%xmm4
        movaps nb130nf_iy(%esp),%xmm5
        movaps nb130nf_iz(%esp),%xmm6

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

        movss nb130nf_krf(%esp),%xmm7
        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movss %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb130nf_three(%esp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb130nf_half(%esp),%xmm0
        mulss  %xmm4,%xmm7      ## xmm7=krsq 
        subss %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv    
        movaps %xmm0,%xmm1
        mulss %xmm0,%xmm3
        addss  nb130nf_vctot(%esp),%xmm3
        movss %xmm3,nb130nf_vctot(%esp)

        ## LJ table
        mulss  %xmm1,%xmm4 ## r
        mulss  nb130nf_tsc(%esp),%xmm4   ## rtab

        movaps %xmm1,%xmm0 ## copy of rinv
        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movl nb130nf_VFtab(%ebp),%esi
        movd %mm6,%eax

        ## dispersion 
        movlps (%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss  nb130nf_c6(%esp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addss  nb130nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb130nf_Vvdwtot(%esp)

        ## repulsion 
        movlps 16(%esi,%eax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%esi,%eax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss  nb130nf_c12(%esp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb130nf_Vvdwtot(%esp),%xmm5
        movss %xmm5,nb130nf_Vvdwtot(%esp)


_nb_kernel130nf_ia32_sse.nb130nf_updateouterdata: 

        ## get n from stack
        movl nb130nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb130nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb130nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb130nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb130nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb130nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb130nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel130nf_ia32_sse.nb130nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb130nf_n(%esp)
        jmp _nb_kernel130nf_ia32_sse.nb130nf_outer
_nb_kernel130nf_ia32_sse.nb130nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb130nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel130nf_ia32_sse.nb130nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel130nf_ia32_sse.nb130nf_threadloop
_nb_kernel130nf_ia32_sse.nb130nf_end: 
        emms

        movl nb130nf_nouter(%esp),%eax
        movl nb130nf_ninner(%esp),%ebx
        movl nb130nf_outeriter(%ebp),%ecx
        movl nb130nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb130nf_salign(%esp),%eax
        addl %eax,%esp
        addl $304,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






