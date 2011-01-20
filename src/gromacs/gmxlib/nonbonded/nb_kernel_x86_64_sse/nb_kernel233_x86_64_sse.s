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





.globl nb_kernel233_x86_64_sse
.globl _nb_kernel233_x86_64_sse
nb_kernel233_x86_64_sse:        
_nb_kernel233_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb233_fshift, 16
.set nb233_gid, 24
.set nb233_pos, 32
.set nb233_faction, 40
.set nb233_charge, 48
.set nb233_p_facel, 56
.set nb233_argkrf, 64
.set nb233_argcrf, 72
.set nb233_Vc, 80
.set nb233_type, 88
.set nb233_p_ntype, 96
.set nb233_vdwparam, 104
.set nb233_Vvdw, 112
.set nb233_p_tabscale, 120
.set nb233_VFtab, 128
.set nb233_invsqrta, 136
.set nb233_dvda, 144
.set nb233_p_gbtabscale, 152
.set nb233_GBtab, 160
.set nb233_p_nthreads, 168
.set nb233_count, 176
.set nb233_mtx, 184
.set nb233_outeriter, 192
.set nb233_inneriter, 200
.set nb233_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb233_ixO, 0
.set nb233_iyO, 16
.set nb233_izO, 32
.set nb233_ixH1, 48
.set nb233_iyH1, 64
.set nb233_izH1, 80
.set nb233_ixH2, 96
.set nb233_iyH2, 112
.set nb233_izH2, 128
.set nb233_ixM, 144
.set nb233_iyM, 160
.set nb233_izM, 176
.set nb233_iqM, 192
.set nb233_iqH, 208
.set nb233_dxO, 224
.set nb233_dyO, 240
.set nb233_dzO, 256
.set nb233_dxH1, 272
.set nb233_dyH1, 288
.set nb233_dzH1, 304
.set nb233_dxH2, 320
.set nb233_dyH2, 336
.set nb233_dzH2, 352
.set nb233_dxM, 368
.set nb233_dyM, 384
.set nb233_dzM, 400
.set nb233_qqM, 416
.set nb233_qqH, 432
.set nb233_rinvH1, 448
.set nb233_rinvH2, 464
.set nb233_rinvM, 480
.set nb233_two, 496
.set nb233_c6, 512
.set nb233_c12, 528
.set nb233_tsc, 544
.set nb233_fstmp, 560
.set nb233_krf, 576
.set nb233_crf, 592
.set nb233_krsqH1, 608
.set nb233_krsqH2, 624
.set nb233_krsqM, 640
.set nb233_vctot, 656
.set nb233_Vvdwtot, 672
.set nb233_fixO, 688
.set nb233_fiyO, 704
.set nb233_fizO, 720
.set nb233_fixH1, 736
.set nb233_fiyH1, 752
.set nb233_fizH1, 768
.set nb233_fixH2, 784
.set nb233_fiyH2, 800
.set nb233_fizH2, 816
.set nb233_fixM, 832
.set nb233_fiyM, 848
.set nb233_fizM, 864
.set nb233_fjx, 880
.set nb233_fjy, 896
.set nb233_fjz, 912
.set nb233_half, 928
.set nb233_three, 944
.set nb233_rsqOO, 960
.set nb233_facel, 976
.set nb233_iinr, 984
.set nb233_jindex, 992
.set nb233_jjnr, 1000
.set nb233_shift, 1008
.set nb233_shiftvec, 1016
.set nb233_innerjjnr, 1032
.set nb233_is3, 1040
.set nb233_ii3, 1044
.set nb233_nri, 1048
.set nb233_ntia, 1052
.set nb233_innerk, 1056
.set nb233_n, 1060
.set nb233_nn1, 1064
.set nb233_nouter, 1068
.set nb233_ninner, 1072

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1096,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb233_nouter(%rsp)
        movl %eax,nb233_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb233_nri(%rsp)
        movq %rsi,nb233_iinr(%rsp)
        movq %rdx,nb233_jindex(%rsp)
        movq %rcx,nb233_jjnr(%rsp)
        movq %r8,nb233_shift(%rsp)
        movq %r9,nb233_shiftvec(%rsp)
        movq nb233_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb233_facel(%rsp)

        movq nb233_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb233_tsc(%rsp)

        movq nb233_argkrf(%rbp),%rsi
        movq nb233_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb233_krf(%rsp)
        movaps %xmm2,nb233_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb233_half(%rsp)
        movss nb233_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb233_half(%rsp)
        movaps %xmm2,nb233_two(%rsp)
        movaps %xmm3,nb233_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb233_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb233_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb233_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb233_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb233_iqM(%rsp)
        movaps %xmm4,nb233_iqH(%rsp)

        movq  nb233_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb233_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233_ntia(%rsp)
_nb_kernel233_x86_64_sse.nb233_threadloop: 
        movq  nb233_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel233_x86_64_sse.nb233_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233_x86_64_sse.nb233_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233_n(%rsp)
        movl %ebx,nb233_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233_x86_64_sse.nb233_outerstart
        jmp _nb_kernel233_x86_64_sse.nb233_end

_nb_kernel233_x86_64_sse.nb233_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233_nouter(%rsp),%ebx
        movl %ebx,nb233_nouter(%rsp)

_nb_kernel233_x86_64_sse.nb233_outer: 
        movq  nb233_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb233_is3(%rsp)      ## store is3 

        movq  nb233_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb233_iinr(%rsp),%rcx             ## ecx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb233_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb233_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb233_ixO(%rsp)
        movaps %xmm4,nb233_iyO(%rsp)
        movaps %xmm5,nb233_izO(%rsp)
        movaps %xmm6,nb233_ixH1(%rsp)
        movaps %xmm7,nb233_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb233_izH1(%rsp)
        movaps %xmm0,nb233_ixH2(%rsp)
        movaps %xmm1,nb233_iyH2(%rsp)
        movaps %xmm2,nb233_izH2(%rsp)
        movaps %xmm3,nb233_ixM(%rsp)
        movaps %xmm4,nb233_iyM(%rsp)
        movaps %xmm5,nb233_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb233_vctot(%rsp)
        movaps %xmm4,nb233_Vvdwtot(%rsp)
        movaps %xmm4,nb233_fixO(%rsp)
        movaps %xmm4,nb233_fiyO(%rsp)
        movaps %xmm4,nb233_fizO(%rsp)
        movaps %xmm4,nb233_fixH1(%rsp)
        movaps %xmm4,nb233_fiyH1(%rsp)
        movaps %xmm4,nb233_fizH1(%rsp)
        movaps %xmm4,nb233_fixH2(%rsp)
        movaps %xmm4,nb233_fiyH2(%rsp)
        movaps %xmm4,nb233_fizH2(%rsp)
        movaps %xmm4,nb233_fixM(%rsp)
        movaps %xmm4,nb233_fiyM(%rsp)
        movaps %xmm4,nb233_fizM(%rsp)

        movq  nb233_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb233_pos(%rbp),%rsi
        movq  nb233_faction(%rbp),%rdi
        movq  nb233_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb233_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb233_ninner(%rsp),%ecx
        movl  %ecx,nb233_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb233_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel233_x86_64_sse.nb233_unroll_loop
        jmp   _nb_kernel233_x86_64_sse.nb233_odd_inner
_nb_kernel233_x86_64_sse.nb233_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb233_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb233_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb233_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb233_iqM(%rsp),%xmm3
        mulps  nb233_iqH(%rsp),%xmm4

        movaps  %xmm3,nb233_qqM(%rsp)
        movaps  %xmm4,nb233_qqH(%rsp)

        movq nb233_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb233_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb233_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        movlps (%rsi,%r8,4),%xmm6
        movlps (%rsi,%r10,4),%xmm7
        movhps (%rsi,%r9,4),%xmm6
        movhps (%rsi,%r11,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movaps %xmm4,nb233_c6(%rsp)
        movaps %xmm6,nb233_c12(%rsp)

        movq nb233_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    ## O interaction
    ## copy to xmm3-xmm5
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5

    subps nb233_ixO(%rsp),%xmm3
    subps nb233_iyO(%rsp),%xmm4
    subps nb233_izO(%rsp),%xmm5

    movaps %xmm3,nb233_dxO(%rsp)
    movaps %xmm4,nb233_dyO(%rsp)
    movaps %xmm5,nb233_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    ## xmm3=rsq

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm3,%xmm5
        movaps %xmm5,%xmm15
        mulps %xmm5,%xmm5
        movaps nb233_three(%rsp),%xmm4
        mulps %xmm3,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm15,%xmm4
        mulps nb233_half(%rsp),%xmm4
        movaps %xmm4,%xmm15
        mulps  %xmm4,%xmm3
    ## xmm15=rinv
    ## xmm3=r

    mulps nb233_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm3,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 8
    pslld   $3,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm3   ## xmm3=eps

    ## move to integer registers
    movhlps %xmm5,%xmm6
    movd    %xmm5,%r8d
    movd    %xmm6,%r10d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm6,%xmm6
    movd    %xmm5,%r9d
    movd    %xmm6,%r11d
    ## xmm3=eps
    ## xmm15=rinv

        movq nb233_VFtab(%rbp),%rsi
    ## calculate LJ table
    movlps (%rsi,%r8,4),%xmm5
        movlps 16(%rsi,%r8,4),%xmm9

        movlps (%rsi,%r10,4),%xmm7
        movlps 16(%rsi,%r10,4),%xmm11

        movhps (%rsi,%r9,4),%xmm5
        movhps 16(%rsi,%r9,4),%xmm9

        movhps (%rsi,%r11,4),%xmm7
        movhps 16(%rsi,%r11,4),%xmm11

    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 8(%rsi,%r8,4),%xmm7
        movlps 24(%rsi,%r8,4),%xmm11

        movlps 8(%rsi,%r10,4),%xmm13
        movlps 24(%rsi,%r10,4),%xmm14

        movhps 8(%rsi,%r9,4),%xmm7
        movhps 24(%rsi,%r9,4),%xmm11

        movhps 8(%rsi,%r11,4),%xmm13
        movhps 24(%rsi,%r11,4),%xmm14

    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm13,%xmm6 ## 10001000
        shufps $136,%xmm14,%xmm10 ## 10001000
        shufps $221,%xmm13,%xmm7 ## 11011101
        shufps $221,%xmm14,%xmm11 ## 11011101
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulps  %xmm3,%xmm7   ## Heps
    mulps  %xmm3,%xmm11
    mulps  %xmm3,%xmm6  ## Geps
    mulps  %xmm3,%xmm10
    mulps  %xmm3,%xmm7  ## Heps2
    mulps  %xmm3,%xmm11
    addps  %xmm6,%xmm5 ## F+Geps
    addps  %xmm10,%xmm9
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm11,%xmm9
    addps  %xmm7,%xmm7   ## 2*Heps2
    addps  %xmm11,%xmm11
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm10,%xmm11

    addps  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm9,%xmm11
    mulps  %xmm3,%xmm5 ## eps*Fp
    mulps  %xmm3,%xmm9
    movaps nb233_c6(%rsp),%xmm12
    movaps nb233_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb233_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb233_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7

    mulps  nb233_tsc(%rsp),%xmm7
    mulps  %xmm15,%xmm7  ## -fscal
    xorps  %xmm9,%xmm9

    subps  %xmm7,%xmm9    ## fscal
    movaps %xmm9,%xmm10
    movaps %xmm9,%xmm11

    mulps  nb233_dxO(%rsp),%xmm9    ## fx/fy/fz
    mulps  nb233_dyO(%rsp),%xmm10
    mulps  nb233_dzO(%rsp),%xmm11

    ## save j force temporarily
    movaps %xmm9,nb233_fjx(%rsp)
    movaps %xmm10,nb233_fjy(%rsp)
    movaps %xmm11,nb233_fjz(%rsp)

    ## increment i O force
    addps nb233_fixO(%rsp),%xmm9
    addps nb233_fiyO(%rsp),%xmm10
    addps nb233_fizO(%rsp),%xmm11
    movaps %xmm9,nb233_fixO(%rsp)
    movaps %xmm10,nb233_fiyO(%rsp)
    movaps %xmm11,nb233_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.                
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb233_ixH1(%rsp),%xmm0
    subps nb233_iyH1(%rsp),%xmm1
    subps nb233_izH1(%rsp),%xmm2
    subps nb233_ixH2(%rsp),%xmm3
    subps nb233_iyH2(%rsp),%xmm4
    subps nb233_izH2(%rsp),%xmm5
    subps nb233_ixM(%rsp),%xmm6
    subps nb233_iyM(%rsp),%xmm7
    subps nb233_izM(%rsp),%xmm8

        movaps %xmm0,nb233_dxH1(%rsp)
        movaps %xmm1,nb233_dyH1(%rsp)
        movaps %xmm2,nb233_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb233_dxH2(%rsp)
        movaps %xmm4,nb233_dyH2(%rsp)
        movaps %xmm5,nb233_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb233_dxM(%rsp)
        movaps %xmm7,nb233_dyM(%rsp)
        movaps %xmm8,nb233_dzM(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb233_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  nb233_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1 
        mulps   %xmm4,%xmm10 ## rinvH2
    mulps   %xmm4,%xmm11 ## rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb233_krf(%rsp),%xmm2
    mulps  %xmm9,%xmm9  ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  %xmm2,%xmm0 ## k*rsq
    mulps  %xmm2,%xmm3
    mulps  %xmm2,%xmm6
    movaps %xmm0,%xmm2 ## copy of k*rsq
    movaps %xmm3,%xmm5
    movaps %xmm6,%xmm8
    addps  %xmm1,%xmm2 ## rinv+krsq
    addps  %xmm4,%xmm5
    addps  %xmm7,%xmm8
    movaps nb233_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb233_qqH(%rsp),%xmm12
    movaps nb233_qqM(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm12,%xmm4
    mulps  %xmm13,%xmm7
    addps  nb233_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,%xmm15

    movaps %xmm2,nb233_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j forces to local temp variables 
        movq nb233_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm9 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm10 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm9 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm10 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm11    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm12    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm6     ## jzb
    movss  8(%rdi,%rdx,4),%xmm8     ## jzd
    movlhps %xmm6,%xmm11 ## jza  -  jzb  -
    movlhps %xmm8,%xmm12 ## jzc  -  jzd -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jza jzb jzc jzd

    ## xmm9: jxa jya jxb jyb 
    ## xmm10: jxc jyc jxd jyd
    ## xmm11: jza jzb jzc jzd

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb233_dxH1(%rsp),%xmm0
        mulps nb233_dyH1(%rsp),%xmm1
        mulps nb233_dzH1(%rsp),%xmm2
        mulps nb233_dxH2(%rsp),%xmm3
        mulps nb233_dyH2(%rsp),%xmm4
        mulps nb233_dzH2(%rsp),%xmm5
        mulps nb233_dxM(%rsp),%xmm6
        mulps nb233_dyM(%rsp),%xmm7
        mulps nb233_dzM(%rsp),%xmm8

    ## fetch forces from O interaction
    movaps nb233_fjx(%rsp),%xmm13
    movaps nb233_fjy(%rsp),%xmm14
    addps  nb233_fjz(%rsp),%xmm11

    addps %xmm0,%xmm13
    addps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb233_fixH1(%rsp),%xmm0
    addps nb233_fiyH1(%rsp),%xmm1
    addps nb233_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb233_fixH2(%rsp),%xmm3
    addps nb233_fiyH2(%rsp),%xmm4
    addps nb233_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb233_fixM(%rsp),%xmm6
    addps nb233_fiyM(%rsp),%xmm7
    addps nb233_fizM(%rsp),%xmm8

    movaps %xmm0,nb233_fixH1(%rsp)
    movaps %xmm1,nb233_fiyH1(%rsp)
    movaps %xmm2,nb233_fizH1(%rsp)
    movaps %xmm3,nb233_fixH2(%rsp)
    movaps %xmm4,nb233_fiyH2(%rsp)
    movaps %xmm5,nb233_fizH2(%rsp)
    movaps %xmm6,nb233_fixM(%rsp)
    movaps %xmm7,nb233_fiyM(%rsp)
    movaps %xmm8,nb233_fizM(%rsp)

    ## xmm9 = fjx
    ## xmm10 = fjy
    ## xmm11 = fjz
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12 ## fjzc fjzd

    movlps %xmm9,(%rdi,%rax,4)
    movhps %xmm9,(%rdi,%rbx,4)
    movlps %xmm10,(%rdi,%rcx,4)
    movhps %xmm10,(%rdi,%rdx,4)
    movss  %xmm11,8(%rdi,%rax,4)
    movss  %xmm12,8(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,8(%rdi,%rbx,4)
    movss  %xmm12,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb233_innerk(%rsp)
        jl    _nb_kernel233_x86_64_sse.nb233_odd_inner
        jmp   _nb_kernel233_x86_64_sse.nb233_unroll_loop
_nb_kernel233_x86_64_sse.nb233_odd_inner: 
        addl $4,nb233_innerk(%rsp)
        jnz   _nb_kernel233_x86_64_sse.nb233_odd_loop
        jmp   _nb_kernel233_x86_64_sse.nb233_updateouterdata
_nb_kernel233_x86_64_sse.nb233_odd_loop: 
        movq  nb233_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb233_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb233_iqM(%rsp),%xmm4
        movq nb233_charge(%rbp),%rsi
        movhps nb233_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb233_qqM(%rsp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb233_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb233_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb233_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb233_c6(%rsp)
        movaps %xmm7,nb233_c12(%rsp)

        movq nb233_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb233_ixO(%rsp),%xmm0
        movss nb233_iyO(%rsp),%xmm1
        movss nb233_izO(%rsp),%xmm2
        movss nb233_ixH1(%rsp),%xmm3
        movss nb233_iyH1(%rsp),%xmm4
        movss nb233_izH1(%rsp),%xmm5
        unpcklps nb233_ixH2(%rsp),%xmm0         ## ixO ixH2 - -
        unpcklps nb233_iyH2(%rsp),%xmm1         ## iyO iyH2 - -
        unpcklps nb233_izH2(%rsp),%xmm2         ## izO izH2 - -
        unpcklps nb233_ixM(%rsp),%xmm3          ## ixH1 ixM - -
        unpcklps nb233_iyM(%rsp),%xmm4          ## iyH1 iyM - -
        unpcklps nb233_izM(%rsp),%xmm5          ## izH1 izM - -
        unpcklps %xmm3,%xmm0    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm4,%xmm1    ## same for y
        unpcklps %xmm5,%xmm2    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use O distances for storage
        movaps %xmm3,nb233_dxO(%rsp)
        movaps %xmm4,nb233_dyO(%rsp)
        movaps %xmm5,nb233_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 
        movaps %xmm4,%xmm0
        mulps nb233_krf(%rsp),%xmm0
        movaps %xmm0,nb233_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb233_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb233_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb233_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movd %eax,%mm0

        movq nb233_VFtab(%rbp),%rsi
        movd %mm6,%eax

        ## dispersion 
        movlps (%rsi,%rax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb233_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb233_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb233_tsc(%rsp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb233_Vvdwtot(%rsp),%xmm5
        movss %xmm7,nb233_fstmp(%rsp)
        movss %xmm5,nb233_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%rax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb233_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb233_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb233_tsc(%rsp),%xmm7
        addss  nb233_fstmp(%rsp),%xmm7
        movss %xmm7,nb233_fstmp(%rsp)
        addss  nb233_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb233_Vvdwtot(%rsp)

        movd %mm0,%eax

        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb233_krsqM(%rsp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb233_crf(%rsp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb233_two(%rsp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq
        movaps %xmm5,%xmm2
        mulps  nb233_qqM(%rsp),%xmm2    ## xmm2=vcoul 
        mulps  nb233_qqM(%rsp),%xmm1    ## xmm1=coul part of fs 
        mulps  %xmm0,%xmm1
        subss  nb233_fstmp(%rsp),%xmm1
        mulps  %xmm0,%xmm1
        movaps %xmm1,%xmm4

        addps  nb233_vctot(%rsp),%xmm2
        movaps %xmm2,nb233_vctot(%rsp)


        movaps nb233_dxO(%rsp),%xmm0
        movaps nb233_dyO(%rsp),%xmm1
        movaps nb233_dzO(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb233_fixO(%rsp),%xmm3
        movss  nb233_fiyO(%rsp),%xmm4
        movss  nb233_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb233_fixO(%rsp)
        movss  %xmm4,nb233_fiyO(%rsp)
        movss  %xmm5,nb233_fizO(%rsp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixH1(%rsp),%xmm3
        addss  nb233_fiyH1(%rsp),%xmm4
        addss  nb233_fizH1(%rsp),%xmm5
        movss  %xmm3,nb233_fixH1(%rsp)
        movss  %xmm4,nb233_fiyH1(%rsp)
        movss  %xmm5,nb233_fizH1(%rsp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixH2(%rsp),%xmm3
        addss  nb233_fiyH2(%rsp),%xmm4
        addss  nb233_fizH2(%rsp),%xmm5
        movss  %xmm3,nb233_fixH2(%rsp)
        movss  %xmm4,nb233_fiyH2(%rsp)
        movss  %xmm5,nb233_fizH2(%rsp)          ## updated the H2 force 

        movq nb233_faction(%rbp),%rdi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb233_fixM(%rsp),%xmm3
        addss  nb233_fiyM(%rsp),%xmm4
        addss  nb233_fizM(%rsp),%xmm5
        movss  %xmm3,nb233_fixM(%rsp)
        movss  %xmm4,nb233_fiyM(%rsp)
        movss  %xmm5,nb233_fizM(%rsp)   ## updated the M force 

        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,%xmm0
        movaps  %xmm4,%xmm1
        movaps  %xmm5,%xmm2

        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2
        unpcklps %xmm1,%xmm0    ## x,y sum in xmm0, z sum in xmm2

        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb233_innerk(%rsp)
        jz    _nb_kernel233_x86_64_sse.nb233_updateouterdata
        jmp   _nb_kernel233_x86_64_sse.nb233_odd_loop
_nb_kernel233_x86_64_sse.nb233_updateouterdata: 
        movl  nb233_ii3(%rsp),%ecx
        movq  nb233_faction(%rbp),%rdi
        movq  nb233_fshift(%rbp),%rsi
        movl  nb233_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb233_fixO(%rsp),%xmm0
        movaps nb233_fiyO(%rsp),%xmm1
        movaps nb233_fizO(%rsp),%xmm2

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
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb233_fixH1(%rsp),%xmm0
        movaps nb233_fiyH1(%rsp),%xmm1
        movaps nb233_fizH1(%rsp),%xmm2

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
        movss  12(%rdi,%rcx,4),%xmm3
        movss  16(%rdi,%rcx,4),%xmm4
        movss  20(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,12(%rdi,%rcx,4)
        movss  %xmm4,16(%rdi,%rcx,4)
        movss  %xmm5,20(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb233_fixH2(%rsp),%xmm0
        movaps nb233_fiyH2(%rsp),%xmm1
        movaps nb233_fizH2(%rsp),%xmm2

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
        movss  24(%rdi,%rcx,4),%xmm3
        movss  28(%rdi,%rcx,4),%xmm4
        movss  32(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,24(%rdi,%rcx,4)
        movss  %xmm4,28(%rdi,%rcx,4)
        movss  %xmm5,32(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb233_fixM(%rsp),%xmm0
        movaps nb233_fiyM(%rsp),%xmm1
        movaps nb233_fizM(%rsp),%xmm2

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
        movss  36(%rdi,%rcx,4),%xmm3
        movss  40(%rdi,%rcx,4),%xmm4
        movss  44(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,36(%rdi,%rcx,4)
        movss  %xmm4,40(%rdi,%rcx,4)
        movss  %xmm5,44(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb233_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb233_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb233_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb233_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb233_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb233_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb233_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233_x86_64_sse.nb233_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233_n(%rsp)
        jmp _nb_kernel233_x86_64_sse.nb233_outer
_nb_kernel233_x86_64_sse.nb233_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233_x86_64_sse.nb233_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233_x86_64_sse.nb233_threadloop
_nb_kernel233_x86_64_sse.nb233_end: 
        movl nb233_nouter(%rsp),%eax
        movl nb233_ninner(%rsp),%ebx
        movq nb233_outeriter(%rbp),%rcx
        movq nb233_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1096,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel233nf_x86_64_sse
.globl _nb_kernel233nf_x86_64_sse
nb_kernel233nf_x86_64_sse:      
_nb_kernel233nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb233nf_fshift, 16
.set nb233nf_gid, 24
.set nb233nf_pos, 32
.set nb233nf_faction, 40
.set nb233nf_charge, 48
.set nb233nf_p_facel, 56
.set nb233nf_argkrf, 64
.set nb233nf_argcrf, 72
.set nb233nf_Vc, 80
.set nb233nf_type, 88
.set nb233nf_p_ntype, 96
.set nb233nf_vdwparam, 104
.set nb233nf_Vvdw, 112
.set nb233nf_p_tabscale, 120
.set nb233nf_VFtab, 128
.set nb233nf_invsqrta, 136
.set nb233nf_dvda, 144
.set nb233nf_p_gbtabscale, 152
.set nb233nf_GBtab, 160
.set nb233nf_p_nthreads, 168
.set nb233nf_count, 176
.set nb233nf_mtx, 184
.set nb233nf_outeriter, 192
.set nb233nf_inneriter, 200
.set nb233nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb233nf_ixO, 0
.set nb233nf_iyO, 16
.set nb233nf_izO, 32
.set nb233nf_ixH1, 48
.set nb233nf_iyH1, 64
.set nb233nf_izH1, 80
.set nb233nf_ixH2, 96
.set nb233nf_iyH2, 112
.set nb233nf_izH2, 128
.set nb233nf_ixM, 144
.set nb233nf_iyM, 160
.set nb233nf_izM, 176
.set nb233nf_iqM, 192
.set nb233nf_iqH, 208
.set nb233nf_qqM, 224
.set nb233nf_qqH, 240
.set nb233nf_rinvH1, 256
.set nb233nf_rinvH2, 272
.set nb233nf_rinvM, 288
.set nb233nf_tsc, 304
.set nb233nf_c6, 320
.set nb233nf_c12, 336
.set nb233nf_krf, 352
.set nb233nf_crf, 368
.set nb233nf_krsqH1, 384
.set nb233nf_krsqH2, 400
.set nb233nf_krsqM, 416
.set nb233nf_vctot, 432
.set nb233nf_Vvdwtot, 448
.set nb233nf_half, 464
.set nb233nf_three, 480
.set nb233nf_rsqO, 496
.set nb233nf_facel, 512
.set nb233nf_iinr, 520
.set nb233nf_jindex, 528
.set nb233nf_jjnr, 536
.set nb233nf_shift, 544
.set nb233nf_shiftvec, 552
.set nb233nf_innerjjnr, 560
.set nb233nf_nri, 568
.set nb233nf_is3, 572
.set nb233nf_ii3, 576
.set nb233nf_ntia, 580
.set nb233nf_innerk, 584
.set nb233nf_n, 588
.set nb233nf_nn1, 592
.set nb233nf_nouter, 596
.set nb233nf_ninner, 600


        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $616,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb233nf_nouter(%rsp)
        movl %eax,nb233nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb233nf_nri(%rsp)
        movq %rsi,nb233nf_iinr(%rsp)
        movq %rdx,nb233nf_jindex(%rsp)
        movq %rcx,nb233nf_jjnr(%rsp)
        movq %r8,nb233nf_shift(%rsp)
        movq %r9,nb233nf_shiftvec(%rsp)
        movq nb233nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb233nf_facel(%rsp)

        movq nb233nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb233nf_tsc(%rsp)

        movq nb233nf_argkrf(%rbp),%rsi
        movq nb233nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb233nf_krf(%rsp)
        movaps %xmm2,nb233nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb233nf_half(%rsp)
        movss nb233nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb233nf_half(%rsp)
        movaps %xmm3,nb233nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb233nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb233nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb233nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb233nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb233nf_iqM(%rsp)
        movaps %xmm4,nb233nf_iqH(%rsp)

        movq  nb233nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb233nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb233nf_ntia(%rsp)

_nb_kernel233nf_x86_64_sse.nb233nf_threadloop: 
        movq  nb233nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel233nf_x86_64_sse.nb233nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel233nf_x86_64_sse.nb233nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb233nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb233nf_n(%rsp)
        movl %ebx,nb233nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel233nf_x86_64_sse.nb233nf_outerstart
        jmp _nb_kernel233nf_x86_64_sse.nb233nf_end

_nb_kernel233nf_x86_64_sse.nb233nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb233nf_nouter(%rsp),%ebx
        movl %ebx,nb233nf_nouter(%rsp)

_nb_kernel233nf_x86_64_sse.nb233nf_outer: 
        movq  nb233nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb233nf_is3(%rsp)            ## store is3 

        movq  nb233nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb233nf_iinr(%rsp),%rcx           ## ecx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb233nf_pos(%rbp),%rax    ## eax = base of pos[]  
        movl  %ebx,nb233nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb233nf_ixO(%rsp)
        movaps %xmm4,nb233nf_iyO(%rsp)
        movaps %xmm5,nb233nf_izO(%rsp)
        movaps %xmm6,nb233nf_ixH1(%rsp)
        movaps %xmm7,nb233nf_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb233nf_izH1(%rsp)
        movaps %xmm0,nb233nf_ixH2(%rsp)
        movaps %xmm1,nb233nf_iyH2(%rsp)
        movaps %xmm2,nb233nf_izH2(%rsp)
        movaps %xmm3,nb233nf_ixM(%rsp)
        movaps %xmm4,nb233nf_iyM(%rsp)
        movaps %xmm5,nb233nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb233nf_vctot(%rsp)
        movaps %xmm4,nb233nf_Vvdwtot(%rsp)

        movq  nb233nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb233nf_pos(%rbp),%rsi
        movq  nb233nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb233nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb233nf_ninner(%rsp),%ecx
        movl  %ecx,nb233nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb233nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel233nf_x86_64_sse.nb233nf_unroll_loop
        jmp   _nb_kernel233nf_x86_64_sse.nb233nf_odd_inner
_nb_kernel233nf_x86_64_sse.nb233nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb233nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb233nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb233nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb233nf_iqM(%rsp),%xmm3
        mulps  nb233nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb233nf_qqM(%rsp)
        movaps  %xmm4,nb233nf_qqH(%rsp)

        movq nb233nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb233nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb233nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm6 ## constant 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb233nf_c6(%rsp)
        movaps %xmm6,nb233nf_c12(%rsp)

        movq nb233nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## constant 10001000

        shufps $136,%xmm5,%xmm0 ## constant 10001000
        shufps $221,%xmm5,%xmm1 ## constant 11011101            

        ## move ixO-izO to xmm4-xmm6 
        movaps nb233nf_ixO(%rsp),%xmm4
        movaps nb233nf_iyO(%rsp),%xmm5
        movaps nb233nf_izO(%rsp),%xmm6

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
        movaps %xmm4,%xmm7
        ## rsqO in xmm7

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb233nf_ixH1(%rsp),%xmm4
        movaps nb233nf_iyH1(%rsp),%xmm5
        movaps nb233nf_izH1(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb233nf_ixH2(%rsp),%xmm3
        movaps nb233nf_iyH2(%rsp),%xmm4
        movaps nb233nf_izH2(%rsp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movaps nb233nf_iyM(%rsp),%xmm3
        movaps nb233nf_izM(%rsp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb233nf_ixM(%rsp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movaps %xmm4,%xmm0
        movaps %xmm5,%xmm1
        movaps %xmm6,%xmm2
        mulps  nb233nf_krf(%rsp),%xmm0
        mulps  nb233nf_krf(%rsp),%xmm1
        mulps  nb233nf_krf(%rsp),%xmm2
        movaps %xmm0,nb233nf_krsqM(%rsp)
        movaps %xmm1,nb233nf_krsqH2(%rsp)
        movaps %xmm2,nb233nf_krsqH1(%rsp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%rsp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%rsp),%xmm0
        movaps  %xmm0,nb233nf_rinvH1(%rsp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%rsp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%rsp),%xmm0
        movaps  %xmm0,nb233nf_rinvH2(%rsp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%rsp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%rsp),%xmm0
        movaps  %xmm0,nb233nf_rinvM(%rsp)

        ## Do the O LJ-only interaction directly.       
        ## rsqO is in xmm7
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb233nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb233nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm0
        ## xmm0=rinvO

        mulps %xmm0,%xmm7
        mulps nb233nf_tsc(%rsp),%xmm7   ## rtab

        movhlps %xmm7,%xmm5
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6
        pslld $3,%mm7

        movq nb233nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## dispersion 
        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half dispersion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of dispersion table 
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

        movaps nb233nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb233nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb233nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%rax,4),%xmm5
        movlps 16(%rsi,%rcx,4),%xmm7
        movhps 16(%rsi,%rbx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm7    ## got half repulsion table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movlps 24(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rbx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm3    ## other half of repulsion table 
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

        movaps nb233nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb233nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb233nf_Vvdwtot(%rsp)

        ## Do H1 interaction
        movaps  nb233nf_rinvH1(%rsp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqH1(%rsp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqH(%rsp),%xmm0

        addps  nb233nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb233nf_vctot(%rsp)

        ## Done with H1, do H2 interactions
        movaps  nb233nf_rinvH2(%rsp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqH2(%rsp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqH(%rsp),%xmm0

        addps  nb233nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb233nf_vctot(%rsp)

        ## Done with H2, do M interactions
        movaps  nb233nf_rinvM(%rsp),%xmm7
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq
        movaps %xmm7,%xmm0
        movaps nb233nf_krsqM(%rsp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb233nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqM(%rsp),%xmm0

        addps  nb233nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb233nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb233nf_innerk(%rsp)
        jl    _nb_kernel233nf_x86_64_sse.nb233nf_odd_inner
        jmp   _nb_kernel233nf_x86_64_sse.nb233nf_unroll_loop
_nb_kernel233nf_x86_64_sse.nb233nf_odd_inner: 
        addl $4,nb233nf_innerk(%rsp)
        jnz   _nb_kernel233nf_x86_64_sse.nb233nf_odd_loop
        jmp   _nb_kernel233nf_x86_64_sse.nb233nf_updateouterdata
_nb_kernel233nf_x86_64_sse.nb233nf_odd_loop: 
        movq  nb233nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb233nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb233nf_iqM(%rsp),%xmm4
        movq nb233nf_charge(%rbp),%rsi
        movhps nb233nf_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb233nf_qqM(%rsp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb233nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb233nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb233nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb233nf_c6(%rsp)
        movaps %xmm7,nb233nf_c12(%rsp)

        movq nb233nf_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb233nf_ixO(%rsp),%xmm3
        movss nb233nf_iyO(%rsp),%xmm4
        movss nb233nf_izO(%rsp),%xmm5
        movss nb233nf_ixH1(%rsp),%xmm0
        movss nb233nf_iyH1(%rsp),%xmm1
        movss nb233nf_izH1(%rsp),%xmm2
        unpcklps nb233nf_ixH2(%rsp),%xmm3       ## ixO ixH2 - -
        unpcklps nb233nf_iyH2(%rsp),%xmm4       ## iyO iyH2 - -
        unpcklps nb233nf_izH2(%rsp),%xmm5       ## izO izH2 - -
        unpcklps nb233nf_ixM(%rsp),%xmm0        ## ixH1 ixM - -
        unpcklps nb233nf_iyM(%rsp),%xmm1        ## iyH1 iyM - -
        unpcklps nb233nf_izM(%rsp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 
        movaps %xmm4,%xmm0
        mulps nb233nf_krf(%rsp),%xmm0
        movaps %xmm0,nb233nf_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb233nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb233nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv, xmm4=rsq

        mulps %xmm0,%xmm4
        mulps  nb233nf_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb233nf_VFtab(%rbp),%rsi
        movd %mm6,%eax

        ## dispersion 
        movlps (%rsi,%rax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%rsi,%rax,4),%xmm7
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

        movss nb233nf_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb233nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb233nf_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%rax,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%rsi,%rax,4),%xmm7
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

        movss nb233nf_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb233nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb233nf_Vvdwtot(%rsp)

        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb233nf_krsqM(%rsp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb233nf_crf(%rsp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb233nf_qqM(%rsp),%xmm5          ## xmm2=vcoul 

        addps  nb233nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb233nf_vctot(%rsp)

        decl nb233nf_innerk(%rsp)
        jz    _nb_kernel233nf_x86_64_sse.nb233nf_updateouterdata
        jmp   _nb_kernel233nf_x86_64_sse.nb233nf_odd_loop
_nb_kernel233nf_x86_64_sse.nb233nf_updateouterdata: 
        ## get n from stack
        movl nb233nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb233nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb233nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb233nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb233nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb233nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb233nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel233nf_x86_64_sse.nb233nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb233nf_n(%rsp)
        jmp _nb_kernel233nf_x86_64_sse.nb233nf_outer
_nb_kernel233nf_x86_64_sse.nb233nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb233nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel233nf_x86_64_sse.nb233nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel233nf_x86_64_sse.nb233nf_threadloop
_nb_kernel233nf_x86_64_sse.nb233nf_end: 

        movl nb233nf_nouter(%rsp),%eax
        movl nb233nf_ninner(%rsp),%ebx
        movq nb233nf_outeriter(%rbp),%rcx
        movq nb233nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $616,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




