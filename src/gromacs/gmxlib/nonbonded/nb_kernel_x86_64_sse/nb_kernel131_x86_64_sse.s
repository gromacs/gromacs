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





.globl nb_kernel131_x86_64_sse
.globl _nb_kernel131_x86_64_sse
nb_kernel131_x86_64_sse:        
_nb_kernel131_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb131_fshift, 16
.set nb131_gid, 24
.set nb131_pos, 32
.set nb131_faction, 40
.set nb131_charge, 48
.set nb131_p_facel, 56
.set nb131_argkrf, 64
.set nb131_argcrf, 72
.set nb131_Vc, 80
.set nb131_type, 88
.set nb131_p_ntype, 96
.set nb131_vdwparam, 104
.set nb131_Vvdw, 112
.set nb131_p_tabscale, 120
.set nb131_VFtab, 128
.set nb131_invsqrta, 136
.set nb131_dvda, 144
.set nb131_p_gbtabscale, 152
.set nb131_GBtab, 160
.set nb131_p_nthreads, 168
.set nb131_count, 176
.set nb131_mtx, 184
.set nb131_outeriter, 192
.set nb131_inneriter, 200
.set nb131_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb131_ixO, 0
.set nb131_iyO, 16
.set nb131_izO, 32
.set nb131_ixH1, 48
.set nb131_iyH1, 64
.set nb131_izH1, 80
.set nb131_ixH2, 96
.set nb131_iyH2, 112
.set nb131_izH2, 128
.set nb131_iqO, 144
.set nb131_iqH, 160
.set nb131_dxO, 176
.set nb131_dyO, 192
.set nb131_dzO, 208
.set nb131_dxH1, 224
.set nb131_dyH1, 240
.set nb131_dzH1, 256
.set nb131_dxH2, 272
.set nb131_dyH2, 288
.set nb131_dzH2, 304
.set nb131_qqO, 320
.set nb131_qqH, 336
.set nb131_c6, 352
.set nb131_c12, 368
.set nb131_tsc, 384
.set nb131_fstmp, 400
.set nb131_vctot, 416
.set nb131_Vvdwtot, 432
.set nb131_fixO, 448
.set nb131_fiyO, 464
.set nb131_fizO, 480
.set nb131_fixH1, 496
.set nb131_fiyH1, 512
.set nb131_fizH1, 528
.set nb131_fixH2, 544
.set nb131_fiyH2, 560
.set nb131_fizH2, 576
.set nb131_fjx, 592
.set nb131_fjy, 608
.set nb131_fjz, 624
.set nb131_half, 640
.set nb131_three, 656
.set nb131_two, 672
.set nb131_krf, 688
.set nb131_crf, 704
.set nb131_rsqO, 720
.set nb131_rsqH1, 736
.set nb131_rsqH2, 752
.set nb131_rinvO, 768
.set nb131_rinvH1, 784
.set nb131_rinvH2, 800
.set nb131_facel, 816
.set nb131_iinr, 824
.set nb131_jindex, 832
.set nb131_jjnr, 840
.set nb131_shift, 848
.set nb131_shiftvec, 856
.set nb131_innerjjnr, 864
.set nb131_nri, 872
.set nb131_is3, 876
.set nb131_ii3, 880
.set nb131_ntia, 884
.set nb131_innerk, 888
.set nb131_n, 892
.set nb131_nn1, 896
.set nb131_nouter, 900
.set nb131_ninner, 904

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $920,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb131_nouter(%rsp)
        movl %eax,nb131_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb131_nri(%rsp)
        movq %rsi,nb131_iinr(%rsp)
        movq %rdx,nb131_jindex(%rsp)
        movq %rcx,nb131_jjnr(%rsp)
        movq %r8,nb131_shift(%rsp)
        movq %r9,nb131_shiftvec(%rsp)
        movq nb131_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb131_facel(%rsp)

        movq nb131_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb131_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb131_half(%rsp)
        movss nb131_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb131_half(%rsp)
        movaps %xmm2,nb131_two(%rsp)
        movaps %xmm3,nb131_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb131_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb131_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb131_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb131_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb131_iqO(%rsp)
        movaps %xmm4,nb131_iqH(%rsp)

        movq  nb131_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb131_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131_ntia(%rsp)

_nb_kernel131_x86_64_sse.nb131_threadloop: 
        movq  nb131_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel131_x86_64_sse.nb131_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131_x86_64_sse.nb131_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131_n(%rsp)
        movl %ebx,nb131_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131_x86_64_sse.nb131_outerstart
        jmp _nb_kernel131_x86_64_sse.nb131_end

_nb_kernel131_x86_64_sse.nb131_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131_nouter(%rsp),%ebx
        movl %ebx,nb131_nouter(%rsp)

_nb_kernel131_x86_64_sse.nb131_outer: 
        movq  nb131_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb131_is3(%rsp)      ## store is3 

        movq  nb131_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb131_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb131_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb131_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb131_ixO(%rsp)
        movaps %xmm4,nb131_iyO(%rsp)
        movaps %xmm5,nb131_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb131_ixH1(%rsp)
        movaps %xmm1,nb131_iyH1(%rsp)
        movaps %xmm2,nb131_izH1(%rsp)
        movaps %xmm3,nb131_ixH2(%rsp)
        movaps %xmm4,nb131_iyH2(%rsp)
        movaps %xmm5,nb131_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb131_vctot(%rsp)
        movaps %xmm4,nb131_Vvdwtot(%rsp)
        movaps %xmm4,nb131_fixO(%rsp)
        movaps %xmm4,nb131_fiyO(%rsp)
        movaps %xmm4,nb131_fizO(%rsp)
        movaps %xmm4,nb131_fixH1(%rsp)
        movaps %xmm4,nb131_fiyH1(%rsp)
        movaps %xmm4,nb131_fizH1(%rsp)
        movaps %xmm4,nb131_fixH2(%rsp)
        movaps %xmm4,nb131_fiyH2(%rsp)
        movaps %xmm4,nb131_fizH2(%rsp)

        movq  nb131_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb131_pos(%rbp),%rsi
        movq  nb131_faction(%rbp),%rdi
        movq  nb131_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb131_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb131_ninner(%rsp),%ecx
        movl  %ecx,nb131_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb131_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel131_x86_64_sse.nb131_unroll_loop
        jmp   _nb_kernel131_x86_64_sse.nb131_odd_inner
_nb_kernel131_x86_64_sse.nb131_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb131_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb131_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb131_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb131_iqO(%rsp),%xmm3
        mulps  nb131_iqH(%rsp),%xmm4

        movaps  %xmm3,nb131_qqO(%rsp)
        movaps  %xmm4,nb131_qqH(%rsp)

        movq nb131_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb131_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb131_ntia(%rsp),%edi
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

        movaps %xmm4,nb131_c6(%rsp)
        movaps %xmm6,nb131_c12(%rsp)

        movq nb131_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four j coordinates to xmm0-xmm2         
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

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb131_ixO(%rsp),%xmm0
    subps nb131_iyO(%rsp),%xmm1
    subps nb131_izO(%rsp),%xmm2
    subps nb131_ixH1(%rsp),%xmm3
    subps nb131_iyH1(%rsp),%xmm4
    subps nb131_izH1(%rsp),%xmm5
    subps nb131_ixH2(%rsp),%xmm6
    subps nb131_iyH2(%rsp),%xmm7
    subps nb131_izH2(%rsp),%xmm8

        movaps %xmm0,nb131_dxO(%rsp)
        movaps %xmm1,nb131_dyO(%rsp)
        movaps %xmm2,nb131_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb131_dxH1(%rsp)
        movaps %xmm4,nb131_dyH1(%rsp)
        movaps %xmm5,nb131_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb131_dxH2(%rsp)
        movaps %xmm7,nb131_dyH2(%rsp)
        movaps %xmm8,nb131_dzH2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb131_three(%rsp),%xmm9
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

        movaps  nb131_half(%rsp),%xmm2
        mulps   %xmm2,%xmm9 ## rinvO 
        mulps   %xmm2,%xmm10 ## rinvH1
    mulps   %xmm2,%xmm11 ## rinvH2

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movaps %xmm0,nb131_rsqO(%rsp)
    movaps %xmm3,nb131_rsqH1(%rsp)
    movaps %xmm6,nb131_rsqH2(%rsp)
    movaps %xmm9,nb131_rinvO(%rsp)
    movaps %xmm10,nb131_rinvH1(%rsp)
    movaps %xmm11,nb131_rinvH2(%rsp)

    ## table LJ interaction
    mulps  %xmm9,%xmm0
    mulps  nb131_tsc(%rsp),%xmm0   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm0,%xmm1

    ## convert back to float
    cvtdq2ps  %xmm1,%xmm2

    ## multiply by 8
    pslld   $3,%xmm1

    ## move to integer registers
    movhlps %xmm1,%xmm13
    movd    %xmm1,%r8d
    movd    %xmm13,%r10d
    shufps $1,%xmm1,%xmm1
    shufps $1,%xmm13,%xmm13
    movd    %xmm1,%r9d
    movd    %xmm13,%r11d

    ## calculate eps
    subps     %xmm2,%xmm0
    movq nb131_VFtab(%rbp),%rsi

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

    mulps  %xmm0,%xmm7   ## Heps
    mulps  %xmm0,%xmm11
    mulps  %xmm0,%xmm6  ## Geps
    mulps  %xmm0,%xmm10
    mulps  %xmm0,%xmm7  ## Heps2
    mulps  %xmm0,%xmm11
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
    mulps  %xmm0,%xmm5 ## eps*Fp
    mulps  %xmm0,%xmm9
    movaps nb131_c6(%rsp),%xmm12
    movaps nb131_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb131_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb131_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7
    mulps  nb131_tsc(%rsp),%xmm7

    movaps nb131_rinvO(%rsp),%xmm9
    movaps nb131_rinvH1(%rsp),%xmm10
    movaps nb131_rinvH2(%rsp),%xmm11
    movaps %xmm9,%xmm0  ## rinv
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2

    mulps  %xmm10,%xmm10 ## rinvsq
    mulps  %xmm11,%xmm11
    mulps  nb131_qqO(%rsp),%xmm0
    mulps  nb131_qqH(%rsp),%xmm1
    mulps  nb131_qqH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    subps  %xmm7,%xmm9
    mulps  nb131_rinvO(%rsp),%xmm9

    addps nb131_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb131_vctot(%rsp)

        ## move j  forces to local temp variables 
        movq nb131_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm0 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm2    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm4    ## jzb
    movss  8(%rdi,%rdx,4),%xmm5    ## jzd
    movlhps %xmm4,%xmm2 ## jza  -  jzb  -
    movlhps %xmm5,%xmm3 ## jzc  -  jzd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jza jzb jzc jzd

    ## xmm0: jxa jya jxb jyb 
    ## xmm1: jxc jyc jxd jyd
    ## xmm2: jza jzb jzc jzd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb131_dxO(%rsp),%xmm7
        mulps nb131_dyO(%rsp),%xmm8
        mulps nb131_dzO(%rsp),%xmm9
        mulps nb131_dxH1(%rsp),%xmm10
        mulps nb131_dyH1(%rsp),%xmm11
        mulps nb131_dzH1(%rsp),%xmm12
        mulps nb131_dxH2(%rsp),%xmm13
        mulps nb131_dyH2(%rsp),%xmm14
        mulps nb131_dzH2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb131_fixO(%rsp),%xmm7
    addps nb131_fiyO(%rsp),%xmm8
    addps nb131_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb131_fixH1(%rsp),%xmm10
    addps nb131_fiyH1(%rsp),%xmm11
    addps nb131_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb131_fixH2(%rsp),%xmm13
    addps nb131_fiyH2(%rsp),%xmm14
    addps nb131_fizH2(%rsp),%xmm15

    movaps %xmm7,nb131_fixO(%rsp)
    movaps %xmm8,nb131_fiyO(%rsp)
    movaps %xmm9,nb131_fizO(%rsp)
    movaps %xmm10,nb131_fixH1(%rsp)
    movaps %xmm11,nb131_fiyH1(%rsp)
    movaps %xmm12,nb131_fizH1(%rsp)
    movaps %xmm13,nb131_fixH2(%rsp)
    movaps %xmm14,nb131_fiyH2(%rsp)
    movaps %xmm15,nb131_fizH2(%rsp)

    ## xmm0 = fjx
    ## xmm1 = fjy
    ## xmm2 = fjz
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3

    movlps %xmm0,(%rdi,%rax,4)
    movhps %xmm0,(%rdi,%rbx,4)
    movlps %xmm1,(%rdi,%rcx,4)
    movhps %xmm1,(%rdi,%rdx,4)
    movss  %xmm2,8(%rdi,%rax,4)
    movss  %xmm3,8(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,8(%rdi,%rbx,4)
    movss  %xmm3,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb131_innerk(%rsp)
        jl    _nb_kernel131_x86_64_sse.nb131_odd_inner
        jmp   _nb_kernel131_x86_64_sse.nb131_unroll_loop
_nb_kernel131_x86_64_sse.nb131_odd_inner: 
        addl $4,nb131_innerk(%rsp)
        jnz   _nb_kernel131_x86_64_sse.nb131_odd_loop
        jmp   _nb_kernel131_x86_64_sse.nb131_updateouterdata
_nb_kernel131_x86_64_sse.nb131_odd_loop: 
        movq  nb131_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb131_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb131_iqO(%rsp),%xmm4
        movq nb131_charge(%rbp),%rsi
        movhps nb131_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb131_qqO(%rsp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb131_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb131_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb131_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb131_c6(%rsp)
        movaps %xmm7,nb131_c12(%rsp)

        movq nb131_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb131_ixO(%rsp),%xmm0
        movss nb131_iyO(%rsp),%xmm1
        movss nb131_izO(%rsp),%xmm2

        movlps nb131_ixH1(%rsp),%xmm6
        movlps nb131_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb131_iyH1(%rsp),%xmm6
        movlps nb131_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb131_izH1(%rsp),%xmm6
        movlps nb131_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb131_dxO(%rsp)
        movaps %xmm4,nb131_dyO(%rsp)
        movaps %xmm5,nb131_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb131_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb131_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        ## rsq still in xmm4, rinv in xmm0.
        mulps  %xmm0,%xmm4
        mulps  nb131_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb131_VFtab(%rbp),%rsi
        movd %mm6,%r12d

        ## dispersion 
        movlps (%rsi,%r12,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%rsi,%r12,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb131_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb131_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb131_tsc(%rsp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb131_Vvdwtot(%rsp),%xmm5
        movss %xmm7,nb131_fstmp(%rsp)
        movss %xmm5,nb131_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%r12,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%rsi,%r12,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb131_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb131_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb131_tsc(%rsp),%xmm7
        addss  nb131_fstmp(%rsp),%xmm7
        movss %xmm7,nb131_fstmp(%rsp)

        addss  nb131_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb131_Vvdwtot(%rsp)

        ## coulomb in analytical form. xmm0=rinv
        movaps %xmm0,%xmm2              ## rinv
        mulps  nb131_qqO(%rsp),%xmm2    ## vcoulO
        movaps %xmm2,%xmm3              ## vcoulO

        addps  nb131_vctot(%rsp),%xmm2
        movaps %xmm2,nb131_vctot(%rsp)

        mulps  %xmm0,%xmm3
        subss  nb131_fstmp(%rsp),%xmm3
        mulps  %xmm0,%xmm3

        movaps nb131_dxO(%rsp),%xmm0
        movaps nb131_dyO(%rsp),%xmm1
        movaps nb131_dzO(%rsp),%xmm2

        mulps  %xmm3,%xmm0
        mulps  %xmm3,%xmm1
        mulps  %xmm3,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb131_fixO(%rsp),%xmm3
        movss  nb131_fiyO(%rsp),%xmm4
        movss  nb131_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb131_fixO(%rsp)
        movss  %xmm4,nb131_fiyO(%rsp)
        movss  %xmm5,nb131_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb131_fixH1(%rsp),%xmm3
        addss  nb131_fiyH1(%rsp),%xmm4
        addss  nb131_fizH1(%rsp),%xmm5
        movss  %xmm3,nb131_fixH1(%rsp)
        movss  %xmm4,nb131_fiyH1(%rsp)
        movss  %xmm5,nb131_fizH1(%rsp)          ## updated the H1 force 

        movq nb131_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb131_fixH2(%rsp),%xmm3
        addss  nb131_fiyH2(%rsp),%xmm4
        addss  nb131_fizH2(%rsp),%xmm5
        movss  %xmm3,nb131_fixH2(%rsp)
        movss  %xmm4,nb131_fiyH2(%rsp)
        movss  %xmm5,nb131_fizH2(%rsp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2   ## z sum in xmm2 
        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb131_innerk(%rsp)
        jz    _nb_kernel131_x86_64_sse.nb131_updateouterdata
        jmp   _nb_kernel131_x86_64_sse.nb131_odd_loop
_nb_kernel131_x86_64_sse.nb131_updateouterdata: 
        movl  nb131_ii3(%rsp),%ecx
        movq  nb131_faction(%rbp),%rdi
        movq  nb131_fshift(%rbp),%rsi
        movl  nb131_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb131_fixO(%rsp),%xmm0
        movaps nb131_fiyO(%rsp),%xmm1
        movaps nb131_fizO(%rsp),%xmm2

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
        movaps nb131_fixH1(%rsp),%xmm0
        movaps nb131_fiyH1(%rsp),%xmm1
        movaps nb131_fizH1(%rsp),%xmm2

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
        movaps nb131_fixH2(%rsp),%xmm0
        movaps nb131_fiyH2(%rsp),%xmm1
        movaps nb131_fizH2(%rsp),%xmm2

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

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb131_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb131_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb131_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb131_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb131_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb131_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb131_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131_x86_64_sse.nb131_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131_n(%rsp)
        jmp _nb_kernel131_x86_64_sse.nb131_outer
_nb_kernel131_x86_64_sse.nb131_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131_x86_64_sse.nb131_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131_x86_64_sse.nb131_threadloop
_nb_kernel131_x86_64_sse.nb131_end: 
        movl nb131_nouter(%rsp),%eax
        movl nb131_ninner(%rsp),%ebx
        movq nb131_outeriter(%rbp),%rcx
        movq nb131_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $920,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel131nf_x86_64_sse
.globl _nb_kernel131nf_x86_64_sse
nb_kernel131nf_x86_64_sse:      
_nb_kernel131nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb131nf_fshift, 16
.set nb131nf_gid, 24
.set nb131nf_pos, 32
.set nb131nf_faction, 40
.set nb131nf_charge, 48
.set nb131nf_p_facel, 56
.set nb131nf_argkrf, 64
.set nb131nf_argcrf, 72
.set nb131nf_Vc, 80
.set nb131nf_type, 88
.set nb131nf_p_ntype, 96
.set nb131nf_vdwparam, 104
.set nb131nf_Vvdw, 112
.set nb131nf_p_tabscale, 120
.set nb131nf_VFtab, 128
.set nb131nf_invsqrta, 136
.set nb131nf_dvda, 144
.set nb131nf_p_gbtabscale, 152
.set nb131nf_GBtab, 160
.set nb131nf_p_nthreads, 168
.set nb131nf_count, 176
.set nb131nf_mtx, 184
.set nb131nf_outeriter, 192
.set nb131nf_inneriter, 200
.set nb131nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb131nf_ixO, 0
.set nb131nf_iyO, 16
.set nb131nf_izO, 32
.set nb131nf_ixH1, 48
.set nb131nf_iyH1, 64
.set nb131nf_izH1, 80
.set nb131nf_ixH2, 96
.set nb131nf_iyH2, 112
.set nb131nf_izH2, 128
.set nb131nf_iqO, 144
.set nb131nf_iqH, 160
.set nb131nf_qqO, 176
.set nb131nf_qqH, 192
.set nb131nf_c6, 208
.set nb131nf_c12, 224
.set nb131nf_vctot, 240
.set nb131nf_Vvdwtot, 256
.set nb131nf_half, 272
.set nb131nf_three, 288
.set nb131nf_rinvO, 304
.set nb131nf_rinvH1, 320
.set nb131nf_rinvH2, 336
.set nb131nf_krsqO, 352
.set nb131nf_krsqH1, 368
.set nb131nf_krsqH2, 384
.set nb131nf_tsc, 400
.set nb131nf_facel, 416
.set nb131nf_iinr, 424
.set nb131nf_jindex, 432
.set nb131nf_jjnr, 440
.set nb131nf_shift, 448
.set nb131nf_shiftvec, 456
.set nb131nf_innerjjnr, 464
.set nb131nf_nri, 472
.set nb131nf_ntia, 476
.set nb131nf_is3, 480
.set nb131nf_ii3, 484
.set nb131nf_innerk, 488
.set nb131nf_n, 492
.set nb131nf_nn1, 496
.set nb131nf_nouter, 500
.set nb131nf_ninner, 504

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $520,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb131nf_nouter(%rsp)
        movl %eax,nb131nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb131nf_nri(%rsp)
        movq %rsi,nb131nf_iinr(%rsp)
        movq %rdx,nb131nf_jindex(%rsp)
        movq %rcx,nb131nf_jjnr(%rsp)
        movq %r8,nb131nf_shift(%rsp)
        movq %r9,nb131nf_shiftvec(%rsp)
        movq nb131nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb131nf_facel(%rsp)

        movq nb131nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb131nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb131nf_half(%rsp)
        movss nb131nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb131nf_half(%rsp)
        movaps %xmm3,nb131nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb131nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb131nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb131nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb131nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb131nf_iqO(%rsp)
        movaps %xmm4,nb131nf_iqH(%rsp)

        movq  nb131nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb131nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb131nf_ntia(%rsp)

_nb_kernel131nf_x86_64_sse.nb131nf_threadloop: 
        movq  nb131nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel131nf_x86_64_sse.nb131nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel131nf_x86_64_sse.nb131nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb131nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb131nf_n(%rsp)
        movl %ebx,nb131nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel131nf_x86_64_sse.nb131nf_outerstart
        jmp _nb_kernel131nf_x86_64_sse.nb131nf_end

_nb_kernel131nf_x86_64_sse.nb131nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb131nf_nouter(%rsp),%ebx
        movl %ebx,nb131nf_nouter(%rsp)

_nb_kernel131nf_x86_64_sse.nb131nf_outer: 
        movq  nb131nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb131nf_is3(%rsp)            ## store is3 

        movq  nb131nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb131nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb131nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb131nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb131nf_ixO(%rsp)
        movaps %xmm4,nb131nf_iyO(%rsp)
        movaps %xmm5,nb131nf_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb131nf_ixH1(%rsp)
        movaps %xmm1,nb131nf_iyH1(%rsp)
        movaps %xmm2,nb131nf_izH1(%rsp)
        movaps %xmm3,nb131nf_ixH2(%rsp)
        movaps %xmm4,nb131nf_iyH2(%rsp)
        movaps %xmm5,nb131nf_izH2(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb131nf_vctot(%rsp)
        movaps %xmm4,nb131nf_Vvdwtot(%rsp)

        movq  nb131nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb131nf_pos(%rbp),%rsi
        movq  nb131nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb131nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb131nf_ninner(%rsp),%ecx
        movl  %ecx,nb131nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb131nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel131nf_x86_64_sse.nb131nf_unroll_loop
        jmp   _nb_kernel131nf_x86_64_sse.nb131nf_odd_inner
_nb_kernel131nf_x86_64_sse.nb131nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb131nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb131nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb131nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb131nf_iqO(%rsp),%xmm3
        mulps  nb131nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb131nf_qqO(%rsp)
        movaps  %xmm4,nb131nf_qqH(%rsp)

        movq nb131nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb131nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb131nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb131nf_c6(%rsp)
        movaps %xmm6,nb131nf_c12(%rsp)

        movq nb131nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
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
        movaps nb131nf_ixO(%rsp),%xmm4
        movaps nb131nf_iyO(%rsp),%xmm5
        movaps nb131nf_izO(%rsp),%xmm6

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
        movaps nb131nf_ixH1(%rsp),%xmm4
        movaps nb131nf_iyH1(%rsp),%xmm5
        movaps nb131nf_izH1(%rsp),%xmm6

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
        movaps nb131nf_ixH2(%rsp),%xmm3
        movaps nb131nf_iyH2(%rsp),%xmm4
        movaps nb131nf_izH2(%rsp),%xmm5

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
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

        ## start with rsqO (still in xmm7) - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb131nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb131nf_half(%rsp),%xmm4
        movaps  %xmm4,nb131nf_rinvO(%rsp)       ## rinvO in xmm0, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb131nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb131nf_half(%rsp),%xmm4
        movaps  %xmm4,nb131nf_rinvH1(%rsp)      ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb131nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb131nf_half(%rsp),%xmm4
        movaps  %xmm4,nb131nf_rinvH2(%rsp)      ## rinvH2 in xmm5 


        ## do O table interactions - rsqO in xmm7.
        mulps nb131nf_rinvO(%rsp),%xmm7
        mulps nb131nf_tsc(%rsp),%xmm7   ## rtab

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

        movq nb131nf_VFtab(%rbp),%rsi
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

        movaps nb131nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb131nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb131nf_Vvdwtot(%rsp)

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

        movaps nb131nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb131nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb131nf_Vvdwtot(%rsp)

        movaps nb131nf_rinvO(%rsp),%xmm2
        mulps  nb131nf_qqO(%rsp),%xmm2

        ## H1 & H2 interactions 
        movaps nb131nf_rinvH1(%rsp),%xmm3
        addps  nb131nf_rinvH2(%rsp),%xmm3

        mulps  nb131nf_qqH(%rsp),%xmm3

        addps  %xmm3,%xmm2
        addps  nb131nf_vctot(%rsp),%xmm2
        movaps %xmm2,nb131nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb131nf_innerk(%rsp)
        jl    _nb_kernel131nf_x86_64_sse.nb131nf_odd_inner
        jmp   _nb_kernel131nf_x86_64_sse.nb131nf_unroll_loop
_nb_kernel131nf_x86_64_sse.nb131nf_odd_inner: 
        addl $4,nb131nf_innerk(%rsp)
        jnz   _nb_kernel131nf_x86_64_sse.nb131nf_odd_loop
        jmp   _nb_kernel131nf_x86_64_sse.nb131nf_updateouterdata
_nb_kernel131nf_x86_64_sse.nb131nf_odd_loop: 
        movq  nb131nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb131nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb131nf_iqO(%rsp),%xmm4
        movq nb131nf_charge(%rbp),%rsi
        movhps nb131nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb131nf_qqO(%rsp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb131nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb131nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb131nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb131nf_c6(%rsp)
        movaps %xmm7,nb131nf_c12(%rsp)

        movq nb131nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb131nf_ixO(%rsp),%xmm3
        movss nb131nf_iyO(%rsp),%xmm4
        movss nb131nf_izO(%rsp),%xmm5

        movlps nb131nf_ixH1(%rsp),%xmm6
        movlps nb131nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb131nf_iyH1(%rsp),%xmm6
        movlps nb131nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb131nf_izH1(%rsp),%xmm6
        movlps nb131nf_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb131nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb131nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## constant 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## constant 11100000

        ## rsq still in xmm4, rinv in xmm0.
        mulps  %xmm0,%xmm4
        mulps  nb131nf_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb131nf_VFtab(%rbp),%rsi
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

        movss nb131nf_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        addss  nb131nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb131nf_Vvdwtot(%rsp)

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

        movss nb131nf_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  nb131nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb131nf_Vvdwtot(%rsp)

        ## xmm0=rinv
        mulps  nb131nf_qqO(%rsp),%xmm0          ## xmm0=vcoul 

        addps  nb131nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb131nf_vctot(%rsp)

        decl nb131nf_innerk(%rsp)
        jz    _nb_kernel131nf_x86_64_sse.nb131nf_updateouterdata
        jmp   _nb_kernel131nf_x86_64_sse.nb131nf_odd_loop
_nb_kernel131nf_x86_64_sse.nb131nf_updateouterdata: 
        ## get n from stack
        movl nb131nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb131nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb131nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb131nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb131nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb131nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb131nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel131nf_x86_64_sse.nb131nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb131nf_n(%rsp)
        jmp _nb_kernel131nf_x86_64_sse.nb131nf_outer
_nb_kernel131nf_x86_64_sse.nb131nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb131nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel131nf_x86_64_sse.nb131nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel131nf_x86_64_sse.nb131nf_threadloop
_nb_kernel131nf_x86_64_sse.nb131nf_end: 
        movl nb131nf_nouter(%rsp),%eax
        movl nb131nf_ninner(%rsp),%ebx
        movq nb131nf_outeriter(%rbp),%rcx
        movq nb131nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $520,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



