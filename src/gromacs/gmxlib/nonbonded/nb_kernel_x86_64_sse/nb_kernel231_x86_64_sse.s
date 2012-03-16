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






.globl nb_kernel231_x86_64_sse
.globl _nb_kernel231_x86_64_sse
nb_kernel231_x86_64_sse:        
_nb_kernel231_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb231_fshift, 16
.set nb231_gid, 24
.set nb231_pos, 32
.set nb231_faction, 40
.set nb231_charge, 48
.set nb231_p_facel, 56
.set nb231_argkrf, 64
.set nb231_argcrf, 72
.set nb231_Vc, 80
.set nb231_type, 88
.set nb231_p_ntype, 96
.set nb231_vdwparam, 104
.set nb231_Vvdw, 112
.set nb231_p_tabscale, 120
.set nb231_VFtab, 128
.set nb231_invsqrta, 136
.set nb231_dvda, 144
.set nb231_p_gbtabscale, 152
.set nb231_GBtab, 160
.set nb231_p_nthreads, 168
.set nb231_count, 176
.set nb231_mtx, 184
.set nb231_outeriter, 192
.set nb231_inneriter, 200
.set nb231_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb231_ixO, 0
.set nb231_iyO, 16
.set nb231_izO, 32
.set nb231_ixH1, 48
.set nb231_iyH1, 64
.set nb231_izH1, 80
.set nb231_ixH2, 96
.set nb231_iyH2, 112
.set nb231_izH2, 128
.set nb231_iqO, 144
.set nb231_iqH, 160
.set nb231_dxO, 176
.set nb231_dyO, 192
.set nb231_dzO, 208
.set nb231_dxH1, 224
.set nb231_dyH1, 240
.set nb231_dzH1, 256
.set nb231_dxH2, 272
.set nb231_dyH2, 288
.set nb231_dzH2, 304
.set nb231_qqO, 320
.set nb231_qqH, 336
.set nb231_c6, 352
.set nb231_c12, 368
.set nb231_tsc, 384
.set nb231_fstmp, 400
.set nb231_vctot, 416
.set nb231_Vvdwtot, 432
.set nb231_fixO, 448
.set nb231_fiyO, 464
.set nb231_fizO, 480
.set nb231_fixH1, 496
.set nb231_fiyH1, 512
.set nb231_fizH1, 528
.set nb231_fixH2, 544
.set nb231_fiyH2, 560
.set nb231_fizH2, 576
.set nb231_fjx, 592
.set nb231_fjy, 608
.set nb231_fjz, 624
.set nb231_half, 640
.set nb231_three, 656
.set nb231_two, 672
.set nb231_krf, 688
.set nb231_crf, 704
.set nb231_rsqO, 720
.set nb231_rsqH1, 736
.set nb231_rsqH2, 752
.set nb231_rinvO, 768
.set nb231_rinvH1, 784
.set nb231_rinvH2, 800
.set nb231_facel, 816
.set nb231_iinr, 824
.set nb231_jindex, 832
.set nb231_jjnr, 840
.set nb231_shift, 848
.set nb231_shiftvec, 856
.set nb231_innerjjnr, 864
.set nb231_nri, 872
.set nb231_is3, 876
.set nb231_ii3, 880
.set nb231_ntia, 884
.set nb231_innerk, 888
.set nb231_n, 892
.set nb231_nn1, 896
.set nb231_nouter, 900
.set nb231_ninner, 904

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $920,%rsp  # # local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb231_nouter(%rsp)
        movl %eax,nb231_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb231_nri(%rsp)
        movq %rsi,nb231_iinr(%rsp)
        movq %rdx,nb231_jindex(%rsp)
        movq %rcx,nb231_jjnr(%rsp)
        movq %r8,nb231_shift(%rsp)
        movq %r9,nb231_shiftvec(%rsp)
        movq nb231_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb231_facel(%rsp)

        movq nb231_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb231_tsc(%rsp)

        movq nb231_argkrf(%rbp),%rsi
        movq nb231_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb231_krf(%rsp)
        movaps %xmm2,nb231_crf(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb231_half(%rsp)
        movss nb231_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb231_half(%rsp)
        movaps %xmm2,nb231_two(%rsp)
        movaps %xmm3,nb231_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb231_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb231_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb231_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb231_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb231_iqO(%rsp)
        movaps %xmm4,nb231_iqH(%rsp)

        movq  nb231_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb231_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231_ntia(%rsp)

_nb_kernel231_x86_64_sse.nb231_threadloop: 
        movq  nb231_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel231_x86_64_sse.nb231_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231_x86_64_sse.nb231_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231_n(%rsp)
        movl %ebx,nb231_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231_x86_64_sse.nb231_outerstart
        jmp _nb_kernel231_x86_64_sse.nb231_end

_nb_kernel231_x86_64_sse.nb231_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231_nouter(%rsp),%ebx
        movl %ebx,nb231_nouter(%rsp)

_nb_kernel231_x86_64_sse.nb231_outer: 
        movq  nb231_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb231_is3(%rsp)      ## store is3 

        movq  nb231_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb231_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb231_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb231_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb231_ixO(%rsp)
        movaps %xmm4,nb231_iyO(%rsp)
        movaps %xmm5,nb231_izO(%rsp)

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
        movaps %xmm0,nb231_ixH1(%rsp)
        movaps %xmm1,nb231_iyH1(%rsp)
        movaps %xmm2,nb231_izH1(%rsp)
        movaps %xmm3,nb231_ixH2(%rsp)
        movaps %xmm4,nb231_iyH2(%rsp)
        movaps %xmm5,nb231_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb231_vctot(%rsp)
        movaps %xmm4,nb231_Vvdwtot(%rsp)
        movaps %xmm4,nb231_fixO(%rsp)
        movaps %xmm4,nb231_fiyO(%rsp)
        movaps %xmm4,nb231_fizO(%rsp)
        movaps %xmm4,nb231_fixH1(%rsp)
        movaps %xmm4,nb231_fiyH1(%rsp)
        movaps %xmm4,nb231_fizH1(%rsp)
        movaps %xmm4,nb231_fixH2(%rsp)
        movaps %xmm4,nb231_fiyH2(%rsp)
        movaps %xmm4,nb231_fizH2(%rsp)

        movq  nb231_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb231_pos(%rbp),%rsi
        movq  nb231_faction(%rbp),%rdi
        movq  nb231_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb231_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb231_ninner(%rsp),%ecx
        movl  %ecx,nb231_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb231_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel231_x86_64_sse.nb231_unroll_loop
        jmp   _nb_kernel231_x86_64_sse.nb231_odd_inner
_nb_kernel231_x86_64_sse.nb231_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb231_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb231_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb231_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb231_iqO(%rsp),%xmm3
        mulps  nb231_iqH(%rsp),%xmm4

        movaps  %xmm3,nb231_qqO(%rsp)
        movaps  %xmm4,nb231_qqH(%rsp)

        movq nb231_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb231_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb231_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        movlps (%rsi,%r8,4),%xmm6
        movlps (%rsi,%r10,4),%xmm7
        movhps (%rsi,%r9,4),%xmm6
        movhps (%rsi,%r11,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movaps %xmm4,nb231_c6(%rsp)
        movaps %xmm6,nb231_c12(%rsp)

        movq nb231_pos(%rbp),%rsi        ## base of pos[] 

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

    subps nb231_ixO(%rsp),%xmm0
    subps nb231_iyO(%rsp),%xmm1
    subps nb231_izO(%rsp),%xmm2
    subps nb231_ixH1(%rsp),%xmm3
    subps nb231_iyH1(%rsp),%xmm4
    subps nb231_izH1(%rsp),%xmm5
    subps nb231_ixH2(%rsp),%xmm6
    subps nb231_iyH2(%rsp),%xmm7
    subps nb231_izH2(%rsp),%xmm8

        movaps %xmm0,nb231_dxO(%rsp)
        movaps %xmm1,nb231_dyO(%rsp)
        movaps %xmm2,nb231_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb231_dxH1(%rsp)
        movaps %xmm4,nb231_dyH1(%rsp)
        movaps %xmm5,nb231_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb231_dxH2(%rsp)
        movaps %xmm7,nb231_dyH2(%rsp)
        movaps %xmm8,nb231_dzH2(%rsp)
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

        movaps  nb231_three(%rsp),%xmm9
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

        movaps  nb231_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvO
        mulps   %xmm4,%xmm10 ## rinvH1
    mulps   %xmm4,%xmm11 ## rinvH2

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movaps %xmm0,nb231_rsqO(%rsp)
    movaps %xmm3,nb231_rsqH1(%rsp)
    movaps %xmm6,nb231_rsqH2(%rsp)
    movaps %xmm9,nb231_rinvO(%rsp)
    movaps %xmm10,nb231_rinvH1(%rsp)
    movaps %xmm11,nb231_rinvH2(%rsp)

    ## table LJ interaction
    mulps  %xmm9,%xmm0
    mulps  nb231_tsc(%rsp),%xmm0   ## rtab

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
    movq nb231_VFtab(%rbp),%rsi

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
    movaps nb231_c6(%rsp),%xmm12
    movaps nb231_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb231_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb231_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7
    mulps  nb231_tsc(%rsp),%xmm7
    movaps %xmm7,nb231_fstmp(%rsp)

    ## Coulomb reaction-field interaction
    movaps nb231_rsqO(%rsp),%xmm0
    movaps nb231_rsqH1(%rsp),%xmm3
    movaps nb231_rsqH2(%rsp),%xmm6
    movaps nb231_rinvO(%rsp),%xmm9
    movaps nb231_rinvH1(%rsp),%xmm10
    movaps nb231_rinvH2(%rsp),%xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb231_krf(%rsp),%xmm2
    mulps  %xmm10,%xmm10 ## rinvsq
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
    subps  nb231_crf(%rsp),%xmm2     ## rinv+krsq-crf
    subps  nb231_crf(%rsp),%xmm5
    subps  nb231_crf(%rsp),%xmm8
    mulps  nb231_qqO(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf)
    mulps  nb231_qqH(%rsp),%xmm5   ## voul=qq*(rinv+ krsq-crf)
    mulps  nb231_qqH(%rsp),%xmm8   ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  nb231_qqO(%rsp),%xmm1     ## (rinv-2*krsq)*qq
    mulps  nb231_qqH(%rsp),%xmm4
    mulps  nb231_qqH(%rsp),%xmm7
    addps  nb231_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,%xmm15
    mulps  %xmm9,%xmm1  ## fijC
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7
    ## xmm1, xmm4, xmm7 contains fscal coul
    ## xmm15 contains vctot
    subps nb231_fstmp(%rsp),%xmm1
    mulps %xmm9,%xmm1

    movaps %xmm2,nb231_vctot(%rsp)

        ## move j  forces to local temp variables 
        movq nb231_faction(%rbp),%rdi
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

        mulps nb231_dxO(%rsp),%xmm0
        mulps nb231_dyO(%rsp),%xmm1
        mulps nb231_dzO(%rsp),%xmm2
        mulps nb231_dxH1(%rsp),%xmm3
        mulps nb231_dyH1(%rsp),%xmm4
        mulps nb231_dzH1(%rsp),%xmm5
        mulps nb231_dxH2(%rsp),%xmm6
        mulps nb231_dyH2(%rsp),%xmm7
        mulps nb231_dzH2(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb231_fixO(%rsp),%xmm0
    addps nb231_fiyO(%rsp),%xmm1
    addps nb231_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb231_fixH1(%rsp),%xmm3
    addps nb231_fiyH1(%rsp),%xmm4
    addps nb231_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb231_fixH2(%rsp),%xmm6
    addps nb231_fiyH2(%rsp),%xmm7
    addps nb231_fizH2(%rsp),%xmm8

    movaps %xmm0,nb231_fixO(%rsp)
    movaps %xmm1,nb231_fiyO(%rsp)
    movaps %xmm2,nb231_fizO(%rsp)
    movaps %xmm3,nb231_fixH1(%rsp)
    movaps %xmm4,nb231_fiyH1(%rsp)
    movaps %xmm5,nb231_fizH1(%rsp)
    movaps %xmm6,nb231_fixH2(%rsp)
    movaps %xmm7,nb231_fiyH2(%rsp)
    movaps %xmm8,nb231_fizH2(%rsp)

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
        subl $4,nb231_innerk(%rsp)
        jl    _nb_kernel231_x86_64_sse.nb231_odd_inner
        jmp   _nb_kernel231_x86_64_sse.nb231_unroll_loop
_nb_kernel231_x86_64_sse.nb231_odd_inner: 
        addl $4,nb231_innerk(%rsp)
        jnz   _nb_kernel231_x86_64_sse.nb231_odd_loop
        jmp   _nb_kernel231_x86_64_sse.nb231_updateouterdata
_nb_kernel231_x86_64_sse.nb231_odd_loop: 
        movq  nb231_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb231_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb231_iqO(%rsp),%xmm4
        movq nb231_charge(%rbp),%rsi
        movhps nb231_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb231_qqO(%rsp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb231_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb231_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb231_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb231_c6(%rsp)
        movaps %xmm7,nb231_c12(%rsp)

        movq nb231_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb231_ixO(%rsp),%xmm0
        movss nb231_iyO(%rsp),%xmm1
        movss nb231_izO(%rsp),%xmm2

        movlps nb231_ixH1(%rsp),%xmm6
        movlps nb231_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb231_iyH1(%rsp),%xmm6
        movlps nb231_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb231_izH1(%rsp),%xmm6
        movlps nb231_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb231_dxO(%rsp)
        movaps %xmm4,nb231_dyO(%rsp)
        movaps %xmm5,nb231_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        movaps %xmm0,nb231_rsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb231_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb231_half(%rsp),%xmm0
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
        mulps  nb231_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb231_VFtab(%rbp),%rsi
        movd %mm6,%r8d

        ## dispersion 
        movlps (%rsi,%r8,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 8(%rsi,%r8,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## dispersion table ready, in xmm4-xmm7         

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb231_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm7       ## fijD 
        mulss  %xmm4,%xmm5       ## Vvdw6 
        mulss  nb231_tsc(%rsp),%xmm7
        ## put scalar force on stack Update Vvdwtot directly 
        addss  nb231_Vvdwtot(%rsp),%xmm5
        movss %xmm7,nb231_fstmp(%rsp)
        movss %xmm5,nb231_Vvdwtot(%rsp)

        ## repulsion 
        movlps 16(%rsi,%r8,4),%xmm5
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## constant 10001000
        shufps $221,%xmm7,%xmm5 ## constant 11011101

        movlps 24(%rsi,%r8,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## constant 10001000
        shufps $221,%xmm3,%xmm7 ## constant 11011101
        ## table ready, in xmm4-xmm7    
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        mulss  nb231_two(%rsp),%xmm7    ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movss nb231_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        mulss  nb231_tsc(%rsp),%xmm7
        addss  nb231_fstmp(%rsp),%xmm7
        movss %xmm7,nb231_fstmp(%rsp)

        addss  nb231_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb231_Vvdwtot(%rsp)

        movaps nb231_rsqO(%rsp),%xmm3
        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps %xmm0,%xmm4
        mulps nb231_krf(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        mulps  nb231_two(%rsp),%xmm3
        subps  nb231_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb231_qqO(%rsp),%xmm0    ## xmm0=vcoul 
        mulps  nb231_qqO(%rsp),%xmm1    ## xmm1=coul part of fs 

        mulps  %xmm4,%xmm1
        subss  nb231_fstmp(%rsp),%xmm1
        mulps  %xmm1,%xmm4

        addps  nb231_vctot(%rsp),%xmm0
        movaps %xmm0,nb231_vctot(%rsp)

        movaps nb231_dxO(%rsp),%xmm0
        movaps nb231_dyO(%rsp),%xmm1
        movaps nb231_dzO(%rsp),%xmm2


        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb231_fixO(%rsp),%xmm3
        movss  nb231_fiyO(%rsp),%xmm4
        movss  nb231_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb231_fixO(%rsp)
        movss  %xmm4,nb231_fiyO(%rsp)
        movss  %xmm5,nb231_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## constant 11100110     ;# shift right 
        shufps $230,%xmm4,%xmm4 ## constant 11100110
        shufps $230,%xmm5,%xmm5 ## constant 11100110
        addss  nb231_fixH1(%rsp),%xmm3
        addss  nb231_fiyH1(%rsp),%xmm4
        addss  nb231_fizH1(%rsp),%xmm5
        movss  %xmm3,nb231_fixH1(%rsp)
        movss  %xmm4,nb231_fiyH1(%rsp)
        movss  %xmm5,nb231_fizH1(%rsp)          ## updated the H1 force 

        movq nb231_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## constant 11100111     ;# shift right 
        shufps $231,%xmm4,%xmm4 ## constant 11100111
        shufps $231,%xmm5,%xmm5 ## constant 11100111
        addss  nb231_fixH2(%rsp),%xmm3
        addss  nb231_fiyH2(%rsp),%xmm4
        addss  nb231_fizH2(%rsp),%xmm5
        movss  %xmm3,nb231_fixH2(%rsp)
        movss  %xmm4,nb231_fiyH2(%rsp)
        movss  %xmm5,nb231_fizH2(%rsp)          ## updated the H2 force 

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

        decl nb231_innerk(%rsp)
        jz    _nb_kernel231_x86_64_sse.nb231_updateouterdata
        jmp   _nb_kernel231_x86_64_sse.nb231_odd_loop
_nb_kernel231_x86_64_sse.nb231_updateouterdata: 
        movl  nb231_ii3(%rsp),%ecx
        movq  nb231_faction(%rbp),%rdi
        movq  nb231_fshift(%rbp),%rsi
        movl  nb231_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb231_fixO(%rsp),%xmm0
        movaps nb231_fiyO(%rsp),%xmm1
        movaps nb231_fizO(%rsp),%xmm2

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
        movaps nb231_fixH1(%rsp),%xmm0
        movaps nb231_fiyH1(%rsp),%xmm1
        movaps nb231_fizH1(%rsp),%xmm2

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
        movaps nb231_fixH2(%rsp),%xmm0
        movaps nb231_fiyH2(%rsp),%xmm1
        movaps nb231_fizH2(%rsp),%xmm2

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
        movl nb231_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb231_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb231_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb231_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb231_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb231_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb231_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231_x86_64_sse.nb231_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231_n(%rsp)
        jmp _nb_kernel231_x86_64_sse.nb231_outer
_nb_kernel231_x86_64_sse.nb231_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231_x86_64_sse.nb231_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231_x86_64_sse.nb231_threadloop
_nb_kernel231_x86_64_sse.nb231_end: 
        movl nb231_nouter(%rsp),%eax
        movl nb231_ninner(%rsp),%ebx
        movq nb231_outeriter(%rbp),%rcx
        movq nb231_inneriter(%rbp),%rdx
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




.globl nb_kernel231nf_x86_64_sse
.globl _nb_kernel231nf_x86_64_sse
nb_kernel231nf_x86_64_sse:      
_nb_kernel231nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb231nf_fshift, 16
.set nb231nf_gid, 24
.set nb231nf_pos, 32
.set nb231nf_faction, 40
.set nb231nf_charge, 48
.set nb231nf_p_facel, 56
.set nb231nf_argkrf, 64
.set nb231nf_argcrf, 72
.set nb231nf_Vc, 80
.set nb231nf_type, 88
.set nb231nf_p_ntype, 96
.set nb231nf_vdwparam, 104
.set nb231nf_Vvdw, 112
.set nb231nf_p_tabscale, 120
.set nb231nf_VFtab, 128
.set nb231nf_invsqrta, 136
.set nb231nf_dvda, 144
.set nb231nf_p_gbtabscale, 152
.set nb231nf_GBtab, 160
.set nb231nf_p_nthreads, 168
.set nb231nf_count, 176
.set nb231nf_mtx, 184
.set nb231nf_outeriter, 192
.set nb231nf_inneriter, 200
.set nb231nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb231nf_ixO, 0
.set nb231nf_iyO, 16
.set nb231nf_izO, 32
.set nb231nf_ixH1, 48
.set nb231nf_iyH1, 64
.set nb231nf_izH1, 80
.set nb231nf_ixH2, 96
.set nb231nf_iyH2, 112
.set nb231nf_izH2, 128
.set nb231nf_iqO, 144
.set nb231nf_iqH, 160
.set nb231nf_qqO, 176
.set nb231nf_qqH, 192
.set nb231nf_c6, 208
.set nb231nf_c12, 224
.set nb231nf_vctot, 240
.set nb231nf_Vvdwtot, 256
.set nb231nf_half, 272
.set nb231nf_three, 288
.set nb231nf_krf, 304
.set nb231nf_crf, 320
.set nb231nf_rinvO, 336
.set nb231nf_rinvH1, 352
.set nb231nf_rinvH2, 368
.set nb231nf_krsqO, 384
.set nb231nf_krsqH1, 400
.set nb231nf_krsqH2, 416
.set nb231nf_tsc, 432
.set nb231nf_facel, 448
.set nb231nf_iinr, 456
.set nb231nf_jindex, 464
.set nb231nf_jjnr, 472
.set nb231nf_shift, 480
.set nb231nf_shiftvec, 488
.set nb231nf_innerjjnr, 496
.set nb231nf_nri, 504
.set nb231nf_ntia, 508
.set nb231nf_is3, 512
.set nb231nf_ii3, 516
.set nb231nf_innerk, 520
.set nb231nf_n, 524
.set nb231nf_nn1, 528
.set nb231nf_nouter, 532
.set nb231nf_ninner, 536

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $552,%rsp       # # local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb231nf_nouter(%rsp)
        movl %eax,nb231nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb231nf_nri(%rsp)
        movq %rsi,nb231nf_iinr(%rsp)
        movq %rdx,nb231nf_jindex(%rsp)
        movq %rcx,nb231nf_jjnr(%rsp)
        movq %r8,nb231nf_shift(%rsp)
        movq %r9,nb231nf_shiftvec(%rsp)
        movq nb231nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb231nf_facel(%rsp)

        movq nb231nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb231nf_tsc(%rsp)

        movq nb231nf_argkrf(%rbp),%rsi
        movq nb231nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb231nf_krf(%rsp)
        movaps %xmm2,nb231nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb231nf_half(%rsp)
        movss nb231nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb231nf_half(%rsp)
        movaps %xmm3,nb231nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb231nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb231nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb231nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb231nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb231nf_iqO(%rsp)
        movaps %xmm4,nb231nf_iqH(%rsp)

        movq  nb231nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb231nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb231nf_ntia(%rsp)

_nb_kernel231nf_x86_64_sse.nb231nf_threadloop: 
        movq  nb231nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel231nf_x86_64_sse.nb231nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel231nf_x86_64_sse.nb231nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb231nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb231nf_n(%rsp)
        movl %ebx,nb231nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel231nf_x86_64_sse.nb231nf_outerstart
        jmp _nb_kernel231nf_x86_64_sse.nb231nf_end

_nb_kernel231nf_x86_64_sse.nb231nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb231nf_nouter(%rsp),%ebx
        movl %ebx,nb231nf_nouter(%rsp)

_nb_kernel231nf_x86_64_sse.nb231nf_outer: 
        movq  nb231nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb231nf_is3(%rsp)            ## store is3 

        movq  nb231nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb231nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb231nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb231nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb231nf_ixO(%rsp)
        movaps %xmm4,nb231nf_iyO(%rsp)
        movaps %xmm5,nb231nf_izO(%rsp)

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
        movaps %xmm0,nb231nf_ixH1(%rsp)
        movaps %xmm1,nb231nf_iyH1(%rsp)
        movaps %xmm2,nb231nf_izH1(%rsp)
        movaps %xmm3,nb231nf_ixH2(%rsp)
        movaps %xmm4,nb231nf_iyH2(%rsp)
        movaps %xmm5,nb231nf_izH2(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb231nf_vctot(%rsp)
        movaps %xmm4,nb231nf_Vvdwtot(%rsp)

        movq  nb231nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb231nf_pos(%rbp),%rsi
        movq  nb231nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb231nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb231nf_ninner(%rsp),%ecx
        movl  %ecx,nb231nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb231nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel231nf_x86_64_sse.nb231nf_unroll_loop
        jmp   _nb_kernel231nf_x86_64_sse.nb231nf_odd_inner
_nb_kernel231nf_x86_64_sse.nb231nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb231nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb231nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb231nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## constant 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb231nf_iqO(%rsp),%xmm3
        mulps  nb231nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb231nf_qqO(%rsp)
        movaps  %xmm4,nb231nf_qqH(%rsp)

        movq nb231nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb231nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb231nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb231nf_c6(%rsp)
        movaps %xmm6,nb231nf_c12(%rsp)

        movq nb231nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb231nf_ixO(%rsp),%xmm4
        movaps nb231nf_iyO(%rsp),%xmm5
        movaps nb231nf_izO(%rsp),%xmm6

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
        movaps nb231nf_ixH1(%rsp),%xmm4
        movaps nb231nf_iyH1(%rsp),%xmm5
        movaps nb231nf_izH1(%rsp),%xmm6

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
        movaps nb231nf_ixH2(%rsp),%xmm3
        movaps nb231nf_iyH2(%rsp),%xmm4
        movaps nb231nf_izH2(%rsp),%xmm5

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


        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb231nf_krf(%rsp),%xmm0
        mulps  nb231nf_krf(%rsp),%xmm1
        mulps  nb231nf_krf(%rsp),%xmm2

        movaps %xmm0,nb231nf_krsqH2(%rsp)
        movaps %xmm1,nb231nf_krsqH1(%rsp)
        movaps %xmm2,nb231nf_krsqO(%rsp)

        ## start with rsqO (still in xmm7) - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%rsp),%xmm4
        movaps  %xmm4,nb231nf_rinvO(%rsp)       ## rinvO in xmm0, rsqO in xmm7

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%rsp),%xmm4
        movaps  %xmm4,nb231nf_rinvH1(%rsp)      ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb231nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## constant 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb231nf_half(%rsp),%xmm4
        movaps  %xmm4,nb231nf_rinvH2(%rsp)      ## rinvH2 in xmm5 


        ## do O table interactions - rsqO in xmm7.
        mulps nb231nf_rinvO(%rsp),%xmm7
        mulps nb231nf_tsc(%rsp),%xmm7   ## rtab

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

        movq nb231nf_VFtab(%rbp),%rsi
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

        movaps nb231nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5       ## Vvdw6 

        addps  nb231nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb231nf_Vvdwtot(%rsp)

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

        movaps nb231nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb231nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb231nf_Vvdwtot(%rsp)

        movaps nb231nf_rinvO(%rsp),%xmm2
        movaps nb231nf_krsqO(%rsp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%rsp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqO(%rsp),%xmm2

        addps  nb231nf_vctot(%rsp),%xmm2
        movaps %xmm2,nb231nf_vctot(%rsp)

        ## H1 interactions 
        movaps nb231nf_rinvH1(%rsp),%xmm2
        movaps nb231nf_krsqH1(%rsp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%rsp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqH(%rsp),%xmm2

        addps  nb231nf_vctot(%rsp),%xmm2
        movaps %xmm2,nb231nf_vctot(%rsp)

        ## H2 interactions 
        movaps nb231nf_rinvH2(%rsp),%xmm2
        movaps nb231nf_krsqH2(%rsp),%xmm1
        addps  %xmm1,%xmm2 ## rinv+krsq
        subps  nb231nf_crf(%rsp),%xmm2   ## rinv+krsq-crf
        mulps  nb231nf_qqH(%rsp),%xmm2
        addps  nb231nf_vctot(%rsp),%xmm2
        movaps %xmm2,nb231nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb231nf_innerk(%rsp)
        jl    _nb_kernel231nf_x86_64_sse.nb231nf_odd_inner
        jmp   _nb_kernel231nf_x86_64_sse.nb231nf_unroll_loop
_nb_kernel231nf_x86_64_sse.nb231nf_odd_inner: 
        addl $4,nb231nf_innerk(%rsp)
        jnz   _nb_kernel231nf_x86_64_sse.nb231nf_odd_loop
        jmp   _nb_kernel231nf_x86_64_sse.nb231nf_updateouterdata
_nb_kernel231nf_x86_64_sse.nb231nf_odd_loop: 
        movq  nb231nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb231nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb231nf_iqO(%rsp),%xmm4
        movq nb231nf_charge(%rbp),%rsi
        movhps nb231nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb231nf_qqO(%rsp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb231nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb231nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb231nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## constant 11111100
        shufps $253,%xmm7,%xmm7 ## constant 11111101
        movaps %xmm6,nb231nf_c6(%rsp)
        movaps %xmm7,nb231nf_c12(%rsp)

        movq nb231nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb231nf_ixO(%rsp),%xmm3
        movss nb231nf_iyO(%rsp),%xmm4
        movss nb231nf_izO(%rsp),%xmm5

        movlps nb231nf_ixH1(%rsp),%xmm6
        movlps nb231nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb231nf_iyH1(%rsp),%xmm6
        movlps nb231nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb231nf_izH1(%rsp),%xmm6
        movlps nb231nf_izH2(%rsp),%xmm7
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

        movaps %xmm4,%xmm0
        mulps nb231nf_krf(%rsp),%xmm0
        movaps %xmm0,nb231nf_krsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb231nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb231nf_half(%rsp),%xmm0
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
        mulps  nb231nf_tsc(%rsp),%xmm4   ## rtab

        cvttps2pi %xmm4,%mm6
        cvtpi2ps %mm6,%xmm6
        subss  %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $3,%mm6

        movq nb231nf_VFtab(%rbp),%rsi
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

        movss nb231nf_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm5       ## Vvdw6 

        addss  nb231nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb231nf_Vvdwtot(%rsp)

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

        movss nb231nf_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  nb231nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb231nf_Vvdwtot(%rsp)

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps %xmm0,%xmm4
        movaps nb231nf_krsqO(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 

        subps  nb231nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 

        mulps  nb231nf_qqO(%rsp),%xmm0          ## xmm0=vcoul 

        addps  nb231nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb231nf_vctot(%rsp)

        decl nb231nf_innerk(%rsp)
        jz    _nb_kernel231nf_x86_64_sse.nb231nf_updateouterdata
        jmp   _nb_kernel231nf_x86_64_sse.nb231nf_odd_loop
_nb_kernel231nf_x86_64_sse.nb231nf_updateouterdata: 
        ## get n from stack
        movl nb231nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb231nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb231nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb231nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb231nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb231nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb231nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel231nf_x86_64_sse.nb231nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb231nf_n(%rsp)
        jmp _nb_kernel231nf_x86_64_sse.nb231nf_outer
_nb_kernel231nf_x86_64_sse.nb231nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb231nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel231nf_x86_64_sse.nb231nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel231nf_x86_64_sse.nb231nf_threadloop
_nb_kernel231nf_x86_64_sse.nb231nf_end: 
        movl nb231nf_nouter(%rsp),%eax
        movl nb231nf_ninner(%rsp),%ebx
        movq nb231nf_outeriter(%rbp),%rcx
        movq nb231nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $552,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



