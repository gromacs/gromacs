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







.globl nb_kernel311_x86_64_sse
.globl _nb_kernel311_x86_64_sse
nb_kernel311_x86_64_sse:        
_nb_kernel311_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb311_fshift, 16
.set nb311_gid, 24
.set nb311_pos, 32
.set nb311_faction, 40
.set nb311_charge, 48
.set nb311_p_facel, 56
.set nb311_argkrf, 64
.set nb311_argcrf, 72
.set nb311_Vc, 80
.set nb311_type, 88
.set nb311_p_ntype, 96
.set nb311_vdwparam, 104
.set nb311_Vvdw, 112
.set nb311_p_tabscale, 120
.set nb311_VFtab, 128
.set nb311_invsqrta, 136
.set nb311_dvda, 144
.set nb311_p_gbtabscale, 152
.set nb311_GBtab, 160
.set nb311_p_nthreads, 168
.set nb311_count, 176
.set nb311_mtx, 184
.set nb311_outeriter, 192
.set nb311_inneriter, 200
.set nb311_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb311_ixO, 0
.set nb311_iyO, 16
.set nb311_izO, 32
.set nb311_ixH1, 48
.set nb311_iyH1, 64
.set nb311_izH1, 80
.set nb311_ixH2, 96
.set nb311_iyH2, 112
.set nb311_izH2, 128
.set nb311_iqO, 144
.set nb311_iqH, 160
.set nb311_dxO, 176
.set nb311_dyO, 192
.set nb311_dzO, 208
.set nb311_dxH1, 224
.set nb311_dyH1, 240
.set nb311_dzH1, 256
.set nb311_dxH2, 272
.set nb311_dyH2, 288
.set nb311_dzH2, 304
.set nb311_qqO, 320
.set nb311_qqH, 336
.set nb311_rinvO, 352
.set nb311_rinvH1, 368
.set nb311_rinvH2, 384
.set nb311_rO, 400
.set nb311_rH1, 416
.set nb311_rH2, 432
.set nb311_tsc, 448
.set nb311_fstmp, 464
.set nb311_c6, 480
.set nb311_c12, 496
.set nb311_six, 512
.set nb311_twelve, 528
.set nb311_vctot, 544
.set nb311_Vvdwtot, 560
.set nb311_fixO, 576
.set nb311_fiyO, 592
.set nb311_fizO, 608
.set nb311_fixH1, 624
.set nb311_fiyH1, 640
.set nb311_fizH1, 656
.set nb311_fixH2, 672
.set nb311_fiyH2, 688
.set nb311_fizH2, 704
.set nb311_epsO, 720
.set nb311_epsH1, 736
.set nb311_epsH2, 752
.set nb311_half, 768
.set nb311_three, 784
.set nb311_is3, 800
.set nb311_ii3, 804
.set nb311_nri, 808
.set nb311_iinr, 816
.set nb311_jindex, 824
.set nb311_jjnr, 832
.set nb311_shift, 840
.set nb311_shiftvec, 848
.set nb311_facel, 856
.set nb311_innerjjnr, 864
.set nb311_ntia, 872
.set nb311_innerk, 876
.set nb311_n, 880
.set nb311_nn1, 884
.set nb311_nouter, 888
.set nb311_ninner, 892
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $904,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb311_nouter(%rsp)
        movl %eax,nb311_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb311_nri(%rsp)
        movq %rsi,nb311_iinr(%rsp)
        movq %rdx,nb311_jindex(%rsp)
        movq %rcx,nb311_jjnr(%rsp)
        movq %r8,nb311_shift(%rsp)
        movq %r9,nb311_shiftvec(%rsp)
        movq nb311_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb311_facel(%rsp)

        movq nb311_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb311_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb311_half(%rsp)
        movss nb311_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb311_half(%rsp)
        movaps %xmm3,nb311_three(%rsp)
        movaps %xmm4,nb311_six(%rsp)
        movaps %xmm5,nb311_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb311_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb311_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb311_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb311_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb311_iqO(%rsp)
        movaps %xmm4,nb311_iqH(%rsp)

        movq  nb311_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb311_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311_ntia(%rsp)

_nb_kernel311_x86_64_sse.nb311_threadloop: 
        movq  nb311_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel311_x86_64_sse.nb311_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311_x86_64_sse.nb311_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311_n(%rsp)
        movl %ebx,nb311_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311_x86_64_sse.nb311_outerstart
        jmp _nb_kernel311_x86_64_sse.nb311_end

_nb_kernel311_x86_64_sse.nb311_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311_nouter(%rsp),%ebx
        movl %ebx,nb311_nouter(%rsp)

_nb_kernel311_x86_64_sse.nb311_outer: 
        movq  nb311_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb311_is3(%rsp)      ## store is3 

        movq  nb311_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb311_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb311_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb311_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb311_ixO(%rsp)
        movaps %xmm4,nb311_iyO(%rsp)
        movaps %xmm5,nb311_izO(%rsp)

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
        movaps %xmm0,nb311_ixH1(%rsp)
        movaps %xmm1,nb311_iyH1(%rsp)
        movaps %xmm2,nb311_izH1(%rsp)
        movaps %xmm3,nb311_ixH2(%rsp)
        movaps %xmm4,nb311_iyH2(%rsp)
        movaps %xmm5,nb311_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb311_vctot(%rsp)
        movaps %xmm4,nb311_Vvdwtot(%rsp)
        movaps %xmm4,nb311_fixO(%rsp)
        movaps %xmm4,nb311_fiyO(%rsp)
        movaps %xmm4,nb311_fizO(%rsp)
        movaps %xmm4,nb311_fixH1(%rsp)
        movaps %xmm4,nb311_fiyH1(%rsp)
        movaps %xmm4,nb311_fizH1(%rsp)
        movaps %xmm4,nb311_fixH2(%rsp)
        movaps %xmm4,nb311_fiyH2(%rsp)
        movaps %xmm4,nb311_fizH2(%rsp)

        movq  nb311_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb311_pos(%rbp),%rsi
        movq  nb311_faction(%rbp),%rdi
        movq  nb311_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb311_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb311_ninner(%rsp),%ecx
        movl  %ecx,nb311_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb311_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel311_x86_64_sse.nb311_unroll_loop
        jmp   _nb_kernel311_x86_64_sse.nb311_odd_inner
_nb_kernel311_x86_64_sse.nb311_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb311_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb311_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb311_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb311_iqO(%rsp),%xmm3
        mulps  nb311_iqH(%rsp),%xmm4

        movaps  %xmm3,nb311_qqO(%rsp)
        movaps  %xmm4,nb311_qqH(%rsp)

        movq nb311_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb311_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb311_ntia(%rsp),%edi
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

        movaps %xmm4,nb311_c6(%rsp)
        movaps %xmm6,nb311_c12(%rsp)

        movq nb311_pos(%rbp),%rsi        ## base of pos[] 

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

    movd %eax,%mm0
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    subps nb311_ixO(%rsp),%xmm0
    subps nb311_iyO(%rsp),%xmm1
    subps nb311_izO(%rsp),%xmm2
    subps nb311_ixH1(%rsp),%xmm3
    subps nb311_iyH1(%rsp),%xmm4
    subps nb311_izH1(%rsp),%xmm5
    subps nb311_ixH2(%rsp),%xmm6
    subps nb311_iyH2(%rsp),%xmm7
    subps nb311_izH2(%rsp),%xmm8

        movaps %xmm0,nb311_dxO(%rsp)
        movaps %xmm1,nb311_dyO(%rsp)
        movaps %xmm2,nb311_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb311_dxH1(%rsp)
        movaps %xmm4,nb311_dyH1(%rsp)
        movaps %xmm5,nb311_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb311_dxH2(%rsp)
        movaps %xmm7,nb311_dyH2(%rsp)
        movaps %xmm8,nb311_dzH2(%rsp)
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

        movaps  nb311_three(%rsp),%xmm9
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

        movaps  nb311_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvO
        mulps   %xmm4,%xmm10 ## rinvH1
    mulps   %xmm4,%xmm11 ## rinvH2

        movaps  %xmm9,nb311_rinvO(%rsp)
        movaps  %xmm10,nb311_rinvH1(%rsp)
        movaps  %xmm11,nb311_rinvH2(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb311_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttps2dq %xmm0,%xmm1
    cvttps2dq %xmm3,%xmm4
    cvttps2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2ps  %xmm1,%xmm2
    cvtdq2ps  %xmm4,%xmm5
    cvtdq2ps  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    movhlps %xmm1,%xmm13
    movhlps %xmm4,%xmm14
    movhlps %xmm7,%xmm15
    movd    %xmm1,%eax
    movd    %xmm4,%r8d
    movd    %xmm7,%r12d
    movd    %xmm13,%ecx
    movd    %xmm14,%r10d
    movd    %xmm15,%r14d
    pshufd $1,%xmm1,%xmm1
    pshufd $1,%xmm4,%xmm4
    pshufd $1,%xmm7,%xmm7
    pshufd $1,%xmm13,%xmm13
    pshufd $1,%xmm14,%xmm14
    pshufd $1,%xmm15,%xmm15
    movd    %xmm1,%ebx
    movd    %xmm4,%r9d
    movd    %xmm7,%r13d
    movd    %xmm13,%edx
    movd    %xmm14,%r11d
    movd    %xmm15,%r15d

    movq nb311_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb311_epsO(%rsp)
    movaps    %xmm3,nb311_epsH1(%rsp)
    movaps    %xmm6,nb311_epsH2(%rsp)

    ## Load LOTS of table data
        movlps (%rsi,%rax,4),%xmm1
        movlps (%rsi,%r8,4),%xmm5
        movlps (%rsi,%r12,4),%xmm9

        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%r10,4),%xmm7
        movlps (%rsi,%r14,4),%xmm11

        movhps (%rsi,%rbx,4),%xmm1
        movhps (%rsi,%r9,4),%xmm5
        movhps (%rsi,%r13,4),%xmm9

        movhps (%rsi,%rdx,4),%xmm3
        movhps (%rsi,%r11,4),%xmm7
        movhps (%rsi,%r15,4),%xmm11

    movaps %xmm1,%xmm0
    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm3,%xmm0 ## 10001000
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm3,%xmm1 ## 11011101
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm3
        movlps 8(%rsi,%r8,4),%xmm7
        movlps 8(%rsi,%r12,4),%xmm11

        movlps 8(%rsi,%rcx,4),%xmm12
        movlps 8(%rsi,%r10,4),%xmm13
        movlps 8(%rsi,%r14,4),%xmm14

        movhps 8(%rsi,%rbx,4),%xmm3
        movhps 8(%rsi,%r9,4),%xmm7
        movhps 8(%rsi,%r13,4),%xmm11

        movhps 8(%rsi,%rdx,4),%xmm12
        movhps 8(%rsi,%r11,4),%xmm13
        movhps 8(%rsi,%r15,4),%xmm14

    movaps %xmm3,%xmm2
    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm12,%xmm2 ## 10001000
        shufps $136,%xmm13,%xmm6 ## 10001000
        shufps $136,%xmm14,%xmm10 ## 10001000
        shufps $221,%xmm12,%xmm3 ## 11011101
        shufps $221,%xmm13,%xmm7 ## 11011101
        shufps $221,%xmm14,%xmm11 ## 11011101
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movaps nb311_epsO(%rsp),%xmm12
    movaps nb311_epsH1(%rsp),%xmm13
    movaps nb311_epsH2(%rsp),%xmm14

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm13,%xmm7
    mulps  %xmm14,%xmm11
    mulps  %xmm12,%xmm2  ## Geps
    mulps  %xmm13,%xmm6
    mulps  %xmm14,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm13,%xmm7
    mulps  %xmm14,%xmm11

    addps  %xmm2,%xmm1  ## F+Geps
    addps  %xmm6,%xmm5
    addps  %xmm10,%xmm9
    addps  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm5
    addps  %xmm11,%xmm9
    addps  %xmm3,%xmm3   ## 2*Heps2
    addps  %xmm7,%xmm7
    addps  %xmm11,%xmm11
    addps  %xmm2,%xmm3   ## 2*Heps2+Geps
    addps  %xmm6,%xmm7
    addps  %xmm10,%xmm11
    addps  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm5,%xmm7
    addps  %xmm9,%xmm11
    mulps  %xmm12,%xmm1  ## eps*Fp
    mulps  %xmm13,%xmm5
    mulps  %xmm14,%xmm9
    movaps nb311_qqO(%rsp),%xmm12
    movaps nb311_qqH(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    ## calculate LJ
    movaps nb311_rinvO(%rsp),%xmm12
    mulps  %xmm12,%xmm12 ## rinvsq
    movaps %xmm12,%xmm13 ## rinvsq
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm13,%xmm12 ## rinv6
    movaps %xmm12,%xmm13 ## rinv6
    mulps  %xmm12,%xmm12 ## rinv12
        mulps  nb311_c6(%rsp),%xmm13
        mulps  nb311_c12(%rsp),%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14

        addps  nb311_Vvdwtot(%rsp),%xmm14
        mulps  nb311_six(%rsp),%xmm13
        mulps  nb311_twelve(%rsp),%xmm12
        movaps %xmm14,nb311_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal    
    mulps  nb311_rinvO(%rsp),%xmm12
    movaps %xmm12,nb311_fstmp(%rsp)


    ## accumulate vctot
    addps  nb311_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb311_vctot(%rsp)

    movaps nb311_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

    movq nb311_faction(%rbp),%rdi
        ## move j forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm11 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm12 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm11 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm12 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm13    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm14    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm2     ## jzb
    movss  8(%rdi,%rdx,4),%xmm5     ## jzd
    movlhps %xmm2,%xmm13 ## jza  -  jzb  -
    movlhps %xmm5,%xmm14 ## jzc  -  jzd -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jza jzb jzc jzd

    ## xmm11: jxa jya jxb jyb 
    ## xmm12: jxc jyc jxd jyd
    ## xmm13: jza jzb jzc jzd

    movaps nb311_fstmp(%rsp),%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    mulps  nb311_rinvO(%rsp),%xmm0
    mulps  nb311_rinvH1(%rsp),%xmm4
    mulps  nb311_rinvH2(%rsp),%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb311_dxO(%rsp),%xmm0
        mulps nb311_dyO(%rsp),%xmm1
        mulps nb311_dzO(%rsp),%xmm2
        mulps nb311_dxH1(%rsp),%xmm3
        mulps nb311_dyH1(%rsp),%xmm4
        mulps nb311_dzH1(%rsp),%xmm5
        mulps nb311_dxH2(%rsp),%xmm6
        mulps nb311_dyH2(%rsp),%xmm7
        mulps nb311_dzH2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb311_fixO(%rsp),%xmm0
    addps nb311_fiyO(%rsp),%xmm1
    addps nb311_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb311_fixH1(%rsp),%xmm3
    addps nb311_fiyH1(%rsp),%xmm4
    addps nb311_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb311_fixH2(%rsp),%xmm6
    addps nb311_fiyH2(%rsp),%xmm7
    addps nb311_fizH2(%rsp),%xmm8

    movaps %xmm0,nb311_fixO(%rsp)
    movaps %xmm1,nb311_fiyO(%rsp)
    movaps %xmm2,nb311_fizO(%rsp)
    movaps %xmm3,nb311_fixH1(%rsp)
    movaps %xmm4,nb311_fiyH1(%rsp)
    movaps %xmm5,nb311_fizH1(%rsp)
    movaps %xmm6,nb311_fixH2(%rsp)
    movaps %xmm7,nb311_fiyH2(%rsp)
    movaps %xmm8,nb311_fizH2(%rsp)

    ## xmm14 = fjx
    ## xmm15 = fjy
    ## xmm13 = fjz
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fjzc fjzd

    movlps %xmm11,(%rdi,%rax,4)
    movhps %xmm11,(%rdi,%rbx,4)
    movlps %xmm12,(%rdi,%rcx,4)
    movhps %xmm12,(%rdi,%rdx,4)
    movss  %xmm13,8(%rdi,%rax,4)
    movss  %xmm14,8(%rdi,%rcx,4)
    shufps $1,%xmm13,%xmm13
    shufps $1,%xmm14,%xmm14
    movss  %xmm13,8(%rdi,%rbx,4)
    movss  %xmm14,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb311_innerk(%rsp)
        jl    _nb_kernel311_x86_64_sse.nb311_odd_inner
        jmp   _nb_kernel311_x86_64_sse.nb311_unroll_loop
_nb_kernel311_x86_64_sse.nb311_odd_inner: 
        addl $4,nb311_innerk(%rsp)
        jnz   _nb_kernel311_x86_64_sse.nb311_odd_loop
        jmp   _nb_kernel311_x86_64_sse.nb311_updateouterdata
_nb_kernel311_x86_64_sse.nb311_odd_loop: 
        movq  nb311_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb311_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb311_iqO(%rsp),%xmm4
        movq nb311_charge(%rbp),%rsi
        movhps nb311_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb311_qqO(%rsp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb311_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb311_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb311_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb311_c6(%rsp)
        movaps %xmm7,nb311_c12(%rsp)

        movq nb311_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb311_ixO(%rsp),%xmm0
        movss nb311_iyO(%rsp),%xmm1
        movss nb311_izO(%rsp),%xmm2

        movlps nb311_ixH1(%rsp),%xmm6
        movlps nb311_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb311_iyH1(%rsp),%xmm6
        movlps nb311_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb311_izH1(%rsp),%xmm6
        movlps nb311_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb311_dxO(%rsp)
        movaps %xmm4,nb311_dyO(%rsp)
        movaps %xmm5,nb311_dzO(%rsp)

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
        movaps nb311_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb311_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000
        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb311_rinvO(%rsp)

        mulps nb311_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7   ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm7
    movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %eax,%mm0
    movd %ecx,%mm1
    movd %edx,%mm2

    movq nb311_VFtab(%rbp),%rsi
    movd %mm6,%eax
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx

    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## coulomb table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    addps  %xmm7,%xmm7      ## two*Heps2 
    movaps nb311_qqO(%rsp),%xmm0
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm0 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm0 fijC 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb311_vctot(%rsp),%xmm5
    movaps %xmm5,nb311_vctot(%rsp)

        ## do nontable L-J 
        movaps nb311_rinvO(%rsp),%xmm2
        mulps  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb311_c6(%rsp),%xmm1
        mulps  nb311_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subps  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulps  nb311_six(%rsp),%xmm1
        mulps  nb311_twelve(%rsp),%xmm4
        subps  %xmm1,%xmm4
        addps  nb311_Vvdwtot(%rsp),%xmm3
        mulps  nb311_rinvO(%rsp),%xmm4
        mulps  nb311_tsc(%rsp),%xmm0
        subps  %xmm0,%xmm4
        movaps %xmm3,nb311_Vvdwtot(%rsp)
        mulps  nb311_rinvO(%rsp),%xmm4

    movd %mm0,%eax
    movd %mm1,%ecx
    movd %mm2,%edx

        movaps nb311_dxO(%rsp),%xmm0
        movaps nb311_dyO(%rsp),%xmm1
        movaps nb311_dzO(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force) 
        movss  nb311_fixO(%rsp),%xmm3
        movss  nb311_fiyO(%rsp),%xmm4
        movss  nb311_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb311_fixO(%rsp)
        movss  %xmm4,nb311_fiyO(%rsp)
        movss  %xmm5,nb311_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb311_fixH1(%rsp),%xmm3
        addss  nb311_fiyH1(%rsp),%xmm4
        addss  nb311_fizH1(%rsp),%xmm5
        movss  %xmm3,nb311_fixH1(%rsp)
        movss  %xmm4,nb311_fiyH1(%rsp)
        movss  %xmm5,nb311_fizH1(%rsp)          ## updated the H1 force 

        movq nb311_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb311_fixH2(%rsp),%xmm3
        addss  nb311_fiyH2(%rsp),%xmm4
        addss  nb311_fizH2(%rsp),%xmm5
        movss  %xmm3,nb311_fixH2(%rsp)
        movss  %xmm4,nb311_fiyH2(%rsp)
        movss  %xmm5,nb311_fizH2(%rsp)          ## updated the H2 force 

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

        decl nb311_innerk(%rsp)
        jz    _nb_kernel311_x86_64_sse.nb311_updateouterdata
        jmp   _nb_kernel311_x86_64_sse.nb311_odd_loop
_nb_kernel311_x86_64_sse.nb311_updateouterdata: 
        movl  nb311_ii3(%rsp),%ecx
        movq  nb311_faction(%rbp),%rdi
        movq  nb311_fshift(%rbp),%rsi
        movl  nb311_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb311_fixO(%rsp),%xmm0
        movaps nb311_fiyO(%rsp),%xmm1
        movaps nb311_fizO(%rsp),%xmm2

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
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb311_fixH1(%rsp),%xmm0
        movaps nb311_fiyH1(%rsp),%xmm1
        movaps nb311_fizH1(%rsp),%xmm2

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
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb311_fixH2(%rsp),%xmm0
        movaps nb311_fiyH2(%rsp),%xmm1
        movaps nb311_fizH2(%rsp),%xmm2

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
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb311_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb311_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb311_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb311_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb311_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb311_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb311_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311_x86_64_sse.nb311_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311_n(%rsp)
        jmp _nb_kernel311_x86_64_sse.nb311_outer
_nb_kernel311_x86_64_sse.nb311_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311_x86_64_sse.nb311_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311_x86_64_sse.nb311_threadloop
_nb_kernel311_x86_64_sse.nb311_end: 
        movl nb311_nouter(%rsp),%eax
        movl nb311_ninner(%rsp),%ebx
        movq nb311_outeriter(%rbp),%rcx
        movq nb311_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $904,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel311nf_x86_64_sse
.globl _nb_kernel311nf_x86_64_sse
nb_kernel311nf_x86_64_sse:      
_nb_kernel311nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb311nf_fshift, 16
.set nb311nf_gid, 24
.set nb311nf_pos, 32
.set nb311nf_faction, 40
.set nb311nf_charge, 48
.set nb311nf_p_facel, 56
.set nb311nf_argkrf, 64
.set nb311nf_argcrf, 72
.set nb311nf_Vc, 80
.set nb311nf_type, 88
.set nb311nf_p_ntype, 96
.set nb311nf_vdwparam, 104
.set nb311nf_Vvdw, 112
.set nb311nf_p_tabscale, 120
.set nb311nf_VFtab, 128
.set nb311nf_invsqrta, 136
.set nb311nf_dvda, 144
.set nb311nf_p_gbtabscale, 152
.set nb311nf_GBtab, 160
.set nb311nf_p_nthreads, 168
.set nb311nf_count, 176
.set nb311nf_mtx, 184
.set nb311nf_outeriter, 192
.set nb311nf_inneriter, 200
.set nb311nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb311nf_ixO, 0
.set nb311nf_iyO, 16
.set nb311nf_izO, 32
.set nb311nf_ixH1, 48
.set nb311nf_iyH1, 64
.set nb311nf_izH1, 80
.set nb311nf_ixH2, 96
.set nb311nf_iyH2, 112
.set nb311nf_izH2, 128
.set nb311nf_iqO, 144
.set nb311nf_iqH, 160
.set nb311nf_qqO, 176
.set nb311nf_qqH, 192
.set nb311nf_rinvO, 208
.set nb311nf_rinvH1, 224
.set nb311nf_rinvH2, 240
.set nb311nf_rO, 256
.set nb311nf_rH1, 272
.set nb311nf_rH2, 288
.set nb311nf_tsc, 304
.set nb311nf_c6, 320
.set nb311nf_c12, 336
.set nb311nf_vctot, 352
.set nb311nf_Vvdwtot, 368
.set nb311nf_half, 384
.set nb311nf_three, 400
.set nb311nf_is3, 416
.set nb311nf_ii3, 420
.set nb311nf_nri, 424
.set nb311nf_iinr, 432
.set nb311nf_jindex, 440
.set nb311nf_jjnr, 448
.set nb311nf_shift, 456
.set nb311nf_shiftvec, 464
.set nb311nf_facel, 472
.set nb311nf_innerjjnr, 480
.set nb311nf_ntia, 488
.set nb311nf_innerk, 492
.set nb311nf_n, 496
.set nb311nf_nn1, 500
.set nb311nf_nouter, 504
.set nb311nf_ninner, 508
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
        movl %eax,nb311nf_nouter(%rsp)
        movl %eax,nb311nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb311nf_nri(%rsp)
        movq %rsi,nb311nf_iinr(%rsp)
        movq %rdx,nb311nf_jindex(%rsp)
        movq %rcx,nb311nf_jjnr(%rsp)
        movq %r8,nb311nf_shift(%rsp)
        movq %r9,nb311nf_shiftvec(%rsp)
        movq nb311nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb311nf_facel(%rsp)

        movq nb311nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb311nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb311nf_half(%rsp)
        movss nb311nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb311nf_half(%rsp)
        movaps %xmm3,nb311nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb311nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb311nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb311nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb311nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb311nf_iqO(%rsp)
        movaps %xmm4,nb311nf_iqH(%rsp)

        movq  nb311nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb311nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb311nf_ntia(%rsp)

_nb_kernel311nf_x86_64_sse.nb311nf_threadloop: 
        movq  nb311nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel311nf_x86_64_sse.nb311nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel311nf_x86_64_sse.nb311nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb311nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb311nf_n(%rsp)
        movl %ebx,nb311nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel311nf_x86_64_sse.nb311nf_outerstart
        jmp _nb_kernel311nf_x86_64_sse.nb311nf_end

_nb_kernel311nf_x86_64_sse.nb311nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb311nf_nouter(%rsp),%ebx
        movl %ebx,nb311nf_nouter(%rsp)

_nb_kernel311nf_x86_64_sse.nb311nf_outer: 
        movq  nb311nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb311nf_is3(%rsp)            ## store is3 

        movq  nb311nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb311nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb311nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb311nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb311nf_ixO(%rsp)
        movaps %xmm4,nb311nf_iyO(%rsp)
        movaps %xmm5,nb311nf_izO(%rsp)

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
        movaps %xmm0,nb311nf_ixH1(%rsp)
        movaps %xmm1,nb311nf_iyH1(%rsp)
        movaps %xmm2,nb311nf_izH1(%rsp)
        movaps %xmm3,nb311nf_ixH2(%rsp)
        movaps %xmm4,nb311nf_iyH2(%rsp)
        movaps %xmm5,nb311nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb311nf_vctot(%rsp)
        movaps %xmm4,nb311nf_Vvdwtot(%rsp)

        movq  nb311nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb311nf_pos(%rbp),%rsi
        movq  nb311nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb311nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb311nf_ninner(%rsp),%ecx
        movl  %ecx,nb311nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb311nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel311nf_x86_64_sse.nb311nf_unroll_loop
        jmp   _nb_kernel311nf_x86_64_sse.nb311nf_odd_inner
_nb_kernel311nf_x86_64_sse.nb311nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb311nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb311nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb311nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb311nf_iqO(%rsp),%xmm3
        mulps  nb311nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb311nf_qqO(%rsp)
        movaps  %xmm4,nb311nf_qqH(%rsp)

        movq nb311nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb311nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb311nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb311nf_c6(%rsp)
        movaps %xmm6,nb311nf_c12(%rsp)

        movq nb311nf_pos(%rbp),%rsi        ## base of pos[] 

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

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixO-izO to xmm4-xmm6 
        movaps nb311nf_ixO(%rsp),%xmm4
        movaps nb311nf_iyO(%rsp),%xmm5
        movaps nb311nf_izO(%rsp),%xmm6

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
        movaps nb311nf_ixH1(%rsp),%xmm4
        movaps nb311nf_iyH1(%rsp),%xmm5
        movaps nb311nf_izH1(%rsp),%xmm6

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
        movaps nb311nf_ixH2(%rsp),%xmm3
        movaps nb311nf_iyH2(%rsp),%xmm4
        movaps nb311nf_izH2(%rsp),%xmm5

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

        ## start with rsqO - seed to xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb311nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb311nf_half(%rsp),%xmm4
        movaps  %xmm4,nb311nf_rinvO(%rsp)       ## rinvO in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb311nf_rO(%rsp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb311nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb311nf_half(%rsp),%xmm4
        movaps  %xmm4,nb311nf_rinvH1(%rsp)      ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb311nf_rH1(%rsp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb311nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb311nf_half(%rsp),%xmm4
        movaps  %xmm4,nb311nf_rinvH2(%rsp)      ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb311nf_rH2(%rsp)

        ## do O interactions 
        ## rO is still in xmm7 
        mulps   nb311nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm4
    movlhps %xmm4,%xmm3

    subps %xmm3,%xmm7

        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %eax,%mm0
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    movq nb311nf_VFtab(%rbp),%rsi
    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rbx,4),%xmm5
    movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rbx,4),%xmm7
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb311nf_qqO(%rsp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV 

        ## do nontable L-J 
        movaps nb311nf_rinvO(%rsp),%xmm2
        mulps  %xmm2,%xmm2

    ## at this point mm5 contains vcoul and 
    ## increment vcoul - then we can get rid of mm5 
    addps  nb311nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb311nf_vctot(%rsp)

        movaps %xmm2,%xmm1
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb311nf_c6(%rsp),%xmm1
        mulps  nb311nf_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subps  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addps  nb311nf_Vvdwtot(%rsp),%xmm3
        movaps %xmm3,nb311nf_Vvdwtot(%rsp)

        ## Done with O interactions - now H1! 
        movaps nb311nf_rH1(%rsp),%xmm7
        mulps   nb311nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm4
    movlhps %xmm4,%xmm3

    subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rbx,4),%xmm5
    movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rbx,4),%xmm7
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp        
    movaps nb311nf_qqH(%rsp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addps  nb311nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb311nf_vctot(%rsp)

        ## Done with H1, finally we do H2 interactions 
        movaps nb311nf_rH2(%rsp),%xmm7
        mulps   nb311nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7   ## mm6/mm7 contain lu indices 

    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm4
    movlhps %xmm4,%xmm3

    subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rbx,4),%xmm5
    movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rbx,4),%xmm7
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## coulomb table ready, in xmm4-xmm7      

    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    movaps nb311nf_qqH(%rsp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul 
    addps  nb311nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb311nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb311nf_innerk(%rsp)
        jl    _nb_kernel311nf_x86_64_sse.nb311nf_odd_inner
        jmp   _nb_kernel311nf_x86_64_sse.nb311nf_unroll_loop
_nb_kernel311nf_x86_64_sse.nb311nf_odd_inner: 
        addl $4,nb311nf_innerk(%rsp)
        jnz   _nb_kernel311nf_x86_64_sse.nb311nf_odd_loop
        jmp   _nb_kernel311nf_x86_64_sse.nb311nf_updateouterdata
_nb_kernel311nf_x86_64_sse.nb311nf_odd_loop: 
        movq  nb311nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb311nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb311nf_iqO(%rsp),%xmm4
        movq nb311nf_charge(%rbp),%rsi
        movhps nb311nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb311nf_qqO(%rsp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb311nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb311nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb311nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb311nf_c6(%rsp)
        movaps %xmm7,nb311nf_c12(%rsp)

        movq nb311nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movss nb311nf_ixO(%rsp),%xmm3
        movss nb311nf_iyO(%rsp),%xmm4
        movss nb311nf_izO(%rsp),%xmm5

        movlps nb311nf_ixH1(%rsp),%xmm6
        movlps nb311nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb311nf_iyH1(%rsp),%xmm6
        movlps nb311nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb311nf_izH1(%rsp),%xmm6
        movlps nb311nf_izH2(%rsp),%xmm7
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
        movaps nb311nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb311nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb311nf_rinvO(%rsp)

        mulps nb311nf_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7   ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm7
    movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %eax,%mm0
    movd %ecx,%mm1
    movd %edx,%mm2

    movq nb311nf_VFtab(%rbp),%rsi
    movd %mm6,%eax
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx

    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## coulomb table ready, in xmm4-xmm7      
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp      
    movaps nb311nf_qqO(%rsp),%xmm0
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm0,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul  
    ## increment vcoul - then we can get rid of mm5 
    addps  nb311nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb311nf_vctot(%rsp)

        ## do nontable L-J 
        movaps nb311nf_rinvO(%rsp),%xmm2
        mulps  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb311nf_c6(%rsp),%xmm1
        mulps  nb311nf_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subps  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addps  nb311nf_Vvdwtot(%rsp),%xmm3
        movaps %xmm3,nb311nf_Vvdwtot(%rsp)

        decl nb311nf_innerk(%rsp)
        jz    _nb_kernel311nf_x86_64_sse.nb311nf_updateouterdata
        jmp   _nb_kernel311nf_x86_64_sse.nb311nf_odd_loop
_nb_kernel311nf_x86_64_sse.nb311nf_updateouterdata: 
        ## get n from stack
        movl nb311nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb311nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb311nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb311nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb311nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb311nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb311nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel311nf_x86_64_sse.nb311nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb311nf_n(%rsp)
        jmp _nb_kernel311nf_x86_64_sse.nb311nf_outer
_nb_kernel311nf_x86_64_sse.nb311nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb311nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel311nf_x86_64_sse.nb311nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel311nf_x86_64_sse.nb311nf_threadloop
_nb_kernel311nf_x86_64_sse.nb311nf_end: 
        movl nb311nf_nouter(%rsp),%eax
        movl nb311nf_ninner(%rsp),%ebx
        movq nb311nf_outeriter(%rbp),%rcx
        movq nb311nf_inneriter(%rbp),%rdx
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



