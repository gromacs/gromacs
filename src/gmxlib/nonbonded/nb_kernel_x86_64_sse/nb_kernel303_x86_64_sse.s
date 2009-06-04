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






.globl nb_kernel303_x86_64_sse
.globl _nb_kernel303_x86_64_sse
nb_kernel303_x86_64_sse:        
_nb_kernel303_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb303_fshift, 16
.set nb303_gid, 24
.set nb303_pos, 32
.set nb303_faction, 40
.set nb303_charge, 48
.set nb303_p_facel, 56
.set nb303_argkrf, 64
.set nb303_argcrf, 72
.set nb303_Vc, 80
.set nb303_type, 88
.set nb303_p_ntype, 96
.set nb303_vdwparam, 104
.set nb303_Vvdw, 112
.set nb303_p_tabscale, 120
.set nb303_VFtab, 128
.set nb303_invsqrta, 136
.set nb303_dvda, 144
.set nb303_p_gbtabscale, 152
.set nb303_GBtab, 160
.set nb303_p_nthreads, 168
.set nb303_count, 176
.set nb303_mtx, 184
.set nb303_outeriter, 192
.set nb303_inneriter, 200
.set nb303_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb303_ixH1, 0
.set nb303_iyH1, 16
.set nb303_izH1, 32
.set nb303_ixH2, 48
.set nb303_iyH2, 64
.set nb303_izH2, 80
.set nb303_ixM, 96
.set nb303_iyM, 112
.set nb303_izM, 128
.set nb303_iqH, 144
.set nb303_iqM, 160
.set nb303_dxH1, 176
.set nb303_dyH1, 192
.set nb303_dzH1, 208
.set nb303_dxH2, 224
.set nb303_dyH2, 240
.set nb303_dzH2, 256
.set nb303_dxM, 272
.set nb303_dyM, 288
.set nb303_dzM, 304
.set nb303_qqH, 320
.set nb303_qqM, 336
.set nb303_rinvH1, 352
.set nb303_rinvH2, 368
.set nb303_rinvM, 384
.set nb303_rH1, 400
.set nb303_rH2, 416
.set nb303_rM, 432
.set nb303_tsc, 448
.set nb303_two, 464
.set nb303_vctot, 480
.set nb303_fixH1, 496
.set nb303_fiyH1, 512
.set nb303_fizH1, 528
.set nb303_fixH2, 544
.set nb303_fiyH2, 560
.set nb303_fizH2, 576
.set nb303_fixM, 592
.set nb303_fiyM, 608
.set nb303_fizM, 624
.set nb303_epsH1, 640
.set nb303_epsH2, 656
.set nb303_epsM, 672
.set nb303_half, 688
.set nb303_three, 704
.set nb303_is3, 720
.set nb303_ii3, 724
.set nb303_nri, 728
.set nb303_iinr, 736
.set nb303_jindex, 744
.set nb303_jjnr, 752
.set nb303_shift, 760
.set nb303_shiftvec, 768
.set nb303_facel, 776
.set nb303_innerjjnr, 784
.set nb303_innerk, 792
.set nb303_n, 796
.set nb303_nn1, 800
.set nb303_nouter, 804
.set nb303_ninner, 808

        push %rbp
        movq %rsp,%rbp
        push %rbx


        push %r12
        push %r13
        push %r14
        push %r15

        subq $824,%rsp          ## local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb303_nouter(%rsp)
        movl %eax,nb303_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb303_nri(%rsp)
        movq %rsi,nb303_iinr(%rsp)
        movq %rdx,nb303_jindex(%rsp)
        movq %rcx,nb303_jjnr(%rsp)
        movq %r8,nb303_shift(%rsp)
        movq %r9,nb303_shiftvec(%rsp)
        movq nb303_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb303_facel(%rsp)

        movq nb303_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb303_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb303_half(%rsp)
        movss nb303_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb303_half(%rsp)
        movaps %xmm2,nb303_two(%rsp)
        movaps %xmm3,nb303_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb303_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb303_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb303_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb303_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb303_iqH(%rsp)
        movaps %xmm4,nb303_iqM(%rsp)

_nb_kernel303_x86_64_sse.nb303_threadloop: 
        movq  nb303_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel303_x86_64_sse.nb303_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303_x86_64_sse.nb303_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303_n(%rsp)
        movl %ebx,nb303_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303_x86_64_sse.nb303_outerstart
        jmp _nb_kernel303_x86_64_sse.nb303_end

_nb_kernel303_x86_64_sse.nb303_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303_nouter(%rsp),%ebx
        movl %ebx,nb303_nouter(%rsp)

_nb_kernel303_x86_64_sse.nb303_outer: 
        movq  nb303_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb303_is3(%rsp)      ## store is3 

        movq  nb303_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb303_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb303_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb303_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb303_ixH1(%rsp)
        movaps %xmm4,nb303_iyH1(%rsp)
        movaps %xmm5,nb303_izH1(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%rax,%rbx,4),%xmm0
        addss 28(%rax,%rbx,4),%xmm1
        addss 32(%rax,%rbx,4),%xmm2
        addss 36(%rax,%rbx,4),%xmm3
        addss 40(%rax,%rbx,4),%xmm4
        addss 44(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb303_ixH2(%rsp)
        movaps %xmm1,nb303_iyH2(%rsp)
        movaps %xmm2,nb303_izH2(%rsp)
        movaps %xmm3,nb303_ixM(%rsp)
        movaps %xmm4,nb303_iyM(%rsp)
        movaps %xmm5,nb303_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb303_vctot(%rsp)
        movaps %xmm4,nb303_fixH1(%rsp)
        movaps %xmm4,nb303_fiyH1(%rsp)
        movaps %xmm4,nb303_fizH1(%rsp)
        movaps %xmm4,nb303_fixH2(%rsp)
        movaps %xmm4,nb303_fiyH2(%rsp)
        movaps %xmm4,nb303_fizH2(%rsp)
        movaps %xmm4,nb303_fixM(%rsp)
        movaps %xmm4,nb303_fiyM(%rsp)
        movaps %xmm4,nb303_fizM(%rsp)

        movq  nb303_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb303_pos(%rbp),%rsi
        movq  nb303_faction(%rbp),%rdi
        movq  nb303_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb303_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb303_ninner(%rsp),%ecx
        movl  %ecx,nb303_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb303_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel303_x86_64_sse.nb303_unroll_loop
        jmp   _nb_kernel303_x86_64_sse.nb303_odd_inner
_nb_kernel303_x86_64_sse.nb303_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb303_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb303_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb303_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb303_iqH(%rsp),%xmm3
        mulps  nb303_iqM(%rsp),%xmm4

        movaps  %xmm3,nb303_qqH(%rsp)
        movaps  %xmm4,nb303_qqM(%rsp)

        movq nb303_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
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

    subps nb303_ixH1(%rsp),%xmm0
    subps nb303_iyH1(%rsp),%xmm1
    subps nb303_izH1(%rsp),%xmm2
    subps nb303_ixH2(%rsp),%xmm3
    subps nb303_iyH2(%rsp),%xmm4
    subps nb303_izH2(%rsp),%xmm5
    subps nb303_ixM(%rsp),%xmm6
    subps nb303_iyM(%rsp),%xmm7
    subps nb303_izM(%rsp),%xmm8

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm0,nb303_dxH1(%rsp)
        movaps %xmm1,nb303_dyH1(%rsp)
        movaps %xmm2,nb303_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb303_dxH2(%rsp)
        movaps %xmm4,nb303_dyH2(%rsp)
        movaps %xmm5,nb303_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb303_dxM(%rsp)
        movaps %xmm7,nb303_dyM(%rsp)
        movaps %xmm8,nb303_dzM(%rsp)
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

        movaps  nb303_three(%rsp),%xmm9
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

        movaps  nb303_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1 
        mulps   %xmm4,%xmm10 ## rinvH2
    mulps   %xmm4,%xmm11 ## rinvM

        movaps  %xmm9,nb303_rinvH1(%rsp)
        movaps  %xmm10,nb303_rinvH2(%rsp)
        movaps  %xmm11,nb303_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb303_tsc(%rsp),%xmm1
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

    movq nb303_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb303_epsH1(%rsp)
    movaps    %xmm3,nb303_epsH2(%rsp)
    movaps    %xmm6,nb303_epsM(%rsp)

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

    movaps nb303_epsH1(%rsp),%xmm12
    movaps nb303_epsH2(%rsp),%xmm13
    movaps nb303_epsM(%rsp),%xmm14

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
    movaps nb303_qqH(%rsp),%xmm12
    movaps nb303_qqM(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm12,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm12,%xmm7
    mulps  %xmm13,%xmm11

    ## accumulate vctot
    addps  nb303_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb303_vctot(%rsp)

    movaps nb303_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

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

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb303_rinvH1(%rsp),%xmm3
    mulps  nb303_rinvH2(%rsp),%xmm7
    mulps  nb303_rinvM(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb303_dxH1(%rsp),%xmm0
        mulps nb303_dyH1(%rsp),%xmm1
        mulps nb303_dzH1(%rsp),%xmm2
        mulps nb303_dxH2(%rsp),%xmm3
        mulps nb303_dyH2(%rsp),%xmm4
        mulps nb303_dzH2(%rsp),%xmm5
        mulps nb303_dxM(%rsp),%xmm6
        mulps nb303_dyM(%rsp),%xmm7
        mulps nb303_dzM(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb303_fixH1(%rsp),%xmm0
    addps nb303_fiyH1(%rsp),%xmm1
    addps nb303_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb303_fixH2(%rsp),%xmm3
    addps nb303_fiyH2(%rsp),%xmm4
    addps nb303_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb303_fixM(%rsp),%xmm6
    addps nb303_fiyM(%rsp),%xmm7
    addps nb303_fizM(%rsp),%xmm8

    movaps %xmm0,nb303_fixH1(%rsp)
    movaps %xmm1,nb303_fiyH1(%rsp)
    movaps %xmm2,nb303_fizH1(%rsp)
    movaps %xmm3,nb303_fixH2(%rsp)
    movaps %xmm4,nb303_fiyH2(%rsp)
    movaps %xmm5,nb303_fizH2(%rsp)
    movaps %xmm6,nb303_fixM(%rsp)
    movaps %xmm7,nb303_fiyM(%rsp)
    movaps %xmm8,nb303_fizM(%rsp)

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
        subl $4,nb303_innerk(%rsp)
        jl    _nb_kernel303_x86_64_sse.nb303_odd_inner
        jmp   _nb_kernel303_x86_64_sse.nb303_unroll_loop
_nb_kernel303_x86_64_sse.nb303_odd_inner: 
        addl $4,nb303_innerk(%rsp)
        jnz   _nb_kernel303_x86_64_sse.nb303_odd_loop
        jmp   _nb_kernel303_x86_64_sse.nb303_updateouterdata
_nb_kernel303_x86_64_sse.nb303_odd_loop: 
        movq  nb303_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb303_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb303_iqM(%rsp),%xmm4
        movq nb303_charge(%rbp),%rsi
        movhps nb303_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb303_qqM(%rsp)    ## use dummy qq for storage 

        movq nb303_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb303_ixM(%rsp),%xmm0
        movss nb303_iyM(%rsp),%xmm1
        movss nb303_izM(%rsp),%xmm2

        movlps nb303_ixH1(%rsp),%xmm6
        movlps nb303_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb303_iyH1(%rsp),%xmm6
        movlps nb303_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb303_izH1(%rsp),%xmm6
        movlps nb303_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use dummy dx for storage
        movaps %xmm3,nb303_dxM(%rsp)
        movaps %xmm4,nb303_dyM(%rsp)
        movaps %xmm5,nb303_dzM(%rsp)

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
        movaps nb303_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb303_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb303_rinvM(%rsp)

        mulps nb303_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
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

        movq nb303_VFtab(%rbp),%rsi
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
        mulps  nb303_two(%rsp),%xmm7            ## two*Heps2 
        movaps nb303_qqM(%rsp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        addps  nb303_vctot(%rsp),%xmm5
        movaps %xmm5,nb303_vctot(%rsp)

        xorps %xmm4,%xmm4
        mulps  nb303_tsc(%rsp),%xmm0
        mulps  nb303_rinvM(%rsp),%xmm0
        subps  %xmm0,%xmm4

        movd %mm0,%eax
        movd %mm1,%ecx
        movd %mm2,%edx

        movaps nb303_dxM(%rsp),%xmm0
        movaps nb303_dyM(%rsp),%xmm1
        movaps nb303_dzM(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force) 
        movss  nb303_fixM(%rsp),%xmm3
        movss  nb303_fiyM(%rsp),%xmm4
        movss  nb303_fizM(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb303_fixM(%rsp)
        movss  %xmm4,nb303_fiyM(%rsp)
        movss  %xmm5,nb303_fizM(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb303_fixH1(%rsp),%xmm3
        addss  nb303_fiyH1(%rsp),%xmm4
        addss  nb303_fizH1(%rsp),%xmm5
        movss  %xmm3,nb303_fixH1(%rsp)
        movss  %xmm4,nb303_fiyH1(%rsp)
        movss  %xmm5,nb303_fizH1(%rsp)          ## updated the H1 force 

        movq nb303_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb303_fixH2(%rsp),%xmm3
        addss  nb303_fiyH2(%rsp),%xmm4
        addss  nb303_fizH2(%rsp),%xmm5
        movss  %xmm3,nb303_fixH2(%rsp)
        movss  %xmm4,nb303_fiyH2(%rsp)
        movss  %xmm5,nb303_fizH2(%rsp)          ## updated the H2 force 

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
        addss    %xmm1,%xmm2    ## z sum in xmm2 
        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb303_innerk(%rsp)
        jz    _nb_kernel303_x86_64_sse.nb303_updateouterdata
        jmp   _nb_kernel303_x86_64_sse.nb303_odd_loop
_nb_kernel303_x86_64_sse.nb303_updateouterdata: 
        movl  nb303_ii3(%rsp),%ecx
        movq  nb303_faction(%rbp),%rdi
        movq  nb303_fshift(%rbp),%rsi
        movl  nb303_is3(%rsp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixH1(%rsp),%xmm0
        movaps nb303_fiyH1(%rsp),%xmm1
        movaps nb303_fizH1(%rsp),%xmm2

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
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixH2(%rsp),%xmm0
        movaps nb303_fiyH2(%rsp),%xmm1
        movaps nb303_fizH2(%rsp),%xmm2

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

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb303_fixM(%rsp),%xmm0
        movaps nb303_fiyM(%rsp),%xmm1
        movaps nb303_fizM(%rsp),%xmm2

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
        movl nb303_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb303_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb303_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb303_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb303_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303_x86_64_sse.nb303_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303_n(%rsp)
        jmp _nb_kernel303_x86_64_sse.nb303_outer
_nb_kernel303_x86_64_sse.nb303_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303_x86_64_sse.nb303_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303_x86_64_sse.nb303_threadloop
_nb_kernel303_x86_64_sse.nb303_end: 
        movl nb303_nouter(%rsp),%eax
        movl nb303_ninner(%rsp),%ebx
        movq nb303_outeriter(%rbp),%rcx
        movq nb303_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $824,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel303nf_x86_64_sse
.globl _nb_kernel303nf_x86_64_sse
nb_kernel303nf_x86_64_sse:      
_nb_kernel303nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb303nf_fshift, 16
.set nb303nf_gid, 24
.set nb303nf_pos, 32
.set nb303nf_faction, 40
.set nb303nf_charge, 48
.set nb303nf_p_facel, 56
.set nb303nf_argkrf, 64
.set nb303nf_argcrf, 72
.set nb303nf_Vc, 80
.set nb303nf_type, 88
.set nb303nf_p_ntype, 96
.set nb303nf_vdwparam, 104
.set nb303nf_Vvdw, 112
.set nb303nf_p_tabscale, 120
.set nb303nf_VFtab, 128
.set nb303nf_invsqrta, 136
.set nb303nf_dvda, 144
.set nb303nf_p_gbtabscale, 152
.set nb303nf_GBtab, 160
.set nb303nf_p_nthreads, 168
.set nb303nf_count, 176
.set nb303nf_mtx, 184
.set nb303nf_outeriter, 192
.set nb303nf_inneriter, 200
.set nb303nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb303nf_ixH1, 0
.set nb303nf_iyH1, 16
.set nb303nf_izH1, 32
.set nb303nf_ixH2, 48
.set nb303nf_iyH2, 64
.set nb303nf_izH2, 80
.set nb303nf_ixM, 96
.set nb303nf_iyM, 112
.set nb303nf_izM, 128
.set nb303nf_iqH, 144
.set nb303nf_iqM, 160
.set nb303nf_qqH, 176
.set nb303nf_qqM, 192
.set nb303nf_rinvH1, 208
.set nb303nf_rinvH2, 224
.set nb303nf_rinvM, 240
.set nb303nf_rH1, 256
.set nb303nf_rH2, 272
.set nb303nf_rM, 288
.set nb303nf_tsc, 304
.set nb303nf_vctot, 320
.set nb303nf_half, 336
.set nb303nf_three, 352
.set nb303nf_is3, 368
.set nb303nf_ii3, 372
.set nb303nf_nri, 376
.set nb303nf_iinr, 384
.set nb303nf_jindex, 392
.set nb303nf_jjnr, 400
.set nb303nf_shift, 408
.set nb303nf_shiftvec, 416
.set nb303nf_facel, 424
.set nb303nf_innerjjnr, 432
.set nb303nf_innerk, 440
.set nb303nf_n, 444
.set nb303nf_nn1, 448
.set nb303nf_nouter, 452
.set nb303nf_ninner, 456
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $472,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb303nf_nouter(%rsp)
        movl %eax,nb303nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb303nf_nri(%rsp)
        movq %rsi,nb303nf_iinr(%rsp)
        movq %rdx,nb303nf_jindex(%rsp)
        movq %rcx,nb303nf_jjnr(%rsp)
        movq %r8,nb303nf_shift(%rsp)
        movq %r9,nb303nf_shiftvec(%rsp)
        movq nb303nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb303nf_facel(%rsp)

        movq nb303nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb303nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb303nf_half(%rsp)
        movss nb303nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb303nf_half(%rsp)
        movaps %xmm3,nb303nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb303nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb303nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb303nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb303nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb303nf_iqH(%rsp)
        movaps %xmm4,nb303nf_iqM(%rsp)

_nb_kernel303nf_x86_64_sse.nb303nf_threadloop: 
        movq  nb303nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel303nf_x86_64_sse.nb303nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel303nf_x86_64_sse.nb303nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb303nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb303nf_n(%rsp)
        movl %ebx,nb303nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel303nf_x86_64_sse.nb303nf_outerstart
        jmp _nb_kernel303nf_x86_64_sse.nb303nf_end

_nb_kernel303nf_x86_64_sse.nb303nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb303nf_nouter(%rsp),%ebx
        movl %ebx,nb303nf_nouter(%rsp)

_nb_kernel303nf_x86_64_sse.nb303nf_outer: 
        movq  nb303nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb303nf_is3(%rsp)            ## store is3 

        movq  nb303nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb303nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb303nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb303nf_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb303nf_ixH1(%rsp)
        movaps %xmm4,nb303nf_iyH1(%rsp)
        movaps %xmm5,nb303nf_izH1(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%rax,%rbx,4),%xmm0
        addss 28(%rax,%rbx,4),%xmm1
        addss 32(%rax,%rbx,4),%xmm2
        addss 36(%rax,%rbx,4),%xmm3
        addss 40(%rax,%rbx,4),%xmm4
        addss 44(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb303nf_ixH2(%rsp)
        movaps %xmm1,nb303nf_iyH2(%rsp)
        movaps %xmm2,nb303nf_izH2(%rsp)
        movaps %xmm3,nb303nf_ixM(%rsp)
        movaps %xmm4,nb303nf_iyM(%rsp)
        movaps %xmm5,nb303nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb303nf_vctot(%rsp)

        movq  nb303nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb303nf_pos(%rbp),%rsi
        movq  nb303nf_faction(%rbp),%rdi
        movq  nb303nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb303nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb303nf_ninner(%rsp),%ecx
        movl  %ecx,nb303nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb303nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel303nf_x86_64_sse.nb303nf_unroll_loop
        jmp   _nb_kernel303nf_x86_64_sse.nb303nf_odd_inner
_nb_kernel303nf_x86_64_sse.nb303nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb303nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb303nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb303nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb303nf_iqH(%rsp),%xmm3
        mulps  nb303nf_iqM(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb303nf_qqH(%rsp)
        movaps  %xmm4,nb303nf_qqM(%rsp)

        movq nb303nf_pos(%rbp),%rsi     ## base of pos[] 

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

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb303nf_ixH1(%rsp),%xmm4
        movaps nb303nf_iyH1(%rsp),%xmm5
        movaps nb303nf_izH1(%rsp),%xmm6

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
        ## rsqH1 in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb303nf_ixH2(%rsp),%xmm4
        movaps nb303nf_iyH2(%rsp),%xmm5
        movaps nb303nf_izH2(%rsp),%xmm6

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
        ## rsqH2 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb303nf_ixM(%rsp),%xmm3
        movaps nb303nf_iyM(%rsp),%xmm4
        movaps nb303nf_izM(%rsp),%xmm5

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
        ## rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

        ## start with rsqH1 - seed to xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%rsp),%xmm4
        movaps  %xmm4,nb303nf_rinvH1(%rsp)      ## rinvH1 in xmm4 
        mulps   %xmm4,%xmm7
        movaps  %xmm7,nb303nf_rH1(%rsp)

        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%rsp),%xmm4
        movaps  %xmm4,nb303nf_rinvH2(%rsp)      ## rinvH2 in xmm4 
        mulps   %xmm4,%xmm6
        movaps  %xmm6,nb303nf_rH2(%rsp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb303nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb303nf_half(%rsp),%xmm4
        movaps  %xmm4,nb303nf_rinvM(%rsp)       ## rinvM in xmm4 
        mulps   %xmm4,%xmm5
        movaps  %xmm5,nb303nf_rM(%rsp)

        ## do H1 interactions 
        ## rH1 is still in xmm7 
        mulps   nb303nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

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

        movq nb303nf_VFtab(%rbp),%rsi
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
        movaps nb303nf_qqH(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  

        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb303nf_vctot(%rsp)

        ## Done with H1 interactions - now H2! 
        movaps nb303nf_rH2(%rsp),%xmm7
        mulps   nb303nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

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
        movaps nb303nf_qqH(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb303nf_vctot(%rsp)

        ## Done with H2, finally we do M interactions 
        movaps nb303nf_rM(%rsp),%xmm7
        mulps   nb303nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

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
        movaps nb303nf_qqM(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb303nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb303nf_innerk(%rsp)
        jl    _nb_kernel303nf_x86_64_sse.nb303nf_odd_inner
        jmp   _nb_kernel303nf_x86_64_sse.nb303nf_unroll_loop
_nb_kernel303nf_x86_64_sse.nb303nf_odd_inner: 
        addl $4,nb303nf_innerk(%rsp)
        jnz   _nb_kernel303nf_x86_64_sse.nb303nf_odd_loop
        jmp   _nb_kernel303nf_x86_64_sse.nb303nf_updateouterdata
_nb_kernel303nf_x86_64_sse.nb303nf_odd_loop: 
        movq  nb303nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb303nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb303nf_iqM(%rsp),%xmm4
        movq nb303nf_charge(%rbp),%rsi
        movhps nb303nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb303nf_qqM(%rsp)          ## use dummy qq for storage 

        movq nb303nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb303nf_ixM(%rsp),%xmm3
        movss nb303nf_iyM(%rsp),%xmm4
        movss nb303nf_izM(%rsp),%xmm5

        movlps nb303nf_ixH1(%rsp),%xmm6
        movlps nb303nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb303nf_iyH1(%rsp),%xmm6
        movlps nb303nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb303nf_izH1(%rsp),%xmm6
        movlps nb303nf_izH2(%rsp),%xmm7
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
        movaps nb303nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb303nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        mulps %xmm0,%xmm4       ## xmm4=r 
        movaps %xmm0,nb303nf_rinvM(%rsp)

        mulps nb303nf_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
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

        movq nb303nf_VFtab(%rbp),%rsi
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
        movaps nb303nf_qqM(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb303nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb303nf_vctot(%rsp)

        decl nb303nf_innerk(%rsp)
        jz    _nb_kernel303nf_x86_64_sse.nb303nf_updateouterdata
        jmp   _nb_kernel303nf_x86_64_sse.nb303nf_odd_loop
_nb_kernel303nf_x86_64_sse.nb303nf_updateouterdata: 
        ## get n from stack
        movl nb303nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb303nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb303nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb303nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb303nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel303nf_x86_64_sse.nb303nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb303nf_n(%rsp)
        jmp _nb_kernel303nf_x86_64_sse.nb303nf_outer
_nb_kernel303nf_x86_64_sse.nb303nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb303nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel303nf_x86_64_sse.nb303nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel303nf_x86_64_sse.nb303nf_threadloop
_nb_kernel303nf_x86_64_sse.nb303nf_end: 
        movl nb303nf_nouter(%rsp),%eax
        movl nb303nf_ninner(%rsp),%ebx
        movq nb303nf_outeriter(%rbp),%rcx
        movq nb303nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $472,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

