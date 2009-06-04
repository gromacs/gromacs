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







.globl nb_kernel304_x86_64_sse
.globl _nb_kernel304_x86_64_sse
nb_kernel304_x86_64_sse:        
_nb_kernel304_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb304_fshift, 16
.set nb304_gid, 24
.set nb304_pos, 32
.set nb304_faction, 40
.set nb304_charge, 48
.set nb304_p_facel, 56
.set nb304_argkrf, 64
.set nb304_argcrf, 72
.set nb304_Vc, 80
.set nb304_type, 88
.set nb304_p_ntype, 96
.set nb304_vdwparam, 104
.set nb304_Vvdw, 112
.set nb304_p_tabscale, 120
.set nb304_VFtab, 128
.set nb304_invsqrta, 136
.set nb304_dvda, 144
.set nb304_p_gbtabscale, 152
.set nb304_GBtab, 160
.set nb304_p_nthreads, 168
.set nb304_count, 176
.set nb304_mtx, 184
.set nb304_outeriter, 192
.set nb304_inneriter, 200
.set nb304_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb304_ixH1, 0
.set nb304_iyH1, 16
.set nb304_izH1, 32
.set nb304_ixH2, 48
.set nb304_iyH2, 64
.set nb304_izH2, 80
.set nb304_ixM, 96
.set nb304_iyM, 112
.set nb304_izM, 128
.set nb304_jxH1, 144
.set nb304_jyH1, 160
.set nb304_jzH1, 176
.set nb304_jxH2, 192
.set nb304_jyH2, 208
.set nb304_jzH2, 224
.set nb304_jxM, 240
.set nb304_jyM, 256
.set nb304_jzM, 272
.set nb304_dxH1H1, 288
.set nb304_dyH1H1, 304
.set nb304_dzH1H1, 320
.set nb304_dxH1H2, 336
.set nb304_dyH1H2, 352
.set nb304_dzH1H2, 368
.set nb304_dxH1M, 384
.set nb304_dyH1M, 400
.set nb304_dzH1M, 416
.set nb304_dxH2H1, 432
.set nb304_dyH2H1, 448
.set nb304_dzH2H1, 464
.set nb304_dxH2H2, 480
.set nb304_dyH2H2, 496
.set nb304_dzH2H2, 512
.set nb304_dxH2M, 528
.set nb304_dyH2M, 544
.set nb304_dzH2M, 560
.set nb304_dxMH1, 576
.set nb304_dyMH1, 592
.set nb304_dzMH1, 608
.set nb304_dxMH2, 624
.set nb304_dyMH2, 640
.set nb304_dzMH2, 656
.set nb304_dxMM, 672
.set nb304_dyMM, 688
.set nb304_dzMM, 704
.set nb304_qqHH, 720
.set nb304_qqMH, 736
.set nb304_qqMM, 752
.set nb304_two, 768
.set nb304_tsc, 784
.set nb304_vctot, 800
.set nb304_fixH1, 816
.set nb304_fiyH1, 832
.set nb304_fizH1, 848
.set nb304_fixH2, 864
.set nb304_fiyH2, 880
.set nb304_fizH2, 896
.set nb304_fixM, 912
.set nb304_fiyM, 928
.set nb304_fizM, 944
.set nb304_epsH1, 960
.set nb304_epsH2, 976
.set nb304_epsM, 992
.set nb304_fjxH2, 1008
.set nb304_fjyH2, 1024
.set nb304_fjzH2, 1040
.set nb304_fjxM, 1056
.set nb304_fjyM, 1072
.set nb304_fjzM, 1088
.set nb304_half, 1104
.set nb304_three, 1120
.set nb304_rsqH1H1, 1136
.set nb304_rsqH1H2, 1152
.set nb304_rsqH1M, 1168
.set nb304_rsqH2H1, 1184
.set nb304_rsqH2H2, 1200
.set nb304_rsqH2M, 1216
.set nb304_rsqMH1, 1232
.set nb304_rsqMH2, 1248
.set nb304_rsqMM, 1264
.set nb304_rinvH1H1, 1280
.set nb304_rinvH1H2, 1296
.set nb304_rinvH1M, 1312
.set nb304_rinvH2H1, 1328
.set nb304_rinvH2H2, 1344
.set nb304_rinvH2M, 1360
.set nb304_rinvMH1, 1376
.set nb304_rinvMH2, 1392
.set nb304_rinvMM, 1408
.set nb304_is3, 1424
.set nb304_ii3, 1428
.set nb304_nri, 1432
.set nb304_iinr, 1440
.set nb304_jindex, 1448
.set nb304_jjnr, 1456
.set nb304_shift, 1464
.set nb304_shiftvec, 1472
.set nb304_facel, 1480
.set nb304_innerjjnr, 1488
.set nb304_innerk, 1496
.set nb304_n, 1500
.set nb304_nn1, 1504
.set nb304_nouter, 1508
.set nb304_ninner, 1512
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1528,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb304_nouter(%rsp)
        movl %eax,nb304_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb304_nri(%rsp)
        movq %rsi,nb304_iinr(%rsp)
        movq %rdx,nb304_jindex(%rsp)
        movq %rcx,nb304_jjnr(%rsp)
        movq %r8,nb304_shift(%rsp)
        movq %r9,nb304_shiftvec(%rsp)
        movq nb304_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb304_facel(%rsp)

        movq nb304_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb304_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb304_half(%rsp)
        movss nb304_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb304_half(%rsp)
        movaps %xmm2,nb304_two(%rsp)
        movaps %xmm3,nb304_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb304_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb304_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb304_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb304_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304_qqHH(%rsp)
        movaps %xmm4,nb304_qqMH(%rsp)
        movaps %xmm5,nb304_qqMM(%rsp)

_nb_kernel304_x86_64_sse.nb304_threadloop: 
        movq  nb304_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel304_x86_64_sse.nb304_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304_x86_64_sse.nb304_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304_n(%rsp)
        movl %ebx,nb304_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304_x86_64_sse.nb304_outerstart
        jmp _nb_kernel304_x86_64_sse.nb304_end

_nb_kernel304_x86_64_sse.nb304_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304_nouter(%rsp),%ebx
        movl %ebx,nb304_nouter(%rsp)

_nb_kernel304_x86_64_sse.nb304_outer: 
        movq  nb304_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb304_is3(%rsp)      ## store is3 

        movq  nb304_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb304_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb304_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb304_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304_ixH1(%rsp)
        movaps %xmm4,nb304_iyH1(%rsp)
        movaps %xmm5,nb304_izH1(%rsp)

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
        movaps %xmm0,nb304_ixH2(%rsp)
        movaps %xmm1,nb304_iyH2(%rsp)
        movaps %xmm2,nb304_izH2(%rsp)
        movaps %xmm3,nb304_ixM(%rsp)
        movaps %xmm4,nb304_iyM(%rsp)
        movaps %xmm5,nb304_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb304_vctot(%rsp)
        movaps %xmm4,nb304_fixH1(%rsp)
        movaps %xmm4,nb304_fiyH1(%rsp)
        movaps %xmm4,nb304_fizH1(%rsp)
        movaps %xmm4,nb304_fixH2(%rsp)
        movaps %xmm4,nb304_fiyH2(%rsp)
        movaps %xmm4,nb304_fizH2(%rsp)
        movaps %xmm4,nb304_fixM(%rsp)
        movaps %xmm4,nb304_fiyM(%rsp)
        movaps %xmm4,nb304_fizM(%rsp)

        movq  nb304_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb304_pos(%rbp),%rsi
        movq  nb304_faction(%rbp),%rdi
        movq  nb304_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb304_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb304_ninner(%rsp),%ecx
        movl  %ecx,nb304_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb304_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel304_x86_64_sse.nb304_unroll_loop
        jmp   _nb_kernel304_x86_64_sse.nb304_single_check
_nb_kernel304_x86_64_sse.nb304_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb304_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb304_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb304_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j H1 coordinates to local temp variables 
    movlps 12(%rsi,%rax,4),%xmm0    ## jxH1a jyH1a  -   -
    movlps 12(%rsi,%rcx,4),%xmm1    ## jxH1c jyH1c  -   -
    movhps 12(%rsi,%rbx,4),%xmm0    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rsi,%rdx,4),%xmm1    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rsi,%rax,4),%xmm2    ## jzH1a  -  -  -
    movss  20(%rsi,%rcx,4),%xmm3    ## jzH1c  -  -  -
    movhps 20(%rsi,%rbx,4),%xmm2    ## jzH1a  -  jzH1b  -
    movhps 20(%rsi,%rdx,4),%xmm3    ## jzH1c  -  jzH1d -

    movd %eax,%mm0 ## save j3 in mm0-mm3
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH1a jxH1c jyH1a jyH1c        
    unpckhps %xmm1,%xmm4 ## jxH1b jxH1d jyH1b jyH1d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps  $136,%xmm3,%xmm2  ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb304_ixH1(%rsp),%xmm0
    subps nb304_iyH1(%rsp),%xmm1
    subps nb304_izH1(%rsp),%xmm2
    subps nb304_ixH2(%rsp),%xmm3
    subps nb304_iyH2(%rsp),%xmm4
    subps nb304_izH2(%rsp),%xmm5
    subps nb304_ixM(%rsp),%xmm6
    subps nb304_iyM(%rsp),%xmm7
    subps nb304_izM(%rsp),%xmm8

        movaps %xmm0,nb304_dxH1H1(%rsp)
        movaps %xmm1,nb304_dyH1H1(%rsp)
        movaps %xmm2,nb304_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH2H1(%rsp)
        movaps %xmm4,nb304_dyH2H1(%rsp)
        movaps %xmm5,nb304_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb304_dxMH1(%rsp)
        movaps %xmm7,nb304_dyMH1(%rsp)
        movaps %xmm8,nb304_dzMH1(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb304_three(%rsp),%xmm9
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

        movaps  nb304_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H1 
        mulps   %xmm4,%xmm10 ## rinvH2H1
    mulps   %xmm4,%xmm11 ## rinvMH1

        movaps  %xmm9,nb304_rinvH1H1(%rsp)
        movaps  %xmm10,nb304_rinvH2H1(%rsp)
        movaps  %xmm11,nb304_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb304_tsc(%rsp),%xmm1
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

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb304_epsH1(%rsp)
    movaps    %xmm3,nb304_epsH2(%rsp)
    movaps    %xmm6,nb304_epsM(%rsp)

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

    movaps nb304_epsH1(%rsp),%xmm12
    movaps nb304_epsH2(%rsp),%xmm13
    movaps nb304_epsM(%rsp),%xmm14

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
    movaps nb304_qqHH(%rsp),%xmm12
    movaps nb304_qqMH(%rsp),%xmm13
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
    addps  nb304_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb304_vctot(%rsp)

    movaps nb304_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j H1 forces to local temp variables 
    movlps 12(%rdi,%rax,4),%xmm11    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm12    ## jxH1c jyH1c  -   -
    movhps 12(%rdi,%rbx,4),%xmm11    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm12    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rdi,%rax,4),%xmm13    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm14    ## jzH1c  -  -  -
    movhps 20(%rdi,%rbx,4),%xmm13    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm14    ## jzH1c  -  jzH1d -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm11: jxH1a jyH1a jxH1b jyH1b 
    ## xmm12: jxH1c jyH1c jxH1d jyH1d
    ## xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb304_rinvH1H1(%rsp),%xmm3
    mulps  nb304_rinvH2H1(%rsp),%xmm7
    mulps  nb304_rinvMH1(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb304_dxH1H1(%rsp),%xmm0
        mulps nb304_dyH1H1(%rsp),%xmm1
        mulps nb304_dzH1H1(%rsp),%xmm2
        mulps nb304_dxH2H1(%rsp),%xmm3
        mulps nb304_dyH2H1(%rsp),%xmm4
        mulps nb304_dzH2H1(%rsp),%xmm5
        mulps nb304_dxMH1(%rsp),%xmm6
        mulps nb304_dyMH1(%rsp),%xmm7
        mulps nb304_dzMH1(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb304_fixH1(%rsp),%xmm0
    addps nb304_fiyH1(%rsp),%xmm1
    addps nb304_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb304_fixH2(%rsp),%xmm3
    addps nb304_fiyH2(%rsp),%xmm4
    addps nb304_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb304_fixM(%rsp),%xmm6
    addps nb304_fiyM(%rsp),%xmm7
    addps nb304_fizM(%rsp),%xmm8

    movaps %xmm0,nb304_fixH1(%rsp)
    movaps %xmm1,nb304_fiyH1(%rsp)
    movaps %xmm2,nb304_fizH1(%rsp)
    movaps %xmm3,nb304_fixH2(%rsp)
    movaps %xmm4,nb304_fiyH2(%rsp)
    movaps %xmm5,nb304_fizH2(%rsp)
    movaps %xmm6,nb304_fixM(%rsp)
    movaps %xmm7,nb304_fiyM(%rsp)
    movaps %xmm8,nb304_fizM(%rsp)

    ## xmm14 = fH1x
    ## xmm15 = fH1y
    ## xmm13 = fH1z
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fH1zc fH1zd

    movlps %xmm11,12(%rdi,%rax,4)
    movhps %xmm11,12(%rdi,%rbx,4)
    movlps %xmm12,12(%rdi,%rcx,4)
    movhps %xmm12,12(%rdi,%rdx,4)
    movss  %xmm13,20(%rdi,%rax,4)
    movss  %xmm14,20(%rdi,%rcx,4)
    shufps $1,%xmm13,%xmm13
    shufps $1,%xmm14,%xmm14
    movss  %xmm13,20(%rdi,%rbx,4)
    movss  %xmm14,20(%rdi,%rdx,4)

        ## move j H2 coordinates to local temp variables 
        movq  nb304_pos(%rbp),%rsi
    movlps 24(%rsi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rsi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rsi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rsi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rsi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rsi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movhps 32(%rsi,%rbx,4),%xmm2    ## jzH2a  -  jzH2b  -
    movhps 32(%rsi,%rdx,4),%xmm3    ## jzH2c  -  jzH2d -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH2a jxH2c jyH2a jyH2c        
    unpckhps %xmm1,%xmm4 ## jxH2b jxH2d jyH2b jyH2d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps  $136,%xmm3,%xmm2  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb304_ixH1(%rsp),%xmm0
    subps nb304_iyH1(%rsp),%xmm1
    subps nb304_izH1(%rsp),%xmm2
    subps nb304_ixH2(%rsp),%xmm3
    subps nb304_iyH2(%rsp),%xmm4
    subps nb304_izH2(%rsp),%xmm5
    subps nb304_ixM(%rsp),%xmm6
    subps nb304_iyM(%rsp),%xmm7
    subps nb304_izM(%rsp),%xmm8

        movaps %xmm0,nb304_dxH1H2(%rsp)
        movaps %xmm1,nb304_dyH1H2(%rsp)
        movaps %xmm2,nb304_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH2H2(%rsp)
        movaps %xmm4,nb304_dyH2H2(%rsp)
        movaps %xmm5,nb304_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb304_dxMH2(%rsp)
        movaps %xmm7,nb304_dyMH2(%rsp)
        movaps %xmm8,nb304_dzMH2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb304_three(%rsp),%xmm9
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

        movaps  nb304_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H2
        mulps   %xmm4,%xmm10 ## rinvH2H2
    mulps   %xmm4,%xmm11 ## rinvMH2

        movaps  %xmm9,nb304_rinvH1H2(%rsp)
        movaps  %xmm10,nb304_rinvH2H2(%rsp)
        movaps  %xmm11,nb304_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb304_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    movq nb304_VFtab(%rbp),%rsi

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

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb304_epsH1(%rsp)
    movaps    %xmm3,nb304_epsH2(%rsp)
    movaps    %xmm6,nb304_epsM(%rsp)


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

    movaps nb304_epsH1(%rsp),%xmm12
    movaps nb304_epsH2(%rsp),%xmm13
    movaps nb304_epsM(%rsp),%xmm14

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
    movaps nb304_qqHH(%rsp),%xmm12
    movaps nb304_qqMH(%rsp),%xmm13
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
    addps  nb304_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb304_vctot(%rsp)

    movaps nb304_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm11    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm12    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm11    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm12    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm13    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm14    ## jzH2c  -  -  -
    movhps 32(%rdi,%rbx,4),%xmm13    ## jzH2a  -  jzH2b  -
    movhps 32(%rdi,%rdx,4),%xmm14    ## jzH2c  -  jzH2d -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm11: jxH2a jyH2a jxH2b jyH2b 
    ## xmm12: jxH2c jyH2c jxH2d jyH2d
    ## xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb304_rinvH1H2(%rsp),%xmm3
    mulps  nb304_rinvH2H2(%rsp),%xmm7
    mulps  nb304_rinvMH2(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb304_dxH1H2(%rsp),%xmm0
        mulps nb304_dyH1H2(%rsp),%xmm1
        mulps nb304_dzH1H2(%rsp),%xmm2
        mulps nb304_dxH2H2(%rsp),%xmm3
        mulps nb304_dyH2H2(%rsp),%xmm4
        mulps nb304_dzH2H2(%rsp),%xmm5
        mulps nb304_dxMH2(%rsp),%xmm6
        mulps nb304_dyMH2(%rsp),%xmm7
        mulps nb304_dzMH2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb304_fixH1(%rsp),%xmm0
    addps nb304_fiyH1(%rsp),%xmm1
    addps nb304_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb304_fixH2(%rsp),%xmm3
    addps nb304_fiyH2(%rsp),%xmm4
    addps nb304_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb304_fixM(%rsp),%xmm6
    addps nb304_fiyM(%rsp),%xmm7
    addps nb304_fizM(%rsp),%xmm8

    movaps %xmm0,nb304_fixH1(%rsp)
    movaps %xmm1,nb304_fiyH1(%rsp)
    movaps %xmm2,nb304_fizH1(%rsp)
    movaps %xmm3,nb304_fixH2(%rsp)
    movaps %xmm4,nb304_fiyH2(%rsp)
    movaps %xmm5,nb304_fizH2(%rsp)
    movaps %xmm6,nb304_fixM(%rsp)
    movaps %xmm7,nb304_fiyM(%rsp)
    movaps %xmm8,nb304_fizM(%rsp)

    ## xmm14 = fH2x
    ## xmm15 = fH2y
    ## xmm13 = fH2z
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fH2zc fH2zd

    movlps %xmm11,24(%rdi,%rax,4)
    movhps %xmm11,24(%rdi,%rbx,4)
    movlps %xmm12,24(%rdi,%rcx,4)
    movhps %xmm12,24(%rdi,%rdx,4)
    movss  %xmm13,32(%rdi,%rax,4)
    movss  %xmm14,32(%rdi,%rcx,4)
    shufps $1,%xmm13,%xmm13
    shufps $1,%xmm14,%xmm14
    movss  %xmm13,32(%rdi,%rbx,4)
    movss  %xmm14,32(%rdi,%rdx,4)

        movq  nb304_pos(%rbp),%rsi
        ## move j M coordinates to local temp variables 
    movlps 36(%rsi,%rax,4),%xmm0    ## jxMa jyMa  -   -
    movlps 36(%rsi,%rcx,4),%xmm1    ## jxMc jyMc  -   -
    movhps 36(%rsi,%rbx,4),%xmm0    ## jxMa jyMa jxMb jyMb 
    movhps 36(%rsi,%rdx,4),%xmm1    ## jxMc jyMc jxMd jyMd 

    movss  44(%rsi,%rax,4),%xmm2    ## jzMa  -  -  -
    movss  44(%rsi,%rcx,4),%xmm3    ## jzMc  -  -  -
    movss  44(%rsi,%rbx,4),%xmm5    ## jzMb  -  -  -
    movss  44(%rsi,%rdx,4),%xmm6    ## jzMd  -  -  -
    movlhps %xmm5,%xmm2 ## jzMa  -  jzMb  -
    movlhps %xmm6,%xmm3 ## jzMc  -  jzMd -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxMa jxMc jyMa jyMc        
    unpckhps %xmm1,%xmm4 ## jxMb jxMd jyMb jyMd
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps  $136,%xmm3,%xmm2  ## 10001000 => jzMa jzMb jzMc jzMd

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb304_ixH1(%rsp),%xmm0
    subps nb304_iyH1(%rsp),%xmm1
    subps nb304_izH1(%rsp),%xmm2
    subps nb304_ixH2(%rsp),%xmm3
    subps nb304_iyH2(%rsp),%xmm4
    subps nb304_izH2(%rsp),%xmm5
    subps nb304_ixM(%rsp),%xmm6
    subps nb304_iyM(%rsp),%xmm7
    subps nb304_izM(%rsp),%xmm8

        movaps %xmm0,nb304_dxH1M(%rsp)
        movaps %xmm1,nb304_dyH1M(%rsp)
        movaps %xmm2,nb304_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb304_dxH2M(%rsp)
        movaps %xmm4,nb304_dyH2M(%rsp)
        movaps %xmm5,nb304_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb304_dxMM(%rsp)
        movaps %xmm7,nb304_dyMM(%rsp)
        movaps %xmm8,nb304_dzMM(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb304_three(%rsp),%xmm9
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

        movaps  nb304_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1M
        mulps   %xmm4,%xmm10 ## rinvH2M
    mulps   %xmm4,%xmm11 ## rinvMM

        movaps  %xmm9,nb304_rinvH1M(%rsp)
        movaps  %xmm10,nb304_rinvH2M(%rsp)
        movaps  %xmm11,nb304_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb304_tsc(%rsp),%xmm1
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

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb304_epsH1(%rsp)
    movaps    %xmm3,nb304_epsH2(%rsp)
    movaps    %xmm6,nb304_epsM(%rsp)

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

    movaps nb304_epsH1(%rsp),%xmm12
    movaps nb304_epsH2(%rsp),%xmm13
    movaps nb304_epsM(%rsp),%xmm14

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
    movaps nb304_qqMH(%rsp),%xmm12
    movaps nb304_qqMM(%rsp),%xmm13
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
    addps  nb304_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb304_vctot(%rsp)

    movaps nb304_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j M forces to local temp variables 
    movlps 36(%rdi,%rax,4),%xmm11    ## jxMa jyMa  -   -
    movlps 36(%rdi,%rcx,4),%xmm12    ## jxMc jyMc  -   -
    movhps 36(%rdi,%rbx,4),%xmm11    ## jxMa jyMa jxMb jyMb 
    movhps 36(%rdi,%rdx,4),%xmm12    ## jxMc jyMc jxMd jyMd 

    movss  44(%rdi,%rax,4),%xmm13    ## jzMa  -  -  -
    movss  44(%rdi,%rcx,4),%xmm14    ## jzMc  -  -  -
    movss  44(%rdi,%rbx,4),%xmm1     ## jzMb  -  -  -
    movss  44(%rdi,%rdx,4),%xmm2     ## jzMd  -  -  -
    movlhps %xmm1,%xmm13 ## jzMa  -  jzMb  -
    movlhps %xmm2,%xmm14 ## jzMc  -  jzMd -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jzMa jzMb jzMc jzMd

    ## xmm11: jxMa jyMa jxMb jyMb 
    ## xmm12: jxMc jyMc jxMd jyMd
    ## xmm13: jzMa jzMb jzMc jzMd

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb304_rinvH1M(%rsp),%xmm3
    mulps  nb304_rinvH2M(%rsp),%xmm7
    mulps  nb304_rinvMM(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb304_dxH1M(%rsp),%xmm0
        mulps nb304_dyH1M(%rsp),%xmm1
        mulps nb304_dzH1M(%rsp),%xmm2
        mulps nb304_dxH2M(%rsp),%xmm3
        mulps nb304_dyH2M(%rsp),%xmm4
        mulps nb304_dzH2M(%rsp),%xmm5
        mulps nb304_dxMM(%rsp),%xmm6
        mulps nb304_dyMM(%rsp),%xmm7
        mulps nb304_dzMM(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb304_fixH1(%rsp),%xmm0
    addps nb304_fiyH1(%rsp),%xmm1
    addps nb304_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb304_fixH2(%rsp),%xmm3
    addps nb304_fiyH2(%rsp),%xmm4
    addps nb304_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb304_fixM(%rsp),%xmm6
    addps nb304_fiyM(%rsp),%xmm7
    addps nb304_fizM(%rsp),%xmm8

    movaps %xmm0,nb304_fixH1(%rsp)
    movaps %xmm1,nb304_fiyH1(%rsp)
    movaps %xmm2,nb304_fizH1(%rsp)
    movaps %xmm3,nb304_fixH2(%rsp)
    movaps %xmm4,nb304_fiyH2(%rsp)
    movaps %xmm5,nb304_fizH2(%rsp)
    movaps %xmm6,nb304_fixM(%rsp)
    movaps %xmm7,nb304_fiyM(%rsp)
    movaps %xmm8,nb304_fizM(%rsp)

    ## xmm14 = fMx
    ## xmm15 = fMy
    ## xmm13 = fMz
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fMzc fMzd

    movlps %xmm11,36(%rdi,%rax,4)
    movhps %xmm11,36(%rdi,%rbx,4)
    movlps %xmm12,36(%rdi,%rcx,4)
    movhps %xmm12,36(%rdi,%rdx,4)
    movss  %xmm13,44(%rdi,%rax,4)
    movss  %xmm14,44(%rdi,%rcx,4)
    shufps $1,%xmm13,%xmm13
    shufps $1,%xmm14,%xmm14
    movss  %xmm13,44(%rdi,%rbx,4)
    movss  %xmm14,44(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb304_innerk(%rsp)
        jl    _nb_kernel304_x86_64_sse.nb304_single_check
        jmp   _nb_kernel304_x86_64_sse.nb304_unroll_loop
_nb_kernel304_x86_64_sse.nb304_single_check: 
        addl $4,nb304_innerk(%rsp)
        jnz   _nb_kernel304_x86_64_sse.nb304_single_loop
        jmp   _nb_kernel304_x86_64_sse.nb304_updateouterdata
_nb_kernel304_x86_64_sse.nb304_single_loop: 
        movq  nb304_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb304_innerjjnr(%rsp)

        movq nb304_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm0,%xmm0
        xorps %xmm1,%xmm1
        xorps %xmm2,%xmm2
        movss 36(%rsi,%rax,4),%xmm0             ## jxM  -  -  -
        movss 40(%rsi,%rax,4),%xmm1             ## jyM  -  -  -
        movss 44(%rsi,%rax,4),%xmm2             ## jzM  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm5            ## xmm5 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm5,%xmm7                    ## xmm7 = jzH1 jzH2   -    -

        movlhps %xmm6,%xmm0                     ## xmm0 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm1 ## 11100100     ;# xmm1 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm2 ## 01000100     ;# xmm2 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM 
        movaps %xmm0,nb304_jxM(%rsp)
        movaps %xmm1,nb304_jyM(%rsp)
        movaps %xmm2,nb304_jzM(%rsp)
        subps  nb304_ixM(%rsp),%xmm0
        subps  nb304_iyM(%rsp),%xmm1
        subps  nb304_izM(%rsp),%xmm2
        movaps %xmm0,nb304_dxMM(%rsp)
        movaps %xmm1,nb304_dyMM(%rsp)
        movaps %xmm2,nb304_dzMM(%rsp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1  ## rinv
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb304_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movq nb304_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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
        mulps  nb304_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMM(%rsp),%xmm3
        movhps  nb304_qqMH(%rsp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb304_vctot(%rsp),%xmm5
        movaps %xmm5,nb304_vctot(%rsp)
        xorps  %xmm2,%xmm2
        mulps  nb304_tsc(%rsp),%xmm3

        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb304_dxMM(%rsp),%xmm0
        mulps   nb304_dyMM(%rsp),%xmm1
        mulps   nb304_dzMM(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb304_fjxM(%rsp)
        movaps  %xmm4,nb304_fjyM(%rsp)
        movaps  %xmm5,nb304_fjzM(%rsp)
        addps   nb304_fixM(%rsp),%xmm0
        addps   nb304_fiyM(%rsp),%xmm1
        addps   nb304_fizM(%rsp),%xmm2
        movaps  %xmm0,nb304_fixM(%rsp)
        movaps  %xmm1,nb304_fiyM(%rsp)
        movaps  %xmm2,nb304_fizM(%rsp)


        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb304_jxM(%rsp),%xmm0
    movaps  nb304_jyM(%rsp),%xmm1
    movaps  nb304_jzM(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5

        subps   nb304_ixH1(%rsp),%xmm0
        subps   nb304_iyH1(%rsp),%xmm1
        subps   nb304_izH1(%rsp),%xmm2
        subps   nb304_ixH2(%rsp),%xmm3
        subps   nb304_iyH2(%rsp),%xmm4
        subps   nb304_izH2(%rsp),%xmm5

        movaps %xmm0,nb304_dxH1M(%rsp)
        movaps %xmm1,nb304_dyH1M(%rsp)
        movaps %xmm2,nb304_dzH1M(%rsp)
        movaps %xmm3,nb304_dxH2M(%rsp)
        movaps %xmm4,nb304_dyH2M(%rsp)
        movaps %xmm5,nb304_dzH2M(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## start with H1, save H2 data 
        movaps %xmm4,nb304_rsqH2M(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb304_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb304_rinvH2M(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb304_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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
        mulps  nb304_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMH(%rsp),%xmm3
        movhps  nb304_qqHH(%rsp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb304_vctot(%rsp),%xmm5
        movaps %xmm5,nb304_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb304_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb304_dxH1M(%rsp),%xmm0
        mulps   nb304_dyH1M(%rsp),%xmm1
        mulps   nb304_dzH1M(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb304_fjxM(%rsp),%xmm3
        movaps  nb304_fjyM(%rsp),%xmm4
        movaps  nb304_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb304_fjxM(%rsp)
        movaps  %xmm4,nb304_fjyM(%rsp)
        movaps  %xmm5,nb304_fjzM(%rsp)
        addps   nb304_fixH1(%rsp),%xmm0
        addps   nb304_fiyH1(%rsp),%xmm1
        addps   nb304_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb304_fixH1(%rsp)
        movaps  %xmm1,nb304_fiyH1(%rsp)
        movaps  %xmm2,nb304_fizH1(%rsp)
        ## do table for H2 - j water interaction 
        movaps nb304_rinvH2M(%rsp),%xmm0
        movaps nb304_rsqH2M(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb304_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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
        mulps  nb304_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304_qqMH(%rsp),%xmm3
        movhps  nb304_qqHH(%rsp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb304_vctot(%rsp),%xmm5
        movaps %xmm5,nb304_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb304_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb304_dxH2M(%rsp),%xmm0
        mulps   nb304_dyH2M(%rsp),%xmm1
        mulps   nb304_dzH2M(%rsp),%xmm2
        movaps  nb304_fjxM(%rsp),%xmm3
        movaps  nb304_fjyM(%rsp),%xmm4
        movaps  nb304_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb304_faction(%rbp),%rsi
        movaps  %xmm3,nb304_fjxM(%rsp)
        movaps  %xmm4,nb304_fjyM(%rsp)
        movaps  %xmm5,nb304_fjzM(%rsp)
        addps   nb304_fixH2(%rsp),%xmm0
        addps   nb304_fiyH2(%rsp),%xmm1
        addps   nb304_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb304_fixH2(%rsp)
        movaps  %xmm1,nb304_fiyH2(%rsp)
        movaps  %xmm2,nb304_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  36(%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb304_fjxM(%rsp),%xmm3
        movaps  nb304_fjyM(%rsp),%xmm4
        movaps  nb304_fjzM(%rsp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## 00000010
        shufps $3,%xmm7,%xmm7 ## 00000011
        addss   44(%rsi,%rax,4),%xmm5
        addss   20(%rsi,%rax,4),%xmm6
        addss   32(%rsi,%rax,4),%xmm7
        movss   %xmm5,44(%rsi,%rax,4)
        movss   %xmm6,20(%rsi,%rax,4)
        movss   %xmm7,32(%rsi,%rax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,36(%rsi,%rax,4)
        movlps  %xmm1,12(%rsi,%rax,4)
        movhps  %xmm1,24(%rsi,%rax,4)

        decl nb304_innerk(%rsp)
        jz    _nb_kernel304_x86_64_sse.nb304_updateouterdata
        jmp   _nb_kernel304_x86_64_sse.nb304_single_loop
_nb_kernel304_x86_64_sse.nb304_updateouterdata: 
        movl  nb304_ii3(%rsp),%ecx
        movq  nb304_faction(%rbp),%rdi
        movq  nb304_fshift(%rbp),%rsi
        movl  nb304_is3(%rsp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb304_fixH1(%rsp),%xmm0
        movaps nb304_fiyH1(%rsp),%xmm1
        movaps nb304_fizH1(%rsp),%xmm2

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

        ## accumulate H2 i forces in xmm0, xmm1, xmm2 
        movaps nb304_fixH2(%rsp),%xmm0
        movaps nb304_fiyH2(%rsp),%xmm1
        movaps nb304_fizH2(%rsp),%xmm2

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
        movaps nb304_fixM(%rsp),%xmm0
        movaps nb304_fiyM(%rsp),%xmm1
        movaps nb304_fizM(%rsp),%xmm2

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
        movl nb304_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb304_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb304_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb304_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb304_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304_x86_64_sse.nb304_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304_n(%rsp)
        jmp _nb_kernel304_x86_64_sse.nb304_outer
_nb_kernel304_x86_64_sse.nb304_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304_x86_64_sse.nb304_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304_x86_64_sse.nb304_threadloop
_nb_kernel304_x86_64_sse.nb304_end: 
        movl nb304_nouter(%rsp),%eax
        movl nb304_ninner(%rsp),%ebx
        movq nb304_outeriter(%rbp),%rcx
        movq nb304_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1528,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel304nf_x86_64_sse
.globl _nb_kernel304nf_x86_64_sse
nb_kernel304nf_x86_64_sse:      
_nb_kernel304nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb304nf_fshift, 16
.set nb304nf_gid, 24
.set nb304nf_pos, 32
.set nb304nf_faction, 40
.set nb304nf_charge, 48
.set nb304nf_p_facel, 56
.set nb304nf_argkrf, 64
.set nb304nf_argcrf, 72
.set nb304nf_Vc, 80
.set nb304nf_type, 88
.set nb304nf_p_ntype, 96
.set nb304nf_vdwparam, 104
.set nb304nf_Vvdw, 112
.set nb304nf_p_tabscale, 120
.set nb304nf_VFtab, 128
.set nb304nf_invsqrta, 136
.set nb304nf_dvda, 144
.set nb304nf_p_gbtabscale, 152
.set nb304nf_GBtab, 160
.set nb304nf_p_nthreads, 168
.set nb304nf_count, 176
.set nb304nf_mtx, 184
.set nb304nf_outeriter, 192
.set nb304nf_inneriter, 200
.set nb304nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb304nf_ixH1, 0
.set nb304nf_iyH1, 16
.set nb304nf_izH1, 32
.set nb304nf_ixH2, 48
.set nb304nf_iyH2, 64
.set nb304nf_izH2, 80
.set nb304nf_ixM, 96
.set nb304nf_iyM, 112
.set nb304nf_izM, 128
.set nb304nf_jxH1, 144
.set nb304nf_jyH1, 160
.set nb304nf_jzH1, 176
.set nb304nf_jxH2, 192
.set nb304nf_jyH2, 208
.set nb304nf_jzH2, 224
.set nb304nf_jxM, 240
.set nb304nf_jyM, 256
.set nb304nf_jzM, 272
.set nb304nf_qqHH, 288
.set nb304nf_qqMH, 304
.set nb304nf_qqMM, 320
.set nb304nf_tsc, 336
.set nb304nf_vctot, 352
.set nb304nf_half, 368
.set nb304nf_three, 384
.set nb304nf_rsqH1H1, 400
.set nb304nf_rsqH1H2, 416
.set nb304nf_rsqH1M, 432
.set nb304nf_rsqH2H1, 448
.set nb304nf_rsqH2H2, 464
.set nb304nf_rsqH2M, 480
.set nb304nf_rsqMH1, 496
.set nb304nf_rsqMH2, 512
.set nb304nf_rsqMM, 528
.set nb304nf_rinvH1H1, 544
.set nb304nf_rinvH1H2, 560
.set nb304nf_rinvH1M, 576
.set nb304nf_rinvH2H1, 592
.set nb304nf_rinvH2H2, 608
.set nb304nf_rinvH2M, 624
.set nb304nf_rinvMH1, 640
.set nb304nf_rinvMH2, 656
.set nb304nf_rinvMM, 672
.set nb304nf_is3, 688
.set nb304nf_ii3, 692
.set nb304nf_nri, 696
.set nb304nf_iinr, 704
.set nb304nf_jindex, 712
.set nb304nf_jjnr, 720
.set nb304nf_shift, 728
.set nb304nf_shiftvec, 736
.set nb304nf_facel, 744
.set nb304nf_innerjjnr, 752
.set nb304nf_innerk, 760
.set nb304nf_n, 764
.set nb304nf_nn1, 768
.set nb304nf_nouter, 772
.set nb304nf_ninner, 776
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $792,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb304nf_nouter(%rsp)
        movl %eax,nb304nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb304nf_nri(%rsp)
        movq %rsi,nb304nf_iinr(%rsp)
        movq %rdx,nb304nf_jindex(%rsp)
        movq %rcx,nb304nf_jjnr(%rsp)
        movq %r8,nb304nf_shift(%rsp)
        movq %r9,nb304nf_shiftvec(%rsp)
        movq nb304nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb304nf_facel(%rsp)

        movq nb304nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb304nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb304nf_half(%rsp)
        movss nb304nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb304nf_half(%rsp)
        movaps %xmm3,nb304nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb304nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb304nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb304nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb304nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304nf_qqHH(%rsp)
        movaps %xmm4,nb304nf_qqMH(%rsp)
        movaps %xmm5,nb304nf_qqMM(%rsp)

_nb_kernel304nf_x86_64_sse.nb304nf_threadloop: 
        movq  nb304nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel304nf_x86_64_sse.nb304nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304nf_x86_64_sse.nb304nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304nf_n(%rsp)
        movl %ebx,nb304nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304nf_x86_64_sse.nb304nf_outerstart
        jmp _nb_kernel304nf_x86_64_sse.nb304nf_end

_nb_kernel304nf_x86_64_sse.nb304nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304nf_nouter(%rsp),%ebx
        movl %ebx,nb304nf_nouter(%rsp)

_nb_kernel304nf_x86_64_sse.nb304nf_outer: 
        movq  nb304nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb304nf_is3(%rsp)            ## store is3 

        movq  nb304nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb304nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb304nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb304nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb304nf_ixH1(%rsp)
        movaps %xmm4,nb304nf_iyH1(%rsp)
        movaps %xmm5,nb304nf_izH1(%rsp)

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
        movaps %xmm0,nb304nf_ixH2(%rsp)
        movaps %xmm1,nb304nf_iyH2(%rsp)
        movaps %xmm2,nb304nf_izH2(%rsp)
        movaps %xmm3,nb304nf_ixM(%rsp)
        movaps %xmm4,nb304nf_iyM(%rsp)
        movaps %xmm5,nb304nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb304nf_vctot(%rsp)

        movq  nb304nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb304nf_pos(%rbp),%rsi
        movq  nb304nf_faction(%rbp),%rdi
        movq  nb304nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb304nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb304nf_ninner(%rsp),%ecx
        movl  %ecx,nb304nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb304nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel304nf_x86_64_sse.nb304nf_unroll_loop
        jmp   _nb_kernel304nf_x86_64_sse.nb304nf_single_check
_nb_kernel304nf_x86_64_sse.nb304nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb304nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb304nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb304nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j coordinates to local temp variables 
        movlps 12(%rsi,%rax,4),%xmm2
        movlps 24(%rsi,%rax,4),%xmm3
        movlps 36(%rsi,%rax,4),%xmm4

        movlps 12(%rsi,%rbx,4),%xmm5
        movlps 24(%rsi,%rbx,4),%xmm6
        movlps 36(%rsi,%rbx,4),%xmm7

        movhps 12(%rsi,%rcx,4),%xmm2
        movhps 24(%rsi,%rcx,4),%xmm3
        movhps 36(%rsi,%rcx,4),%xmm4

        movhps 12(%rsi,%rdx,4),%xmm5
        movhps 24(%rsi,%rdx,4),%xmm6
        movhps 36(%rsi,%rdx,4),%xmm7

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb304nf_jxH1(%rsp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb304nf_jyH1(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb304nf_jxH2(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb304nf_jyH2(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb304nf_jxM(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb304nf_jyM(%rsp)

        movss  20(%rsi,%rax,4),%xmm0
        movss  32(%rsi,%rax,4),%xmm1
        movss  44(%rsi,%rax,4),%xmm2

        movss  20(%rsi,%rcx,4),%xmm3
        movss  32(%rsi,%rcx,4),%xmm4
        movss  44(%rsi,%rcx,4),%xmm5

        movhps 16(%rsi,%rbx,4),%xmm0
        movhps 28(%rsi,%rbx,4),%xmm1
        movhps 40(%rsi,%rbx,4),%xmm2

        movhps 16(%rsi,%rdx,4),%xmm3
        movhps 28(%rsi,%rdx,4),%xmm4
        movhps 40(%rsi,%rdx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## 11001100
        shufps $204,%xmm4,%xmm1 ## 11001100
        shufps $204,%xmm5,%xmm2 ## 11001100
        movaps %xmm0,nb304nf_jzH1(%rsp)
        movaps %xmm1,nb304nf_jzH2(%rsp)
        movaps %xmm2,nb304nf_jzM(%rsp)

        movaps nb304nf_ixH1(%rsp),%xmm0
        movaps nb304nf_iyH1(%rsp),%xmm1
        movaps nb304nf_izH1(%rsp),%xmm2
        movaps nb304nf_ixH1(%rsp),%xmm3
        movaps nb304nf_iyH1(%rsp),%xmm4
        movaps nb304nf_izH1(%rsp),%xmm5
        subps  nb304nf_jxH1(%rsp),%xmm0
        subps  nb304nf_jyH1(%rsp),%xmm1
        subps  nb304nf_jzH1(%rsp),%xmm2
        subps  nb304nf_jxH2(%rsp),%xmm3
        subps  nb304nf_jyH2(%rsp),%xmm4
        subps  nb304nf_jzH2(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304nf_rsqH1H1(%rsp)
        movaps %xmm3,nb304nf_rsqH1H2(%rsp)

        movaps nb304nf_ixH1(%rsp),%xmm0
        movaps nb304nf_iyH1(%rsp),%xmm1
        movaps nb304nf_izH1(%rsp),%xmm2
        movaps nb304nf_ixH2(%rsp),%xmm3
        movaps nb304nf_iyH2(%rsp),%xmm4
        movaps nb304nf_izH2(%rsp),%xmm5
        subps  nb304nf_jxM(%rsp),%xmm0
        subps  nb304nf_jyM(%rsp),%xmm1
        subps  nb304nf_jzM(%rsp),%xmm2
        subps  nb304nf_jxH1(%rsp),%xmm3
        subps  nb304nf_jyH1(%rsp),%xmm4
        subps  nb304nf_jzH1(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304nf_rsqH1M(%rsp)
        movaps %xmm3,nb304nf_rsqH2H1(%rsp)

        movaps nb304nf_ixH2(%rsp),%xmm0
        movaps nb304nf_iyH2(%rsp),%xmm1
        movaps nb304nf_izH2(%rsp),%xmm2
        movaps nb304nf_ixH2(%rsp),%xmm3
        movaps nb304nf_iyH2(%rsp),%xmm4
        movaps nb304nf_izH2(%rsp),%xmm5
        subps  nb304nf_jxH2(%rsp),%xmm0
        subps  nb304nf_jyH2(%rsp),%xmm1
        subps  nb304nf_jzH2(%rsp),%xmm2
        subps  nb304nf_jxM(%rsp),%xmm3
        subps  nb304nf_jyM(%rsp),%xmm4
        subps  nb304nf_jzM(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb304nf_rsqH2H2(%rsp)
        movaps %xmm3,nb304nf_rsqH2M(%rsp)

        movaps nb304nf_ixM(%rsp),%xmm0
        movaps nb304nf_iyM(%rsp),%xmm1
        movaps nb304nf_izM(%rsp),%xmm2
        movaps nb304nf_ixM(%rsp),%xmm3
        movaps nb304nf_iyM(%rsp),%xmm4
        movaps nb304nf_izM(%rsp),%xmm5
        subps  nb304nf_jxH1(%rsp),%xmm0
        subps  nb304nf_jyH1(%rsp),%xmm1
        subps  nb304nf_jzH1(%rsp),%xmm2
        subps  nb304nf_jxH2(%rsp),%xmm3
        subps  nb304nf_jyH2(%rsp),%xmm4
        subps  nb304nf_jzH2(%rsp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb304nf_rsqMH1(%rsp)
        movaps %xmm4,nb304nf_rsqMH2(%rsp)

        movaps nb304nf_ixM(%rsp),%xmm0
        movaps nb304nf_iyM(%rsp),%xmm1
        movaps nb304nf_izM(%rsp),%xmm2
        subps  nb304nf_jxM(%rsp),%xmm0
        subps  nb304nf_jyM(%rsp),%xmm1
        subps  nb304nf_jzM(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb304nf_rsqMM(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%rsp),%xmm3   ## rinvMM
        mulps   nb304nf_half(%rsp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb304nf_rinvMM(%rsp)
        movaps  %xmm7,nb304nf_rinvMH2(%rsp)

        rsqrtps nb304nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb304nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH1H1(%rsp),%xmm1
        mulps   nb304nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%rsp),%xmm3
        mulps   nb304nf_half(%rsp),%xmm7
        movaps  %xmm3,nb304nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb304nf_rinvH1H2(%rsp)

        rsqrtps nb304nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb304nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH1M(%rsp),%xmm1
        mulps   nb304nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%rsp),%xmm3
        mulps   nb304nf_half(%rsp),%xmm7
        movaps  %xmm3,nb304nf_rinvH1M(%rsp)
        movaps  %xmm7,nb304nf_rinvH2H1(%rsp)

        rsqrtps nb304nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb304nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb304nf_rsqH2H2(%rsp),%xmm1
        mulps   nb304nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%rsp),%xmm3
        mulps   nb304nf_half(%rsp),%xmm7
        movaps  %xmm3,nb304nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb304nf_rinvH2M(%rsp)

        rsqrtps nb304nf_rsqMH1(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304nf_three(%rsp),%xmm3
        mulps   nb304nf_rsqMH1(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304nf_half(%rsp),%xmm3
        movaps  %xmm3,nb304nf_rinvMH1(%rsp)

        ## start with H1-H1 interaction 
        movaps nb304nf_rinvH1H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb304nf_VFtab(%rbp),%rsi
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
        movaps nb304nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## H1-H2 interaction 
        movaps nb304nf_rinvH1H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## H1-M interaction  
        movaps nb304nf_rinvH1M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## H2-H1 interaction 
        movaps nb304nf_rinvH2H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## H2-H2 interaction 
        movaps nb304nf_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## H2-M interaction 
        movaps nb304nf_rinvH2M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## M-H1 interaction 
        movaps nb304nf_rinvMH1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## M-H2 interaction 
        movaps nb304nf_rinvMH2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## M-M interaction 
        movaps nb304nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb304nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld   $2,%mm6
        pslld   $2,%mm7

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
        movaps nb304nf_qqMM(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb304nf_innerk(%rsp)
        jl    _nb_kernel304nf_x86_64_sse.nb304nf_single_check
        jmp   _nb_kernel304nf_x86_64_sse.nb304nf_unroll_loop
_nb_kernel304nf_x86_64_sse.nb304nf_single_check: 
        addl $4,nb304nf_innerk(%rsp)
        jnz   _nb_kernel304nf_x86_64_sse.nb304nf_single_loop
        jmp   _nb_kernel304nf_x86_64_sse.nb304nf_updateouterdata
_nb_kernel304nf_x86_64_sse.nb304nf_single_loop: 
        movq  nb304nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb304nf_innerjjnr(%rsp)

        movq nb304nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%rsi,%rax,4),%xmm3             ## jxM  -  -  -
        movss 40(%rsi,%rax,4),%xmm4             ## jyM  -  -  -
        movss 44(%rsi,%rax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb304nf_ixM(%rsp),%xmm0
        movaps  nb304nf_iyM(%rsp),%xmm1
        movaps  nb304nf_izM(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM 
        movaps %xmm3,nb304nf_jxM(%rsp)
        movaps %xmm4,nb304nf_jyM(%rsp)
        movaps %xmm5,nb304nf_jzM(%rsp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb304nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb304nf_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb304nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movq nb304nf_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMM(%rsp),%xmm3
        movhps  nb304nf_qqMH(%rsp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb304nf_ixH1(%rsp),%xmm0
        movaps  nb304nf_iyH1(%rsp),%xmm1
        movaps  nb304nf_izH1(%rsp),%xmm2
        movaps  nb304nf_ixH2(%rsp),%xmm3
        movaps  nb304nf_iyH2(%rsp),%xmm4
        movaps  nb304nf_izH2(%rsp),%xmm5
        subps   nb304nf_jxM(%rsp),%xmm0
        subps   nb304nf_jyM(%rsp),%xmm1
        subps   nb304nf_jzM(%rsp),%xmm2
        subps   nb304nf_jxM(%rsp),%xmm3
        subps   nb304nf_jyM(%rsp),%xmm4
        subps   nb304nf_jzM(%rsp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## start with H1, save H2 data 
        movaps %xmm4,nb304nf_rsqH2M(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb304nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb304nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb304nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb304nf_rinvH2M(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb304nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMH(%rsp),%xmm3
        movhps  nb304nf_qqHH(%rsp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        ## do table for H2 - j water interaction 
        movaps nb304nf_rinvH2M(%rsp),%xmm0
        movaps nb304nf_rsqH2M(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb304nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
        cvttps2pi %xmm1,%mm6
        cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6
        pslld   $2,%mm7

        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        movlps (%rsi,%rbx,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rbx,4),%xmm7
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

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb304nf_qqMH(%rsp),%xmm3
        movhps  nb304nf_qqHH(%rsp),%xmm3

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb304nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb304nf_vctot(%rsp)

        decl nb304nf_innerk(%rsp)
        jz    _nb_kernel304nf_x86_64_sse.nb304nf_updateouterdata
        jmp   _nb_kernel304nf_x86_64_sse.nb304nf_single_loop
_nb_kernel304nf_x86_64_sse.nb304nf_updateouterdata: 
        ## get n from stack
        movl nb304nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb304nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb304nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb304nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb304nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304nf_x86_64_sse.nb304nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304nf_n(%rsp)
        jmp _nb_kernel304nf_x86_64_sse.nb304nf_outer
_nb_kernel304nf_x86_64_sse.nb304nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304nf_x86_64_sse.nb304nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304nf_x86_64_sse.nb304nf_threadloop
_nb_kernel304nf_x86_64_sse.nb304nf_end: 
        movl nb304nf_nouter(%rsp),%eax
        movl nb304nf_ninner(%rsp),%ebx
        movq nb304nf_outeriter(%rbp),%rcx
        movq nb304nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $792,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





