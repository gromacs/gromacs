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





.globl nb_kernel332_x86_64_sse
.globl _nb_kernel332_x86_64_sse
nb_kernel332_x86_64_sse:        
_nb_kernel332_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb332_fshift, 16
.set nb332_gid, 24
.set nb332_pos, 32
.set nb332_faction, 40
.set nb332_charge, 48
.set nb332_p_facel, 56
.set nb332_argkrf, 64
.set nb332_argcrf, 72
.set nb332_Vc, 80
.set nb332_type, 88
.set nb332_p_ntype, 96
.set nb332_vdwparam, 104
.set nb332_Vvdw, 112
.set nb332_p_tabscale, 120
.set nb332_VFtab, 128
.set nb332_invsqrta, 136
.set nb332_dvda, 144
.set nb332_p_gbtabscale, 152
.set nb332_GBtab, 160
.set nb332_p_nthreads, 168
.set nb332_count, 176
.set nb332_mtx, 184
.set nb332_outeriter, 192
.set nb332_inneriter, 200
.set nb332_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb332_ixO, 0
.set nb332_iyO, 16
.set nb332_izO, 32
.set nb332_ixH1, 48
.set nb332_iyH1, 64
.set nb332_izH1, 80
.set nb332_ixH2, 96
.set nb332_iyH2, 112
.set nb332_izH2, 128
.set nb332_jxO, 144
.set nb332_jyO, 160
.set nb332_jzO, 176
.set nb332_jxH1, 192
.set nb332_jyH1, 208
.set nb332_jzH1, 224
.set nb332_jxH2, 240
.set nb332_jyH2, 256
.set nb332_jzH2, 272
.set nb332_dxOO, 288
.set nb332_dyOO, 304
.set nb332_dzOO, 320
.set nb332_dxOH1, 336
.set nb332_dyOH1, 352
.set nb332_dzOH1, 368
.set nb332_dxOH2, 384
.set nb332_dyOH2, 400
.set nb332_dzOH2, 416
.set nb332_dxH1O, 432
.set nb332_dyH1O, 448
.set nb332_dzH1O, 464
.set nb332_dxH1H1, 480
.set nb332_dyH1H1, 496
.set nb332_dzH1H1, 512
.set nb332_dxH1H2, 528
.set nb332_dyH1H2, 544
.set nb332_dzH1H2, 560
.set nb332_dxH2O, 576
.set nb332_dyH2O, 592
.set nb332_dzH2O, 608
.set nb332_dxH2H1, 624
.set nb332_dyH2H1, 640
.set nb332_dzH2H1, 656
.set nb332_dxH2H2, 672
.set nb332_dyH2H2, 688
.set nb332_dzH2H2, 704
.set nb332_qqOO, 720
.set nb332_qqOH, 736
.set nb332_qqHH, 752
.set nb332_two, 768
.set nb332_tsc, 784
.set nb332_c6, 800
.set nb332_c12, 816
.set nb332_vctot, 832
.set nb332_Vvdwtot, 848
.set nb332_fixO, 864
.set nb332_fiyO, 880
.set nb332_fizO, 896
.set nb332_fixH1, 912
.set nb332_fiyH1, 928
.set nb332_fizH1, 944
.set nb332_fixH2, 960
.set nb332_fiyH2, 976
.set nb332_fizH2, 992
.set nb332_fjxO, 1008
.set nb332_fjyO, 1024
.set nb332_fjzO, 1040
.set nb332_fjxH1, 1056
.set nb332_fjyH1, 1072
.set nb332_fjzH1, 1088
.set nb332_fjxH2, 1104
.set nb332_fjyH2, 1120
.set nb332_fjzH2, 1136
.set nb332_half, 1152
.set nb332_three, 1168
.set nb332_epsO, 1184
.set nb332_epsH1, 1200
.set nb332_epsH2, 1216
.set nb332_rsqH1O, 1232
.set nb332_rsqOO, 1248
.set nb332_rsqH1H2, 1264
.set nb332_rsqH2O, 1280
.set nb332_fstmpH1, 1296
.set nb332_fstmpH2, 1312
.set nb332_rinvOO, 1328
.set nb332_rinvOH1, 1344
.set nb332_rinvOH2, 1360
.set nb332_rinvH1O, 1376
.set nb332_rinvH1H1, 1392
.set nb332_rinvH1H2, 1408
.set nb332_rinvH2O, 1424
.set nb332_rinvH2H1, 1440
.set nb332_rinvH2H2, 1456
.set nb332_fstmp, 1472
.set nb332_is3, 1488
.set nb332_ii3, 1492
.set nb332_nri, 1496
.set nb332_iinr, 1504
.set nb332_jindex, 1512
.set nb332_jjnr, 1520
.set nb332_shift, 1528
.set nb332_shiftvec, 1536
.set nb332_facel, 1544
.set nb332_innerjjnr, 1552
.set nb332_innerk, 1560
.set nb332_n, 1564
.set nb332_nn1, 1568
.set nb332_nouter, 1572
.set nb332_ninner, 1576
        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
        subq $1592,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb332_nouter(%rsp)
        movl %eax,nb332_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb332_nri(%rsp)
        movq %rsi,nb332_iinr(%rsp)
        movq %rdx,nb332_jindex(%rsp)
        movq %rcx,nb332_jjnr(%rsp)
        movq %r8,nb332_shift(%rsp)
        movq %r9,nb332_shiftvec(%rsp)
        movq nb332_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb332_facel(%rsp)

        movq nb332_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb332_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb332_half(%rsp)
        movss nb332_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb332_half(%rsp)
        movaps %xmm2,nb332_two(%rsp)
        movaps %xmm3,nb332_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb332_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb332_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb332_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb332_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb332_qqOO(%rsp)
        movaps %xmm4,nb332_qqOH(%rsp)
        movaps %xmm5,nb332_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb332_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb332_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb332_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb332_c6(%rsp)
        movaps %xmm1,nb332_c12(%rsp)

_nb_kernel332_x86_64_sse.nb332_threadloop: 
        movq  nb332_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel332_x86_64_sse.nb332_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel332_x86_64_sse.nb332_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb332_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb332_n(%rsp)
        movl %ebx,nb332_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel332_x86_64_sse.nb332_outerstart
        jmp _nb_kernel332_x86_64_sse.nb332_end

_nb_kernel332_x86_64_sse.nb332_outerstart: 
        ## ebx contains number of outer iterations
        addl nb332_nouter(%rsp),%ebx
        movl %ebx,nb332_nouter(%rsp)

_nb_kernel332_x86_64_sse.nb332_outer: 
        movq  nb332_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb332_is3(%rsp)      ## store is3 

        movq  nb332_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb332_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb332_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb332_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb332_ixO(%rsp)
        movaps %xmm4,nb332_iyO(%rsp)
        movaps %xmm5,nb332_izO(%rsp)

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
        movaps %xmm0,nb332_ixH1(%rsp)
        movaps %xmm1,nb332_iyH1(%rsp)
        movaps %xmm2,nb332_izH1(%rsp)
        movaps %xmm3,nb332_ixH2(%rsp)
        movaps %xmm4,nb332_iyH2(%rsp)
        movaps %xmm5,nb332_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb332_vctot(%rsp)
        movaps %xmm4,nb332_Vvdwtot(%rsp)
        movaps %xmm4,nb332_fixO(%rsp)
        movaps %xmm4,nb332_fiyO(%rsp)
        movaps %xmm4,nb332_fizO(%rsp)
        movaps %xmm4,nb332_fixH1(%rsp)
        movaps %xmm4,nb332_fiyH1(%rsp)
        movaps %xmm4,nb332_fizH1(%rsp)
        movaps %xmm4,nb332_fixH2(%rsp)
        movaps %xmm4,nb332_fiyH2(%rsp)
        movaps %xmm4,nb332_fizH2(%rsp)

        movq  nb332_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb332_pos(%rbp),%rsi
        movq  nb332_faction(%rbp),%rdi
        movq  nb332_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb332_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb332_ninner(%rsp),%ecx
        movl  %ecx,nb332_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb332_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel332_x86_64_sse.nb332_unroll_loop
        jmp   _nb_kernel332_x86_64_sse.nb332_single_check
_nb_kernel332_x86_64_sse.nb332_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb332_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb332_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb332_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j O coordinates to local temp variables 
    movlps (%rsi,%rax,4),%xmm0 ## jxOa jyOa  -   -
    movlps (%rsi,%rcx,4),%xmm1 ## jxOc jyOc  -   -
    movhps (%rsi,%rbx,4),%xmm0 ## jxOa jyOa jxOb jyOb 
    movhps (%rsi,%rdx,4),%xmm1 ## jxOc jyOc jxOd jyOd 

    movss  8(%rsi,%rax,4),%xmm2    ## jzOa  -  -  -
    movss  8(%rsi,%rcx,4),%xmm3    ## jzOc  -  -  -
    movhps 8(%rsi,%rbx,4),%xmm2    ## jzOa  -  jzOb  -
    movhps 8(%rsi,%rdx,4),%xmm3    ## jzOc  -  jzOd -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxOa jxOc jyOa jyOc        
    unpckhps %xmm1,%xmm4 ## jxOb jxOd jyOb jyOd
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb332_ixO(%rsp),%xmm0
    subps nb332_iyO(%rsp),%xmm1
    subps nb332_izO(%rsp),%xmm2
    subps nb332_ixH1(%rsp),%xmm3
    subps nb332_iyH1(%rsp),%xmm4
    subps nb332_izH1(%rsp),%xmm5
    subps nb332_ixH2(%rsp),%xmm6
    subps nb332_iyH2(%rsp),%xmm7
    subps nb332_izH2(%rsp),%xmm8

    movd %eax,%mm0 ## save j3 in mm0-mm3
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

        movaps %xmm0,nb332_dxOO(%rsp)
        movaps %xmm1,nb332_dyOO(%rsp)
        movaps %xmm2,nb332_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb332_dxH1O(%rsp)
        movaps %xmm4,nb332_dyH1O(%rsp)
        movaps %xmm5,nb332_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb332_dxH2O(%rsp)
        movaps %xmm7,nb332_dyH2O(%rsp)
        movaps %xmm8,nb332_dzH2O(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb332_three(%rsp),%xmm9
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

        movaps  nb332_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOO 
        mulps   %xmm4,%xmm10 ## rinvH1O
    mulps   %xmm4,%xmm11 ## rinvH2O

        movaps  %xmm9,nb332_rinvOO(%rsp)
        movaps  %xmm10,nb332_rinvH1O(%rsp)
        movaps  %xmm11,nb332_rinvH2O(%rsp)

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movaps %xmm9,nb332_rinvOO(%rsp)
    movaps nb332_tsc(%rsp),%xmm1

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

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm1,%xmm10
    movaps  %xmm4,%xmm11
    movaps  %xmm7,%xmm12
    pslld   $1,%xmm1
    pslld   $1,%xmm4
    pslld   $1,%xmm7
    paddd   %xmm10,%xmm1
    paddd   %xmm11,%xmm4
    paddd   %xmm12,%xmm7

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

    movq nb332_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb332_epsO(%rsp)
    movaps    %xmm3,nb332_epsH1(%rsp)
    movaps    %xmm6,nb332_epsH2(%rsp)

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

    mulps  nb332_epsO(%rsp),%xmm3     ## Heps
    mulps  nb332_epsH1(%rsp),%xmm7
    mulps  nb332_epsH2(%rsp),%xmm11
    mulps  nb332_epsO(%rsp),%xmm2     ## Geps
    mulps  nb332_epsH1(%rsp),%xmm6
    mulps  nb332_epsH2(%rsp),%xmm10
    mulps  nb332_epsO(%rsp),%xmm3     ## Heps2
    mulps  nb332_epsH1(%rsp),%xmm7
    mulps  nb332_epsH2(%rsp),%xmm11

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
    mulps  nb332_epsO(%rsp),%xmm1     ## eps*Fp
    mulps  nb332_epsH1(%rsp),%xmm5
    mulps  nb332_epsH2(%rsp),%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb332_qqOO(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb332_qqOH(%rsp),%xmm5
    mulps  nb332_qqOH(%rsp),%xmm9
    mulps  nb332_qqOO(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb332_qqOH(%rsp),%xmm7
    mulps  nb332_qqOH(%rsp),%xmm11

    ## accumulate vctot
    addps  nb332_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb332_vctot(%rsp)

    movaps %xmm7,%xmm2
    movaps %xmm11,%xmm1

    ## fij coul in xmm3, xmm2, xmm1    

    ## calculate LJ table
    movlps 16(%rsi,%rax,4),%xmm5
        movlps 32(%rsi,%rax,4),%xmm9

        movlps 16(%rsi,%rcx,4),%xmm7
        movlps 32(%rsi,%rcx,4),%xmm11

        movhps 16(%rsi,%rbx,4),%xmm5
        movhps 32(%rsi,%rbx,4),%xmm9

        movhps 16(%rsi,%rdx,4),%xmm7
        movhps 32(%rsi,%rdx,4),%xmm11

    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movlps 40(%rsi,%rax,4),%xmm11

        movlps 24(%rsi,%rcx,4),%xmm13
        movlps 40(%rsi,%rcx,4),%xmm14

        movhps 24(%rsi,%rbx,4),%xmm7
        movhps 40(%rsi,%rbx,4),%xmm11

        movhps 24(%rsi,%rdx,4),%xmm13
        movhps 40(%rsi,%rdx,4),%xmm14

    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm13,%xmm6 ## 10001000
        shufps $136,%xmm14,%xmm10 ## 10001000
        shufps $221,%xmm13,%xmm7 ## 11011101
        shufps $221,%xmm14,%xmm11 ## 11011101
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    movaps nb332_epsO(%rsp),%xmm0

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
    movaps nb332_c6(%rsp),%xmm12
    movaps nb332_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb332_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7

    addps  %xmm7,%xmm3
    movaps nb332_tsc(%rsp),%xmm10

    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm2
    mulps  %xmm10,%xmm1

        ## move j O forces to local temp variables 
    movq nb332_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm11 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm12 ## jxOc jyOc  -   -

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    movhps (%rdi,%rbx,4),%xmm11 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm12 ## jxOc jyOc jxOd jyOd 

    subps  %xmm3,%xmm0
    subps  %xmm2,%xmm4
    subps  %xmm1,%xmm8

    movss  8(%rdi,%rax,4),%xmm13    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm14    ## jzOc  -  -  -
    mulps  nb332_rinvOO(%rsp),%xmm0
    mulps  nb332_rinvH1O(%rsp),%xmm4
    mulps  nb332_rinvH2O(%rsp),%xmm8
    movhps 8(%rdi,%rbx,4),%xmm13    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm14    ## jzOc  -  jzOd -    
    shufps $136,%xmm14,%xmm13 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm11: jxOa jyOa jxOb jyOb 
    ## xmm12: jxOc jyOc jxOd jyOd
    ## xmm13: jzOa jzOb jzOc jzOd

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb332_dxOO(%rsp),%xmm0
        mulps nb332_dyOO(%rsp),%xmm1
        mulps nb332_dzOO(%rsp),%xmm2
        mulps nb332_dxH1O(%rsp),%xmm3
        mulps nb332_dyH1O(%rsp),%xmm4
        mulps nb332_dzH1O(%rsp),%xmm5
        mulps nb332_dxH2O(%rsp),%xmm6
        mulps nb332_dyH2O(%rsp),%xmm7
        mulps nb332_dzH2O(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb332_fixO(%rsp),%xmm0
    addps nb332_fiyO(%rsp),%xmm1
    addps nb332_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb332_fixH1(%rsp),%xmm3
    addps nb332_fiyH1(%rsp),%xmm4
    addps nb332_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb332_fixH2(%rsp),%xmm6
    addps nb332_fiyH2(%rsp),%xmm7
    addps nb332_fizH2(%rsp),%xmm8

    movaps %xmm0,nb332_fixO(%rsp)
    movaps %xmm1,nb332_fiyO(%rsp)
    movaps %xmm2,nb332_fizO(%rsp)
    movaps %xmm3,nb332_fixH1(%rsp)
    movaps %xmm4,nb332_fiyH1(%rsp)
    movaps %xmm5,nb332_fizH1(%rsp)
    movaps %xmm6,nb332_fixH2(%rsp)
    movaps %xmm7,nb332_fiyH2(%rsp)
    movaps %xmm8,nb332_fizH2(%rsp)

    ## xmm11 = fOx
    ## xmm12 = fOy
    ## xmm13 = fOz
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fOzc fOzd

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

        ## move j H1 coordinates to local temp variables 
        movq  nb332_pos(%rbp),%rsi
    movlps 12(%rsi,%rax,4),%xmm0    ## jxH1a jyH1a  -   -
    movlps 12(%rsi,%rcx,4),%xmm1    ## jxH1c jyH1c  -   -
    movhps 12(%rsi,%rbx,4),%xmm0    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rsi,%rdx,4),%xmm1    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rsi,%rax,4),%xmm2    ## jzH1a  -  -  -
    movss  20(%rsi,%rcx,4),%xmm3    ## jzH1c  -  -  -
    movhps 20(%rsi,%rbx,4),%xmm2    ## jzH1a  -  jzH1b  -
    movhps 20(%rsi,%rdx,4),%xmm3    ## jzH1c  -  jzH1d -

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

    subps nb332_ixO(%rsp),%xmm0
    subps nb332_iyO(%rsp),%xmm1
    subps nb332_izO(%rsp),%xmm2
    subps nb332_ixH1(%rsp),%xmm3
    subps nb332_iyH1(%rsp),%xmm4
    subps nb332_izH1(%rsp),%xmm5
    subps nb332_ixH2(%rsp),%xmm6
    subps nb332_iyH2(%rsp),%xmm7
    subps nb332_izH2(%rsp),%xmm8

        movaps %xmm0,nb332_dxOH1(%rsp)
        movaps %xmm1,nb332_dyOH1(%rsp)
        movaps %xmm2,nb332_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb332_dxH1H1(%rsp)
        movaps %xmm4,nb332_dyH1H1(%rsp)
        movaps %xmm5,nb332_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb332_dxH2H1(%rsp)
        movaps %xmm7,nb332_dyH2H1(%rsp)
        movaps %xmm8,nb332_dzH2H1(%rsp)
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

        movaps  nb332_three(%rsp),%xmm9
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

        movaps  nb332_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH1
        mulps   %xmm4,%xmm10 ## rinvH1H1
    mulps   %xmm4,%xmm11 ## rinvH2H1

        movaps  %xmm9,nb332_rinvOH1(%rsp)
        movaps  %xmm10,nb332_rinvH1H1(%rsp)
        movaps  %xmm11,nb332_rinvH2H1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb332_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    movq nb332_VFtab(%rbp),%rsi

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

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm1,%xmm10
    movaps  %xmm4,%xmm11
    movaps  %xmm7,%xmm12
    pslld   $1,%xmm1
    pslld   $1,%xmm4
    pslld   $1,%xmm7
    paddd   %xmm10,%xmm1
    paddd   %xmm11,%xmm4
    paddd   %xmm12,%xmm7

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

    movaps    %xmm0,nb332_epsO(%rsp)
    movaps    %xmm3,nb332_epsH1(%rsp)
    movaps    %xmm6,nb332_epsH2(%rsp)


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

    movaps nb332_epsO(%rsp),%xmm12
    movaps nb332_epsH1(%rsp),%xmm13
    movaps nb332_epsH2(%rsp),%xmm14

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
    movaps nb332_qqOH(%rsp),%xmm12
    movaps nb332_qqHH(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    ## accumulate vctot
    addps  nb332_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb332_vctot(%rsp)

    movaps nb332_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j H1 forces to local temp variables 
    movq nb332_faction(%rbp),%rdi
    movlps 12(%rdi,%rax,4),%xmm11    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm12    ## jxH1c jyH1c  -   -

    ## xmm11: jxH1a jyH1a jxH1b jyH1b 
    ## xmm12: jxH1c jyH1c jxH1d jyH1d
    ## xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    movhps 12(%rdi,%rbx,4),%xmm11    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm12    ## jxH1c jyH1c jxH1d jyH1d 

    mulps  nb332_rinvOH1(%rsp),%xmm3
    mulps  nb332_rinvH1H1(%rsp),%xmm7
    mulps  nb332_rinvH2H1(%rsp),%xmm10

    movss  20(%rdi,%rax,4),%xmm13    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm14    ## jzH1c  -  -  -
    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movhps 20(%rdi,%rbx,4),%xmm13    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm14    ## jzH1c  -  jzH1d -

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

        mulps nb332_dxOH1(%rsp),%xmm0
        mulps nb332_dyOH1(%rsp),%xmm1
        mulps nb332_dzOH1(%rsp),%xmm2
        mulps nb332_dxH1H1(%rsp),%xmm3
        mulps nb332_dyH1H1(%rsp),%xmm4
        mulps nb332_dzH1H1(%rsp),%xmm5
        mulps nb332_dxH2H1(%rsp),%xmm6
        mulps nb332_dyH2H1(%rsp),%xmm7
        mulps nb332_dzH2H1(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb332_fixO(%rsp),%xmm0
    addps nb332_fiyO(%rsp),%xmm1
    addps nb332_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb332_fixH1(%rsp),%xmm3
    addps nb332_fiyH1(%rsp),%xmm4
    addps nb332_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb332_fixH2(%rsp),%xmm6
    addps nb332_fiyH2(%rsp),%xmm7
    addps nb332_fizH2(%rsp),%xmm8

    movaps %xmm0,nb332_fixO(%rsp)
    movaps %xmm1,nb332_fiyO(%rsp)
    movaps %xmm2,nb332_fizO(%rsp)
    movaps %xmm3,nb332_fixH1(%rsp)
    movaps %xmm4,nb332_fiyH1(%rsp)
    movaps %xmm5,nb332_fizH1(%rsp)
    movaps %xmm6,nb332_fixH2(%rsp)
    movaps %xmm7,nb332_fiyH2(%rsp)
    movaps %xmm8,nb332_fizH2(%rsp)

    ## xmm11 = fH1x
    ## xmm12 = fH1y
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

        movq  nb332_pos(%rbp),%rsi
        ## move j H2 coordinates to local temp variables 
    movlps 24(%rsi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rsi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rsi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rsi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rsi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rsi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movss  32(%rsi,%rbx,4),%xmm5    ## jzH2b  -  -  -
    movss  32(%rsi,%rdx,4),%xmm6    ## jzH2d  -  -  -
    movlhps %xmm5,%xmm2 ## jzH2a  -  jzH2b  -
    movlhps %xmm6,%xmm3 ## jzH2c  -  jzH2d -

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

    subps nb332_ixO(%rsp),%xmm0
    subps nb332_iyO(%rsp),%xmm1
    subps nb332_izO(%rsp),%xmm2
    subps nb332_ixH1(%rsp),%xmm3
    subps nb332_iyH1(%rsp),%xmm4
    subps nb332_izH1(%rsp),%xmm5
    subps nb332_ixH2(%rsp),%xmm6
    subps nb332_iyH2(%rsp),%xmm7
    subps nb332_izH2(%rsp),%xmm8

        movaps %xmm0,nb332_dxOH2(%rsp)
        movaps %xmm1,nb332_dyOH2(%rsp)
        movaps %xmm2,nb332_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb332_dxH1H2(%rsp)
        movaps %xmm4,nb332_dyH1H2(%rsp)
        movaps %xmm5,nb332_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb332_dxH2H2(%rsp)
        movaps %xmm7,nb332_dyH2H2(%rsp)
        movaps %xmm8,nb332_dzH2H2(%rsp)
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

        movaps  nb332_three(%rsp),%xmm9
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

        movaps  nb332_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH2
        mulps   %xmm4,%xmm10 ## rinvH1H2
    mulps   %xmm4,%xmm11 ## rinvH2H2

        movaps  %xmm9,nb332_rinvOH2(%rsp)
        movaps  %xmm10,nb332_rinvH1H2(%rsp)
        movaps  %xmm11,nb332_rinvH2H2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb332_tsc(%rsp),%xmm1
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

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm1,%xmm10
    movaps  %xmm4,%xmm11
    movaps  %xmm7,%xmm12
    pslld   $1,%xmm1
    pslld   $1,%xmm4
    pslld   $1,%xmm7
    paddd   %xmm10,%xmm1
    paddd   %xmm11,%xmm4
    paddd   %xmm12,%xmm7

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

    movq nb332_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb332_epsO(%rsp)
    movaps    %xmm3,nb332_epsH1(%rsp)
    movaps    %xmm6,nb332_epsH2(%rsp)

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

    movaps nb332_epsO(%rsp),%xmm12
    movaps nb332_epsH1(%rsp),%xmm13
    movaps nb332_epsH2(%rsp),%xmm14

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
    movaps nb332_qqOH(%rsp),%xmm12
    movaps nb332_qqHH(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    ## accumulate vctot
    addps  nb332_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb332_vctot(%rsp)

    movaps nb332_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx


        ## move j H2 forces to local temp variables 
    movq nb332_faction(%rbp),%rdi
    movlps 24(%rdi,%rax,4),%xmm11    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm12    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm11    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm12    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm13    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm14    ## jzH2c  -  -  -
    movss  32(%rdi,%rbx,4),%xmm1    ## jzH2b  -  -  -
    movss  32(%rdi,%rdx,4),%xmm2    ## jzH2d  -  -  -
    movlhps %xmm1,%xmm13 ## jzH2a  -  jzH2b  -
    movlhps %xmm2,%xmm14 ## jzH2c  -  jzH2d -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm11: jxH2a jyH2a jxH2b jyH2b 
    ## xmm12: jxH2c jyH2c jxH2d jyH2d
    ## xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb332_rinvOH2(%rsp),%xmm3
    mulps  nb332_rinvH1H2(%rsp),%xmm7
    mulps  nb332_rinvH2H2(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb332_dxOH2(%rsp),%xmm0
        mulps nb332_dyOH2(%rsp),%xmm1
        mulps nb332_dzOH2(%rsp),%xmm2
        mulps nb332_dxH1H2(%rsp),%xmm3
        mulps nb332_dyH1H2(%rsp),%xmm4
        mulps nb332_dzH1H2(%rsp),%xmm5
        mulps nb332_dxH2H2(%rsp),%xmm6
        mulps nb332_dyH2H2(%rsp),%xmm7
        mulps nb332_dzH2H2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb332_fixO(%rsp),%xmm0
    addps nb332_fiyO(%rsp),%xmm1
    addps nb332_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb332_fixH1(%rsp),%xmm3
    addps nb332_fiyH1(%rsp),%xmm4
    addps nb332_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb332_fixH2(%rsp),%xmm6
    addps nb332_fiyH2(%rsp),%xmm7
    addps nb332_fizH2(%rsp),%xmm8

    movaps %xmm0,nb332_fixO(%rsp)
    movaps %xmm1,nb332_fiyO(%rsp)
    movaps %xmm2,nb332_fizO(%rsp)
    movaps %xmm3,nb332_fixH1(%rsp)
    movaps %xmm4,nb332_fiyH1(%rsp)
    movaps %xmm5,nb332_fizH1(%rsp)
    movaps %xmm6,nb332_fixH2(%rsp)
    movaps %xmm7,nb332_fiyH2(%rsp)
    movaps %xmm8,nb332_fizH2(%rsp)

    ## xmm11 = fH2x
    ## xmm12 = fH2y
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

        ## should we do one more iteration? 
        subl $4,nb332_innerk(%rsp)
        jl    _nb_kernel332_x86_64_sse.nb332_single_check
        jmp   _nb_kernel332_x86_64_sse.nb332_unroll_loop
_nb_kernel332_x86_64_sse.nb332_single_check: 
        addl $4,nb332_innerk(%rsp)
        jnz   _nb_kernel332_x86_64_sse.nb332_single_loop
        jmp   _nb_kernel332_x86_64_sse.nb332_updateouterdata
_nb_kernel332_x86_64_sse.nb332_single_loop: 
        movq  nb332_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb332_innerjjnr(%rsp)

        movq nb332_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm0,%xmm0
        xorps %xmm1,%xmm1
        xorps %xmm2,%xmm2

        movss (%rsi,%rax,4),%xmm0               ## jxO  -  -  -
        movss 4(%rsi,%rax,4),%xmm1              ## jyO  -  -  -
        movss 8(%rsi,%rax,4),%xmm2              ## jzO  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm5            ## xmm5 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm5,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movlhps %xmm6,%xmm0                     ## xmm0 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm1 ## 11100100     ;# xmm1 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm2 ## 01000100     ;# xmm2 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm0,nb332_jxO(%rsp)
        movaps %xmm1,nb332_jyO(%rsp)
        movaps %xmm2,nb332_jzO(%rsp)
        subps  nb332_ixO(%rsp),%xmm0
        subps  nb332_iyO(%rsp),%xmm1
        subps  nb332_izO(%rsp),%xmm2
        movaps %xmm0,nb332_dxOO(%rsp)
        movaps %xmm1,nb332_dyOO(%rsp)
        movaps %xmm2,nb332_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb332_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb332_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb332_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

        movq nb332_VFtab(%rbp),%rsi

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    mulps  nb332_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb332_qqOO(%rsp),%xmm3
        movhps  nb332_qqOH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 

    addps  nb332_vctot(%rsp),%xmm5
    movaps %xmm5,nb332_vctot(%rsp)
    ## put scalar force on stack temporarily 
    movaps %xmm3,nb332_fstmp(%rsp)

    ## dispersion 
        movss  16(%rsi,%rbx,4),%xmm4
        movss  20(%rsi,%rbx,4),%xmm5
        movss  24(%rsi,%rbx,4),%xmm6
        movss  28(%rsi,%rbx,4),%xmm7
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb332_two(%rsp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 
        xorps  %xmm4,%xmm4
    movss  nb332_c6(%rsp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    addps  nb332_fstmp(%rsp),%xmm7   ## add to fscal 

    ## put scalar force on stack Update Vvdwtot directly 
    addps  nb332_Vvdwtot(%rsp),%xmm5
    movaps %xmm7,nb332_fstmp(%rsp)
    movaps %xmm5,nb332_Vvdwtot(%rsp)

    ## repulsion 
        movss  32(%rsi,%rbx,4),%xmm4
        movss  36(%rsi,%rbx,4),%xmm5
        movss  40(%rsi,%rbx,4),%xmm6
        movss  44(%rsi,%rbx,4),%xmm7
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb332_two(%rsp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

        xorps  %xmm4,%xmm4
    movss  nb332_c12(%rsp),%xmm4
    mulps  %xmm4,%xmm7 ## fijR 
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb332_fstmp(%rsp),%xmm7

    addps  nb332_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332_Vvdwtot(%rsp)
    xorps  %xmm1,%xmm1

    mulps nb332_tsc(%rsp),%xmm7
    mulps %xmm0,%xmm7
    subps  %xmm7,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        mulps   nb332_dxOO(%rsp),%xmm0
        mulps   nb332_dyOO(%rsp),%xmm1
        mulps   nb332_dzOO(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb332_fjxO(%rsp)
        movaps  %xmm4,nb332_fjyO(%rsp)
        movaps  %xmm5,nb332_fjzO(%rsp)
        addps   nb332_fixO(%rsp),%xmm0
        addps   nb332_fiyO(%rsp),%xmm1
        addps   nb332_fizO(%rsp),%xmm2
        movaps  %xmm0,nb332_fixO(%rsp)
        movaps  %xmm1,nb332_fiyO(%rsp)
        movaps  %xmm2,nb332_fizO(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb332_jxO(%rsp),%xmm0
    movaps  nb332_jyO(%rsp),%xmm1
    movaps  nb332_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  nb332_ixH1(%rsp),%xmm0
        subps  nb332_iyH1(%rsp),%xmm1
        subps  nb332_izH1(%rsp),%xmm2
        subps  nb332_ixH2(%rsp),%xmm3
        subps  nb332_iyH2(%rsp),%xmm4
        subps  nb332_izH2(%rsp),%xmm5
    movaps %xmm0,nb332_dxH1O(%rsp)
        movaps %xmm1,nb332_dyH1O(%rsp)
        movaps %xmm2,nb332_dzH1O(%rsp)
        movaps %xmm3,nb332_dxH2O(%rsp)
        movaps %xmm4,nb332_dyH2O(%rsp)
        movaps %xmm5,nb332_dzH2O(%rsp)
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
        movaps %xmm4,nb332_rsqH2O(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb332_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb332_rinvH2O(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb332_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    mulps  nb332_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb332_qqOH(%rsp),%xmm3
        movhps  nb332_qqHH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 
    addps  nb332_vctot(%rsp),%xmm5
    movaps %xmm5,nb332_vctot(%rsp)

    xorps  %xmm1,%xmm1

    mulps nb332_tsc(%rsp),%xmm3
    mulps %xmm0,%xmm3
    subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb332_dxH1O(%rsp),%xmm0
        mulps   nb332_dyH1O(%rsp),%xmm1
        mulps   nb332_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb332_fjxO(%rsp),%xmm3
        movaps  nb332_fjyO(%rsp),%xmm4
        movaps  nb332_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb332_fjxO(%rsp)
        movaps  %xmm4,nb332_fjyO(%rsp)
        movaps  %xmm5,nb332_fjzO(%rsp)
        addps   nb332_fixH1(%rsp),%xmm0
        addps   nb332_fiyH1(%rsp),%xmm1
        addps   nb332_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb332_fixH1(%rsp)
        movaps  %xmm1,nb332_fiyH1(%rsp)
        movaps  %xmm2,nb332_fizH1(%rsp)
        ## do table for H2 - j water interaction 
        movaps nb332_rinvH2O(%rsp),%xmm0
        movaps nb332_rsqH2O(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb332_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    mulps  nb332_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb332_qqOH(%rsp),%xmm3
        movhps  nb332_qqHH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 
    addps  nb332_vctot(%rsp),%xmm5
    movaps %xmm5,nb332_vctot(%rsp)

    xorps  %xmm1,%xmm1

    mulps nb332_tsc(%rsp),%xmm3
    mulps %xmm0,%xmm3
    subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb332_dxH2O(%rsp),%xmm0
        mulps   nb332_dyH2O(%rsp),%xmm1
        mulps   nb332_dzH2O(%rsp),%xmm2
        movaps  nb332_fjxO(%rsp),%xmm3
        movaps  nb332_fjyO(%rsp),%xmm4
        movaps  nb332_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb332_faction(%rbp),%rsi
        movaps  %xmm3,nb332_fjxO(%rsp)
        movaps  %xmm4,nb332_fjyO(%rsp)
        movaps  %xmm5,nb332_fjzO(%rsp)
        addps   nb332_fixH2(%rsp),%xmm0
        addps   nb332_fiyH2(%rsp),%xmm1
        addps   nb332_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb332_fixH2(%rsp)
        movaps  %xmm1,nb332_fiyH2(%rsp)
        movaps  %xmm2,nb332_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb332_fjxO(%rsp),%xmm3
        movaps  nb332_fjyO(%rsp),%xmm4
        movaps  nb332_fjzO(%rsp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## 00000010
        shufps $3,%xmm7,%xmm7 ## 00000011
        addss   8(%rsi,%rax,4),%xmm5
        addss   20(%rsi,%rax,4),%xmm6
        addss   32(%rsi,%rax,4),%xmm7
        movss   %xmm5,8(%rsi,%rax,4)
        movss   %xmm6,20(%rsi,%rax,4)
        movss   %xmm7,32(%rsi,%rax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%rsi,%rax,4)
        movlps  %xmm1,12(%rsi,%rax,4)
        movhps  %xmm1,24(%rsi,%rax,4)

        decl nb332_innerk(%rsp)
        jz    _nb_kernel332_x86_64_sse.nb332_updateouterdata
        jmp   _nb_kernel332_x86_64_sse.nb332_single_loop
_nb_kernel332_x86_64_sse.nb332_updateouterdata: 
        movl  nb332_ii3(%rsp),%ecx
        movq  nb332_faction(%rbp),%rdi
        movq  nb332_fshift(%rbp),%rsi
        movl  nb332_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb332_fixO(%rsp),%xmm0
        movaps nb332_fiyO(%rsp),%xmm1
        movaps nb332_fizO(%rsp),%xmm2

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
        movaps nb332_fixH1(%rsp),%xmm0
        movaps nb332_fiyH1(%rsp),%xmm1
        movaps nb332_fizH1(%rsp),%xmm2

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
        movaps nb332_fixH2(%rsp),%xmm0
        movaps nb332_fiyH2(%rsp),%xmm1
        movaps nb332_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 i xmm0-xmm2 

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
        movl nb332_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb332_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb332_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb332_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb332_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb332_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb332_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel332_x86_64_sse.nb332_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb332_n(%rsp)
        jmp _nb_kernel332_x86_64_sse.nb332_outer
_nb_kernel332_x86_64_sse.nb332_outerend: 
        ## check if more outer neighborlists remain
        movl  nb332_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel332_x86_64_sse.nb332_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel332_x86_64_sse.nb332_threadloop
_nb_kernel332_x86_64_sse.nb332_end: 
        movl nb332_nouter(%rsp),%eax
        movl nb332_ninner(%rsp),%ebx
        movq nb332_outeriter(%rbp),%rcx
        movq nb332_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1592,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel332nf_x86_64_sse
.globl _nb_kernel332nf_x86_64_sse
nb_kernel332nf_x86_64_sse:      
_nb_kernel332nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb332nf_fshift, 16
.set nb332nf_gid, 24
.set nb332nf_pos, 32
.set nb332nf_faction, 40
.set nb332nf_charge, 48
.set nb332nf_p_facel, 56
.set nb332nf_argkrf, 64
.set nb332nf_argcrf, 72
.set nb332nf_Vc, 80
.set nb332nf_type, 88
.set nb332nf_p_ntype, 96
.set nb332nf_vdwparam, 104
.set nb332nf_Vvdw, 112
.set nb332nf_p_tabscale, 120
.set nb332nf_VFtab, 128
.set nb332nf_invsqrta, 136
.set nb332nf_dvda, 144
.set nb332nf_p_gbtabscale, 152
.set nb332nf_GBtab, 160
.set nb332nf_p_nthreads, 168
.set nb332nf_count, 176
.set nb332nf_mtx, 184
.set nb332nf_outeriter, 192
.set nb332nf_inneriter, 200
.set nb332nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb332nf_ixO, 0
.set nb332nf_iyO, 16
.set nb332nf_izO, 32
.set nb332nf_ixH1, 48
.set nb332nf_iyH1, 64
.set nb332nf_izH1, 80
.set nb332nf_ixH2, 96
.set nb332nf_iyH2, 112
.set nb332nf_izH2, 128
.set nb332nf_jxO, 144
.set nb332nf_jyO, 160
.set nb332nf_jzO, 176
.set nb332nf_jxH1, 192
.set nb332nf_jyH1, 208
.set nb332nf_jzH1, 224
.set nb332nf_jxH2, 240
.set nb332nf_jyH2, 256
.set nb332nf_jzH2, 272
.set nb332nf_qqOO, 288
.set nb332nf_qqOH, 304
.set nb332nf_qqHH, 320
.set nb332nf_tsc, 336
.set nb332nf_c6, 352
.set nb332nf_c12, 368
.set nb332nf_vctot, 384
.set nb332nf_Vvdwtot, 400
.set nb332nf_half, 416
.set nb332nf_three, 432
.set nb332nf_rsqOO, 448
.set nb332nf_rsqOH1, 464
.set nb332nf_rsqOH2, 480
.set nb332nf_rsqH1O, 496
.set nb332nf_rsqH1H1, 512
.set nb332nf_rsqH1H2, 528
.set nb332nf_rsqH2O, 544
.set nb332nf_rsqH2H1, 560
.set nb332nf_rsqH2H2, 576
.set nb332nf_rinvOO, 592
.set nb332nf_rinvOH1, 608
.set nb332nf_rinvOH2, 624
.set nb332nf_rinvH1O, 640
.set nb332nf_rinvH1H1, 656
.set nb332nf_rinvH1H2, 672
.set nb332nf_rinvH2O, 688
.set nb332nf_rinvH2H1, 704
.set nb332nf_rinvH2H2, 720
.set nb332nf_is3, 736
.set nb332nf_ii3, 740
.set nb332nf_nri, 744
.set nb332nf_iinr, 752
.set nb332nf_jindex, 760
.set nb332nf_jjnr, 768
.set nb332nf_shift, 776
.set nb332nf_shiftvec, 784
.set nb332nf_facel, 792
.set nb332nf_innerjjnr, 800
.set nb332nf_innerk, 808
.set nb332nf_n, 812
.set nb332nf_nn1, 816
.set nb332nf_nouter, 820
.set nb332nf_ninner, 824
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms
        subq $840,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb332nf_nouter(%rsp)
        movl %eax,nb332nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb332nf_nri(%rsp)
        movq %rsi,nb332nf_iinr(%rsp)
        movq %rdx,nb332nf_jindex(%rsp)
        movq %rcx,nb332nf_jjnr(%rsp)
        movq %r8,nb332nf_shift(%rsp)
        movq %r9,nb332nf_shiftvec(%rsp)
        movq nb332nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb332nf_facel(%rsp)

        movq nb332nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb332nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb332nf_half(%rsp)
        movss nb332nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb332nf_half(%rsp)
        movaps %xmm3,nb332nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb332nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb332nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb332nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb332nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb332nf_qqOO(%rsp)
        movaps %xmm4,nb332nf_qqOH(%rsp)
        movaps %xmm5,nb332nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb332nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb332nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb332nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb332nf_c6(%rsp)
        movaps %xmm1,nb332nf_c12(%rsp)

_nb_kernel332nf_x86_64_sse.nb332nf_threadloop: 
        movq  nb332nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel332nf_x86_64_sse.nb332nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel332nf_x86_64_sse.nb332nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb332nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb332nf_n(%rsp)
        movl %ebx,nb332nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel332nf_x86_64_sse.nb332nf_outerstart
        jmp _nb_kernel332nf_x86_64_sse.nb332nf_end

_nb_kernel332nf_x86_64_sse.nb332nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb332nf_nouter(%rsp),%ebx
        movl %ebx,nb332nf_nouter(%rsp)

_nb_kernel332nf_x86_64_sse.nb332nf_outer: 
        movq  nb332nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb332nf_is3(%rsp)            ## store is3 

        movq  nb332nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb332nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb332nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb332nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb332nf_ixO(%rsp)
        movaps %xmm4,nb332nf_iyO(%rsp)
        movaps %xmm5,nb332nf_izO(%rsp)

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
        movaps %xmm0,nb332nf_ixH1(%rsp)
        movaps %xmm1,nb332nf_iyH1(%rsp)
        movaps %xmm2,nb332nf_izH1(%rsp)
        movaps %xmm3,nb332nf_ixH2(%rsp)
        movaps %xmm4,nb332nf_iyH2(%rsp)
        movaps %xmm5,nb332nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb332nf_vctot(%rsp)
        movaps %xmm4,nb332nf_Vvdwtot(%rsp)

        movq  nb332nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb332nf_pos(%rbp),%rsi
        movq  nb332nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb332nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb332nf_ninner(%rsp),%ecx
        movl  %ecx,nb332nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb332nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel332nf_x86_64_sse.nb332nf_unroll_loop
        jmp   _nb_kernel332nf_x86_64_sse.nb332nf_single_check
_nb_kernel332nf_x86_64_sse.nb332nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb332nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb332nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb332nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j coordinates to local temp variables 
        movlps (%rsi,%rax,4),%xmm2
        movlps 12(%rsi,%rax,4),%xmm3
        movlps 24(%rsi,%rax,4),%xmm4

        movlps (%rsi,%rbx,4),%xmm5
        movlps 12(%rsi,%rbx,4),%xmm6
        movlps 24(%rsi,%rbx,4),%xmm7

        movhps (%rsi,%rcx,4),%xmm2
        movhps 12(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rcx,4),%xmm4

        movhps (%rsi,%rdx,4),%xmm5
        movhps 12(%rsi,%rdx,4),%xmm6
        movhps 24(%rsi,%rdx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd  
        movaps %xmm0,nb332nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb332nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb332nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb332nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb332nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb332nf_jyH2(%rsp)

        movss  8(%rsi,%rax,4),%xmm0
        movss  20(%rsi,%rax,4),%xmm1
        movss  32(%rsi,%rax,4),%xmm2

        movss  8(%rsi,%rcx,4),%xmm3
        movss  20(%rsi,%rcx,4),%xmm4
        movss  32(%rsi,%rcx,4),%xmm5

        movhps 4(%rsi,%rbx,4),%xmm0
        movhps 16(%rsi,%rbx,4),%xmm1
        movhps 28(%rsi,%rbx,4),%xmm2

        movhps 4(%rsi,%rdx,4),%xmm3
        movhps 16(%rsi,%rdx,4),%xmm4
        movhps 28(%rsi,%rdx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## 11001100
        shufps $204,%xmm4,%xmm1 ## 11001100
        shufps $204,%xmm5,%xmm2 ## 11001100
        movaps %xmm0,nb332nf_jzO(%rsp)
        movaps %xmm1,nb332nf_jzH1(%rsp)
        movaps %xmm2,nb332nf_jzH2(%rsp)

        movaps nb332nf_ixO(%rsp),%xmm0
        movaps nb332nf_iyO(%rsp),%xmm1
        movaps nb332nf_izO(%rsp),%xmm2
        movaps nb332nf_ixO(%rsp),%xmm3
        movaps nb332nf_iyO(%rsp),%xmm4
        movaps nb332nf_izO(%rsp),%xmm5
        subps  nb332nf_jxO(%rsp),%xmm0
        subps  nb332nf_jyO(%rsp),%xmm1
        subps  nb332nf_jzO(%rsp),%xmm2
        subps  nb332nf_jxH1(%rsp),%xmm3
        subps  nb332nf_jyH1(%rsp),%xmm4
        subps  nb332nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb332nf_rsqOO(%rsp)
        movaps %xmm3,nb332nf_rsqOH1(%rsp)

        movaps nb332nf_ixO(%rsp),%xmm0
        movaps nb332nf_iyO(%rsp),%xmm1
        movaps nb332nf_izO(%rsp),%xmm2
        movaps nb332nf_ixH1(%rsp),%xmm3
        movaps nb332nf_iyH1(%rsp),%xmm4
        movaps nb332nf_izH1(%rsp),%xmm5
        subps  nb332nf_jxH2(%rsp),%xmm0
        subps  nb332nf_jyH2(%rsp),%xmm1
        subps  nb332nf_jzH2(%rsp),%xmm2
        subps  nb332nf_jxO(%rsp),%xmm3
        subps  nb332nf_jyO(%rsp),%xmm4
        subps  nb332nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb332nf_rsqOH2(%rsp)
        movaps %xmm3,nb332nf_rsqH1O(%rsp)

        movaps nb332nf_ixH1(%rsp),%xmm0
        movaps nb332nf_iyH1(%rsp),%xmm1
        movaps nb332nf_izH1(%rsp),%xmm2
        movaps nb332nf_ixH1(%rsp),%xmm3
        movaps nb332nf_iyH1(%rsp),%xmm4
        movaps nb332nf_izH1(%rsp),%xmm5
        subps  nb332nf_jxH1(%rsp),%xmm0
        subps  nb332nf_jyH1(%rsp),%xmm1
        subps  nb332nf_jzH1(%rsp),%xmm2
        subps  nb332nf_jxH2(%rsp),%xmm3
        subps  nb332nf_jyH2(%rsp),%xmm4
        subps  nb332nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb332nf_rsqH1H1(%rsp)
        movaps %xmm3,nb332nf_rsqH1H2(%rsp)

        movaps nb332nf_ixH2(%rsp),%xmm0
        movaps nb332nf_iyH2(%rsp),%xmm1
        movaps nb332nf_izH2(%rsp),%xmm2
        movaps nb332nf_ixH2(%rsp),%xmm3
        movaps nb332nf_iyH2(%rsp),%xmm4
        movaps nb332nf_izH2(%rsp),%xmm5
        subps  nb332nf_jxO(%rsp),%xmm0
        subps  nb332nf_jyO(%rsp),%xmm1
        subps  nb332nf_jzO(%rsp),%xmm2
        subps  nb332nf_jxH1(%rsp),%xmm3
        subps  nb332nf_jyH1(%rsp),%xmm4
        subps  nb332nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb332nf_rsqH2O(%rsp)
        movaps %xmm4,nb332nf_rsqH2H1(%rsp)

        movaps nb332nf_ixH2(%rsp),%xmm0
        movaps nb332nf_iyH2(%rsp),%xmm1
        movaps nb332nf_izH2(%rsp),%xmm2
        subps  nb332nf_jxH2(%rsp),%xmm0
        subps  nb332nf_jyH2(%rsp),%xmm1
        subps  nb332nf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb332nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb332nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb332nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb332nf_rinvH2H1(%rsp)

        rsqrtps nb332nf_rsqOO(%rsp),%xmm1
        rsqrtps nb332nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb332nf_rsqOO(%rsp),%xmm1
        mulps   nb332nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332nf_half(%rsp),%xmm3
        mulps   nb332nf_half(%rsp),%xmm7
        movaps  %xmm3,nb332nf_rinvOO(%rsp)
        movaps  %xmm7,nb332nf_rinvOH1(%rsp)

        rsqrtps nb332nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb332nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb332nf_rsqOH2(%rsp),%xmm1
        mulps   nb332nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332nf_half(%rsp),%xmm3
        mulps   nb332nf_half(%rsp),%xmm7
        movaps  %xmm3,nb332nf_rinvOH2(%rsp)
        movaps  %xmm7,nb332nf_rinvH1O(%rsp)

        rsqrtps nb332nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb332nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb332nf_rsqH1H1(%rsp),%xmm1
        mulps   nb332nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332nf_half(%rsp),%xmm3
        mulps   nb332nf_half(%rsp),%xmm7
        movaps  %xmm3,nb332nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb332nf_rinvH1H2(%rsp)

        rsqrtps nb332nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb332nf_three(%rsp),%xmm3
        mulps   nb332nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb332nf_half(%rsp),%xmm3
        movaps  %xmm3,nb332nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movaps nb332nf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7

    movd %eax,%mm0
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    movq nb332nf_VFtab(%rbp),%rsi
    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqOO(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

    ## dispersion 
    movlps 16(%rsi,%rax,4),%xmm5
    movlps 16(%rsi,%rcx,4),%xmm7
    movhps 16(%rsi,%rbx,4),%xmm5
    movhps 16(%rsi,%rdx,4),%xmm7    ## got half dispersion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 24(%rsi,%rax,4),%xmm7
    movlps 24(%rsi,%rcx,4),%xmm3
    movhps 24(%rsi,%rbx,4),%xmm7
    movhps 24(%rsi,%rdx,4),%xmm3    ## other half of dispersion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb332nf_c6(%rsp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 

    ## put scalar force on stack Update Vvdwtot directly 
    addps  nb332nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332nf_Vvdwtot(%rsp)

    ## repulsion 
    movlps 32(%rsi,%rax,4),%xmm5
    movlps 32(%rsi,%rcx,4),%xmm7
    movhps 32(%rsi,%rbx,4),%xmm5
    movhps 32(%rsi,%rdx,4),%xmm7    ## got half repulsion table 
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## 10001000
    shufps $221,%xmm7,%xmm5 ## 11011101

    movlps 40(%rsi,%rax,4),%xmm7
    movlps 40(%rsi,%rcx,4),%xmm3
    movhps 40(%rsi,%rbx,4),%xmm7
    movhps 40(%rsi,%rdx,4),%xmm3    ## other half of repulsion table 
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## 10001000
    shufps $221,%xmm3,%xmm7 ## 11011101
    ## table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb332nf_c12(%rsp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb332nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332nf_Vvdwtot(%rsp)

        ## O-H1 interaction 
        movaps nb332nf_rinvOH1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqOH1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## O-H2 interaction  
        movaps nb332nf_rinvOH2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqOH2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H1-O interaction 
        movaps nb332nf_rinvH1O(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH1O(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H1-H1 interaction 
        movaps nb332nf_rinvH1H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H1-H2 interaction 
        movaps nb332nf_rinvH1H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H2-O interaction 
        movaps nb332nf_rinvH2O(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH2O(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H2-H1 interaction 
        movaps nb332nf_rinvH2H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## H2-H2 interaction 
        movaps nb332nf_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb332nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1
        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
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

    lea  (%rax,%rax,2),%rax
    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
    movaps nb332nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb332nf_innerk(%rsp)
        jl    _nb_kernel332nf_x86_64_sse.nb332nf_single_check
        jmp   _nb_kernel332nf_x86_64_sse.nb332nf_unroll_loop
_nb_kernel332nf_x86_64_sse.nb332nf_single_check: 
        addl $4,nb332nf_innerk(%rsp)
        jnz   _nb_kernel332nf_x86_64_sse.nb332nf_single_loop
        jmp   _nb_kernel332nf_x86_64_sse.nb332nf_updateouterdata
_nb_kernel332nf_x86_64_sse.nb332nf_single_loop: 
        movq  nb332nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb332nf_innerjjnr(%rsp)

        movq nb332nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%rsi,%rax,4),%xmm3               ## jxO  -  -  -
        movss 4(%rsi,%rax,4),%xmm4              ## jyO  -  -  -
        movss 8(%rsi,%rax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%rsi,%rax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%rsi,%rax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%rsi,%rax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## 11011000      ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb332nf_ixO(%rsp),%xmm0
        movaps  nb332nf_iyO(%rsp),%xmm1
        movaps  nb332nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb332nf_jxO(%rsp)
        movaps %xmm4,nb332nf_jyO(%rsp)
        movaps %xmm5,nb332nf_jzO(%rsp)
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
        movaps  nb332nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb332nf_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb332nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

        movq nb332nf_VFtab(%rbp),%rsi

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
        movss   nb332nf_qqOO(%rsp),%xmm3
        movhps  nb332nf_qqOH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point xmm5 contains vcoul 

    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

    ## dispersion 
        movss  16(%rsi,%rbx,4),%xmm4
        movss  20(%rsi,%rbx,4),%xmm5
        movss  24(%rsi,%rbx,4),%xmm6
        movss  28(%rsi,%rbx,4),%xmm7
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 
        xorps  %xmm4,%xmm4
    movss  nb332nf_c6(%rsp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 
    ## put scalar force on stack 
    addps  nb332nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332nf_Vvdwtot(%rsp)

    ## repulsion 
        movss  32(%rsi,%rbx,4),%xmm4
        movss  36(%rsi,%rbx,4),%xmm5
        movss  40(%rsi,%rbx,4),%xmm6
        movss  44(%rsi,%rbx,4),%xmm7
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

        xorps  %xmm4,%xmm4
    movss  nb332nf_c12(%rsp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb332nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb332nf_Vvdwtot(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb332nf_ixH1(%rsp),%xmm0
        movaps  nb332nf_iyH1(%rsp),%xmm1
        movaps  nb332nf_izH1(%rsp),%xmm2
        movaps  nb332nf_ixH2(%rsp),%xmm3
        movaps  nb332nf_iyH2(%rsp),%xmm4
        movaps  nb332nf_izH2(%rsp),%xmm5
        subps   nb332nf_jxO(%rsp),%xmm0
        subps   nb332nf_jyO(%rsp),%xmm1
        subps   nb332nf_jzO(%rsp),%xmm2
        subps   nb332nf_jxO(%rsp),%xmm3
        subps   nb332nf_jyO(%rsp),%xmm4
        subps   nb332nf_jzO(%rsp),%xmm5
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
        movaps %xmm4,nb332nf_rsqH2O(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb332nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb332nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb332nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb332nf_rinvH2O(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb332nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
        movss   nb332nf_qqOH(%rsp),%xmm3
        movhps  nb332nf_qqHH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
    ## at this point xmm5 contains vcoul 
    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)


        ## do table for H2 - j water interaction 
        movaps nb332nf_rinvH2O(%rsp),%xmm0
        movaps nb332nf_rsqH2O(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb332nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $2,%mm6
    pslld $2,%mm7
    movd %mm6,%ebx
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

    lea  (%rbx,%rbx,2),%rbx
    lea  (%rcx,%rcx,2),%rcx
    lea  (%rdx,%rdx,2),%rdx

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
        movss   nb332nf_qqOH(%rsp),%xmm3
        movhps  nb332nf_qqHH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point xmm5 contains vcoul 
    addps  nb332nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb332nf_vctot(%rsp)

        decl nb332nf_innerk(%rsp)
        jz    _nb_kernel332nf_x86_64_sse.nb332nf_updateouterdata
        jmp   _nb_kernel332nf_x86_64_sse.nb332nf_single_loop
_nb_kernel332nf_x86_64_sse.nb332nf_updateouterdata: 
        ## get n from stack
        movl nb332nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb332nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb332nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb332nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb332nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb332nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb332nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel332nf_x86_64_sse.nb332nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb332nf_n(%rsp)
        jmp _nb_kernel332nf_x86_64_sse.nb332nf_outer
_nb_kernel332nf_x86_64_sse.nb332nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb332nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel332nf_x86_64_sse.nb332nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel332nf_x86_64_sse.nb332nf_threadloop
_nb_kernel332nf_x86_64_sse.nb332nf_end: 
        movl nb332nf_nouter(%rsp),%eax
        movl nb332nf_ninner(%rsp),%ebx
        movq nb332nf_outeriter(%rbp),%rcx
        movq nb332nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $840,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret

