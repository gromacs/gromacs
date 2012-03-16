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





.globl nb_kernel132_x86_64_sse
.globl _nb_kernel132_x86_64_sse
nb_kernel132_x86_64_sse:        
_nb_kernel132_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb132_fshift, 16
.set nb132_gid, 24
.set nb132_pos, 32
.set nb132_faction, 40
.set nb132_charge, 48
.set nb132_p_facel, 56
.set nb132_argkrf, 64
.set nb132_argcrf, 72
.set nb132_Vc, 80
.set nb132_type, 88
.set nb132_p_ntype, 96
.set nb132_vdwparam, 104
.set nb132_Vvdw, 112
.set nb132_p_tabscale, 120
.set nb132_VFtab, 128
.set nb132_invsqrta, 136
.set nb132_dvda, 144
.set nb132_p_gbtabscale, 152
.set nb132_GBtab, 160
.set nb132_p_nthreads, 168
.set nb132_count, 176
.set nb132_mtx, 184
.set nb132_outeriter, 192
.set nb132_inneriter, 200
.set nb132_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb132_ixO, 0
.set nb132_iyO, 16
.set nb132_izO, 32
.set nb132_ixH1, 48
.set nb132_iyH1, 64
.set nb132_izH1, 80
.set nb132_ixH2, 96
.set nb132_iyH2, 112
.set nb132_izH2, 128
.set nb132_jxO, 144
.set nb132_jyO, 160
.set nb132_jzO, 176
.set nb132_jxH1, 192
.set nb132_jyH1, 208
.set nb132_jzH1, 224
.set nb132_jxH2, 240
.set nb132_jyH2, 256
.set nb132_jzH2, 272
.set nb132_dxOO, 288
.set nb132_dyOO, 304
.set nb132_dzOO, 320
.set nb132_dxOH1, 336
.set nb132_dyOH1, 352
.set nb132_dzOH1, 368
.set nb132_dxOH2, 384
.set nb132_dyOH2, 400
.set nb132_dzOH2, 416
.set nb132_dxH1O, 432
.set nb132_dyH1O, 448
.set nb132_dzH1O, 464
.set nb132_dxH1H1, 480
.set nb132_dyH1H1, 496
.set nb132_dzH1H1, 512
.set nb132_dxH1H2, 528
.set nb132_dyH1H2, 544
.set nb132_dzH1H2, 560
.set nb132_dxH2O, 576
.set nb132_dyH2O, 592
.set nb132_dzH2O, 608
.set nb132_dxH2H1, 624
.set nb132_dyH2H1, 640
.set nb132_dzH2H1, 656
.set nb132_dxH2H2, 672
.set nb132_dyH2H2, 688
.set nb132_dzH2H2, 704
.set nb132_qqOO, 720
.set nb132_qqOH, 736
.set nb132_qqHH, 752
.set nb132_c6, 768
.set nb132_c12, 784
.set nb132_tsc, 800
.set nb132_fstmp, 816
.set nb132_vctot, 832
.set nb132_Vvdwtot, 848
.set nb132_fixO, 864
.set nb132_fiyO, 880
.set nb132_fizO, 896
.set nb132_fixH1, 912
.set nb132_fiyH1, 928
.set nb132_fizH1, 944
.set nb132_fixH2, 960
.set nb132_fiyH2, 976
.set nb132_fizH2, 992
.set nb132_fjxO, 1008
.set nb132_fjyO, 1024
.set nb132_fjzO, 1040
.set nb132_fjxH1, 1056
.set nb132_fjyH1, 1072
.set nb132_fjzH1, 1088
.set nb132_fjxH2, 1104
.set nb132_fjyH2, 1120
.set nb132_fjzH2, 1136
.set nb132_half, 1152
.set nb132_three, 1168
.set nb132_rsqOO, 1184
.set nb132_rsqOH1, 1200
.set nb132_rsqOH2, 1216
.set nb132_rsqH1O, 1232
.set nb132_rsqH1H1, 1248
.set nb132_rsqH1H2, 1264
.set nb132_rsqH2O, 1280
.set nb132_rsqH2H1, 1296
.set nb132_rsqH2H2, 1312
.set nb132_rinvOO, 1328
.set nb132_rinvOH1, 1344
.set nb132_rinvOH2, 1360
.set nb132_rinvH1O, 1376
.set nb132_rinvH1H1, 1392
.set nb132_rinvH1H2, 1408
.set nb132_rinvH2O, 1424
.set nb132_rinvH2H1, 1440
.set nb132_rinvH2H2, 1456
.set nb132_two, 1472
.set nb132_krf, 1488
.set nb132_crf, 1504
.set nb132_nri, 1520
.set nb132_iinr, 1528
.set nb132_jindex, 1536
.set nb132_jjnr, 1544
.set nb132_shift, 1552
.set nb132_shiftvec, 1560
.set nb132_facel, 1568
.set nb132_innerjjnr, 1576
.set nb132_is3, 1584
.set nb132_ii3, 1588
.set nb132_innerk, 1592
.set nb132_n, 1596
.set nb132_nn1, 1600
.set nb132_nouter, 1604
.set nb132_ninner, 1608

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1624,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb132_nouter(%rsp)
        movl %eax,nb132_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb132_nri(%rsp)
        movq %rsi,nb132_iinr(%rsp)
        movq %rdx,nb132_jindex(%rsp)
        movq %rcx,nb132_jjnr(%rsp)
        movq %r8,nb132_shift(%rsp)
        movq %r9,nb132_shiftvec(%rsp)
        movq nb132_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb132_facel(%rsp)

        movq nb132_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb132_tsc(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb132_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb132_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb132_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb132_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132_qqOO(%rsp)
        movaps %xmm4,nb132_qqOH(%rsp)
        movaps %xmm5,nb132_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb132_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb132_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addq  %rcx,%rdx
        movq  nb132_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb132_c6(%rsp)
        movaps %xmm1,nb132_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb132_half(%rsp)
        movss nb132_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb132_half(%rsp)
        movaps %xmm2,nb132_two(%rsp)
        movaps %xmm3,nb132_three(%rsp)

_nb_kernel132_x86_64_sse.nb132_threadloop: 
        movq  nb132_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel132_x86_64_sse.nb132_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel132_x86_64_sse.nb132_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb132_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb132_n(%rsp)
        movl %ebx,nb132_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel132_x86_64_sse.nb132_outerstart
        jmp _nb_kernel132_x86_64_sse.nb132_end

_nb_kernel132_x86_64_sse.nb132_outerstart: 
        ## ebx contains number of outer iterations
        addl nb132_nouter(%rsp),%ebx
        movl %ebx,nb132_nouter(%rsp)

_nb_kernel132_x86_64_sse.nb132_outer: 
        movq  nb132_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb132_is3(%rsp)      ## store is3 

        movq  nb132_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb132_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb132_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb132_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132_ixO(%rsp)
        movaps %xmm4,nb132_iyO(%rsp)
        movaps %xmm5,nb132_izO(%rsp)

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
        movaps %xmm0,nb132_ixH1(%rsp)
        movaps %xmm1,nb132_iyH1(%rsp)
        movaps %xmm2,nb132_izH1(%rsp)
        movaps %xmm3,nb132_ixH2(%rsp)
        movaps %xmm4,nb132_iyH2(%rsp)
        movaps %xmm5,nb132_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb132_vctot(%rsp)
        movaps %xmm4,nb132_Vvdwtot(%rsp)
        movaps %xmm4,nb132_fixO(%rsp)
        movaps %xmm4,nb132_fiyO(%rsp)
        movaps %xmm4,nb132_fizO(%rsp)
        movaps %xmm4,nb132_fixH1(%rsp)
        movaps %xmm4,nb132_fiyH1(%rsp)
        movaps %xmm4,nb132_fizH1(%rsp)
        movaps %xmm4,nb132_fixH2(%rsp)
        movaps %xmm4,nb132_fiyH2(%rsp)
        movaps %xmm4,nb132_fizH2(%rsp)

        movq  nb132_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb132_pos(%rbp),%rsi
        movq  nb132_faction(%rbp),%rdi
        movq  nb132_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb132_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb132_ninner(%rsp),%ecx
        movl  %ecx,nb132_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb132_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel132_x86_64_sse.nb132_unroll_loop
        jmp   _nb_kernel132_x86_64_sse.nb132_single_check
_nb_kernel132_x86_64_sse.nb132_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb132_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb132_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb132_pos(%rbp),%rsi        ## base of pos[] 

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

    subps nb132_ixO(%rsp),%xmm0
    subps nb132_iyO(%rsp),%xmm1
    subps nb132_izO(%rsp),%xmm2
    subps nb132_ixH1(%rsp),%xmm3
    subps nb132_iyH1(%rsp),%xmm4
    subps nb132_izH1(%rsp),%xmm5
    subps nb132_ixH2(%rsp),%xmm6
    subps nb132_iyH2(%rsp),%xmm7
    subps nb132_izH2(%rsp),%xmm8

        movaps %xmm0,nb132_dxOO(%rsp)
        movaps %xmm1,nb132_dyOO(%rsp)
        movaps %xmm2,nb132_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH1O(%rsp)
        movaps %xmm4,nb132_dyH1O(%rsp)
        movaps %xmm5,nb132_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb132_dxH2O(%rsp)
        movaps %xmm7,nb132_dyH2O(%rsp)
        movaps %xmm8,nb132_dzH2O(%rsp)
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

        movaps  nb132_three(%rsp),%xmm9
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

        movaps  nb132_half(%rsp),%xmm2
        mulps   %xmm2,%xmm9 ## rinvOO 
        mulps   %xmm2,%xmm10 ## rinvH1O
    mulps   %xmm2,%xmm11 ## rinvH2O

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movaps %xmm0,nb132_rsqOO(%rsp)
    movaps %xmm3,nb132_rsqOH1(%rsp)
    movaps %xmm6,nb132_rsqOH2(%rsp)
    movaps %xmm9,nb132_rinvOO(%rsp)
    movaps %xmm10,nb132_rinvOH1(%rsp)
    movaps %xmm11,nb132_rinvOH2(%rsp)

    ## table LJ interaction
    mulps  %xmm9,%xmm0
    mulps  nb132_tsc(%rsp),%xmm0   ## rtab

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
    movq nb132_VFtab(%rbp),%rsi

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
    movaps nb132_c6(%rsp),%xmm12
    movaps nb132_c12(%rsp),%xmm13
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulps  %xmm13,%xmm9 ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb132_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb132_Vvdwtot(%rsp)

    mulps  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulps  %xmm13,%xmm11  ## FF*c12  = fnb12
    addps  %xmm11,%xmm7
    mulps  nb132_tsc(%rsp),%xmm7

    movaps nb132_rinvOO(%rsp),%xmm9
    movaps nb132_rinvOH1(%rsp),%xmm10
    movaps nb132_rinvOH2(%rsp),%xmm11
    movaps %xmm9,%xmm0  ## rinv
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2

    mulps  %xmm10,%xmm10 ## rinvsq
    mulps  %xmm11,%xmm11
    mulps  nb132_qqOO(%rsp),%xmm0
    mulps  nb132_qqOH(%rsp),%xmm1
    mulps  nb132_qqOH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    subps  %xmm7,%xmm9
    mulps  nb132_rinvOO(%rsp),%xmm9

    addps nb132_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb132_vctot(%rsp)

        ## move j O forces to local temp variables 
        movq nb132_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm0 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm2    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm2    ## jzOa  -  jzOb -
    movhps 8(%rdi,%rdx,4),%xmm3    ## jzOc  -  jzOd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm0: jxOa jyOa jxOb jyOb 
    ## xmm1: jxOc jyOc jxOd jyOd
    ## xmm2: jzOa jzOb jzOc jzOd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb132_dxOO(%rsp),%xmm7
        mulps nb132_dyOO(%rsp),%xmm8
        mulps nb132_dzOO(%rsp),%xmm9
        mulps nb132_dxH1O(%rsp),%xmm10
        mulps nb132_dyH1O(%rsp),%xmm11
        mulps nb132_dzH1O(%rsp),%xmm12
        mulps nb132_dxH2O(%rsp),%xmm13
        mulps nb132_dyH2O(%rsp),%xmm14
        mulps nb132_dzH2O(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb132_fixO(%rsp),%xmm7
    addps nb132_fiyO(%rsp),%xmm8
    addps nb132_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb132_fixH1(%rsp),%xmm10
    addps nb132_fiyH1(%rsp),%xmm11
    addps nb132_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb132_fixH2(%rsp),%xmm13
    addps nb132_fiyH2(%rsp),%xmm14
    addps nb132_fizH2(%rsp),%xmm15

    movaps %xmm7,nb132_fixO(%rsp)
    movaps %xmm8,nb132_fiyO(%rsp)
    movaps %xmm9,nb132_fizO(%rsp)
    movaps %xmm10,nb132_fixH1(%rsp)
    movaps %xmm11,nb132_fiyH1(%rsp)
    movaps %xmm12,nb132_fizH1(%rsp)
    movaps %xmm13,nb132_fixH2(%rsp)
    movaps %xmm14,nb132_fiyH2(%rsp)
    movaps %xmm15,nb132_fizH2(%rsp)

    ## xmm0 = fOx
    ## xmm1 = fOy
    ## xmm2 = fOz
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fOzc fOzd

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


        ## move j H1 coordinates to local temp variables 
        movq nb132_pos(%rbp),%rsi        ## base of pos[] 
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

    subps nb132_ixO(%rsp),%xmm0
    subps nb132_iyO(%rsp),%xmm1
    subps nb132_izO(%rsp),%xmm2
    subps nb132_ixH1(%rsp),%xmm3
    subps nb132_iyH1(%rsp),%xmm4
    subps nb132_izH1(%rsp),%xmm5
    subps nb132_ixH2(%rsp),%xmm6
    subps nb132_iyH2(%rsp),%xmm7
    subps nb132_izH2(%rsp),%xmm8

        movaps %xmm0,nb132_dxOH1(%rsp)
        movaps %xmm1,nb132_dyOH1(%rsp)
        movaps %xmm2,nb132_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH1H1(%rsp)
        movaps %xmm4,nb132_dyH1H1(%rsp)
        movaps %xmm5,nb132_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb132_dxH2H1(%rsp)
        movaps %xmm7,nb132_dyH2H1(%rsp)
        movaps %xmm8,nb132_dzH2H1(%rsp)
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

        movaps  nb132_three(%rsp),%xmm9
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

        movaps  nb132_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOH1
        mulps   %xmm0,%xmm10 ## rinvH1H1
    mulps   %xmm0,%xmm11 ## rinvH2H1

        ## H1 interactions 
    movaps %xmm9,%xmm0  ## rinv
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9   ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb132_qqOH(%rsp),%xmm0      ## vcoul
    mulps  nb132_qqHH(%rsp),%xmm1
    mulps  nb132_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb132_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb132_vctot(%rsp)

        ## move j H1 forces to local temp variables 
        movq nb132_faction(%rbp),%rdi        ## base of faction[] 
    movlps 12(%rdi,%rax,4),%xmm0    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm1    ## jxH1c jyH1c  -   -
    movhps 12(%rdi,%rbx,4),%xmm0    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm1    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rdi,%rax,4),%xmm2    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm3    ## jzH1c  -  -  -
    movhps 20(%rdi,%rbx,4),%xmm2    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm3    ## jzH1c  -  jzH1d -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm0: jxH1a jyH1a jxH1b jyH1b 
    ## xmm1: jxH1c jyH1c jxH1d jyH1d
    ## xmm2: jzH1a jzH1b jzH1c jzH1d

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb132_dxOH1(%rsp),%xmm7
        mulps nb132_dyOH1(%rsp),%xmm8
        mulps nb132_dzOH1(%rsp),%xmm9
        mulps nb132_dxH1H1(%rsp),%xmm10
        mulps nb132_dyH1H1(%rsp),%xmm11
        mulps nb132_dzH1H1(%rsp),%xmm12
        mulps nb132_dxH2H1(%rsp),%xmm13
        mulps nb132_dyH2H1(%rsp),%xmm14
        mulps nb132_dzH2H1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb132_fixO(%rsp),%xmm7
    addps nb132_fiyO(%rsp),%xmm8
    addps nb132_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb132_fixH1(%rsp),%xmm10
    addps nb132_fiyH1(%rsp),%xmm11
    addps nb132_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb132_fixH2(%rsp),%xmm13
    addps nb132_fiyH2(%rsp),%xmm14
    addps nb132_fizH2(%rsp),%xmm15

    movaps %xmm7,nb132_fixO(%rsp)
    movaps %xmm8,nb132_fiyO(%rsp)
    movaps %xmm9,nb132_fizO(%rsp)
    movaps %xmm10,nb132_fixH1(%rsp)
    movaps %xmm11,nb132_fiyH1(%rsp)
    movaps %xmm12,nb132_fizH1(%rsp)
    movaps %xmm13,nb132_fixH2(%rsp)
    movaps %xmm14,nb132_fiyH2(%rsp)
    movaps %xmm15,nb132_fizH2(%rsp)

    ## xmm0 = fH1x
    ## xmm1 = fH1y
    ## xmm2 = fH1z
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fH1zc fH1zd

    movlps %xmm0,12(%rdi,%rax,4)
    movhps %xmm0,12(%rdi,%rbx,4)
    movlps %xmm1,12(%rdi,%rcx,4)
    movhps %xmm1,12(%rdi,%rdx,4)
    movss  %xmm2,20(%rdi,%rax,4)
    movss  %xmm3,20(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,20(%rdi,%rbx,4)
    movss  %xmm3,20(%rdi,%rdx,4)


        ## move j H2 coordinates to local temp variables 
        movq nb132_pos(%rbp),%rsi        ## base of pos[] 
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

    subps nb132_ixO(%rsp),%xmm0
    subps nb132_iyO(%rsp),%xmm1
    subps nb132_izO(%rsp),%xmm2
    subps nb132_ixH1(%rsp),%xmm3
    subps nb132_iyH1(%rsp),%xmm4
    subps nb132_izH1(%rsp),%xmm5
    subps nb132_ixH2(%rsp),%xmm6
    subps nb132_iyH2(%rsp),%xmm7
    subps nb132_izH2(%rsp),%xmm8

        movaps %xmm0,nb132_dxOH2(%rsp)
        movaps %xmm1,nb132_dyOH2(%rsp)
        movaps %xmm2,nb132_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH1H2(%rsp)
        movaps %xmm4,nb132_dyH1H2(%rsp)
        movaps %xmm5,nb132_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb132_dxH2H2(%rsp)
        movaps %xmm7,nb132_dyH2H2(%rsp)
        movaps %xmm8,nb132_dzH2H2(%rsp)
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

        movaps  nb132_three(%rsp),%xmm9
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

        movaps  nb132_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOH2 
        mulps   %xmm0,%xmm10 ## rinvH1H2
    mulps   %xmm0,%xmm11 ## rinvH2H2

        ## H2 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb132_qqOH(%rsp),%xmm0
    mulps  nb132_qqHH(%rsp),%xmm1
    mulps  nb132_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb132_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb132_vctot(%rsp)

        ## move j H2 forces to local temp variables 
        movq nb132_faction(%rbp),%rdi        ## base of faction[] 
    movlps 24(%rdi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movss  32(%rdi,%rbx,4),%xmm7    ## jzH2b  -  -  -
    movss  32(%rdi,%rdx,4),%xmm8    ## jzH2d  -  -  -
    movlhps %xmm7,%xmm2 ## jzH2a  -  jzH2b  -
    movlhps %xmm8,%xmm3 ## jzH2c  -  jzH2d -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm0: jxH2a jyH2a jxH2b jyH2b 
    ## xmm1: jxH2c jyH2c jxH2d jyH2d
    ## xmm2: jzH2a jzH2b jzH2c jzH2d

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb132_dxOH2(%rsp),%xmm7
        mulps nb132_dyOH2(%rsp),%xmm8
        mulps nb132_dzOH2(%rsp),%xmm9
        mulps nb132_dxH1H2(%rsp),%xmm10
        mulps nb132_dyH1H2(%rsp),%xmm11
        mulps nb132_dzH1H2(%rsp),%xmm12
        mulps nb132_dxH2H2(%rsp),%xmm13
        mulps nb132_dyH2H2(%rsp),%xmm14
        mulps nb132_dzH2H2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb132_fixO(%rsp),%xmm7
    addps nb132_fiyO(%rsp),%xmm8
    addps nb132_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb132_fixH1(%rsp),%xmm10
    addps nb132_fiyH1(%rsp),%xmm11
    addps nb132_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb132_fixH2(%rsp),%xmm13
    addps nb132_fiyH2(%rsp),%xmm14
    addps nb132_fizH2(%rsp),%xmm15

    movaps %xmm7,nb132_fixO(%rsp)
    movaps %xmm8,nb132_fiyO(%rsp)
    movaps %xmm9,nb132_fizO(%rsp)
    movaps %xmm10,nb132_fixH1(%rsp)
    movaps %xmm11,nb132_fiyH1(%rsp)
    movaps %xmm12,nb132_fizH1(%rsp)
    movaps %xmm13,nb132_fixH2(%rsp)
    movaps %xmm14,nb132_fiyH2(%rsp)
    movaps %xmm15,nb132_fizH2(%rsp)

    ## xmm0 = fH2x
    ## xmm1 = fH2y
    ## xmm2 = fH2z
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fH2zc fH2zd

    movlps %xmm0,24(%rdi,%rax,4)
    movhps %xmm0,24(%rdi,%rbx,4)
    movlps %xmm1,24(%rdi,%rcx,4)
    movhps %xmm1,24(%rdi,%rdx,4)
    movss  %xmm2,32(%rdi,%rax,4)
    movss  %xmm3,32(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,32(%rdi,%rbx,4)
    movss  %xmm3,32(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb132_innerk(%rsp)
        jl    _nb_kernel132_x86_64_sse.nb132_single_check
        jmp   _nb_kernel132_x86_64_sse.nb132_unroll_loop
_nb_kernel132_x86_64_sse.nb132_single_check: 
        addl $4,nb132_innerk(%rsp)
        jnz   _nb_kernel132_x86_64_sse.nb132_single_loop
        jmp   _nb_kernel132_x86_64_sse.nb132_updateouterdata
_nb_kernel132_x86_64_sse.nb132_single_loop: 
        movq  nb132_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb132_innerjjnr(%rsp)

        movq nb132_pos(%rbp),%rsi
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
        movaps %xmm0,nb132_jxO(%rsp)
        movaps %xmm1,nb132_jyO(%rsp)
        movaps %xmm2,nb132_jzO(%rsp)
        subps  nb132_ixO(%rsp),%xmm0
        subps  nb132_iyO(%rsp),%xmm1
        subps  nb132_izO(%rsp),%xmm2
        movaps %xmm0,nb132_dxOO(%rsp)
        movaps %xmm1,nb132_dyOO(%rsp)
        movaps %xmm2,nb132_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        movaps  %xmm0,nb132_rsqOO(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132_half(%rsp),%xmm3   ## rinv iO - j water in xmm3
        movaps  %xmm3,nb132_rinvOO(%rsp)

        movaps  %xmm3,%xmm0 ## rinv

        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb132_qqOO(%rsp),%xmm4
        movhps  nb132_qqOH(%rsp),%xmm4

        mulps  %xmm4,%xmm3  ## vcoul
        movaps %xmm3,%xmm6

        mulps %xmm0,%xmm3
        movaps %xmm3,nb132_fstmp(%rsp)   ## save it
        addps  nb132_vctot(%rsp),%xmm6
    movaps %xmm6,nb132_vctot(%rsp)

        movss %xmm0,%xmm1  ## rinv
        mulss  nb132_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulss  nb132_tsc(%rsp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movq nb132_VFtab(%rbp),%rsi
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
    addss  %xmm7,%xmm7      ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb132_c6(%rsp),%xmm4
    mulss  %xmm4,%xmm7   ## fijD 
    mulss  %xmm4,%xmm5   ## Vvdw6 
        movss  nb132_fstmp(%rsp),%xmm3
        mulps  nb132_tsc(%rsp),%xmm7
        subss  %xmm7,%xmm3
        movss  %xmm3,nb132_fstmp(%rsp)

    addss  nb132_Vvdwtot(%rsp),%xmm5
    movss %xmm5,nb132_Vvdwtot(%rsp)

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
    addss  %xmm7,%xmm7      ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb132_c12(%rsp),%xmm4
    mulss  %xmm4,%xmm7 ## fijR 
    mulss  %xmm4,%xmm5 ## Vvdw12 
        movaps nb132_fstmp(%rsp),%xmm3
        mulss  nb132_tsc(%rsp),%xmm7
        subss  %xmm7,%xmm3

    addss  nb132_Vvdwtot(%rsp),%xmm5
    movss %xmm5,nb132_Vvdwtot(%rsp)

        mulps  %xmm3,%xmm0
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2

        mulps   nb132_dxOO(%rsp),%xmm0
        mulps   nb132_dyOO(%rsp),%xmm1
        mulps   nb132_dzOO(%rsp),%xmm2

        ## initial update for j forces 
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        addps   nb132_fixO(%rsp),%xmm0
        addps   nb132_fiyO(%rsp),%xmm1
        addps   nb132_fizO(%rsp),%xmm2
        movaps  %xmm0,nb132_fixO(%rsp)
        movaps  %xmm1,nb132_fiyO(%rsp)
        movaps  %xmm2,nb132_fizO(%rsp)
        movaps  %xmm3,nb132_fjxO(%rsp)
        movaps  %xmm4,nb132_fjyO(%rsp)
        movaps  %xmm5,nb132_fjzO(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb132_jxO(%rsp),%xmm0
    movaps  nb132_jyO(%rsp),%xmm1
    movaps  nb132_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  nb132_ixH1(%rsp),%xmm0
        subps  nb132_iyH1(%rsp),%xmm1
        subps  nb132_izH1(%rsp),%xmm2
        subps  nb132_ixH2(%rsp),%xmm3
        subps  nb132_iyH2(%rsp),%xmm4
        subps  nb132_izH2(%rsp),%xmm5
    movaps %xmm0,nb132_dxH1O(%rsp)
        movaps %xmm1,nb132_dyH1O(%rsp)
        movaps %xmm2,nb132_dzH1O(%rsp)
        movaps %xmm3,nb132_dxH2O(%rsp)
        movaps %xmm4,nb132_dyH2O(%rsp)
        movaps %xmm5,nb132_dzH2O(%rsp)
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

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb132_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb132_qqOH(%rsp),%xmm6
        movhps  nb132_qqHH(%rsp),%xmm6

        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm1
        mulps   %xmm6,%xmm3 ## vcoul
        mulps   %xmm6,%xmm7 ## vcoul
        mulps   %xmm0,%xmm0
        mulps   %xmm1,%xmm1
        mulps   %xmm3,%xmm0
        addps   %xmm7,%xmm3
        mulps   %xmm1,%xmm7
        addps   nb132_vctot(%rsp),%xmm3
        movaps  %xmm3,nb132_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb132_dxH1O(%rsp),%xmm0
        mulps   nb132_dyH1O(%rsp),%xmm1
        mulps   nb132_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        addps   nb132_fixH1(%rsp),%xmm0
        addps   nb132_fiyH1(%rsp),%xmm1
        addps   nb132_fizH1(%rsp),%xmm2
        addps   nb132_fjxO(%rsp),%xmm3
        addps   nb132_fjyO(%rsp),%xmm4
        addps   nb132_fjzO(%rsp),%xmm5
        movaps  %xmm0,nb132_fixH1(%rsp)
        movaps  %xmm1,nb132_fiyH1(%rsp)
        movaps  %xmm2,nb132_fizH1(%rsp)
        movaps  %xmm3,nb132_fjxO(%rsp)
        movaps  %xmm4,nb132_fjyO(%rsp)
        movaps  %xmm5,nb132_fjzO(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb132_dxH2O(%rsp),%xmm0
        mulps   nb132_dyH2O(%rsp),%xmm1
        mulps   nb132_dzH2O(%rsp),%xmm2

    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        addps   nb132_fixH2(%rsp),%xmm0
        addps   nb132_fiyH2(%rsp),%xmm1
        addps   nb132_fizH2(%rsp),%xmm2
        addps   nb132_fjxO(%rsp),%xmm3
        addps   nb132_fjyO(%rsp),%xmm4
        addps   nb132_fjzO(%rsp),%xmm5
        movaps  %xmm0,nb132_fixH2(%rsp)
        movaps  %xmm1,nb132_fiyH2(%rsp)
        movaps  %xmm2,nb132_fizH2(%rsp)
    ## fjx/fjy/fjz in xmm3-xmm5

        movq    nb132_faction(%rbp),%rsi

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0        ## 
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1

        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
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

        decl nb132_innerk(%rsp)
        jz    _nb_kernel132_x86_64_sse.nb132_updateouterdata
        jmp   _nb_kernel132_x86_64_sse.nb132_single_loop
_nb_kernel132_x86_64_sse.nb132_updateouterdata: 
        movl  nb132_ii3(%rsp),%ecx
        movq  nb132_faction(%rbp),%rdi
        movq  nb132_fshift(%rbp),%rsi
        movl  nb132_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb132_fixO(%rsp),%xmm0
        movaps nb132_fiyO(%rsp),%xmm1
        movaps nb132_fizO(%rsp),%xmm2

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
        movaps nb132_fixH1(%rsp),%xmm0
        movaps nb132_fiyH1(%rsp),%xmm1
        movaps nb132_fizH1(%rsp),%xmm2

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
        movaps nb132_fixH2(%rsp),%xmm0
        movaps nb132_fiyH2(%rsp),%xmm1
        movaps nb132_fizH2(%rsp),%xmm2

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
        movl nb132_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb132_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb132_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb132_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb132_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb132_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb132_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel132_x86_64_sse.nb132_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb132_n(%rsp)
        jmp _nb_kernel132_x86_64_sse.nb132_outer
_nb_kernel132_x86_64_sse.nb132_outerend: 
        ## check if more outer neighborlists remain
        movl  nb132_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel132_x86_64_sse.nb132_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel132_x86_64_sse.nb132_threadloop
_nb_kernel132_x86_64_sse.nb132_end: 
        movl nb132_nouter(%rsp),%eax
        movl nb132_ninner(%rsp),%ebx
        movq nb132_outeriter(%rbp),%rcx
        movq nb132_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1624,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel132nf_x86_64_sse
.globl _nb_kernel132nf_x86_64_sse
nb_kernel132nf_x86_64_sse:      
_nb_kernel132nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb132nf_fshift, 16
.set nb132nf_gid, 24
.set nb132nf_pos, 32
.set nb132nf_faction, 40
.set nb132nf_charge, 48
.set nb132nf_p_facel, 56
.set nb132nf_argkrf, 64
.set nb132nf_argcrf, 72
.set nb132nf_Vc, 80
.set nb132nf_type, 88
.set nb132nf_p_ntype, 96
.set nb132nf_vdwparam, 104
.set nb132nf_Vvdw, 112
.set nb132nf_p_tabscale, 120
.set nb132nf_VFtab, 128
.set nb132nf_invsqrta, 136
.set nb132nf_dvda, 144
.set nb132nf_p_gbtabscale, 152
.set nb132nf_GBtab, 160
.set nb132nf_p_nthreads, 168
.set nb132nf_count, 176
.set nb132nf_mtx, 184
.set nb132nf_outeriter, 192
.set nb132nf_inneriter, 200
.set nb132nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb132nf_ixO, 0
.set nb132nf_iyO, 16
.set nb132nf_izO, 32
.set nb132nf_ixH1, 48
.set nb132nf_iyH1, 64
.set nb132nf_izH1, 80
.set nb132nf_ixH2, 96
.set nb132nf_iyH2, 112
.set nb132nf_izH2, 128
.set nb132nf_jxO, 144
.set nb132nf_jyO, 160
.set nb132nf_jzO, 176
.set nb132nf_jxH1, 192
.set nb132nf_jyH1, 208
.set nb132nf_jzH1, 224
.set nb132nf_jxH2, 240
.set nb132nf_jyH2, 256
.set nb132nf_jzH2, 272
.set nb132nf_qqOO, 288
.set nb132nf_qqOH, 304
.set nb132nf_qqHH, 320
.set nb132nf_c6, 336
.set nb132nf_c12, 352
.set nb132nf_vctot, 368
.set nb132nf_Vvdwtot, 384
.set nb132nf_half, 400
.set nb132nf_three, 416
.set nb132nf_rsqOO, 432
.set nb132nf_rsqOH1, 448
.set nb132nf_rsqOH2, 464
.set nb132nf_rsqH1O, 480
.set nb132nf_rsqH1H1, 496
.set nb132nf_rsqH1H2, 512
.set nb132nf_rsqH2O, 528
.set nb132nf_rsqH2H1, 544
.set nb132nf_rsqH2H2, 560
.set nb132nf_rinvOO, 576
.set nb132nf_rinvOH1, 592
.set nb132nf_rinvOH2, 608
.set nb132nf_rinvH1O, 624
.set nb132nf_rinvH1H1, 640
.set nb132nf_rinvH1H2, 656
.set nb132nf_rinvH2O, 672
.set nb132nf_rinvH2H1, 688
.set nb132nf_rinvH2H2, 704
.set nb132nf_krf, 720
.set nb132nf_crf, 736
.set nb132nf_tsc, 752
.set nb132nf_nri, 768
.set nb132nf_iinr, 776
.set nb132nf_jindex, 784
.set nb132nf_jjnr, 792
.set nb132nf_shift, 800
.set nb132nf_shiftvec, 808
.set nb132nf_facel, 816
.set nb132nf_innerjjnr, 824
.set nb132nf_is3, 832
.set nb132nf_ii3, 836
.set nb132nf_innerk, 840
.set nb132nf_n, 844
.set nb132nf_nn1, 848
.set nb132nf_nouter, 852
.set nb132nf_ninner, 856
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $872,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb132nf_nouter(%rsp)
        movl %eax,nb132nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb132nf_nri(%rsp)
        movq %rsi,nb132nf_iinr(%rsp)
        movq %rdx,nb132nf_jindex(%rsp)
        movq %rcx,nb132nf_jjnr(%rsp)
        movq %r8,nb132nf_shift(%rsp)
        movq %r9,nb132nf_shiftvec(%rsp)
        movq nb132nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb132nf_facel(%rsp)

        movq nb132nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb132nf_tsc(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb132nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb132nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb132nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb132nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132nf_qqOO(%rsp)
        movaps %xmm4,nb132nf_qqOH(%rsp)
        movaps %xmm5,nb132nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb132nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb132nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addq  %rcx,%rdx
        movq  nb132nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb132nf_c6(%rsp)
        movaps %xmm1,nb132nf_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb132nf_half(%rsp)
        movss nb132nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb132nf_half(%rsp)
        movaps %xmm3,nb132nf_three(%rsp)

_nb_kernel132nf_x86_64_sse.nb132nf_threadloop: 
        movq  nb132nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel132nf_x86_64_sse.nb132nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel132nf_x86_64_sse.nb132nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb132nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb132nf_n(%rsp)
        movl %ebx,nb132nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel132nf_x86_64_sse.nb132nf_outerstart
        jmp _nb_kernel132nf_x86_64_sse.nb132nf_end

_nb_kernel132nf_x86_64_sse.nb132nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb132nf_nouter(%rsp),%ebx
        movl %ebx,nb132nf_nouter(%rsp)

_nb_kernel132nf_x86_64_sse.nb132nf_outer: 
        movq  nb132nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb132nf_is3(%rsp)            ## store is3 

        movq  nb132nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb132nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb132nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb132nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132nf_ixO(%rsp)
        movaps %xmm4,nb132nf_iyO(%rsp)
        movaps %xmm5,nb132nf_izO(%rsp)

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
        movaps %xmm0,nb132nf_ixH1(%rsp)
        movaps %xmm1,nb132nf_iyH1(%rsp)
        movaps %xmm2,nb132nf_izH1(%rsp)
        movaps %xmm3,nb132nf_ixH2(%rsp)
        movaps %xmm4,nb132nf_iyH2(%rsp)
        movaps %xmm5,nb132nf_izH2(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb132nf_vctot(%rsp)
        movaps %xmm4,nb132nf_Vvdwtot(%rsp)

        movq  nb132nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb132nf_pos(%rbp),%rsi
        movq  nb132nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb132nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb132nf_ninner(%rsp),%ecx
        movl  %ecx,nb132nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb132nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel132nf_x86_64_sse.nb132nf_unroll_loop
        jmp   _nb_kernel132nf_x86_64_sse.nb132nf_single_check
_nb_kernel132nf_x86_64_sse.nb132nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb132nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb132nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb132nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps %xmm0,nb132nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb132nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb132nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb132nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb132nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb132nf_jyH2(%rsp)

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

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb132nf_jzO(%rsp)
        movaps %xmm1,nb132nf_jzH1(%rsp)
        movaps %xmm2,nb132nf_jzH2(%rsp)

        movaps nb132nf_ixO(%rsp),%xmm0
        movaps nb132nf_iyO(%rsp),%xmm1
        movaps nb132nf_izO(%rsp),%xmm2
        movaps nb132nf_ixO(%rsp),%xmm3
        movaps nb132nf_iyO(%rsp),%xmm4
        movaps nb132nf_izO(%rsp),%xmm5
        subps  nb132nf_jxO(%rsp),%xmm0
        subps  nb132nf_jyO(%rsp),%xmm1
        subps  nb132nf_jzO(%rsp),%xmm2
        subps  nb132nf_jxH1(%rsp),%xmm3
        subps  nb132nf_jyH1(%rsp),%xmm4
        subps  nb132nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb132nf_rsqOO(%rsp)
        movaps %xmm3,nb132nf_rsqOH1(%rsp)

        movaps nb132nf_ixO(%rsp),%xmm0
        movaps nb132nf_iyO(%rsp),%xmm1
        movaps nb132nf_izO(%rsp),%xmm2
        movaps nb132nf_ixH1(%rsp),%xmm3
        movaps nb132nf_iyH1(%rsp),%xmm4
        movaps nb132nf_izH1(%rsp),%xmm5
        subps  nb132nf_jxH2(%rsp),%xmm0
        subps  nb132nf_jyH2(%rsp),%xmm1
        subps  nb132nf_jzH2(%rsp),%xmm2
        subps  nb132nf_jxO(%rsp),%xmm3
        subps  nb132nf_jyO(%rsp),%xmm4
        subps  nb132nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb132nf_rsqOH2(%rsp)
        movaps %xmm3,nb132nf_rsqH1O(%rsp)

        movaps nb132nf_ixH1(%rsp),%xmm0
        movaps nb132nf_iyH1(%rsp),%xmm1
        movaps nb132nf_izH1(%rsp),%xmm2
        movaps nb132nf_ixH1(%rsp),%xmm3
        movaps nb132nf_iyH1(%rsp),%xmm4
        movaps nb132nf_izH1(%rsp),%xmm5
        subps  nb132nf_jxH1(%rsp),%xmm0
        subps  nb132nf_jyH1(%rsp),%xmm1
        subps  nb132nf_jzH1(%rsp),%xmm2
        subps  nb132nf_jxH2(%rsp),%xmm3
        subps  nb132nf_jyH2(%rsp),%xmm4
        subps  nb132nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb132nf_rsqH1H1(%rsp)
        movaps %xmm3,nb132nf_rsqH1H2(%rsp)

        movaps nb132nf_ixH2(%rsp),%xmm0
        movaps nb132nf_iyH2(%rsp),%xmm1
        movaps nb132nf_izH2(%rsp),%xmm2
        movaps nb132nf_ixH2(%rsp),%xmm3
        movaps nb132nf_iyH2(%rsp),%xmm4
        movaps nb132nf_izH2(%rsp),%xmm5
        subps  nb132nf_jxO(%rsp),%xmm0
        subps  nb132nf_jyO(%rsp),%xmm1
        subps  nb132nf_jzO(%rsp),%xmm2
        subps  nb132nf_jxH1(%rsp),%xmm3
        subps  nb132nf_jyH1(%rsp),%xmm4
        subps  nb132nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb132nf_rsqH2O(%rsp)
        movaps %xmm4,nb132nf_rsqH2H1(%rsp)

        movaps nb132nf_ixH2(%rsp),%xmm0
        movaps nb132nf_iyH2(%rsp),%xmm1
        movaps nb132nf_izH2(%rsp),%xmm2
        subps  nb132nf_jxH2(%rsp),%xmm0
        subps  nb132nf_jyH2(%rsp),%xmm1
        subps  nb132nf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb132nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb132nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb132nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb132nf_rinvH2H1(%rsp)

        rsqrtps nb132nf_rsqOO(%rsp),%xmm1
        rsqrtps nb132nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqOO(%rsp),%xmm1
        mulps   nb132nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%rsp),%xmm3
        mulps   nb132nf_half(%rsp),%xmm7
        movaps  %xmm3,nb132nf_rinvOO(%rsp)
        movaps  %xmm7,nb132nf_rinvOH1(%rsp)

        rsqrtps nb132nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb132nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqOH2(%rsp),%xmm1
        mulps   nb132nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%rsp),%xmm3
        mulps   nb132nf_half(%rsp),%xmm7
        movaps  %xmm3,nb132nf_rinvOH2(%rsp)
        movaps  %xmm7,nb132nf_rinvH1O(%rsp)

        rsqrtps nb132nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb132nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqH1H1(%rsp),%xmm1
        mulps   nb132nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%rsp),%xmm3
        mulps   nb132nf_half(%rsp),%xmm7
        movaps  %xmm3,nb132nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb132nf_rinvH1H2(%rsp)

        rsqrtps nb132nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132nf_three(%rsp),%xmm3
        mulps   nb132nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132nf_half(%rsp),%xmm3
        movaps  %xmm3,nb132nf_rinvH2O(%rsp)

        ## start with OO interaction - first the table LJ part
        movaps nb132nf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb132nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulps  nb132nf_tsc(%rsp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6
    pslld $3,%mm7

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movq nb132nf_VFtab(%rbp),%rsi

    ## dispersion 
    movlps (%rsi,%rax,4),%xmm5
    movlps (%rsi,%rcx,4),%xmm7
    movhps (%rsi,%rbx,4),%xmm5
    movhps (%rsi,%rdx,4),%xmm7 ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%rsi,%rax,4),%xmm7
    movlps 8(%rsi,%rcx,4),%xmm3
    movhps 8(%rsi,%rbx,4),%xmm7
    movhps 8(%rsi,%rdx,4),%xmm3    ## other half of table  
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

    movaps nb132nf_c6(%rsp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 

    addps  nb132nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb132nf_Vvdwtot(%rsp)

    ## repulsion 
    movlps 16(%rsi,%rax,4),%xmm5
    movlps 16(%rsi,%rcx,4),%xmm7
    movhps 16(%rsi,%rbx,4),%xmm5
    movhps 16(%rsi,%rdx,4),%xmm7    ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%rsi,%rax,4),%xmm7
    movlps 24(%rsi,%rcx,4),%xmm3
    movhps 24(%rsi,%rbx,4),%xmm7
    movhps 24(%rsi,%rdx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb132nf_c12(%rsp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 

    addps  nb132nf_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb132nf_Vvdwtot(%rsp)

        ## Coulomb interactions 
        movaps nb132nf_rinvOH1(%rsp),%xmm0
        movaps nb132nf_rinvH1H1(%rsp),%xmm1
        movaps nb132nf_rinvOO(%rsp),%xmm2
        addps  nb132nf_rinvOH2(%rsp),%xmm0
        addps  nb132nf_rinvH1H2(%rsp),%xmm1
        addps  nb132nf_rinvH1O(%rsp),%xmm0
        addps  nb132nf_rinvH2H1(%rsp),%xmm1
        addps  nb132nf_rinvH2O(%rsp),%xmm0
        addps  nb132nf_rinvH2H2(%rsp),%xmm1

        mulps  nb132nf_qqOH(%rsp),%xmm0
        mulps  nb132nf_qqHH(%rsp),%xmm1
        mulps  nb132nf_qqOO(%rsp),%xmm2

        addps  %xmm1,%xmm0
        addps  nb132nf_vctot(%rsp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb132nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb132nf_innerk(%rsp)
        jl    _nb_kernel132nf_x86_64_sse.nb132nf_single_check
        jmp   _nb_kernel132nf_x86_64_sse.nb132nf_unroll_loop
_nb_kernel132nf_x86_64_sse.nb132nf_single_check: 
        addl $4,nb132nf_innerk(%rsp)
        jnz   _nb_kernel132nf_x86_64_sse.nb132nf_single_loop
        jmp   _nb_kernel132nf_x86_64_sse.nb132nf_updateouterdata
_nb_kernel132nf_x86_64_sse.nb132nf_single_loop: 
        movq  nb132nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb132nf_innerjjnr(%rsp)

        movq nb132nf_pos(%rbp),%rsi
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

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb132nf_ixO(%rsp),%xmm0
        movaps  nb132nf_iyO(%rsp),%xmm1
        movaps  nb132nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb132nf_jxO(%rsp)
        movaps %xmm4,nb132nf_jyO(%rsp)
        movaps %xmm5,nb132nf_jzO(%rsp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        movaps %xmm0,nb132nf_rsqOO(%rsp)

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132nf_half(%rsp),%xmm3   ## rinv iO - j water in xmm3
        movaps  %xmm3,nb132nf_rinvOO(%rsp)


        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb132nf_qqOO(%rsp),%xmm4
        movhps  nb132nf_qqOH(%rsp),%xmm4

        mulps %xmm4,%xmm3       ## vcoul  
        addps  nb132nf_vctot(%rsp),%xmm3
    movaps %xmm3,nb132nf_vctot(%rsp)


        movaps nb132nf_rinvOO(%rsp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb132nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulss  nb132nf_tsc(%rsp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movq nb132nf_VFtab(%rbp),%rsi
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

    movss nb132nf_c6(%rsp),%xmm4
    mulss  %xmm4,%xmm5   ## Vvdw6 

    addss  nb132nf_Vvdwtot(%rsp),%xmm5
    movss %xmm5,nb132nf_Vvdwtot(%rsp)

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

    movss nb132nf_c12(%rsp),%xmm4
    mulss  %xmm4,%xmm5 ## Vvdw12 
    addss  nb132nf_Vvdwtot(%rsp),%xmm5
    movss %xmm5,nb132nf_Vvdwtot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb132nf_ixH1(%rsp),%xmm0
        movaps  nb132nf_iyH1(%rsp),%xmm1
        movaps  nb132nf_izH1(%rsp),%xmm2
        movaps  nb132nf_ixH2(%rsp),%xmm3
        movaps  nb132nf_iyH2(%rsp),%xmm4
        movaps  nb132nf_izH2(%rsp),%xmm5
        subps   nb132nf_jxO(%rsp),%xmm0
        subps   nb132nf_jyO(%rsp),%xmm1
        subps   nb132nf_jzO(%rsp),%xmm2
        subps   nb132nf_jxO(%rsp),%xmm3
        subps   nb132nf_jyO(%rsp),%xmm4
        subps   nb132nf_jzO(%rsp),%xmm5
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

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb132nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb132nf_qqOH(%rsp),%xmm6
        movhps  nb132nf_qqHH(%rsp),%xmm6

        addps   %xmm7,%xmm3
        mulps   %xmm6,%xmm3

        addps   nb132nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb132nf_vctot(%rsp)

        decl nb132nf_innerk(%rsp)
        jz    _nb_kernel132nf_x86_64_sse.nb132nf_updateouterdata
        jmp   _nb_kernel132nf_x86_64_sse.nb132nf_single_loop
_nb_kernel132nf_x86_64_sse.nb132nf_updateouterdata: 
        ## get n from stack
        movl nb132nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb132nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb132nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb132nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb132nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb132nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb132nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel132nf_x86_64_sse.nb132nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb132nf_n(%rsp)
        jmp _nb_kernel132nf_x86_64_sse.nb132nf_outer
_nb_kernel132nf_x86_64_sse.nb132nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb132nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel132nf_x86_64_sse.nb132nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel132nf_x86_64_sse.nb132nf_threadloop
_nb_kernel132nf_x86_64_sse.nb132nf_end: 
        movl nb132nf_nouter(%rsp),%eax
        movl nb132nf_ninner(%rsp),%ebx
        movq nb132nf_outeriter(%rbp),%rcx
        movq nb132nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $872,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


