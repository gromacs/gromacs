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






.globl nb_kernel112_x86_64_sse
.globl _nb_kernel112_x86_64_sse
nb_kernel112_x86_64_sse:        
_nb_kernel112_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb112_fshift, 16
.set nb112_gid, 24
.set nb112_pos, 32
.set nb112_faction, 40
.set nb112_charge, 48
.set nb112_p_facel, 56
.set nb112_argkrf, 64
.set nb112_argcrf, 72
.set nb112_Vc, 80
.set nb112_type, 88
.set nb112_p_ntype, 96
.set nb112_vdwparam, 104
.set nb112_Vvdw, 112
.set nb112_p_tabscale, 120
.set nb112_VFtab, 128
.set nb112_invsqrta, 136
.set nb112_dvda, 144
.set nb112_p_gbtabscale, 152
.set nb112_GBtab, 160
.set nb112_p_nthreads, 168
.set nb112_count, 176
.set nb112_mtx, 184
.set nb112_outeriter, 192
.set nb112_inneriter, 200
.set nb112_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb112_ixO, 0
.set nb112_iyO, 16
.set nb112_izO, 32
.set nb112_ixH1, 48
.set nb112_iyH1, 64
.set nb112_izH1, 80
.set nb112_ixH2, 96
.set nb112_iyH2, 112
.set nb112_izH2, 128
.set nb112_jxO, 144
.set nb112_jyO, 160
.set nb112_jzO, 176
.set nb112_jxH1, 192
.set nb112_jyH1, 208
.set nb112_jzH1, 224
.set nb112_jxH2, 240
.set nb112_jyH2, 256
.set nb112_jzH2, 272
.set nb112_dxOO, 288
.set nb112_dyOO, 304
.set nb112_dzOO, 320
.set nb112_dxOH1, 336
.set nb112_dyOH1, 352
.set nb112_dzOH1, 368
.set nb112_dxOH2, 384
.set nb112_dyOH2, 400
.set nb112_dzOH2, 416
.set nb112_dxH1O, 432
.set nb112_dyH1O, 448
.set nb112_dzH1O, 464
.set nb112_dxH1H1, 480
.set nb112_dyH1H1, 496
.set nb112_dzH1H1, 512
.set nb112_dxH1H2, 528
.set nb112_dyH1H2, 544
.set nb112_dzH1H2, 560
.set nb112_dxH2O, 576
.set nb112_dyH2O, 592
.set nb112_dzH2O, 608
.set nb112_dxH2H1, 624
.set nb112_dyH2H1, 640
.set nb112_dzH2H1, 656
.set nb112_dxH2H2, 672
.set nb112_dyH2H2, 688
.set nb112_dzH2H2, 704
.set nb112_qqOO, 720
.set nb112_qqOH, 736
.set nb112_qqHH, 752
.set nb112_c6, 768
.set nb112_c12, 784
.set nb112_six, 800
.set nb112_twelve, 816
.set nb112_vctot, 832
.set nb112_Vvdwtot, 848
.set nb112_fixO, 864
.set nb112_fiyO, 880
.set nb112_fizO, 896
.set nb112_fixH1, 912
.set nb112_fiyH1, 928
.set nb112_fizH1, 944
.set nb112_fixH2, 960
.set nb112_fiyH2, 976
.set nb112_fizH2, 992
.set nb112_fjxO, 1008
.set nb112_fjyO, 1024
.set nb112_fjzO, 1040
.set nb112_fjxH1, 1056
.set nb112_fjyH1, 1072
.set nb112_fjzH1, 1088
.set nb112_fjxH2, 1104
.set nb112_fjyH2, 1120
.set nb112_fjzH2, 1136
.set nb112_half, 1152
.set nb112_three, 1168
.set nb112_rsqOO, 1184
.set nb112_rsqOH1, 1200
.set nb112_rsqOH2, 1216
.set nb112_rsqH1O, 1232
.set nb112_rsqH1H1, 1248
.set nb112_rsqH1H2, 1264
.set nb112_rsqH2O, 1280
.set nb112_rsqH2H1, 1296
.set nb112_rsqH2H2, 1312
.set nb112_rinvOO, 1328
.set nb112_rinvOH1, 1344
.set nb112_rinvOH2, 1360
.set nb112_rinvH1O, 1376
.set nb112_rinvH1H1, 1392
.set nb112_rinvH1H2, 1408
.set nb112_rinvH2O, 1424
.set nb112_rinvH2H1, 1440
.set nb112_rinvH2H2, 1456
.set nb112_is3, 1472
.set nb112_ii3, 1476
.set nb112_nri, 1492
.set nb112_iinr, 1500
.set nb112_jindex, 1508
.set nb112_jjnr, 1516
.set nb112_shift, 1524
.set nb112_shiftvec, 1532
.set nb112_facel, 1540
.set nb112_innerjjnr, 1548
.set nb112_innerk, 1556
.set nb112_n, 1560
.set nb112_nn1, 1564
.set nb112_nouter, 1568
.set nb112_ninner, 1572

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1592,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb112_nouter(%rsp)
        movl %eax,nb112_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb112_nri(%rsp)
        movq %rsi,nb112_iinr(%rsp)
        movq %rdx,nb112_jindex(%rsp)
        movq %rcx,nb112_jjnr(%rsp)
        movq %r8,nb112_shift(%rsp)
        movq %r9,nb112_shiftvec(%rsp)
        movq nb112_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb112_facel(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb112_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb112_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb112_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb112_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb112_qqOO(%rsp)
        movaps %xmm4,nb112_qqOH(%rsp)
        movaps %xmm5,nb112_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb112_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb112_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb112_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb112_c6(%rsp)
        movaps %xmm1,nb112_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb112_half(%rsp)
        movss nb112_half(%rsp),%xmm1
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
        movaps %xmm1,nb112_half(%rsp)
        movaps %xmm3,nb112_three(%rsp)
        movaps %xmm4,nb112_six(%rsp)
        movaps %xmm5,nb112_twelve(%rsp)

_nb_kernel112_x86_64_sse.nb112_threadloop: 
        movq  nb112_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel112_x86_64_sse.nb112_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112_x86_64_sse.nb112_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112_n(%rsp)
        movl %ebx,nb112_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112_x86_64_sse.nb112_outerstart
        jmp _nb_kernel112_x86_64_sse.nb112_end

_nb_kernel112_x86_64_sse.nb112_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112_nouter(%rsp),%ebx
        movl %ebx,nb112_nouter(%rsp)

_nb_kernel112_x86_64_sse.nb112_outer: 
        movq  nb112_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb112_is3(%rsp)      ## store is3 

        movq  nb112_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb112_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb112_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb112_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb112_ixO(%rsp)
        movaps %xmm4,nb112_iyO(%rsp)
        movaps %xmm5,nb112_izO(%rsp)

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
        movaps %xmm0,nb112_ixH1(%rsp)
        movaps %xmm1,nb112_iyH1(%rsp)
        movaps %xmm2,nb112_izH1(%rsp)
        movaps %xmm3,nb112_ixH2(%rsp)
        movaps %xmm4,nb112_iyH2(%rsp)
        movaps %xmm5,nb112_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb112_vctot(%rsp)
        movaps %xmm4,nb112_Vvdwtot(%rsp)
        movaps %xmm4,nb112_fixO(%rsp)
        movaps %xmm4,nb112_fiyO(%rsp)
        movaps %xmm4,nb112_fizO(%rsp)
        movaps %xmm4,nb112_fixH1(%rsp)
        movaps %xmm4,nb112_fiyH1(%rsp)
        movaps %xmm4,nb112_fizH1(%rsp)
        movaps %xmm4,nb112_fixH2(%rsp)
        movaps %xmm4,nb112_fiyH2(%rsp)
        movaps %xmm4,nb112_fizH2(%rsp)

        movq  nb112_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb112_pos(%rbp),%rsi
        movq  nb112_faction(%rbp),%rdi
        movq  nb112_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb112_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb112_ninner(%rsp),%ecx
        movl  %ecx,nb112_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb112_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel112_x86_64_sse.nb112_unroll_loop
        jmp   _nb_kernel112_x86_64_sse.nb112_single_check
_nb_kernel112_x86_64_sse.nb112_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb112_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb112_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb112_pos(%rbp),%rsi        ## base of pos[] 

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


    subps nb112_ixO(%rsp),%xmm0
    subps nb112_iyO(%rsp),%xmm1
    subps nb112_izO(%rsp),%xmm2
    subps nb112_ixH1(%rsp),%xmm3
    subps nb112_iyH1(%rsp),%xmm4
    subps nb112_izH1(%rsp),%xmm5
    subps nb112_ixH2(%rsp),%xmm6
    subps nb112_iyH2(%rsp),%xmm7
    subps nb112_izH2(%rsp),%xmm8

        movaps %xmm0,nb112_dxOO(%rsp)
        movaps %xmm1,nb112_dyOO(%rsp)
        movaps %xmm2,nb112_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb112_dxH1O(%rsp)
        movaps %xmm4,nb112_dyH1O(%rsp)
        movaps %xmm5,nb112_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb112_dxH2O(%rsp)
        movaps %xmm7,nb112_dyH2O(%rsp)
        movaps %xmm8,nb112_dzH2O(%rsp)
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

        movaps  nb112_three(%rsp),%xmm9
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

        movaps  nb112_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOO 
        mulps   %xmm0,%xmm10 ## rinvH1O
    mulps   %xmm0,%xmm11 ## rinvH2O

        ## O interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9   ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    movaps %xmm9,%xmm12
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm9,%xmm12 ## rinv6
    mulps  nb112_qqOO(%rsp),%xmm0
    mulps  nb112_qqOH(%rsp),%xmm1
    mulps  nb112_qqOH(%rsp),%xmm2
    movaps %xmm12,%xmm13 ## rinv6
    mulps %xmm12,%xmm12 ## rinv12
        mulps  nb112_c6(%rsp),%xmm13
        mulps  nb112_c12(%rsp),%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14

        addps  nb112_Vvdwtot(%rsp),%xmm14
        mulps  nb112_six(%rsp),%xmm13
        mulps  nb112_twelve(%rsp),%xmm12
        movaps %xmm14,nb112_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal        

    addps  %xmm0,%xmm12

    mulps  %xmm12,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb112_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb112_vctot(%rsp)

        ## move j O forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm0 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm2    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm2    ## jzOa  -  jzOb  -
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

        mulps nb112_dxOO(%rsp),%xmm7
        mulps nb112_dyOO(%rsp),%xmm8
        mulps nb112_dzOO(%rsp),%xmm9
        mulps nb112_dxH1O(%rsp),%xmm10
        mulps nb112_dyH1O(%rsp),%xmm11
        mulps nb112_dzH1O(%rsp),%xmm12
        mulps nb112_dxH2O(%rsp),%xmm13
        mulps nb112_dyH2O(%rsp),%xmm14
        mulps nb112_dzH2O(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb112_fixO(%rsp),%xmm7
    addps nb112_fiyO(%rsp),%xmm8
    addps nb112_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb112_fixH1(%rsp),%xmm10
    addps nb112_fiyH1(%rsp),%xmm11
    addps nb112_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb112_fixH2(%rsp),%xmm13
    addps nb112_fiyH2(%rsp),%xmm14
    addps nb112_fizH2(%rsp),%xmm15

    movaps %xmm7,nb112_fixO(%rsp)
    movaps %xmm8,nb112_fiyO(%rsp)
    movaps %xmm9,nb112_fizO(%rsp)
    movaps %xmm10,nb112_fixH1(%rsp)
    movaps %xmm11,nb112_fiyH1(%rsp)
    movaps %xmm12,nb112_fizH1(%rsp)
    movaps %xmm13,nb112_fixH2(%rsp)
    movaps %xmm14,nb112_fiyH2(%rsp)
    movaps %xmm15,nb112_fizH2(%rsp)

    ## xmm3 = fOx , xmm4 = fOy, xmm5=fOz
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

    subps nb112_ixO(%rsp),%xmm0
    subps nb112_iyO(%rsp),%xmm1
    subps nb112_izO(%rsp),%xmm2
    subps nb112_ixH1(%rsp),%xmm3
    subps nb112_iyH1(%rsp),%xmm4
    subps nb112_izH1(%rsp),%xmm5
    subps nb112_ixH2(%rsp),%xmm6
    subps nb112_iyH2(%rsp),%xmm7
    subps nb112_izH2(%rsp),%xmm8

        movaps %xmm0,nb112_dxOH1(%rsp)
        movaps %xmm1,nb112_dyOH1(%rsp)
        movaps %xmm2,nb112_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb112_dxH1H1(%rsp)
        movaps %xmm4,nb112_dyH1H1(%rsp)
        movaps %xmm5,nb112_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb112_dxH2H1(%rsp)
        movaps %xmm7,nb112_dyH2H1(%rsp)
        movaps %xmm8,nb112_dzH2H1(%rsp)
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

        movaps  nb112_three(%rsp),%xmm9
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

        movaps  nb112_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOH1
        mulps   %xmm0,%xmm10 ## rinvH1H1
    mulps   %xmm0,%xmm11 ## rinvH2H1

        ## H1 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb112_qqOH(%rsp),%xmm0
    mulps  nb112_qqHH(%rsp),%xmm1
    mulps  nb112_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb112_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb112_vctot(%rsp)

        ## move j H1 forces to local temp variables 
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

        mulps nb112_dxOH1(%rsp),%xmm7
        mulps nb112_dyOH1(%rsp),%xmm8
        mulps nb112_dzOH1(%rsp),%xmm9
        mulps nb112_dxH1H1(%rsp),%xmm10
        mulps nb112_dyH1H1(%rsp),%xmm11
        mulps nb112_dzH1H1(%rsp),%xmm12
        mulps nb112_dxH2H1(%rsp),%xmm13
        mulps nb112_dyH2H1(%rsp),%xmm14
        mulps nb112_dzH2H1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb112_fixO(%rsp),%xmm7
    addps nb112_fiyO(%rsp),%xmm8
    addps nb112_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb112_fixH1(%rsp),%xmm10
    addps nb112_fiyH1(%rsp),%xmm11
    addps nb112_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb112_fixH2(%rsp),%xmm13
    addps nb112_fiyH2(%rsp),%xmm14
    addps nb112_fizH2(%rsp),%xmm15

    movaps %xmm7,nb112_fixO(%rsp)
    movaps %xmm8,nb112_fiyO(%rsp)
    movaps %xmm9,nb112_fizO(%rsp)
    movaps %xmm10,nb112_fixH1(%rsp)
    movaps %xmm11,nb112_fiyH1(%rsp)
    movaps %xmm12,nb112_fizH1(%rsp)
    movaps %xmm13,nb112_fixH2(%rsp)
    movaps %xmm14,nb112_fiyH2(%rsp)
    movaps %xmm15,nb112_fizH2(%rsp)

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

    subps nb112_ixO(%rsp),%xmm0
    subps nb112_iyO(%rsp),%xmm1
    subps nb112_izO(%rsp),%xmm2
    subps nb112_ixH1(%rsp),%xmm3
    subps nb112_iyH1(%rsp),%xmm4
    subps nb112_izH1(%rsp),%xmm5
    subps nb112_ixH2(%rsp),%xmm6
    subps nb112_iyH2(%rsp),%xmm7
    subps nb112_izH2(%rsp),%xmm8

        movaps %xmm0,nb112_dxOH2(%rsp)
        movaps %xmm1,nb112_dyOH2(%rsp)
        movaps %xmm2,nb112_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb112_dxH1H2(%rsp)
        movaps %xmm4,nb112_dyH1H2(%rsp)
        movaps %xmm5,nb112_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb112_dxH2H2(%rsp)
        movaps %xmm7,nb112_dyH2H2(%rsp)
        movaps %xmm8,nb112_dzH2H2(%rsp)
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

        movaps  nb112_three(%rsp),%xmm9
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

        movaps  nb112_half(%rsp),%xmm0
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
    mulps  nb112_qqOH(%rsp),%xmm0
    mulps  nb112_qqHH(%rsp),%xmm1
    mulps  nb112_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb112_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb112_vctot(%rsp)

        ## move j H2 forces to local temp variables 
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

        mulps nb112_dxOH2(%rsp),%xmm7
        mulps nb112_dyOH2(%rsp),%xmm8
        mulps nb112_dzOH2(%rsp),%xmm9
        mulps nb112_dxH1H2(%rsp),%xmm10
        mulps nb112_dyH1H2(%rsp),%xmm11
        mulps nb112_dzH1H2(%rsp),%xmm12
        mulps nb112_dxH2H2(%rsp),%xmm13
        mulps nb112_dyH2H2(%rsp),%xmm14
        mulps nb112_dzH2H2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb112_fixO(%rsp),%xmm7
    addps nb112_fiyO(%rsp),%xmm8
    addps nb112_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb112_fixH1(%rsp),%xmm10
    addps nb112_fiyH1(%rsp),%xmm11
    addps nb112_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb112_fixH2(%rsp),%xmm13
    addps nb112_fiyH2(%rsp),%xmm14
    addps nb112_fizH2(%rsp),%xmm15

    movaps %xmm7,nb112_fixO(%rsp)
    movaps %xmm8,nb112_fiyO(%rsp)
    movaps %xmm9,nb112_fizO(%rsp)
    movaps %xmm10,nb112_fixH1(%rsp)
    movaps %xmm11,nb112_fiyH1(%rsp)
    movaps %xmm12,nb112_fizH1(%rsp)
    movaps %xmm13,nb112_fixH2(%rsp)
    movaps %xmm14,nb112_fiyH2(%rsp)
    movaps %xmm15,nb112_fizH2(%rsp)

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
        subl $4,nb112_innerk(%rsp)
        jl    _nb_kernel112_x86_64_sse.nb112_single_check
        jmp   _nb_kernel112_x86_64_sse.nb112_unroll_loop
_nb_kernel112_x86_64_sse.nb112_single_check: 
        addl $4,nb112_innerk(%rsp)
        jnz   _nb_kernel112_x86_64_sse.nb112_single_loop
        jmp   _nb_kernel112_x86_64_sse.nb112_updateouterdata
_nb_kernel112_x86_64_sse.nb112_single_loop: 
        movq  nb112_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb112_innerjjnr(%rsp)

        movq nb112_pos(%rbp),%rsi
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
        movaps %xmm0,nb112_jxO(%rsp)
        movaps %xmm1,nb112_jyO(%rsp)
        movaps %xmm2,nb112_jzO(%rsp)
        subps  nb112_ixO(%rsp),%xmm0
        subps  nb112_iyO(%rsp),%xmm1
        subps  nb112_izO(%rsp),%xmm2
        movaps %xmm0,nb112_dxOO(%rsp)
        movaps %xmm1,nb112_dyOO(%rsp)
        movaps %xmm2,nb112_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb112_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb112_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb112_qqOO(%rsp),%xmm4
        movss   %xmm0,%xmm1
        movhps  nb112_qqOH(%rsp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   nb112_c6(%rsp),%xmm1
        mulss   nb112_c12(%rsp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   nb112_Vvdwtot(%rsp),%xmm4
        mulss   nb112_six(%rsp),%xmm1
        mulss   nb112_twelve(%rsp),%xmm2
        movaps  %xmm4,nb112_Vvdwtot(%rsp)
        subss   %xmm1,%xmm2     ## fsD+ fsR 
        addps   %xmm3,%xmm2     ## fsC+ fsD+ fsR 

        addps   nb112_vctot(%rsp),%xmm3
        mulps   %xmm2,%xmm0     ## total fscal 
        movaps  %xmm3,nb112_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb112_dxOO(%rsp),%xmm0
        mulps   nb112_dyOO(%rsp),%xmm1
        mulps   nb112_dzOO(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb112_fjxO(%rsp)
        movaps  %xmm4,nb112_fjyO(%rsp)
        movaps  %xmm5,nb112_fjzO(%rsp)
        addps   nb112_fixO(%rsp),%xmm0
        addps   nb112_fiyO(%rsp),%xmm1
        addps   nb112_fizO(%rsp),%xmm2
        movaps  %xmm0,nb112_fixO(%rsp)
        movaps  %xmm1,nb112_fiyO(%rsp)
        movaps  %xmm2,nb112_fizO(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb112_jxO(%rsp),%xmm0
    movaps  nb112_jyO(%rsp),%xmm1
    movaps  nb112_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  nb112_ixH1(%rsp),%xmm0
        subps  nb112_iyH1(%rsp),%xmm1
        subps  nb112_izH1(%rsp),%xmm2
        subps  nb112_ixH2(%rsp),%xmm3
        subps  nb112_iyH2(%rsp),%xmm4
        subps  nb112_izH2(%rsp),%xmm5
    movaps %xmm0,nb112_dxH1O(%rsp)
        movaps %xmm1,nb112_dyH1O(%rsp)
        movaps %xmm2,nb112_dzH1O(%rsp)
        movaps %xmm3,nb112_dxH2O(%rsp)
        movaps %xmm4,nb112_dyH2O(%rsp)
        movaps %xmm5,nb112_dzH2O(%rsp)
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
        movaps  %xmm1,%xmm2  ## do coulomb interaction 
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb112_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb112_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb112_qqHH(%rsp),%xmm6
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   nb112_vctot(%rsp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,nb112_vctot(%rsp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb112_dxH1O(%rsp),%xmm0
        mulps   nb112_dyH1O(%rsp),%xmm1
        mulps   nb112_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb112_fjxO(%rsp),%xmm3
        movaps  nb112_fjyO(%rsp),%xmm4
        movaps  nb112_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb112_fjxO(%rsp)
        movaps  %xmm4,nb112_fjyO(%rsp)
        movaps  %xmm5,nb112_fjzO(%rsp)
        addps   nb112_fixH1(%rsp),%xmm0
        addps   nb112_fiyH1(%rsp),%xmm1
        addps   nb112_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb112_fixH1(%rsp)
        movaps  %xmm1,nb112_fiyH1(%rsp)
        movaps  %xmm2,nb112_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb112_dxH2O(%rsp),%xmm0
        mulps   nb112_dyH2O(%rsp),%xmm1
        mulps   nb112_dzH2O(%rsp),%xmm2
        movaps  nb112_fjxO(%rsp),%xmm3
        movaps  nb112_fjyO(%rsp),%xmm4
        movaps  nb112_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb112_faction(%rbp),%rsi
        movaps  %xmm3,nb112_fjxO(%rsp)
        movaps  %xmm4,nb112_fjyO(%rsp)
        movaps  %xmm5,nb112_fjzO(%rsp)
        addps   nb112_fixH2(%rsp),%xmm0
        addps   nb112_fiyH2(%rsp),%xmm1
        addps   nb112_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb112_fixH2(%rsp)
        movaps  %xmm1,nb112_fiyH2(%rsp)
        movaps  %xmm2,nb112_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb112_fjxO(%rsp),%xmm3
        movaps  nb112_fjyO(%rsp),%xmm4
        movaps  nb112_fjzO(%rsp),%xmm5
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

        decl nb112_innerk(%rsp)
        jz    _nb_kernel112_x86_64_sse.nb112_updateouterdata
        jmp   _nb_kernel112_x86_64_sse.nb112_single_loop
_nb_kernel112_x86_64_sse.nb112_updateouterdata: 
        movl  nb112_ii3(%rsp),%ecx
        movq  nb112_faction(%rbp),%rdi
        movq  nb112_fshift(%rbp),%rsi
        movl  nb112_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb112_fixO(%rsp),%xmm0
        movaps nb112_fiyO(%rsp),%xmm1
        movaps nb112_fizO(%rsp),%xmm2

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
        movaps nb112_fixH1(%rsp),%xmm0
        movaps nb112_fiyH1(%rsp),%xmm1
        movaps nb112_fizH1(%rsp),%xmm2

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
        movaps nb112_fixH2(%rsp),%xmm0
        movaps nb112_fiyH2(%rsp),%xmm1
        movaps nb112_fizH2(%rsp),%xmm2

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
        movl nb112_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb112_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb112_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb112_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb112_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb112_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb112_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112_x86_64_sse.nb112_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112_n(%rsp)
        jmp _nb_kernel112_x86_64_sse.nb112_outer
_nb_kernel112_x86_64_sse.nb112_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112_x86_64_sse.nb112_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112_x86_64_sse.nb112_threadloop
_nb_kernel112_x86_64_sse.nb112_end: 


        movl nb112_nouter(%rsp),%eax
        movl nb112_ninner(%rsp),%ebx
        movq nb112_outeriter(%rbp),%rcx
        movq nb112_inneriter(%rbp),%rdx
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





.globl nb_kernel112nf_x86_64_sse
.globl _nb_kernel112nf_x86_64_sse
nb_kernel112nf_x86_64_sse:      
_nb_kernel112nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb112nf_fshift, 16
.set nb112nf_gid, 24
.set nb112nf_pos, 32
.set nb112nf_faction, 40
.set nb112nf_charge, 48
.set nb112nf_p_facel, 56
.set nb112nf_argkrf, 64
.set nb112nf_argcrf, 72
.set nb112nf_Vc, 80
.set nb112nf_type, 88
.set nb112nf_p_ntype, 96
.set nb112nf_vdwparam, 104
.set nb112nf_Vvdw, 112
.set nb112nf_p_tabscale, 120
.set nb112nf_VFtab, 128
.set nb112nf_invsqrta, 136
.set nb112nf_dvda, 144
.set nb112nf_p_gbtabscale, 152
.set nb112nf_GBtab, 160
.set nb112nf_p_nthreads, 168
.set nb112nf_count, 176
.set nb112nf_mtx, 184
.set nb112nf_outeriter, 192
.set nb112nf_inneriter, 200
.set nb112nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb112nf_ixO, 0
.set nb112nf_iyO, 16
.set nb112nf_izO, 32
.set nb112nf_ixH1, 48
.set nb112nf_iyH1, 64
.set nb112nf_izH1, 80
.set nb112nf_ixH2, 96
.set nb112nf_iyH2, 112
.set nb112nf_izH2, 128
.set nb112nf_jxO, 144
.set nb112nf_jyO, 160
.set nb112nf_jzO, 176
.set nb112nf_jxH1, 192
.set nb112nf_jyH1, 208
.set nb112nf_jzH1, 224
.set nb112nf_jxH2, 240
.set nb112nf_jyH2, 256
.set nb112nf_jzH2, 272
.set nb112nf_qqOO, 288
.set nb112nf_qqOH, 304
.set nb112nf_qqHH, 320
.set nb112nf_c6, 336
.set nb112nf_c12, 352
.set nb112nf_vctot, 368
.set nb112nf_Vvdwtot, 384
.set nb112nf_half, 400
.set nb112nf_three, 416
.set nb112nf_rsqOO, 432
.set nb112nf_rsqOH1, 448
.set nb112nf_rsqOH2, 464
.set nb112nf_rsqH1O, 480
.set nb112nf_rsqH1H1, 496
.set nb112nf_rsqH1H2, 512
.set nb112nf_rsqH2O, 528
.set nb112nf_rsqH2H1, 544
.set nb112nf_rsqH2H2, 560
.set nb112nf_rinvOO, 576
.set nb112nf_rinvOH1, 592
.set nb112nf_rinvOH2, 608
.set nb112nf_rinvH1O, 624
.set nb112nf_rinvH1H1, 640
.set nb112nf_rinvH1H2, 656
.set nb112nf_rinvH2O, 672
.set nb112nf_rinvH2H1, 688
.set nb112nf_rinvH2H2, 704
.set nb112nf_is3, 720
.set nb112nf_ii3, 724
.set nb112nf_nri, 740
.set nb112nf_iinr, 748
.set nb112nf_jindex, 756
.set nb112nf_jjnr, 764
.set nb112nf_shift, 772
.set nb112nf_shiftvec, 780
.set nb112nf_facel, 788
.set nb112nf_innerjjnr, 796
.set nb112nf_innerk, 804
.set nb112nf_n, 812
.set nb112nf_nn1, 816
.set nb112nf_nouter, 820
.set nb112nf_ninner, 824

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $840,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb112nf_nouter(%rsp)
        movl %eax,nb112nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb112nf_nri(%rsp)
        movq %rsi,nb112nf_iinr(%rsp)
        movq %rdx,nb112nf_jindex(%rsp)
        movq %rcx,nb112nf_jjnr(%rsp)
        movq %r8,nb112nf_shift(%rsp)
        movq %r9,nb112nf_shiftvec(%rsp)
        movq nb112nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb112nf_facel(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  nb112nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb112nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb112nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb112nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb112nf_qqOO(%rsp)
        movaps %xmm4,nb112nf_qqOH(%rsp)
        movaps %xmm5,nb112nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb112nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb112nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb112nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb112nf_c6(%rsp)
        movaps %xmm1,nb112nf_c12(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb112nf_half(%rsp)
        movss nb112nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb112nf_half(%rsp)
        movaps %xmm3,nb112nf_three(%rsp)

_nb_kernel112nf_x86_64_sse.nb112nf_threadloop: 
        movq  nb112nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel112nf_x86_64_sse.nb112nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel112nf_x86_64_sse.nb112nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb112nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb112nf_n(%rsp)
        movl %ebx,nb112nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel112nf_x86_64_sse.nb112nf_outerstart
        jmp _nb_kernel112nf_x86_64_sse.nb112nf_end

_nb_kernel112nf_x86_64_sse.nb112nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb112nf_nouter(%rsp),%ebx
        movl %ebx,nb112nf_nouter(%rsp)

_nb_kernel112nf_x86_64_sse.nb112nf_outer: 
        movq  nb112nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb112nf_is3(%rsp)            ## store is3 

        movq  nb112nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb112nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb112nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb112nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb112nf_ixO(%rsp)
        movaps %xmm4,nb112nf_iyO(%rsp)
        movaps %xmm5,nb112nf_izO(%rsp)

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
        movaps %xmm0,nb112nf_ixH1(%rsp)
        movaps %xmm1,nb112nf_iyH1(%rsp)
        movaps %xmm2,nb112nf_izH1(%rsp)
        movaps %xmm3,nb112nf_ixH2(%rsp)
        movaps %xmm4,nb112nf_iyH2(%rsp)
        movaps %xmm5,nb112nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb112nf_vctot(%rsp)
        movaps %xmm4,nb112nf_Vvdwtot(%rsp)

        movq  nb112nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb112nf_pos(%rbp),%rsi
        movq  nb112nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb112nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb112nf_ninner(%rsp),%ecx
        movl  %ecx,nb112nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb112nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel112nf_x86_64_sse.nb112nf_unroll_loop
        jmp   _nb_kernel112nf_x86_64_sse.nb112nf_single_check
_nb_kernel112nf_x86_64_sse.nb112nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb112nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb112nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb112nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps %xmm0,nb112nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb112nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb112nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb112nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb112nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb112nf_jyH2(%rsp)

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
        movaps %xmm0,nb112nf_jzO(%rsp)
        movaps %xmm1,nb112nf_jzH1(%rsp)
        movaps %xmm2,nb112nf_jzH2(%rsp)

        movaps nb112nf_ixO(%rsp),%xmm0
        movaps nb112nf_iyO(%rsp),%xmm1
        movaps nb112nf_izO(%rsp),%xmm2
        movaps nb112nf_ixO(%rsp),%xmm3
        movaps nb112nf_iyO(%rsp),%xmm4
        movaps nb112nf_izO(%rsp),%xmm5
        subps  nb112nf_jxO(%rsp),%xmm0
        subps  nb112nf_jyO(%rsp),%xmm1
        subps  nb112nf_jzO(%rsp),%xmm2
        subps  nb112nf_jxH1(%rsp),%xmm3
        subps  nb112nf_jyH1(%rsp),%xmm4
        subps  nb112nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb112nf_rsqOO(%rsp)
        movaps %xmm3,nb112nf_rsqOH1(%rsp)

        movaps nb112nf_ixO(%rsp),%xmm0
        movaps nb112nf_iyO(%rsp),%xmm1
        movaps nb112nf_izO(%rsp),%xmm2
        movaps nb112nf_ixH1(%rsp),%xmm3
        movaps nb112nf_iyH1(%rsp),%xmm4
        movaps nb112nf_izH1(%rsp),%xmm5
        subps  nb112nf_jxH2(%rsp),%xmm0
        subps  nb112nf_jyH2(%rsp),%xmm1
        subps  nb112nf_jzH2(%rsp),%xmm2
        subps  nb112nf_jxO(%rsp),%xmm3
        subps  nb112nf_jyO(%rsp),%xmm4
        subps  nb112nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb112nf_rsqOH2(%rsp)
        movaps %xmm3,nb112nf_rsqH1O(%rsp)

        movaps nb112nf_ixH1(%rsp),%xmm0
        movaps nb112nf_iyH1(%rsp),%xmm1
        movaps nb112nf_izH1(%rsp),%xmm2
        movaps nb112nf_ixH1(%rsp),%xmm3
        movaps nb112nf_iyH1(%rsp),%xmm4
        movaps nb112nf_izH1(%rsp),%xmm5
        subps  nb112nf_jxH1(%rsp),%xmm0
        subps  nb112nf_jyH1(%rsp),%xmm1
        subps  nb112nf_jzH1(%rsp),%xmm2
        subps  nb112nf_jxH2(%rsp),%xmm3
        subps  nb112nf_jyH2(%rsp),%xmm4
        subps  nb112nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb112nf_rsqH1H1(%rsp)
        movaps %xmm3,nb112nf_rsqH1H2(%rsp)

        movaps nb112nf_ixH2(%rsp),%xmm0
        movaps nb112nf_iyH2(%rsp),%xmm1
        movaps nb112nf_izH2(%rsp),%xmm2
        movaps nb112nf_ixH2(%rsp),%xmm3
        movaps nb112nf_iyH2(%rsp),%xmm4
        movaps nb112nf_izH2(%rsp),%xmm5
        subps  nb112nf_jxO(%rsp),%xmm0
        subps  nb112nf_jyO(%rsp),%xmm1
        subps  nb112nf_jzO(%rsp),%xmm2
        subps  nb112nf_jxH1(%rsp),%xmm3
        subps  nb112nf_jyH1(%rsp),%xmm4
        subps  nb112nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb112nf_rsqH2O(%rsp)
        movaps %xmm4,nb112nf_rsqH2H1(%rsp)

        movaps nb112nf_ixH2(%rsp),%xmm0
        movaps nb112nf_iyH2(%rsp),%xmm1
        movaps nb112nf_izH2(%rsp),%xmm2
        subps  nb112nf_jxH2(%rsp),%xmm0
        subps  nb112nf_jyH2(%rsp),%xmm1
        subps  nb112nf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb112nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb112nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb112nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb112nf_rinvH2H1(%rsp)

        rsqrtps nb112nf_rsqOO(%rsp),%xmm1
        rsqrtps nb112nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb112nf_rsqOO(%rsp),%xmm1
        mulps   nb112nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112nf_half(%rsp),%xmm3
        mulps   nb112nf_half(%rsp),%xmm7
        movaps  %xmm3,nb112nf_rinvOO(%rsp)
        movaps  %xmm7,nb112nf_rinvOH1(%rsp)

        rsqrtps nb112nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb112nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb112nf_rsqOH2(%rsp),%xmm1
        mulps   nb112nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112nf_half(%rsp),%xmm3
        mulps   nb112nf_half(%rsp),%xmm7
        movaps  %xmm3,nb112nf_rinvOH2(%rsp)
        movaps  %xmm7,nb112nf_rinvH1O(%rsp)

        rsqrtps nb112nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb112nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb112nf_rsqH1H1(%rsp),%xmm1
        mulps   nb112nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112nf_half(%rsp),%xmm3
        mulps   nb112nf_half(%rsp),%xmm7
        movaps  %xmm3,nb112nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb112nf_rinvH1H2(%rsp)

        rsqrtps nb112nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb112nf_three(%rsp),%xmm3
        mulps   nb112nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb112nf_half(%rsp),%xmm3
        movaps  %xmm3,nb112nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movaps nb112nf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm7
        mulps  %xmm0,%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulps  nb112nf_qqOO(%rsp),%xmm7
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb112nf_c6(%rsp),%xmm1
        mulps  nb112nf_c12(%rsp),%xmm2
        subps  %xmm1,%xmm2      ## xmm3=Vvdw12-Vvdw6 
        addps  nb112nf_Vvdwtot(%rsp),%xmm2
        movaps %xmm2,nb112nf_Vvdwtot(%rsp)
        addps  nb112nf_vctot(%rsp),%xmm7

        ## all other interaction 
        movaps nb112nf_rinvOH1(%rsp),%xmm0
        movaps nb112nf_rinvH1H1(%rsp),%xmm1
        addps  nb112nf_rinvOH2(%rsp),%xmm0
        addps  nb112nf_rinvH1H2(%rsp),%xmm1
        addps  nb112nf_rinvH1O(%rsp),%xmm0
        addps  nb112nf_rinvH2H1(%rsp),%xmm1
        addps  nb112nf_rinvH2O(%rsp),%xmm0
        addps  nb112nf_rinvH2H2(%rsp),%xmm1

        mulps nb112nf_qqOH(%rsp),%xmm0
        mulps nb112nf_qqHH(%rsp),%xmm1
        addps %xmm0,%xmm7
        addps %xmm1,%xmm7
        movaps %xmm7,nb112nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb112nf_innerk(%rsp)
        jl    _nb_kernel112nf_x86_64_sse.nb112nf_single_check
        jmp   _nb_kernel112nf_x86_64_sse.nb112nf_unroll_loop
_nb_kernel112nf_x86_64_sse.nb112nf_single_check: 
        addl $4,nb112nf_innerk(%rsp)
        jnz   _nb_kernel112nf_x86_64_sse.nb112nf_single_loop
        jmp   _nb_kernel112nf_x86_64_sse.nb112nf_updateouterdata
_nb_kernel112nf_x86_64_sse.nb112nf_single_loop: 
        movq  nb112nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb112nf_innerjjnr(%rsp)

        movq nb112nf_pos(%rbp),%rsi
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
        movaps  nb112nf_ixO(%rsp),%xmm0
        movaps  nb112nf_iyO(%rsp),%xmm1
        movaps  nb112nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb112nf_jxO(%rsp)
        movaps %xmm4,nb112nf_jyO(%rsp)
        movaps %xmm5,nb112nf_jzO(%rsp)
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
        movaps  nb112nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb112nf_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb112nf_qqOO(%rsp),%xmm4
        movss   %xmm0,%xmm1
        movhps  nb112nf_qqOH(%rsp),%xmm4
        mulss   %xmm0,%xmm1
        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulss   %xmm0,%xmm1     ## xmm1(0)=rinvsix 
        movaps  %xmm1,%xmm2     ## zero everything else in xmm2 
        mulss   %xmm2,%xmm2     ## xmm2=rinvtwelve 

        mulss   nb112nf_c6(%rsp),%xmm1
        mulss   nb112nf_c12(%rsp),%xmm2
        movaps  %xmm2,%xmm4
        subss   %xmm1,%xmm4     ## Vvdwtot=Vvdw12-Vvdw6 
        addps   nb112nf_Vvdwtot(%rsp),%xmm4
        movaps  %xmm4,nb112nf_Vvdwtot(%rsp)

        addps   nb112nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb112nf_vctot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb112nf_ixH1(%rsp),%xmm0
        movaps  nb112nf_iyH1(%rsp),%xmm1
        movaps  nb112nf_izH1(%rsp),%xmm2
        movaps  nb112nf_ixH2(%rsp),%xmm3
        movaps  nb112nf_iyH2(%rsp),%xmm4
        movaps  nb112nf_izH2(%rsp),%xmm5
        subps   nb112nf_jxO(%rsp),%xmm0
        subps   nb112nf_jyO(%rsp),%xmm1
        subps   nb112nf_jzO(%rsp),%xmm2
        subps   nb112nf_jxO(%rsp),%xmm3
        subps   nb112nf_jyO(%rsp),%xmm4
        subps   nb112nf_jzO(%rsp),%xmm5
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
        movaps  %xmm1,%xmm2  ## do coulomb interaction 
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb112nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb112nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb112nf_half(%rsp),%xmm7   ## rinv H2 - j water  
        addps   %xmm7,%xmm3
        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb112nf_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb112nf_qqHH(%rsp),%xmm6
        mulps   %xmm6,%xmm3     ## total vcoul 

        addps   nb112nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb112nf_vctot(%rsp)

        decl nb112nf_innerk(%rsp)
        jz    _nb_kernel112nf_x86_64_sse.nb112nf_updateouterdata
        jmp   _nb_kernel112nf_x86_64_sse.nb112nf_single_loop
_nb_kernel112nf_x86_64_sse.nb112nf_updateouterdata: 
        ## get n from stack
        movl nb112nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb112nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb112nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb112nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb112nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb112nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb112nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel112nf_x86_64_sse.nb112nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb112nf_n(%rsp)
        jmp _nb_kernel112nf_x86_64_sse.nb112nf_outer
_nb_kernel112nf_x86_64_sse.nb112nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb112nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel112nf_x86_64_sse.nb112nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel112nf_x86_64_sse.nb112nf_threadloop
_nb_kernel112nf_x86_64_sse.nb112nf_end: 



        movl nb112nf_nouter(%rsp),%eax
        movl nb112nf_ninner(%rsp),%ebx
        movq nb112nf_outeriter(%rbp),%rcx
        movq nb112nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $840,%rsp
        emms
        pop %rbx
        pop    %rbp
        ret

