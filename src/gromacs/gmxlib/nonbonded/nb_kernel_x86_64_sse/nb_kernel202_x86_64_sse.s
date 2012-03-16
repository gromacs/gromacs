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






.globl nb_kernel202_x86_64_sse
.globl _nb_kernel202_x86_64_sse
nb_kernel202_x86_64_sse:        
_nb_kernel202_x86_64_sse:       
.set nb202_fshift, 16
.set nb202_gid, 24
.set nb202_pos, 32
.set nb202_faction, 40
.set nb202_charge, 48
.set nb202_p_facel, 56
.set nb202_argkrf, 64
.set nb202_argcrf, 72
.set nb202_Vc, 80
.set nb202_type, 88
.set nb202_p_ntype, 96
.set nb202_vdwparam, 104
.set nb202_Vvdw, 112
.set nb202_p_tabscale, 120
.set nb202_VFtab, 128
.set nb202_invsqrta, 136
.set nb202_dvda, 144
.set nb202_p_gbtabscale, 152
.set nb202_GBtab, 160
.set nb202_p_nthreads, 168
.set nb202_count, 176
.set nb202_mtx, 184
.set nb202_outeriter, 192
.set nb202_inneriter, 200
.set nb202_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb202_ixO, 0
.set nb202_iyO, 16
.set nb202_izO, 32
.set nb202_ixH1, 48
.set nb202_iyH1, 64
.set nb202_izH1, 80
.set nb202_ixH2, 96
.set nb202_iyH2, 112
.set nb202_izH2, 128
.set nb202_jxO, 144
.set nb202_jyO, 160
.set nb202_jzO, 176
.set nb202_jxH1, 192
.set nb202_jyH1, 208
.set nb202_jzH1, 224
.set nb202_jxH2, 240
.set nb202_jyH2, 256
.set nb202_jzH2, 272
.set nb202_dxOO, 288
.set nb202_dyOO, 304
.set nb202_dzOO, 320
.set nb202_dxOH1, 336
.set nb202_dyOH1, 352
.set nb202_dzOH1, 368
.set nb202_dxOH2, 384
.set nb202_dyOH2, 400
.set nb202_dzOH2, 416
.set nb202_dxH1O, 432
.set nb202_dyH1O, 448
.set nb202_dzH1O, 464
.set nb202_dxH1H1, 480
.set nb202_dyH1H1, 496
.set nb202_dzH1H1, 512
.set nb202_dxH1H2, 528
.set nb202_dyH1H2, 544
.set nb202_dzH1H2, 560
.set nb202_dxH2O, 576
.set nb202_dyH2O, 592
.set nb202_dzH2O, 608
.set nb202_dxH2H1, 624
.set nb202_dyH2H1, 640
.set nb202_dzH2H1, 656
.set nb202_dxH2H2, 672
.set nb202_dyH2H2, 688
.set nb202_dzH2H2, 704
.set nb202_qqOO, 720
.set nb202_qqOH, 736
.set nb202_qqHH, 752
.set nb202_vctot, 768
.set nb202_fixO, 784
.set nb202_fiyO, 800
.set nb202_fizO, 816
.set nb202_fixH1, 832
.set nb202_fiyH1, 848
.set nb202_fizH1, 864
.set nb202_fixH2, 880
.set nb202_fiyH2, 896
.set nb202_fizH2, 912
.set nb202_fjxO, 928
.set nb202_fjyO, 944
.set nb202_fjzO, 960
.set nb202_fjxH1, 976
.set nb202_fjyH1, 992
.set nb202_fjzH1, 1008
.set nb202_fjxH2, 1024
.set nb202_fjyH2, 1040
.set nb202_fjzH2, 1056
.set nb202_half, 1072
.set nb202_three, 1088
.set nb202_rsqOO, 1104
.set nb202_rsqOH1, 1120
.set nb202_rsqOH2, 1136
.set nb202_rsqH1O, 1152
.set nb202_rsqH1H1, 1168
.set nb202_rsqH1H2, 1184
.set nb202_rsqH2O, 1200
.set nb202_rsqH2H1, 1216
.set nb202_rsqH2H2, 1232
.set nb202_rinvOO, 1248
.set nb202_rinvOH1, 1264
.set nb202_rinvOH2, 1280
.set nb202_rinvH1O, 1296
.set nb202_rinvH1H1, 1312
.set nb202_rinvH1H2, 1328
.set nb202_rinvH2O, 1344
.set nb202_rinvH2H1, 1360
.set nb202_rinvH2H2, 1376
.set nb202_two, 1392
.set nb202_krf, 1408
.set nb202_crf, 1424
.set nb202_innerjjnr, 1440
.set nb202_nri, 1448
.set nb202_iinr, 1456
.set nb202_jindex, 1464
.set nb202_jjnr, 1472
.set nb202_shift, 1480
.set nb202_shiftvec, 1488
.set nb202_facel, 1496
.set nb202_is3, 1504
.set nb202_ii3, 1508
.set nb202_innerk, 1512
.set nb202_n, 1516
.set nb202_nn1, 1520
.set nb202_nouter, 1524
.set nb202_ninner, 1528

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15


        subq $1544,%rsp         ## local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb202_nouter(%rsp)
        movl %eax,nb202_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb202_nri(%rsp)
        movq %rsi,nb202_iinr(%rsp)
        movq %rdx,nb202_jindex(%rsp)
        movq %rcx,nb202_jjnr(%rsp)
        movq %r8,nb202_shift(%rsp)
        movq %r9,nb202_shiftvec(%rsp)
        movq nb202_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb202_facel(%rsp)


        movq nb202_argkrf(%rbp),%rsi
        movq nb202_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb202_krf(%rsp)
        movaps %xmm2,nb202_crf(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb202_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb202_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb202_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb202_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202_qqOO(%rsp)
        movaps %xmm4,nb202_qqOH(%rsp)
        movaps %xmm5,nb202_qqHH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb202_half(%rsp)
        movss nb202_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb202_half(%rsp)
        movaps %xmm2,nb202_two(%rsp)
        movaps %xmm3,nb202_three(%rsp)

_nb_kernel202_x86_64_sse.nb202_threadloop: 
        movq  nb202_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel202_x86_64_sse.nb202_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202_x86_64_sse.nb202_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb202_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb202_n(%rsp)
        movl %ebx,nb202_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel202_x86_64_sse.nb202_outerstart
        jmp _nb_kernel202_x86_64_sse.nb202_end

_nb_kernel202_x86_64_sse.nb202_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202_nouter(%rsp),%ebx
        movl %ebx,nb202_nouter(%rsp)

_nb_kernel202_x86_64_sse.nb202_outer: 
        movq  nb202_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb202_is3(%rsp)      ## store is3 

        movq  nb202_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb202_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb202_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb202_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202_ixO(%rsp)
        movaps %xmm4,nb202_iyO(%rsp)
        movaps %xmm5,nb202_izO(%rsp)

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
        movaps %xmm0,nb202_ixH1(%rsp)
        movaps %xmm1,nb202_iyH1(%rsp)
        movaps %xmm2,nb202_izH1(%rsp)
        movaps %xmm3,nb202_ixH2(%rsp)
        movaps %xmm4,nb202_iyH2(%rsp)
        movaps %xmm5,nb202_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb202_vctot(%rsp)
        movaps %xmm4,nb202_fixO(%rsp)
        movaps %xmm4,nb202_fiyO(%rsp)
        movaps %xmm4,nb202_fizO(%rsp)
        movaps %xmm4,nb202_fixH1(%rsp)
        movaps %xmm4,nb202_fiyH1(%rsp)
        movaps %xmm4,nb202_fizH1(%rsp)
        movaps %xmm4,nb202_fixH2(%rsp)
        movaps %xmm4,nb202_fiyH2(%rsp)
        movaps %xmm4,nb202_fizH2(%rsp)

        movq  nb202_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movq  nb202_pos(%rbp),%rsi
        movq  nb202_faction(%rbp),%rdi
        movq  nb202_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb202_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb202_ninner(%rsp),%ecx
        movl  %ecx,nb202_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb202_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel202_x86_64_sse.nb202_unroll_loop
        jmp   _nb_kernel202_x86_64_sse.nb202_single_check
_nb_kernel202_x86_64_sse.nb202_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb202_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb202_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb202_pos(%rbp),%rsi        ## base of pos[] 

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

    subps nb202_ixO(%rsp),%xmm0
    subps nb202_iyO(%rsp),%xmm1
    subps nb202_izO(%rsp),%xmm2
    subps nb202_ixH1(%rsp),%xmm3
    subps nb202_iyH1(%rsp),%xmm4
    subps nb202_izH1(%rsp),%xmm5
    subps nb202_ixH2(%rsp),%xmm6
    subps nb202_iyH2(%rsp),%xmm7
    subps nb202_izH2(%rsp),%xmm8

        movaps %xmm0,nb202_dxOO(%rsp)
        movaps %xmm1,nb202_dyOO(%rsp)
        movaps %xmm2,nb202_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH1O(%rsp)
        movaps %xmm4,nb202_dyH1O(%rsp)
        movaps %xmm5,nb202_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb202_dxH2O(%rsp)
        movaps %xmm7,nb202_dyH2O(%rsp)
        movaps %xmm8,nb202_dzH2O(%rsp)
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

        movaps  nb202_three(%rsp),%xmm9
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

        movaps  nb202_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOO 
        mulps   %xmm4,%xmm10 ## rinvH1O
    mulps   %xmm4,%xmm11 ## rinvH2O

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb202_krf(%rsp),%xmm2
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
    movaps nb202_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb202_qqOO(%rsp),%xmm12
    movaps nb202_qqOH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm13,%xmm4
    mulps  %xmm13,%xmm7
    addps  nb202_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,%xmm15

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j O forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm9 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm10 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm9 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm10 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm11    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm12    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm11    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm12    ## jzOc  -  jzOd -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm9: jxOa jyOa jxOb jyOb 
    ## xmm10: jxOc jyOc jxOd jyOd
    ## xmm11: jzOa jzOb jzOc jzOd

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb202_dxOO(%rsp),%xmm0
        mulps nb202_dyOO(%rsp),%xmm1
        mulps nb202_dzOO(%rsp),%xmm2
        mulps nb202_dxH1O(%rsp),%xmm3
        mulps nb202_dyH1O(%rsp),%xmm4
        mulps nb202_dzH1O(%rsp),%xmm5
        mulps nb202_dxH2O(%rsp),%xmm6
        mulps nb202_dyH2O(%rsp),%xmm7
        mulps nb202_dzH2O(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb202_fixO(%rsp),%xmm0
    addps nb202_fiyO(%rsp),%xmm1
    addps nb202_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb202_fixH1(%rsp),%xmm3
    addps nb202_fiyH1(%rsp),%xmm4
    addps nb202_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb202_fixH2(%rsp),%xmm6
    addps nb202_fiyH2(%rsp),%xmm7
    addps nb202_fizH2(%rsp),%xmm8

    movaps %xmm0,nb202_fixO(%rsp)
    movaps %xmm1,nb202_fiyO(%rsp)
    movaps %xmm2,nb202_fizO(%rsp)
    movaps %xmm3,nb202_fixH1(%rsp)
    movaps %xmm4,nb202_fiyH1(%rsp)
    movaps %xmm5,nb202_fizH1(%rsp)
    movaps %xmm6,nb202_fixH2(%rsp)
    movaps %xmm7,nb202_fiyH2(%rsp)
    movaps %xmm8,nb202_fizH2(%rsp)

    ## xmm9 = fOx
    ## xmm10 = fOy
    ## xmm11 = fOz
    movaps %xmm13,%xmm0
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm0

    addps %xmm13,%xmm9
    addps %xmm0,%xmm10

    movhlps  %xmm11,%xmm12 ## fOzc fOzd

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


    subps nb202_ixO(%rsp),%xmm0
    subps nb202_iyO(%rsp),%xmm1
    subps nb202_izO(%rsp),%xmm2
    subps nb202_ixH1(%rsp),%xmm3
    subps nb202_iyH1(%rsp),%xmm4
    subps nb202_izH1(%rsp),%xmm5
    subps nb202_ixH2(%rsp),%xmm6
    subps nb202_iyH2(%rsp),%xmm7
    subps nb202_izH2(%rsp),%xmm8

        movaps %xmm0,nb202_dxOH1(%rsp)
        movaps %xmm1,nb202_dyOH1(%rsp)
        movaps %xmm2,nb202_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH1H1(%rsp)
        movaps %xmm4,nb202_dyH1H1(%rsp)
        movaps %xmm5,nb202_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb202_dxH2H1(%rsp)
        movaps %xmm7,nb202_dyH2H1(%rsp)
        movaps %xmm8,nb202_dzH2H1(%rsp)
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

        movaps  nb202_three(%rsp),%xmm9
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

        movaps  nb202_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH1
        mulps   %xmm4,%xmm10 ## rinvH1H1
    mulps   %xmm4,%xmm11 ## rinvH2H1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb202_krf(%rsp),%xmm2
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
    movaps nb202_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb202_qqOH(%rsp),%xmm12
    movaps nb202_qqHH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm13,%xmm4
    mulps  %xmm13,%xmm7
    addps  %xmm2,%xmm15
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm15

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j H1 forces to local temp variables 
    movlps 12(%rdi,%rax,4),%xmm9    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm10    ## jxH1c jyH1c  -   -
    movhps 12(%rdi,%rbx,4),%xmm9    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm10    ## jxH1c jyH1c jxH1d jyH1d 

    movss  20(%rdi,%rax,4),%xmm11    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm12    ## jzH1c  -  -  -
    movhps 20(%rdi,%rbx,4),%xmm11    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm12    ## jzH1c  -  jzH1d -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm9: jxH1a jyH1a jxH1b jyH1b 
    ## xmm10: jxH1c jyH1c jxH1d jyH1d
    ## xmm11: jzH1a jzH1b jzH1c jzH1d

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb202_dxOH1(%rsp),%xmm0
        mulps nb202_dyOH1(%rsp),%xmm1
        mulps nb202_dzOH1(%rsp),%xmm2
        mulps nb202_dxH1H1(%rsp),%xmm3
        mulps nb202_dyH1H1(%rsp),%xmm4
        mulps nb202_dzH1H1(%rsp),%xmm5
        mulps nb202_dxH2H1(%rsp),%xmm6
        mulps nb202_dyH2H1(%rsp),%xmm7
        mulps nb202_dzH2H1(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb202_fixO(%rsp),%xmm0
    addps nb202_fiyO(%rsp),%xmm1
    addps nb202_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb202_fixH1(%rsp),%xmm3
    addps nb202_fiyH1(%rsp),%xmm4
    addps nb202_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb202_fixH2(%rsp),%xmm6
    addps nb202_fiyH2(%rsp),%xmm7
    addps nb202_fizH2(%rsp),%xmm8

    movaps %xmm0,nb202_fixO(%rsp)
    movaps %xmm1,nb202_fiyO(%rsp)
    movaps %xmm2,nb202_fizO(%rsp)
    movaps %xmm3,nb202_fixH1(%rsp)
    movaps %xmm4,nb202_fiyH1(%rsp)
    movaps %xmm5,nb202_fizH1(%rsp)
    movaps %xmm6,nb202_fixH2(%rsp)
    movaps %xmm7,nb202_fiyH2(%rsp)
    movaps %xmm8,nb202_fizH2(%rsp)

    ## xmm9  = fH1x
    ## xmm10 = fH1y
    ## xmm11 = fH1z
    movaps %xmm13,%xmm0
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm0

    addps %xmm13,%xmm9
    addps %xmm0,%xmm10

    movhlps  %xmm11,%xmm12 ## fH1zc fH1zd

    movlps %xmm9,12(%rdi,%rax,4)
    movhps %xmm9,12(%rdi,%rbx,4)
    movlps %xmm10,12(%rdi,%rcx,4)
    movhps %xmm10,12(%rdi,%rdx,4)
    movss  %xmm11,20(%rdi,%rax,4)
    movss  %xmm12,20(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,20(%rdi,%rbx,4)
    movss  %xmm12,20(%rdi,%rdx,4)

        ## move j H2 coordinates to local temp variables 
    movlps 24(%rsi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rsi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rsi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rsi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rsi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rsi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movss  32(%rsi,%rbx,4),%xmm5    ## jzOb  -  -  -
    movss  32(%rsi,%rdx,4),%xmm6    ## jzOd  -  -  -
    movlhps %xmm5,%xmm2 ## jzOa  -  jzOb  -
    movlhps %xmm6,%xmm3 ## jzOc  -  jzOd -

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH2a jxH2c jyH2a jyH2c        
    unpckhps %xmm1,%xmm4 ## jxH2b jxH2d jyH2b jyH2d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb202_ixO(%rsp),%xmm0
    subps nb202_iyO(%rsp),%xmm1
    subps nb202_izO(%rsp),%xmm2
    subps nb202_ixH1(%rsp),%xmm3
    subps nb202_iyH1(%rsp),%xmm4
    subps nb202_izH1(%rsp),%xmm5
    subps nb202_ixH2(%rsp),%xmm6
    subps nb202_iyH2(%rsp),%xmm7
    subps nb202_izH2(%rsp),%xmm8

        movaps %xmm0,nb202_dxOH2(%rsp)
        movaps %xmm1,nb202_dyOH2(%rsp)
        movaps %xmm2,nb202_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH1H2(%rsp)
        movaps %xmm4,nb202_dyH1H2(%rsp)
        movaps %xmm5,nb202_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb202_dxH2H2(%rsp)
        movaps %xmm7,nb202_dyH2H2(%rsp)
        movaps %xmm8,nb202_dzH2H2(%rsp)
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

        movaps  nb202_three(%rsp),%xmm9
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

        movaps  nb202_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH2
        mulps   %xmm4,%xmm10 ## rinvH1H2
    mulps   %xmm4,%xmm11 ## rinvH2H2

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb202_krf(%rsp),%xmm2
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
    movaps nb202_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb202_qqOH(%rsp),%xmm12
    movaps nb202_qqHH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm13,%xmm4
    mulps  %xmm13,%xmm7
    addps  %xmm8,%xmm5
    addps  %xmm15,%xmm2
    addps  %xmm5,%xmm2
    movaps %xmm2,nb202_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm9    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm10    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm9    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm10    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm11    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm12    ## jzH2c  -  -  -
    movss  32(%rdi,%rbx,4),%xmm2    ## jzH2b  -  -  -
    movss  32(%rdi,%rdx,4),%xmm3    ## jzH2d  -  -  -
    movlhps %xmm2,%xmm11 ## jzH2a  -  jzH2b  -
    movlhps %xmm3,%xmm12 ## jzH2c  -  jzH2d -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## xmm9: jxH2a jyH2a jxH2b jyH2b 
    ## xmm10: jxH2c jyH2c jxH2d jyH2d
    ## xmm11: jzH2a jzH2b jzH2c jzH2d

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb202_dxOH2(%rsp),%xmm0
        mulps nb202_dyOH2(%rsp),%xmm1
        mulps nb202_dzOH2(%rsp),%xmm2
        mulps nb202_dxH1H2(%rsp),%xmm3
        mulps nb202_dyH1H2(%rsp),%xmm4
        mulps nb202_dzH1H2(%rsp),%xmm5
        mulps nb202_dxH2H2(%rsp),%xmm6
        mulps nb202_dyH2H2(%rsp),%xmm7
        mulps nb202_dzH2H2(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb202_fixO(%rsp),%xmm0
    addps nb202_fiyO(%rsp),%xmm1
    addps nb202_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb202_fixH1(%rsp),%xmm3
    addps nb202_fiyH1(%rsp),%xmm4
    addps nb202_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb202_fixH2(%rsp),%xmm6
    addps nb202_fiyH2(%rsp),%xmm7
    addps nb202_fizH2(%rsp),%xmm8

    movaps %xmm0,nb202_fixO(%rsp)
    movaps %xmm1,nb202_fiyO(%rsp)
    movaps %xmm2,nb202_fizO(%rsp)
    movaps %xmm3,nb202_fixH1(%rsp)
    movaps %xmm4,nb202_fiyH1(%rsp)
    movaps %xmm5,nb202_fizH1(%rsp)
    movaps %xmm6,nb202_fixH2(%rsp)
    movaps %xmm7,nb202_fiyH2(%rsp)
    movaps %xmm8,nb202_fizH2(%rsp)

    ## xmm0 = fH2x
    ## xmm1 = fH2y
    ## xmm2 = fH2z
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12 ## fH2zc fH2zd

    movlps %xmm9,24(%rdi,%rax,4)
    movhps %xmm9,24(%rdi,%rbx,4)
    movlps %xmm10,24(%rdi,%rcx,4)
    movhps %xmm10,24(%rdi,%rdx,4)
    movss  %xmm11,32(%rdi,%rax,4)
    movss  %xmm12,32(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,32(%rdi,%rbx,4)
    movss  %xmm12,32(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb202_innerk(%rsp)
        jl    _nb_kernel202_x86_64_sse.nb202_single_check
        jmp   _nb_kernel202_x86_64_sse.nb202_unroll_loop
_nb_kernel202_x86_64_sse.nb202_single_check: 
        addl $4,nb202_innerk(%rsp)
        jnz   _nb_kernel202_x86_64_sse.nb202_single_loop
        jmp   _nb_kernel202_x86_64_sse.nb202_updateouterdata
_nb_kernel202_x86_64_sse.nb202_single_loop: 
        movq  nb202_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb202_innerjjnr(%rsp)

        movq nb202_pos(%rbp),%rsi
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

        movaps %xmm0,nb202_jxO(%rsp)
        movaps %xmm1,nb202_jyO(%rsp)
        movaps %xmm2,nb202_jzO(%rsp)
        subps  nb202_ixO(%rsp),%xmm0
        subps  nb202_iyO(%rsp),%xmm1
        subps  nb202_izO(%rsp),%xmm2
        movaps %xmm0,nb202_dxOO(%rsp)
        movaps %xmm1,nb202_dyOO(%rsp)
        movaps %xmm2,nb202_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb202_krf(%rsp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        movaps  %xmm6,%xmm7    ## xmm7=krsq 
        mulps   %xmm1,%xmm1
        movaps  nb202_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202_half(%rsp),%xmm3   ## rinv iO - j water 



        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        mulps   nb202_two(%rsp),%xmm7
        subps   nb202_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb202_qqOO(%rsp),%xmm4
        movhps  nb202_qqOH(%rsp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  
        mulps %xmm4,%xmm3       ## coul part of fs  


        addps   nb202_vctot(%rsp),%xmm6
        mulps   %xmm3,%xmm0     ## total fscal 
        movaps  %xmm6,nb202_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb202_dxOO(%rsp),%xmm0
        mulps   nb202_dyOO(%rsp),%xmm1
        mulps   nb202_dzOO(%rsp),%xmm2

        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb202_fjxO(%rsp)
        movaps  %xmm4,nb202_fjyO(%rsp)
        movaps  %xmm5,nb202_fjzO(%rsp)
        addps   nb202_fixO(%rsp),%xmm0
        addps   nb202_fiyO(%rsp),%xmm1
        addps   nb202_fizO(%rsp),%xmm2
        movaps  %xmm0,nb202_fixO(%rsp)
        movaps  %xmm1,nb202_fiyO(%rsp)
        movaps  %xmm2,nb202_fizO(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb202_jxO(%rsp),%xmm0
    movaps  nb202_jyO(%rsp),%xmm1
    movaps  nb202_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
    subps   nb202_ixH1(%rsp),%xmm0
    subps   nb202_iyH1(%rsp),%xmm1
    subps   nb202_izH1(%rsp),%xmm2
    subps   nb202_ixH2(%rsp),%xmm3
    subps   nb202_iyH2(%rsp),%xmm4
    subps   nb202_izH2(%rsp),%xmm5

        movaps %xmm0,nb202_dxH1O(%rsp)
        movaps %xmm1,nb202_dyH1O(%rsp)
        movaps %xmm2,nb202_dzH1O(%rsp)
        movaps %xmm3,nb202_dxH2O(%rsp)
        movaps %xmm4,nb202_dyH2O(%rsp)
        movaps %xmm5,nb202_dzH2O(%rsp)
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
        movaps  nb202_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb202_half(%rsp),%xmm7   ## rinv H2 - j water  

        mulps nb202_krf(%rsp),%xmm0   ## krsq 
        mulps nb202_krf(%rsp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb202_qqOH(%rsp),%xmm6
        movhps  nb202_qqHH(%rsp),%xmm6
        movaps  %xmm0,%xmm1
        movaps  %xmm4,%xmm5
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb202_crf(%rsp),%xmm0
        subps   nb202_crf(%rsp),%xmm4
        mulps   nb202_two(%rsp),%xmm1
        mulps   nb202_two(%rsp),%xmm5
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb202_vctot(%rsp),%xmm4
        movaps  %xmm4,nb202_vctot(%rsp)
        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm4
        mulps   %xmm3,%xmm3
        mulps   %xmm7,%xmm7
        subps   %xmm1,%xmm0
        subps   %xmm5,%xmm4
        mulps   %xmm6,%xmm0
        mulps   %xmm6,%xmm4
        mulps   %xmm3,%xmm0     ## fscal 
        mulps   %xmm4,%xmm7     ## fscal 

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb202_dxH1O(%rsp),%xmm0
        mulps   nb202_dyH1O(%rsp),%xmm1
        mulps   nb202_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb202_fjxO(%rsp),%xmm3
        movaps  nb202_fjyO(%rsp),%xmm4
        movaps  nb202_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb202_fjxO(%rsp)
        movaps  %xmm4,nb202_fjyO(%rsp)
        movaps  %xmm5,nb202_fjzO(%rsp)
        addps   nb202_fixH1(%rsp),%xmm0
        addps   nb202_fiyH1(%rsp),%xmm1
        addps   nb202_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb202_fixH1(%rsp)
        movaps  %xmm1,nb202_fiyH1(%rsp)
        movaps  %xmm2,nb202_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb202_dxH2O(%rsp),%xmm0
        mulps   nb202_dyH2O(%rsp),%xmm1
        mulps   nb202_dzH2O(%rsp),%xmm2
        movaps  nb202_fjxO(%rsp),%xmm3
        movaps  nb202_fjyO(%rsp),%xmm4
        movaps  nb202_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb202_faction(%rbp),%rsi
        movaps  %xmm3,nb202_fjxO(%rsp)
        movaps  %xmm4,nb202_fjyO(%rsp)
        movaps  %xmm5,nb202_fjzO(%rsp)
        addps   nb202_fixH2(%rsp),%xmm0
        addps   nb202_fiyH2(%rsp),%xmm1
        addps   nb202_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb202_fixH2(%rsp)
        movaps  %xmm1,nb202_fiyH2(%rsp)
        movaps  %xmm2,nb202_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb202_fjxO(%rsp),%xmm3
        movaps  nb202_fjyO(%rsp),%xmm4
        movaps  nb202_fjzO(%rsp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## 00000010
        shufps $3,%xmm7,%xmm7 # # 00000011
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

        decl nb202_innerk(%rsp)
        jz    _nb_kernel202_x86_64_sse.nb202_updateouterdata
        jmp   _nb_kernel202_x86_64_sse.nb202_single_loop
_nb_kernel202_x86_64_sse.nb202_updateouterdata: 
        movl  nb202_ii3(%rsp),%ecx
        movq  nb202_faction(%rbp),%rdi
        movq  nb202_fshift(%rbp),%rsi
        movl  nb202_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb202_fixO(%rsp),%xmm0
        movaps nb202_fiyO(%rsp),%xmm1
        movaps nb202_fizO(%rsp),%xmm2

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
        movaps nb202_fixH1(%rsp),%xmm0
        movaps nb202_fiyH1(%rsp),%xmm1
        movaps nb202_fizH1(%rsp),%xmm2

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
        movaps nb202_fixH2(%rsp),%xmm0
        movaps nb202_fiyH2(%rsp),%xmm1
        movaps nb202_fizH2(%rsp),%xmm2

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
        movl nb202_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb202_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb202_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb202_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb202_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202_x86_64_sse.nb202_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202_n(%rsp)
        jmp _nb_kernel202_x86_64_sse.nb202_outer
_nb_kernel202_x86_64_sse.nb202_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202_x86_64_sse.nb202_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202_x86_64_sse.nb202_threadloop
_nb_kernel202_x86_64_sse.nb202_end: 

        movl nb202_nouter(%rsp),%eax
        movl nb202_ninner(%rsp),%ebx
        movq nb202_outeriter(%rbp),%rcx
        movq nb202_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1544,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel202nf_x86_64_sse
.globl _nb_kernel202nf_x86_64_sse
nb_kernel202nf_x86_64_sse:      
_nb_kernel202nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb202nf_fshift, 16
.set nb202nf_gid, 24
.set nb202nf_pos, 32
.set nb202nf_faction, 40
.set nb202nf_charge, 48
.set nb202nf_p_facel, 56
.set nb202nf_argkrf, 64
.set nb202nf_argcrf, 72
.set nb202nf_Vc, 80
.set nb202nf_type, 88
.set nb202nf_p_ntype, 96
.set nb202nf_vdwparam, 104
.set nb202nf_Vvdw, 112
.set nb202nf_p_tabscale, 120
.set nb202nf_VFtab, 128
.set nb202nf_invsqrta, 136
.set nb202nf_dvda, 144
.set nb202nf_p_gbtabscale, 152
.set nb202nf_GBtab, 160
.set nb202nf_p_nthreads, 168
.set nb202nf_count, 176
.set nb202nf_mtx, 184
.set nb202nf_outeriter, 192
.set nb202nf_inneriter, 200
.set nb202nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb202nf_ixO, 0
.set nb202nf_iyO, 16
.set nb202nf_izO, 32
.set nb202nf_ixH1, 48
.set nb202nf_iyH1, 64
.set nb202nf_izH1, 80
.set nb202nf_ixH2, 96
.set nb202nf_iyH2, 112
.set nb202nf_izH2, 128
.set nb202nf_jxO, 144
.set nb202nf_jyO, 160
.set nb202nf_jzO, 176
.set nb202nf_jxH1, 192
.set nb202nf_jyH1, 208
.set nb202nf_jzH1, 224
.set nb202nf_jxH2, 240
.set nb202nf_jyH2, 256
.set nb202nf_jzH2, 272
.set nb202nf_qqOO, 288
.set nb202nf_qqOH, 304
.set nb202nf_qqHH, 320
.set nb202nf_vctot, 336
.set nb202nf_half, 352
.set nb202nf_three, 368
.set nb202nf_rsqOO, 384
.set nb202nf_rsqOH1, 400
.set nb202nf_rsqOH2, 416
.set nb202nf_rsqH1O, 432
.set nb202nf_rsqH1H1, 448
.set nb202nf_rsqH1H2, 464
.set nb202nf_rsqH2O, 480
.set nb202nf_rsqH2H1, 496
.set nb202nf_rsqH2H2, 512
.set nb202nf_rinvOO, 528
.set nb202nf_rinvOH1, 544
.set nb202nf_rinvOH2, 560
.set nb202nf_rinvH1O, 576
.set nb202nf_rinvH1H1, 592
.set nb202nf_rinvH1H2, 608
.set nb202nf_rinvH2O, 624
.set nb202nf_rinvH2H1, 640
.set nb202nf_rinvH2H2, 656
.set nb202nf_krf, 672
.set nb202nf_crf, 688
.set nb202nf_is3, 704
.set nb202nf_ii3, 708
.set nb202nf_innerjjnr, 712
.set nb202nf_nri, 720
.set nb202nf_iinr, 728
.set nb202nf_jindex, 736
.set nb202nf_jjnr, 744
.set nb202nf_shift, 752
.set nb202nf_shiftvec, 760
.set nb202nf_facel, 768
.set nb202nf_innerk, 776
.set nb202nf_n, 780
.set nb202nf_nn1, 784
.set nb202nf_nouter, 788
.set nb202nf_ninner, 792

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms
        subq $808,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb202nf_nouter(%rsp)
        movl %eax,nb202nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb202nf_nri(%rsp)
        movq %rsi,nb202nf_iinr(%rsp)
        movq %rdx,nb202nf_jindex(%rsp)
        movq %rcx,nb202nf_jjnr(%rsp)
        movq %r8,nb202nf_shift(%rsp)
        movq %r9,nb202nf_shiftvec(%rsp)
        movq nb202nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb202nf_facel(%rsp)

        movq nb202nf_argkrf(%rbp),%rsi
        movq nb202nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb202nf_krf(%rsp)
        movaps %xmm2,nb202nf_crf(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  nb202nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb202nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb202nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb202nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202nf_qqOO(%rsp)
        movaps %xmm4,nb202nf_qqOH(%rsp)
        movaps %xmm5,nb202nf_qqHH(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb202nf_half(%rsp)
        movss nb202nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb202nf_half(%rsp)
        movaps %xmm3,nb202nf_three(%rsp)

_nb_kernel202nf_x86_64_sse.nb202nf_threadloop: 
        movq  nb202nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel202nf_x86_64_sse.nb202nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202nf_x86_64_sse.nb202nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb202nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb202nf_n(%rsp)
        movl %ebx,nb202nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel202nf_x86_64_sse.nb202nf_outerstart
        jmp _nb_kernel202nf_x86_64_sse.nb202nf_end

_nb_kernel202nf_x86_64_sse.nb202nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202nf_nouter(%rsp),%ebx
        movl %ebx,nb202nf_nouter(%rsp)

_nb_kernel202nf_x86_64_sse.nb202nf_outer: 
        movq  nb202nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb202nf_is3(%rsp)            ## store is3 

        movq  nb202nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb202nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb202nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb202nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202nf_ixO(%rsp)
        movaps %xmm4,nb202nf_iyO(%rsp)
        movaps %xmm5,nb202nf_izO(%rsp)

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
        movaps %xmm0,nb202nf_ixH1(%rsp)
        movaps %xmm1,nb202nf_iyH1(%rsp)
        movaps %xmm2,nb202nf_izH1(%rsp)
        movaps %xmm3,nb202nf_ixH2(%rsp)
        movaps %xmm4,nb202nf_iyH2(%rsp)
        movaps %xmm5,nb202nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb202nf_vctot(%rsp)

        movq  nb202nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb202nf_pos(%rbp),%rsi
        movq  nb202nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb202nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb202nf_ninner(%rsp),%ecx
        movl  %ecx,nb202nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb202nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel202nf_x86_64_sse.nb202nf_unroll_loop
        jmp   _nb_kernel202nf_x86_64_sse.nb202nf_single_check
_nb_kernel202nf_x86_64_sse.nb202nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb202nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb202nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb202nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps %xmm0,nb202nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb202nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb202nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb202nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb202nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb202nf_jyH2(%rsp)

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
        movaps %xmm0,nb202nf_jzO(%rsp)
        movaps %xmm1,nb202nf_jzH1(%rsp)
        movaps %xmm2,nb202nf_jzH2(%rsp)

        movaps nb202nf_ixO(%rsp),%xmm0
        movaps nb202nf_iyO(%rsp),%xmm1
        movaps nb202nf_izO(%rsp),%xmm2
        movaps nb202nf_ixO(%rsp),%xmm3
        movaps nb202nf_iyO(%rsp),%xmm4
        movaps nb202nf_izO(%rsp),%xmm5
        subps  nb202nf_jxO(%rsp),%xmm0
        subps  nb202nf_jyO(%rsp),%xmm1
        subps  nb202nf_jzO(%rsp),%xmm2
        subps  nb202nf_jxH1(%rsp),%xmm3
        subps  nb202nf_jyH1(%rsp),%xmm4
        subps  nb202nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb202nf_rsqOO(%rsp)
        movaps %xmm3,nb202nf_rsqOH1(%rsp)

        movaps nb202nf_ixO(%rsp),%xmm0
        movaps nb202nf_iyO(%rsp),%xmm1
        movaps nb202nf_izO(%rsp),%xmm2
        movaps nb202nf_ixH1(%rsp),%xmm3
        movaps nb202nf_iyH1(%rsp),%xmm4
        movaps nb202nf_izH1(%rsp),%xmm5
        subps  nb202nf_jxH2(%rsp),%xmm0
        subps  nb202nf_jyH2(%rsp),%xmm1
        subps  nb202nf_jzH2(%rsp),%xmm2
        subps  nb202nf_jxO(%rsp),%xmm3
        subps  nb202nf_jyO(%rsp),%xmm4
        subps  nb202nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb202nf_rsqOH2(%rsp)
        movaps %xmm3,nb202nf_rsqH1O(%rsp)

        movaps nb202nf_ixH1(%rsp),%xmm0
        movaps nb202nf_iyH1(%rsp),%xmm1
        movaps nb202nf_izH1(%rsp),%xmm2
        movaps nb202nf_ixH1(%rsp),%xmm3
        movaps nb202nf_iyH1(%rsp),%xmm4
        movaps nb202nf_izH1(%rsp),%xmm5
        subps  nb202nf_jxH1(%rsp),%xmm0
        subps  nb202nf_jyH1(%rsp),%xmm1
        subps  nb202nf_jzH1(%rsp),%xmm2
        subps  nb202nf_jxH2(%rsp),%xmm3
        subps  nb202nf_jyH2(%rsp),%xmm4
        subps  nb202nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb202nf_rsqH1H1(%rsp)
        movaps %xmm3,nb202nf_rsqH1H2(%rsp)

        movaps nb202nf_ixH2(%rsp),%xmm0
        movaps nb202nf_iyH2(%rsp),%xmm1
        movaps nb202nf_izH2(%rsp),%xmm2
        movaps nb202nf_ixH2(%rsp),%xmm3
        movaps nb202nf_iyH2(%rsp),%xmm4
        movaps nb202nf_izH2(%rsp),%xmm5
        subps  nb202nf_jxO(%rsp),%xmm0
        subps  nb202nf_jyO(%rsp),%xmm1
        subps  nb202nf_jzO(%rsp),%xmm2
        subps  nb202nf_jxH1(%rsp),%xmm3
        subps  nb202nf_jyH1(%rsp),%xmm4
        subps  nb202nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb202nf_rsqH2O(%rsp)
        movaps %xmm4,nb202nf_rsqH2H1(%rsp)

        movaps nb202nf_ixH2(%rsp),%xmm0
        movaps nb202nf_iyH2(%rsp),%xmm1
        movaps nb202nf_izH2(%rsp),%xmm2
        subps  nb202nf_jxH2(%rsp),%xmm0
        subps  nb202nf_jyH2(%rsp),%xmm1
        subps  nb202nf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb202nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb202nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb202nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb202nf_rinvH2H1(%rsp)

        rsqrtps nb202nf_rsqOO(%rsp),%xmm1
        rsqrtps nb202nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqOO(%rsp),%xmm1
        mulps   nb202nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%rsp),%xmm3
        mulps   nb202nf_half(%rsp),%xmm7
        movaps  %xmm3,nb202nf_rinvOO(%rsp)
        movaps  %xmm7,nb202nf_rinvOH1(%rsp)

        rsqrtps nb202nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb202nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqOH2(%rsp),%xmm1
        mulps   nb202nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%rsp),%xmm3
        mulps   nb202nf_half(%rsp),%xmm7
        movaps  %xmm3,nb202nf_rinvOH2(%rsp)
        movaps  %xmm7,nb202nf_rinvH1O(%rsp)

        rsqrtps nb202nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb202nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqH1H1(%rsp),%xmm1
        mulps   nb202nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%rsp),%xmm3
        mulps   nb202nf_half(%rsp),%xmm7
        movaps  %xmm3,nb202nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb202nf_rinvH1H2(%rsp)

        rsqrtps nb202nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb202nf_three(%rsp),%xmm3
        mulps   nb202nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202nf_half(%rsp),%xmm3
        movaps  %xmm3,nb202nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movaps nb202nf_krf(%rsp),%xmm6
        mulps  nb202nf_rsqOO(%rsp),%xmm6   ## xmm5=krsq 
        addps  nb202nf_rinvOO(%rsp),%xmm6       ## xmm6=rinv+ krsq 
        subps  nb202nf_crf(%rsp),%xmm6
        mulps  nb202nf_qqOO(%rsp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addps  nb202nf_vctot(%rsp),%xmm6   ## local vctot summation variable 

        ## O-H interactions 
        movaps nb202nf_krf(%rsp),%xmm0
        movaps nb202nf_krf(%rsp),%xmm1
        movaps nb202nf_krf(%rsp),%xmm2
        movaps nb202nf_krf(%rsp),%xmm3
        mulps  nb202nf_rsqOH1(%rsp),%xmm0   ## krsq 
        mulps  nb202nf_rsqOH2(%rsp),%xmm1   ## krsq 
        mulps  nb202nf_rsqH1O(%rsp),%xmm2   ## krsq 
        mulps  nb202nf_rsqH2O(%rsp),%xmm3   ## krsq 
        addps  nb202nf_rinvOH1(%rsp),%xmm0      ## rinv+ krsq 
        addps  nb202nf_rinvOH2(%rsp),%xmm1      ## rinv+ krsq 
        addps  nb202nf_rinvH1O(%rsp),%xmm2      ## rinv+ krsq 
        addps  nb202nf_rinvH2O(%rsp),%xmm3      ## rinv+ krsq 
        subps  nb202nf_crf(%rsp),%xmm0
        subps  nb202nf_crf(%rsp),%xmm1
        subps  nb202nf_crf(%rsp),%xmm2
        subps  nb202nf_crf(%rsp),%xmm3
        mulps  nb202nf_qqOH(%rsp),%xmm0   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%rsp),%xmm1   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%rsp),%xmm3   ## voul=qq*(rinv+ krsq-crf) 
        addps %xmm0,%xmm6
        addps %xmm2,%xmm1
        addps %xmm3,%xmm6
        addps %xmm1,%xmm6

        ## H-H interactions 
        movaps nb202nf_krf(%rsp),%xmm0
        movaps nb202nf_krf(%rsp),%xmm1
        movaps nb202nf_krf(%rsp),%xmm2
        movaps nb202nf_krf(%rsp),%xmm3
        mulps  nb202nf_rsqH1H1(%rsp),%xmm0   ## krsq 
        mulps  nb202nf_rsqH1H2(%rsp),%xmm1   ## krsq 
        mulps  nb202nf_rsqH2H1(%rsp),%xmm2   ## krsq 
        mulps  nb202nf_rsqH2H2(%rsp),%xmm3   ## krsq 
        addps  nb202nf_rinvH1H1(%rsp),%xmm0     ## rinv+ krsq 
        addps  nb202nf_rinvH1H2(%rsp),%xmm1     ## rinv+ krsq 
        addps  nb202nf_rinvH2H1(%rsp),%xmm2     ## rinv+ krsq 
        addps  nb202nf_rinvH2H2(%rsp),%xmm3     ## rinv+ krsq 
        subps  nb202nf_crf(%rsp),%xmm0
        subps  nb202nf_crf(%rsp),%xmm1
        subps  nb202nf_crf(%rsp),%xmm2
        subps  nb202nf_crf(%rsp),%xmm3
        mulps  nb202nf_qqHH(%rsp),%xmm0   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%rsp),%xmm1   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%rsp),%xmm3   ## voul=qq*(rinv+ krsq-crf) 
        addps %xmm0,%xmm6
        addps %xmm2,%xmm1
        addps %xmm3,%xmm6
        addps %xmm1,%xmm6
        movaps %xmm6,nb202nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb202nf_innerk(%rsp)
        jl    _nb_kernel202nf_x86_64_sse.nb202nf_single_check
        jmp   _nb_kernel202nf_x86_64_sse.nb202nf_unroll_loop
_nb_kernel202nf_x86_64_sse.nb202nf_single_check: 
        addl $4,nb202nf_innerk(%rsp)
        jnz   _nb_kernel202nf_x86_64_sse.nb202nf_single_loop
        jmp   _nb_kernel202nf_x86_64_sse.nb202nf_updateouterdata
_nb_kernel202nf_x86_64_sse.nb202nf_single_loop: 
        movq  nb202nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb202nf_innerjjnr(%rsp)

        movq nb202nf_pos(%rbp),%rsi
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
        movaps  nb202nf_ixO(%rsp),%xmm0
        movaps  nb202nf_iyO(%rsp),%xmm1
        movaps  nb202nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb202nf_jxO(%rsp)
        movaps %xmm4,nb202nf_jyO(%rsp)
        movaps %xmm5,nb202nf_jzO(%rsp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb202nf_krf(%rsp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb202nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202nf_half(%rsp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb202nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb202nf_qqOO(%rsp),%xmm4
        movhps  nb202nf_qqOH(%rsp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  

        addps   nb202nf_vctot(%rsp),%xmm6
        movaps  %xmm6,nb202nf_vctot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb202nf_ixH1(%rsp),%xmm0
        movaps  nb202nf_iyH1(%rsp),%xmm1
        movaps  nb202nf_izH1(%rsp),%xmm2
        movaps  nb202nf_ixH2(%rsp),%xmm3
        movaps  nb202nf_iyH2(%rsp),%xmm4
        movaps  nb202nf_izH2(%rsp),%xmm5
        subps   nb202nf_jxO(%rsp),%xmm0
        subps   nb202nf_jyO(%rsp),%xmm1
        subps   nb202nf_jzO(%rsp),%xmm2
        subps   nb202nf_jxO(%rsp),%xmm3
        subps   nb202nf_jyO(%rsp),%xmm4
        subps   nb202nf_jzO(%rsp),%xmm5
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
        movaps  nb202nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb202nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        mulps nb202nf_krf(%rsp),%xmm0   ## krsq 
        mulps nb202nf_krf(%rsp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb202nf_qqOH(%rsp),%xmm6
        movhps  nb202nf_qqHH(%rsp),%xmm6
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb202nf_crf(%rsp),%xmm0
        subps   nb202nf_crf(%rsp),%xmm4
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb202nf_vctot(%rsp),%xmm4
        movaps  %xmm4,nb202nf_vctot(%rsp)

        decl nb202nf_innerk(%rsp)
        jz    _nb_kernel202nf_x86_64_sse.nb202nf_updateouterdata
        jmp   _nb_kernel202nf_x86_64_sse.nb202nf_single_loop
_nb_kernel202nf_x86_64_sse.nb202nf_updateouterdata: 
        ## get n from stack
        movl nb202nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb202nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb202nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb202nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb202nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202nf_x86_64_sse.nb202nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202nf_n(%rsp)
        jmp _nb_kernel202nf_x86_64_sse.nb202nf_outer
_nb_kernel202nf_x86_64_sse.nb202nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202nf_x86_64_sse.nb202nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202nf_x86_64_sse.nb202nf_threadloop
_nb_kernel202nf_x86_64_sse.nb202nf_end: 

        movl nb202nf_nouter(%rsp),%eax
        movl nb202nf_ninner(%rsp),%ebx
        movq nb202nf_outeriter(%rbp),%rcx
        movq nb202nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $808,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret




