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







.globl nb_kernel102_x86_64_sse
.globl _nb_kernel102_x86_64_sse
nb_kernel102_x86_64_sse:        
_nb_kernel102_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb102_fshift, 16
.set nb102_gid, 24
.set nb102_pos, 32
.set nb102_faction, 40
.set nb102_charge, 48
.set nb102_p_facel, 56
.set nb102_argkrf, 64
.set nb102_argcrf, 72
.set nb102_Vc, 80
.set nb102_type, 88
.set nb102_p_ntype, 96
.set nb102_vdwparam, 104
.set nb102_Vvdw, 112
.set nb102_p_tabscale, 120
.set nb102_VFtab, 128
.set nb102_invsqrta, 136
.set nb102_dvda, 144
.set nb102_p_gbtabscale, 152
.set nb102_GBtab, 160
.set nb102_p_nthreads, 168
.set nb102_count, 176
.set nb102_mtx, 184
.set nb102_outeriter, 192
.set nb102_inneriter, 200
.set nb102_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb102_ixO, 0
.set nb102_iyO, 16
.set nb102_izO, 32
.set nb102_ixH1, 48
.set nb102_iyH1, 64
.set nb102_izH1, 80
.set nb102_ixH2, 96
.set nb102_iyH2, 112
.set nb102_izH2, 128
.set nb102_jxO, 144
.set nb102_jyO, 160
.set nb102_jzO, 176
.set nb102_jxH1, 192
.set nb102_jyH1, 208
.set nb102_jzH1, 224
.set nb102_jxH2, 240
.set nb102_jyH2, 256
.set nb102_jzH2, 272
.set nb102_dxOO, 288
.set nb102_dyOO, 304
.set nb102_dzOO, 320
.set nb102_dxOH1, 336
.set nb102_dyOH1, 352
.set nb102_dzOH1, 368
.set nb102_dxOH2, 384
.set nb102_dyOH2, 400
.set nb102_dzOH2, 416
.set nb102_dxH1O, 432
.set nb102_dyH1O, 448
.set nb102_dzH1O, 464
.set nb102_dxH1H1, 480
.set nb102_dyH1H1, 496
.set nb102_dzH1H1, 512
.set nb102_dxH1H2, 528
.set nb102_dyH1H2, 544
.set nb102_dzH1H2, 560
.set nb102_dxH2O, 576
.set nb102_dyH2O, 592
.set nb102_dzH2O, 608
.set nb102_dxH2H1, 624
.set nb102_dyH2H1, 640
.set nb102_dzH2H1, 656
.set nb102_dxH2H2, 672
.set nb102_dyH2H2, 688
.set nb102_dzH2H2, 704
.set nb102_qqOO, 720
.set nb102_qqOH, 736
.set nb102_qqHH, 752
.set nb102_vctot, 768
.set nb102_fixO, 784
.set nb102_fiyO, 800
.set nb102_fizO, 816
.set nb102_fixH1, 832
.set nb102_fiyH1, 848
.set nb102_fizH1, 864
.set nb102_fixH2, 880
.set nb102_fiyH2, 896
.set nb102_fizH2, 912
.set nb102_fjxO, 928
.set nb102_fjyO, 944
.set nb102_fjzO, 960
.set nb102_fjxH1, 976
.set nb102_fjyH1, 992
.set nb102_fjzH1, 1008
.set nb102_fjxH2, 1024
.set nb102_fjyH2, 1040
.set nb102_fjzH2, 1056
.set nb102_half, 1072
.set nb102_three, 1088
.set nb102_rsqOO, 1104
.set nb102_rsqOH1, 1120
.set nb102_rsqOH2, 1136
.set nb102_rsqH1O, 1152
.set nb102_rsqH1H1, 1168
.set nb102_rsqH1H2, 1184
.set nb102_rsqH2O, 1200
.set nb102_rsqH2H1, 1216
.set nb102_rsqH2H2, 1232
.set nb102_rinvOO, 1248
.set nb102_rinvOH1, 1264
.set nb102_rinvOH2, 1280
.set nb102_rinvH1O, 1296
.set nb102_rinvH1H1, 1312
.set nb102_rinvH1H2, 1328
.set nb102_rinvH2O, 1344
.set nb102_rinvH2H1, 1360
.set nb102_rinvH2H2, 1376
.set nb102_is3, 1392
.set nb102_ii3, 1396
.set nb102_innerjjnr, 1400
.set nb102_nri, 1408
.set nb102_iinr, 1416
.set nb102_jindex, 1424
.set nb102_jjnr, 1432
.set nb102_shift, 1440
.set nb102_shiftvec, 1448
.set nb102_facel, 1456
.set nb102_innerk, 1464
.set nb102_n, 1472
.set nb102_nn1, 1480
.set nb102_nouter, 1484
.set nb102_ninner, 1488
.set nb102_salign, 1492

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
        subq $1512,%rsp
        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb102_nouter(%rsp)
        movl %eax,nb102_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb102_nri(%rsp)
        movq %rsi,nb102_iinr(%rsp)
        movq %rdx,nb102_jindex(%rsp)
        movq %rcx,nb102_jjnr(%rsp)
        movq %r8,nb102_shift(%rsp)
        movq %r9,nb102_shiftvec(%rsp)
        movq nb102_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb102_facel(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb102_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb102_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb102_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb102_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102_qqOO(%rsp)
        movaps %xmm4,nb102_qqOH(%rsp)
        movaps %xmm5,nb102_qqHH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb102_half(%rsp)
        movss nb102_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb102_half(%rsp)
        movaps %xmm3,nb102_three(%rsp)
        movq nb102_faction(%rbp),%rdi


_nb_kernel102_x86_64_sse.nb102_threadloop: 
        movq  nb102_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel102_x86_64_sse.nb102_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102_x86_64_sse.nb102_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102_n(%rsp)
        movl %ebx,nb102_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel102_x86_64_sse.nb102_outerstart
        jmp _nb_kernel102_x86_64_sse.nb102_end

_nb_kernel102_x86_64_sse.nb102_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102_nouter(%rsp),%ebx
        movl %ebx,nb102_nouter(%rsp)

_nb_kernel102_x86_64_sse.nb102_outer: 
        movq  nb102_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb102_is3(%rsp)      ## store is3 

        movq  nb102_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb102_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb102_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb102_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102_ixO(%rsp)
        movaps %xmm4,nb102_iyO(%rsp)
        movaps %xmm5,nb102_izO(%rsp)

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
        movaps %xmm0,nb102_ixH1(%rsp)
        movaps %xmm1,nb102_iyH1(%rsp)
        movaps %xmm2,nb102_izH1(%rsp)
        movaps %xmm3,nb102_ixH2(%rsp)
        movaps %xmm4,nb102_iyH2(%rsp)
        movaps %xmm5,nb102_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb102_vctot(%rsp)
        movaps %xmm4,nb102_fixO(%rsp)
        movaps %xmm4,nb102_fiyO(%rsp)
        movaps %xmm4,nb102_fizO(%rsp)
        movaps %xmm4,nb102_fixH1(%rsp)
        movaps %xmm4,nb102_fiyH1(%rsp)
        movaps %xmm4,nb102_fizH1(%rsp)
        movaps %xmm4,nb102_fixH2(%rsp)
        movaps %xmm4,nb102_fiyH2(%rsp)
        movaps %xmm4,nb102_fizH2(%rsp)

        movq  nb102_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movq  nb102_pos(%rbp),%rsi
        movq  nb102_faction(%rbp),%rdi
        movq  nb102_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb102_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb102_ninner(%rsp),%ecx
        movl  %ecx,nb102_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb102_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel102_x86_64_sse.nb102_unroll_loop
        jmp   _nb_kernel102_x86_64_sse.nb102_single_check
_nb_kernel102_x86_64_sse.nb102_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb102_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb102_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb102_pos(%rbp),%rsi        ## base of pos[] 

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


    subps nb102_ixO(%rsp),%xmm0
    subps nb102_iyO(%rsp),%xmm1
    subps nb102_izO(%rsp),%xmm2
    subps nb102_ixH1(%rsp),%xmm3
    subps nb102_iyH1(%rsp),%xmm4
    subps nb102_izH1(%rsp),%xmm5
    subps nb102_ixH2(%rsp),%xmm6
    subps nb102_iyH2(%rsp),%xmm7
    subps nb102_izH2(%rsp),%xmm8

        movaps %xmm0,nb102_dxOO(%rsp)
        movaps %xmm1,nb102_dyOO(%rsp)
        movaps %xmm2,nb102_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH1O(%rsp)
        movaps %xmm4,nb102_dyH1O(%rsp)
        movaps %xmm5,nb102_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb102_dxH2O(%rsp)
        movaps %xmm7,nb102_dyH2O(%rsp)
        movaps %xmm8,nb102_dzH2O(%rsp)
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

        movaps  nb102_three(%rsp),%xmm9
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

        movaps  nb102_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvOO 
        mulps   %xmm0,%xmm10 ## rinvH1O
    mulps   %xmm0,%xmm11 ## rinvH2O

        ## O interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb102_qqOO(%rsp),%xmm0
    mulps  nb102_qqOH(%rsp),%xmm1
    mulps  nb102_qqOH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb102_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb102_vctot(%rsp)

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

        mulps nb102_dxOO(%rsp),%xmm7
        mulps nb102_dyOO(%rsp),%xmm8
        mulps nb102_dzOO(%rsp),%xmm9
        mulps nb102_dxH1O(%rsp),%xmm10
        mulps nb102_dyH1O(%rsp),%xmm11
        mulps nb102_dzH1O(%rsp),%xmm12
        mulps nb102_dxH2O(%rsp),%xmm13
        mulps nb102_dyH2O(%rsp),%xmm14
        mulps nb102_dzH2O(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb102_fixO(%rsp),%xmm7
    addps nb102_fiyO(%rsp),%xmm8
    addps nb102_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb102_fixH1(%rsp),%xmm10
    addps nb102_fiyH1(%rsp),%xmm11
    addps nb102_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb102_fixH2(%rsp),%xmm13
    addps nb102_fiyH2(%rsp),%xmm14
    addps nb102_fizH2(%rsp),%xmm15

    movaps %xmm7,nb102_fixO(%rsp)
    movaps %xmm8,nb102_fiyO(%rsp)
    movaps %xmm9,nb102_fizO(%rsp)
    movaps %xmm10,nb102_fixH1(%rsp)
    movaps %xmm11,nb102_fiyH1(%rsp)
    movaps %xmm12,nb102_fizH1(%rsp)
    movaps %xmm13,nb102_fixH2(%rsp)
    movaps %xmm14,nb102_fiyH2(%rsp)
    movaps %xmm15,nb102_fizH2(%rsp)

    ## xmm3 = fOx , xmm4 = fOy
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

    subps nb102_ixO(%rsp),%xmm0
    subps nb102_iyO(%rsp),%xmm1
    subps nb102_izO(%rsp),%xmm2
    subps nb102_ixH1(%rsp),%xmm3
    subps nb102_iyH1(%rsp),%xmm4
    subps nb102_izH1(%rsp),%xmm5
    subps nb102_ixH2(%rsp),%xmm6
    subps nb102_iyH2(%rsp),%xmm7
    subps nb102_izH2(%rsp),%xmm8

        movaps %xmm0,nb102_dxOH1(%rsp)
        movaps %xmm1,nb102_dyOH1(%rsp)
        movaps %xmm2,nb102_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH1H1(%rsp)
        movaps %xmm4,nb102_dyH1H1(%rsp)
        movaps %xmm5,nb102_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb102_dxH2H1(%rsp)
        movaps %xmm7,nb102_dyH2H1(%rsp)
        movaps %xmm8,nb102_dzH2H1(%rsp)
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

        movaps  nb102_three(%rsp),%xmm9
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

        movaps  nb102_half(%rsp),%xmm0
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
    mulps  nb102_qqOH(%rsp),%xmm0
    mulps  nb102_qqHH(%rsp),%xmm1
    mulps  nb102_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb102_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb102_vctot(%rsp)

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

        mulps nb102_dxOH1(%rsp),%xmm7
        mulps nb102_dyOH1(%rsp),%xmm8
        mulps nb102_dzOH1(%rsp),%xmm9
        mulps nb102_dxH1H1(%rsp),%xmm10
        mulps nb102_dyH1H1(%rsp),%xmm11
        mulps nb102_dzH1H1(%rsp),%xmm12
        mulps nb102_dxH2H1(%rsp),%xmm13
        mulps nb102_dyH2H1(%rsp),%xmm14
        mulps nb102_dzH2H1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb102_fixO(%rsp),%xmm7
    addps nb102_fiyO(%rsp),%xmm8
    addps nb102_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb102_fixH1(%rsp),%xmm10
    addps nb102_fiyH1(%rsp),%xmm11
    addps nb102_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb102_fixH2(%rsp),%xmm13
    addps nb102_fiyH2(%rsp),%xmm14
    addps nb102_fizH2(%rsp),%xmm15

    movaps %xmm7,nb102_fixO(%rsp)
    movaps %xmm8,nb102_fiyO(%rsp)
    movaps %xmm9,nb102_fizO(%rsp)
    movaps %xmm10,nb102_fixH1(%rsp)
    movaps %xmm11,nb102_fiyH1(%rsp)
    movaps %xmm12,nb102_fizH1(%rsp)
    movaps %xmm13,nb102_fixH2(%rsp)
    movaps %xmm14,nb102_fiyH2(%rsp)
    movaps %xmm15,nb102_fizH2(%rsp)

    ## xmm3 = fH1x , xmm4 = fH1y
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

    subps nb102_ixO(%rsp),%xmm0
    subps nb102_iyO(%rsp),%xmm1
    subps nb102_izO(%rsp),%xmm2
    subps nb102_ixH1(%rsp),%xmm3
    subps nb102_iyH1(%rsp),%xmm4
    subps nb102_izH1(%rsp),%xmm5
    subps nb102_ixH2(%rsp),%xmm6
    subps nb102_iyH2(%rsp),%xmm7
    subps nb102_izH2(%rsp),%xmm8

        movaps %xmm0,nb102_dxOH2(%rsp)
        movaps %xmm1,nb102_dyOH2(%rsp)
        movaps %xmm2,nb102_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb102_dxH1H2(%rsp)
        movaps %xmm4,nb102_dyH1H2(%rsp)
        movaps %xmm5,nb102_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb102_dxH2H2(%rsp)
        movaps %xmm7,nb102_dyH2H2(%rsp)
        movaps %xmm8,nb102_dzH2H2(%rsp)
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

        movaps  nb102_three(%rsp),%xmm9
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

        movaps  nb102_half(%rsp),%xmm0
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
    mulps  nb102_qqOH(%rsp),%xmm0
    mulps  nb102_qqHH(%rsp),%xmm1
    mulps  nb102_qqHH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb102_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb102_vctot(%rsp)

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm2     ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm3     ## jzH2c  -  -  -
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

        mulps nb102_dxOH2(%rsp),%xmm7
        mulps nb102_dyOH2(%rsp),%xmm8
        mulps nb102_dzOH2(%rsp),%xmm9
        mulps nb102_dxH1H2(%rsp),%xmm10
        mulps nb102_dyH1H2(%rsp),%xmm11
        mulps nb102_dzH1H2(%rsp),%xmm12
        mulps nb102_dxH2H2(%rsp),%xmm13
        mulps nb102_dyH2H2(%rsp),%xmm14
        mulps nb102_dzH2H2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb102_fixO(%rsp),%xmm7
    addps nb102_fiyO(%rsp),%xmm8
    addps nb102_fizO(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb102_fixH1(%rsp),%xmm10
    addps nb102_fiyH1(%rsp),%xmm11
    addps nb102_fizH1(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb102_fixH2(%rsp),%xmm13
    addps nb102_fiyH2(%rsp),%xmm14
    addps nb102_fizH2(%rsp),%xmm15

    movaps %xmm7,nb102_fixO(%rsp)
    movaps %xmm8,nb102_fiyO(%rsp)
    movaps %xmm9,nb102_fizO(%rsp)
    movaps %xmm10,nb102_fixH1(%rsp)
    movaps %xmm11,nb102_fiyH1(%rsp)
    movaps %xmm12,nb102_fizH1(%rsp)
    movaps %xmm13,nb102_fixH2(%rsp)
    movaps %xmm14,nb102_fiyH2(%rsp)
    movaps %xmm15,nb102_fizH2(%rsp)

    ## xmm3 = fH2x , xmm4 = fH2y
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
        subl $4,nb102_innerk(%rsp)
        jl    _nb_kernel102_x86_64_sse.nb102_single_check
        jmp   _nb_kernel102_x86_64_sse.nb102_unroll_loop

_nb_kernel102_x86_64_sse.nb102_single_check: 
        addl $4,nb102_innerk(%rsp)
        jnz   _nb_kernel102_x86_64_sse.nb102_single_loop
        jmp   _nb_kernel102_x86_64_sse.nb102_updateouterdata
_nb_kernel102_x86_64_sse.nb102_single_loop: 
        movq  nb102_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax                      ## index is 32 bits 
        addq $4,nb102_innerjjnr(%rsp)

        movq nb102_pos(%rbp),%rsi
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

        movlhps %xmm6,%xmm0                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm1 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm2 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  

        movaps %xmm0,nb102_jxO(%rsp)
        movaps %xmm1,nb102_jyO(%rsp)
        movaps %xmm2,nb102_jzO(%rsp)
        subps  nb102_ixO(%rsp),%xmm0
        subps  nb102_iyO(%rsp),%xmm1
        subps  nb102_izO(%rsp),%xmm2
        movaps %xmm0,nb102_dxOO(%rsp)
        movaps %xmm1,nb102_dyOO(%rsp)
        movaps %xmm2,nb102_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb102_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb102_qqOO(%rsp),%xmm4

        movhps  nb102_qqOH(%rsp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulps   %xmm3,%xmm0     ## total fscal 
        addps   nb102_vctot(%rsp),%xmm3
        movaps  %xmm3,nb102_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb102_dxOO(%rsp),%xmm0
        mulps   nb102_dyOO(%rsp),%xmm1
        mulps   nb102_dzOO(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb102_fjxO(%rsp)
        movaps  %xmm4,nb102_fjyO(%rsp)
        movaps  %xmm5,nb102_fjzO(%rsp)
        addps   nb102_fixO(%rsp),%xmm0
        addps   nb102_fiyO(%rsp),%xmm1
        addps   nb102_fizO(%rsp),%xmm2
        movaps  %xmm0,nb102_fixO(%rsp)
        movaps  %xmm1,nb102_fiyO(%rsp)
        movaps  %xmm2,nb102_fizO(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb102_jxO(%rsp),%xmm0
    movaps  nb102_jyO(%rsp),%xmm1
    movaps  nb102_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  nb102_ixH1(%rsp),%xmm0
        subps  nb102_iyH1(%rsp),%xmm1
        subps  nb102_izH1(%rsp),%xmm2
        subps  nb102_ixH2(%rsp),%xmm3
        subps  nb102_iyH2(%rsp),%xmm4
        subps  nb102_izH2(%rsp),%xmm5

        movaps %xmm0,nb102_dxH1O(%rsp)
        movaps %xmm1,nb102_dyH1O(%rsp)
        movaps %xmm2,nb102_dzH1O(%rsp)
        movaps %xmm3,nb102_dxH2O(%rsp)
        movaps %xmm4,nb102_dyH2O(%rsp)
        movaps %xmm5,nb102_dzH2O(%rsp)
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
        movaps  nb102_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb102_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb102_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb102_qqHH(%rsp),%xmm6
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   nb102_vctot(%rsp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,nb102_vctot(%rsp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb102_dxH1O(%rsp),%xmm0
        mulps   nb102_dyH1O(%rsp),%xmm1
        mulps   nb102_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb102_fjxO(%rsp),%xmm3
        movaps  nb102_fjyO(%rsp),%xmm4
        movaps  nb102_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb102_fjxO(%rsp)
        movaps  %xmm4,nb102_fjyO(%rsp)
        movaps  %xmm5,nb102_fjzO(%rsp)
        addps   nb102_fixH1(%rsp),%xmm0
        addps   nb102_fiyH1(%rsp),%xmm1
        addps   nb102_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb102_fixH1(%rsp)
        movaps  %xmm1,nb102_fiyH1(%rsp)
        movaps  %xmm2,nb102_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb102_dxH2O(%rsp),%xmm0
        mulps   nb102_dyH2O(%rsp),%xmm1
        mulps   nb102_dzH2O(%rsp),%xmm2
        movaps  nb102_fjxO(%rsp),%xmm3
        movaps  nb102_fjyO(%rsp),%xmm4
        movaps  nb102_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb102_faction(%rbp),%rsi
        movaps  %xmm3,nb102_fjxO(%rsp)
        movaps  %xmm4,nb102_fjyO(%rsp)
        movaps  %xmm5,nb102_fjzO(%rsp)
        addps   nb102_fixH2(%rsp),%xmm0
        addps   nb102_fiyH2(%rsp),%xmm1
        addps   nb102_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb102_fixH2(%rsp)
        movaps  %xmm1,nb102_fiyH2(%rsp)
        movaps  %xmm2,nb102_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb102_fjxO(%rsp),%xmm3
        movaps  nb102_fjyO(%rsp),%xmm4
        movaps  nb102_fjzO(%rsp),%xmm5
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

        decl  nb102_innerk(%rsp)
        jz    _nb_kernel102_x86_64_sse.nb102_updateouterdata
        jmp   _nb_kernel102_x86_64_sse.nb102_single_loop
_nb_kernel102_x86_64_sse.nb102_updateouterdata: 
        movl  nb102_ii3(%rsp),%ecx
        movq  nb102_faction(%rbp),%rdi
        movq  nb102_fshift(%rbp),%rsi
        movl  nb102_is3(%rsp),%edx

        ## accumulate Oi forces in xmm0, xmm1, xmm2 
        movaps nb102_fixO(%rsp),%xmm0
        movaps nb102_fiyO(%rsp),%xmm1
        movaps nb102_fizO(%rsp),%xmm2

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
        movaps nb102_fixH1(%rsp),%xmm0
        movaps nb102_fiyH1(%rsp),%xmm1
        movaps nb102_fizH1(%rsp),%xmm2

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
        movaps nb102_fixH2(%rsp),%xmm0
        movaps nb102_fiyH2(%rsp),%xmm1
        movaps nb102_fizH2(%rsp),%xmm2

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
        movl nb102_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb102_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb102_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb102_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb102_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102_x86_64_sse.nb102_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102_n(%rsp)
        jmp _nb_kernel102_x86_64_sse.nb102_outer
_nb_kernel102_x86_64_sse.nb102_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102_x86_64_sse.nb102_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102_x86_64_sse.nb102_threadloop
_nb_kernel102_x86_64_sse.nb102_end: 

        movl nb102_nouter(%rsp),%eax
        movl nb102_ninner(%rsp),%ebx
        movq nb102_outeriter(%rbp),%rcx
        movq nb102_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1512,%rsp
        emms

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel102nf_x86_64_sse
.globl _nb_kernel102nf_x86_64_sse
nb_kernel102nf_x86_64_sse:      
_nb_kernel102nf_x86_64_sse:     
.set nb102nf_fshift, 16
.set nb102nf_gid, 24
.set nb102nf_pos, 32
.set nb102nf_faction, 40
.set nb102nf_charge, 48
.set nb102nf_p_facel, 56
.set nb102nf_argkrf, 64
.set nb102nf_argcrf, 72
.set nb102nf_Vc, 80
.set nb102nf_type, 88
.set nb102nf_p_ntype, 96
.set nb102nf_vdwparam, 104
.set nb102nf_Vvdw, 112
.set nb102nf_p_tabscale, 120
.set nb102nf_VFtab, 128
.set nb102nf_invsqrta, 136
.set nb102nf_dvda, 144
.set nb102nf_p_gbtabscale, 152
.set nb102nf_GBtab, 160
.set nb102nf_p_nthreads, 168
.set nb102nf_count, 176
.set nb102nf_mtx, 184
.set nb102nf_outeriter, 192
.set nb102nf_inneriter, 200
.set nb102nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb102nf_ixO, 0
.set nb102nf_iyO, 16
.set nb102nf_izO, 32
.set nb102nf_ixH1, 48
.set nb102nf_iyH1, 64
.set nb102nf_izH1, 80
.set nb102nf_ixH2, 96
.set nb102nf_iyH2, 112
.set nb102nf_izH2, 128
.set nb102nf_jxO, 144
.set nb102nf_jyO, 160
.set nb102nf_jzO, 176
.set nb102nf_jxH1, 192
.set nb102nf_jyH1, 208
.set nb102nf_jzH1, 224
.set nb102nf_jxH2, 240
.set nb102nf_jyH2, 256
.set nb102nf_jzH2, 272
.set nb102nf_qqOO, 288
.set nb102nf_qqOH, 304
.set nb102nf_qqHH, 320
.set nb102nf_vctot, 336
.set nb102nf_half, 352
.set nb102nf_three, 368
.set nb102nf_rsqOO, 384
.set nb102nf_rsqOH1, 400
.set nb102nf_rsqOH2, 416
.set nb102nf_rsqH1O, 432
.set nb102nf_rsqH1H1, 448
.set nb102nf_rsqH1H2, 464
.set nb102nf_rsqH2O, 480
.set nb102nf_rsqH2H1, 496
.set nb102nf_rsqH2H2, 512
.set nb102nf_rinvOO, 528
.set nb102nf_rinvOH1, 544
.set nb102nf_rinvOH2, 560
.set nb102nf_rinvH1O, 576
.set nb102nf_rinvH1H1, 592
.set nb102nf_rinvH1H2, 608
.set nb102nf_rinvH2O, 624
.set nb102nf_rinvH2H1, 640
.set nb102nf_rinvH2H2, 656
.set nb102nf_is3, 672
.set nb102nf_ii3, 676
.set nb102nf_nri, 680
.set nb102nf_iinr, 688
.set nb102nf_jindex, 696
.set nb102nf_jjnr, 704
.set nb102nf_shift, 712
.set nb102nf_shiftvec, 720
.set nb102nf_facel, 728
.set nb102nf_innerjjnr, 736
.set nb102nf_innerk, 744
.set nb102nf_n, 748
.set nb102nf_nn1, 752
.set nb102nf_nouter, 756
.set nb102nf_ninner, 760

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $776,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb102nf_nouter(%rsp)
        movl %eax,nb102nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb102nf_nri(%rsp)
        movq %rsi,nb102nf_iinr(%rsp)
        movq %rdx,nb102nf_jindex(%rsp)
        movq %rcx,nb102nf_jjnr(%rsp)
        movq %r8,nb102nf_shift(%rsp)
        movq %r9,nb102nf_shiftvec(%rsp)
        movq nb102nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb102nf_facel(%rsp)


        ## assume we have at least one i particle - start directly 
        movq  nb102nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb102nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb102nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb102nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102nf_qqOO(%rsp)
        movaps %xmm4,nb102nf_qqOH(%rsp)
        movaps %xmm5,nb102nf_qqHH(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb102nf_half(%rsp)
        movss nb102nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb102nf_half(%rsp)
        movaps %xmm3,nb102nf_three(%rsp)

_nb_kernel102nf_x86_64_sse.nb102nf_threadloop: 
        movq  nb102nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel102nf_x86_64_sse.nb102nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel102nf_x86_64_sse.nb102nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb102nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb102nf_n(%rsp)
        movl %ebx,nb102nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel102nf_x86_64_sse.nb102nf_outerstart
        jmp _nb_kernel102nf_x86_64_sse.nb102nf_end

_nb_kernel102nf_x86_64_sse.nb102nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb102nf_nouter(%rsp),%ebx
        movl %ebx,nb102nf_nouter(%rsp)

_nb_kernel102nf_x86_64_sse.nb102nf_outer: 
        movq  nb102nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb102nf_is3(%rsp)            ## store is3 

        movq  nb102nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb102nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb102nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb102nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb102nf_ixO(%rsp)
        movaps %xmm4,nb102nf_iyO(%rsp)
        movaps %xmm5,nb102nf_izO(%rsp)

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
        movaps %xmm0,nb102nf_ixH1(%rsp)
        movaps %xmm1,nb102nf_iyH1(%rsp)
        movaps %xmm2,nb102nf_izH1(%rsp)
        movaps %xmm3,nb102nf_ixH2(%rsp)
        movaps %xmm4,nb102nf_iyH2(%rsp)
        movaps %xmm5,nb102nf_izH2(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb102nf_vctot(%rsp)

        movq  nb102nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movq  nb102nf_pos(%rbp),%rsi
        movq  nb102nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb102nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb102nf_ninner(%rsp),%ecx
        movl  %ecx,nb102nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb102nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel102nf_x86_64_sse.nb102nf_unroll_loop
        jmp   _nb_kernel102nf_x86_64_sse.nb102nf_single_check
_nb_kernel102nf_x86_64_sse.nb102nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb102nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb102nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb102nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps %xmm0,nb102nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb102nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb102nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb102nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb102nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb102nf_jyH2(%rsp)

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
        movaps %xmm0,nb102nf_jzO(%rsp)
        movaps %xmm1,nb102nf_jzH1(%rsp)
        movaps %xmm2,nb102nf_jzH2(%rsp)

        movaps nb102nf_ixO(%rsp),%xmm0
        movaps nb102nf_iyO(%rsp),%xmm1
        movaps nb102nf_izO(%rsp),%xmm2
        movaps nb102nf_ixO(%rsp),%xmm3
        movaps nb102nf_iyO(%rsp),%xmm4
        movaps nb102nf_izO(%rsp),%xmm5
        subps  nb102nf_jxO(%rsp),%xmm0
        subps  nb102nf_jyO(%rsp),%xmm1
        subps  nb102nf_jzO(%rsp),%xmm2
        subps  nb102nf_jxH1(%rsp),%xmm3
        subps  nb102nf_jyH1(%rsp),%xmm4
        subps  nb102nf_jzH1(%rsp),%xmm5

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
        movaps %xmm0,nb102nf_rsqOO(%rsp)
        movaps %xmm3,nb102nf_rsqOH1(%rsp)

        movaps nb102nf_ixO(%rsp),%xmm0
        movaps nb102nf_iyO(%rsp),%xmm1
        movaps nb102nf_izO(%rsp),%xmm2
        movaps nb102nf_ixH1(%rsp),%xmm3
        movaps nb102nf_iyH1(%rsp),%xmm4
        movaps nb102nf_izH1(%rsp),%xmm5
        subps  nb102nf_jxH2(%rsp),%xmm0
        subps  nb102nf_jyH2(%rsp),%xmm1
        subps  nb102nf_jzH2(%rsp),%xmm2
        subps  nb102nf_jxO(%rsp),%xmm3
        subps  nb102nf_jyO(%rsp),%xmm4
        subps  nb102nf_jzO(%rsp),%xmm5

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
        movaps %xmm0,nb102nf_rsqOH2(%rsp)
        movaps %xmm3,nb102nf_rsqH1O(%rsp)

        movaps nb102nf_ixH1(%rsp),%xmm0
        movaps nb102nf_iyH1(%rsp),%xmm1
        movaps nb102nf_izH1(%rsp),%xmm2
        movaps nb102nf_ixH1(%rsp),%xmm3
        movaps nb102nf_iyH1(%rsp),%xmm4
        movaps nb102nf_izH1(%rsp),%xmm5
        subps  nb102nf_jxH1(%rsp),%xmm0
        subps  nb102nf_jyH1(%rsp),%xmm1
        subps  nb102nf_jzH1(%rsp),%xmm2
        subps  nb102nf_jxH2(%rsp),%xmm3
        subps  nb102nf_jyH2(%rsp),%xmm4
        subps  nb102nf_jzH2(%rsp),%xmm5

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
        movaps %xmm0,nb102nf_rsqH1H1(%rsp)
        movaps %xmm3,nb102nf_rsqH1H2(%rsp)

        movaps nb102nf_ixH2(%rsp),%xmm0
        movaps nb102nf_iyH2(%rsp),%xmm1
        movaps nb102nf_izH2(%rsp),%xmm2
        movaps nb102nf_ixH2(%rsp),%xmm3
        movaps nb102nf_iyH2(%rsp),%xmm4
        movaps nb102nf_izH2(%rsp),%xmm5
        subps  nb102nf_jxO(%rsp),%xmm0
        subps  nb102nf_jyO(%rsp),%xmm1
        subps  nb102nf_jzO(%rsp),%xmm2
        subps  nb102nf_jxH1(%rsp),%xmm3
        subps  nb102nf_jyH1(%rsp),%xmm4
        subps  nb102nf_jzH1(%rsp),%xmm5

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
        movaps %xmm0,nb102nf_rsqH2O(%rsp)
        movaps %xmm4,nb102nf_rsqH2H1(%rsp)

        movaps nb102nf_ixH2(%rsp),%xmm0
        movaps nb102nf_iyH2(%rsp),%xmm1
        movaps nb102nf_izH2(%rsp),%xmm2
        subps  nb102nf_jxH2(%rsp),%xmm0
        subps  nb102nf_jyH2(%rsp),%xmm1
        subps  nb102nf_jzH2(%rsp),%xmm2

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb102nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb102nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb102nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb102nf_rinvH2H1(%rsp)

        rsqrtps nb102nf_rsqOO(%rsp),%xmm1
        rsqrtps nb102nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqOO(%rsp),%xmm1
        mulps   nb102nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%rsp),%xmm3
        mulps   nb102nf_half(%rsp),%xmm7
        movaps  %xmm3,nb102nf_rinvOO(%rsp)
        movaps  %xmm7,nb102nf_rinvOH1(%rsp)

        rsqrtps nb102nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb102nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqOH2(%rsp),%xmm1
        mulps   nb102nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%rsp),%xmm3
        mulps   nb102nf_half(%rsp),%xmm7
        movaps  %xmm3,nb102nf_rinvOH2(%rsp)
        movaps  %xmm7,nb102nf_rinvH1O(%rsp)

        rsqrtps nb102nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb102nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb102nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb102nf_rsqH1H1(%rsp),%xmm1
        mulps   nb102nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%rsp),%xmm3
        mulps   nb102nf_half(%rsp),%xmm7
        movaps  %xmm3,nb102nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb102nf_rinvH1H2(%rsp)

        rsqrtps nb102nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb102nf_three(%rsp),%xmm3
        mulps   nb102nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102nf_half(%rsp),%xmm3
        movaps  %xmm3,nb102nf_rinvH2O(%rsp)

        ## sum OO pot in xmm0, OH in xmm1 HH in xmm2 
        movaps nb102nf_rinvOO(%rsp),%xmm0
        movaps nb102nf_rinvOH1(%rsp),%xmm1
        movaps nb102nf_rinvH1H1(%rsp),%xmm2
        addps  nb102nf_rinvOH2(%rsp),%xmm1
        addps  nb102nf_rinvH1H2(%rsp),%xmm2
        addps  nb102nf_rinvH1O(%rsp),%xmm1
        addps  nb102nf_rinvH2H1(%rsp),%xmm2
        addps  nb102nf_rinvH2O(%rsp),%xmm1
        addps  nb102nf_rinvH2H2(%rsp),%xmm2

        mulps  nb102nf_qqOO(%rsp),%xmm0
        mulps  nb102nf_qqOH(%rsp),%xmm1
        mulps  nb102nf_qqHH(%rsp),%xmm2
        addps  nb102nf_vctot(%rsp),%xmm0
        addps  %xmm2,%xmm1
        addps  %xmm1,%xmm0
        movaps  %xmm0,nb102nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb102nf_innerk(%rsp)
        jl    _nb_kernel102nf_x86_64_sse.nb102nf_single_check
        jmp   _nb_kernel102nf_x86_64_sse.nb102nf_unroll_loop
_nb_kernel102nf_x86_64_sse.nb102nf_single_check: 
        addl $4,nb102nf_innerk(%rsp)
        jnz   _nb_kernel102nf_x86_64_sse.nb102nf_single_loop
        jmp   _nb_kernel102nf_x86_64_sse.nb102nf_updateouterdata
_nb_kernel102nf_x86_64_sse.nb102nf_single_loop: 
        movq  nb102nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb102nf_innerjjnr(%rsp)

        movq nb102nf_pos(%rbp),%rsi
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
        movaps  nb102nf_ixO(%rsp),%xmm0
        movaps  nb102nf_iyO(%rsp),%xmm1
        movaps  nb102nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb102nf_jxO(%rsp)
        movaps %xmm4,nb102nf_jyO(%rsp)
        movaps %xmm5,nb102nf_jzO(%rsp)
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
        movaps  nb102nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb102nf_half(%rsp),%xmm3   ## rinv iO - j water 

        xorps   %xmm1,%xmm1

        xorps   %xmm4,%xmm4

        ## fetch charges to xmm4 (temporary) 
        movss   nb102nf_qqOO(%rsp),%xmm4

        movhps  nb102nf_qqOH(%rsp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 

        addps   nb102nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb102nf_vctot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously: 
        movaps  nb102nf_ixH1(%rsp),%xmm0
        movaps  nb102nf_iyH1(%rsp),%xmm1
        movaps  nb102nf_izH1(%rsp),%xmm2
        movaps  nb102nf_ixH2(%rsp),%xmm3
        movaps  nb102nf_iyH2(%rsp),%xmm4
        movaps  nb102nf_izH2(%rsp),%xmm5
        subps   nb102nf_jxO(%rsp),%xmm0
        subps   nb102nf_jyO(%rsp),%xmm1
        subps   nb102nf_jzO(%rsp),%xmm2
        subps   nb102nf_jxO(%rsp),%xmm3
        subps   nb102nf_jyO(%rsp),%xmm4
        subps   nb102nf_jzO(%rsp),%xmm5
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
        movaps  nb102nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb102nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb102nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movss   nb102nf_qqOH(%rsp),%xmm6
        movaps  %xmm7,%xmm4
        movhps  nb102nf_qqHH(%rsp),%xmm6
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        addps   %xmm7,%xmm3     ## total vcoul 
        addps   nb102nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb102nf_vctot(%rsp)

        decl  nb102nf_innerk(%rsp)
        jz    _nb_kernel102nf_x86_64_sse.nb102nf_updateouterdata
        jmp   _nb_kernel102nf_x86_64_sse.nb102nf_single_loop
_nb_kernel102nf_x86_64_sse.nb102nf_updateouterdata: 
        ## get n from stack
        movl nb102nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb102nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb102nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb102nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb102nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel102nf_x86_64_sse.nb102nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb102nf_n(%rsp)
        jmp _nb_kernel102nf_x86_64_sse.nb102nf_outer
_nb_kernel102nf_x86_64_sse.nb102nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb102nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel102nf_x86_64_sse.nb102nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel102nf_x86_64_sse.nb102nf_threadloop
_nb_kernel102nf_x86_64_sse.nb102nf_end: 


        movl nb102nf_nouter(%rsp),%eax
        movl nb102nf_ninner(%rsp),%ebx
        movq nb102nf_outeriter(%rbp),%rcx
        movq nb102nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $776,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret




