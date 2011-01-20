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





.globl nb_kernel312_x86_64_sse
.globl _nb_kernel312_x86_64_sse
nb_kernel312_x86_64_sse:        
_nb_kernel312_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb312_fshift, 16
.set nb312_gid, 24
.set nb312_pos, 32
.set nb312_faction, 40
.set nb312_charge, 48
.set nb312_p_facel, 56
.set nb312_argkrf, 64
.set nb312_argcrf, 72
.set nb312_Vc, 80
.set nb312_type, 88
.set nb312_p_ntype, 96
.set nb312_vdwparam, 104
.set nb312_Vvdw, 112
.set nb312_p_tabscale, 120
.set nb312_VFtab, 128
.set nb312_invsqrta, 136
.set nb312_dvda, 144
.set nb312_p_gbtabscale, 152
.set nb312_GBtab, 160
.set nb312_p_nthreads, 168
.set nb312_count, 176
.set nb312_mtx, 184
.set nb312_outeriter, 192
.set nb312_inneriter, 200
.set nb312_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb312_ixO, 0
.set nb312_iyO, 16
.set nb312_izO, 32
.set nb312_ixH1, 48
.set nb312_iyH1, 64
.set nb312_izH1, 80
.set nb312_ixH2, 96
.set nb312_iyH2, 112
.set nb312_izH2, 128
.set nb312_jxO, 144
.set nb312_jyO, 160
.set nb312_jzO, 176
.set nb312_jxH1, 192
.set nb312_jyH1, 208
.set nb312_jzH1, 224
.set nb312_jxH2, 240
.set nb312_jyH2, 256
.set nb312_jzH2, 272
.set nb312_dxOO, 288
.set nb312_dyOO, 304
.set nb312_dzOO, 320
.set nb312_dxOH1, 336
.set nb312_dyOH1, 352
.set nb312_dzOH1, 368
.set nb312_dxOH2, 384
.set nb312_dyOH2, 400
.set nb312_dzOH2, 416
.set nb312_dxH1O, 432
.set nb312_dyH1O, 448
.set nb312_dzH1O, 464
.set nb312_dxH1H1, 480
.set nb312_dyH1H1, 496
.set nb312_dzH1H1, 512
.set nb312_dxH1H2, 528
.set nb312_dyH1H2, 544
.set nb312_dzH1H2, 560
.set nb312_dxH2O, 576
.set nb312_dyH2O, 592
.set nb312_dzH2O, 608
.set nb312_dxH2H1, 624
.set nb312_dyH2H1, 640
.set nb312_dzH2H1, 656
.set nb312_dxH2H2, 672
.set nb312_dyH2H2, 688
.set nb312_dzH2H2, 704
.set nb312_qqOO, 720
.set nb312_qqOH, 736
.set nb312_qqHH, 752
.set nb312_two, 768
.set nb312_tsc, 784
.set nb312_c6, 800
.set nb312_c12, 816
.set nb312_six, 832
.set nb312_twelve, 848
.set nb312_vctot, 864
.set nb312_Vvdwtot, 880
.set nb312_fixO, 896
.set nb312_fiyO, 912
.set nb312_fizO, 928
.set nb312_fixH1, 944
.set nb312_fiyH1, 960
.set nb312_fizH1, 976
.set nb312_fixH2, 992
.set nb312_fiyH2, 1008
.set nb312_fizH2, 1024
.set nb312_fjxO, 1040
.set nb312_fjyO, 1056
.set nb312_fjzO, 1072
.set nb312_fjxH1, 1088
.set nb312_fjyH1, 1104
.set nb312_fjzH1, 1120
.set nb312_fjxH2, 1136
.set nb312_fjyH2, 1152
.set nb312_fjzH2, 1168
.set nb312_half, 1184
.set nb312_three, 1200
.set nb312_epsO, 1216
.set nb312_epsH1, 1232
.set nb312_epsH2, 1248
.set nb312_rsqH1O, 1264
.set nb312_rsqH1H1, 1280
.set nb312_rsqH1H2, 1296
.set nb312_rsqH2O, 1312
.set nb312_rsqH2H1, 1328
.set nb312_rsqH2H2, 1344
.set nb312_rinvOO, 1360
.set nb312_rinvOH1, 1376
.set nb312_rinvOH2, 1392
.set nb312_rinvH1O, 1408
.set nb312_rinvH1H1, 1424
.set nb312_rinvH1H2, 1440
.set nb312_rinvH2O, 1456
.set nb312_rinvH2H1, 1472
.set nb312_rinvH2H2, 1488
.set nb312_fstmp, 1504
.set nb312_is3, 1520
.set nb312_ii3, 1524
.set nb312_nri, 1528
.set nb312_iinr, 1536
.set nb312_jindex, 1544
.set nb312_jjnr, 1552
.set nb312_shift, 1560
.set nb312_shiftvec, 1568
.set nb312_facel, 1576
.set nb312_innerjjnr, 1584
.set nb312_innerk, 1592
.set nb312_n, 1596
.set nb312_nn1, 1600
.set nb312_nouter, 1604
.set nb312_ninner, 1608
        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
        subq $1624,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb312_nouter(%rsp)
        movl %eax,nb312_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb312_nri(%rsp)
        movq %rsi,nb312_iinr(%rsp)
        movq %rdx,nb312_jindex(%rsp)
        movq %rcx,nb312_jjnr(%rsp)
        movq %r8,nb312_shift(%rsp)
        movq %r9,nb312_shiftvec(%rsp)
        movq nb312_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb312_facel(%rsp)

        movq nb312_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb312_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb312_half(%rsp)
        movss nb312_half(%rsp),%xmm1
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
        movaps %xmm1,nb312_half(%rsp)
        movaps %xmm2,nb312_two(%rsp)
        movaps %xmm3,nb312_three(%rsp)
        movaps %xmm4,nb312_six(%rsp)
        movaps %xmm5,nb312_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb312_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb312_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb312_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb312_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb312_qqOO(%rsp)
        movaps %xmm4,nb312_qqOH(%rsp)
        movaps %xmm5,nb312_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb312_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb312_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb312_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb312_c6(%rsp)
        movaps %xmm1,nb312_c12(%rsp)

_nb_kernel312_x86_64_sse.nb312_threadloop: 
        movq  nb312_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel312_x86_64_sse.nb312_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel312_x86_64_sse.nb312_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb312_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb312_n(%rsp)
        movl %ebx,nb312_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel312_x86_64_sse.nb312_outerstart
        jmp _nb_kernel312_x86_64_sse.nb312_end
_nb_kernel312_x86_64_sse.nb312_outerstart: 
        ## ebx contains number of outer iterations
        addl nb312_nouter(%rsp),%ebx
        movl %ebx,nb312_nouter(%rsp)

_nb_kernel312_x86_64_sse.nb312_outer: 
        movq  nb312_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb312_is3(%rsp)      ## store is3 

        movq  nb312_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb312_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb312_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb312_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb312_ixO(%rsp)
        movaps %xmm4,nb312_iyO(%rsp)
        movaps %xmm5,nb312_izO(%rsp)

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
        movaps %xmm0,nb312_ixH1(%rsp)
        movaps %xmm1,nb312_iyH1(%rsp)
        movaps %xmm2,nb312_izH1(%rsp)
        movaps %xmm3,nb312_ixH2(%rsp)
        movaps %xmm4,nb312_iyH2(%rsp)
        movaps %xmm5,nb312_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb312_vctot(%rsp)
        movaps %xmm4,nb312_Vvdwtot(%rsp)
        movaps %xmm4,nb312_fixO(%rsp)
        movaps %xmm4,nb312_fiyO(%rsp)
        movaps %xmm4,nb312_fizO(%rsp)
        movaps %xmm4,nb312_fixH1(%rsp)
        movaps %xmm4,nb312_fiyH1(%rsp)
        movaps %xmm4,nb312_fizH1(%rsp)
        movaps %xmm4,nb312_fixH2(%rsp)
        movaps %xmm4,nb312_fiyH2(%rsp)
        movaps %xmm4,nb312_fizH2(%rsp)

        movq  nb312_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb312_pos(%rbp),%rsi
        movq  nb312_faction(%rbp),%rdi
        movq  nb312_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb312_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb312_ninner(%rsp),%ecx
        movl  %ecx,nb312_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb312_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel312_x86_64_sse.nb312_unroll_loop
        jmp   _nb_kernel312_x86_64_sse.nb312_single_check
_nb_kernel312_x86_64_sse.nb312_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb312_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb312_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb312_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## load j O coordinates
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

    subps nb312_ixO(%rsp),%xmm0
    subps nb312_iyO(%rsp),%xmm1
    subps nb312_izO(%rsp),%xmm2
    subps nb312_ixH1(%rsp),%xmm3
    subps nb312_iyH1(%rsp),%xmm4
    subps nb312_izH1(%rsp),%xmm5
    subps nb312_ixH2(%rsp),%xmm6
    subps nb312_iyH2(%rsp),%xmm7
    subps nb312_izH2(%rsp),%xmm8

        movaps %xmm0,nb312_dxOO(%rsp)
        movaps %xmm1,nb312_dyOO(%rsp)
        movaps %xmm2,nb312_dzOO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb312_dxH1O(%rsp)
        movaps %xmm4,nb312_dyH1O(%rsp)
        movaps %xmm5,nb312_dzH1O(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb312_dxH2O(%rsp)
        movaps %xmm7,nb312_dyH2O(%rsp)
        movaps %xmm8,nb312_dzH2O(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

    movd %eax,%mm0 ## save j3 in mm0-mm3
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

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

        movaps  nb312_three(%rsp),%xmm9
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

        movaps  nb312_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOO 
        mulps   %xmm4,%xmm10 ## rinvH1O
    mulps   %xmm4,%xmm11 ## rinvH2O

        movaps  %xmm9,nb312_rinvOO(%rsp)
        movaps  %xmm10,nb312_rinvH1O(%rsp)
        movaps  %xmm11,nb312_rinvH2O(%rsp)

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb312_tsc(%rsp),%xmm1
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

    movq nb312_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb312_epsO(%rsp)
    movaps    %xmm3,nb312_epsH1(%rsp)
    movaps    %xmm6,nb312_epsH2(%rsp)

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

    mulps  nb312_epsO(%rsp),%xmm3     ## Heps
    mulps  nb312_epsH1(%rsp),%xmm7
    mulps  nb312_epsH2(%rsp),%xmm11
    mulps  nb312_epsO(%rsp),%xmm2     ## Geps
    mulps  nb312_epsH1(%rsp),%xmm6
    mulps  nb312_epsH2(%rsp),%xmm10
    mulps  nb312_epsO(%rsp),%xmm3     ## Heps2
    mulps  nb312_epsH1(%rsp),%xmm7
    mulps  nb312_epsH2(%rsp),%xmm11

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
    mulps  nb312_epsO(%rsp),%xmm1     ## eps*Fp
    mulps  nb312_epsH1(%rsp),%xmm5
    mulps  nb312_epsH2(%rsp),%xmm9
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  nb312_qqOO(%rsp),%xmm1     ## VV*qq = vcoul
    mulps  nb312_qqOH(%rsp),%xmm5
    mulps  nb312_qqOH(%rsp),%xmm9
    mulps  nb312_qqOO(%rsp),%xmm3      ## FF*qq = fij
    mulps  nb312_qqOH(%rsp),%xmm7
    mulps  nb312_qqOH(%rsp),%xmm11

    ## calculate LJ
    movaps nb312_rinvOO(%rsp),%xmm12
    mulps  %xmm12,%xmm12 ## rinvsq
    movaps %xmm12,%xmm13 ## rinvsq
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm13,%xmm12 ## rinv6
    movaps %xmm12,%xmm13 ## rinv6
    mulps  %xmm12,%xmm12 ## rinv12
        mulps  nb312_c6(%rsp),%xmm13
        mulps  nb312_c12(%rsp),%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14

        addps  nb312_Vvdwtot(%rsp),%xmm14
        mulps  nb312_six(%rsp),%xmm13
        mulps  nb312_twelve(%rsp),%xmm12
        movaps %xmm14,nb312_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal    
    mulps  nb312_rinvOO(%rsp),%xmm12
    movaps %xmm12,nb312_fstmp(%rsp)

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

    ## accumulate vctot
    addps  nb312_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb312_vctot(%rsp)

    movaps nb312_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

        ## move j O forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm11 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm12 ## jxOc jyOc  -   -

    ## xmm11: jxOa jyOa jxOb jyOb 
    ## xmm12: jxOc jyOc jxOd jyOd
    ## xmm13: jzOa jzOb jzOc jzOd


    movaps nb312_fstmp(%rsp),%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8
    movhps (%rdi,%rbx,4),%xmm11 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm12 ## jxOc jyOc jxOd jyOd 

    mulps  nb312_rinvOO(%rsp),%xmm0
    mulps  nb312_rinvH1O(%rsp),%xmm4
    mulps  nb312_rinvH2O(%rsp),%xmm8

    movss  8(%rdi,%rax,4),%xmm13    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm14    ## jzOc  -  -  -

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

    movhps 8(%rdi,%rbx,4),%xmm13    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm14    ## jzOc  -  jzOd -    
    shufps $136,%xmm14,%xmm13 ## 10001000 => jzOa jzOb jzOc jzOd

        mulps nb312_dxOO(%rsp),%xmm0
        mulps nb312_dyOO(%rsp),%xmm1
        mulps nb312_dzOO(%rsp),%xmm2
        mulps nb312_dxH1O(%rsp),%xmm3
        mulps nb312_dyH1O(%rsp),%xmm4
        mulps nb312_dzH1O(%rsp),%xmm5
        mulps nb312_dxH2O(%rsp),%xmm6
        mulps nb312_dyH2O(%rsp),%xmm7
        mulps nb312_dzH2O(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb312_fixO(%rsp),%xmm0
    addps nb312_fiyO(%rsp),%xmm1
    addps nb312_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb312_fixH1(%rsp),%xmm3
    addps nb312_fiyH1(%rsp),%xmm4
    addps nb312_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb312_fixH2(%rsp),%xmm6
    addps nb312_fiyH2(%rsp),%xmm7
    addps nb312_fizH2(%rsp),%xmm8

    movaps %xmm0,nb312_fixO(%rsp)
    movaps %xmm1,nb312_fiyO(%rsp)
    movaps %xmm2,nb312_fizO(%rsp)
    movaps %xmm3,nb312_fixH1(%rsp)
    movaps %xmm4,nb312_fiyH1(%rsp)
    movaps %xmm5,nb312_fizH1(%rsp)
    movaps %xmm6,nb312_fixH2(%rsp)
    movaps %xmm7,nb312_fiyH2(%rsp)
    movaps %xmm8,nb312_fizH2(%rsp)

    ## xmm11 = fOx
    ## xmm12 = fOy
    ## xmm13 = fOz
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fOzc fOzd

    pshufd $1,%xmm13,%xmm3
    pshufd $3,%xmm13,%xmm4

    movlps %xmm11,(%rdi,%rax,4)
    movhps %xmm11,(%rdi,%rbx,4)
    movlps %xmm12,(%rdi,%rcx,4)
    movhps %xmm12,(%rdi,%rdx,4)
    movss  %xmm13,8(%rdi,%rax,4)
    movss  %xmm14,8(%rdi,%rcx,4)
    movss  %xmm3,8(%rdi,%rbx,4)
    movss  %xmm4,8(%rdi,%rdx,4)

        ## move j H1 coordinates to local temp variables 
        movq  nb312_pos(%rbp),%rsi
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

    subps nb312_ixO(%rsp),%xmm0
    subps nb312_iyO(%rsp),%xmm1
    subps nb312_izO(%rsp),%xmm2
    subps nb312_ixH1(%rsp),%xmm3
    subps nb312_iyH1(%rsp),%xmm4
    subps nb312_izH1(%rsp),%xmm5
    subps nb312_ixH2(%rsp),%xmm6
    subps nb312_iyH2(%rsp),%xmm7
    subps nb312_izH2(%rsp),%xmm8

        movaps %xmm0,nb312_dxOH1(%rsp)
        movaps %xmm1,nb312_dyOH1(%rsp)
        movaps %xmm2,nb312_dzOH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb312_dxH1H1(%rsp)
        movaps %xmm4,nb312_dyH1H1(%rsp)
        movaps %xmm5,nb312_dzH1H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb312_dxH2H1(%rsp)
        movaps %xmm7,nb312_dyH2H1(%rsp)
        movaps %xmm8,nb312_dzH2H1(%rsp)
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

        movaps  nb312_three(%rsp),%xmm9
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

        movaps  nb312_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH1
        mulps   %xmm4,%xmm10 ## rinvH1H1
    mulps   %xmm4,%xmm11 ## rinvH2H1

        movaps  %xmm9,nb312_rinvOH1(%rsp)
        movaps  %xmm10,nb312_rinvH1H1(%rsp)
        movaps  %xmm11,nb312_rinvH2H1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb312_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    movq nb312_VFtab(%rbp),%rsi

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

    movaps    %xmm0,nb312_epsO(%rsp)
    movaps    %xmm3,nb312_epsH1(%rsp)
    movaps    %xmm6,nb312_epsH2(%rsp)


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

    movaps nb312_epsO(%rsp),%xmm12
    movaps nb312_epsH1(%rsp),%xmm13
    movaps nb312_epsH2(%rsp),%xmm14

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
    movaps nb312_qqOH(%rsp),%xmm12
    movaps nb312_qqHH(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

    ## accumulate vctot
    addps  nb312_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb312_vctot(%rsp)

    movaps nb312_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

        ## move j H1 forces to local temp variables 
    movlps 12(%rdi,%rax,4),%xmm11    ## jxH1a jyH1a  -   -
    movlps 12(%rdi,%rcx,4),%xmm12    ## jxH1c jyH1c  -   -

    ## xmm11: jxH1a jyH1a jxH1b jyH1b 
    ## xmm12: jxH1c jyH1c jxH1d jyH1d
    ## xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb312_rinvOH1(%rsp),%xmm3
    mulps  nb312_rinvH1H1(%rsp),%xmm7
    mulps  nb312_rinvH2H1(%rsp),%xmm10

    movhps 12(%rdi,%rbx,4),%xmm11    ## jxH1a jyH1a jxH1b jyH1b 
    movhps 12(%rdi,%rdx,4),%xmm12    ## jxH1c jyH1c jxH1d jyH1d 

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movss  20(%rdi,%rax,4),%xmm13    ## jzH1a  -  -  -
    movss  20(%rdi,%rcx,4),%xmm14    ## jzH1c  -  -  -

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

    movhps 20(%rdi,%rbx,4),%xmm13    ## jzH1a  -  jzH1b  -
    movhps 20(%rdi,%rdx,4),%xmm14    ## jzH1c  -  jzH1d -    
    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

        mulps nb312_dxOH1(%rsp),%xmm0
        mulps nb312_dyOH1(%rsp),%xmm1
        mulps nb312_dzOH1(%rsp),%xmm2
        mulps nb312_dxH1H1(%rsp),%xmm3
        mulps nb312_dyH1H1(%rsp),%xmm4
        mulps nb312_dzH1H1(%rsp),%xmm5
        mulps nb312_dxH2H1(%rsp),%xmm6
        mulps nb312_dyH2H1(%rsp),%xmm7
        mulps nb312_dzH2H1(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb312_fixO(%rsp),%xmm0
    addps nb312_fiyO(%rsp),%xmm1
    addps nb312_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb312_fixH1(%rsp),%xmm3
    addps nb312_fiyH1(%rsp),%xmm4
    addps nb312_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb312_fixH2(%rsp),%xmm6
    addps nb312_fiyH2(%rsp),%xmm7
    addps nb312_fizH2(%rsp),%xmm8

    movaps %xmm0,nb312_fixO(%rsp)
    movaps %xmm1,nb312_fiyO(%rsp)
    movaps %xmm2,nb312_fizO(%rsp)
    movaps %xmm3,nb312_fixH1(%rsp)
    movaps %xmm4,nb312_fiyH1(%rsp)
    movaps %xmm5,nb312_fizH1(%rsp)
    movaps %xmm6,nb312_fixH2(%rsp)
    movaps %xmm7,nb312_fiyH2(%rsp)
    movaps %xmm8,nb312_fizH2(%rsp)

    ## xmm11 = fH1x
    ## xmm12 = fH1y
    ## xmm13 = fH1z
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fH1zc fH1zd

    pshufd $1,%xmm13,%xmm3
    pshufd $3,%xmm13,%xmm4

    movlps %xmm11,12(%rdi,%rax,4)
    movhps %xmm11,12(%rdi,%rbx,4)
    movlps %xmm12,12(%rdi,%rcx,4)
    movhps %xmm12,12(%rdi,%rdx,4)
    movss  %xmm13,20(%rdi,%rax,4)
    movss  %xmm14,20(%rdi,%rcx,4)
    movss  %xmm3,20(%rdi,%rbx,4)
    movss  %xmm4,20(%rdi,%rdx,4)

        movq  nb312_pos(%rbp),%rsi
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

    subps nb312_ixO(%rsp),%xmm0
    subps nb312_iyO(%rsp),%xmm1
    subps nb312_izO(%rsp),%xmm2
    subps nb312_ixH1(%rsp),%xmm3
    subps nb312_iyH1(%rsp),%xmm4
    subps nb312_izH1(%rsp),%xmm5
    subps nb312_ixH2(%rsp),%xmm6
    subps nb312_iyH2(%rsp),%xmm7
    subps nb312_izH2(%rsp),%xmm8

        movaps %xmm0,nb312_dxOH2(%rsp)
        movaps %xmm1,nb312_dyOH2(%rsp)
        movaps %xmm2,nb312_dzOH2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb312_dxH1H2(%rsp)
        movaps %xmm4,nb312_dyH1H2(%rsp)
        movaps %xmm5,nb312_dzH1H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb312_dxH2H2(%rsp)
        movaps %xmm7,nb312_dyH2H2(%rsp)
        movaps %xmm8,nb312_dzH2H2(%rsp)
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

        movaps  nb312_three(%rsp),%xmm9
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

        movaps  nb312_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvOH2
        mulps   %xmm4,%xmm10 ## rinvH1H2
    mulps   %xmm4,%xmm11 ## rinvH2H2

        movaps  %xmm9,nb312_rinvOH2(%rsp)
        movaps  %xmm10,nb312_rinvH1H2(%rsp)
        movaps  %xmm11,nb312_rinvH2H2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb312_tsc(%rsp),%xmm1
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

    movq nb312_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb312_epsO(%rsp)
    movaps    %xmm3,nb312_epsH1(%rsp)
    movaps    %xmm6,nb312_epsH2(%rsp)

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

    movaps nb312_epsO(%rsp),%xmm12
    movaps nb312_epsH1(%rsp),%xmm13
    movaps nb312_epsH2(%rsp),%xmm14

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
    movaps nb312_qqOH(%rsp),%xmm12
    movaps nb312_qqHH(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm13,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm13,%xmm7
    mulps  %xmm13,%xmm11

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

    ## accumulate vctot
    addps  nb312_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb312_vctot(%rsp)

    movaps nb312_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm11    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm12    ## jxH2c jyH2c  -   -

    ## xmm11: jxH2a jyH2a jxH2b jyH2b 
    ## xmm12: jxH2c jyH2c jxH2d jyH2d
    ## xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    movhps 24(%rdi,%rbx,4),%xmm11    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm12    ## jxH2c jyH2c jxH2d jyH2d 

    mulps  nb312_rinvOH2(%rsp),%xmm3
    mulps  nb312_rinvH1H2(%rsp),%xmm7
    mulps  nb312_rinvH2H2(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movss  32(%rdi,%rax,4),%xmm13    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm14    ## jzH2c  -  -  -

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

    movss  32(%rdi,%rbx,4),%xmm15    ## jzH2b  -  -  -
    movss  32(%rdi,%rdx,4),%xmm9    ## jzH2d  -  -  -
    movlhps %xmm15,%xmm13 ## jzH2a  -  jzH2b  -
    movlhps %xmm9,%xmm14 ## jzH2c  -  jzH2d -
    shufps $136,%xmm14,%xmm13 ## 10001000 => jzH2a jzH2b jzH2c jzH2d

        mulps nb312_dxOH2(%rsp),%xmm0
        mulps nb312_dyOH2(%rsp),%xmm1
        mulps nb312_dzOH2(%rsp),%xmm2
        mulps nb312_dxH1H2(%rsp),%xmm3
        mulps nb312_dyH1H2(%rsp),%xmm4
        mulps nb312_dzH1H2(%rsp),%xmm5
        mulps nb312_dxH2H2(%rsp),%xmm6
        mulps nb312_dyH2H2(%rsp),%xmm7
        mulps nb312_dzH2H2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb312_fixO(%rsp),%xmm0
    addps nb312_fiyO(%rsp),%xmm1
    addps nb312_fizO(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb312_fixH1(%rsp),%xmm3
    addps nb312_fiyH1(%rsp),%xmm4
    addps nb312_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb312_fixH2(%rsp),%xmm6
    addps nb312_fiyH2(%rsp),%xmm7
    addps nb312_fizH2(%rsp),%xmm8

    movaps %xmm0,nb312_fixO(%rsp)
    movaps %xmm1,nb312_fiyO(%rsp)
    movaps %xmm2,nb312_fizO(%rsp)
    movaps %xmm3,nb312_fixH1(%rsp)
    movaps %xmm4,nb312_fiyH1(%rsp)
    movaps %xmm5,nb312_fizH1(%rsp)
    movaps %xmm6,nb312_fixH2(%rsp)
    movaps %xmm7,nb312_fiyH2(%rsp)
    movaps %xmm8,nb312_fizH2(%rsp)

    ## xmm11 = fH2x
    ## xmm12 = fH2y
    ## xmm13 = fH2z
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fH2zc fH2zd

    pshufd $1,%xmm13,%xmm3
    pshufd $3,%xmm13,%xmm4

    movlps %xmm11,24(%rdi,%rax,4)
    movhps %xmm11,24(%rdi,%rbx,4)
    movlps %xmm12,24(%rdi,%rcx,4)
    movhps %xmm12,24(%rdi,%rdx,4)
    movss  %xmm13,32(%rdi,%rax,4)
    movss  %xmm14,32(%rdi,%rcx,4)
    movss  %xmm3,32(%rdi,%rbx,4)
    movss  %xmm4,32(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb312_innerk(%rsp)
        jl    _nb_kernel312_x86_64_sse.nb312_single_check
        jmp   _nb_kernel312_x86_64_sse.nb312_unroll_loop
_nb_kernel312_x86_64_sse.nb312_single_check: 
        addl $4,nb312_innerk(%rsp)
        jnz   _nb_kernel312_x86_64_sse.nb312_single_loop
        jmp   _nb_kernel312_x86_64_sse.nb312_updateouterdata
_nb_kernel312_x86_64_sse.nb312_single_loop: 
        movq  nb312_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb312_innerjjnr(%rsp)

        movq nb312_pos(%rbp),%rsi
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
        movaps %xmm0,nb312_jxO(%rsp)
        movaps %xmm1,nb312_jyO(%rsp)
        movaps %xmm2,nb312_jzO(%rsp)
        subps  nb312_ixO(%rsp),%xmm0
        subps  nb312_iyO(%rsp),%xmm1
        subps  nb312_izO(%rsp),%xmm2
        movaps %xmm0,nb312_dxOO(%rsp)
        movaps %xmm1,nb312_dyOO(%rsp)
        movaps %xmm2,nb312_dzOO(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb312_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb312_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb312_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

        movq nb312_VFtab(%rbp),%rsi

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
    mulps  nb312_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb312_qqOO(%rsp),%xmm3
        movhps  nb312_qqOH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 

    addps  nb312_vctot(%rsp),%xmm5
    movaps %xmm5,nb312_vctot(%rsp)

        mulps  nb312_tsc(%rsp),%xmm3

        ## start doing lj 
        xorps  %xmm2,%xmm2
        movss  %xmm0,%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm2,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulss  nb312_c6(%rsp),%xmm1
        mulss  nb312_c12(%rsp),%xmm2
        movaps %xmm2,%xmm4
        subss  %xmm1,%xmm4
        addps  nb312_Vvdwtot(%rsp),%xmm4
        mulss  nb312_six(%rsp),%xmm1
        mulss  nb312_twelve(%rsp),%xmm2
        movaps %xmm4,nb312_Vvdwtot(%rsp)
        subss  %xmm1,%xmm2
        mulss  %xmm0,%xmm2

        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb312_dxOO(%rsp),%xmm0
        mulps   nb312_dyOO(%rsp),%xmm1
        mulps   nb312_dzOO(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb312_fjxO(%rsp)
        movaps  %xmm4,nb312_fjyO(%rsp)
        movaps  %xmm5,nb312_fjzO(%rsp)
        addps   nb312_fixO(%rsp),%xmm0
        addps   nb312_fiyO(%rsp),%xmm1
        addps   nb312_fizO(%rsp),%xmm2
        movaps  %xmm0,nb312_fixO(%rsp)
        movaps  %xmm1,nb312_fiyO(%rsp)
        movaps  %xmm2,nb312_fizO(%rsp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb312_jxO(%rsp),%xmm0
    movaps  nb312_jyO(%rsp),%xmm1
    movaps  nb312_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps  nb312_ixH1(%rsp),%xmm0
        subps  nb312_iyH1(%rsp),%xmm1
        subps  nb312_izH1(%rsp),%xmm2
        subps  nb312_ixH2(%rsp),%xmm3
        subps  nb312_iyH2(%rsp),%xmm4
        subps  nb312_izH2(%rsp),%xmm5
    movaps %xmm0,nb312_dxH1O(%rsp)
        movaps %xmm1,nb312_dyH1O(%rsp)
        movaps %xmm2,nb312_dzH1O(%rsp)
        movaps %xmm3,nb312_dxH2O(%rsp)
        movaps %xmm4,nb312_dyH2O(%rsp)
        movaps %xmm5,nb312_dzH2O(%rsp)
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
        movaps %xmm4,nb312_rsqH2O(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb312_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb312_rinvH2O(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb312_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

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
    mulps  nb312_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb312_qqOH(%rsp),%xmm3
        movhps  nb312_qqHH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 
    addps  nb312_vctot(%rsp),%xmm5
    movaps %xmm5,nb312_vctot(%rsp)

    xorps  %xmm1,%xmm1

    mulps nb312_tsc(%rsp),%xmm3
    mulps %xmm0,%xmm3
    subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb312_dxH1O(%rsp),%xmm0
        mulps   nb312_dyH1O(%rsp),%xmm1
        mulps   nb312_dzH1O(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb312_fjxO(%rsp),%xmm3
        movaps  nb312_fjyO(%rsp),%xmm4
        movaps  nb312_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb312_fjxO(%rsp)
        movaps  %xmm4,nb312_fjyO(%rsp)
        movaps  %xmm5,nb312_fjzO(%rsp)
        addps   nb312_fixH1(%rsp),%xmm0
        addps   nb312_fiyH1(%rsp),%xmm1
        addps   nb312_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb312_fixH1(%rsp)
        movaps  %xmm1,nb312_fiyH1(%rsp)
        movaps  %xmm2,nb312_fizH1(%rsp)
        ## do table for H2 - j water interaction 
        movaps nb312_rinvH2O(%rsp),%xmm0
        movaps nb312_rsqH2O(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb312_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

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
    mulps  nb312_two(%rsp),%xmm7         ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb312_qqOH(%rsp),%xmm3
        movhps  nb312_qqHH(%rsp),%xmm3

    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    mulps  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point xmm5 contains vcoul and xmm3 fijC 
    addps  nb312_vctot(%rsp),%xmm5
    movaps %xmm5,nb312_vctot(%rsp)

    xorps  %xmm1,%xmm1

    mulps nb312_tsc(%rsp),%xmm3
    mulps %xmm0,%xmm3
    subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb312_dxH2O(%rsp),%xmm0
        mulps   nb312_dyH2O(%rsp),%xmm1
        mulps   nb312_dzH2O(%rsp),%xmm2
        movaps  nb312_fjxO(%rsp),%xmm3
        movaps  nb312_fjyO(%rsp),%xmm4
        movaps  nb312_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb312_faction(%rbp),%rsi
        movaps  %xmm3,nb312_fjxO(%rsp)
        movaps  %xmm4,nb312_fjyO(%rsp)
        movaps  %xmm5,nb312_fjzO(%rsp)
        addps   nb312_fixH2(%rsp),%xmm0
        addps   nb312_fiyH2(%rsp),%xmm1
        addps   nb312_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb312_fixH2(%rsp)
        movaps  %xmm1,nb312_fiyH2(%rsp)
        movaps  %xmm2,nb312_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  (%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb312_fjxO(%rsp),%xmm3
        movaps  nb312_fjyO(%rsp),%xmm4
        movaps  nb312_fjzO(%rsp),%xmm5
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

        decl nb312_innerk(%rsp)
        jz    _nb_kernel312_x86_64_sse.nb312_updateouterdata
        jmp   _nb_kernel312_x86_64_sse.nb312_single_loop
_nb_kernel312_x86_64_sse.nb312_updateouterdata: 
        movl  nb312_ii3(%rsp),%ecx
        movq  nb312_faction(%rbp),%rdi
        movq  nb312_fshift(%rbp),%rsi
        movl  nb312_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb312_fixO(%rsp),%xmm0
        movaps nb312_fiyO(%rsp),%xmm1
        movaps nb312_fizO(%rsp),%xmm2

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
        movaps nb312_fixH1(%rsp),%xmm0
        movaps nb312_fiyH1(%rsp),%xmm1
        movaps nb312_fizH1(%rsp),%xmm2

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
        movaps nb312_fixH2(%rsp),%xmm0
        movaps nb312_fiyH2(%rsp),%xmm1
        movaps nb312_fizH2(%rsp),%xmm2

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
        movl nb312_n(%rsp),%esi
    ## get group index for i particle 
    movq  nb312_gid(%rbp),%rdx          ## base of gid[]
    movl  (%rdx,%rsi,4),%edx            ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb312_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb312_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb312_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb312_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb312_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel312_x86_64_sse.nb312_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb312_n(%rsp)
        jmp _nb_kernel312_x86_64_sse.nb312_outer
_nb_kernel312_x86_64_sse.nb312_outerend: 
        ## check if more outer neighborlists remain
        movl  nb312_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel312_x86_64_sse.nb312_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel312_x86_64_sse.nb312_threadloop
_nb_kernel312_x86_64_sse.nb312_end: 
        movl nb312_nouter(%rsp),%eax
        movl nb312_ninner(%rsp),%ebx
        movq nb312_outeriter(%rbp),%rcx
        movq nb312_inneriter(%rbp),%rdx
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





.globl nb_kernel312nf_x86_64_sse
.globl _nb_kernel312nf_x86_64_sse
nb_kernel312nf_x86_64_sse:      
_nb_kernel312nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb312nf_fshift, 16
.set nb312nf_gid, 24
.set nb312nf_pos, 32
.set nb312nf_faction, 40
.set nb312nf_charge, 48
.set nb312nf_p_facel, 56
.set nb312nf_argkrf, 64
.set nb312nf_argcrf, 72
.set nb312nf_Vc, 80
.set nb312nf_type, 88
.set nb312nf_p_ntype, 96
.set nb312nf_vdwparam, 104
.set nb312nf_Vvdw, 112
.set nb312nf_p_tabscale, 120
.set nb312nf_VFtab, 128
.set nb312nf_invsqrta, 136
.set nb312nf_dvda, 144
.set nb312nf_p_gbtabscale, 152
.set nb312nf_GBtab, 160
.set nb312nf_p_nthreads, 168
.set nb312nf_count, 176
.set nb312nf_mtx, 184
.set nb312nf_outeriter, 192
.set nb312nf_inneriter, 200
.set nb312nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb312nf_ixO, 0
.set nb312nf_iyO, 16
.set nb312nf_izO, 32
.set nb312nf_ixH1, 48
.set nb312nf_iyH1, 64
.set nb312nf_izH1, 80
.set nb312nf_ixH2, 96
.set nb312nf_iyH2, 112
.set nb312nf_izH2, 128
.set nb312nf_jxO, 144
.set nb312nf_jyO, 160
.set nb312nf_jzO, 176
.set nb312nf_jxH1, 192
.set nb312nf_jyH1, 208
.set nb312nf_jzH1, 224
.set nb312nf_jxH2, 240
.set nb312nf_jyH2, 256
.set nb312nf_jzH2, 272
.set nb312nf_qqOO, 288
.set nb312nf_qqOH, 304
.set nb312nf_qqHH, 320
.set nb312nf_tsc, 336
.set nb312nf_c6, 352
.set nb312nf_c12, 368
.set nb312nf_vctot, 384
.set nb312nf_Vvdwtot, 400
.set nb312nf_half, 416
.set nb312nf_three, 432
.set nb312nf_rsqOO, 448
.set nb312nf_rsqOH1, 464
.set nb312nf_rsqOH2, 480
.set nb312nf_rsqH1O, 496
.set nb312nf_rsqH1H1, 512
.set nb312nf_rsqH1H2, 528
.set nb312nf_rsqH2O, 544
.set nb312nf_rsqH2H1, 560
.set nb312nf_rsqH2H2, 576
.set nb312nf_rinvOO, 592
.set nb312nf_rinvOH1, 608
.set nb312nf_rinvOH2, 624
.set nb312nf_rinvH1O, 640
.set nb312nf_rinvH1H1, 656
.set nb312nf_rinvH1H2, 672
.set nb312nf_rinvH2O, 688
.set nb312nf_rinvH2H1, 704
.set nb312nf_rinvH2H2, 720
.set nb312nf_is3, 736
.set nb312nf_ii3, 740
.set nb312nf_nri, 744
.set nb312nf_iinr, 752
.set nb312nf_jindex, 760
.set nb312nf_jjnr, 768
.set nb312nf_shift, 776
.set nb312nf_shiftvec, 784
.set nb312nf_facel, 792
.set nb312nf_innerjjnr, 800
.set nb312nf_innerk, 808
.set nb312nf_n, 812
.set nb312nf_nn1, 816
.set nb312nf_nouter, 820
.set nb312nf_ninner, 824
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms
        subq $840,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb312nf_nouter(%rsp)
        movl %eax,nb312nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb312nf_nri(%rsp)
        movq %rsi,nb312nf_iinr(%rsp)
        movq %rdx,nb312nf_jindex(%rsp)
        movq %rcx,nb312nf_jjnr(%rsp)
        movq %r8,nb312nf_shift(%rsp)
        movq %r9,nb312nf_shiftvec(%rsp)
        movq nb312nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb312nf_facel(%rsp)

        movq nb312nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb312nf_tsc(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb312nf_half(%rsp)
        movss nb312nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb312nf_half(%rsp)
        movaps %xmm3,nb312nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb312nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb312nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%rdx,%rbx,4),%xmm5
        movq nb312nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb312nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb312nf_qqOO(%rsp)
        movaps %xmm4,nb312nf_qqOH(%rsp)
        movaps %xmm5,nb312nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb312nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb312nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb312nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## 01010101
        movaps %xmm0,nb312nf_c6(%rsp)
        movaps %xmm1,nb312nf_c12(%rsp)

_nb_kernel312nf_x86_64_sse.nb312nf_threadloop: 
        movq  nb312nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel312nf_x86_64_sse.nb312nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel312nf_x86_64_sse.nb312nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb312nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb312nf_n(%rsp)
        movl %ebx,nb312nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel312nf_x86_64_sse.nb312nf_outerstart
        jmp _nb_kernel312nf_x86_64_sse.nb312nf_end
_nb_kernel312nf_x86_64_sse.nb312nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb312nf_nouter(%rsp),%ebx
        movl %ebx,nb312nf_nouter(%rsp)

_nb_kernel312nf_x86_64_sse.nb312nf_outer: 
        movq  nb312nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb312nf_is3(%rsp)            ## store is3 

        movq  nb312nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb312nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb312nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb312nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb312nf_ixO(%rsp)
        movaps %xmm4,nb312nf_iyO(%rsp)
        movaps %xmm5,nb312nf_izO(%rsp)

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
        movaps %xmm0,nb312nf_ixH1(%rsp)
        movaps %xmm1,nb312nf_iyH1(%rsp)
        movaps %xmm2,nb312nf_izH1(%rsp)
        movaps %xmm3,nb312nf_ixH2(%rsp)
        movaps %xmm4,nb312nf_iyH2(%rsp)
        movaps %xmm5,nb312nf_izH2(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb312nf_vctot(%rsp)
        movaps %xmm4,nb312nf_Vvdwtot(%rsp)

        movq  nb312nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb312nf_pos(%rbp),%rsi
        movq  nb312nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb312nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb312nf_ninner(%rsp),%ecx
        movl  %ecx,nb312nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb312nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel312nf_x86_64_sse.nb312nf_unroll_loop
        jmp   _nb_kernel312nf_x86_64_sse.nb312nf_single_check
_nb_kernel312nf_x86_64_sse.nb312nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb312nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb312nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb312nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps %xmm0,nb312nf_jxO(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb312nf_jyO(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb312nf_jxH1(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb312nf_jyH1(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb312nf_jxH2(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb312nf_jyH2(%rsp)

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
        movaps %xmm0,nb312nf_jzO(%rsp)
        movaps %xmm1,nb312nf_jzH1(%rsp)
        movaps %xmm2,nb312nf_jzH2(%rsp)

        movaps nb312nf_ixO(%rsp),%xmm0
        movaps nb312nf_iyO(%rsp),%xmm1
        movaps nb312nf_izO(%rsp),%xmm2
        movaps nb312nf_ixO(%rsp),%xmm3
        movaps nb312nf_iyO(%rsp),%xmm4
        movaps nb312nf_izO(%rsp),%xmm5
        subps  nb312nf_jxO(%rsp),%xmm0
        subps  nb312nf_jyO(%rsp),%xmm1
        subps  nb312nf_jzO(%rsp),%xmm2
        subps  nb312nf_jxH1(%rsp),%xmm3
        subps  nb312nf_jyH1(%rsp),%xmm4
        subps  nb312nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb312nf_rsqOO(%rsp)
        movaps %xmm3,nb312nf_rsqOH1(%rsp)

        movaps nb312nf_ixO(%rsp),%xmm0
        movaps nb312nf_iyO(%rsp),%xmm1
        movaps nb312nf_izO(%rsp),%xmm2
        movaps nb312nf_ixH1(%rsp),%xmm3
        movaps nb312nf_iyH1(%rsp),%xmm4
        movaps nb312nf_izH1(%rsp),%xmm5
        subps  nb312nf_jxH2(%rsp),%xmm0
        subps  nb312nf_jyH2(%rsp),%xmm1
        subps  nb312nf_jzH2(%rsp),%xmm2
        subps  nb312nf_jxO(%rsp),%xmm3
        subps  nb312nf_jyO(%rsp),%xmm4
        subps  nb312nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb312nf_rsqOH2(%rsp)
        movaps %xmm3,nb312nf_rsqH1O(%rsp)

        movaps nb312nf_ixH1(%rsp),%xmm0
        movaps nb312nf_iyH1(%rsp),%xmm1
        movaps nb312nf_izH1(%rsp),%xmm2
        movaps nb312nf_ixH1(%rsp),%xmm3
        movaps nb312nf_iyH1(%rsp),%xmm4
        movaps nb312nf_izH1(%rsp),%xmm5
        subps  nb312nf_jxH1(%rsp),%xmm0
        subps  nb312nf_jyH1(%rsp),%xmm1
        subps  nb312nf_jzH1(%rsp),%xmm2
        subps  nb312nf_jxH2(%rsp),%xmm3
        subps  nb312nf_jyH2(%rsp),%xmm4
        subps  nb312nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb312nf_rsqH1H1(%rsp)
        movaps %xmm3,nb312nf_rsqH1H2(%rsp)

        movaps nb312nf_ixH2(%rsp),%xmm0
        movaps nb312nf_iyH2(%rsp),%xmm1
        movaps nb312nf_izH2(%rsp),%xmm2
        movaps nb312nf_ixH2(%rsp),%xmm3
        movaps nb312nf_iyH2(%rsp),%xmm4
        movaps nb312nf_izH2(%rsp),%xmm5
        subps  nb312nf_jxO(%rsp),%xmm0
        subps  nb312nf_jyO(%rsp),%xmm1
        subps  nb312nf_jzO(%rsp),%xmm2
        subps  nb312nf_jxH1(%rsp),%xmm3
        subps  nb312nf_jyH1(%rsp),%xmm4
        subps  nb312nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb312nf_rsqH2O(%rsp)
        movaps %xmm4,nb312nf_rsqH2H1(%rsp)

        movaps nb312nf_ixH2(%rsp),%xmm0
        movaps nb312nf_iyH2(%rsp),%xmm1
        movaps nb312nf_izH2(%rsp),%xmm2
        subps  nb312nf_jxH2(%rsp),%xmm0
        subps  nb312nf_jyH2(%rsp),%xmm1
        subps  nb312nf_jzH2(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb312nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb312nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb312nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb312nf_rinvH2H1(%rsp)

        rsqrtps nb312nf_rsqOO(%rsp),%xmm1
        rsqrtps nb312nf_rsqOH1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb312nf_rsqOO(%rsp),%xmm1
        mulps   nb312nf_rsqOH1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312nf_half(%rsp),%xmm3
        mulps   nb312nf_half(%rsp),%xmm7
        movaps  %xmm3,nb312nf_rinvOO(%rsp)
        movaps  %xmm7,nb312nf_rinvOH1(%rsp)

        rsqrtps nb312nf_rsqOH2(%rsp),%xmm1
        rsqrtps nb312nf_rsqH1O(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb312nf_rsqOH2(%rsp),%xmm1
        mulps   nb312nf_rsqH1O(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312nf_half(%rsp),%xmm3
        mulps   nb312nf_half(%rsp),%xmm7
        movaps  %xmm3,nb312nf_rinvOH2(%rsp)
        movaps  %xmm7,nb312nf_rinvH1O(%rsp)

        rsqrtps nb312nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb312nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb312nf_rsqH1H1(%rsp),%xmm1
        mulps   nb312nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312nf_half(%rsp),%xmm3
        mulps   nb312nf_half(%rsp),%xmm7
        movaps  %xmm3,nb312nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb312nf_rinvH1H2(%rsp)

        rsqrtps nb312nf_rsqH2O(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb312nf_three(%rsp),%xmm3
        mulps   nb312nf_rsqH2O(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb312nf_half(%rsp),%xmm3
        movaps  %xmm3,nb312nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movaps nb312nf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1

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

    movq nb312nf_VFtab(%rbp),%rsi
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
    movaps nb312nf_qqOO(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## start doing lj 
        movaps %xmm0,%xmm2
        mulps  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulps  %xmm2,%xmm1
        mulps  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb312nf_c6(%rsp),%xmm1
        mulps  nb312nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb312nf_Vvdwtot(%rsp),%xmm4
        movaps %xmm4,nb312nf_Vvdwtot(%rsp)

        ## O-H1 interaction 
        movaps nb312nf_rinvOH1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqOH1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## O-H2 interaction  
        movaps nb312nf_rinvOH2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqOH2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H1-O interaction 
        movaps nb312nf_rinvH1O(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH1O(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H1-H1 interaction 
        movaps nb312nf_rinvH1H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H1-H2 interaction 
        movaps nb312nf_rinvH1H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H2-O interaction 
        movaps nb312nf_rinvH2O(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH2O(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqOH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H2-H1 interaction 
        movaps nb312nf_rinvH2H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## H2-H2 interaction 
        movaps nb312nf_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb312nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1
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
    movaps nb312nf_qqHH(%rsp),%xmm3
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb312nf_innerk(%rsp)
        jl    _nb_kernel312nf_x86_64_sse.nb312nf_single_check
        jmp   _nb_kernel312nf_x86_64_sse.nb312nf_unroll_loop
_nb_kernel312nf_x86_64_sse.nb312nf_single_check: 
        addl $4,nb312nf_innerk(%rsp)
        jnz   _nb_kernel312nf_x86_64_sse.nb312nf_single_loop
        jmp   _nb_kernel312nf_x86_64_sse.nb312nf_updateouterdata
_nb_kernel312nf_x86_64_sse.nb312nf_single_loop: 
        movq  nb312nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb312nf_innerjjnr(%rsp)

        movq nb312nf_pos(%rbp),%rsi
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
        movaps  nb312nf_ixO(%rsp),%xmm0
        movaps  nb312nf_iyO(%rsp),%xmm1
        movaps  nb312nf_izO(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb312nf_jxO(%rsp)
        movaps %xmm4,nb312nf_jyO(%rsp)
        movaps %xmm5,nb312nf_jzO(%rsp)
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
        movaps  nb312nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb312nf_half(%rsp),%xmm3   ## rinv iO - j water 

        movaps  %xmm3,%xmm1
        mulps   %xmm0,%xmm1     ## xmm1=r 
        movaps  %xmm3,%xmm0     ## xmm0=rinv 
        mulps  nb312nf_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

        movq nb312nf_VFtab(%rbp),%rsi

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
        movss   nb312nf_qqOO(%rsp),%xmm3
        movhps  nb312nf_qqOH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point xmm5 contains vcoul 

    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## start doing lj 
        xorps  %xmm2,%xmm2
        movss  %xmm0,%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm2,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulss  nb312nf_c6(%rsp),%xmm1
        mulss  nb312nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm4
        subss  %xmm1,%xmm4
        addps  nb312nf_Vvdwtot(%rsp),%xmm4
        movaps %xmm4,nb312nf_Vvdwtot(%rsp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb312nf_ixH1(%rsp),%xmm0
        movaps  nb312nf_iyH1(%rsp),%xmm1
        movaps  nb312nf_izH1(%rsp),%xmm2
        movaps  nb312nf_ixH2(%rsp),%xmm3
        movaps  nb312nf_iyH2(%rsp),%xmm4
        movaps  nb312nf_izH2(%rsp),%xmm5
        subps   nb312nf_jxO(%rsp),%xmm0
        subps   nb312nf_jyO(%rsp),%xmm1
        subps   nb312nf_jzO(%rsp),%xmm2
        subps   nb312nf_jxO(%rsp),%xmm3
        subps   nb312nf_jyO(%rsp),%xmm4
        subps   nb312nf_jzO(%rsp),%xmm5
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
        movaps %xmm4,nb312nf_rsqH2O(%rsp)

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb312nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb312nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb312nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## start with H1, save H2 data 
        movaps %xmm7,nb312nf_rinvH2O(%rsp)

        movaps %xmm3,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb312nf_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

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
        movss   nb312nf_qqOH(%rsp),%xmm3
        movhps  nb312nf_qqHH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point xmm5 contains vcoul 
    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        ## do table for H2 - j water interaction 
        movaps nb312nf_rinvH2O(%rsp),%xmm0
        movaps nb312nf_rsqH2O(%rsp),%xmm1
        mulps  %xmm0,%xmm1      ## xmm0=rinv, xmm1=r 
        mulps  nb312nf_tsc(%rsp),%xmm1

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
    movd %mm7,%edx              ## table indices in ebx,ecx,edx 

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
        movss   nb312nf_qqOH(%rsp),%xmm3
        movhps  nb312nf_qqHH(%rsp),%xmm3

    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 
    mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point xmm5 contains vcoul 
    addps  nb312nf_vctot(%rsp),%xmm5
    movaps %xmm5,nb312nf_vctot(%rsp)

        decl nb312nf_innerk(%rsp)
        jz    _nb_kernel312nf_x86_64_sse.nb312nf_updateouterdata
        jmp   _nb_kernel312nf_x86_64_sse.nb312nf_single_loop
_nb_kernel312nf_x86_64_sse.nb312nf_updateouterdata: 
        ## get n from stack
        movl nb312nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb312nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb312nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb312nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb312nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb312nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb312nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel312nf_x86_64_sse.nb312nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb312nf_n(%rsp)
        jmp _nb_kernel312nf_x86_64_sse.nb312nf_outer
_nb_kernel312nf_x86_64_sse.nb312nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb312nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel312nf_x86_64_sse.nb312nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel312nf_x86_64_sse.nb312nf_threadloop
_nb_kernel312nf_x86_64_sse.nb312nf_end: 
        movl nb312nf_nouter(%rsp),%eax
        movl nb312nf_ninner(%rsp),%ebx
        movq nb312nf_outeriter(%rbp),%rcx
        movq nb312nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $840,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret



