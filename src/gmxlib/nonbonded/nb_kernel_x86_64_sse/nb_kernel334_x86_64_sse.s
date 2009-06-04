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





.globl nb_kernel334_x86_64_sse
.globl _nb_kernel334_x86_64_sse
nb_kernel334_x86_64_sse:        
_nb_kernel334_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb334_fshift, 16
.set nb334_gid, 24
.set nb334_pos, 32
.set nb334_faction, 40
.set nb334_charge, 48
.set nb334_p_facel, 56
.set nb334_argkrf, 64
.set nb334_argcrf, 72
.set nb334_Vc, 80
.set nb334_type, 88
.set nb334_p_ntype, 96
.set nb334_vdwparam, 104
.set nb334_Vvdw, 112
.set nb334_p_tabscale, 120
.set nb334_VFtab, 128
.set nb334_invsqrta, 136
.set nb334_dvda, 144
.set nb334_p_gbtabscale, 152
.set nb334_GBtab, 160
.set nb334_p_nthreads, 168
.set nb334_count, 176
.set nb334_mtx, 184
.set nb334_outeriter, 192
.set nb334_inneriter, 200
.set nb334_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb334_ixO, 0
.set nb334_iyO, 16
.set nb334_izO, 32
.set nb334_ixH1, 48
.set nb334_iyH1, 64
.set nb334_izH1, 80
.set nb334_ixH2, 96
.set nb334_iyH2, 112
.set nb334_izH2, 128
.set nb334_ixM, 144
.set nb334_iyM, 160
.set nb334_izM, 176
.set nb334_jxO, 192
.set nb334_jyO, 208
.set nb334_jzO, 224
.set nb334_jxH1, 240
.set nb334_jyH1, 256
.set nb334_jzH1, 272
.set nb334_jxH2, 288
.set nb334_jyH2, 304
.set nb334_jzH2, 320
.set nb334_jxM, 336
.set nb334_jyM, 352
.set nb334_jzM, 368
.set nb334_dxOO, 384
.set nb334_dyOO, 400
.set nb334_dzOO, 416
.set nb334_dxH1H1, 432
.set nb334_dyH1H1, 448
.set nb334_dzH1H1, 464
.set nb334_dxH1H2, 480
.set nb334_dyH1H2, 496
.set nb334_dzH1H2, 512
.set nb334_dxH1M, 528
.set nb334_dyH1M, 544
.set nb334_dzH1M, 560
.set nb334_dxH2H1, 576
.set nb334_dyH2H1, 592
.set nb334_dzH2H1, 608
.set nb334_dxH2H2, 624
.set nb334_dyH2H2, 640
.set nb334_dzH2H2, 656
.set nb334_dxH2M, 672
.set nb334_dyH2M, 688
.set nb334_dzH2M, 704
.set nb334_dxMH1, 720
.set nb334_dyMH1, 736
.set nb334_dzMH1, 752
.set nb334_dxMH2, 768
.set nb334_dyMH2, 784
.set nb334_dzMH2, 800
.set nb334_dxMM, 816
.set nb334_dyMM, 832
.set nb334_dzMM, 848
.set nb334_qqMM, 864
.set nb334_qqMH, 880
.set nb334_qqHH, 896
.set nb334_two, 912
.set nb334_tsc, 928
.set nb334_c6, 944
.set nb334_c12, 960
.set nb334_vctot, 976
.set nb334_Vvdwtot, 992
.set nb334_fixO, 1008
.set nb334_fiyO, 1024
.set nb334_fizO, 1040
.set nb334_fixH1, 1056
.set nb334_fiyH1, 1072
.set nb334_fizH1, 1088
.set nb334_fixH2, 1104
.set nb334_fiyH2, 1120
.set nb334_fizH2, 1136
.set nb334_fixM, 1152
.set nb334_fiyM, 1168
.set nb334_fizM, 1184
.set nb334_fjxO, 1200
.set nb334_fjyO, 1216
.set nb334_fjzO, 1232
.set nb334_epsH1, 1248
.set nb334_epsH2, 1264
.set nb334_epsM, 1280
.set nb334_fjxH2, 1296
.set nb334_fjyH2, 1312
.set nb334_fjzH2, 1328
.set nb334_fjxM, 1344
.set nb334_fjyM, 1360
.set nb334_fjzM, 1376
.set nb334_half, 1392
.set nb334_three, 1408
.set nb334_rsqOO, 1424
.set nb334_rsqH1H1, 1440
.set nb334_rsqH1H2, 1456
.set nb334_rsqH1M, 1472
.set nb334_rsqH2H1, 1488
.set nb334_rsqH2H2, 1504
.set nb334_rsqH2M, 1520
.set nb334_rsqMH1, 1536
.set nb334_rsqMH2, 1552
.set nb334_rsqMM, 1568
.set nb334_rinvOO, 1584
.set nb334_rinvH1H1, 1600
.set nb334_rinvH1H2, 1616
.set nb334_rinvH1M, 1632
.set nb334_rinvH2H1, 1648
.set nb334_rinvH2H2, 1664
.set nb334_rinvH2M, 1680
.set nb334_rinvMH1, 1696
.set nb334_rinvMH2, 1712
.set nb334_rinvMM, 1728
.set nb334_fstmp, 1744
.set nb334_is3, 1760
.set nb334_ii3, 1764
.set nb334_nri, 1768
.set nb334_iinr, 1776
.set nb334_jindex, 1784
.set nb334_jjnr, 1792
.set nb334_shift, 1800
.set nb334_shiftvec, 1808
.set nb334_facel, 1816
.set nb334_innerjjnr, 1824
.set nb334_innerk, 1832
.set nb334_n, 1836
.set nb334_nn1, 1840
.set nb334_nouter, 1844
.set nb334_ninner, 1848
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1864,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb334_nouter(%rsp)
        movl %eax,nb334_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb334_nri(%rsp)
        movq %rsi,nb334_iinr(%rsp)
        movq %rdx,nb334_jindex(%rsp)
        movq %rcx,nb334_jjnr(%rsp)
        movq %r8,nb334_shift(%rsp)
        movq %r9,nb334_shiftvec(%rsp)
        movq nb334_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb334_facel(%rsp)

        movq nb334_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb334_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb334_half(%rsp)
        movss nb334_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb334_half(%rsp)
        movaps %xmm2,nb334_two(%rsp)
        movaps %xmm3,nb334_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb334_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb334_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb334_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb334_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb334_qqMM(%rsp)
        movaps %xmm4,nb334_qqMH(%rsp)
        movaps %xmm5,nb334_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb334_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb334_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx       ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb334_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb334_c6(%rsp)
        movaps %xmm1,nb334_c12(%rsp)

_nb_kernel334_x86_64_sse.nb334_threadloop: 
        movq  nb334_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel334_x86_64_sse.nb334_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334_x86_64_sse.nb334_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb334_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb334_n(%rsp)
        movl %ebx,nb334_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel334_x86_64_sse.nb334_outerstart
        jmp _nb_kernel334_x86_64_sse.nb334_end

_nb_kernel334_x86_64_sse.nb334_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334_nouter(%rsp),%ebx
        movl %ebx,nb334_nouter(%rsp)

_nb_kernel334_x86_64_sse.nb334_outer: 
        movq  nb334_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb334_is3(%rsp)      ## store is3 

        movq  nb334_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb334_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb334_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb334_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb334_ixO(%rsp)
        movaps %xmm4,nb334_iyO(%rsp)
        movaps %xmm5,nb334_izO(%rsp)
        movaps %xmm6,nb334_ixH1(%rsp)
        movaps %xmm7,nb334_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb334_izH1(%rsp)
        movaps %xmm0,nb334_ixH2(%rsp)
        movaps %xmm1,nb334_iyH2(%rsp)
        movaps %xmm2,nb334_izH2(%rsp)
        movaps %xmm3,nb334_ixM(%rsp)
        movaps %xmm4,nb334_iyM(%rsp)
        movaps %xmm5,nb334_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb334_vctot(%rsp)
        movaps %xmm4,nb334_Vvdwtot(%rsp)
        movaps %xmm4,nb334_fixO(%rsp)
        movaps %xmm4,nb334_fiyO(%rsp)
        movaps %xmm4,nb334_fizO(%rsp)
        movaps %xmm4,nb334_fixH1(%rsp)
        movaps %xmm4,nb334_fiyH1(%rsp)
        movaps %xmm4,nb334_fizH1(%rsp)
        movaps %xmm4,nb334_fixH2(%rsp)
        movaps %xmm4,nb334_fiyH2(%rsp)
        movaps %xmm4,nb334_fizH2(%rsp)
        movaps %xmm4,nb334_fixM(%rsp)
        movaps %xmm4,nb334_fiyM(%rsp)
        movaps %xmm4,nb334_fizM(%rsp)

        movq  nb334_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb334_pos(%rbp),%rsi
        movq  nb334_faction(%rbp),%rdi
        movq  nb334_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb334_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb334_ninner(%rsp),%ecx
        movl  %ecx,nb334_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb334_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel334_x86_64_sse.nb334_unroll_loop
        jmp   _nb_kernel334_x86_64_sse.nb334_single_check
_nb_kernel334_x86_64_sse.nb334_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb334_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb334_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb334_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
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

    subps nb334_ixO(%rsp),%xmm0
    subps nb334_iyO(%rsp),%xmm1
    subps nb334_izO(%rsp),%xmm2

    ## store dx/dy/dz
    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    movaps %xmm2,%xmm15

    ## square it
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2

        addps  %xmm0,%xmm1
        addps  %xmm2,%xmm1
    ## rsq in xmm1

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb334_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4
        mulps nb334_half(%rsp),%xmm4
        movaps %xmm4,%xmm2
        mulps  %xmm4,%xmm1
    ## xmm2=rinv
    ## xmm1=r

    mulps nb334_tsc(%rsp),%xmm1   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm1,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm4

    ## multiply by 4
    pslld   $2,%xmm5

    ## multiply by three (copy, mult. by two, add back)
    movaps  %xmm5,%xmm6
    pslld   $1,%xmm5
    paddd   %xmm6,%xmm5

    ## calculate eps
    subps     %xmm4,%xmm1

    ## move to integer registers
    movhlps %xmm5,%xmm6
    movd    %xmm5,%r8d
    movd    %xmm6,%r10d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm6,%xmm6
    movd    %xmm5,%r9d
    movd    %xmm6,%r11d
    ## table indices in r8-r11

    ## xmm1=eps
    ## xmm2=rinv

        movq nb334_VFtab(%rbp),%rsi
    ## load LJ dispersion and repulsion in parallel 
    ## NB: We are using a combined (LJ+coul) table, 
    ## so the LJ table data is offset 4*4 = 16 bytes
    movlps 16(%rsi,%r8,4),%xmm5
        movlps 32(%rsi,%r8,4),%xmm9
        movlps 16(%rsi,%r10,4),%xmm7
        movlps 32(%rsi,%r10,4),%xmm11
        movhps 16(%rsi,%r9,4),%xmm5
        movhps 32(%rsi,%r9,4),%xmm9
        movhps 16(%rsi,%r11,4),%xmm7
        movhps 32(%rsi,%r11,4),%xmm11

    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 24(%rsi,%r8,4),%xmm7
        movlps 40(%rsi,%r8,4),%xmm11
        movlps 24(%rsi,%r10,4),%xmm0
        movlps 40(%rsi,%r10,4),%xmm3
        movhps 24(%rsi,%r9,4),%xmm7
        movhps 40(%rsi,%r9,4),%xmm11
        movhps 24(%rsi,%r11,4),%xmm0
        movhps 40(%rsi,%r11,4),%xmm3

    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm0,%xmm6 ## 10001000
        shufps $136,%xmm3,%xmm10 ## 10001000
        shufps $221,%xmm0,%xmm7 ## 11011101
        shufps $221,%xmm3,%xmm11 ## 11011101
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulps  %xmm1,%xmm7   ## Heps
    mulps  %xmm1,%xmm11
    mulps  %xmm1,%xmm6  ## Geps
    mulps  %xmm1,%xmm10
    mulps  %xmm1,%xmm7  ## Heps2
    mulps  %xmm1,%xmm11
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
    mulps  %xmm1,%xmm5 ## eps*Fp
    mulps  %xmm1,%xmm9
    addps  %xmm4,%xmm5 ## VV
    addps  %xmm8,%xmm9

    mulps  nb334_c6(%rsp),%xmm5    ## VV*c6 = vnb6
    mulps  nb334_c12(%rsp),%xmm9    ## VV*c12 = vnb12
    addps  %xmm9,%xmm5
    addps  nb334_Vvdwtot(%rsp),%xmm5
    movaps %xmm5,nb334_Vvdwtot(%rsp)

    mulps  nb334_c6(%rsp),%xmm7     ## FF*c6 = fnb6
    mulps  nb334_c12(%rsp),%xmm11     ## FF*c12  = fnb12
    addps  %xmm11,%xmm7

    mulps  nb334_tsc(%rsp),%xmm7
    mulps  %xmm2,%xmm7
    xorps  %xmm9,%xmm9

    subps  %xmm7,%xmm9

    ## fx/fy/fz
    mulps  %xmm9,%xmm13
    mulps  %xmm9,%xmm14
    mulps  %xmm9,%xmm15

    ## increment i force
    movaps nb334_fixO(%rsp),%xmm0
    movaps nb334_fiyO(%rsp),%xmm1
    movaps nb334_fizO(%rsp),%xmm2
    addps  %xmm13,%xmm0
    addps  %xmm14,%xmm1
    addps  %xmm15,%xmm2
    movaps %xmm0,nb334_fixO(%rsp)
    movaps %xmm1,nb334_fiyO(%rsp)
    movaps %xmm2,nb334_fizO(%rsp)

        ## move j O forces to local temp variables 
    movq nb334_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm4 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm5 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm4 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm5 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm6    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm7    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm6    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm7    ## jzOc  -  jzOd -

    shufps $136,%xmm7,%xmm6 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm4: jxOa jyOa jxOb jyOb 
    ## xmm5: jxOc jyOc jxOd jyOd
    ## xmm6: jzOa jzOb jzOc jzOd

    ## update O forces
    movaps %xmm13,%xmm12
    unpcklps %xmm14,%xmm13  ## (local) fjx1 fjx1 fjy1 fjy2
    unpckhps %xmm14,%xmm12  ## (local) fjx3 fjx4 fjy3 fjy4

    addps %xmm13,%xmm4
    addps %xmm12,%xmm5
    addps %xmm15,%xmm6

    movhlps  %xmm6,%xmm7 ## fH1zc fH1zd

    movlps %xmm4,(%rdi,%rax,4)
    movhps %xmm4,(%rdi,%rbx,4)
    movlps %xmm5,(%rdi,%rcx,4)
    movhps %xmm5,(%rdi,%rdx,4)
    movss  %xmm6,8(%rdi,%rax,4)
    movss  %xmm7,8(%rdi,%rcx,4)
    shufps $1,%xmm6,%xmm6
    shufps $1,%xmm7,%xmm7
    movss  %xmm6,8(%rdi,%rbx,4)
    movss  %xmm7,8(%rdi,%rdx,4)
    ## done with OO interaction

    ## move j H1 coordinates to local temp variables 
    movq nb334_pos(%rbp),%rsi
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

    subps nb334_ixH1(%rsp),%xmm0
    subps nb334_iyH1(%rsp),%xmm1
    subps nb334_izH1(%rsp),%xmm2
    subps nb334_ixH2(%rsp),%xmm3
    subps nb334_iyH2(%rsp),%xmm4
    subps nb334_izH2(%rsp),%xmm5
    subps nb334_ixM(%rsp),%xmm6
    subps nb334_iyM(%rsp),%xmm7
    subps nb334_izM(%rsp),%xmm8

        movaps %xmm0,nb334_dxH1H1(%rsp)
        movaps %xmm1,nb334_dyH1H1(%rsp)
        movaps %xmm2,nb334_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH2H1(%rsp)
        movaps %xmm4,nb334_dyH2H1(%rsp)
        movaps %xmm5,nb334_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb334_dxMH1(%rsp)
        movaps %xmm7,nb334_dyMH1(%rsp)
        movaps %xmm8,nb334_dzMH1(%rsp)
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

        movaps  nb334_three(%rsp),%xmm9
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

        movaps  nb334_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H1 
        mulps   %xmm4,%xmm10 ## rinvH2H1
    mulps   %xmm4,%xmm11 ## rinvMH1

        movaps  %xmm9,nb334_rinvH1H1(%rsp)
        movaps  %xmm10,nb334_rinvH2H1(%rsp)
        movaps  %xmm11,nb334_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb334_tsc(%rsp),%xmm1
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

    movq nb334_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb334_epsH1(%rsp)
    movaps    %xmm3,nb334_epsH2(%rsp)
    movaps    %xmm6,nb334_epsM(%rsp)

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

    movaps nb334_epsH1(%rsp),%xmm12
    movaps nb334_epsH2(%rsp),%xmm13
    movaps nb334_epsM(%rsp),%xmm14

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
    movaps nb334_qqHH(%rsp),%xmm12
    movaps nb334_qqMH(%rsp),%xmm13
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
    addps  nb334_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb334_vctot(%rsp)

    movaps nb334_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j H1 forces to local temp variables 
    movq nb334_faction(%rbp),%rdi
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

    mulps  nb334_rinvH1H1(%rsp),%xmm3
    mulps  nb334_rinvH2H1(%rsp),%xmm7
    mulps  nb334_rinvMH1(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb334_dxH1H1(%rsp),%xmm0
        mulps nb334_dyH1H1(%rsp),%xmm1
        mulps nb334_dzH1H1(%rsp),%xmm2
        mulps nb334_dxH2H1(%rsp),%xmm3
        mulps nb334_dyH2H1(%rsp),%xmm4
        mulps nb334_dzH2H1(%rsp),%xmm5
        mulps nb334_dxMH1(%rsp),%xmm6
        mulps nb334_dyMH1(%rsp),%xmm7
        mulps nb334_dzMH1(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb334_fixH1(%rsp),%xmm0
    addps nb334_fiyH1(%rsp),%xmm1
    addps nb334_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb334_fixH2(%rsp),%xmm3
    addps nb334_fiyH2(%rsp),%xmm4
    addps nb334_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb334_fixM(%rsp),%xmm6
    addps nb334_fiyM(%rsp),%xmm7
    addps nb334_fizM(%rsp),%xmm8

    movaps %xmm0,nb334_fixH1(%rsp)
    movaps %xmm1,nb334_fiyH1(%rsp)
    movaps %xmm2,nb334_fizH1(%rsp)
    movaps %xmm3,nb334_fixH2(%rsp)
    movaps %xmm4,nb334_fiyH2(%rsp)
    movaps %xmm5,nb334_fizH2(%rsp)
    movaps %xmm6,nb334_fixM(%rsp)
    movaps %xmm7,nb334_fiyM(%rsp)
    movaps %xmm8,nb334_fizM(%rsp)

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
        movq  nb334_pos(%rbp),%rsi
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

    subps nb334_ixH1(%rsp),%xmm0
    subps nb334_iyH1(%rsp),%xmm1
    subps nb334_izH1(%rsp),%xmm2
    subps nb334_ixH2(%rsp),%xmm3
    subps nb334_iyH2(%rsp),%xmm4
    subps nb334_izH2(%rsp),%xmm5
    subps nb334_ixM(%rsp),%xmm6
    subps nb334_iyM(%rsp),%xmm7
    subps nb334_izM(%rsp),%xmm8

        movaps %xmm0,nb334_dxH1H2(%rsp)
        movaps %xmm1,nb334_dyH1H2(%rsp)
        movaps %xmm2,nb334_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH2H2(%rsp)
        movaps %xmm4,nb334_dyH2H2(%rsp)
        movaps %xmm5,nb334_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb334_dxMH2(%rsp)
        movaps %xmm7,nb334_dyMH2(%rsp)
        movaps %xmm8,nb334_dzMH2(%rsp)
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

        movaps  nb334_three(%rsp),%xmm9
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

        movaps  nb334_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H2
        mulps   %xmm4,%xmm10 ## rinvH2H2
    mulps   %xmm4,%xmm11 ## rinvMH2

        movaps  %xmm9,nb334_rinvH1H2(%rsp)
        movaps  %xmm10,nb334_rinvH2H2(%rsp)
        movaps  %xmm11,nb334_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb334_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    movq nb334_VFtab(%rbp),%rsi

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

    movaps    %xmm0,nb334_epsH1(%rsp)
    movaps    %xmm3,nb334_epsH2(%rsp)
    movaps    %xmm6,nb334_epsM(%rsp)


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

    movaps nb334_epsH1(%rsp),%xmm12
    movaps nb334_epsH2(%rsp),%xmm13
    movaps nb334_epsM(%rsp),%xmm14

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
    movaps nb334_qqHH(%rsp),%xmm12
    movaps nb334_qqMH(%rsp),%xmm13
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
    addps  nb334_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb334_vctot(%rsp)

    movaps nb334_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j H2 forces to local temp variables 
    movq nb334_faction(%rbp),%rdi
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

    mulps  nb334_rinvH1H2(%rsp),%xmm3
    mulps  nb334_rinvH2H2(%rsp),%xmm7
    mulps  nb334_rinvMH2(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb334_dxH1H2(%rsp),%xmm0
        mulps nb334_dyH1H2(%rsp),%xmm1
        mulps nb334_dzH1H2(%rsp),%xmm2
        mulps nb334_dxH2H2(%rsp),%xmm3
        mulps nb334_dyH2H2(%rsp),%xmm4
        mulps nb334_dzH2H2(%rsp),%xmm5
        mulps nb334_dxMH2(%rsp),%xmm6
        mulps nb334_dyMH2(%rsp),%xmm7
        mulps nb334_dzMH2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb334_fixH1(%rsp),%xmm0
    addps nb334_fiyH1(%rsp),%xmm1
    addps nb334_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb334_fixH2(%rsp),%xmm3
    addps nb334_fiyH2(%rsp),%xmm4
    addps nb334_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb334_fixM(%rsp),%xmm6
    addps nb334_fiyM(%rsp),%xmm7
    addps nb334_fizM(%rsp),%xmm8

    movaps %xmm0,nb334_fixH1(%rsp)
    movaps %xmm1,nb334_fiyH1(%rsp)
    movaps %xmm2,nb334_fizH1(%rsp)
    movaps %xmm3,nb334_fixH2(%rsp)
    movaps %xmm4,nb334_fiyH2(%rsp)
    movaps %xmm5,nb334_fizH2(%rsp)
    movaps %xmm6,nb334_fixM(%rsp)
    movaps %xmm7,nb334_fiyM(%rsp)
    movaps %xmm8,nb334_fizM(%rsp)

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

        movq  nb334_pos(%rbp),%rsi
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

    subps nb334_ixH1(%rsp),%xmm0
    subps nb334_iyH1(%rsp),%xmm1
    subps nb334_izH1(%rsp),%xmm2
    subps nb334_ixH2(%rsp),%xmm3
    subps nb334_iyH2(%rsp),%xmm4
    subps nb334_izH2(%rsp),%xmm5
    subps nb334_ixM(%rsp),%xmm6
    subps nb334_iyM(%rsp),%xmm7
    subps nb334_izM(%rsp),%xmm8

        movaps %xmm0,nb334_dxH1M(%rsp)
        movaps %xmm1,nb334_dyH1M(%rsp)
        movaps %xmm2,nb334_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb334_dxH2M(%rsp)
        movaps %xmm4,nb334_dyH2M(%rsp)
        movaps %xmm5,nb334_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb334_dxMM(%rsp)
        movaps %xmm7,nb334_dyMM(%rsp)
        movaps %xmm8,nb334_dzMM(%rsp)
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

        movaps  nb334_three(%rsp),%xmm9
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

        movaps  nb334_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1M
        mulps   %xmm4,%xmm10 ## rinvH2M
    mulps   %xmm4,%xmm11 ## rinvMM

        movaps  %xmm9,nb334_rinvH1M(%rsp)
        movaps  %xmm10,nb334_rinvH2M(%rsp)
        movaps  %xmm11,nb334_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb334_tsc(%rsp),%xmm1
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

    movq nb334_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb334_epsH1(%rsp)
    movaps    %xmm3,nb334_epsH2(%rsp)
    movaps    %xmm6,nb334_epsM(%rsp)

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

    movaps nb334_epsH1(%rsp),%xmm12
    movaps nb334_epsH2(%rsp),%xmm13
    movaps nb334_epsM(%rsp),%xmm14

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
    movaps nb334_qqMH(%rsp),%xmm12
    movaps nb334_qqMM(%rsp),%xmm13
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
    addps  nb334_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb334_vctot(%rsp)

    movaps nb334_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax ## restore j3 from mm0-mm3
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j M forces to local temp variables 
    movq nb334_faction(%rbp),%rdi
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

    mulps  nb334_rinvH1M(%rsp),%xmm3
    mulps  nb334_rinvH2M(%rsp),%xmm7
    mulps  nb334_rinvMM(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb334_dxH1M(%rsp),%xmm0
        mulps nb334_dyH1M(%rsp),%xmm1
        mulps nb334_dzH1M(%rsp),%xmm2
        mulps nb334_dxH2M(%rsp),%xmm3
        mulps nb334_dyH2M(%rsp),%xmm4
        mulps nb334_dzH2M(%rsp),%xmm5
        mulps nb334_dxMM(%rsp),%xmm6
        mulps nb334_dyMM(%rsp),%xmm7
        mulps nb334_dzMM(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb334_fixH1(%rsp),%xmm0
    addps nb334_fiyH1(%rsp),%xmm1
    addps nb334_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb334_fixH2(%rsp),%xmm3
    addps nb334_fiyH2(%rsp),%xmm4
    addps nb334_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb334_fixM(%rsp),%xmm6
    addps nb334_fiyM(%rsp),%xmm7
    addps nb334_fizM(%rsp),%xmm8

    movaps %xmm0,nb334_fixH1(%rsp)
    movaps %xmm1,nb334_fiyH1(%rsp)
    movaps %xmm2,nb334_fizH1(%rsp)
    movaps %xmm3,nb334_fixH2(%rsp)
    movaps %xmm4,nb334_fiyH2(%rsp)
    movaps %xmm5,nb334_fizH2(%rsp)
    movaps %xmm6,nb334_fixM(%rsp)
    movaps %xmm7,nb334_fiyM(%rsp)
    movaps %xmm8,nb334_fizM(%rsp)

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
        subl $4,nb334_innerk(%rsp)
        jl    _nb_kernel334_x86_64_sse.nb334_single_check
        jmp   _nb_kernel334_x86_64_sse.nb334_unroll_loop
_nb_kernel334_x86_64_sse.nb334_single_check: 
        addl $4,nb334_innerk(%rsp)
        jnz   _nb_kernel334_x86_64_sse.nb334_single_loop
        jmp   _nb_kernel334_x86_64_sse.nb334_updateouterdata
_nb_kernel334_x86_64_sse.nb334_single_loop: 
        movq  nb334_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb334_innerjjnr(%rsp)

        movq nb334_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates
        movlps (%rsi,%rax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%rsi,%rax,4),%xmm4            ## H1y H1z 
        movlps 32(%rsi,%rax,4),%xmm5            ## H2z  Mx 
        movhps 8(%rsi,%rax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%rsi,%rax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%rsi,%rax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb334_jxO(%rsp)
        movaps %xmm3,nb334_jyO(%rsp)
        movaps %xmm1,nb334_jzO(%rsp)

    movaps %xmm6,%xmm0  ## jxO
    movaps %xmm1,%xmm2  ## jzO
    movaps %xmm3,%xmm1  ## jyO
    movaps %xmm3,%xmm4  ## jyO
    movaps %xmm6,%xmm3  ## jxO
    movaps %xmm2,%xmm5  ## jzO

        ## do O and H1 in parallel
        subps  nb334_ixO(%rsp),%xmm0
        subps  nb334_iyO(%rsp),%xmm1
        subps  nb334_izO(%rsp),%xmm2
        subps  nb334_ixH1(%rsp),%xmm3
        subps  nb334_iyH1(%rsp),%xmm4
        subps  nb334_izH1(%rsp),%xmm5

        movaps %xmm0,nb334_dxOO(%rsp)
        movaps %xmm1,nb334_dyOO(%rsp)
        movaps %xmm2,nb334_dzOO(%rsp)
        movaps %xmm3,nb334_dxH1H1(%rsp)
        movaps %xmm4,nb334_dyH1H1(%rsp)
        movaps %xmm5,nb334_dzH1H1(%rsp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        movaps %xmm0,nb334_rsqOO(%rsp)
        movaps %xmm4,nb334_rsqH1H1(%rsp)

        ## do 1/sqrt(x) for O and  H1
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%rsp),%xmm3   ## rinv O - j water 
        mulps   nb334_half(%rsp),%xmm7   ## rinv H1 - j water  

        movaps %xmm3,nb334_rinvOO(%rsp)
        movaps %xmm7,nb334_rinvH1H1(%rsp)

        movq nb334_VFtab(%rbp),%rsi

        ## do O table LJ interaction
        movaps %xmm3,%xmm0
        movaps %xmm0,%xmm1
        mulss  nb334_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulss  nb334_tsc(%rsp),%xmm1

        cvttps2pi %xmm1,%mm6
        cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6

        movd %mm6,%ebx
        lea  (%rbx,%rbx,2),%rbx

        ## load dispersion table data into xmm4
        movlps 16(%rsi,%rbx,4),%xmm4
        movlps 24(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb334_two(%rsp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm7      ## fijD 
        mulss  %xmm4,%xmm5      ## Vvdw6 

        ## save scalar force in xmm3. Update Vvdwtot directly 
        addss  nb334_Vvdwtot(%rsp),%xmm5
        movaps %xmm7,%xmm3 ## fscal 
        movss %xmm5,nb334_Vvdwtot(%rsp)

        ## load repulsion table data into xmm4
        movlps 32(%rsi,%rbx,4),%xmm4
        movlps 40(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  nb334_two(%rsp),%xmm7            ## two*Heps2 
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm7 ## fijR 
        mulss  %xmm4,%xmm5 ## Vvdw12 
        addss  %xmm3,%xmm7

        addss  nb334_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb334_Vvdwtot(%rsp)

        xorps  %xmm1,%xmm1
        mulss nb334_tsc(%rsp),%xmm7
        mulss %xmm0,%xmm7
        subss  %xmm7,%xmm1

        movaps %xmm1,%xmm0
        movaps %xmm1,%xmm2

        mulss  nb334_dxOO(%rsp),%xmm0
        mulss  nb334_dyOO(%rsp),%xmm1
        mulss  nb334_dzOO(%rsp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addss   %xmm0,%xmm3
        addss   %xmm1,%xmm4
        addss   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%rsp)
        movaps  %xmm4,nb334_fjyO(%rsp)
        movaps  %xmm5,nb334_fjzO(%rsp)
        addss   nb334_fixO(%rsp),%xmm0
        addss   nb334_fiyO(%rsp),%xmm1
        addss   nb334_fizO(%rsp),%xmm2
        movss  %xmm0,nb334_fixO(%rsp)
        movss  %xmm1,nb334_fiyO(%rsp)
        movss  %xmm2,nb334_fizO(%rsp)

        ## do  H1 coulomb interaction
        movaps nb334_rinvH1H1(%rsp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb334_rsqH1H1(%rsp),%xmm1        ## r
        mulps nb334_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movq nb334_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqHH(%rsp),%xmm3
        movhps  nb334_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb334_vctot(%rsp),%xmm5
        movaps %xmm5,nb334_vctot(%rsp)

        mulps  nb334_tsc(%rsp),%xmm3
        xorps  %xmm2,%xmm2
        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb334_dxH1H1(%rsp),%xmm0
        mulps   nb334_dyH1H1(%rsp),%xmm1
        mulps   nb334_dzH1H1(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb334_fjxO(%rsp),%xmm3
        movaps  nb334_fjyO(%rsp),%xmm4
        movaps  nb334_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%rsp)
        movaps  %xmm4,nb334_fjyO(%rsp)
        movaps  %xmm5,nb334_fjzO(%rsp)
        addps   nb334_fixH1(%rsp),%xmm0
        addps   nb334_fiyH1(%rsp),%xmm1
        addps   nb334_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb334_fixH1(%rsp)
        movaps  %xmm1,nb334_fiyH1(%rsp)
        movaps  %xmm2,nb334_fizH1(%rsp)

        ## i H2 & M simultaneously first get i particle coords: 
    movaps  nb334_jxO(%rsp),%xmm0
    movaps  nb334_jyO(%rsp),%xmm1
    movaps  nb334_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps   nb334_ixH2(%rsp),%xmm0
        subps   nb334_iyH2(%rsp),%xmm1
        subps   nb334_izH2(%rsp),%xmm2
        subps   nb334_ixM(%rsp),%xmm3
        subps   nb334_iyM(%rsp),%xmm4
        subps   nb334_izM(%rsp),%xmm5
        movaps %xmm0,nb334_dxH2H2(%rsp)
        movaps %xmm1,nb334_dyH2H2(%rsp)
        movaps %xmm2,nb334_dzH2H2(%rsp)
        movaps %xmm3,nb334_dxMM(%rsp)
        movaps %xmm4,nb334_dyMM(%rsp)
        movaps %xmm5,nb334_dzMM(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb334_rsqH2H2(%rsp)
        movaps %xmm4,nb334_rsqMM(%rsp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334_half(%rsp),%xmm3   ## rinv H2 - j water 
        mulps   nb334_half(%rsp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb334_rinvH2H2(%rsp)
        movaps %xmm7,nb334_rinvMM(%rsp)

        movaps %xmm3,%xmm1
        mulps  nb334_rsqH2H2(%rsp),%xmm1        ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb334_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqHH(%rsp),%xmm3
        movhps  nb334_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb334_vctot(%rsp),%xmm5
        movaps %xmm5,nb334_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb334_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb334_dxH2H2(%rsp),%xmm0
        mulps   nb334_dyH2H2(%rsp),%xmm1
        mulps   nb334_dzH2H2(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb334_fjxO(%rsp),%xmm3
        movaps  nb334_fjyO(%rsp),%xmm4
        movaps  nb334_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb334_fjxO(%rsp)
        movaps  %xmm4,nb334_fjyO(%rsp)
        movaps  %xmm5,nb334_fjzO(%rsp)
        addps   nb334_fixH2(%rsp),%xmm0
        addps   nb334_fiyH2(%rsp),%xmm1
        addps   nb334_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb334_fixH2(%rsp)
        movaps  %xmm1,nb334_fiyH2(%rsp)
        movaps  %xmm2,nb334_fizH2(%rsp)

        ## do table for i M - j water interaction 
        movaps nb334_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334_rsqMM(%rsp),%xmm1          ## xmm0=rinv, xmm1=r 
        mulps  nb334_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movq nb334_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## # coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  nb334_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334_qqMH(%rsp),%xmm3
        movhps  nb334_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb334_vctot(%rsp),%xmm5
        movaps %xmm5,nb334_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb334_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb334_dxMM(%rsp),%xmm0
        mulps   nb334_dyMM(%rsp),%xmm1
        mulps   nb334_dzMM(%rsp),%xmm2
        movaps  nb334_fjxO(%rsp),%xmm3
        movaps  nb334_fjyO(%rsp),%xmm4
        movaps  nb334_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb334_faction(%rbp),%rsi
        movaps  %xmm3,nb334_fjxO(%rsp)
        movaps  %xmm4,nb334_fjyO(%rsp)
        movaps  %xmm5,nb334_fjzO(%rsp)
        addps   nb334_fixM(%rsp),%xmm0
        addps   nb334_fiyM(%rsp),%xmm1
        addps   nb334_fizM(%rsp),%xmm2
        movaps  %xmm0,nb334_fixM(%rsp)
        movaps  %xmm1,nb334_fiyM(%rsp)
        movaps  %xmm2,nb334_fizM(%rsp)

        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb334_fjxO(%rsp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb334_fjyO(%rsp),%xmm1   ## Oy H1y H2y My
        movaps  nb334_fjzO(%rsp),%xmm2   ## Oz H1z H2z Mz

        movaps  %xmm0,%xmm3
        movaps  %xmm0,%xmm4
        unpcklps %xmm1,%xmm3            ## Ox Oy - -
        shufps $0x1,%xmm2,%xmm4        ## h1x - Oz -
        movaps  %xmm1,%xmm5
        movaps  %xmm0,%xmm6
        unpcklps %xmm2,%xmm5            ## - - H1y H1z
        unpckhps %xmm1,%xmm6            ## h2x h2y - - 
        unpckhps %xmm2,%xmm1            ## - - My Mz

        shufps  $0x32,%xmm0,%xmm2 ## (00110010) h2z - Mx -
        shufps  $36,%xmm4,%xmm3 ## 00100100 ;# Ox Oy Oz H1x 
        shufps  $78,%xmm6,%xmm5 ## 01001110 ;# h1y h1z h2x h2y
        shufps  $232,%xmm1,%xmm2 ## 11101000 ;# h2z mx my mz

        movlps  (%rsi,%rax,4),%xmm0
        movlps  16(%rsi,%rax,4),%xmm1
        movlps  32(%rsi,%rax,4),%xmm4
        movhps  8(%rsi,%rax,4),%xmm0
        movhps  24(%rsi,%rax,4),%xmm1
        movhps  40(%rsi,%rax,4),%xmm4
        addps   %xmm3,%xmm0
        addps   %xmm5,%xmm1
        addps   %xmm2,%xmm4
        movlps   %xmm0,(%rsi,%rax,4)
        movlps   %xmm1,16(%rsi,%rax,4)
        movlps   %xmm4,32(%rsi,%rax,4)
        movhps   %xmm0,8(%rsi,%rax,4)
        movhps   %xmm1,24(%rsi,%rax,4)
        movhps   %xmm4,40(%rsi,%rax,4)

        decl nb334_innerk(%rsp)
        jz    _nb_kernel334_x86_64_sse.nb334_updateouterdata
        jmp   _nb_kernel334_x86_64_sse.nb334_single_loop
_nb_kernel334_x86_64_sse.nb334_updateouterdata: 
        movl  nb334_ii3(%rsp),%ecx
        movq  nb334_faction(%rbp),%rdi
        movq  nb334_fshift(%rbp),%rsi
        movl  nb334_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb334_fixO(%rsp),%xmm0
        movaps nb334_fiyO(%rsp),%xmm1
        movaps nb334_fizO(%rsp),%xmm2

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
        movaps nb334_fixH1(%rsp),%xmm0
        movaps nb334_fiyH1(%rsp),%xmm1
        movaps nb334_fizH1(%rsp),%xmm2

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
        movaps nb334_fixH2(%rsp),%xmm0
        movaps nb334_fiyH2(%rsp),%xmm1
        movaps nb334_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb334_fixM(%rsp),%xmm0
        movaps nb334_fiyM(%rsp),%xmm1
        movaps nb334_fizM(%rsp),%xmm2

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
        movl nb334_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb334_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb334_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb334_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb334_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb334_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb334_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334_x86_64_sse.nb334_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334_n(%rsp)
        jmp _nb_kernel334_x86_64_sse.nb334_outer
_nb_kernel334_x86_64_sse.nb334_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334_x86_64_sse.nb334_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334_x86_64_sse.nb334_threadloop
_nb_kernel334_x86_64_sse.nb334_end: 
        movl nb334_nouter(%rsp),%eax
        movl nb334_ninner(%rsp),%ebx
        movq nb334_outeriter(%rbp),%rcx
        movq nb334_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1864,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel334nf_x86_64_sse
.globl _nb_kernel334nf_x86_64_sse
nb_kernel334nf_x86_64_sse:      
_nb_kernel334nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb334nf_fshift, 16
.set nb334nf_gid, 24
.set nb334nf_pos, 32
.set nb334nf_faction, 40
.set nb334nf_charge, 48
.set nb334nf_p_facel, 56
.set nb334nf_argkrf, 64
.set nb334nf_argcrf, 72
.set nb334nf_Vc, 80
.set nb334nf_type, 88
.set nb334nf_p_ntype, 96
.set nb334nf_vdwparam, 104
.set nb334nf_Vvdw, 112
.set nb334nf_p_tabscale, 120
.set nb334nf_VFtab, 128
.set nb334nf_invsqrta, 136
.set nb334nf_dvda, 144
.set nb334nf_p_gbtabscale, 152
.set nb334nf_GBtab, 160
.set nb334nf_p_nthreads, 168
.set nb334nf_count, 176
.set nb334nf_mtx, 184
.set nb334nf_outeriter, 192
.set nb334nf_inneriter, 200
.set nb334nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb334nf_ixO, 0
.set nb334nf_iyO, 16
.set nb334nf_izO, 32
.set nb334nf_ixH1, 48
.set nb334nf_iyH1, 64
.set nb334nf_izH1, 80
.set nb334nf_ixH2, 96
.set nb334nf_iyH2, 112
.set nb334nf_izH2, 128
.set nb334nf_ixM, 144
.set nb334nf_iyM, 160
.set nb334nf_izM, 176
.set nb334nf_jxO, 192
.set nb334nf_jyO, 208
.set nb334nf_jzO, 224
.set nb334nf_jxH1, 240
.set nb334nf_jyH1, 256
.set nb334nf_jzH1, 272
.set nb334nf_jxH2, 288
.set nb334nf_jyH2, 304
.set nb334nf_jzH2, 320
.set nb334nf_jxM, 336
.set nb334nf_jyM, 352
.set nb334nf_jzM, 368
.set nb334nf_qqMM, 384
.set nb334nf_qqMH, 400
.set nb334nf_qqHH, 416
.set nb334nf_tsc, 432
.set nb334nf_c6, 448
.set nb334nf_c12, 464
.set nb334nf_vctot, 480
.set nb334nf_Vvdwtot, 496
.set nb334nf_half, 512
.set nb334nf_three, 528
.set nb334nf_rsqOO, 544
.set nb334nf_rsqH1H1, 560
.set nb334nf_rsqH1H2, 576
.set nb334nf_rsqH1M, 592
.set nb334nf_rsqH2H1, 608
.set nb334nf_rsqH2H2, 704
.set nb334nf_rsqH2M, 720
.set nb334nf_rsqMH1, 736
.set nb334nf_rsqMH2, 752
.set nb334nf_rsqMM, 768
.set nb334nf_rinvOO, 784
.set nb334nf_rinvH1H1, 800
.set nb334nf_rinvH1H2, 816
.set nb334nf_rinvH1M, 832
.set nb334nf_rinvH2H1, 848
.set nb334nf_rinvH2H2, 864
.set nb334nf_rinvH2M, 880
.set nb334nf_rinvMH1, 896
.set nb334nf_rinvMH2, 912
.set nb334nf_rinvMM, 928
.set nb334nf_is3, 944
.set nb334nf_ii3, 948
.set nb334nf_nri, 952
.set nb334nf_iinr, 960
.set nb334nf_jindex, 968
.set nb334nf_jjnr, 976
.set nb334nf_shift, 984
.set nb334nf_shiftvec, 992
.set nb334nf_facel, 1000
.set nb334nf_innerjjnr, 1008
.set nb334nf_innerk, 1016
.set nb334nf_n, 1020
.set nb334nf_nn1, 1024
.set nb334nf_nouter, 1028
.set nb334nf_ninner, 1032
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1048,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb334nf_nouter(%rsp)
        movl %eax,nb334nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb334nf_nri(%rsp)
        movq %rsi,nb334nf_iinr(%rsp)
        movq %rdx,nb334nf_jindex(%rsp)
        movq %rcx,nb334nf_jjnr(%rsp)
        movq %r8,nb334nf_shift(%rsp)
        movq %r9,nb334nf_shiftvec(%rsp)
        movq nb334nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb334nf_facel(%rsp)

        movq nb334nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb334nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb334nf_half(%rsp)
        movss nb334nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb334nf_half(%rsp)
        movaps %xmm3,nb334nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb334nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb334nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb334nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb334nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb334nf_qqMM(%rsp)
        movaps %xmm4,nb334nf_qqMH(%rsp)
        movaps %xmm5,nb334nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb334nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb334nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx       ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb334nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb334nf_c6(%rsp)
        movaps %xmm1,nb334nf_c12(%rsp)

_nb_kernel334nf_x86_64_sse.nb334nf_threadloop: 
        movq  nb334nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel334nf_x86_64_sse.nb334nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel334nf_x86_64_sse.nb334nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb334nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb334nf_n(%rsp)
        movl %ebx,nb334nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel334nf_x86_64_sse.nb334nf_outerstart
        jmp _nb_kernel334nf_x86_64_sse.nb334nf_end

_nb_kernel334nf_x86_64_sse.nb334nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb334nf_nouter(%rsp),%ebx
        movl %ebx,nb334nf_nouter(%rsp)

_nb_kernel334nf_x86_64_sse.nb334nf_outer: 
        movq  nb334nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb334nf_is3(%rsp)            ## store is3 

        movq  nb334nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb334nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb334nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb334nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb334nf_ixO(%rsp)
        movaps %xmm4,nb334nf_iyO(%rsp)
        movaps %xmm5,nb334nf_izO(%rsp)
        movaps %xmm6,nb334nf_ixH1(%rsp)
        movaps %xmm7,nb334nf_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb334nf_izH1(%rsp)
        movaps %xmm0,nb334nf_ixH2(%rsp)
        movaps %xmm1,nb334nf_iyH2(%rsp)
        movaps %xmm2,nb334nf_izH2(%rsp)
        movaps %xmm3,nb334nf_ixM(%rsp)
        movaps %xmm4,nb334nf_iyM(%rsp)
        movaps %xmm5,nb334nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb334nf_vctot(%rsp)
        movaps %xmm4,nb334nf_Vvdwtot(%rsp)

        movq  nb334nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb334nf_pos(%rbp),%rsi
        movq  nb334nf_faction(%rbp),%rdi
        movq  nb334nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb334nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb334nf_ninner(%rsp),%ecx
        movl  %ecx,nb334nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb334nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel334nf_x86_64_sse.nb334nf_unroll_loop
        jmp   _nb_kernel334nf_x86_64_sse.nb334nf_single_check
_nb_kernel334nf_x86_64_sse.nb334nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb334nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb334nf_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb334nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%rsi,%rax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%rsi,%rcx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%rsi,%rbx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%rsi,%rdx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%rsi,%rax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%rsi,%rcx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%rsi,%rbx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%rsi,%rdx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jxO(%rsp)
        movaps %xmm1,nb334nf_jyO(%rsp)
        movaps %xmm2,nb334nf_jzO(%rsp)
        movaps %xmm3,nb334nf_jxH1(%rsp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%rsi,%rax,4),%xmm1
        movlps 16(%rsi,%rcx,4),%xmm4
        movhps 16(%rsi,%rbx,4),%xmm1
        movhps 16(%rsi,%rdx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%rsi,%rax,4),%xmm3
        movlps 24(%rsi,%rcx,4),%xmm5
        movhps 24(%rsi,%rbx,4),%xmm3
        movhps 24(%rsi,%rdx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jyH1(%rsp)
        movaps %xmm1,nb334nf_jzH1(%rsp)
        movaps %xmm2,nb334nf_jxH2(%rsp)
        movaps %xmm3,nb334nf_jyH2(%rsp)

        ## Load H2z Mx My Mz 
        movlps 32(%rsi,%rax,4),%xmm1
        movlps 32(%rsi,%rcx,4),%xmm4
        movhps 32(%rsi,%rbx,4),%xmm1
        movhps 32(%rsi,%rdx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%rsi,%rax,4),%xmm3
        movlps 40(%rsi,%rcx,4),%xmm5
        movhps 40(%rsi,%rbx,4),%xmm3
        movhps 40(%rsi,%rdx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb334nf_jzH2(%rsp)
        movaps %xmm1,nb334nf_jxM(%rsp)
        movaps %xmm2,nb334nf_jyM(%rsp)
        movaps %xmm3,nb334nf_jzM(%rsp)

        ## start calculating pairwise distances
        movaps nb334nf_ixO(%rsp),%xmm0
        movaps nb334nf_iyO(%rsp),%xmm1
        movaps nb334nf_izO(%rsp),%xmm2
        movaps nb334nf_ixH1(%rsp),%xmm3
        movaps nb334nf_iyH1(%rsp),%xmm4
        movaps nb334nf_izH1(%rsp),%xmm5
        subps  nb334nf_jxO(%rsp),%xmm0
        subps  nb334nf_jyO(%rsp),%xmm1
        subps  nb334nf_jzO(%rsp),%xmm2
        subps  nb334nf_jxH1(%rsp),%xmm3
        subps  nb334nf_jyH1(%rsp),%xmm4
        subps  nb334nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb334nf_rsqOO(%rsp)
        movaps %xmm3,nb334nf_rsqH1H1(%rsp)

        movaps nb334nf_ixH1(%rsp),%xmm0
        movaps nb334nf_iyH1(%rsp),%xmm1
        movaps nb334nf_izH1(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH2(%rsp),%xmm0
        subps  nb334nf_jyH2(%rsp),%xmm1
        subps  nb334nf_jzH2(%rsp),%xmm2
        subps  nb334nf_jxM(%rsp),%xmm3
        subps  nb334nf_jyM(%rsp),%xmm4
        subps  nb334nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb334nf_rsqH1H2(%rsp)
        movaps %xmm3,nb334nf_rsqH1M(%rsp)

        movaps nb334nf_ixH2(%rsp),%xmm0
        movaps nb334nf_iyH2(%rsp),%xmm1
        movaps nb334nf_izH2(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH1(%rsp),%xmm0
        subps  nb334nf_jyH1(%rsp),%xmm1
        subps  nb334nf_jzH1(%rsp),%xmm2
        subps  nb334nf_jxH2(%rsp),%xmm3
        subps  nb334nf_jyH2(%rsp),%xmm4
        subps  nb334nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb334nf_rsqH2H1(%rsp)
        movaps %xmm3,nb334nf_rsqH2H2(%rsp)

        movaps nb334nf_ixH2(%rsp),%xmm0
        movaps nb334nf_iyH2(%rsp),%xmm1
        movaps nb334nf_izH2(%rsp),%xmm2
        movaps nb334nf_ixM(%rsp),%xmm3
        movaps nb334nf_iyM(%rsp),%xmm4
        movaps nb334nf_izM(%rsp),%xmm5
        subps  nb334nf_jxM(%rsp),%xmm0
        subps  nb334nf_jyM(%rsp),%xmm1
        subps  nb334nf_jzM(%rsp),%xmm2
        subps  nb334nf_jxH1(%rsp),%xmm3
        subps  nb334nf_jyH1(%rsp),%xmm4
        subps  nb334nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb334nf_rsqH2M(%rsp)
        movaps %xmm4,nb334nf_rsqMH1(%rsp)

        movaps nb334nf_ixM(%rsp),%xmm0
        movaps nb334nf_iyM(%rsp),%xmm1
        movaps nb334nf_izM(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb334nf_jxH2(%rsp),%xmm0
        subps  nb334nf_jyH2(%rsp),%xmm1
        subps  nb334nf_jzH2(%rsp),%xmm2
        subps  nb334nf_jxM(%rsp),%xmm3
        subps  nb334nf_jyM(%rsp),%xmm4
        subps  nb334nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb334nf_rsqMH2(%rsp)
        movaps %xmm4,nb334nf_rsqMM(%rsp)

        ## Invsqrt for O-O
        rsqrtps  nb334nf_rsqOO(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334nf_three(%rsp),%xmm3
        mulps   nb334nf_rsqOO(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334nf_half(%rsp),%xmm3   ## rinvOO
        movaps %xmm3,nb334nf_rinvOO(%rsp)

        ## Invsqrt for H1-H1 and H1-H2
        rsqrtps nb334nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb334nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH1H1(%rsp),%xmm1
        mulps   nb334nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3   ## rinvH1H1 
        mulps   nb334nf_half(%rsp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb334nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb334nf_rinvH1H2(%rsp)

        rsqrtps nb334nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb334nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH1M(%rsp),%xmm1
        mulps   nb334nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3
        mulps   nb334nf_half(%rsp),%xmm7
        movaps  %xmm3,nb334nf_rinvH1M(%rsp)
        movaps  %xmm7,nb334nf_rinvH2H1(%rsp)

        rsqrtps nb334nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb334nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqH2H2(%rsp),%xmm1
        mulps   nb334nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3
        mulps   nb334nf_half(%rsp),%xmm7
        movaps  %xmm3,nb334nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb334nf_rinvH2M(%rsp)

        rsqrtps nb334nf_rsqMH1(%rsp),%xmm1
        rsqrtps nb334nf_rsqMH2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb334nf_rsqMH1(%rsp),%xmm1
        mulps   nb334nf_rsqMH2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3
        mulps   nb334nf_half(%rsp),%xmm7
        movaps  %xmm3,nb334nf_rinvMH1(%rsp)
        movaps  %xmm7,nb334nf_rinvMH2(%rsp)

        rsqrtps nb334nf_rsqMM(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb334nf_three(%rsp),%xmm3
        mulps   nb334nf_rsqMM(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb334nf_half(%rsp),%xmm3
        movaps  %xmm3,nb334nf_rinvMM(%rsp)

        ## start with OO table interaction
        movaps nb334nf_rinvOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqOO(%rsp),%xmm1   ## xmm1=r
        mulps  nb334nf_tsc(%rsp),%xmm1

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

        movq nb334nf_VFtab(%rbp),%rsi
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

        ## load dispersion table data into xmm4-xmm7
        movlps 16(%rsi,%rax,4),%xmm5
        movlps 16(%rsi,%rcx,4),%xmm7
        movhps 16(%rsi,%rbx,4),%xmm5
        movhps 16(%rsi,%rdx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 24(%rsi,%rax,4),%xmm7
        movlps 24(%rsi,%rcx,4),%xmm3
        movhps 24(%rsi,%rbx,4),%xmm7
        movhps 24(%rsi,%rdx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## dispersion table YFGH ready in xmm4-xmm7
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c6(%rsp),%xmm4
        mulps  %xmm4,%xmm5      ## Vvdw6 

        ## Update Vvdwtot directly 
        addps  nb334nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb334nf_Vvdwtot(%rsp)

        ## load repulsion table data into xmm4-xmm7
        movlps 32(%rsi,%rax,4),%xmm5
        movlps 32(%rsi,%rcx,4),%xmm7
        movhps 32(%rsi,%rbx,4),%xmm5
        movhps 32(%rsi,%rdx,4),%xmm7    ## got half table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 40(%rsi,%rax,4),%xmm7
        movlps 40(%rsi,%rcx,4),%xmm3
        movhps 40(%rsi,%rbx,4),%xmm7
        movhps 40(%rsi,%rdx,4),%xmm3    ## other half of table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## repulsion table YFGH ready in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c12(%rsp),%xmm4
        mulps  %xmm4,%xmm5 ## Vvdw12 

        addps  nb334nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb334nf_Vvdwtot(%rsp)

        ## Coulomb interactions - first H1H1
        movaps nb334nf_rinvH1H1(%rsp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1

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

        movq nb334nf_VFtab(%rbp),%rsi
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
        movaps nb334nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## H1-H2 interaction 
        movaps nb334nf_rinvH1H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## H1-M interaction  
        movaps nb334nf_rinvH1M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## H2-H1 interaction 
        movaps nb334nf_rinvH2H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## H2-H2 interaction 
        movaps nb334nf_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## H2-M interaction 
        movaps nb334nf_rinvH2M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## M-H1 interaction 
        movaps nb334nf_rinvMH1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## M-H2 interaction 
        movaps nb334nf_rinvMH2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## M-M interaction 
        movaps nb334nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1
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
        movaps nb334nf_qqMM(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)
        ## should we do one more iteration? 
        subl $4,nb334nf_innerk(%rsp)
        jl    _nb_kernel334nf_x86_64_sse.nb334nf_single_check
        jmp   _nb_kernel334nf_x86_64_sse.nb334nf_unroll_loop
_nb_kernel334nf_x86_64_sse.nb334nf_single_check: 
        addl $4,nb334nf_innerk(%rsp)
        jnz   _nb_kernel334nf_x86_64_sse.nb334nf_single_loop
        jmp   _nb_kernel334nf_x86_64_sse.nb334nf_updateouterdata
_nb_kernel334nf_x86_64_sse.nb334nf_single_loop: 
        movq  nb334nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb334nf_innerjjnr(%rsp)

        movq nb334nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates
        movlps (%rsi,%rax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%rsi,%rax,4),%xmm4            ## H1y H1z 
        movlps 32(%rsi,%rax,4),%xmm5            ## H2z  Mx 
        movhps 8(%rsi,%rax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%rsi,%rax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%rsi,%rax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb334nf_jxO(%rsp)
        movaps %xmm3,nb334nf_jyO(%rsp)
        movaps %xmm1,nb334nf_jzO(%rsp)

        ## do O and H1 in parallel
        movaps nb334nf_ixO(%rsp),%xmm0
        movaps nb334nf_iyO(%rsp),%xmm1
        movaps nb334nf_izO(%rsp),%xmm2
        movaps nb334nf_ixH1(%rsp),%xmm3
        movaps nb334nf_iyH1(%rsp),%xmm4
        movaps nb334nf_izH1(%rsp),%xmm5
        subps  nb334nf_jxO(%rsp),%xmm0
        subps  nb334nf_jyO(%rsp),%xmm1
        subps  nb334nf_jzO(%rsp),%xmm2
        subps  nb334nf_jxO(%rsp),%xmm3
        subps  nb334nf_jyO(%rsp),%xmm4
        subps  nb334nf_jzO(%rsp),%xmm5

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        movaps %xmm0,nb334nf_rsqOO(%rsp)
        movaps %xmm4,nb334nf_rsqH1H1(%rsp)

        ## do 1/sqrt(x) for O and  H1
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3   ## rinv O - j water 
        mulps   nb334nf_half(%rsp),%xmm7   ## rinv H1 - j water  

        movaps %xmm3,nb334nf_rinvOO(%rsp)
        movaps %xmm7,nb334nf_rinvH1H1(%rsp)

        movq nb334nf_VFtab(%rbp),%rsi

        ## do O table LJ interaction
        movaps %xmm3,%xmm0
        movaps %xmm0,%xmm1
        mulss  nb334nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulss  nb334nf_tsc(%rsp),%xmm1

        cvttps2pi %xmm1,%mm6
        cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 
        pslld   $2,%mm6

        movd %mm6,%ebx
        lea  (%rbx,%rbx,2),%rbx

        ## load dispersion table data into xmm4
        movlps 16(%rsi,%rbx,4),%xmm4
        movlps 24(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## dispersion table YFGH ready in xmm4-xmm7
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c6(%rsp),%xmm4
        mulss  %xmm4,%xmm5      ## Vvdw6 

        ## Update Vvdwtot directly 
        addss  nb334nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb334nf_Vvdwtot(%rsp)

        ## load repulsion table data into xmm4
        movlps 32(%rsi,%rbx,4),%xmm4
        movlps 40(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm7,%xmm7
        ## repulsion table YFGH ready in xmm4-xmm7

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 

        movaps nb334nf_c12(%rsp),%xmm4
        mulss  %xmm4,%xmm5 ## Vvdw12 

        addss  nb334nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb334nf_Vvdwtot(%rsp)

        ## do  H1 coulomb interaction
        movaps nb334nf_rinvH1H1(%rsp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqH1H1(%rsp),%xmm1      ## r
        mulps nb334nf_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movq nb334_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqHH(%rsp),%xmm3
        movhps  nb334nf_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb334nf_ixH2(%rsp),%xmm0
        movaps  nb334nf_iyH2(%rsp),%xmm1
        movaps  nb334nf_izH2(%rsp),%xmm2
        movaps  nb334nf_ixM(%rsp),%xmm3
        movaps  nb334nf_iyM(%rsp),%xmm4
        movaps  nb334nf_izM(%rsp),%xmm5
        subps   nb334nf_jxO(%rsp),%xmm0
        subps   nb334nf_jyO(%rsp),%xmm1
        subps   nb334nf_jzO(%rsp),%xmm2
        subps   nb334nf_jxO(%rsp),%xmm3
        subps   nb334nf_jyO(%rsp),%xmm4
        subps   nb334nf_jzO(%rsp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH2 in xmm0 
        addps %xmm5,%xmm4       ## have rsqM in xmm4 

        ## start with H2, save data 
        movaps %xmm0,nb334nf_rsqH2H2(%rsp)
        movaps %xmm4,nb334nf_rsqMM(%rsp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb334nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb334nf_half(%rsp),%xmm3   ## rinv H2 - j water 
        mulps   nb334nf_half(%rsp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb334nf_rinvH2H2(%rsp)
        movaps %xmm7,nb334nf_rinvMM(%rsp)

        movaps %xmm3,%xmm1
        mulps  nb334nf_rsqH2H2(%rsp),%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb334nf_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqHH(%rsp),%xmm3
        movhps  nb334nf_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        ## do table for i M - j water interaction 
        movaps nb334nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb334nf_rsqMM(%rsp),%xmm1        ## xmm0=rinv, xmm1=r 
        mulps  nb334nf_tsc(%rsp),%xmm1

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

        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx          ## table indices in ebx,ecx,edx 

        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx
        lea  (%rdx,%rdx,2),%rdx

        movq nb334_VFtab(%rbp),%rsi

        movlps (%rsi,%rbx,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%rdx,4),%xmm7
        movhps 8(%rsi,%rbx,4),%xmm4
        movhps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm3,%xmm6
        unpcklps %xmm7,%xmm6
        unpckhps %xmm7,%xmm3
        movaps %xmm4,%xmm5
        movaps %xmm4,%xmm7
        shufps $0x40,%xmm6,%xmm4
        shufps $0xE4,%xmm6,%xmm5
        movaps %xmm7,%xmm6
        shufps $0x48,%xmm3,%xmm6
        shufps $0xEC,%xmm3,%xmm7
        ## # coulomb table ready, in xmm4-xmm7

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb334nf_qqMH(%rsp),%xmm3
        movhps  nb334nf_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul
        addps  nb334nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb334nf_vctot(%rsp)

        decl nb334nf_innerk(%rsp)
        jz    _nb_kernel334nf_x86_64_sse.nb334nf_updateouterdata
        jmp   _nb_kernel334nf_x86_64_sse.nb334nf_single_loop
_nb_kernel334nf_x86_64_sse.nb334nf_updateouterdata: 
        ## get n from stack
        movl nb334nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb334nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb334nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb334nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb334nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb334nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb334nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel334nf_x86_64_sse.nb334nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb334nf_n(%rsp)
        jmp _nb_kernel334nf_x86_64_sse.nb334nf_outer
_nb_kernel334nf_x86_64_sse.nb334nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb334nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel334nf_x86_64_sse.nb334nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel334nf_x86_64_sse.nb334nf_threadloop
_nb_kernel334nf_x86_64_sse.nb334nf_end: 
        movl nb334nf_nouter(%rsp),%eax
        movl nb334nf_ninner(%rsp),%ebx
        movq nb334nf_outeriter(%rbp),%rcx
        movq nb334nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1048,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


