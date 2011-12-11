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






.globl nb_kernel314_x86_64_sse
.globl _nb_kernel314_x86_64_sse
nb_kernel314_x86_64_sse:        
_nb_kernel314_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb314_fshift, 16
.set nb314_gid, 24
.set nb314_pos, 32
.set nb314_faction, 40
.set nb314_charge, 48
.set nb314_p_facel, 56
.set nb314_argkrf, 64
.set nb314_argcrf, 72
.set nb314_Vc, 80
.set nb314_type, 88
.set nb314_p_ntype, 96
.set nb314_vdwparam, 104
.set nb314_Vvdw, 112
.set nb314_p_tabscale, 120
.set nb314_VFtab, 128
.set nb314_invsqrta, 136
.set nb314_dvda, 144
.set nb314_p_gbtabscale, 152
.set nb314_GBtab, 160
.set nb314_p_nthreads, 168
.set nb314_count, 176
.set nb314_mtx, 184
.set nb314_outeriter, 192
.set nb314_inneriter, 200
.set nb314_work, 208
        ## bottom of stack is cache-aligned for sse use 
.set nb314_ixO, 0
.set nb314_iyO, 16
.set nb314_izO, 32
.set nb314_ixH1, 48
.set nb314_iyH1, 64
.set nb314_izH1, 80
.set nb314_ixH2, 96
.set nb314_iyH2, 112
.set nb314_izH2, 128
.set nb314_ixM, 144
.set nb314_iyM, 160
.set nb314_izM, 176
.set nb314_jxO, 192
.set nb314_jyO, 208
.set nb314_jzO, 224
.set nb314_jxH1, 240
.set nb314_jyH1, 256
.set nb314_jzH1, 272
.set nb314_jxH2, 288
.set nb314_jyH2, 304
.set nb314_jzH2, 320
.set nb314_jxM, 336
.set nb314_jyM, 352
.set nb314_jzM, 368
.set nb314_dxOO, 384
.set nb314_dyOO, 400
.set nb314_dzOO, 416
.set nb314_dxH1H1, 432
.set nb314_dyH1H1, 448
.set nb314_dzH1H1, 464
.set nb314_dxH1H2, 480
.set nb314_dyH1H2, 496
.set nb314_dzH1H2, 512
.set nb314_dxH1M, 528
.set nb314_dyH1M, 544
.set nb314_dzH1M, 560
.set nb314_dxH2H1, 576
.set nb314_dyH2H1, 592
.set nb314_dzH2H1, 608
.set nb314_dxH2H2, 624
.set nb314_dyH2H2, 640
.set nb314_dzH2H2, 656
.set nb314_dxH2M, 672
.set nb314_dyH2M, 688
.set nb314_dzH2M, 704
.set nb314_dxMH1, 720
.set nb314_dyMH1, 736
.set nb314_dzMH1, 752
.set nb314_dxMH2, 768
.set nb314_dyMH2, 784
.set nb314_dzMH2, 800
.set nb314_dxMM, 816
.set nb314_dyMM, 832
.set nb314_dzMM, 848
.set nb314_qqMM, 864
.set nb314_qqMH, 880
.set nb314_qqHH, 896
.set nb314_two, 912
.set nb314_tsc, 928
.set nb314_c6, 944
.set nb314_c12, 960
.set nb314_six, 976
.set nb314_twelve, 992
.set nb314_vctot, 1008
.set nb314_Vvdwtot, 1024
.set nb314_fixO, 1040
.set nb314_fiyO, 1056
.set nb314_fizO, 1072
.set nb314_fixH1, 1088
.set nb314_fiyH1, 1104
.set nb314_fizH1, 1120
.set nb314_fixH2, 1136
.set nb314_fiyH2, 1152
.set nb314_fizH2, 1168
.set nb314_fixM, 1184
.set nb314_fiyM, 1200
.set nb314_fizM, 1216
.set nb314_fjxO, 1232
.set nb314_fjyO, 1248
.set nb314_fjzO, 1264
.set nb314_fjxH1, 1280
.set nb314_fjyH1, 1296
.set nb314_fjzH1, 1312
.set nb314_fjxH2, 1328
.set nb314_fjyH2, 1344
.set nb314_fjzH2, 1360
.set nb314_epsH1, 1376
.set nb314_epsH2, 1392
.set nb314_epsM, 1408
.set nb314_half, 1424
.set nb314_three, 1440
.set nb314_rsqOO, 1456
.set nb314_rsqH1H1, 1472
.set nb314_rsqH1H2, 1488
.set nb314_rsqH1M, 1504
.set nb314_rsqH2H1, 1520
.set nb314_rsqH2H2, 1536
.set nb314_rsqH2M, 1552
.set nb314_rsqMH1, 1568
.set nb314_rsqMH2, 1584
.set nb314_rsqMM, 1600
.set nb314_rinvsqOO, 1616
.set nb314_rinvH1H1, 1632
.set nb314_rinvH1H2, 1648
.set nb314_rinvH1M, 1664
.set nb314_rinvH2H1, 1680
.set nb314_rinvH2H2, 1696
.set nb314_rinvH2M, 1712
.set nb314_rinvMH1, 1728
.set nb314_rinvMH2, 1744
.set nb314_rinvMM, 1760
.set nb314_fstmp, 1776
.set nb314_is3, 1792
.set nb314_ii3, 1796
.set nb314_nri, 1800
.set nb314_iinr, 1808
.set nb314_jindex, 1816
.set nb314_jjnr, 1824
.set nb314_shift, 1832
.set nb314_shiftvec, 1840
.set nb314_facel, 1848
.set nb314_innerjjnr, 1856
.set nb314_innerk, 1864
.set nb314_n, 1868
.set nb314_nn1, 1872
.set nb314_nouter, 1876
.set nb314_ninner, 1880
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1896,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb314_nouter(%rsp)
        movl %eax,nb314_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb314_nri(%rsp)
        movq %rsi,nb314_iinr(%rsp)
        movq %rdx,nb314_jindex(%rsp)
        movq %rcx,nb314_jjnr(%rsp)
        movq %r8,nb314_shift(%rsp)
        movq %r9,nb314_shiftvec(%rsp)
        movq nb314_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb314_facel(%rsp)

        movq nb314_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb314_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb314_half(%rsp)
        movss nb314_half(%rsp),%xmm1
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
        movaps %xmm1,nb314_half(%rsp)
        movaps %xmm2,nb314_two(%rsp)
        movaps %xmm3,nb314_three(%rsp)
        movaps %xmm4,nb314_six(%rsp)
        movaps %xmm5,nb314_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb314_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb314_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb314_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb314_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb314_qqMM(%rsp)
        movaps %xmm4,nb314_qqMH(%rsp)
        movaps %xmm5,nb314_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb314_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb314_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx       ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb314_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb314_c6(%rsp)
        movaps %xmm1,nb314_c12(%rsp)

_nb_kernel314_x86_64_sse.nb314_threadloop: 
        movq  nb314_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel314_x86_64_sse.nb314_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314_x86_64_sse.nb314_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314_n(%rsp)
        movl %ebx,nb314_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314_x86_64_sse.nb314_outerstart
        jmp _nb_kernel314_x86_64_sse.nb314_end

_nb_kernel314_x86_64_sse.nb314_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314_nouter(%rsp),%ebx
        movl %ebx,nb314_nouter(%rsp)

_nb_kernel314_x86_64_sse.nb314_outer: 
        movq  nb314_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb314_is3(%rsp)      ## store is3 

        movq  nb314_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb314_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb314_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb314_ii3(%rsp)

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
        movaps %xmm3,nb314_ixO(%rsp)
        movaps %xmm4,nb314_iyO(%rsp)
        movaps %xmm5,nb314_izO(%rsp)
        movaps %xmm6,nb314_ixH1(%rsp)
        movaps %xmm7,nb314_iyH1(%rsp)

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
        movaps %xmm6,nb314_izH1(%rsp)
        movaps %xmm0,nb314_ixH2(%rsp)
        movaps %xmm1,nb314_iyH2(%rsp)
        movaps %xmm2,nb314_izH2(%rsp)
        movaps %xmm3,nb314_ixM(%rsp)
        movaps %xmm4,nb314_iyM(%rsp)
        movaps %xmm5,nb314_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb314_vctot(%rsp)
        movaps %xmm4,nb314_Vvdwtot(%rsp)
        movaps %xmm4,nb314_fixO(%rsp)
        movaps %xmm4,nb314_fiyO(%rsp)
        movaps %xmm4,nb314_fizO(%rsp)
        movaps %xmm4,nb314_fixH1(%rsp)
        movaps %xmm4,nb314_fiyH1(%rsp)
        movaps %xmm4,nb314_fizH1(%rsp)
        movaps %xmm4,nb314_fixH2(%rsp)
        movaps %xmm4,nb314_fiyH2(%rsp)
        movaps %xmm4,nb314_fizH2(%rsp)
        movaps %xmm4,nb314_fixM(%rsp)
        movaps %xmm4,nb314_fiyM(%rsp)
        movaps %xmm4,nb314_fizM(%rsp)

        movq  nb314_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb314_pos(%rbp),%rsi
        movq  nb314_faction(%rbp),%rdi
        movq  nb314_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb314_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb314_ninner(%rsp),%ecx
        movl  %ecx,nb314_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb314_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel314_x86_64_sse.nb314_unroll_loop
        jmp   _nb_kernel314_x86_64_sse.nb314_single_check
_nb_kernel314_x86_64_sse.nb314_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb314_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb314_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb314_pos(%rbp),%rsi       ## base of pos[] 

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

    subps nb314_ixO(%rsp),%xmm0
    subps nb314_iyO(%rsp),%xmm1
    subps nb314_izO(%rsp),%xmm2

    movaps %xmm0,%xmm4
    movaps %xmm1,%xmm5
    movaps %xmm2,%xmm6

    ## square it
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2

        addps  %xmm0,%xmm1
        addps  %xmm2,%xmm1
    ## rsq in xmm1

        ## move j O forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm10 ## jxOa jyOa  -   -
    movlps (%rdi,%rcx,4),%xmm11 ## jxOc jyOc  -   -
    movhps (%rdi,%rbx,4),%xmm10 ## jxOa jyOa jxOb jyOb 
    movhps (%rdi,%rdx,4),%xmm11 ## jxOc jyOc jxOd jyOd 

    movss  8(%rdi,%rax,4),%xmm12    ## jzOa  -  -  -
    movss  8(%rdi,%rcx,4),%xmm13    ## jzOc  -  -  -
    movhps 8(%rdi,%rbx,4),%xmm12    ## jzOa  -  jzOb  -
    movhps 8(%rdi,%rdx,4),%xmm13    ## jzOc  -  jzOd -

    shufps $136,%xmm13,%xmm12 ## 10001000 => jzOa jzOb jzOc jzOd

    ## xmm10: jxOa jyOa jxOb jyOb 
    ## xmm11: jxOc jyOc jxOd jyOd
    ## xmm12: jzOa jzOb jzOc jzOd

    ## calc rinvsq=1/rsq
        rcpps %xmm1,%xmm2
        movaps nb314_two(%rsp),%xmm0
        mulps %xmm2,%xmm1
        subps %xmm1,%xmm0
        mulps %xmm2,%xmm0       ## xmm0=rinvsq

        movaps %xmm0,%xmm1  ## rinvsq

        mulps  %xmm0,%xmm0  ## rinv4
        mulps  %xmm1,%xmm0      ## rinv6
        movaps %xmm0,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinv12

        mulps  nb314_c6(%rsp),%xmm0
        mulps  nb314_c12(%rsp),%xmm2
        movaps %xmm2,%xmm8
    subps  %xmm0,%xmm2  ## Vvdw=Vvdw12-Vvdw6 
        mulps  nb314_six(%rsp),%xmm0
        mulps  nb314_twelve(%rsp),%xmm8
        subps  %xmm0,%xmm8
        mulps  %xmm1,%xmm8      ## xmm8=total fscal 

    ## add potential to Vvdwtot
        addps  nb314_Vvdwtot(%rsp),%xmm2
    movaps %xmm2,nb314_Vvdwtot(%rsp)

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm8,%xmm4
        mulps  %xmm8,%xmm5
        mulps  %xmm8,%xmm6

    ## increment i force
    movaps nb314_fixO(%rsp),%xmm0
    movaps nb314_fiyO(%rsp),%xmm1
    movaps nb314_fizO(%rsp),%xmm2
    addps  %xmm4,%xmm0
    addps  %xmm5,%xmm1
    addps  %xmm6,%xmm2
    movaps %xmm0,nb314_fixO(%rsp)
    movaps %xmm1,nb314_fiyO(%rsp)
    movaps %xmm2,nb314_fizO(%rsp)

    ## update O forces
    ## xmm3 = fH1x , xmm4 = fH1y
    movaps %xmm4,%xmm3
    unpcklps %xmm5,%xmm4  ## fjx1 fjx1 fjy1 fjy2
    unpckhps %xmm5,%xmm3  ## fjx3 fjx4 fjy3 fjy4

    addps %xmm4,%xmm10
    addps %xmm3,%xmm11
    addps %xmm6,%xmm12

    movhlps  %xmm12,%xmm13 ## fH1zc fH1zd

    movlps %xmm10,(%rdi,%rax,4)
    movhps %xmm10,(%rdi,%rbx,4)
    movlps %xmm11,(%rdi,%rcx,4)
    movhps %xmm11,(%rdi,%rdx,4)
    movss  %xmm12,8(%rdi,%rax,4)
    movss  %xmm13,8(%rdi,%rcx,4)
    shufps $1,%xmm12,%xmm12
    shufps $1,%xmm13,%xmm13
    movss  %xmm12,8(%rdi,%rbx,4)
    movss  %xmm13,8(%rdi,%rdx,4)
    ## done with OO interaction.

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

    subps nb314_ixH1(%rsp),%xmm0
    subps nb314_iyH1(%rsp),%xmm1
    subps nb314_izH1(%rsp),%xmm2
    subps nb314_ixH2(%rsp),%xmm3
    subps nb314_iyH2(%rsp),%xmm4
    subps nb314_izH2(%rsp),%xmm5
    subps nb314_ixM(%rsp),%xmm6
    subps nb314_iyM(%rsp),%xmm7
    subps nb314_izM(%rsp),%xmm8

        movaps %xmm0,nb314_dxH1H1(%rsp)
        movaps %xmm1,nb314_dyH1H1(%rsp)
        movaps %xmm2,nb314_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH2H1(%rsp)
        movaps %xmm4,nb314_dyH2H1(%rsp)
        movaps %xmm5,nb314_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb314_dxMH1(%rsp)
        movaps %xmm7,nb314_dyMH1(%rsp)
        movaps %xmm8,nb314_dzMH1(%rsp)
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

        movaps  nb314_three(%rsp),%xmm9
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

        movaps  nb314_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H1 
        mulps   %xmm4,%xmm10 ## rinvH2H1
    mulps   %xmm4,%xmm11 ## rinvMH1

        movaps  %xmm9,nb314_rinvH1H1(%rsp)
        movaps  %xmm10,nb314_rinvH2H1(%rsp)
        movaps  %xmm11,nb314_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb314_tsc(%rsp),%xmm1
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

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb314_epsH1(%rsp)
    movaps    %xmm3,nb314_epsH2(%rsp)
    movaps    %xmm6,nb314_epsM(%rsp)

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

    movaps nb314_epsH1(%rsp),%xmm12
    movaps nb314_epsH2(%rsp),%xmm13
    movaps nb314_epsM(%rsp),%xmm14

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
    movaps nb314_qqHH(%rsp),%xmm12
    movaps nb314_qqMH(%rsp),%xmm13
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
    addps  nb314_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb314_vctot(%rsp)

    movaps nb314_tsc(%rsp),%xmm10
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

    mulps  nb314_rinvH1H1(%rsp),%xmm3
    mulps  nb314_rinvH2H1(%rsp),%xmm7
    mulps  nb314_rinvMH1(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb314_dxH1H1(%rsp),%xmm0
        mulps nb314_dyH1H1(%rsp),%xmm1
        mulps nb314_dzH1H1(%rsp),%xmm2
        mulps nb314_dxH2H1(%rsp),%xmm3
        mulps nb314_dyH2H1(%rsp),%xmm4
        mulps nb314_dzH2H1(%rsp),%xmm5
        mulps nb314_dxMH1(%rsp),%xmm6
        mulps nb314_dyMH1(%rsp),%xmm7
        mulps nb314_dzMH1(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb314_fixH1(%rsp),%xmm0
    addps nb314_fiyH1(%rsp),%xmm1
    addps nb314_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb314_fixH2(%rsp),%xmm3
    addps nb314_fiyH2(%rsp),%xmm4
    addps nb314_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb314_fixM(%rsp),%xmm6
    addps nb314_fiyM(%rsp),%xmm7
    addps nb314_fizM(%rsp),%xmm8

    movaps %xmm0,nb314_fixH1(%rsp)
    movaps %xmm1,nb314_fiyH1(%rsp)
    movaps %xmm2,nb314_fizH1(%rsp)
    movaps %xmm3,nb314_fixH2(%rsp)
    movaps %xmm4,nb314_fiyH2(%rsp)
    movaps %xmm5,nb314_fizH2(%rsp)
    movaps %xmm6,nb314_fixM(%rsp)
    movaps %xmm7,nb314_fiyM(%rsp)
    movaps %xmm8,nb314_fizM(%rsp)

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
        movq  nb314_pos(%rbp),%rsi
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

    subps nb314_ixH1(%rsp),%xmm0
    subps nb314_iyH1(%rsp),%xmm1
    subps nb314_izH1(%rsp),%xmm2
    subps nb314_ixH2(%rsp),%xmm3
    subps nb314_iyH2(%rsp),%xmm4
    subps nb314_izH2(%rsp),%xmm5
    subps nb314_ixM(%rsp),%xmm6
    subps nb314_iyM(%rsp),%xmm7
    subps nb314_izM(%rsp),%xmm8

        movaps %xmm0,nb314_dxH1H2(%rsp)
        movaps %xmm1,nb314_dyH1H2(%rsp)
        movaps %xmm2,nb314_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH2H2(%rsp)
        movaps %xmm4,nb314_dyH2H2(%rsp)
        movaps %xmm5,nb314_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb314_dxMH2(%rsp)
        movaps %xmm7,nb314_dyMH2(%rsp)
        movaps %xmm8,nb314_dzMH2(%rsp)
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

        movaps  nb314_three(%rsp),%xmm9
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

        movaps  nb314_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H2
        mulps   %xmm4,%xmm10 ## rinvH2H2
    mulps   %xmm4,%xmm11 ## rinvMH2

        movaps  %xmm9,nb314_rinvH1H2(%rsp)
        movaps  %xmm10,nb314_rinvH2H2(%rsp)
        movaps  %xmm11,nb314_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb314_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    movq nb314_VFtab(%rbp),%rsi

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

    movaps    %xmm0,nb314_epsH1(%rsp)
    movaps    %xmm3,nb314_epsH2(%rsp)
    movaps    %xmm6,nb314_epsM(%rsp)


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

    movaps nb314_epsH1(%rsp),%xmm12
    movaps nb314_epsH2(%rsp),%xmm13
    movaps nb314_epsM(%rsp),%xmm14

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
    movaps nb314_qqHH(%rsp),%xmm12
    movaps nb314_qqMH(%rsp),%xmm13
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
    addps  nb314_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb314_vctot(%rsp)

    movaps nb314_tsc(%rsp),%xmm10
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

    mulps  nb314_rinvH1H2(%rsp),%xmm3
    mulps  nb314_rinvH2H2(%rsp),%xmm7
    mulps  nb314_rinvMH2(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb314_dxH1H2(%rsp),%xmm0
        mulps nb314_dyH1H2(%rsp),%xmm1
        mulps nb314_dzH1H2(%rsp),%xmm2
        mulps nb314_dxH2H2(%rsp),%xmm3
        mulps nb314_dyH2H2(%rsp),%xmm4
        mulps nb314_dzH2H2(%rsp),%xmm5
        mulps nb314_dxMH2(%rsp),%xmm6
        mulps nb314_dyMH2(%rsp),%xmm7
        mulps nb314_dzMH2(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb314_fixH1(%rsp),%xmm0
    addps nb314_fiyH1(%rsp),%xmm1
    addps nb314_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb314_fixH2(%rsp),%xmm3
    addps nb314_fiyH2(%rsp),%xmm4
    addps nb314_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb314_fixM(%rsp),%xmm6
    addps nb314_fiyM(%rsp),%xmm7
    addps nb314_fizM(%rsp),%xmm8

    movaps %xmm0,nb314_fixH1(%rsp)
    movaps %xmm1,nb314_fiyH1(%rsp)
    movaps %xmm2,nb314_fizH1(%rsp)
    movaps %xmm3,nb314_fixH2(%rsp)
    movaps %xmm4,nb314_fiyH2(%rsp)
    movaps %xmm5,nb314_fizH2(%rsp)
    movaps %xmm6,nb314_fixM(%rsp)
    movaps %xmm7,nb314_fiyM(%rsp)
    movaps %xmm8,nb314_fizM(%rsp)

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

        movq  nb314_pos(%rbp),%rsi
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

    subps nb314_ixH1(%rsp),%xmm0
    subps nb314_iyH1(%rsp),%xmm1
    subps nb314_izH1(%rsp),%xmm2
    subps nb314_ixH2(%rsp),%xmm3
    subps nb314_iyH2(%rsp),%xmm4
    subps nb314_izH2(%rsp),%xmm5
    subps nb314_ixM(%rsp),%xmm6
    subps nb314_iyM(%rsp),%xmm7
    subps nb314_izM(%rsp),%xmm8

        movaps %xmm0,nb314_dxH1M(%rsp)
        movaps %xmm1,nb314_dyH1M(%rsp)
        movaps %xmm2,nb314_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb314_dxH2M(%rsp)
        movaps %xmm4,nb314_dyH2M(%rsp)
        movaps %xmm5,nb314_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb314_dxMM(%rsp)
        movaps %xmm7,nb314_dyMM(%rsp)
        movaps %xmm8,nb314_dzMM(%rsp)
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

        movaps  nb314_three(%rsp),%xmm9
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

        movaps  nb314_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1M
        mulps   %xmm4,%xmm10 ## rinvH2M
    mulps   %xmm4,%xmm11 ## rinvMM

        movaps  %xmm9,nb314_rinvH1M(%rsp)
        movaps  %xmm10,nb314_rinvH2M(%rsp)
        movaps  %xmm11,nb314_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb314_tsc(%rsp),%xmm1
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

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb314_epsH1(%rsp)
    movaps    %xmm3,nb314_epsH2(%rsp)
    movaps    %xmm6,nb314_epsM(%rsp)

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

    movaps nb314_epsH1(%rsp),%xmm12
    movaps nb314_epsH2(%rsp),%xmm13
    movaps nb314_epsM(%rsp),%xmm14

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
    movaps nb314_qqMH(%rsp),%xmm12
    movaps nb314_qqMM(%rsp),%xmm13
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
    addps  nb314_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb314_vctot(%rsp)

    movaps nb314_tsc(%rsp),%xmm10
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

    mulps  nb314_rinvH1M(%rsp),%xmm3
    mulps  nb314_rinvH2M(%rsp),%xmm7
    mulps  nb314_rinvMM(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb314_dxH1M(%rsp),%xmm0
        mulps nb314_dyH1M(%rsp),%xmm1
        mulps nb314_dzH1M(%rsp),%xmm2
        mulps nb314_dxH2M(%rsp),%xmm3
        mulps nb314_dyH2M(%rsp),%xmm4
        mulps nb314_dzH2M(%rsp),%xmm5
        mulps nb314_dxMM(%rsp),%xmm6
        mulps nb314_dyMM(%rsp),%xmm7
        mulps nb314_dzMM(%rsp),%xmm8

    movaps %xmm0,%xmm14
    movaps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb314_fixH1(%rsp),%xmm0
    addps nb314_fiyH1(%rsp),%xmm1
    addps nb314_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb314_fixH2(%rsp),%xmm3
    addps nb314_fiyH2(%rsp),%xmm4
    addps nb314_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb314_fixM(%rsp),%xmm6
    addps nb314_fiyM(%rsp),%xmm7
    addps nb314_fizM(%rsp),%xmm8

    movaps %xmm0,nb314_fixH1(%rsp)
    movaps %xmm1,nb314_fiyH1(%rsp)
    movaps %xmm2,nb314_fizH1(%rsp)
    movaps %xmm3,nb314_fixH2(%rsp)
    movaps %xmm4,nb314_fiyH2(%rsp)
    movaps %xmm5,nb314_fizH2(%rsp)
    movaps %xmm6,nb314_fixM(%rsp)
    movaps %xmm7,nb314_fiyM(%rsp)
    movaps %xmm8,nb314_fizM(%rsp)

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
        subl $4,nb314_innerk(%rsp)
        jl    _nb_kernel314_x86_64_sse.nb314_single_check
        jmp   _nb_kernel314_x86_64_sse.nb314_unroll_loop
_nb_kernel314_x86_64_sse.nb314_single_check: 
        addl $4,nb314_innerk(%rsp)
        jnz   _nb_kernel314_x86_64_sse.nb314_single_loop
        jmp   _nb_kernel314_x86_64_sse.nb314_updateouterdata
_nb_kernel314_x86_64_sse.nb314_single_loop: 
        movq  nb314_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb314_innerjjnr(%rsp)

        movq nb314_pos(%rbp),%rsi
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
        movaps %xmm6,nb314_jxO(%rsp)
        movaps %xmm3,nb314_jyO(%rsp)
        movaps %xmm1,nb314_jzO(%rsp)

    movaps %xmm6,%xmm0  ## jxO
    movaps %xmm1,%xmm2  ## jzO
    movaps %xmm3,%xmm1  ## jyO
    movaps %xmm3,%xmm4  ## jyO
    movaps %xmm6,%xmm3  ## jxO
    movaps %xmm2,%xmm5  ## jzO

        ## do O and H1 in parallel
        subps  nb314_ixO(%rsp),%xmm0
        subps  nb314_iyO(%rsp),%xmm1
        subps  nb314_izO(%rsp),%xmm2
        subps  nb314_ixH1(%rsp),%xmm3
        subps  nb314_iyH1(%rsp),%xmm4
        subps  nb314_izH1(%rsp),%xmm5

        movaps %xmm0,nb314_dxOO(%rsp)
        movaps %xmm1,nb314_dyOO(%rsp)
        movaps %xmm2,nb314_dzOO(%rsp)
        movaps %xmm3,nb314_dxH1H1(%rsp)
        movaps %xmm4,nb314_dyH1H1(%rsp)
        movaps %xmm5,nb314_dzH1H1(%rsp)

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
        ## Save H1 data in H1H1 
        movaps %xmm4,nb314_rsqH1H1(%rsp)

        ## do 1/x for O and 1/sqrt(x) for H1
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb314_two(%rsp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb314_three(%rsp),%xmm7
        mulss  %xmm1,%xmm2      ## 1/r2


        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## 1/r6
        mulps   nb314_half(%rsp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb314_rinvH1H1(%rsp)

        mulss  %xmm2,%xmm2      ## 1/r12
        mulss  nb314_c6(%rsp),%xmm1
        mulss  nb314_c12(%rsp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb314_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb314_Vvdwtot(%rsp)
        mulss  nb314_six(%rsp),%xmm1
        mulss  nb314_twelve(%rsp),%xmm2
        subss  %xmm1,%xmm2
        mulss  %xmm2,%xmm0      ## fscal
        movss  %xmm0,%xmm1
        movss  %xmm0,%xmm2
        mulss  nb314_dxOO(%rsp),%xmm0
        mulss  nb314_dyOO(%rsp),%xmm1
        mulss  nb314_dzOO(%rsp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addss   %xmm0,%xmm3
        addss   %xmm1,%xmm4
        addss   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%rsp)
        movaps  %xmm4,nb314_fjyO(%rsp)
        movaps  %xmm5,nb314_fjzO(%rsp)
        addss   nb314_fixO(%rsp),%xmm0
        addss   nb314_fiyO(%rsp),%xmm1
        addss   nb314_fizO(%rsp),%xmm2
        movss  %xmm0,nb314_fixO(%rsp)
        movss  %xmm1,nb314_fiyO(%rsp)
        movss  %xmm2,nb314_fizO(%rsp)

        ## do  H1 coulomb interaction
        movaps nb314_rinvH1H1(%rsp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb314_rsqH1H1(%rsp),%xmm1        ## r
        mulps nb314_tsc(%rsp),%xmm1

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

        movq nb314_VFtab(%rbp),%rsi

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
        addps  %xmm7,%xmm7      ## two*Heps2 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb314_qqHH(%rsp),%xmm3
        movhps  nb314_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 

        addps  nb314_vctot(%rsp),%xmm5
        movaps %xmm5,nb314_vctot(%rsp)

        mulps  nb314_tsc(%rsp),%xmm3
        xorps  %xmm2,%xmm2
        subps  %xmm3,%xmm2
        mulps  %xmm2,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb314_dxH1H1(%rsp),%xmm0
        mulps   nb314_dyH1H1(%rsp),%xmm1
        mulps   nb314_dzH1H1(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb314_fjxO(%rsp),%xmm3
        movaps  nb314_fjyO(%rsp),%xmm4
        movaps  nb314_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%rsp)
        movaps  %xmm4,nb314_fjyO(%rsp)
        movaps  %xmm5,nb314_fjzO(%rsp)
        addps   nb314_fixH1(%rsp),%xmm0
        addps   nb314_fiyH1(%rsp),%xmm1
        addps   nb314_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb314_fixH1(%rsp)
        movaps  %xmm1,nb314_fiyH1(%rsp)
        movaps  %xmm2,nb314_fizH1(%rsp)

        ## i H2 & M simultaneously first get i particle coords: 
    movaps  nb314_jxO(%rsp),%xmm0
    movaps  nb314_jyO(%rsp),%xmm1
    movaps  nb314_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps   nb314_ixH2(%rsp),%xmm0
        subps   nb314_iyH2(%rsp),%xmm1
        subps   nb314_izH2(%rsp),%xmm2
        subps   nb314_ixM(%rsp),%xmm3
        subps   nb314_iyM(%rsp),%xmm4
        subps   nb314_izM(%rsp),%xmm5
        movaps %xmm0,nb314_dxH2H2(%rsp)
        movaps %xmm1,nb314_dyH2H2(%rsp)
        movaps %xmm2,nb314_dzH2H2(%rsp)
        movaps %xmm3,nb314_dxMM(%rsp)
        movaps %xmm4,nb314_dyMM(%rsp)
        movaps %xmm5,nb314_dzMM(%rsp)
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
        movaps %xmm0,nb314_rsqH2H2(%rsp)
        movaps %xmm4,nb314_rsqMM(%rsp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314_half(%rsp),%xmm3   ## rinv H2 - j water 
        mulps   nb314_half(%rsp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb314_rinvH2H2(%rsp)
        movaps %xmm7,nb314_rinvMM(%rsp)

        movaps %xmm3,%xmm1
        mulps  nb314_rsqH2H2(%rsp),%xmm1        ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb314_tsc(%rsp),%xmm1

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
        mulps  nb314_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3

        ## fetch charges to xmm3 (temporary) 
        movss   nb314_qqHH(%rsp),%xmm3
        movhps  nb314_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb314_vctot(%rsp),%xmm5
        movaps %xmm5,nb314_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb314_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2
        mulps   nb314_dxH2H2(%rsp),%xmm0
        mulps   nb314_dyH2H2(%rsp),%xmm1
        mulps   nb314_dzH2H2(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb314_fjxO(%rsp),%xmm3
        movaps  nb314_fjyO(%rsp),%xmm4
        movaps  nb314_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb314_fjxO(%rsp)
        movaps  %xmm4,nb314_fjyO(%rsp)
        movaps  %xmm5,nb314_fjzO(%rsp)
        addps   nb314_fixH2(%rsp),%xmm0
        addps   nb314_fiyH2(%rsp),%xmm1
        addps   nb314_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb314_fixH2(%rsp)
        movaps  %xmm1,nb314_fiyH2(%rsp)
        movaps  %xmm2,nb314_fizH2(%rsp)

        ## do table for i M - j water interaction 
        movaps nb314_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314_rsqMM(%rsp),%xmm1          ## xmm0=rinv, xmm1=r 
        mulps  nb314_tsc(%rsp),%xmm1

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
        mulps  nb314_two(%rsp),%xmm7            ## two*Heps2 

        xorps  %xmm3,%xmm3
        ## fetch charges to xmm3 (temporary) 
        movss   nb314_qqMH(%rsp),%xmm3
        movhps  nb314_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC 
        addps  nb314_vctot(%rsp),%xmm5
        movaps %xmm5,nb314_vctot(%rsp)

        xorps  %xmm1,%xmm1

        mulps nb314_tsc(%rsp),%xmm3
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm1

        movaps  %xmm1,%xmm0
        movaps  %xmm1,%xmm2

        mulps   nb314_dxMM(%rsp),%xmm0
        mulps   nb314_dyMM(%rsp),%xmm1
        mulps   nb314_dzMM(%rsp),%xmm2
        movaps  nb314_fjxO(%rsp),%xmm3
        movaps  nb314_fjyO(%rsp),%xmm4
        movaps  nb314_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb314_faction(%rbp),%rsi
        movaps  %xmm3,nb314_fjxO(%rsp)
        movaps  %xmm4,nb314_fjyO(%rsp)
        movaps  %xmm5,nb314_fjzO(%rsp)
        addps   nb314_fixM(%rsp),%xmm0
        addps   nb314_fiyM(%rsp),%xmm1
        addps   nb314_fizM(%rsp),%xmm2
        movaps  %xmm0,nb314_fixM(%rsp)
        movaps  %xmm1,nb314_fiyM(%rsp)
        movaps  %xmm2,nb314_fizM(%rsp)

        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb314_fjxO(%rsp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb314_fjyO(%rsp),%xmm1   ## Oy H1y H2y My
        movaps  nb314_fjzO(%rsp),%xmm2   ## Oz H1z H2z Mz

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

        decl nb314_innerk(%rsp)
        jz    _nb_kernel314_x86_64_sse.nb314_updateouterdata
        jmp   _nb_kernel314_x86_64_sse.nb314_single_loop
_nb_kernel314_x86_64_sse.nb314_updateouterdata: 
        movl  nb314_ii3(%rsp),%ecx
        movq  nb314_faction(%rbp),%rdi
        movq  nb314_fshift(%rbp),%rsi
        movl  nb314_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb314_fixO(%rsp),%xmm0
        movaps nb314_fiyO(%rsp),%xmm1
        movaps nb314_fizO(%rsp),%xmm2

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
        movaps nb314_fixH1(%rsp),%xmm0
        movaps nb314_fiyH1(%rsp),%xmm1
        movaps nb314_fizH1(%rsp),%xmm2

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
        movaps nb314_fixH2(%rsp),%xmm0
        movaps nb314_fiyH2(%rsp),%xmm1
        movaps nb314_fizH2(%rsp),%xmm2

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
        movaps nb314_fixM(%rsp),%xmm0
        movaps nb314_fiyM(%rsp),%xmm1
        movaps nb314_fizM(%rsp),%xmm2

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
        movl nb314_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb314_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb314_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb314_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb314_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb314_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb314_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314_x86_64_sse.nb314_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314_n(%rsp)
        jmp _nb_kernel314_x86_64_sse.nb314_outer
_nb_kernel314_x86_64_sse.nb314_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314_x86_64_sse.nb314_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314_x86_64_sse.nb314_threadloop
_nb_kernel314_x86_64_sse.nb314_end: 
        movl nb314_nouter(%rsp),%eax
        movl nb314_ninner(%rsp),%ebx
        movq nb314_outeriter(%rbp),%rcx
        movq nb314_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1896,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel314nf_x86_64_sse
.globl _nb_kernel314nf_x86_64_sse
nb_kernel314nf_x86_64_sse:      
_nb_kernel314nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb314nf_fshift, 16
.set nb314nf_gid, 24
.set nb314nf_pos, 32
.set nb314nf_faction, 40
.set nb314nf_charge, 48
.set nb314nf_p_facel, 56
.set nb314nf_argkrf, 64
.set nb314nf_argcrf, 72
.set nb314nf_Vc, 80
.set nb314nf_type, 88
.set nb314nf_p_ntype, 96
.set nb314nf_vdwparam, 104
.set nb314nf_Vvdw, 112
.set nb314nf_p_tabscale, 120
.set nb314nf_VFtab, 128
.set nb314nf_invsqrta, 136
.set nb314nf_dvda, 144
.set nb314nf_p_gbtabscale, 152
.set nb314nf_GBtab, 160
.set nb314nf_p_nthreads, 168
.set nb314nf_count, 176
.set nb314nf_mtx, 184
.set nb314nf_outeriter, 192
.set nb314nf_inneriter, 200
.set nb314nf_work, 208
        ## bottom of stack is cache-aligned for sse use 
.set nb314nf_ixO, 0
.set nb314nf_iyO, 16
.set nb314nf_izO, 32
.set nb314nf_ixH1, 48
.set nb314nf_iyH1, 64
.set nb314nf_izH1, 80
.set nb314nf_ixH2, 96
.set nb314nf_iyH2, 112
.set nb314nf_izH2, 128
.set nb314nf_ixM, 144
.set nb314nf_iyM, 160
.set nb314nf_izM, 176
.set nb314nf_jxO, 192
.set nb314nf_jyO, 208
.set nb314nf_jzO, 224
.set nb314nf_jxH1, 240
.set nb314nf_jyH1, 256
.set nb314nf_jzH1, 272
.set nb314nf_jxH2, 288
.set nb314nf_jyH2, 304
.set nb314nf_jzH2, 320
.set nb314nf_jxM, 336
.set nb314nf_jyM, 352
.set nb314nf_jzM, 368
.set nb314nf_qqMM, 384
.set nb314nf_qqMH, 400
.set nb314nf_qqHH, 416
.set nb314nf_two, 432
.set nb314nf_tsc, 448
.set nb314nf_c6, 464
.set nb314nf_c12, 480
.set nb314nf_vctot, 496
.set nb314nf_Vvdwtot, 512
.set nb314nf_half, 528
.set nb314nf_three, 544
.set nb314nf_rsqOO, 560
.set nb314nf_rsqH1H1, 576
.set nb314nf_rsqH1H2, 592
.set nb314nf_rsqH1M, 608
.set nb314nf_rsqH2H1, 624
.set nb314nf_rsqH2H2, 640
.set nb314nf_rsqH2M, 656
.set nb314nf_rsqMH1, 672
.set nb314nf_rsqMH2, 688
.set nb314nf_rsqMM, 704
.set nb314nf_rinvsqOO, 720
.set nb314nf_rinvH1H1, 736
.set nb314nf_rinvH1H2, 752
.set nb314nf_rinvH1M, 768
.set nb314nf_rinvH2H1, 784
.set nb314nf_rinvH2H2, 800
.set nb314nf_rinvH2M, 816
.set nb314nf_rinvMH1, 832
.set nb314nf_rinvMH2, 848
.set nb314nf_rinvMM, 864
.set nb314nf_is3, 880
.set nb314nf_ii3, 884
.set nb314nf_nri, 888
.set nb314nf_iinr, 896
.set nb314nf_jindex, 904
.set nb314nf_jjnr, 912
.set nb314nf_shift, 920
.set nb314nf_shiftvec, 928
.set nb314nf_facel, 936
.set nb314nf_innerjjnr, 944
.set nb314nf_innerk, 952
.set nb314nf_n, 956
.set nb314nf_nn1, 960
.set nb314nf_nouter, 964
.set nb314nf_ninner, 968
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $984,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb314nf_nouter(%rsp)
        movl %eax,nb314nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb314nf_nri(%rsp)
        movq %rsi,nb314nf_iinr(%rsp)
        movq %rdx,nb314nf_jindex(%rsp)
        movq %rcx,nb314nf_jjnr(%rsp)
        movq %r8,nb314nf_shift(%rsp)
        movq %r9,nb314nf_shiftvec(%rsp)
        movq nb314nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb314nf_facel(%rsp)

        movq nb314nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb314nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb314nf_half(%rsp)
        movss nb314nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb314nf_half(%rsp)
        movaps %xmm2,nb314nf_two(%rsp)
        movaps %xmm3,nb314nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb314nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb314nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb314nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb314nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb314nf_qqMM(%rsp)
        movaps %xmm4,nb314nf_qqMH(%rsp)
        movaps %xmm5,nb314nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb314nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb314nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx       ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb314nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb314nf_c6(%rsp)
        movaps %xmm1,nb314nf_c12(%rsp)

_nb_kernel314nf_x86_64_sse.nb314nf_threadloop: 
        movq  nb314nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel314nf_x86_64_sse.nb314nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314nf_x86_64_sse.nb314nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314nf_n(%rsp)
        movl %ebx,nb314nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314nf_x86_64_sse.nb314nf_outerstart
        jmp _nb_kernel314nf_x86_64_sse.nb314nf_end

_nb_kernel314nf_x86_64_sse.nb314nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314nf_nouter(%rsp),%ebx
        movl %ebx,nb314nf_nouter(%rsp)

_nb_kernel314nf_x86_64_sse.nb314nf_outer: 
        movq  nb314nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb314nf_is3(%rsp)            ## store is3 

        movq  nb314nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb314nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb314nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb314nf_ii3(%rsp)

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
        movaps %xmm3,nb314nf_ixO(%rsp)
        movaps %xmm4,nb314nf_iyO(%rsp)
        movaps %xmm5,nb314nf_izO(%rsp)
        movaps %xmm6,nb314nf_ixH1(%rsp)
        movaps %xmm7,nb314nf_iyH1(%rsp)

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
        movaps %xmm6,nb314nf_izH1(%rsp)
        movaps %xmm0,nb314nf_ixH2(%rsp)
        movaps %xmm1,nb314nf_iyH2(%rsp)
        movaps %xmm2,nb314nf_izH2(%rsp)
        movaps %xmm3,nb314nf_ixM(%rsp)
        movaps %xmm4,nb314nf_iyM(%rsp)
        movaps %xmm5,nb314nf_izM(%rsp)

        ## clear vctot  
        xorps %xmm4,%xmm4
        movaps %xmm4,nb314nf_vctot(%rsp)
        movaps %xmm4,nb314nf_Vvdwtot(%rsp)

        movq  nb314nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb314nf_pos(%rbp),%rsi
        movq  nb314nf_faction(%rbp),%rdi
        movq  nb314nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb314nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb314nf_ninner(%rsp),%ecx
        movl  %ecx,nb314nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb314nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel314nf_x86_64_sse.nb314nf_unroll_loop
        jmp   _nb_kernel314nf_x86_64_sse.nb314nf_single_check
_nb_kernel314nf_x86_64_sse.nb314nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb314nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb314nf_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb314nf_pos(%rbp),%rsi     ## base of pos[] 

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
        movaps %xmm0,nb314nf_jxO(%rsp)
        movaps %xmm1,nb314nf_jyO(%rsp)
        movaps %xmm2,nb314nf_jzO(%rsp)
        movaps %xmm3,nb314nf_jxH1(%rsp)

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
        movaps %xmm0,nb314nf_jyH1(%rsp)
        movaps %xmm1,nb314nf_jzH1(%rsp)
        movaps %xmm2,nb314nf_jxH2(%rsp)
        movaps %xmm3,nb314nf_jyH2(%rsp)

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
        movaps %xmm0,nb314nf_jzH2(%rsp)
        movaps %xmm1,nb314nf_jxM(%rsp)
        movaps %xmm2,nb314nf_jyM(%rsp)
        movaps %xmm3,nb314nf_jzM(%rsp)

        ## start calculating pairwise distances
        movaps nb314nf_ixO(%rsp),%xmm0
        movaps nb314nf_iyO(%rsp),%xmm1
        movaps nb314nf_izO(%rsp),%xmm2
        movaps nb314nf_ixH1(%rsp),%xmm3
        movaps nb314nf_iyH1(%rsp),%xmm4
        movaps nb314nf_izH1(%rsp),%xmm5
        subps  nb314nf_jxO(%rsp),%xmm0
        subps  nb314nf_jyO(%rsp),%xmm1
        subps  nb314nf_jzO(%rsp),%xmm2
        subps  nb314nf_jxH1(%rsp),%xmm3
        subps  nb314nf_jyH1(%rsp),%xmm4
        subps  nb314nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqOO(%rsp)
        movaps %xmm3,nb314nf_rsqH1H1(%rsp)

        movaps nb314nf_ixH1(%rsp),%xmm0
        movaps nb314nf_iyH1(%rsp),%xmm1
        movaps nb314nf_izH1(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH2(%rsp),%xmm0
        subps  nb314nf_jyH2(%rsp),%xmm1
        subps  nb314nf_jzH2(%rsp),%xmm2
        subps  nb314nf_jxM(%rsp),%xmm3
        subps  nb314nf_jyM(%rsp),%xmm4
        subps  nb314nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqH1H2(%rsp)
        movaps %xmm3,nb314nf_rsqH1M(%rsp)

        movaps nb314nf_ixH2(%rsp),%xmm0
        movaps nb314nf_iyH2(%rsp),%xmm1
        movaps nb314nf_izH2(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH1(%rsp),%xmm0
        subps  nb314nf_jyH1(%rsp),%xmm1
        subps  nb314nf_jzH1(%rsp),%xmm2
        subps  nb314nf_jxH2(%rsp),%xmm3
        subps  nb314nf_jyH2(%rsp),%xmm4
        subps  nb314nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqH2H1(%rsp)
        movaps %xmm3,nb314nf_rsqH2H2(%rsp)

        movaps nb314nf_ixH2(%rsp),%xmm0
        movaps nb314nf_iyH2(%rsp),%xmm1
        movaps nb314nf_izH2(%rsp),%xmm2
        movaps nb314nf_ixM(%rsp),%xmm3
        movaps nb314nf_iyM(%rsp),%xmm4
        movaps nb314nf_izM(%rsp),%xmm5
        subps  nb314nf_jxM(%rsp),%xmm0
        subps  nb314nf_jyM(%rsp),%xmm1
        subps  nb314nf_jzM(%rsp),%xmm2
        subps  nb314nf_jxH1(%rsp),%xmm3
        subps  nb314nf_jyH1(%rsp),%xmm4
        subps  nb314nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqH2M(%rsp)
        movaps %xmm4,nb314nf_rsqMH1(%rsp)

        movaps nb314nf_ixM(%rsp),%xmm0
        movaps nb314nf_iyM(%rsp),%xmm1
        movaps nb314nf_izM(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb314nf_jxH2(%rsp),%xmm0
        subps  nb314nf_jyH2(%rsp),%xmm1
        subps  nb314nf_jzH2(%rsp),%xmm2
        subps  nb314nf_jxM(%rsp),%xmm3
        subps  nb314nf_jyM(%rsp),%xmm4
        subps  nb314nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqMH2(%rsp)
        movaps %xmm4,nb314nf_rsqMM(%rsp)

        ## start by doing reciprocal for OO
        movaps  nb314nf_rsqOO(%rsp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb314nf_two(%rsp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb314nf_rinvsqOO(%rsp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb314nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb314nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH1H1(%rsp),%xmm1
        mulps   nb314nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%rsp),%xmm3   ## rinvH1H1 
        mulps   nb314nf_half(%rsp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb314nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb314nf_rinvH1H2(%rsp)

        rsqrtps nb314nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb314nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH1M(%rsp),%xmm1
        mulps   nb314nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%rsp),%xmm3
        mulps   nb314nf_half(%rsp),%xmm7
        movaps  %xmm3,nb314nf_rinvH1M(%rsp)
        movaps  %xmm7,nb314nf_rinvH2H1(%rsp)

        rsqrtps nb314nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb314nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqH2H2(%rsp),%xmm1
        mulps   nb314nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%rsp),%xmm3
        mulps   nb314nf_half(%rsp),%xmm7
        movaps  %xmm3,nb314nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb314nf_rinvH2M(%rsp)

        rsqrtps nb314nf_rsqMH1(%rsp),%xmm1
        rsqrtps nb314nf_rsqMH2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb314nf_rsqMH1(%rsp),%xmm1
        mulps   nb314nf_rsqMH2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%rsp),%xmm3
        mulps   nb314nf_half(%rsp),%xmm7
        movaps  %xmm3,nb314nf_rinvMH1(%rsp)
        movaps  %xmm7,nb314nf_rinvMH2(%rsp)

        rsqrtps nb314nf_rsqMM(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb314nf_three(%rsp),%xmm3
        mulps   nb314nf_rsqMM(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb314nf_half(%rsp),%xmm3
        movaps  %xmm3,nb314nf_rinvMM(%rsp)

        ## start with OO LJ interaction
        movaps nb314nf_rinvsqOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb314nf_c6(%rsp),%xmm1
        mulps  nb314nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb314nf_Vvdwtot(%rsp),%xmm4
        movaps %xmm4,nb314nf_Vvdwtot(%rsp)

        ## Coulomb interactions - first H1H1
        movaps nb314nf_rinvH1H1(%rsp),%xmm0

        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1

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

        movq nb314nf_VFtab(%rbp),%rsi
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
        movaps nb314nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## update vctot 
        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## H1-H2 interaction 
        movaps nb314nf_rinvH1H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## H1-M interaction  
        movaps nb314nf_rinvH1M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## H2-H1 interaction 
        movaps nb314nf_rinvH2H1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## H2-H2 interaction 
        movaps nb314nf_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqHH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## H2-M interaction 
        movaps nb314nf_rinvH2M(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## M-H1 interaction 
        movaps nb314nf_rinvMH1(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## M-H2 interaction 
        movaps nb314nf_rinvMH2(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqMH(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## M-M interaction 
        movaps nb314nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1
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
        movaps nb314nf_qqMM(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb314nf_innerk(%rsp)
        jl    _nb_kernel314nf_x86_64_sse.nb314nf_single_check
        jmp   _nb_kernel314nf_x86_64_sse.nb314nf_unroll_loop
_nb_kernel314nf_x86_64_sse.nb314nf_single_check: 
        addl $4,nb314nf_innerk(%rsp)
        jnz   _nb_kernel314nf_x86_64_sse.nb314nf_single_loop
        jmp   _nb_kernel314nf_x86_64_sse.nb314nf_updateouterdata
_nb_kernel314nf_x86_64_sse.nb314nf_single_loop: 
        movq  nb314nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb314nf_innerjjnr(%rsp)

        movq nb314nf_pos(%rbp),%rsi
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
        movaps %xmm6,nb314nf_jxO(%rsp)
        movaps %xmm3,nb314nf_jyO(%rsp)
        movaps %xmm1,nb314nf_jzO(%rsp)

        ## do O and H1 in parallel
        movaps nb314nf_ixO(%rsp),%xmm0
        movaps nb314nf_iyO(%rsp),%xmm1
        movaps nb314nf_izO(%rsp),%xmm2
        movaps nb314nf_ixH1(%rsp),%xmm3
        movaps nb314nf_iyH1(%rsp),%xmm4
        movaps nb314nf_izH1(%rsp),%xmm5
        subps  nb314nf_jxO(%rsp),%xmm0
        subps  nb314nf_jyO(%rsp),%xmm1
        subps  nb314nf_jzO(%rsp),%xmm2
        subps  nb314nf_jxO(%rsp),%xmm3
        subps  nb314nf_jyO(%rsp),%xmm4
        subps  nb314nf_jzO(%rsp),%xmm5

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
        ## Save H1 data in H1H1 
        movaps %xmm4,nb314nf_rsqH1H1(%rsp)

        ## do 1/x for O and 1/sqrt(x) for H1
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb314nf_two(%rsp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb314nf_three(%rsp),%xmm7
        mulss  %xmm1,%xmm2      ## 1/r2


        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## 1/r6
        mulps   nb314nf_half(%rsp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb314nf_rinvH1H1(%rsp)

        mulss  %xmm2,%xmm2      ## 1/r12
        mulss  nb314nf_c6(%rsp),%xmm1
        mulss  nb314nf_c12(%rsp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb314nf_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb314nf_Vvdwtot(%rsp)

        ## do  H1 coulomb interaction
        movaps nb314nf_rinvH1H1(%rsp),%xmm0   ## rinv 
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqH1H1(%rsp),%xmm1      ## r
        mulps nb314nf_tsc(%rsp),%xmm1

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

        movq nb314nf_VFtab(%rbp),%rsi

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

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb314nf_qqHH(%rsp),%xmm3
        movhps  nb314nf_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 

        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## i H2 & M simultaneously first get i particle coords: 
        movaps  nb314nf_ixH2(%rsp),%xmm0
        movaps  nb314nf_iyH2(%rsp),%xmm1
        movaps  nb314nf_izH2(%rsp),%xmm2
        movaps  nb314nf_ixM(%rsp),%xmm3
        movaps  nb314nf_iyM(%rsp),%xmm4
        movaps  nb314nf_izM(%rsp),%xmm5
        subps   nb314nf_jxO(%rsp),%xmm0
        subps   nb314nf_jyO(%rsp),%xmm1
        subps   nb314nf_jzO(%rsp),%xmm2
        subps   nb314nf_jxO(%rsp),%xmm3
        subps   nb314nf_jyO(%rsp),%xmm4
        subps   nb314nf_jzO(%rsp),%xmm5
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
        movaps %xmm0,nb314nf_rsqH2H2(%rsp)
        movaps %xmm4,nb314nf_rsqMM(%rsp)
        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb314nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb314nf_half(%rsp),%xmm3   ## rinv H2 - j water 
        mulps   nb314nf_half(%rsp),%xmm7   ## rinv M - j water  

        movaps %xmm3,nb314nf_rinvH2H2(%rsp)
        movaps %xmm7,nb314nf_rinvMM(%rsp)

        movaps %xmm3,%xmm1
        mulps  nb314nf_rsqH2H2(%rsp),%xmm1      ## xmm1=r 
        movaps %xmm3,%xmm0      ## xmm0=rinv 
        mulps  nb314nf_tsc(%rsp),%xmm1

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
        movss   nb314nf_qqHH(%rsp),%xmm3
        movhps  nb314nf_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        ## do table for i M - j water interaction 
        movaps nb314nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb314nf_rsqMM(%rsp),%xmm1        ## xmm0=rinv, xmm1=r 
        mulps  nb314nf_tsc(%rsp),%xmm1

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
        movss   nb314nf_qqMH(%rsp),%xmm3
        movhps  nb314nf_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001

        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point xmm5 contains vcoul 
        addps  nb314nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb314nf_vctot(%rsp)

        decl nb314nf_innerk(%rsp)
        jz    _nb_kernel314nf_x86_64_sse.nb314nf_updateouterdata
        jmp   _nb_kernel314nf_x86_64_sse.nb314nf_single_loop
_nb_kernel314nf_x86_64_sse.nb314nf_updateouterdata: 
        ## get n from stack
        movl nb314nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb314nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb314nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb314nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb314nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb314nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb314nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314nf_x86_64_sse.nb314nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314nf_n(%rsp)
        jmp _nb_kernel314nf_x86_64_sse.nb314nf_outer
_nb_kernel314nf_x86_64_sse.nb314nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314nf_x86_64_sse.nb314nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314nf_x86_64_sse.nb314nf_threadloop
_nb_kernel314nf_x86_64_sse.nb314nf_end: 
        movl nb314nf_nouter(%rsp),%eax
        movl nb314nf_ninner(%rsp),%ebx
        movq nb314nf_outeriter(%rbp),%rcx
        movq nb314nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $984,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



