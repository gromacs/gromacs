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






.globl nb_kernel114_x86_64_sse
.globl _nb_kernel114_x86_64_sse
nb_kernel114_x86_64_sse:        
_nb_kernel114_x86_64_sse:       
.set nb114_fshift, 16
.set nb114_gid, 24
.set nb114_pos, 32
.set nb114_faction, 40
.set nb114_charge, 48
.set nb114_p_facel, 56
.set nb114_argkrf, 64
.set nb114_argcrf, 72
.set nb114_Vc, 80
.set nb114_type, 88
.set nb114_p_ntype, 96
.set nb114_vdwparam, 104
.set nb114_Vvdw, 112
.set nb114_p_tabscale, 120
.set nb114_VFtab, 128
.set nb114_invsqrta, 136
.set nb114_dvda, 144
.set nb114_p_gbtabscale, 152
.set nb114_GBtab, 160
.set nb114_p_nthreads, 168
.set nb114_count, 176
.set nb114_mtx, 184
.set nb114_outeriter, 192
.set nb114_inneriter, 200
.set nb114_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb114_ixO, 0
.set nb114_iyO, 16
.set nb114_izO, 32
.set nb114_ixH1, 48
.set nb114_iyH1, 64
.set nb114_izH1, 80
.set nb114_ixH2, 96
.set nb114_iyH2, 112
.set nb114_izH2, 128
.set nb114_ixM, 144
.set nb114_iyM, 160
.set nb114_izM, 176
.set nb114_jxO, 192
.set nb114_jyO, 208
.set nb114_jzO, 224
.set nb114_jxH1, 240
.set nb114_jyH1, 256
.set nb114_jzH1, 272
.set nb114_jxH2, 288
.set nb114_jyH2, 304
.set nb114_jzH2, 320
.set nb114_jxM, 336
.set nb114_jyM, 352
.set nb114_jzM, 368
.set nb114_dxOO, 384
.set nb114_dyOO, 400
.set nb114_dzOO, 416
.set nb114_dxH1H1, 432
.set nb114_dyH1H1, 448
.set nb114_dzH1H1, 464
.set nb114_dxH1H2, 480
.set nb114_dyH1H2, 496
.set nb114_dzH1H2, 512
.set nb114_dxH1M, 528
.set nb114_dyH1M, 544
.set nb114_dzH1M, 560
.set nb114_dxH2H1, 576
.set nb114_dyH2H1, 592
.set nb114_dzH2H1, 608
.set nb114_dxH2H2, 624
.set nb114_dyH2H2, 640
.set nb114_dzH2H2, 656
.set nb114_dxH2M, 672
.set nb114_dyH2M, 688
.set nb114_dzH2M, 704
.set nb114_dxMH1, 720
.set nb114_dyMH1, 736
.set nb114_dzMH1, 752
.set nb114_dxMH2, 768
.set nb114_dyMH2, 784
.set nb114_dzMH2, 800
.set nb114_dxMM, 816
.set nb114_dyMM, 832
.set nb114_dzMM, 848
.set nb114_qqMM, 864
.set nb114_qqMH, 880
.set nb114_qqHH, 896
.set nb114_two, 912
.set nb114_c6, 928
.set nb114_c12, 944
.set nb114_six, 960
.set nb114_twelve, 976
.set nb114_vctot, 992
.set nb114_Vvdwtot, 1008
.set nb114_fixO, 1024
.set nb114_fiyO, 1040
.set nb114_fizO, 1056
.set nb114_fixH1, 1072
.set nb114_fiyH1, 1088
.set nb114_fizH1, 1104
.set nb114_fixH2, 1120
.set nb114_fiyH2, 1136
.set nb114_fizH2, 1152
.set nb114_fixM, 1168
.set nb114_fiyM, 1184
.set nb114_fizM, 1200
.set nb114_fjxO, 1216
.set nb114_fjyO, 1232
.set nb114_fjzO, 1248
.set nb114_fjxH1, 1264
.set nb114_fjyH1, 1280
.set nb114_fjzH1, 1296
.set nb114_fjxH2, 1312
.set nb114_fjyH2, 1328
.set nb114_fjzH2, 1344
.set nb114_fjxM, 1360
.set nb114_fjyM, 1376
.set nb114_fjzM, 1392
.set nb114_half, 1408
.set nb114_three, 1424
.set nb114_rsqOO, 1440
.set nb114_rsqH1H1, 1456
.set nb114_rsqH1H2, 1472
.set nb114_rsqH1M, 1488
.set nb114_rsqH2H1, 1504
.set nb114_rsqH2H2, 1520
.set nb114_rsqH2M, 1536
.set nb114_rsqMH1, 1552
.set nb114_rsqMH2, 1568
.set nb114_rsqMM, 1584
.set nb114_rinvsqOO, 1600
.set nb114_rinvH1H1, 1616
.set nb114_rinvH1H2, 1632
.set nb114_rinvH1M, 1648
.set nb114_rinvH2H1, 1664
.set nb114_rinvH2H2, 1680
.set nb114_rinvH2M, 1696
.set nb114_rinvMH1, 1712
.set nb114_rinvMH2, 1728
.set nb114_rinvMM, 1744
.set nb114_fstmp, 1760
.set nb114_nri, 1776
.set nb114_iinr, 1784
.set nb114_jindex, 1792
.set nb114_jjnr, 1800
.set nb114_shift, 1808
.set nb114_shiftvec, 1816
.set nb114_facel, 1824
.set nb114_innerjjnr, 1832
.set nb114_is3, 1840
.set nb114_ii3, 1844
.set nb114_innerk, 1848
.set nb114_n, 1856
.set nb114_nn1, 1860
.set nb114_nouter, 1864
.set nb114_ninner, 1868

        push %rbp
        movq %rsp,%rbp
        push %rbx


        push %r12
        push %r13
        push %r14
        push %r15

        subq $1880,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb114_nouter(%rsp)
        movl %eax,nb114_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb114_nri(%rsp)
        movq %rsi,nb114_iinr(%rsp)
        movq %rdx,nb114_jindex(%rsp)
        movq %rcx,nb114_jjnr(%rsp)
        movq %r8,nb114_shift(%rsp)
        movq %r9,nb114_shiftvec(%rsp)
        movq nb114_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb114_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb114_half(%rsp)
        movss nb114_half(%rsp),%xmm1
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
        movaps %xmm1,nb114_half(%rsp)
        movaps %xmm2,nb114_two(%rsp)
        movaps %xmm3,nb114_three(%rsp)
        movaps %xmm4,nb114_six(%rsp)
        movaps %xmm5,nb114_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb114_iinr(%rsp),%rcx     ## rcx = pointer into iinr[]
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb114_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb114_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb114_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb114_qqMM(%rsp)
        movaps %xmm4,nb114_qqMH(%rsp)
        movaps %xmm5,nb114_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb114_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb114_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx ## rcx = ntia = 2*ntype*type[ii0] 
        addq  %rcx,%rdx
        movq  nb114_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb114_c6(%rsp)
        movaps %xmm1,nb114_c12(%rsp)

_nb_kernel114_x86_64_sse.nb114_threadloop: 
        movq  nb114_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel114_x86_64_sse.nb114_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114_x86_64_sse.nb114_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114_n(%rsp)
        movl %ebx,nb114_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114_x86_64_sse.nb114_outerstart
        jmp _nb_kernel114_x86_64_sse.nb114_end

_nb_kernel114_x86_64_sse.nb114_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114_nouter(%rsp),%ebx
        movl %ebx,nb114_nouter(%rsp)

_nb_kernel114_x86_64_sse.nb114_outer: 
        movq  nb114_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb114_is3(%rsp)      ## store is3 

        movq  nb114_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb114_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb114_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb114_ii3(%rsp)

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
        movaps %xmm3,nb114_ixO(%rsp)
        movaps %xmm4,nb114_iyO(%rsp)
        movaps %xmm5,nb114_izO(%rsp)
        movaps %xmm6,nb114_ixH1(%rsp)
        movaps %xmm7,nb114_iyH1(%rsp)

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
        movaps %xmm6,nb114_izH1(%rsp)
        movaps %xmm0,nb114_ixH2(%rsp)
        movaps %xmm1,nb114_iyH2(%rsp)
        movaps %xmm2,nb114_izH2(%rsp)
        movaps %xmm3,nb114_ixM(%rsp)
        movaps %xmm4,nb114_iyM(%rsp)
        movaps %xmm5,nb114_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb114_vctot(%rsp)
        movaps %xmm4,nb114_Vvdwtot(%rsp)
        movaps %xmm4,nb114_fixO(%rsp)
        movaps %xmm4,nb114_fiyO(%rsp)
        movaps %xmm4,nb114_fizO(%rsp)
        movaps %xmm4,nb114_fixH1(%rsp)
        movaps %xmm4,nb114_fiyH1(%rsp)
        movaps %xmm4,nb114_fizH1(%rsp)
        movaps %xmm4,nb114_fixH2(%rsp)
        movaps %xmm4,nb114_fiyH2(%rsp)
        movaps %xmm4,nb114_fizH2(%rsp)
        movaps %xmm4,nb114_fixM(%rsp)
        movaps %xmm4,nb114_fiyM(%rsp)
        movaps %xmm4,nb114_fizM(%rsp)

        movq  nb114_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb114_pos(%rbp),%rsi
        movq  nb114_faction(%rbp),%rdi
        movq  nb114_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb114_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb114_ninner(%rsp),%ecx
        movl  %ecx,nb114_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb114_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel114_x86_64_sse.nb114_unroll_loop
        jmp   _nb_kernel114_x86_64_sse.nb114_single_check
_nb_kernel114_x86_64_sse.nb114_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb114_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb114_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb114_pos(%rbp),%rsi       ## base of pos[] 

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

    subps nb114_ixO(%rsp),%xmm0
    subps nb114_iyO(%rsp),%xmm1
    subps nb114_izO(%rsp),%xmm2

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
        movaps nb114_two(%rsp),%xmm0
        mulps %xmm2,%xmm1
        subps %xmm1,%xmm0
        mulps %xmm2,%xmm0       ## xmm0=rinvsq

        movaps %xmm0,%xmm1  ## rinvsq

        mulps  %xmm0,%xmm0  ## rinv4
        mulps  %xmm1,%xmm0      ## rinv6
        movaps %xmm0,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinv12

        mulps  nb114_c6(%rsp),%xmm0
        mulps  nb114_c12(%rsp),%xmm2
        movaps %xmm2,%xmm8
    subps  %xmm0,%xmm2  ## Vvdw=Vvdw12-Vvdw6 
        mulps  nb114_six(%rsp),%xmm0
        mulps  nb114_twelve(%rsp),%xmm8
        subps  %xmm0,%xmm8
        mulps  %xmm1,%xmm8      ## xmm8=total fscal 

    ## add potential to Vvdwtot
        addps  nb114_Vvdwtot(%rsp),%xmm2
    movaps %xmm2,nb114_Vvdwtot(%rsp)

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm8,%xmm4
        mulps  %xmm8,%xmm5
        mulps  %xmm8,%xmm6

    ## increment i force
    movaps nb114_fixO(%rsp),%xmm0
    movaps nb114_fiyO(%rsp),%xmm1
    movaps nb114_fizO(%rsp),%xmm2
    addps  %xmm4,%xmm0
    addps  %xmm5,%xmm1
    addps  %xmm6,%xmm2
    movaps %xmm0,nb114_fixO(%rsp)
    movaps %xmm1,nb114_fiyO(%rsp)
    movaps %xmm2,nb114_fizO(%rsp)

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

    movaps %xmm0,%xmm4
    unpcklps %xmm1,%xmm0 ## jxH1a jxH1c jyH1a jyH1c        
    unpckhps %xmm1,%xmm4 ## jxH1b jxH1d jyH1b jyH1d
    movaps %xmm0,%xmm1
    unpcklps %xmm4,%xmm0 ## x
    unpckhps %xmm4,%xmm1 ## y

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzH1a jzH1b jzH1c jzH1d

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb114_ixH1(%rsp),%xmm0
    subps nb114_iyH1(%rsp),%xmm1
    subps nb114_izH1(%rsp),%xmm2
    subps nb114_ixH2(%rsp),%xmm3
    subps nb114_iyH2(%rsp),%xmm4
    subps nb114_izH2(%rsp),%xmm5
    subps nb114_ixM(%rsp),%xmm6
    subps nb114_iyM(%rsp),%xmm7
    subps nb114_izM(%rsp),%xmm8

        movaps %xmm0,nb114_dxH1H1(%rsp)
        movaps %xmm1,nb114_dyH1H1(%rsp)
        movaps %xmm2,nb114_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH2H1(%rsp)
        movaps %xmm4,nb114_dyH2H1(%rsp)
        movaps %xmm5,nb114_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb114_dxMH1(%rsp)
        movaps %xmm7,nb114_dyMH1(%rsp)
        movaps %xmm8,nb114_dzMH1(%rsp)
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

        movaps  nb114_three(%rsp),%xmm9
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

        movaps  nb114_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvH1H1
        mulps   %xmm0,%xmm10 ## rinvH2H1
    mulps   %xmm0,%xmm11 ## rinvMH1

        ## H1 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb114_qqHH(%rsp),%xmm0
    mulps  nb114_qqHH(%rsp),%xmm1
    mulps  nb114_qqMH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb114_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb114_vctot(%rsp)

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

        mulps nb114_dxH1H1(%rsp),%xmm7
        mulps nb114_dyH1H1(%rsp),%xmm8
        mulps nb114_dzH1H1(%rsp),%xmm9
        mulps nb114_dxH2H1(%rsp),%xmm10
        mulps nb114_dyH2H1(%rsp),%xmm11
        mulps nb114_dzH2H1(%rsp),%xmm12
        mulps nb114_dxMH1(%rsp),%xmm13
        mulps nb114_dyMH1(%rsp),%xmm14
        mulps nb114_dzMH1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb114_fixH1(%rsp),%xmm7
    addps nb114_fiyH1(%rsp),%xmm8
    addps nb114_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb114_fixH2(%rsp),%xmm10
    addps nb114_fiyH2(%rsp),%xmm11
    addps nb114_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb114_fixM(%rsp),%xmm13
    addps nb114_fiyM(%rsp),%xmm14
    addps nb114_fizM(%rsp),%xmm15

    movaps %xmm7,nb114_fixH1(%rsp)
    movaps %xmm8,nb114_fiyH1(%rsp)
    movaps %xmm9,nb114_fizH1(%rsp)
    movaps %xmm10,nb114_fixH2(%rsp)
    movaps %xmm11,nb114_fiyH2(%rsp)
    movaps %xmm12,nb114_fizH2(%rsp)
    movaps %xmm13,nb114_fixM(%rsp)
    movaps %xmm14,nb114_fiyM(%rsp)
    movaps %xmm15,nb114_fizM(%rsp)

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

    subps nb114_ixH1(%rsp),%xmm0
    subps nb114_iyH1(%rsp),%xmm1
    subps nb114_izH1(%rsp),%xmm2
    subps nb114_ixH2(%rsp),%xmm3
    subps nb114_iyH2(%rsp),%xmm4
    subps nb114_izH2(%rsp),%xmm5
    subps nb114_ixM(%rsp),%xmm6
    subps nb114_iyM(%rsp),%xmm7
    subps nb114_izM(%rsp),%xmm8

        movaps %xmm0,nb114_dxH1H2(%rsp)
        movaps %xmm1,nb114_dyH1H2(%rsp)
        movaps %xmm2,nb114_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH2H2(%rsp)
        movaps %xmm4,nb114_dyH2H2(%rsp)
        movaps %xmm5,nb114_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb114_dxMH2(%rsp)
        movaps %xmm7,nb114_dyMH2(%rsp)
        movaps %xmm8,nb114_dzMH2(%rsp)
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

        movaps  nb114_three(%rsp),%xmm9
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

        movaps  nb114_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvH1H2
        mulps   %xmm0,%xmm10 ## rinvH2H2
    mulps   %xmm0,%xmm11 ## rinvMH2

        ## H2 interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb114_qqHH(%rsp),%xmm0
    mulps  nb114_qqHH(%rsp),%xmm1
    mulps  nb114_qqMH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb114_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb114_vctot(%rsp)

        ## move j H2 forces to local temp variables 
    movlps 24(%rdi,%rax,4),%xmm0    ## jxH2a jyH2a  -   -
    movlps 24(%rdi,%rcx,4),%xmm1    ## jxH2c jyH2c  -   -
    movhps 24(%rdi,%rbx,4),%xmm0    ## jxH2a jyH2a jxH2b jyH2b 
    movhps 24(%rdi,%rdx,4),%xmm1    ## jxH2c jyH2c jxH2d jyH2d 

    movss  32(%rdi,%rax,4),%xmm2    ## jzH2a  -  -  -
    movss  32(%rdi,%rcx,4),%xmm3    ## jzH2c  -  -  -
    movhps 32(%rdi,%rbx,4),%xmm2    ## jzH2a  -  jzH2b  -
    movhps 32(%rdi,%rdx,4),%xmm3    ## jzH2c  -  jzH2d -

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

        mulps nb114_dxH1H2(%rsp),%xmm7
        mulps nb114_dyH1H2(%rsp),%xmm8
        mulps nb114_dzH1H2(%rsp),%xmm9
        mulps nb114_dxH2H2(%rsp),%xmm10
        mulps nb114_dyH2H2(%rsp),%xmm11
        mulps nb114_dzH2H2(%rsp),%xmm12
        mulps nb114_dxMH2(%rsp),%xmm13
        mulps nb114_dyMH2(%rsp),%xmm14
        mulps nb114_dzMH2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb114_fixH1(%rsp),%xmm7
    addps nb114_fiyH1(%rsp),%xmm8
    addps nb114_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb114_fixH2(%rsp),%xmm10
    addps nb114_fiyH2(%rsp),%xmm11
    addps nb114_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb114_fixM(%rsp),%xmm13
    addps nb114_fiyM(%rsp),%xmm14
    addps nb114_fizM(%rsp),%xmm15

    movaps %xmm7,nb114_fixH1(%rsp)
    movaps %xmm8,nb114_fiyH1(%rsp)
    movaps %xmm9,nb114_fizH1(%rsp)
    movaps %xmm10,nb114_fixH2(%rsp)
    movaps %xmm11,nb114_fiyH2(%rsp)
    movaps %xmm12,nb114_fizH2(%rsp)
    movaps %xmm13,nb114_fixM(%rsp)
    movaps %xmm14,nb114_fiyM(%rsp)
    movaps %xmm15,nb114_fizM(%rsp)

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

    subps nb114_ixH1(%rsp),%xmm0
    subps nb114_iyH1(%rsp),%xmm1
    subps nb114_izH1(%rsp),%xmm2
    subps nb114_ixH2(%rsp),%xmm3
    subps nb114_iyH2(%rsp),%xmm4
    subps nb114_izH2(%rsp),%xmm5
    subps nb114_ixM(%rsp),%xmm6
    subps nb114_iyM(%rsp),%xmm7
    subps nb114_izM(%rsp),%xmm8

        movaps %xmm0,nb114_dxH1M(%rsp)
        movaps %xmm1,nb114_dyH1M(%rsp)
        movaps %xmm2,nb114_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb114_dxH2M(%rsp)
        movaps %xmm4,nb114_dyH2M(%rsp)
        movaps %xmm5,nb114_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb114_dxMM(%rsp)
        movaps %xmm7,nb114_dyMM(%rsp)
        movaps %xmm8,nb114_dzMM(%rsp)
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

        movaps  nb114_three(%rsp),%xmm9
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

        movaps  nb114_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvH1M 
        mulps   %xmm0,%xmm10 ## rinvH2M
    mulps   %xmm0,%xmm11 ## rinvMM

        ## M interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb114_qqMH(%rsp),%xmm0
    mulps  nb114_qqMH(%rsp),%xmm1
    mulps  nb114_qqMM(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb114_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb114_vctot(%rsp)

        ## move j M forces to local temp variables 
    movlps 36(%rdi,%rax,4),%xmm0    ## jxMa jyMa  -   -
    movlps 36(%rdi,%rcx,4),%xmm1    ## jxMc jyMc  -   -
    movhps 36(%rdi,%rbx,4),%xmm0    ## jxMa jyMa jxMb jyMb 
    movhps 36(%rdi,%rdx,4),%xmm1    ## jxMc jyMc jxMd jyMd 

    movss  44(%rdi,%rax,4),%xmm2    ## jzMa  -  -  -
    movss  44(%rdi,%rcx,4),%xmm3    ## jzMc  -  -  -
    movss  44(%rdi,%rbx,4),%xmm7    ## jzMb  -  -  -
    movss  44(%rdi,%rdx,4),%xmm8    ## jzMd  -  -  -
    movlhps %xmm7,%xmm2 ## jzMa  -  jzMb  -
    movlhps %xmm8,%xmm3 ## jzMc  -  jzMd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzMa jzMb jzMc jzMd

    ## xmm0: jxMa jyMa jxMb jyMb 
    ## xmm1: jxMc jyMc jxMd jyMd
    ## xmm2: jzMa jzMb jzMc jzMd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb114_dxH1M(%rsp),%xmm7
        mulps nb114_dyH1M(%rsp),%xmm8
        mulps nb114_dzH1M(%rsp),%xmm9
        mulps nb114_dxH2M(%rsp),%xmm10
        mulps nb114_dyH2M(%rsp),%xmm11
        mulps nb114_dzH2M(%rsp),%xmm12
        mulps nb114_dxMM(%rsp),%xmm13
        mulps nb114_dyMM(%rsp),%xmm14
        mulps nb114_dzMM(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb114_fixH1(%rsp),%xmm7
    addps nb114_fiyH1(%rsp),%xmm8
    addps nb114_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb114_fixH2(%rsp),%xmm10
    addps nb114_fiyH2(%rsp),%xmm11
    addps nb114_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb114_fixM(%rsp),%xmm13
    addps nb114_fiyM(%rsp),%xmm14
    addps nb114_fizM(%rsp),%xmm15

    movaps %xmm7,nb114_fixH1(%rsp)
    movaps %xmm8,nb114_fiyH1(%rsp)
    movaps %xmm9,nb114_fizH1(%rsp)
    movaps %xmm10,nb114_fixH2(%rsp)
    movaps %xmm11,nb114_fiyH2(%rsp)
    movaps %xmm12,nb114_fizH2(%rsp)
    movaps %xmm13,nb114_fixM(%rsp)
    movaps %xmm14,nb114_fiyM(%rsp)
    movaps %xmm15,nb114_fizM(%rsp)

    ## xmm3 = fMx , xmm4 = fMy
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fMzc fMzd

    movlps %xmm0,36(%rdi,%rax,4)
    movhps %xmm0,36(%rdi,%rbx,4)
    movlps %xmm1,36(%rdi,%rcx,4)
    movhps %xmm1,36(%rdi,%rdx,4)
    movss  %xmm2,44(%rdi,%rax,4)
    movss  %xmm3,44(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,44(%rdi,%rbx,4)
    movss  %xmm3,44(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb114_innerk(%rsp)
        jl    _nb_kernel114_x86_64_sse.nb114_single_check
        jmp   _nb_kernel114_x86_64_sse.nb114_unroll_loop
_nb_kernel114_x86_64_sse.nb114_single_check: 
        addl $4,nb114_innerk(%rsp)
        jnz   _nb_kernel114_x86_64_sse.nb114_single_loop
        jmp   _nb_kernel114_x86_64_sse.nb114_updateouterdata
_nb_kernel114_x86_64_sse.nb114_single_loop: 
        movq  nb114_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb114_innerjjnr(%rsp)

        movq nb114_pos(%rbp),%rsi
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
        movaps %xmm6,nb114_jxO(%rsp)
        movaps %xmm3,nb114_jyO(%rsp)
        movaps %xmm1,nb114_jzO(%rsp)

    movaps %xmm6,%xmm0  ## jxO
    movaps %xmm1,%xmm2  ## jzO
    movaps %xmm3,%xmm1  ## jyO
    movaps %xmm3,%xmm4  ## jyO
    movaps %xmm6,%xmm3  ## jxO
    movaps %xmm2,%xmm5  ## jzO

        ## do O and M in parallel
        subps  nb114_ixO(%rsp),%xmm0
        subps  nb114_iyO(%rsp),%xmm1
        subps  nb114_izO(%rsp),%xmm2
        subps  nb114_ixM(%rsp),%xmm3
        subps  nb114_iyM(%rsp),%xmm4
        subps  nb114_izM(%rsp),%xmm5

        movaps %xmm0,nb114_dxOO(%rsp)
        movaps %xmm1,nb114_dyOO(%rsp)
        movaps %xmm2,nb114_dzOO(%rsp)
        movaps %xmm3,nb114_dxMM(%rsp)
        movaps %xmm4,nb114_dyMM(%rsp)
        movaps %xmm5,nb114_dzMM(%rsp)

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
        ## Save M data 
        movaps %xmm4,nb114_rsqMM(%rsp)

        ## do 1/x for O and 1/sqrt(x) for M
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb114_two(%rsp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb114_three(%rsp),%xmm7
        mulss  %xmm1,%xmm2      ## 1/r2

        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## 1/r6
        mulps   nb114_half(%rsp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb114_rinvMM(%rsp)

        mulss  %xmm2,%xmm2      ## 1/r12
        mulss  nb114_c6(%rsp),%xmm1
        mulss  nb114_c12(%rsp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb114_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb114_Vvdwtot(%rsp)
        mulss  nb114_six(%rsp),%xmm1
        mulss  nb114_twelve(%rsp),%xmm2
        subss  %xmm1,%xmm2
        mulss  %xmm2,%xmm0      ## fscal
        movss  %xmm0,%xmm1
        movss  %xmm0,%xmm2
        mulss  nb114_dxOO(%rsp),%xmm0
        mulss  nb114_dyOO(%rsp),%xmm1
        mulss  nb114_dzOO(%rsp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addss   %xmm0,%xmm3
        addss   %xmm1,%xmm4
        addss   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%rsp)
        movaps  %xmm4,nb114_fjyO(%rsp)
        movaps  %xmm5,nb114_fjzO(%rsp)
        addss   nb114_fixO(%rsp),%xmm0
        addss   nb114_fiyO(%rsp),%xmm1
        addss   nb114_fizO(%rsp),%xmm2
        movss  %xmm0,nb114_fixO(%rsp)
        movss  %xmm1,nb114_fiyO(%rsp)
        movss  %xmm2,nb114_fizO(%rsp)

        ## do  M coulomb interaction
        movaps nb114_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqMH(%rsp),%xmm3
        movhps  nb114_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul    
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%rsp),%xmm4
        movaps %xmm4,nb114_vctot(%rsp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxMM(%rsp),%xmm0
        mulps   nb114_dyMM(%rsp),%xmm1
        mulps   nb114_dzMM(%rsp),%xmm2
        ## update forces M - j water 
        movaps  nb114_fjxO(%rsp),%xmm3
        movaps  nb114_fjyO(%rsp),%xmm4
        movaps  nb114_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%rsp)
        movaps  %xmm4,nb114_fjyO(%rsp)
        movaps  %xmm5,nb114_fjzO(%rsp)
        addps   nb114_fixM(%rsp),%xmm0
        addps   nb114_fiyM(%rsp),%xmm1
        addps   nb114_fizM(%rsp),%xmm2
        movaps  %xmm0,nb114_fixM(%rsp)
        movaps  %xmm1,nb114_fiyM(%rsp)
        movaps  %xmm2,nb114_fizM(%rsp)

        ## i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb114_jxO(%rsp),%xmm0
    movaps  nb114_jyO(%rsp),%xmm1
    movaps  nb114_jzO(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5
        subps   nb114_ixH1(%rsp),%xmm0
        subps   nb114_iyH1(%rsp),%xmm1
        subps   nb114_izH1(%rsp),%xmm2
        subps   nb114_ixH2(%rsp),%xmm3
        subps   nb114_iyH2(%rsp),%xmm4
        subps   nb114_izH2(%rsp),%xmm5
        movaps %xmm0,nb114_dxH1H1(%rsp)
        movaps %xmm1,nb114_dyH1H1(%rsp)
        movaps %xmm2,nb114_dzH1H1(%rsp)
        movaps %xmm3,nb114_dxH2H2(%rsp)
        movaps %xmm4,nb114_dyH2H2(%rsp)
        movaps %xmm5,nb114_dzH2H2(%rsp)
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
        movaps  %xmm0,nb114_rsqH1H1(%rsp)
        movaps  %xmm4,nb114_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114_half(%rsp),%xmm3   ## rinvH1H1
        mulps   nb114_half(%rsp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb114_rinvH1H1(%rsp)
        movaps  %xmm7,nb114_rinvH2H2(%rsp)

        ## Do H1 coulomb interaction
        movaps nb114_rinvH1H1(%rsp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqHH(%rsp),%xmm3
        movhps  nb114_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%rsp),%xmm4
        movaps %xmm4,nb114_vctot(%rsp)


        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxH1H1(%rsp),%xmm0
        mulps   nb114_dyH1H1(%rsp),%xmm1
        mulps   nb114_dzH1H1(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb114_fjxO(%rsp),%xmm3
        movaps  nb114_fjyO(%rsp),%xmm4
        movaps  nb114_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%rsp)
        movaps  %xmm4,nb114_fjyO(%rsp)
        movaps  %xmm5,nb114_fjzO(%rsp)
        addps   nb114_fixH1(%rsp),%xmm0
        addps   nb114_fiyH1(%rsp),%xmm1
        addps   nb114_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb114_fixH1(%rsp)
        movaps  %xmm1,nb114_fiyH1(%rsp)
        movaps  %xmm2,nb114_fizH1(%rsp)

        ## H2 Coulomb
        movaps nb114_rinvH2H2(%rsp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv 
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114_qqHH(%rsp),%xmm3
        movhps  nb114_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul
        mulps  %xmm4,%xmm0
        addps  nb114_vctot(%rsp),%xmm4   ## local vctot summation variable
        movaps %xmm4,nb114_vctot(%rsp)

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb114_dxH2H2(%rsp),%xmm0
        mulps   nb114_dyH2H2(%rsp),%xmm1
        mulps   nb114_dzH2H2(%rsp),%xmm2
        ## update forces H2 - j water 
        movaps  nb114_fjxO(%rsp),%xmm3
        movaps  nb114_fjyO(%rsp),%xmm4
        movaps  nb114_fjzO(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb114_fjxO(%rsp)
        movaps  %xmm4,nb114_fjyO(%rsp)
        movaps  %xmm5,nb114_fjzO(%rsp)
        addps   nb114_fixH2(%rsp),%xmm0
        addps   nb114_fiyH2(%rsp),%xmm1
        addps   nb114_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb114_fixH2(%rsp)
        movaps  %xmm1,nb114_fiyH2(%rsp)
        movaps  %xmm2,nb114_fizH2(%rsp)

        movq    nb114_faction(%rbp),%rsi
        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb114_fjxO(%rsp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb114_fjyO(%rsp),%xmm1   ## Oy H1y H2y My
        movaps  nb114_fjzO(%rsp),%xmm2   ## Oz H1z H2z Mz

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

        decl nb114_innerk(%rsp)
        jz    _nb_kernel114_x86_64_sse.nb114_updateouterdata
        jmp   _nb_kernel114_x86_64_sse.nb114_single_loop
_nb_kernel114_x86_64_sse.nb114_updateouterdata: 
        movl  nb114_ii3(%rsp),%ecx
        movq  nb114_faction(%rbp),%rdi
        movq  nb114_fshift(%rbp),%rsi
        movl  nb114_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb114_fixO(%rsp),%xmm0
        movaps nb114_fiyO(%rsp),%xmm1
        movaps nb114_fizO(%rsp),%xmm2

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
        movaps nb114_fixH1(%rsp),%xmm0
        movaps nb114_fiyH1(%rsp),%xmm1
        movaps nb114_fizH1(%rsp),%xmm2

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
        movaps nb114_fixH2(%rsp),%xmm0
        movaps nb114_fiyH2(%rsp),%xmm1
        movaps nb114_fizH2(%rsp),%xmm2

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
        movaps nb114_fixM(%rsp),%xmm0
        movaps nb114_fiyM(%rsp),%xmm1
        movaps nb114_fizM(%rsp),%xmm2

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
        movl nb114_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb114_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb114_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb114_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb114_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb114_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb114_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114_x86_64_sse.nb114_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114_n(%rsp)
        jmp _nb_kernel114_x86_64_sse.nb114_outer
_nb_kernel114_x86_64_sse.nb114_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114_x86_64_sse.nb114_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114_x86_64_sse.nb114_threadloop
_nb_kernel114_x86_64_sse.nb114_end: 


        emms

        movl nb114_nouter(%rsp),%eax
        movl nb114_ninner(%rsp),%ebx
        movq nb114_outeriter(%rbp),%rcx
        movq nb114_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1880,%rsp

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel114nf_x86_64_sse
.globl _nb_kernel114nf_x86_64_sse
nb_kernel114nf_x86_64_sse:      
_nb_kernel114nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb114nf_fshift, 16
.set nb114nf_gid, 24
.set nb114nf_pos, 32
.set nb114nf_faction, 40
.set nb114nf_charge, 48
.set nb114nf_p_facel, 56
.set nb114nf_argkrf, 64
.set nb114nf_argcrf, 72
.set nb114nf_Vc, 80
.set nb114nf_type, 88
.set nb114nf_p_ntype, 96
.set nb114nf_vdwparam, 104
.set nb114nf_Vvdw, 112
.set nb114nf_p_tabscale, 120
.set nb114nf_VFtab, 128
.set nb114nf_invsqrta, 136
.set nb114nf_dvda, 144
.set nb114nf_p_gbtabscale, 152
.set nb114nf_GBtab, 160
.set nb114nf_p_nthreads, 168
.set nb114nf_count, 176
.set nb114nf_mtx, 184
.set nb114nf_outeriter, 192
.set nb114nf_inneriter, 200
.set nb114nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb114nf_ixO, 0
.set nb114nf_iyO, 16
.set nb114nf_izO, 32
.set nb114nf_ixH1, 48
.set nb114nf_iyH1, 64
.set nb114nf_izH1, 80
.set nb114nf_ixH2, 96
.set nb114nf_iyH2, 112
.set nb114nf_izH2, 128
.set nb114nf_ixM, 144
.set nb114nf_iyM, 160
.set nb114nf_izM, 176
.set nb114nf_jxO, 192
.set nb114nf_jyO, 208
.set nb114nf_jzO, 224
.set nb114nf_jxH1, 240
.set nb114nf_jyH1, 256
.set nb114nf_jzH1, 272
.set nb114nf_jxH2, 288
.set nb114nf_jyH2, 304
.set nb114nf_jzH2, 320
.set nb114nf_jxM, 336
.set nb114nf_jyM, 352
.set nb114nf_jzM, 368
.set nb114nf_qqMM, 384
.set nb114nf_qqMH, 400
.set nb114nf_qqHH, 416
.set nb114nf_two, 432
.set nb114nf_c6, 448
.set nb114nf_c12, 464
.set nb114nf_vctot, 480
.set nb114nf_Vvdwtot, 496
.set nb114nf_half, 512
.set nb114nf_three, 528
.set nb114nf_rsqOO, 544
.set nb114nf_rsqH1H1, 560
.set nb114nf_rsqH1H2, 576
.set nb114nf_rsqH1M, 592
.set nb114nf_rsqH2H1, 608
.set nb114nf_rsqH2H2, 624
.set nb114nf_rsqH2M, 640
.set nb114nf_rsqMH1, 656
.set nb114nf_rsqMH2, 672
.set nb114nf_rsqMM, 688
.set nb114nf_rinvsqOO, 704
.set nb114nf_rinvH1H1, 720
.set nb114nf_rinvH1H2, 736
.set nb114nf_rinvH1M, 752
.set nb114nf_rinvH2H1, 768
.set nb114nf_rinvH2H2, 784
.set nb114nf_rinvH2M, 800
.set nb114nf_rinvMH1, 816
.set nb114nf_rinvMH2, 832
.set nb114nf_rinvMM, 848
.set nb114nf_nri, 864
.set nb114nf_iinr, 872
.set nb114nf_jindex, 880
.set nb114nf_jjnr, 888
.set nb114nf_shift, 896
.set nb114nf_shiftvec, 904
.set nb114nf_facel, 912
.set nb114nf_innerjjnr, 920
.set nb114nf_innerk, 928
.set nb114nf_is3, 932
.set nb114nf_ii3, 936
.set nb114nf_n, 940
.set nb114nf_nn1, 944
.set nb114nf_nouter, 948
.set nb114nf_ninner, 952
.set nb114nf_salign, 956

        push %rbp
        movq %rsp,%rbp
        push %rbx


        push %r12
        push %r13
        push %r14
        push %r15

        subq $968,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb114nf_nouter(%rsp)
        movl %eax,nb114nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb114nf_nri(%rsp)
        movq %rsi,nb114nf_iinr(%rsp)
        movq %rdx,nb114nf_jindex(%rsp)
        movq %rcx,nb114nf_jjnr(%rsp)
        movq %r8,nb114nf_shift(%rsp)
        movq %r9,nb114nf_shiftvec(%rsp)
        movq nb114nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb114nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb114nf_half(%rsp)
        movss nb114nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb114nf_half(%rsp)
        movaps %xmm2,nb114nf_two(%rsp)
        movaps %xmm3,nb114nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb114nf_iinr(%rsp),%rcx     ## rcx = pointer into iinr[]
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb114nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm5
        movss 12(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movq nb114nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb114nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb114nf_qqMM(%rsp)
        movaps %xmm4,nb114nf_qqMH(%rsp)
        movaps %xmm5,nb114nf_qqHH(%rsp)

        xorps %xmm0,%xmm0
        movq  nb114nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb114nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx ## rcx = ntia = 2*ntype*type[ii0] 
        addq  %rcx,%rdx
        movq  nb114nf_vdwparam(%rbp),%rax
        movlps (%rax,%rdx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb114nf_c6(%rsp)
        movaps %xmm1,nb114nf_c12(%rsp)

_nb_kernel114nf_x86_64_sse.nb114nf_threadloop: 
        movq  nb114nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel114nf_x86_64_sse.nb114nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114nf_x86_64_sse.nb114nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114nf_n(%rsp)
        movl %ebx,nb114nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114nf_x86_64_sse.nb114nf_outerstart
        jmp _nb_kernel114nf_x86_64_sse.nb114nf_end

_nb_kernel114nf_x86_64_sse.nb114nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114nf_nouter(%rsp),%ebx
        movl %ebx,nb114nf_nouter(%rsp)

_nb_kernel114nf_x86_64_sse.nb114nf_outer: 
        movq  nb114nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb114nf_is3(%rsp)            ## store is3 

        movq  nb114nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb114nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb114nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb114nf_ii3(%rsp)

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
        movaps %xmm3,nb114nf_ixO(%rsp)
        movaps %xmm4,nb114nf_iyO(%rsp)
        movaps %xmm5,nb114nf_izO(%rsp)
        movaps %xmm6,nb114nf_ixH1(%rsp)
        movaps %xmm7,nb114nf_iyH1(%rsp)

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
        movaps %xmm6,nb114nf_izH1(%rsp)
        movaps %xmm0,nb114nf_ixH2(%rsp)
        movaps %xmm1,nb114nf_iyH2(%rsp)
        movaps %xmm2,nb114nf_izH2(%rsp)
        movaps %xmm3,nb114nf_ixM(%rsp)
        movaps %xmm4,nb114nf_iyM(%rsp)
        movaps %xmm5,nb114nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb114nf_vctot(%rsp)
        movaps %xmm4,nb114nf_Vvdwtot(%rsp)

        movq  nb114nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb114nf_pos(%rbp),%rsi
        movq  nb114nf_faction(%rbp),%rdi
        movq  nb114nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb114nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb114nf_ninner(%rsp),%ecx
        movl  %ecx,nb114nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb114nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel114nf_x86_64_sse.nb114nf_unroll_loop
        jmp   _nb_kernel114nf_x86_64_sse.nb114nf_single_check
_nb_kernel114nf_x86_64_sse.nb114nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb114nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb114nf_innerjjnr(%rsp)             ## advance pointer (unroll 4) 

        movq nb114nf_pos(%rbp),%rsi     ## base of pos[] 

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
        movaps %xmm0,nb114nf_jxO(%rsp)
        movaps %xmm1,nb114nf_jyO(%rsp)
        movaps %xmm2,nb114nf_jzO(%rsp)
        movaps %xmm3,nb114nf_jxH1(%rsp)

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
        movaps %xmm0,nb114nf_jyH1(%rsp)
        movaps %xmm1,nb114nf_jzH1(%rsp)
        movaps %xmm2,nb114nf_jxH2(%rsp)
        movaps %xmm3,nb114nf_jyH2(%rsp)

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
        movaps %xmm0,nb114nf_jzH2(%rsp)
        movaps %xmm1,nb114nf_jxM(%rsp)
        movaps %xmm2,nb114nf_jyM(%rsp)
        movaps %xmm3,nb114nf_jzM(%rsp)

        ## start calculating pairwise distances
        movaps nb114nf_ixO(%rsp),%xmm0
        movaps nb114nf_iyO(%rsp),%xmm1
        movaps nb114nf_izO(%rsp),%xmm2
        movaps nb114nf_ixH1(%rsp),%xmm3
        movaps nb114nf_iyH1(%rsp),%xmm4
        movaps nb114nf_izH1(%rsp),%xmm5
        subps  nb114nf_jxO(%rsp),%xmm0
        subps  nb114nf_jyO(%rsp),%xmm1
        subps  nb114nf_jzO(%rsp),%xmm2
        subps  nb114nf_jxH1(%rsp),%xmm3
        subps  nb114nf_jyH1(%rsp),%xmm4
        subps  nb114nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb114nf_rsqOO(%rsp)
        movaps %xmm3,nb114nf_rsqH1H1(%rsp)

        movaps nb114nf_ixH1(%rsp),%xmm0
        movaps nb114nf_iyH1(%rsp),%xmm1
        movaps nb114nf_izH1(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH2(%rsp),%xmm0
        subps  nb114nf_jyH2(%rsp),%xmm1
        subps  nb114nf_jzH2(%rsp),%xmm2
        subps  nb114nf_jxM(%rsp),%xmm3
        subps  nb114nf_jyM(%rsp),%xmm4
        subps  nb114nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH1H2(%rsp)
        movaps %xmm3,nb114nf_rsqH1M(%rsp)

        movaps nb114nf_ixH2(%rsp),%xmm0
        movaps nb114nf_iyH2(%rsp),%xmm1
        movaps nb114nf_izH2(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH1(%rsp),%xmm0
        subps  nb114nf_jyH1(%rsp),%xmm1
        subps  nb114nf_jzH1(%rsp),%xmm2
        subps  nb114nf_jxH2(%rsp),%xmm3
        subps  nb114nf_jyH2(%rsp),%xmm4
        subps  nb114nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH2H1(%rsp)
        movaps %xmm3,nb114nf_rsqH2H2(%rsp)

        movaps nb114nf_ixH2(%rsp),%xmm0
        movaps nb114nf_iyH2(%rsp),%xmm1
        movaps nb114nf_izH2(%rsp),%xmm2
        movaps nb114nf_ixM(%rsp),%xmm3
        movaps nb114nf_iyM(%rsp),%xmm4
        movaps nb114nf_izM(%rsp),%xmm5
        subps  nb114nf_jxM(%rsp),%xmm0
        subps  nb114nf_jyM(%rsp),%xmm1
        subps  nb114nf_jzM(%rsp),%xmm2
        subps  nb114nf_jxH1(%rsp),%xmm3
        subps  nb114nf_jyH1(%rsp),%xmm4
        subps  nb114nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb114nf_rsqH2M(%rsp)
        movaps %xmm4,nb114nf_rsqMH1(%rsp)

        movaps nb114nf_ixM(%rsp),%xmm0
        movaps nb114nf_iyM(%rsp),%xmm1
        movaps nb114nf_izM(%rsp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb114nf_jxH2(%rsp),%xmm0
        subps  nb114nf_jyH2(%rsp),%xmm1
        subps  nb114nf_jzH2(%rsp),%xmm2
        subps  nb114nf_jxM(%rsp),%xmm3
        subps  nb114nf_jyM(%rsp),%xmm4
        subps  nb114nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb114nf_rsqMH2(%rsp)
        movaps %xmm4,nb114nf_rsqMM(%rsp)

        ## start by doing reciprocal for OO
        movaps  nb114nf_rsqOO(%rsp),%xmm7
        rcpps   %xmm7,%xmm2
        movaps  nb114nf_two(%rsp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps %xmm2,nb114nf_rinvsqOO(%rsp)

        ## next step is invsqrt - do two at a time.
        rsqrtps nb114nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb114nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH1H1(%rsp),%xmm1
        mulps   nb114nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%rsp),%xmm3   ## rinvH1H1 
        mulps   nb114nf_half(%rsp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb114nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb114nf_rinvH1H2(%rsp)

        rsqrtps nb114nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb114nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH1M(%rsp),%xmm1
        mulps   nb114nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%rsp),%xmm3
        mulps   nb114nf_half(%rsp),%xmm7
        movaps  %xmm3,nb114nf_rinvH1M(%rsp)
        movaps  %xmm7,nb114nf_rinvH2H1(%rsp)

        rsqrtps nb114nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb114nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqH2H2(%rsp),%xmm1
        mulps   nb114nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%rsp),%xmm3
        mulps   nb114nf_half(%rsp),%xmm7
        movaps  %xmm3,nb114nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb114nf_rinvH2M(%rsp)

        rsqrtps nb114nf_rsqMH1(%rsp),%xmm1
        rsqrtps nb114nf_rsqMH2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb114nf_rsqMH1(%rsp),%xmm1
        mulps   nb114nf_rsqMH2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%rsp),%xmm3
        mulps   nb114nf_half(%rsp),%xmm7
        movaps  %xmm3,nb114nf_rinvMH1(%rsp)
        movaps  %xmm7,nb114nf_rinvMH2(%rsp)

        rsqrtps nb114nf_rsqMM(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb114nf_three(%rsp),%xmm3
        mulps   nb114nf_rsqMM(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb114nf_half(%rsp),%xmm3
        movaps  %xmm3,nb114nf_rinvMM(%rsp)

        ## start with OO LJ interaction
        movaps nb114nf_rinvsqOO(%rsp),%xmm0
        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1      ## rinv4
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb114nf_c6(%rsp),%xmm1
        mulps  nb114nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm4
        subps  %xmm1,%xmm4
        addps  nb114nf_Vvdwtot(%rsp),%xmm4
        movaps %xmm4,nb114nf_Vvdwtot(%rsp)

        ## Coulomb interactions 
        ## all H-H interactions
        movaps nb114nf_rinvH1H1(%rsp),%xmm0
        addps  nb114nf_rinvH1H2(%rsp),%xmm0
        addps  nb114nf_rinvH2H1(%rsp),%xmm0
        addps  nb114nf_rinvH2H2(%rsp),%xmm0
        mulps  nb114nf_qqHH(%rsp),%xmm0
        ## all M-H interactions
        movaps nb114nf_rinvH1M(%rsp),%xmm1
        addps  nb114nf_rinvH2M(%rsp),%xmm1
        addps  nb114nf_rinvMH1(%rsp),%xmm1
        addps  nb114nf_rinvMH2(%rsp),%xmm1
        mulps  nb114nf_qqMH(%rsp),%xmm1
        ## The M-M interaction
        movaps nb114nf_rinvMM(%rsp),%xmm2
        mulps  nb114nf_qqMM(%rsp),%xmm2
        addps  %xmm1,%xmm0
        addps  nb114nf_vctot(%rsp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb114nf_vctot(%rsp)
        ## should we do one more iteration? 
        subl $4,nb114nf_innerk(%rsp)
        jl    _nb_kernel114nf_x86_64_sse.nb114nf_single_check
        jmp   _nb_kernel114nf_x86_64_sse.nb114nf_unroll_loop
_nb_kernel114nf_x86_64_sse.nb114nf_single_check: 
        addl $4,nb114nf_innerk(%rsp)
        jnz   _nb_kernel114nf_x86_64_sse.nb114nf_single_loop
        jmp   _nb_kernel114nf_x86_64_sse.nb114nf_updateouterdata
_nb_kernel114nf_x86_64_sse.nb114nf_single_loop: 
        movq  nb114nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb114nf_innerjjnr(%rsp)

        movq nb114nf_pos(%rbp),%rsi
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
        movaps %xmm6,nb114nf_jxO(%rsp)
        movaps %xmm3,nb114nf_jyO(%rsp)
        movaps %xmm1,nb114nf_jzO(%rsp)

        ## do O and M in parallel
        movaps nb114nf_ixO(%rsp),%xmm0
        movaps nb114nf_iyO(%rsp),%xmm1
        movaps nb114nf_izO(%rsp),%xmm2
        movaps nb114nf_ixM(%rsp),%xmm3
        movaps nb114nf_iyM(%rsp),%xmm4
        movaps nb114nf_izM(%rsp),%xmm5
        subps  nb114nf_jxO(%rsp),%xmm0
        subps  nb114nf_jyO(%rsp),%xmm1
        subps  nb114nf_jzO(%rsp),%xmm2
        subps  nb114nf_jxO(%rsp),%xmm3
        subps  nb114nf_jyO(%rsp),%xmm4
        subps  nb114nf_jzO(%rsp),%xmm5

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
        ## Save M data 
        movaps %xmm4,nb114nf_rsqMM(%rsp)

        ## do 1/x for O and 1/sqrt(x) for M
        rcpss  %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movss  nb114nf_two(%rsp),%xmm2
        movaps  %xmm5,%xmm6
        mulss  %xmm1,%xmm0
        mulps   %xmm5,%xmm5
        subss  %xmm0,%xmm2
        movaps  nb114nf_three(%rsp),%xmm7
        mulss  %xmm1,%xmm2      ## 1/r2

        mulps   %xmm4,%xmm5
        movss  %xmm2,%xmm0
        subps   %xmm5,%xmm7
        mulss  %xmm2,%xmm2
        mulps   %xmm6,%xmm7
        mulss  %xmm0,%xmm2      ## 1/r6
        mulps   nb114nf_half(%rsp),%xmm7   ## rinv iH1 - j water 
        movss  %xmm2,%xmm1
        movaps %xmm7,nb114nf_rinvMM(%rsp)

        mulss  %xmm2,%xmm2      ## 1/r12
        mulss  nb114nf_c6(%rsp),%xmm1
        mulss  nb114nf_c12(%rsp),%xmm2
        movss  %xmm2,%xmm3
        subss  %xmm1,%xmm3
        addss  nb114nf_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb114nf_Vvdwtot(%rsp)

        ## do  M coulomb interaction
        movaps nb114nf_rinvMM(%rsp),%xmm0
        movaps %xmm0,%xmm4      ## xmm4=rinv

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114nf_qqMH(%rsp),%xmm3
        movhps  nb114nf_qqMM(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm3,%xmm4      ## xmm4=voul    
        addps  nb114nf_vctot(%rsp),%xmm4
        movaps %xmm4,nb114nf_vctot(%rsp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb114nf_ixH1(%rsp),%xmm0
        movaps  nb114nf_iyH1(%rsp),%xmm1
        movaps  nb114nf_izH1(%rsp),%xmm2
        movaps  nb114nf_ixH2(%rsp),%xmm3
        movaps  nb114nf_iyH2(%rsp),%xmm4
        movaps  nb114nf_izH2(%rsp),%xmm5
        subps   nb114nf_jxO(%rsp),%xmm0
        subps   nb114nf_jyO(%rsp),%xmm1
        subps   nb114nf_jzO(%rsp),%xmm2
        subps   nb114nf_jxO(%rsp),%xmm3
        subps   nb114nf_jyO(%rsp),%xmm4
        subps   nb114nf_jzO(%rsp),%xmm5
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
        movaps  %xmm0,nb114nf_rsqH1H1(%rsp)
        movaps  %xmm4,nb114nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb114nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb114nf_half(%rsp),%xmm3   ## rinvH1H1
        mulps   nb114nf_half(%rsp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb114nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb114nf_rinvH2H2(%rsp)

        ## Do H1 & H2 coulomb interaction
        movaps nb114nf_rinvH1H1(%rsp),%xmm0
        addps  nb114nf_rinvH2H2(%rsp),%xmm0

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb114nf_qqHH(%rsp),%xmm3
        movhps  nb114nf_qqMH(%rsp),%xmm3
        shufps $193,%xmm3,%xmm3 ## 11000001 

        mulps  %xmm3,%xmm0      ## xmm4=voul
        addps  nb114nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb114nf_vctot(%rsp)

        decl nb114nf_innerk(%rsp)
        jz    _nb_kernel114nf_x86_64_sse.nb114nf_updateouterdata
        jmp   _nb_kernel114nf_x86_64_sse.nb114nf_single_loop
_nb_kernel114nf_x86_64_sse.nb114nf_updateouterdata: 
        ## get n from stack
        movl nb114nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb114nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb114nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb114nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb114nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb114nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb114nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114nf_x86_64_sse.nb114nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114nf_n(%rsp)
        jmp _nb_kernel114nf_x86_64_sse.nb114nf_outer
_nb_kernel114nf_x86_64_sse.nb114nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114nf_x86_64_sse.nb114nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114nf_x86_64_sse.nb114nf_threadloop
_nb_kernel114nf_x86_64_sse.nb114nf_end: 


        emms

        movl nb114nf_nouter(%rsp),%eax
        movl nb114nf_ninner(%rsp),%ebx
        movq nb114nf_outeriter(%rbp),%rcx
        movq nb114nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $968,%rsp

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


