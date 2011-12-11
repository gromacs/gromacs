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





.globl nb_kernel104_x86_64_sse
.globl _nb_kernel104_x86_64_sse
nb_kernel104_x86_64_sse:        
_nb_kernel104_x86_64_sse:       
.set nb104_fshift, 16
.set nb104_gid, 24
.set nb104_pos, 32
.set nb104_faction, 40
.set nb104_charge, 48
.set nb104_p_facel, 56
.set nb104_argkrf, 64
.set nb104_argcrf, 72
.set nb104_Vc, 80
.set nb104_type, 88
.set nb104_p_ntype, 96
.set nb104_vdwparam, 104
.set nb104_Vvdw, 112
.set nb104_p_tabscale, 120
.set nb104_VFtab, 128
.set nb104_invsqrta, 136
.set nb104_dvda, 144
.set nb104_p_gbtabscale, 152
.set nb104_GBtab, 160
.set nb104_p_nthreads, 168
.set nb104_count, 176
.set nb104_mtx, 184
.set nb104_outeriter, 192
.set nb104_inneriter, 200
.set nb104_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb104_ixH1, 0
.set nb104_iyH1, 16
.set nb104_izH1, 32
.set nb104_ixH2, 48
.set nb104_iyH2, 64
.set nb104_izH2, 80
.set nb104_ixM, 96
.set nb104_iyM, 112
.set nb104_izM, 128
.set nb104_jxH1, 144
.set nb104_jyH1, 160
.set nb104_jzH1, 176
.set nb104_jxH2, 192
.set nb104_jyH2, 208
.set nb104_jzH2, 224
.set nb104_jxM, 240
.set nb104_jyM, 256
.set nb104_jzM, 272
.set nb104_dxH1H1, 288
.set nb104_dyH1H1, 304
.set nb104_dzH1H1, 320
.set nb104_dxH1H2, 336
.set nb104_dyH1H2, 352
.set nb104_dzH1H2, 368
.set nb104_dxH1M, 384
.set nb104_dyH1M, 400
.set nb104_dzH1M, 416
.set nb104_dxH2H1, 432
.set nb104_dyH2H1, 448
.set nb104_dzH2H1, 464
.set nb104_dxH2H2, 480
.set nb104_dyH2H2, 496
.set nb104_dzH2H2, 512
.set nb104_dxH2M, 528
.set nb104_dyH2M, 544
.set nb104_dzH2M, 560
.set nb104_dxMH1, 576
.set nb104_dyMH1, 592
.set nb104_dzMH1, 608
.set nb104_dxMH2, 624
.set nb104_dyMH2, 640
.set nb104_dzMH2, 656
.set nb104_dxMM, 672
.set nb104_dyMM, 688
.set nb104_dzMM, 704
.set nb104_qqHH, 720
.set nb104_qqMH, 736
.set nb104_qqMM, 752
.set nb104_vctot, 768
.set nb104_fixH1, 784
.set nb104_fiyH1, 800
.set nb104_fizH1, 816
.set nb104_fixH2, 832
.set nb104_fiyH2, 848
.set nb104_fizH2, 864
.set nb104_fixM, 880
.set nb104_fiyM, 896
.set nb104_fizM, 912
.set nb104_fjxH1, 928
.set nb104_fjyH1, 944
.set nb104_fjzH1, 960
.set nb104_fjxH2, 976
.set nb104_fjyH2, 992
.set nb104_fjzH2, 1008
.set nb104_fjxM, 1024
.set nb104_fjyM, 1040
.set nb104_fjzM, 1056
.set nb104_half, 1072
.set nb104_three, 1088
.set nb104_rsqH1H1, 1104
.set nb104_rsqH1H2, 1120
.set nb104_rsqH1M, 1136
.set nb104_rsqH2H1, 1152
.set nb104_rsqH2H2, 1168
.set nb104_rsqH2M, 1184
.set nb104_rsqMH1, 1200
.set nb104_rsqMH2, 1216
.set nb104_rsqMM, 1232
.set nb104_rinvH1H1, 1248
.set nb104_rinvH1H2, 1264
.set nb104_rinvH1M, 1280
.set nb104_rinvH2H1, 1296
.set nb104_rinvH2H2, 1312
.set nb104_rinvH2M, 1328
.set nb104_rinvMH1, 1344
.set nb104_rinvMH2, 1360
.set nb104_rinvMM, 1376
.set nb104_nri, 1392
.set nb104_innerjjnr, 1400
.set nb104_iinr, 1408
.set nb104_jindex, 1416
.set nb104_jjnr, 1424
.set nb104_shift, 1432
.set nb104_shiftvec, 1440
.set nb104_facel, 1448
.set nb104_is3, 1456
.set nb104_ii3, 1464
.set nb104_innerk, 1472
.set nb104_n, 1480
.set nb104_nn1, 1484
.set nb104_nouter, 1488
.set nb104_ninner, 1492

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1512,%rsp
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb104_nouter(%rsp)
        movl %eax,nb104_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb104_nri(%rsp)
        movq %rsi,nb104_iinr(%rsp)
        movq %rdx,nb104_jindex(%rsp)
        movq %rcx,nb104_jjnr(%rsp)
        movq %r8,nb104_shift(%rsp)
        movq %r9,nb104_shiftvec(%rsp)
        movq nb104_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb104_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb104_half(%rsp)
        movss nb104_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb104_half(%rsp)
        movaps %xmm3,nb104_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb104_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb104_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb104_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb104_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104_qqHH(%rsp)
        movaps %xmm4,nb104_qqMH(%rsp)
        movaps %xmm5,nb104_qqMM(%rsp)

_nb_kernel104_x86_64_sse.nb104_threadloop: 
        movq  nb104_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel104_x86_64_sse.nb104_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104_x86_64_sse.nb104_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104_n(%rsp)
        movl %ebx,nb104_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104_x86_64_sse.nb104_outerstart
        jmp _nb_kernel104_x86_64_sse.nb104_end

_nb_kernel104_x86_64_sse.nb104_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104_nouter(%rsp),%ebx
        movl %ebx,nb104_nouter(%rsp)

_nb_kernel104_x86_64_sse.nb104_outer: 
        movq  nb104_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb104_is3(%rsp)      ## store is3 

        movq  nb104_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb104_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb104_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb104_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104_ixH1(%rsp)
        movaps %xmm4,nb104_iyH1(%rsp)
        movaps %xmm5,nb104_izH1(%rsp)

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
        movaps %xmm0,nb104_ixH2(%rsp)
        movaps %xmm1,nb104_iyH2(%rsp)
        movaps %xmm2,nb104_izH2(%rsp)
        movaps %xmm3,nb104_ixM(%rsp)
        movaps %xmm4,nb104_iyM(%rsp)
        movaps %xmm5,nb104_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb104_vctot(%rsp)
        movaps %xmm4,nb104_fixH1(%rsp)
        movaps %xmm4,nb104_fiyH1(%rsp)
        movaps %xmm4,nb104_fizH1(%rsp)
        movaps %xmm4,nb104_fixH2(%rsp)
        movaps %xmm4,nb104_fiyH2(%rsp)
        movaps %xmm4,nb104_fizH2(%rsp)
        movaps %xmm4,nb104_fixM(%rsp)
        movaps %xmm4,nb104_fiyM(%rsp)
        movaps %xmm4,nb104_fizM(%rsp)

        movq  nb104_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb104_pos(%rbp),%rsi
        movq  nb104_faction(%rbp),%rdi
        movq  nb104_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb104_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb104_ninner(%rsp),%ecx
        movl  %ecx,nb104_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb104_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel104_x86_64_sse.nb104_unroll_loop
        jmp   _nb_kernel104_x86_64_sse.nb104_single_check
_nb_kernel104_x86_64_sse.nb104_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb104_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb104_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb104_pos(%rbp),%rsi       ## base of pos[] 

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


    subps nb104_ixH1(%rsp),%xmm0
    subps nb104_iyH1(%rsp),%xmm1
    subps nb104_izH1(%rsp),%xmm2
    subps nb104_ixH2(%rsp),%xmm3
    subps nb104_iyH2(%rsp),%xmm4
    subps nb104_izH2(%rsp),%xmm5
    subps nb104_ixM(%rsp),%xmm6
    subps nb104_iyM(%rsp),%xmm7
    subps nb104_izM(%rsp),%xmm8

        movaps %xmm0,nb104_dxH1H1(%rsp)
        movaps %xmm1,nb104_dyH1H1(%rsp)
        movaps %xmm2,nb104_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH2H1(%rsp)
        movaps %xmm4,nb104_dyH2H1(%rsp)
        movaps %xmm5,nb104_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb104_dxMH1(%rsp)
        movaps %xmm7,nb104_dyMH1(%rsp)
        movaps %xmm8,nb104_dzMH1(%rsp)
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

        movaps  nb104_three(%rsp),%xmm9
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

        movaps  nb104_half(%rsp),%xmm0
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
    mulps  nb104_qqHH(%rsp),%xmm0
    mulps  nb104_qqHH(%rsp),%xmm1
    mulps  nb104_qqMH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb104_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb104_vctot(%rsp)

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

        mulps nb104_dxH1H1(%rsp),%xmm7
        mulps nb104_dyH1H1(%rsp),%xmm8
        mulps nb104_dzH1H1(%rsp),%xmm9
        mulps nb104_dxH2H1(%rsp),%xmm10
        mulps nb104_dyH2H1(%rsp),%xmm11
        mulps nb104_dzH2H1(%rsp),%xmm12
        mulps nb104_dxMH1(%rsp),%xmm13
        mulps nb104_dyMH1(%rsp),%xmm14
        mulps nb104_dzMH1(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb104_fixH1(%rsp),%xmm7
    addps nb104_fiyH1(%rsp),%xmm8
    addps nb104_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb104_fixH2(%rsp),%xmm10
    addps nb104_fiyH2(%rsp),%xmm11
    addps nb104_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb104_fixM(%rsp),%xmm13
    addps nb104_fiyM(%rsp),%xmm14
    addps nb104_fizM(%rsp),%xmm15

    movaps %xmm7,nb104_fixH1(%rsp)
    movaps %xmm8,nb104_fiyH1(%rsp)
    movaps %xmm9,nb104_fizH1(%rsp)
    movaps %xmm10,nb104_fixH2(%rsp)
    movaps %xmm11,nb104_fiyH2(%rsp)
    movaps %xmm12,nb104_fizH2(%rsp)
    movaps %xmm13,nb104_fixM(%rsp)
    movaps %xmm14,nb104_fiyM(%rsp)
    movaps %xmm15,nb104_fizM(%rsp)

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

    subps nb104_ixH1(%rsp),%xmm0
    subps nb104_iyH1(%rsp),%xmm1
    subps nb104_izH1(%rsp),%xmm2
    subps nb104_ixH2(%rsp),%xmm3
    subps nb104_iyH2(%rsp),%xmm4
    subps nb104_izH2(%rsp),%xmm5
    subps nb104_ixM(%rsp),%xmm6
    subps nb104_iyM(%rsp),%xmm7
    subps nb104_izM(%rsp),%xmm8

        movaps %xmm0,nb104_dxH1H2(%rsp)
        movaps %xmm1,nb104_dyH1H2(%rsp)
        movaps %xmm2,nb104_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH2H2(%rsp)
        movaps %xmm4,nb104_dyH2H2(%rsp)
        movaps %xmm5,nb104_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb104_dxMH2(%rsp)
        movaps %xmm7,nb104_dyMH2(%rsp)
        movaps %xmm8,nb104_dzMH2(%rsp)
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

        movaps  nb104_three(%rsp),%xmm9
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

        movaps  nb104_half(%rsp),%xmm0
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
    mulps  nb104_qqHH(%rsp),%xmm0
    mulps  nb104_qqHH(%rsp),%xmm1
    mulps  nb104_qqMH(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb104_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb104_vctot(%rsp)

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

        mulps nb104_dxH1H2(%rsp),%xmm7
        mulps nb104_dyH1H2(%rsp),%xmm8
        mulps nb104_dzH1H2(%rsp),%xmm9
        mulps nb104_dxH2H2(%rsp),%xmm10
        mulps nb104_dyH2H2(%rsp),%xmm11
        mulps nb104_dzH2H2(%rsp),%xmm12
        mulps nb104_dxMH2(%rsp),%xmm13
        mulps nb104_dyMH2(%rsp),%xmm14
        mulps nb104_dzMH2(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb104_fixH1(%rsp),%xmm7
    addps nb104_fiyH1(%rsp),%xmm8
    addps nb104_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb104_fixH2(%rsp),%xmm10
    addps nb104_fiyH2(%rsp),%xmm11
    addps nb104_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb104_fixM(%rsp),%xmm13
    addps nb104_fiyM(%rsp),%xmm14
    addps nb104_fizM(%rsp),%xmm15

    movaps %xmm7,nb104_fixH1(%rsp)
    movaps %xmm8,nb104_fiyH1(%rsp)
    movaps %xmm9,nb104_fizH1(%rsp)
    movaps %xmm10,nb104_fixH2(%rsp)
    movaps %xmm11,nb104_fiyH2(%rsp)
    movaps %xmm12,nb104_fizH2(%rsp)
    movaps %xmm13,nb104_fixM(%rsp)
    movaps %xmm14,nb104_fiyM(%rsp)
    movaps %xmm15,nb104_fizM(%rsp)

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

    subps nb104_ixH1(%rsp),%xmm0
    subps nb104_iyH1(%rsp),%xmm1
    subps nb104_izH1(%rsp),%xmm2
    subps nb104_ixH2(%rsp),%xmm3
    subps nb104_iyH2(%rsp),%xmm4
    subps nb104_izH2(%rsp),%xmm5
    subps nb104_ixM(%rsp),%xmm6
    subps nb104_iyM(%rsp),%xmm7
    subps nb104_izM(%rsp),%xmm8

        movaps %xmm0,nb104_dxH1M(%rsp)
        movaps %xmm1,nb104_dyH1M(%rsp)
        movaps %xmm2,nb104_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb104_dxH2M(%rsp)
        movaps %xmm4,nb104_dyH2M(%rsp)
        movaps %xmm5,nb104_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb104_dxMM(%rsp)
        movaps %xmm7,nb104_dyMM(%rsp)
        movaps %xmm8,nb104_dzMM(%rsp)
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

        movaps  nb104_three(%rsp),%xmm9
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

        movaps  nb104_half(%rsp),%xmm0
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
    mulps  nb104_qqMH(%rsp),%xmm0
    mulps  nb104_qqMH(%rsp),%xmm1
    mulps  nb104_qqMM(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb104_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb104_vctot(%rsp)

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

        mulps nb104_dxH1M(%rsp),%xmm7
        mulps nb104_dyH1M(%rsp),%xmm8
        mulps nb104_dzH1M(%rsp),%xmm9
        mulps nb104_dxH2M(%rsp),%xmm10
        mulps nb104_dyH2M(%rsp),%xmm11
        mulps nb104_dzH2M(%rsp),%xmm12
        mulps nb104_dxMM(%rsp),%xmm13
        mulps nb104_dyMM(%rsp),%xmm14
        mulps nb104_dzMM(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb104_fixH1(%rsp),%xmm7
    addps nb104_fiyH1(%rsp),%xmm8
    addps nb104_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb104_fixH2(%rsp),%xmm10
    addps nb104_fiyH2(%rsp),%xmm11
    addps nb104_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb104_fixM(%rsp),%xmm13
    addps nb104_fiyM(%rsp),%xmm14
    addps nb104_fizM(%rsp),%xmm15

    movaps %xmm7,nb104_fixH1(%rsp)
    movaps %xmm8,nb104_fiyH1(%rsp)
    movaps %xmm9,nb104_fizH1(%rsp)
    movaps %xmm10,nb104_fixH2(%rsp)
    movaps %xmm11,nb104_fiyH2(%rsp)
    movaps %xmm12,nb104_fizH2(%rsp)
    movaps %xmm13,nb104_fixM(%rsp)
    movaps %xmm14,nb104_fiyM(%rsp)
    movaps %xmm15,nb104_fizM(%rsp)

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
        subl $4,nb104_innerk(%rsp)
        jl    _nb_kernel104_x86_64_sse.nb104_single_check
        jmp   _nb_kernel104_x86_64_sse.nb104_unroll_loop
_nb_kernel104_x86_64_sse.nb104_single_check: 
        addl $4,nb104_innerk(%rsp)
        jnz   _nb_kernel104_x86_64_sse.nb104_single_loop
        jmp   _nb_kernel104_x86_64_sse.nb104_updateouterdata
_nb_kernel104_x86_64_sse.nb104_single_loop: 
        movq  nb104_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb104_innerjjnr(%rsp)

        movq nb104_pos(%rbp),%rsi
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
        movaps %xmm0,nb104_jxM(%rsp)
        movaps %xmm1,nb104_jyM(%rsp)
        movaps %xmm2,nb104_jzM(%rsp)
        subps  nb104_ixM(%rsp),%xmm0
        subps  nb104_iyM(%rsp),%xmm1
        subps  nb104_izM(%rsp),%xmm2
        movaps %xmm0,nb104_dxMM(%rsp)
        movaps %xmm1,nb104_dyMM(%rsp)
        movaps %xmm2,nb104_dzMM(%rsp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104_half(%rsp),%xmm3   ## rinv iM- j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq

        ## fetch charges to xmm4
        movss   nb104_qqMM(%rsp),%xmm4
        movhps  nb104_qqMH(%rsp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        mulps   %xmm3,%xmm0     ## total fscal 
        addps   nb104_vctot(%rsp),%xmm3
        movaps  %xmm3,nb104_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb104_dxMM(%rsp),%xmm0
        mulps   nb104_dyMM(%rsp),%xmm1
        mulps   nb104_dzMM(%rsp),%xmm2
        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb104_fjxM(%rsp)
        movaps  %xmm4,nb104_fjyM(%rsp)
        movaps  %xmm5,nb104_fjzM(%rsp)
        addps   nb104_fixM(%rsp),%xmm0
        addps   nb104_fiyM(%rsp),%xmm1
        addps   nb104_fizM(%rsp),%xmm2
        movaps  %xmm0,nb104_fixM(%rsp)
        movaps  %xmm1,nb104_fiyM(%rsp)
        movaps  %xmm2,nb104_fizM(%rsp)


        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb104_jxM(%rsp),%xmm0
    movaps  nb104_jyM(%rsp),%xmm1
    movaps  nb104_jzM(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5

        subps   nb104_ixH1(%rsp),%xmm0
        subps   nb104_iyH1(%rsp),%xmm1
        subps   nb104_izH1(%rsp),%xmm2
        subps   nb104_ixH2(%rsp),%xmm3
        subps   nb104_iyH2(%rsp),%xmm4
        subps   nb104_izH2(%rsp),%xmm5

        movaps %xmm0,nb104_dxH1M(%rsp)
        movaps %xmm1,nb104_dyH1M(%rsp)
        movaps %xmm2,nb104_dzH1M(%rsp)
        movaps %xmm3,nb104_dxH2M(%rsp)
        movaps %xmm4,nb104_dyH2M(%rsp)
        movaps %xmm5,nb104_dzH2M(%rsp)
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
        movaps  nb104_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb104_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb104_qqMH(%rsp),%xmm6
        movhps  nb104_qqHH(%rsp),%xmm6

        ## do coulomb interaction 
        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm4
        mulps   %xmm0,%xmm0     ## rinvsq 
        mulps   %xmm4,%xmm4     ## rinvsq 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        movaps  %xmm3,%xmm2
        addps   %xmm7,%xmm2     ## total vcoul 
        mulps   %xmm3,%xmm0     ## fscal 

        addps   nb104_vctot(%rsp),%xmm2
        mulps   %xmm4,%xmm7     ## fscal 
        movaps  %xmm2,nb104_vctot(%rsp)
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb104_dxH1M(%rsp),%xmm0
        mulps   nb104_dyH1M(%rsp),%xmm1
        mulps   nb104_dzH1M(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb104_fjxM(%rsp),%xmm3
        movaps  nb104_fjyM(%rsp),%xmm4
        movaps  nb104_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb104_fjxM(%rsp)
        movaps  %xmm4,nb104_fjyM(%rsp)
        movaps  %xmm5,nb104_fjzM(%rsp)
        addps   nb104_fixH1(%rsp),%xmm0
        addps   nb104_fiyH1(%rsp),%xmm1
        addps   nb104_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb104_fixH1(%rsp)
        movaps  %xmm1,nb104_fiyH1(%rsp)
        movaps  %xmm2,nb104_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb104_dxH2M(%rsp),%xmm0
        mulps   nb104_dyH2M(%rsp),%xmm1
        mulps   nb104_dzH2M(%rsp),%xmm2
        movaps  nb104_fjxM(%rsp),%xmm3
        movaps  nb104_fjyM(%rsp),%xmm4
        movaps  nb104_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb104_faction(%rbp),%rsi
        movaps  %xmm3,nb104_fjxM(%rsp)
        movaps  %xmm4,nb104_fjyM(%rsp)
        movaps  %xmm5,nb104_fjzM(%rsp)
        addps   nb104_fixH2(%rsp),%xmm0
        addps   nb104_fiyH2(%rsp),%xmm1
        addps   nb104_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb104_fixH2(%rsp)
        movaps  %xmm1,nb104_fiyH2(%rsp)
        movaps  %xmm2,nb104_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  36(%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb104_fjxM(%rsp),%xmm3
        movaps  nb104_fjyM(%rsp),%xmm4
        movaps  nb104_fjzM(%rsp),%xmm5
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

        decl  nb104_innerk(%rsp)
        jz    _nb_kernel104_x86_64_sse.nb104_updateouterdata
        jmp   _nb_kernel104_x86_64_sse.nb104_single_loop
_nb_kernel104_x86_64_sse.nb104_updateouterdata: 
        movl  nb104_ii3(%rsp),%ecx
        movq  nb104_faction(%rbp),%rdi
        movq  nb104_fshift(%rbp),%rsi
        movl  nb104_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb104_fixH1(%rsp),%xmm0
        movaps nb104_fiyH1(%rsp),%xmm1
        movaps nb104_fizH1(%rsp),%xmm2

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

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb104_fixH2(%rsp),%xmm0
        movaps nb104_fiyH2(%rsp),%xmm1
        movaps nb104_fizH2(%rsp),%xmm2

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
        movaps nb104_fixM(%rsp),%xmm0
        movaps nb104_fiyM(%rsp),%xmm1
        movaps nb104_fizM(%rsp),%xmm2

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
        movl nb104_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb104_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb104_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb104_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb104_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104_x86_64_sse.nb104_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104_n(%rsp)
        jmp _nb_kernel104_x86_64_sse.nb104_outer
_nb_kernel104_x86_64_sse.nb104_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104_x86_64_sse.nb104_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104_x86_64_sse.nb104_threadloop
_nb_kernel104_x86_64_sse.nb104_end: 


        movl nb104_nouter(%rsp),%eax
        movl nb104_ninner(%rsp),%ebx
        movq nb104_outeriter(%rbp),%rcx
        movq nb104_inneriter(%rbp),%rdx
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







.globl nb_kernel104nf_x86_64_sse
.globl _nb_kernel104nf_x86_64_sse
nb_kernel104nf_x86_64_sse:      
_nb_kernel104nf_x86_64_sse:     
.set nb104nf_fshift, 16
.set nb104nf_gid, 24
.set nb104nf_pos, 32
.set nb104nf_faction, 40
.set nb104nf_charge, 48
.set nb104nf_p_facel, 56
.set nb104nf_argkrf, 64
.set nb104nf_argcrf, 72
.set nb104nf_Vc, 80
.set nb104nf_type, 88
.set nb104nf_p_ntype, 96
.set nb104nf_vdwparam, 104
.set nb104nf_Vvdw, 112
.set nb104nf_p_tabscale, 120
.set nb104nf_VFtab, 128
.set nb104nf_invsqrta, 136
.set nb104nf_dvda, 144
.set nb104nf_p_gbtabscale, 152
.set nb104nf_GBtab, 160
.set nb104nf_p_nthreads, 168
.set nb104nf_count, 176
.set nb104nf_mtx, 184
.set nb104nf_outeriter, 192
.set nb104nf_inneriter, 200
.set nb104nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use         
.set nb104nf_ixH1, 0
.set nb104nf_iyH1, 16
.set nb104nf_izH1, 32
.set nb104nf_ixH2, 48
.set nb104nf_iyH2, 64
.set nb104nf_izH2, 80
.set nb104nf_ixM, 96
.set nb104nf_iyM, 112
.set nb104nf_izM, 128
.set nb104nf_jxH1, 144
.set nb104nf_jyH1, 160
.set nb104nf_jzH1, 176
.set nb104nf_jxH2, 192
.set nb104nf_jyH2, 208
.set nb104nf_jzH2, 224
.set nb104nf_jxM, 240
.set nb104nf_jyM, 256
.set nb104nf_jzM, 272
.set nb104nf_dxMM, 288
.set nb104nf_dyMM, 304
.set nb104nf_dzMM, 320
.set nb104nf_qqHH, 336
.set nb104nf_qqMH, 352
.set nb104nf_qqMM, 368
.set nb104nf_vctot, 384
.set nb104nf_half, 400
.set nb104nf_three, 416
.set nb104nf_rsqH1H1, 432
.set nb104nf_rsqH1H2, 448
.set nb104nf_rsqH1M, 464
.set nb104nf_rsqH2H1, 480
.set nb104nf_rsqH2H2, 496
.set nb104nf_rsqH2M, 512
.set nb104nf_rsqMH1, 528
.set nb104nf_rsqMH2, 544
.set nb104nf_rsqMM, 560
.set nb104nf_rinvH1H1, 576
.set nb104nf_rinvH1H2, 592
.set nb104nf_rinvH1M, 608
.set nb104nf_rinvH2H1, 624
.set nb104nf_rinvH2H2, 640
.set nb104nf_rinvH2M, 656
.set nb104nf_rinvMH1, 672
.set nb104nf_rinvMH2, 688
.set nb104nf_rinvMM, 704
.set nb104nf_is3, 720
.set nb104nf_ii3, 724
.set nb104nf_nri, 728
.set nb104nf_iinr, 736
.set nb104nf_jindex, 744
.set nb104nf_jjnr, 752
.set nb104nf_shift, 760
.set nb104nf_shiftvec, 768
.set nb104nf_facel, 776
.set nb104nf_innerjjnr, 784
.set nb104nf_innerk, 792
.set nb104nf_n, 800
.set nb104nf_nn1, 804
.set nb104nf_nouter, 808
.set nb104nf_ninner, 812

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms
        subq $824,%rsp
        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb104nf_nouter(%rsp)
        movl %eax,nb104nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb104nf_nri(%rsp)
        movq %rsi,nb104nf_iinr(%rsp)
        movq %rdx,nb104nf_jindex(%rsp)
        movq %rcx,nb104nf_jjnr(%rsp)
        movq %r8,nb104nf_shift(%rsp)
        movq %r9,nb104nf_shiftvec(%rsp)
        movq nb104nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb104nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb104nf_half(%rsp)
        movss nb104nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb104nf_half(%rsp)
        movaps %xmm3,nb104nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb104nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb104nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb104nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb104nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104nf_qqHH(%rsp)
        movaps %xmm4,nb104nf_qqMH(%rsp)
        movaps %xmm5,nb104nf_qqMM(%rsp)

_nb_kernel104nf_x86_64_sse.nb104nf_threadloop: 
        movq  nb104nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel104nf_x86_64_sse.nb104nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104nf_x86_64_sse.nb104nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104nf_n(%rsp)
        movl %ebx,nb104nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104nf_x86_64_sse.nb104nf_outerstart
        jmp _nb_kernel104nf_x86_64_sse.nb104nf_end

_nb_kernel104nf_x86_64_sse.nb104nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104nf_nouter(%rsp),%ebx
        movl %ebx,nb104nf_nouter(%rsp)

_nb_kernel104nf_x86_64_sse.nb104nf_outer: 
        movq  nb104nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb104nf_is3(%rsp)            ## store is3 

        movq  nb104nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb104nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb104nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb104nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb104nf_ixH1(%rsp)
        movaps %xmm4,nb104nf_iyH1(%rsp)
        movaps %xmm5,nb104nf_izH1(%rsp)

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
        movaps %xmm0,nb104nf_ixH2(%rsp)
        movaps %xmm1,nb104nf_iyH2(%rsp)
        movaps %xmm2,nb104nf_izH2(%rsp)
        movaps %xmm3,nb104nf_ixM(%rsp)
        movaps %xmm4,nb104nf_iyM(%rsp)
        movaps %xmm5,nb104nf_izM(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb104nf_vctot(%rsp)

        movq  nb104nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb104nf_pos(%rbp),%rsi
        movq  nb104nf_faction(%rbp),%rdi
        movq  nb104nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb104nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb104nf_ninner(%rsp),%ecx
        movl  %ecx,nb104nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb104nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel104nf_x86_64_sse.nb104nf_unroll_loop
        jmp   _nb_kernel104nf_x86_64_sse.nb104nf_single_check
_nb_kernel104nf_x86_64_sse.nb104nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb104nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb104nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb104nf_pos(%rbp),%rsi     ## base of pos[] 

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

        ## current state:       
        ## xmm2= jxh1a  jyH1a  jxH1c  jyH1c 
        ## xmm3= jxH2a jyH2a jxH2c jyH2c 
        ## xmm4= jxMa jyMa jxMc jyMc 
        ## xmm5= jxH1b  jyH1b  jxH1d  jyH1d 
        ## xmm6= jxH2b jyH2b jxH2d jyH2d 
        ## xmm7= jxMb jyMb jxMd jyMd 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxH1a  jxH1b  jyH1a  jyH1b 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH2a jxH2b jyH2a jyH2b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxH1c  jxH1d  jyH1c  jyH1d 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH2c jxH2d jyH2c jyH2d  
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxMa jxMb jyMa jyMb            
        unpckhps %xmm7,%xmm5    ## xmm5= jxMc jxMd jyMc jyMd     
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxH1a  jxH1b  jxH1c  jxH1d  
        movaps %xmm0,nb104nf_jxH1(%rsp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyH1a  jyH1b  jyH1c  jyH1d 
        movaps %xmm2,nb104nf_jyH1(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb104nf_jxH2(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb104nf_jyH2(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb104nf_jxM(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb104nf_jyM(%rsp)

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
        movaps %xmm0,nb104nf_jzH1(%rsp)
        movaps %xmm1,nb104nf_jzH2(%rsp)
        movaps %xmm2,nb104nf_jzM(%rsp)

        movaps nb104nf_ixH1(%rsp),%xmm0
        movaps nb104nf_iyH1(%rsp),%xmm1
        movaps nb104nf_izH1(%rsp),%xmm2
        movaps nb104nf_ixH1(%rsp),%xmm3
        movaps nb104nf_iyH1(%rsp),%xmm4
        movaps nb104nf_izH1(%rsp),%xmm5
        subps  nb104nf_jxH1(%rsp),%xmm0
        subps  nb104nf_jyH1(%rsp),%xmm1
        subps  nb104nf_jzH1(%rsp),%xmm2
        subps  nb104nf_jxH2(%rsp),%xmm3
        subps  nb104nf_jyH2(%rsp),%xmm4
        subps  nb104nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH1H1(%rsp)
        movaps %xmm3,nb104nf_rsqH1H2(%rsp)

        movaps nb104nf_ixH1(%rsp),%xmm0
        movaps nb104nf_iyH1(%rsp),%xmm1
        movaps nb104nf_izH1(%rsp),%xmm2
        movaps nb104nf_ixH2(%rsp),%xmm3
        movaps nb104nf_iyH2(%rsp),%xmm4
        movaps nb104nf_izH2(%rsp),%xmm5
        subps  nb104nf_jxM(%rsp),%xmm0
        subps  nb104nf_jyM(%rsp),%xmm1
        subps  nb104nf_jzM(%rsp),%xmm2
        subps  nb104nf_jxH1(%rsp),%xmm3
        subps  nb104nf_jyH1(%rsp),%xmm4
        subps  nb104nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH1M(%rsp)
        movaps %xmm3,nb104nf_rsqH2H1(%rsp)

        movaps nb104nf_ixH2(%rsp),%xmm0
        movaps nb104nf_iyH2(%rsp),%xmm1
        movaps nb104nf_izH2(%rsp),%xmm2
        movaps nb104nf_ixH2(%rsp),%xmm3
        movaps nb104nf_iyH2(%rsp),%xmm4
        movaps nb104nf_izH2(%rsp),%xmm5
        subps  nb104nf_jxH2(%rsp),%xmm0
        subps  nb104nf_jyH2(%rsp),%xmm1
        subps  nb104nf_jzH2(%rsp),%xmm2
        subps  nb104nf_jxM(%rsp),%xmm3
        subps  nb104nf_jyM(%rsp),%xmm4
        subps  nb104nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb104nf_rsqH2H2(%rsp)
        movaps %xmm3,nb104nf_rsqH2M(%rsp)

        movaps nb104nf_ixM(%rsp),%xmm0
        movaps nb104nf_iyM(%rsp),%xmm1
        movaps nb104nf_izM(%rsp),%xmm2
        movaps nb104nf_ixM(%rsp),%xmm3
        movaps nb104nf_iyM(%rsp),%xmm4
        movaps nb104nf_izM(%rsp),%xmm5
        subps  nb104nf_jxH1(%rsp),%xmm0
        subps  nb104nf_jyH1(%rsp),%xmm1
        subps  nb104nf_jzH1(%rsp),%xmm2
        subps  nb104nf_jxH2(%rsp),%xmm3
        subps  nb104nf_jyH2(%rsp),%xmm4
        subps  nb104nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb104nf_rsqMH1(%rsp)
        movaps %xmm4,nb104nf_rsqMH2(%rsp)

        movaps nb104nf_ixM(%rsp),%xmm0
        movaps nb104nf_iyM(%rsp),%xmm1
        movaps nb104nf_izM(%rsp),%xmm2
        subps  nb104nf_jxM(%rsp),%xmm0
        subps  nb104nf_jyM(%rsp),%xmm1
        subps  nb104nf_jzM(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb104nf_rsqMM(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%rsp),%xmm3   ## rinvMM 
        mulps   nb104nf_half(%rsp),%xmm7   ## rinvMH2 
        movaps  %xmm3,nb104nf_rinvMM(%rsp)
        movaps  %xmm7,nb104nf_rinvMH2(%rsp)

        rsqrtps nb104nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb104nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH1H1(%rsp),%xmm1
        mulps   nb104nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%rsp),%xmm3
        mulps   nb104nf_half(%rsp),%xmm7
        movaps  %xmm3,nb104nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb104nf_rinvH1H2(%rsp)

        rsqrtps nb104nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb104nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH1M(%rsp),%xmm1
        mulps   nb104nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%rsp),%xmm3
        mulps   nb104nf_half(%rsp),%xmm7
        movaps  %xmm3,nb104nf_rinvH1M(%rsp)
        movaps  %xmm7,nb104nf_rinvH2H1(%rsp)

        rsqrtps nb104nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb104nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb104nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb104nf_rsqH2H2(%rsp),%xmm1
        mulps   nb104nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%rsp),%xmm3
        mulps   nb104nf_half(%rsp),%xmm7
        movaps  %xmm3,nb104nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb104nf_rinvH2M(%rsp)

        rsqrtps nb104nf_rsqMH1(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb104nf_three(%rsp),%xmm3
        mulps   nb104nf_rsqMH1(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104nf_half(%rsp),%xmm3
        movaps  %xmm3,nb104nf_rinvMH1(%rsp)

        ## all H-H interactions
        movaps nb104nf_rinvH1H1(%rsp),%xmm0
        addps  nb104nf_rinvH1H2(%rsp),%xmm0
        addps  nb104nf_rinvH2H1(%rsp),%xmm0
        addps  nb104nf_rinvH2H2(%rsp),%xmm0
        mulps  nb104nf_qqHH(%rsp),%xmm0
        ## all M-H interactions
        movaps nb104nf_rinvH1M(%rsp),%xmm1
        addps  nb104nf_rinvH2M(%rsp),%xmm1
        addps  nb104nf_rinvMH1(%rsp),%xmm1
        addps  nb104nf_rinvMH2(%rsp),%xmm1
        mulps  nb104nf_qqMH(%rsp),%xmm1
        ## The M-M interaction
        movaps nb104nf_rinvMM(%rsp),%xmm2
        mulps  nb104nf_qqMM(%rsp),%xmm2
        addps  %xmm1,%xmm0
        addps  nb104nf_vctot(%rsp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb104nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb104nf_innerk(%rsp)
        jl    _nb_kernel104nf_x86_64_sse.nb104nf_single_check
        jmp   _nb_kernel104nf_x86_64_sse.nb104nf_unroll_loop
_nb_kernel104nf_x86_64_sse.nb104nf_single_check: 
        addl $4,nb104nf_innerk(%rsp)
        jnz   _nb_kernel104nf_x86_64_sse.nb104nf_single_loop
        jmp   _nb_kernel104nf_x86_64_sse.nb104nf_updateouterdata
_nb_kernel104nf_x86_64_sse.nb104nf_single_loop: 
        movq  nb104nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb104nf_innerjjnr(%rsp)

        movq nb104nf_pos(%rbp),%rsi
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
        movaps  nb104nf_ixM(%rsp),%xmm0
        movaps  nb104nf_iyM(%rsp),%xmm1
        movaps  nb104nf_izM(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM 
        movaps %xmm3,nb104nf_jxM(%rsp)
        movaps %xmm4,nb104nf_jyM(%rsp)
        movaps %xmm5,nb104nf_jzM(%rsp)
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
        movaps  nb104nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb104nf_half(%rsp),%xmm3   ## rinv iM- j water 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq

        ## fetch charges to xmm4
        movss   nb104nf_qqMM(%rsp),%xmm4
        movhps  nb104nf_qqMH(%rsp),%xmm4

        mulps   %xmm4,%xmm3     ## xmm3=vcoul 
        addps   nb104nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb104nf_vctot(%rsp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb104nf_ixH1(%rsp),%xmm0
        movaps  nb104nf_iyH1(%rsp),%xmm1
        movaps  nb104nf_izH1(%rsp),%xmm2
        movaps  nb104nf_ixH2(%rsp),%xmm3
        movaps  nb104nf_iyH2(%rsp),%xmm4
        movaps  nb104nf_izH2(%rsp),%xmm5
        subps   nb104nf_jxM(%rsp),%xmm0
        subps   nb104nf_jyM(%rsp),%xmm1
        subps   nb104nf_jzM(%rsp),%xmm2
        subps   nb104nf_jxM(%rsp),%xmm3
        subps   nb104nf_jyM(%rsp),%xmm4
        subps   nb104nf_jzM(%rsp),%xmm5
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
        movaps  nb104nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb104nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb104nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb104nf_qqMH(%rsp),%xmm6
        movhps  nb104nf_qqHH(%rsp),%xmm6

        ## do coulomb interaction 
        mulps   %xmm6,%xmm3     ## vcoul 
        mulps   %xmm6,%xmm7     ## vcoul 
        addps   %xmm7,%xmm3     ## total vcoul 
        addps   nb104nf_vctot(%rsp),%xmm3
        movaps  %xmm3,nb104nf_vctot(%rsp)

        decl  nb104nf_innerk(%rsp)
        jz    _nb_kernel104nf_x86_64_sse.nb104nf_updateouterdata
        jmp   _nb_kernel104nf_x86_64_sse.nb104nf_single_loop
_nb_kernel104nf_x86_64_sse.nb104nf_updateouterdata: 
        ## get n from stack
        movl nb104nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb104nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb104nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb104nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb104nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104nf_x86_64_sse.nb104nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104nf_n(%rsp)
        jmp _nb_kernel104nf_x86_64_sse.nb104nf_outer
_nb_kernel104nf_x86_64_sse.nb104nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104nf_x86_64_sse.nb104nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104nf_x86_64_sse.nb104nf_threadloop
_nb_kernel104nf_x86_64_sse.nb104nf_end: 


        movl nb104nf_nouter(%rsp),%eax
        movl nb104nf_ninner(%rsp),%ebx
        movq nb104nf_outeriter(%rbp),%rcx
        movq nb104nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $824,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret





