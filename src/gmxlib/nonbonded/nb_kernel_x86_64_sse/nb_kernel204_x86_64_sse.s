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

		
.globl nb_kernel204_x86_64_sse
.globl _nb_kernel204_x86_64_sse
nb_kernel204_x86_64_sse:        
_nb_kernel204_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb204_fshift, 16
.set nb204_gid, 24
.set nb204_pos, 32
.set nb204_faction, 40
.set nb204_charge, 48
.set nb204_p_facel, 56
.set nb204_argkrf, 64
.set nb204_argcrf, 72
.set nb204_Vc, 80
.set nb204_type, 88
.set nb204_p_ntype, 96
.set nb204_vdwparam, 104
.set nb204_Vvdw, 112
.set nb204_p_tabscale, 120
.set nb204_VFtab, 128
.set nb204_invsqrta, 136
.set nb204_dvda, 144
.set nb204_p_gbtabscale, 152
.set nb204_GBtab, 160
.set nb204_p_nthreads, 168
.set nb204_count, 176
.set nb204_mtx, 184
.set nb204_outeriter, 192
.set nb204_inneriter, 200
.set nb204_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb204_ixH1, 0
.set nb204_iyH1, 16
.set nb204_izH1, 32
.set nb204_ixH2, 48
.set nb204_iyH2, 64
.set nb204_izH2, 80
.set nb204_ixM, 96
.set nb204_iyM, 112
.set nb204_izM, 128
.set nb204_jxH1, 144
.set nb204_jyH1, 160
.set nb204_jzH1, 176
.set nb204_jxH2, 192
.set nb204_jyH2, 208
.set nb204_jzH2, 224
.set nb204_jxM, 240
.set nb204_jyM, 256
.set nb204_jzM, 272
.set nb204_dxH1H1, 288
.set nb204_dyH1H1, 304
.set nb204_dzH1H1, 320
.set nb204_dxH1H2, 336
.set nb204_dyH1H2, 352
.set nb204_dzH1H2, 368
.set nb204_dxH1M, 384
.set nb204_dyH1M, 400
.set nb204_dzH1M, 416
.set nb204_dxH2H1, 432
.set nb204_dyH2H1, 448
.set nb204_dzH2H1, 464
.set nb204_dxH2H2, 480
.set nb204_dyH2H2, 496
.set nb204_dzH2H2, 512
.set nb204_dxH2M, 528
.set nb204_dyH2M, 544
.set nb204_dzH2M, 560
.set nb204_dxMH1, 576
.set nb204_dyMH1, 592
.set nb204_dzMH1, 608
.set nb204_dxMH2, 624
.set nb204_dyMH2, 640
.set nb204_dzMH2, 656
.set nb204_dxMM, 672
.set nb204_dyMM, 688
.set nb204_dzMM, 704
.set nb204_qqHH, 720
.set nb204_qqMH, 736
.set nb204_qqMM, 752
.set nb204_vctot, 768
.set nb204_fixH1, 784
.set nb204_fiyH1, 800
.set nb204_fizH1, 816
.set nb204_fixH2, 832
.set nb204_fiyH2, 848
.set nb204_fizH2, 864
.set nb204_fixM, 880
.set nb204_fiyM, 896
.set nb204_fizM, 912
.set nb204_fjxH1, 928
.set nb204_fjyH1, 944
.set nb204_fjzH1, 960
.set nb204_fjxH2, 976
.set nb204_fjyH2, 992
.set nb204_fjzH2, 1008
.set nb204_fjxM, 1024
.set nb204_fjyM, 1040
.set nb204_fjzM, 1056
.set nb204_half, 1072
.set nb204_three, 1088
.set nb204_rsqH1H1, 1104
.set nb204_rsqH1H2, 1120
.set nb204_rsqH1M, 1136
.set nb204_rsqH2H1, 1152
.set nb204_rsqH2H2, 1168
.set nb204_rsqH2M, 1184
.set nb204_rsqMH1, 1200
.set nb204_rsqMH2, 1216
.set nb204_rsqMM, 1232
.set nb204_rinvH1H1, 1248
.set nb204_rinvH1H2, 1264
.set nb204_rinvH1M, 1280
.set nb204_rinvH2H1, 1296
.set nb204_rinvH2H2, 1312
.set nb204_rinvH2M, 1328
.set nb204_rinvMH1, 1344
.set nb204_rinvMH2, 1360
.set nb204_rinvMM, 1376
.set nb204_two, 1392
.set nb204_krf, 1408
.set nb204_crf, 1424
.set nb204_is3, 1440
.set nb204_ii3, 1444
.set nb204_innerjjnr, 1448
.set nb204_nri, 1456
.set nb204_iinr, 1464
.set nb204_jindex, 1472
.set nb204_jjnr, 1480
.set nb204_shift, 1488
.set nb204_shiftvec, 1496
.set nb204_facel, 1504
.set nb204_innerk, 1512
.set nb204_n, 1516
.set nb204_nn1, 1520
.set nb204_nouter, 1524
.set nb204_ninner, 1528

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1544,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb204_nouter(%rsp)
        movl %eax,nb204_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb204_nri(%rsp)
        movq %rsi,nb204_iinr(%rsp)
        movq %rdx,nb204_jindex(%rsp)
        movq %rcx,nb204_jjnr(%rsp)
        movq %r8,nb204_shift(%rsp)
        movq %r9,nb204_shiftvec(%rsp)
        movq nb204_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb204_facel(%rsp)


        movq nb204_argkrf(%rbp),%rsi
        movq nb204_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb204_krf(%rsp)
        movaps %xmm2,nb204_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb204_half(%rsp)
        movss nb204_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb204_half(%rsp)
        movaps %xmm2,nb204_two(%rsp)
        movaps %xmm3,nb204_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb204_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb204_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb204_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb204_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204_qqHH(%rsp)
        movaps %xmm4,nb204_qqMH(%rsp)
        movaps %xmm5,nb204_qqMM(%rsp)

_nb_kernel204_x86_64_sse.nb204_threadloop: 
        movq  nb204_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel204_x86_64_sse.nb204_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204_x86_64_sse.nb204_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204_n(%rsp)
        movl %ebx,nb204_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204_x86_64_sse.nb204_outerstart
        jmp _nb_kernel204_x86_64_sse.nb204_end

_nb_kernel204_x86_64_sse.nb204_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204_nouter(%rsp),%ebx
        movl %ebx,nb204_nouter(%rsp)

_nb_kernel204_x86_64_sse.nb204_outer: 
        movq  nb204_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb204_is3(%rsp)      ## store is3 

        movq  nb204_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb204_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb204_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb204_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204_ixH1(%rsp)
        movaps %xmm4,nb204_iyH1(%rsp)
        movaps %xmm5,nb204_izH1(%rsp)

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
        movaps %xmm0,nb204_ixH2(%rsp)
        movaps %xmm1,nb204_iyH2(%rsp)
        movaps %xmm2,nb204_izH2(%rsp)
        movaps %xmm3,nb204_ixM(%rsp)
        movaps %xmm4,nb204_iyM(%rsp)
        movaps %xmm5,nb204_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb204_vctot(%rsp)
        movaps %xmm4,nb204_fixH1(%rsp)
        movaps %xmm4,nb204_fiyH1(%rsp)
        movaps %xmm4,nb204_fizH1(%rsp)
        movaps %xmm4,nb204_fixH2(%rsp)
        movaps %xmm4,nb204_fiyH2(%rsp)
        movaps %xmm4,nb204_fizH2(%rsp)
        movaps %xmm4,nb204_fixM(%rsp)
        movaps %xmm4,nb204_fiyM(%rsp)
        movaps %xmm4,nb204_fizM(%rsp)

        movq  nb204_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb204_pos(%rbp),%rsi
        movq  nb204_faction(%rbp),%rdi
        movq  nb204_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb204_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb204_ninner(%rsp),%ecx
        movl  %ecx,nb204_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb204_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel204_x86_64_sse.nb204_unroll_loop
        jmp   _nb_kernel204_x86_64_sse.nb204_single_check
_nb_kernel204_x86_64_sse.nb204_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb204_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb204_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb204_pos(%rbp),%rsi       ## base of pos[] 

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

    subps nb204_ixH1(%rsp),%xmm0
    subps nb204_iyH1(%rsp),%xmm1
    subps nb204_izH1(%rsp),%xmm2
    subps nb204_ixH2(%rsp),%xmm3
    subps nb204_iyH2(%rsp),%xmm4
    subps nb204_izH2(%rsp),%xmm5
    subps nb204_ixM(%rsp),%xmm6
    subps nb204_iyM(%rsp),%xmm7
    subps nb204_izM(%rsp),%xmm8

        movaps %xmm0,nb204_dxH1H1(%rsp)
        movaps %xmm1,nb204_dyH1H1(%rsp)
        movaps %xmm2,nb204_dzH1H1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH2H1(%rsp)
        movaps %xmm4,nb204_dyH2H1(%rsp)
        movaps %xmm5,nb204_dzH2H1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb204_dxMH1(%rsp)
        movaps %xmm7,nb204_dyMH1(%rsp)
        movaps %xmm8,nb204_dzMH1(%rsp)
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

        movaps  nb204_three(%rsp),%xmm9
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

        movaps  nb204_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H1 
        mulps   %xmm4,%xmm10 ## rinvH2H1
    mulps   %xmm4,%xmm11 ## rinvMH1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb204_krf(%rsp),%xmm2
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
    movaps nb204_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb204_qqHH(%rsp),%xmm12
    movaps nb204_qqMH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm12,%xmm4
    mulps  %xmm13,%xmm7
    addps  nb204_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,%xmm15

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

        mulps nb204_dxH1H1(%rsp),%xmm0
        mulps nb204_dyH1H1(%rsp),%xmm1
        mulps nb204_dzH1H1(%rsp),%xmm2
        mulps nb204_dxH2H1(%rsp),%xmm3
        mulps nb204_dyH2H1(%rsp),%xmm4
        mulps nb204_dzH2H1(%rsp),%xmm5
        mulps nb204_dxMH1(%rsp),%xmm6
        mulps nb204_dyMH1(%rsp),%xmm7
        mulps nb204_dzMH1(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb204_fixH1(%rsp),%xmm0
    addps nb204_fiyH1(%rsp),%xmm1
    addps nb204_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb204_fixH2(%rsp),%xmm3
    addps nb204_fiyH2(%rsp),%xmm4
    addps nb204_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb204_fixM(%rsp),%xmm6
    addps nb204_fiyM(%rsp),%xmm7
    addps nb204_fizM(%rsp),%xmm8

    movaps %xmm0,nb204_fixH1(%rsp)
    movaps %xmm1,nb204_fiyH1(%rsp)
    movaps %xmm2,nb204_fizH1(%rsp)
    movaps %xmm3,nb204_fixH2(%rsp)
    movaps %xmm4,nb204_fiyH2(%rsp)
    movaps %xmm5,nb204_fizH2(%rsp)
    movaps %xmm6,nb204_fixM(%rsp)
    movaps %xmm7,nb204_fiyM(%rsp)
    movaps %xmm8,nb204_fizM(%rsp)

    ## xmm9 = fH1x
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


    subps nb204_ixH1(%rsp),%xmm0
    subps nb204_iyH1(%rsp),%xmm1
    subps nb204_izH1(%rsp),%xmm2
    subps nb204_ixH2(%rsp),%xmm3
    subps nb204_iyH2(%rsp),%xmm4
    subps nb204_izH2(%rsp),%xmm5
    subps nb204_ixM(%rsp),%xmm6
    subps nb204_iyM(%rsp),%xmm7
    subps nb204_izM(%rsp),%xmm8

        movaps %xmm0,nb204_dxH1H2(%rsp)
        movaps %xmm1,nb204_dyH1H2(%rsp)
        movaps %xmm2,nb204_dzH1H2(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH2H2(%rsp)
        movaps %xmm4,nb204_dyH2H2(%rsp)
        movaps %xmm5,nb204_dzH2H2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb204_dxMH2(%rsp)
        movaps %xmm7,nb204_dyMH2(%rsp)
        movaps %xmm8,nb204_dzMH2(%rsp)
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

        movaps  nb204_three(%rsp),%xmm9
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

        movaps  nb204_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1H2
        mulps   %xmm4,%xmm10 ## rinvH2H2
    mulps   %xmm4,%xmm11 ## rinvMH2

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb204_krf(%rsp),%xmm2
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
    movaps nb204_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb204_qqHH(%rsp),%xmm12
    movaps nb204_qqMH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm12,%xmm4
    mulps  %xmm13,%xmm7
    addps  %xmm2,%xmm15
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm15

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
    movhps 32(%rdi,%rbx,4),%xmm11    ## jzH2a  -  jzH2b  -
    movhps 32(%rdi,%rdx,4),%xmm12    ## jzH2c  -  jzH2d -

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

        mulps nb204_dxH1H2(%rsp),%xmm0
        mulps nb204_dyH1H2(%rsp),%xmm1
        mulps nb204_dzH1H2(%rsp),%xmm2
        mulps nb204_dxH2H2(%rsp),%xmm3
        mulps nb204_dyH2H2(%rsp),%xmm4
        mulps nb204_dzH2H2(%rsp),%xmm5
        mulps nb204_dxMH2(%rsp),%xmm6
        mulps nb204_dyMH2(%rsp),%xmm7
        mulps nb204_dzMH2(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb204_fixH1(%rsp),%xmm0
    addps nb204_fiyH1(%rsp),%xmm1
    addps nb204_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb204_fixH2(%rsp),%xmm3
    addps nb204_fiyH2(%rsp),%xmm4
    addps nb204_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb204_fixM(%rsp),%xmm6
    addps nb204_fiyM(%rsp),%xmm7
    addps nb204_fizM(%rsp),%xmm8

    movaps %xmm0,nb204_fixH1(%rsp)
    movaps %xmm1,nb204_fiyH1(%rsp)
    movaps %xmm2,nb204_fizH1(%rsp)
    movaps %xmm3,nb204_fixH2(%rsp)
    movaps %xmm4,nb204_fiyH2(%rsp)
    movaps %xmm5,nb204_fizH2(%rsp)
    movaps %xmm6,nb204_fixM(%rsp)
    movaps %xmm7,nb204_fiyM(%rsp)
    movaps %xmm8,nb204_fizM(%rsp)

    ## xmm9  = fH2x
    ## xmm10 = fH2y
    ## xmm11 = fH2z
    movaps %xmm13,%xmm0
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm0

    addps %xmm13,%xmm9
    addps %xmm0,%xmm10

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

    shufps $136,%xmm3,%xmm2 ## 10001000 => jzMa jzMb jzMc jzMd

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb204_ixH1(%rsp),%xmm0
    subps nb204_iyH1(%rsp),%xmm1
    subps nb204_izH1(%rsp),%xmm2
    subps nb204_ixH2(%rsp),%xmm3
    subps nb204_iyH2(%rsp),%xmm4
    subps nb204_izH2(%rsp),%xmm5
    subps nb204_ixM(%rsp),%xmm6
    subps nb204_iyM(%rsp),%xmm7
    subps nb204_izM(%rsp),%xmm8

        movaps %xmm0,nb204_dxH1M(%rsp)
        movaps %xmm1,nb204_dyH1M(%rsp)
        movaps %xmm2,nb204_dzH1M(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH2M(%rsp)
        movaps %xmm4,nb204_dyH2M(%rsp)
        movaps %xmm5,nb204_dzH2M(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb204_dxMM(%rsp)
        movaps %xmm7,nb204_dyMM(%rsp)
        movaps %xmm8,nb204_dzMM(%rsp)
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

        movaps  nb204_three(%rsp),%xmm9
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

        movaps  nb204_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1M
        mulps   %xmm4,%xmm10 ## rinvH2M
    mulps   %xmm4,%xmm11 ## rinvMM

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb204_krf(%rsp),%xmm2
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
    movaps nb204_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb204_qqMH(%rsp),%xmm12
    movaps nb204_qqMM(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm12,%xmm4
    mulps  %xmm13,%xmm7
    addps  %xmm8,%xmm5
    addps  %xmm15,%xmm2
    addps  %xmm5,%xmm2
    movaps %xmm2,nb204_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j M forces to local temp variables 
    movlps 36(%rdi,%rax,4),%xmm9    ## jxMa jyMa  -   -
    movlps 36(%rdi,%rcx,4),%xmm10    ## jxMc jyMc  -   -
    movhps 36(%rdi,%rbx,4),%xmm9    ## jxMa jyMa jxMb jyMb 
    movhps 36(%rdi,%rdx,4),%xmm10    ## jxMc jyMc jxMd jyMd 

    movss  44(%rdi,%rax,4),%xmm11    ## jzMa  -  -  -
    movss  44(%rdi,%rcx,4),%xmm12    ## jzMc  -  -  -
    movss  44(%rdi,%rbx,4),%xmm2     ## jzMb  -  -  -
    movss  44(%rdi,%rdx,4),%xmm3     ## jzMd  -  -  -
    movlhps %xmm2,%xmm11 ## jzMa  -  jzMb  -
    movlhps %xmm3,%xmm12 ## jzMc  -  jzMd -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jzMa jzMb jzMc jzMd

    ## xmm9: jxMa jyMa jxMb jyMb 
    ## xmm10: jxMc jyMc jxMd jyMd
    ## xmm11: jzMa jzMb jzMc jzMd

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb204_dxH1M(%rsp),%xmm0
        mulps nb204_dyH1M(%rsp),%xmm1
        mulps nb204_dzH1M(%rsp),%xmm2
        mulps nb204_dxH2M(%rsp),%xmm3
        mulps nb204_dyH2M(%rsp),%xmm4
        mulps nb204_dzH2M(%rsp),%xmm5
        mulps nb204_dxMM(%rsp),%xmm6
        mulps nb204_dyMM(%rsp),%xmm7
        mulps nb204_dzMM(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb204_fixH1(%rsp),%xmm0
    addps nb204_fiyH1(%rsp),%xmm1
    addps nb204_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb204_fixH2(%rsp),%xmm3
    addps nb204_fiyH2(%rsp),%xmm4
    addps nb204_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb204_fixM(%rsp),%xmm6
    addps nb204_fiyM(%rsp),%xmm7
    addps nb204_fizM(%rsp),%xmm8

    movaps %xmm0,nb204_fixH1(%rsp)
    movaps %xmm1,nb204_fiyH1(%rsp)
    movaps %xmm2,nb204_fizH1(%rsp)
    movaps %xmm3,nb204_fixH2(%rsp)
    movaps %xmm4,nb204_fiyH2(%rsp)
    movaps %xmm5,nb204_fizH2(%rsp)
    movaps %xmm6,nb204_fixM(%rsp)
    movaps %xmm7,nb204_fiyM(%rsp)
    movaps %xmm8,nb204_fizM(%rsp)

    ## xmm0 = fMx
    ## xmm1 = fMy
    ## xmm2 = fMz
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12 ## fMzc fMzd

    movlps %xmm9,36(%rdi,%rax,4)
    movhps %xmm9,36(%rdi,%rbx,4)
    movlps %xmm10,36(%rdi,%rcx,4)
    movhps %xmm10,36(%rdi,%rdx,4)
    movss  %xmm11,44(%rdi,%rax,4)
    movss  %xmm12,44(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,44(%rdi,%rbx,4)
    movss  %xmm12,44(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb204_innerk(%rsp)
        jl    _nb_kernel204_x86_64_sse.nb204_single_check
        jmp   _nb_kernel204_x86_64_sse.nb204_unroll_loop
_nb_kernel204_x86_64_sse.nb204_single_check: 
        addl $4,nb204_innerk(%rsp)
        jnz   _nb_kernel204_x86_64_sse.nb204_single_loop
        jmp   _nb_kernel204_x86_64_sse.nb204_updateouterdata
_nb_kernel204_x86_64_sse.nb204_single_loop: 
        movq  nb204_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb204_innerjjnr(%rsp)

        movq nb204_pos(%rbp),%rsi
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
        movaps %xmm0,nb204_jxM(%rsp)
        movaps %xmm1,nb204_jyM(%rsp)
        movaps %xmm2,nb204_jzM(%rsp)
        subps  nb204_ixM(%rsp),%xmm0
        subps  nb204_iyM(%rsp),%xmm1
        subps  nb204_izM(%rsp),%xmm2
        movaps %xmm0,nb204_dxMM(%rsp)
        movaps %xmm1,nb204_dyMM(%rsp)
        movaps %xmm2,nb204_dzMM(%rsp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb204_krf(%rsp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        movaps  %xmm6,%xmm7     ## xmm7=krsq 
        mulps   %xmm1,%xmm1
        movaps  nb204_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204_half(%rsp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        mulps   nb204_two(%rsp),%xmm7
        subps   nb204_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb204_qqMM(%rsp),%xmm4
        movhps  nb204_qqMH(%rsp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  
        mulps %xmm4,%xmm3       ## coul part of fs  


        addps   nb204_vctot(%rsp),%xmm6
        mulps   %xmm3,%xmm0     ## total fscal 
        movaps  %xmm6,nb204_vctot(%rsp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb204_dxMM(%rsp),%xmm0
        mulps   nb204_dyMM(%rsp),%xmm1
        mulps   nb204_dzMM(%rsp),%xmm2

        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb204_fjxM(%rsp)
        movaps  %xmm4,nb204_fjyM(%rsp)
        movaps  %xmm5,nb204_fjzM(%rsp)
        addps   nb204_fixM(%rsp),%xmm0
        addps   nb204_fiyM(%rsp),%xmm1
        addps   nb204_fizM(%rsp),%xmm2
        movaps  %xmm0,nb204_fixM(%rsp)
        movaps  %xmm1,nb204_fiyM(%rsp)
        movaps  %xmm2,nb204_fizM(%rsp)


        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  nb204_jxM(%rsp),%xmm0
    movaps  nb204_jyM(%rsp),%xmm1
    movaps  nb204_jzM(%rsp),%xmm2
    movaps  %xmm0,%xmm3
    movaps  %xmm1,%xmm4
    movaps  %xmm2,%xmm5

        subps   nb204_ixH1(%rsp),%xmm0
        subps   nb204_iyH1(%rsp),%xmm1
        subps   nb204_izH1(%rsp),%xmm2
        subps   nb204_ixH2(%rsp),%xmm3
        subps   nb204_iyH2(%rsp),%xmm4
        subps   nb204_izH2(%rsp),%xmm5

        movaps %xmm0,nb204_dxH1M(%rsp)
        movaps %xmm1,nb204_dyH1M(%rsp)
        movaps %xmm2,nb204_dzH1M(%rsp)
        movaps %xmm3,nb204_dxH2M(%rsp)
        movaps %xmm4,nb204_dyH2M(%rsp)
        movaps %xmm5,nb204_dzH2M(%rsp)
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
        movaps  nb204_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb204_half(%rsp),%xmm7   ## rinv H2 - j water  

        mulps nb204_krf(%rsp),%xmm0   ## krsq 
        mulps nb204_krf(%rsp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb204_qqMH(%rsp),%xmm6
        movhps  nb204_qqHH(%rsp),%xmm6
        movaps  %xmm0,%xmm1
        movaps  %xmm4,%xmm5
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb204_crf(%rsp),%xmm0
        subps   nb204_crf(%rsp),%xmm4
        mulps   nb204_two(%rsp),%xmm1
        mulps   nb204_two(%rsp),%xmm5
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb204_vctot(%rsp),%xmm4
        movaps  %xmm4,nb204_vctot(%rsp)
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
        mulps   nb204_dxH1M(%rsp),%xmm0
        mulps   nb204_dyH1M(%rsp),%xmm1
        mulps   nb204_dzH1M(%rsp),%xmm2
        ## update forces H1 - j water 
        movaps  nb204_fjxM(%rsp),%xmm3
        movaps  nb204_fjyM(%rsp),%xmm4
        movaps  nb204_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,nb204_fjxM(%rsp)
        movaps  %xmm4,nb204_fjyM(%rsp)
        movaps  %xmm5,nb204_fjzM(%rsp)
        addps   nb204_fixH1(%rsp),%xmm0
        addps   nb204_fiyH1(%rsp),%xmm1
        addps   nb204_fizH1(%rsp),%xmm2
        movaps  %xmm0,nb204_fixH1(%rsp)
        movaps  %xmm1,nb204_fiyH1(%rsp)
        movaps  %xmm2,nb204_fizH1(%rsp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb204_dxH2M(%rsp),%xmm0
        mulps   nb204_dyH2M(%rsp),%xmm1
        mulps   nb204_dzH2M(%rsp),%xmm2
        movaps  nb204_fjxM(%rsp),%xmm3
        movaps  nb204_fjyM(%rsp),%xmm4
        movaps  nb204_fjzM(%rsp),%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movq    nb204_faction(%rbp),%rsi
        movaps  %xmm3,nb204_fjxM(%rsp)
        movaps  %xmm4,nb204_fjyM(%rsp)
        movaps  %xmm5,nb204_fjzM(%rsp)
        addps   nb204_fixH2(%rsp),%xmm0
        addps   nb204_fiyH2(%rsp),%xmm1
        addps   nb204_fizH2(%rsp),%xmm2
        movaps  %xmm0,nb204_fixH2(%rsp)
        movaps  %xmm1,nb204_fiyH2(%rsp)
        movaps  %xmm2,nb204_fizH2(%rsp)

        ## update j water forces from local variables 
        movlps  36(%rsi,%rax,4),%xmm0
        movlps  12(%rsi,%rax,4),%xmm1
        movhps  24(%rsi,%rax,4),%xmm1
        movaps  nb204_fjxM(%rsp),%xmm3
        movaps  nb204_fjyM(%rsp),%xmm4
        movaps  nb204_fjzM(%rsp),%xmm5
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

        decl nb204_innerk(%rsp)
        jz    _nb_kernel204_x86_64_sse.nb204_updateouterdata
        jmp   _nb_kernel204_x86_64_sse.nb204_single_loop
_nb_kernel204_x86_64_sse.nb204_updateouterdata: 
        movl  nb204_ii3(%rsp),%ecx
        movq  nb204_faction(%rbp),%rdi
        movq  nb204_fshift(%rbp),%rsi
        movl  nb204_is3(%rsp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb204_fixH1(%rsp),%xmm0
        movaps nb204_fiyH1(%rsp),%xmm1
        movaps nb204_fizH1(%rsp),%xmm2

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
        movaps nb204_fixH2(%rsp),%xmm0
        movaps nb204_fiyH2(%rsp),%xmm1
        movaps nb204_fizH2(%rsp),%xmm2

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
        movaps nb204_fixM(%rsp),%xmm0
        movaps nb204_fiyM(%rsp),%xmm1
        movaps nb204_fizM(%rsp),%xmm2

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
        movl nb204_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb204_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb204_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb204_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb204_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204_x86_64_sse.nb204_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204_n(%rsp)
        jmp _nb_kernel204_x86_64_sse.nb204_outer
_nb_kernel204_x86_64_sse.nb204_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204_x86_64_sse.nb204_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204_x86_64_sse.nb204_threadloop
_nb_kernel204_x86_64_sse.nb204_end: 

        movl nb204_nouter(%rsp),%eax
        movl nb204_ninner(%rsp),%ebx
        movq nb204_outeriter(%rbp),%rcx
        movq nb204_inneriter(%rbp),%rdx
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








.globl nb_kernel204nf_x86_64_sse
.globl _nb_kernel204nf_x86_64_sse
nb_kernel204nf_x86_64_sse:      
_nb_kernel204nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb204nf_fshift, 16
.set nb204nf_gid, 24
.set nb204nf_pos, 32
.set nb204nf_faction, 40
.set nb204nf_charge, 48
.set nb204nf_p_facel, 56
.set nb204nf_argkrf, 64
.set nb204nf_argcrf, 72
.set nb204nf_Vc, 80
.set nb204nf_type, 88
.set nb204nf_p_ntype, 96
.set nb204nf_vdwparam, 104
.set nb204nf_Vvdw, 112
.set nb204nf_p_tabscale, 120
.set nb204nf_VFtab, 128
.set nb204nf_invsqrta, 136
.set nb204nf_dvda, 144
.set nb204nf_p_gbtabscale, 152
.set nb204nf_GBtab, 160
.set nb204nf_p_nthreads, 168
.set nb204nf_count, 176
.set nb204nf_mtx, 184
.set nb204nf_outeriter, 192
.set nb204nf_inneriter, 200
.set nb204nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb204nf_ixH1, 0
.set nb204nf_iyH1, 16
.set nb204nf_izH1, 32
.set nb204nf_ixH2, 48
.set nb204nf_iyH2, 64
.set nb204nf_izH2, 80
.set nb204nf_ixM, 96
.set nb204nf_iyM, 112
.set nb204nf_izM, 128
.set nb204nf_jxH1, 144
.set nb204nf_jyH1, 160
.set nb204nf_jzH1, 176
.set nb204nf_jxH2, 192
.set nb204nf_jyH2, 208
.set nb204nf_jzH2, 224
.set nb204nf_jxM, 240
.set nb204nf_jyM, 256
.set nb204nf_jzM, 272
.set nb204nf_qqHH, 288
.set nb204nf_qqMH, 304
.set nb204nf_qqMM, 320
.set nb204nf_vctot, 336
.set nb204nf_half, 352
.set nb204nf_three, 368
.set nb204nf_rsqH1H1, 384
.set nb204nf_rsqH1H2, 400
.set nb204nf_rsqH1M, 416
.set nb204nf_rsqH2H1, 432
.set nb204nf_rsqH2H2, 448
.set nb204nf_rsqH2M, 464
.set nb204nf_rsqMH1, 480
.set nb204nf_rsqMH2, 496
.set nb204nf_rsqMM, 512
.set nb204nf_rinvH1H1, 528
.set nb204nf_rinvH1H2, 544
.set nb204nf_rinvH1M, 560
.set nb204nf_rinvH2H1, 576
.set nb204nf_rinvH2H2, 592
.set nb204nf_rinvH2M, 608
.set nb204nf_rinvMH1, 624
.set nb204nf_rinvMH2, 640
.set nb204nf_rinvMM, 656
.set nb204nf_krf, 672
.set nb204nf_crf, 688
.set nb204nf_is3, 704
.set nb204nf_ii3, 708
.set nb204nf_innerjjnr, 712
.set nb204nf_nri, 720
.set nb204nf_iinr, 728
.set nb204nf_jindex, 736
.set nb204nf_jjnr, 744
.set nb204nf_shift, 752
.set nb204nf_shiftvec, 760
.set nb204nf_facel, 768
.set nb204nf_innerk, 776
.set nb204nf_n, 780
.set nb204nf_nn1, 784
.set nb204nf_nouter, 788
.set nb204nf_ninner, 792

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $808,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb204nf_nouter(%rsp)
        movl %eax,nb204nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb204nf_nri(%rsp)
        movq %rsi,nb204nf_iinr(%rsp)
        movq %rdx,nb204nf_jindex(%rsp)
        movq %rcx,nb204nf_jjnr(%rsp)
        movq %r8,nb204nf_shift(%rsp)
        movq %r9,nb204nf_shiftvec(%rsp)
        movq nb204nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb204nf_facel(%rsp)

        movq nb204nf_argkrf(%rbp),%rsi
        movq nb204nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb204nf_krf(%rsp)
        movaps %xmm2,nb204nf_crf(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb204nf_half(%rsp)
        movss nb204nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb204nf_half(%rsp)
        movaps %xmm3,nb204nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb204nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb204nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%rdx,%rbx,4),%xmm5
        movq nb204nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb204nf_facel(%rsp),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204nf_qqHH(%rsp)
        movaps %xmm4,nb204nf_qqMH(%rsp)
        movaps %xmm5,nb204nf_qqMM(%rsp)

_nb_kernel204nf_x86_64_sse.nb204nf_threadloop: 
        movq  nb204nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel204nf_x86_64_sse.nb204nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204nf_x86_64_sse.nb204nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204nf_n(%rsp)
        movl %ebx,nb204nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204nf_x86_64_sse.nb204nf_outerstart
        jmp _nb_kernel204nf_x86_64_sse.nb204nf_end

_nb_kernel204nf_x86_64_sse.nb204nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204nf_nouter(%rsp),%ebx
        movl %ebx,nb204nf_nouter(%rsp)

_nb_kernel204nf_x86_64_sse.nb204nf_outer: 
        movq  nb204nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb204nf_is3(%rsp)            ## store is3 

        movq  nb204nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb204nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb204nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb204nf_ii3(%rsp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204nf_ixH1(%rsp)
        movaps %xmm4,nb204nf_iyH1(%rsp)
        movaps %xmm5,nb204nf_izH1(%rsp)

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
        movaps %xmm0,nb204nf_ixH2(%rsp)
        movaps %xmm1,nb204nf_iyH2(%rsp)
        movaps %xmm2,nb204nf_izH2(%rsp)
        movaps %xmm3,nb204nf_ixM(%rsp)
        movaps %xmm4,nb204nf_iyM(%rsp)
        movaps %xmm5,nb204nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb204nf_vctot(%rsp)

        movq  nb204nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb204nf_pos(%rbp),%rsi
        movq  nb204nf_faction(%rbp),%rdi
        movq  nb204nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb204nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb204nf_ninner(%rsp),%ecx
        movl  %ecx,nb204nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb204nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel204nf_x86_64_sse.nb204nf_unroll_loop
        jmp   _nb_kernel204nf_x86_64_sse.nb204nf_single_check
_nb_kernel204nf_x86_64_sse.nb204nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb204nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 

        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb204nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb204nf_pos(%rbp),%rsi     ## base of pos[] 

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

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb204nf_jxH1(%rsp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb204nf_jyH1(%rsp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb204nf_jxH2(%rsp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb204nf_jyH2(%rsp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb204nf_jxM(%rsp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb204nf_jyM(%rsp)

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
        movaps %xmm0,nb204nf_jzH1(%rsp)
        movaps %xmm1,nb204nf_jzH2(%rsp)
        movaps %xmm2,nb204nf_jzM(%rsp)

        movaps nb204nf_ixH1(%rsp),%xmm0
        movaps nb204nf_iyH1(%rsp),%xmm1
        movaps nb204nf_izH1(%rsp),%xmm2
        movaps nb204nf_ixH1(%rsp),%xmm3
        movaps nb204nf_iyH1(%rsp),%xmm4
        movaps nb204nf_izH1(%rsp),%xmm5
        subps  nb204nf_jxH1(%rsp),%xmm0
        subps  nb204nf_jyH1(%rsp),%xmm1
        subps  nb204nf_jzH1(%rsp),%xmm2
        subps  nb204nf_jxH2(%rsp),%xmm3
        subps  nb204nf_jyH2(%rsp),%xmm4
        subps  nb204nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb204nf_rsqH1H1(%rsp)
        movaps %xmm3,nb204nf_rsqH1H2(%rsp)

        movaps nb204nf_ixH1(%rsp),%xmm0
        movaps nb204nf_iyH1(%rsp),%xmm1
        movaps nb204nf_izH1(%rsp),%xmm2
        movaps nb204nf_ixH2(%rsp),%xmm3
        movaps nb204nf_iyH2(%rsp),%xmm4
        movaps nb204nf_izH2(%rsp),%xmm5
        subps  nb204nf_jxM(%rsp),%xmm0
        subps  nb204nf_jyM(%rsp),%xmm1
        subps  nb204nf_jzM(%rsp),%xmm2
        subps  nb204nf_jxH1(%rsp),%xmm3
        subps  nb204nf_jyH1(%rsp),%xmm4
        subps  nb204nf_jzH1(%rsp),%xmm5
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
        movaps %xmm0,nb204nf_rsqH1M(%rsp)
        movaps %xmm3,nb204nf_rsqH2H1(%rsp)

        movaps nb204nf_ixH2(%rsp),%xmm0
        movaps nb204nf_iyH2(%rsp),%xmm1
        movaps nb204nf_izH2(%rsp),%xmm2
        movaps nb204nf_ixH2(%rsp),%xmm3
        movaps nb204nf_iyH2(%rsp),%xmm4
        movaps nb204nf_izH2(%rsp),%xmm5
        subps  nb204nf_jxH2(%rsp),%xmm0
        subps  nb204nf_jyH2(%rsp),%xmm1
        subps  nb204nf_jzH2(%rsp),%xmm2
        subps  nb204nf_jxM(%rsp),%xmm3
        subps  nb204nf_jyM(%rsp),%xmm4
        subps  nb204nf_jzM(%rsp),%xmm5
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
        movaps %xmm0,nb204nf_rsqH2H2(%rsp)
        movaps %xmm3,nb204nf_rsqH2M(%rsp)

        movaps nb204nf_ixM(%rsp),%xmm0
        movaps nb204nf_iyM(%rsp),%xmm1
        movaps nb204nf_izM(%rsp),%xmm2
        movaps nb204nf_ixM(%rsp),%xmm3
        movaps nb204nf_iyM(%rsp),%xmm4
        movaps nb204nf_izM(%rsp),%xmm5
        subps  nb204nf_jxH1(%rsp),%xmm0
        subps  nb204nf_jyH1(%rsp),%xmm1
        subps  nb204nf_jzH1(%rsp),%xmm2
        subps  nb204nf_jxH2(%rsp),%xmm3
        subps  nb204nf_jyH2(%rsp),%xmm4
        subps  nb204nf_jzH2(%rsp),%xmm5
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
        movaps %xmm0,nb204nf_rsqMH1(%rsp)
        movaps %xmm4,nb204nf_rsqMH2(%rsp)

        movaps nb204nf_ixM(%rsp),%xmm0
        movaps nb204nf_iyM(%rsp),%xmm1
        movaps nb204nf_izM(%rsp),%xmm2
        subps  nb204nf_jxM(%rsp),%xmm0
        subps  nb204nf_jyM(%rsp),%xmm1
        subps  nb204nf_jzM(%rsp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb204nf_rsqMM(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%rsp),%xmm3   ## rinvH2H2 
        mulps   nb204nf_half(%rsp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb204nf_rinvMM(%rsp)
        movaps  %xmm7,nb204nf_rinvMH2(%rsp)

        rsqrtps nb204nf_rsqH1H1(%rsp),%xmm1
        rsqrtps nb204nf_rsqH1H2(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH1H1(%rsp),%xmm1
        mulps   nb204nf_rsqH1H2(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%rsp),%xmm3
        mulps   nb204nf_half(%rsp),%xmm7
        movaps  %xmm3,nb204nf_rinvH1H1(%rsp)
        movaps  %xmm7,nb204nf_rinvH1H2(%rsp)

        rsqrtps nb204nf_rsqH1M(%rsp),%xmm1
        rsqrtps nb204nf_rsqH2H1(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH1M(%rsp),%xmm1
        mulps   nb204nf_rsqH2H1(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%rsp),%xmm3
        mulps   nb204nf_half(%rsp),%xmm7
        movaps  %xmm3,nb204nf_rinvH1M(%rsp)
        movaps  %xmm7,nb204nf_rinvH2H1(%rsp)

        rsqrtps nb204nf_rsqH2H2(%rsp),%xmm1
        rsqrtps nb204nf_rsqH2M(%rsp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH2H2(%rsp),%xmm1
        mulps   nb204nf_rsqH2M(%rsp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%rsp),%xmm3
        mulps   nb204nf_half(%rsp),%xmm7
        movaps  %xmm3,nb204nf_rinvH2H2(%rsp)
        movaps  %xmm7,nb204nf_rinvH2M(%rsp)

        rsqrtps nb204nf_rsqMH1(%rsp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb204nf_three(%rsp),%xmm3
        mulps   nb204nf_rsqMH1(%rsp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204nf_half(%rsp),%xmm3
        movaps  %xmm3,nb204nf_rinvMH1(%rsp)

        ## Coulomb interactions 
        ## add all H-H rsq in xmm2, H-M rsq xmm4
        ## H-H rinv in xmm0, H-M in xmm1
        movaps nb204nf_rinvH1H1(%rsp),%xmm0
        movaps nb204nf_rinvH1M(%rsp),%xmm1
        movaps nb204nf_rsqH1H1(%rsp),%xmm2
        movaps nb204nf_rsqH1M(%rsp),%xmm4
        addps  nb204nf_rinvH1H2(%rsp),%xmm0
        addps  nb204nf_rinvH2M(%rsp),%xmm1
        addps  nb204nf_rsqH1H2(%rsp),%xmm2
        addps  nb204nf_rsqH2M(%rsp),%xmm4
        addps  nb204nf_rinvH2H1(%rsp),%xmm0
        addps  nb204nf_rinvMH1(%rsp),%xmm1
        addps  nb204nf_rsqH2H1(%rsp),%xmm2
        addps  nb204nf_rsqMH1(%rsp),%xmm4
        addps  nb204nf_rinvH2H2(%rsp),%xmm0
        addps  nb204nf_rinvMH2(%rsp),%xmm1
        addps  nb204nf_rsqH2H2(%rsp),%xmm2
        addps  nb204nf_rsqMH2(%rsp),%xmm4
        movaps nb204nf_krf(%rsp),%xmm5
        movaps nb204nf_crf(%rsp),%xmm6

        ## calc 4*crf in xmm7
        movaps %xmm6,%xmm7
        addps  %xmm7,%xmm7
        addps  %xmm7,%xmm7
        mulps  %xmm5,%xmm2 ## H-H krsq
        mulps  %xmm5,%xmm4 ## H-M krsq
        addps  %xmm2,%xmm0 ## H-H rinv+krsq
        addps  %xmm4,%xmm1 ## H-M rinv+krsq
        subps  %xmm7,%xmm0 ## H-H rinv+krsq-crf
        subps  %xmm7,%xmm1 ## H-M rinv+krsq-crf
        mulps  nb204nf_qqHH(%rsp),%xmm0
        mulps  nb204nf_qqMH(%rsp),%xmm1
        addps  %xmm1,%xmm0
        addps  nb204nf_vctot(%rsp),%xmm0
        ## M-M interaction
        movaps nb204nf_rinvMM(%rsp),%xmm4
        mulps  nb204nf_rsqMM(%rsp),%xmm5   ## krsq
        addps  %xmm4,%xmm5                 ## rinv+krsq
        subps nb204nf_crf(%rsp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb204nf_qqMM(%rsp),%xmm5
        addps %xmm0,%xmm5
        movaps %xmm5,nb204nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb204nf_innerk(%rsp)
        jl    _nb_kernel204nf_x86_64_sse.nb204nf_single_check
        jmp   _nb_kernel204nf_x86_64_sse.nb204nf_unroll_loop
_nb_kernel204nf_x86_64_sse.nb204nf_single_check: 
        addl $4,nb204nf_innerk(%rsp)
        jnz   _nb_kernel204nf_x86_64_sse.nb204nf_single_loop
        jmp   _nb_kernel204nf_x86_64_sse.nb204nf_updateouterdata
_nb_kernel204nf_x86_64_sse.nb204nf_single_loop: 
        movq  nb204nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb204nf_innerjjnr(%rsp)

        movq nb204nf_pos(%rbp),%rsi
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
        movaps  nb204nf_ixM(%rsp),%xmm0
        movaps  nb204nf_iyM(%rsp),%xmm1
        movaps  nb204nf_izM(%rsp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## 11100100     ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## 01000100     ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM
        movaps %xmm3,nb204nf_jxM(%rsp)
        movaps %xmm4,nb204nf_jyM(%rsp)
        movaps %xmm5,nb204nf_jzM(%rsp)
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
        mulps   nb204nf_krf(%rsp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb204nf_three(%rsp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204nf_half(%rsp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb204nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb204nf_qqMM(%rsp),%xmm4
        movhps  nb204nf_qqMH(%rsp),%xmm4
        mulps %xmm4,%xmm6       ## vcoul  
        addps   nb204nf_vctot(%rsp),%xmm6
        movaps  %xmm6,nb204nf_vctot(%rsp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb204nf_ixH1(%rsp),%xmm0
        movaps  nb204nf_iyH1(%rsp),%xmm1
        movaps  nb204nf_izH1(%rsp),%xmm2
        movaps  nb204nf_ixH2(%rsp),%xmm3
        movaps  nb204nf_iyH2(%rsp),%xmm4
        movaps  nb204nf_izH2(%rsp),%xmm5
        subps   nb204nf_jxM(%rsp),%xmm0
        subps   nb204nf_jyM(%rsp),%xmm1
        subps   nb204nf_jzM(%rsp),%xmm2
        subps   nb204nf_jxM(%rsp),%xmm3
        subps   nb204nf_jyM(%rsp),%xmm4
        subps   nb204nf_jzM(%rsp),%xmm5
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
        movaps  nb204nf_three(%rsp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%rsp),%xmm3   ## rinv H1 - j water 
        mulps   nb204nf_half(%rsp),%xmm7   ## rinv H2 - j water  

        mulps nb204nf_krf(%rsp),%xmm0   ## krsq 
        mulps nb204nf_krf(%rsp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb204nf_qqMH(%rsp),%xmm6
        movhps  nb204nf_qqHH(%rsp),%xmm6
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb204nf_crf(%rsp),%xmm0
        subps   nb204nf_crf(%rsp),%xmm4
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb204nf_vctot(%rsp),%xmm4
        movaps  %xmm4,nb204nf_vctot(%rsp)
        decl nb204nf_innerk(%rsp)
        jz    _nb_kernel204nf_x86_64_sse.nb204nf_updateouterdata
        jmp   _nb_kernel204nf_x86_64_sse.nb204nf_single_loop
_nb_kernel204nf_x86_64_sse.nb204nf_updateouterdata: 
        ## get n from stack
        movl nb204nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb204nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb204nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb204nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb204nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204nf_x86_64_sse.nb204nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204nf_n(%rsp)
        jmp _nb_kernel204nf_x86_64_sse.nb204nf_outer
_nb_kernel204nf_x86_64_sse.nb204nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204nf_x86_64_sse.nb204nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204nf_x86_64_sse.nb204nf_threadloop
_nb_kernel204nf_x86_64_sse.nb204nf_end: 

        movl nb204nf_nouter(%rsp),%eax
        movl nb204nf_ninner(%rsp),%ebx
        movq nb204nf_outeriter(%rbp),%rcx
        movq nb204nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $808,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

