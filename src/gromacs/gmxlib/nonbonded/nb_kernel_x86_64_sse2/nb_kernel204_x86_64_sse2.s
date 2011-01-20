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






.globl nb_kernel204_x86_64_sse2
.globl _nb_kernel204_x86_64_sse2
nb_kernel204_x86_64_sse2:       
_nb_kernel204_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
.set nb204_ixM, 0
.set nb204_iyM, 16
.set nb204_izM, 32
.set nb204_ixH1, 48
.set nb204_iyH1, 64
.set nb204_izH1, 80
.set nb204_ixH2, 96
.set nb204_iyH2, 112
.set nb204_izH2, 128
.set nb204_jxM, 144
.set nb204_jyM, 160
.set nb204_jzM, 176
.set nb204_jxH1, 192
.set nb204_jyH1, 208
.set nb204_jzH1, 224
.set nb204_jxH2, 240
.set nb204_jyH2, 256
.set nb204_jzH2, 272
.set nb204_dxMM, 288
.set nb204_dyMM, 304
.set nb204_dzMM, 320
.set nb204_dxMH1, 336
.set nb204_dyMH1, 352
.set nb204_dzMH1, 368
.set nb204_dxMH2, 384
.set nb204_dyMH2, 400
.set nb204_dzMH2, 416
.set nb204_dxH1M, 432
.set nb204_dyH1M, 448
.set nb204_dzH1M, 464
.set nb204_dxH1H1, 480
.set nb204_dyH1H1, 496
.set nb204_dzH1H1, 512
.set nb204_dxH1H2, 528
.set nb204_dyH1H2, 544
.set nb204_dzH1H2, 560
.set nb204_dxH2M, 576
.set nb204_dyH2M, 592
.set nb204_dzH2M, 608
.set nb204_dxH2H1, 624
.set nb204_dyH2H1, 640
.set nb204_dzH2H1, 656
.set nb204_dxH2H2, 672
.set nb204_dyH2H2, 688
.set nb204_dzH2H2, 704
.set nb204_qqMM, 720
.set nb204_qqMH, 736
.set nb204_qqHH, 752
.set nb204_vctot, 768
.set nb204_fixM, 784
.set nb204_fiyM, 800
.set nb204_fizM, 816
.set nb204_fixH1, 832
.set nb204_fiyH1, 848
.set nb204_fizH1, 864
.set nb204_fixH2, 880
.set nb204_fiyH2, 896
.set nb204_fizH2, 912
.set nb204_fjxM, 928
.set nb204_fjyM, 944
.set nb204_fjzM, 960
.set nb204_fjxH1, 976
.set nb204_fjyH1, 992
.set nb204_fjzH1, 1008
.set nb204_fjxH2, 1024
.set nb204_fjyH2, 1040
.set nb204_fjzH2, 1056
.set nb204_half, 1072
.set nb204_three, 1088
.set nb204_rsqMM, 1104
.set nb204_rsqMH1, 1120
.set nb204_rsqMH2, 1136
.set nb204_rsqH1M, 1152
.set nb204_rsqH1H1, 1168
.set nb204_rsqH1H2, 1184
.set nb204_rsqH2M, 1200
.set nb204_rsqH2H1, 1216
.set nb204_rsqH2H2, 1232
.set nb204_rinvMM, 1248
.set nb204_rinvMH1, 1264
.set nb204_rinvMH2, 1280
.set nb204_rinvH1M, 1296
.set nb204_rinvH1H1, 1312
.set nb204_rinvH1H2, 1328
.set nb204_rinvH2M, 1344
.set nb204_rinvH2H1, 1360
.set nb204_rinvH2H2, 1376
.set nb204_two, 1392
.set nb204_krf, 1408
.set nb204_crf, 1424
.set nb204_is3, 1440
.set nb204_ii3, 1444
.set nb204_nri, 1448
.set nb204_iinr, 1456
.set nb204_jindex, 1464
.set nb204_jjnr, 1472
.set nb204_shift, 1480
.set nb204_shiftvec, 1488
.set nb204_facel, 1496
.set nb204_innerjjnr, 1504
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb204_facel(%rsp)

        movq nb204_argkrf(%rbp),%rsi
        movq nb204_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb204_krf(%rsp)
        movapd %xmm2,nb204_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb204_half(%rsp)
        movl %ebx,nb204_half+4(%rsp)
        movsd nb204_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb204_half(%rsp)
        movapd %xmm2,nb204_two(%rsp)
        movapd %xmm3,nb204_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb204_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb204_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb204_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204_qqMM(%rsp)
        movapd %xmm4,nb204_qqMH(%rsp)
        movapd %xmm5,nb204_qqHH(%rsp)

_nb_kernel204_x86_64_sse2.nb204_threadloop: 
        movq  nb204_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel204_x86_64_sse2.nb204_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204_x86_64_sse2.nb204_spinlock

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
        jg  _nb_kernel204_x86_64_sse2.nb204_outerstart
        jmp _nb_kernel204_x86_64_sse2.nb204_end

_nb_kernel204_x86_64_sse2.nb204_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204_nouter(%rsp),%ebx
        movl %ebx,nb204_nouter(%rsp)

_nb_kernel204_x86_64_sse2.nb204_outer: 
        movq  nb204_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb204_is3(%rsp)      ## store is3 

        movq  nb204_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb204_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb204_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb204_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204_ixH1(%rsp)
        movapd %xmm4,nb204_iyH1(%rsp)
        movapd %xmm5,nb204_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb204_ixH2(%rsp)
        movapd %xmm1,nb204_iyH2(%rsp)
        movapd %xmm2,nb204_izH2(%rsp)
        movapd %xmm3,nb204_ixM(%rsp)
        movapd %xmm4,nb204_iyM(%rsp)
        movapd %xmm5,nb204_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb204_vctot(%rsp)
        movapd %xmm4,nb204_fixM(%rsp)
        movapd %xmm4,nb204_fiyM(%rsp)
        movapd %xmm4,nb204_fizM(%rsp)
        movapd %xmm4,nb204_fixH1(%rsp)
        movapd %xmm4,nb204_fiyH1(%rsp)
        movapd %xmm4,nb204_fizH1(%rsp)
        movapd %xmm4,nb204_fixH2(%rsp)
        movapd %xmm4,nb204_fiyH2(%rsp)
        movapd %xmm4,nb204_fizH2(%rsp)

        movq  nb204_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb204_pos(%rbp),%rsi
        movq  nb204_faction(%rbp),%rdi
        movq  nb204_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb204_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb204_ninner(%rsp),%ecx
        movl  %ecx,nb204_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb204_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel204_x86_64_sse2.nb204_unroll_loop
        jmp   _nb_kernel204_x86_64_sse2.nb204_checksingle
_nb_kernel204_x86_64_sse2.nb204_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb204_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb204_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb204_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j H1 coordinates to local temp variables 
    movlpd 24(%rsi,%rax,8),%xmm0
    movlpd 32(%rsi,%rax,8),%xmm1
    movlpd 40(%rsi,%rax,8),%xmm2
    movhpd 24(%rsi,%rbx,8),%xmm0
    movhpd 32(%rsi,%rbx,8),%xmm1
    movhpd 40(%rsi,%rbx,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb204_ixH1(%rsp),%xmm0
    subpd nb204_iyH1(%rsp),%xmm1
    subpd nb204_izH1(%rsp),%xmm2
    subpd nb204_ixH2(%rsp),%xmm3
    subpd nb204_iyH2(%rsp),%xmm4
    subpd nb204_izH2(%rsp),%xmm5
    subpd nb204_ixM(%rsp),%xmm6
    subpd nb204_iyM(%rsp),%xmm7
    subpd nb204_izM(%rsp),%xmm8

        movapd %xmm0,nb204_dxH1H1(%rsp)
        movapd %xmm1,nb204_dyH1H1(%rsp)
        movapd %xmm2,nb204_dzH1H1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH2H1(%rsp)
        movapd %xmm4,nb204_dyH2H1(%rsp)
        movapd %xmm5,nb204_dzH2H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb204_dxMH1(%rsp)
        movapd %xmm7,nb204_dyMH1(%rsp)
        movapd %xmm8,nb204_dzMH1(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb204_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb204_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H1 
        mulpd   %xmm15,%xmm10 ##   rinvH2H1
    mulpd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb204_krf(%rsp),%xmm2
    mulpd  %xmm9,%xmm9  ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  %xmm2,%xmm0 ## k*rsq
    mulpd  %xmm2,%xmm3
    mulpd  %xmm2,%xmm6
    movapd %xmm0,%xmm2 ## copy of k*rsq
    movapd %xmm3,%xmm5
    movapd %xmm6,%xmm8
    addpd  %xmm1,%xmm2 ## rinv+krsq
    addpd  %xmm4,%xmm5
    addpd  %xmm7,%xmm8
    movapd nb204_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb204_qqHH(%rsp),%xmm12
    movapd nb204_qqMH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb204_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb204_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb204_dxH1H1(%rsp),%xmm7
        mulpd nb204_dyH1H1(%rsp),%xmm8
        mulpd nb204_dzH1H1(%rsp),%xmm9
        mulpd nb204_dxH2H1(%rsp),%xmm10
        mulpd nb204_dyH2H1(%rsp),%xmm11
        mulpd nb204_dzH2H1(%rsp),%xmm12
        mulpd nb204_dxMH1(%rsp),%xmm13
        mulpd nb204_dyMH1(%rsp),%xmm14
        mulpd nb204_dzMH1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb204_fixH1(%rsp),%xmm7
    addpd nb204_fiyH1(%rsp),%xmm8
    addpd nb204_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb204_fixH2(%rsp),%xmm10
    addpd nb204_fiyH2(%rsp),%xmm11
    addpd nb204_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb204_fixM(%rsp),%xmm13
    addpd nb204_fiyM(%rsp),%xmm14
    addpd nb204_fizM(%rsp),%xmm15

    movapd %xmm7,nb204_fixH1(%rsp)
    movapd %xmm8,nb204_fiyH1(%rsp)
    movapd %xmm9,nb204_fizH1(%rsp)
    movapd %xmm10,nb204_fixH2(%rsp)
    movapd %xmm11,nb204_fiyH2(%rsp)
    movapd %xmm12,nb204_fizH2(%rsp)
    movapd %xmm13,nb204_fixM(%rsp)
    movapd %xmm14,nb204_fiyM(%rsp)
    movapd %xmm15,nb204_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movlpd 48(%rsi,%rax,8),%xmm0
    movlpd 56(%rsi,%rax,8),%xmm1
    movlpd 64(%rsi,%rax,8),%xmm2
    movhpd 48(%rsi,%rbx,8),%xmm0
    movhpd 56(%rsi,%rbx,8),%xmm1
    movhpd 64(%rsi,%rbx,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb204_ixH1(%rsp),%xmm0
    subpd nb204_iyH1(%rsp),%xmm1
    subpd nb204_izH1(%rsp),%xmm2
    subpd nb204_ixH2(%rsp),%xmm3
    subpd nb204_iyH2(%rsp),%xmm4
    subpd nb204_izH2(%rsp),%xmm5
    subpd nb204_ixM(%rsp),%xmm6
    subpd nb204_iyM(%rsp),%xmm7
    subpd nb204_izM(%rsp),%xmm8

        movapd %xmm0,nb204_dxH1H2(%rsp)
        movapd %xmm1,nb204_dyH1H2(%rsp)
        movapd %xmm2,nb204_dzH1H2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH2H2(%rsp)
        movapd %xmm4,nb204_dyH2H2(%rsp)
        movapd %xmm5,nb204_dzH2H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb204_dxMH2(%rsp)
        movapd %xmm7,nb204_dyMH2(%rsp)
        movapd %xmm8,nb204_dzMH2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb204_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb204_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H2
        mulpd   %xmm15,%xmm10 ##   rinvH2H2
    mulpd   %xmm15,%xmm11 ##   rinvMH2


        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb204_krf(%rsp),%xmm2
    mulpd  %xmm9,%xmm9  ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  %xmm2,%xmm0 ## k*rsq
    mulpd  %xmm2,%xmm3
    mulpd  %xmm2,%xmm6
    movapd %xmm0,%xmm2 ## copy of k*rsq
    movapd %xmm3,%xmm5
    movapd %xmm6,%xmm8
    addpd  %xmm1,%xmm2 ## rinv+krsq
    addpd  %xmm4,%xmm5
    addpd  %xmm7,%xmm8
    movapd nb204_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb204_qqHH(%rsp),%xmm12
    movapd nb204_qqMH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb204_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd  %xmm2,nb204_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb204_dxH1H2(%rsp),%xmm7
        mulpd nb204_dyH1H2(%rsp),%xmm8
        mulpd nb204_dzH1H2(%rsp),%xmm9
        mulpd nb204_dxH2H2(%rsp),%xmm10
        mulpd nb204_dyH2H2(%rsp),%xmm11
        mulpd nb204_dzH2H2(%rsp),%xmm12
        mulpd nb204_dxMH2(%rsp),%xmm13
        mulpd nb204_dyMH2(%rsp),%xmm14
        mulpd nb204_dzMH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb204_fixH1(%rsp),%xmm7
    addpd nb204_fiyH1(%rsp),%xmm8
    addpd nb204_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb204_fixH2(%rsp),%xmm10
    addpd nb204_fiyH2(%rsp),%xmm11
    addpd nb204_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb204_fixM(%rsp),%xmm13
    addpd nb204_fiyM(%rsp),%xmm14
    addpd nb204_fizM(%rsp),%xmm15

    movapd %xmm7,nb204_fixH1(%rsp)
    movapd %xmm8,nb204_fiyH1(%rsp)
    movapd %xmm9,nb204_fizH1(%rsp)
    movapd %xmm10,nb204_fixH2(%rsp)
    movapd %xmm11,nb204_fiyH2(%rsp)
    movapd %xmm12,nb204_fizH2(%rsp)
    movapd %xmm13,nb204_fixM(%rsp)
    movapd %xmm14,nb204_fiyM(%rsp)
    movapd %xmm15,nb204_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## move j M coordinates to local temp variables 
    movlpd 72(%rsi,%rax,8),%xmm0
    movlpd 80(%rsi,%rax,8),%xmm1
    movlpd 88(%rsi,%rax,8),%xmm2
    movhpd 72(%rsi,%rbx,8),%xmm0
    movhpd 80(%rsi,%rbx,8),%xmm1
    movhpd 88(%rsi,%rbx,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb204_ixH1(%rsp),%xmm0
    subpd nb204_iyH1(%rsp),%xmm1
    subpd nb204_izH1(%rsp),%xmm2
    subpd nb204_ixH2(%rsp),%xmm3
    subpd nb204_iyH2(%rsp),%xmm4
    subpd nb204_izH2(%rsp),%xmm5
    subpd nb204_ixM(%rsp),%xmm6
    subpd nb204_iyM(%rsp),%xmm7
    subpd nb204_izM(%rsp),%xmm8

        movapd %xmm0,nb204_dxH1M(%rsp)
        movapd %xmm1,nb204_dyH1M(%rsp)
        movapd %xmm2,nb204_dzH1M(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH2M(%rsp)
        movapd %xmm4,nb204_dyH2M(%rsp)
        movapd %xmm5,nb204_dzH2M(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb204_dxMM(%rsp)
        movapd %xmm7,nb204_dyMM(%rsp)
        movapd %xmm8,nb204_dzMM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb204_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb204_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb204_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1M
        mulpd   %xmm15,%xmm10 ##   rinvH2M
    mulpd   %xmm15,%xmm11 ##   rinvMM


        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb204_krf(%rsp),%xmm2
    mulpd  %xmm9,%xmm9  ## rinvsq
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  %xmm2,%xmm0 ## k*rsq
    mulpd  %xmm2,%xmm3
    mulpd  %xmm2,%xmm6
    movapd %xmm0,%xmm2 ## copy of k*rsq
    movapd %xmm3,%xmm5
    movapd %xmm6,%xmm8
    addpd  %xmm1,%xmm2 ## rinv+krsq
    addpd  %xmm4,%xmm5
    addpd  %xmm7,%xmm8
    movapd nb204_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb204_qqMH(%rsp),%xmm12
    movapd nb204_qqMM(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm12,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb204_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd  %xmm2,nb204_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j M forces to xmm0-xmm2
        movlpd 72(%rdi,%rax,8),%xmm0
        movlpd 80(%rdi,%rax,8),%xmm1
        movlpd 88(%rdi,%rax,8),%xmm2
        movhpd 72(%rdi,%rbx,8),%xmm0
        movhpd 80(%rdi,%rbx,8),%xmm1
        movhpd 88(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb204_dxH1M(%rsp),%xmm7
        mulpd nb204_dyH1M(%rsp),%xmm8
        mulpd nb204_dzH1M(%rsp),%xmm9
        mulpd nb204_dxH2M(%rsp),%xmm10
        mulpd nb204_dyH2M(%rsp),%xmm11
        mulpd nb204_dzH2M(%rsp),%xmm12
        mulpd nb204_dxMM(%rsp),%xmm13
        mulpd nb204_dyMM(%rsp),%xmm14
        mulpd nb204_dzMM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb204_fixH1(%rsp),%xmm7
    addpd nb204_fiyH1(%rsp),%xmm8
    addpd nb204_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb204_fixH2(%rsp),%xmm10
    addpd nb204_fiyH2(%rsp),%xmm11
    addpd nb204_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb204_fixM(%rsp),%xmm13
    addpd nb204_fiyM(%rsp),%xmm14
    addpd nb204_fizM(%rsp),%xmm15

    movapd %xmm7,nb204_fixH1(%rsp)
    movapd %xmm8,nb204_fiyH1(%rsp)
    movapd %xmm9,nb204_fizH1(%rsp)
    movapd %xmm10,nb204_fixH2(%rsp)
    movapd %xmm11,nb204_fiyH2(%rsp)
    movapd %xmm12,nb204_fizH2(%rsp)
    movapd %xmm13,nb204_fixM(%rsp)
    movapd %xmm14,nb204_fiyM(%rsp)
    movapd %xmm15,nb204_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movlpd %xmm0,72(%rdi,%rax,8)
        movlpd %xmm1,80(%rdi,%rax,8)
        movlpd %xmm2,88(%rdi,%rax,8)
        movhpd %xmm0,72(%rdi,%rbx,8)
        movhpd %xmm1,80(%rdi,%rbx,8)
        movhpd %xmm2,88(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb204_innerk(%rsp)
        jl    _nb_kernel204_x86_64_sse2.nb204_checksingle
        jmp   _nb_kernel204_x86_64_sse2.nb204_unroll_loop
_nb_kernel204_x86_64_sse2.nb204_checksingle: 
        movl  nb204_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel204_x86_64_sse2.nb204_dosingle
        jmp   _nb_kernel204_x86_64_sse2.nb204_updateouterdata
_nb_kernel204_x86_64_sse2.nb204_dosingle: 
        movq  nb204_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb204_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax


        ## move j H1 coordinates to local temp variables 
    movsd 24(%rsi,%rax,8),%xmm0
    movsd 32(%rsi,%rax,8),%xmm1
    movsd 40(%rsi,%rax,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb204_ixH1(%rsp),%xmm0
    subsd nb204_iyH1(%rsp),%xmm1
    subsd nb204_izH1(%rsp),%xmm2
    subsd nb204_ixH2(%rsp),%xmm3
    subsd nb204_iyH2(%rsp),%xmm4
    subsd nb204_izH2(%rsp),%xmm5
    subsd nb204_ixM(%rsp),%xmm6
    subsd nb204_iyM(%rsp),%xmm7
    subsd nb204_izM(%rsp),%xmm8

        movsd %xmm0,nb204_dxH1H1(%rsp)
        movsd %xmm1,nb204_dyH1H1(%rsp)
        movsd %xmm2,nb204_dzH1H1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb204_dxH2H1(%rsp)
        movsd %xmm4,nb204_dyH2H1(%rsp)
        movsd %xmm5,nb204_dzH2H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb204_dxMH1(%rsp)
        movsd %xmm7,nb204_dyMH1(%rsp)
        movsd %xmm8,nb204_dzMH1(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb204_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb204_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H1 
        mulsd   %xmm15,%xmm10 ##   rinvH2H1
    mulsd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb204_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd  %xmm0,%xmm2 ## copy of k*rsq
    movsd  %xmm3,%xmm5
    movsd  %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd  nb204_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd  nb204_qqHH(%rsp),%xmm12
    movsd  nb204_qqMH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb204_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb204_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb204_dxH1H1(%rsp),%xmm7
        mulsd nb204_dyH1H1(%rsp),%xmm8
        mulsd nb204_dzH1H1(%rsp),%xmm9
        mulsd nb204_dxH2H1(%rsp),%xmm10
        mulsd nb204_dyH2H1(%rsp),%xmm11
        mulsd nb204_dzH2H1(%rsp),%xmm12
        mulsd nb204_dxMH1(%rsp),%xmm13
        mulsd nb204_dyMH1(%rsp),%xmm14
        mulsd nb204_dzMH1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb204_fixH1(%rsp),%xmm7
    addsd nb204_fiyH1(%rsp),%xmm8
    addsd nb204_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb204_fixH2(%rsp),%xmm10
    addsd nb204_fiyH2(%rsp),%xmm11
    addsd nb204_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb204_fixM(%rsp),%xmm13
    addsd nb204_fiyM(%rsp),%xmm14
    addsd nb204_fizM(%rsp),%xmm15

    movsd %xmm7,nb204_fixH1(%rsp)
    movsd %xmm8,nb204_fiyH1(%rsp)
    movsd %xmm9,nb204_fizH1(%rsp)
    movsd %xmm10,nb204_fixH2(%rsp)
    movsd %xmm11,nb204_fiyH2(%rsp)
    movsd %xmm12,nb204_fizH2(%rsp)
    movsd %xmm13,nb204_fixM(%rsp)
    movsd %xmm14,nb204_fiyM(%rsp)
    movsd %xmm15,nb204_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movsd 48(%rsi,%rax,8),%xmm0
    movsd 56(%rsi,%rax,8),%xmm1
    movsd 64(%rsi,%rax,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb204_ixH1(%rsp),%xmm0
    subsd nb204_iyH1(%rsp),%xmm1
    subsd nb204_izH1(%rsp),%xmm2
    subsd nb204_ixH2(%rsp),%xmm3
    subsd nb204_iyH2(%rsp),%xmm4
    subsd nb204_izH2(%rsp),%xmm5
    subsd nb204_ixM(%rsp),%xmm6
    subsd nb204_iyM(%rsp),%xmm7
    subsd nb204_izM(%rsp),%xmm8

        movsd %xmm0,nb204_dxH1H2(%rsp)
        movsd %xmm1,nb204_dyH1H2(%rsp)
        movsd %xmm2,nb204_dzH1H2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb204_dxH2H2(%rsp)
        movsd %xmm4,nb204_dyH2H2(%rsp)
        movsd %xmm5,nb204_dzH2H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb204_dxMH2(%rsp)
        movsd %xmm7,nb204_dyMH2(%rsp)
        movsd %xmm8,nb204_dzMH2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb204_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb204_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H2
        mulsd   %xmm15,%xmm10 ##   rinvH2H2
    mulsd   %xmm15,%xmm11 ##   rinvMH2

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb204_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd %xmm0,%xmm2 ## copy of k*rsq
    movsd %xmm3,%xmm5
    movsd %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd nb204_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb204_qqHH(%rsp),%xmm12
    movsd nb204_qqMH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb204_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb204_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb204_dxH1H2(%rsp),%xmm7
        mulsd nb204_dyH1H2(%rsp),%xmm8
        mulsd nb204_dzH1H2(%rsp),%xmm9
        mulsd nb204_dxH2H2(%rsp),%xmm10
        mulsd nb204_dyH2H2(%rsp),%xmm11
        mulsd nb204_dzH2H2(%rsp),%xmm12
        mulsd nb204_dxMH2(%rsp),%xmm13
        mulsd nb204_dyMH2(%rsp),%xmm14
        mulsd nb204_dzMH2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb204_fixH1(%rsp),%xmm7
    addsd nb204_fiyH1(%rsp),%xmm8
    addsd nb204_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb204_fixH2(%rsp),%xmm10
    addsd nb204_fiyH2(%rsp),%xmm11
    addsd nb204_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb204_fixM(%rsp),%xmm13
    addsd nb204_fiyM(%rsp),%xmm14
    addsd nb204_fizM(%rsp),%xmm15

    movsd %xmm7,nb204_fixH1(%rsp)
    movsd %xmm8,nb204_fiyH1(%rsp)
    movsd %xmm9,nb204_fizH1(%rsp)
    movsd %xmm10,nb204_fixH2(%rsp)
    movsd %xmm11,nb204_fiyH2(%rsp)
    movsd %xmm12,nb204_fizH2(%rsp)
    movsd %xmm13,nb204_fixM(%rsp)
    movsd %xmm14,nb204_fiyM(%rsp)
    movsd %xmm15,nb204_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

        ## move j M coordinates to local temp variables 
    movsd 72(%rsi,%rax,8),%xmm0
    movsd 80(%rsi,%rax,8),%xmm1
    movsd 88(%rsi,%rax,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb204_ixH1(%rsp),%xmm0
    subsd nb204_iyH1(%rsp),%xmm1
    subsd nb204_izH1(%rsp),%xmm2
    subsd nb204_ixH2(%rsp),%xmm3
    subsd nb204_iyH2(%rsp),%xmm4
    subsd nb204_izH2(%rsp),%xmm5
    subsd nb204_ixM(%rsp),%xmm6
    subsd nb204_iyM(%rsp),%xmm7
    subsd nb204_izM(%rsp),%xmm8

        movsd %xmm0,nb204_dxH1M(%rsp)
        movsd %xmm1,nb204_dyH1M(%rsp)
        movsd %xmm2,nb204_dzH1M(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb204_dxH2M(%rsp)
        movsd %xmm4,nb204_dyH2M(%rsp)
        movsd %xmm5,nb204_dzH2M(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb204_dxMM(%rsp)
        movsd %xmm7,nb204_dyMM(%rsp)
        movsd %xmm8,nb204_dzMM(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb204_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb204_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb204_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1M
        mulsd   %xmm15,%xmm10 ##   rinvH2M
    mulsd   %xmm15,%xmm11 ##   rinvMM

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb204_krf(%rsp),%xmm2
    mulsd  %xmm9,%xmm9  ## rinvsq
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  %xmm2,%xmm0 ## k*rsq
    mulsd  %xmm2,%xmm3
    mulsd  %xmm2,%xmm6
    movsd %xmm0,%xmm2 ## copy of k*rsq
    movsd %xmm3,%xmm5
    movsd %xmm6,%xmm8
    addsd  %xmm1,%xmm2 ## rinv+krsq
    addsd  %xmm4,%xmm5
    addsd  %xmm7,%xmm8
    movsd nb204_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb204_qqMH(%rsp),%xmm12
    movsd nb204_qqMM(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm12,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm12,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb204_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd  %xmm2,nb204_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j M forces to xmm0-xmm2
        movsd 72(%rdi,%rax,8),%xmm0
        movsd 80(%rdi,%rax,8),%xmm1
        movsd 88(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb204_dxH1M(%rsp),%xmm7
        mulsd nb204_dyH1M(%rsp),%xmm8
        mulsd nb204_dzH1M(%rsp),%xmm9
        mulsd nb204_dxH2M(%rsp),%xmm10
        mulsd nb204_dyH2M(%rsp),%xmm11
        mulsd nb204_dzH2M(%rsp),%xmm12
        mulsd nb204_dxMM(%rsp),%xmm13
        mulsd nb204_dyMM(%rsp),%xmm14
        mulsd nb204_dzMM(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb204_fixH1(%rsp),%xmm7
    addsd nb204_fiyH1(%rsp),%xmm8
    addsd nb204_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb204_fixH2(%rsp),%xmm10
    addsd nb204_fiyH2(%rsp),%xmm11
    addsd nb204_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb204_fixM(%rsp),%xmm13
    addsd nb204_fiyM(%rsp),%xmm14
    addsd nb204_fizM(%rsp),%xmm15

    movsd %xmm7,nb204_fixH1(%rsp)
    movsd %xmm8,nb204_fiyH1(%rsp)
    movsd %xmm9,nb204_fizH1(%rsp)
    movsd %xmm10,nb204_fixH2(%rsp)
    movsd %xmm11,nb204_fiyH2(%rsp)
    movsd %xmm12,nb204_fizH2(%rsp)
    movsd %xmm13,nb204_fixM(%rsp)
    movsd %xmm14,nb204_fiyM(%rsp)
    movsd %xmm15,nb204_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movsd %xmm0,72(%rdi,%rax,8)
        movsd %xmm1,80(%rdi,%rax,8)
        movsd %xmm2,88(%rdi,%rax,8)

_nb_kernel204_x86_64_sse2.nb204_updateouterdata: 
        movl  nb204_ii3(%rsp),%ecx
        movq  nb204_faction(%rbp),%rdi
        movq  nb204_fshift(%rbp),%rsi
        movl  nb204_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixH1(%rsp),%xmm0
        movapd nb204_fiyH1(%rsp),%xmm1
        movapd nb204_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixH2(%rsp),%xmm0
        movapd nb204_fiyH2(%rsp),%xmm1
        movapd nb204_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixM(%rsp),%xmm0
        movapd nb204_fiyM(%rsp),%xmm1
        movapd nb204_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  72(%rdi,%rcx,8),%xmm3
        movsd  80(%rdi,%rcx,8),%xmm4
        movsd  88(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,72(%rdi,%rcx,8)
        movsd  %xmm4,80(%rdi,%rcx,8)
        movsd  %xmm5,88(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb204_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb204_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb204_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb204_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb204_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204_x86_64_sse2.nb204_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204_n(%rsp)
        jmp _nb_kernel204_x86_64_sse2.nb204_outer
_nb_kernel204_x86_64_sse2.nb204_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204_x86_64_sse2.nb204_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204_x86_64_sse2.nb204_threadloop
_nb_kernel204_x86_64_sse2.nb204_end: 
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





.globl nb_kernel204nf_x86_64_sse2
.globl _nb_kernel204nf_x86_64_sse2
nb_kernel204nf_x86_64_sse2:     
_nb_kernel204nf_x86_64_sse2:    
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
.set nb204nf_ixM, 0
.set nb204nf_iyM, 16
.set nb204nf_izM, 32
.set nb204nf_ixH1, 48
.set nb204nf_iyH1, 64
.set nb204nf_izH1, 80
.set nb204nf_ixH2, 96
.set nb204nf_iyH2, 112
.set nb204nf_izH2, 128
.set nb204nf_jxM, 144
.set nb204nf_jyM, 160
.set nb204nf_jzM, 176
.set nb204nf_jxH1, 192
.set nb204nf_jyH1, 208
.set nb204nf_jzH1, 224
.set nb204nf_jxH2, 240
.set nb204nf_jyH2, 256
.set nb204nf_jzH2, 272
.set nb204nf_qqMM, 288
.set nb204nf_qqMH, 304
.set nb204nf_qqHH, 320
.set nb204nf_vctot, 336
.set nb204nf_half, 352
.set nb204nf_three, 368
.set nb204nf_rsqMM, 384
.set nb204nf_rsqMH1, 400
.set nb204nf_rsqMH2, 416
.set nb204nf_rsqH1M, 432
.set nb204nf_rsqH1H1, 448
.set nb204nf_rsqH1H2, 464
.set nb204nf_rsqH2M, 480
.set nb204nf_rsqH2H1, 496
.set nb204nf_rsqH2H2, 512
.set nb204nf_rinvMM, 528
.set nb204nf_rinvMH1, 544
.set nb204nf_rinvMH2, 560
.set nb204nf_rinvH1M, 576
.set nb204nf_rinvH1H1, 592
.set nb204nf_rinvH1H2, 608
.set nb204nf_rinvH2M, 624
.set nb204nf_rinvH2H1, 640
.set nb204nf_rinvH2H2, 656
.set nb204nf_krf, 672
.set nb204nf_crf, 688
.set nb204nf_is3, 704
.set nb204nf_ii3, 708
.set nb204nf_nri, 712
.set nb204nf_iinr, 720
.set nb204nf_jindex, 728
.set nb204nf_jjnr, 736
.set nb204nf_shift, 744
.set nb204nf_shiftvec, 752
.set nb204nf_facel, 760
.set nb204nf_innerjjnr, 768
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb204nf_facel(%rsp)

        movq nb204nf_argkrf(%rbp),%rsi
        movq nb204nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb204nf_krf(%rsp)
        movapd %xmm2,nb204nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb204nf_half(%rsp)
        movl %ebx,nb204nf_half+4(%rsp)
        movsd nb204nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb204nf_half(%rsp)
        movapd %xmm3,nb204nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb204nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb204nf_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb204nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204nf_qqMM(%rsp)
        movapd %xmm4,nb204nf_qqMH(%rsp)
        movapd %xmm5,nb204nf_qqHH(%rsp)

_nb_kernel204nf_x86_64_sse2.nb204nf_threadloop: 
        movq  nb204nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel204nf_x86_64_sse2.nb204nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204nf_x86_64_sse2.nb204nf_spinlock

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
        jg  _nb_kernel204nf_x86_64_sse2.nb204nf_outerstart
        jmp _nb_kernel204nf_x86_64_sse2.nb204nf_end

_nb_kernel204nf_x86_64_sse2.nb204nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204nf_nouter(%rsp),%ebx
        movl %ebx,nb204nf_nouter(%rsp)

_nb_kernel204nf_x86_64_sse2.nb204nf_outer: 
        movq  nb204nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb204nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb204nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb204nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb204nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204nf_ixH1(%rsp)
        movapd %xmm4,nb204nf_iyH1(%rsp)
        movapd %xmm5,nb204nf_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb204nf_ixH2(%rsp)
        movapd %xmm1,nb204nf_iyH2(%rsp)
        movapd %xmm2,nb204nf_izH2(%rsp)
        movapd %xmm3,nb204nf_ixM(%rsp)
        movapd %xmm4,nb204nf_iyM(%rsp)
        movapd %xmm5,nb204nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb204nf_vctot(%rsp)

        movq  nb204nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb204nf_pos(%rbp),%rsi
        movq  nb204nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb204nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb204nf_ninner(%rsp),%ecx
        movl  %ecx,nb204nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb204nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel204nf_x86_64_sse2.nb204nf_unroll_loop
        jmp   _nb_kernel204nf_x86_64_sse2.nb204nf_checksingle
_nb_kernel204nf_x86_64_sse2.nb204nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb204nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb204nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb204nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        movlpd 24(%rsi,%rax,8),%xmm2
        movlpd 32(%rsi,%rax,8),%xmm3
        movlpd 40(%rsi,%rax,8),%xmm4
        movlpd 48(%rsi,%rax,8),%xmm5
        movlpd 56(%rsi,%rax,8),%xmm6
        movlpd 64(%rsi,%rax,8),%xmm7
        movhpd 24(%rsi,%rbx,8),%xmm2
        movhpd 32(%rsi,%rbx,8),%xmm3
        movhpd 40(%rsi,%rbx,8),%xmm4
        movhpd 48(%rsi,%rbx,8),%xmm5
        movhpd 56(%rsi,%rbx,8),%xmm6
        movhpd 64(%rsi,%rbx,8),%xmm7
        movapd  %xmm2,nb204nf_jxH1(%rsp)
        movapd  %xmm3,nb204nf_jyH1(%rsp)
        movapd  %xmm4,nb204nf_jzH1(%rsp)
        movapd  %xmm5,nb204nf_jxH2(%rsp)
        movapd  %xmm6,nb204nf_jyH2(%rsp)
        movapd  %xmm7,nb204nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movhpd 72(%rsi,%rbx,8),%xmm2
        movhpd 80(%rsi,%rbx,8),%xmm3
        movhpd 88(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb204nf_jxM(%rsp)
        movapd  %xmm3,nb204nf_jyM(%rsp)
        movapd  %xmm4,nb204nf_jzM(%rsp)

        movapd nb204nf_ixM(%rsp),%xmm0
        movapd nb204nf_iyM(%rsp),%xmm1
        movapd nb204nf_izM(%rsp),%xmm2
        movapd nb204nf_ixM(%rsp),%xmm3
        movapd nb204nf_iyM(%rsp),%xmm4
        movapd nb204nf_izM(%rsp),%xmm5
        subpd  nb204nf_jxM(%rsp),%xmm0
        subpd  nb204nf_jyM(%rsp),%xmm1
        subpd  nb204nf_jzM(%rsp),%xmm2
        subpd  nb204nf_jxH1(%rsp),%xmm3
        subpd  nb204nf_jyH1(%rsp),%xmm4
        subpd  nb204nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqMM(%rsp)
        movapd %xmm3,nb204nf_rsqMH1(%rsp)

        movapd nb204nf_ixM(%rsp),%xmm0
        movapd nb204nf_iyM(%rsp),%xmm1
        movapd nb204nf_izM(%rsp),%xmm2
        movapd nb204nf_ixH1(%rsp),%xmm3
        movapd nb204nf_iyH1(%rsp),%xmm4
        movapd nb204nf_izH1(%rsp),%xmm5
        subpd  nb204nf_jxH2(%rsp),%xmm0
        subpd  nb204nf_jyH2(%rsp),%xmm1
        subpd  nb204nf_jzH2(%rsp),%xmm2
        subpd  nb204nf_jxM(%rsp),%xmm3
        subpd  nb204nf_jyM(%rsp),%xmm4
        subpd  nb204nf_jzM(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqMH2(%rsp)
        movapd %xmm3,nb204nf_rsqH1M(%rsp)

        movapd nb204nf_ixH1(%rsp),%xmm0
        movapd nb204nf_iyH1(%rsp),%xmm1
        movapd nb204nf_izH1(%rsp),%xmm2
        movapd nb204nf_ixH1(%rsp),%xmm3
        movapd nb204nf_iyH1(%rsp),%xmm4
        movapd nb204nf_izH1(%rsp),%xmm5
        subpd  nb204nf_jxH1(%rsp),%xmm0
        subpd  nb204nf_jyH1(%rsp),%xmm1
        subpd  nb204nf_jzH1(%rsp),%xmm2
        subpd  nb204nf_jxH2(%rsp),%xmm3
        subpd  nb204nf_jyH2(%rsp),%xmm4
        subpd  nb204nf_jzH2(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqH1H1(%rsp)
        movapd %xmm3,nb204nf_rsqH1H2(%rsp)

        movapd nb204nf_ixH2(%rsp),%xmm0
        movapd nb204nf_iyH2(%rsp),%xmm1
        movapd nb204nf_izH2(%rsp),%xmm2
        movapd nb204nf_ixH2(%rsp),%xmm3
        movapd nb204nf_iyH2(%rsp),%xmm4
        movapd nb204nf_izH2(%rsp),%xmm5
        subpd  nb204nf_jxM(%rsp),%xmm0
        subpd  nb204nf_jyM(%rsp),%xmm1
        subpd  nb204nf_jzM(%rsp),%xmm2
        subpd  nb204nf_jxH1(%rsp),%xmm3
        subpd  nb204nf_jyH1(%rsp),%xmm4
        subpd  nb204nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb204nf_rsqH2M(%rsp)
        movapd %xmm4,nb204nf_rsqH2H1(%rsp)

        movapd nb204nf_ixH2(%rsp),%xmm0
        movapd nb204nf_iyH2(%rsp),%xmm1
        movapd nb204nf_izH2(%rsp),%xmm2
        subpd  nb204nf_jxH2(%rsp),%xmm0
        subpd  nb204nf_jyH2(%rsp),%xmm1
        subpd  nb204nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb204nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb204nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH2H2(%rsp)
        movapd %xmm5,nb204nf_rinvH2H1(%rsp)

        movapd nb204nf_rsqMM(%rsp),%xmm0
        movapd nb204nf_rsqMH1(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb204nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb204nf_rinvMM(%rsp)
        movapd %xmm5,nb204nf_rinvMH1(%rsp)

        movapd nb204nf_rsqMH2(%rsp),%xmm0
        movapd nb204nf_rsqH1M(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb204nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvMH2(%rsp)
        movapd %xmm5,nb204nf_rinvH1M(%rsp)

        movapd nb204nf_rsqH1H1(%rsp),%xmm0
        movapd nb204nf_rsqH1H2(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb204nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH1H1(%rsp)
        movapd %xmm5,nb204nf_rinvH1H2(%rsp)

        movapd nb204nf_rsqH2M(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb204nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb204nf_krf(%rsp),%xmm6
        mulpd  nb204nf_rsqMM(%rsp),%xmm6        ## xmm5=krsq 
        addpd  nb204nf_rinvMM(%rsp),%xmm6       ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm6

        mulpd  nb204nf_qqMM(%rsp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  nb204nf_vctot(%rsp),%xmm6   ## local vctot summation variable 

        ## M-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulpd  nb204nf_rsqMH1(%rsp),%xmm5       ## xmm5=krsq 
        addpd  nb204nf_rinvMH1(%rsp),%xmm5      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm5

        mulpd  nb204nf_qqMH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## M-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulpd  nb204nf_rsqMH2(%rsp),%xmm7       ## xmm5=krsq 
        addpd  nb204nf_rinvMH2(%rsp),%xmm7      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm7

        mulpd  nb204nf_qqMH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-M interaction 
        movapd nb204nf_krf(%rsp),%xmm4
        mulpd  nb204nf_rsqH1M(%rsp),%xmm4       ## xmm5=krsq 
        addpd  nb204nf_rinvH1M(%rsp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm4

        mulpd  nb204nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulpd  nb204nf_rsqH1H1(%rsp),%xmm5      ## xmm5=krsq 
        addpd  nb204nf_rinvH1H1(%rsp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm5

        mulpd  nb204nf_qqHH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulpd  nb204nf_rsqH1H2(%rsp),%xmm7      ## xmm5=krsq 
        addpd  nb204nf_rinvH1H2(%rsp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm7

        mulpd  nb204nf_qqHH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-M interaction 
        movapd nb204nf_krf(%rsp),%xmm4
        mulpd  nb204nf_rsqH2M(%rsp),%xmm4       ## xmm5=krsq 
        addpd  nb204nf_rinvH2M(%rsp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm4

        mulpd  nb204nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulpd  nb204nf_rsqH2H1(%rsp),%xmm5      ## xmm5=krsq 
        addpd  nb204nf_rinvH2H1(%rsp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm5

        mulpd  nb204nf_qqHH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulpd  nb204nf_rsqH2H2(%rsp),%xmm7      ## xmm5=krsq 
        addpd  nb204nf_rinvH2H2(%rsp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%rsp),%xmm7

        mulpd  nb204nf_qqHH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 
        movapd %xmm6,nb204nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb204nf_innerk(%rsp)
        jl    _nb_kernel204nf_x86_64_sse2.nb204nf_checksingle
        jmp   _nb_kernel204nf_x86_64_sse2.nb204nf_unroll_loop
_nb_kernel204nf_x86_64_sse2.nb204nf_checksingle: 
        movl  nb204nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel204nf_x86_64_sse2.nb204nf_dosingle
        jmp   _nb_kernel204nf_x86_64_sse2.nb204nf_updateouterdata
_nb_kernel204nf_x86_64_sse2.nb204nf_dosingle: 
        movq  nb204nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb204nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coordinates to local temp variables 
        movlpd 24(%rsi,%rax,8),%xmm2
        movlpd 32(%rsi,%rax,8),%xmm3
        movlpd 40(%rsi,%rax,8),%xmm4
        movlpd 48(%rsi,%rax,8),%xmm5
        movlpd 56(%rsi,%rax,8),%xmm6
        movlpd 64(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb204nf_jxH1(%rsp)
        movapd  %xmm3,nb204nf_jyH1(%rsp)
        movapd  %xmm4,nb204nf_jzH1(%rsp)
        movapd  %xmm5,nb204nf_jxH2(%rsp)
        movapd  %xmm6,nb204nf_jyH2(%rsp)
        movapd  %xmm7,nb204nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb204nf_jxM(%rsp)
        movapd  %xmm3,nb204nf_jyM(%rsp)
        movapd  %xmm4,nb204nf_jzM(%rsp)

        movapd nb204nf_ixM(%rsp),%xmm0
        movapd nb204nf_iyM(%rsp),%xmm1
        movapd nb204nf_izM(%rsp),%xmm2
        movapd nb204nf_ixM(%rsp),%xmm3
        movapd nb204nf_iyM(%rsp),%xmm4
        movapd nb204nf_izM(%rsp),%xmm5
        subsd  nb204nf_jxM(%rsp),%xmm0
        subsd  nb204nf_jyM(%rsp),%xmm1
        subsd  nb204nf_jzM(%rsp),%xmm2
        subsd  nb204nf_jxH1(%rsp),%xmm3
        subsd  nb204nf_jyH1(%rsp),%xmm4
        subsd  nb204nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqMM(%rsp)
        movapd %xmm3,nb204nf_rsqMH1(%rsp)

        movapd nb204nf_ixM(%rsp),%xmm0
        movapd nb204nf_iyM(%rsp),%xmm1
        movapd nb204nf_izM(%rsp),%xmm2
        movapd nb204nf_ixH1(%rsp),%xmm3
        movapd nb204nf_iyH1(%rsp),%xmm4
        movapd nb204nf_izH1(%rsp),%xmm5
        subsd  nb204nf_jxH2(%rsp),%xmm0
        subsd  nb204nf_jyH2(%rsp),%xmm1
        subsd  nb204nf_jzH2(%rsp),%xmm2
        subsd  nb204nf_jxM(%rsp),%xmm3
        subsd  nb204nf_jyM(%rsp),%xmm4
        subsd  nb204nf_jzM(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqMH2(%rsp)
        movapd %xmm3,nb204nf_rsqH1M(%rsp)

        movapd nb204nf_ixH1(%rsp),%xmm0
        movapd nb204nf_iyH1(%rsp),%xmm1
        movapd nb204nf_izH1(%rsp),%xmm2
        movapd nb204nf_ixH1(%rsp),%xmm3
        movapd nb204nf_iyH1(%rsp),%xmm4
        movapd nb204nf_izH1(%rsp),%xmm5
        subsd  nb204nf_jxH1(%rsp),%xmm0
        subsd  nb204nf_jyH1(%rsp),%xmm1
        subsd  nb204nf_jzH1(%rsp),%xmm2
        subsd  nb204nf_jxH2(%rsp),%xmm3
        subsd  nb204nf_jyH2(%rsp),%xmm4
        subsd  nb204nf_jzH2(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204nf_rsqH1H1(%rsp)
        movapd %xmm3,nb204nf_rsqH1H2(%rsp)

        movapd nb204nf_ixH2(%rsp),%xmm0
        movapd nb204nf_iyH2(%rsp),%xmm1
        movapd nb204nf_izH2(%rsp),%xmm2
        movapd nb204nf_ixH2(%rsp),%xmm3
        movapd nb204nf_iyH2(%rsp),%xmm4
        movapd nb204nf_izH2(%rsp),%xmm5
        subsd  nb204nf_jxM(%rsp),%xmm0
        subsd  nb204nf_jyM(%rsp),%xmm1
        subsd  nb204nf_jzM(%rsp),%xmm2
        subsd  nb204nf_jxH1(%rsp),%xmm3
        subsd  nb204nf_jyH1(%rsp),%xmm4
        subsd  nb204nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb204nf_rsqH2M(%rsp)
        movapd %xmm4,nb204nf_rsqH2H1(%rsp)

        movapd nb204nf_ixH2(%rsp),%xmm0
        movapd nb204nf_iyH2(%rsp),%xmm1
        movapd nb204nf_izH2(%rsp),%xmm2
        subsd  nb204nf_jxH2(%rsp),%xmm0
        subsd  nb204nf_jyH2(%rsp),%xmm1
        subsd  nb204nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb204nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb204nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH2H2(%rsp)
        movapd %xmm5,nb204nf_rinvH2H1(%rsp)

        movapd nb204nf_rsqMM(%rsp),%xmm0
        movapd nb204nf_rsqMH1(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb204nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb204nf_rinvMM(%rsp)
        movapd %xmm5,nb204nf_rinvMH1(%rsp)

        movapd nb204nf_rsqMH2(%rsp),%xmm0
        movapd nb204nf_rsqH1M(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb204nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvMH2(%rsp)
        movapd %xmm5,nb204nf_rinvH1M(%rsp)

        movapd nb204nf_rsqH1H1(%rsp),%xmm0
        movapd nb204nf_rsqH1H2(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb204nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH1H1(%rsp)
        movapd %xmm5,nb204nf_rinvH1H2(%rsp)

        movapd nb204nf_rsqH2M(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb204nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb204nf_krf(%rsp),%xmm6
        mulsd  nb204nf_rsqMM(%rsp),%xmm6        ## xmm5=krsq 
        addsd  nb204nf_rinvMM(%rsp),%xmm6       ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm6

        mulsd  nb204nf_qqMM(%rsp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  nb204nf_vctot(%rsp),%xmm6   ## local vctot summation variable 

        ## M-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulsd  nb204nf_rsqMH1(%rsp),%xmm5       ## xmm5=krsq 
        addsd  nb204nf_rinvMH1(%rsp),%xmm5      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm5

        mulsd  nb204nf_qqMH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## M-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulsd  nb204nf_rsqMH2(%rsp),%xmm7       ## xmm5=krsq 
        addsd  nb204nf_rinvMH2(%rsp),%xmm7      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm7

        mulsd  nb204nf_qqMH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-M interaction 
        movapd nb204nf_krf(%rsp),%xmm4
        mulsd  nb204nf_rsqH1M(%rsp),%xmm4       ## xmm5=krsq 
        addsd  nb204nf_rinvH1M(%rsp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm4

        mulsd  nb204nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulsd  nb204nf_rsqH1H1(%rsp),%xmm5      ## xmm5=krsq 
        addsd  nb204nf_rinvH1H1(%rsp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm5

        mulsd  nb204nf_qqHH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulsd  nb204nf_rsqH1H2(%rsp),%xmm7      ## xmm5=krsq 
        addsd  nb204nf_rinvH1H2(%rsp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm7

        mulsd  nb204nf_qqHH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-M interaction 
        movapd nb204nf_krf(%rsp),%xmm4
        mulsd  nb204nf_rsqH2M(%rsp),%xmm4       ## xmm5=krsq 
        addsd  nb204nf_rinvH2M(%rsp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm4

        mulsd  nb204nf_qqMH(%rsp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb204nf_krf(%rsp),%xmm5
        mulsd  nb204nf_rsqH2H1(%rsp),%xmm5      ## xmm5=krsq 
        addsd  nb204nf_rinvH2H1(%rsp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm5

        mulsd  nb204nf_qqHH(%rsp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb204nf_krf(%rsp),%xmm7
        mulsd  nb204nf_rsqH2H2(%rsp),%xmm7      ## xmm5=krsq 
        addsd  nb204nf_rinvH2H2(%rsp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%rsp),%xmm7

        mulsd  nb204nf_qqHH(%rsp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 
        movlpd %xmm6,nb204nf_vctot(%rsp)

_nb_kernel204nf_x86_64_sse2.nb204nf_updateouterdata: 
        ## get n from stack
        movl nb204nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb204nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb204nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb204nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb204nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204nf_x86_64_sse2.nb204nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204nf_n(%rsp)
        jmp _nb_kernel204nf_x86_64_sse2.nb204nf_outer
_nb_kernel204nf_x86_64_sse2.nb204nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204nf_x86_64_sse2.nb204nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204nf_x86_64_sse2.nb204nf_threadloop
_nb_kernel204nf_x86_64_sse2.nb204nf_end: 
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


