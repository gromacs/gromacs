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






.globl nb_kernel304_x86_64_sse2
.globl _nb_kernel304_x86_64_sse2
nb_kernel304_x86_64_sse2:       
_nb_kernel304_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb304_fshift, 16
.set nb304_gid, 24
.set nb304_pos, 32
.set nb304_faction, 40
.set nb304_charge, 48
.set nb304_p_facel, 56
.set nb304_argkrf, 64
.set nb304_argcrf, 72
.set nb304_Vc, 80
.set nb304_type, 88
.set nb304_p_ntype, 96
.set nb304_vdwparam, 104
.set nb304_Vvdw, 112
.set nb304_p_tabscale, 120
.set nb304_VFtab, 128
.set nb304_invsqrta, 136
.set nb304_dvda, 144
.set nb304_p_gbtabscale, 152
.set nb304_GBtab, 160
.set nb304_p_nthreads, 168
.set nb304_count, 176
.set nb304_mtx, 184
.set nb304_outeriter, 192
.set nb304_inneriter, 200
.set nb304_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb304_ixM, 0
.set nb304_iyM, 16
.set nb304_izM, 32
.set nb304_ixH1, 48
.set nb304_iyH1, 64
.set nb304_izH1, 80
.set nb304_ixH2, 96
.set nb304_iyH2, 112
.set nb304_izH2, 128
.set nb304_jxM, 144
.set nb304_jyM, 160
.set nb304_jzM, 176
.set nb304_jxH1, 192
.set nb304_jyH1, 208
.set nb304_jzH1, 224
.set nb304_jxH2, 240
.set nb304_jyH2, 256
.set nb304_jzH2, 272
.set nb304_dxMM, 288
.set nb304_dyMM, 304
.set nb304_dzMM, 320
.set nb304_dxMH1, 336
.set nb304_dyMH1, 352
.set nb304_dzMH1, 368
.set nb304_dxMH2, 384
.set nb304_dyMH2, 400
.set nb304_dzMH2, 416
.set nb304_dxH1M, 432
.set nb304_dyH1M, 448
.set nb304_dzH1M, 464
.set nb304_dxH1H1, 480
.set nb304_dyH1H1, 496
.set nb304_dzH1H1, 512
.set nb304_dxH1H2, 528
.set nb304_dyH1H2, 544
.set nb304_dzH1H2, 560
.set nb304_dxH2M, 576
.set nb304_dyH2M, 592
.set nb304_dzH2M, 608
.set nb304_dxH2H1, 624
.set nb304_dyH2H1, 640
.set nb304_dzH2H1, 656
.set nb304_dxH2H2, 672
.set nb304_dyH2H2, 688
.set nb304_dzH2H2, 704
.set nb304_qqMM, 720
.set nb304_qqMH, 736
.set nb304_qqHH, 752
.set nb304_two, 768
.set nb304_tsc, 784
.set nb304_vctot, 800
.set nb304_fixM, 816
.set nb304_fiyM, 832
.set nb304_fizM, 848
.set nb304_fixH1, 864
.set nb304_fiyH1, 880
.set nb304_fizH1, 896
.set nb304_fixH2, 912
.set nb304_fiyH2, 928
.set nb304_fizH2, 944
.set nb304_epsH1, 960
.set nb304_epsH2, 976
.set nb304_epsM, 992
.set nb304_fjxH1, 1008
.set nb304_fjyH1, 1024
.set nb304_fjzH1, 1040
.set nb304_fjxH2, 1056
.set nb304_fjyH2, 1072
.set nb304_fjzH2, 1088
.set nb304_half, 1104
.set nb304_three, 1120
.set nb304_rsqMM, 1136
.set nb304_rsqMH1, 1152
.set nb304_rsqMH2, 1168
.set nb304_rsqH1M, 1184
.set nb304_rsqH1H1, 1200
.set nb304_rsqH1H2, 1216
.set nb304_rsqH2M, 1232
.set nb304_rsqH2H1, 1248
.set nb304_rsqH2H2, 1264
.set nb304_rinvMM, 1280
.set nb304_rinvMH1, 1296
.set nb304_rinvMH2, 1312
.set nb304_rinvH1M, 1328
.set nb304_rinvH1H1, 1344
.set nb304_rinvH1H2, 1360
.set nb304_rinvH2M, 1376
.set nb304_rinvH2H1, 1392
.set nb304_rinvH2H2, 1408
.set nb304_is3, 1424
.set nb304_ii3, 1428
.set nb304_nri, 1432
.set nb304_iinr, 1440
.set nb304_jindex, 1448
.set nb304_jjnr, 1456
.set nb304_shift, 1464
.set nb304_shiftvec, 1472
.set nb304_facel, 1480
.set nb304_innerjjnr, 1488
.set nb304_innerk, 1496
.set nb304_n, 1500
.set nb304_nn1, 1504
.set nb304_nouter, 1508
.set nb304_ninner, 1512
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1528,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb304_nouter(%rsp)
        movl %eax,nb304_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb304_nri(%rsp)
        movq %rsi,nb304_iinr(%rsp)
        movq %rdx,nb304_jindex(%rsp)
        movq %rcx,nb304_jjnr(%rsp)
        movq %r8,nb304_shift(%rsp)
        movq %r9,nb304_shiftvec(%rsp)
        movq nb304_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb304_facel(%rsp)

        movq nb304_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb304_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb304_half(%rsp)
        movl %ebx,nb304_half+4(%rsp)
        movsd nb304_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb304_half(%rsp)
        movapd %xmm2,nb304_two(%rsp)
        movapd %xmm3,nb304_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb304_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb304_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb304_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb304_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304_qqMM(%rsp)
        movapd %xmm4,nb304_qqMH(%rsp)
        movapd %xmm5,nb304_qqHH(%rsp)

_nb_kernel304_x86_64_sse2.nb304_threadloop: 
        movq  nb304_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel304_x86_64_sse2.nb304_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304_x86_64_sse2.nb304_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304_n(%rsp)
        movl %ebx,nb304_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304_x86_64_sse2.nb304_outerstart
        jmp _nb_kernel304_x86_64_sse2.nb304_end

_nb_kernel304_x86_64_sse2.nb304_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304_nouter(%rsp),%ebx
        movl %ebx,nb304_nouter(%rsp)

_nb_kernel304_x86_64_sse2.nb304_outer: 
        movq  nb304_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb304_is3(%rsp)      ## store is3 

        movq  nb304_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb304_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb304_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb304_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304_ixH1(%rsp)
        movapd %xmm4,nb304_iyH1(%rsp)
        movapd %xmm5,nb304_izH1(%rsp)

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
        movapd %xmm0,nb304_ixH2(%rsp)
        movapd %xmm1,nb304_iyH2(%rsp)
        movapd %xmm2,nb304_izH2(%rsp)
        movapd %xmm3,nb304_ixM(%rsp)
        movapd %xmm4,nb304_iyM(%rsp)
        movapd %xmm5,nb304_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb304_vctot(%rsp)
        movapd %xmm4,nb304_fixM(%rsp)
        movapd %xmm4,nb304_fiyM(%rsp)
        movapd %xmm4,nb304_fizM(%rsp)
        movapd %xmm4,nb304_fixH1(%rsp)
        movapd %xmm4,nb304_fiyH1(%rsp)
        movapd %xmm4,nb304_fizH1(%rsp)
        movapd %xmm4,nb304_fixH2(%rsp)
        movapd %xmm4,nb304_fiyH2(%rsp)
        movapd %xmm4,nb304_fizH2(%rsp)

        movq  nb304_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb304_pos(%rbp),%rsi
        movq  nb304_faction(%rbp),%rdi
        movq  nb304_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb304_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb304_ninner(%rsp),%ecx
        movl  %ecx,nb304_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb304_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel304_x86_64_sse2.nb304_unroll_loop
        jmp   _nb_kernel304_x86_64_sse2.nb304_checksingle
_nb_kernel304_x86_64_sse2.nb304_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb304_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb304_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb304_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb304_ixH1(%rsp),%xmm0
    subpd nb304_iyH1(%rsp),%xmm1
    subpd nb304_izH1(%rsp),%xmm2
    subpd nb304_ixH2(%rsp),%xmm3
    subpd nb304_iyH2(%rsp),%xmm4
    subpd nb304_izH2(%rsp),%xmm5
    subpd nb304_ixM(%rsp),%xmm6
    subpd nb304_iyM(%rsp),%xmm7
    subpd nb304_izM(%rsp),%xmm8

        movapd %xmm0,nb304_dxH1H1(%rsp)
        movapd %xmm1,nb304_dyH1H1(%rsp)
        movapd %xmm2,nb304_dzH1H1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH2H1(%rsp)
        movapd %xmm4,nb304_dyH2H1(%rsp)
        movapd %xmm5,nb304_dzH2H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb304_dxMH1(%rsp)
        movapd %xmm7,nb304_dyMH1(%rsp)
        movapd %xmm8,nb304_dzMH1(%rsp)
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

        movapd  nb304_three(%rsp),%xmm9
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

        movapd  nb304_half(%rsp),%xmm15
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

        movapd  nb304_three(%rsp),%xmm1
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

        movapd  nb304_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H1 
        mulpd   %xmm15,%xmm10 ##   rinvH2H1
    mulpd   %xmm15,%xmm11 ##   rinvMH1

        movapd  %xmm9,nb304_rinvH1H1(%rsp)
        movapd  %xmm10,nb304_rinvH2H1(%rsp)
        movapd  %xmm11,nb304_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb304_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb304_epsH1(%rsp)
    movapd    %xmm3,nb304_epsH2(%rsp)
    movapd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb304_epsH1(%rsp),%xmm12
    movapd nb304_epsH2(%rsp),%xmm13
    movapd nb304_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb304_qqHH(%rsp),%xmm12
    movapd nb304_qqMH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb304_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb304_vctot(%rsp)

    movapd nb304_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb304_rinvH1H1(%rsp),%xmm4
    mulpd nb304_rinvH2H1(%rsp),%xmm8
    mulpd nb304_rinvMH1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb304_dxH1H1(%rsp),%xmm3
        mulpd nb304_dyH1H1(%rsp),%xmm4
        mulpd nb304_dzH1H1(%rsp),%xmm5
        mulpd nb304_dxH2H1(%rsp),%xmm7
        mulpd nb304_dyH2H1(%rsp),%xmm8
        mulpd nb304_dzH2H1(%rsp),%xmm9
        mulpd nb304_dxMH1(%rsp),%xmm10
        mulpd nb304_dyMH1(%rsp),%xmm11
        mulpd nb304_dzMH1(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb304_fixH1(%rsp),%xmm3
    addpd nb304_fiyH1(%rsp),%xmm4
    addpd nb304_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb304_fixH2(%rsp),%xmm7
    addpd nb304_fiyH2(%rsp),%xmm8
    addpd nb304_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb304_fixM(%rsp),%xmm10
    addpd nb304_fiyM(%rsp),%xmm11
    addpd nb304_fizM(%rsp),%xmm12

    movapd %xmm3,nb304_fixH1(%rsp)
    movapd %xmm4,nb304_fiyH1(%rsp)
    movapd %xmm5,nb304_fizH1(%rsp)
    movapd %xmm7,nb304_fixH2(%rsp)
    movapd %xmm8,nb304_fiyH2(%rsp)
    movapd %xmm9,nb304_fizH2(%rsp)
    movapd %xmm10,nb304_fixM(%rsp)
    movapd %xmm11,nb304_fiyM(%rsp)
    movapd %xmm12,nb304_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb304_pos(%rbp),%rsi
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

    subpd nb304_ixH1(%rsp),%xmm0
    subpd nb304_iyH1(%rsp),%xmm1
    subpd nb304_izH1(%rsp),%xmm2
    subpd nb304_ixH2(%rsp),%xmm3
    subpd nb304_iyH2(%rsp),%xmm4
    subpd nb304_izH2(%rsp),%xmm5
    subpd nb304_ixM(%rsp),%xmm6
    subpd nb304_iyM(%rsp),%xmm7
    subpd nb304_izM(%rsp),%xmm8

        movapd %xmm0,nb304_dxH1H2(%rsp)
        movapd %xmm1,nb304_dyH1H2(%rsp)
        movapd %xmm2,nb304_dzH1H2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH2H2(%rsp)
        movapd %xmm4,nb304_dyH2H2(%rsp)
        movapd %xmm5,nb304_dzH2H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb304_dxMH2(%rsp)
        movapd %xmm7,nb304_dyMH2(%rsp)
        movapd %xmm8,nb304_dzMH2(%rsp)
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

        movapd  nb304_three(%rsp),%xmm9
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

        movapd  nb304_half(%rsp),%xmm15
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

        movapd  nb304_three(%rsp),%xmm1
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

        movapd  nb304_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H2
        mulpd   %xmm15,%xmm10 ##   rinvH2H2
    mulpd   %xmm15,%xmm11 ##   rinvMH2

        movapd  %xmm9,nb304_rinvH1H2(%rsp)
        movapd  %xmm10,nb304_rinvH2H2(%rsp)
        movapd  %xmm11,nb304_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb304_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb304_epsH1(%rsp)
    movapd    %xmm3,nb304_epsH2(%rsp)
    movapd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb304_epsH1(%rsp),%xmm12
    movapd nb304_epsH2(%rsp),%xmm13
    movapd nb304_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb304_qqHH(%rsp),%xmm12
    movapd nb304_qqMH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb304_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb304_vctot(%rsp)

    movapd nb304_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb304_rinvH1H2(%rsp),%xmm4
    mulpd nb304_rinvH2H2(%rsp),%xmm8
    mulpd nb304_rinvMH2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb304_dxH1H2(%rsp),%xmm3
        mulpd nb304_dyH1H2(%rsp),%xmm4
        mulpd nb304_dzH1H2(%rsp),%xmm5
        mulpd nb304_dxH2H2(%rsp),%xmm7
        mulpd nb304_dyH2H2(%rsp),%xmm8
        mulpd nb304_dzH2H2(%rsp),%xmm9
        mulpd nb304_dxMH2(%rsp),%xmm10
        mulpd nb304_dyMH2(%rsp),%xmm11
        mulpd nb304_dzMH2(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb304_fixH1(%rsp),%xmm3
    addpd nb304_fiyH1(%rsp),%xmm4
    addpd nb304_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb304_fixH2(%rsp),%xmm7
    addpd nb304_fiyH2(%rsp),%xmm8
    addpd nb304_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb304_fixM(%rsp),%xmm10
    addpd nb304_fiyM(%rsp),%xmm11
    addpd nb304_fizM(%rsp),%xmm12

    movapd %xmm3,nb304_fixH1(%rsp)
    movapd %xmm4,nb304_fiyH1(%rsp)
    movapd %xmm5,nb304_fizH1(%rsp)
    movapd %xmm7,nb304_fixH2(%rsp)
    movapd %xmm8,nb304_fiyH2(%rsp)
    movapd %xmm9,nb304_fizH2(%rsp)
    movapd %xmm10,nb304_fixM(%rsp)
    movapd %xmm11,nb304_fiyM(%rsp)
    movapd %xmm12,nb304_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## move j M coordinates to local temp variables 
    movq nb304_pos(%rbp),%rsi
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

    subpd nb304_ixH1(%rsp),%xmm0
    subpd nb304_iyH1(%rsp),%xmm1
    subpd nb304_izH1(%rsp),%xmm2
    subpd nb304_ixH2(%rsp),%xmm3
    subpd nb304_iyH2(%rsp),%xmm4
    subpd nb304_izH2(%rsp),%xmm5
    subpd nb304_ixM(%rsp),%xmm6
    subpd nb304_iyM(%rsp),%xmm7
    subpd nb304_izM(%rsp),%xmm8

        movapd %xmm0,nb304_dxH1M(%rsp)
        movapd %xmm1,nb304_dyH1M(%rsp)
        movapd %xmm2,nb304_dzH1M(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH2M(%rsp)
        movapd %xmm4,nb304_dyH2M(%rsp)
        movapd %xmm5,nb304_dzH2M(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb304_dxMM(%rsp)
        movapd %xmm7,nb304_dyMM(%rsp)
        movapd %xmm8,nb304_dzMM(%rsp)
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

        movapd  nb304_three(%rsp),%xmm9
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

        movapd  nb304_half(%rsp),%xmm15
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

        movapd  nb304_three(%rsp),%xmm1
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

        movapd  nb304_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1M
        mulpd   %xmm15,%xmm10 ##   rinvH2M
    mulpd   %xmm15,%xmm11 ##   rinvMM


        movapd  %xmm9,nb304_rinvH1M(%rsp)
        movapd  %xmm10,nb304_rinvH2M(%rsp)
        movapd  %xmm11,nb304_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb304_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb304_epsH1(%rsp)
    movapd    %xmm3,nb304_epsH2(%rsp)
    movapd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb304_epsH1(%rsp),%xmm12
    movapd nb304_epsH2(%rsp),%xmm13
    movapd nb304_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb304_qqMH(%rsp),%xmm12
    movapd nb304_qqMM(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb304_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb304_vctot(%rsp)

    movapd nb304_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb304_rinvH1M(%rsp),%xmm4
    mulpd nb304_rinvH2M(%rsp),%xmm8
    mulpd nb304_rinvMM(%rsp),%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movlpd 72(%rdi,%rax,8),%xmm0
        movlpd 80(%rdi,%rax,8),%xmm1
        movlpd 88(%rdi,%rax,8),%xmm2
        movhpd 72(%rdi,%rbx,8),%xmm0
        movhpd 80(%rdi,%rbx,8),%xmm1
        movhpd 88(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb304_dxH1M(%rsp),%xmm3
        mulpd nb304_dyH1M(%rsp),%xmm4
        mulpd nb304_dzH1M(%rsp),%xmm5
        mulpd nb304_dxH2M(%rsp),%xmm7
        mulpd nb304_dyH2M(%rsp),%xmm8
        mulpd nb304_dzH2M(%rsp),%xmm9
        mulpd nb304_dxMM(%rsp),%xmm10
        mulpd nb304_dyMM(%rsp),%xmm11
        mulpd nb304_dzMM(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb304_fixH1(%rsp),%xmm3
    addpd nb304_fiyH1(%rsp),%xmm4
    addpd nb304_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb304_fixH2(%rsp),%xmm7
    addpd nb304_fiyH2(%rsp),%xmm8
    addpd nb304_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb304_fixM(%rsp),%xmm10
    addpd nb304_fiyM(%rsp),%xmm11
    addpd nb304_fizM(%rsp),%xmm12

    movapd %xmm3,nb304_fixH1(%rsp)
    movapd %xmm4,nb304_fiyH1(%rsp)
    movapd %xmm5,nb304_fizH1(%rsp)
    movapd %xmm7,nb304_fixH2(%rsp)
    movapd %xmm8,nb304_fiyH2(%rsp)
    movapd %xmm9,nb304_fizH2(%rsp)
    movapd %xmm10,nb304_fixM(%rsp)
    movapd %xmm11,nb304_fiyM(%rsp)
    movapd %xmm12,nb304_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movlpd %xmm0,72(%rdi,%rax,8)
        movlpd %xmm1,80(%rdi,%rax,8)
        movlpd %xmm2,88(%rdi,%rax,8)
        movhpd %xmm0,72(%rdi,%rbx,8)
        movhpd %xmm1,80(%rdi,%rbx,8)
        movhpd %xmm2,88(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb304_innerk(%rsp)
        jl    _nb_kernel304_x86_64_sse2.nb304_checksingle
        jmp   _nb_kernel304_x86_64_sse2.nb304_unroll_loop
_nb_kernel304_x86_64_sse2.nb304_checksingle: 
        movl  nb304_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel304_x86_64_sse2.nb304_dosingle
        jmp   _nb_kernel304_x86_64_sse2.nb304_updateouterdata
_nb_kernel304_x86_64_sse2.nb304_dosingle: 
        movq  nb304_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb304_pos(%rbp),%rsi
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

    subsd nb304_ixH1(%rsp),%xmm0
    subsd nb304_iyH1(%rsp),%xmm1
    subsd nb304_izH1(%rsp),%xmm2
    subsd nb304_ixH2(%rsp),%xmm3
    subsd nb304_iyH2(%rsp),%xmm4
    subsd nb304_izH2(%rsp),%xmm5
    subsd nb304_ixM(%rsp),%xmm6
    subsd nb304_iyM(%rsp),%xmm7
    subsd nb304_izM(%rsp),%xmm8

        movsd %xmm0,nb304_dxH1H1(%rsp)
        movsd %xmm1,nb304_dyH1H1(%rsp)
        movsd %xmm2,nb304_dzH1H1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb304_dxH2H1(%rsp)
        movsd %xmm4,nb304_dyH2H1(%rsp)
        movsd %xmm5,nb304_dzH2H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb304_dxMH1(%rsp)
        movsd %xmm7,nb304_dyMH1(%rsp)
        movsd %xmm8,nb304_dzMH1(%rsp)
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

        movsd  nb304_three(%rsp),%xmm9
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

        movsd  nb304_half(%rsp),%xmm15
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

        movsd  nb304_three(%rsp),%xmm1
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

        movsd  nb304_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H1 
        mulsd   %xmm15,%xmm10 ##   rinvH2H1
    mulsd   %xmm15,%xmm11 ##   rinvMH1

        movsd  %xmm9,nb304_rinvH1H1(%rsp)
        movsd  %xmm10,nb304_rinvH2H1(%rsp)
        movsd  %xmm11,nb304_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb304_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb304_epsH1(%rsp)
    movsd    %xmm3,nb304_epsH2(%rsp)
    movsd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb304_epsH1(%rsp),%xmm12
    movsd nb304_epsH2(%rsp),%xmm13
    movsd nb304_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb304_qqHH(%rsp),%xmm12
    movsd nb304_qqMH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb304_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb304_vctot(%rsp)

    movsd nb304_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb304_rinvH1H1(%rsp),%xmm4
    mulsd nb304_rinvH2H1(%rsp),%xmm8
    mulsd nb304_rinvMH1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb304_dxH1H1(%rsp),%xmm3
        mulsd nb304_dyH1H1(%rsp),%xmm4
        mulsd nb304_dzH1H1(%rsp),%xmm5
        mulsd nb304_dxH2H1(%rsp),%xmm7
        mulsd nb304_dyH2H1(%rsp),%xmm8
        mulsd nb304_dzH2H1(%rsp),%xmm9
        mulsd nb304_dxMH1(%rsp),%xmm10
        mulsd nb304_dyMH1(%rsp),%xmm11
        mulsd nb304_dzMH1(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb304_fixH1(%rsp),%xmm3
    addsd nb304_fiyH1(%rsp),%xmm4
    addsd nb304_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb304_fixH2(%rsp),%xmm7
    addsd nb304_fiyH2(%rsp),%xmm8
    addsd nb304_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb304_fixM(%rsp),%xmm10
    addsd nb304_fiyM(%rsp),%xmm11
    addsd nb304_fizM(%rsp),%xmm12

    movsd %xmm3,nb304_fixH1(%rsp)
    movsd %xmm4,nb304_fiyH1(%rsp)
    movsd %xmm5,nb304_fizH1(%rsp)
    movsd %xmm7,nb304_fixH2(%rsp)
    movsd %xmm8,nb304_fiyH2(%rsp)
    movsd %xmm9,nb304_fizH2(%rsp)
    movsd %xmm10,nb304_fixM(%rsp)
    movsd %xmm11,nb304_fiyM(%rsp)
    movsd %xmm12,nb304_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movq nb304_pos(%rbp),%rsi
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

    subsd nb304_ixH1(%rsp),%xmm0
    subsd nb304_iyH1(%rsp),%xmm1
    subsd nb304_izH1(%rsp),%xmm2
    subsd nb304_ixH2(%rsp),%xmm3
    subsd nb304_iyH2(%rsp),%xmm4
    subsd nb304_izH2(%rsp),%xmm5
    subsd nb304_ixM(%rsp),%xmm6
    subsd nb304_iyM(%rsp),%xmm7
    subsd nb304_izM(%rsp),%xmm8

        movsd %xmm0,nb304_dxH1H2(%rsp)
        movsd %xmm1,nb304_dyH1H2(%rsp)
        movsd %xmm2,nb304_dzH1H2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb304_dxH2H2(%rsp)
        movsd %xmm4,nb304_dyH2H2(%rsp)
        movsd %xmm5,nb304_dzH2H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb304_dxMH2(%rsp)
        movsd %xmm7,nb304_dyMH2(%rsp)
        movsd %xmm8,nb304_dzMH2(%rsp)
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

        movsd  nb304_three(%rsp),%xmm9
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

        movsd  nb304_half(%rsp),%xmm15
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

        movsd  nb304_three(%rsp),%xmm1
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

        movsd  nb304_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H2
        mulsd   %xmm15,%xmm10 ##   rinvH2H2
    mulsd   %xmm15,%xmm11 ##   rinvMH2

        movsd  %xmm9,nb304_rinvH1H2(%rsp)
        movsd  %xmm10,nb304_rinvH2H2(%rsp)
        movsd  %xmm11,nb304_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb304_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb304_epsH1(%rsp)
    movsd    %xmm3,nb304_epsH2(%rsp)
    movsd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb304_epsH1(%rsp),%xmm12
    movsd nb304_epsH2(%rsp),%xmm13
    movsd nb304_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb304_qqHH(%rsp),%xmm12
    movsd nb304_qqMH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb304_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb304_vctot(%rsp)

    movsd nb304_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb304_rinvH1H2(%rsp),%xmm4
    mulsd nb304_rinvH2H2(%rsp),%xmm8
    mulsd nb304_rinvMH2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb304_dxH1H2(%rsp),%xmm3
        mulsd nb304_dyH1H2(%rsp),%xmm4
        mulsd nb304_dzH1H2(%rsp),%xmm5
        mulsd nb304_dxH2H2(%rsp),%xmm7
        mulsd nb304_dyH2H2(%rsp),%xmm8
        mulsd nb304_dzH2H2(%rsp),%xmm9
        mulsd nb304_dxMH2(%rsp),%xmm10
        mulsd nb304_dyMH2(%rsp),%xmm11
        mulsd nb304_dzMH2(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb304_fixH1(%rsp),%xmm3
    addsd nb304_fiyH1(%rsp),%xmm4
    addsd nb304_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb304_fixH2(%rsp),%xmm7
    addsd nb304_fiyH2(%rsp),%xmm8
    addsd nb304_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb304_fixM(%rsp),%xmm10
    addsd nb304_fiyM(%rsp),%xmm11
    addsd nb304_fizM(%rsp),%xmm12

    movsd %xmm3,nb304_fixH1(%rsp)
    movsd %xmm4,nb304_fiyH1(%rsp)
    movsd %xmm5,nb304_fizH1(%rsp)
    movsd %xmm7,nb304_fixH2(%rsp)
    movsd %xmm8,nb304_fiyH2(%rsp)
    movsd %xmm9,nb304_fizH2(%rsp)
    movsd %xmm10,nb304_fixM(%rsp)
    movsd %xmm11,nb304_fiyM(%rsp)
    movsd %xmm12,nb304_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

        ## move j M coordinates to local temp variables 
    movq nb304_pos(%rbp),%rsi
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

    subsd nb304_ixH1(%rsp),%xmm0
    subsd nb304_iyH1(%rsp),%xmm1
    subsd nb304_izH1(%rsp),%xmm2
    subsd nb304_ixH2(%rsp),%xmm3
    subsd nb304_iyH2(%rsp),%xmm4
    subsd nb304_izH2(%rsp),%xmm5
    subsd nb304_ixM(%rsp),%xmm6
    subsd nb304_iyM(%rsp),%xmm7
    subsd nb304_izM(%rsp),%xmm8

        movsd %xmm0,nb304_dxH1M(%rsp)
        movsd %xmm1,nb304_dyH1M(%rsp)
        movsd %xmm2,nb304_dzH1M(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb304_dxH2M(%rsp)
        movsd %xmm4,nb304_dyH2M(%rsp)
        movsd %xmm5,nb304_dzH2M(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb304_dxMM(%rsp)
        movsd %xmm7,nb304_dyMM(%rsp)
        movsd %xmm8,nb304_dzMM(%rsp)
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

        movsd  nb304_three(%rsp),%xmm9
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

        movsd  nb304_half(%rsp),%xmm15
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

        movsd  nb304_three(%rsp),%xmm1
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

        movsd  nb304_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1M
        mulsd   %xmm15,%xmm10 ##   rinvH2M
    mulsd   %xmm15,%xmm11 ##   rinvMM

        movsd  %xmm9,nb304_rinvH1M(%rsp)
        movsd  %xmm10,nb304_rinvH2M(%rsp)
        movsd  %xmm11,nb304_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb304_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb304_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb304_epsH1(%rsp)
    movsd    %xmm3,nb304_epsH2(%rsp)
    movsd    %xmm6,nb304_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb304_epsH1(%rsp),%xmm12
    movsd nb304_epsH2(%rsp),%xmm13
    movsd nb304_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb304_qqMH(%rsp),%xmm12
    movsd nb304_qqMM(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb304_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb304_vctot(%rsp)

    movsd nb304_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb304_rinvH1M(%rsp),%xmm4
    mulsd nb304_rinvH2M(%rsp),%xmm8
    mulsd nb304_rinvMM(%rsp),%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb304_faction(%rbp),%rdi
        movsd 72(%rdi,%rax,8),%xmm0
        movsd 80(%rdi,%rax,8),%xmm1
        movsd 88(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb304_dxH1M(%rsp),%xmm3
        mulsd nb304_dyH1M(%rsp),%xmm4
        mulsd nb304_dzH1M(%rsp),%xmm5
        mulsd nb304_dxH2M(%rsp),%xmm7
        mulsd nb304_dyH2M(%rsp),%xmm8
        mulsd nb304_dzH2M(%rsp),%xmm9
        mulsd nb304_dxMM(%rsp),%xmm10
        mulsd nb304_dyMM(%rsp),%xmm11
        mulsd nb304_dzMM(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb304_fixH1(%rsp),%xmm3
    addsd nb304_fiyH1(%rsp),%xmm4
    addsd nb304_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb304_fixH2(%rsp),%xmm7
    addsd nb304_fiyH2(%rsp),%xmm8
    addsd nb304_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb304_fixM(%rsp),%xmm10
    addsd nb304_fiyM(%rsp),%xmm11
    addsd nb304_fizM(%rsp),%xmm12

    movsd %xmm3,nb304_fixH1(%rsp)
    movsd %xmm4,nb304_fiyH1(%rsp)
    movsd %xmm5,nb304_fizH1(%rsp)
    movsd %xmm7,nb304_fixH2(%rsp)
    movsd %xmm8,nb304_fiyH2(%rsp)
    movsd %xmm9,nb304_fizH2(%rsp)
    movsd %xmm10,nb304_fixM(%rsp)
    movsd %xmm11,nb304_fiyM(%rsp)
    movsd %xmm12,nb304_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movsd %xmm0,72(%rdi,%rax,8)
        movsd %xmm1,80(%rdi,%rax,8)
        movsd %xmm2,88(%rdi,%rax,8)

_nb_kernel304_x86_64_sse2.nb304_updateouterdata: 
        movl  nb304_ii3(%rsp),%ecx
        movq  nb304_faction(%rbp),%rdi
        movq  nb304_fshift(%rbp),%rsi
        movl  nb304_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb304_fixH1(%rsp),%xmm0
        movapd nb304_fiyH1(%rsp),%xmm1
        movapd nb304_fizH1(%rsp),%xmm2

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
        movapd nb304_fixH2(%rsp),%xmm0
        movapd nb304_fiyH2(%rsp),%xmm1
        movapd nb304_fizH2(%rsp),%xmm2

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
        movapd nb304_fixM(%rsp),%xmm0
        movapd nb304_fiyM(%rsp),%xmm1
        movapd nb304_fizM(%rsp),%xmm2

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
        movl nb304_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb304_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb304_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb304_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb304_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304_x86_64_sse2.nb304_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304_n(%rsp)
        jmp _nb_kernel304_x86_64_sse2.nb304_outer
_nb_kernel304_x86_64_sse2.nb304_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304_x86_64_sse2.nb304_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304_x86_64_sse2.nb304_threadloop
_nb_kernel304_x86_64_sse2.nb304_end: 
        movl nb304_nouter(%rsp),%eax
        movl nb304_ninner(%rsp),%ebx
        movq nb304_outeriter(%rbp),%rcx
        movq nb304_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1528,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel304nf_x86_64_sse2
.globl _nb_kernel304nf_x86_64_sse2
nb_kernel304nf_x86_64_sse2:     
_nb_kernel304nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb304nf_fshift, 16
.set nb304nf_gid, 24
.set nb304nf_pos, 32
.set nb304nf_faction, 40
.set nb304nf_charge, 48
.set nb304nf_p_facel, 56
.set nb304nf_argkrf, 64
.set nb304nf_argcrf, 72
.set nb304nf_Vc, 80
.set nb304nf_type, 88
.set nb304nf_p_ntype, 96
.set nb304nf_vdwparam, 104
.set nb304nf_Vvdw, 112
.set nb304nf_p_tabscale, 120
.set nb304nf_VFtab, 128
.set nb304nf_invsqrta, 136
.set nb304nf_dvda, 144
.set nb304nf_p_gbtabscale, 152
.set nb304nf_GBtab, 160
.set nb304nf_p_nthreads, 168
.set nb304nf_count, 176
.set nb304nf_mtx, 184
.set nb304nf_outeriter, 192
.set nb304nf_inneriter, 200
.set nb304nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb304nf_ixM, 0
.set nb304nf_iyM, 16
.set nb304nf_izM, 32
.set nb304nf_ixH1, 48
.set nb304nf_iyH1, 64
.set nb304nf_izH1, 80
.set nb304nf_ixH2, 96
.set nb304nf_iyH2, 112
.set nb304nf_izH2, 128
.set nb304nf_jxM, 144
.set nb304nf_jyM, 160
.set nb304nf_jzM, 176
.set nb304nf_jxH1, 192
.set nb304nf_jyH1, 208
.set nb304nf_jzH1, 224
.set nb304nf_jxH2, 240
.set nb304nf_jyH2, 256
.set nb304nf_jzH2, 272
.set nb304nf_qqMM, 288
.set nb304nf_qqMH, 304
.set nb304nf_qqHH, 320
.set nb304nf_tsc, 336
.set nb304nf_vctot, 352
.set nb304nf_half, 368
.set nb304nf_three, 384
.set nb304nf_rsqMM, 400
.set nb304nf_rsqMH1, 416
.set nb304nf_rsqMH2, 432
.set nb304nf_rsqH1M, 448
.set nb304nf_rsqH1H1, 464
.set nb304nf_rsqH1H2, 480
.set nb304nf_rsqH2M, 496
.set nb304nf_rsqH2H1, 512
.set nb304nf_rsqH2H2, 528
.set nb304nf_rinvMM, 544
.set nb304nf_rinvMH1, 560
.set nb304nf_rinvMH2, 576
.set nb304nf_rinvH1M, 592
.set nb304nf_rinvH1H1, 608
.set nb304nf_rinvH1H2, 624
.set nb304nf_rinvH2M, 640
.set nb304nf_rinvH2H1, 656
.set nb304nf_rinvH2H2, 672
.set nb304nf_is3, 688
.set nb304nf_ii3, 692
.set nb304nf_nri, 696
.set nb304nf_iinr, 704
.set nb304nf_jindex, 712
.set nb304nf_jjnr, 720
.set nb304nf_shift, 728
.set nb304nf_shiftvec, 736
.set nb304nf_facel, 744
.set nb304nf_innerjjnr, 752
.set nb304nf_innerk, 760
.set nb304nf_n, 764
.set nb304nf_nn1, 768
.set nb304nf_nouter, 772
.set nb304nf_ninner, 776
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $792,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb304nf_nouter(%rsp)
        movl %eax,nb304nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb304nf_nri(%rsp)
        movq %rsi,nb304nf_iinr(%rsp)
        movq %rdx,nb304nf_jindex(%rsp)
        movq %rcx,nb304nf_jjnr(%rsp)
        movq %r8,nb304nf_shift(%rsp)
        movq %r9,nb304nf_shiftvec(%rsp)
        movq nb304nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb304nf_facel(%rsp)

        movq nb304nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb304nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb304nf_half(%rsp)
        movl %ebx,nb304nf_half+4(%rsp)
        movsd nb304nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb304nf_half(%rsp)
        movapd %xmm3,nb304nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb304nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb304nf_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb304nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb304nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304nf_qqMM(%rsp)
        movapd %xmm4,nb304nf_qqMH(%rsp)
        movapd %xmm5,nb304nf_qqHH(%rsp)

_nb_kernel304nf_x86_64_sse2.nb304nf_threadloop: 
        movq  nb304nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel304nf_x86_64_sse2.nb304nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304nf_x86_64_sse2.nb304nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304nf_n(%rsp)
        movl %ebx,nb304nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304nf_x86_64_sse2.nb304nf_outerstart
        jmp _nb_kernel304nf_x86_64_sse2.nb304nf_end

_nb_kernel304nf_x86_64_sse2.nb304nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304nf_nouter(%rsp),%ebx
        movl %ebx,nb304nf_nouter(%rsp)

_nb_kernel304nf_x86_64_sse2.nb304nf_outer: 
        movq  nb304nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb304nf_is3(%rsp)            ## store is3 

        movq  nb304nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb304nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb304nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb304nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304nf_ixH1(%rsp)
        movapd %xmm4,nb304nf_iyH1(%rsp)
        movapd %xmm5,nb304nf_izH1(%rsp)

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
        movapd %xmm0,nb304nf_ixH2(%rsp)
        movapd %xmm1,nb304nf_iyH2(%rsp)
        movapd %xmm2,nb304nf_izH2(%rsp)
        movapd %xmm3,nb304nf_ixM(%rsp)
        movapd %xmm4,nb304nf_iyM(%rsp)
        movapd %xmm5,nb304nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb304nf_vctot(%rsp)

        movq  nb304nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb304nf_pos(%rbp),%rsi
        movq  nb304nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb304nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb304nf_ninner(%rsp),%ecx
        movl  %ecx,nb304nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb304nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel304nf_x86_64_sse2.nb304nf_unroll_loop
        jmp   _nb_kernel304nf_x86_64_sse2.nb304nf_checksingle
_nb_kernel304nf_x86_64_sse2.nb304nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb304nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb304nf_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb304nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd  %xmm2,nb304nf_jxH1(%rsp)
        movapd  %xmm3,nb304nf_jyH1(%rsp)
        movapd  %xmm4,nb304nf_jzH1(%rsp)
        movapd  %xmm5,nb304nf_jxH2(%rsp)
        movapd  %xmm6,nb304nf_jyH2(%rsp)
        movapd  %xmm7,nb304nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movhpd 72(%rsi,%rbx,8),%xmm2
        movhpd 80(%rsi,%rbx,8),%xmm3
        movhpd 88(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb304nf_jxM(%rsp)
        movapd  %xmm3,nb304nf_jyM(%rsp)
        movapd  %xmm4,nb304nf_jzM(%rsp)

        movapd nb304nf_ixM(%rsp),%xmm0
        movapd nb304nf_iyM(%rsp),%xmm1
        movapd nb304nf_izM(%rsp),%xmm2
        movapd nb304nf_ixM(%rsp),%xmm3
        movapd nb304nf_iyM(%rsp),%xmm4
        movapd nb304nf_izM(%rsp),%xmm5
        subpd  nb304nf_jxM(%rsp),%xmm0
        subpd  nb304nf_jyM(%rsp),%xmm1
        subpd  nb304nf_jzM(%rsp),%xmm2
        subpd  nb304nf_jxH1(%rsp),%xmm3
        subpd  nb304nf_jyH1(%rsp),%xmm4
        subpd  nb304nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMM(%rsp)
        movapd %xmm3,nb304nf_rsqMH1(%rsp)

        movapd nb304nf_ixM(%rsp),%xmm0
        movapd nb304nf_iyM(%rsp),%xmm1
        movapd nb304nf_izM(%rsp),%xmm2
        movapd nb304nf_ixH1(%rsp),%xmm3
        movapd nb304nf_iyH1(%rsp),%xmm4
        movapd nb304nf_izH1(%rsp),%xmm5
        subpd  nb304nf_jxH2(%rsp),%xmm0
        subpd  nb304nf_jyH2(%rsp),%xmm1
        subpd  nb304nf_jzH2(%rsp),%xmm2
        subpd  nb304nf_jxM(%rsp),%xmm3
        subpd  nb304nf_jyM(%rsp),%xmm4
        subpd  nb304nf_jzM(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMH2(%rsp)
        movapd %xmm3,nb304nf_rsqH1M(%rsp)

        movapd nb304nf_ixH1(%rsp),%xmm0
        movapd nb304nf_iyH1(%rsp),%xmm1
        movapd nb304nf_izH1(%rsp),%xmm2
        movapd nb304nf_ixH1(%rsp),%xmm3
        movapd nb304nf_iyH1(%rsp),%xmm4
        movapd nb304nf_izH1(%rsp),%xmm5
        subpd  nb304nf_jxH1(%rsp),%xmm0
        subpd  nb304nf_jyH1(%rsp),%xmm1
        subpd  nb304nf_jzH1(%rsp),%xmm2
        subpd  nb304nf_jxH2(%rsp),%xmm3
        subpd  nb304nf_jyH2(%rsp),%xmm4
        subpd  nb304nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH1H1(%rsp)
        movapd %xmm3,nb304nf_rsqH1H2(%rsp)

        movapd nb304nf_ixH2(%rsp),%xmm0
        movapd nb304nf_iyH2(%rsp),%xmm1
        movapd nb304nf_izH2(%rsp),%xmm2
        movapd nb304nf_ixH2(%rsp),%xmm3
        movapd nb304nf_iyH2(%rsp),%xmm4
        movapd nb304nf_izH2(%rsp),%xmm5
        subpd  nb304nf_jxM(%rsp),%xmm0
        subpd  nb304nf_jyM(%rsp),%xmm1
        subpd  nb304nf_jzM(%rsp),%xmm2
        subpd  nb304nf_jxH1(%rsp),%xmm3
        subpd  nb304nf_jyH1(%rsp),%xmm4
        subpd  nb304nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH2M(%rsp)
        movapd %xmm4,nb304nf_rsqH2H1(%rsp)

        movapd nb304nf_ixH2(%rsp),%xmm0
        movapd nb304nf_iyH2(%rsp),%xmm1
        movapd nb304nf_izH2(%rsp),%xmm2
        subpd  nb304nf_jxH2(%rsp),%xmm0
        subpd  nb304nf_jyH2(%rsp),%xmm1
        subpd  nb304nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb304nf_rsqH2H2(%rsp)

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
        movapd  nb304nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb304nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH2H2(%rsp)
        movapd %xmm5,nb304nf_rinvH2H1(%rsp)

        movapd nb304nf_rsqMM(%rsp),%xmm0
        movapd nb304nf_rsqMH1(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb304nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb304nf_rinvMM(%rsp)
        movapd %xmm5,nb304nf_rinvMH1(%rsp)

        movapd nb304nf_rsqMH2(%rsp),%xmm0
        movapd nb304nf_rsqH1M(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb304nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvMH2(%rsp)
        movapd %xmm5,nb304nf_rinvH1M(%rsp)

        movapd nb304nf_rsqH1H1(%rsp),%xmm0
        movapd nb304nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb304nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH1H1(%rsp)
        movapd %xmm5,nb304nf_rinvH1H2(%rsp)

        movapd nb304nf_rsqH2M(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb304nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb304nf_rinvMM(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMM(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## M-H1 interaction 
        movapd nb304nf_rinvMH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## M-H2 interaction  
        movapd nb304nf_rinvMH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H1-M interaction 
        movapd nb304nf_rinvH1M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H1-H1 interaction 
        movapd nb304nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H1-H2 interaction 
        movapd nb304nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H2-M interaction 
        movapd nb304nf_rinvH2M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb304nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb304nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb304nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb304nf_innerk(%rsp)
        jl    _nb_kernel304nf_x86_64_sse2.nb304nf_checksingle
        jmp   _nb_kernel304nf_x86_64_sse2.nb304nf_unroll_loop
_nb_kernel304nf_x86_64_sse2.nb304nf_checksingle: 
        movl  nb304nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel304nf_x86_64_sse2.nb304nf_dosingle
        jmp   _nb_kernel304nf_x86_64_sse2.nb304nf_updateouterdata
_nb_kernel304nf_x86_64_sse2.nb304nf_dosingle: 
        movq  nb304nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb304nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coordinates to local temp variables 
        movlpd 24(%rsi,%rax,8),%xmm2
        movlpd 32(%rsi,%rax,8),%xmm3
        movlpd 40(%rsi,%rax,8),%xmm4
        movlpd 48(%rsi,%rax,8),%xmm5
        movlpd 56(%rsi,%rax,8),%xmm6
        movlpd 64(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb304nf_jxH1(%rsp)
        movapd  %xmm3,nb304nf_jyH1(%rsp)
        movapd  %xmm4,nb304nf_jzH1(%rsp)
        movapd  %xmm5,nb304nf_jxH2(%rsp)
        movapd  %xmm6,nb304nf_jyH2(%rsp)
        movapd  %xmm7,nb304nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb304nf_jxM(%rsp)
        movapd  %xmm3,nb304nf_jyM(%rsp)
        movapd  %xmm4,nb304nf_jzM(%rsp)

        movapd nb304nf_ixM(%rsp),%xmm0
        movapd nb304nf_iyM(%rsp),%xmm1
        movapd nb304nf_izM(%rsp),%xmm2
        movapd nb304nf_ixM(%rsp),%xmm3
        movapd nb304nf_iyM(%rsp),%xmm4
        movapd nb304nf_izM(%rsp),%xmm5
        subsd  nb304nf_jxM(%rsp),%xmm0
        subsd  nb304nf_jyM(%rsp),%xmm1
        subsd  nb304nf_jzM(%rsp),%xmm2
        subsd  nb304nf_jxH1(%rsp),%xmm3
        subsd  nb304nf_jyH1(%rsp),%xmm4
        subsd  nb304nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMM(%rsp)
        movapd %xmm3,nb304nf_rsqMH1(%rsp)

        movapd nb304nf_ixM(%rsp),%xmm0
        movapd nb304nf_iyM(%rsp),%xmm1
        movapd nb304nf_izM(%rsp),%xmm2
        movapd nb304nf_ixH1(%rsp),%xmm3
        movapd nb304nf_iyH1(%rsp),%xmm4
        movapd nb304nf_izH1(%rsp),%xmm5
        subsd  nb304nf_jxH2(%rsp),%xmm0
        subsd  nb304nf_jyH2(%rsp),%xmm1
        subsd  nb304nf_jzH2(%rsp),%xmm2
        subsd  nb304nf_jxM(%rsp),%xmm3
        subsd  nb304nf_jyM(%rsp),%xmm4
        subsd  nb304nf_jzM(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMH2(%rsp)
        movapd %xmm3,nb304nf_rsqH1M(%rsp)

        movapd nb304nf_ixH1(%rsp),%xmm0
        movapd nb304nf_iyH1(%rsp),%xmm1
        movapd nb304nf_izH1(%rsp),%xmm2
        movapd nb304nf_ixH1(%rsp),%xmm3
        movapd nb304nf_iyH1(%rsp),%xmm4
        movapd nb304nf_izH1(%rsp),%xmm5
        subsd  nb304nf_jxH1(%rsp),%xmm0
        subsd  nb304nf_jyH1(%rsp),%xmm1
        subsd  nb304nf_jzH1(%rsp),%xmm2
        subsd  nb304nf_jxH2(%rsp),%xmm3
        subsd  nb304nf_jyH2(%rsp),%xmm4
        subsd  nb304nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH1H1(%rsp)
        movapd %xmm3,nb304nf_rsqH1H2(%rsp)

        movapd nb304nf_ixH2(%rsp),%xmm0
        movapd nb304nf_iyH2(%rsp),%xmm1
        movapd nb304nf_izH2(%rsp),%xmm2
        movapd nb304nf_ixH2(%rsp),%xmm3
        movapd nb304nf_iyH2(%rsp),%xmm4
        movapd nb304nf_izH2(%rsp),%xmm5
        subsd  nb304nf_jxM(%rsp),%xmm0
        subsd  nb304nf_jyM(%rsp),%xmm1
        subsd  nb304nf_jzM(%rsp),%xmm2
        subsd  nb304nf_jxH1(%rsp),%xmm3
        subsd  nb304nf_jyH1(%rsp),%xmm4
        subsd  nb304nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH2M(%rsp)
        movapd %xmm4,nb304nf_rsqH2H1(%rsp)

        movapd nb304nf_ixH2(%rsp),%xmm0
        movapd nb304nf_iyH2(%rsp),%xmm1
        movapd nb304nf_izH2(%rsp),%xmm2
        subsd  nb304nf_jxH2(%rsp),%xmm0
        subsd  nb304nf_jyH2(%rsp),%xmm1
        subsd  nb304nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb304nf_rsqH2H2(%rsp)

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
        movapd  nb304nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb304nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH2H2(%rsp)
        movapd %xmm5,nb304nf_rinvH2H1(%rsp)

        movapd nb304nf_rsqMM(%rsp),%xmm0
        movapd nb304nf_rsqMH1(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb304nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb304nf_rinvMM(%rsp)
        movapd %xmm5,nb304nf_rinvMH1(%rsp)

        movapd nb304nf_rsqMH2(%rsp),%xmm0
        movapd nb304nf_rsqH1M(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb304nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvMH2(%rsp)
        movapd %xmm5,nb304nf_rinvH1M(%rsp)

        movapd nb304nf_rsqH1H1(%rsp),%xmm0
        movapd nb304nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb304nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb304nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH1H1(%rsp)
        movapd %xmm5,nb304nf_rinvH1H2(%rsp)

        movapd nb304nf_rsqH2M(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb304nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb304nf_rinvMM(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMM(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## M-H1 interaction 
        movapd nb304nf_rinvMH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## M-H2 interaction  
        movapd nb304nf_rinvMH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H1-M interaction 
        movapd nb304nf_rinvH1M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H1-H1 interaction 
        movapd nb304nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H1-H2 interaction 
        movapd nb304nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H2-M interaction 
        movapd nb304nf_rinvH2M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb304nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb304nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb304nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%rsp)

_nb_kernel304nf_x86_64_sse2.nb304nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb304nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb304nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb304nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb304nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb304nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304nf_x86_64_sse2.nb304nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304nf_n(%rsp)
        jmp _nb_kernel304nf_x86_64_sse2.nb304nf_outer
_nb_kernel304nf_x86_64_sse2.nb304nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304nf_x86_64_sse2.nb304nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304nf_x86_64_sse2.nb304nf_threadloop
_nb_kernel304nf_x86_64_sse2.nb304nf_end: 
        movl nb304nf_nouter(%rsp),%eax
        movl nb304nf_ninner(%rsp),%ebx
        movq nb304nf_outeriter(%rbp),%rcx
        movq nb304nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $792,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

