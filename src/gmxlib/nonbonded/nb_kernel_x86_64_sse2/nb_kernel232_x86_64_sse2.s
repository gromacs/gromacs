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




.globl nb_kernel232_x86_64_sse2
.globl _nb_kernel232_x86_64_sse2
nb_kernel232_x86_64_sse2:       
_nb_kernel232_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb232_fshift, 16
.set nb232_gid, 24
.set nb232_pos, 32
.set nb232_faction, 40
.set nb232_charge, 48
.set nb232_p_facel, 56
.set nb232_argkrf, 64
.set nb232_argcrf, 72
.set nb232_Vc, 80
.set nb232_type, 88
.set nb232_p_ntype, 96
.set nb232_vdwparam, 104
.set nb232_Vvdw, 112
.set nb232_p_tabscale, 120
.set nb232_VFtab, 128
.set nb232_invsqrta, 136
.set nb232_dvda, 144
.set nb232_p_gbtabscale, 152
.set nb232_GBtab, 160
.set nb232_p_nthreads, 168
.set nb232_count, 176
.set nb232_mtx, 184
.set nb232_outeriter, 192
.set nb232_inneriter, 200
.set nb232_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb232_ixO, 0
.set nb232_iyO, 16
.set nb232_izO, 32
.set nb232_ixH1, 48
.set nb232_iyH1, 64
.set nb232_izH1, 80
.set nb232_ixH2, 96
.set nb232_iyH2, 112
.set nb232_izH2, 128
.set nb232_jxO, 144
.set nb232_jyO, 160
.set nb232_jzO, 176
.set nb232_jxH1, 192
.set nb232_jyH1, 208
.set nb232_jzH1, 224
.set nb232_jxH2, 240
.set nb232_jyH2, 256
.set nb232_jzH2, 272
.set nb232_dxOO, 288
.set nb232_dyOO, 304
.set nb232_dzOO, 320
.set nb232_dxOH1, 336
.set nb232_dyOH1, 352
.set nb232_dzOH1, 368
.set nb232_dxOH2, 384
.set nb232_dyOH2, 400
.set nb232_dzOH2, 416
.set nb232_dxH1O, 432
.set nb232_dyH1O, 448
.set nb232_dzH1O, 464
.set nb232_dxH1H1, 480
.set nb232_dyH1H1, 496
.set nb232_dzH1H1, 512
.set nb232_dxH1H2, 528
.set nb232_dyH1H2, 544
.set nb232_dzH1H2, 560
.set nb232_dxH2O, 576
.set nb232_dyH2O, 592
.set nb232_dzH2O, 608
.set nb232_dxH2H1, 624
.set nb232_dyH2H1, 640
.set nb232_dzH2H1, 656
.set nb232_dxH2H2, 672
.set nb232_dyH2H2, 688
.set nb232_dzH2H2, 704
.set nb232_qqOO, 720
.set nb232_qqOH, 736
.set nb232_qqHH, 752
.set nb232_c6, 768
.set nb232_c12, 784
.set nb232_tsc, 800
.set nb232_fstmp, 816
.set nb232_vctot, 832
.set nb232_Vvdwtot, 848
.set nb232_fixO, 864
.set nb232_fiyO, 880
.set nb232_fizO, 896
.set nb232_fixH1, 912
.set nb232_fiyH1, 928
.set nb232_fizH1, 944
.set nb232_fixH2, 960
.set nb232_fiyH2, 976
.set nb232_fizH2, 992
.set nb232_fjxO, 1008
.set nb232_fjyO, 1024
.set nb232_fjzO, 1040
.set nb232_fjxH1, 1056
.set nb232_fjyH1, 1072
.set nb232_fjzH1, 1088
.set nb232_fjxH2, 1104
.set nb232_fjyH2, 1120
.set nb232_fjzH2, 1136
.set nb232_half, 1152
.set nb232_three, 1168
.set nb232_rsqOO, 1184
.set nb232_rsqOH1, 1200
.set nb232_rsqOH2, 1216
.set nb232_rsqH1O, 1232
.set nb232_rsqH1H1, 1248
.set nb232_rsqH1H2, 1264
.set nb232_rsqH2O, 1280
.set nb232_rsqH2H1, 1296
.set nb232_rsqH2H2, 1312
.set nb232_rinvOO, 1328
.set nb232_rinvOH1, 1344
.set nb232_rinvOH2, 1360
.set nb232_rinvH1O, 1376
.set nb232_rinvH1H1, 1392
.set nb232_rinvH1H2, 1408
.set nb232_rinvH2O, 1424
.set nb232_rinvH2H1, 1440
.set nb232_rinvH2H2, 1456
.set nb232_two, 1472
.set nb232_krf, 1488
.set nb232_crf, 1504
.set nb232_nri, 1520
.set nb232_iinr, 1528
.set nb232_jindex, 1536
.set nb232_jjnr, 1544
.set nb232_shift, 1552
.set nb232_shiftvec, 1560
.set nb232_facel, 1568
.set nb232_innerjjnr, 1576
.set nb232_is3, 1584
.set nb232_ii3, 1588
.set nb232_innerk, 1592
.set nb232_n, 1596
.set nb232_nn1, 1600
.set nb232_nouter, 1604
.set nb232_ninner, 1608
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1624,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb232_nouter(%rsp)
        movl %eax,nb232_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb232_nri(%rsp)
        movq %rsi,nb232_iinr(%rsp)
        movq %rdx,nb232_jindex(%rsp)
        movq %rcx,nb232_jjnr(%rsp)
        movq %r8,nb232_shift(%rsp)
        movq %r9,nb232_shiftvec(%rsp)
        movq nb232_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb232_facel(%rsp)

        movq nb232_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb232_tsc(%rsp)

        movq nb232_argkrf(%rbp),%rsi
        movq nb232_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb232_krf(%rsp)
        movapd %xmm2,nb232_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb232_half(%rsp)
        movl %ebx,nb232_half+4(%rsp)
        movsd nb232_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb232_half(%rsp)
        movapd %xmm2,nb232_two(%rsp)
        movapd %xmm3,nb232_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb232_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb232_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb232_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb232_qqOO(%rsp)
        movapd %xmm4,nb232_qqOH(%rsp)
        movapd %xmm5,nb232_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb232_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb232_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb232_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb232_c6(%rsp)
        movapd %xmm1,nb232_c12(%rsp)

_nb_kernel232_x86_64_sse2.nb232_threadloop: 
        movq  nb232_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel232_x86_64_sse2.nb232_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel232_x86_64_sse2.nb232_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb232_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb232_n(%rsp)
        movl %ebx,nb232_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel232_x86_64_sse2.nb232_outerstart
        jmp _nb_kernel232_x86_64_sse2.nb232_end

_nb_kernel232_x86_64_sse2.nb232_outerstart: 
        ## ebx contains number of outer iterations
        addl nb232_nouter(%rsp),%ebx
        movl %ebx,nb232_nouter(%rsp)

_nb_kernel232_x86_64_sse2.nb232_outer: 
        movq  nb232_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb232_is3(%rsp)      ## store is3 

        movq  nb232_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movlpd (%rax,%rbx,8),%xmm0
        movlpd 8(%rax,%rbx,8),%xmm1
        movlpd 16(%rax,%rbx,8),%xmm2

        movq  nb232_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb232_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb232_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb232_ixO(%rsp)
        movapd %xmm4,nb232_iyO(%rsp)
        movapd %xmm5,nb232_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb232_ixH1(%rsp)
        movapd %xmm1,nb232_iyH1(%rsp)
        movapd %xmm2,nb232_izH1(%rsp)
        movapd %xmm3,nb232_ixH2(%rsp)
        movapd %xmm4,nb232_iyH2(%rsp)
        movapd %xmm5,nb232_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb232_vctot(%rsp)
        movapd %xmm4,nb232_Vvdwtot(%rsp)
        movapd %xmm4,nb232_fixO(%rsp)
        movapd %xmm4,nb232_fiyO(%rsp)
        movapd %xmm4,nb232_fizO(%rsp)
        movapd %xmm4,nb232_fixH1(%rsp)
        movapd %xmm4,nb232_fiyH1(%rsp)
        movapd %xmm4,nb232_fizH1(%rsp)
        movapd %xmm4,nb232_fixH2(%rsp)
        movapd %xmm4,nb232_fiyH2(%rsp)
        movapd %xmm4,nb232_fizH2(%rsp)

        movq  nb232_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb232_pos(%rbp),%rsi
        movq  nb232_faction(%rbp),%rdi
        movq  nb232_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb232_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb232_ninner(%rsp),%ecx
        movl  %ecx,nb232_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb232_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel232_x86_64_sse2.nb232_unroll_loop
        jmp   _nb_kernel232_x86_64_sse2.nb232_checksingle
_nb_kernel232_x86_64_sse2.nb232_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb232_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb232_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb232_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j O coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb232_ixO(%rsp),%xmm0
    subpd nb232_iyO(%rsp),%xmm1
    subpd nb232_izO(%rsp),%xmm2
    subpd nb232_ixH1(%rsp),%xmm3
    subpd nb232_iyH1(%rsp),%xmm4
    subpd nb232_izH1(%rsp),%xmm5
    subpd nb232_ixH2(%rsp),%xmm6
    subpd nb232_iyH2(%rsp),%xmm7
    subpd nb232_izH2(%rsp),%xmm8

        movapd %xmm0,nb232_dxOO(%rsp)
        movapd %xmm1,nb232_dyOO(%rsp)
        movapd %xmm2,nb232_dzOO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb232_dxH1O(%rsp)
        movapd %xmm4,nb232_dyH1O(%rsp)
        movapd %xmm5,nb232_dzH1O(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb232_dxH2O(%rsp)
        movapd %xmm7,nb232_dyH2O(%rsp)
        movapd %xmm8,nb232_dzH2O(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
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

        movapd  nb232_three(%rsp),%xmm9
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb232_three(%rsp),%xmm1
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOO 
        mulpd   %xmm15,%xmm10 ##   rinvH1O
    mulpd   %xmm15,%xmm11 ##   rinvH2O

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm0,nb232_rsqOO(%rsp)
    movapd %xmm3,nb232_rsqOH1(%rsp)
    movapd %xmm6,nb232_rsqOH2(%rsp)
    movapd %xmm9,nb232_rinvOO(%rsp)
    movapd %xmm10,nb232_rinvOH1(%rsp)
    movapd %xmm11,nb232_rinvOH2(%rsp)

    ## table LJ interaction
    mulpd  %xmm9,%xmm0
    mulpd  nb232_tsc(%rsp),%xmm0   ## rtab

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2

    ## multiply by 8
    pslld   $3,%xmm1

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    movd    %xmm1,%r8d
    movd    %xmm13,%r10d

    ## calculate eps
    subpd     %xmm2,%xmm0
    movq nb232_VFtab(%rbp),%rsi

    movlpd (%rsi,%r8,8),%xmm4
        movlpd 8(%rsi,%r8,8),%xmm5
        movlpd 16(%rsi,%r8,8),%xmm6
        movlpd 24(%rsi,%r8,8),%xmm7
    movlpd 32(%rsi,%r8,8),%xmm8
        movlpd 40(%rsi,%r8,8),%xmm9
        movlpd 48(%rsi,%r8,8),%xmm10
        movlpd 56(%rsi,%r8,8),%xmm11

    movhpd (%rsi,%r10,8),%xmm4
        movhpd 8(%rsi,%r10,8),%xmm5
        movhpd 16(%rsi,%r10,8),%xmm6
        movhpd 24(%rsi,%r10,8),%xmm7
    movhpd 32(%rsi,%r10,8),%xmm8
        movhpd 40(%rsi,%r10,8),%xmm9
        movhpd 48(%rsi,%r10,8),%xmm10
        movhpd 56(%rsi,%r10,8),%xmm11
    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulpd  %xmm0,%xmm7   ## Heps
    mulpd  %xmm0,%xmm11
    mulpd  %xmm0,%xmm6  ## Geps
    mulpd  %xmm0,%xmm10
    mulpd  %xmm0,%xmm7  ## Heps2
    mulpd  %xmm0,%xmm11
    addpd  %xmm6,%xmm5 ## F+Geps
    addpd  %xmm10,%xmm9
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm11,%xmm9
    addpd  %xmm7,%xmm7   ## 2*Heps2
    addpd  %xmm11,%xmm11
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm10,%xmm11

    addpd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm9,%xmm11
    mulpd  %xmm0,%xmm5 ## eps*Fp
    mulpd  %xmm0,%xmm9
    movapd nb232_c6(%rsp),%xmm12
    movapd nb232_c12(%rsp),%xmm13
    addpd  %xmm4,%xmm5 ## VV
    addpd  %xmm8,%xmm9

    mulpd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulpd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addpd  %xmm9,%xmm5
    addpd  nb232_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb232_Vvdwtot(%rsp)

    mulpd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulpd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addpd  %xmm11,%xmm7
    mulpd  nb232_tsc(%rsp),%xmm7
    movapd %xmm7,nb232_fstmp(%rsp)

    ## Coulomb reaction-field interaction
    movapd nb232_rsqOO(%rsp),%xmm0
    movapd nb232_rsqOH1(%rsp),%xmm3
    movapd nb232_rsqOH2(%rsp),%xmm6
    movapd nb232_rinvOO(%rsp),%xmm9
    movapd nb232_rinvOH1(%rsp),%xmm10
    movapd nb232_rinvOH2(%rsp),%xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb232_krf(%rsp),%xmm2
    mulpd  %xmm10,%xmm10 ## rinvsq
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
    subpd  nb232_crf(%rsp),%xmm2     ## rinv+krsq-crf
    subpd  nb232_crf(%rsp),%xmm5
    subpd  nb232_crf(%rsp),%xmm8
    mulpd  nb232_qqOO(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf)
    mulpd  nb232_qqOH(%rsp),%xmm5   ## voul=qq*(rinv+ krsq-crf)
    mulpd  nb232_qqOH(%rsp),%xmm8   ## voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  nb232_qqOO(%rsp),%xmm1     ## (rinv-2*krsq)*qq
    mulpd  nb232_qqOH(%rsp),%xmm4
    mulpd  nb232_qqOH(%rsp),%xmm7
    addpd  nb232_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb232_vctot(%rsp)
    mulpd  %xmm9,%xmm1  ## fijC
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11
    subpd nb232_fstmp(%rsp),%xmm1
    mulpd %xmm1,%xmm9

    ## move j O forces to xmm0-xmm2
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb232_dxOO(%rsp),%xmm7
        mulpd nb232_dyOO(%rsp),%xmm8
        mulpd nb232_dzOO(%rsp),%xmm9
        mulpd nb232_dxH1O(%rsp),%xmm10
        mulpd nb232_dyH1O(%rsp),%xmm11
        mulpd nb232_dzH1O(%rsp),%xmm12
        mulpd nb232_dxH2O(%rsp),%xmm13
        mulpd nb232_dyH2O(%rsp),%xmm14
        mulpd nb232_dzH2O(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb232_fixO(%rsp),%xmm7
    addpd nb232_fiyO(%rsp),%xmm8
    addpd nb232_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb232_fixH1(%rsp),%xmm10
    addpd nb232_fiyH1(%rsp),%xmm11
    addpd nb232_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb232_fixH2(%rsp),%xmm13
    addpd nb232_fiyH2(%rsp),%xmm14
    addpd nb232_fizH2(%rsp),%xmm15

    movapd %xmm7,nb232_fixO(%rsp)
    movapd %xmm8,nb232_fiyO(%rsp)
    movapd %xmm9,nb232_fizO(%rsp)
    movapd %xmm10,nb232_fixH1(%rsp)
    movapd %xmm11,nb232_fiyH1(%rsp)
    movapd %xmm12,nb232_fizH1(%rsp)
    movapd %xmm13,nb232_fixH2(%rsp)
    movapd %xmm14,nb232_fiyH2(%rsp)
    movapd %xmm15,nb232_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## move j H1 coordinates to local temp variables 
    movq nb232_pos(%rbp),%rsi
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

    subpd nb232_ixO(%rsp),%xmm0
    subpd nb232_iyO(%rsp),%xmm1
    subpd nb232_izO(%rsp),%xmm2
    subpd nb232_ixH1(%rsp),%xmm3
    subpd nb232_iyH1(%rsp),%xmm4
    subpd nb232_izH1(%rsp),%xmm5
    subpd nb232_ixH2(%rsp),%xmm6
    subpd nb232_iyH2(%rsp),%xmm7
    subpd nb232_izH2(%rsp),%xmm8

        movapd %xmm0,nb232_dxOH1(%rsp)
        movapd %xmm1,nb232_dyOH1(%rsp)
        movapd %xmm2,nb232_dzOH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb232_dxH1H1(%rsp)
        movapd %xmm4,nb232_dyH1H1(%rsp)
        movapd %xmm5,nb232_dzH1H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb232_dxH2H1(%rsp)
        movapd %xmm7,nb232_dyH2H1(%rsp)
        movapd %xmm8,nb232_dzH2H1(%rsp)
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

        movapd  nb232_three(%rsp),%xmm9
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb232_three(%rsp),%xmm1
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH1
        mulpd   %xmm15,%xmm10 ##   rinvH1H1
    mulpd   %xmm15,%xmm11 ##   rinvH2H1


        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb232_krf(%rsp),%xmm2
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
    movapd nb232_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb232_qqOH(%rsp),%xmm12
    movapd nb232_qqHH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm13,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb232_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb232_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb232_faction(%rbp),%rdi
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

        mulpd nb232_dxOH1(%rsp),%xmm7
        mulpd nb232_dyOH1(%rsp),%xmm8
        mulpd nb232_dzOH1(%rsp),%xmm9
        mulpd nb232_dxH1H1(%rsp),%xmm10
        mulpd nb232_dyH1H1(%rsp),%xmm11
        mulpd nb232_dzH1H1(%rsp),%xmm12
        mulpd nb232_dxH2H1(%rsp),%xmm13
        mulpd nb232_dyH2H1(%rsp),%xmm14
        mulpd nb232_dzH2H1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb232_fixO(%rsp),%xmm7
    addpd nb232_fiyO(%rsp),%xmm8
    addpd nb232_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb232_fixH1(%rsp),%xmm10
    addpd nb232_fiyH1(%rsp),%xmm11
    addpd nb232_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb232_fixH2(%rsp),%xmm13
    addpd nb232_fiyH2(%rsp),%xmm14
    addpd nb232_fizH2(%rsp),%xmm15

    movapd %xmm7,nb232_fixO(%rsp)
    movapd %xmm8,nb232_fiyO(%rsp)
    movapd %xmm9,nb232_fizO(%rsp)
    movapd %xmm10,nb232_fixH1(%rsp)
    movapd %xmm11,nb232_fiyH1(%rsp)
    movapd %xmm12,nb232_fizH1(%rsp)
    movapd %xmm13,nb232_fixH2(%rsp)
    movapd %xmm14,nb232_fiyH2(%rsp)
    movapd %xmm15,nb232_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb232_pos(%rbp),%rsi
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

    subpd nb232_ixO(%rsp),%xmm0
    subpd nb232_iyO(%rsp),%xmm1
    subpd nb232_izO(%rsp),%xmm2
    subpd nb232_ixH1(%rsp),%xmm3
    subpd nb232_iyH1(%rsp),%xmm4
    subpd nb232_izH1(%rsp),%xmm5
    subpd nb232_ixH2(%rsp),%xmm6
    subpd nb232_iyH2(%rsp),%xmm7
    subpd nb232_izH2(%rsp),%xmm8

        movapd %xmm0,nb232_dxOH2(%rsp)
        movapd %xmm1,nb232_dyOH2(%rsp)
        movapd %xmm2,nb232_dzOH2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb232_dxH1H2(%rsp)
        movapd %xmm4,nb232_dyH1H2(%rsp)
        movapd %xmm5,nb232_dzH1H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb232_dxH2H2(%rsp)
        movapd %xmm7,nb232_dyH2H2(%rsp)
        movapd %xmm8,nb232_dzH2H2(%rsp)
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

        movapd  nb232_three(%rsp),%xmm9
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb232_three(%rsp),%xmm1
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

        movapd  nb232_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH2
        mulpd   %xmm15,%xmm10 ##   rinvH1H2
    mulpd   %xmm15,%xmm11 ##   rinvH2H2


        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb232_krf(%rsp),%xmm2
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
    movapd nb232_crf(%rsp),%xmm14
    subpd  %xmm14,%xmm2  ## rinv+krsq-crf
    subpd  %xmm14,%xmm5
    subpd  %xmm14,%xmm8
    movapd nb232_qqOH(%rsp),%xmm12
    movapd nb232_qqHH(%rsp),%xmm13
    mulpd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulpd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addpd  %xmm0,%xmm0 ## 2*krsq
    addpd  %xmm3,%xmm3
    addpd  %xmm6,%xmm6
    subpd  %xmm0,%xmm1 ## rinv-2*krsq
    subpd  %xmm3,%xmm4
    subpd  %xmm6,%xmm7
    mulpd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulpd  %xmm13,%xmm4
    mulpd  %xmm13,%xmm7
    addpd  nb232_vctot(%rsp),%xmm2
    addpd  %xmm8,%xmm5
    addpd  %xmm5,%xmm2
    movapd %xmm2,nb232_vctot(%rsp)

    mulpd  %xmm1,%xmm9  ## fscal
    mulpd  %xmm4,%xmm10
    mulpd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb232_faction(%rbp),%rdi
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

        mulpd nb232_dxOH2(%rsp),%xmm7
        mulpd nb232_dyOH2(%rsp),%xmm8
        mulpd nb232_dzOH2(%rsp),%xmm9
        mulpd nb232_dxH1H2(%rsp),%xmm10
        mulpd nb232_dyH1H2(%rsp),%xmm11
        mulpd nb232_dzH1H2(%rsp),%xmm12
        mulpd nb232_dxH2H2(%rsp),%xmm13
        mulpd nb232_dyH2H2(%rsp),%xmm14
        mulpd nb232_dzH2H2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb232_fixO(%rsp),%xmm7
    addpd nb232_fiyO(%rsp),%xmm8
    addpd nb232_fizO(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb232_fixH1(%rsp),%xmm10
    addpd nb232_fiyH1(%rsp),%xmm11
    addpd nb232_fizH1(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb232_fixH2(%rsp),%xmm13
    addpd nb232_fiyH2(%rsp),%xmm14
    addpd nb232_fizH2(%rsp),%xmm15

    movapd %xmm7,nb232_fixO(%rsp)
    movapd %xmm8,nb232_fiyO(%rsp)
    movapd %xmm9,nb232_fizO(%rsp)
    movapd %xmm10,nb232_fixH1(%rsp)
    movapd %xmm11,nb232_fiyH1(%rsp)
    movapd %xmm12,nb232_fizH1(%rsp)
    movapd %xmm13,nb232_fixH2(%rsp)
    movapd %xmm14,nb232_fiyH2(%rsp)
    movapd %xmm15,nb232_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb232_innerk(%rsp)
        jl    _nb_kernel232_x86_64_sse2.nb232_checksingle
        jmp   _nb_kernel232_x86_64_sse2.nb232_unroll_loop
_nb_kernel232_x86_64_sse2.nb232_checksingle: 
        movl  nb232_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel232_x86_64_sse2.nb232_dosingle
        jmp   _nb_kernel232_x86_64_sse2.nb232_updateouterdata
_nb_kernel232_x86_64_sse2.nb232_dosingle: 
        movq  nb232_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb232_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j O coordinates to local temp variables 
    movsd (%rsi,%rax,8),%xmm0
    movsd 8(%rsi,%rax,8),%xmm1
    movsd 16(%rsi,%rax,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb232_ixO(%rsp),%xmm0
    subsd nb232_iyO(%rsp),%xmm1
    subsd nb232_izO(%rsp),%xmm2
    subsd nb232_ixH1(%rsp),%xmm3
    subsd nb232_iyH1(%rsp),%xmm4
    subsd nb232_izH1(%rsp),%xmm5
    subsd nb232_ixH2(%rsp),%xmm6
    subsd nb232_iyH2(%rsp),%xmm7
    subsd nb232_izH2(%rsp),%xmm8

        movsd %xmm0,nb232_dxOO(%rsp)
        movsd %xmm1,nb232_dyOO(%rsp)
        movsd %xmm2,nb232_dzOO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb232_dxH1O(%rsp)
        movsd %xmm4,nb232_dyH1O(%rsp)
        movsd %xmm5,nb232_dzH1O(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb232_dxH2O(%rsp)
        movsd %xmm7,nb232_dyH2O(%rsp)
        movsd %xmm8,nb232_dzH2O(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movapd  nb232_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb232_three(%rsp),%xmm1
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

        movapd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOO 
        mulsd   %xmm15,%xmm10 ##   rinvH1O
    mulsd   %xmm15,%xmm11 ##   rinvH2O

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11
    movapd %xmm0,nb232_rsqOO(%rsp)
    movapd %xmm3,nb232_rsqOH1(%rsp)
    movapd %xmm6,nb232_rsqOH2(%rsp)
    movapd %xmm9,nb232_rinvOO(%rsp)
    movapd %xmm10,nb232_rinvOH1(%rsp)
    movapd %xmm11,nb232_rinvOH2(%rsp)

    ## table LJ interaction
    mulsd  %xmm9,%xmm0
    mulsd  nb232_tsc(%rsp),%xmm0   ## rtab

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2

    ## multiply by 8
    shll  $3,%r8d

    ## calculate eps
    subsd     %xmm2,%xmm0
    movq nb232_VFtab(%rbp),%rsi

    movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm5
        movsd 16(%rsi,%r8,8),%xmm6
        movsd 24(%rsi,%r8,8),%xmm7
    movsd 32(%rsi,%r8,8),%xmm8
        movsd 40(%rsi,%r8,8),%xmm9
        movsd 48(%rsi,%r8,8),%xmm10
        movsd 56(%rsi,%r8,8),%xmm11

    ## dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    mulsd  %xmm0,%xmm7   ## Heps
    mulsd  %xmm0,%xmm11
    mulsd  %xmm0,%xmm6  ## Geps
    mulsd  %xmm0,%xmm10
    mulsd  %xmm0,%xmm7  ## Heps2
    mulsd  %xmm0,%xmm11
    addsd  %xmm6,%xmm5 ## F+Geps
    addsd  %xmm10,%xmm9
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm11,%xmm9
    addsd  %xmm7,%xmm7   ## 2*Heps2
    addsd  %xmm11,%xmm11
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm10,%xmm11

    addsd  %xmm5,%xmm7 ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm9,%xmm11
    mulsd  %xmm0,%xmm5 ## eps*Fp
    mulsd  %xmm0,%xmm9
    movsd nb232_c6(%rsp),%xmm12
    movsd nb232_c12(%rsp),%xmm13
    addsd  %xmm4,%xmm5 ## VV
    addsd  %xmm8,%xmm9

    mulsd  %xmm12,%xmm5 ## VV*c6 = vnb6
    mulsd  %xmm13,%xmm9 ## VV*c12 = vnb12
    addsd  %xmm9,%xmm5
    addsd  nb232_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb232_Vvdwtot(%rsp)

    mulsd  %xmm12,%xmm7  ## FF*c6 = fnb6
    mulsd  %xmm13,%xmm11  ## FF*c12  = fnb12
    addsd  %xmm11,%xmm7
    mulsd  nb232_tsc(%rsp),%xmm7
    movsd %xmm7,nb232_fstmp(%rsp)

    ## Coulomb reaction-field interaction
    movapd nb232_rsqOO(%rsp),%xmm0
    movapd nb232_rsqOH1(%rsp),%xmm3
    movapd nb232_rsqOH2(%rsp),%xmm6
    movapd nb232_rinvOO(%rsp),%xmm9
    movapd nb232_rinvOH1(%rsp),%xmm10
    movapd nb232_rinvOH2(%rsp),%xmm11

    movapd %xmm9,%xmm1 ## copy of rinv
    movapd %xmm10,%xmm4
    movapd %xmm11,%xmm7
    movapd nb232_krf(%rsp),%xmm2
    mulsd  %xmm10,%xmm10 ## rinvsq
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
    subsd  nb232_crf(%rsp),%xmm2     ## rinv+krsq-crf
    subsd  nb232_crf(%rsp),%xmm5
    subsd  nb232_crf(%rsp),%xmm8
    mulsd  nb232_qqOO(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf)
    mulsd  nb232_qqOH(%rsp),%xmm5   ## voul=qq*(rinv+ krsq-crf)
    mulsd  nb232_qqOH(%rsp),%xmm8   ## voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  nb232_qqOO(%rsp),%xmm1     ## (rinv-2*krsq)*qq
    mulsd  nb232_qqOH(%rsp),%xmm4
    mulsd  nb232_qqOH(%rsp),%xmm7
    addsd  nb232_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd %xmm2,nb232_vctot(%rsp)
    mulsd  %xmm9,%xmm1  ## fijC
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    subsd nb232_fstmp(%rsp),%xmm1
    mulsd %xmm1,%xmm9

    ## move j O forces to xmm0-xmm2
    movq nb232_faction(%rbp),%rdi
        movsd (%rdi,%rax,8),%xmm0
        movsd 8(%rdi,%rax,8),%xmm1
        movsd 16(%rdi,%rax,8),%xmm2

    movsd %xmm9,%xmm7
    movsd %xmm9,%xmm8
    movsd %xmm11,%xmm13
    movsd %xmm11,%xmm14
    movsd %xmm11,%xmm15
    movsd %xmm10,%xmm11
    movsd %xmm10,%xmm12

        mulsd nb232_dxOO(%rsp),%xmm7
        mulsd nb232_dyOO(%rsp),%xmm8
        mulsd nb232_dzOO(%rsp),%xmm9
        mulsd nb232_dxH1O(%rsp),%xmm10
        mulsd nb232_dyH1O(%rsp),%xmm11
        mulsd nb232_dzH1O(%rsp),%xmm12
        mulsd nb232_dxH2O(%rsp),%xmm13
        mulsd nb232_dyH2O(%rsp),%xmm14
        mulsd nb232_dzH2O(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb232_fixO(%rsp),%xmm7
    addsd nb232_fiyO(%rsp),%xmm8
    addsd nb232_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb232_fixH1(%rsp),%xmm10
    addsd nb232_fiyH1(%rsp),%xmm11
    addsd nb232_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb232_fixH2(%rsp),%xmm13
    addsd nb232_fiyH2(%rsp),%xmm14
    addsd nb232_fizH2(%rsp),%xmm15

    movsd %xmm7,nb232_fixO(%rsp)
    movsd %xmm8,nb232_fiyO(%rsp)
    movsd %xmm9,nb232_fizO(%rsp)
    movsd %xmm10,nb232_fixH1(%rsp)
    movsd %xmm11,nb232_fiyH1(%rsp)
    movsd %xmm12,nb232_fizH1(%rsp)
    movsd %xmm13,nb232_fixH2(%rsp)
    movsd %xmm14,nb232_fiyH2(%rsp)
    movsd %xmm15,nb232_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movsd %xmm0,(%rdi,%rax,8)
        movsd %xmm1,8(%rdi,%rax,8)
        movsd %xmm2,16(%rdi,%rax,8)

        ## move j H1 coordinates to local temp variables 
    movq nb232_pos(%rbp),%rsi
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

    subsd nb232_ixO(%rsp),%xmm0
    subsd nb232_iyO(%rsp),%xmm1
    subsd nb232_izO(%rsp),%xmm2
    subsd nb232_ixH1(%rsp),%xmm3
    subsd nb232_iyH1(%rsp),%xmm4
    subsd nb232_izH1(%rsp),%xmm5
    subsd nb232_ixH2(%rsp),%xmm6
    subsd nb232_iyH2(%rsp),%xmm7
    subsd nb232_izH2(%rsp),%xmm8

        movsd %xmm0,nb232_dxOH1(%rsp)
        movsd %xmm1,nb232_dyOH1(%rsp)
        movsd %xmm2,nb232_dzOH1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb232_dxH1H1(%rsp)
        movsd %xmm4,nb232_dyH1H1(%rsp)
        movsd %xmm5,nb232_dzH1H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb232_dxH2H1(%rsp)
        movsd %xmm7,nb232_dyH2H1(%rsp)
        movsd %xmm8,nb232_dzH2H1(%rsp)
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

        movsd  nb232_three(%rsp),%xmm9
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

        movsd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb232_three(%rsp),%xmm1
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

        movsd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH1
        mulsd   %xmm15,%xmm10 ##   rinvH1H1
    mulsd   %xmm15,%xmm11 ##   rinvH2H1

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb232_krf(%rsp),%xmm2
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
    movsd nb232_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb232_qqOH(%rsp),%xmm12
    movsd nb232_qqHH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm13,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb232_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd %xmm2,nb232_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb232_faction(%rbp),%rdi
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

        mulsd nb232_dxOH1(%rsp),%xmm7
        mulsd nb232_dyOH1(%rsp),%xmm8
        mulsd nb232_dzOH1(%rsp),%xmm9
        mulsd nb232_dxH1H1(%rsp),%xmm10
        mulsd nb232_dyH1H1(%rsp),%xmm11
        mulsd nb232_dzH1H1(%rsp),%xmm12
        mulsd nb232_dxH2H1(%rsp),%xmm13
        mulsd nb232_dyH2H1(%rsp),%xmm14
        mulsd nb232_dzH2H1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb232_fixO(%rsp),%xmm7
    addsd nb232_fiyO(%rsp),%xmm8
    addsd nb232_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb232_fixH1(%rsp),%xmm10
    addsd nb232_fiyH1(%rsp),%xmm11
    addsd nb232_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb232_fixH2(%rsp),%xmm13
    addsd nb232_fiyH2(%rsp),%xmm14
    addsd nb232_fizH2(%rsp),%xmm15

    movsd %xmm7,nb232_fixO(%rsp)
    movsd %xmm8,nb232_fiyO(%rsp)
    movsd %xmm9,nb232_fizO(%rsp)
    movsd %xmm10,nb232_fixH1(%rsp)
    movsd %xmm11,nb232_fiyH1(%rsp)
    movsd %xmm12,nb232_fizH1(%rsp)
    movsd %xmm13,nb232_fixH2(%rsp)
    movsd %xmm14,nb232_fiyH2(%rsp)
    movsd %xmm15,nb232_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movq nb232_pos(%rbp),%rsi
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

    subsd nb232_ixO(%rsp),%xmm0
    subsd nb232_iyO(%rsp),%xmm1
    subsd nb232_izO(%rsp),%xmm2
    subsd nb232_ixH1(%rsp),%xmm3
    subsd nb232_iyH1(%rsp),%xmm4
    subsd nb232_izH1(%rsp),%xmm5
    subsd nb232_ixH2(%rsp),%xmm6
    subsd nb232_iyH2(%rsp),%xmm7
    subsd nb232_izH2(%rsp),%xmm8

        movsd %xmm0,nb232_dxOH2(%rsp)
        movsd %xmm1,nb232_dyOH2(%rsp)
        movsd %xmm2,nb232_dzOH2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb232_dxH1H2(%rsp)
        movsd %xmm4,nb232_dyH1H2(%rsp)
        movsd %xmm5,nb232_dzH1H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb232_dxH2H2(%rsp)
        movsd %xmm7,nb232_dyH2H2(%rsp)
        movsd %xmm8,nb232_dzH2H2(%rsp)
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

        movsd  nb232_three(%rsp),%xmm9
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

        movsd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb232_three(%rsp),%xmm1
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

        movsd  nb232_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH2
        mulsd   %xmm15,%xmm10 ##   rinvH1H2
    mulsd   %xmm15,%xmm11 ##   rinvH2H2

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd %xmm9,%xmm1 ## copy of rinv
    movsd %xmm10,%xmm4
    movsd %xmm11,%xmm7
    movsd nb232_krf(%rsp),%xmm2
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
    movsd nb232_crf(%rsp),%xmm14
    subsd  %xmm14,%xmm2  ## rinv+krsq-crf
    subsd  %xmm14,%xmm5
    subsd  %xmm14,%xmm8
    movsd nb232_qqOH(%rsp),%xmm12
    movsd nb232_qqHH(%rsp),%xmm13
    mulsd  %xmm12,%xmm2 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm5 ## xmm6=voul=qq*(rinv+ krsq-crf)
    mulsd  %xmm13,%xmm8 ## xmm6=voul=qq*(rinv+ krsq-crf)
    addsd  %xmm0,%xmm0 ## 2*krsq
    addsd  %xmm3,%xmm3
    addsd  %xmm6,%xmm6
    subsd  %xmm0,%xmm1 ## rinv-2*krsq
    subsd  %xmm3,%xmm4
    subsd  %xmm6,%xmm7
    mulsd  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulsd  %xmm13,%xmm4
    mulsd  %xmm13,%xmm7
    addsd  nb232_vctot(%rsp),%xmm2
    addsd  %xmm8,%xmm5
    addsd  %xmm5,%xmm2
    movsd %xmm2,nb232_vctot(%rsp)

    mulsd  %xmm1,%xmm9  ## fscal
    mulsd  %xmm4,%xmm10
    mulsd  %xmm7,%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb232_faction(%rbp),%rdi
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

        mulsd nb232_dxOH2(%rsp),%xmm7
        mulsd nb232_dyOH2(%rsp),%xmm8
        mulsd nb232_dzOH2(%rsp),%xmm9
        mulsd nb232_dxH1H2(%rsp),%xmm10
        mulsd nb232_dyH1H2(%rsp),%xmm11
        mulsd nb232_dzH1H2(%rsp),%xmm12
        mulsd nb232_dxH2H2(%rsp),%xmm13
        mulsd nb232_dyH2H2(%rsp),%xmm14
        mulsd nb232_dzH2H2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb232_fixO(%rsp),%xmm7
    addsd nb232_fiyO(%rsp),%xmm8
    addsd nb232_fizO(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb232_fixH1(%rsp),%xmm10
    addsd nb232_fiyH1(%rsp),%xmm11
    addsd nb232_fizH1(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb232_fixH2(%rsp),%xmm13
    addsd nb232_fiyH2(%rsp),%xmm14
    addsd nb232_fizH2(%rsp),%xmm15

    movsd %xmm7,nb232_fixO(%rsp)
    movsd %xmm8,nb232_fiyO(%rsp)
    movsd %xmm9,nb232_fizO(%rsp)
    movsd %xmm10,nb232_fixH1(%rsp)
    movsd %xmm11,nb232_fiyH1(%rsp)
    movsd %xmm12,nb232_fizH1(%rsp)
    movsd %xmm13,nb232_fixH2(%rsp)
    movsd %xmm14,nb232_fiyH2(%rsp)
    movsd %xmm15,nb232_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

_nb_kernel232_x86_64_sse2.nb232_updateouterdata: 
        movl  nb232_ii3(%rsp),%ecx
        movq  nb232_faction(%rbp),%rdi
        movq  nb232_fshift(%rbp),%rsi
        movl  nb232_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb232_fixO(%rsp),%xmm0
        movapd nb232_fiyO(%rsp),%xmm1
        movapd nb232_fizO(%rsp),%xmm2

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
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb232_fixH1(%rsp),%xmm0
        movapd nb232_fiyH1(%rsp),%xmm1
        movapd nb232_fizH1(%rsp),%xmm2

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
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb232_fixH2(%rsp),%xmm0
        movapd nb232_fiyH2(%rsp),%xmm1
        movapd nb232_fizH2(%rsp),%xmm2

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
        movl nb232_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb232_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb232_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb232_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb232_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb232_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb232_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel232_x86_64_sse2.nb232_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb232_n(%rsp)
        jmp _nb_kernel232_x86_64_sse2.nb232_outer
_nb_kernel232_x86_64_sse2.nb232_outerend: 
        ## check if more outer neighborlists remain
        movl  nb232_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel232_x86_64_sse2.nb232_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel232_x86_64_sse2.nb232_threadloop
_nb_kernel232_x86_64_sse2.nb232_end: 
        movl nb232_nouter(%rsp),%eax
        movl nb232_ninner(%rsp),%ebx
        movq nb232_outeriter(%rbp),%rcx
        movq nb232_inneriter(%rbp),%rdx
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




.globl nb_kernel232nf_x86_64_sse2
.globl _nb_kernel232nf_x86_64_sse2
nb_kernel232nf_x86_64_sse2:     
_nb_kernel232nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb232nf_fshift, 16
.set nb232nf_gid, 24
.set nb232nf_pos, 32
.set nb232nf_faction, 40
.set nb232nf_charge, 48
.set nb232nf_p_facel, 56
.set nb232nf_argkrf, 64
.set nb232nf_argcrf, 72
.set nb232nf_Vc, 80
.set nb232nf_type, 88
.set nb232nf_p_ntype, 96
.set nb232nf_vdwparam, 104
.set nb232nf_Vvdw, 112
.set nb232nf_p_tabscale, 120
.set nb232nf_VFtab, 128
.set nb232nf_invsqrta, 136
.set nb232nf_dvda, 144
.set nb232nf_p_gbtabscale, 152
.set nb232nf_GBtab, 160
.set nb232nf_p_nthreads, 168
.set nb232nf_count, 176
.set nb232nf_mtx, 184
.set nb232nf_outeriter, 192
.set nb232nf_inneriter, 200
.set nb232nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb232nf_ixO, 0
.set nb232nf_iyO, 16
.set nb232nf_izO, 32
.set nb232nf_ixH1, 48
.set nb232nf_iyH1, 64
.set nb232nf_izH1, 80
.set nb232nf_ixH2, 96
.set nb232nf_iyH2, 112
.set nb232nf_izH2, 128
.set nb232nf_jxO, 144
.set nb232nf_jyO, 160
.set nb232nf_jzO, 176
.set nb232nf_jxH1, 192
.set nb232nf_jyH1, 208
.set nb232nf_jzH1, 224
.set nb232nf_jxH2, 240
.set nb232nf_jyH2, 256
.set nb232nf_jzH2, 272
.set nb232nf_qqOO, 288
.set nb232nf_qqOH, 304
.set nb232nf_qqHH, 320
.set nb232nf_c6, 336
.set nb232nf_c12, 352
.set nb232nf_vctot, 368
.set nb232nf_Vvdwtot, 384
.set nb232nf_half, 400
.set nb232nf_three, 416
.set nb232nf_rsqOO, 432
.set nb232nf_rsqOH1, 448
.set nb232nf_rsqOH2, 464
.set nb232nf_rsqH1O, 480
.set nb232nf_rsqH1H1, 496
.set nb232nf_rsqH1H2, 512
.set nb232nf_rsqH2O, 528
.set nb232nf_rsqH2H1, 544
.set nb232nf_rsqH2H2, 560
.set nb232nf_rinvOO, 576
.set nb232nf_rinvOH1, 592
.set nb232nf_rinvOH2, 608
.set nb232nf_rinvH1O, 624
.set nb232nf_rinvH1H1, 640
.set nb232nf_rinvH1H2, 656
.set nb232nf_rinvH2O, 672
.set nb232nf_rinvH2H1, 688
.set nb232nf_rinvH2H2, 704
.set nb232nf_krf, 720
.set nb232nf_crf, 736
.set nb232nf_tsc, 752
.set nb232nf_nri, 768
.set nb232nf_iinr, 776
.set nb232nf_jindex, 784
.set nb232nf_jjnr, 792
.set nb232nf_shift, 800
.set nb232nf_shiftvec, 808
.set nb232nf_facel, 816
.set nb232nf_innerjjnr, 824
.set nb232nf_is3, 832
.set nb232nf_ii3, 836
.set nb232nf_innerk, 840
.set nb232nf_n, 844
.set nb232nf_nn1, 848
.set nb232nf_nouter, 852
.set nb232nf_ninner, 856
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $872,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb232nf_nouter(%rsp)
        movl %eax,nb232nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb232nf_nri(%rsp)
        movq %rsi,nb232nf_iinr(%rsp)
        movq %rdx,nb232nf_jindex(%rsp)
        movq %rcx,nb232nf_jjnr(%rsp)
        movq %r8,nb232nf_shift(%rsp)
        movq %r9,nb232nf_shiftvec(%rsp)
        movq nb232nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb232nf_facel(%rsp)

        movq nb232nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb232nf_tsc(%rsp)

        movq nb232nf_argkrf(%rbp),%rsi
        movq nb232nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb232nf_krf(%rsp)
        movapd %xmm2,nb232nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb232nf_half(%rsp)
        movl %ebx,nb232nf_half+4(%rsp)
        movsd nb232nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb232nf_half(%rsp)
        movapd %xmm3,nb232nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb232nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb232nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5

        movsd nb232nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb232nf_qqOO(%rsp)
        movapd %xmm4,nb232nf_qqOH(%rsp)
        movapd %xmm5,nb232nf_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb232nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb232nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb232nf_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb232nf_c6(%rsp)
        movapd %xmm1,nb232nf_c12(%rsp)

_nb_kernel232nf_x86_64_sse2.nb232nf_threadloop: 
        movq  nb232nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel232nf_x86_64_sse2.nb232nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel232nf_x86_64_sse2.nb232nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb232nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb232nf_n(%rsp)
        movl %ebx,nb232nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel232nf_x86_64_sse2.nb232nf_outerstart
        jmp _nb_kernel232nf_x86_64_sse2.nb232nf_end

_nb_kernel232nf_x86_64_sse2.nb232nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb232nf_nouter(%rsp),%ebx
        movl %ebx,nb232nf_nouter(%rsp)

_nb_kernel232nf_x86_64_sse2.nb232nf_outer: 
        movq  nb232nf_shift(%rsp),%rax        ## eax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb232nf_is3(%rsp)            ## store is3 

        movq  nb232nf_shiftvec(%rsp),%rax     ## eax = base of shiftvec[] 

        movlpd (%rax,%rbx,8),%xmm0
        movlpd 8(%rax,%rbx,8),%xmm1
        movlpd 16(%rax,%rbx,8),%xmm2

        movq  nb232nf_iinr(%rsp),%rcx         ## ecx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb232nf_pos(%rbp),%rax      ## eax = base of pos[]  
        movl  %ebx,nb232nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb232nf_ixO(%rsp)
        movapd %xmm4,nb232nf_iyO(%rsp)
        movapd %xmm5,nb232nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb232nf_ixH1(%rsp)
        movapd %xmm1,nb232nf_iyH1(%rsp)
        movapd %xmm2,nb232nf_izH1(%rsp)
        movapd %xmm3,nb232nf_ixH2(%rsp)
        movapd %xmm4,nb232nf_iyH2(%rsp)
        movapd %xmm5,nb232nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb232nf_vctot(%rsp)
        movapd %xmm4,nb232nf_Vvdwtot(%rsp)

        movq  nb232nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb232nf_pos(%rbp),%rsi
        movq  nb232nf_faction(%rbp),%rdi
        movq  nb232nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb232nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb232nf_ninner(%rsp),%ecx
        movl  %ecx,nb232nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb232nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel232nf_x86_64_sse2.nb232nf_unroll_loop
        jmp   _nb_kernel232nf_x86_64_sse2.nb232nf_checksingle
_nb_kernel232nf_x86_64_sse2.nb232nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb232nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb232nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb232nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movhpd (%rsi,%rbx,8),%xmm2
        movhpd 8(%rsi,%rbx,8),%xmm3
        movhpd 16(%rsi,%rbx,8),%xmm4
        movhpd 24(%rsi,%rbx,8),%xmm5
        movhpd 32(%rsi,%rbx,8),%xmm6
        movhpd 40(%rsi,%rbx,8),%xmm7
        movapd  %xmm2,nb232nf_jxO(%rsp)
        movapd  %xmm3,nb232nf_jyO(%rsp)
        movapd  %xmm4,nb232nf_jzO(%rsp)
        movapd  %xmm5,nb232nf_jxH1(%rsp)
        movapd  %xmm6,nb232nf_jyH1(%rsp)
        movapd  %xmm7,nb232nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 48(%rsi,%rbx,8),%xmm2
        movhpd 56(%rsi,%rbx,8),%xmm3
        movhpd 64(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb232nf_jxH2(%rsp)
        movapd  %xmm3,nb232nf_jyH2(%rsp)
        movapd  %xmm4,nb232nf_jzH2(%rsp)

        movapd nb232nf_ixO(%rsp),%xmm0
        movapd nb232nf_iyO(%rsp),%xmm1
        movapd nb232nf_izO(%rsp),%xmm2
        movapd nb232nf_ixO(%rsp),%xmm3
        movapd nb232nf_iyO(%rsp),%xmm4
        movapd nb232nf_izO(%rsp),%xmm5
        subpd  nb232nf_jxO(%rsp),%xmm0
        subpd  nb232nf_jyO(%rsp),%xmm1
        subpd  nb232nf_jzO(%rsp),%xmm2
        subpd  nb232nf_jxH1(%rsp),%xmm3
        subpd  nb232nf_jyH1(%rsp),%xmm4
        subpd  nb232nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqOO(%rsp)
        movapd %xmm3,nb232nf_rsqOH1(%rsp)

        movapd nb232nf_ixO(%rsp),%xmm0
        movapd nb232nf_iyO(%rsp),%xmm1
        movapd nb232nf_izO(%rsp),%xmm2
        movapd nb232nf_ixH1(%rsp),%xmm3
        movapd nb232nf_iyH1(%rsp),%xmm4
        movapd nb232nf_izH1(%rsp),%xmm5
        subpd  nb232nf_jxH2(%rsp),%xmm0
        subpd  nb232nf_jyH2(%rsp),%xmm1
        subpd  nb232nf_jzH2(%rsp),%xmm2
        subpd  nb232nf_jxO(%rsp),%xmm3
        subpd  nb232nf_jyO(%rsp),%xmm4
        subpd  nb232nf_jzO(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqOH2(%rsp)
        movapd %xmm3,nb232nf_rsqH1O(%rsp)

        movapd nb232nf_ixH1(%rsp),%xmm0
        movapd nb232nf_iyH1(%rsp),%xmm1
        movapd nb232nf_izH1(%rsp),%xmm2
        movapd nb232nf_ixH1(%rsp),%xmm3
        movapd nb232nf_iyH1(%rsp),%xmm4
        movapd nb232nf_izH1(%rsp),%xmm5
        subpd  nb232nf_jxH1(%rsp),%xmm0
        subpd  nb232nf_jyH1(%rsp),%xmm1
        subpd  nb232nf_jzH1(%rsp),%xmm2
        subpd  nb232nf_jxH2(%rsp),%xmm3
        subpd  nb232nf_jyH2(%rsp),%xmm4
        subpd  nb232nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqH1H1(%rsp)
        movapd %xmm3,nb232nf_rsqH1H2(%rsp)

        movapd nb232nf_ixH2(%rsp),%xmm0
        movapd nb232nf_iyH2(%rsp),%xmm1
        movapd nb232nf_izH2(%rsp),%xmm2
        movapd nb232nf_ixH2(%rsp),%xmm3
        movapd nb232nf_iyH2(%rsp),%xmm4
        movapd nb232nf_izH2(%rsp),%xmm5
        subpd  nb232nf_jxO(%rsp),%xmm0
        subpd  nb232nf_jyO(%rsp),%xmm1
        subpd  nb232nf_jzO(%rsp),%xmm2
        subpd  nb232nf_jxH1(%rsp),%xmm3
        subpd  nb232nf_jyH1(%rsp),%xmm4
        subpd  nb232nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqH2O(%rsp)
        movapd %xmm4,nb232nf_rsqH2H1(%rsp)

        movapd nb232nf_ixH2(%rsp),%xmm0
        movapd nb232nf_iyH2(%rsp),%xmm1
        movapd nb232nf_izH2(%rsp),%xmm2
        subpd  nb232nf_jxH2(%rsp),%xmm0
        subpd  nb232nf_jyH2(%rsp),%xmm1
        subpd  nb232nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb232nf_rsqH2H2(%rsp)

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
        movapd  nb232nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb232nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvH2H2(%rsp)
        movapd %xmm5,nb232nf_rinvH2H1(%rsp)

        movapd nb232nf_rsqOO(%rsp),%xmm0
        movapd nb232nf_rsqOH1(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb232nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb232nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb232nf_rinvOO(%rsp)
        movapd %xmm5,nb232nf_rinvOH1(%rsp)

        movapd nb232nf_rsqOH2(%rsp),%xmm0
        movapd nb232nf_rsqH1O(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb232nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvOH2(%rsp)
        movapd %xmm5,nb232nf_rinvH1O(%rsp)

        movapd nb232nf_rsqH1H1(%rsp),%xmm0
        movapd nb232nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb232nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvH1H1(%rsp)
        movapd %xmm5,nb232nf_rinvH1H2(%rsp)

        movapd nb232nf_rsqH2O(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb232nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb232nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb232nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb232nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb232nf_rinvH2O(%rsp)

        ## start with OO interaction    
        movapd nb232nf_rinvOO(%rsp),%xmm0
        movapd %xmm0,%xmm7              ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        mulpd  nb232nf_rsqOO(%rsp),%xmm5        ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6              ## xmm6=rinv+ krsq 
        subpd  nb232nf_crf(%rsp),%xmm6          ## rinv+krsq-crf        
        mulpd  nb232nf_qqOO(%rsp),%xmm6         ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb232nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb232nf_vctot(%rsp)

        ## LJ table interaction
        movapd nb232nf_rsqOO(%rsp),%xmm4

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb232nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8 

        movq nb232nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx

        ## dispersion 
        movlpd (%rsi,%rax,8),%xmm4      ## Y1   
        movlpd (%rsi,%rbx,8),%xmm3      ## Y2 
        movhpd 8(%rsi,%rax,8),%xmm4     ## Y1 F1        
        movhpd 8(%rsi,%rbx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rax,8),%xmm6    ## G1
        movlpd 16(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 24(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 24(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb232nf_c6(%rsp),%xmm4
        mulpd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addpd  nb232nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb232nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rax,8),%xmm4    ## Y1   
        movlpd 32(%rsi,%rbx,8),%xmm3    ## Y2 
        movhpd 40(%rsi,%rax,8),%xmm4    ## Y1 F1        
        movhpd 40(%rsi,%rbx,8),%xmm3    ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rax,8),%xmm6    ## G1
        movlpd 48(%rsi,%rbx,8),%xmm3    ## G2
        movhpd 56(%rsi,%rax,8),%xmm6    ## G1 H1        
        movhpd 56(%rsi,%rbx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb232nf_c12(%rsp),%xmm4
        mulpd  %xmm4,%xmm5

        addpd  nb232nf_Vvdwtot(%rsp),%xmm5
        movapd %xmm5,nb232nf_Vvdwtot(%rsp)

        ## O-H1 interaction 
        movapd nb232nf_rinvOH1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqOH1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd  %xmm0,%xmm0
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 

        movapd nb232nf_vctot(%rsp),%xmm6
        addpd  %xmm4,%xmm6

        ## O-H2 interaction  
        movapd nb232nf_rinvOH2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqOH2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-O interaction 
        movapd nb232nf_rinvH1O(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH1O(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H1 interaction 
        movapd nb232nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH1H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H2 interaction 
        movapd nb232nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH1H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-O interaction 
        movapd nb232nf_rinvH2O(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH2O(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H1 interaction 
        movapd nb232nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH2H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H2 interaction 
        movapd nb232nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb232nf_rsqH2H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb232nf_crf(%rsp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        movapd %xmm6,nb232nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb232nf_innerk(%rsp)
        jl    _nb_kernel232nf_x86_64_sse2.nb232nf_checksingle
        jmp   _nb_kernel232nf_x86_64_sse2.nb232nf_unroll_loop
_nb_kernel232nf_x86_64_sse2.nb232nf_checksingle: 
        movl  nb232nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel232nf_x86_64_sse2.nb232nf_dosingle
        jmp   _nb_kernel232nf_x86_64_sse2.nb232nf_updateouterdata
_nb_kernel232nf_x86_64_sse2.nb232nf_dosingle: 
        movq  nb232nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb232nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb232nf_jxO(%rsp)
        movapd  %xmm3,nb232nf_jyO(%rsp)
        movapd  %xmm4,nb232nf_jzO(%rsp)
        movapd  %xmm5,nb232nf_jxH1(%rsp)
        movapd  %xmm6,nb232nf_jyH1(%rsp)
        movapd  %xmm7,nb232nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb232nf_jxH2(%rsp)
        movapd  %xmm3,nb232nf_jyH2(%rsp)
        movapd  %xmm4,nb232nf_jzH2(%rsp)

        movapd nb232nf_ixO(%rsp),%xmm0
        movapd nb232nf_iyO(%rsp),%xmm1
        movapd nb232nf_izO(%rsp),%xmm2
        movapd nb232nf_ixO(%rsp),%xmm3
        movapd nb232nf_iyO(%rsp),%xmm4
        movapd nb232nf_izO(%rsp),%xmm5
        subsd  nb232nf_jxO(%rsp),%xmm0
        subsd  nb232nf_jyO(%rsp),%xmm1
        subsd  nb232nf_jzO(%rsp),%xmm2
        subsd  nb232nf_jxH1(%rsp),%xmm3
        subsd  nb232nf_jyH1(%rsp),%xmm4
        subsd  nb232nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqOO(%rsp)
        movapd %xmm3,nb232nf_rsqOH1(%rsp)

        movapd nb232nf_ixO(%rsp),%xmm0
        movapd nb232nf_iyO(%rsp),%xmm1
        movapd nb232nf_izO(%rsp),%xmm2
        movapd nb232nf_ixH1(%rsp),%xmm3
        movapd nb232nf_iyH1(%rsp),%xmm4
        movapd nb232nf_izH1(%rsp),%xmm5
        subsd  nb232nf_jxH2(%rsp),%xmm0
        subsd  nb232nf_jyH2(%rsp),%xmm1
        subsd  nb232nf_jzH2(%rsp),%xmm2
        subsd  nb232nf_jxO(%rsp),%xmm3
        subsd  nb232nf_jyO(%rsp),%xmm4
        subsd  nb232nf_jzO(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqOH2(%rsp)
        movapd %xmm3,nb232nf_rsqH1O(%rsp)

        movapd nb232nf_ixH1(%rsp),%xmm0
        movapd nb232nf_iyH1(%rsp),%xmm1
        movapd nb232nf_izH1(%rsp),%xmm2
        movapd nb232nf_ixH1(%rsp),%xmm3
        movapd nb232nf_iyH1(%rsp),%xmm4
        movapd nb232nf_izH1(%rsp),%xmm5
        subsd  nb232nf_jxH1(%rsp),%xmm0
        subsd  nb232nf_jyH1(%rsp),%xmm1
        subsd  nb232nf_jzH1(%rsp),%xmm2
        subsd  nb232nf_jxH2(%rsp),%xmm3
        subsd  nb232nf_jyH2(%rsp),%xmm4
        subsd  nb232nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqH1H1(%rsp)
        movapd %xmm3,nb232nf_rsqH1H2(%rsp)

        movapd nb232nf_ixH2(%rsp),%xmm0
        movapd nb232nf_iyH2(%rsp),%xmm1
        movapd nb232nf_izH2(%rsp),%xmm2
        movapd nb232nf_ixH2(%rsp),%xmm3
        movapd nb232nf_iyH2(%rsp),%xmm4
        movapd nb232nf_izH2(%rsp),%xmm5
        subsd  nb232nf_jxO(%rsp),%xmm0
        subsd  nb232nf_jyO(%rsp),%xmm1
        subsd  nb232nf_jzO(%rsp),%xmm2
        subsd  nb232nf_jxH1(%rsp),%xmm3
        subsd  nb232nf_jyH1(%rsp),%xmm4
        subsd  nb232nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb232nf_rsqH2O(%rsp)
        movapd %xmm4,nb232nf_rsqH2H1(%rsp)

        movapd nb232nf_ixH2(%rsp),%xmm0
        movapd nb232nf_iyH2(%rsp),%xmm1
        movapd nb232nf_izH2(%rsp),%xmm2
        subsd  nb232nf_jxH2(%rsp),%xmm0
        subsd  nb232nf_jyH2(%rsp),%xmm1
        subsd  nb232nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb232nf_rsqH2H2(%rsp)

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
        movapd  nb232nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb232nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvH2H2(%rsp)
        movapd %xmm5,nb232nf_rinvH2H1(%rsp)

        movapd nb232nf_rsqOO(%rsp),%xmm0
        movapd nb232nf_rsqOH1(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb232nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb232nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb232nf_rinvOO(%rsp)
        movapd %xmm5,nb232nf_rinvOH1(%rsp)

        movapd nb232nf_rsqOH2(%rsp),%xmm0
        movapd nb232nf_rsqH1O(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb232nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvOH2(%rsp)
        movapd %xmm5,nb232nf_rinvH1O(%rsp)

        movapd nb232nf_rsqH1H1(%rsp),%xmm0
        movapd nb232nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb232nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb232nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb232nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb232nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb232nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb232nf_rinvH1H1(%rsp)
        movapd %xmm5,nb232nf_rinvH1H2(%rsp)

        movapd nb232nf_rsqH2O(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb232nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb232nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb232nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb232nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb232nf_rinvH2O(%rsp)

        ## start with OO interaction    
        movsd nb232nf_rinvOO(%rsp),%xmm0
        movsd %xmm0,%xmm7               ## xmm7=rinv 
        movsd nb232nf_krf(%rsp),%xmm5
        mulsd  nb232nf_rsqOO(%rsp),%xmm5        ## xmm5=krsq 
        movsd %xmm5,%xmm6               ## krsq 
        addsd  %xmm7,%xmm6              ## xmm6=rinv+ krsq 
        subsd  nb232nf_crf(%rsp),%xmm6          ## rinv+krsq-crf        
        mulsd  nb232nf_qqOO(%rsp),%xmm6         ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb232nf_vctot(%rsp),%xmm6
        movsd %xmm6,nb232nf_vctot(%rsp)

        ## LJ table interaction
        movsd nb232nf_rsqOO(%rsp),%xmm4

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb232nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%ebx    ## mm6 = lu idx 
        cvtsi2sd %ebx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%ebx

        movq nb232nf_VFtab(%rbp),%rsi

        ## dispersion 
        movlpd (%rsi,%rbx,8),%xmm4      ## Y1   
        movhpd 8(%rsi,%rbx,8),%xmm4     ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 24(%rsi,%rbx,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## dispersion table ready, in xmm4-xmm7         
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb232nf_c6(%rsp),%xmm4
        mulsd  %xmm4,%xmm5       ## Vvdw6 

        ## Update Vvdwtot directly 
        addsd  nb232nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb232nf_Vvdwtot(%rsp)

        ## repulsion 
        movlpd 32(%rsi,%rbx,8),%xmm4    ## Y1   
        movhpd 40(%rsi,%rbx,8),%xmm4    ## Y1 F1        

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 48(%rsi,%rbx,8),%xmm6    ## G1
        movhpd 56(%rsi,%rbx,8),%xmm6    ## G1 H1        

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 

        ## table ready, in xmm4-xmm7    
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movsd nb232nf_c12(%rsp),%xmm4
        mulsd  %xmm4,%xmm5

        addsd  nb232nf_Vvdwtot(%rsp),%xmm5
        movsd %xmm5,nb232nf_Vvdwtot(%rsp)

        ## O-H1 interaction 
        movapd nb232nf_rinvOH1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqOH1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 

        movsd  nb232nf_vctot(%rsp),%xmm6
        addsd  %xmm4,%xmm6

        ## O-H2 interaction  
        movapd nb232nf_rinvOH2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqOH2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-O interaction 
        movapd nb232nf_rinvH1O(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH1O(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H1 interaction 
        movapd nb232nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH1H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H2 interaction 
        movapd nb232nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH1H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-O interaction 
        movapd nb232nf_rinvH2O(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH2O(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb232nf_qqOH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H1 interaction 
        movapd nb232nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH2H1(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H2 interaction 
        movapd nb232nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb232nf_krf(%rsp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb232nf_rsqH2H2(%rsp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb232nf_crf(%rsp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb232nf_qqHH(%rsp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        movsd %xmm6,nb232nf_vctot(%rsp)

_nb_kernel232nf_x86_64_sse2.nb232nf_updateouterdata: 
        ## get n from stack
        movl nb232nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb232nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb232nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb232nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb232nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb232nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb232nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel232nf_x86_64_sse2.nb232nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb232nf_n(%rsp)
        jmp _nb_kernel232nf_x86_64_sse2.nb232nf_outer
_nb_kernel232nf_x86_64_sse2.nb232nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb232nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel232nf_x86_64_sse2.nb232nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel232nf_x86_64_sse2.nb232nf_threadloop
_nb_kernel232nf_x86_64_sse2.nb232nf_end: 
        movl nb232nf_nouter(%rsp),%eax
        movl nb232nf_ninner(%rsp),%ebx
        movq nb232nf_outeriter(%rbp),%rcx
        movq nb232nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $872,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





