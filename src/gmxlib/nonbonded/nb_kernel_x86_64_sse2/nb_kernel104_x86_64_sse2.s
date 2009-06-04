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






.globl nb_kernel104_x86_64_sse2
.globl _nb_kernel104_x86_64_sse2
nb_kernel104_x86_64_sse2:       
_nb_kernel104_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use        
.set nb104_ixM, 0
.set nb104_iyM, 16
.set nb104_izM, 32
.set nb104_ixH1, 48
.set nb104_iyH1, 64
.set nb104_izH1, 80
.set nb104_ixH2, 96
.set nb104_iyH2, 112
.set nb104_izH2, 128
.set nb104_jxM, 144
.set nb104_jyM, 160
.set nb104_jzM, 176
.set nb104_jxH1, 192
.set nb104_jyH1, 208
.set nb104_jzH1, 224
.set nb104_jxH2, 240
.set nb104_jyH2, 256
.set nb104_jzH2, 272
.set nb104_dxMM, 288
.set nb104_dyMM, 304
.set nb104_dzMM, 320
.set nb104_dxMH1, 336
.set nb104_dyMH1, 352
.set nb104_dzMH1, 368
.set nb104_dxMH2, 384
.set nb104_dyMH2, 400
.set nb104_dzMH2, 416
.set nb104_dxH1M, 432
.set nb104_dyH1M, 448
.set nb104_dzH1M, 464
.set nb104_dxH1H1, 480
.set nb104_dyH1H1, 496
.set nb104_dzH1H1, 512
.set nb104_dxH1H2, 528
.set nb104_dyH1H2, 544
.set nb104_dzH1H2, 560
.set nb104_dxH2M, 576
.set nb104_dyH2M, 592
.set nb104_dzH2M, 608
.set nb104_dxH2H1, 624
.set nb104_dyH2H1, 640
.set nb104_dzH2H1, 656
.set nb104_dxH2H2, 672
.set nb104_dyH2H2, 688
.set nb104_dzH2H2, 704
.set nb104_qqMM, 720
.set nb104_qqMH, 736
.set nb104_qqHH, 752
.set nb104_vctot, 768
.set nb104_fixM, 784
.set nb104_fiyM, 800
.set nb104_fizM, 816
.set nb104_fixH1, 832
.set nb104_fiyH1, 848
.set nb104_fizH1, 864
.set nb104_fixH2, 880
.set nb104_fiyH2, 896
.set nb104_fizH2, 912
.set nb104_fjxM, 928
.set nb104_fjyM, 944
.set nb104_fjzM, 960
.set nb104_fjxH1, 976
.set nb104_fjyH1, 992
.set nb104_fjzH1, 1008
.set nb104_fjxH2, 1024
.set nb104_fjyH2, 1040
.set nb104_fjzH2, 1056
.set nb104_half, 1072
.set nb104_three, 1088
.set nb104_rsqMM, 1104
.set nb104_rsqMH1, 1120
.set nb104_rsqMH2, 1136
.set nb104_rsqH1M, 1152
.set nb104_rsqH1H1, 1168
.set nb104_rsqH1H2, 1184
.set nb104_rsqH2M, 1200
.set nb104_rsqH2H1, 1216
.set nb104_rsqH2H2, 1232
.set nb104_rinvMM, 1248
.set nb104_rinvMH1, 1264
.set nb104_rinvMH2, 1280
.set nb104_rinvH1M, 1296
.set nb104_rinvH1H1, 1312
.set nb104_rinvH1H2, 1328
.set nb104_rinvH2M, 1344
.set nb104_rinvH2H1, 1360
.set nb104_rinvH2H2, 1376
.set nb104_nri, 1392
.set nb104_innerjjnr, 1400
.set nb104_iinr, 1408
.set nb104_jindex, 1416
.set nb104_jjnr, 1424
.set nb104_shift, 1432
.set nb104_shiftvec, 1440
.set nb104_facel, 1448
.set nb104_is3, 1456
.set nb104_ii3, 1460
.set nb104_innerk, 1464
.set nb104_n, 1468
.set nb104_nn1, 1472
.set nb104_nouter, 1476
.set nb104_ninner, 1480
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1496,%rsp         ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb104_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb104_half(%rsp)
        movl %ebx,nb104_half+4(%rsp)
        movsd nb104_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb104_half(%rsp)
        movapd %xmm3,nb104_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb104_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb104_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3     ## qM 
        movsd %xmm3,%xmm4               ## qM 
        movsd 8(%rdx,%rbx,8),%xmm5      ## qH 
        movq nb104_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb104_facel(%rsp),%xmm6   ## facel 
        mulsd  %xmm3,%xmm3              ## qM*qM 
        mulsd  %xmm5,%xmm4              ## qM*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104_qqMM(%rsp)
        movapd %xmm4,nb104_qqMH(%rsp)
        movapd %xmm5,nb104_qqHH(%rsp)

_nb_kernel104_x86_64_sse2.nb104_threadloop: 
        movq  nb104_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel104_x86_64_sse2.nb104_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104_x86_64_sse2.nb104_spinlock

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
        jg  _nb_kernel104_x86_64_sse2.nb104_outerstart
        jmp _nb_kernel104_x86_64_sse2.nb104_end

_nb_kernel104_x86_64_sse2.nb104_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104_nouter(%rsp),%ebx
        movl %ebx,nb104_nouter(%rsp)

_nb_kernel104_x86_64_sse2.nb104_outer: 
        movq  nb104_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb104_is3(%rsp)      ## store is3 

        movq  nb104_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb104_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb104_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb104_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104_ixH1(%rsp)
        movapd %xmm4,nb104_iyH1(%rsp)
        movapd %xmm5,nb104_izH1(%rsp)

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
        movapd %xmm0,nb104_ixH2(%rsp)
        movapd %xmm1,nb104_iyH2(%rsp)
        movapd %xmm2,nb104_izH2(%rsp)
        movapd %xmm3,nb104_ixM(%rsp)
        movapd %xmm4,nb104_iyM(%rsp)
        movapd %xmm5,nb104_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb104_vctot(%rsp)
        movapd %xmm4,nb104_fixM(%rsp)
        movapd %xmm4,nb104_fiyM(%rsp)
        movapd %xmm4,nb104_fizM(%rsp)
        movapd %xmm4,nb104_fixH1(%rsp)
        movapd %xmm4,nb104_fiyH1(%rsp)
        movapd %xmm4,nb104_fizH1(%rsp)
        movapd %xmm4,nb104_fixH2(%rsp)
        movapd %xmm4,nb104_fiyH2(%rsp)
        movapd %xmm4,nb104_fizH2(%rsp)

        movq  nb104_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb104_pos(%rbp),%rsi
        movq  nb104_faction(%rbp),%rdi
        movq  nb104_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb104_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb104_ninner(%rsp),%ecx
        movl  %ecx,nb104_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb104_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel104_x86_64_sse2.nb104_unroll_loop
        jmp   _nb_kernel104_x86_64_sse2.nb104_checksingle
_nb_kernel104_x86_64_sse2.nb104_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb104_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb104_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb104_pos(%rbp),%rsi        ## base of pos[] 

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

    subpd nb104_ixH1(%rsp),%xmm0
    subpd nb104_iyH1(%rsp),%xmm1
    subpd nb104_izH1(%rsp),%xmm2
    subpd nb104_ixH2(%rsp),%xmm3
    subpd nb104_iyH2(%rsp),%xmm4
    subpd nb104_izH2(%rsp),%xmm5
    subpd nb104_ixM(%rsp),%xmm6
    subpd nb104_iyM(%rsp),%xmm7
    subpd nb104_izM(%rsp),%xmm8

        movapd %xmm0,nb104_dxH1H1(%rsp)
        movapd %xmm1,nb104_dyH1H1(%rsp)
        movapd %xmm2,nb104_dzH1H1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH2H1(%rsp)
        movapd %xmm4,nb104_dyH2H1(%rsp)
        movapd %xmm5,nb104_dzH2H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb104_dxMH1(%rsp)
        movapd %xmm7,nb104_dyMH1(%rsp)
        movapd %xmm8,nb104_dzMH1(%rsp)
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

        movapd  nb104_three(%rsp),%xmm9
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

        movapd  nb104_half(%rsp),%xmm15
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

        movapd  nb104_three(%rsp),%xmm1
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

        movapd  nb104_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H1 
        mulpd   %xmm15,%xmm10 ##   rinvH2H1
    mulpd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb104_qqHH(%rsp),%xmm0
    mulpd  nb104_qqHH(%rsp),%xmm1
    mulpd  nb104_qqMH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb104_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb104_vctot(%rsp)

    ## move j H1 forces to xmm0-xmm2
        movq  nb104_faction(%rbp),%rdi
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

        mulpd nb104_dxH1H1(%rsp),%xmm7
        mulpd nb104_dyH1H1(%rsp),%xmm8
        mulpd nb104_dzH1H1(%rsp),%xmm9
        mulpd nb104_dxH2H1(%rsp),%xmm10
        mulpd nb104_dyH2H1(%rsp),%xmm11
        mulpd nb104_dzH2H1(%rsp),%xmm12
        mulpd nb104_dxMH1(%rsp),%xmm13
        mulpd nb104_dyMH1(%rsp),%xmm14
        mulpd nb104_dzMH1(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb104_fixH1(%rsp),%xmm7
    addpd nb104_fiyH1(%rsp),%xmm8
    addpd nb104_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb104_fixH2(%rsp),%xmm10
    addpd nb104_fiyH2(%rsp),%xmm11
    addpd nb104_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb104_fixM(%rsp),%xmm13
    addpd nb104_fiyM(%rsp),%xmm14
    addpd nb104_fizM(%rsp),%xmm15

    movapd %xmm7,nb104_fixH1(%rsp)
    movapd %xmm8,nb104_fiyH1(%rsp)
    movapd %xmm9,nb104_fizH1(%rsp)
    movapd %xmm10,nb104_fixH2(%rsp)
    movapd %xmm11,nb104_fiyH2(%rsp)
    movapd %xmm12,nb104_fizH2(%rsp)
    movapd %xmm13,nb104_fixM(%rsp)
    movapd %xmm14,nb104_fiyM(%rsp)
    movapd %xmm15,nb104_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
        movq  nb104_pos(%rbp),%rsi
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

    subpd nb104_ixH1(%rsp),%xmm0
    subpd nb104_iyH1(%rsp),%xmm1
    subpd nb104_izH1(%rsp),%xmm2
    subpd nb104_ixH2(%rsp),%xmm3
    subpd nb104_iyH2(%rsp),%xmm4
    subpd nb104_izH2(%rsp),%xmm5
    subpd nb104_ixM(%rsp),%xmm6
    subpd nb104_iyM(%rsp),%xmm7
    subpd nb104_izM(%rsp),%xmm8

        movapd %xmm0,nb104_dxH1H2(%rsp)
        movapd %xmm1,nb104_dyH1H2(%rsp)
        movapd %xmm2,nb104_dzH1H2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH2H2(%rsp)
        movapd %xmm4,nb104_dyH2H2(%rsp)
        movapd %xmm5,nb104_dzH2H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb104_dxMH2(%rsp)
        movapd %xmm7,nb104_dyMH2(%rsp)
        movapd %xmm8,nb104_dzMH2(%rsp)
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

        movapd  nb104_three(%rsp),%xmm9
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

        movapd  nb104_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH1H2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb104_three(%rsp),%xmm1
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

        movapd  nb104_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H2
        mulpd   %xmm15,%xmm10 ##   rinvH2H2
    mulpd   %xmm15,%xmm11 ##   rinvMH2

        ## H2 interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb104_qqHH(%rsp),%xmm0
    mulpd  nb104_qqHH(%rsp),%xmm1
    mulpd  nb104_qqMH(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb104_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb104_vctot(%rsp)

    ## move j H2 forces to xmm0-xmm2
        movq  nb104_faction(%rbp),%rdi
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

        mulpd nb104_dxH1H2(%rsp),%xmm7
        mulpd nb104_dyH1H2(%rsp),%xmm8
        mulpd nb104_dzH1H2(%rsp),%xmm9
        mulpd nb104_dxH2H2(%rsp),%xmm10
        mulpd nb104_dyH2H2(%rsp),%xmm11
        mulpd nb104_dzH2H2(%rsp),%xmm12
        mulpd nb104_dxMH2(%rsp),%xmm13
        mulpd nb104_dyMH2(%rsp),%xmm14
        mulpd nb104_dzMH2(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb104_fixH1(%rsp),%xmm7
    addpd nb104_fiyH1(%rsp),%xmm8
    addpd nb104_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb104_fixH2(%rsp),%xmm10
    addpd nb104_fiyH2(%rsp),%xmm11
    addpd nb104_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb104_fixM(%rsp),%xmm13
    addpd nb104_fiyM(%rsp),%xmm14
    addpd nb104_fizM(%rsp),%xmm15

    movapd %xmm7,nb104_fixH1(%rsp)
    movapd %xmm8,nb104_fiyH1(%rsp)
    movapd %xmm9,nb104_fizH1(%rsp)
    movapd %xmm10,nb104_fixH2(%rsp)
    movapd %xmm11,nb104_fiyH2(%rsp)
    movapd %xmm12,nb104_fizH2(%rsp)
    movapd %xmm13,nb104_fixM(%rsp)
    movapd %xmm14,nb104_fiyM(%rsp)
    movapd %xmm15,nb104_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## move j M coordinates to local temp variables 
        movq  nb104_pos(%rbp),%rsi
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

    subpd nb104_ixH1(%rsp),%xmm0
    subpd nb104_iyH1(%rsp),%xmm1
    subpd nb104_izH1(%rsp),%xmm2
    subpd nb104_ixH2(%rsp),%xmm3
    subpd nb104_iyH2(%rsp),%xmm4
    subpd nb104_izH2(%rsp),%xmm5
    subpd nb104_ixM(%rsp),%xmm6
    subpd nb104_iyM(%rsp),%xmm7
    subpd nb104_izM(%rsp),%xmm8

        movapd %xmm0,nb104_dxH1M(%rsp)
        movapd %xmm1,nb104_dyH1M(%rsp)
        movapd %xmm2,nb104_dzH1M(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH2M(%rsp)
        movapd %xmm4,nb104_dyH2M(%rsp)
        movapd %xmm5,nb104_dzH2M(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb104_dxMM(%rsp)
        movapd %xmm7,nb104_dyMM(%rsp)
        movapd %xmm8,nb104_dzMM(%rsp)
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

        movapd  nb104_three(%rsp),%xmm9
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

        movapd  nb104_half(%rsp),%xmm15
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

        movapd  nb104_three(%rsp),%xmm1
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

        movapd  nb104_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1M
        mulpd   %xmm15,%xmm10 ##   rinvH2M
    mulpd   %xmm15,%xmm11 ##   rinvMM

        ## M interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb104_qqMH(%rsp),%xmm0
    mulpd  nb104_qqMH(%rsp),%xmm1
    mulpd  nb104_qqMM(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb104_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb104_vctot(%rsp)

    ## move j M forces to xmm0-xmm2
        movq  nb104_faction(%rbp),%rdi
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

        mulpd nb104_dxH1M(%rsp),%xmm7
        mulpd nb104_dyH1M(%rsp),%xmm8
        mulpd nb104_dzH1M(%rsp),%xmm9
        mulpd nb104_dxH2M(%rsp),%xmm10
        mulpd nb104_dyH2M(%rsp),%xmm11
        mulpd nb104_dzH2M(%rsp),%xmm12
        mulpd nb104_dxMM(%rsp),%xmm13
        mulpd nb104_dyMM(%rsp),%xmm14
        mulpd nb104_dzMM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb104_fixH1(%rsp),%xmm7
    addpd nb104_fiyH1(%rsp),%xmm8
    addpd nb104_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb104_fixH2(%rsp),%xmm10
    addpd nb104_fiyH2(%rsp),%xmm11
    addpd nb104_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb104_fixM(%rsp),%xmm13
    addpd nb104_fiyM(%rsp),%xmm14
    addpd nb104_fizM(%rsp),%xmm15

    movapd %xmm7,nb104_fixH1(%rsp)
    movapd %xmm8,nb104_fiyH1(%rsp)
    movapd %xmm9,nb104_fizH1(%rsp)
    movapd %xmm10,nb104_fixH2(%rsp)
    movapd %xmm11,nb104_fiyH2(%rsp)
    movapd %xmm12,nb104_fizH2(%rsp)
    movapd %xmm13,nb104_fixM(%rsp)
    movapd %xmm14,nb104_fiyM(%rsp)
    movapd %xmm15,nb104_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movlpd %xmm0,72(%rdi,%rax,8)
        movlpd %xmm1,80(%rdi,%rax,8)
        movlpd %xmm2,88(%rdi,%rax,8)
        movhpd %xmm0,72(%rdi,%rbx,8)
        movhpd %xmm1,80(%rdi,%rbx,8)
        movhpd %xmm2,88(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb104_innerk(%rsp)
        jl    _nb_kernel104_x86_64_sse2.nb104_checksingle
        jmp   _nb_kernel104_x86_64_sse2.nb104_unroll_loop
_nb_kernel104_x86_64_sse2.nb104_checksingle: 
        movl  nb104_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel104_x86_64_sse2.nb104_dosingle
        jmp   _nb_kernel104_x86_64_sse2.nb104_updateouterdata
_nb_kernel104_x86_64_sse2.nb104_dosingle: 
        movq  nb104_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb104_pos(%rbp),%rsi
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

    subsd nb104_ixH1(%rsp),%xmm0
    subsd nb104_iyH1(%rsp),%xmm1
    subsd nb104_izH1(%rsp),%xmm2
    subsd nb104_ixH2(%rsp),%xmm3
    subsd nb104_iyH2(%rsp),%xmm4
    subsd nb104_izH2(%rsp),%xmm5
    subsd nb104_ixM(%rsp),%xmm6
    subsd nb104_iyM(%rsp),%xmm7
    subsd nb104_izM(%rsp),%xmm8

        movsd %xmm0,nb104_dxH1H1(%rsp)
        movsd %xmm1,nb104_dyH1H1(%rsp)
        movsd %xmm2,nb104_dzH1H1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb104_dxH2H1(%rsp)
        movsd %xmm4,nb104_dyH2H1(%rsp)
        movsd %xmm5,nb104_dzH2H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb104_dxMH1(%rsp)
        movsd %xmm7,nb104_dyMH1(%rsp)
        movsd %xmm8,nb104_dzMH1(%rsp)
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

        movsd  nb104_three(%rsp),%xmm9
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

        movsd  nb104_half(%rsp),%xmm15
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

        movsd  nb104_three(%rsp),%xmm1
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

        movsd  nb104_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H1 
        mulsd   %xmm15,%xmm10 ##   rinvH2H1
    mulsd   %xmm15,%xmm11 ##   rinvMH1

        ## H1 interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb104_qqHH(%rsp),%xmm0
    mulsd  nb104_qqHH(%rsp),%xmm1
    mulsd  nb104_qqMH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb104_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb104_vctot(%rsp)

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

        mulsd nb104_dxH1H1(%rsp),%xmm7
        mulsd nb104_dyH1H1(%rsp),%xmm8
        mulsd nb104_dzH1H1(%rsp),%xmm9
        mulsd nb104_dxH2H1(%rsp),%xmm10
        mulsd nb104_dyH2H1(%rsp),%xmm11
        mulsd nb104_dzH2H1(%rsp),%xmm12
        mulsd nb104_dxMH1(%rsp),%xmm13
        mulsd nb104_dyMH1(%rsp),%xmm14
        mulsd nb104_dzMH1(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb104_fixH1(%rsp),%xmm7
    addsd nb104_fiyH1(%rsp),%xmm8
    addsd nb104_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb104_fixH2(%rsp),%xmm10
    addsd nb104_fiyH2(%rsp),%xmm11
    addsd nb104_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb104_fixM(%rsp),%xmm13
    addsd nb104_fiyM(%rsp),%xmm14
    addsd nb104_fizM(%rsp),%xmm15

    movsd %xmm7,nb104_fixH1(%rsp)
    movsd %xmm8,nb104_fiyH1(%rsp)
    movsd %xmm9,nb104_fizH1(%rsp)
    movsd %xmm10,nb104_fixH2(%rsp)
    movsd %xmm11,nb104_fiyH2(%rsp)
    movsd %xmm12,nb104_fizH2(%rsp)
    movsd %xmm13,nb104_fixM(%rsp)
    movsd %xmm14,nb104_fiyM(%rsp)
    movsd %xmm15,nb104_fizM(%rsp)

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

    subsd nb104_ixH1(%rsp),%xmm0
    subsd nb104_iyH1(%rsp),%xmm1
    subsd nb104_izH1(%rsp),%xmm2
    subsd nb104_ixH2(%rsp),%xmm3
    subsd nb104_iyH2(%rsp),%xmm4
    subsd nb104_izH2(%rsp),%xmm5
    subsd nb104_ixM(%rsp),%xmm6
    subsd nb104_iyM(%rsp),%xmm7
    subsd nb104_izM(%rsp),%xmm8

        movsd %xmm0,nb104_dxH1H2(%rsp)
        movsd %xmm1,nb104_dyH1H2(%rsp)
        movsd %xmm2,nb104_dzH1H2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb104_dxH2H2(%rsp)
        movsd %xmm4,nb104_dyH2H2(%rsp)
        movsd %xmm5,nb104_dzH2H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb104_dxMH2(%rsp)
        movsd %xmm7,nb104_dyMH2(%rsp)
        movsd %xmm8,nb104_dzMH2(%rsp)
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

        movsd  nb104_three(%rsp),%xmm9
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

        movsd  nb104_half(%rsp),%xmm15
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

        movsd  nb104_three(%rsp),%xmm1
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

        movsd  nb104_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H2
        mulsd   %xmm15,%xmm10 ##   rinvH2H2
    mulsd   %xmm15,%xmm11 ##   rinvMH2

        ## H2 interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb104_qqHH(%rsp),%xmm0
    mulsd  nb104_qqHH(%rsp),%xmm1
    mulsd  nb104_qqMH(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb104_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb104_vctot(%rsp)

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

        mulsd nb104_dxH1H2(%rsp),%xmm7
        mulsd nb104_dyH1H2(%rsp),%xmm8
        mulsd nb104_dzH1H2(%rsp),%xmm9
        mulsd nb104_dxH2H2(%rsp),%xmm10
        mulsd nb104_dyH2H2(%rsp),%xmm11
        mulsd nb104_dzH2H2(%rsp),%xmm12
        mulsd nb104_dxMH2(%rsp),%xmm13
        mulsd nb104_dyMH2(%rsp),%xmm14
        mulsd nb104_dzMH2(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb104_fixH1(%rsp),%xmm7
    addsd nb104_fiyH1(%rsp),%xmm8
    addsd nb104_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb104_fixH2(%rsp),%xmm10
    addsd nb104_fiyH2(%rsp),%xmm11
    addsd nb104_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb104_fixM(%rsp),%xmm13
    addsd nb104_fiyM(%rsp),%xmm14
    addsd nb104_fizM(%rsp),%xmm15

    movsd %xmm7,nb104_fixH1(%rsp)
    movsd %xmm8,nb104_fiyH1(%rsp)
    movsd %xmm9,nb104_fizH1(%rsp)
    movsd %xmm10,nb104_fixH2(%rsp)
    movsd %xmm11,nb104_fiyH2(%rsp)
    movsd %xmm12,nb104_fizH2(%rsp)
    movsd %xmm13,nb104_fixM(%rsp)
    movsd %xmm14,nb104_fiyM(%rsp)
    movsd %xmm15,nb104_fizM(%rsp)

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

    subsd nb104_ixH1(%rsp),%xmm0
    subsd nb104_iyH1(%rsp),%xmm1
    subsd nb104_izH1(%rsp),%xmm2
    subsd nb104_ixH2(%rsp),%xmm3
    subsd nb104_iyH2(%rsp),%xmm4
    subsd nb104_izH2(%rsp),%xmm5
    subsd nb104_ixM(%rsp),%xmm6
    subsd nb104_iyM(%rsp),%xmm7
    subsd nb104_izM(%rsp),%xmm8

        movsd %xmm0,nb104_dxH1M(%rsp)
        movsd %xmm1,nb104_dyH1M(%rsp)
        movsd %xmm2,nb104_dzH1M(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb104_dxH2M(%rsp)
        movsd %xmm4,nb104_dyH2M(%rsp)
        movsd %xmm5,nb104_dzH2M(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb104_dxMM(%rsp)
        movsd %xmm7,nb104_dyMM(%rsp)
        movsd %xmm8,nb104_dzMM(%rsp)
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

        movsd  nb104_three(%rsp),%xmm9
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

        movsd  nb104_half(%rsp),%xmm15
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

        movsd  nb104_three(%rsp),%xmm1
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

        movsd  nb104_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1M
        mulsd   %xmm15,%xmm10 ##   rinvH2M
    mulsd   %xmm15,%xmm11 ##   rinvMM

        ## M interactions 
    movsd %xmm9,%xmm0
    movsd %xmm10,%xmm1
    movsd %xmm11,%xmm2
    mulsd  %xmm9,%xmm9
    mulsd  %xmm10,%xmm10
    mulsd  %xmm11,%xmm11
    mulsd  nb104_qqMH(%rsp),%xmm0
    mulsd  nb104_qqMH(%rsp),%xmm1
    mulsd  nb104_qqMM(%rsp),%xmm2
    mulsd  %xmm0,%xmm9
    mulsd  %xmm1,%xmm10
    mulsd  %xmm2,%xmm11

    addsd nb104_vctot(%rsp),%xmm0
    addsd %xmm2,%xmm1
    addsd %xmm1,%xmm0
    movsd %xmm0,nb104_vctot(%rsp)

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

        mulsd nb104_dxH1M(%rsp),%xmm7
        mulsd nb104_dyH1M(%rsp),%xmm8
        mulsd nb104_dzH1M(%rsp),%xmm9
        mulsd nb104_dxH2M(%rsp),%xmm10
        mulsd nb104_dyH2M(%rsp),%xmm11
        mulsd nb104_dzH2M(%rsp),%xmm12
        mulsd nb104_dxMM(%rsp),%xmm13
        mulsd nb104_dyMM(%rsp),%xmm14
        mulsd nb104_dzMM(%rsp),%xmm15

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb104_fixH1(%rsp),%xmm7
    addsd nb104_fiyH1(%rsp),%xmm8
    addsd nb104_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb104_fixH2(%rsp),%xmm10
    addsd nb104_fiyH2(%rsp),%xmm11
    addsd nb104_fizH2(%rsp),%xmm12

    addsd %xmm13,%xmm0
    addsd %xmm14,%xmm1
    addsd %xmm15,%xmm2
    addsd nb104_fixM(%rsp),%xmm13
    addsd nb104_fiyM(%rsp),%xmm14
    addsd nb104_fizM(%rsp),%xmm15

    movsd %xmm7,nb104_fixH1(%rsp)
    movsd %xmm8,nb104_fiyH1(%rsp)
    movsd %xmm9,nb104_fizH1(%rsp)
    movsd %xmm10,nb104_fixH2(%rsp)
    movsd %xmm11,nb104_fiyH2(%rsp)
    movsd %xmm12,nb104_fizH2(%rsp)
    movsd %xmm13,nb104_fixM(%rsp)
    movsd %xmm14,nb104_fiyM(%rsp)
    movsd %xmm15,nb104_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movsd %xmm0,72(%rdi,%rax,8)
        movsd %xmm1,80(%rdi,%rax,8)
        movsd %xmm2,88(%rdi,%rax,8)

_nb_kernel104_x86_64_sse2.nb104_updateouterdata: 
        movl  nb104_ii3(%rsp),%ecx
        movq  nb104_faction(%rbp),%rdi
        movq  nb104_fshift(%rbp),%rsi
        movl  nb104_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb104_fixH1(%rsp),%xmm0
        movapd nb104_fiyH1(%rsp),%xmm1
        movapd nb104_fizH1(%rsp),%xmm2

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
        movapd nb104_fixH2(%rsp),%xmm0
        movapd nb104_fiyH2(%rsp),%xmm1
        movapd nb104_fizH2(%rsp),%xmm2

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
        movapd nb104_fixM(%rsp),%xmm0
        movapd nb104_fiyM(%rsp),%xmm1
        movapd nb104_fizM(%rsp),%xmm2

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
        movl nb104_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb104_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb104_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb104_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb104_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104_x86_64_sse2.nb104_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104_n(%rsp)
        jmp _nb_kernel104_x86_64_sse2.nb104_outer
_nb_kernel104_x86_64_sse2.nb104_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104_x86_64_sse2.nb104_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104_x86_64_sse2.nb104_threadloop
_nb_kernel104_x86_64_sse2.nb104_end: 
        movl nb104_nouter(%rsp),%eax
        movl nb104_ninner(%rsp),%ebx
        movq nb104_outeriter(%rbp),%rcx
        movq nb104_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1496,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel104nf_x86_64_sse2
.globl _nb_kernel104nf_x86_64_sse2
nb_kernel104nf_x86_64_sse2:     
_nb_kernel104nf_x86_64_sse2:    
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
        ## bottom of stack is cache-aligned for sse2 use 
.set nb104nf_ixM, 0
.set nb104nf_iyM, 16
.set nb104nf_izM, 32
.set nb104nf_ixH1, 48
.set nb104nf_iyH1, 64
.set nb104nf_izH1, 80
.set nb104nf_ixH2, 96
.set nb104nf_iyH2, 112
.set nb104nf_izH2, 128
.set nb104nf_jxM, 144
.set nb104nf_jyM, 160
.set nb104nf_jzM, 176
.set nb104nf_jxH1, 192
.set nb104nf_jyH1, 208
.set nb104nf_jzH1, 224
.set nb104nf_jxH2, 240
.set nb104nf_jyH2, 256
.set nb104nf_jzH2, 272
.set nb104nf_qqMM, 288
.set nb104nf_qqMH, 304
.set nb104nf_qqHH, 320
.set nb104nf_vctot, 336
.set nb104nf_half, 352
.set nb104nf_three, 368
.set nb104nf_rsqMM, 384
.set nb104nf_rsqMH1, 400
.set nb104nf_rsqMH2, 416
.set nb104nf_rsqH1M, 432
.set nb104nf_rsqH1H1, 448
.set nb104nf_rsqH1H2, 464
.set nb104nf_rsqH2M, 480
.set nb104nf_rsqH2H1, 496
.set nb104nf_rsqH2H2, 512
.set nb104nf_rinvMM, 528
.set nb104nf_rinvMH1, 544
.set nb104nf_rinvMH2, 560
.set nb104nf_rinvH1M, 576
.set nb104nf_rinvH1H1, 592
.set nb104nf_rinvH1H2, 608
.set nb104nf_rinvH2M, 624
.set nb104nf_rinvH2H1, 640
.set nb104nf_rinvH2H2, 656
.set nb104nf_is3, 672
.set nb104nf_ii3, 676
.set nb104nf_nri, 692
.set nb104nf_iinr, 700
.set nb104nf_jindex, 708
.set nb104nf_jjnr, 716
.set nb104nf_shift, 724
.set nb104nf_shiftvec, 732
.set nb104nf_facel, 740
.set nb104nf_innerjjnr, 748
.set nb104nf_innerk, 756
.set nb104nf_n, 760
.set nb104nf_nn1, 764
.set nb104nf_nouter, 768
.set nb104nf_ninner, 772
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb104nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb104nf_half(%rsp)
        movl %ebx,nb104nf_half+4(%rsp)
        movsd nb104nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb104nf_half(%rsp)
        movapd %xmm3,nb104nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb104nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb104nf_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3     ## qM 
        movsd %xmm3,%xmm4               ## qM 
        movsd 8(%rdx,%rbx,8),%xmm5      ## qH 
        movq nb104nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb104nf_facel(%rsp),%xmm6         ## facel 
        mulsd  %xmm3,%xmm3              ## qM*qM 
        mulsd  %xmm5,%xmm4              ## qM*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104nf_qqMM(%rsp)
        movapd %xmm4,nb104nf_qqMH(%rsp)
        movapd %xmm5,nb104nf_qqHH(%rsp)

_nb_kernel104nf_x86_64_sse2.nb104nf_threadloop: 
        movq  nb104nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel104nf_x86_64_sse2.nb104nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104nf_x86_64_sse2.nb104nf_spinlock

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
        jg  _nb_kernel104nf_x86_64_sse2.nb104nf_outerstart
        jmp _nb_kernel104nf_x86_64_sse2.nb104nf_end

_nb_kernel104nf_x86_64_sse2.nb104nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104nf_nouter(%rsp),%ebx
        movl %ebx,nb104nf_nouter(%rsp)

_nb_kernel104nf_x86_64_sse2.nb104nf_outer: 
        movq  nb104nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb104nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb104nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb104nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb104nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104nf_ixH1(%rsp)
        movapd %xmm4,nb104nf_iyH1(%rsp)
        movapd %xmm5,nb104nf_izH1(%rsp)

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
        movapd %xmm0,nb104nf_ixH2(%rsp)
        movapd %xmm1,nb104nf_iyH2(%rsp)
        movapd %xmm2,nb104nf_izH2(%rsp)
        movapd %xmm3,nb104nf_ixM(%rsp)
        movapd %xmm4,nb104nf_iyM(%rsp)
        movapd %xmm5,nb104nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb104nf_vctot(%rsp)

        movq  nb104nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb104nf_pos(%rbp),%rsi
        movq  nb104nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb104nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb104nf_ninner(%rsp),%ecx
        movl  %ecx,nb104nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb104nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel104nf_x86_64_sse2.nb104nf_unroll_loop
        jmp   _nb_kernel104nf_x86_64_sse2.nb104nf_checksingle
_nb_kernel104nf_x86_64_sse2.nb104nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb104nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb104nf_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb104nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd  %xmm2,nb104nf_jxH1(%rsp)
        movapd  %xmm3,nb104nf_jyH1(%rsp)
        movapd  %xmm4,nb104nf_jzH1(%rsp)
        movapd  %xmm5,nb104nf_jxH2(%rsp)
        movapd  %xmm6,nb104nf_jyH2(%rsp)
        movapd  %xmm7,nb104nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movhpd 72(%rsi,%rbx,8),%xmm2
        movhpd 80(%rsi,%rbx,8),%xmm3
        movhpd 88(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb104nf_jxM(%rsp)
        movapd  %xmm3,nb104nf_jyM(%rsp)
        movapd  %xmm4,nb104nf_jzM(%rsp)

        movapd nb104nf_ixM(%rsp),%xmm0
        movapd nb104nf_iyM(%rsp),%xmm1
        movapd nb104nf_izM(%rsp),%xmm2
        movapd nb104nf_ixM(%rsp),%xmm3
        movapd nb104nf_iyM(%rsp),%xmm4
        movapd nb104nf_izM(%rsp),%xmm5
        subpd  nb104nf_jxM(%rsp),%xmm0
        subpd  nb104nf_jyM(%rsp),%xmm1
        subpd  nb104nf_jzM(%rsp),%xmm2
        subpd  nb104nf_jxH1(%rsp),%xmm3
        subpd  nb104nf_jyH1(%rsp),%xmm4
        subpd  nb104nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMM(%rsp)
        movapd %xmm3,nb104nf_rsqMH1(%rsp)

        movapd nb104nf_ixM(%rsp),%xmm0
        movapd nb104nf_iyM(%rsp),%xmm1
        movapd nb104nf_izM(%rsp),%xmm2
        movapd nb104nf_ixH1(%rsp),%xmm3
        movapd nb104nf_iyH1(%rsp),%xmm4
        movapd nb104nf_izH1(%rsp),%xmm5
        subpd  nb104nf_jxH2(%rsp),%xmm0
        subpd  nb104nf_jyH2(%rsp),%xmm1
        subpd  nb104nf_jzH2(%rsp),%xmm2
        subpd  nb104nf_jxM(%rsp),%xmm3
        subpd  nb104nf_jyM(%rsp),%xmm4
        subpd  nb104nf_jzM(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMH2(%rsp)
        movapd %xmm3,nb104nf_rsqH1M(%rsp)

        movapd nb104nf_ixH1(%rsp),%xmm0
        movapd nb104nf_iyH1(%rsp),%xmm1
        movapd nb104nf_izH1(%rsp),%xmm2
        movapd nb104nf_ixH1(%rsp),%xmm3
        movapd nb104nf_iyH1(%rsp),%xmm4
        movapd nb104nf_izH1(%rsp),%xmm5
        subpd  nb104nf_jxH1(%rsp),%xmm0
        subpd  nb104nf_jyH1(%rsp),%xmm1
        subpd  nb104nf_jzH1(%rsp),%xmm2
        subpd  nb104nf_jxH2(%rsp),%xmm3
        subpd  nb104nf_jyH2(%rsp),%xmm4
        subpd  nb104nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH1H1(%rsp)
        movapd %xmm3,nb104nf_rsqH1H2(%rsp)

        movapd nb104nf_ixH2(%rsp),%xmm0
        movapd nb104nf_iyH2(%rsp),%xmm1
        movapd nb104nf_izH2(%rsp),%xmm2
        movapd nb104nf_ixH2(%rsp),%xmm3
        movapd nb104nf_iyH2(%rsp),%xmm4
        movapd nb104nf_izH2(%rsp),%xmm5
        subpd  nb104nf_jxM(%rsp),%xmm0
        subpd  nb104nf_jyM(%rsp),%xmm1
        subpd  nb104nf_jzM(%rsp),%xmm2
        subpd  nb104nf_jxH1(%rsp),%xmm3
        subpd  nb104nf_jyH1(%rsp),%xmm4
        subpd  nb104nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH2M(%rsp)
        movapd %xmm4,nb104nf_rsqH2H1(%rsp)

        movapd nb104nf_ixH2(%rsp),%xmm0
        movapd nb104nf_iyH2(%rsp),%xmm1
        movapd nb104nf_izH2(%rsp),%xmm2
        subpd  nb104nf_jxH2(%rsp),%xmm0
        subpd  nb104nf_jyH2(%rsp),%xmm1
        subpd  nb104nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb104nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0 (h2h2) , xmm4 (h2h1) 
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb104nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH2H2(%rsp)
        movapd %xmm5,nb104nf_rinvH2H1(%rsp)

        movapd nb104nf_rsqMM(%rsp),%xmm0
        movapd nb104nf_rsqMH1(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb104nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb104nf_rinvMM(%rsp)
        movapd %xmm5,nb104nf_rinvMH1(%rsp)

        movapd nb104nf_rsqMH2(%rsp),%xmm0
        movapd nb104nf_rsqH1M(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb104nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvMH2(%rsp)
        movapd %xmm5,nb104nf_rinvH1M(%rsp)

        movapd nb104nf_rsqH1H1(%rsp),%xmm0
        movapd nb104nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb104nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH1H1(%rsp)
        movapd %xmm5,nb104nf_rinvH1H2(%rsp)

        movapd nb104nf_rsqH2M(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb104nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb104nf_rinvMM(%rsp),%xmm0
        mulpd  nb104nf_qqMM(%rsp),%xmm0
        addpd  nb104nf_vctot(%rsp),%xmm0

        ## other interactions 
        movapd nb104nf_rinvMH1(%rsp),%xmm1
        movapd nb104nf_rinvH1H1(%rsp),%xmm2

        addpd nb104nf_rinvMH2(%rsp),%xmm1
        addpd nb104nf_rinvH1H2(%rsp),%xmm2

        addpd nb104nf_rinvH1M(%rsp),%xmm1
        addpd nb104nf_rinvH2H1(%rsp),%xmm2

        addpd nb104nf_rinvH2M(%rsp),%xmm1
        addpd nb104nf_rinvH2H2(%rsp),%xmm2

        mulpd nb104nf_qqMH(%rsp),%xmm1
        mulpd nb104nf_qqHH(%rsp),%xmm2

        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0

        movapd %xmm0,nb104nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb104nf_innerk(%rsp)
        jl    _nb_kernel104nf_x86_64_sse2.nb104nf_checksingle
        jmp   _nb_kernel104nf_x86_64_sse2.nb104nf_unroll_loop
_nb_kernel104nf_x86_64_sse2.nb104nf_checksingle: 
        movl  nb104nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel104nf_x86_64_sse2.nb104nf_dosingle
        jmp   _nb_kernel104nf_x86_64_sse2.nb104nf_updateouterdata
_nb_kernel104nf_x86_64_sse2.nb104nf_dosingle: 
        movq  nb104nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb104nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coordinates to local temp variables 
        movlpd 24(%rsi,%rax,8),%xmm2
        movlpd 32(%rsi,%rax,8),%xmm3
        movlpd 40(%rsi,%rax,8),%xmm4
        movlpd 48(%rsi,%rax,8),%xmm5
        movlpd 56(%rsi,%rax,8),%xmm6
        movlpd 64(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb104nf_jxH1(%rsp)
        movapd  %xmm3,nb104nf_jyH1(%rsp)
        movapd  %xmm4,nb104nf_jzH1(%rsp)
        movapd  %xmm5,nb104nf_jxH2(%rsp)
        movapd  %xmm6,nb104nf_jyH2(%rsp)
        movapd  %xmm7,nb104nf_jzH2(%rsp)
        movlpd 72(%rsi,%rax,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 88(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb104nf_jxM(%rsp)
        movapd  %xmm3,nb104nf_jyM(%rsp)
        movapd  %xmm4,nb104nf_jzM(%rsp)

        movapd nb104nf_ixM(%rsp),%xmm0
        movapd nb104nf_iyM(%rsp),%xmm1
        movapd nb104nf_izM(%rsp),%xmm2
        movapd nb104nf_ixM(%rsp),%xmm3
        movapd nb104nf_iyM(%rsp),%xmm4
        movapd nb104nf_izM(%rsp),%xmm5
        subsd  nb104nf_jxM(%rsp),%xmm0
        subsd  nb104nf_jyM(%rsp),%xmm1
        subsd  nb104nf_jzM(%rsp),%xmm2
        subsd  nb104nf_jxH1(%rsp),%xmm3
        subsd  nb104nf_jyH1(%rsp),%xmm4
        subsd  nb104nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMM(%rsp)
        movapd %xmm3,nb104nf_rsqMH1(%rsp)

        movapd nb104nf_ixM(%rsp),%xmm0
        movapd nb104nf_iyM(%rsp),%xmm1
        movapd nb104nf_izM(%rsp),%xmm2
        movapd nb104nf_ixH1(%rsp),%xmm3
        movapd nb104nf_iyH1(%rsp),%xmm4
        movapd nb104nf_izH1(%rsp),%xmm5
        subsd  nb104nf_jxH2(%rsp),%xmm0
        subsd  nb104nf_jyH2(%rsp),%xmm1
        subsd  nb104nf_jzH2(%rsp),%xmm2
        subsd  nb104nf_jxM(%rsp),%xmm3
        subsd  nb104nf_jyM(%rsp),%xmm4
        subsd  nb104nf_jzM(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMH2(%rsp)
        movapd %xmm3,nb104nf_rsqH1M(%rsp)

        movapd nb104nf_ixH1(%rsp),%xmm0
        movapd nb104nf_iyH1(%rsp),%xmm1
        movapd nb104nf_izH1(%rsp),%xmm2
        movapd nb104nf_ixH1(%rsp),%xmm3
        movapd nb104nf_iyH1(%rsp),%xmm4
        movapd nb104nf_izH1(%rsp),%xmm5
        subsd  nb104nf_jxH1(%rsp),%xmm0
        subsd  nb104nf_jyH1(%rsp),%xmm1
        subsd  nb104nf_jzH1(%rsp),%xmm2
        subsd  nb104nf_jxH2(%rsp),%xmm3
        subsd  nb104nf_jyH2(%rsp),%xmm4
        subsd  nb104nf_jzH2(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH1H1(%rsp)
        movapd %xmm3,nb104nf_rsqH1H2(%rsp)

        movapd nb104nf_ixH2(%rsp),%xmm0
        movapd nb104nf_iyH2(%rsp),%xmm1
        movapd nb104nf_izH2(%rsp),%xmm2
        movapd nb104nf_ixH2(%rsp),%xmm3
        movapd nb104nf_iyH2(%rsp),%xmm4
        movapd nb104nf_izH2(%rsp),%xmm5
        subsd  nb104nf_jxM(%rsp),%xmm0
        subsd  nb104nf_jyM(%rsp),%xmm1
        subsd  nb104nf_jzM(%rsp),%xmm2
        subsd  nb104nf_jxH1(%rsp),%xmm3
        subsd  nb104nf_jyH1(%rsp),%xmm4
        subsd  nb104nf_jzH1(%rsp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH2M(%rsp)
        movapd %xmm4,nb104nf_rsqH2H1(%rsp)

        movapd nb104nf_ixH2(%rsp),%xmm0
        movapd nb104nf_iyH2(%rsp),%xmm1
        movapd nb104nf_izH2(%rsp),%xmm2
        subsd  nb104nf_jxH2(%rsp),%xmm0
        subsd  nb104nf_jyH2(%rsp),%xmm1
        subsd  nb104nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb104nf_rsqH2H2(%rsp)

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
        movapd  nb104nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb104nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH2H2(%rsp)
        movapd %xmm5,nb104nf_rinvH2H1(%rsp)

        movapd nb104nf_rsqMM(%rsp),%xmm0
        movapd nb104nf_rsqMH1(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb104nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb104nf_rinvMM(%rsp)
        movapd %xmm5,nb104nf_rinvMH1(%rsp)

        movapd nb104nf_rsqMH2(%rsp),%xmm0
        movapd nb104nf_rsqH1M(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb104nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvMH2(%rsp)
        movapd %xmm5,nb104nf_rinvH1M(%rsp)

        movapd nb104nf_rsqH1H1(%rsp),%xmm0
        movapd nb104nf_rsqH1H2(%rsp),%xmm4
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
        movapd  nb104nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb104nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH1H1(%rsp)
        movapd %xmm5,nb104nf_rinvH1H2(%rsp)

        movapd nb104nf_rsqH2M(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb104nf_rinvH2M(%rsp)

        ## start with MM interaction 
        movapd nb104nf_rinvMM(%rsp),%xmm0
        mulpd  nb104nf_qqMM(%rsp),%xmm0
        addpd  nb104nf_vctot(%rsp),%xmm0

        ## other interactions 
        movapd nb104nf_rinvMH1(%rsp),%xmm1
        movapd nb104nf_rinvH1H1(%rsp),%xmm2

        addsd nb104nf_rinvMH2(%rsp),%xmm1
        addsd nb104nf_rinvH1H2(%rsp),%xmm2

        addsd nb104nf_rinvH1M(%rsp),%xmm1
        addsd nb104nf_rinvH2H1(%rsp),%xmm2

        addsd nb104nf_rinvH2M(%rsp),%xmm1
        addsd nb104nf_rinvH2H2(%rsp),%xmm2

        mulsd nb104nf_qqMH(%rsp),%xmm1
        mulsd nb104nf_qqHH(%rsp),%xmm2

        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0

        movlpd %xmm0,nb104nf_vctot(%rsp)

_nb_kernel104nf_x86_64_sse2.nb104nf_updateouterdata: 
        ## get n from stack
        movl nb104nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb104nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        movapd nb104nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb104nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb104nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104nf_x86_64_sse2.nb104nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104nf_n(%rsp)
        jmp _nb_kernel104nf_x86_64_sse2.nb104nf_outer
_nb_kernel104nf_x86_64_sse2.nb104nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104nf_x86_64_sse2.nb104nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104nf_x86_64_sse2.nb104nf_threadloop
_nb_kernel104nf_x86_64_sse2.nb104nf_end: 
        movl nb104nf_nouter(%rsp),%eax
        movl nb104nf_ninner(%rsp),%ebx
        movq nb104nf_outeriter(%rbp),%rcx
        movq nb104nf_inneriter(%rbp),%rdx
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

