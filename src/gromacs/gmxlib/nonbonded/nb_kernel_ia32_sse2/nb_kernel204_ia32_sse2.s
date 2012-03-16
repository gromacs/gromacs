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



.globl nb_kernel204_ia32_sse2
.globl _nb_kernel204_ia32_sse2
nb_kernel204_ia32_sse2: 
_nb_kernel204_ia32_sse2:        
.set nb204_p_nri, 8
.set nb204_iinr, 12
.set nb204_jindex, 16
.set nb204_jjnr, 20
.set nb204_shift, 24
.set nb204_shiftvec, 28
.set nb204_fshift, 32
.set nb204_gid, 36
.set nb204_pos, 40
.set nb204_faction, 44
.set nb204_charge, 48
.set nb204_p_facel, 52
.set nb204_argkrf, 56
.set nb204_argcrf, 60
.set nb204_Vc, 64
.set nb204_type, 68
.set nb204_p_ntype, 72
.set nb204_vdwparam, 76
.set nb204_Vvdw, 80
.set nb204_p_tabscale, 84
.set nb204_VFtab, 88
.set nb204_invsqrta, 92
.set nb204_dvda, 96
.set nb204_p_gbtabscale, 100
.set nb204_GBtab, 104
.set nb204_p_nthreads, 108
.set nb204_count, 112
.set nb204_mtx, 116
.set nb204_outeriter, 120
.set nb204_inneriter, 124
.set nb204_work, 128
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
.set nb204_innerjjnr, 1448
.set nb204_innerk, 1452
.set nb204_n, 1456
.set nb204_nn1, 1460
.set nb204_nri, 1464
.set nb204_nouter, 1468
.set nb204_ninner, 1472
.set nb204_salign, 1476
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1480,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb204_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb204_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb204_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb204_nouter(%esp)
        movl %eax,nb204_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb204_half(%esp)
        movl %ebx,nb204_half+4(%esp)
        movsd nb204_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb204_half(%esp)
        movapd %xmm2,nb204_two(%esp)
        movapd %xmm3,nb204_three(%esp)

        movl nb204_argkrf(%ebp),%esi
        movl nb204_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb204_krf(%esp)
        movapd %xmm6,nb204_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb204_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb204_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb204_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204_qqMM(%esp)
        movapd %xmm4,nb204_qqMH(%esp)
        movapd %xmm5,nb204_qqHH(%esp)

_nb_kernel204_ia32_sse2.nb204_threadloop: 
        movl  nb204_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel204_ia32_sse2.nb204_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204_ia32_sse2.nb204_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204_n(%esp)
        movl %ebx,nb204_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204_ia32_sse2.nb204_outerstart
        jmp _nb_kernel204_ia32_sse2.nb204_end

_nb_kernel204_ia32_sse2.nb204_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204_nouter(%esp),%ebx
        movl %ebx,nb204_nouter(%esp)

_nb_kernel204_ia32_sse2.nb204_outer: 
        movl  nb204_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb204_is3(%esp)      ## store is3 

        movl  nb204_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb204_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb204_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb204_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204_ixH1(%esp)
        movapd %xmm4,nb204_iyH1(%esp)
        movapd %xmm5,nb204_izH1(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%eax,%ebx,8),%xmm0
        addsd 56(%eax,%ebx,8),%xmm1
        addsd 64(%eax,%ebx,8),%xmm2
        addsd 72(%eax,%ebx,8),%xmm3
        addsd 80(%eax,%ebx,8),%xmm4
        addsd 88(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb204_ixH2(%esp)
        movapd %xmm1,nb204_iyH2(%esp)
        movapd %xmm2,nb204_izH2(%esp)
        movapd %xmm3,nb204_ixM(%esp)
        movapd %xmm4,nb204_iyM(%esp)
        movapd %xmm5,nb204_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb204_vctot(%esp)
        movapd %xmm4,nb204_fixM(%esp)
        movapd %xmm4,nb204_fiyM(%esp)
        movapd %xmm4,nb204_fizM(%esp)
        movapd %xmm4,nb204_fixH1(%esp)
        movapd %xmm4,nb204_fiyH1(%esp)
        movapd %xmm4,nb204_fizH1(%esp)
        movapd %xmm4,nb204_fixH2(%esp)
        movapd %xmm4,nb204_fiyH2(%esp)
        movapd %xmm4,nb204_fizH2(%esp)

        movl  nb204_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb204_pos(%ebp),%esi
        movl  nb204_faction(%ebp),%edi
        movl  nb204_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb204_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb204_ninner(%esp),%ecx
        movl  %ecx,nb204_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb204_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel204_ia32_sse2.nb204_unroll_loop
        jmp   _nb_kernel204_ia32_sse2.nb204_checksingle
_nb_kernel204_ia32_sse2.nb204_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb204_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb204_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb204_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movhpd 24(%esi,%ebx,8),%xmm2
        movhpd 32(%esi,%ebx,8),%xmm3
        movhpd 40(%esi,%ebx,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%ebx,8),%xmm6
        movhpd 64(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb204_jxH1(%esp)
        movapd  %xmm3,nb204_jyH1(%esp)
        movapd  %xmm4,nb204_jzH1(%esp)
        movapd  %xmm5,nb204_jxH2(%esp)
        movapd  %xmm6,nb204_jyH2(%esp)
        movapd  %xmm7,nb204_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb204_jxM(%esp)
        movapd  %xmm3,nb204_jyM(%esp)
        movapd  %xmm4,nb204_jzM(%esp)

        movapd nb204_ixM(%esp),%xmm0
        movapd nb204_iyM(%esp),%xmm1
        movapd nb204_izM(%esp),%xmm2
        movapd nb204_ixM(%esp),%xmm3
        movapd nb204_iyM(%esp),%xmm4
        movapd nb204_izM(%esp),%xmm5
        subpd  nb204_jxM(%esp),%xmm0
        subpd  nb204_jyM(%esp),%xmm1
        subpd  nb204_jzM(%esp),%xmm2
        subpd  nb204_jxH1(%esp),%xmm3
        subpd  nb204_jyH1(%esp),%xmm4
        subpd  nb204_jzH1(%esp),%xmm5
        movapd %xmm0,nb204_dxMM(%esp)
        movapd %xmm1,nb204_dyMM(%esp)
        movapd %xmm2,nb204_dzMM(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxMH1(%esp)
        movapd %xmm4,nb204_dyMH1(%esp)
        movapd %xmm5,nb204_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqMM(%esp)
        movapd %xmm3,nb204_rsqMH1(%esp)

        movapd nb204_ixM(%esp),%xmm0
        movapd nb204_iyM(%esp),%xmm1
        movapd nb204_izM(%esp),%xmm2
        movapd nb204_ixH1(%esp),%xmm3
        movapd nb204_iyH1(%esp),%xmm4
        movapd nb204_izH1(%esp),%xmm5
        subpd  nb204_jxH2(%esp),%xmm0
        subpd  nb204_jyH2(%esp),%xmm1
        subpd  nb204_jzH2(%esp),%xmm2
        subpd  nb204_jxM(%esp),%xmm3
        subpd  nb204_jyM(%esp),%xmm4
        subpd  nb204_jzM(%esp),%xmm5
        movapd %xmm0,nb204_dxMH2(%esp)
        movapd %xmm1,nb204_dyMH2(%esp)
        movapd %xmm2,nb204_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH1M(%esp)
        movapd %xmm4,nb204_dyH1M(%esp)
        movapd %xmm5,nb204_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqMH2(%esp)
        movapd %xmm3,nb204_rsqH1M(%esp)

        movapd nb204_ixH1(%esp),%xmm0
        movapd nb204_iyH1(%esp),%xmm1
        movapd nb204_izH1(%esp),%xmm2
        movapd nb204_ixH1(%esp),%xmm3
        movapd nb204_iyH1(%esp),%xmm4
        movapd nb204_izH1(%esp),%xmm5
        subpd  nb204_jxH1(%esp),%xmm0
        subpd  nb204_jyH1(%esp),%xmm1
        subpd  nb204_jzH1(%esp),%xmm2
        subpd  nb204_jxH2(%esp),%xmm3
        subpd  nb204_jyH2(%esp),%xmm4
        subpd  nb204_jzH2(%esp),%xmm5
        movapd %xmm0,nb204_dxH1H1(%esp)
        movapd %xmm1,nb204_dyH1H1(%esp)
        movapd %xmm2,nb204_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH1H2(%esp)
        movapd %xmm4,nb204_dyH1H2(%esp)
        movapd %xmm5,nb204_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqH1H1(%esp)
        movapd %xmm3,nb204_rsqH1H2(%esp)

        movapd nb204_ixH2(%esp),%xmm0
        movapd nb204_iyH2(%esp),%xmm1
        movapd nb204_izH2(%esp),%xmm2
        movapd nb204_ixH2(%esp),%xmm3
        movapd nb204_iyH2(%esp),%xmm4
        movapd nb204_izH2(%esp),%xmm5
        subpd  nb204_jxM(%esp),%xmm0
        subpd  nb204_jyM(%esp),%xmm1
        subpd  nb204_jzM(%esp),%xmm2
        subpd  nb204_jxH1(%esp),%xmm3
        subpd  nb204_jyH1(%esp),%xmm4
        subpd  nb204_jzH1(%esp),%xmm5
        movapd %xmm0,nb204_dxH2M(%esp)
        movapd %xmm1,nb204_dyH2M(%esp)
        movapd %xmm2,nb204_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH2H1(%esp)
        movapd %xmm4,nb204_dyH2H1(%esp)
        movapd %xmm5,nb204_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb204_rsqH2M(%esp)
        movapd %xmm4,nb204_rsqH2H1(%esp)

        movapd nb204_ixH2(%esp),%xmm0
        movapd nb204_iyH2(%esp),%xmm1
        movapd nb204_izH2(%esp),%xmm2
        subpd  nb204_jxH2(%esp),%xmm0
        subpd  nb204_jyH2(%esp),%xmm1
        subpd  nb204_jzH2(%esp),%xmm2
        movapd %xmm0,nb204_dxH2H2(%esp)
        movapd %xmm1,nb204_dyH2H2(%esp)
        movapd %xmm2,nb204_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb204_rsqH2H2(%esp)

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
        movapd  nb204_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm3   ## iter1 
        mulpd   nb204_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm1   ## rinv 
        mulpd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvH2H2(%esp)
        movapd %xmm5,nb204_rinvH2H1(%esp)

        movapd nb204_rsqMM(%esp),%xmm0
        movapd nb204_rsqMH1(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb204_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm1   ## rinv 
        mulpd   nb204_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb204_rinvMM(%esp)
        movapd %xmm5,nb204_rinvMH1(%esp)

        movapd nb204_rsqMH2(%esp),%xmm0
        movapd nb204_rsqH1M(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm3   ## iter1 
        mulpd   nb204_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm1   ## rinv 
        mulpd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvMH2(%esp)
        movapd %xmm5,nb204_rinvH1M(%esp)

        movapd nb204_rsqH1H1(%esp),%xmm0
        movapd nb204_rsqH1H2(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm3   ## iter1a 
        mulpd   nb204_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204_half(%esp),%xmm1   ## rinv 
        mulpd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvH1H1(%esp)
        movapd %xmm5,nb204_rinvH1H2(%esp)

        movapd nb204_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb204_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb204_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## rinvsq 
        mulpd  nb204_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb204_crf(%esp),%xmm6

        mulpd  nb204_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb204_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb204_dxMM(%esp),%xmm0
        mulpd nb204_dyMM(%esp),%xmm1
        mulpd nb204_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixM(%esp),%xmm0
        addpd nb204_fiyM(%esp),%xmm1
        addpd nb204_fizM(%esp),%xmm2
        movapd %xmm3,nb204_fjxM(%esp)
        movapd %xmm4,nb204_fjyM(%esp)
        movapd %xmm5,nb204_fjzM(%esp)
        movapd %xmm0,nb204_fixM(%esp)
        movapd %xmm1,nb204_fiyM(%esp)
        movapd %xmm2,nb204_fizM(%esp)

        ## M-H1 interaction 
        movapd nb204_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd  %xmm0,%xmm0
        subpd  nb204_crf(%esp),%xmm4
        mulpd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb204_dxMH1(%esp),%xmm0
        mulpd nb204_dyMH1(%esp),%xmm1
        mulpd nb204_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixM(%esp),%xmm0
        addpd nb204_fiyM(%esp),%xmm1
        addpd nb204_fizM(%esp),%xmm2
        movapd %xmm3,nb204_fjxH1(%esp)
        movapd %xmm4,nb204_fjyH1(%esp)
        movapd %xmm5,nb204_fjzH1(%esp)
        movapd %xmm0,nb204_fixM(%esp)
        movapd %xmm1,nb204_fiyM(%esp)
        movapd %xmm2,nb204_fizM(%esp)

        ## M-H2 interaction  
        movapd nb204_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb204_crf(%esp),%xmm4
        mulpd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb204_dxMH2(%esp),%xmm0
        mulpd nb204_dyMH2(%esp),%xmm1
        mulpd nb204_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixM(%esp),%xmm0
        addpd nb204_fiyM(%esp),%xmm1
        addpd nb204_fizM(%esp),%xmm2
        movapd %xmm3,nb204_fjxH2(%esp)
        movapd %xmm4,nb204_fjyH2(%esp)
        movapd %xmm5,nb204_fjzH2(%esp)
        movapd %xmm0,nb204_fixM(%esp)
        movapd %xmm1,nb204_fiyM(%esp)
        movapd %xmm2,nb204_fizM(%esp)

        ## H1-M interaction 
        movapd nb204_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb204_crf(%esp),%xmm4
        mulpd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxM(%esp),%xmm3
        movapd nb204_fjyM(%esp),%xmm4
        movapd nb204_fjzM(%esp),%xmm5
        mulpd nb204_dxH1M(%esp),%xmm0
        mulpd nb204_dyH1M(%esp),%xmm1
        mulpd nb204_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH1(%esp),%xmm0
        addpd nb204_fiyH1(%esp),%xmm1
        addpd nb204_fizH1(%esp),%xmm2
        movapd %xmm3,nb204_fjxM(%esp)
        movapd %xmm4,nb204_fjyM(%esp)
        movapd %xmm5,nb204_fjzM(%esp)
        movapd %xmm0,nb204_fixH1(%esp)
        movapd %xmm1,nb204_fiyH1(%esp)
        movapd %xmm2,nb204_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb204_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb204_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH1(%esp),%xmm3
        movapd nb204_fjyH1(%esp),%xmm4
        movapd nb204_fjzH1(%esp),%xmm5
        mulpd nb204_dxH1H1(%esp),%xmm0
        mulpd nb204_dyH1H1(%esp),%xmm1
        mulpd nb204_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH1(%esp),%xmm0
        addpd nb204_fiyH1(%esp),%xmm1
        addpd nb204_fizH1(%esp),%xmm2
        movapd %xmm3,nb204_fjxH1(%esp)
        movapd %xmm4,nb204_fjyH1(%esp)
        movapd %xmm5,nb204_fjzH1(%esp)
        movapd %xmm0,nb204_fixH1(%esp)
        movapd %xmm1,nb204_fiyH1(%esp)
        movapd %xmm2,nb204_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb204_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb204_crf(%esp),%xmm4
        mulpd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH2(%esp),%xmm3
        movapd nb204_fjyH2(%esp),%xmm4
        movapd nb204_fjzH2(%esp),%xmm5
        mulpd nb204_dxH1H2(%esp),%xmm0
        mulpd nb204_dyH1H2(%esp),%xmm1
        mulpd nb204_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH1(%esp),%xmm0
        addpd nb204_fiyH1(%esp),%xmm1
        addpd nb204_fizH1(%esp),%xmm2
        movapd %xmm3,nb204_fjxH2(%esp)
        movapd %xmm4,nb204_fjyH2(%esp)
        movapd %xmm5,nb204_fjzH2(%esp)
        movapd %xmm0,nb204_fixH1(%esp)
        movapd %xmm1,nb204_fiyH1(%esp)
        movapd %xmm2,nb204_fizH1(%esp)

        ## H2-M interaction 
        movapd nb204_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb204_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxM(%esp),%xmm3
        movapd nb204_fjyM(%esp),%xmm4
        movapd nb204_fjzM(%esp),%xmm5
        mulpd nb204_dxH2M(%esp),%xmm0
        mulpd nb204_dyH2M(%esp),%xmm1
        mulpd nb204_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH2(%esp),%xmm0
        addpd nb204_fiyH2(%esp),%xmm1
        addpd nb204_fizH2(%esp),%xmm2
        movapd %xmm3,nb204_fjxM(%esp)
        movapd %xmm4,nb204_fjyM(%esp)
        movapd %xmm5,nb204_fjzM(%esp)
        movapd %xmm0,nb204_fixH2(%esp)
        movapd %xmm1,nb204_fiyH2(%esp)
        movapd %xmm2,nb204_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb204_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb204_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH1(%esp),%xmm3
        movapd nb204_fjyH1(%esp),%xmm4
        movapd nb204_fjzH1(%esp),%xmm5
        mulpd nb204_dxH2H1(%esp),%xmm0
        mulpd nb204_dyH2H1(%esp),%xmm1
        mulpd nb204_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH2(%esp),%xmm0
        addpd nb204_fiyH2(%esp),%xmm1
        addpd nb204_fizH2(%esp),%xmm2
        movapd %xmm3,nb204_fjxH1(%esp)
        movapd %xmm4,nb204_fjyH1(%esp)
        movapd %xmm5,nb204_fjzH1(%esp)
        movapd %xmm0,nb204_fixH2(%esp)
        movapd %xmm1,nb204_fiyH2(%esp)
        movapd %xmm2,nb204_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb204_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb204_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb204_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb204_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movapd %xmm6,nb204_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb204_fjxH2(%esp),%xmm3
        movapd nb204_fjyH2(%esp),%xmm4
        movapd nb204_fjzH2(%esp),%xmm5
        mulpd nb204_dxH2H2(%esp),%xmm0
        mulpd nb204_dyH2H2(%esp),%xmm1
        mulpd nb204_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb204_fixH2(%esp),%xmm0
        addpd nb204_fiyH2(%esp),%xmm1
        addpd nb204_fizH2(%esp),%xmm2
        movapd %xmm3,nb204_fjxH2(%esp)
        movapd %xmm4,nb204_fjyH2(%esp)
        movapd %xmm5,nb204_fjzH2(%esp)
        movapd %xmm0,nb204_fixH2(%esp)
        movapd %xmm1,nb204_fiyH2(%esp)
        movapd %xmm2,nb204_fizH2(%esp)

        movl nb204_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        movlpd 24(%edi,%eax,8),%xmm0
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 40(%edi,%eax,8),%xmm2
        movlpd 48(%edi,%eax,8),%xmm3
        movlpd 56(%edi,%eax,8),%xmm4
        movlpd 64(%edi,%eax,8),%xmm5
        movlpd 72(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%eax,8),%xmm7
        movhpd 24(%edi,%ebx,8),%xmm0
        movhpd 32(%edi,%ebx,8),%xmm1
        movhpd 40(%edi,%ebx,8),%xmm2
        movhpd 48(%edi,%ebx,8),%xmm3
        movhpd 56(%edi,%ebx,8),%xmm4
        movhpd 64(%edi,%ebx,8),%xmm5
        movhpd 72(%edi,%ebx,8),%xmm6
        movhpd 80(%edi,%ebx,8),%xmm7
        addpd nb204_fjxH1(%esp),%xmm0
        addpd nb204_fjyH1(%esp),%xmm1
        addpd nb204_fjzH1(%esp),%xmm2
        addpd nb204_fjxH2(%esp),%xmm3
        addpd nb204_fjyH2(%esp),%xmm4
        addpd nb204_fjzH2(%esp),%xmm5
        addpd nb204_fjxM(%esp),%xmm6
        addpd nb204_fjyM(%esp),%xmm7
        movlpd %xmm0,24(%edi,%eax,8)
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movlpd %xmm4,56(%edi,%eax,8)
        movlpd %xmm5,64(%edi,%eax,8)
        movlpd %xmm6,72(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%eax,8)
        movhpd %xmm0,24(%edi,%ebx,8)
        movhpd %xmm1,32(%edi,%ebx,8)
        movhpd %xmm2,40(%edi,%ebx,8)
        movhpd %xmm3,48(%edi,%ebx,8)
        movhpd %xmm4,56(%edi,%ebx,8)
        movhpd %xmm5,64(%edi,%ebx,8)
        movhpd %xmm6,72(%edi,%ebx,8)
        movhpd %xmm7,80(%edi,%ebx,8)

        movlpd 88(%edi,%eax,8),%xmm0
        movhpd 88(%edi,%ebx,8),%xmm0
        addpd nb204_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)
        movhpd %xmm0,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb204_innerk(%esp)
        jl    _nb_kernel204_ia32_sse2.nb204_checksingle
        jmp   _nb_kernel204_ia32_sse2.nb204_unroll_loop
_nb_kernel204_ia32_sse2.nb204_checksingle: 
        movl  nb204_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel204_ia32_sse2.nb204_dosingle
        jmp   _nb_kernel204_ia32_sse2.nb204_updateouterdata
_nb_kernel204_ia32_sse2.nb204_dosingle: 
        movl  nb204_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb204_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb204_jxH1(%esp)
        movapd  %xmm3,nb204_jyH1(%esp)
        movapd  %xmm4,nb204_jzH1(%esp)
        movapd  %xmm5,nb204_jxH2(%esp)
        movapd  %xmm6,nb204_jyH2(%esp)
        movapd  %xmm7,nb204_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb204_jxM(%esp)
        movapd  %xmm3,nb204_jyM(%esp)
        movapd  %xmm4,nb204_jzM(%esp)

        movapd nb204_ixM(%esp),%xmm0
        movapd nb204_iyM(%esp),%xmm1
        movapd nb204_izM(%esp),%xmm2
        movapd nb204_ixM(%esp),%xmm3
        movapd nb204_iyM(%esp),%xmm4
        movapd nb204_izM(%esp),%xmm5
        subsd  nb204_jxM(%esp),%xmm0
        subsd  nb204_jyM(%esp),%xmm1
        subsd  nb204_jzM(%esp),%xmm2
        subsd  nb204_jxH1(%esp),%xmm3
        subsd  nb204_jyH1(%esp),%xmm4
        subsd  nb204_jzH1(%esp),%xmm5
        movapd %xmm0,nb204_dxMM(%esp)
        movapd %xmm1,nb204_dyMM(%esp)
        movapd %xmm2,nb204_dzMM(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxMH1(%esp)
        movapd %xmm4,nb204_dyMH1(%esp)
        movapd %xmm5,nb204_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqMM(%esp)
        movapd %xmm3,nb204_rsqMH1(%esp)

        movapd nb204_ixM(%esp),%xmm0
        movapd nb204_iyM(%esp),%xmm1
        movapd nb204_izM(%esp),%xmm2
        movapd nb204_ixH1(%esp),%xmm3
        movapd nb204_iyH1(%esp),%xmm4
        movapd nb204_izH1(%esp),%xmm5
        subsd  nb204_jxH2(%esp),%xmm0
        subsd  nb204_jyH2(%esp),%xmm1
        subsd  nb204_jzH2(%esp),%xmm2
        subsd  nb204_jxM(%esp),%xmm3
        subsd  nb204_jyM(%esp),%xmm4
        subsd  nb204_jzM(%esp),%xmm5
        movapd %xmm0,nb204_dxMH2(%esp)
        movapd %xmm1,nb204_dyMH2(%esp)
        movapd %xmm2,nb204_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH1M(%esp)
        movapd %xmm4,nb204_dyH1M(%esp)
        movapd %xmm5,nb204_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqMH2(%esp)
        movapd %xmm3,nb204_rsqH1M(%esp)

        movapd nb204_ixH1(%esp),%xmm0
        movapd nb204_iyH1(%esp),%xmm1
        movapd nb204_izH1(%esp),%xmm2
        movapd nb204_ixH1(%esp),%xmm3
        movapd nb204_iyH1(%esp),%xmm4
        movapd nb204_izH1(%esp),%xmm5
        subsd  nb204_jxH1(%esp),%xmm0
        subsd  nb204_jyH1(%esp),%xmm1
        subsd  nb204_jzH1(%esp),%xmm2
        subsd  nb204_jxH2(%esp),%xmm3
        subsd  nb204_jyH2(%esp),%xmm4
        subsd  nb204_jzH2(%esp),%xmm5
        movapd %xmm0,nb204_dxH1H1(%esp)
        movapd %xmm1,nb204_dyH1H1(%esp)
        movapd %xmm2,nb204_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH1H2(%esp)
        movapd %xmm4,nb204_dyH1H2(%esp)
        movapd %xmm5,nb204_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb204_rsqH1H1(%esp)
        movapd %xmm3,nb204_rsqH1H2(%esp)

        movapd nb204_ixH2(%esp),%xmm0
        movapd nb204_iyH2(%esp),%xmm1
        movapd nb204_izH2(%esp),%xmm2
        movapd nb204_ixH2(%esp),%xmm3
        movapd nb204_iyH2(%esp),%xmm4
        movapd nb204_izH2(%esp),%xmm5
        subsd  nb204_jxM(%esp),%xmm0
        subsd  nb204_jyM(%esp),%xmm1
        subsd  nb204_jzM(%esp),%xmm2
        subsd  nb204_jxH1(%esp),%xmm3
        subsd  nb204_jyH1(%esp),%xmm4
        subsd  nb204_jzH1(%esp),%xmm5
        movapd %xmm0,nb204_dxH2M(%esp)
        movapd %xmm1,nb204_dyH2M(%esp)
        movapd %xmm2,nb204_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb204_dxH2H1(%esp)
        movapd %xmm4,nb204_dyH2H1(%esp)
        movapd %xmm5,nb204_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb204_rsqH2M(%esp)
        movapd %xmm4,nb204_rsqH2H1(%esp)

        movapd nb204_ixH2(%esp),%xmm0
        movapd nb204_iyH2(%esp),%xmm1
        movapd nb204_izH2(%esp),%xmm2
        subsd  nb204_jxH2(%esp),%xmm0
        subsd  nb204_jyH2(%esp),%xmm1
        subsd  nb204_jzH2(%esp),%xmm2
        movapd %xmm0,nb204_dxH2H2(%esp)
        movapd %xmm1,nb204_dyH2H2(%esp)
        movapd %xmm2,nb204_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb204_rsqH2H2(%esp)

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
        movapd  nb204_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm3   ## iter1 
        mulsd   nb204_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm1   ## rinv 
        mulsd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvH2H2(%esp)
        movapd %xmm5,nb204_rinvH2H1(%esp)

        movapd nb204_rsqMM(%esp),%xmm0
        movapd nb204_rsqMH1(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb204_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm1   ## rinv 
        mulsd   nb204_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb204_rinvMM(%esp)
        movapd %xmm5,nb204_rinvMH1(%esp)

        movapd nb204_rsqMH2(%esp),%xmm0
        movapd nb204_rsqH1M(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm3   ## iter1 
        mulsd   nb204_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm1   ## rinv 
        mulsd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvMH2(%esp)
        movapd %xmm5,nb204_rinvH1M(%esp)

        movapd nb204_rsqH1H1(%esp),%xmm0
        movapd nb204_rsqH1H2(%esp),%xmm4
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
        movapd  nb204_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm3   ## iter1a 
        mulsd   nb204_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204_half(%esp),%xmm1   ## rinv 
        mulsd   nb204_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204_rinvH1H1(%esp)
        movapd %xmm5,nb204_rinvH1H2(%esp)

        movapd nb204_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb204_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb204_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb204_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb204_crf(%esp),%xmm6

        mulsd  nb204_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb204_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb204_dxMM(%esp),%xmm0
        mulsd nb204_dyMM(%esp),%xmm1
        mulsd nb204_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixM(%esp),%xmm0
        addsd nb204_fiyM(%esp),%xmm1
        addsd nb204_fizM(%esp),%xmm2
        movlpd %xmm3,nb204_fjxM(%esp)
        movlpd %xmm4,nb204_fjyM(%esp)
        movlpd %xmm5,nb204_fjzM(%esp)
        movlpd %xmm0,nb204_fixM(%esp)
        movlpd %xmm1,nb204_fiyM(%esp)
        movlpd %xmm2,nb204_fizM(%esp)

        ## M-H1 interaction 
        movapd nb204_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb204_crf(%esp),%xmm4
        mulsd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb204_dxMH1(%esp),%xmm0
        mulsd nb204_dyMH1(%esp),%xmm1
        mulsd nb204_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixM(%esp),%xmm0
        addsd nb204_fiyM(%esp),%xmm1
        addsd nb204_fizM(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH1(%esp)
        movlpd %xmm4,nb204_fjyH1(%esp)
        movlpd %xmm5,nb204_fjzH1(%esp)
        movlpd %xmm0,nb204_fixM(%esp)
        movlpd %xmm1,nb204_fiyM(%esp)
        movlpd %xmm2,nb204_fizM(%esp)

        ## M-H2 interaction  
        movapd nb204_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb204_crf(%esp),%xmm4
        mulsd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb204_dxMH2(%esp),%xmm0
        mulsd nb204_dyMH2(%esp),%xmm1
        mulsd nb204_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixM(%esp),%xmm0
        addsd nb204_fiyM(%esp),%xmm1
        addsd nb204_fizM(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH2(%esp)
        movlpd %xmm4,nb204_fjyH2(%esp)
        movlpd %xmm5,nb204_fjzH2(%esp)
        movlpd %xmm0,nb204_fixM(%esp)
        movlpd %xmm1,nb204_fiyM(%esp)
        movlpd %xmm2,nb204_fizM(%esp)

        ## H1-M interaction 
        movapd nb204_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb204_crf(%esp),%xmm4
        mulsd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxM(%esp),%xmm3
        movapd nb204_fjyM(%esp),%xmm4
        movapd nb204_fjzM(%esp),%xmm5
        mulsd nb204_dxH1M(%esp),%xmm0
        mulsd nb204_dyH1M(%esp),%xmm1
        mulsd nb204_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH1(%esp),%xmm0
        addsd nb204_fiyH1(%esp),%xmm1
        addsd nb204_fizH1(%esp),%xmm2
        movlpd %xmm3,nb204_fjxM(%esp)
        movlpd %xmm4,nb204_fjyM(%esp)
        movlpd %xmm5,nb204_fjzM(%esp)
        movlpd %xmm0,nb204_fixH1(%esp)
        movlpd %xmm1,nb204_fiyH1(%esp)
        movlpd %xmm2,nb204_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb204_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb204_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH1(%esp),%xmm3
        movapd nb204_fjyH1(%esp),%xmm4
        movapd nb204_fjzH1(%esp),%xmm5
        mulsd nb204_dxH1H1(%esp),%xmm0
        mulsd nb204_dyH1H1(%esp),%xmm1
        mulsd nb204_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH1(%esp),%xmm0
        addsd nb204_fiyH1(%esp),%xmm1
        addsd nb204_fizH1(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH1(%esp)
        movlpd %xmm4,nb204_fjyH1(%esp)
        movlpd %xmm5,nb204_fjzH1(%esp)
        movlpd %xmm0,nb204_fixH1(%esp)
        movlpd %xmm1,nb204_fiyH1(%esp)
        movlpd %xmm2,nb204_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb204_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb204_crf(%esp),%xmm4
        mulsd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH2(%esp),%xmm3
        movapd nb204_fjyH2(%esp),%xmm4
        movapd nb204_fjzH2(%esp),%xmm5
        mulsd nb204_dxH1H2(%esp),%xmm0
        mulsd nb204_dyH1H2(%esp),%xmm1
        mulsd nb204_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH1(%esp),%xmm0
        addsd nb204_fiyH1(%esp),%xmm1
        addsd nb204_fizH1(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH2(%esp)
        movlpd %xmm4,nb204_fjyH2(%esp)
        movlpd %xmm5,nb204_fjzH2(%esp)
        movlpd %xmm0,nb204_fixH1(%esp)
        movlpd %xmm1,nb204_fiyH1(%esp)
        movlpd %xmm2,nb204_fizH1(%esp)

        ## H2-M interaction 
        movapd nb204_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb204_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxM(%esp),%xmm3
        movapd nb204_fjyM(%esp),%xmm4
        movapd nb204_fjzM(%esp),%xmm5
        mulsd nb204_dxH2M(%esp),%xmm0
        mulsd nb204_dyH2M(%esp),%xmm1
        mulsd nb204_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH2(%esp),%xmm0
        addsd nb204_fiyH2(%esp),%xmm1
        addsd nb204_fizH2(%esp),%xmm2
        movlpd %xmm3,nb204_fjxM(%esp)
        movlpd %xmm4,nb204_fjyM(%esp)
        movlpd %xmm5,nb204_fjzM(%esp)
        movlpd %xmm0,nb204_fixH2(%esp)
        movlpd %xmm1,nb204_fiyH2(%esp)
        movlpd %xmm2,nb204_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb204_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb204_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb204_fjxH1(%esp),%xmm3
        movapd nb204_fjyH1(%esp),%xmm4
        movapd nb204_fjzH1(%esp),%xmm5
        mulsd nb204_dxH2H1(%esp),%xmm0
        mulsd nb204_dyH2H1(%esp),%xmm1
        mulsd nb204_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH2(%esp),%xmm0
        addsd nb204_fiyH2(%esp),%xmm1
        addsd nb204_fizH2(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH1(%esp)
        movlpd %xmm4,nb204_fjyH1(%esp)
        movlpd %xmm5,nb204_fjzH1(%esp)
        movlpd %xmm0,nb204_fixH2(%esp)
        movlpd %xmm1,nb204_fiyH2(%esp)
        movlpd %xmm2,nb204_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb204_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb204_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb204_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb204_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb204_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsMH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movlpd %xmm6,nb204_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb204_fjxH2(%esp),%xmm3
        movapd nb204_fjyH2(%esp),%xmm4
        movapd nb204_fjzH2(%esp),%xmm5
        mulsd nb204_dxH2H2(%esp),%xmm0
        mulsd nb204_dyH2H2(%esp),%xmm1
        mulsd nb204_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb204_fixH2(%esp),%xmm0
        addsd nb204_fiyH2(%esp),%xmm1
        addsd nb204_fizH2(%esp),%xmm2
        movlpd %xmm3,nb204_fjxH2(%esp)
        movlpd %xmm4,nb204_fjyH2(%esp)
        movlpd %xmm5,nb204_fjzH2(%esp)
        movlpd %xmm0,nb204_fixH2(%esp)
        movlpd %xmm1,nb204_fiyH2(%esp)
        movlpd %xmm2,nb204_fizH2(%esp)

        movl nb204_faction(%ebp),%edi
        ## Did all interactions - now update j forces 
        movlpd 24(%edi,%eax,8),%xmm0
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 40(%edi,%eax,8),%xmm2
        movlpd 48(%edi,%eax,8),%xmm3
        movlpd 56(%edi,%eax,8),%xmm4
        movlpd 64(%edi,%eax,8),%xmm5
        movlpd 72(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%eax,8),%xmm7
        addsd nb204_fjxH1(%esp),%xmm0
        addsd nb204_fjyH1(%esp),%xmm1
        addsd nb204_fjzH1(%esp),%xmm2
        addsd nb204_fjxH2(%esp),%xmm3
        addsd nb204_fjyH2(%esp),%xmm4
        addsd nb204_fjzH2(%esp),%xmm5
        addsd nb204_fjxM(%esp),%xmm6
        addsd nb204_fjyM(%esp),%xmm7
        movlpd %xmm0,24(%edi,%eax,8)
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movlpd %xmm4,56(%edi,%eax,8)
        movlpd %xmm5,64(%edi,%eax,8)
        movlpd %xmm6,72(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%eax,8)

        movlpd 88(%edi,%eax,8),%xmm0
        addsd nb204_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)

_nb_kernel204_ia32_sse2.nb204_updateouterdata: 
        movl  nb204_ii3(%esp),%ecx
        movl  nb204_faction(%ebp),%edi
        movl  nb204_fshift(%ebp),%esi
        movl  nb204_is3(%esp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixH1(%esp),%xmm0
        movapd nb204_fiyH1(%esp),%xmm1
        movapd nb204_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixH2(%esp),%xmm0
        movapd nb204_fiyH2(%esp),%xmm1
        movapd nb204_fizH2(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb204_fixM(%esp),%xmm0
        movapd nb204_fiyM(%esp),%xmm1
        movapd nb204_fizM(%esp),%xmm2

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
        movsd  72(%edi,%ecx,8),%xmm3
        movsd  80(%edi,%ecx,8),%xmm4
        movsd  88(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,72(%edi,%ecx,8)
        movsd  %xmm4,80(%edi,%ecx,8)
        movsd  %xmm5,88(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb204_n(%esp),%esi
        ## get group index for i particle 
        movl  nb204_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb204_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb204_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb204_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204_ia32_sse2.nb204_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204_n(%esp)
        jmp _nb_kernel204_ia32_sse2.nb204_outer
_nb_kernel204_ia32_sse2.nb204_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204_ia32_sse2.nb204_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204_ia32_sse2.nb204_threadloop
_nb_kernel204_ia32_sse2.nb204_end: 
        emms

        movl nb204_nouter(%esp),%eax
        movl nb204_ninner(%esp),%ebx
        movl nb204_outeriter(%ebp),%ecx
        movl nb204_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb204_salign(%esp),%eax
        addl %eax,%esp
        addl $1480,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret




.globl nb_kernel204nf_ia32_sse2
.globl _nb_kernel204nf_ia32_sse2
nb_kernel204nf_ia32_sse2:       
_nb_kernel204nf_ia32_sse2:      
.set nb204nf_p_nri, 8
.set nb204nf_iinr, 12
.set nb204nf_jindex, 16
.set nb204nf_jjnr, 20
.set nb204nf_shift, 24
.set nb204nf_shiftvec, 28
.set nb204nf_fshift, 32
.set nb204nf_gid, 36
.set nb204nf_pos, 40
.set nb204nf_faction, 44
.set nb204nf_charge, 48
.set nb204nf_p_facel, 52
.set nb204nf_argkrf, 56
.set nb204nf_argcrf, 60
.set nb204nf_Vc, 64
.set nb204nf_type, 68
.set nb204nf_p_ntype, 72
.set nb204nf_vdwparam, 76
.set nb204nf_Vvdw, 80
.set nb204nf_p_tabscale, 84
.set nb204nf_VFtab, 88
.set nb204nf_invsqrta, 92
.set nb204nf_dvda, 96
.set nb204nf_p_gbtabscale, 100
.set nb204nf_GBtab, 104
.set nb204nf_p_nthreads, 108
.set nb204nf_count, 112
.set nb204nf_mtx, 116
.set nb204nf_outeriter, 120
.set nb204nf_inneriter, 124
.set nb204nf_work, 128
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
.set nb204nf_innerjjnr, 712
.set nb204nf_innerk, 716
.set nb204nf_n, 720
.set nb204nf_nn1, 724
.set nb204nf_nri, 728
.set nb204nf_nouter, 732
.set nb204nf_ninner, 736
.set nb204nf_salign, 740
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $744,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb204nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb204nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb204nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb204nf_nouter(%esp)
        movl %eax,nb204nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb204nf_half(%esp)
        movl %ebx,nb204nf_half+4(%esp)
        movsd nb204nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb204nf_half(%esp)
        movapd %xmm3,nb204nf_three(%esp)

        movl nb204nf_argkrf(%ebp),%esi
        movl nb204nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb204nf_krf(%esp)
        movapd %xmm6,nb204nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb204nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb204nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb204nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204nf_qqMM(%esp)
        movapd %xmm4,nb204nf_qqMH(%esp)
        movapd %xmm5,nb204nf_qqHH(%esp)

_nb_kernel204nf_ia32_sse2.nb204nf_threadloop: 
        movl  nb204nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel204nf_ia32_sse2.nb204nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204nf_ia32_sse2.nb204nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204nf_n(%esp)
        movl %ebx,nb204nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204nf_ia32_sse2.nb204nf_outerstart
        jmp _nb_kernel204nf_ia32_sse2.nb204nf_end

_nb_kernel204nf_ia32_sse2.nb204nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204nf_nouter(%esp),%ebx
        movl %ebx,nb204nf_nouter(%esp)

_nb_kernel204nf_ia32_sse2.nb204nf_outer: 
        movl  nb204nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb204nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb204nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb204nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb204nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb204nf_ixH1(%esp)
        movapd %xmm4,nb204nf_iyH1(%esp)
        movapd %xmm5,nb204nf_izH1(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%eax,%ebx,8),%xmm0
        addsd 56(%eax,%ebx,8),%xmm1
        addsd 64(%eax,%ebx,8),%xmm2
        addsd 72(%eax,%ebx,8),%xmm3
        addsd 80(%eax,%ebx,8),%xmm4
        addsd 88(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb204nf_ixH2(%esp)
        movapd %xmm1,nb204nf_iyH2(%esp)
        movapd %xmm2,nb204nf_izH2(%esp)
        movapd %xmm3,nb204nf_ixM(%esp)
        movapd %xmm4,nb204nf_iyM(%esp)
        movapd %xmm5,nb204nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb204nf_vctot(%esp)

        movl  nb204nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb204nf_pos(%ebp),%esi
        movl  nb204nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb204nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb204nf_ninner(%esp),%ecx
        movl  %ecx,nb204nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb204nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel204nf_ia32_sse2.nb204nf_unroll_loop
        jmp   _nb_kernel204nf_ia32_sse2.nb204nf_checksingle
_nb_kernel204nf_ia32_sse2.nb204nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb204nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb204nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb204nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movhpd 24(%esi,%ebx,8),%xmm2
        movhpd 32(%esi,%ebx,8),%xmm3
        movhpd 40(%esi,%ebx,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%ebx,8),%xmm6
        movhpd 64(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb204nf_jxH1(%esp)
        movapd  %xmm3,nb204nf_jyH1(%esp)
        movapd  %xmm4,nb204nf_jzH1(%esp)
        movapd  %xmm5,nb204nf_jxH2(%esp)
        movapd  %xmm6,nb204nf_jyH2(%esp)
        movapd  %xmm7,nb204nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb204nf_jxM(%esp)
        movapd  %xmm3,nb204nf_jyM(%esp)
        movapd  %xmm4,nb204nf_jzM(%esp)

        movapd nb204nf_ixM(%esp),%xmm0
        movapd nb204nf_iyM(%esp),%xmm1
        movapd nb204nf_izM(%esp),%xmm2
        movapd nb204nf_ixM(%esp),%xmm3
        movapd nb204nf_iyM(%esp),%xmm4
        movapd nb204nf_izM(%esp),%xmm5
        subpd  nb204nf_jxM(%esp),%xmm0
        subpd  nb204nf_jyM(%esp),%xmm1
        subpd  nb204nf_jzM(%esp),%xmm2
        subpd  nb204nf_jxH1(%esp),%xmm3
        subpd  nb204nf_jyH1(%esp),%xmm4
        subpd  nb204nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqMM(%esp)
        movapd %xmm3,nb204nf_rsqMH1(%esp)

        movapd nb204nf_ixM(%esp),%xmm0
        movapd nb204nf_iyM(%esp),%xmm1
        movapd nb204nf_izM(%esp),%xmm2
        movapd nb204nf_ixH1(%esp),%xmm3
        movapd nb204nf_iyH1(%esp),%xmm4
        movapd nb204nf_izH1(%esp),%xmm5
        subpd  nb204nf_jxH2(%esp),%xmm0
        subpd  nb204nf_jyH2(%esp),%xmm1
        subpd  nb204nf_jzH2(%esp),%xmm2
        subpd  nb204nf_jxM(%esp),%xmm3
        subpd  nb204nf_jyM(%esp),%xmm4
        subpd  nb204nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqMH2(%esp)
        movapd %xmm3,nb204nf_rsqH1M(%esp)

        movapd nb204nf_ixH1(%esp),%xmm0
        movapd nb204nf_iyH1(%esp),%xmm1
        movapd nb204nf_izH1(%esp),%xmm2
        movapd nb204nf_ixH1(%esp),%xmm3
        movapd nb204nf_iyH1(%esp),%xmm4
        movapd nb204nf_izH1(%esp),%xmm5
        subpd  nb204nf_jxH1(%esp),%xmm0
        subpd  nb204nf_jyH1(%esp),%xmm1
        subpd  nb204nf_jzH1(%esp),%xmm2
        subpd  nb204nf_jxH2(%esp),%xmm3
        subpd  nb204nf_jyH2(%esp),%xmm4
        subpd  nb204nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqH1H1(%esp)
        movapd %xmm3,nb204nf_rsqH1H2(%esp)

        movapd nb204nf_ixH2(%esp),%xmm0
        movapd nb204nf_iyH2(%esp),%xmm1
        movapd nb204nf_izH2(%esp),%xmm2
        movapd nb204nf_ixH2(%esp),%xmm3
        movapd nb204nf_iyH2(%esp),%xmm4
        movapd nb204nf_izH2(%esp),%xmm5
        subpd  nb204nf_jxM(%esp),%xmm0
        subpd  nb204nf_jyM(%esp),%xmm1
        subpd  nb204nf_jzM(%esp),%xmm2
        subpd  nb204nf_jxH1(%esp),%xmm3
        subpd  nb204nf_jyH1(%esp),%xmm4
        subpd  nb204nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqH2M(%esp)
        movapd %xmm4,nb204nf_rsqH2H1(%esp)

        movapd nb204nf_ixH2(%esp),%xmm0
        movapd nb204nf_iyH2(%esp),%xmm1
        movapd nb204nf_izH2(%esp),%xmm2
        subpd  nb204nf_jxH2(%esp),%xmm0
        subpd  nb204nf_jyH2(%esp),%xmm1
        subpd  nb204nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb204nf_rsqH2H2(%esp)

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
        movapd  nb204nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb204nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH2H2(%esp)
        movapd %xmm5,nb204nf_rinvH2H1(%esp)

        movapd nb204nf_rsqMM(%esp),%xmm0
        movapd nb204nf_rsqMH1(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb204nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb204nf_rinvMM(%esp)
        movapd %xmm5,nb204nf_rinvMH1(%esp)

        movapd nb204nf_rsqMH2(%esp),%xmm0
        movapd nb204nf_rsqH1M(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb204nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvMH2(%esp)
        movapd %xmm5,nb204nf_rinvH1M(%esp)

        movapd nb204nf_rsqH1H1(%esp),%xmm0
        movapd nb204nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb204nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH1H1(%esp)
        movapd %xmm5,nb204nf_rinvH1H2(%esp)

        movapd nb204nf_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb204nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb204nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb204nf_krf(%esp),%xmm6
        mulpd  nb204nf_rsqMM(%esp),%xmm6        ## xmm5=krsq 
        addpd  nb204nf_rinvMM(%esp),%xmm6       ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm6

        mulpd  nb204nf_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  nb204nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## M-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulpd  nb204nf_rsqMH1(%esp),%xmm5       ## xmm5=krsq 
        addpd  nb204nf_rinvMH1(%esp),%xmm5      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm5

        mulpd  nb204nf_qqMH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## M-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulpd  nb204nf_rsqMH2(%esp),%xmm7       ## xmm5=krsq 
        addpd  nb204nf_rinvMH2(%esp),%xmm7      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm7

        mulpd  nb204nf_qqMH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-M interaction 
        movapd nb204nf_krf(%esp),%xmm4
        mulpd  nb204nf_rsqH1M(%esp),%xmm4       ## xmm5=krsq 
        addpd  nb204nf_rinvH1M(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm4

        mulpd  nb204nf_qqMH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulpd  nb204nf_rsqH1H1(%esp),%xmm5      ## xmm5=krsq 
        addpd  nb204nf_rinvH1H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm5

        mulpd  nb204nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulpd  nb204nf_rsqH1H2(%esp),%xmm7      ## xmm5=krsq 
        addpd  nb204nf_rinvH1H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm7

        mulpd  nb204nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-M interaction 
        movapd nb204nf_krf(%esp),%xmm4
        mulpd  nb204nf_rsqH2M(%esp),%xmm4       ## xmm5=krsq 
        addpd  nb204nf_rinvH2M(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm4

        mulpd  nb204nf_qqMH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulpd  nb204nf_rsqH2H1(%esp),%xmm5      ## xmm5=krsq 
        addpd  nb204nf_rinvH2H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm5

        mulpd  nb204nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulpd  nb204nf_rsqH2H2(%esp),%xmm7      ## xmm5=krsq 
        addpd  nb204nf_rinvH2H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb204nf_crf(%esp),%xmm7

        mulpd  nb204nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 
        movapd %xmm6,nb204nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb204nf_innerk(%esp)
        jl    _nb_kernel204nf_ia32_sse2.nb204nf_checksingle
        jmp   _nb_kernel204nf_ia32_sse2.nb204nf_unroll_loop
_nb_kernel204nf_ia32_sse2.nb204nf_checksingle: 
        movl  nb204nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel204nf_ia32_sse2.nb204nf_dosingle
        jmp   _nb_kernel204nf_ia32_sse2.nb204nf_updateouterdata
_nb_kernel204nf_ia32_sse2.nb204nf_dosingle: 
        movl  nb204nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb204nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb204nf_jxH1(%esp)
        movapd  %xmm3,nb204nf_jyH1(%esp)
        movapd  %xmm4,nb204nf_jzH1(%esp)
        movapd  %xmm5,nb204nf_jxH2(%esp)
        movapd  %xmm6,nb204nf_jyH2(%esp)
        movapd  %xmm7,nb204nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb204nf_jxM(%esp)
        movapd  %xmm3,nb204nf_jyM(%esp)
        movapd  %xmm4,nb204nf_jzM(%esp)

        movapd nb204nf_ixM(%esp),%xmm0
        movapd nb204nf_iyM(%esp),%xmm1
        movapd nb204nf_izM(%esp),%xmm2
        movapd nb204nf_ixM(%esp),%xmm3
        movapd nb204nf_iyM(%esp),%xmm4
        movapd nb204nf_izM(%esp),%xmm5
        subsd  nb204nf_jxM(%esp),%xmm0
        subsd  nb204nf_jyM(%esp),%xmm1
        subsd  nb204nf_jzM(%esp),%xmm2
        subsd  nb204nf_jxH1(%esp),%xmm3
        subsd  nb204nf_jyH1(%esp),%xmm4
        subsd  nb204nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqMM(%esp)
        movapd %xmm3,nb204nf_rsqMH1(%esp)

        movapd nb204nf_ixM(%esp),%xmm0
        movapd nb204nf_iyM(%esp),%xmm1
        movapd nb204nf_izM(%esp),%xmm2
        movapd nb204nf_ixH1(%esp),%xmm3
        movapd nb204nf_iyH1(%esp),%xmm4
        movapd nb204nf_izH1(%esp),%xmm5
        subsd  nb204nf_jxH2(%esp),%xmm0
        subsd  nb204nf_jyH2(%esp),%xmm1
        subsd  nb204nf_jzH2(%esp),%xmm2
        subsd  nb204nf_jxM(%esp),%xmm3
        subsd  nb204nf_jyM(%esp),%xmm4
        subsd  nb204nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqMH2(%esp)
        movapd %xmm3,nb204nf_rsqH1M(%esp)

        movapd nb204nf_ixH1(%esp),%xmm0
        movapd nb204nf_iyH1(%esp),%xmm1
        movapd nb204nf_izH1(%esp),%xmm2
        movapd nb204nf_ixH1(%esp),%xmm3
        movapd nb204nf_iyH1(%esp),%xmm4
        movapd nb204nf_izH1(%esp),%xmm5
        subsd  nb204nf_jxH1(%esp),%xmm0
        subsd  nb204nf_jyH1(%esp),%xmm1
        subsd  nb204nf_jzH1(%esp),%xmm2
        subsd  nb204nf_jxH2(%esp),%xmm3
        subsd  nb204nf_jyH2(%esp),%xmm4
        subsd  nb204nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqH1H1(%esp)
        movapd %xmm3,nb204nf_rsqH1H2(%esp)

        movapd nb204nf_ixH2(%esp),%xmm0
        movapd nb204nf_iyH2(%esp),%xmm1
        movapd nb204nf_izH2(%esp),%xmm2
        movapd nb204nf_ixH2(%esp),%xmm3
        movapd nb204nf_iyH2(%esp),%xmm4
        movapd nb204nf_izH2(%esp),%xmm5
        subsd  nb204nf_jxM(%esp),%xmm0
        subsd  nb204nf_jyM(%esp),%xmm1
        subsd  nb204nf_jzM(%esp),%xmm2
        subsd  nb204nf_jxH1(%esp),%xmm3
        subsd  nb204nf_jyH1(%esp),%xmm4
        subsd  nb204nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb204nf_rsqH2M(%esp)
        movapd %xmm4,nb204nf_rsqH2H1(%esp)

        movapd nb204nf_ixH2(%esp),%xmm0
        movapd nb204nf_iyH2(%esp),%xmm1
        movapd nb204nf_izH2(%esp),%xmm2
        subsd  nb204nf_jxH2(%esp),%xmm0
        subsd  nb204nf_jyH2(%esp),%xmm1
        subsd  nb204nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb204nf_rsqH2H2(%esp)

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
        movapd  nb204nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb204nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH2H2(%esp)
        movapd %xmm5,nb204nf_rinvH2H1(%esp)

        movapd nb204nf_rsqMM(%esp),%xmm0
        movapd nb204nf_rsqMH1(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb204nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb204nf_rinvMM(%esp)
        movapd %xmm5,nb204nf_rinvMH1(%esp)

        movapd nb204nf_rsqMH2(%esp),%xmm0
        movapd nb204nf_rsqH1M(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb204nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvMH2(%esp)
        movapd %xmm5,nb204nf_rinvH1M(%esp)

        movapd nb204nf_rsqH1H1(%esp),%xmm0
        movapd nb204nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb204nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb204nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb204nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb204nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb204nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb204nf_rinvH1H1(%esp)
        movapd %xmm5,nb204nf_rinvH1H2(%esp)

        movapd nb204nf_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb204nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb204nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb204nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb204nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb204nf_krf(%esp),%xmm6
        mulsd  nb204nf_rsqMM(%esp),%xmm6        ## xmm5=krsq 
        addsd  nb204nf_rinvMM(%esp),%xmm6       ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm6

        mulsd  nb204nf_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  nb204nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## M-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulsd  nb204nf_rsqMH1(%esp),%xmm5       ## xmm5=krsq 
        addsd  nb204nf_rinvMH1(%esp),%xmm5      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm5

        mulsd  nb204nf_qqMH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## M-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulsd  nb204nf_rsqMH2(%esp),%xmm7       ## xmm5=krsq 
        addsd  nb204nf_rinvMH2(%esp),%xmm7      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm7

        mulsd  nb204nf_qqMH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-M interaction 
        movapd nb204nf_krf(%esp),%xmm4
        mulsd  nb204nf_rsqH1M(%esp),%xmm4       ## xmm5=krsq 
        addsd  nb204nf_rinvH1M(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm4

        mulsd  nb204nf_qqMH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulsd  nb204nf_rsqH1H1(%esp),%xmm5      ## xmm5=krsq 
        addsd  nb204nf_rinvH1H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm5

        mulsd  nb204nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulsd  nb204nf_rsqH1H2(%esp),%xmm7      ## xmm5=krsq 
        addsd  nb204nf_rinvH1H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm7

        mulsd  nb204nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-M interaction 
        movapd nb204nf_krf(%esp),%xmm4
        mulsd  nb204nf_rsqH2M(%esp),%xmm4       ## xmm5=krsq 
        addsd  nb204nf_rinvH2M(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm4

        mulsd  nb204nf_qqMH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb204nf_krf(%esp),%xmm5
        mulsd  nb204nf_rsqH2H1(%esp),%xmm5      ## xmm5=krsq 
        addsd  nb204nf_rinvH2H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm5

        mulsd  nb204nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb204nf_krf(%esp),%xmm7
        mulsd  nb204nf_rsqH2H2(%esp),%xmm7      ## xmm5=krsq 
        addsd  nb204nf_rinvH2H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb204nf_crf(%esp),%xmm7

        mulsd  nb204nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 
        movlpd %xmm6,nb204nf_vctot(%esp)

_nb_kernel204nf_ia32_sse2.nb204nf_updateouterdata: 
        ## get n from stack
        movl nb204nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb204nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb204nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb204nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb204nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204nf_ia32_sse2.nb204nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204nf_n(%esp)
        jmp _nb_kernel204nf_ia32_sse2.nb204nf_outer
_nb_kernel204nf_ia32_sse2.nb204nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204nf_ia32_sse2.nb204nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204nf_ia32_sse2.nb204nf_threadloop
_nb_kernel204nf_ia32_sse2.nb204nf_end: 
        emms

        movl nb204nf_nouter(%esp),%eax
        movl nb204nf_ninner(%esp),%ebx
        movl nb204nf_outeriter(%ebp),%ecx
        movl nb204nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb204nf_salign(%esp),%eax
        addl %eax,%esp
        addl $744,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



