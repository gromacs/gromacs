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



.globl nb_kernel214_ia32_sse2
.globl _nb_kernel214_ia32_sse2
nb_kernel214_ia32_sse2: 
_nb_kernel214_ia32_sse2:        
.set nb214_p_nri, 8
.set nb214_iinr, 12
.set nb214_jindex, 16
.set nb214_jjnr, 20
.set nb214_shift, 24
.set nb214_shiftvec, 28
.set nb214_fshift, 32
.set nb214_gid, 36
.set nb214_pos, 40
.set nb214_faction, 44
.set nb214_charge, 48
.set nb214_p_facel, 52
.set nb214_argkrf, 56
.set nb214_argcrf, 60
.set nb214_Vc, 64
.set nb214_type, 68
.set nb214_p_ntype, 72
.set nb214_vdwparam, 76
.set nb214_Vvdw, 80
.set nb214_p_tabscale, 84
.set nb214_VFtab, 88
.set nb214_invsqrta, 92
.set nb214_dvda, 96
.set nb214_p_gbtabscale, 100
.set nb214_GBtab, 104
.set nb214_p_nthreads, 108
.set nb214_count, 112
.set nb214_mtx, 116
.set nb214_outeriter, 120
.set nb214_inneriter, 124
.set nb214_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb214_ixO, 0
.set nb214_iyO, 16
.set nb214_izO, 32
.set nb214_ixH1, 48
.set nb214_iyH1, 64
.set nb214_izH1, 80
.set nb214_ixH2, 96
.set nb214_iyH2, 112
.set nb214_izH2, 128
.set nb214_ixM, 144
.set nb214_iyM, 160
.set nb214_izM, 176
.set nb214_jxO, 192
.set nb214_jyO, 208
.set nb214_jzO, 224
.set nb214_jxH1, 240
.set nb214_jyH1, 256
.set nb214_jzH1, 272
.set nb214_jxH2, 288
.set nb214_jyH2, 304
.set nb214_jzH2, 320
.set nb214_jxM, 336
.set nb214_jyM, 352
.set nb214_jzM, 368
.set nb214_dxOO, 384
.set nb214_dyOO, 400
.set nb214_dzOO, 416
.set nb214_dxH1H1, 432
.set nb214_dyH1H1, 448
.set nb214_dzH1H1, 464
.set nb214_dxH1H2, 480
.set nb214_dyH1H2, 496
.set nb214_dzH1H2, 512
.set nb214_dxH1M, 528
.set nb214_dyH1M, 544
.set nb214_dzH1M, 560
.set nb214_dxH2H1, 576
.set nb214_dyH2H1, 592
.set nb214_dzH2H1, 608
.set nb214_dxH2H2, 624
.set nb214_dyH2H2, 640
.set nb214_dzH2H2, 656
.set nb214_dxH2M, 672
.set nb214_dyH2M, 688
.set nb214_dzH2M, 704
.set nb214_dxMH1, 720
.set nb214_dyMH1, 736
.set nb214_dzMH1, 752
.set nb214_dxMH2, 768
.set nb214_dyMH2, 784
.set nb214_dzMH2, 800
.set nb214_dxMM, 816
.set nb214_dyMM, 832
.set nb214_dzMM, 848
.set nb214_qqMM, 864
.set nb214_qqMH, 880
.set nb214_qqHH, 896
.set nb214_two, 912
.set nb214_c6, 944
.set nb214_c12, 960
.set nb214_vctot, 976
.set nb214_Vvdwtot, 992
.set nb214_fixO, 1008
.set nb214_fiyO, 1024
.set nb214_fizO, 1040
.set nb214_fixH1, 1056
.set nb214_fiyH1, 1072
.set nb214_fizH1, 1088
.set nb214_fixH2, 1104
.set nb214_fiyH2, 1120
.set nb214_fizH2, 1136
.set nb214_fixM, 1152
.set nb214_fiyM, 1168
.set nb214_fizM, 1184
.set nb214_fjxO, 1200
.set nb214_fjyO, 1216
.set nb214_fjzO, 1232
.set nb214_fjxH1, 1248
.set nb214_fjyH1, 1264
.set nb214_fjzH1, 1280
.set nb214_fjxH2, 1296
.set nb214_fjyH2, 1312
.set nb214_fjzH2, 1328
.set nb214_fjxM, 1344
.set nb214_fjyM, 1360
.set nb214_fjzM, 1376
.set nb214_half, 1392
.set nb214_three, 1408
.set nb214_six, 1424
.set nb214_twelve, 1440
.set nb214_rsqOO, 1456
.set nb214_rsqH1H1, 1472
.set nb214_rsqH1H2, 1488
.set nb214_rsqH1M, 1504
.set nb214_rsqH2H1, 1520
.set nb214_rsqH2H2, 1536
.set nb214_rsqH2M, 1552
.set nb214_rsqMH1, 1568
.set nb214_rsqMH2, 1584
.set nb214_rsqMM, 1600
.set nb214_rinvsqOO, 1616
.set nb214_rinvH1H1, 1632
.set nb214_rinvH1H2, 1648
.set nb214_rinvH1M, 1664
.set nb214_rinvH2H1, 1680
.set nb214_rinvH2H2, 1696
.set nb214_rinvH2M, 1712
.set nb214_rinvMH1, 1728
.set nb214_rinvMH2, 1744
.set nb214_rinvMM, 1760
.set nb214_krf, 1776
.set nb214_crf, 1792
.set nb214_is3, 1808
.set nb214_ii3, 1812
.set nb214_innerjjnr, 1816
.set nb214_innerk, 1820
.set nb214_n, 1824
.set nb214_nn1, 1828
.set nb214_nri, 1832
.set nb214_nouter, 1836
.set nb214_ninner, 1840
.set nb214_salign, 1844
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1848,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb214_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb214_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb214_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb214_nouter(%esp)
        movl %eax,nb214_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb214_half(%esp)
        movl %ebx,nb214_half+4(%esp)
        movsd nb214_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb214_half(%esp)
        movapd %xmm2,nb214_two(%esp)
        movapd %xmm3,nb214_three(%esp)
        movapd %xmm4,nb214_six(%esp)
        movapd %xmm5,nb214_twelve(%esp)

        movl nb214_argkrf(%ebp),%esi
        movl nb214_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb214_krf(%esp)
        movapd %xmm6,nb214_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb214_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb214_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb214_p_facel(%ebp),%esi
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
        movapd %xmm3,nb214_qqMM(%esp)
        movapd %xmm4,nb214_qqMH(%esp)
        movapd %xmm5,nb214_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb214_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb214_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb214_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb214_c6(%esp)
        movapd %xmm1,nb214_c12(%esp)

_nb_kernel214_ia32_sse2.nb214_threadloop: 
        movl  nb214_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel214_ia32_sse2.nb214_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel214_ia32_sse2.nb214_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb214_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb214_n(%esp)
        movl %ebx,nb214_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel214_ia32_sse2.nb214_outerstart
        jmp _nb_kernel214_ia32_sse2.nb214_end

_nb_kernel214_ia32_sse2.nb214_outerstart: 
        ## ebx contains number of outer iterations
        addl nb214_nouter(%esp),%ebx
        movl %ebx,nb214_nouter(%esp)

_nb_kernel214_ia32_sse2.nb214_outer: 
        movl  nb214_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb214_is3(%esp)      ## store is3 

        movl  nb214_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb214_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb214_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb214_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb214_ixO(%esp)
        movapd %xmm4,nb214_iyO(%esp)
        movapd %xmm5,nb214_izO(%esp)
        movapd %xmm6,nb214_ixH1(%esp)
        movapd %xmm7,nb214_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb214_izH1(%esp)
        movapd %xmm0,nb214_ixH2(%esp)
        movapd %xmm1,nb214_iyH2(%esp)
        movapd %xmm2,nb214_izH2(%esp)
        movapd %xmm3,nb214_ixM(%esp)
        movapd %xmm4,nb214_iyM(%esp)
        movapd %xmm5,nb214_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb214_vctot(%esp)
        movapd %xmm4,nb214_Vvdwtot(%esp)
        movapd %xmm4,nb214_fixO(%esp)
        movapd %xmm4,nb214_fiyO(%esp)
        movapd %xmm4,nb214_fizO(%esp)
        movapd %xmm4,nb214_fixH1(%esp)
        movapd %xmm4,nb214_fiyH1(%esp)
        movapd %xmm4,nb214_fizH1(%esp)
        movapd %xmm4,nb214_fixH2(%esp)
        movapd %xmm4,nb214_fiyH2(%esp)
        movapd %xmm4,nb214_fizH2(%esp)
        movapd %xmm4,nb214_fixM(%esp)
        movapd %xmm4,nb214_fiyM(%esp)
        movapd %xmm4,nb214_fizM(%esp)

        movl  nb214_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb214_pos(%ebp),%esi
        movl  nb214_faction(%ebp),%edi
        movl  nb214_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb214_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb214_ninner(%esp),%ecx
        movl  %ecx,nb214_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb214_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel214_ia32_sse2.nb214_unroll_loop
        jmp   _nb_kernel214_ia32_sse2.nb214_checksingle
_nb_kernel214_ia32_sse2.nb214_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb214_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb214_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb214_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb214_jxO(%esp)
        movapd  %xmm1,nb214_jyO(%esp)
        movapd  %xmm3,nb214_jzO(%esp)
        movapd  %xmm4,nb214_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb214_jyH1(%esp)
        movapd  %xmm1,nb214_jzH1(%esp)
        movapd  %xmm3,nb214_jxH2(%esp)
        movapd  %xmm4,nb214_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb214_jzH2(%esp)
        movapd  %xmm1,nb214_jxM(%esp)
        movapd  %xmm3,nb214_jyM(%esp)
        movapd  %xmm4,nb214_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb214_ixO(%esp),%xmm0
        movapd nb214_iyO(%esp),%xmm1
        movapd nb214_izO(%esp),%xmm2
        movapd nb214_ixH1(%esp),%xmm3
        movapd nb214_iyH1(%esp),%xmm4
        movapd nb214_izH1(%esp),%xmm5
        subpd  nb214_jxO(%esp),%xmm0
        subpd  nb214_jyO(%esp),%xmm1
        subpd  nb214_jzO(%esp),%xmm2
        subpd  nb214_jxH1(%esp),%xmm3
        subpd  nb214_jyH1(%esp),%xmm4
        subpd  nb214_jzH1(%esp),%xmm5
        movapd %xmm0,nb214_dxOO(%esp)
        movapd %xmm1,nb214_dyOO(%esp)
        movapd %xmm2,nb214_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH1H1(%esp)
        movapd %xmm4,nb214_dyH1H1(%esp)
        movapd %xmm5,nb214_dzH1H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqOO(%esp)
        movapd %xmm3,nb214_rsqH1H1(%esp)

        movapd nb214_ixH1(%esp),%xmm0
        movapd nb214_iyH1(%esp),%xmm1
        movapd nb214_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214_jxH2(%esp),%xmm0
        subpd  nb214_jyH2(%esp),%xmm1
        subpd  nb214_jzH2(%esp),%xmm2
        subpd  nb214_jxM(%esp),%xmm3
        subpd  nb214_jyM(%esp),%xmm4
        subpd  nb214_jzM(%esp),%xmm5
        movapd %xmm0,nb214_dxH1H2(%esp)
        movapd %xmm1,nb214_dyH1H2(%esp)
        movapd %xmm2,nb214_dzH1H2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH1M(%esp)
        movapd %xmm4,nb214_dyH1M(%esp)
        movapd %xmm5,nb214_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqH1H2(%esp)
        movapd %xmm3,nb214_rsqH1M(%esp)

        movapd nb214_ixH2(%esp),%xmm0
        movapd nb214_iyH2(%esp),%xmm1
        movapd nb214_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214_jxH1(%esp),%xmm0
        subpd  nb214_jyH1(%esp),%xmm1
        subpd  nb214_jzH1(%esp),%xmm2
        subpd  nb214_jxH2(%esp),%xmm3
        subpd  nb214_jyH2(%esp),%xmm4
        subpd  nb214_jzH2(%esp),%xmm5
        movapd %xmm0,nb214_dxH2H1(%esp)
        movapd %xmm1,nb214_dyH2H1(%esp)
        movapd %xmm2,nb214_dzH2H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH2H2(%esp)
        movapd %xmm4,nb214_dyH2H2(%esp)
        movapd %xmm5,nb214_dzH2H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqH2H1(%esp)
        movapd %xmm3,nb214_rsqH2H2(%esp)

        movapd nb214_ixH2(%esp),%xmm0
        movapd nb214_iyH2(%esp),%xmm1
        movapd nb214_izH2(%esp),%xmm2
        movapd nb214_ixM(%esp),%xmm3
        movapd nb214_iyM(%esp),%xmm4
        movapd nb214_izM(%esp),%xmm5
        subpd  nb214_jxM(%esp),%xmm0
        subpd  nb214_jyM(%esp),%xmm1
        subpd  nb214_jzM(%esp),%xmm2
        subpd  nb214_jxH1(%esp),%xmm3
        subpd  nb214_jyH1(%esp),%xmm4
        subpd  nb214_jzH1(%esp),%xmm5
        movapd %xmm0,nb214_dxH2M(%esp)
        movapd %xmm1,nb214_dyH2M(%esp)
        movapd %xmm2,nb214_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxMH1(%esp)
        movapd %xmm4,nb214_dyMH1(%esp)
        movapd %xmm5,nb214_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb214_rsqH2M(%esp)
        movapd %xmm4,nb214_rsqMH1(%esp)

        movapd nb214_ixM(%esp),%xmm0
        movapd nb214_iyM(%esp),%xmm1
        movapd nb214_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214_jxH2(%esp),%xmm0
        subpd  nb214_jyH2(%esp),%xmm1
        subpd  nb214_jzH2(%esp),%xmm2
        subpd  nb214_jxM(%esp),%xmm3
        subpd  nb214_jyM(%esp),%xmm4
        subpd  nb214_jzM(%esp),%xmm5
        movapd %xmm0,nb214_dxMH2(%esp)
        movapd %xmm1,nb214_dyMH2(%esp)
        movapd %xmm2,nb214_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxMM(%esp)
        movapd %xmm4,nb214_dyMM(%esp)
        movapd %xmm5,nb214_dzMM(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb214_rsqMH2(%esp)
        movapd %xmm4,nb214_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb214_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm3   ## iter1 
        mulpd   nb214_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm1   ## rinv 
        mulpd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvMH2(%esp)
        movapd %xmm5,nb214_rinvMM(%esp)

        movapd nb214_rsqOO(%esp),%xmm0
        movapd nb214_rsqH1H1(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb214_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm1   ## rinv 
        mulpd   nb214_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb214_rinvsqOO(%esp)
        movapd %xmm5,nb214_rinvH1H1(%esp)

        movapd nb214_rsqH1H2(%esp),%xmm0
        movapd nb214_rsqH1M(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm3   ## iter1 
        mulpd   nb214_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm1   ## rinv 
        mulpd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvH1H2(%esp)
        movapd %xmm5,nb214_rinvH1M(%esp)

        movapd nb214_rsqH2H1(%esp),%xmm0
        movapd nb214_rsqH2H2(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm3   ## iter1a 
        mulpd   nb214_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm1   ## rinv 
        mulpd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvH2H1(%esp)
        movapd %xmm5,nb214_rinvH2H2(%esp)

        movapd nb214_rsqMH1(%esp),%xmm0
        movapd nb214_rsqH2M(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm3   ## iter1a 
        mulpd   nb214_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214_half(%esp),%xmm1   ## rinv 
        mulpd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvMH1(%esp)
        movapd %xmm5,nb214_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb214_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb214_c6(%esp),%xmm1
        mulpd  nb214_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb214_Vvdwtot(%esp),%xmm3
        mulpd  nb214_six(%esp),%xmm1
        mulpd  nb214_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm3,nb214_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1 ## fscal

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb214_dxOO(%esp),%xmm0
        mulpd nb214_dyOO(%esp),%xmm1
        mulpd nb214_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixO(%esp),%xmm0
        addpd nb214_fiyO(%esp),%xmm1
        addpd nb214_fizO(%esp),%xmm2
        movapd %xmm3,nb214_fjxO(%esp)
        movapd %xmm4,nb214_fjyO(%esp)
        movapd %xmm5,nb214_fjzO(%esp)
        movapd %xmm0,nb214_fixO(%esp)
        movapd %xmm1,nb214_fiyO(%esp)
        movapd %xmm2,nb214_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb214_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb214_dxH1H1(%esp),%xmm0
        mulpd nb214_dyH1H1(%esp),%xmm1
        mulpd nb214_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH1(%esp),%xmm0
        addpd nb214_fiyH1(%esp),%xmm1
        addpd nb214_fizH1(%esp),%xmm2
        movapd %xmm3,nb214_fjxH1(%esp)
        movapd %xmm4,nb214_fjyH1(%esp)
        movapd %xmm5,nb214_fjzH1(%esp)
        movapd %xmm0,nb214_fixH1(%esp)
        movapd %xmm1,nb214_fiyH1(%esp)
        movapd %xmm2,nb214_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb214_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb214_dxH1H2(%esp),%xmm0
        mulpd nb214_dyH1H2(%esp),%xmm1
        mulpd nb214_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH1(%esp),%xmm0
        addpd nb214_fiyH1(%esp),%xmm1
        addpd nb214_fizH1(%esp),%xmm2
        movapd %xmm3,nb214_fjxH2(%esp)
        movapd %xmm4,nb214_fjyH2(%esp)
        movapd %xmm5,nb214_fjzH2(%esp)
        movapd %xmm0,nb214_fixH1(%esp)
        movapd %xmm1,nb214_fiyH1(%esp)
        movapd %xmm2,nb214_fizH1(%esp)

        ## H1-M interaction 
        movapd nb214_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb214_dxH1M(%esp),%xmm0
        mulpd nb214_dyH1M(%esp),%xmm1
        mulpd nb214_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH1(%esp),%xmm0
        addpd nb214_fiyH1(%esp),%xmm1
        addpd nb214_fizH1(%esp),%xmm2
        movapd %xmm3,nb214_fjxM(%esp)
        movapd %xmm4,nb214_fjyM(%esp)
        movapd %xmm5,nb214_fjzM(%esp)
        movapd %xmm0,nb214_fixH1(%esp)
        movapd %xmm1,nb214_fiyH1(%esp)
        movapd %xmm2,nb214_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb214_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH1(%esp),%xmm3
        movapd nb214_fjyH1(%esp),%xmm4
        movapd nb214_fjzH1(%esp),%xmm5
        mulpd nb214_dxH2H1(%esp),%xmm0
        mulpd nb214_dyH2H1(%esp),%xmm1
        mulpd nb214_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH2(%esp),%xmm0
        addpd nb214_fiyH2(%esp),%xmm1
        addpd nb214_fizH2(%esp),%xmm2
        movapd %xmm3,nb214_fjxH1(%esp)
        movapd %xmm4,nb214_fjyH1(%esp)
        movapd %xmm5,nb214_fjzH1(%esp)
        movapd %xmm0,nb214_fixH2(%esp)
        movapd %xmm1,nb214_fiyH2(%esp)
        movapd %xmm2,nb214_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb214_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH2(%esp),%xmm3
        movapd nb214_fjyH2(%esp),%xmm4
        movapd nb214_fjzH2(%esp),%xmm5
        mulpd nb214_dxH2H2(%esp),%xmm0
        mulpd nb214_dyH2H2(%esp),%xmm1
        mulpd nb214_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH2(%esp),%xmm0
        addpd nb214_fiyH2(%esp),%xmm1
        addpd nb214_fizH2(%esp),%xmm2
        movapd %xmm3,nb214_fjxH2(%esp)
        movapd %xmm4,nb214_fjyH2(%esp)
        movapd %xmm5,nb214_fjzH2(%esp)
        movapd %xmm0,nb214_fixH2(%esp)
        movapd %xmm1,nb214_fiyH2(%esp)
        movapd %xmm2,nb214_fizH2(%esp)

        ## H2-M interaction 
        movapd nb214_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxM(%esp),%xmm3
        movapd nb214_fjyM(%esp),%xmm4
        movapd nb214_fjzM(%esp),%xmm5
        mulpd nb214_dxH2M(%esp),%xmm0
        mulpd nb214_dyH2M(%esp),%xmm1
        mulpd nb214_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixH2(%esp),%xmm0
        addpd nb214_fiyH2(%esp),%xmm1
        addpd nb214_fizH2(%esp),%xmm2
        movapd %xmm3,nb214_fjxM(%esp)
        movapd %xmm4,nb214_fjyM(%esp)
        movapd %xmm5,nb214_fjzM(%esp)
        movapd %xmm0,nb214_fixH2(%esp)
        movapd %xmm1,nb214_fiyH2(%esp)
        movapd %xmm2,nb214_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb214_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH1(%esp),%xmm3
        movapd nb214_fjyH1(%esp),%xmm4
        movapd nb214_fjzH1(%esp),%xmm5
        mulpd nb214_dxMH1(%esp),%xmm0
        mulpd nb214_dyMH1(%esp),%xmm1
        mulpd nb214_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixM(%esp),%xmm0
        addpd nb214_fiyM(%esp),%xmm1
        addpd nb214_fizM(%esp),%xmm2
        movapd %xmm3,nb214_fjxH1(%esp)
        movapd %xmm4,nb214_fjyH1(%esp)
        movapd %xmm5,nb214_fjzH1(%esp)
        movapd %xmm0,nb214_fixM(%esp)
        movapd %xmm1,nb214_fiyM(%esp)
        movapd %xmm2,nb214_fizM(%esp)

        ## M-H2 interaction 
        movapd nb214_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH2(%esp),%xmm3
        movapd nb214_fjyH2(%esp),%xmm4
        movapd nb214_fjzH2(%esp),%xmm5
        mulpd nb214_dxMH2(%esp),%xmm0
        mulpd nb214_dyMH2(%esp),%xmm1
        mulpd nb214_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixM(%esp),%xmm0
        addpd nb214_fiyM(%esp),%xmm1
        addpd nb214_fizM(%esp),%xmm2
        movapd %xmm3,nb214_fjxH2(%esp)
        movapd %xmm4,nb214_fjyH2(%esp)
        movapd %xmm5,nb214_fjzH2(%esp)
        movapd %xmm0,nb214_fixM(%esp)
        movapd %xmm1,nb214_fiyM(%esp)
        movapd %xmm2,nb214_fizM(%esp)

        ## M-M interaction 
        movapd nb214_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb214_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulpd  nb214_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214_crf(%esp),%xmm6
        mulpd  nb214_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb214_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb214_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0
        movapd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxM(%esp),%xmm3
        movapd nb214_fjyM(%esp),%xmm4
        movapd nb214_fjzM(%esp),%xmm5
        mulpd nb214_dxMM(%esp),%xmm0
        mulpd nb214_dyMM(%esp),%xmm1
        mulpd nb214_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb214_fixM(%esp),%xmm0
        addpd nb214_fiyM(%esp),%xmm1
        addpd nb214_fizM(%esp),%xmm2
        movapd %xmm3,nb214_fjxM(%esp)
        movapd %xmm4,nb214_fjyM(%esp)
        movapd %xmm5,nb214_fjzM(%esp)
        movapd %xmm0,nb214_fixM(%esp)
        movapd %xmm1,nb214_fiyM(%esp)
        movapd %xmm2,nb214_fizM(%esp)

        movl nb214_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - transpose fjxO, fjyO and fjzO, fjxH1
        movapd nb214_fjxO(%esp),%xmm0
        movapd nb214_fjyO(%esp),%xmm1
        movapd nb214_fjzO(%esp),%xmm2
        movapd nb214_fjxH1(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd (%edi,%eax,8),%xmm1
        movlpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%eax,8),%xmm1
        movhpd 8(%edi,%ebx,8),%xmm3
        movlpd 16(%edi,%eax,8),%xmm6
        movlpd 16(%edi,%ebx,8),%xmm7
        movhpd 24(%edi,%eax,8),%xmm6
        movhpd 24(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,(%edi,%eax,8)
        movlpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%eax,8)
        movhpd %xmm3,8(%edi,%ebx,8)
        movlpd %xmm6,16(%edi,%eax,8)
        movlpd %xmm7,16(%edi,%ebx,8)
        movhpd %xmm6,24(%edi,%eax,8)
        movhpd %xmm7,24(%edi,%ebx,8)

        ## Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
        movapd nb214_fjyH1(%esp),%xmm0
        movapd nb214_fjzH1(%esp),%xmm1
        movapd nb214_fjxH2(%esp),%xmm2
        movapd nb214_fjyH2(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 32(%edi,%ebx,8),%xmm3
        movhpd 40(%edi,%eax,8),%xmm1
        movhpd 40(%edi,%ebx,8),%xmm3
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 48(%edi,%ebx,8),%xmm7
        movhpd 56(%edi,%eax,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm3,32(%edi,%ebx,8)
        movhpd %xmm1,40(%edi,%eax,8)
        movhpd %xmm3,40(%edi,%ebx,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,48(%edi,%ebx,8)
        movhpd %xmm6,56(%edi,%eax,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        ## Step3 - transpose fjzH2, fjxM and fjyM, fjzM
        movapd nb214_fjzH2(%esp),%xmm0
        movapd nb214_fjxM(%esp),%xmm1
        movapd nb214_fjyM(%esp),%xmm2
        movapd nb214_fjzM(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 64(%edi,%eax,8),%xmm1
        movlpd 64(%edi,%ebx,8),%xmm3
        movhpd 72(%edi,%eax,8),%xmm1
        movhpd 72(%edi,%ebx,8),%xmm3
        movlpd 80(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%ebx,8),%xmm7
        movhpd 88(%edi,%eax,8),%xmm6
        movhpd 88(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,64(%edi,%eax,8)
        movlpd %xmm3,64(%edi,%ebx,8)
        movhpd %xmm1,72(%edi,%eax,8)
        movhpd %xmm3,72(%edi,%ebx,8)
        movlpd %xmm6,80(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%ebx,8)
        movhpd %xmm6,88(%edi,%eax,8)
        movhpd %xmm7,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb214_innerk(%esp)
        jl    _nb_kernel214_ia32_sse2.nb214_checksingle
        jmp   _nb_kernel214_ia32_sse2.nb214_unroll_loop
_nb_kernel214_ia32_sse2.nb214_checksingle: 
        movl  nb214_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel214_ia32_sse2.nb214_dosingle
        jmp   _nb_kernel214_ia32_sse2.nb214_updateouterdata
_nb_kernel214_ia32_sse2.nb214_dosingle: 
        movl  nb214_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb214_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb214_jxO(%esp)
        movsd  %xmm1,nb214_jzO(%esp)
        movsd  %xmm2,nb214_jyH1(%esp)
        movsd  %xmm3,nb214_jxH2(%esp)
        movsd  %xmm4,nb214_jzH2(%esp)
        movsd  %xmm5,nb214_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb214_jyO(%esp)
        movsd  %xmm1,nb214_jxH1(%esp)
        movsd  %xmm2,nb214_jzH1(%esp)
        movsd  %xmm3,nb214_jyH2(%esp)
        movsd  %xmm4,nb214_jxM(%esp)
        movsd  %xmm5,nb214_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb214_ixO(%esp),%xmm0
        movapd nb214_iyO(%esp),%xmm1
        movapd nb214_izO(%esp),%xmm2
        movapd nb214_ixH1(%esp),%xmm3
        movapd nb214_iyH1(%esp),%xmm4
        movapd nb214_izH1(%esp),%xmm5
        subsd  nb214_jxO(%esp),%xmm0
        subsd  nb214_jyO(%esp),%xmm1
        subsd  nb214_jzO(%esp),%xmm2
        subsd  nb214_jxH1(%esp),%xmm3
        subsd  nb214_jyH1(%esp),%xmm4
        subsd  nb214_jzH1(%esp),%xmm5
        movapd %xmm0,nb214_dxOO(%esp)
        movapd %xmm1,nb214_dyOO(%esp)
        movapd %xmm2,nb214_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH1H1(%esp)
        movapd %xmm4,nb214_dyH1H1(%esp)
        movapd %xmm5,nb214_dzH1H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqOO(%esp)
        movapd %xmm3,nb214_rsqH1H1(%esp)

        movapd nb214_ixH1(%esp),%xmm0
        movapd nb214_iyH1(%esp),%xmm1
        movapd nb214_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214_jxH2(%esp),%xmm0
        subsd  nb214_jyH2(%esp),%xmm1
        subsd  nb214_jzH2(%esp),%xmm2
        subsd  nb214_jxM(%esp),%xmm3
        subsd  nb214_jyM(%esp),%xmm4
        subsd  nb214_jzM(%esp),%xmm5
        movapd %xmm0,nb214_dxH1H2(%esp)
        movapd %xmm1,nb214_dyH1H2(%esp)
        movapd %xmm2,nb214_dzH1H2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH1M(%esp)
        movapd %xmm4,nb214_dyH1M(%esp)
        movapd %xmm5,nb214_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqH1H2(%esp)
        movapd %xmm3,nb214_rsqH1M(%esp)

        movapd nb214_ixH2(%esp),%xmm0
        movapd nb214_iyH2(%esp),%xmm1
        movapd nb214_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214_jxH1(%esp),%xmm0
        subsd  nb214_jyH1(%esp),%xmm1
        subsd  nb214_jzH1(%esp),%xmm2
        subsd  nb214_jxH2(%esp),%xmm3
        subsd  nb214_jyH2(%esp),%xmm4
        subsd  nb214_jzH2(%esp),%xmm5
        movapd %xmm0,nb214_dxH2H1(%esp)
        movapd %xmm1,nb214_dyH2H1(%esp)
        movapd %xmm2,nb214_dzH2H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxH2H2(%esp)
        movapd %xmm4,nb214_dyH2H2(%esp)
        movapd %xmm5,nb214_dzH2H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb214_rsqH2H1(%esp)
        movapd %xmm3,nb214_rsqH2H2(%esp)

        movapd nb214_ixH2(%esp),%xmm0
        movapd nb214_iyH2(%esp),%xmm1
        movapd nb214_izH2(%esp),%xmm2
        movapd nb214_ixM(%esp),%xmm3
        movapd nb214_iyM(%esp),%xmm4
        movapd nb214_izM(%esp),%xmm5
        subsd  nb214_jxM(%esp),%xmm0
        subsd  nb214_jyM(%esp),%xmm1
        subsd  nb214_jzM(%esp),%xmm2
        subsd  nb214_jxH1(%esp),%xmm3
        subsd  nb214_jyH1(%esp),%xmm4
        subsd  nb214_jzH1(%esp),%xmm5
        movapd %xmm0,nb214_dxH2M(%esp)
        movapd %xmm1,nb214_dyH2M(%esp)
        movapd %xmm2,nb214_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxMH1(%esp)
        movapd %xmm4,nb214_dyMH1(%esp)
        movapd %xmm5,nb214_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb214_rsqH2M(%esp)
        movapd %xmm4,nb214_rsqMH1(%esp)

        movapd nb214_ixM(%esp),%xmm0
        movapd nb214_iyM(%esp),%xmm1
        movapd nb214_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214_jxH2(%esp),%xmm0
        subsd  nb214_jyH2(%esp),%xmm1
        subsd  nb214_jzH2(%esp),%xmm2
        subsd  nb214_jxM(%esp),%xmm3
        subsd  nb214_jyM(%esp),%xmm4
        subsd  nb214_jzM(%esp),%xmm5
        movapd %xmm0,nb214_dxMH2(%esp)
        movapd %xmm1,nb214_dyMH2(%esp)
        movapd %xmm2,nb214_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb214_dxMM(%esp)
        movapd %xmm4,nb214_dyMM(%esp)
        movapd %xmm5,nb214_dzMM(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb214_rsqMH2(%esp)
        movapd %xmm4,nb214_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb214_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm3   ## iter1 
        mulsd   nb214_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm1   ## rinv 
        mulsd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvMH2(%esp)
        movapd %xmm5,nb214_rinvMM(%esp)

        movapd nb214_rsqOO(%esp),%xmm0
        movapd nb214_rsqH1H1(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb214_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm1   ## rinv 
        mulsd   nb214_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb214_rinvsqOO(%esp)
        movapd %xmm5,nb214_rinvH1H1(%esp)

        movapd nb214_rsqH1H2(%esp),%xmm0
        movapd nb214_rsqH1M(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm3   ## iter1 
        mulsd   nb214_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm1   ## rinv 
        mulsd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvH1H2(%esp)
        movapd %xmm5,nb214_rinvH1M(%esp)

        movapd nb214_rsqH2H1(%esp),%xmm0
        movapd nb214_rsqH2H2(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm3   ## iter1a 
        mulsd   nb214_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm1   ## rinv 
        mulsd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvH2H1(%esp)
        movapd %xmm5,nb214_rinvH2H2(%esp)

        movapd nb214_rsqMH1(%esp),%xmm0
        movapd nb214_rsqH2M(%esp),%xmm4
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
        movapd  nb214_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm3   ## iter1a 
        mulsd   nb214_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214_half(%esp),%xmm1   ## rinv 
        mulsd   nb214_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214_rinvMH1(%esp)
        movapd %xmm5,nb214_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb214_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb214_c6(%esp),%xmm1
        mulsd  nb214_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb214_Vvdwtot(%esp),%xmm3
        mulsd  nb214_six(%esp),%xmm1
        mulsd  nb214_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movsd %xmm3,nb214_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1


        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb214_dxOO(%esp),%xmm0
        mulsd nb214_dyOO(%esp),%xmm1
        mulsd nb214_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixO(%esp),%xmm0
        addsd nb214_fiyO(%esp),%xmm1
        addsd nb214_fizO(%esp),%xmm2
        movsd %xmm3,nb214_fjxO(%esp)
        movsd %xmm4,nb214_fjyO(%esp)
        movsd %xmm5,nb214_fjzO(%esp)
        movsd %xmm0,nb214_fixO(%esp)
        movsd %xmm1,nb214_fiyO(%esp)
        movsd %xmm2,nb214_fizO(%esp)

        ## H1-H1 interaction 
        movsd nb214_rinvH1H1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb214_dxH1H1(%esp),%xmm0
        mulsd nb214_dyH1H1(%esp),%xmm1
        mulsd nb214_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH1(%esp),%xmm0
        addsd nb214_fiyH1(%esp),%xmm1
        addsd nb214_fizH1(%esp),%xmm2
        movsd %xmm3,nb214_fjxH1(%esp)
        movsd %xmm4,nb214_fjyH1(%esp)
        movsd %xmm5,nb214_fjzH1(%esp)
        movsd %xmm0,nb214_fixH1(%esp)
        movsd %xmm1,nb214_fiyH1(%esp)
        movsd %xmm2,nb214_fizH1(%esp)

        ## H1-H2 interaction  
        movsd nb214_rinvH1H2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb214_dxH1H2(%esp),%xmm0
        mulsd nb214_dyH1H2(%esp),%xmm1
        mulsd nb214_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH1(%esp),%xmm0
        addsd nb214_fiyH1(%esp),%xmm1
        addsd nb214_fizH1(%esp),%xmm2
        movsd %xmm3,nb214_fjxH2(%esp)
        movsd %xmm4,nb214_fjyH2(%esp)
        movsd %xmm5,nb214_fjzH2(%esp)
        movsd %xmm0,nb214_fixH1(%esp)
        movsd %xmm1,nb214_fiyH1(%esp)
        movsd %xmm2,nb214_fizH1(%esp)

        ## H1-M interaction 
        movsd nb214_rinvH1M(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb214_dxH1M(%esp),%xmm0
        mulsd nb214_dyH1M(%esp),%xmm1
        mulsd nb214_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH1(%esp),%xmm0
        addsd nb214_fiyH1(%esp),%xmm1
        addsd nb214_fizH1(%esp),%xmm2
        movsd %xmm3,nb214_fjxM(%esp)
        movsd %xmm4,nb214_fjyM(%esp)
        movsd %xmm5,nb214_fjzM(%esp)
        movsd %xmm0,nb214_fixH1(%esp)
        movsd %xmm1,nb214_fiyH1(%esp)
        movsd %xmm2,nb214_fizH1(%esp)

        ## H2-H1 interaction 
        movsd nb214_rinvH2H1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH1(%esp),%xmm3
        movapd nb214_fjyH1(%esp),%xmm4
        movapd nb214_fjzH1(%esp),%xmm5
        mulsd nb214_dxH2H1(%esp),%xmm0
        mulsd nb214_dyH2H1(%esp),%xmm1
        mulsd nb214_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH2(%esp),%xmm0
        addsd nb214_fiyH2(%esp),%xmm1
        addsd nb214_fizH2(%esp),%xmm2
        movsd %xmm3,nb214_fjxH1(%esp)
        movsd %xmm4,nb214_fjyH1(%esp)
        movsd %xmm5,nb214_fjzH1(%esp)
        movsd %xmm0,nb214_fixH2(%esp)
        movsd %xmm1,nb214_fiyH2(%esp)
        movsd %xmm2,nb214_fizH2(%esp)

        ## H2-H2 interaction 
        movsd nb214_rinvH2H2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movsd nb214_fjxH2(%esp),%xmm3
        movsd nb214_fjyH2(%esp),%xmm4
        movsd nb214_fjzH2(%esp),%xmm5
        mulsd nb214_dxH2H2(%esp),%xmm0
        mulsd nb214_dyH2H2(%esp),%xmm1
        mulsd nb214_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH2(%esp),%xmm0
        addsd nb214_fiyH2(%esp),%xmm1
        addsd nb214_fizH2(%esp),%xmm2
        movsd %xmm3,nb214_fjxH2(%esp)
        movsd %xmm4,nb214_fjyH2(%esp)
        movsd %xmm5,nb214_fjzH2(%esp)
        movsd %xmm0,nb214_fixH2(%esp)
        movsd %xmm1,nb214_fiyH2(%esp)
        movsd %xmm2,nb214_fizH2(%esp)

        ## H2-M interaction 
        movsd nb214_rinvH2M(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxM(%esp),%xmm3
        movapd nb214_fjyM(%esp),%xmm4
        movapd nb214_fjzM(%esp),%xmm5
        mulsd nb214_dxH2M(%esp),%xmm0
        mulsd nb214_dyH2M(%esp),%xmm1
        mulsd nb214_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixH2(%esp),%xmm0
        addsd nb214_fiyH2(%esp),%xmm1
        addsd nb214_fizH2(%esp),%xmm2
        movsd %xmm3,nb214_fjxM(%esp)
        movsd %xmm4,nb214_fjyM(%esp)
        movsd %xmm5,nb214_fjzM(%esp)
        movsd %xmm0,nb214_fixH2(%esp)
        movsd %xmm1,nb214_fiyH2(%esp)
        movsd %xmm2,nb214_fizH2(%esp)

        ## M-H1 interaction 
        movsd nb214_rinvMH1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH1(%esp),%xmm3
        movapd nb214_fjyH1(%esp),%xmm4
        movapd nb214_fjzH1(%esp),%xmm5
        mulsd nb214_dxMH1(%esp),%xmm0
        mulsd nb214_dyMH1(%esp),%xmm1
        mulsd nb214_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixM(%esp),%xmm0
        addsd nb214_fiyM(%esp),%xmm1
        addsd nb214_fizM(%esp),%xmm2
        movsd %xmm3,nb214_fjxH1(%esp)
        movsd %xmm4,nb214_fjyH1(%esp)
        movsd %xmm5,nb214_fjzH1(%esp)
        movsd %xmm0,nb214_fixM(%esp)
        movsd %xmm1,nb214_fiyM(%esp)
        movsd %xmm2,nb214_fizM(%esp)

        ## M-H2 interaction 
        movsd nb214_rinvMH2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxH2(%esp),%xmm3
        movapd nb214_fjyH2(%esp),%xmm4
        movapd nb214_fjzH2(%esp),%xmm5
        mulsd nb214_dxMH2(%esp),%xmm0
        mulsd nb214_dyMH2(%esp),%xmm1
        mulsd nb214_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixM(%esp),%xmm0
        addsd nb214_fiyM(%esp),%xmm1
        addsd nb214_fizM(%esp),%xmm2
        movsd %xmm3,nb214_fjxH2(%esp)
        movsd %xmm4,nb214_fjyH2(%esp)
        movsd %xmm5,nb214_fjzH2(%esp)
        movsd %xmm0,nb214_fixM(%esp)
        movsd %xmm1,nb214_fiyM(%esp)
        movsd %xmm2,nb214_fizM(%esp)

        ## M-M interaction 
        movsd nb214_rinvMM(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        movsd nb214_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb214_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214_crf(%esp),%xmm6
        mulsd  nb214_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb214_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb214_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb214_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0
        movsd %xmm6,nb214_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb214_fjxM(%esp),%xmm3
        movapd nb214_fjyM(%esp),%xmm4
        movapd nb214_fjzM(%esp),%xmm5
        mulsd nb214_dxMM(%esp),%xmm0
        mulsd nb214_dyMM(%esp),%xmm1
        mulsd nb214_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb214_fixM(%esp),%xmm0
        addsd nb214_fiyM(%esp),%xmm1
        addsd nb214_fizM(%esp),%xmm2
        movsd %xmm3,nb214_fjxM(%esp)
        movsd %xmm4,nb214_fjyM(%esp)
        movsd %xmm5,nb214_fjzM(%esp)
        movsd %xmm0,nb214_fixM(%esp)
        movsd %xmm1,nb214_fiyM(%esp)
        movsd %xmm2,nb214_fizM(%esp)

        movl nb214_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - merge forces
        movlpd nb214_fjxO(%esp),%xmm0
        movlpd nb214_fjzO(%esp),%xmm1
        movlpd nb214_fjyH1(%esp),%xmm2
        movlpd nb214_fjxH2(%esp),%xmm3
        movlpd nb214_fjzH2(%esp),%xmm4
        movlpd nb214_fjyM(%esp),%xmm5

        movhpd nb214_fjyO(%esp),%xmm0
        movhpd nb214_fjxH1(%esp),%xmm1
        movhpd nb214_fjzH1(%esp),%xmm2
        movhpd nb214_fjyH2(%esp),%xmm3
        movhpd nb214_fjxM(%esp),%xmm4
        movhpd nb214_fjzM(%esp),%xmm5

        movlps (%edi,%eax,8),%xmm6
        movhps 8(%edi,%eax,8),%xmm6
        movlps 16(%edi,%eax,8),%xmm7
        movhps 24(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm0
        addpd  %xmm7,%xmm1
        movlps 32(%edi,%eax,8),%xmm6
        movhps 40(%edi,%eax,8),%xmm6
        movlps 48(%edi,%eax,8),%xmm7
        movhps 56(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm2
        addpd  %xmm7,%xmm3
        movlps 64(%edi,%eax,8),%xmm6
        movhps 72(%edi,%eax,8),%xmm6
        movlps 80(%edi,%eax,8),%xmm7
        movhps 88(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm4
        addpd  %xmm7,%xmm5

        movlpd %xmm0,(%edi,%eax,8)
        movhpd %xmm0,8(%edi,%eax,8)
        movlpd %xmm1,16(%edi,%eax,8)
        movhpd %xmm1,24(%edi,%eax,8)
        movlpd %xmm2,32(%edi,%eax,8)
        movhpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movhpd %xmm3,56(%edi,%eax,8)
        movlpd %xmm4,64(%edi,%eax,8)
        movhpd %xmm4,72(%edi,%eax,8)
        movlpd %xmm5,80(%edi,%eax,8)
        movhpd %xmm5,88(%edi,%eax,8)

_nb_kernel214_ia32_sse2.nb214_updateouterdata: 
        movl  nb214_ii3(%esp),%ecx
        movl  nb214_faction(%ebp),%edi
        movl  nb214_fshift(%ebp),%esi
        movl  nb214_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb214_fixO(%esp),%xmm0
        movapd nb214_fiyO(%esp),%xmm1
        movapd nb214_fizO(%esp),%xmm2

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
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb214_fixH1(%esp),%xmm0
        movapd nb214_fiyH1(%esp),%xmm1
        movapd nb214_fizH1(%esp),%xmm2

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
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb214_fixH2(%esp),%xmm0
        movapd nb214_fiyH2(%esp),%xmm1
        movapd nb214_fizH2(%esp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb214_fixM(%esp),%xmm0
        movapd nb214_fiyM(%esp),%xmm1
        movapd nb214_fizM(%esp),%xmm2

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
        movl nb214_n(%esp),%esi
        ## get group index for i particle 
        movl  nb214_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb214_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb214_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb214_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb214_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb214_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel214_ia32_sse2.nb214_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb214_n(%esp)
        jmp _nb_kernel214_ia32_sse2.nb214_outer
_nb_kernel214_ia32_sse2.nb214_outerend: 
        ## check if more outer neighborlists remain
        movl  nb214_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel214_ia32_sse2.nb214_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel214_ia32_sse2.nb214_threadloop
_nb_kernel214_ia32_sse2.nb214_end: 
        emms

        movl nb214_nouter(%esp),%eax
        movl nb214_ninner(%esp),%ebx
        movl nb214_outeriter(%ebp),%ecx
        movl nb214_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb214_salign(%esp),%eax
        addl %eax,%esp
        addl $1848,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel214nf_ia32_sse2
.globl _nb_kernel214nf_ia32_sse2
nb_kernel214nf_ia32_sse2:       
_nb_kernel214nf_ia32_sse2:      
.set nb214nf_p_nri, 8
.set nb214nf_iinr, 12
.set nb214nf_jindex, 16
.set nb214nf_jjnr, 20
.set nb214nf_shift, 24
.set nb214nf_shiftvec, 28
.set nb214nf_fshift, 32
.set nb214nf_gid, 36
.set nb214nf_pos, 40
.set nb214nf_faction, 44
.set nb214nf_charge, 48
.set nb214nf_p_facel, 52
.set nb214nf_argkrf, 56
.set nb214nf_argcrf, 60
.set nb214nf_Vc, 64
.set nb214nf_type, 68
.set nb214nf_p_ntype, 72
.set nb214nf_vdwparam, 76
.set nb214nf_Vvdw, 80
.set nb214nf_p_tabscale, 84
.set nb214nf_VFtab, 88
.set nb214nf_invsqrta, 92
.set nb214nf_dvda, 96
.set nb214nf_p_gbtabscale, 100
.set nb214nf_GBtab, 104
.set nb214nf_p_nthreads, 108
.set nb214nf_count, 112
.set nb214nf_mtx, 116
.set nb214nf_outeriter, 120
.set nb214nf_inneriter, 124
.set nb214nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb214nf_ixO, 0
.set nb214nf_iyO, 16
.set nb214nf_izO, 32
.set nb214nf_ixH1, 48
.set nb214nf_iyH1, 64
.set nb214nf_izH1, 80
.set nb214nf_ixH2, 96
.set nb214nf_iyH2, 112
.set nb214nf_izH2, 128
.set nb214nf_ixM, 144
.set nb214nf_iyM, 160
.set nb214nf_izM, 176
.set nb214nf_jxO, 192
.set nb214nf_jyO, 208
.set nb214nf_jzO, 224
.set nb214nf_jxH1, 240
.set nb214nf_jyH1, 256
.set nb214nf_jzH1, 272
.set nb214nf_jxH2, 288
.set nb214nf_jyH2, 304
.set nb214nf_jzH2, 320
.set nb214nf_jxM, 336
.set nb214nf_jyM, 352
.set nb214nf_jzM, 368
.set nb214nf_qqMM, 384
.set nb214nf_qqMH, 400
.set nb214nf_qqHH, 416
.set nb214nf_two, 432
.set nb214nf_c6, 448
.set nb214nf_c12, 464
.set nb214nf_vctot, 480
.set nb214nf_Vvdwtot, 496
.set nb214nf_half, 512
.set nb214nf_three, 528
.set nb214nf_rsqOO, 544
.set nb214nf_rsqH1H1, 560
.set nb214nf_rsqH1H2, 576
.set nb214nf_rsqH1M, 592
.set nb214nf_rsqH2H1, 608
.set nb214nf_rsqH2H2, 624
.set nb214nf_rsqH2M, 640
.set nb214nf_rsqMH1, 656
.set nb214nf_rsqMH2, 672
.set nb214nf_rsqMM, 688
.set nb214nf_rinvsqOO, 704
.set nb214nf_rinvH1H1, 720
.set nb214nf_rinvH1H2, 736
.set nb214nf_rinvH1M, 752
.set nb214nf_rinvH2H1, 768
.set nb214nf_rinvH2H2, 784
.set nb214nf_rinvH2M, 800
.set nb214nf_rinvMH1, 816
.set nb214nf_rinvMH2, 832
.set nb214nf_rinvMM, 848
.set nb214nf_krf, 864
.set nb214nf_crf, 880
.set nb214nf_is3, 896
.set nb214nf_ii3, 900
.set nb214nf_innerjjnr, 904
.set nb214nf_innerk, 908
.set nb214nf_n, 912
.set nb214nf_nn1, 916
.set nb214nf_nri, 920
.set nb214nf_nouter, 924
.set nb214nf_ninner, 928
.set nb214nf_salign, 932
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $936,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb214nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb214nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb214nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb214nf_nouter(%esp)
        movl %eax,nb214nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb214nf_half(%esp)
        movl %ebx,nb214nf_half+4(%esp)
        movsd nb214nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb214nf_half(%esp)
        movapd %xmm2,nb214nf_two(%esp)
        movapd %xmm3,nb214nf_three(%esp)

        movl nb214nf_argkrf(%ebp),%esi
        movl nb214nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb214nf_krf(%esp)
        movapd %xmm6,nb214nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb214nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb214nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb214nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb214nf_qqMM(%esp)
        movapd %xmm4,nb214nf_qqMH(%esp)
        movapd %xmm5,nb214nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb214nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb214nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb214nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb214nf_c6(%esp)
        movapd %xmm1,nb214nf_c12(%esp)

_nb_kernel214nf_ia32_sse2.nb214nf_threadloop: 
        movl  nb214nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel214nf_ia32_sse2.nb214nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel214nf_ia32_sse2.nb214nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb214nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb214nf_n(%esp)
        movl %ebx,nb214nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel214nf_ia32_sse2.nb214nf_outerstart
        jmp _nb_kernel214nf_ia32_sse2.nb214nf_end

_nb_kernel214nf_ia32_sse2.nb214nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb214nf_nouter(%esp),%ebx
        movl %ebx,nb214nf_nouter(%esp)

_nb_kernel214nf_ia32_sse2.nb214nf_outer: 
        movl  nb214nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb214nf_is3(%esp)            ## store is3 

        movl  nb214nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb214nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb214nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb214nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb214nf_ixO(%esp)
        movapd %xmm4,nb214nf_iyO(%esp)
        movapd %xmm5,nb214nf_izO(%esp)
        movapd %xmm6,nb214nf_ixH1(%esp)
        movapd %xmm7,nb214nf_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb214nf_izH1(%esp)
        movapd %xmm0,nb214nf_ixH2(%esp)
        movapd %xmm1,nb214nf_iyH2(%esp)
        movapd %xmm2,nb214nf_izH2(%esp)
        movapd %xmm3,nb214nf_ixM(%esp)
        movapd %xmm4,nb214nf_iyM(%esp)
        movapd %xmm5,nb214nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb214nf_vctot(%esp)
        movapd %xmm4,nb214nf_Vvdwtot(%esp)

        movl  nb214nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb214nf_pos(%ebp),%esi
        movl  nb214nf_faction(%ebp),%edi
        movl  nb214nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb214nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb214nf_ninner(%esp),%ecx
        movl  %ecx,nb214nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb214nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel214nf_ia32_sse2.nb214nf_unroll_loop
        jmp   _nb_kernel214nf_ia32_sse2.nb214nf_checksingle
_nb_kernel214nf_ia32_sse2.nb214nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb214nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb214nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb214nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb214nf_jxO(%esp)
        movapd  %xmm1,nb214nf_jyO(%esp)
        movapd  %xmm3,nb214nf_jzO(%esp)
        movapd  %xmm4,nb214nf_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb214nf_jyH1(%esp)
        movapd  %xmm1,nb214nf_jzH1(%esp)
        movapd  %xmm3,nb214nf_jxH2(%esp)
        movapd  %xmm4,nb214nf_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb214nf_jzH2(%esp)
        movapd  %xmm1,nb214nf_jxM(%esp)
        movapd  %xmm3,nb214nf_jyM(%esp)
        movapd  %xmm4,nb214nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb214nf_ixO(%esp),%xmm0
        movapd nb214nf_iyO(%esp),%xmm1
        movapd nb214nf_izO(%esp),%xmm2
        movapd nb214nf_ixH1(%esp),%xmm3
        movapd nb214nf_iyH1(%esp),%xmm4
        movapd nb214nf_izH1(%esp),%xmm5
        subpd  nb214nf_jxO(%esp),%xmm0
        subpd  nb214nf_jyO(%esp),%xmm1
        subpd  nb214nf_jzO(%esp),%xmm2
        subpd  nb214nf_jxH1(%esp),%xmm3
        subpd  nb214nf_jyH1(%esp),%xmm4
        subpd  nb214nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqOO(%esp)
        movapd %xmm3,nb214nf_rsqH1H1(%esp)

        movapd nb214nf_ixH1(%esp),%xmm0
        movapd nb214nf_iyH1(%esp),%xmm1
        movapd nb214nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214nf_jxH2(%esp),%xmm0
        subpd  nb214nf_jyH2(%esp),%xmm1
        subpd  nb214nf_jzH2(%esp),%xmm2
        subpd  nb214nf_jxM(%esp),%xmm3
        subpd  nb214nf_jyM(%esp),%xmm4
        subpd  nb214nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH1H2(%esp)
        movapd %xmm3,nb214nf_rsqH1M(%esp)

        movapd nb214nf_ixH2(%esp),%xmm0
        movapd nb214nf_iyH2(%esp),%xmm1
        movapd nb214nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214nf_jxH1(%esp),%xmm0
        subpd  nb214nf_jyH1(%esp),%xmm1
        subpd  nb214nf_jzH1(%esp),%xmm2
        subpd  nb214nf_jxH2(%esp),%xmm3
        subpd  nb214nf_jyH2(%esp),%xmm4
        subpd  nb214nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH2H1(%esp)
        movapd %xmm3,nb214nf_rsqH2H2(%esp)

        movapd nb214nf_ixH2(%esp),%xmm0
        movapd nb214nf_iyH2(%esp),%xmm1
        movapd nb214nf_izH2(%esp),%xmm2
        movapd nb214nf_ixM(%esp),%xmm3
        movapd nb214nf_iyM(%esp),%xmm4
        movapd nb214nf_izM(%esp),%xmm5
        subpd  nb214nf_jxM(%esp),%xmm0
        subpd  nb214nf_jyM(%esp),%xmm1
        subpd  nb214nf_jzM(%esp),%xmm2
        subpd  nb214nf_jxH1(%esp),%xmm3
        subpd  nb214nf_jyH1(%esp),%xmm4
        subpd  nb214nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH2M(%esp)
        movapd %xmm4,nb214nf_rsqMH1(%esp)

        movapd nb214nf_ixM(%esp),%xmm0
        movapd nb214nf_iyM(%esp),%xmm1
        movapd nb214nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb214nf_jxH2(%esp),%xmm0
        subpd  nb214nf_jyH2(%esp),%xmm1
        subpd  nb214nf_jzH2(%esp),%xmm2
        subpd  nb214nf_jxM(%esp),%xmm3
        subpd  nb214nf_jyM(%esp),%xmm4
        subpd  nb214nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqMH2(%esp)
        movapd %xmm4,nb214nf_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb214nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvMH2(%esp)
        movapd %xmm5,nb214nf_rinvMM(%esp)

        movapd nb214nf_rsqOO(%esp),%xmm0
        movapd nb214nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb214nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb214nf_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb214nf_rinvsqOO(%esp)
        movapd %xmm5,nb214nf_rinvH1H1(%esp)

        movapd nb214nf_rsqH1H2(%esp),%xmm0
        movapd nb214nf_rsqH1M(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb214nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvH1H2(%esp)
        movapd %xmm5,nb214nf_rinvH1M(%esp)

        movapd nb214nf_rsqH2H1(%esp),%xmm0
        movapd nb214nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb214nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvH2H1(%esp)
        movapd %xmm5,nb214nf_rinvH2H2(%esp)

        movapd nb214nf_rsqMH1(%esp),%xmm0
        movapd nb214nf_rsqH2M(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb214nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvMH1(%esp)
        movapd %xmm5,nb214nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb214nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb214nf_c6(%esp),%xmm1
        mulpd  nb214nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb214nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb214nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movapd nb214nf_rinvH1H1(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        addpd  %xmm0,%xmm5      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm5
        mulpd  nb214nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm5   ## local vctot summation variable 
        movapd %xmm5,nb214nf_vctot(%esp)

        ## H1-H2 interaction  
        movapd nb214nf_rinvH1H2(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        addpd  %xmm0,%xmm5      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm5
        mulpd  nb214nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm5   ## local vctot summation variable 
        movapd %xmm5,nb214nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb214nf_rinvH1M(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb214nf_rinvH2H1(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb214nf_rinvH2H2(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb214nf_rinvH2M(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb214nf_rinvMH1(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## M-H2 interaction 
        movapd nb214nf_rinvMH2(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## M-M interaction 
        movapd nb214nf_rinvMM(%esp),%xmm0
        movapd nb214nf_krf(%esp),%xmm5

        mulpd  nb214nf_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb214nf_crf(%esp),%xmm6
        mulpd  nb214nf_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addpd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movapd %xmm6,nb214nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb214nf_innerk(%esp)
        jl    _nb_kernel214nf_ia32_sse2.nb214nf_checksingle
        jmp   _nb_kernel214nf_ia32_sse2.nb214nf_unroll_loop
_nb_kernel214nf_ia32_sse2.nb214nf_checksingle: 
        movl  nb214nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel214nf_ia32_sse2.nb214nf_dosingle
        jmp   _nb_kernel214nf_ia32_sse2.nb214nf_updateouterdata
_nb_kernel214nf_ia32_sse2.nb214nf_dosingle: 
        movl  nb214nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb214nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb214nf_jxO(%esp)
        movsd  %xmm1,nb214nf_jzO(%esp)
        movsd  %xmm2,nb214nf_jyH1(%esp)
        movsd  %xmm3,nb214nf_jxH2(%esp)
        movsd  %xmm4,nb214nf_jzH2(%esp)
        movsd  %xmm5,nb214nf_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb214nf_jyO(%esp)
        movsd  %xmm1,nb214nf_jxH1(%esp)
        movsd  %xmm2,nb214nf_jzH1(%esp)
        movsd  %xmm3,nb214nf_jyH2(%esp)
        movsd  %xmm4,nb214nf_jxM(%esp)
        movsd  %xmm5,nb214nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb214nf_ixO(%esp),%xmm0
        movapd nb214nf_iyO(%esp),%xmm1
        movapd nb214nf_izO(%esp),%xmm2
        movapd nb214nf_ixH1(%esp),%xmm3
        movapd nb214nf_iyH1(%esp),%xmm4
        movapd nb214nf_izH1(%esp),%xmm5
        subsd  nb214nf_jxO(%esp),%xmm0
        subsd  nb214nf_jyO(%esp),%xmm1
        subsd  nb214nf_jzO(%esp),%xmm2
        subsd  nb214nf_jxH1(%esp),%xmm3
        subsd  nb214nf_jyH1(%esp),%xmm4
        subsd  nb214nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqOO(%esp)
        movapd %xmm3,nb214nf_rsqH1H1(%esp)

        movapd nb214nf_ixH1(%esp),%xmm0
        movapd nb214nf_iyH1(%esp),%xmm1
        movapd nb214nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214nf_jxH2(%esp),%xmm0
        subsd  nb214nf_jyH2(%esp),%xmm1
        subsd  nb214nf_jzH2(%esp),%xmm2
        subsd  nb214nf_jxM(%esp),%xmm3
        subsd  nb214nf_jyM(%esp),%xmm4
        subsd  nb214nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH1H2(%esp)
        movapd %xmm3,nb214nf_rsqH1M(%esp)

        movapd nb214nf_ixH2(%esp),%xmm0
        movapd nb214nf_iyH2(%esp),%xmm1
        movapd nb214nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214nf_jxH1(%esp),%xmm0
        subsd  nb214nf_jyH1(%esp),%xmm1
        subsd  nb214nf_jzH1(%esp),%xmm2
        subsd  nb214nf_jxH2(%esp),%xmm3
        subsd  nb214nf_jyH2(%esp),%xmm4
        subsd  nb214nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH2H1(%esp)
        movapd %xmm3,nb214nf_rsqH2H2(%esp)

        movapd nb214nf_ixH2(%esp),%xmm0
        movapd nb214nf_iyH2(%esp),%xmm1
        movapd nb214nf_izH2(%esp),%xmm2
        movapd nb214nf_ixM(%esp),%xmm3
        movapd nb214nf_iyM(%esp),%xmm4
        movapd nb214nf_izM(%esp),%xmm5
        subsd  nb214nf_jxM(%esp),%xmm0
        subsd  nb214nf_jyM(%esp),%xmm1
        subsd  nb214nf_jzM(%esp),%xmm2
        subsd  nb214nf_jxH1(%esp),%xmm3
        subsd  nb214nf_jyH1(%esp),%xmm4
        subsd  nb214nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqH2M(%esp)
        movapd %xmm4,nb214nf_rsqMH1(%esp)

        movapd nb214nf_ixM(%esp),%xmm0
        movapd nb214nf_iyM(%esp),%xmm1
        movapd nb214nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb214nf_jxH2(%esp),%xmm0
        subsd  nb214nf_jyH2(%esp),%xmm1
        subsd  nb214nf_jzH2(%esp),%xmm2
        subsd  nb214nf_jxM(%esp),%xmm3
        subsd  nb214nf_jyM(%esp),%xmm4
        subsd  nb214nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb214nf_rsqMH2(%esp)
        movapd %xmm4,nb214nf_rsqMM(%esp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb214nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvMH2(%esp)
        movapd %xmm5,nb214nf_rinvMM(%esp)

        movapd nb214nf_rsqOO(%esp),%xmm0
        movapd nb214nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb214nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb214nf_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb214nf_rinvsqOO(%esp)
        movapd %xmm5,nb214nf_rinvH1H1(%esp)

        movapd nb214nf_rsqH1H2(%esp),%xmm0
        movapd nb214nf_rsqH1M(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb214nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvH1H2(%esp)
        movapd %xmm5,nb214nf_rinvH1M(%esp)

        movapd nb214nf_rsqH2H1(%esp),%xmm0
        movapd nb214nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb214nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvH2H1(%esp)
        movapd %xmm5,nb214nf_rinvH2H2(%esp)

        movapd nb214nf_rsqMH1(%esp),%xmm0
        movapd nb214nf_rsqH2M(%esp),%xmm4
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
        movapd  nb214nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb214nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb214nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb214nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb214nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb214nf_rinvMH1(%esp)
        movapd %xmm5,nb214nf_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb214nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb214nf_c6(%esp),%xmm1
        mulsd  nb214nf_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb214nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb214nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movsd nb214nf_rinvH1H1(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## H1-H2 interaction  
        movsd nb214nf_rinvH1H2(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## H1-M interaction 
        movsd nb214nf_rinvH1M(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## H2-H1 interaction 
        movsd nb214nf_rinvH2H1(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## H2-H2 interaction 
        movsd nb214nf_rinvH2H2(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## H2-M interaction 
        movsd nb214nf_rinvH2M(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## M-H1 interaction 
        movsd nb214nf_rinvMH1(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## M-H2 interaction 
        movsd nb214nf_rinvMH2(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqMH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

        ## M-M interaction 
        movsd nb214nf_rinvMM(%esp),%xmm0
        movsd nb214nf_krf(%esp),%xmm5

        mulsd  nb214nf_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movsd %xmm5,%xmm6
        addsd  %xmm0,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb214nf_crf(%esp),%xmm6
        mulsd  nb214nf_qqMM(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addsd  nb214nf_vctot(%esp),%xmm6   ## local vctot summation variable 
        movsd %xmm6,nb214nf_vctot(%esp)

_nb_kernel214nf_ia32_sse2.nb214nf_updateouterdata: 
        ## get n from stack
        movl nb214nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb214nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb214nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb214nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb214nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb214nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb214nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel214nf_ia32_sse2.nb214nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb214nf_n(%esp)
        jmp _nb_kernel214nf_ia32_sse2.nb214nf_outer
_nb_kernel214nf_ia32_sse2.nb214nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb214nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel214nf_ia32_sse2.nb214nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel214nf_ia32_sse2.nb214nf_threadloop
_nb_kernel214nf_ia32_sse2.nb214nf_end: 
        emms

        movl nb214nf_nouter(%esp),%eax
        movl nb214nf_ninner(%esp),%ebx
        movl nb214nf_outeriter(%ebp),%ecx
        movl nb214nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb214nf_salign(%esp),%eax
        addl %eax,%esp
        addl $936,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret

