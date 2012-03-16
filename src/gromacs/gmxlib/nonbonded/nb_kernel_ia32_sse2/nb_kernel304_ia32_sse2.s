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



.globl nb_kernel304_ia32_sse2
.globl _nb_kernel304_ia32_sse2
nb_kernel304_ia32_sse2: 
_nb_kernel304_ia32_sse2:        
.set nb304_p_nri, 8
.set nb304_iinr, 12
.set nb304_jindex, 16
.set nb304_jjnr, 20
.set nb304_shift, 24
.set nb304_shiftvec, 28
.set nb304_fshift, 32
.set nb304_gid, 36
.set nb304_pos, 40
.set nb304_faction, 44
.set nb304_charge, 48
.set nb304_p_facel, 52
.set nb304_argkrf, 56
.set nb304_argcrf, 60
.set nb304_Vc, 64
.set nb304_type, 68
.set nb304_p_ntype, 72
.set nb304_vdwparam, 76
.set nb304_Vvdw, 80
.set nb304_p_tabscale, 84
.set nb304_VFtab, 88
.set nb304_invsqrta, 92
.set nb304_dvda, 96
.set nb304_p_gbtabscale, 100
.set nb304_GBtab, 104
.set nb304_p_nthreads, 108
.set nb304_count, 112
.set nb304_mtx, 116
.set nb304_outeriter, 120
.set nb304_inneriter, 124
.set nb304_work, 128
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
.set nb304_fjxM, 960
.set nb304_fjyM, 976
.set nb304_fjzM, 992
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
.set nb304_innerjjnr, 1432
.set nb304_innerk, 1436
.set nb304_n, 1440
.set nb304_nn1, 1444
.set nb304_nri, 1448
.set nb304_nouter, 1452
.set nb304_ninner, 1456
.set nb304_salign, 1460
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1464,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb304_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb304_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb304_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb304_nouter(%esp)
        movl %eax,nb304_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb304_half(%esp)
        movl %ebx,nb304_half+4(%esp)
        movsd nb304_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb304_half(%esp)
        movapd %xmm2,nb304_two(%esp)
        movapd %xmm3,nb304_three(%esp)

        movl nb304_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb304_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb304_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb304_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb304_p_facel(%ebp),%esi
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
        movapd %xmm3,nb304_qqMM(%esp)
        movapd %xmm4,nb304_qqMH(%esp)
        movapd %xmm5,nb304_qqHH(%esp)

_nb_kernel304_ia32_sse2.nb304_threadloop: 
        movl  nb304_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel304_ia32_sse2.nb304_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304_ia32_sse2.nb304_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304_n(%esp)
        movl %ebx,nb304_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304_ia32_sse2.nb304_outerstart
        jmp _nb_kernel304_ia32_sse2.nb304_end

_nb_kernel304_ia32_sse2.nb304_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304_nouter(%esp),%ebx
        movl %ebx,nb304_nouter(%esp)

_nb_kernel304_ia32_sse2.nb304_outer: 
        movl  nb304_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb304_is3(%esp)      ## store is3 

        movl  nb304_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb304_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb304_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb304_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304_ixH1(%esp)
        movapd %xmm4,nb304_iyH1(%esp)
        movapd %xmm5,nb304_izH1(%esp)

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
        movapd %xmm0,nb304_ixH2(%esp)
        movapd %xmm1,nb304_iyH2(%esp)
        movapd %xmm2,nb304_izH2(%esp)
        movapd %xmm3,nb304_ixM(%esp)
        movapd %xmm4,nb304_iyM(%esp)
        movapd %xmm5,nb304_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb304_vctot(%esp)
        movapd %xmm4,nb304_fixM(%esp)
        movapd %xmm4,nb304_fiyM(%esp)
        movapd %xmm4,nb304_fizM(%esp)
        movapd %xmm4,nb304_fixH1(%esp)
        movapd %xmm4,nb304_fiyH1(%esp)
        movapd %xmm4,nb304_fizH1(%esp)
        movapd %xmm4,nb304_fixH2(%esp)
        movapd %xmm4,nb304_fiyH2(%esp)
        movapd %xmm4,nb304_fizH2(%esp)

        movl  nb304_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb304_pos(%ebp),%esi
        movl  nb304_faction(%ebp),%edi
        movl  nb304_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb304_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb304_ninner(%esp),%ecx
        movl  %ecx,nb304_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb304_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel304_ia32_sse2.nb304_unroll_loop
        jmp   _nb_kernel304_ia32_sse2.nb304_checksingle
_nb_kernel304_ia32_sse2.nb304_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb304_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb304_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb304_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb304_jxH1(%esp)
        movapd  %xmm3,nb304_jyH1(%esp)
        movapd  %xmm4,nb304_jzH1(%esp)
        movapd  %xmm5,nb304_jxH2(%esp)
        movapd  %xmm6,nb304_jyH2(%esp)
        movapd  %xmm7,nb304_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb304_jxM(%esp)
        movapd  %xmm3,nb304_jyM(%esp)
        movapd  %xmm4,nb304_jzM(%esp)

        movapd nb304_ixM(%esp),%xmm0
        movapd nb304_iyM(%esp),%xmm1
        movapd nb304_izM(%esp),%xmm2
        movapd nb304_ixM(%esp),%xmm3
        movapd nb304_iyM(%esp),%xmm4
        movapd nb304_izM(%esp),%xmm5
        subpd  nb304_jxM(%esp),%xmm0
        subpd  nb304_jyM(%esp),%xmm1
        subpd  nb304_jzM(%esp),%xmm2
        subpd  nb304_jxH1(%esp),%xmm3
        subpd  nb304_jyH1(%esp),%xmm4
        subpd  nb304_jzH1(%esp),%xmm5
        movapd %xmm0,nb304_dxMM(%esp)
        movapd %xmm1,nb304_dyMM(%esp)
        movapd %xmm2,nb304_dzMM(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxMH1(%esp)
        movapd %xmm4,nb304_dyMH1(%esp)
        movapd %xmm5,nb304_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqMM(%esp)
        movapd %xmm3,nb304_rsqMH1(%esp)

        movapd nb304_ixM(%esp),%xmm0
        movapd nb304_iyM(%esp),%xmm1
        movapd nb304_izM(%esp),%xmm2
        movapd nb304_ixH1(%esp),%xmm3
        movapd nb304_iyH1(%esp),%xmm4
        movapd nb304_izH1(%esp),%xmm5
        subpd  nb304_jxH2(%esp),%xmm0
        subpd  nb304_jyH2(%esp),%xmm1
        subpd  nb304_jzH2(%esp),%xmm2
        subpd  nb304_jxM(%esp),%xmm3
        subpd  nb304_jyM(%esp),%xmm4
        subpd  nb304_jzM(%esp),%xmm5
        movapd %xmm0,nb304_dxMH2(%esp)
        movapd %xmm1,nb304_dyMH2(%esp)
        movapd %xmm2,nb304_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH1M(%esp)
        movapd %xmm4,nb304_dyH1M(%esp)
        movapd %xmm5,nb304_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqMH2(%esp)
        movapd %xmm3,nb304_rsqH1M(%esp)

        movapd nb304_ixH1(%esp),%xmm0
        movapd nb304_iyH1(%esp),%xmm1
        movapd nb304_izH1(%esp),%xmm2
        movapd nb304_ixH1(%esp),%xmm3
        movapd nb304_iyH1(%esp),%xmm4
        movapd nb304_izH1(%esp),%xmm5
        subpd  nb304_jxH1(%esp),%xmm0
        subpd  nb304_jyH1(%esp),%xmm1
        subpd  nb304_jzH1(%esp),%xmm2
        subpd  nb304_jxH2(%esp),%xmm3
        subpd  nb304_jyH2(%esp),%xmm4
        subpd  nb304_jzH2(%esp),%xmm5
        movapd %xmm0,nb304_dxH1H1(%esp)
        movapd %xmm1,nb304_dyH1H1(%esp)
        movapd %xmm2,nb304_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH1H2(%esp)
        movapd %xmm4,nb304_dyH1H2(%esp)
        movapd %xmm5,nb304_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqH1H1(%esp)
        movapd %xmm3,nb304_rsqH1H2(%esp)

        movapd nb304_ixH2(%esp),%xmm0
        movapd nb304_iyH2(%esp),%xmm1
        movapd nb304_izH2(%esp),%xmm2
        movapd nb304_ixH2(%esp),%xmm3
        movapd nb304_iyH2(%esp),%xmm4
        movapd nb304_izH2(%esp),%xmm5
        subpd  nb304_jxM(%esp),%xmm0
        subpd  nb304_jyM(%esp),%xmm1
        subpd  nb304_jzM(%esp),%xmm2
        subpd  nb304_jxH1(%esp),%xmm3
        subpd  nb304_jyH1(%esp),%xmm4
        subpd  nb304_jzH1(%esp),%xmm5
        movapd %xmm0,nb304_dxH2M(%esp)
        movapd %xmm1,nb304_dyH2M(%esp)
        movapd %xmm2,nb304_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH2H1(%esp)
        movapd %xmm4,nb304_dyH2H1(%esp)
        movapd %xmm5,nb304_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb304_rsqH2M(%esp)
        movapd %xmm4,nb304_rsqH2H1(%esp)

        movapd nb304_ixH2(%esp),%xmm0
        movapd nb304_iyH2(%esp),%xmm1
        movapd nb304_izH2(%esp),%xmm2
        subpd  nb304_jxH2(%esp),%xmm0
        subpd  nb304_jyH2(%esp),%xmm1
        subpd  nb304_jzH2(%esp),%xmm2
        movapd %xmm0,nb304_dxH2H2(%esp)
        movapd %xmm1,nb304_dyH2H2(%esp)
        movapd %xmm2,nb304_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb304_rsqH2H2(%esp)

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
        movapd  nb304_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm3   ## iter1 
        mulpd   nb304_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm1   ## rinv 
        mulpd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvH2H2(%esp)
        movapd %xmm5,nb304_rinvH2H1(%esp)

        movapd nb304_rsqMM(%esp),%xmm0
        movapd nb304_rsqMH1(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb304_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm1   ## rinv 
        mulpd   nb304_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb304_rinvMM(%esp)
        movapd %xmm5,nb304_rinvMH1(%esp)

        movapd nb304_rsqMH2(%esp),%xmm0
        movapd nb304_rsqH1M(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm3   ## iter1 
        mulpd   nb304_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm1   ## rinv 
        mulpd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvMH2(%esp)
        movapd %xmm5,nb304_rinvH1M(%esp)

        movapd nb304_rsqH1H1(%esp),%xmm0
        movapd nb304_rsqH1H2(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm3   ## iter1a 
        mulpd   nb304_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304_half(%esp),%xmm1   ## rinv 
        mulpd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvH1H1(%esp)
        movapd %xmm5,nb304_rinvH1H2(%esp)

        movapd nb304_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb304_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb304_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMM(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb304_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movapd %xmm5,nb304_vctot(%esp)
        mulpd  nb304_tsc(%esp),%xmm3

        subpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm0      ## mult by rinv 

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb304_dxMM(%esp),%xmm0
        mulpd nb304_dyMM(%esp),%xmm1
        mulpd nb304_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixM(%esp),%xmm0
        addpd nb304_fiyM(%esp),%xmm1
        addpd nb304_fizM(%esp),%xmm2
        movapd %xmm3,nb304_fjxM(%esp)
        movapd %xmm4,nb304_fjyM(%esp)
        movapd %xmm5,nb304_fjzM(%esp)
        movapd %xmm0,nb304_fixM(%esp)
        movapd %xmm1,nb304_fiyM(%esp)
        movapd %xmm2,nb304_fizM(%esp)

        ## M-H1 interaction 
        movapd nb304_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb304_dxMH1(%esp),%xmm0
        mulpd nb304_dyMH1(%esp),%xmm1
        mulpd nb304_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixM(%esp),%xmm0
        addpd nb304_fiyM(%esp),%xmm1
        addpd nb304_fizM(%esp),%xmm2
        movapd %xmm3,nb304_fjxH1(%esp)
        movapd %xmm4,nb304_fjyH1(%esp)
        movapd %xmm5,nb304_fjzH1(%esp)
        movapd %xmm0,nb304_fixM(%esp)
        movapd %xmm1,nb304_fiyM(%esp)
        movapd %xmm2,nb304_fizM(%esp)

        ## M-H2 interaction  
        movapd nb304_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb304_dxMH2(%esp),%xmm0
        mulpd nb304_dyMH2(%esp),%xmm1
        mulpd nb304_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixM(%esp),%xmm0
        addpd nb304_fiyM(%esp),%xmm1
        addpd nb304_fizM(%esp),%xmm2
        movapd %xmm3,nb304_fjxH2(%esp)
        movapd %xmm4,nb304_fjyH2(%esp)
        movapd %xmm5,nb304_fjzH2(%esp)
        movapd %xmm0,nb304_fixM(%esp)
        movapd %xmm1,nb304_fiyM(%esp)
        movapd %xmm2,nb304_fizM(%esp)

        ## H1-M interaction 
        movapd nb304_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxM(%esp),%xmm3
        movapd nb304_fjyM(%esp),%xmm4
        movapd nb304_fjzM(%esp),%xmm5
        mulpd nb304_dxH1M(%esp),%xmm0
        mulpd nb304_dyH1M(%esp),%xmm1
        mulpd nb304_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH1(%esp),%xmm0
        addpd nb304_fiyH1(%esp),%xmm1
        addpd nb304_fizH1(%esp),%xmm2
        movapd %xmm3,nb304_fjxM(%esp)
        movapd %xmm4,nb304_fjyM(%esp)
        movapd %xmm5,nb304_fjzM(%esp)
        movapd %xmm0,nb304_fixH1(%esp)
        movapd %xmm1,nb304_fiyH1(%esp)
        movapd %xmm2,nb304_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb304_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH1(%esp),%xmm3
        movapd nb304_fjyH1(%esp),%xmm4
        movapd nb304_fjzH1(%esp),%xmm5
        mulpd nb304_dxH1H1(%esp),%xmm0
        mulpd nb304_dyH1H1(%esp),%xmm1
        mulpd nb304_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH1(%esp),%xmm0
        addpd nb304_fiyH1(%esp),%xmm1
        addpd nb304_fizH1(%esp),%xmm2
        movapd %xmm3,nb304_fjxH1(%esp)
        movapd %xmm4,nb304_fjyH1(%esp)
        movapd %xmm5,nb304_fjzH1(%esp)
        movapd %xmm0,nb304_fixH1(%esp)
        movapd %xmm1,nb304_fiyH1(%esp)
        movapd %xmm2,nb304_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb304_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH2(%esp),%xmm3
        movapd nb304_fjyH2(%esp),%xmm4
        movapd nb304_fjzH2(%esp),%xmm5
        mulpd nb304_dxH1H2(%esp),%xmm0
        mulpd nb304_dyH1H2(%esp),%xmm1
        mulpd nb304_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH1(%esp),%xmm0
        addpd nb304_fiyH1(%esp),%xmm1
        addpd nb304_fizH1(%esp),%xmm2
        movapd %xmm3,nb304_fjxH2(%esp)
        movapd %xmm4,nb304_fjyH2(%esp)
        movapd %xmm5,nb304_fjzH2(%esp)
        movapd %xmm0,nb304_fixH1(%esp)
        movapd %xmm1,nb304_fiyH1(%esp)
        movapd %xmm2,nb304_fizH1(%esp)

        ## H2-M interaction 
        movapd nb304_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxM(%esp),%xmm3
        movapd nb304_fjyM(%esp),%xmm4
        movapd nb304_fjzM(%esp),%xmm5
        mulpd nb304_dxH2M(%esp),%xmm0
        mulpd nb304_dyH2M(%esp),%xmm1
        mulpd nb304_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH2(%esp),%xmm0
        addpd nb304_fiyH2(%esp),%xmm1
        addpd nb304_fizH2(%esp),%xmm2
        movapd %xmm3,nb304_fjxM(%esp)
        movapd %xmm4,nb304_fjyM(%esp)
        movapd %xmm5,nb304_fjzM(%esp)
        movapd %xmm0,nb304_fixH2(%esp)
        movapd %xmm1,nb304_fiyH2(%esp)
        movapd %xmm2,nb304_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb304_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH1(%esp),%xmm3
        movapd nb304_fjyH1(%esp),%xmm4
        movapd nb304_fjzH1(%esp),%xmm5
        mulpd nb304_dxH2H1(%esp),%xmm0
        mulpd nb304_dyH2H1(%esp),%xmm1
        mulpd nb304_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH2(%esp),%xmm0
        addpd nb304_fiyH2(%esp),%xmm1
        addpd nb304_fizH2(%esp),%xmm2
        movapd %xmm3,nb304_fjxH1(%esp)
        movapd %xmm4,nb304_fjyH1(%esp)
        movapd %xmm5,nb304_fjzH1(%esp)
        movapd %xmm0,nb304_fixH2(%esp)
        movapd %xmm1,nb304_fiyH2(%esp)
        movapd %xmm2,nb304_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb304_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb304_vctot(%esp),%xmm5
    movapd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb304_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH2(%esp),%xmm3
        movapd nb304_fjyH2(%esp),%xmm4
        movapd nb304_fjzH2(%esp),%xmm5
        mulpd nb304_dxH2H2(%esp),%xmm0
        mulpd nb304_dyH2H2(%esp),%xmm1
        mulpd nb304_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb304_fixH2(%esp),%xmm0
        addpd nb304_fiyH2(%esp),%xmm1
        addpd nb304_fizH2(%esp),%xmm2
        movapd %xmm3,nb304_fjxH2(%esp)
        movapd %xmm4,nb304_fjyH2(%esp)
        movapd %xmm5,nb304_fjzH2(%esp)
        movapd %xmm0,nb304_fixH2(%esp)
        movapd %xmm1,nb304_fiyH2(%esp)
        movapd %xmm2,nb304_fizH2(%esp)

        movl nb304_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx

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
        addpd nb304_fjxH1(%esp),%xmm0
        addpd nb304_fjyH1(%esp),%xmm1
        addpd nb304_fjzH1(%esp),%xmm2
        addpd nb304_fjxH2(%esp),%xmm3
        addpd nb304_fjyH2(%esp),%xmm4
        addpd nb304_fjzH2(%esp),%xmm5
        addpd nb304_fjxM(%esp),%xmm6
        addpd nb304_fjyM(%esp),%xmm7
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
        addpd nb304_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)
        movhpd %xmm0,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb304_innerk(%esp)
        jl    _nb_kernel304_ia32_sse2.nb304_checksingle
        jmp   _nb_kernel304_ia32_sse2.nb304_unroll_loop
_nb_kernel304_ia32_sse2.nb304_checksingle: 
        movl  nb304_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel304_ia32_sse2.nb304_dosingle
        jmp   _nb_kernel304_ia32_sse2.nb304_updateouterdata
_nb_kernel304_ia32_sse2.nb304_dosingle: 
        movl  nb304_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb304_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb304_jxH1(%esp)
        movapd  %xmm3,nb304_jyH1(%esp)
        movapd  %xmm4,nb304_jzH1(%esp)
        movapd  %xmm5,nb304_jxH2(%esp)
        movapd  %xmm6,nb304_jyH2(%esp)
        movapd  %xmm7,nb304_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb304_jxM(%esp)
        movapd  %xmm3,nb304_jyM(%esp)
        movapd  %xmm4,nb304_jzM(%esp)

        movapd nb304_ixM(%esp),%xmm0
        movapd nb304_iyM(%esp),%xmm1
        movapd nb304_izM(%esp),%xmm2
        movapd nb304_ixM(%esp),%xmm3
        movapd nb304_iyM(%esp),%xmm4
        movapd nb304_izM(%esp),%xmm5
        subsd  nb304_jxM(%esp),%xmm0
        subsd  nb304_jyM(%esp),%xmm1
        subsd  nb304_jzM(%esp),%xmm2
        subsd  nb304_jxH1(%esp),%xmm3
        subsd  nb304_jyH1(%esp),%xmm4
        subsd  nb304_jzH1(%esp),%xmm5
        movapd %xmm0,nb304_dxMM(%esp)
        movapd %xmm1,nb304_dyMM(%esp)
        movapd %xmm2,nb304_dzMM(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxMH1(%esp)
        movapd %xmm4,nb304_dyMH1(%esp)
        movapd %xmm5,nb304_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqMM(%esp)
        movapd %xmm3,nb304_rsqMH1(%esp)

        movapd nb304_ixM(%esp),%xmm0
        movapd nb304_iyM(%esp),%xmm1
        movapd nb304_izM(%esp),%xmm2
        movapd nb304_ixH1(%esp),%xmm3
        movapd nb304_iyH1(%esp),%xmm4
        movapd nb304_izH1(%esp),%xmm5
        subsd  nb304_jxH2(%esp),%xmm0
        subsd  nb304_jyH2(%esp),%xmm1
        subsd  nb304_jzH2(%esp),%xmm2
        subsd  nb304_jxM(%esp),%xmm3
        subsd  nb304_jyM(%esp),%xmm4
        subsd  nb304_jzM(%esp),%xmm5
        movapd %xmm0,nb304_dxMH2(%esp)
        movapd %xmm1,nb304_dyMH2(%esp)
        movapd %xmm2,nb304_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH1M(%esp)
        movapd %xmm4,nb304_dyH1M(%esp)
        movapd %xmm5,nb304_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqMH2(%esp)
        movapd %xmm3,nb304_rsqH1M(%esp)

        movapd nb304_ixH1(%esp),%xmm0
        movapd nb304_iyH1(%esp),%xmm1
        movapd nb304_izH1(%esp),%xmm2
        movapd nb304_ixH1(%esp),%xmm3
        movapd nb304_iyH1(%esp),%xmm4
        movapd nb304_izH1(%esp),%xmm5
        subsd  nb304_jxH1(%esp),%xmm0
        subsd  nb304_jyH1(%esp),%xmm1
        subsd  nb304_jzH1(%esp),%xmm2
        subsd  nb304_jxH2(%esp),%xmm3
        subsd  nb304_jyH2(%esp),%xmm4
        subsd  nb304_jzH2(%esp),%xmm5
        movapd %xmm0,nb304_dxH1H1(%esp)
        movapd %xmm1,nb304_dyH1H1(%esp)
        movapd %xmm2,nb304_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH1H2(%esp)
        movapd %xmm4,nb304_dyH1H2(%esp)
        movapd %xmm5,nb304_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb304_rsqH1H1(%esp)
        movapd %xmm3,nb304_rsqH1H2(%esp)

        movapd nb304_ixH2(%esp),%xmm0
        movapd nb304_iyH2(%esp),%xmm1
        movapd nb304_izH2(%esp),%xmm2
        movapd nb304_ixH2(%esp),%xmm3
        movapd nb304_iyH2(%esp),%xmm4
        movapd nb304_izH2(%esp),%xmm5
        subsd  nb304_jxM(%esp),%xmm0
        subsd  nb304_jyM(%esp),%xmm1
        subsd  nb304_jzM(%esp),%xmm2
        subsd  nb304_jxH1(%esp),%xmm3
        subsd  nb304_jyH1(%esp),%xmm4
        subsd  nb304_jzH1(%esp),%xmm5
        movapd %xmm0,nb304_dxH2M(%esp)
        movapd %xmm1,nb304_dyH2M(%esp)
        movapd %xmm2,nb304_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb304_dxH2H1(%esp)
        movapd %xmm4,nb304_dyH2H1(%esp)
        movapd %xmm5,nb304_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb304_rsqH2M(%esp)
        movapd %xmm4,nb304_rsqH2H1(%esp)

        movapd nb304_ixH2(%esp),%xmm0
        movapd nb304_iyH2(%esp),%xmm1
        movapd nb304_izH2(%esp),%xmm2
        subsd  nb304_jxH2(%esp),%xmm0
        subsd  nb304_jyH2(%esp),%xmm1
        subsd  nb304_jzH2(%esp),%xmm2
        movapd %xmm0,nb304_dxH2H2(%esp)
        movapd %xmm1,nb304_dyH2H2(%esp)
        movapd %xmm2,nb304_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb304_rsqH2H2(%esp)

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
        movapd  nb304_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm3   ## iter1 
        mulsd   nb304_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm1   ## rinv 
        mulsd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvH2H2(%esp)
        movapd %xmm5,nb304_rinvH2H1(%esp)

        movapd nb304_rsqMM(%esp),%xmm0
        movapd nb304_rsqMH1(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb304_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm1   ## rinv 
        mulsd   nb304_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb304_rinvMM(%esp)
        movapd %xmm5,nb304_rinvMH1(%esp)

        movapd nb304_rsqMH2(%esp),%xmm0
        movapd nb304_rsqH1M(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm3   ## iter1 
        mulsd   nb304_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm1   ## rinv 
        mulsd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvMH2(%esp)
        movapd %xmm5,nb304_rinvH1M(%esp)

        movapd nb304_rsqH1H1(%esp),%xmm0
        movapd nb304_rsqH1H2(%esp),%xmm4
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
        movapd  nb304_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm3   ## iter1a 
        mulsd   nb304_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304_half(%esp),%xmm1   ## rinv 
        mulsd   nb304_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304_rinvH1H1(%esp)
        movapd %xmm5,nb304_rinvH1H2(%esp)

        movapd nb304_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb304_rinvH2M(%esp)

        movd %eax,%mm0
        ## start with MM interaction 
        movapd nb304_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMM(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb304_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movlpd %xmm5,nb304_vctot(%esp)
        mulsd  nb304_tsc(%esp),%xmm3

        subsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb304_dxMM(%esp),%xmm0
        mulsd nb304_dyMM(%esp),%xmm1
        mulsd nb304_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixM(%esp),%xmm0
        addsd nb304_fiyM(%esp),%xmm1
        addsd nb304_fizM(%esp),%xmm2
        movlpd %xmm3,nb304_fjxM(%esp)
        movlpd %xmm4,nb304_fjyM(%esp)
        movlpd %xmm5,nb304_fjzM(%esp)
        movlpd %xmm0,nb304_fixM(%esp)
        movlpd %xmm1,nb304_fiyM(%esp)
        movlpd %xmm2,nb304_fizM(%esp)

        ## M-H1 interaction 
        movapd nb304_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb304_dxMH1(%esp),%xmm0
        mulsd nb304_dyMH1(%esp),%xmm1
        mulsd nb304_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixM(%esp),%xmm0
        addsd nb304_fiyM(%esp),%xmm1
        addsd nb304_fizM(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH1(%esp)
        movlpd %xmm4,nb304_fjyH1(%esp)
        movlpd %xmm5,nb304_fjzH1(%esp)
        movlpd %xmm0,nb304_fixM(%esp)
        movlpd %xmm1,nb304_fiyM(%esp)
        movlpd %xmm2,nb304_fizM(%esp)

        ## M-H2 interaction  
        movapd nb304_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb304_dxMH2(%esp),%xmm0
        mulsd nb304_dyMH2(%esp),%xmm1
        mulsd nb304_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixM(%esp),%xmm0
        addsd nb304_fiyM(%esp),%xmm1
        addsd nb304_fizM(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH2(%esp)
        movlpd %xmm4,nb304_fjyH2(%esp)
        movlpd %xmm5,nb304_fjzH2(%esp)
        movlpd %xmm0,nb304_fixM(%esp)
        movlpd %xmm1,nb304_fiyM(%esp)
        movlpd %xmm2,nb304_fizM(%esp)

        ## H1-M interaction 
        movapd nb304_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxM(%esp),%xmm3
        movapd nb304_fjyM(%esp),%xmm4
        movapd nb304_fjzM(%esp),%xmm5
        mulsd nb304_dxH1M(%esp),%xmm0
        mulsd nb304_dyH1M(%esp),%xmm1
        mulsd nb304_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH1(%esp),%xmm0
        addsd nb304_fiyH1(%esp),%xmm1
        addsd nb304_fizH1(%esp),%xmm2
        movlpd %xmm3,nb304_fjxM(%esp)
        movlpd %xmm4,nb304_fjyM(%esp)
        movlpd %xmm5,nb304_fjzM(%esp)
        movlpd %xmm0,nb304_fixH1(%esp)
        movlpd %xmm1,nb304_fiyH1(%esp)
        movlpd %xmm2,nb304_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb304_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH1(%esp),%xmm3
        movapd nb304_fjyH1(%esp),%xmm4
        movapd nb304_fjzH1(%esp),%xmm5
        mulsd nb304_dxH1H1(%esp),%xmm0
        mulsd nb304_dyH1H1(%esp),%xmm1
        mulsd nb304_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH1(%esp),%xmm0
        addsd nb304_fiyH1(%esp),%xmm1
        addsd nb304_fizH1(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH1(%esp)
        movlpd %xmm4,nb304_fjyH1(%esp)
        movlpd %xmm5,nb304_fjzH1(%esp)
        movlpd %xmm0,nb304_fixH1(%esp)
        movlpd %xmm1,nb304_fiyH1(%esp)
        movlpd %xmm2,nb304_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb304_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH2(%esp),%xmm3
        movapd nb304_fjyH2(%esp),%xmm4
        movapd nb304_fjzH2(%esp),%xmm5
        mulsd nb304_dxH1H2(%esp),%xmm0
        mulsd nb304_dyH1H2(%esp),%xmm1
        mulsd nb304_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH1(%esp),%xmm0
        addsd nb304_fiyH1(%esp),%xmm1
        addsd nb304_fizH1(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH2(%esp)
        movlpd %xmm4,nb304_fjyH2(%esp)
        movlpd %xmm5,nb304_fjzH2(%esp)
        movlpd %xmm0,nb304_fixH1(%esp)
        movlpd %xmm1,nb304_fiyH1(%esp)
        movlpd %xmm2,nb304_fizH1(%esp)

        ## H2-M interaction 
        movapd nb304_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxM(%esp),%xmm3
        movapd nb304_fjyM(%esp),%xmm4
        movapd nb304_fjzM(%esp),%xmm5
        mulsd nb304_dxH2M(%esp),%xmm0
        mulsd nb304_dyH2M(%esp),%xmm1
        mulsd nb304_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH2(%esp),%xmm0
        addsd nb304_fiyH2(%esp),%xmm1
        addsd nb304_fizH2(%esp),%xmm2
        movlpd %xmm3,nb304_fjxM(%esp)
        movlpd %xmm4,nb304_fjyM(%esp)
        movlpd %xmm5,nb304_fjzM(%esp)
        movlpd %xmm0,nb304_fixH2(%esp)
        movlpd %xmm1,nb304_fiyH2(%esp)
        movlpd %xmm2,nb304_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb304_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH1(%esp),%xmm3
        movapd nb304_fjyH1(%esp),%xmm4
        movapd nb304_fjzH1(%esp),%xmm5
        mulsd nb304_dxH2H1(%esp),%xmm0
        mulsd nb304_dyH2H1(%esp),%xmm1
        mulsd nb304_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH2(%esp),%xmm0
        addsd nb304_fiyH2(%esp),%xmm1
        addsd nb304_fizH2(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH1(%esp)
        movlpd %xmm4,nb304_fjyH1(%esp)
        movlpd %xmm5,nb304_fjzH1(%esp)
        movlpd %xmm0,nb304_fixH2(%esp)
        movlpd %xmm1,nb304_fiyH2(%esp)
        movlpd %xmm2,nb304_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb304_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb304_two(%esp),%xmm7    ## two*Heps2 
        movapd nb304_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb304_vctot(%esp),%xmm5
    movlpd %xmm5,nb304_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb304_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb304_fjxH2(%esp),%xmm3
        movapd nb304_fjyH2(%esp),%xmm4
        movapd nb304_fjzH2(%esp),%xmm5
        mulsd nb304_dxH2H2(%esp),%xmm0
        mulsd nb304_dyH2H2(%esp),%xmm1
        mulsd nb304_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb304_fixH2(%esp),%xmm0
        addsd nb304_fiyH2(%esp),%xmm1
        addsd nb304_fizH2(%esp),%xmm2
        movlpd %xmm3,nb304_fjxH2(%esp)
        movlpd %xmm4,nb304_fjyH2(%esp)
        movlpd %xmm5,nb304_fjzH2(%esp)
        movlpd %xmm0,nb304_fixH2(%esp)
        movlpd %xmm1,nb304_fiyH2(%esp)
        movlpd %xmm2,nb304_fizH2(%esp)

        movl nb304_faction(%ebp),%edi

        movd %mm0,%eax

        ## Did all interactions - now update j forces 
        movlpd 24(%edi,%eax,8),%xmm0
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 40(%edi,%eax,8),%xmm2
        movlpd 48(%edi,%eax,8),%xmm3
        movlpd 56(%edi,%eax,8),%xmm4
        movlpd 64(%edi,%eax,8),%xmm5
        movlpd 72(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%eax,8),%xmm7
        addsd nb304_fjxH1(%esp),%xmm0
        addsd nb304_fjyH1(%esp),%xmm1
        addsd nb304_fjzH1(%esp),%xmm2
        addsd nb304_fjxH2(%esp),%xmm3
        addsd nb304_fjyH2(%esp),%xmm4
        addsd nb304_fjzH2(%esp),%xmm5
        addsd nb304_fjxM(%esp),%xmm6
        addsd nb304_fjyM(%esp),%xmm7
        movlpd %xmm0,24(%edi,%eax,8)
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movlpd %xmm4,56(%edi,%eax,8)
        movlpd %xmm5,64(%edi,%eax,8)
        movlpd %xmm6,72(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%eax,8)

        movlpd 88(%edi,%eax,8),%xmm0
        addsd nb304_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)

_nb_kernel304_ia32_sse2.nb304_updateouterdata: 
        movl  nb304_ii3(%esp),%ecx
        movl  nb304_faction(%ebp),%edi
        movl  nb304_fshift(%ebp),%esi
        movl  nb304_is3(%esp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb304_fixH1(%esp),%xmm0
        movapd nb304_fiyH1(%esp),%xmm1
        movapd nb304_fizH1(%esp),%xmm2

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
        movapd nb304_fixH2(%esp),%xmm0
        movapd nb304_fiyH2(%esp),%xmm1
        movapd nb304_fizH2(%esp),%xmm2

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
        movapd nb304_fixM(%esp),%xmm0
        movapd nb304_fiyM(%esp),%xmm1
        movapd nb304_fizM(%esp),%xmm2

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
        movl nb304_n(%esp),%esi
        ## get group index for i particle 
        movl  nb304_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb304_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb304_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb304_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304_ia32_sse2.nb304_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304_n(%esp)
        jmp _nb_kernel304_ia32_sse2.nb304_outer
_nb_kernel304_ia32_sse2.nb304_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304_ia32_sse2.nb304_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304_ia32_sse2.nb304_threadloop
_nb_kernel304_ia32_sse2.nb304_end: 
        emms

        movl nb304_nouter(%esp),%eax
        movl nb304_ninner(%esp),%ebx
        movl nb304_outeriter(%ebp),%ecx
        movl nb304_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb304_salign(%esp),%eax
        addl %eax,%esp
        addl $1464,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret


.globl nb_kernel304nf_ia32_sse2
.globl _nb_kernel304nf_ia32_sse2
nb_kernel304nf_ia32_sse2:       
_nb_kernel304nf_ia32_sse2:      
.set nb304nf_p_nri, 8
.set nb304nf_iinr, 12
.set nb304nf_jindex, 16
.set nb304nf_jjnr, 20
.set nb304nf_shift, 24
.set nb304nf_shiftvec, 28
.set nb304nf_fshift, 32
.set nb304nf_gid, 36
.set nb304nf_pos, 40
.set nb304nf_faction, 44
.set nb304nf_charge, 48
.set nb304nf_p_facel, 52
.set nb304nf_argkrf, 56
.set nb304nf_argcrf, 60
.set nb304nf_Vc, 64
.set nb304nf_type, 68
.set nb304nf_p_ntype, 72
.set nb304nf_vdwparam, 76
.set nb304nf_Vvdw, 80
.set nb304nf_p_tabscale, 84
.set nb304nf_VFtab, 88
.set nb304nf_invsqrta, 92
.set nb304nf_dvda, 96
.set nb304nf_p_gbtabscale, 100
.set nb304nf_GBtab, 104
.set nb304nf_p_nthreads, 108
.set nb304nf_count, 112
.set nb304nf_mtx, 116
.set nb304nf_outeriter, 120
.set nb304nf_inneriter, 124
.set nb304nf_work, 128
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
.set nb304nf_innerjjnr, 696
.set nb304nf_innerk, 700
.set nb304nf_n, 704
.set nb304nf_nn1, 708
.set nb304nf_nri, 712
.set nb304nf_nouter, 716
.set nb304nf_ninner, 720
.set nb304nf_salign, 724
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $728,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb304nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb304nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb304nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb304nf_nouter(%esp)
        movl %eax,nb304nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb304nf_half(%esp)
        movl %ebx,nb304nf_half+4(%esp)
        movsd nb304nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb304nf_half(%esp)
        movapd %xmm3,nb304nf_three(%esp)

        movl nb304nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb304nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb304nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb304nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb304nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb304nf_qqMM(%esp)
        movapd %xmm4,nb304nf_qqMH(%esp)
        movapd %xmm5,nb304nf_qqHH(%esp)

_nb_kernel304nf_ia32_sse2.nb304nf_threadloop: 
        movl  nb304nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel304nf_ia32_sse2.nb304nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel304nf_ia32_sse2.nb304nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb304nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb304nf_n(%esp)
        movl %ebx,nb304nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel304nf_ia32_sse2.nb304nf_outerstart
        jmp _nb_kernel304nf_ia32_sse2.nb304nf_end

_nb_kernel304nf_ia32_sse2.nb304nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb304nf_nouter(%esp),%ebx
        movl %ebx,nb304nf_nouter(%esp)

_nb_kernel304nf_ia32_sse2.nb304nf_outer: 
        movl  nb304nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb304nf_is3(%esp)            ## store is3 

        movl  nb304nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb304nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb304nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb304nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb304nf_ixH1(%esp)
        movapd %xmm4,nb304nf_iyH1(%esp)
        movapd %xmm5,nb304nf_izH1(%esp)

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
        movapd %xmm0,nb304nf_ixH2(%esp)
        movapd %xmm1,nb304nf_iyH2(%esp)
        movapd %xmm2,nb304nf_izH2(%esp)
        movapd %xmm3,nb304nf_ixM(%esp)
        movapd %xmm4,nb304nf_iyM(%esp)
        movapd %xmm5,nb304nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb304nf_vctot(%esp)

        movl  nb304nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb304nf_pos(%ebp),%esi
        movl  nb304nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb304nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb304nf_ninner(%esp),%ecx
        movl  %ecx,nb304nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb304nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel304nf_ia32_sse2.nb304nf_unroll_loop
        jmp   _nb_kernel304nf_ia32_sse2.nb304nf_checksingle
_nb_kernel304nf_ia32_sse2.nb304nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb304nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb304nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb304nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb304nf_jxH1(%esp)
        movapd  %xmm3,nb304nf_jyH1(%esp)
        movapd  %xmm4,nb304nf_jzH1(%esp)
        movapd  %xmm5,nb304nf_jxH2(%esp)
        movapd  %xmm6,nb304nf_jyH2(%esp)
        movapd  %xmm7,nb304nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb304nf_jxM(%esp)
        movapd  %xmm3,nb304nf_jyM(%esp)
        movapd  %xmm4,nb304nf_jzM(%esp)

        movapd nb304nf_ixM(%esp),%xmm0
        movapd nb304nf_iyM(%esp),%xmm1
        movapd nb304nf_izM(%esp),%xmm2
        movapd nb304nf_ixM(%esp),%xmm3
        movapd nb304nf_iyM(%esp),%xmm4
        movapd nb304nf_izM(%esp),%xmm5
        subpd  nb304nf_jxM(%esp),%xmm0
        subpd  nb304nf_jyM(%esp),%xmm1
        subpd  nb304nf_jzM(%esp),%xmm2
        subpd  nb304nf_jxH1(%esp),%xmm3
        subpd  nb304nf_jyH1(%esp),%xmm4
        subpd  nb304nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMM(%esp)
        movapd %xmm3,nb304nf_rsqMH1(%esp)

        movapd nb304nf_ixM(%esp),%xmm0
        movapd nb304nf_iyM(%esp),%xmm1
        movapd nb304nf_izM(%esp),%xmm2
        movapd nb304nf_ixH1(%esp),%xmm3
        movapd nb304nf_iyH1(%esp),%xmm4
        movapd nb304nf_izH1(%esp),%xmm5
        subpd  nb304nf_jxH2(%esp),%xmm0
        subpd  nb304nf_jyH2(%esp),%xmm1
        subpd  nb304nf_jzH2(%esp),%xmm2
        subpd  nb304nf_jxM(%esp),%xmm3
        subpd  nb304nf_jyM(%esp),%xmm4
        subpd  nb304nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMH2(%esp)
        movapd %xmm3,nb304nf_rsqH1M(%esp)

        movapd nb304nf_ixH1(%esp),%xmm0
        movapd nb304nf_iyH1(%esp),%xmm1
        movapd nb304nf_izH1(%esp),%xmm2
        movapd nb304nf_ixH1(%esp),%xmm3
        movapd nb304nf_iyH1(%esp),%xmm4
        movapd nb304nf_izH1(%esp),%xmm5
        subpd  nb304nf_jxH1(%esp),%xmm0
        subpd  nb304nf_jyH1(%esp),%xmm1
        subpd  nb304nf_jzH1(%esp),%xmm2
        subpd  nb304nf_jxH2(%esp),%xmm3
        subpd  nb304nf_jyH2(%esp),%xmm4
        subpd  nb304nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH1H1(%esp)
        movapd %xmm3,nb304nf_rsqH1H2(%esp)

        movapd nb304nf_ixH2(%esp),%xmm0
        movapd nb304nf_iyH2(%esp),%xmm1
        movapd nb304nf_izH2(%esp),%xmm2
        movapd nb304nf_ixH2(%esp),%xmm3
        movapd nb304nf_iyH2(%esp),%xmm4
        movapd nb304nf_izH2(%esp),%xmm5
        subpd  nb304nf_jxM(%esp),%xmm0
        subpd  nb304nf_jyM(%esp),%xmm1
        subpd  nb304nf_jzM(%esp),%xmm2
        subpd  nb304nf_jxH1(%esp),%xmm3
        subpd  nb304nf_jyH1(%esp),%xmm4
        subpd  nb304nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH2M(%esp)
        movapd %xmm4,nb304nf_rsqH2H1(%esp)

        movapd nb304nf_ixH2(%esp),%xmm0
        movapd nb304nf_iyH2(%esp),%xmm1
        movapd nb304nf_izH2(%esp),%xmm2
        subpd  nb304nf_jxH2(%esp),%xmm0
        subpd  nb304nf_jyH2(%esp),%xmm1
        subpd  nb304nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb304nf_rsqH2H2(%esp)

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
        movapd  nb304nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb304nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH2H2(%esp)
        movapd %xmm5,nb304nf_rinvH2H1(%esp)

        movapd nb304nf_rsqMM(%esp),%xmm0
        movapd nb304nf_rsqMH1(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb304nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb304nf_rinvMM(%esp)
        movapd %xmm5,nb304nf_rinvMH1(%esp)

        movapd nb304nf_rsqMH2(%esp),%xmm0
        movapd nb304nf_rsqH1M(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb304nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvMH2(%esp)
        movapd %xmm5,nb304nf_rinvH1M(%esp)

        movapd nb304nf_rsqH1H1(%esp),%xmm0
        movapd nb304nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb304nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH1H1(%esp)
        movapd %xmm5,nb304nf_rinvH1H2(%esp)

        movapd nb304nf_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb304nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb304nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb304nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMM(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb304nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## M-H2 interaction  
        movapd nb304nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb304nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb304nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb304nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb304nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb304nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb304nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb304nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb304nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb304nf_vctot(%esp),%xmm5
    movapd %xmm5,nb304nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb304nf_innerk(%esp)
        jl    _nb_kernel304nf_ia32_sse2.nb304nf_checksingle
        jmp   _nb_kernel304nf_ia32_sse2.nb304nf_unroll_loop
_nb_kernel304nf_ia32_sse2.nb304nf_checksingle: 
        movl  nb304nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel304nf_ia32_sse2.nb304nf_dosingle
        jmp   _nb_kernel304nf_ia32_sse2.nb304nf_updateouterdata
_nb_kernel304nf_ia32_sse2.nb304nf_dosingle: 
        movl  nb304nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb304nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb304nf_jxH1(%esp)
        movapd  %xmm3,nb304nf_jyH1(%esp)
        movapd  %xmm4,nb304nf_jzH1(%esp)
        movapd  %xmm5,nb304nf_jxH2(%esp)
        movapd  %xmm6,nb304nf_jyH2(%esp)
        movapd  %xmm7,nb304nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb304nf_jxM(%esp)
        movapd  %xmm3,nb304nf_jyM(%esp)
        movapd  %xmm4,nb304nf_jzM(%esp)

        movapd nb304nf_ixM(%esp),%xmm0
        movapd nb304nf_iyM(%esp),%xmm1
        movapd nb304nf_izM(%esp),%xmm2
        movapd nb304nf_ixM(%esp),%xmm3
        movapd nb304nf_iyM(%esp),%xmm4
        movapd nb304nf_izM(%esp),%xmm5
        subsd  nb304nf_jxM(%esp),%xmm0
        subsd  nb304nf_jyM(%esp),%xmm1
        subsd  nb304nf_jzM(%esp),%xmm2
        subsd  nb304nf_jxH1(%esp),%xmm3
        subsd  nb304nf_jyH1(%esp),%xmm4
        subsd  nb304nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMM(%esp)
        movapd %xmm3,nb304nf_rsqMH1(%esp)

        movapd nb304nf_ixM(%esp),%xmm0
        movapd nb304nf_iyM(%esp),%xmm1
        movapd nb304nf_izM(%esp),%xmm2
        movapd nb304nf_ixH1(%esp),%xmm3
        movapd nb304nf_iyH1(%esp),%xmm4
        movapd nb304nf_izH1(%esp),%xmm5
        subsd  nb304nf_jxH2(%esp),%xmm0
        subsd  nb304nf_jyH2(%esp),%xmm1
        subsd  nb304nf_jzH2(%esp),%xmm2
        subsd  nb304nf_jxM(%esp),%xmm3
        subsd  nb304nf_jyM(%esp),%xmm4
        subsd  nb304nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqMH2(%esp)
        movapd %xmm3,nb304nf_rsqH1M(%esp)

        movapd nb304nf_ixH1(%esp),%xmm0
        movapd nb304nf_iyH1(%esp),%xmm1
        movapd nb304nf_izH1(%esp),%xmm2
        movapd nb304nf_ixH1(%esp),%xmm3
        movapd nb304nf_iyH1(%esp),%xmm4
        movapd nb304nf_izH1(%esp),%xmm5
        subsd  nb304nf_jxH1(%esp),%xmm0
        subsd  nb304nf_jyH1(%esp),%xmm1
        subsd  nb304nf_jzH1(%esp),%xmm2
        subsd  nb304nf_jxH2(%esp),%xmm3
        subsd  nb304nf_jyH2(%esp),%xmm4
        subsd  nb304nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH1H1(%esp)
        movapd %xmm3,nb304nf_rsqH1H2(%esp)

        movapd nb304nf_ixH2(%esp),%xmm0
        movapd nb304nf_iyH2(%esp),%xmm1
        movapd nb304nf_izH2(%esp),%xmm2
        movapd nb304nf_ixH2(%esp),%xmm3
        movapd nb304nf_iyH2(%esp),%xmm4
        movapd nb304nf_izH2(%esp),%xmm5
        subsd  nb304nf_jxM(%esp),%xmm0
        subsd  nb304nf_jyM(%esp),%xmm1
        subsd  nb304nf_jzM(%esp),%xmm2
        subsd  nb304nf_jxH1(%esp),%xmm3
        subsd  nb304nf_jyH1(%esp),%xmm4
        subsd  nb304nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb304nf_rsqH2M(%esp)
        movapd %xmm4,nb304nf_rsqH2H1(%esp)

        movapd nb304nf_ixH2(%esp),%xmm0
        movapd nb304nf_iyH2(%esp),%xmm1
        movapd nb304nf_izH2(%esp),%xmm2
        subsd  nb304nf_jxH2(%esp),%xmm0
        subsd  nb304nf_jyH2(%esp),%xmm1
        subsd  nb304nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb304nf_rsqH2H2(%esp)

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
        movapd  nb304nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb304nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH2H2(%esp)
        movapd %xmm5,nb304nf_rinvH2H1(%esp)

        movapd nb304nf_rsqMM(%esp),%xmm0
        movapd nb304nf_rsqMH1(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb304nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb304nf_rinvMM(%esp)
        movapd %xmm5,nb304nf_rinvMH1(%esp)

        movapd nb304nf_rsqMH2(%esp),%xmm0
        movapd nb304nf_rsqH1M(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb304nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvMH2(%esp)
        movapd %xmm5,nb304nf_rinvH1M(%esp)

        movapd nb304nf_rsqH1H1(%esp),%xmm0
        movapd nb304nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb304nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb304nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb304nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb304nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb304nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb304nf_rinvH1H1(%esp)
        movapd %xmm5,nb304nf_rinvH1H2(%esp)

        movapd nb304nf_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb304nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb304nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb304nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb304nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb304nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMM(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb304nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## M-H2 interaction  
        movapd nb304nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb304nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb304nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb304nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb304nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb304nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb304nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb304nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb304nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb304nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb304nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb304nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb304nf_vctot(%esp)

_nb_kernel304nf_ia32_sse2.nb304nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb304nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb304nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb304nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb304nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb304nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel304nf_ia32_sse2.nb304nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb304nf_n(%esp)
        jmp _nb_kernel304nf_ia32_sse2.nb304nf_outer
_nb_kernel304nf_ia32_sse2.nb304nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb304nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel304nf_ia32_sse2.nb304nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel304nf_ia32_sse2.nb304nf_threadloop
_nb_kernel304nf_ia32_sse2.nb304nf_end: 
        emms

        movl nb304nf_nouter(%esp),%eax
        movl nb304nf_ninner(%esp),%ebx
        movl nb304nf_outeriter(%ebp),%ecx
        movl nb304nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb304nf_salign(%esp),%eax
        addl %eax,%esp
        addl $728,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



