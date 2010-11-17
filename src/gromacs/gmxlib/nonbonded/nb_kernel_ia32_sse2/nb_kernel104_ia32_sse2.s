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



.globl nb_kernel104_ia32_sse2
.globl _nb_kernel104_ia32_sse2
nb_kernel104_ia32_sse2: 
_nb_kernel104_ia32_sse2:        
.set nb104_p_nri, 8
.set nb104_iinr, 12
.set nb104_jindex, 16
.set nb104_jjnr, 20
.set nb104_shift, 24
.set nb104_shiftvec, 28
.set nb104_fshift, 32
.set nb104_gid, 36
.set nb104_pos, 40
.set nb104_faction, 44
.set nb104_charge, 48
.set nb104_p_facel, 52
.set nb104_argkrf, 56
.set nb104_argcrf, 60
.set nb104_Vc, 64
.set nb104_type, 68
.set nb104_p_ntype, 72
.set nb104_vdwparam, 76
.set nb104_Vvdw, 80
.set nb104_p_tabscale, 84
.set nb104_VFtab, 88
.set nb104_invsqrta, 92
.set nb104_dvda, 96
.set nb104_p_gbtabscale, 100
.set nb104_GBtab, 104
.set nb104_p_nthreads, 108
.set nb104_count, 112
.set nb104_mtx, 116
.set nb104_outeriter, 120
.set nb104_inneriter, 124
.set nb104_work, 128
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
.set nb104_is3, 1392
.set nb104_ii3, 1396
.set nb104_innerjjnr, 1400
.set nb104_innerk, 1404
.set nb104_n, 1408
.set nb104_nn1, 1412
.set nb104_nri, 1416
.set nb104_nouter, 1420
.set nb104_ninner, 1424
.set nb104_salign, 1428
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1432,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb104_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb104_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb104_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb104_nouter(%esp)
        movl %eax,nb104_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb104_half(%esp)
        movl %ebx,nb104_half+4(%esp)
        movsd nb104_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb104_half(%esp)
        movapd %xmm3,nb104_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb104_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb104_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3     ## qM 
        movsd %xmm3,%xmm4               ## qM 
        movsd 8(%edx,%ebx,8),%xmm5      ## qH 
        movl nb104_p_facel(%ebp),%esi
        movsd (%esi),%xmm6      ## facel 
        mulsd  %xmm3,%xmm3              ## qM*qM 
        mulsd  %xmm5,%xmm4              ## qM*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104_qqMM(%esp)
        movapd %xmm4,nb104_qqMH(%esp)
        movapd %xmm5,nb104_qqHH(%esp)

_nb_kernel104_ia32_sse2.nb104_threadloop: 
        movl  nb104_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel104_ia32_sse2.nb104_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104_ia32_sse2.nb104_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104_n(%esp)
        movl %ebx,nb104_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104_ia32_sse2.nb104_outerstart
        jmp _nb_kernel104_ia32_sse2.nb104_end

_nb_kernel104_ia32_sse2.nb104_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104_nouter(%esp),%ebx
        movl %ebx,nb104_nouter(%esp)

_nb_kernel104_ia32_sse2.nb104_outer: 
        movl  nb104_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb104_is3(%esp)      ## store is3 

        movl  nb104_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb104_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb104_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb104_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104_ixH1(%esp)
        movapd %xmm4,nb104_iyH1(%esp)
        movapd %xmm5,nb104_izH1(%esp)

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
        movapd %xmm0,nb104_ixH2(%esp)
        movapd %xmm1,nb104_iyH2(%esp)
        movapd %xmm2,nb104_izH2(%esp)
        movapd %xmm3,nb104_ixM(%esp)
        movapd %xmm4,nb104_iyM(%esp)
        movapd %xmm5,nb104_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb104_vctot(%esp)
        movapd %xmm4,nb104_fixM(%esp)
        movapd %xmm4,nb104_fiyM(%esp)
        movapd %xmm4,nb104_fizM(%esp)
        movapd %xmm4,nb104_fixH1(%esp)
        movapd %xmm4,nb104_fiyH1(%esp)
        movapd %xmm4,nb104_fizH1(%esp)
        movapd %xmm4,nb104_fixH2(%esp)
        movapd %xmm4,nb104_fiyH2(%esp)
        movapd %xmm4,nb104_fizH2(%esp)

        movl  nb104_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx     ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb104_pos(%ebp),%esi
        movl  nb104_faction(%ebp),%edi
        movl  nb104_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb104_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb104_ninner(%esp),%ecx
        movl  %ecx,nb104_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb104_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel104_ia32_sse2.nb104_unroll_loop
        jmp   _nb_kernel104_ia32_sse2.nb104_checksingle
_nb_kernel104_ia32_sse2.nb104_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb104_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb104_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb104_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb104_jxH1(%esp)
        movapd  %xmm3,nb104_jyH1(%esp)
        movapd  %xmm4,nb104_jzH1(%esp)
        movapd  %xmm5,nb104_jxH2(%esp)
        movapd  %xmm6,nb104_jyH2(%esp)
        movapd  %xmm7,nb104_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb104_jxM(%esp)
        movapd  %xmm3,nb104_jyM(%esp)
        movapd  %xmm4,nb104_jzM(%esp)

        movapd nb104_ixM(%esp),%xmm0
        movapd nb104_iyM(%esp),%xmm1
        movapd nb104_izM(%esp),%xmm2
        movapd nb104_ixM(%esp),%xmm3
        movapd nb104_iyM(%esp),%xmm4
        movapd nb104_izM(%esp),%xmm5
        subpd  nb104_jxM(%esp),%xmm0
        subpd  nb104_jyM(%esp),%xmm1
        subpd  nb104_jzM(%esp),%xmm2
        subpd  nb104_jxH1(%esp),%xmm3
        subpd  nb104_jyH1(%esp),%xmm4
        subpd  nb104_jzH1(%esp),%xmm5
        movapd %xmm0,nb104_dxMM(%esp)
        movapd %xmm1,nb104_dyMM(%esp)
        movapd %xmm2,nb104_dzMM(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxMH1(%esp)
        movapd %xmm4,nb104_dyMH1(%esp)
        movapd %xmm5,nb104_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqMM(%esp)
        movapd %xmm3,nb104_rsqMH1(%esp)

        movapd nb104_ixM(%esp),%xmm0
        movapd nb104_iyM(%esp),%xmm1
        movapd nb104_izM(%esp),%xmm2
        movapd nb104_ixH1(%esp),%xmm3
        movapd nb104_iyH1(%esp),%xmm4
        movapd nb104_izH1(%esp),%xmm5
        subpd  nb104_jxH2(%esp),%xmm0
        subpd  nb104_jyH2(%esp),%xmm1
        subpd  nb104_jzH2(%esp),%xmm2
        subpd  nb104_jxM(%esp),%xmm3
        subpd  nb104_jyM(%esp),%xmm4
        subpd  nb104_jzM(%esp),%xmm5
        movapd %xmm0,nb104_dxMH2(%esp)
        movapd %xmm1,nb104_dyMH2(%esp)
        movapd %xmm2,nb104_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH1M(%esp)
        movapd %xmm4,nb104_dyH1M(%esp)
        movapd %xmm5,nb104_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqMH2(%esp)
        movapd %xmm3,nb104_rsqH1M(%esp)

        movapd nb104_ixH1(%esp),%xmm0
        movapd nb104_iyH1(%esp),%xmm1
        movapd nb104_izH1(%esp),%xmm2
        movapd nb104_ixH1(%esp),%xmm3
        movapd nb104_iyH1(%esp),%xmm4
        movapd nb104_izH1(%esp),%xmm5
        subpd  nb104_jxH1(%esp),%xmm0
        subpd  nb104_jyH1(%esp),%xmm1
        subpd  nb104_jzH1(%esp),%xmm2
        subpd  nb104_jxH2(%esp),%xmm3
        subpd  nb104_jyH2(%esp),%xmm4
        subpd  nb104_jzH2(%esp),%xmm5
        movapd %xmm0,nb104_dxH1H1(%esp)
        movapd %xmm1,nb104_dyH1H1(%esp)
        movapd %xmm2,nb104_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH1H2(%esp)
        movapd %xmm4,nb104_dyH1H2(%esp)
        movapd %xmm5,nb104_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqH1H1(%esp)
        movapd %xmm3,nb104_rsqH1H2(%esp)

        movapd nb104_ixH2(%esp),%xmm0
        movapd nb104_iyH2(%esp),%xmm1
        movapd nb104_izH2(%esp),%xmm2
        movapd nb104_ixH2(%esp),%xmm3
        movapd nb104_iyH2(%esp),%xmm4
        movapd nb104_izH2(%esp),%xmm5
        subpd  nb104_jxM(%esp),%xmm0
        subpd  nb104_jyM(%esp),%xmm1
        subpd  nb104_jzM(%esp),%xmm2
        subpd  nb104_jxH1(%esp),%xmm3
        subpd  nb104_jyH1(%esp),%xmm4
        subpd  nb104_jzH1(%esp),%xmm5
        movapd %xmm0,nb104_dxH2M(%esp)
        movapd %xmm1,nb104_dyH2M(%esp)
        movapd %xmm2,nb104_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH2H1(%esp)
        movapd %xmm4,nb104_dyH2H1(%esp)
        movapd %xmm5,nb104_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb104_rsqH2M(%esp)
        movapd %xmm4,nb104_rsqH2H1(%esp)

        movapd nb104_ixH2(%esp),%xmm0
        movapd nb104_iyH2(%esp),%xmm1
        movapd nb104_izH2(%esp),%xmm2
        subpd  nb104_jxH2(%esp),%xmm0
        subpd  nb104_jyH2(%esp),%xmm1
        subpd  nb104_jzH2(%esp),%xmm2
        movapd %xmm0,nb104_dxH2H2(%esp)
        movapd %xmm1,nb104_dyH2H2(%esp)
        movapd %xmm2,nb104_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb104_rsqH2H2(%esp)

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
        movapd  nb104_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm3   ## iter1 
        mulpd   nb104_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm1   ## rinv 
        mulpd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvH2H2(%esp)
        movapd %xmm5,nb104_rinvH2H1(%esp)

        movapd nb104_rsqMM(%esp),%xmm0
        movapd nb104_rsqMH1(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb104_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm1   ## rinv 
        mulpd   nb104_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb104_rinvMM(%esp)
        movapd %xmm5,nb104_rinvMH1(%esp)

        movapd nb104_rsqMH2(%esp),%xmm0
        movapd nb104_rsqH1M(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm3   ## iter1 
        mulpd   nb104_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm1   ## rinv 
        mulpd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvMH2(%esp)
        movapd %xmm5,nb104_rinvH1M(%esp)

        movapd nb104_rsqH1H1(%esp),%xmm0
        movapd nb104_rsqH1H2(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm3   ## iter1a 
        mulpd   nb104_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104_half(%esp),%xmm1   ## rinv 
        mulpd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvH1H1(%esp)
        movapd %xmm5,nb104_rinvH1H2(%esp)

        movapd nb104_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb104_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb104_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulpd  %xmm0,%xmm0              ## rinvsq 
        mulpd  nb104_qqMM(%esp),%xmm7
        mulpd  %xmm7,%xmm0
        addpd  nb104_vctot(%esp),%xmm7
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb104_dxMM(%esp),%xmm0
        mulpd nb104_dyMM(%esp),%xmm1
        mulpd nb104_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixM(%esp),%xmm0
        addpd nb104_fiyM(%esp),%xmm1
        addpd nb104_fizM(%esp),%xmm2
        movapd %xmm3,nb104_fjxM(%esp)
        movapd %xmm4,nb104_fjyM(%esp)
        movapd %xmm5,nb104_fjzM(%esp)
        movapd %xmm0,nb104_fixM(%esp)
        movapd %xmm1,nb104_fiyM(%esp)
        movapd %xmm2,nb104_fizM(%esp)

        ## M-H1 interaction 
        movapd nb104_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqMH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsMH1  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb104_dxMH1(%esp),%xmm0
        mulpd nb104_dyMH1(%esp),%xmm1
        mulpd nb104_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixM(%esp),%xmm0
        addpd nb104_fiyM(%esp),%xmm1
        addpd nb104_fizM(%esp),%xmm2
        movapd %xmm3,nb104_fjxH1(%esp)
        movapd %xmm4,nb104_fjyH1(%esp)
        movapd %xmm5,nb104_fjzH1(%esp)
        movapd %xmm0,nb104_fixM(%esp)
        movapd %xmm1,nb104_fiyM(%esp)
        movapd %xmm2,nb104_fizM(%esp)

        ## M-H2 interaction  
        movapd nb104_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqMH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsMH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb104_dxMH2(%esp),%xmm0
        mulpd nb104_dyMH2(%esp),%xmm1
        mulpd nb104_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixM(%esp),%xmm0
        addpd nb104_fiyM(%esp),%xmm1
        addpd nb104_fizM(%esp),%xmm2
        movapd %xmm3,nb104_fjxH2(%esp)
        movapd %xmm4,nb104_fjyH2(%esp)
        movapd %xmm5,nb104_fjzH2(%esp)
        movapd %xmm0,nb104_fixM(%esp)
        movapd %xmm1,nb104_fiyM(%esp)
        movapd %xmm2,nb104_fizM(%esp)

        ## H1-M interaction 
        movapd nb104_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqMH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1M 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxM(%esp),%xmm3
        movapd nb104_fjyM(%esp),%xmm4
        movapd nb104_fjzM(%esp),%xmm5
        mulpd nb104_dxH1M(%esp),%xmm0
        mulpd nb104_dyH1M(%esp),%xmm1
        mulpd nb104_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH1(%esp),%xmm0
        addpd nb104_fiyH1(%esp),%xmm1
        addpd nb104_fizH1(%esp),%xmm2
        movapd %xmm3,nb104_fjxM(%esp)
        movapd %xmm4,nb104_fjyM(%esp)
        movapd %xmm5,nb104_fjzM(%esp)
        movapd %xmm0,nb104_fixH1(%esp)
        movapd %xmm1,nb104_fiyH1(%esp)
        movapd %xmm2,nb104_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb104_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH1H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH1(%esp),%xmm3
        movapd nb104_fjyH1(%esp),%xmm4
        movapd nb104_fjzH1(%esp),%xmm5
        mulpd nb104_dxH1H1(%esp),%xmm0
        mulpd nb104_dyH1H1(%esp),%xmm1
        mulpd nb104_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH1(%esp),%xmm0
        addpd nb104_fiyH1(%esp),%xmm1
        addpd nb104_fizH1(%esp),%xmm2
        movapd %xmm3,nb104_fjxH1(%esp)
        movapd %xmm4,nb104_fjyH1(%esp)
        movapd %xmm5,nb104_fjzH1(%esp)
        movapd %xmm0,nb104_fixH1(%esp)
        movapd %xmm1,nb104_fiyH1(%esp)
        movapd %xmm2,nb104_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb104_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsMH2  
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH2(%esp),%xmm3
        movapd nb104_fjyH2(%esp),%xmm4
        movapd nb104_fjzH2(%esp),%xmm5
        mulpd nb104_dxH1H2(%esp),%xmm0
        mulpd nb104_dyH1H2(%esp),%xmm1
        mulpd nb104_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH1(%esp),%xmm0
        addpd nb104_fiyH1(%esp),%xmm1
        addpd nb104_fizH1(%esp),%xmm2
        movapd %xmm3,nb104_fjxH2(%esp)
        movapd %xmm4,nb104_fjyH2(%esp)
        movapd %xmm5,nb104_fjzH2(%esp)
        movapd %xmm0,nb104_fixH1(%esp)
        movapd %xmm1,nb104_fiyH1(%esp)
        movapd %xmm2,nb104_fizH1(%esp)

        ## H2-M interaction 
        movapd nb104_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqMH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2M 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxM(%esp),%xmm3
        movapd nb104_fjyM(%esp),%xmm4
        movapd nb104_fjzM(%esp),%xmm5
        mulpd nb104_dxH2M(%esp),%xmm0
        mulpd nb104_dyH2M(%esp),%xmm1
        mulpd nb104_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH2(%esp),%xmm0
        addpd nb104_fiyH2(%esp),%xmm1
        addpd nb104_fizH2(%esp),%xmm2
        movapd %xmm3,nb104_fjxM(%esp)
        movapd %xmm4,nb104_fjyM(%esp)
        movapd %xmm5,nb104_fjzM(%esp)
        movapd %xmm0,nb104_fixH2(%esp)
        movapd %xmm1,nb104_fiyH2(%esp)
        movapd %xmm2,nb104_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb104_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H1 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH1(%esp),%xmm3
        movapd nb104_fjyH1(%esp),%xmm4
        movapd nb104_fjzH1(%esp),%xmm5
        mulpd nb104_dxH2H1(%esp),%xmm0
        mulpd nb104_dyH2H1(%esp),%xmm1
        mulpd nb104_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH2(%esp),%xmm0
        addpd nb104_fiyH2(%esp),%xmm1
        addpd nb104_fizH2(%esp),%xmm2
        movapd %xmm3,nb104_fjxH1(%esp)
        movapd %xmm4,nb104_fjyH1(%esp)
        movapd %xmm5,nb104_fjzH1(%esp)
        movapd %xmm0,nb104_fixH2(%esp)
        movapd %xmm1,nb104_fiyH2(%esp)
        movapd %xmm2,nb104_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb104_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd %xmm0,%xmm0
        mulpd nb104_qqHH(%esp),%xmm1
        mulpd %xmm1,%xmm0       ## fsH2H2 
        addpd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm7,nb104_vctot(%esp)
        movapd %xmm0,%xmm2
        movapd nb104_fjxH2(%esp),%xmm3
        movapd nb104_fjyH2(%esp),%xmm4
        movapd nb104_fjzH2(%esp),%xmm5
        mulpd nb104_dxH2H2(%esp),%xmm0
        mulpd nb104_dyH2H2(%esp),%xmm1
        mulpd nb104_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb104_fixH2(%esp),%xmm0
        addpd nb104_fiyH2(%esp),%xmm1
        addpd nb104_fizH2(%esp),%xmm2
        movapd %xmm3,nb104_fjxH2(%esp)
        movapd %xmm4,nb104_fjyH2(%esp)
        movapd %xmm5,nb104_fjzH2(%esp)
        movapd %xmm0,nb104_fixH2(%esp)
        movapd %xmm1,nb104_fiyH2(%esp)
        movapd %xmm2,nb104_fizH2(%esp)

        movl nb104_faction(%ebp),%edi

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
        addpd nb104_fjxH1(%esp),%xmm0
        addpd nb104_fjyH1(%esp),%xmm1
        addpd nb104_fjzH1(%esp),%xmm2
        addpd nb104_fjxH2(%esp),%xmm3
        addpd nb104_fjyH2(%esp),%xmm4
        addpd nb104_fjzH2(%esp),%xmm5
        addpd nb104_fjxM(%esp),%xmm6
        addpd nb104_fjyM(%esp),%xmm7
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
        addpd nb104_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)
        movhpd %xmm0,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb104_innerk(%esp)
        jl    _nb_kernel104_ia32_sse2.nb104_checksingle
        jmp   _nb_kernel104_ia32_sse2.nb104_unroll_loop
_nb_kernel104_ia32_sse2.nb104_checksingle: 
        movl  nb104_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel104_ia32_sse2.nb104_dosingle
        jmp   _nb_kernel104_ia32_sse2.nb104_updateouterdata
_nb_kernel104_ia32_sse2.nb104_dosingle: 
        movl  nb104_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb104_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb104_jxH1(%esp)
        movapd  %xmm3,nb104_jyH1(%esp)
        movapd  %xmm4,nb104_jzH1(%esp)
        movapd  %xmm5,nb104_jxH2(%esp)
        movapd  %xmm6,nb104_jyH2(%esp)
        movapd  %xmm7,nb104_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb104_jxM(%esp)
        movapd  %xmm3,nb104_jyM(%esp)
        movapd  %xmm4,nb104_jzM(%esp)

        movapd nb104_ixM(%esp),%xmm0
        movapd nb104_iyM(%esp),%xmm1
        movapd nb104_izM(%esp),%xmm2
        movapd nb104_ixM(%esp),%xmm3
        movapd nb104_iyM(%esp),%xmm4
        movapd nb104_izM(%esp),%xmm5
        subsd  nb104_jxM(%esp),%xmm0
        subsd  nb104_jyM(%esp),%xmm1
        subsd  nb104_jzM(%esp),%xmm2
        subsd  nb104_jxH1(%esp),%xmm3
        subsd  nb104_jyH1(%esp),%xmm4
        subsd  nb104_jzH1(%esp),%xmm5
        movapd %xmm0,nb104_dxMM(%esp)
        movapd %xmm1,nb104_dyMM(%esp)
        movapd %xmm2,nb104_dzMM(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxMH1(%esp)
        movapd %xmm4,nb104_dyMH1(%esp)
        movapd %xmm5,nb104_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqMM(%esp)
        movapd %xmm3,nb104_rsqMH1(%esp)

        movapd nb104_ixM(%esp),%xmm0
        movapd nb104_iyM(%esp),%xmm1
        movapd nb104_izM(%esp),%xmm2
        movapd nb104_ixH1(%esp),%xmm3
        movapd nb104_iyH1(%esp),%xmm4
        movapd nb104_izH1(%esp),%xmm5
        subsd  nb104_jxH2(%esp),%xmm0
        subsd  nb104_jyH2(%esp),%xmm1
        subsd  nb104_jzH2(%esp),%xmm2
        subsd  nb104_jxM(%esp),%xmm3
        subsd  nb104_jyM(%esp),%xmm4
        subsd  nb104_jzM(%esp),%xmm5
        movapd %xmm0,nb104_dxMH2(%esp)
        movapd %xmm1,nb104_dyMH2(%esp)
        movapd %xmm2,nb104_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH1M(%esp)
        movapd %xmm4,nb104_dyH1M(%esp)
        movapd %xmm5,nb104_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqMH2(%esp)
        movapd %xmm3,nb104_rsqH1M(%esp)

        movapd nb104_ixH1(%esp),%xmm0
        movapd nb104_iyH1(%esp),%xmm1
        movapd nb104_izH1(%esp),%xmm2
        movapd nb104_ixH1(%esp),%xmm3
        movapd nb104_iyH1(%esp),%xmm4
        movapd nb104_izH1(%esp),%xmm5
        subsd  nb104_jxH1(%esp),%xmm0
        subsd  nb104_jyH1(%esp),%xmm1
        subsd  nb104_jzH1(%esp),%xmm2
        subsd  nb104_jxH2(%esp),%xmm3
        subsd  nb104_jyH2(%esp),%xmm4
        subsd  nb104_jzH2(%esp),%xmm5
        movapd %xmm0,nb104_dxH1H1(%esp)
        movapd %xmm1,nb104_dyH1H1(%esp)
        movapd %xmm2,nb104_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH1H2(%esp)
        movapd %xmm4,nb104_dyH1H2(%esp)
        movapd %xmm5,nb104_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb104_rsqH1H1(%esp)
        movapd %xmm3,nb104_rsqH1H2(%esp)

        movapd nb104_ixH2(%esp),%xmm0
        movapd nb104_iyH2(%esp),%xmm1
        movapd nb104_izH2(%esp),%xmm2
        movapd nb104_ixH2(%esp),%xmm3
        movapd nb104_iyH2(%esp),%xmm4
        movapd nb104_izH2(%esp),%xmm5
        subsd  nb104_jxM(%esp),%xmm0
        subsd  nb104_jyM(%esp),%xmm1
        subsd  nb104_jzM(%esp),%xmm2
        subsd  nb104_jxH1(%esp),%xmm3
        subsd  nb104_jyH1(%esp),%xmm4
        subsd  nb104_jzH1(%esp),%xmm5
        movapd %xmm0,nb104_dxH2M(%esp)
        movapd %xmm1,nb104_dyH2M(%esp)
        movapd %xmm2,nb104_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb104_dxH2H1(%esp)
        movapd %xmm4,nb104_dyH2H1(%esp)
        movapd %xmm5,nb104_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb104_rsqH2M(%esp)
        movapd %xmm4,nb104_rsqH2H1(%esp)

        movapd nb104_ixH2(%esp),%xmm0
        movapd nb104_iyH2(%esp),%xmm1
        movapd nb104_izH2(%esp),%xmm2
        subsd  nb104_jxH2(%esp),%xmm0
        subsd  nb104_jyH2(%esp),%xmm1
        subsd  nb104_jzH2(%esp),%xmm2
        movapd %xmm0,nb104_dxH2H2(%esp)
        movapd %xmm1,nb104_dyH2H2(%esp)
        movapd %xmm2,nb104_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb104_rsqH2H2(%esp)

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
        movapd  nb104_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm3   ## iter1 
        mulsd   nb104_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm1   ## rinv 
        mulsd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvH2H2(%esp)
        movapd %xmm5,nb104_rinvH2H1(%esp)

        movapd nb104_rsqMM(%esp),%xmm0
        movapd nb104_rsqMH1(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb104_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm1   ## rinv 
        mulsd   nb104_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb104_rinvMM(%esp)
        movapd %xmm5,nb104_rinvMH1(%esp)

        movapd nb104_rsqMH2(%esp),%xmm0
        movapd nb104_rsqH1M(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm3   ## iter1 
        mulsd   nb104_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm1   ## rinv 
        mulsd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvMH2(%esp)
        movapd %xmm5,nb104_rinvH1M(%esp)

        movapd nb104_rsqH1H1(%esp),%xmm0
        movapd nb104_rsqH1H2(%esp),%xmm4
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
        movapd  nb104_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm3   ## iter1a 
        mulsd   nb104_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104_half(%esp),%xmm1   ## rinv 
        mulsd   nb104_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104_rinvH1H1(%esp)
        movapd %xmm5,nb104_rinvH1H2(%esp)

        movapd nb104_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb104_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb104_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7
        mulsd  %xmm0,%xmm0
        mulsd  nb104_qqMM(%esp),%xmm7
        mulsd  %xmm7,%xmm0
        addsd  nb104_vctot(%esp),%xmm7
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb104_dxMM(%esp),%xmm0
        mulsd nb104_dyMM(%esp),%xmm1
        mulsd nb104_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixM(%esp),%xmm0
        addsd nb104_fiyM(%esp),%xmm1
        addsd nb104_fizM(%esp),%xmm2
        movlpd %xmm3,nb104_fjxM(%esp)
        movlpd %xmm4,nb104_fjyM(%esp)
        movlpd %xmm5,nb104_fjzM(%esp)
        movlpd %xmm0,nb104_fixM(%esp)
        movlpd %xmm1,nb104_fiyM(%esp)
        movlpd %xmm2,nb104_fizM(%esp)

        ## M-H1 interaction 
        movapd nb104_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqMH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsMH1  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb104_dxMH1(%esp),%xmm0
        mulsd nb104_dyMH1(%esp),%xmm1
        mulsd nb104_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixM(%esp),%xmm0
        addsd nb104_fiyM(%esp),%xmm1
        addsd nb104_fizM(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH1(%esp)
        movlpd %xmm4,nb104_fjyH1(%esp)
        movlpd %xmm5,nb104_fjzH1(%esp)
        movlpd %xmm0,nb104_fixM(%esp)
        movlpd %xmm1,nb104_fiyM(%esp)
        movlpd %xmm2,nb104_fizM(%esp)

        ## M-H2 interaction  
        movapd nb104_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqMH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsMH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb104_dxMH2(%esp),%xmm0
        mulsd nb104_dyMH2(%esp),%xmm1
        mulsd nb104_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixM(%esp),%xmm0
        addsd nb104_fiyM(%esp),%xmm1
        addsd nb104_fizM(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH2(%esp)
        movlpd %xmm4,nb104_fjyH2(%esp)
        movlpd %xmm5,nb104_fjzH2(%esp)
        movlpd %xmm0,nb104_fixM(%esp)
        movlpd %xmm1,nb104_fiyM(%esp)
        movlpd %xmm2,nb104_fizM(%esp)

        ## H1-M interaction 
        movapd nb104_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqMH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1M 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxM(%esp),%xmm3
        movapd nb104_fjyM(%esp),%xmm4
        movapd nb104_fjzM(%esp),%xmm5
        mulsd nb104_dxH1M(%esp),%xmm0
        mulsd nb104_dyH1M(%esp),%xmm1
        mulsd nb104_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH1(%esp),%xmm0
        addsd nb104_fiyH1(%esp),%xmm1
        addsd nb104_fizH1(%esp),%xmm2
        movlpd %xmm3,nb104_fjxM(%esp)
        movlpd %xmm4,nb104_fjyM(%esp)
        movlpd %xmm5,nb104_fjzM(%esp)
        movlpd %xmm0,nb104_fixH1(%esp)
        movlpd %xmm1,nb104_fiyH1(%esp)
        movlpd %xmm2,nb104_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb104_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH1H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH1(%esp),%xmm3
        movapd nb104_fjyH1(%esp),%xmm4
        movapd nb104_fjzH1(%esp),%xmm5
        mulsd nb104_dxH1H1(%esp),%xmm0
        mulsd nb104_dyH1H1(%esp),%xmm1
        mulsd nb104_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH1(%esp),%xmm0
        addsd nb104_fiyH1(%esp),%xmm1
        addsd nb104_fizH1(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH1(%esp)
        movlpd %xmm4,nb104_fjyH1(%esp)
        movlpd %xmm5,nb104_fjzH1(%esp)
        movlpd %xmm0,nb104_fixH1(%esp)
        movlpd %xmm1,nb104_fiyH1(%esp)
        movlpd %xmm2,nb104_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb104_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsMH2  
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH2(%esp),%xmm3
        movapd nb104_fjyH2(%esp),%xmm4
        movapd nb104_fjzH2(%esp),%xmm5
        mulsd nb104_dxH1H2(%esp),%xmm0
        mulsd nb104_dyH1H2(%esp),%xmm1
        mulsd nb104_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH1(%esp),%xmm0
        addsd nb104_fiyH1(%esp),%xmm1
        addsd nb104_fizH1(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH2(%esp)
        movlpd %xmm4,nb104_fjyH2(%esp)
        movlpd %xmm5,nb104_fjzH2(%esp)
        movlpd %xmm0,nb104_fixH1(%esp)
        movlpd %xmm1,nb104_fiyH1(%esp)
        movlpd %xmm2,nb104_fizH1(%esp)

        ## H2-M interaction 
        movapd nb104_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqMH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2M 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxM(%esp),%xmm3
        movapd nb104_fjyM(%esp),%xmm4
        movapd nb104_fjzM(%esp),%xmm5
        mulsd nb104_dxH2M(%esp),%xmm0
        mulsd nb104_dyH2M(%esp),%xmm1
        mulsd nb104_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH2(%esp),%xmm0
        addsd nb104_fiyH2(%esp),%xmm1
        addsd nb104_fizH2(%esp),%xmm2
        movlpd %xmm3,nb104_fjxM(%esp)
        movlpd %xmm4,nb104_fjyM(%esp)
        movlpd %xmm5,nb104_fjzM(%esp)
        movlpd %xmm0,nb104_fixH2(%esp)
        movlpd %xmm1,nb104_fiyH2(%esp)
        movlpd %xmm2,nb104_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb104_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H1 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2
        movapd nb104_fjxH1(%esp),%xmm3
        movapd nb104_fjyH1(%esp),%xmm4
        movapd nb104_fjzH1(%esp),%xmm5
        mulsd nb104_dxH2H1(%esp),%xmm0
        mulsd nb104_dyH2H1(%esp),%xmm1
        mulsd nb104_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH2(%esp),%xmm0
        addsd nb104_fiyH2(%esp),%xmm1
        addsd nb104_fizH2(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH1(%esp)
        movlpd %xmm4,nb104_fjyH1(%esp)
        movlpd %xmm5,nb104_fjzH1(%esp)
        movlpd %xmm0,nb104_fixH2(%esp)
        movlpd %xmm1,nb104_fiyH2(%esp)
        movlpd %xmm2,nb104_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb104_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd %xmm0,%xmm0
        mulsd nb104_qqHH(%esp),%xmm1
        mulsd %xmm1,%xmm0       ## fsH2H2 
        addsd %xmm1,%xmm7       ## add to local vctot 
        movapd %xmm0,%xmm1
        movsd %xmm7,nb104_vctot(%esp)
        movapd %xmm0,%xmm2
        movapd nb104_fjxH2(%esp),%xmm3
        movapd nb104_fjyH2(%esp),%xmm4
        movapd nb104_fjzH2(%esp),%xmm5
        mulsd nb104_dxH2H2(%esp),%xmm0
        mulsd nb104_dyH2H2(%esp),%xmm1
        mulsd nb104_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb104_fixH2(%esp),%xmm0
        addsd nb104_fiyH2(%esp),%xmm1
        addsd nb104_fizH2(%esp),%xmm2
        movlpd %xmm3,nb104_fjxH2(%esp)
        movlpd %xmm4,nb104_fjyH2(%esp)
        movlpd %xmm5,nb104_fjzH2(%esp)
        movlpd %xmm0,nb104_fixH2(%esp)
        movlpd %xmm1,nb104_fiyH2(%esp)
        movlpd %xmm2,nb104_fizH2(%esp)

        movl nb104_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        movlpd 24(%edi,%eax,8),%xmm0
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 40(%edi,%eax,8),%xmm2
        movlpd 48(%edi,%eax,8),%xmm3
        movlpd 56(%edi,%eax,8),%xmm4
        movlpd 64(%edi,%eax,8),%xmm5
        movlpd 72(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%eax,8),%xmm7
        addsd nb104_fjxH1(%esp),%xmm0
        addsd nb104_fjyH1(%esp),%xmm1
        addsd nb104_fjzH1(%esp),%xmm2
        addsd nb104_fjxH2(%esp),%xmm3
        addsd nb104_fjyH2(%esp),%xmm4
        addsd nb104_fjzH2(%esp),%xmm5
        addsd nb104_fjxM(%esp),%xmm6
        addsd nb104_fjyM(%esp),%xmm7
        movlpd %xmm0,24(%edi,%eax,8)
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movlpd %xmm4,56(%edi,%eax,8)
        movlpd %xmm5,64(%edi,%eax,8)
        movlpd %xmm6,72(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%eax,8)

        movlpd 88(%edi,%eax,8),%xmm0
        addsd nb104_fjzM(%esp),%xmm0
        movlpd %xmm0,88(%edi,%eax,8)

_nb_kernel104_ia32_sse2.nb104_updateouterdata: 
        movl  nb104_ii3(%esp),%ecx
        movl  nb104_faction(%ebp),%edi
        movl  nb104_fshift(%ebp),%esi
        movl  nb104_is3(%esp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb104_fixH1(%esp),%xmm0
        movapd nb104_fiyH1(%esp),%xmm1
        movapd nb104_fizH1(%esp),%xmm2

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
        movapd nb104_fixH2(%esp),%xmm0
        movapd nb104_fiyH2(%esp),%xmm1
        movapd nb104_fizH2(%esp),%xmm2

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
        movapd nb104_fixM(%esp),%xmm0
        movapd nb104_fiyM(%esp),%xmm1
        movapd nb104_fizM(%esp),%xmm2

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
        movl nb104_n(%esp),%esi
        ## get group index for i particle 
        movl  nb104_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb104_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb104_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb104_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104_ia32_sse2.nb104_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104_n(%esp)
        jmp _nb_kernel104_ia32_sse2.nb104_outer
_nb_kernel104_ia32_sse2.nb104_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104_ia32_sse2.nb104_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104_ia32_sse2.nb104_threadloop
_nb_kernel104_ia32_sse2.nb104_end: 
        emms

        movl nb104_nouter(%esp),%eax
        movl nb104_ninner(%esp),%ebx
        movl nb104_outeriter(%ebp),%ecx
        movl nb104_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb104_salign(%esp),%eax
        addl %eax,%esp
        addl $1432,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



.globl nb_kernel104nf_ia32_sse2
.globl _nb_kernel104nf_ia32_sse2
nb_kernel104nf_ia32_sse2:       
_nb_kernel104nf_ia32_sse2:      
.set nb104nf_p_nri, 8
.set nb104nf_iinr, 12
.set nb104nf_jindex, 16
.set nb104nf_jjnr, 20
.set nb104nf_shift, 24
.set nb104nf_shiftvec, 28
.set nb104nf_fshift, 32
.set nb104nf_gid, 36
.set nb104nf_pos, 40
.set nb104nf_faction, 44
.set nb104nf_charge, 48
.set nb104nf_p_facel, 52
.set nb104nf_argkrf, 56
.set nb104nf_argcrf, 60
.set nb104nf_Vc, 64
.set nb104nf_type, 68
.set nb104nf_p_ntype, 72
.set nb104nf_vdwparam, 76
.set nb104nf_Vvdw, 80
.set nb104nf_p_tabscale, 84
.set nb104nf_VFtab, 88
.set nb104nf_invsqrta, 92
.set nb104nf_dvda, 96
.set nb104nf_p_gbtabscale, 100
.set nb104nf_GBtab, 104
.set nb104nf_p_nthreads, 108
.set nb104nf_count, 112
.set nb104nf_mtx, 116
.set nb104nf_outeriter, 120
.set nb104nf_inneriter, 124
.set nb104nf_work, 128
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
.set nb104nf_innerjjnr, 680
.set nb104nf_innerk, 684
.set nb104nf_n, 688
.set nb104nf_nn1, 692
.set nb104nf_nri, 696
.set nb104nf_nouter, 700
.set nb104nf_ninner, 704
.set nb104nf_salign, 708
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $712,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb104nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb104nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb104nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb104nf_nouter(%esp)
        movl %eax,nb104nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb104nf_half(%esp)
        movl %ebx,nb104nf_half+4(%esp)
        movsd nb104nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb104nf_half(%esp)
        movapd %xmm3,nb104nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb104nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb104nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3     ## qM 
        movsd %xmm3,%xmm4               ## qM 
        movsd 8(%edx,%ebx,8),%xmm5      ## qH 
        movl nb104nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6      ## facel 
        mulsd  %xmm3,%xmm3              ## qM*qM 
        mulsd  %xmm5,%xmm4              ## qM*qH 
        mulsd  %xmm5,%xmm5              ## qH*qH 
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104nf_qqMM(%esp)
        movapd %xmm4,nb104nf_qqMH(%esp)
        movapd %xmm5,nb104nf_qqHH(%esp)

_nb_kernel104nf_ia32_sse2.nb104nf_threadloop: 
        movl  nb104nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel104nf_ia32_sse2.nb104nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel104nf_ia32_sse2.nb104nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb104nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb104nf_n(%esp)
        movl %ebx,nb104nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel104nf_ia32_sse2.nb104nf_outerstart
        jmp _nb_kernel104nf_ia32_sse2.nb104nf_end

_nb_kernel104nf_ia32_sse2.nb104nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb104nf_nouter(%esp),%ebx
        movl %ebx,nb104nf_nouter(%esp)

_nb_kernel104nf_ia32_sse2.nb104nf_outer: 
        movl  nb104nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb104nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb104nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb104nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb104nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm3
        addsd 32(%eax,%ebx,8),%xmm4
        addsd 40(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb104nf_ixH1(%esp)
        movapd %xmm4,nb104nf_iyH1(%esp)
        movapd %xmm5,nb104nf_izH1(%esp)

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
        movapd %xmm0,nb104nf_ixH2(%esp)
        movapd %xmm1,nb104nf_iyH2(%esp)
        movapd %xmm2,nb104nf_izH2(%esp)
        movapd %xmm3,nb104nf_ixM(%esp)
        movapd %xmm4,nb104nf_iyM(%esp)
        movapd %xmm5,nb104nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb104nf_vctot(%esp)

        movl  nb104nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb104nf_pos(%ebp),%esi
        movl  nb104nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb104nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb104nf_ninner(%esp),%ecx
        movl  %ecx,nb104nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb104nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel104nf_ia32_sse2.nb104nf_unroll_loop
        jmp   _nb_kernel104nf_ia32_sse2.nb104nf_checksingle
_nb_kernel104nf_ia32_sse2.nb104nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb104nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb104nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb104nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb104nf_jxH1(%esp)
        movapd  %xmm3,nb104nf_jyH1(%esp)
        movapd  %xmm4,nb104nf_jzH1(%esp)
        movapd  %xmm5,nb104nf_jxH2(%esp)
        movapd  %xmm6,nb104nf_jyH2(%esp)
        movapd  %xmm7,nb104nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%ebx,8),%xmm2
        movhpd 80(%esi,%ebx,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb104nf_jxM(%esp)
        movapd  %xmm3,nb104nf_jyM(%esp)
        movapd  %xmm4,nb104nf_jzM(%esp)

        movapd nb104nf_ixM(%esp),%xmm0
        movapd nb104nf_iyM(%esp),%xmm1
        movapd nb104nf_izM(%esp),%xmm2
        movapd nb104nf_ixM(%esp),%xmm3
        movapd nb104nf_iyM(%esp),%xmm4
        movapd nb104nf_izM(%esp),%xmm5
        subpd  nb104nf_jxM(%esp),%xmm0
        subpd  nb104nf_jyM(%esp),%xmm1
        subpd  nb104nf_jzM(%esp),%xmm2
        subpd  nb104nf_jxH1(%esp),%xmm3
        subpd  nb104nf_jyH1(%esp),%xmm4
        subpd  nb104nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMM(%esp)
        movapd %xmm3,nb104nf_rsqMH1(%esp)

        movapd nb104nf_ixM(%esp),%xmm0
        movapd nb104nf_iyM(%esp),%xmm1
        movapd nb104nf_izM(%esp),%xmm2
        movapd nb104nf_ixH1(%esp),%xmm3
        movapd nb104nf_iyH1(%esp),%xmm4
        movapd nb104nf_izH1(%esp),%xmm5
        subpd  nb104nf_jxH2(%esp),%xmm0
        subpd  nb104nf_jyH2(%esp),%xmm1
        subpd  nb104nf_jzH2(%esp),%xmm2
        subpd  nb104nf_jxM(%esp),%xmm3
        subpd  nb104nf_jyM(%esp),%xmm4
        subpd  nb104nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMH2(%esp)
        movapd %xmm3,nb104nf_rsqH1M(%esp)

        movapd nb104nf_ixH1(%esp),%xmm0
        movapd nb104nf_iyH1(%esp),%xmm1
        movapd nb104nf_izH1(%esp),%xmm2
        movapd nb104nf_ixH1(%esp),%xmm3
        movapd nb104nf_iyH1(%esp),%xmm4
        movapd nb104nf_izH1(%esp),%xmm5
        subpd  nb104nf_jxH1(%esp),%xmm0
        subpd  nb104nf_jyH1(%esp),%xmm1
        subpd  nb104nf_jzH1(%esp),%xmm2
        subpd  nb104nf_jxH2(%esp),%xmm3
        subpd  nb104nf_jyH2(%esp),%xmm4
        subpd  nb104nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH1H1(%esp)
        movapd %xmm3,nb104nf_rsqH1H2(%esp)

        movapd nb104nf_ixH2(%esp),%xmm0
        movapd nb104nf_iyH2(%esp),%xmm1
        movapd nb104nf_izH2(%esp),%xmm2
        movapd nb104nf_ixH2(%esp),%xmm3
        movapd nb104nf_iyH2(%esp),%xmm4
        movapd nb104nf_izH2(%esp),%xmm5
        subpd  nb104nf_jxM(%esp),%xmm0
        subpd  nb104nf_jyM(%esp),%xmm1
        subpd  nb104nf_jzM(%esp),%xmm2
        subpd  nb104nf_jxH1(%esp),%xmm3
        subpd  nb104nf_jyH1(%esp),%xmm4
        subpd  nb104nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH2M(%esp)
        movapd %xmm4,nb104nf_rsqH2H1(%esp)

        movapd nb104nf_ixH2(%esp),%xmm0
        movapd nb104nf_iyH2(%esp),%xmm1
        movapd nb104nf_izH2(%esp),%xmm2
        subpd  nb104nf_jxH2(%esp),%xmm0
        subpd  nb104nf_jyH2(%esp),%xmm1
        subpd  nb104nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb104nf_rsqH2H2(%esp)

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
        movapd  nb104nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb104nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH2H2(%esp)
        movapd %xmm5,nb104nf_rinvH2H1(%esp)

        movapd nb104nf_rsqMM(%esp),%xmm0
        movapd nb104nf_rsqMH1(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb104nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb104nf_rinvMM(%esp)
        movapd %xmm5,nb104nf_rinvMH1(%esp)

        movapd nb104nf_rsqMH2(%esp),%xmm0
        movapd nb104nf_rsqH1M(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb104nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvMH2(%esp)
        movapd %xmm5,nb104nf_rinvH1M(%esp)

        movapd nb104nf_rsqH1H1(%esp),%xmm0
        movapd nb104nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb104nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH1H1(%esp)
        movapd %xmm5,nb104nf_rinvH1H2(%esp)

        movapd nb104nf_rsqH2M(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb104nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb104nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb104nf_rinvMM(%esp),%xmm0
        mulpd  nb104nf_qqMM(%esp),%xmm0
        addpd  nb104nf_vctot(%esp),%xmm0

        ## other interactions 
        movapd nb104nf_rinvMH1(%esp),%xmm1
        movapd nb104nf_rinvH1H1(%esp),%xmm2

        addpd nb104nf_rinvMH2(%esp),%xmm1
        addpd nb104nf_rinvH1H2(%esp),%xmm2

        addpd nb104nf_rinvH1M(%esp),%xmm1
        addpd nb104nf_rinvH2H1(%esp),%xmm2

        addpd nb104nf_rinvH2M(%esp),%xmm1
        addpd nb104nf_rinvH2H2(%esp),%xmm2

        mulpd nb104nf_qqMH(%esp),%xmm1
        mulpd nb104nf_qqHH(%esp),%xmm2

        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0

        movapd %xmm0,nb104nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb104nf_innerk(%esp)
        jl    _nb_kernel104nf_ia32_sse2.nb104nf_checksingle
        jmp   _nb_kernel104nf_ia32_sse2.nb104nf_unroll_loop
_nb_kernel104nf_ia32_sse2.nb104nf_checksingle: 
        movl  nb104nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel104nf_ia32_sse2.nb104nf_dosingle
        jmp   _nb_kernel104nf_ia32_sse2.nb104nf_updateouterdata
_nb_kernel104nf_ia32_sse2.nb104nf_dosingle: 
        movl  nb104nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb104nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## move j coordinates to local temp variables 
        movlpd 24(%esi,%eax,8),%xmm2
        movlpd 32(%esi,%eax,8),%xmm3
        movlpd 40(%esi,%eax,8),%xmm4
        movlpd 48(%esi,%eax,8),%xmm5
        movlpd 56(%esi,%eax,8),%xmm6
        movlpd 64(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb104nf_jxH1(%esp)
        movapd  %xmm3,nb104nf_jyH1(%esp)
        movapd  %xmm4,nb104nf_jzH1(%esp)
        movapd  %xmm5,nb104nf_jxH2(%esp)
        movapd  %xmm6,nb104nf_jyH2(%esp)
        movapd  %xmm7,nb104nf_jzH2(%esp)
        movlpd 72(%esi,%eax,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 88(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb104nf_jxM(%esp)
        movapd  %xmm3,nb104nf_jyM(%esp)
        movapd  %xmm4,nb104nf_jzM(%esp)

        movapd nb104nf_ixM(%esp),%xmm0
        movapd nb104nf_iyM(%esp),%xmm1
        movapd nb104nf_izM(%esp),%xmm2
        movapd nb104nf_ixM(%esp),%xmm3
        movapd nb104nf_iyM(%esp),%xmm4
        movapd nb104nf_izM(%esp),%xmm5
        subsd  nb104nf_jxM(%esp),%xmm0
        subsd  nb104nf_jyM(%esp),%xmm1
        subsd  nb104nf_jzM(%esp),%xmm2
        subsd  nb104nf_jxH1(%esp),%xmm3
        subsd  nb104nf_jyH1(%esp),%xmm4
        subsd  nb104nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMM(%esp)
        movapd %xmm3,nb104nf_rsqMH1(%esp)

        movapd nb104nf_ixM(%esp),%xmm0
        movapd nb104nf_iyM(%esp),%xmm1
        movapd nb104nf_izM(%esp),%xmm2
        movapd nb104nf_ixH1(%esp),%xmm3
        movapd nb104nf_iyH1(%esp),%xmm4
        movapd nb104nf_izH1(%esp),%xmm5
        subsd  nb104nf_jxH2(%esp),%xmm0
        subsd  nb104nf_jyH2(%esp),%xmm1
        subsd  nb104nf_jzH2(%esp),%xmm2
        subsd  nb104nf_jxM(%esp),%xmm3
        subsd  nb104nf_jyM(%esp),%xmm4
        subsd  nb104nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqMH2(%esp)
        movapd %xmm3,nb104nf_rsqH1M(%esp)

        movapd nb104nf_ixH1(%esp),%xmm0
        movapd nb104nf_iyH1(%esp),%xmm1
        movapd nb104nf_izH1(%esp),%xmm2
        movapd nb104nf_ixH1(%esp),%xmm3
        movapd nb104nf_iyH1(%esp),%xmm4
        movapd nb104nf_izH1(%esp),%xmm5
        subsd  nb104nf_jxH1(%esp),%xmm0
        subsd  nb104nf_jyH1(%esp),%xmm1
        subsd  nb104nf_jzH1(%esp),%xmm2
        subsd  nb104nf_jxH2(%esp),%xmm3
        subsd  nb104nf_jyH2(%esp),%xmm4
        subsd  nb104nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH1H1(%esp)
        movapd %xmm3,nb104nf_rsqH1H2(%esp)

        movapd nb104nf_ixH2(%esp),%xmm0
        movapd nb104nf_iyH2(%esp),%xmm1
        movapd nb104nf_izH2(%esp),%xmm2
        movapd nb104nf_ixH2(%esp),%xmm3
        movapd nb104nf_iyH2(%esp),%xmm4
        movapd nb104nf_izH2(%esp),%xmm5
        subsd  nb104nf_jxM(%esp),%xmm0
        subsd  nb104nf_jyM(%esp),%xmm1
        subsd  nb104nf_jzM(%esp),%xmm2
        subsd  nb104nf_jxH1(%esp),%xmm3
        subsd  nb104nf_jyH1(%esp),%xmm4
        subsd  nb104nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb104nf_rsqH2M(%esp)
        movapd %xmm4,nb104nf_rsqH2H1(%esp)

        movapd nb104nf_ixH2(%esp),%xmm0
        movapd nb104nf_iyH2(%esp),%xmm1
        movapd nb104nf_izH2(%esp),%xmm2
        subsd  nb104nf_jxH2(%esp),%xmm0
        subsd  nb104nf_jyH2(%esp),%xmm1
        subsd  nb104nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb104nf_rsqH2H2(%esp)

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
        movapd  nb104nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb104nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH2H2(%esp)
        movapd %xmm5,nb104nf_rinvH2H1(%esp)

        movapd nb104nf_rsqMM(%esp),%xmm0
        movapd nb104nf_rsqMH1(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb104nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb104nf_rinvMM(%esp)
        movapd %xmm5,nb104nf_rinvMH1(%esp)

        movapd nb104nf_rsqMH2(%esp),%xmm0
        movapd nb104nf_rsqH1M(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb104nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvMH2(%esp)
        movapd %xmm5,nb104nf_rinvH1M(%esp)

        movapd nb104nf_rsqH1H1(%esp),%xmm0
        movapd nb104nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb104nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb104nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb104nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb104nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb104nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb104nf_rinvH1H1(%esp)
        movapd %xmm5,nb104nf_rinvH1H2(%esp)

        movapd nb104nf_rsqH2M(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb104nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb104nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb104nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb104nf_rinvH2M(%esp)

        ## start with MM interaction 
        movapd nb104nf_rinvMM(%esp),%xmm0
        mulpd  nb104nf_qqMM(%esp),%xmm0
        addpd  nb104nf_vctot(%esp),%xmm0

        ## other interactions 
        movapd nb104nf_rinvMH1(%esp),%xmm1
        movapd nb104nf_rinvH1H1(%esp),%xmm2

        addsd nb104nf_rinvMH2(%esp),%xmm1
        addsd nb104nf_rinvH1H2(%esp),%xmm2

        addsd nb104nf_rinvH1M(%esp),%xmm1
        addsd nb104nf_rinvH2H1(%esp),%xmm2

        addsd nb104nf_rinvH2M(%esp),%xmm1
        addsd nb104nf_rinvH2H2(%esp),%xmm2

        mulsd nb104nf_qqMH(%esp),%xmm1
        mulsd nb104nf_qqHH(%esp),%xmm2

        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0

        movlpd %xmm0,nb104nf_vctot(%esp)

_nb_kernel104nf_ia32_sse2.nb104nf_updateouterdata: 
        ## get n from stack
        movl nb104nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb104nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        movapd nb104nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb104nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb104nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel104nf_ia32_sse2.nb104nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb104nf_n(%esp)
        jmp _nb_kernel104nf_ia32_sse2.nb104nf_outer
_nb_kernel104nf_ia32_sse2.nb104nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb104nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel104nf_ia32_sse2.nb104nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel104nf_ia32_sse2.nb104nf_threadloop
_nb_kernel104nf_ia32_sse2.nb104nf_end: 
        emms

        movl nb104nf_nouter(%esp),%eax
        movl nb104nf_ninner(%esp),%ebx
        movl nb104nf_outeriter(%ebp),%ecx
        movl nb104nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb104nf_salign(%esp),%eax
        addl %eax,%esp
        addl $712,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



