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


.globl nb_kernel234_ia32_sse
.globl _nb_kernel234_ia32_sse
nb_kernel234_ia32_sse:  
_nb_kernel234_ia32_sse: 
.set nb234_p_nri, 8
.set nb234_iinr, 12
.set nb234_jindex, 16
.set nb234_jjnr, 20
.set nb234_shift, 24
.set nb234_shiftvec, 28
.set nb234_fshift, 32
.set nb234_gid, 36
.set nb234_pos, 40
.set nb234_faction, 44
.set nb234_charge, 48
.set nb234_p_facel, 52
.set nb234_argkrf, 56
.set nb234_argcrf, 60
.set nb234_Vc, 64
.set nb234_type, 68
.set nb234_p_ntype, 72
.set nb234_vdwparam, 76
.set nb234_Vvdw, 80
.set nb234_p_tabscale, 84
.set nb234_VFtab, 88
.set nb234_invsqrta, 92
.set nb234_dvda, 96
.set nb234_p_gbtabscale, 100
.set nb234_GBtab, 104
.set nb234_p_nthreads, 108
.set nb234_count, 112
.set nb234_mtx, 116
.set nb234_outeriter, 120
.set nb234_inneriter, 124
.set nb234_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb234_ixO, 0
.set nb234_iyO, 16
.set nb234_izO, 32
.set nb234_ixH1, 48
.set nb234_iyH1, 64
.set nb234_izH1, 80
.set nb234_ixH2, 96
.set nb234_iyH2, 112
.set nb234_izH2, 128
.set nb234_ixM, 144
.set nb234_iyM, 160
.set nb234_izM, 176
.set nb234_jxO, 192
.set nb234_jyO, 208
.set nb234_jzO, 224
.set nb234_jxH1, 240
.set nb234_jyH1, 256
.set nb234_jzH1, 272
.set nb234_jxH2, 288
.set nb234_jyH2, 304
.set nb234_jzH2, 320
.set nb234_jxM, 336
.set nb234_jyM, 352
.set nb234_jzM, 368
.set nb234_dxOO, 384
.set nb234_dyOO, 400
.set nb234_dzOO, 416
.set nb234_dxH1H1, 432
.set nb234_dyH1H1, 448
.set nb234_dzH1H1, 464
.set nb234_dxH1H2, 480
.set nb234_dyH1H2, 496
.set nb234_dzH1H2, 512
.set nb234_dxH1M, 528
.set nb234_dyH1M, 544
.set nb234_dzH1M, 560
.set nb234_dxH2H1, 576
.set nb234_dyH2H1, 592
.set nb234_dzH2H1, 608
.set nb234_dxH2H2, 624
.set nb234_dyH2H2, 640
.set nb234_dzH2H2, 656
.set nb234_dxH2M, 672
.set nb234_dyH2M, 688
.set nb234_dzH2M, 704
.set nb234_dxMH1, 720
.set nb234_dyMH1, 736
.set nb234_dzMH1, 752
.set nb234_dxMH2, 768
.set nb234_dyMH2, 784
.set nb234_dzMH2, 800
.set nb234_dxMM, 816
.set nb234_dyMM, 832
.set nb234_dzMM, 848
.set nb234_qqMM, 864
.set nb234_qqMH, 880
.set nb234_qqHH, 896
.set nb234_two, 912
.set nb234_c6, 928
.set nb234_c12, 944
.set nb234_tsc, 960
.set nb234_fscal, 976
.set nb234_vctot, 992
.set nb234_Vvdwtot, 1008
.set nb234_fixO, 1024
.set nb234_fiyO, 1040
.set nb234_fizO, 1056
.set nb234_fixH1, 1072
.set nb234_fiyH1, 1088
.set nb234_fizH1, 1104
.set nb234_fixH2, 1120
.set nb234_fiyH2, 1136
.set nb234_fizH2, 1152
.set nb234_fixM, 1168
.set nb234_fiyM, 1184
.set nb234_fizM, 1200
.set nb234_fjxO, 1216
.set nb234_fjyO, 1232
.set nb234_fjzO, 1248
.set nb234_fjxH1, 1264
.set nb234_fjyH1, 1280
.set nb234_fjzH1, 1296
.set nb234_fjxH2, 1312
.set nb234_fjyH2, 1328
.set nb234_fjzH2, 1344
.set nb234_fjxM, 1360
.set nb234_fjyM, 1376
.set nb234_fjzM, 1392
.set nb234_half, 1408
.set nb234_three, 1424
.set nb234_rsqOO, 1440
.set nb234_rsqH1H1, 1456
.set nb234_rsqH1H2, 1472
.set nb234_rsqH1M, 1488
.set nb234_rsqH2H1, 1504
.set nb234_rsqH2H2, 1520
.set nb234_rsqH2M, 1536
.set nb234_rsqMH1, 1552
.set nb234_rsqMH2, 1568
.set nb234_rsqMM, 1584
.set nb234_rinvOO, 1600
.set nb234_rinvH1H1, 1616
.set nb234_rinvH1H2, 1632
.set nb234_rinvH1M, 1648
.set nb234_rinvH2H1, 1664
.set nb234_rinvH2H2, 1680
.set nb234_rinvH2M, 1696
.set nb234_rinvMH1, 1712
.set nb234_rinvMH2, 1728
.set nb234_rinvMM, 1744
.set nb234_fstmp, 1760
.set nb234_krf, 1776
.set nb234_crf, 1792
.set nb234_is3, 1808
.set nb234_ii3, 1812
.set nb234_innerjjnr, 1816
.set nb234_innerk, 1820
.set nb234_n, 1824
.set nb234_nn1, 1828
.set nb234_nri, 1832
.set nb234_nouter, 1836
.set nb234_ninner, 1840
.set nb234_salign, 1844
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
        movl %eax,nb234_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb234_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb234_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb234_nouter(%esp)
        movl %eax,nb234_ninner(%esp)

        movl nb234_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb234_tsc(%esp)

        movl nb234_argkrf(%ebp),%esi
        movl nb234_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb234_krf(%esp)
        movaps %xmm6,nb234_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb234_half(%esp)
        movss nb234_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb234_half(%esp)
        movaps %xmm2,nb234_two(%esp)
        movaps %xmm3,nb234_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb234_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb234_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb234_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb234_qqMM(%esp)
        movaps %xmm4,nb234_qqMH(%esp)
        movaps %xmm5,nb234_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb234_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb234_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb234_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb234_c6(%esp)
        movaps %xmm1,nb234_c12(%esp)

_nb_kernel234_ia32_sse.nb234_threadloop: 
        movl  nb234_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel234_ia32_sse.nb234_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel234_ia32_sse.nb234_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb234_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb234_n(%esp)
        movl %ebx,nb234_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel234_ia32_sse.nb234_outerstart
        jmp _nb_kernel234_ia32_sse.nb234_end

_nb_kernel234_ia32_sse.nb234_outerstart: 
        ## ebx contains number of outer iterations
        addl nb234_nouter(%esp),%ebx
        movl %ebx,nb234_nouter(%esp)

_nb_kernel234_ia32_sse.nb234_outer: 
        movl  nb234_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb234_is3(%esp)      ## store is3 

        movl  nb234_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb234_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb234_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb234_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb234_ixO(%esp)
        movaps %xmm4,nb234_iyO(%esp)
        movaps %xmm5,nb234_izO(%esp)
        movaps %xmm6,nb234_ixH1(%esp)
        movaps %xmm7,nb234_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb234_izH1(%esp)
        movaps %xmm0,nb234_ixH2(%esp)
        movaps %xmm1,nb234_iyH2(%esp)
        movaps %xmm2,nb234_izH2(%esp)
        movaps %xmm3,nb234_ixM(%esp)
        movaps %xmm4,nb234_iyM(%esp)
        movaps %xmm5,nb234_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb234_vctot(%esp)
        movaps %xmm4,nb234_Vvdwtot(%esp)
        movaps %xmm4,nb234_fixO(%esp)
        movaps %xmm4,nb234_fiyO(%esp)
        movaps %xmm4,nb234_fizO(%esp)
        movaps %xmm4,nb234_fixH1(%esp)
        movaps %xmm4,nb234_fiyH1(%esp)
        movaps %xmm4,nb234_fizH1(%esp)
        movaps %xmm4,nb234_fixH2(%esp)
        movaps %xmm4,nb234_fiyH2(%esp)
        movaps %xmm4,nb234_fizH2(%esp)
        movaps %xmm4,nb234_fixM(%esp)
        movaps %xmm4,nb234_fiyM(%esp)
        movaps %xmm4,nb234_fizM(%esp)

        movl  nb234_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb234_pos(%ebp),%esi
        movl  nb234_faction(%ebp),%edi
        movl  nb234_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb234_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb234_ninner(%esp),%ecx
        movl  %ecx,nb234_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb234_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel234_ia32_sse.nb234_unroll_loop
        jmp   _nb_kernel234_ia32_sse.nb234_single_check
_nb_kernel234_ia32_sse.nb234_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb234_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb234_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb234_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234_jxO(%esp)
        movaps %xmm1,nb234_jyO(%esp)
        movaps %xmm2,nb234_jzO(%esp)
        movaps %xmm3,nb234_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234_jyH1(%esp)
        movaps %xmm1,nb234_jzH1(%esp)
        movaps %xmm2,nb234_jxH2(%esp)
        movaps %xmm3,nb234_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234_jzH2(%esp)
        movaps %xmm1,nb234_jxM(%esp)
        movaps %xmm2,nb234_jyM(%esp)
        movaps %xmm3,nb234_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb234_ixO(%esp),%xmm0
        movaps nb234_iyO(%esp),%xmm1
        movaps nb234_izO(%esp),%xmm2
        movaps nb234_ixH1(%esp),%xmm3
        movaps nb234_iyH1(%esp),%xmm4
        movaps nb234_izH1(%esp),%xmm5
        subps  nb234_jxO(%esp),%xmm0
        subps  nb234_jyO(%esp),%xmm1
        subps  nb234_jzO(%esp),%xmm2
        subps  nb234_jxH1(%esp),%xmm3
        subps  nb234_jyH1(%esp),%xmm4
        subps  nb234_jzH1(%esp),%xmm5
        movaps %xmm0,nb234_dxOO(%esp)
        movaps %xmm1,nb234_dyOO(%esp)
        movaps %xmm2,nb234_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb234_dxH1H1(%esp)
        movaps %xmm4,nb234_dyH1H1(%esp)
        movaps %xmm5,nb234_dzH1H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb234_rsqOO(%esp)
        movaps %xmm3,nb234_rsqH1H1(%esp)

        movaps nb234_ixH1(%esp),%xmm0
        movaps nb234_iyH1(%esp),%xmm1
        movaps nb234_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234_jxH2(%esp),%xmm0
        subps  nb234_jyH2(%esp),%xmm1
        subps  nb234_jzH2(%esp),%xmm2
        subps  nb234_jxM(%esp),%xmm3
        subps  nb234_jyM(%esp),%xmm4
        subps  nb234_jzM(%esp),%xmm5
        movaps %xmm0,nb234_dxH1H2(%esp)
        movaps %xmm1,nb234_dyH1H2(%esp)
        movaps %xmm2,nb234_dzH1H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb234_dxH1M(%esp)
        movaps %xmm4,nb234_dyH1M(%esp)
        movaps %xmm5,nb234_dzH1M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb234_rsqH1H2(%esp)
        movaps %xmm3,nb234_rsqH1M(%esp)

        movaps nb234_ixH2(%esp),%xmm0
        movaps nb234_iyH2(%esp),%xmm1
        movaps nb234_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234_jxH1(%esp),%xmm0
        subps  nb234_jyH1(%esp),%xmm1
        subps  nb234_jzH1(%esp),%xmm2
        subps  nb234_jxH2(%esp),%xmm3
        subps  nb234_jyH2(%esp),%xmm4
        subps  nb234_jzH2(%esp),%xmm5
        movaps %xmm0,nb234_dxH2H1(%esp)
        movaps %xmm1,nb234_dyH2H1(%esp)
        movaps %xmm2,nb234_dzH2H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb234_dxH2H2(%esp)
        movaps %xmm4,nb234_dyH2H2(%esp)
        movaps %xmm5,nb234_dzH2H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb234_rsqH2H1(%esp)
        movaps %xmm3,nb234_rsqH2H2(%esp)

        movaps nb234_ixH2(%esp),%xmm0
        movaps nb234_iyH2(%esp),%xmm1
        movaps nb234_izH2(%esp),%xmm2
        movaps nb234_ixM(%esp),%xmm3
        movaps nb234_iyM(%esp),%xmm4
        movaps nb234_izM(%esp),%xmm5
        subps  nb234_jxM(%esp),%xmm0
        subps  nb234_jyM(%esp),%xmm1
        subps  nb234_jzM(%esp),%xmm2
        subps  nb234_jxH1(%esp),%xmm3
        subps  nb234_jyH1(%esp),%xmm4
        subps  nb234_jzH1(%esp),%xmm5
        movaps %xmm0,nb234_dxH2M(%esp)
        movaps %xmm1,nb234_dyH2M(%esp)
        movaps %xmm2,nb234_dzH2M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb234_dxMH1(%esp)
        movaps %xmm4,nb234_dyMH1(%esp)
        movaps %xmm5,nb234_dzMH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb234_rsqH2M(%esp)
        movaps %xmm4,nb234_rsqMH1(%esp)

        movaps nb234_ixM(%esp),%xmm0
        movaps nb234_iyM(%esp),%xmm1
        movaps nb234_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234_jxH2(%esp),%xmm0
        subps  nb234_jyH2(%esp),%xmm1
        subps  nb234_jzH2(%esp),%xmm2
        subps  nb234_jxM(%esp),%xmm3
        subps  nb234_jyM(%esp),%xmm4
        subps  nb234_jzM(%esp),%xmm5
        movaps %xmm0,nb234_dxMH2(%esp)
        movaps %xmm1,nb234_dyMH2(%esp)
        movaps %xmm2,nb234_dzMH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb234_dxMM(%esp)
        movaps %xmm4,nb234_dyMM(%esp)
        movaps %xmm5,nb234_dzMM(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb234_rsqMH2(%esp)
        movaps %xmm4,nb234_rsqMM(%esp)

        ## start by doing invsqrt for OO
        rsqrtps nb234_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234_three(%esp),%xmm3
        mulps   nb234_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb234_half(%esp),%xmm3
        movaps  %xmm3,nb234_rinvOO(%esp)

        ## more invsqrt ops - do two at a time.
        rsqrtps nb234_rsqH1H1(%esp),%xmm1
        rsqrtps nb234_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234_rsqH1H1(%esp),%xmm1
        mulps   nb234_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb234_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb234_rinvH1H1(%esp)
        movaps  %xmm7,nb234_rinvH1H2(%esp)

        rsqrtps nb234_rsqH1M(%esp),%xmm1
        rsqrtps nb234_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234_rsqH1M(%esp),%xmm1
        mulps   nb234_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm3
        mulps   nb234_half(%esp),%xmm7
        movaps  %xmm3,nb234_rinvH1M(%esp)
        movaps  %xmm7,nb234_rinvH2H1(%esp)

        rsqrtps nb234_rsqH2H2(%esp),%xmm1
        rsqrtps nb234_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234_rsqH2H2(%esp),%xmm1
        mulps   nb234_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm3
        mulps   nb234_half(%esp),%xmm7
        movaps  %xmm3,nb234_rinvH2H2(%esp)
        movaps  %xmm7,nb234_rinvH2M(%esp)

        rsqrtps nb234_rsqMH1(%esp),%xmm1
        rsqrtps nb234_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234_rsqMH1(%esp),%xmm1
        mulps   nb234_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm3
        mulps   nb234_half(%esp),%xmm7
        movaps  %xmm3,nb234_rinvMH1(%esp)
        movaps  %xmm7,nb234_rinvMH2(%esp)

        rsqrtps nb234_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234_three(%esp),%xmm3
        mulps   nb234_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb234_half(%esp),%xmm3
        movaps  %xmm3,nb234_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb234_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb234_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb234_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6
    pslld $3,%mm7

    movd %eax,%mm0
    movd %ebx,%mm1
    movd %ecx,%mm2
    movd %edx,%mm3

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movl nb234_VFtab(%ebp),%esi

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  nb234_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb234_c6(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    movaps  %xmm7,nb234_fstmp(%esp)

    addps  nb234_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb234_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  nb234_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb234_c12(%esp),%xmm4
    mulps  %xmm4,%xmm7 ## fijR 
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb234_fstmp(%esp),%xmm7
    mulps nb234_tsc(%esp),%xmm7

    addps  nb234_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb234_Vvdwtot(%esp)

        xorps %xmm2,%xmm2
        mulps  %xmm0,%xmm7
        subps %xmm7,%xmm2

        movaps %xmm2,%xmm0
        movaps %xmm2,%xmm1

    movd  %mm0,%eax
    movd  %mm1,%ebx
    movd  %mm2,%ecx
    movd  %mm3,%edx

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb234_dxOO(%esp),%xmm0
        mulps nb234_dyOO(%esp),%xmm1
        mulps nb234_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixO(%esp),%xmm0
        addps nb234_fiyO(%esp),%xmm1
        addps nb234_fizO(%esp),%xmm2
        movaps %xmm3,nb234_fjxO(%esp)
        movaps %xmm4,nb234_fjyO(%esp)
        movaps %xmm5,nb234_fjzO(%esp)
        movaps %xmm0,nb234_fixO(%esp)
        movaps %xmm1,nb234_fiyO(%esp)
        movaps %xmm2,nb234_fizO(%esp)

        ## Coulomb interactions 
        ## start with H1-H1 interaction 
        movaps nb234_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulps  nb234_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234_crf(%esp),%xmm6
        mulps  nb234_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addps  nb234_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb234_dxH1H1(%esp),%xmm0
        mulps nb234_dyH1H1(%esp),%xmm1
        mulps nb234_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH1(%esp),%xmm0
        addps nb234_fiyH1(%esp),%xmm1
        addps nb234_fizH1(%esp),%xmm2
        movaps %xmm3,nb234_fjxH1(%esp)
        movaps %xmm4,nb234_fjyH1(%esp)
        movaps %xmm5,nb234_fjzH1(%esp)
        movaps %xmm0,nb234_fixH1(%esp)
        movaps %xmm1,nb234_fiyH1(%esp)
        movaps %xmm2,nb234_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb234_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps  %xmm0,%xmm0
        mulps  nb234_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH1  
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb234_dxH1H2(%esp),%xmm0
        mulps nb234_dyH1H2(%esp),%xmm1
        mulps nb234_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH1(%esp),%xmm0
        addps nb234_fiyH1(%esp),%xmm1
        addps nb234_fizH1(%esp),%xmm2
        movaps %xmm3,nb234_fjxH2(%esp)
        movaps %xmm4,nb234_fjyH2(%esp)
        movaps %xmm5,nb234_fjzH2(%esp)
        movaps %xmm0,nb234_fixH1(%esp)
        movaps %xmm1,nb234_fiyH1(%esp)
        movaps %xmm2,nb234_fizH1(%esp)

        ## H1-M interaction  
        movaps nb234_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=Rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb234_dxH1M(%esp),%xmm0
        mulps nb234_dyH1M(%esp),%xmm1
        mulps nb234_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH1(%esp),%xmm0
        addps nb234_fiyH1(%esp),%xmm1
        addps nb234_fizH1(%esp),%xmm2
        movaps %xmm3,nb234_fjxM(%esp)
        movaps %xmm4,nb234_fjyM(%esp)
        movaps %xmm5,nb234_fjzM(%esp)
        movaps %xmm0,nb234_fixH1(%esp)
        movaps %xmm1,nb234_fiyH1(%esp)
        movaps %xmm2,nb234_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb234_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb234_fjxH1(%esp),%xmm3
        movaps nb234_fjyH1(%esp),%xmm4
        movaps nb234_fjzH1(%esp),%xmm5
        mulps nb234_dxH2H1(%esp),%xmm0
        mulps nb234_dyH2H1(%esp),%xmm1
        mulps nb234_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH2(%esp),%xmm0
        addps nb234_fiyH2(%esp),%xmm1
        addps nb234_fizH2(%esp),%xmm2
        movaps %xmm3,nb234_fjxH1(%esp)
        movaps %xmm4,nb234_fjyH1(%esp)
        movaps %xmm5,nb234_fjzH1(%esp)
        movaps %xmm0,nb234_fixH2(%esp)
        movaps %xmm1,nb234_fiyH2(%esp)
        movaps %xmm2,nb234_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb234_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb234_fjxH2(%esp),%xmm3
        movaps nb234_fjyH2(%esp),%xmm4
        movaps nb234_fjzH2(%esp),%xmm5
        mulps nb234_dxH2H2(%esp),%xmm0
        mulps nb234_dyH2H2(%esp),%xmm1
        mulps nb234_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH2(%esp),%xmm0
        addps nb234_fiyH2(%esp),%xmm1
        addps nb234_fizH2(%esp),%xmm2
        movaps %xmm3,nb234_fjxH2(%esp)
        movaps %xmm4,nb234_fjyH2(%esp)
        movaps %xmm5,nb234_fjzH2(%esp)
        movaps %xmm0,nb234_fixH2(%esp)
        movaps %xmm1,nb234_fiyH2(%esp)
        movaps %xmm2,nb234_fizH2(%esp)

        ## H2-M interaction 
        movaps nb234_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb234_fjxM(%esp),%xmm3
        movaps nb234_fjyM(%esp),%xmm4
        movaps nb234_fjzM(%esp),%xmm5
        mulps nb234_dxH2M(%esp),%xmm0
        mulps nb234_dyH2M(%esp),%xmm1
        mulps nb234_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixH2(%esp),%xmm0
        addps nb234_fiyH2(%esp),%xmm1
        addps nb234_fizH2(%esp),%xmm2
        movaps %xmm3,nb234_fjxM(%esp)
        movaps %xmm4,nb234_fjyM(%esp)
        movaps %xmm5,nb234_fjzM(%esp)
        movaps %xmm0,nb234_fixH2(%esp)
        movaps %xmm1,nb234_fiyH2(%esp)
        movaps %xmm2,nb234_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb234_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb234_fjxH1(%esp),%xmm3
        movaps nb234_fjyH1(%esp),%xmm4
        movaps nb234_fjzH1(%esp),%xmm5
        mulps nb234_dxMH1(%esp),%xmm0
        mulps nb234_dyMH1(%esp),%xmm1
        mulps nb234_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixM(%esp),%xmm0
        addps nb234_fiyM(%esp),%xmm1
        addps nb234_fizM(%esp),%xmm2
        movaps %xmm3,nb234_fjxH1(%esp)
        movaps %xmm4,nb234_fjyH1(%esp)
        movaps %xmm5,nb234_fjzH1(%esp)
        movaps %xmm0,nb234_fixM(%esp)
        movaps %xmm1,nb234_fiyM(%esp)
        movaps %xmm2,nb234_fizM(%esp)

        ## M-H2 interaction 
        movaps nb234_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb234_fjxH2(%esp),%xmm3
        movaps nb234_fjyH2(%esp),%xmm4
        movaps nb234_fjzH2(%esp),%xmm5
        mulps nb234_dxMH2(%esp),%xmm0
        mulps nb234_dyMH2(%esp),%xmm1
        mulps nb234_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixM(%esp),%xmm0
        addps nb234_fiyM(%esp),%xmm1
        addps nb234_fizM(%esp),%xmm2
        movaps %xmm3,nb234_fjxH2(%esp)
        movaps %xmm4,nb234_fjyH2(%esp)
        movaps %xmm5,nb234_fjzH2(%esp)
        movaps %xmm0,nb234_fixM(%esp)
        movaps %xmm1,nb234_fiyM(%esp)
        movaps %xmm2,nb234_fizM(%esp)

        ## M-M interaction 
        movaps nb234_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234_qqMM(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb234_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps %xmm0,%xmm1
        movaps %xmm6,nb234_vctot(%esp)
        movaps %xmm0,%xmm2

        movaps nb234_fjxM(%esp),%xmm3
        movaps nb234_fjyM(%esp),%xmm4
        movaps nb234_fjzM(%esp),%xmm5
        mulps nb234_dxMM(%esp),%xmm0
        mulps nb234_dyMM(%esp),%xmm1
        mulps nb234_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb234_fixM(%esp),%xmm0
        addps nb234_fiyM(%esp),%xmm1
        addps nb234_fizM(%esp),%xmm2
        movaps %xmm3,nb234_fjxM(%esp)
        movaps %xmm4,nb234_fjyM(%esp)
        movaps %xmm5,nb234_fjzM(%esp)
        movaps %xmm0,nb234_fixM(%esp)
        movaps %xmm1,nb234_fiyM(%esp)
        movaps %xmm2,nb234_fizM(%esp)

        movl nb234_faction(%ebp),%edi
        ## update j forces 
        ## 4 j waters with four atoms each.
        ## step 1 : transpose fjxO, fjyO, fjzO, fjxH1
        movaps nb234_fjxO(%esp),%xmm0
        movaps nb234_fjyO(%esp),%xmm1
        movaps nb234_fjzO(%esp),%xmm2
        movaps nb234_fjxH1(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   (%edi,%eax,4),%xmm1
        movlps   (%edi,%ebx,4),%xmm5
        movlps   (%edi,%ecx,4),%xmm6
        movlps   (%edi,%edx,4),%xmm7
        movhps   8(%edi,%eax,4),%xmm1
        movhps   8(%edi,%ebx,4),%xmm5
        movhps   8(%edi,%ecx,4),%xmm6
        movhps   8(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,(%edi,%eax,4)
        movlps   %xmm5,(%edi,%ebx,4)
        movlps   %xmm6,(%edi,%ecx,4)
        movlps   %xmm7,(%edi,%edx,4)
        movhps   %xmm1,8(%edi,%eax,4)
        movhps   %xmm5,8(%edi,%ebx,4)
        movhps   %xmm6,8(%edi,%ecx,4)
        movhps   %xmm7,8(%edi,%edx,4)

        ## step 2 : transpose fjyH1, fjzH1, fjxH2, fjyH2
        movaps nb234_fjyH1(%esp),%xmm0
        movaps nb234_fjzH1(%esp),%xmm1
        movaps nb234_fjxH2(%esp),%xmm2
        movaps nb234_fjyH2(%esp),%xmm3
        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3
        ## results are now in xmm4, xmm2, xmm0, xmm3
        ## load the corresponding j forces from memory
        movlps   16(%edi,%eax,4),%xmm1
        movlps   16(%edi,%ebx,4),%xmm5
        movlps   16(%edi,%ecx,4),%xmm6
        movlps   16(%edi,%edx,4),%xmm7
        movhps   24(%edi,%eax,4),%xmm1
        movhps   24(%edi,%ebx,4),%xmm5
        movhps   24(%edi,%ecx,4),%xmm6
        movhps   24(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,16(%edi,%eax,4)
        movlps   %xmm5,16(%edi,%ebx,4)
        movlps   %xmm6,16(%edi,%ecx,4)
        movlps   %xmm7,16(%edi,%edx,4)
        movhps   %xmm1,24(%edi,%eax,4)
        movhps   %xmm5,24(%edi,%ebx,4)
        movhps   %xmm6,24(%edi,%ecx,4)
        movhps   %xmm7,24(%edi,%edx,4)

        ## step 3 : transpose fjzH2, fjxM, fjyM, fjzM. xmm4 is scratch
        movaps nb234_fjzH2(%esp),%xmm0
        movaps nb234_fjxM(%esp),%xmm1
        movaps nb234_fjyM(%esp),%xmm2
        movaps nb234_fjzM(%esp),%xmm3

        movaps %xmm0,%xmm4
        movaps %xmm1,%xmm5
        unpcklps %xmm2,%xmm4
        unpcklps %xmm3,%xmm5
        unpckhps %xmm2,%xmm0
        unpckhps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm0,%xmm3
        unpcklps %xmm5,%xmm4
        unpckhps %xmm5,%xmm2

        unpcklps %xmm1,%xmm0
        unpckhps %xmm1,%xmm3

        ## results are now in xmm0, xmm1, xmm2, xmm3
        ## load the corresponding j forces from memory
        movlps   32(%edi,%eax,4),%xmm1
        movlps   32(%edi,%ebx,4),%xmm5
        movlps   32(%edi,%ecx,4),%xmm6
        movlps   32(%edi,%edx,4),%xmm7
        movhps   40(%edi,%eax,4),%xmm1
        movhps   40(%edi,%ebx,4),%xmm5
        movhps   40(%edi,%ecx,4),%xmm6
        movhps   40(%edi,%edx,4),%xmm7
        ## add
        addps    %xmm4,%xmm1
        addps    %xmm2,%xmm5
        addps    %xmm0,%xmm6
        addps    %xmm3,%xmm7
        ## store back
        movlps   %xmm1,32(%edi,%eax,4)
        movlps   %xmm5,32(%edi,%ebx,4)
        movlps   %xmm6,32(%edi,%ecx,4)
        movlps   %xmm7,32(%edi,%edx,4)
        movhps   %xmm1,40(%edi,%eax,4)
        movhps   %xmm5,40(%edi,%ebx,4)
        movhps   %xmm6,40(%edi,%ecx,4)
        movhps   %xmm7,40(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb234_innerk(%esp)
        jl    _nb_kernel234_ia32_sse.nb234_single_check
        jmp   _nb_kernel234_ia32_sse.nb234_unroll_loop
_nb_kernel234_ia32_sse.nb234_single_check: 
        addl $4,nb234_innerk(%esp)
        jnz   _nb_kernel234_ia32_sse.nb234_single_loop
        jmp   _nb_kernel234_ia32_sse.nb234_updateouterdata
_nb_kernel234_ia32_sse.nb234_single_loop: 
        movl  nb234_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb234_innerjjnr(%esp)

        movl nb234_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb234_jxO(%esp)
        movaps %xmm3,nb234_jyO(%esp)
        movaps %xmm1,nb234_jzO(%esp)

        ## do O and M in parallel
        movaps nb234_ixO(%esp),%xmm0
        movaps nb234_iyO(%esp),%xmm1
        movaps nb234_izO(%esp),%xmm2
        movaps nb234_ixM(%esp),%xmm3
        movaps nb234_iyM(%esp),%xmm4
        movaps nb234_izM(%esp),%xmm5
        subps  nb234_jxO(%esp),%xmm0
        subps  nb234_jyO(%esp),%xmm1
        subps  nb234_jzO(%esp),%xmm2
        subps  nb234_jxO(%esp),%xmm3
        subps  nb234_jyO(%esp),%xmm4
        subps  nb234_jzO(%esp),%xmm5

        movaps %xmm0,nb234_dxOO(%esp)
        movaps %xmm1,nb234_dyOO(%esp)
        movaps %xmm2,nb234_dzOO(%esp)
        movaps %xmm3,nb234_dxMM(%esp)
        movaps %xmm4,nb234_dyMM(%esp)
        movaps %xmm5,nb234_dzMM(%esp)

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        ## Save data 
        movaps %xmm0,nb234_rsqOO(%esp)
        movaps %xmm4,nb234_rsqMM(%esp)

        ## do 1/x for O
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234_three(%esp),%xmm3
        mulps   %xmm0,%xmm1     ## rsq*lu*lu 
        subps   %xmm1,%xmm3     ## constant 30-rsq*lu*lu 
        mulps   %xmm2,%xmm3     ## lu*(3-rsq*lu*lu) 
        mulps   nb234_half(%esp),%xmm3
        movaps  %xmm3,nb234_rinvOO(%esp)        ## rinvH2 

        ## 1/sqrt(x) for M
        rsqrtps %xmm4,%xmm5
        movaps  %xmm5,%xmm6
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm7
        mulps   %xmm4,%xmm5
        subps   %xmm5,%xmm7
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm7   ## rinv iH1 - j water 
        movaps %xmm7,nb234_rinvMM(%esp)


        ## LJ table interaction
        movaps nb234_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb234_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb234_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movd %eax,%mm0

    movl nb234_VFtab(%ebp),%esi
    movd %mm6,%eax

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb234_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb234_c6(%esp),%xmm4
    mulss  %xmm4,%xmm7   ## fijD 
    mulss  %xmm4,%xmm5   ## Vvdw6 
        xorps  %xmm3,%xmm3

        mulps  nb234_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3
        movss  %xmm3,nb234_fstmp(%esp)

    addss  nb234_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb234_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  nb234_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb234_c12(%esp),%xmm4
    mulss  %xmm4,%xmm7 ## fijR 
    mulss  %xmm4,%xmm5 ## Vvdw12 
        movaps nb234_fstmp(%esp),%xmm3
        mulss  nb234_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3

    addss  nb234_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb234_Vvdwtot(%esp)

        mulss %xmm3,%xmm0
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movd  %mm0,%eax

        mulss  nb234_dxOO(%esp),%xmm0
        mulss  nb234_dyOO(%esp),%xmm1
        mulss  nb234_dzOO(%esp),%xmm2
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subss   %xmm0,%xmm3
        subss   %xmm1,%xmm4
        subss   %xmm2,%xmm5
        movaps  %xmm3,nb234_fjxO(%esp)
        movaps  %xmm4,nb234_fjyO(%esp)
        movaps  %xmm5,nb234_fjzO(%esp)
        addss   nb234_fixO(%esp),%xmm0
        addss   nb234_fiyO(%esp),%xmm1
        addss   nb234_fizO(%esp),%xmm2
        movss  %xmm0,nb234_fixO(%esp)
        movss  %xmm1,nb234_fiyO(%esp)
        movss  %xmm2,nb234_fizO(%esp)

        ## do  M coulomb interaction
        movaps nb234_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234_qqMH(%esp),%xmm3
        movhps  nb234_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  %xmm3,%xmm7 ## xmm7 = coul part of fscal 

        addps  nb234_vctot(%esp),%xmm6
        movaps %xmm6,nb234_vctot(%esp)
        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb234_dxMM(%esp),%xmm0
        mulps   nb234_dyMM(%esp),%xmm1
        mulps   nb234_dzMM(%esp),%xmm2
        ## update forces M - j water 
        movaps  nb234_fjxO(%esp),%xmm3
        movaps  nb234_fjyO(%esp),%xmm4
        movaps  nb234_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb234_fjxO(%esp)
        movaps  %xmm4,nb234_fjyO(%esp)
        movaps  %xmm5,nb234_fjzO(%esp)
        addps   nb234_fixM(%esp),%xmm0
        addps   nb234_fiyM(%esp),%xmm1
        addps   nb234_fizM(%esp),%xmm2
        movaps  %xmm0,nb234_fixM(%esp)
        movaps  %xmm1,nb234_fiyM(%esp)
        movaps  %xmm2,nb234_fizM(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb234_ixH1(%esp),%xmm0
        movaps  nb234_iyH1(%esp),%xmm1
        movaps  nb234_izH1(%esp),%xmm2
        movaps  nb234_ixH2(%esp),%xmm3
        movaps  nb234_iyH2(%esp),%xmm4
        movaps  nb234_izH2(%esp),%xmm5
        subps   nb234_jxO(%esp),%xmm0
        subps   nb234_jyO(%esp),%xmm1
        subps   nb234_jzO(%esp),%xmm2
        subps   nb234_jxO(%esp),%xmm3
        subps   nb234_jyO(%esp),%xmm4
        subps   nb234_jzO(%esp),%xmm5
        movaps %xmm0,nb234_dxH1H1(%esp)
        movaps %xmm1,nb234_dyH1H1(%esp)
        movaps %xmm2,nb234_dzH1H1(%esp)
        movaps %xmm3,nb234_dxH2H2(%esp)
        movaps %xmm4,nb234_dyH2H2(%esp)
        movaps %xmm5,nb234_dzH2H2(%esp)
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
        movaps  %xmm0,nb234_rsqH1H1(%esp)
        movaps  %xmm4,nb234_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb234_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb234_rinvH1H1(%esp)
        movaps  %xmm7,nb234_rinvH2H2(%esp)

        ## Do H1 coulomb interaction
        movaps nb234_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234_qqHH(%esp),%xmm3
        movhps  nb234_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  %xmm3,%xmm7 ## xmm7 = coul part of fscal 

        addps  nb234_vctot(%esp),%xmm6
        movaps %xmm6,nb234_vctot(%esp)

        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb234_dxH1H1(%esp),%xmm0
        mulps   nb234_dyH1H1(%esp),%xmm1
        mulps   nb234_dzH1H1(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb234_fjxO(%esp),%xmm3
        movaps  nb234_fjyO(%esp),%xmm4
        movaps  nb234_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb234_fjxO(%esp)
        movaps  %xmm4,nb234_fjyO(%esp)
        movaps  %xmm5,nb234_fjzO(%esp)
        addps   nb234_fixH1(%esp),%xmm0
        addps   nb234_fiyH1(%esp),%xmm1
        addps   nb234_fizH1(%esp),%xmm2
        movaps  %xmm0,nb234_fixH1(%esp)
        movaps  %xmm1,nb234_fiyH1(%esp)
        movaps  %xmm2,nb234_fizH1(%esp)

        ## H2 Coulomb
        movaps nb234_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234_qqHH(%esp),%xmm3
        movhps  nb234_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb234_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  %xmm3,%xmm7 ## xmm7 = coul part of fscal 

        addps  nb234_vctot(%esp),%xmm6   ## local vctot summation variable
        movaps %xmm6,nb234_vctot(%esp)
        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        mulps   nb234_dxH2H2(%esp),%xmm0
        mulps   nb234_dyH2H2(%esp),%xmm1
        mulps   nb234_dzH2H2(%esp),%xmm2
        ## update forces H2 - j water 
        movaps  nb234_fjxO(%esp),%xmm3
        movaps  nb234_fjyO(%esp),%xmm4
        movaps  nb234_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb234_fjxO(%esp)
        movaps  %xmm4,nb234_fjyO(%esp)
        movaps  %xmm5,nb234_fjzO(%esp)
        addps   nb234_fixH2(%esp),%xmm0
        addps   nb234_fiyH2(%esp),%xmm1
        addps   nb234_fizH2(%esp),%xmm2
        movaps  %xmm0,nb234_fixH2(%esp)
        movaps  %xmm1,nb234_fiyH2(%esp)
        movaps  %xmm2,nb234_fizH2(%esp)

        movl    nb234_faction(%ebp),%esi
        ## update j water forces from local variables.
        ## transpose back first
        movaps  nb234_fjxO(%esp),%xmm0   ## Ox H1x H2x Mx 
        movaps  nb234_fjyO(%esp),%xmm1   ## Oy H1y H2y My
        movaps  nb234_fjzO(%esp),%xmm2   ## Oz H1z H2z Mz

        movaps  %xmm0,%xmm3
        movaps  %xmm0,%xmm4
        unpcklps %xmm1,%xmm3            ## Ox Oy - -
        shufps $0x1,%xmm2,%xmm4        ## h1x - Oz -
        movaps  %xmm1,%xmm5
        movaps  %xmm0,%xmm6
        unpcklps %xmm2,%xmm5            ## - - H1y H1z
        unpckhps %xmm1,%xmm6            ## h2x h2y - - 
        unpckhps %xmm2,%xmm1            ## - - My Mz

        shufps  $50,%xmm0,%xmm2 ## (00110010) h2z - Mx -
        shufps  $36,%xmm4,%xmm3 ## constant 00100100 ;# Ox Oy Oz H1x 
        shufps  $78,%xmm6,%xmm5 ## constant 01001110 ;# h1y h1z h2x h2y
        shufps  $232,%xmm1,%xmm2 ## constant 11101000 ;# h2z mx my mz

        movlps  (%esi,%eax,4),%xmm0
        movlps  16(%esi,%eax,4),%xmm1
        movlps  32(%esi,%eax,4),%xmm4
        movhps  8(%esi,%eax,4),%xmm0
        movhps  24(%esi,%eax,4),%xmm1
        movhps  40(%esi,%eax,4),%xmm4
        addps   %xmm3,%xmm0
        addps   %xmm5,%xmm1
        addps   %xmm2,%xmm4
        movlps   %xmm0,(%esi,%eax,4)
        movlps   %xmm1,16(%esi,%eax,4)
        movlps   %xmm4,32(%esi,%eax,4)
        movhps   %xmm0,8(%esi,%eax,4)
        movhps   %xmm1,24(%esi,%eax,4)
        movhps   %xmm4,40(%esi,%eax,4)

        decl nb234_innerk(%esp)
        jz    _nb_kernel234_ia32_sse.nb234_updateouterdata
        jmp   _nb_kernel234_ia32_sse.nb234_single_loop
_nb_kernel234_ia32_sse.nb234_updateouterdata: 
        movl  nb234_ii3(%esp),%ecx
        movl  nb234_faction(%ebp),%edi
        movl  nb234_fshift(%ebp),%esi
        movl  nb234_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb234_fixO(%esp),%xmm0
        movaps nb234_fiyO(%esp),%xmm1
        movaps nb234_fizO(%esp),%xmm2

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
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb234_fixH1(%esp),%xmm0
        movaps nb234_fiyH1(%esp),%xmm1
        movaps nb234_fizH1(%esp),%xmm2

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
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb234_fixH2(%esp),%xmm0
        movaps nb234_fiyH2(%esp),%xmm1
        movaps nb234_fizH2(%esp),%xmm2

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
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb234_fixM(%esp),%xmm0
        movaps nb234_fiyM(%esp),%xmm1
        movaps nb234_fizM(%esp),%xmm2

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
        movss  36(%edi,%ecx,4),%xmm3
        movss  40(%edi,%ecx,4),%xmm4
        movss  44(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,36(%edi,%ecx,4)
        movss  %xmm4,40(%edi,%ecx,4)
        movss  %xmm5,44(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb234_n(%esp),%esi
        ## get group index for i particle 
        movl  nb234_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb234_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb234_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb234_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb234_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb234_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel234_ia32_sse.nb234_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb234_n(%esp)
        jmp _nb_kernel234_ia32_sse.nb234_outer
_nb_kernel234_ia32_sse.nb234_outerend: 
        ## check if more outer neighborlists remain
        movl  nb234_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel234_ia32_sse.nb234_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel234_ia32_sse.nb234_threadloop
_nb_kernel234_ia32_sse.nb234_end: 
        emms

        movl nb234_nouter(%esp),%eax
        movl nb234_ninner(%esp),%ebx
        movl nb234_outeriter(%ebp),%ecx
        movl nb234_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb234_salign(%esp),%eax
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



.globl nb_kernel234nf_ia32_sse
.globl _nb_kernel234nf_ia32_sse
nb_kernel234nf_ia32_sse:        
_nb_kernel234nf_ia32_sse:       
.set nb234nf_p_nri, 8
.set nb234nf_iinr, 12
.set nb234nf_jindex, 16
.set nb234nf_jjnr, 20
.set nb234nf_shift, 24
.set nb234nf_shiftvec, 28
.set nb234nf_fshift, 32
.set nb234nf_gid, 36
.set nb234nf_pos, 40
.set nb234nf_faction, 44
.set nb234nf_charge, 48
.set nb234nf_p_facel, 52
.set nb234nf_argkrf, 56
.set nb234nf_argcrf, 60
.set nb234nf_Vc, 64
.set nb234nf_type, 68
.set nb234nf_p_ntype, 72
.set nb234nf_vdwparam, 76
.set nb234nf_Vvdw, 80
.set nb234nf_p_tabscale, 84
.set nb234nf_VFtab, 88
.set nb234nf_invsqrta, 92
.set nb234nf_dvda, 96
.set nb234nf_p_gbtabscale, 100
.set nb234nf_GBtab, 104
.set nb234nf_p_nthreads, 108
.set nb234nf_count, 112
.set nb234nf_mtx, 116
.set nb234nf_outeriter, 120
.set nb234nf_inneriter, 124
.set nb234nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb234nf_ixO, 0
.set nb234nf_iyO, 16
.set nb234nf_izO, 32
.set nb234nf_ixH1, 48
.set nb234nf_iyH1, 64
.set nb234nf_izH1, 80
.set nb234nf_ixH2, 96
.set nb234nf_iyH2, 112
.set nb234nf_izH2, 128
.set nb234nf_ixM, 144
.set nb234nf_iyM, 160
.set nb234nf_izM, 176
.set nb234nf_jxO, 192
.set nb234nf_jyO, 208
.set nb234nf_jzO, 224
.set nb234nf_jxH1, 240
.set nb234nf_jyH1, 256
.set nb234nf_jzH1, 272
.set nb234nf_jxH2, 288
.set nb234nf_jyH2, 304
.set nb234nf_jzH2, 320
.set nb234nf_jxM, 336
.set nb234nf_jyM, 352
.set nb234nf_jzM, 368
.set nb234nf_qqMM, 384
.set nb234nf_qqMH, 400
.set nb234nf_qqHH, 416
.set nb234nf_c6, 432
.set nb234nf_c12, 448
.set nb234nf_tsc, 464
.set nb234nf_vctot, 480
.set nb234nf_Vvdwtot, 496
.set nb234nf_three, 512
.set nb234nf_rsqOO, 528
.set nb234nf_rsqH1H1, 544
.set nb234nf_rsqH1H2, 560
.set nb234nf_rsqH1M, 576
.set nb234nf_rsqH2H1, 592
.set nb234nf_rsqH2H2, 608
.set nb234nf_rsqH2M, 624
.set nb234nf_rsqMH1, 640
.set nb234nf_rsqMH2, 656
.set nb234nf_rsqMM, 672
.set nb234nf_rinvOO, 688
.set nb234nf_rinvH1H1, 704
.set nb234nf_rinvH1H2, 720
.set nb234nf_rinvH1M, 736
.set nb234nf_rinvH2H1, 752
.set nb234nf_rinvH2H2, 768
.set nb234nf_rinvH2M, 784
.set nb234nf_rinvMH1, 800
.set nb234nf_rinvMH2, 816
.set nb234nf_rinvMM, 832
.set nb234nf_krf, 848
.set nb234nf_crf, 864
.set nb234nf_half, 880
.set nb234nf_is3, 896
.set nb234nf_ii3, 900
.set nb234nf_innerjjnr, 904
.set nb234nf_innerk, 908
.set nb234nf_n, 912
.set nb234nf_nn1, 916
.set nb234nf_nri, 920
.set nb234nf_nouter, 924
.set nb234nf_ninner, 928
.set nb234nf_salign, 932
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
        movl %eax,nb234nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb234nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb234nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb234nf_nouter(%esp)
        movl %eax,nb234nf_ninner(%esp)

        movl nb234nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb234nf_tsc(%esp)

        movl nb234nf_argkrf(%ebp),%esi
        movl nb234nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb234nf_krf(%esp)
        movaps %xmm6,nb234nf_crf(%esp)
        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb234nf_half(%esp)
        movss nb234nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb234nf_half(%esp)
        movaps %xmm3,nb234nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb234nf_iinr(%ebp),%ecx     ## ecx = pointer into iinr[]
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb234nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm5
        movss 12(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movl nb234nf_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb234nf_qqMM(%esp)
        movaps %xmm4,nb234nf_qqMH(%esp)
        movaps %xmm5,nb234nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb234nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb234nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb234nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $0x55,%xmm1,%xmm1
        movaps %xmm0,nb234nf_c6(%esp)
        movaps %xmm1,nb234nf_c12(%esp)

_nb_kernel234nf_ia32_sse.nb234nf_threadloop: 
        movl  nb234nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel234nf_ia32_sse.nb234nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel234nf_ia32_sse.nb234nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb234nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb234nf_n(%esp)
        movl %ebx,nb234nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel234nf_ia32_sse.nb234nf_outerstart
        jmp _nb_kernel234nf_ia32_sse.nb234nf_end

_nb_kernel234nf_ia32_sse.nb234nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb234nf_nouter(%esp),%ebx
        movl %ebx,nb234nf_nouter(%esp)

_nb_kernel234nf_ia32_sse.nb234nf_outer: 
        movl  nb234nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb234nf_is3(%esp)            ## store is3 

        movl  nb234nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb234nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb234nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb234nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        addss (%eax,%ebx,4),%xmm3       ## ox
        addss 4(%eax,%ebx,4),%xmm4     ## oy
        addss 8(%eax,%ebx,4),%xmm5     ## oz
        addss 12(%eax,%ebx,4),%xmm6    ## h1x
        addss 16(%eax,%ebx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb234nf_ixO(%esp)
        movaps %xmm4,nb234nf_iyO(%esp)
        movaps %xmm5,nb234nf_izO(%esp)
        movaps %xmm6,nb234nf_ixH1(%esp)
        movaps %xmm7,nb234nf_iyH1(%esp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%eax,%ebx,4),%xmm6    ## h1z
        addss 24(%eax,%ebx,4),%xmm0    ## h2x
        addss 28(%eax,%ebx,4),%xmm1    ## h2y
        addss 32(%eax,%ebx,4),%xmm2    ## h2z
        addss 36(%eax,%ebx,4),%xmm3    ## mx
        addss 40(%eax,%ebx,4),%xmm4    ## my
        addss 44(%eax,%ebx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb234nf_izH1(%esp)
        movaps %xmm0,nb234nf_ixH2(%esp)
        movaps %xmm1,nb234nf_iyH2(%esp)
        movaps %xmm2,nb234nf_izH2(%esp)
        movaps %xmm3,nb234nf_ixM(%esp)
        movaps %xmm4,nb234nf_iyM(%esp)
        movaps %xmm5,nb234nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb234nf_vctot(%esp)
        movaps %xmm4,nb234nf_Vvdwtot(%esp)

        movl  nb234nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb234nf_pos(%ebp),%esi
        movl  nb234nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb234nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb234nf_ninner(%esp),%ecx
        movl  %ecx,nb234nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb234nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel234nf_ia32_sse.nb234nf_unroll_loop
        jmp   _nb_kernel234nf_ia32_sse.nb234nf_single_check
_nb_kernel234nf_ia32_sse.nb234nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb234nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb234nf_innerjjnr(%esp)             ## advance pointer (unroll 4) 

        movl nb234nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables
        ## Load Ox, Oy, Oz, H1x 
        movlps (%esi,%eax,4),%xmm1      ##  Oxa   Oya    -    -
        movlps (%esi,%ecx,4),%xmm4      ##  Oxc   Oyc    -    -
        movhps (%esi,%ebx,4),%xmm1      ##  Oxa   Oya   Oxb   Oyb 
        movhps (%esi,%edx,4),%xmm4      ##  Oxc   Oyc   Oxd   Oyd 
        movaps %xmm1,%xmm0              ##  Oxa   Oya   Oxb   Oyb 
        shufps $0x88,%xmm4,%xmm0       ##  Oxa   Oxb   Oxc   Oxd
        shufps $0xDD,%xmm4,%xmm1       ##  Oya   Oyb   Oyc   Oyd
        movlps 8(%esi,%eax,4),%xmm3     ##  Oza  H1xa    -    -
        movlps 8(%esi,%ecx,4),%xmm5     ##  Ozc  H1xc    -    -
        movhps 8(%esi,%ebx,4),%xmm3     ##  Oza  H1xa   Ozb  H1xb 
        movhps 8(%esi,%edx,4),%xmm5     ##  Ozc  H1xc   Ozd  H1xd 
        movaps %xmm3,%xmm2              ##  Oza  H1xa   Ozb  H1xb 
        shufps $0x88,%xmm5,%xmm2       ##  Oza   Ozb   Ozc   Ozd
        shufps $0xDD,%xmm5,%xmm3       ## H1xa  H1xb  H1xc  H1xd
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234nf_jxO(%esp)
        movaps %xmm1,nb234nf_jyO(%esp)
        movaps %xmm2,nb234nf_jzO(%esp)
        movaps %xmm3,nb234nf_jxH1(%esp)

        ## Load H1y H1z H2x H2y 
        movlps 16(%esi,%eax,4),%xmm1
        movlps 16(%esi,%ecx,4),%xmm4
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 16(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 24(%esi,%eax,4),%xmm3
        movlps 24(%esi,%ecx,4),%xmm5
        movhps 24(%esi,%ebx,4),%xmm3
        movhps 24(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234nf_jyH1(%esp)
        movaps %xmm1,nb234nf_jzH1(%esp)
        movaps %xmm2,nb234nf_jxH2(%esp)
        movaps %xmm3,nb234nf_jyH2(%esp)

        ## Load H2z Mx My Mz 
        movlps 32(%esi,%eax,4),%xmm1
        movlps 32(%esi,%ecx,4),%xmm4
        movhps 32(%esi,%ebx,4),%xmm1
        movhps 32(%esi,%edx,4),%xmm4
        movaps %xmm1,%xmm0
        shufps $0x88,%xmm4,%xmm0
        shufps $0xDD,%xmm4,%xmm1
        movlps 40(%esi,%eax,4),%xmm3
        movlps 40(%esi,%ecx,4),%xmm5
        movhps 40(%esi,%ebx,4),%xmm3
        movhps 40(%esi,%edx,4),%xmm5
        movaps %xmm3,%xmm2
        shufps $0x88,%xmm5,%xmm2
        shufps $0xDD,%xmm5,%xmm3
        ## coordinates in xmm0-xmm3     
        ## store
        movaps %xmm0,nb234nf_jzH2(%esp)
        movaps %xmm1,nb234nf_jxM(%esp)
        movaps %xmm2,nb234nf_jyM(%esp)
        movaps %xmm3,nb234nf_jzM(%esp)

        ## start calculating pairwise distances
        movaps nb234nf_ixO(%esp),%xmm0
        movaps nb234nf_iyO(%esp),%xmm1
        movaps nb234nf_izO(%esp),%xmm2
        movaps nb234nf_ixH1(%esp),%xmm3
        movaps nb234nf_iyH1(%esp),%xmm4
        movaps nb234nf_izH1(%esp),%xmm5
        subps  nb234nf_jxO(%esp),%xmm0
        subps  nb234nf_jyO(%esp),%xmm1
        subps  nb234nf_jzO(%esp),%xmm2
        subps  nb234nf_jxH1(%esp),%xmm3
        subps  nb234nf_jyH1(%esp),%xmm4
        subps  nb234nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb234nf_rsqOO(%esp)
        movaps %xmm3,nb234nf_rsqH1H1(%esp)

        movaps nb234nf_ixH1(%esp),%xmm0
        movaps nb234nf_iyH1(%esp),%xmm1
        movaps nb234nf_izH1(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234nf_jxH2(%esp),%xmm0
        subps  nb234nf_jyH2(%esp),%xmm1
        subps  nb234nf_jzH2(%esp),%xmm2
        subps  nb234nf_jxM(%esp),%xmm3
        subps  nb234nf_jyM(%esp),%xmm4
        subps  nb234nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb234nf_rsqH1H2(%esp)
        movaps %xmm3,nb234nf_rsqH1M(%esp)

        movaps nb234nf_ixH2(%esp),%xmm0
        movaps nb234nf_iyH2(%esp),%xmm1
        movaps nb234nf_izH2(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234nf_jxH1(%esp),%xmm0
        subps  nb234nf_jyH1(%esp),%xmm1
        subps  nb234nf_jzH1(%esp),%xmm2
        subps  nb234nf_jxH2(%esp),%xmm3
        subps  nb234nf_jyH2(%esp),%xmm4
        subps  nb234nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb234nf_rsqH2H1(%esp)
        movaps %xmm3,nb234nf_rsqH2H2(%esp)

        movaps nb234nf_ixH2(%esp),%xmm0
        movaps nb234nf_iyH2(%esp),%xmm1
        movaps nb234nf_izH2(%esp),%xmm2
        movaps nb234nf_ixM(%esp),%xmm3
        movaps nb234nf_iyM(%esp),%xmm4
        movaps nb234nf_izM(%esp),%xmm5
        subps  nb234nf_jxM(%esp),%xmm0
        subps  nb234nf_jyM(%esp),%xmm1
        subps  nb234nf_jzM(%esp),%xmm2
        subps  nb234nf_jxH1(%esp),%xmm3
        subps  nb234nf_jyH1(%esp),%xmm4
        subps  nb234nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb234nf_rsqH2M(%esp)
        movaps %xmm4,nb234nf_rsqMH1(%esp)

        movaps nb234nf_ixM(%esp),%xmm0
        movaps nb234nf_iyM(%esp),%xmm1
        movaps nb234nf_izM(%esp),%xmm2
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        subps  nb234nf_jxH2(%esp),%xmm0
        subps  nb234nf_jyH2(%esp),%xmm1
        subps  nb234nf_jzH2(%esp),%xmm2
        subps  nb234nf_jxM(%esp),%xmm3
        subps  nb234nf_jyM(%esp),%xmm4
        subps  nb234nf_jzM(%esp),%xmm5
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
        movaps %xmm0,nb234nf_rsqMH2(%esp)
        movaps %xmm4,nb234nf_rsqMM(%esp)

        ## start by doing invsqrt for OO
        rsqrtps nb234nf_rsqOO(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234nf_three(%esp),%xmm3
        mulps   nb234nf_rsqOO(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb234nf_half(%esp),%xmm3
        movaps  %xmm3,nb234nf_rinvOO(%esp)

        ## more invsqrt ops - do two at a time.
        rsqrtps nb234nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb234nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234nf_rsqH1H1(%esp),%xmm1
        mulps   nb234nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm3   ## rinvH1H1 
        mulps   nb234nf_half(%esp),%xmm7   ## rinvH1H2 
        movaps  %xmm3,nb234nf_rinvH1H1(%esp)
        movaps  %xmm7,nb234nf_rinvH1H2(%esp)

        rsqrtps nb234nf_rsqH1M(%esp),%xmm1
        rsqrtps nb234nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234nf_rsqH1M(%esp),%xmm1
        mulps   nb234nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm3
        mulps   nb234nf_half(%esp),%xmm7
        movaps  %xmm3,nb234nf_rinvH1M(%esp)
        movaps  %xmm7,nb234nf_rinvH2H1(%esp)

        rsqrtps nb234nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb234nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234nf_rsqH2H2(%esp),%xmm1
        mulps   nb234nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm3
        mulps   nb234nf_half(%esp),%xmm7
        movaps  %xmm3,nb234nf_rinvH2H2(%esp)
        movaps  %xmm7,nb234nf_rinvH2M(%esp)

        rsqrtps nb234nf_rsqMH1(%esp),%xmm1
        rsqrtps nb234nf_rsqMH2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb234nf_rsqMH1(%esp),%xmm1
        mulps   nb234nf_rsqMH2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm3
        mulps   nb234nf_half(%esp),%xmm7
        movaps  %xmm3,nb234nf_rinvMH1(%esp)
        movaps  %xmm7,nb234nf_rinvMH2(%esp)

        rsqrtps nb234nf_rsqMM(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234nf_three(%esp),%xmm3
        mulps   nb234nf_rsqMM(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb234nf_half(%esp),%xmm3
        movaps  %xmm3,nb234nf_rinvMM(%esp)

        ## start with OO LJ interaction
        movaps nb234nf_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb234nf_tsc(%esp),%xmm1

        movhlps %xmm1,%xmm2
    cvttps2pi %xmm1,%mm6
    cvttps2pi %xmm2,%mm7    ## mm6/mm7 contain lu indices 
    cvtpi2ps %mm6,%xmm3
    cvtpi2ps %mm7,%xmm2
        movlhps  %xmm2,%xmm3
        subps    %xmm3,%xmm1    ## xmm1=eps 
    movaps %xmm1,%xmm2
    mulps  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6
    pslld $3,%mm7

    movd %mm6,%eax
    psrlq $32,%mm6
    movd %mm7,%ecx
    psrlq $32,%mm7
    movd %mm6,%ebx
    movd %mm7,%edx

    movl nb234nf_VFtab(%ebp),%esi

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movlps (%esi,%ecx,4),%xmm7
    movhps (%esi,%ebx,4),%xmm5
    movhps (%esi,%edx,4),%xmm7 ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movlps 8(%esi,%ecx,4),%xmm3
    movhps 8(%esi,%ebx,4),%xmm7
    movhps 8(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb234nf_c6(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 

    addps  nb234nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb234nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movlps 16(%esi,%ecx,4),%xmm7
    movhps 16(%esi,%ebx,4),%xmm5
    movhps 16(%esi,%edx,4),%xmm7    ## got half table 

    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movlps 24(%esi,%ecx,4),%xmm3
    movhps 24(%esi,%ebx,4),%xmm7
    movhps 24(%esi,%edx,4),%xmm3    ## other half of table  
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## repulsion table ready, in xmm4-xmm7 
    mulps  %xmm1,%xmm6      ## xmm6=Geps 
    mulps  %xmm2,%xmm7      ## xmm7=Heps2 
    addps  %xmm6,%xmm5
    addps  %xmm7,%xmm5      ## xmm5=Fp 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb234nf_c12(%esp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 

    addps  nb234nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb234nf_Vvdwtot(%esp)

        ## Coulomb interactions 
        ## start with H1-H1 interaction 
        movaps nb234nf_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulps  nb234nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234nf_crf(%esp),%xmm6
        mulps  nb234nf_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addps  nb234nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## H1-H2 interaction 
        movaps nb234nf_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps  %xmm0,%xmm0
        mulps  nb234nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## H1-M interaction  
        movaps nb234nf_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=Rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H1 interaction 
        movaps nb234nf_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H2 interaction 
        movaps nb234nf_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## H2-M interaction 
        movaps nb234nf_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## M-H1 interaction 
        movaps nb234nf_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## M-H2 interaction 
        movaps nb234nf_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 

        ## M-M interaction 
        movaps nb234nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb234nf_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb234nf_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb234nf_qqMM(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm6,nb234nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb234nf_innerk(%esp)
        jl    _nb_kernel234nf_ia32_sse.nb234nf_single_check
        jmp   _nb_kernel234nf_ia32_sse.nb234nf_unroll_loop
_nb_kernel234nf_ia32_sse.nb234nf_single_check: 
        addl $4,nb234nf_innerk(%esp)
        jnz   _nb_kernel234nf_ia32_sse.nb234nf_single_loop
        jmp   _nb_kernel234nf_ia32_sse.nb234nf_updateouterdata
_nb_kernel234nf_ia32_sse.nb234nf_single_loop: 
        movl  nb234nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb234nf_innerjjnr(%esp)

        movl nb234nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates
        movlps (%esi,%eax,4),%xmm3              ##  Ox  Oy  
        movlps 16(%esi,%eax,4),%xmm4            ## H1y H1z 
        movlps 32(%esi,%eax,4),%xmm5            ## H2z  Mx 
        movhps 8(%esi,%eax,4),%xmm3             ##  Ox  Oy  Oz H1x
        movhps 24(%esi,%eax,4),%xmm4            ## H1y H1z H2x H2y
        movhps 40(%esi,%eax,4),%xmm5            ## H2z  Mx  My  Mz
        ## transpose
        movaps %xmm4,%xmm0
        movaps %xmm3,%xmm1
        movaps %xmm4,%xmm2
        movaps %xmm3,%xmm6
        shufps $18,%xmm5,%xmm4 ## (00010010)  h2x - Mx  - 
        shufps $193,%xmm0,%xmm3 ## (11000001)  Oy  - H1y - 
        shufps $35,%xmm5,%xmm2 ## (00100011) H2y - My  - 
        shufps $18,%xmm0,%xmm1 ## (00010010)  Oz  - H1z - 
        ##  xmm6: Ox - - H1x   xmm5: H2z - - Mz 
        shufps $140,%xmm4,%xmm6 ## (10001100) Ox H1x H2x Mx 
        shufps $136,%xmm2,%xmm3 ## (10001000) Oy H1y H2y My 
        shufps $200,%xmm5,%xmm1 ## (11001000) Oz H1z H2z Mz

        ## store all j coordinates in jO  
        movaps %xmm6,nb234nf_jxO(%esp)
        movaps %xmm3,nb234nf_jyO(%esp)
        movaps %xmm1,nb234nf_jzO(%esp)

        ## do O and M in parallel
        movaps nb234nf_ixO(%esp),%xmm0
        movaps nb234nf_iyO(%esp),%xmm1
        movaps nb234nf_izO(%esp),%xmm2
        movaps nb234nf_ixM(%esp),%xmm3
        movaps nb234nf_iyM(%esp),%xmm4
        movaps nb234nf_izM(%esp),%xmm5
        subps  nb234nf_jxO(%esp),%xmm0
        subps  nb234nf_jyO(%esp),%xmm1
        subps  nb234nf_jzO(%esp),%xmm2
        subps  nb234nf_jxO(%esp),%xmm3
        subps  nb234nf_jyO(%esp),%xmm4
        subps  nb234nf_jzO(%esp),%xmm5

        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm4
        addps %xmm5,%xmm4       ## have rsq in xmm4
        ## Save data 
        movaps %xmm0,nb234nf_rsqOO(%esp)
        movaps %xmm4,nb234nf_rsqMM(%esp)

        ## do 1/x for O
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb234nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1     ## rsq*lu*lu 
        subps   %xmm1,%xmm3     ## constant 30-rsq*lu*lu 
        mulps   %xmm2,%xmm3     ## lu*(3-rsq*lu*lu) 
        mulps   nb234nf_half(%esp),%xmm3
        movaps  %xmm3,nb234nf_rinvOO(%esp)      ## rinvH2 

        ## 1/sqrt(x) for M
        rsqrtps %xmm4,%xmm5
        movaps  %xmm5,%xmm6
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm7
        mulps   %xmm4,%xmm5
        subps   %xmm5,%xmm7
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm7   ## rinv iH1 - j water 
        movaps %xmm7,nb234nf_rinvMM(%esp)


        ## LJ table interaction
        movaps nb234nf_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb234nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb234nf_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movl nb234nf_VFtab(%ebp),%esi
    movd %mm6,%eax

    ## dispersion 
    movlps (%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 8(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## dispersion table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb234nf_c6(%esp),%xmm4
    mulss  %xmm4,%xmm5   ## Vvdw6 

    addss  nb234nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb234nf_Vvdwtot(%esp)

    ## repulsion 
    movlps 16(%esi,%eax,4),%xmm5
    movaps %xmm5,%xmm4
    shufps $136,%xmm7,%xmm4 ## constant 10001000
    shufps $221,%xmm7,%xmm5 ## constant 11011101

    movlps 24(%esi,%eax,4),%xmm7
    movaps %xmm7,%xmm6
    shufps $136,%xmm3,%xmm6 ## constant 10001000
    shufps $221,%xmm3,%xmm7 ## constant 11011101
    ## table ready, in xmm4-xmm7 
    mulss  %xmm1,%xmm6      ## xmm6=Geps 
    mulss  %xmm2,%xmm7      ## xmm7=Heps2 
    addss  %xmm6,%xmm5
    addss  %xmm7,%xmm5      ## xmm5=Fp 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb234nf_c12(%esp),%xmm4
    mulss  %xmm4,%xmm5 ## Vvdw12 
    addss  nb234nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb234nf_Vvdwtot(%esp)

        ## do  M coulomb interaction
        movaps nb234nf_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234nf_qqMH(%esp),%xmm3
        movhps  nb234nf_qqMM(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234nf_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234nf_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addps  nb234nf_vctot(%esp),%xmm6
        movaps %xmm6,nb234nf_vctot(%esp)

        ## i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb234nf_ixH1(%esp),%xmm0
        movaps  nb234nf_iyH1(%esp),%xmm1
        movaps  nb234nf_izH1(%esp),%xmm2
        movaps  nb234nf_ixH2(%esp),%xmm3
        movaps  nb234nf_iyH2(%esp),%xmm4
        movaps  nb234nf_izH2(%esp),%xmm5
        subps   nb234nf_jxO(%esp),%xmm0
        subps   nb234nf_jyO(%esp),%xmm1
        subps   nb234nf_jzO(%esp),%xmm2
        subps   nb234nf_jxO(%esp),%xmm3
        subps   nb234nf_jyO(%esp),%xmm4
        subps   nb234nf_jzO(%esp),%xmm5
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
        movaps  %xmm0,nb234nf_rsqH1H1(%esp)
        movaps  %xmm4,nb234nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb234nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb234nf_half(%esp),%xmm3   ## rinvH1H1
        mulps   nb234nf_half(%esp),%xmm7   ## rinvH2H2
        movaps  %xmm3,nb234nf_rinvH1H1(%esp)
        movaps  %xmm7,nb234nf_rinvH2H2(%esp)

        ## Do H1 coulomb interaction
        movaps nb234nf_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234nf_qqHH(%esp),%xmm3
        movhps  nb234nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234nf_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addps  nb234nf_vctot(%esp),%xmm6
        movaps %xmm6,nb234nf_vctot(%esp)

        ## H2 Coulomb
        movaps nb234nf_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb234nf_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        ## fetch charges to xmm3 (temporary) 
        xorps  %xmm3,%xmm3
        movss   nb234nf_qqHH(%esp),%xmm3
        movhps  nb234nf_qqMH(%esp),%xmm3
        shufps $193,%xmm3,%xmm3 ## constant 11000001 

        mulps  nb234nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb234nf_crf(%esp),%xmm6
        mulps  %xmm3,%xmm6 ## xmm6=voul=qq*(rinv+ krsq-crf) 

        addps  nb234nf_vctot(%esp),%xmm6   ## local vctot summation variable
        movaps %xmm6,nb234nf_vctot(%esp)

        decl nb234nf_innerk(%esp)
        jz    _nb_kernel234nf_ia32_sse.nb234nf_updateouterdata
        jmp   _nb_kernel234nf_ia32_sse.nb234nf_single_loop
_nb_kernel234nf_ia32_sse.nb234nf_updateouterdata: 
        ## get n from stack
        movl nb234nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb234nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb234nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb234nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb234nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb234nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb234nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel234nf_ia32_sse.nb234nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb234nf_n(%esp)
        jmp _nb_kernel234nf_ia32_sse.nb234nf_outer
_nb_kernel234nf_ia32_sse.nb234nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb234nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel234nf_ia32_sse.nb234nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel234nf_ia32_sse.nb234nf_threadloop
_nb_kernel234nf_ia32_sse.nb234nf_end: 
        emms

        movl nb234nf_nouter(%esp),%eax
        movl nb234nf_ninner(%esp),%ebx
        movl nb234nf_outeriter(%ebp),%ecx
        movl nb234nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb234nf_salign(%esp),%eax
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


