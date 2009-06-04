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



.globl nb_kernel132_ia32_sse
.globl _nb_kernel132_ia32_sse
nb_kernel132_ia32_sse:  
_nb_kernel132_ia32_sse: 
.set nb132_p_nri, 8
.set nb132_iinr, 12
.set nb132_jindex, 16
.set nb132_jjnr, 20
.set nb132_shift, 24
.set nb132_shiftvec, 28
.set nb132_fshift, 32
.set nb132_gid, 36
.set nb132_pos, 40
.set nb132_faction, 44
.set nb132_charge, 48
.set nb132_p_facel, 52
.set nb132_argkrf, 56
.set nb132_argcrf, 60
.set nb132_Vc, 64
.set nb132_type, 68
.set nb132_p_ntype, 72
.set nb132_vdwparam, 76
.set nb132_Vvdw, 80
.set nb132_p_tabscale, 84
.set nb132_VFtab, 88
.set nb132_invsqrta, 92
.set nb132_dvda, 96
.set nb132_p_gbtabscale, 100
.set nb132_GBtab, 104
.set nb132_p_nthreads, 108
.set nb132_count, 112
.set nb132_mtx, 116
.set nb132_outeriter, 120
.set nb132_inneriter, 124
.set nb132_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb132_ixO, 0
.set nb132_iyO, 16
.set nb132_izO, 32
.set nb132_ixH1, 48
.set nb132_iyH1, 64
.set nb132_izH1, 80
.set nb132_ixH2, 96
.set nb132_iyH2, 112
.set nb132_izH2, 128
.set nb132_jxO, 144
.set nb132_jyO, 160
.set nb132_jzO, 176
.set nb132_jxH1, 192
.set nb132_jyH1, 208
.set nb132_jzH1, 224
.set nb132_jxH2, 240
.set nb132_jyH2, 256
.set nb132_jzH2, 272
.set nb132_dxOO, 288
.set nb132_dyOO, 304
.set nb132_dzOO, 320
.set nb132_dxOH1, 336
.set nb132_dyOH1, 352
.set nb132_dzOH1, 368
.set nb132_dxOH2, 384
.set nb132_dyOH2, 400
.set nb132_dzOH2, 416
.set nb132_dxH1O, 432
.set nb132_dyH1O, 448
.set nb132_dzH1O, 464
.set nb132_dxH1H1, 480
.set nb132_dyH1H1, 496
.set nb132_dzH1H1, 512
.set nb132_dxH1H2, 528
.set nb132_dyH1H2, 544
.set nb132_dzH1H2, 560
.set nb132_dxH2O, 576
.set nb132_dyH2O, 592
.set nb132_dzH2O, 608
.set nb132_dxH2H1, 624
.set nb132_dyH2H1, 640
.set nb132_dzH2H1, 656
.set nb132_dxH2H2, 672
.set nb132_dyH2H2, 688
.set nb132_dzH2H2, 704
.set nb132_qqOO, 720
.set nb132_qqOH, 736
.set nb132_qqHH, 752
.set nb132_c6, 768
.set nb132_c12, 784
.set nb132_tsc, 800
.set nb132_fstmp, 816
.set nb132_vctot, 832
.set nb132_Vvdwtot, 848
.set nb132_fixO, 864
.set nb132_fiyO, 880
.set nb132_fizO, 896
.set nb132_fixH1, 912
.set nb132_fiyH1, 928
.set nb132_fizH1, 944
.set nb132_fixH2, 960
.set nb132_fiyH2, 976
.set nb132_fizH2, 992
.set nb132_fjxO, 1008
.set nb132_fjyO, 1024
.set nb132_fjzO, 1040
.set nb132_fjxH1, 1056
.set nb132_fjyH1, 1072
.set nb132_fjzH1, 1088
.set nb132_fjxH2, 1104
.set nb132_fjyH2, 1120
.set nb132_fjzH2, 1136
.set nb132_fjzH2b, 1140
.set nb132_fjzH2c, 1144
.set nb132_fjzH2d, 1148
.set nb132_half, 1152
.set nb132_three, 1168
.set nb132_rsqOO, 1184
.set nb132_rsqOH1, 1200
.set nb132_rsqOH2, 1216
.set nb132_rsqH1O, 1232
.set nb132_rsqH1H1, 1248
.set nb132_rsqH1H2, 1264
.set nb132_rsqH2O, 1280
.set nb132_rsqH2H1, 1296
.set nb132_rsqH2H2, 1312
.set nb132_rinvOO, 1328
.set nb132_rinvOH1, 1344
.set nb132_rinvOH2, 1360
.set nb132_rinvH1O, 1376
.set nb132_rinvH1H1, 1392
.set nb132_rinvH1H2, 1408
.set nb132_rinvH2O, 1424
.set nb132_rinvH2H1, 1440
.set nb132_rinvH2H2, 1456
.set nb132_two, 1472
.set nb132_is3, 1520
.set nb132_ii3, 1524
.set nb132_innerjjnr, 1528
.set nb132_innerk, 1532
.set nb132_n, 1536
.set nb132_nn1, 1540
.set nb132_nri, 1544
.set nb132_nouter, 1548
.set nb132_ninner, 1552
.set nb132_salign, 1556
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1560,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb132_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb132_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb132_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb132_nouter(%esp)
        movl %eax,nb132_ninner(%esp)

        movl nb132_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb132_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb132_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb132_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb132_p_facel(%ebp),%esi
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
        movaps %xmm3,nb132_qqOO(%esp)
        movaps %xmm4,nb132_qqOH(%esp)
        movaps %xmm5,nb132_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb132_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb132_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb132_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## constant 01010101
        movaps %xmm0,nb132_c6(%esp)
        movaps %xmm1,nb132_c12(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb132_half(%esp)
        movss nb132_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb132_half(%esp)
        movaps %xmm2,nb132_two(%esp)
        movaps %xmm3,nb132_three(%esp)

_nb_kernel132_ia32_sse.nb132_threadloop: 
        movl  nb132_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel132_ia32_sse.nb132_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel132_ia32_sse.nb132_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb132_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb132_n(%esp)
        movl %ebx,nb132_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel132_ia32_sse.nb132_outerstart
        jmp _nb_kernel132_ia32_sse.nb132_end

_nb_kernel132_ia32_sse.nb132_outerstart: 
        ## ebx contains number of outer iterations
        addl nb132_nouter(%esp),%ebx
        movl %ebx,nb132_nouter(%esp)

_nb_kernel132_ia32_sse.nb132_outer: 
        movl  nb132_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb132_is3(%esp)      ## store is3 

        movl  nb132_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb132_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb132_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb132_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132_ixO(%esp)
        movaps %xmm4,nb132_iyO(%esp)
        movaps %xmm5,nb132_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb132_ixH1(%esp)
        movaps %xmm1,nb132_iyH1(%esp)
        movaps %xmm2,nb132_izH1(%esp)
        movaps %xmm3,nb132_ixH2(%esp)
        movaps %xmm4,nb132_iyH2(%esp)
        movaps %xmm5,nb132_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb132_vctot(%esp)
        movaps %xmm4,nb132_Vvdwtot(%esp)
        movaps %xmm4,nb132_fixO(%esp)
        movaps %xmm4,nb132_fiyO(%esp)
        movaps %xmm4,nb132_fizO(%esp)
        movaps %xmm4,nb132_fixH1(%esp)
        movaps %xmm4,nb132_fiyH1(%esp)
        movaps %xmm4,nb132_fizH1(%esp)
        movaps %xmm4,nb132_fixH2(%esp)
        movaps %xmm4,nb132_fiyH2(%esp)
        movaps %xmm4,nb132_fizH2(%esp)

        movl  nb132_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb132_pos(%ebp),%esi
        movl  nb132_faction(%ebp),%edi
        movl  nb132_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb132_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb132_ninner(%esp),%ecx
        movl  %ecx,nb132_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb132_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel132_ia32_sse.nb132_unroll_loop
        jmp   _nb_kernel132_ia32_sse.nb132_single_check
_nb_kernel132_ia32_sse.nb132_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb132_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb132_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb132_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d 
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd 
        movaps %xmm0,nb132_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb132_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb132_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb132_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb132_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb132_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb132_jzO(%esp)
        movaps %xmm1,nb132_jzH1(%esp)
        movaps %xmm2,nb132_jzH2(%esp)

        movaps nb132_ixO(%esp),%xmm0
        movaps nb132_iyO(%esp),%xmm1
        movaps nb132_izO(%esp),%xmm2
        movaps nb132_ixO(%esp),%xmm3
        movaps nb132_iyO(%esp),%xmm4
        movaps nb132_izO(%esp),%xmm5
        subps  nb132_jxO(%esp),%xmm0
        subps  nb132_jyO(%esp),%xmm1
        subps  nb132_jzO(%esp),%xmm2
        subps  nb132_jxH1(%esp),%xmm3
        subps  nb132_jyH1(%esp),%xmm4
        subps  nb132_jzH1(%esp),%xmm5
        movaps %xmm0,nb132_dxOO(%esp)
        movaps %xmm1,nb132_dyOO(%esp)
        movaps %xmm2,nb132_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxOH1(%esp)
        movaps %xmm4,nb132_dyOH1(%esp)
        movaps %xmm5,nb132_dzOH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb132_rsqOO(%esp)
        movaps %xmm3,nb132_rsqOH1(%esp)

        movaps nb132_ixO(%esp),%xmm0
        movaps nb132_iyO(%esp),%xmm1
        movaps nb132_izO(%esp),%xmm2
        movaps nb132_ixH1(%esp),%xmm3
        movaps nb132_iyH1(%esp),%xmm4
        movaps nb132_izH1(%esp),%xmm5
        subps  nb132_jxH2(%esp),%xmm0
        subps  nb132_jyH2(%esp),%xmm1
        subps  nb132_jzH2(%esp),%xmm2
        subps  nb132_jxO(%esp),%xmm3
        subps  nb132_jyO(%esp),%xmm4
        subps  nb132_jzO(%esp),%xmm5
        movaps %xmm0,nb132_dxOH2(%esp)
        movaps %xmm1,nb132_dyOH2(%esp)
        movaps %xmm2,nb132_dzOH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH1O(%esp)
        movaps %xmm4,nb132_dyH1O(%esp)
        movaps %xmm5,nb132_dzH1O(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb132_rsqOH2(%esp)
        movaps %xmm3,nb132_rsqH1O(%esp)

        movaps nb132_ixH1(%esp),%xmm0
        movaps nb132_iyH1(%esp),%xmm1
        movaps nb132_izH1(%esp),%xmm2
        movaps nb132_ixH1(%esp),%xmm3
        movaps nb132_iyH1(%esp),%xmm4
        movaps nb132_izH1(%esp),%xmm5
        subps  nb132_jxH1(%esp),%xmm0
        subps  nb132_jyH1(%esp),%xmm1
        subps  nb132_jzH1(%esp),%xmm2
        subps  nb132_jxH2(%esp),%xmm3
        subps  nb132_jyH2(%esp),%xmm4
        subps  nb132_jzH2(%esp),%xmm5
        movaps %xmm0,nb132_dxH1H1(%esp)
        movaps %xmm1,nb132_dyH1H1(%esp)
        movaps %xmm2,nb132_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH1H2(%esp)
        movaps %xmm4,nb132_dyH1H2(%esp)
        movaps %xmm5,nb132_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb132_rsqH1H1(%esp)
        movaps %xmm3,nb132_rsqH1H2(%esp)

        movaps nb132_ixH2(%esp),%xmm0
        movaps nb132_iyH2(%esp),%xmm1
        movaps nb132_izH2(%esp),%xmm2
        movaps nb132_ixH2(%esp),%xmm3
        movaps nb132_iyH2(%esp),%xmm4
        movaps nb132_izH2(%esp),%xmm5
        subps  nb132_jxO(%esp),%xmm0
        subps  nb132_jyO(%esp),%xmm1
        subps  nb132_jzO(%esp),%xmm2
        subps  nb132_jxH1(%esp),%xmm3
        subps  nb132_jyH1(%esp),%xmm4
        subps  nb132_jzH1(%esp),%xmm5
        movaps %xmm0,nb132_dxH2O(%esp)
        movaps %xmm1,nb132_dyH2O(%esp)
        movaps %xmm2,nb132_dzH2O(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb132_dxH2H1(%esp)
        movaps %xmm4,nb132_dyH2H1(%esp)
        movaps %xmm5,nb132_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb132_rsqH2O(%esp)
        movaps %xmm4,nb132_rsqH2H1(%esp)

        movaps nb132_ixH2(%esp),%xmm0
        movaps nb132_iyH2(%esp),%xmm1
        movaps nb132_izH2(%esp),%xmm2
        subps  nb132_jxH2(%esp),%xmm0
        subps  nb132_jyH2(%esp),%xmm1
        subps  nb132_jzH2(%esp),%xmm2
        movaps %xmm0,nb132_dxH2H2(%esp)
        movaps %xmm1,nb132_dyH2H2(%esp)
        movaps %xmm2,nb132_dzH2H2(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb132_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb132_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb132_rinvH2H2(%esp)
        movaps  %xmm7,nb132_rinvH2H1(%esp)

        rsqrtps nb132_rsqOO(%esp),%xmm1
        rsqrtps nb132_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132_rsqOO(%esp),%xmm1
        mulps   nb132_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%esp),%xmm3
        mulps   nb132_half(%esp),%xmm7
        movaps  %xmm3,nb132_rinvOO(%esp)
        movaps  %xmm7,nb132_rinvOH1(%esp)

        rsqrtps nb132_rsqOH2(%esp),%xmm1
        rsqrtps nb132_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132_rsqOH2(%esp),%xmm1
        mulps   nb132_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%esp),%xmm3
        mulps   nb132_half(%esp),%xmm7
        movaps  %xmm3,nb132_rinvOH2(%esp)
        movaps  %xmm7,nb132_rinvH1O(%esp)

        rsqrtps nb132_rsqH1H1(%esp),%xmm1
        rsqrtps nb132_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132_rsqH1H1(%esp),%xmm1
        mulps   nb132_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%esp),%xmm3
        mulps   nb132_half(%esp),%xmm7
        movaps  %xmm3,nb132_rinvH1H1(%esp)
        movaps  %xmm7,nb132_rinvH1H2(%esp)

        rsqrtps nb132_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132_three(%esp),%xmm3
        mulps   nb132_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132_half(%esp),%xmm3
        movaps  %xmm3,nb132_rinvH2O(%esp)

        ## start with OO interaction - first the table LJ part
        movaps nb132_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb132_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb132_tsc(%esp),%xmm1

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

    movl nb132_VFtab(%ebp),%esi

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
    mulps  nb132_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb132_c6(%esp),%xmm4
    mulps  %xmm4,%xmm7   ## fijD 
    mulps  %xmm4,%xmm5   ## Vvdw6 
    movaps  %xmm7,nb132_fstmp(%esp)

    addps  nb132_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb132_Vvdwtot(%esp)

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
    mulps  nb132_two(%esp),%xmm7         ## two*Heps2 
    addps  %xmm6,%xmm7
    addps  %xmm5,%xmm7 ## xmm7=FF 
    mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addps  %xmm4,%xmm5 ## xmm5=VV 

    movaps nb132_c12(%esp),%xmm4
    mulps  %xmm4,%xmm7 ## fijR 
    mulps  %xmm4,%xmm5 ## Vvdw12 
    addps  nb132_fstmp(%esp),%xmm7
    mulps nb132_tsc(%esp),%xmm7

    addps  nb132_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb132_Vvdwtot(%esp)

        movaps nb132_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm6
        mulps  nb132_qqOO(%esp),%xmm6
        movaps %xmm6,%xmm2

        addps  nb132_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulps  %xmm0,%xmm2
        subps  %xmm7,%xmm2
        mulps  %xmm0,%xmm2

        movaps %xmm2,%xmm0
        movaps %xmm2,%xmm1

    movd  %mm0,%eax
    movd  %mm1,%ebx
    movd  %mm2,%ecx
    movd  %mm3,%edx

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb132_dxOO(%esp),%xmm0
        mulps nb132_dyOO(%esp),%xmm1
        mulps nb132_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixO(%esp),%xmm0
        addps nb132_fiyO(%esp),%xmm1
        addps nb132_fizO(%esp),%xmm2
        movaps %xmm3,nb132_fjxO(%esp)
        movaps %xmm4,nb132_fjyO(%esp)
        movaps %xmm5,nb132_fjzO(%esp)
        movaps %xmm0,nb132_fixO(%esp)
        movaps %xmm1,nb132_fiyO(%esp)
        movaps %xmm2,nb132_fizO(%esp)

        ## O-H1 interaction 
        movaps nb132_rinvOH1(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb132_dxOH1(%esp),%xmm0
        mulps nb132_dyOH1(%esp),%xmm1
        mulps nb132_dzOH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixO(%esp),%xmm0
        addps nb132_fiyO(%esp),%xmm1
        addps nb132_fizO(%esp),%xmm2
        movaps %xmm3,nb132_fjxH1(%esp)
        movaps %xmm4,nb132_fjyH1(%esp)
        movaps %xmm5,nb132_fjzH1(%esp)
        movaps %xmm0,nb132_fixO(%esp)
        movaps %xmm1,nb132_fiyO(%esp)
        movaps %xmm2,nb132_fizO(%esp)

        ## O-H2 interaction  
        movaps nb132_rinvOH2(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb132_dxOH2(%esp),%xmm0
        mulps nb132_dyOH2(%esp),%xmm1
        mulps nb132_dzOH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixO(%esp),%xmm0
        addps nb132_fiyO(%esp),%xmm1
        addps nb132_fizO(%esp),%xmm2
        movaps %xmm3,nb132_fjxH2(%esp)
        movaps %xmm4,nb132_fjyH2(%esp)
        movaps %xmm5,nb132_fjzH2(%esp)
        movaps %xmm0,nb132_fixO(%esp)
        movaps %xmm1,nb132_fiyO(%esp)
        movaps %xmm2,nb132_fizO(%esp)

        ## H1-O interaction 
        movaps nb132_rinvH1O(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb132_fjxO(%esp),%xmm3
        movaps nb132_fjyO(%esp),%xmm4
        movaps nb132_fjzO(%esp),%xmm5
        mulps nb132_dxH1O(%esp),%xmm0
        mulps nb132_dyH1O(%esp),%xmm1
        mulps nb132_dzH1O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH1(%esp),%xmm0
        addps nb132_fiyH1(%esp),%xmm1
        addps nb132_fizH1(%esp),%xmm2
        movaps %xmm3,nb132_fjxO(%esp)
        movaps %xmm4,nb132_fjyO(%esp)
        movaps %xmm5,nb132_fjzO(%esp)
        movaps %xmm0,nb132_fixH1(%esp)
        movaps %xmm1,nb132_fiyH1(%esp)
        movaps %xmm2,nb132_fizH1(%esp)

        ## H1-H1 interaction 
        movaps nb132_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb132_fjxH1(%esp),%xmm3
        movaps nb132_fjyH1(%esp),%xmm4
        movaps nb132_fjzH1(%esp),%xmm5
        mulps nb132_dxH1H1(%esp),%xmm0
        mulps nb132_dyH1H1(%esp),%xmm1
        mulps nb132_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH1(%esp),%xmm0
        addps nb132_fiyH1(%esp),%xmm1
        addps nb132_fizH1(%esp),%xmm2
        movaps %xmm3,nb132_fjxH1(%esp)
        movaps %xmm4,nb132_fjyH1(%esp)
        movaps %xmm5,nb132_fjzH1(%esp)
        movaps %xmm0,nb132_fixH1(%esp)
        movaps %xmm1,nb132_fiyH1(%esp)
        movaps %xmm2,nb132_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb132_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb132_fjxH2(%esp),%xmm3
        movaps nb132_fjyH2(%esp),%xmm4
        movaps nb132_fjzH2(%esp),%xmm5
        mulps nb132_dxH1H2(%esp),%xmm0
        mulps nb132_dyH1H2(%esp),%xmm1
        mulps nb132_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH1(%esp),%xmm0
        addps nb132_fiyH1(%esp),%xmm1
        addps nb132_fizH1(%esp),%xmm2
        movaps %xmm3,nb132_fjxH2(%esp)
        movaps %xmm4,nb132_fjyH2(%esp)
        movaps %xmm5,nb132_fjzH2(%esp)
        movaps %xmm0,nb132_fixH1(%esp)
        movaps %xmm1,nb132_fiyH1(%esp)
        movaps %xmm2,nb132_fizH1(%esp)

        ## H2-O interaction 
        movaps nb132_rinvH2O(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb132_fjxO(%esp),%xmm3
        movaps nb132_fjyO(%esp),%xmm4
        movaps nb132_fjzO(%esp),%xmm5
        mulps nb132_dxH2O(%esp),%xmm0
        mulps nb132_dyH2O(%esp),%xmm1
        mulps nb132_dzH2O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH2(%esp),%xmm0
        addps nb132_fiyH2(%esp),%xmm1
        addps nb132_fizH2(%esp),%xmm2
        movaps %xmm3,nb132_fjxO(%esp)
        movaps %xmm4,nb132_fjyO(%esp)
        movaps %xmm5,nb132_fjzO(%esp)
        movaps %xmm0,nb132_fixH2(%esp)
        movaps %xmm1,nb132_fiyH2(%esp)
        movaps %xmm2,nb132_fizH2(%esp)

        ## H2-H1 interaction 
        movaps nb132_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb132_fjxH1(%esp),%xmm3
        movaps nb132_fjyH1(%esp),%xmm4
        movaps nb132_fjzH1(%esp),%xmm5
        mulps nb132_dxH2H1(%esp),%xmm0
        mulps nb132_dyH2H1(%esp),%xmm1
        mulps nb132_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH2(%esp),%xmm0
        addps nb132_fiyH2(%esp),%xmm1
        addps nb132_fizH2(%esp),%xmm2
        movaps %xmm3,nb132_fjxH1(%esp)
        movaps %xmm4,nb132_fjyH1(%esp)
        movaps %xmm5,nb132_fjzH1(%esp)
        movaps %xmm0,nb132_fixH2(%esp)
        movaps %xmm1,nb132_fiyH2(%esp)
        movaps %xmm2,nb132_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb132_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm4      ## xmm7=rinv 
        mulps  %xmm0,%xmm0  ## rinvsq
        mulps  nb132_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulps  %xmm4,%xmm0

        addps  %xmm4,%xmm6      ## add to local vctot 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps %xmm0,%xmm1
        movaps %xmm6,nb132_vctot(%esp)
        movaps %xmm0,%xmm2

        movaps nb132_fjxH2(%esp),%xmm3
        movaps nb132_fjyH2(%esp),%xmm4
        movaps nb132_fjzH2(%esp),%xmm5
        mulps nb132_dxH2H2(%esp),%xmm0
        mulps nb132_dyH2H2(%esp),%xmm1
        mulps nb132_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb132_fixH2(%esp),%xmm0
        addps nb132_fiyH2(%esp),%xmm1
        addps nb132_fizH2(%esp),%xmm2
        movaps %xmm3,nb132_fjxH2(%esp)
        movaps %xmm4,nb132_fjyH2(%esp)
        movaps %xmm5,nb132_fjzH2(%esp)
        movaps %xmm0,nb132_fixH2(%esp)
        movaps %xmm1,nb132_fiyH2(%esp)
        movaps %xmm2,nb132_fizH2(%esp)

        movl nb132_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxO, fjyO, fjzO, ... , fjzH2.
        ## Each position is a quadruplet of forces for the four 
        ## corresponding j waters, so we need to transpose them before
        ## adding to the memory positions.
        ## 
        ## This _used_ to be a simple transpose, but the resulting high number
        ## of unaligned 128-bit load/stores might trigger a possible hardware 
        ## bug on Athlon and Opteron chips, so I have worked around it
        ## to use 64-bit load/stores instead. The performance hit should be
        ## very modest, since the 128-bit unaligned memory instructions were
        ## slow anyway. 


        ## 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
        movaps nb132_fjxO(%esp),%xmm0   ## xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
        movaps nb132_fjyO(%esp),%xmm2   ## xmm1= fjyOa  fjyOb  fjyOc  fjyOd
        movlps (%edi,%eax,4),%xmm3
        movlps (%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxOa  fjyOa  fjxOb  fjyOb
        unpckhps %xmm2,%xmm1       ## xmm1= fjxOc  fjyOc  fjxOd  fjyOd
        movhps (%edi,%ebx,4),%xmm3
        movhps (%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,(%edi,%eax,4)
        movlps %xmm4,(%edi,%ecx,4)
        movhps %xmm3,(%edi,%ebx,4)
        movhps %xmm4,(%edi,%edx,4)

        ## Oxygen Z & first hydrogen X forces for 4 j particles 
        movaps nb132_fjzO(%esp),%xmm0    ## xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
        movaps nb132_fjxH1(%esp),%xmm2   ## xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
        movlps 8(%edi,%eax,4),%xmm3
        movlps 8(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
        movhps 8(%edi,%ebx,4),%xmm3
        movhps 8(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,8(%edi,%eax,4)
        movlps %xmm4,8(%edi,%ecx,4)
        movhps %xmm3,8(%edi,%ebx,4)
        movhps %xmm4,8(%edi,%edx,4)


        ## First hydrogen Y & Z forces for 4 j particles 
        movaps nb132_fjyH1(%esp),%xmm0    ## xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
        movaps nb132_fjzH1(%esp),%xmm2   ## xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
        movlps 16(%edi,%eax,4),%xmm3
        movlps 16(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
        movhps 16(%edi,%ebx,4),%xmm3
        movhps 16(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,16(%edi,%eax,4)
        movlps %xmm4,16(%edi,%ecx,4)
        movhps %xmm3,16(%edi,%ebx,4)
        movhps %xmm4,16(%edi,%edx,4)


        ## Second hydrogen X & Y forces for 4 j particles 
        movaps nb132_fjxH2(%esp),%xmm0    ## xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
        movaps nb132_fjyH2(%esp),%xmm2   ## xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
        movlps 24(%edi,%eax,4),%xmm3
        movlps 24(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
        movhps 24(%edi,%ebx,4),%xmm3
        movhps 24(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,24(%edi,%eax,4)
        movlps %xmm4,24(%edi,%ecx,4)
        movhps %xmm3,24(%edi,%ebx,4)
        movhps %xmm4,24(%edi,%edx,4)


        ## Second hydrogen Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 32(%edi,%eax,4),%xmm4
        movss 32(%edi,%ebx,4),%xmm5
        movss 32(%edi,%ecx,4),%xmm6
        movss 32(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb132_fjzH2(%esp),%xmm4
        addss nb132_fjzH2b(%esp),%xmm5
        addss nb132_fjzH2c(%esp),%xmm6
        addss nb132_fjzH2d(%esp),%xmm7
        ## store back
        movss %xmm4,32(%edi,%eax,4)
        movss %xmm5,32(%edi,%ebx,4)
        movss %xmm6,32(%edi,%ecx,4)
        movss %xmm7,32(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb132_innerk(%esp)
        jl    _nb_kernel132_ia32_sse.nb132_single_check
        jmp   _nb_kernel132_ia32_sse.nb132_unroll_loop
_nb_kernel132_ia32_sse.nb132_single_check: 
        addl $4,nb132_innerk(%esp)
        jnz   _nb_kernel132_ia32_sse.nb132_single_loop
        jmp   _nb_kernel132_ia32_sse.nb132_updateouterdata
_nb_kernel132_ia32_sse.nb132_single_loop: 
        movl  nb132_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb132_innerjjnr(%esp)

        movl nb132_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb132_ixO(%esp),%xmm0
        movaps  nb132_iyO(%esp),%xmm1
        movaps  nb132_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb132_jxO(%esp)
        movaps %xmm4,nb132_jyO(%esp)
        movaps %xmm5,nb132_jzO(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb132_dxOO(%esp)
        movaps %xmm1,nb132_dyOO(%esp)
        movaps %xmm2,nb132_dzOO(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        movaps %xmm0,nb132_rsqOO(%esp)

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132_half(%esp),%xmm3   ## rinv iO - j water in xmm3
        movaps  %xmm3,nb132_rinvOO(%esp)

        movaps  %xmm3,%xmm0 ## rinv

        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb132_qqOO(%esp),%xmm4
        movhps  nb132_qqOH(%esp),%xmm4

        mulps  %xmm4,%xmm3  ## vcoul
        movaps %xmm3,%xmm6

        mulps %xmm0,%xmm3
        movaps %xmm3,nb132_fstmp(%esp)   ## save it
        addps  nb132_vctot(%esp),%xmm6
    movaps %xmm6,nb132_vctot(%esp)

        movaps nb132_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb132_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb132_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movd %eax,%mm0

    movl nb132_VFtab(%ebp),%esi
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
    mulss  nb132_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb132_c6(%esp),%xmm4
    mulss  %xmm4,%xmm7   ## fijD 
    mulss  %xmm4,%xmm5   ## Vvdw6 
        movss  nb132_fstmp(%esp),%xmm3
        mulps  nb132_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3
        movss  %xmm3,nb132_fstmp(%esp)

    addss  nb132_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb132_Vvdwtot(%esp)

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
    mulss  nb132_two(%esp),%xmm7         ## two*Heps2 
    addss  %xmm6,%xmm7
    addss  %xmm5,%xmm7 ## xmm7=FF 
    mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
    addss  %xmm4,%xmm5 ## xmm5=VV 

    movss nb132_c12(%esp),%xmm4
    mulss  %xmm4,%xmm7 ## fijR 
    mulss  %xmm4,%xmm5 ## Vvdw12 
        movaps nb132_fstmp(%esp),%xmm3
        mulss  nb132_tsc(%esp),%xmm7
        subss  %xmm7,%xmm3

    addss  nb132_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb132_Vvdwtot(%esp)

        mulps  %xmm3,%xmm0
        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2

        mulps   nb132_dxOO(%esp),%xmm0
        mulps   nb132_dyOO(%esp),%xmm1
        mulps   nb132_dzOO(%esp),%xmm2

        movd %mm0,%eax

        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb132_fjxO(%esp)
        movaps  %xmm4,nb132_fjyO(%esp)
        movaps  %xmm5,nb132_fjzO(%esp)
        addps   nb132_fixO(%esp),%xmm0
        addps   nb132_fiyO(%esp),%xmm1
        addps   nb132_fizO(%esp),%xmm2
        movaps  %xmm0,nb132_fixO(%esp)
        movaps  %xmm1,nb132_fiyO(%esp)
        movaps  %xmm2,nb132_fizO(%esp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb132_ixH1(%esp),%xmm0
        movaps  nb132_iyH1(%esp),%xmm1
        movaps  nb132_izH1(%esp),%xmm2
        movaps  nb132_ixH2(%esp),%xmm3
        movaps  nb132_iyH2(%esp),%xmm4
        movaps  nb132_izH2(%esp),%xmm5
        subps   nb132_jxO(%esp),%xmm0
        subps   nb132_jyO(%esp),%xmm1
        subps   nb132_jzO(%esp),%xmm2
        subps   nb132_jxO(%esp),%xmm3
        subps   nb132_jyO(%esp),%xmm4
        subps   nb132_jzO(%esp),%xmm5
        movaps %xmm0,nb132_dxH1O(%esp)
        movaps %xmm1,nb132_dyH1O(%esp)
        movaps %xmm2,nb132_dzH1O(%esp)
        movaps %xmm3,nb132_dxH2O(%esp)
        movaps %xmm4,nb132_dyH2O(%esp)
        movaps %xmm5,nb132_dzH2O(%esp)
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
        movaps  nb132_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb132_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb132_qqOH(%esp),%xmm6
        movhps  nb132_qqHH(%esp),%xmm6

        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm1
        mulps   %xmm6,%xmm3 ## vcoul
        mulps   %xmm6,%xmm7 ## vcoul
        mulps   %xmm0,%xmm0
        mulps   %xmm1,%xmm1
        mulps   %xmm3,%xmm0
        addps   %xmm7,%xmm3
        mulps   %xmm1,%xmm7
        addps   nb132_vctot(%esp),%xmm3
        movaps  %xmm3,nb132_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb132_dxH1O(%esp),%xmm0
        mulps   nb132_dyH1O(%esp),%xmm1
        mulps   nb132_dzH1O(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb132_fjxO(%esp),%xmm3
        movaps  nb132_fjyO(%esp),%xmm4
        movaps  nb132_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb132_fjxO(%esp)
        movaps  %xmm4,nb132_fjyO(%esp)
        movaps  %xmm5,nb132_fjzO(%esp)
        addps   nb132_fixH1(%esp),%xmm0
        addps   nb132_fiyH1(%esp),%xmm1
        addps   nb132_fizH1(%esp),%xmm2
        movaps  %xmm0,nb132_fixH1(%esp)
        movaps  %xmm1,nb132_fiyH1(%esp)
        movaps  %xmm2,nb132_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb132_dxH2O(%esp),%xmm0
        mulps   nb132_dyH2O(%esp),%xmm1
        mulps   nb132_dzH2O(%esp),%xmm2
        movaps  nb132_fjxO(%esp),%xmm3
        movaps  nb132_fjyO(%esp),%xmm4
        movaps  nb132_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb132_faction(%ebp),%esi
        movaps  %xmm3,nb132_fjxO(%esp)
        movaps  %xmm4,nb132_fjyO(%esp)
        movaps  %xmm5,nb132_fjzO(%esp)
        addps   nb132_fixH2(%esp),%xmm0
        addps   nb132_fiyH2(%esp),%xmm1
        addps   nb132_fizH2(%esp),%xmm2
        movaps  %xmm0,nb132_fixH2(%esp)
        movaps  %xmm1,nb132_fiyH2(%esp)
        movaps  %xmm2,nb132_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  (%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb132_fjxO(%esp),%xmm3
        movaps  nb132_fjyO(%esp),%xmm4
        movaps  nb132_fjzO(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   8(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,8(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl nb132_innerk(%esp)
        jz    _nb_kernel132_ia32_sse.nb132_updateouterdata
        jmp   _nb_kernel132_ia32_sse.nb132_single_loop
_nb_kernel132_ia32_sse.nb132_updateouterdata: 
        movl  nb132_ii3(%esp),%ecx
        movl  nb132_faction(%ebp),%edi
        movl  nb132_fshift(%ebp),%esi
        movl  nb132_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb132_fixO(%esp),%xmm0
        movaps nb132_fiyO(%esp),%xmm1
        movaps nb132_fizO(%esp),%xmm2

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
        movaps nb132_fixH1(%esp),%xmm0
        movaps nb132_fiyH1(%esp),%xmm1
        movaps nb132_fizH1(%esp),%xmm2

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
        movaps nb132_fixH2(%esp),%xmm0
        movaps nb132_fiyH2(%esp),%xmm1
        movaps nb132_fizH2(%esp),%xmm2

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

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb132_n(%esp),%esi
        ## get group index for i particle 
        movl  nb132_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb132_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb132_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb132_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb132_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb132_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel132_ia32_sse.nb132_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb132_n(%esp)
        jmp _nb_kernel132_ia32_sse.nb132_outer
_nb_kernel132_ia32_sse.nb132_outerend: 
        ## check if more outer neighborlists remain
        movl  nb132_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel132_ia32_sse.nb132_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel132_ia32_sse.nb132_threadloop
_nb_kernel132_ia32_sse.nb132_end: 
        emms

        movl nb132_nouter(%esp),%eax
        movl nb132_ninner(%esp),%ebx
        movl nb132_outeriter(%ebp),%ecx
        movl nb132_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb132_salign(%esp),%eax
        addl %eax,%esp
        addl $1560,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


.globl nb_kernel132nf_ia32_sse
.globl _nb_kernel132nf_ia32_sse
nb_kernel132nf_ia32_sse:        
_nb_kernel132nf_ia32_sse:       
.set nb132nf_p_nri, 8
.set nb132nf_iinr, 12
.set nb132nf_jindex, 16
.set nb132nf_jjnr, 20
.set nb132nf_shift, 24
.set nb132nf_shiftvec, 28
.set nb132nf_fshift, 32
.set nb132nf_gid, 36
.set nb132nf_pos, 40
.set nb132nf_faction, 44
.set nb132nf_charge, 48
.set nb132nf_p_facel, 52
.set nb132nf_argkrf, 56
.set nb132nf_argcrf, 60
.set nb132nf_Vc, 64
.set nb132nf_type, 68
.set nb132nf_p_ntype, 72
.set nb132nf_vdwparam, 76
.set nb132nf_Vvdw, 80
.set nb132nf_p_tabscale, 84
.set nb132nf_VFtab, 88
.set nb132nf_invsqrta, 92
.set nb132nf_dvda, 96
.set nb132nf_p_gbtabscale, 100
.set nb132nf_GBtab, 104
.set nb132nf_p_nthreads, 108
.set nb132nf_count, 112
.set nb132nf_mtx, 116
.set nb132nf_outeriter, 120
.set nb132nf_inneriter, 124
.set nb132nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb132nf_ixO, 0
.set nb132nf_iyO, 16
.set nb132nf_izO, 32
.set nb132nf_ixH1, 48
.set nb132nf_iyH1, 64
.set nb132nf_izH1, 80
.set nb132nf_ixH2, 96
.set nb132nf_iyH2, 112
.set nb132nf_izH2, 128
.set nb132nf_jxO, 144
.set nb132nf_jyO, 160
.set nb132nf_jzO, 176
.set nb132nf_jxH1, 192
.set nb132nf_jyH1, 208
.set nb132nf_jzH1, 224
.set nb132nf_jxH2, 240
.set nb132nf_jyH2, 256
.set nb132nf_jzH2, 272
.set nb132nf_qqOO, 288
.set nb132nf_qqOH, 304
.set nb132nf_qqHH, 320
.set nb132nf_c6, 336
.set nb132nf_c12, 352
.set nb132nf_tsc, 368
.set nb132nf_vctot, 384
.set nb132nf_Vvdwtot, 400
.set nb132nf_half, 416
.set nb132nf_three, 432
.set nb132nf_rsqOO, 448
.set nb132nf_rsqOH1, 464
.set nb132nf_rsqOH2, 480
.set nb132nf_rsqH1O, 496
.set nb132nf_rsqH1H1, 512
.set nb132nf_rsqH1H2, 528
.set nb132nf_rsqH2O, 544
.set nb132nf_rsqH2H1, 560
.set nb132nf_rsqH2H2, 576
.set nb132nf_rinvOO, 592
.set nb132nf_rinvOH1, 608
.set nb132nf_rinvOH2, 624
.set nb132nf_rinvH1O, 640
.set nb132nf_rinvH1H1, 656
.set nb132nf_rinvH1H2, 672
.set nb132nf_rinvH2O, 688
.set nb132nf_rinvH2H1, 704
.set nb132nf_rinvH2H2, 720
.set nb132nf_is3, 768
.set nb132nf_ii3, 772
.set nb132nf_innerjjnr, 776
.set nb132nf_innerk, 780
.set nb132nf_n, 784
.set nb132nf_nn1, 788
.set nb132nf_nri, 792
.set nb132nf_nouter, 796
.set nb132nf_ninner, 800
.set nb132nf_salign, 804
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $808,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb132nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb132nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb132nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb132nf_nouter(%esp)
        movl %eax,nb132nf_ninner(%esp)

        movl nb132nf_p_tabscale(%ebp),%eax
        movss (%eax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb132nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb132nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb132nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb132nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb132nf_qqOO(%esp)
        movaps %xmm4,nb132nf_qqOH(%esp)
        movaps %xmm5,nb132nf_qqHH(%esp)

        xorps %xmm0,%xmm0
        movl  nb132nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb132nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb132nf_vdwparam(%ebp),%eax
        movlps (%eax,%edx,4),%xmm0
        movaps %xmm0,%xmm1
        shufps $0,%xmm0,%xmm0
        shufps $85,%xmm1,%xmm1 ## constant 01010101
        movaps %xmm0,nb132nf_c6(%esp)
        movaps %xmm1,nb132nf_c12(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb132nf_half(%esp)
        movss nb132nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb132nf_half(%esp)
        movaps %xmm3,nb132nf_three(%esp)

_nb_kernel132nf_ia32_sse.nb132nf_threadloop: 
        movl  nb132nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel132nf_ia32_sse.nb132nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel132nf_ia32_sse.nb132nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb132nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb132nf_n(%esp)
        movl %ebx,nb132nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel132nf_ia32_sse.nb132nf_outerstart
        jmp _nb_kernel132nf_ia32_sse.nb132nf_end

_nb_kernel132nf_ia32_sse.nb132nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb132nf_nouter(%esp),%ebx
        movl %ebx,nb132nf_nouter(%esp)

_nb_kernel132nf_ia32_sse.nb132nf_outer: 
        movl  nb132nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb132nf_is3(%esp)            ## store is3 

        movl  nb132nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb132nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb132nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb132nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb132nf_ixO(%esp)
        movaps %xmm4,nb132nf_iyO(%esp)
        movaps %xmm5,nb132nf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb132nf_ixH1(%esp)
        movaps %xmm1,nb132nf_iyH1(%esp)
        movaps %xmm2,nb132nf_izH1(%esp)
        movaps %xmm3,nb132nf_ixH2(%esp)
        movaps %xmm4,nb132nf_iyH2(%esp)
        movaps %xmm5,nb132nf_izH2(%esp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb132nf_vctot(%esp)
        movaps %xmm4,nb132nf_Vvdwtot(%esp)

        movl  nb132nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb132nf_pos(%ebp),%esi
        movl  nb132nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb132nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb132nf_ninner(%esp),%ecx
        movl  %ecx,nb132nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb132nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel132nf_ia32_sse.nb132nf_unroll_loop
        jmp   _nb_kernel132nf_ia32_sse.nb132nf_single_check
_nb_kernel132nf_ia32_sse.nb132nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb132nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb132nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb132nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d 
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd 
        movaps %xmm0,nb132nf_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb132nf_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb132nf_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb132nf_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb132nf_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb132nf_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb132nf_jzO(%esp)
        movaps %xmm1,nb132nf_jzH1(%esp)
        movaps %xmm2,nb132nf_jzH2(%esp)

        movaps nb132nf_ixO(%esp),%xmm0
        movaps nb132nf_iyO(%esp),%xmm1
        movaps nb132nf_izO(%esp),%xmm2
        movaps nb132nf_ixO(%esp),%xmm3
        movaps nb132nf_iyO(%esp),%xmm4
        movaps nb132nf_izO(%esp),%xmm5
        subps  nb132nf_jxO(%esp),%xmm0
        subps  nb132nf_jyO(%esp),%xmm1
        subps  nb132nf_jzO(%esp),%xmm2
        subps  nb132nf_jxH1(%esp),%xmm3
        subps  nb132nf_jyH1(%esp),%xmm4
        subps  nb132nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb132nf_rsqOO(%esp)
        movaps %xmm3,nb132nf_rsqOH1(%esp)

        movaps nb132nf_ixO(%esp),%xmm0
        movaps nb132nf_iyO(%esp),%xmm1
        movaps nb132nf_izO(%esp),%xmm2
        movaps nb132nf_ixH1(%esp),%xmm3
        movaps nb132nf_iyH1(%esp),%xmm4
        movaps nb132nf_izH1(%esp),%xmm5
        subps  nb132nf_jxH2(%esp),%xmm0
        subps  nb132nf_jyH2(%esp),%xmm1
        subps  nb132nf_jzH2(%esp),%xmm2
        subps  nb132nf_jxO(%esp),%xmm3
        subps  nb132nf_jyO(%esp),%xmm4
        subps  nb132nf_jzO(%esp),%xmm5
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
        movaps %xmm0,nb132nf_rsqOH2(%esp)
        movaps %xmm3,nb132nf_rsqH1O(%esp)

        movaps nb132nf_ixH1(%esp),%xmm0
        movaps nb132nf_iyH1(%esp),%xmm1
        movaps nb132nf_izH1(%esp),%xmm2
        movaps nb132nf_ixH1(%esp),%xmm3
        movaps nb132nf_iyH1(%esp),%xmm4
        movaps nb132nf_izH1(%esp),%xmm5
        subps  nb132nf_jxH1(%esp),%xmm0
        subps  nb132nf_jyH1(%esp),%xmm1
        subps  nb132nf_jzH1(%esp),%xmm2
        subps  nb132nf_jxH2(%esp),%xmm3
        subps  nb132nf_jyH2(%esp),%xmm4
        subps  nb132nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb132nf_rsqH1H1(%esp)
        movaps %xmm3,nb132nf_rsqH1H2(%esp)

        movaps nb132nf_ixH2(%esp),%xmm0
        movaps nb132nf_iyH2(%esp),%xmm1
        movaps nb132nf_izH2(%esp),%xmm2
        movaps nb132nf_ixH2(%esp),%xmm3
        movaps nb132nf_iyH2(%esp),%xmm4
        movaps nb132nf_izH2(%esp),%xmm5
        subps  nb132nf_jxO(%esp),%xmm0
        subps  nb132nf_jyO(%esp),%xmm1
        subps  nb132nf_jzO(%esp),%xmm2
        subps  nb132nf_jxH1(%esp),%xmm3
        subps  nb132nf_jyH1(%esp),%xmm4
        subps  nb132nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb132nf_rsqH2O(%esp)
        movaps %xmm4,nb132nf_rsqH2H1(%esp)

        movaps nb132nf_ixH2(%esp),%xmm0
        movaps nb132nf_iyH2(%esp),%xmm1
        movaps nb132nf_izH2(%esp),%xmm2
        subps  nb132nf_jxH2(%esp),%xmm0
        subps  nb132nf_jyH2(%esp),%xmm1
        subps  nb132nf_jzH2(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb132nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb132nf_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb132nf_rinvH2H2(%esp)
        movaps  %xmm7,nb132nf_rinvH2H1(%esp)

        rsqrtps nb132nf_rsqOO(%esp),%xmm1
        rsqrtps nb132nf_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqOO(%esp),%xmm1
        mulps   nb132nf_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%esp),%xmm3
        mulps   nb132nf_half(%esp),%xmm7
        movaps  %xmm3,nb132nf_rinvOO(%esp)
        movaps  %xmm7,nb132nf_rinvOH1(%esp)

        rsqrtps nb132nf_rsqOH2(%esp),%xmm1
        rsqrtps nb132nf_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqOH2(%esp),%xmm1
        mulps   nb132nf_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%esp),%xmm3
        mulps   nb132nf_half(%esp),%xmm7
        movaps  %xmm3,nb132nf_rinvOH2(%esp)
        movaps  %xmm7,nb132nf_rinvH1O(%esp)

        rsqrtps nb132nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb132nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb132nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb132nf_rsqH1H1(%esp),%xmm1
        mulps   nb132nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%esp),%xmm3
        mulps   nb132nf_half(%esp),%xmm7
        movaps  %xmm3,nb132nf_rinvH1H1(%esp)
        movaps  %xmm7,nb132nf_rinvH1H2(%esp)

        rsqrtps nb132nf_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132nf_three(%esp),%xmm3
        mulps   nb132nf_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132nf_half(%esp),%xmm3
        movaps  %xmm3,nb132nf_rinvH2O(%esp)

        ## start with OO interaction - first the table LJ part
        movaps nb132nf_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm1
        mulps  nb132nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulps  nb132nf_tsc(%esp),%xmm1

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

    movl nb132nf_VFtab(%ebp),%esi

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

    movaps nb132nf_c6(%esp),%xmm4
    mulps  %xmm4,%xmm5   ## Vvdw6 

    addps  nb132nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb132nf_Vvdwtot(%esp)

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

    movaps nb132nf_c12(%esp),%xmm4
    mulps  %xmm4,%xmm5 ## Vvdw12 

    addps  nb132nf_Vvdwtot(%esp),%xmm5
    movaps %xmm5,nb132nf_Vvdwtot(%esp)

        ## Coulomb interactions 
        movaps nb132nf_rinvOH1(%esp),%xmm0
        movaps nb132nf_rinvH1H1(%esp),%xmm1
        movaps nb132nf_rinvOO(%esp),%xmm2
        addps  nb132nf_rinvOH2(%esp),%xmm0
        addps  nb132nf_rinvH1H2(%esp),%xmm1
        addps  nb132nf_rinvH1O(%esp),%xmm0
        addps  nb132nf_rinvH2H1(%esp),%xmm1
        addps  nb132nf_rinvH2O(%esp),%xmm0
        addps  nb132nf_rinvH2H2(%esp),%xmm1

        mulps  nb132nf_qqOH(%esp),%xmm0
        mulps  nb132nf_qqHH(%esp),%xmm1
        mulps  nb132nf_qqOO(%esp),%xmm2

        addps  %xmm1,%xmm0
        addps  nb132nf_vctot(%esp),%xmm2
        addps  %xmm2,%xmm0
        movaps %xmm0,nb132nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb132nf_innerk(%esp)
        jl    _nb_kernel132nf_ia32_sse.nb132nf_single_check
        jmp   _nb_kernel132nf_ia32_sse.nb132nf_unroll_loop
_nb_kernel132nf_ia32_sse.nb132nf_single_check: 
        addl $4,nb132nf_innerk(%esp)
        jnz   _nb_kernel132nf_ia32_sse.nb132nf_single_loop
        jmp   _nb_kernel132nf_ia32_sse.nb132nf_updateouterdata
_nb_kernel132nf_ia32_sse.nb132nf_single_loop: 
        movl  nb132nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb132nf_innerjjnr(%esp)

        movl nb132nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb132nf_ixO(%esp),%xmm0
        movaps  nb132nf_iyO(%esp),%xmm1
        movaps  nb132nf_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb132nf_jxO(%esp)
        movaps %xmm4,nb132nf_jyO(%esp)
        movaps %xmm5,nb132nf_jzO(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 
        movaps %xmm0,nb132nf_rsqOO(%esp)

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb132nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb132nf_half(%esp),%xmm3   ## rinv iO - j water in xmm3
        movaps  %xmm3,nb132nf_rinvOO(%esp)


        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb132nf_qqOO(%esp),%xmm4
        movhps  nb132nf_qqOH(%esp),%xmm4

        mulps %xmm4,%xmm3       ## vcoul  
        addps  nb132nf_vctot(%esp),%xmm3
    movaps %xmm3,nb132nf_vctot(%esp)


        movaps nb132nf_rinvOO(%esp),%xmm0
        movss %xmm0,%xmm1
        mulss  nb132nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulss  nb132nf_tsc(%esp),%xmm1

    cvttps2pi %xmm1,%mm6
    cvtpi2ps %mm6,%xmm3
        subss    %xmm3,%xmm1    ## xmm1=eps 
    movss %xmm1,%xmm2
    mulss  %xmm2,%xmm2      ## xmm2=eps2 
    pslld $3,%mm6

    movl nb132nf_VFtab(%ebp),%esi
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

    movss nb132nf_c6(%esp),%xmm4
    mulss  %xmm4,%xmm5   ## Vvdw6 

    addss  nb132nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb132nf_Vvdwtot(%esp)

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

    movss nb132nf_c12(%esp),%xmm4
    mulss  %xmm4,%xmm5 ## Vvdw12 
    addss  nb132nf_Vvdwtot(%esp),%xmm5
    movss %xmm5,nb132nf_Vvdwtot(%esp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb132nf_ixH1(%esp),%xmm0
        movaps  nb132nf_iyH1(%esp),%xmm1
        movaps  nb132nf_izH1(%esp),%xmm2
        movaps  nb132nf_ixH2(%esp),%xmm3
        movaps  nb132nf_iyH2(%esp),%xmm4
        movaps  nb132nf_izH2(%esp),%xmm5
        subps   nb132nf_jxO(%esp),%xmm0
        subps   nb132nf_jyO(%esp),%xmm1
        subps   nb132nf_jzO(%esp),%xmm2
        subps   nb132nf_jxO(%esp),%xmm3
        subps   nb132nf_jyO(%esp),%xmm4
        subps   nb132nf_jzO(%esp),%xmm5
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
        movaps  nb132nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb132nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb132nf_half(%esp),%xmm7   ## rinv H2 - j water  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb132nf_qqOH(%esp),%xmm6
        movhps  nb132nf_qqHH(%esp),%xmm6

        addps   %xmm7,%xmm3
        mulps   %xmm6,%xmm3

        addps   nb132nf_vctot(%esp),%xmm3
        movaps  %xmm3,nb132nf_vctot(%esp)

        decl nb132nf_innerk(%esp)
        jz    _nb_kernel132nf_ia32_sse.nb132nf_updateouterdata
        jmp   _nb_kernel132nf_ia32_sse.nb132nf_single_loop
_nb_kernel132nf_ia32_sse.nb132nf_updateouterdata: 
        ## get n from stack
        movl nb132nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb132nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb132nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb132nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## accumulate total lj energy and update it 
        movaps nb132nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb132nf_Vvdw(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb132nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel132nf_ia32_sse.nb132nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb132nf_n(%esp)
        jmp _nb_kernel132nf_ia32_sse.nb132nf_outer
_nb_kernel132nf_ia32_sse.nb132nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb132nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel132nf_ia32_sse.nb132nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel132nf_ia32_sse.nb132nf_threadloop
_nb_kernel132nf_ia32_sse.nb132nf_end: 
        emms

        movl nb132nf_nouter(%esp),%eax
        movl nb132nf_ninner(%esp),%ebx
        movl nb132nf_outeriter(%ebp),%ecx
        movl nb132nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb132nf_salign(%esp),%eax
        addl %eax,%esp
        addl $808,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


