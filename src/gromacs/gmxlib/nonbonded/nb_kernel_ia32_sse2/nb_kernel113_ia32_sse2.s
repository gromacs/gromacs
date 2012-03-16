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


.globl nb_kernel113_ia32_sse2
.globl _nb_kernel113_ia32_sse2
nb_kernel113_ia32_sse2: 
_nb_kernel113_ia32_sse2:        
.set nb113_p_nri, 8
.set nb113_iinr, 12
.set nb113_jindex, 16
.set nb113_jjnr, 20
.set nb113_shift, 24
.set nb113_shiftvec, 28
.set nb113_fshift, 32
.set nb113_gid, 36
.set nb113_pos, 40
.set nb113_faction, 44
.set nb113_charge, 48
.set nb113_p_facel, 52
.set nb113_argkrf, 56
.set nb113_argcrf, 60
.set nb113_Vc, 64
.set nb113_type, 68
.set nb113_p_ntype, 72
.set nb113_vdwparam, 76
.set nb113_Vvdw, 80
.set nb113_p_tabscale, 84
.set nb113_VFtab, 88
.set nb113_invsqrta, 92
.set nb113_dvda, 96
.set nb113_p_gbtabscale, 100
.set nb113_GBtab, 104
.set nb113_p_nthreads, 108
.set nb113_count, 112
.set nb113_mtx, 116
.set nb113_outeriter, 120
.set nb113_inneriter, 124
.set nb113_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb113_ixO, 0
.set nb113_iyO, 16
.set nb113_izO, 32
.set nb113_ixH1, 48
.set nb113_iyH1, 64
.set nb113_izH1, 80
.set nb113_ixH2, 96
.set nb113_iyH2, 112
.set nb113_izH2, 128
.set nb113_ixM, 144
.set nb113_iyM, 160
.set nb113_izM, 176
.set nb113_iqH, 192
.set nb113_iqM, 208
.set nb113_dxO, 224
.set nb113_dyO, 240
.set nb113_dzO, 256
.set nb113_dxH1, 272
.set nb113_dyH1, 288
.set nb113_dzH1, 304
.set nb113_dxH2, 320
.set nb113_dyH2, 336
.set nb113_dzH2, 352
.set nb113_dxM, 368
.set nb113_dyM, 384
.set nb113_dzM, 400
.set nb113_qqH, 416
.set nb113_qqM, 432
.set nb113_c6, 448
.set nb113_c12, 464
.set nb113_six, 480
.set nb113_twelve, 496
.set nb113_vctot, 512
.set nb113_Vvdwtot, 528
.set nb113_fixO, 544
.set nb113_fiyO, 560
.set nb113_fizO, 576
.set nb113_fixH1, 592
.set nb113_fiyH1, 608
.set nb113_fizH1, 624
.set nb113_fixH2, 640
.set nb113_fiyH2, 656
.set nb113_fizH2, 672
.set nb113_fixM, 688
.set nb113_fiyM, 704
.set nb113_fizM, 720
.set nb113_fjx, 736
.set nb113_fjy, 752
.set nb113_fjz, 768
.set nb113_half, 784
.set nb113_three, 800
.set nb113_two, 816
.set nb113_rinvH1, 832
.set nb113_rinvH2, 848
.set nb113_rinvM, 864
.set nb113_is3, 880
.set nb113_ii3, 884
.set nb113_ntia, 888
.set nb113_innerjjnr, 892
.set nb113_innerk, 896
.set nb113_n, 900
.set nb113_nn1, 904
.set nb113_nri, 908
.set nb113_nouter, 912
.set nb113_ninner, 916
.set nb113_salign, 920
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $924,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb113_salign(%esp)
        emms

        ## Move args passed by reference to stack
        movl nb113_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb113_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb113_nouter(%esp)
        movl %eax,nb113_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx   ## upper half of 0.5
        movl %eax,nb113_half(%esp)
        movl %ebx,nb113_half+4(%esp)
        movsd nb113_half(%esp),%xmm1
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
        movapd %xmm1,nb113_half(%esp)
        movapd %xmm2,nb113_two(%esp)
        movapd %xmm3,nb113_three(%esp)
        movapd %xmm4,nb113_six(%esp)
        movapd %xmm5,nb113_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb113_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb113_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb113_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb113_iqH(%esp)
        movapd %xmm4,nb113_iqM(%esp)

        movl  nb113_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb113_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113_ntia(%esp)
_nb_kernel113_ia32_sse2.nb113_threadloop: 
        movl  nb113_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel113_ia32_sse2.nb113_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113_ia32_sse2.nb113_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113_n(%esp)
        movl %ebx,nb113_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113_ia32_sse2.nb113_outerstart
        jmp _nb_kernel113_ia32_sse2.nb113_end

_nb_kernel113_ia32_sse2.nb113_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113_nouter(%esp),%ebx
        movl %ebx,nb113_nouter(%esp)

_nb_kernel113_ia32_sse2.nb113_outer: 
        movl  nb113_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb113_is3(%esp)      ## store is3 

        movl  nb113_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb113_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb113_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb113_ii3(%esp)

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
        movapd %xmm3,nb113_ixO(%esp)
        movapd %xmm4,nb113_iyO(%esp)
        movapd %xmm5,nb113_izO(%esp)
        movapd %xmm6,nb113_ixH1(%esp)
        movapd %xmm7,nb113_iyH1(%esp)

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
        movapd %xmm6,nb113_izH1(%esp)
        movapd %xmm0,nb113_ixH2(%esp)
        movapd %xmm1,nb113_iyH2(%esp)
        movapd %xmm2,nb113_izH2(%esp)
        movapd %xmm3,nb113_ixM(%esp)
        movapd %xmm4,nb113_iyM(%esp)
        movapd %xmm5,nb113_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb113_vctot(%esp)
        movapd %xmm4,nb113_Vvdwtot(%esp)
        movapd %xmm4,nb113_fixO(%esp)
        movapd %xmm4,nb113_fiyO(%esp)
        movapd %xmm4,nb113_fizO(%esp)
        movapd %xmm4,nb113_fixH1(%esp)
        movapd %xmm4,nb113_fiyH1(%esp)
        movapd %xmm4,nb113_fizH1(%esp)
        movapd %xmm4,nb113_fixH2(%esp)
        movapd %xmm4,nb113_fiyH2(%esp)
        movapd %xmm4,nb113_fizH2(%esp)
        movapd %xmm4,nb113_fixM(%esp)
        movapd %xmm4,nb113_fiyM(%esp)
        movapd %xmm4,nb113_fizM(%esp)

        movl  nb113_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb113_pos(%ebp),%esi
        movl  nb113_faction(%ebp),%edi
        movl  nb113_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb113_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb113_ninner(%esp),%ecx
        movl  %ecx,nb113_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb113_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel113_ia32_sse2.nb113_unroll_loop
        jmp   _nb_kernel113_ia32_sse2.nb113_checksingle
_nb_kernel113_ia32_sse2.nb113_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb113_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb113_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb113_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113_iqM(%esp),%xmm3
        mulpd  nb113_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb113_qqM(%esp)
        movapd  %xmm4,nb113_qqH(%esp)

        movl nb113_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb113_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb113_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113_c6(%esp)
        movapd %xmm6,nb113_c12(%esp)

        movl nb113_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113_ixO(%esp),%xmm4
        movapd nb113_iyO(%esp),%xmm5
        movapd nb113_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxO(%esp)
        movapd %xmm5,nb113_dyO(%esp)
        movapd %xmm6,nb113_dzO(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113_ixH1(%esp),%xmm4
        movapd nb113_iyH1(%esp),%xmm5
        movapd nb113_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxH1(%esp)
        movapd %xmm5,nb113_dyH1(%esp)
        movapd %xmm6,nb113_dzH1(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113_ixH2(%esp),%xmm3
        movapd nb113_iyH2(%esp),%xmm4
        movapd nb113_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb113_dxH2(%esp)
        movapd %xmm4,nb113_dyH2(%esp)
        movapd %xmm5,nb113_dzH2(%esp)
        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb113_iyM(%esp),%xmm3
        movapd nb113_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb113_dxM(%esp)
        movapd %xmm3,nb113_dyM(%esp)
        movapd %xmm4,nb113_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb113_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb113_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,nb113_rinvM(%esp)

        ## do O interactions directly - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb113_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb113_c6(%esp),%xmm1
        mulpd  nb113_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb113_Vvdwtot(%esp),%xmm3
        mulpd  nb113_six(%esp),%xmm1
        mulpd  nb113_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movapd %xmm3,nb113_Vvdwtot(%esp)

        movapd nb113_dxO(%esp),%xmm0
        movapd nb113_dyO(%esp),%xmm1
        movapd nb113_dzO(%esp),%xmm2
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update O forces 
        movapd nb113_fixO(%esp),%xmm3
        movapd nb113_fiyO(%esp),%xmm4
        movapd nb113_fizO(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb113_fixO(%esp)
        movapd %xmm4,nb113_fiyO(%esp)
        movapd %xmm7,nb113_fizO(%esp)
        ## update j forces with water O 
        movapd %xmm0,nb113_fjx(%esp)
        movapd %xmm1,nb113_fjy(%esp)
        movapd %xmm2,nb113_fjz(%esp)

        ## H1 interactions
        movapd  nb113_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulpd  nb113_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulpd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addpd  nb113_vctot(%esp),%xmm6

        movapd nb113_dxH1(%esp),%xmm0
        movapd nb113_dyH1(%esp),%xmm1
        movapd nb113_dzH1(%esp),%xmm2
        movapd %xmm6,nb113_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb113_fixH1(%esp),%xmm3
        movapd nb113_fiyH1(%esp),%xmm4
        movapd nb113_fizH1(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb113_fixH1(%esp)
        movapd %xmm4,nb113_fiyH1(%esp)
        movapd %xmm7,nb113_fizH1(%esp)
        ## update j forces with water H1 
        addpd  nb113_fjx(%esp),%xmm0
        addpd  nb113_fjy(%esp),%xmm1
        addpd  nb113_fjz(%esp),%xmm2
        movapd %xmm0,nb113_fjx(%esp)
        movapd %xmm1,nb113_fjy(%esp)
        movapd %xmm2,nb113_fjz(%esp)

        ## H2 interactions 
        movapd  nb113_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulpd  nb113_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulpd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addpd  nb113_vctot(%esp),%xmm5

        movapd nb113_dxH2(%esp),%xmm0
        movapd nb113_dyH2(%esp),%xmm1
        movapd nb113_dzH2(%esp),%xmm2
        movapd %xmm5,nb113_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb113_fixH2(%esp),%xmm3
        movapd nb113_fiyH2(%esp),%xmm4
        movapd nb113_fizH2(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb113_fixH2(%esp)
        movapd %xmm4,nb113_fiyH2(%esp)
        movapd %xmm7,nb113_fizH2(%esp)
        ## update j forces with water H2
        addpd  nb113_fjx(%esp),%xmm0
        addpd  nb113_fjy(%esp),%xmm1
        addpd  nb113_fjz(%esp),%xmm2
        movapd %xmm0,nb113_fjx(%esp)
        movapd %xmm1,nb113_fjy(%esp)
        movapd %xmm2,nb113_fjz(%esp)

        ## M interactions 
        movapd  nb113_rinvM(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulpd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulpd  nb113_qqM(%esp),%xmm5    ## xmm5=vcoul 
        mulpd  %xmm5,%xmm4              ## total fsM in xmm4 

        addpd  nb113_vctot(%esp),%xmm5

        movapd nb113_dxM(%esp),%xmm0
        movapd nb113_dyM(%esp),%xmm1
        movapd nb113_dzM(%esp),%xmm2
        movapd %xmm5,nb113_vctot(%esp)
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2

        ## update M forces 
        movapd nb113_fixM(%esp),%xmm3
        movapd nb113_fiyM(%esp),%xmm4
        movapd nb113_fizM(%esp),%xmm7
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm7
        movapd %xmm3,nb113_fixM(%esp)
        movapd %xmm4,nb113_fiyM(%esp)
        movapd %xmm7,nb113_fizM(%esp)

        movl nb113_faction(%ebp),%edi
        ## update j forces 
        addpd  nb113_fjx(%esp),%xmm0
        addpd  nb113_fjy(%esp),%xmm1
        addpd  nb113_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb113_innerk(%esp)
        jl   _nb_kernel113_ia32_sse2.nb113_checksingle
        jmp  _nb_kernel113_ia32_sse2.nb113_unroll_loop
_nb_kernel113_ia32_sse2.nb113_checksingle: 
        movl  nb113_innerk(%esp),%edx
        andl  $1,%edx
        jnz  _nb_kernel113_ia32_sse2.nb113_dosingle
        jmp  _nb_kernel113_ia32_sse2.nb113_updateouterdata
_nb_kernel113_ia32_sse2.nb113_dosingle: 
        movl  nb113_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb113_innerjjnr(%esp)

        movl nb113_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113_iqM(%esp),%xmm3
        mulpd  nb113_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb113_qqM(%esp)
        movapd  %xmm4,nb113_qqH(%esp)

        movl nb113_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb113_vdwparam(%ebp),%esi
        shll %eax
        movl nb113_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113_c6(%esp)
        movapd %xmm6,nb113_c12(%esp)

        movl nb113_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113_ixO(%esp),%xmm4
        movapd nb113_iyO(%esp),%xmm5
        movapd nb113_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxO(%esp)
        movapd %xmm5,nb113_dyO(%esp)
        movapd %xmm6,nb113_dzO(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113_ixH1(%esp),%xmm4
        movapd nb113_iyH1(%esp),%xmm5
        movapd nb113_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxH1(%esp)
        movapd %xmm5,nb113_dyH1(%esp)
        movapd %xmm6,nb113_dzH1(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113_ixH2(%esp),%xmm3
        movapd nb113_iyH2(%esp),%xmm4
        movapd nb113_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## store dr 
        movapd %xmm3,nb113_dxH2(%esp)
        movapd %xmm4,nb113_dyH2(%esp)
        movapd %xmm5,nb113_dzH2(%esp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb113_iyM(%esp),%xmm3
        movapd nb113_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## store dr 
        movapd %xmm2,nb113_dxM(%esp)
        movapd %xmm3,nb113_dyM(%esp)
        movapd %xmm4,nb113_dzM(%esp)
        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvH1(%esp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvH2(%esp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvM(%esp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb113_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb113_c6(%esp),%xmm1
        mulsd  nb113_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb113_Vvdwtot(%esp),%xmm3
        mulsd  nb113_six(%esp),%xmm1
        mulsd  nb113_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb113_Vvdwtot(%esp)

        movapd nb113_dxO(%esp),%xmm0
        movapd nb113_dyO(%esp),%xmm1
        movapd nb113_dzO(%esp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb113_fixO(%esp),%xmm3
        movapd nb113_fiyO(%esp),%xmm4
        movapd nb113_fizO(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixO(%esp)
        movsd %xmm4,nb113_fiyO(%esp)
        movsd %xmm7,nb113_fizO(%esp)
        ## update j forces with water O 
        movsd %xmm0,nb113_fjx(%esp)
        movsd %xmm1,nb113_fjy(%esp)
        movsd %xmm2,nb113_fjz(%esp)

        ## H1 interactions
        movapd  nb113_rinvH1(%esp),%xmm6
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb113_qqH(%esp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%esp),%xmm6

        movapd nb113_dxH1(%esp),%xmm0
        movapd nb113_dyH1(%esp),%xmm1
        movapd nb113_dzH1(%esp),%xmm2
        movsd %xmm6,nb113_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb113_fixH1(%esp),%xmm3
        movapd nb113_fiyH1(%esp),%xmm4
        movapd nb113_fizH1(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixH1(%esp)
        movsd %xmm4,nb113_fiyH1(%esp)
        movsd %xmm7,nb113_fizH1(%esp)
        ## update j forces with water H1 
        addsd  nb113_fjx(%esp),%xmm0
        addsd  nb113_fjy(%esp),%xmm1
        addsd  nb113_fjz(%esp),%xmm2
        movsd %xmm0,nb113_fjx(%esp)
        movsd %xmm1,nb113_fjy(%esp)
        movsd %xmm2,nb113_fjz(%esp)

        ## H2 interactions 
        movapd  nb113_rinvH2(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb113_qqH(%esp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%esp),%xmm5

        movapd nb113_dxH2(%esp),%xmm0
        movapd nb113_dyH2(%esp),%xmm1
        movapd nb113_dzH2(%esp),%xmm2
        movsd %xmm5,nb113_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb113_fixH2(%esp),%xmm3
        movapd nb113_fiyH2(%esp),%xmm4
        movapd nb113_fizH2(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixH2(%esp)
        movsd %xmm4,nb113_fiyH2(%esp)
        movsd %xmm7,nb113_fizH2(%esp)
        ## update j forces with water H2 
        addsd  nb113_fjx(%esp),%xmm0
        addsd  nb113_fjy(%esp),%xmm1
        addsd  nb113_fjz(%esp),%xmm2
        movsd %xmm0,nb113_fjx(%esp)
        movsd %xmm1,nb113_fjy(%esp)
        movsd %xmm2,nb113_fjz(%esp)

        ## M interactions 
        movapd  nb113_rinvM(%esp),%xmm5
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb113_qqM(%esp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%esp),%xmm5

        movapd nb113_dxM(%esp),%xmm0
        movapd nb113_dyM(%esp),%xmm1
        movapd nb113_dzM(%esp),%xmm2
        movsd %xmm5,nb113_vctot(%esp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb113_fixM(%esp),%xmm3
        movapd nb113_fiyM(%esp),%xmm4
        movapd nb113_fizM(%esp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixM(%esp)
        movsd %xmm4,nb113_fiyM(%esp)
        movsd %xmm7,nb113_fizM(%esp)

        movl nb113_faction(%ebp),%edi
        ## update j forces 
        addsd  nb113_fjx(%esp),%xmm0
        addsd  nb113_fjy(%esp),%xmm1
        addsd  nb113_fjz(%esp),%xmm2
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel113_ia32_sse2.nb113_updateouterdata: 
        movl  nb113_ii3(%esp),%ecx
        movl  nb113_faction(%ebp),%edi
        movl  nb113_fshift(%ebp),%esi
        movl  nb113_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb113_fixO(%esp),%xmm0
        movapd nb113_fiyO(%esp),%xmm1
        movapd nb113_fizO(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movapd nb113_fixH1(%esp),%xmm0
        movapd nb113_fiyH1(%esp),%xmm1
        movapd nb113_fizH1(%esp),%xmm2

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
        movapd nb113_fixH2(%esp),%xmm0
        movapd nb113_fiyH2(%esp),%xmm1
        movapd nb113_fizH2(%esp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb113_fixM(%esp),%xmm0
        movapd nb113_fiyM(%esp),%xmm1
        movapd nb113_fizM(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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
        movl nb113_n(%esp),%esi
        ## get group index for i particle 
        movl  nb113_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb113_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb113_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb113_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb113_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb113_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113_ia32_sse2.nb113_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113_n(%esp)
        jmp _nb_kernel113_ia32_sse2.nb113_outer
_nb_kernel113_ia32_sse2.nb113_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113_ia32_sse2.nb113_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113_ia32_sse2.nb113_threadloop
_nb_kernel113_ia32_sse2.nb113_end: 
        emms

        movl nb113_nouter(%esp),%eax
        movl nb113_ninner(%esp),%ebx
        movl nb113_outeriter(%ebp),%ecx
        movl nb113_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb113_salign(%esp),%eax
        addl %eax,%esp
        addl $924,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel113nf_ia32_sse2
.globl _nb_kernel113nf_ia32_sse2
nb_kernel113nf_ia32_sse2:       
_nb_kernel113nf_ia32_sse2:      
.set nb113nf_p_nri, 8
.set nb113nf_iinr, 12
.set nb113nf_jindex, 16
.set nb113nf_jjnr, 20
.set nb113nf_shift, 24
.set nb113nf_shiftvec, 28
.set nb113nf_fshift, 32
.set nb113nf_gid, 36
.set nb113nf_pos, 40
.set nb113nf_faction, 44
.set nb113nf_charge, 48
.set nb113nf_p_facel, 52
.set nb113nf_argkrf, 56
.set nb113nf_argcrf, 60
.set nb113nf_Vc, 64
.set nb113nf_type, 68
.set nb113nf_p_ntype, 72
.set nb113nf_vdwparam, 76
.set nb113nf_Vvdw, 80
.set nb113nf_p_tabscale, 84
.set nb113nf_VFtab, 88
.set nb113nf_invsqrta, 92
.set nb113nf_dvda, 96
.set nb113nf_p_gbtabscale, 100
.set nb113nf_GBtab, 104
.set nb113nf_p_nthreads, 108
.set nb113nf_count, 112
.set nb113nf_mtx, 116
.set nb113nf_outeriter, 120
.set nb113nf_inneriter, 124
.set nb113nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb113nf_ixO, 0
.set nb113nf_iyO, 16
.set nb113nf_izO, 32
.set nb113nf_ixH1, 48
.set nb113nf_iyH1, 64
.set nb113nf_izH1, 80
.set nb113nf_ixH2, 96
.set nb113nf_iyH2, 112
.set nb113nf_izH2, 128
.set nb113nf_ixM, 144
.set nb113nf_iyM, 160
.set nb113nf_izM, 176
.set nb113nf_iqH, 192
.set nb113nf_iqM, 208
.set nb113nf_qqH, 224
.set nb113nf_qqM, 240
.set nb113nf_c6, 256
.set nb113nf_c12, 272
.set nb113nf_vctot, 288
.set nb113nf_Vvdwtot, 304
.set nb113nf_half, 320
.set nb113nf_three, 336
.set nb113nf_two, 352
.set nb113nf_is3, 368
.set nb113nf_ii3, 372
.set nb113nf_ntia, 376
.set nb113nf_innerjjnr, 380
.set nb113nf_innerk, 384
.set nb113nf_n, 388
.set nb113nf_nn1, 392
.set nb113nf_nri, 396
.set nb113nf_nouter, 400
.set nb113nf_ninner, 404
.set nb113nf_salign, 408
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $412,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb113nf_salign(%esp)
        emms

        ## Move args passed by reference to stack
        movl nb113nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb113nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb113nf_nouter(%esp)
        movl %eax,nb113nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb113nf_half(%esp)
        movl %ebx,nb113nf_half+4(%esp)
        movsd nb113nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb113nf_half(%esp)
        movapd %xmm2,nb113nf_two(%esp)
        movapd %xmm3,nb113nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb113nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb113nf_charge(%ebp),%edx
        movsd 8(%edx,%ebx,8),%xmm3
        movsd 24(%edx,%ebx,8),%xmm4
        movl nb113nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb113nf_iqH(%esp)
        movapd %xmm4,nb113nf_iqM(%esp)

        movl  nb113nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl nb113nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113nf_ntia(%esp)

_nb_kernel113nf_ia32_sse2.nb113nf_threadloop: 
        movl  nb113nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel113nf_ia32_sse2.nb113nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113nf_ia32_sse2.nb113nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113nf_n(%esp)
        movl %ebx,nb113nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113nf_ia32_sse2.nb113nf_outerstart
        jmp _nb_kernel113nf_ia32_sse2.nb113nf_end

_nb_kernel113nf_ia32_sse2.nb113nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113nf_nouter(%esp),%ebx
        movl %ebx,nb113nf_nouter(%esp)

_nb_kernel113nf_ia32_sse2.nb113nf_outer: 
        movl  nb113nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb113nf_is3(%esp)            ## store is3 

        movl  nb113nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb113nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb113nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb113nf_ii3(%esp)

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
        movapd %xmm3,nb113nf_ixO(%esp)
        movapd %xmm4,nb113nf_iyO(%esp)
        movapd %xmm5,nb113nf_izO(%esp)
        movapd %xmm6,nb113nf_ixH1(%esp)
        movapd %xmm7,nb113nf_iyH1(%esp)

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
        movapd %xmm6,nb113nf_izH1(%esp)
        movapd %xmm0,nb113nf_ixH2(%esp)
        movapd %xmm1,nb113nf_iyH2(%esp)
        movapd %xmm2,nb113nf_izH2(%esp)
        movapd %xmm3,nb113nf_ixM(%esp)
        movapd %xmm4,nb113nf_iyM(%esp)
        movapd %xmm5,nb113nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb113nf_vctot(%esp)
        movapd %xmm4,nb113nf_Vvdwtot(%esp)

        movl  nb113nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb113nf_pos(%ebp),%esi
        movl  nb113nf_faction(%ebp),%edi
        movl  nb113nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb113nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb113nf_ninner(%esp),%ecx
        movl  %ecx,nb113nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb113nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel113nf_ia32_sse2.nb113nf_unroll_loop
        jmp   _nb_kernel113nf_ia32_sse2.nb113nf_checksingle
_nb_kernel113nf_ia32_sse2.nb113nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb113nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb113nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb113nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113nf_iqM(%esp),%xmm3
        mulpd  nb113nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb113nf_qqM(%esp)
        movapd  %xmm4,nb113nf_qqH(%esp)

        movl nb113nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb113nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb113nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113nf_c6(%esp)
        movapd %xmm6,nb113nf_c12(%esp)

        movl nb113nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113nf_ixO(%esp),%xmm4
        movapd nb113nf_iyO(%esp),%xmm5
        movapd nb113nf_izO(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113nf_ixH1(%esp),%xmm4
        movapd nb113nf_iyH1(%esp),%xmm5
        movapd nb113nf_izH1(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113nf_ixH2(%esp),%xmm3
        movapd nb113nf_iyH2(%esp),%xmm4
        movapd nb113nf_izH2(%esp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb113nf_iyM(%esp),%xmm3
        movapd nb113nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm6     ## rinvH1

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm5     ## rinvH2

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm4     ## rinvM

        ## calculate coulomb potentials from rinv.
        addpd   %xmm5,%xmm6     ## rinvH1+rinvH2
        mulpd   nb113nf_qqM(%esp),%xmm4
        mulpd   nb113nf_qqH(%esp),%xmm6
        addpd   %xmm6,%xmm4
        addpd   nb113nf_vctot(%esp),%xmm4
        movapd  %xmm4,nb113nf_vctot(%esp)

        ## do O interactions - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb113nf_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb113nf_c6(%esp),%xmm1
        mulpd  nb113nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb113nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb113nf_Vvdwtot(%esp)
        ## should we do one more iteration? 
        subl $2,nb113nf_innerk(%esp)
        jl   _nb_kernel113nf_ia32_sse2.nb113nf_checksingle
        jmp  _nb_kernel113nf_ia32_sse2.nb113nf_unroll_loop
_nb_kernel113nf_ia32_sse2.nb113nf_checksingle: 
        addl $2,nb113nf_innerk(%esp)
        jnz  _nb_kernel113nf_ia32_sse2.nb113nf_dosingle
        jmp  _nb_kernel113nf_ia32_sse2.nb113nf_updateouterdata
_nb_kernel113nf_ia32_sse2.nb113nf_dosingle: 
        movl  nb113nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb113nf_innerjjnr(%esp)

        movl nb113nf_charge(%ebp),%esi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113nf_iqM(%esp),%xmm3
        mulpd  nb113nf_iqH(%esp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb113nf_qqM(%esp)
        movapd  %xmm4,nb113nf_qqH(%esp)

        movl nb113nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb113nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb113nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113nf_c6(%esp)
        movapd %xmm6,nb113nf_c12(%esp)

        movl nb113nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113nf_ixO(%esp),%xmm4
        movapd nb113nf_iyO(%esp),%xmm5
        movapd nb113nf_izO(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113nf_ixH1(%esp),%xmm4
        movapd nb113nf_iyH1(%esp),%xmm5
        movapd nb113nf_izH1(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113nf_ixH2(%esp),%xmm3
        movapd nb113nf_iyH2(%esp),%xmm4
        movapd nb113nf_izH2(%esp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb113nf_iyM(%esp),%xmm3
        movapd nb113nf_izM(%esp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113nf_ixM(%esp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,%xmm6

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,%xmm5

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%esp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%esp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113nf_three(%esp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,%xmm4

        ## Calculate coulomb potential
        addsd  %xmm5,%xmm6      ## rinvH1+rinvH2
        mulsd  nb113nf_qqM(%esp),%xmm4
        mulsd  nb113nf_qqH(%esp),%xmm6
        addsd  %xmm6,%xmm4
        addsd nb113nf_vctot(%esp),%xmm4
        movsd %xmm4,nb113nf_vctot(%esp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb113nf_two(%esp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb113nf_c6(%esp),%xmm1
        mulsd  nb113nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb113nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb113nf_Vvdwtot(%esp)

_nb_kernel113nf_ia32_sse2.nb113nf_updateouterdata: 
        ## get n from stack
        movl nb113nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb113nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb113nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb113nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb113nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb113nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb113nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113nf_ia32_sse2.nb113nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113nf_n(%esp)
        jmp _nb_kernel113nf_ia32_sse2.nb113nf_outer
_nb_kernel113nf_ia32_sse2.nb113nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113nf_ia32_sse2.nb113nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113nf_ia32_sse2.nb113nf_threadloop
_nb_kernel113nf_ia32_sse2.nb113nf_end: 
        emms

        movl nb113nf_nouter(%esp),%eax
        movl nb113nf_ninner(%esp),%ebx
        movl nb113nf_outeriter(%ebp),%ecx
        movl nb113nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb113nf_salign(%esp),%eax
        addl %eax,%esp
        addl $412,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




